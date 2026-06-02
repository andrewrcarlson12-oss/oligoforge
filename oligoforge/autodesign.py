"""One-call assay design from an NCBI target.

Flow: fetch sequences for a target query -> design a primer/probe set on the most
complete sequence -> score every candidate for conservation across the WHOLE fetched
set and discrimination against an off-target set -> rank. This turns the manual
fetch/design/validate loop into one call.

It is a SCREEN, not a final answer. It designs on a representative sequence and
validates the result across the set; it does not replace looking at an alignment,
confirming finals in OligoAnalyzer, or bench validation. A bare organism name also
under-specifies the problem: it needs the gene (so the fetched sequences share a
locus) and, for a detection assay, what to discriminate against.
"""
from . import design as D, conservation as C, ncbi as N, thermo as T, profiles as PROF, specificity as SP, structure as STR, refmarkers as RM

# Mitochondrial / ribosomal barcode words: their GenBank records (e.g. MalAvi cytb
# barcodes) are deposited as partial sequences not linked to a Gene record, so a precise
# [Gene Name] fetch would MISS them — keep these as free-text. Protein-coding genes
# (IFNG, IL4, RPL13...) instead get the precise fetch so "interferon gamma" lands on IFNG,
# not IFNGR1 (interferon-gamma receptor 1).
_MARKER_WORDS = {"cytb", "cob", "cytochrome", "oxidase", "nadh", "cox1", "cox2", "cox3", "coi", "co1", "coii", "coiii",
                 "nad1", "nad5", "nad4", "atp6", "18s", "28s", "16s", "23s", "12s",
                 "its", "its1", "its2", "rrna", "ribosomal", "rbcl", "matk", "trnl",
                 "trnh", "psba", "rdrp", "ssu", "lsu", "d-loop",
                 "spacer", "intergenic", "transcribed", "minicircle", "kinetoplast"}


def _is_marker(gene):
    g = (gene or "").lower()
    toks = set(g.replace("-", " ").split())
    return bool(toks & _MARKER_WORDS) or "control region" in g or "ribosomal" in g


def _split_query(q):
    """Best-effort split of a free-text target into (organism, gene): a leading capitalised
    genus (+ lowercase species epithet) is the organism, the remainder is the gene. Used to
    resolve the official gene symbol and to fill the workbench organism/gene fields separately."""
    toks = (q or "").split()
    if not toks:
        return "", ""
    org, i = [], 0
    if toks[0][:1].isupper() and toks[0].isalpha():
        org.append(toks[0]); i = 1
        if i < len(toks) and toks[i].islower() and toks[i].isalpha() and toks[i] not in _MARKER_WORDS:
            org.append(toks[i]); i += 1
    return " ".join(org), " ".join(toks[i:])


def _reference(seqs):
    """Most complete sequence at the transcript scale: longest at or below ~12 kb,
    which admits real mRNAs (immune-gene transcripts run well past 2 kb) while still
    excluding whole mitogenomes (~16 kb+) that would otherwise get sliced in the
    wrong gene. If everything is larger, fall back to the shortest."""
    capped = [s for s in seqs if len(s) <= 12000]
    return max(capped, key=len) if capped else min(seqs, key=len)


def _design_one(template, profile):
    if profile.get("no_probe"):
        fwd, rev = D.enumerate_primers(template, profile)
        pairs = D.pair_primers(fwd, rev, profile)
        if not pairs:
            return None
        p = pairs[0]
        return dict(forward=p["f"], reverse=p["r"], probe=None, amplicon=p["amp"],
                    amplicon_tm=round(T.amplicon_tm(template.upper()[p["fstart"]:p["rend"]]), 1),
                    f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]), probe_info=None,
                    f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
    return D.design_assay(template, profile)


def _candidates(reference, profile, n=3, window=350, step=100):
    """Slide a design window across the WHOLE reference until n distinct candidates
    are found. Previously only two static windows (5' end + exact middle) were tried,
    so a transcript with a bad 5' end and a bad middle could falsely report
    NO ASSAY FOUND while hundreds of valid interior bases went unexamined."""
    L = len(reference)
    if L <= window:
        starts = [0]
    else:
        step_eff = max(step, (L - window) // 29 + 1)   # <= ~30 design windows even on a ~6 kb mitogenome
        starts = list(range(0, L - window + 1, step_eff))
    cands, seen = [], set()
    for s in starts:
        a = _design_one(reference[s:s + window], profile)
        if a and a["forward"] not in seen:
            seen.add(a["forward"]); cands.append(a)
        if len(cands) >= n:
            break
    return cands


def _score(assay, targets, offs, min_ident):
    oligos = {"F": assay["forward"], "R": assay["reverse"]}
    if assay.get("probe"):
        oligos["P"] = assay["probe"]
    cons = {k: C.conservation(s, targets, min_ident) for k, s in oligos.items()}
    disc = {k: C.discrimination(s, offs) for k, s in oligos.items()} if offs else None
    s = (cons["F"]["worst_3prime"] + cons["R"]["worst_3prime"]) / 2
    if "P" in cons:
        s += cons["P"]["mean_ident"]
    if disc:
        # A primer whose 3' end mismatches the off-target won't extend on it, so a single such
        # primer stops the off-target amplicon -- reward the STRONGEST 3'-end blocker, and
        # separately penalize broad off-target similarity across all oligos.
        ids, blocks = [], []
        for k in ("F", "R"):
            dk = disc.get(k)
            if dk and dk.get("n"):
                ids.append(dk.get("max_ident", 0))
                blocks.append(min(dk.get("min_3prime_mismatch", 0), 3))
        dp = disc.get("P")
        if dp and dp.get("n"):
            ids.append(dp.get("max_ident", 0))
        if ids:
            s -= sum(ids) / len(ids)
        if blocks:
            # 9*sum + 5*max equals the old 14*max for a SINGLE blocker (unchanged), and adds 9 per
            # extra 3'-block: a pair that fails the off-target at BOTH 3' ends is far more robust --
            # both terminal mismatches must independently read through to make a false product.
            s += 9.0 * sum(blocks) + 5.0 * max(blocks)
    # Cross-candidate dimer preference: among windows, demote a candidate whose worst self/hetero
    # dimer is worse than a safe floor (-5.5 kcal/mol) so a comparable dimer-clean candidate wins.
    # Candidates milder than the floor get zero penalty, so existing picks are unchanged.
    _wd = min(T.self_dimer(assay["forward"]), T.self_dimer(assay["reverse"]),
              T.hetero_dimer(assay["forward"], assay["reverse"]))
    s -= 4.0 * max(0.0, -5.5 - _wd)
    return round(s, 1), cons, disc


# IUPAC degenerate lookup (reverse of thermo._IUPAC_SETS): frozenset of bases -> code letter
_DEG_CODE = {frozenset(v): k for k, v in T._IUPAC_SETS.items()}


def _seq_quality(seq):
    """Sequence red flags that primer3's per-oligo filters still let through: long
    mononucleotide runs, dinucleotide repeats (slippage / homodimer), extreme GC.
    Returns (flag_strings, score_penalty)."""
    s = (seq or "").upper(); n = len(s)
    if not n:
        return [], 0.0
    flags, pen = [], 0.0
    gc = 100.0 * sum(c in "GC" for c in s) / n
    run = mx = 1
    for i in range(1, n):
        run = run + 1 if s[i] == s[i - 1] else 1
        mx = max(mx, run)
    if mx >= 6:
        flags.append("%d-base run" % mx); pen += 2.0
    elif mx == 5:
        flags.append("5-base run"); pen += 0.5
    best_di = 0
    for i in range(n - 3):
        a, b = s[i], s[i + 1]
        if a == b:
            continue
        j, units = i, 0
        while j + 1 < n and s[j] == a and s[j + 1] == b:
            units += 1; j += 2
        best_di = max(best_di, units)
    if best_di >= 4:
        flags.append("%dx dinucleotide repeat" % best_di); pen += 2.0
    elif best_di == 3:
        flags.append("3x dinucleotide repeat"); pen += 0.8
    if gc > 72:
        flags.append("GC %.0f%% high" % gc); pen += 1.5
    elif gc < 28:
        flags.append("GC %.0f%% low" % gc); pen += 1.0
    return flags, pen


def _best_window(oligo, seq):
    """Lowest-mismatch ungapped placement of oligo in seq; returns (window, mismatches)."""
    L = len(oligo); best, bestmm = None, L + 1
    for i in range(0, len(seq) - L + 1):
        mm = 0
        for a, b in zip(oligo, seq[i:i + L]):
            if a != b:
                mm += 1
                if mm >= bestmm:
                    break
        if mm < bestmm:
            bestmm, best = mm, seq[i:i + L]
            if mm == 0:
                break
    return best, bestmm


def _degenerate(oligo, targets, max_mm_frac=0.25, min_minor=0.18, min_count=2, cap=40):
    """Place the oligo against every target (both orientations), tally the best-matching
    window per target, and emit IUPAC codes where targets disagree -- this is what lets a
    genus design carry a W/R/Y for pan-genus coverage rather than a single consensus base.
    Returns (degenerate_seq, n_degenerate_positions, n_targets_used)."""
    o = (oligo or "").upper(); L = len(o)
    if L < 8:
        return o, 0, 0
    cols = [{} for _ in range(L)]; used = 0
    for t in targets[:cap]:
        t = t.upper()
        w1, m1 = _best_window(o, t)
        try:
            w2, m2 = _best_window(o, T.revcomp(t))
        except Exception:
            w2, m2 = None, L + 1
        w, mm = (w1, m1) if m1 <= m2 else (w2, m2)
        if w is None or mm > L * max_mm_frac:
            continue
        used += 1
        for i, ch in enumerate(w):
            if ch in "ACGT":
                cols[i][ch] = cols[i].get(ch, 0) + 1
    if used < 2:
        return o, 0, used
    out, ndeg = [], 0
    for i, ch0 in enumerate(o):
        c = cols[i]; tot = sum(c.values())
        if not tot:
            out.append(ch0); continue
        maj = max(c, key=c.get)
        keep = {b for b, k in c.items() if k >= min_count and k / tot >= min_minor}
        keep.add(maj)
        out.append(_DEG_CODE.get(frozenset(keep), maj))
        if len(keep) > 1:
            ndeg += 1
    return "".join(out), ndeg, used


def _disc_candidates(reference, profile, offs, want=6, window=350, screen=300):
    """When an off-target set is given, surface primer pairs that DISCRIMINATE it.
    The normal path keeps only the single Tm-optimal pair per window, so a pair whose
    primer 3' end mismatches the off-target -- and would block its amplicon -- is thrown
    away before scoring even when it sits right there in the enumeration. Here we enumerate
    per window, keep every pair whose forward OR reverse primer has a 3'-end mismatch to the
    off-target, rank by 3'-block strength / low off-target identity / Tm-fit, attach a probe,
    and hand them to the scorer next to the Tm-optimal candidates. Best-effort and bounded;
    only runs when offs are provided, so the no-off-target path is unchanged."""
    ref = (reference or "").upper(); L = len(ref)
    offs = [o for o in (offs or []) if o][:12]
    if not offs or L < 40:
        return []
    if L <= window:
        starts = [0]
    else:
        starts = list(range(0, L - window + 1, max(100, (L - window) // 29 + 1)))
    dcache = {}
    def blk3(seq):
        v = dcache.get(seq)
        if v is None:
            if len(dcache) >= 1000:            # bound total alignment work on long, many-window references
                return (0, 100)
            d = C.discrimination(seq, offs)
            v = ((d.get("min_3prime_mismatch", 0), d.get("max_ident", 100)) if d.get("n") else (0, 100))
            dcache[seq] = v
        return v
    ranked = []
    for s in starts:
        win = ref[s:s + window]
        try:
            fwd, rev = D.enumerate_primers(win, profile)
            pairs = D.pair_primers(fwd, rev, profile)
        except Exception:
            continue
        for p in pairs[:screen]:
            fb, fid = blk3(p["f"]); rb, rid = blk3(p["r"])
            blk = max(fb, rb)
            if blk < 1:                       # neither primer 3'-blocks the off-target -> not a discriminator
                continue
            # rank by TOTAL 3'-block first (a pair that fails the off-target at BOTH ends beats one that
            # leans on a single primer), then strongest single block, then low off-target identity, Tm-fit.
            ranked.append((-(fb + rb), -blk, round((fid + rid) / 2.0, 1), p["score"], s, win, p))
    ranked.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    out, seen = [], set()
    no_probe = profile.get("no_probe")
    for _sum, _b, _id, _sc, s, win, p in ranked:
        if p["f"] in seen:
            continue
        if no_probe:
            assay = dict(forward=p["f"], reverse=p["r"], probe=None, amplicon=p["amp"],
                         amplicon_tm=round(T.amplicon_tm(win[p["fstart"]:p["rend"]]), 1),
                         f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]), probe_info=None,
                         f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
        else:
            probe = D.find_probe(win, p["fend"], p["rstart"], p["f"], p["r"], profile)
            if not probe:
                continue
            gb, gs, ge = D.build_gblock(win, p["fstart"], p["rend"])
            assay = dict(forward=p["f"], reverse=p["r"], probe=probe["probe"], amplicon=p["amp"],
                         pair_tm_gap=p.get("gap"), amplicon_tm=round(T.amplicon_tm(win[p["fstart"]:p["rend"]]), 1),
                         f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]), probe_info=probe,
                         gblock=gb, gblock_span=(gs, ge), f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
        seen.add(p["f"]); out.append(assay)
        if len(out) >= want:
            break
    return out


def design_from_sequences(targets, profile, offs=None, min_ident=0.6, n_candidates=5):
    # defensive: de-gap / uppercase / RNA->DNA each sequence and drop anything unusable, so a
    # pasted alignment or RNA sequence is corrected here rather than reaching the Tm engine dirty.
    _clean = []
    for t in (targets or []):
        c, _n, err = T.clean_seq(t) if isinstance(t, str) and t.strip() else ("", None, "empty")
        if c and not err:
            _clean.append(c)
    targets = [t for t in _clean if len(t) > 60]
    if not targets:
        return dict(error="no usable target sequences (need >=1, ideally several)")
    ref = _reference(targets)
    cands = _candidates(ref, profile, n_candidates)
    if offs:                                  # augment with discrimination-aware pairs (3'-blockers the
        try:                                  # Tm-only path never surfaces); best-effort, must not break design
            _dc = _disc_candidates(ref, profile, offs, want=8)
            _have = {c["forward"] for c in cands}
            for _a in _dc:
                if _a["forward"] not in _have:
                    cands.append(_a); _have.add(_a["forward"])
        except Exception:
            pass
    if not cands:
        return dict(error="no primer/probe set met this chemistry's Tm window on the reference. "
                          "AT-rich targets like parasite mtDNA can't reach a high-Tm probe — use the "
                          "Auto setting, or pick a low-Tm / MGB profile. (Very short fetched sequences "
                          "can also lack a full amplicon window.)")
    multi = len(targets) >= 2
    scored = []
    for a in cands:
        sc, cons, disc = _score(a, targets, offs, min_ident)
        qf, pen = _seq_quality(a["forward"])
        rqf, rpen = _seq_quality(a["reverse"]); qf = list(qf) + list(rqf); pen += rpen
        if a.get("probe"):
            pqf, ppen = _seq_quality(a["probe"]); qf += pqf; pen += ppen
        a["quality_flags"] = qf
        if multi:                         # genus / multi-template: add degenerate coverage
            df, ndf, nu = _degenerate(a["forward"], targets)
            dr, ndr, _ = _degenerate(a["reverse"], targets)
            dp, ndp = (None, 0)
            if a.get("probe"):
                dp, ndp, _ = _degenerate(a["probe"], targets)
            if nu >= 2 and (ndf + ndr + ndp) > 0:
                a["forward_deg"], a["reverse_deg"] = df, dr
                if dp is not None:
                    a["probe_deg"] = dp
                a["n_degenerate"] = ndf + ndr + ndp
                a["deg_targets"] = nu
        scored.append(dict(score=round(sc - pen, 1), score_raw=sc, quality_penalty=round(pen, 1),
                           assay=a, conservation=cons, discrimination=disc))
    scored.sort(key=lambda x: -x["score"])
    out = dict(n_targets=len(targets), n_offs=len(offs) if offs else 0,
               reference_len=len(ref), n_candidates=len(scored),
               n_requested=n_candidates, candidates=scored)
    if len(scored) < n_candidates:
        out["constraint_note"] = ("only %d clean set(s) met the %s constraints on this target — "
            "AT-rich or short targets tighten the design space. Try the Auto setting, a low-Tm/MGB "
            "profile, or a wider amplicon window for more options." % (len(scored), profile.get("name", "selected")))
    return out


AUTO_ORDER = ["idt_taqman", "parasite_mtdna", "gc_rich", "parasite_sybr"]   # normal-composition default
_PRETTY = {"idt_taqman": "IDT PrimeTime (ZEN double-quenched probe)",
           "idt_affinity": "IDT Affinity Plus (LNA probe)",
           "thermo_taqman": "Thermo TaqMan (MGB)",
           "parasite_mtdna": "low-Tm TaqMan \u2014 order the probe as IDT Affinity Plus (LNA)",
           "parasite_sybr": "low-Tm SYBR (primers only, no probe)",
           "gc_rich": "GC-rich high-Tm TaqMan (short primers; consider DMSO / 7-deaza-dGTP)",
           "gc_rich_sybr": "GC-rich high-Tm SYBR (no probe)",
           "biorad_probe": "Bio-Rad PrimePCR", "sybr_generic": "SYBR generic"}


def _auto_order(ref):
    """Composition-aware chemistry order for Auto. AT-rich and GC-rich targets each get their
    tuned low-/high-Tm profile tried FIRST (the same first-class treatment), instead of spending
    the first attempt on a generic profile that cannot meet the Tm/GC window \u2014 and, on a GC-rich
    target, failing every profile. 40-62% GC keeps the standard IDT-first order.
    Returns (ordered_profile_keys, gc_percent)."""
    s = (ref or "").upper()
    gc = 100.0 * sum(c in "GC" for c in s) / max(1, len(s))
    if gc < 40.0:
        return ["parasite_mtdna", "idt_taqman", "gc_rich", "parasite_sybr"], gc
    if gc > 62.0:
        return ["gc_rich", "idt_taqman", "parasite_mtdna", "sybr_generic"], gc
    return list(AUTO_ORDER), gc


def _amplicon_on(seq, fwd, amplen):
    """Amplicon start = where the forward primer sits; length is the DESIGNED amplicon
    length. Deliberately does NOT search for the reverse primer: on a long, AT-rich
    whole-mitogenome reference the reverse primer can match far away and yield a multi-kb
    span, and folding that is O(n^3) (minutes). Returns (start, end) or None."""
    if not seq or not amplen:
        return None
    fs = SP._locate(fwd, seq)
    if fs is None:
        return None
    return (fs, min(len(seq), fs + int(amplen)))


def _annotate(out, ref, prefer_junction):
    """Add template-structure (always, offline) and, when prefer_junction is set,
    exon-junction-spanning to each candidate. Structure is folded on the actual designed
    amplicon (from the reference). For the junction call the exon table and the amplicon
    must share a coordinate frame, so junctions are read from whatever transcript the
    gene_table returns and the amplicon is re-located on THAT transcript -- not assumed to
    match the reference (isoform variants differ in exon structure)."""
    cands = out.get("candidates") or []
    junctions = mrna = None
    if prefer_junction and out.get("source_accession"):
        try:
            junctions, jinfo, used = SP.exon_junctions_mrna("", "", out["source_accession"])
        except Exception:
            junctions, jinfo, used = None, "couldn't read exon annotation for the source accession", None
        out["junction_info"] = jinfo
        if junctions and used:
            try:
                recs = N.fetch_accessions([used], "fasta")
                mrna = str(recs[0].seq) if recs else None
            except Exception:
                mrna = None
    for c in cands:
        a = c["assay"]
        rspan = _amplicon_on(ref, a["forward"], a.get("amplicon")) if ref else None
        c["amp_span"] = list(rspan) if rspan else None
        if STR.available() and ref and rspan:
            amp = ref[rspan[0]:rspan[1]]
            f = STR.fold(amp)
            if f:
                fl, rl = len(a["forward"]), len(a["reverse"])
                pp = SP._locate(a["probe"], amp) if a.get("probe") else None
                if a.get("probe") and pp is None:
                    pp = SP._locate(SP._rc_iupac(a["probe"]), amp)
                c["structure"] = dict(
                    mfe=f["mfe"], mfe_per_nt=f["mfe_per_nt"], dna=f["dna_params"],
                    f_paired=STR.site_paired_fraction(f["paired"], 0, fl),
                    r_paired=STR.site_paired_fraction(f["paired"], len(amp) - rl, len(amp)),
                    p_paired=(STR.site_paired_fraction(f["paired"], pp, pp + len(a["probe"]))
                              if (a.get("probe") and pp is not None) else None))
        if prefer_junction:
            jspan = _amplicon_on(mrna, a["forward"], a.get("amplicon")) if (mrna and junctions) else None
            c["spans_junction"] = (any(jspan[0] < j < jspan[1] for j in junctions)
                                   if (junctions and jspan) else None)
    if prefer_junction and junctions:
        cands.sort(key=lambda c: (0 if c.get("spans_junction") else 1, -c["score"]))
    return out


def design_nested(reference, profile, inner_assay, outer_flank_max=600, min_gap=8,
                  outer_amp_max=1400, topk=40):
    """Fully-nested OUTER pair that flanks a given inner (diagnostic) assay on the same
    reference: outer forward upstream of the inner forward, outer reverse downstream of the
    inner reverse, so a second-round reaction only amplifies a correct first-round product
    (the basis of the haemosporidian MalAvi screen). Inner-first, so the diagnostic assay
    keeps its full scoring; outer primer quality is gated by the SAME profile. Returns
    {"outer": {...}} or None when there isn't enough flanking sequence."""
    if not reference or not inner_assay:
        return None
    span = _amplicon_on(reference, inner_assay["forward"], inner_assay.get("amplicon"))
    if not span:
        return None
    I_fs, I_re = span
    L = len(reference)
    up = reference[max(0, I_fs - outer_flank_max):I_fs]
    dn = reference[I_re:min(L, I_re + outer_flank_max)]
    U0 = max(0, I_fs - outer_flank_max)
    if len(up) < profile["len_min"] or len(dn) < profile["len_min"]:
        return None
    fwd, _ = D.enumerate_primers(up, profile)
    _, rev = D.enumerate_primers(dn, profile)
    topt = profile["tm_opt"]
    F = sorted(((U0 + i, U0 + j, w, T.tm(w)) for (i, j, w) in fwd if (U0 + j) <= I_fs - min_gap),
               key=lambda x: abs(x[3] - topt))[:topk]
    R = sorted(((I_re + a, I_re + b, w, T.tm(w)) for (a, b, w) in rev if (I_re + a) >= I_re + min_gap),
               key=lambda x: abs(x[3] - topt))[:topk]
    inner_amp = I_re - I_fs
    best = None
    for (fs, fe, f, ftm) in F:
        for (rs, re, r, rtm) in R:
            amp = re - fs
            if amp < inner_amp + 2 * min_gap or amp > outer_amp_max:
                continue
            gap = abs(ftm - rtm)
            if gap > profile["pair_tm_gap_max"]:
                continue
            if T.hetero_dimer(f, r) <= profile["pair_dimer_min"]:
                continue
            score = abs((ftm + rtm) / 2 - topt) + gap
            if best is None or score < best["score"]:
                best = dict(score=round(score, 2), forward=f, reverse=r, amplicon=amp,
                            f_tm=round(ftm, 1), r_tm=round(rtm, 1), pair_tm_gap=round(gap, 1),
                            f_outside=I_fs - fe, r_outside=rs - I_re)
    return dict(outer=best) if best else None


def design_from_query(target_query, profile_key="auto", off_query=None, n_fetch=20,
                      min_ident=0.6, run_blast=False, blast_mode="remote",
                      blast_db="nt", blast_db_path=None, organism=None, prefer_junction=False,
                      nested=False):
    """Fetch -> design -> (optional) in-silico-PCR the winning pair.

    profile_key may be a specific profile, or "auto": Auto tries the IDT-orderable
    chemistries in order and returns the first that yields a clean assay, so an
    AT-rich parasite target lands on the low-Tm TaqMan instead of just failing."""
    target_query = (target_query or "").strip()
    if not target_query:
        return dict(error='enter a target: an organism plus a gene or marker (e.g. "Plasmodium cytochrome b"), or paste sequences in the Design tab')
    _org, _gene = _split_query(target_query)
    _resolved, _fetch_q = None, target_query
    if _gene and _org and not _is_marker(_gene) and not RM.gold_standard_query(_org, _gene):
        try:
            _resolved = N.resolve_gene(_gene, _org)
        except Exception:
            _resolved = None
        if _resolved and _resolved.get("found"):
            if _resolved.get("in_requested_organism") and _resolved.get("clean"):
                _fetch_q = f'{_resolved["symbol"]}[Gene Name] AND {_org}[Organism]'
            elif _resolved.get("symbol") and not _resolved.get("in_requested_organism"):
                return dict(error=f"{_resolved['symbol']} ({_resolved.get('description','')}) isn't "
                                  f"annotated in {_org} — NCBI has it in "
                                  f"{_resolved.get('organism','another species')}. Design there, or "
                                  f"paste a {_gene} sequence in the Design tab. (Designing on the loosely "
                                  f"matching records in {_org} would risk a paralog like a receptor.)")
    _pairs = N.search_fetch_fasta(_fetch_q, n_fetch)
    if not _pairs and _fetch_q != target_query:
        _pairs = N.search_fetch_fasta(target_query, n_fetch)      # precise found nothing -> free text
    tg = [seq for _, seq in _pairs]
    if not tg:
        msg = f"NCBI returned nothing for: {target_query}"
        if _resolved and _resolved.get("found") and not _resolved.get("in_requested_organism"):
            msg += (f". {_resolved['symbol']} ({_resolved.get('description','')}) is annotated in "
                    f"{_resolved.get('organism','another species')}, not in {_org}.")
        return dict(error=msg)
    _off_pairs = N.search_fetch_fasta(off_query, max(8, n_fetch // 2)) if off_query else []
    off = [seq for _, seq in _off_pairs] if off_query else None

    if profile_key == "auto":
        out, tried = None, []
        _refseq = _reference([t for t in tg if t and len(t) > 60] or tg)
        order, _autogc = _auto_order(_refseq)
        for pk in order:
            tried.append(pk)
            r = design_from_sequences(tg, PROF.PROFILES[pk], off, min_ident)
            if r.get("candidates"):
                out = r
                out["profile_used"] = pk
                break
        if out is None:
            return dict(error="Auto could not place a clean assay in any chemistry for this target "
                              "(tried %s at %.0f%% GC). The region may be too variable, or the fetched "
                              "sequences too short for a full amplicon; try a different gene/region, a "
                              "wider amplicon, or pick a profile manually." % (", ".join(tried), _autogc),
                              tried=tried, target_query=target_query, off_query=off_query)
        out["auto_gc"] = round(_autogc, 1)
    else:
        prof = PROF.PROFILES.get(profile_key)
        if not prof:
            return dict(error=f"unknown profile: {profile_key}")
        out = design_from_sequences(tg, prof, off, min_ident)
        if out.get("error"):
            out["target_query"] = target_query
            out["off_query"] = off_query
            return out
        out["profile_used"] = profile_key

    out["profile_pretty"] = _PRETTY.get(out.get("profile_used"), out.get("profile_used"))
    out["target_query"] = target_query
    out["resolved_gene"] = (_resolved.get("symbol") if (_resolved and _resolved.get("found")) else _gene) or ""
    out["resolved_organism"] = _org or (_resolved.get("organism") if _resolved else "") or ""
    out["off_query"] = off_query
    def _subjects(pairs):
        rows = []
        for desc, _sq in (pairs or [])[:40]:
            parts = (desc or "").split(None, 1)
            rows.append({"acc": parts[0] if parts else "", "title": parts[1] if len(parts) > 1 else "", "len": len(_sq)})
        return rows
    out["target_subjects"] = _subjects(_pairs)   # actual accessions/species behind the conservation call
    out["off_subjects"] = _subjects(_off_pairs)   # ...and the off-target set behind the discrimination call

    if out.get("candidates"):
        try:
            _ref = _reference(tg)
        except Exception:
            _ref = None
        _desc = next((a for a, sq in _pairs if sq == _ref), None) if (_ref and _pairs) else None
        out["source_accession"] = _desc.split()[0] if _desc else None
        _annotate(out, _ref, prefer_junction)
        if out.get("profile_used") in ("idt_affinity", "parasite_mtdna"):
            for c in out["candidates"]:                       # LNA chemistry -> suggest + positions
                pb = (c.get("assay") or {}).get("probe")
                if pb:
                    try:
                        c["assay"]["probe_lna"] = T.suggest_lna(pb)
                    except Exception:
                        pass
        # discrimination-control gBlocks for the lead candidate: positive control (target
        # amplicon + flanks) and the worst-case off-target, so the block can be tested on the bench.
        try:
            _a0 = out["candidates"][0].get("assay") or {}
            _fi = _ref.find(_a0.get("forward", "")) if _ref else -1
            _rcr = T.revcomp(_a0["reverse"]) if _a0.get("reverse") else ""
            _ri = _ref.find(_rcr) if (_ref and _rcr) else -1
            if _ref and _fi >= 0 and _ri >= 0:
                _amp = _ref[_fi:_ri + len(_a0["reverse"])]
                if not _a0.get("gblock"):
                    _gb, _gs, _ge = D.build_gblock(_ref, _fi, _ri + len(_a0["reverse"]))
                    _a0["gblock"] = _gb
                if off:
                    _og = D.build_offtarget_gblock(_amp, off)
                    if _og:
                        _a0["offtarget_gblock"] = _og["seq"]
                        _a0["offtarget_gblock_identity"] = _og["amplicon_identity"]
        except Exception:
            pass
        if nested:
            _prof = PROF.PROFILES.get(out.get("profile_used")) or {}
            try:
                nz = design_nested(_ref, _prof, out["candidates"][0]["assay"]) if (_ref and _prof) else None
            except Exception:
                nz = None
            if nz and nz.get("outer"):
                _inner = out["candidates"][0]["assay"]
                out["nested"] = dict(outer=nz["outer"], inner_amplicon=_inner["amplicon"],
                                     inner_forward=_inner["forward"], inner_reverse=_inner["reverse"])
            else:
                out["nested"] = dict(note="couldn't place a flanking outer pair — the reference needs "
                                          "roughly 150+ bp of usable sequence beyond the inner amplicon "
                                          "on each side. The inner assay above runs fine on its own.")
    if run_blast and out.get("candidates"):
        a = out["candidates"][0]["assay"]
        try:
            out["specificity"] = SP.in_silico_pcr(a["forward"], a["reverse"], mode=blast_mode,
                                                  db=blast_db, db_path=blast_db_path, organism=organism)
        except Exception as e:
            out["specificity"] = dict(error=f"in-silico PCR could not run: {e}")
    return out
