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
from . import design as D, conservation as C, ncbi as N, thermo as T, profiles as PROF, specificity as SP, structure as STR


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
        if "P" in disc and disc["P"].get("n"):
            s -= disc["P"]["max_ident"]
        else:
            s -= (disc["F"].get("max_ident", 0) + disc["R"].get("max_ident", 0)) / 2
    return round(s, 1), cons, disc


def design_from_sequences(targets, profile, offs=None, min_ident=0.6, n_candidates=3):
    targets = [t for t in targets if t and len(t) > 60]
    if not targets:
        return dict(error="no usable target sequences (need >=1, ideally several)")
    ref = _reference(targets)
    cands = _candidates(ref, profile, n_candidates)
    if not cands:
        return dict(error="no primer/probe set met this chemistry's Tm window on the reference. "
                          "AT-rich targets like parasite mtDNA can't reach a high-Tm probe — use the "
                          "Auto setting, or pick a low-Tm / MGB profile. (Very short fetched sequences "
                          "can also lack a full amplicon window.)")
    scored = []
    for a in cands:
        sc, cons, disc = _score(a, targets, offs, min_ident)
        scored.append(dict(score=sc, assay=a, conservation=cons, discrimination=disc))
    scored.sort(key=lambda x: -x["score"])
    return dict(n_targets=len(targets), n_offs=len(offs) if offs else 0,
                reference_len=len(ref), n_candidates=len(scored), candidates=scored)


AUTO_ORDER = ["idt_taqman", "parasite_mtdna", "parasite_sybr"]
_PRETTY = {"idt_taqman": "IDT PrimeTime (ZEN double-quenched probe)",
           "idt_affinity": "IDT Affinity Plus (LNA probe)",
           "thermo_taqman": "Thermo TaqMan (MGB)",
           "parasite_mtdna": "low-Tm TaqMan \u2014 order the probe as IDT Affinity Plus (LNA)",
           "parasite_sybr": "low-Tm SYBR (primers only, no probe)",
           "biorad_probe": "Bio-Rad PrimePCR", "sybr_generic": "SYBR generic"}


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


def design_from_query(target_query, profile_key="auto", off_query=None, n_fetch=20,
                      min_ident=0.6, run_blast=False, blast_mode="remote",
                      blast_db="nt", blast_db_path=None, organism=None, prefer_junction=False):
    """Fetch -> design -> (optional) in-silico-PCR the winning pair.

    profile_key may be a specific profile, or "auto": Auto tries the IDT-orderable
    chemistries in order and returns the first that yields a clean assay, so an
    AT-rich parasite target lands on the low-Tm TaqMan instead of just failing."""
    _pairs = N.search_fetch_fasta(target_query, n_fetch)
    tg = [seq for _, seq in _pairs]
    if not tg:
        return dict(error=f"NCBI returned nothing for: {target_query}")
    off = [seq for _, seq in N.search_fetch_fasta(off_query, max(8, n_fetch // 2))] if off_query else None

    if profile_key == "auto":
        out, tried = None, []
        for pk in AUTO_ORDER:
            tried.append(pk)
            r = design_from_sequences(tg, PROF.PROFILES[pk], off, min_ident)
            if r.get("candidates"):
                out = r
                out["profile_used"] = pk
                break
        if out is None:
            return dict(error="Auto could not place a clean assay in any IDT chemistry for this "
                              "target. The region may be too variable, or the fetched sequences too "
                              "short for a full amplicon; try a different gene/region, or pick MGB "
                              "manually.", tried=tried, target_query=target_query, off_query=off_query)
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
    out["off_query"] = off_query

    if out.get("candidates"):
        try:
            _ref = _reference(tg)
        except Exception:
            _ref = None
        _desc = next((a for a, sq in _pairs if sq == _ref), None) if (_ref and _pairs) else None
        out["source_accession"] = _desc.split()[0] if _desc else None
        _annotate(out, _ref, prefer_junction)
    if run_blast and out.get("candidates"):
        a = out["candidates"][0]["assay"]
        try:
            out["specificity"] = SP.in_silico_pcr(a["forward"], a["reverse"], mode=blast_mode,
                                                  db=blast_db, db_path=blast_db_path, organism=organism)
        except Exception as e:
            out["specificity"] = dict(error=f"in-silico PCR could not run: {e}")
    return out
