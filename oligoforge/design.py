"""Candidate enumeration + scoring, mirroring the manual QC workflow:
enumerate windows, gate each on Tm/GC/3'-end/structure, pair on amplicon size
and Tm match, then search the interior for a clean hydrolysis probe.
"""
from . import thermo as T


def _bump(d, key, n=1):
    if d is not None:
        d[key] = d.get(key, 0) + n


def _primer_rejection(seq, c):
    if not (c["len_min"] <= len(seq) <= c["len_max"]): return "length"
    if not (c["gc_min"] <= T.gc_percent(seq) <= c["gc_max"]): return "gc"
    if not (c["tm_min"] <= T.tm(seq) <= c["tm_max"]): return "tm"
    if c.get("no_three_prime_T") and seq[-1] == "T": return "three_prime_T"
    if T.last5_gc(seq) > c["max_3prime_gc"]: return "three_prime_gc"
    if T.max_run(seq, "G") >= c["max_g_run"]: return "g_run"
    if T.max_run(seq) >= c["max_any_run"]: return "homopolymer"
    ac = c.get("anneal_c", T.ANNEAL_C)
    if T.hairpin_full(seq, ac)[1] <= c["hairpin_min"]: return "hairpin_at_anneal"
    if T.self_dimer_full(seq, ac)[1] <= c["self_dimer_min"]: return "self_dimer_at_anneal"
    return None


def _ok_primer(seq, c):
    return _primer_rejection(seq, c) is None


def enumerate_primers(template, c, audit=None):
    template = template.upper()
    fwd, rev = [], []
    n = len(template)
    if audit is not None:
        audit.update(stage="primer_enumeration", template_length=n,
                     forward_windows=0, reverse_windows=0,
                     forward_retained=0, reverse_retained=0,
                     forward_rejections={}, reverse_rejections={})
    for i in range(n):
        for k in range(c["len_min"], c["len_max"] + 1):
            if i + k > n: break
            w = template[i:i + k]
            rf = _primer_rejection(w, c)
            if audit is not None:
                audit["forward_windows"] += 1
                if rf: _bump(audit["forward_rejections"], rf)
                else: audit["forward_retained"] += 1
            if rf is None: fwd.append((i, i + k, w))
            rv = T.revcomp(w)
            rr = _primer_rejection(rv, c)
            if audit is not None:
                audit["reverse_windows"] += 1
                if rr: _bump(audit["reverse_rejections"], rr)
                else: audit["reverse_retained"] += 1
            if rr is None: rev.append((i, i + k, rv))
    if audit is not None:
        audit["entered"] = audit["forward_windows"] + audit["reverse_windows"]
        audit["retained"] = len(fwd) + len(rev)
        audit["rejected"] = audit["entered"] - audit["retained"]
        audit["hard_gate"] = True
        audit["reversible"] = False
    return fwd, rev


_PAIR_CAP = 300   # rank pairs by Tm-fit first, then dimer-test down the list only this far.
                  # Callers consume the best pair(s); on AT-rich templates hundreds of primers
                  # pass _ok_primer, so the naive O(F x R) hetero-dimer pass cost ~20 s/window.

def pair_primers(fwd, rev, c, audit=None):
    """Best primer pairs, ranked by Tm-fit. Reverse-primer Tms are computed ONCE (were
    recomputed for every f x r pair, twice each); size/Tm-gap-valid pairs are scored cheaply,
    then the expensive hetero-dimer check is applied walking the ranked list until _PAIR_CAP
    dimer-clean pairs are collected. Identical output to the prior filter-then-sort for the
    top pairs every caller actually uses; just orders of magnitude fewer dimer calls."""
    amin, amax = c["amp_min"], c["amp_max"]
    gmin, gapmax, topt = c["min_probe_gap"], c["pair_tm_gap_max"], c["tm_opt"]
    ac = c.get("anneal_c", T.ANNEAL_C)
    rtm = [(rs, re, r, T.tm(r)) for (rs, re, r) in rev]
    prelim = []
    if audit is not None:
        audit.update(stage="primer_pairing", entered=len(fwd) * len(rev),
                     size_rejected=0, interior_gap_rejected=0, tm_gap_rejected=0,
                     dimer_rejected=0, dimer_not_evaluated_due_cap=0,
                     hard_gate=True, reversible=False)
    for (fs, fe, f) in fwd:
        ftm = T.tm(f)
        for (rs, re, r, rt) in rtm:
            amp = re - fs
            if not (amin <= amp <= amax):
                _bump(audit, "size_rejected"); continue
            if rs - fe < gmin:
                _bump(audit, "interior_gap_rejected"); continue
            gap = abs(ftm - rt)
            if gap > gapmax:
                _bump(audit, "tm_gap_rejected"); continue
            score = abs((ftm + rt) / 2 - topt) + gap
            prelim.append((score, fs, fe, f, rs, re, r, amp, gap))
    prelim.sort(key=lambda x: x[0])
    pairs = []
    _sd = {}
    def _self(s):
        v = _sd.get(s)
        if v is None:
            v = T.self_dimer(s); _sd[s] = v
        return v
    for (score, fs, fe, f, rs, re, r, amp, gap) in prelim:
        # REJECT cross-dimers that persist at the annealing temperature; a cross-dimer melted at Ta
        # is not a real interference. The reported/ranked `worst` dimer stays on the 37 C basis so
        # the dimer-floor tie-break (and thus every validated pick) is unchanged.
        if T.hetero_dimer_full(f, r, ac)[1] <= c["pair_dimer_min"]:
            _bump(audit, "dimer_rejected"); continue
        worst = min(_self(f), _self(r), T.hetero_dimer(f, r))
        pairs.append(dict(score=score, fstart=fs, fend=fe, f=f,
                          rstart=rs, rend=re, r=r, amp=amp, gap=gap, dimer=round(worst, 2)))
        if len(pairs) >= _PAIR_CAP: break
    # Final pick is dimer-aware, not Tm-fit alone: a penalty that bites only below a safe dimer
    # floor (-5.5 kcal/mol) lets a clean pair overtake a marginally-better-Tm pair that self- or
    # hetero-dimerizes. Pairs at/above the floor keep their Tm-fit order, so existing designs whose
    # worst dimer is milder than -5.5 are unchanged.
    pairs.sort(key=lambda p: p["score"] + 2.0 * max(0.0, -5.5 - p["dimer"]))
    if audit is not None:
        checked = min(len(prelim), len(pairs) + audit.get("dimer_rejected", 0))
        audit["dimer_not_evaluated_due_cap"] = max(0, len(prelim) - checked)
        audit["cheap_gate_survivors"] = len(prelim)
        audit["retained"] = len(pairs)
        audit["rejected"] = audit["entered"] - len(pairs)
        audit["pair_cap"] = _PAIR_CAP
        audit["candidate_truncation_reversible"] = bool(audit["dimer_not_evaluated_due_cap"])
    return pairs


def enumerate_probe_candidates(template, fend, rstart, fseq, rseq, c, limit=12, audit=None):
    """Return multiple fully screened probes for one primer pair.

    Earlier releases returned only the single locally preferred probe.  That
    greedy choice could discard a stronger complete triplet before target
    coverage, specificity, robustness, or multiplex effects were evaluated.
    The caller now receives a bounded ranked list with exact coordinates.
    """
    template = template.upper()
    inner = template[fend:rstart]
    maxp = max(T.tm(fseq), T.tm(rseq))
    ac = c.get("anneal_c", T.ANNEAL_C)
    out = []
    seen = set()
    if audit is not None:
        audit.update(stage="probe_enumeration", entered=0, retained=0, rejected=0,
                     rejections={}, hard_gate=True, reversible=False,
                     candidate_limit=max(1, int(limit)))
    for a in range(len(inner)):
        for k in range(c["probe_len_min"], c["probe_len_max"] + 1):
            b = a + k
            if b > len(inner):
                break
            genomic = inner[a:b]
            for strand, sub in (("+", genomic), ("-", T.revcomp(genomic))):
                _bump(audit, "entered")
                if sub in seen:
                    if audit is not None: _bump(audit["rejections"], "duplicate_sequence")
                    continue
                if "G" in sub[:3]:
                    if audit is not None: _bump(audit["rejections"], "five_prime_G")
                    continue
                if sub.count("C") < sub.count("G"):
                    if audit is not None: _bump(audit["rejections"], "more_G_than_C")
                    continue
                if T.max_run(sub, "G") >= c["max_g_run"]:
                    if audit is not None: _bump(audit["rejections"], "g_run")
                    continue
                if T.max_run(sub) >= c["max_any_run"]:
                    if audit is not None: _bump(audit["rejections"], "homopolymer")
                    continue
                t = T.tm(sub)
                if not (maxp + c["probe_offset_min"] <= t <= maxp + c["probe_offset_max"]):
                    if audit is not None: _bump(audit["rejections"], "tm_offset")
                    continue
                hdg, hdg_an, htm = T.hairpin_full(sub, ac)
                if hdg_an <= c["probe_hairpin_min"]:
                    if audit is not None: _bump(audit["rejections"], "hairpin_at_anneal")
                    continue
                sd37, sdan, sdtm = T.self_dimer_full(sub, ac)
                if sdan <= c["self_dimer_min"]:
                    if audit is not None: _bump(audit["rejections"], "self_dimer_at_anneal")
                    continue
                df37, dfan, dftm = T.hetero_dimer_full(sub, fseq, ac)
                if dfan <= c["pair_dimer_min"]:
                    if audit is not None: _bump(audit["rejections"], "forward_interaction_at_anneal")
                    continue
                dr37, dran, drtm = T.hetero_dimer_full(sub, rseq, ac)
                if dran <= c["pair_dimer_min"]:
                    if audit is not None: _bump(audit["rejections"], "reverse_interaction_at_anneal")
                    continue
                target = min(max(c["probe_offset_min"], 9.0), c["probe_offset_max"])
                offset = t - maxp
                # Preliminary ordering only.  Full triplet ranking occurs after
                # coverage, specificity, and robustness annotation.
                penalty = (abs(offset - target) +
                           0.35 * max(0.0, -3.0 - hdg) +
                           0.25 * max(0.0, -5.5 - min(sd37, df37, dr37)))
                out.append(dict(probe=sub, strand=strand, tm=t, offset=offset,
                                hairpin_dg=hdg, hairpin_dg_anneal=hdg_an, hairpin_tm=htm,
                                self_dimer=sd37, self_dimer_anneal=sdan, self_dimer_tm=sdtm,
                                dimer_f=df37, dimer_f_anneal=dfan, dimer_f_tm=dftm,
                                dimer_r=dr37, dimer_r_anneal=dran, dimer_r_tm=drtm,
                                start=fend + a, end=fend + b,
                                preliminary_penalty=round(penalty, 5)))
                seen.add(sub)
    out.sort(key=lambda x: (x["preliminary_penalty"], -x["hairpin_dg"],
                            x["start"], x["strand"], x["probe"]))
    # Do not let several nearly identical probe windows consume the complete
    # per-pair beam.  Preserve position, strand and Tm-offset alternatives before
    # the full triplet evidence is available.
    from .candidate_retention import retain_probes_diverse
    selected, diversity_ledger = retain_probes_diverse(out, limit=max(1, int(limit)))
    if audit is not None:
        audit["hard_gate_survivors"] = len(out)
        audit["retained"] = len(selected)
        audit["rejected"] = audit["entered"] - len(selected)
        audit["hard_gate_rejected"] = audit["entered"] - len(out)
        audit["truncated_after_hard_gates"] = max(0, len(out) - len(selected))
        audit["candidate_truncation_reversible"] = bool(audit["truncated_after_hard_gates"])
        audit["diversity_retention"] = diversity_ledger
    return selected


def find_probe(template, fend, rstart, fseq, rseq, c):
    """Backward-compatible best local probe.

    The ranking-truth pipeline calls :func:`enumerate_probe_candidates` and
    evaluates several probes per pair.  Legacy single-assay callers retain the
    first preliminary probe for API compatibility.
    """
    rows = enumerate_probe_candidates(template, fend, rstart, fseq, rseq, c, limit=1)
    return rows[0] if rows else None


def assay_from_pair(template, pair, probe=None):
    """Build one coordinate-exact assay record from a primer pair and optional probe."""
    gb, gs, ge = build_gblock(template, pair["fstart"], pair["rend"])
    out = dict(forward=pair["f"], reverse=pair["r"],
               probe=(probe.get("probe") if probe else None),
               amplicon=pair["amp"], pair_tm_gap=pair["gap"], pair_dimer=pair.get("dimer"),
               amplicon_tm=round(T.amplicon_tm(template.upper()[pair["fstart"]:pair["rend"]]), 1),
               f_tm=T.tm(pair["f"]), r_tm=T.tm(pair["r"]), probe_info=probe,
               gblock=gb, gblock_span=(gs, ge),
               f_xy=(pair["fstart"], pair["fend"]), r_xy=(pair["rstart"], pair["rend"]),
               probe_xy=((probe["start"], probe["end"]) if probe else None),
               preliminary_pair_score=round(float(pair.get("score", 0.0)), 5))
    return out


def generate_assay_candidates(template, c, pair_limit=24, probes_per_pair=4, triplet_limit=96):
    """Jointly enumerate a bounded pool of complete primer/probe triplets.

    Returns ``(assays, ledger)``.  The optimum is explicitly *heuristic bounded*:
    all primer pairs that clear hard gates are ranked, the best ``pair_limit``
    receive multiple probe searches, and the resulting triplets are retained for
    downstream full annotation and structured ranking.
    """
    primer_audit, pair_audit = {}, {}
    fwd, rev = enumerate_primers(template, c, audit=primer_audit)
    all_pairs = pair_primers(fwd, rev, c, audit=pair_audit)
    from .candidate_retention import retain_pairs_diverse
    pairs, pair_retention = retain_pairs_diverse(all_pairs, limit=max(1, int(pair_limit)),
                                                  region_size=max(60, int(c.get("amp_max", 150))),
                                                  amplicon_bin=25, per_near=2)
    ledger = {
        "stage": "window_joint_triplet_search",
        "primer_forward": len(fwd), "primer_reverse": len(rev),
        "pairs_after_hard_gates": len(all_pairs),
        "pairs_fully_explored": len(pairs),
        "pairs_truncated": max(0, len(all_pairs) - len(pairs)),
        "probes_per_pair_limit": max(1, int(probes_per_pair)),
        "triplet_limit": max(1, int(triplet_limit)),
        "hard_gate": False, "reversible": True,
        "search_status": "heuristic_bounded",
        "primer_attrition": primer_audit,
        "pair_attrition": pair_audit,
        "pair_diversity_retention": pair_retention,
        "probe_attrition": [],
    }
    assays = []
    if c.get("no_probe"):
        for pair in pairs:
            assays.append(assay_from_pair(template, pair, None))
    else:
        n_probe_candidates = 0
        pairs_without_probe = 0
        for pair in pairs:
            probe_audit = {}
            probes = enumerate_probe_candidates(template, pair["fend"], pair["rstart"],
                                                 pair["f"], pair["r"], c,
                                                 limit=max(1, int(probes_per_pair)), audit=probe_audit)
            probe_audit["pair_identity"] = [pair["f"], pair["r"]]
            ledger["probe_attrition"].append(probe_audit)
            n_probe_candidates += len(probes)
            if not probes:
                pairs_without_probe += 1
            for probe in probes:
                assays.append(assay_from_pair(template, pair, probe))
        ledger["probe_candidates_retained"] = n_probe_candidates
        ledger["pairs_without_probe"] = pairs_without_probe
    for a in assays:
        a["candidate_rank"] = round(_candidate_rank(a, c), 6)
    assays.sort(key=lambda a: (a["candidate_rank"], a.get("amplicon", 10**9),
                               a.get("f_xy", [0])[0], a.get("forward", ""),
                               a.get("reverse", ""), a.get("probe") or ""))
    before = len(assays)
    from .candidate_retention import retain_diverse
    assays, triplet_retention = retain_diverse(
        assays, limit=max(1, int(triplet_limit)),
        region_size=max(60, int(c.get("amp_max", 150))),
        per_region=max(2, min(8, int(triplet_limit) // 3 or 2)), per_near=2)
    ledger["triplets_generated"] = before
    ledger["triplets_retained"] = len(assays)
    ledger["triplets_truncated"] = max(0, before - len(assays))
    ledger["triplet_diversity_retention"] = triplet_retention
    return assays, ledger

def _rank_single_target_assays(assays, template, c, objective=None, search_ledger=None, annotation_limit=12):
    """Apply the authoritative structured ranker to exact-template candidates.

    This helper intentionally imports the ranking modules lazily: ``ranking`` uses
    the public assay-record schema from this module, while candidate construction
    must remain usable without an import cycle.  The direct designer, batch path
    and sequence viewer all call this same helper.
    """
    if not assays:
        return []
    from . import conservation as C
    from . import ranking as RANK
    from . import ranking_explain as REXPLAIN
    from .candidate_retention import retain_diverse
    target = (template or "").upper()
    annotation_limit = max(1, int(annotation_limit))
    annotation_ledger = None
    if len(assays) > annotation_limit:
        assays, annotation_ledger = retain_diverse(
            assays, limit=annotation_limit,
            region_size=max(60, int(c.get("amp_max", 150))),
            per_region=max(2, min(6, annotation_limit // 3 or 2)), per_near=2)
        if search_ledger is not None:
            search_ledger = dict(search_ledger)
            search_ledger["full_annotation_retention"] = annotation_ledger
            search_ledger["candidates_fully_annotated"] = len(assays)
    scored = []
    for a in assays:
        cons = {"F": C.conservation(a["forward"], [target], 0.6),
                "R": C.conservation(a["reverse"], [target], 0.6)}
        if a.get("probe"):
            cons["P"] = C.conservation(a["probe"], [target], 0.6)
        scored.append(dict(score=0.0, score_raw=a.get("candidate_rank"),
                           quality_penalty=0.0, amplicon_penalty=0.0,
                           assay=a, conservation=cons, discrimination=None))
    objective = objective or ("sybr" if c.get("no_probe") else "balanced")
    ranked, _obj = RANK.rank_candidates(scored, [target], [], c,
                                         objective_name=objective)
    for idx, row in enumerate(ranked):
        competitor = ranked[idx + 1] if idx + 1 < len(ranked) else None
        row["rank_explanation"] = REXPLAIN.explain(row, competitor)
        row["assay"]["ranking_evidence"] = row.get("evidence")
        row["assay"]["rank_trace"] = row.get("rank_trace")
        row["assay"]["rank_explanation"] = row.get("rank_explanation")
        if search_ledger is not None:
            row["assay"]["search_ledger"] = search_ledger
    return ranked


def design_assay(template, c, objective=None):
    """Return the best-supported complete assay on one template.

    Legacy releases selected the first Tm-compatible primer pair and then the
    first acceptable probe.  That greedy path could lose a stronger complete
    triplet.  This implementation constructs a bounded multi-pair/multi-probe
    pool, evaluates every retained triplet with the authoritative hard-gate and
    lexicographic ranker, and returns the first hard-valid result.

    The result remains ``None`` when no hard-valid assay survives so existing API
    contracts stay intact.  ``search_ledger`` and ``rank_trace`` make the bounded
    (heuristic, not exact-global) search auditable.
    """
    template = (template or "").upper()
    if not template:
        return None
    if len(template) <= 650:
        assays, ledger = generate_assay_candidates(
            template, c, pair_limit=8, probes_per_pair=(1 if c.get("no_probe") else 3),
            triplet_limit=24)
    else:
        # Broad targets are searched in target-spanning windows so one dense 5' region
        # cannot consume the complete candidate budget.
        from . import candidate_search as CSEARCH
        assays, ledger = CSEARCH.search(
            template, c, window=440, step=145, budget_s=15.0,
            pair_limit=6, probes_per_pair=(1 if c.get("no_probe") else 2),
            triplets_per_window=12, retained_limit=48, max_windows=14)
    ranked = _rank_single_target_assays(assays, template, c, objective=objective,
                                         search_ledger=ledger, annotation_limit=10)
    winner = next((x for x in ranked if (x.get("evidence") or {}).get("hard_valid")), None)
    return winner["assay"] if winner else None


_GBLOCK_FILLER = ("GCATAGCTGACTGATCAGCTAGTCGATCAGTACGATCGTAGCATCGTAGCTAGCT"
                  "GACGTACGATCGATCAGCTAGCATCGATGCATCGATCGATCGTACGTAGCATGCA")  # GC~50%, no long runs


def build_gblock(template, fstart, rend, flank=40, min_len=125):
    """Amplicon + flanking sequence as a synthesizable standard, at least IDT's
    gBlock minimum (125 bp). Flanks first widen within the template; only if the
    template itself can't reach the minimum is a balanced, low-binding filler
    appended (IDT rejects gene-fragment orders under 125 bp)."""
    s = max(0, fstart - flank)
    e = min(len(template), rend + flank)
    while (e - s) < min_len and (s > 0 or e < len(template)):
        if s > 0:
            s -= 1
        if e < len(template):
            e += 1
    block = template[s:e]
    if len(block) < min_len:
        need = min_len - len(block)
        block += (_GBLOCK_FILLER * (need // len(_GBLOCK_FILLER) + 1))[:need]
    return block, s, e


def build_offtarget_gblock(amplicon, off_seqs, flank=50, min_len=125):
    """Worst-case off-target discrimination control. Finds the off-target sequence whose
    region homologous to the assay's amplicon is the CLOSEST match -- the off-target the
    assay is most likely to fail to reject -- and returns it plus flanks as a synthesizable
    gene fragment. Ordering this next to the positive-control gBlock lets the bench test
    whether off-target rejection holds against the MOST similar off-target, not an average
    one (an average off-target makes a leaky block look clean).
    Returns {seq, off_index, amplicon_identity, span} or None."""
    if not amplicon or not off_seqs:
        return None
    L = len(amplicon)
    if L < 20:
        return None
    best = None  # (identity, off_index, oriented_seq, pos)
    for oi, raw in enumerate(off_seqs):
        for s in (raw, T.revcomp(raw)):
            if len(s) < L:
                continue
            for i in range(len(s) - L + 1):
                cutoff = L if best is None else int(L * (1 - best[0] / 100.0))
                mm = 0
                for a, b in zip(amplicon, s[i:i + L]):
                    if a != b:
                        mm += 1
                        if mm > cutoff:
                            break
                else:
                    ident = 100.0 * (L - mm) / L
                    if best is None or ident > best[0]:
                        best = (ident, oi, s, i)
    if best is None:
        return None
    ident, oi, s, pos = best
    a = max(0, pos - flank)
    e = min(len(s), pos + L + flank)
    while (e - a) < min_len and (a > 0 or e < len(s)):
        if a > 0:
            a -= 1
        if e < len(s):
            e += 1
    block = s[a:e]
    if len(block) < min_len:
        need = min_len - len(block)
        block += (_GBLOCK_FILLER * (need // len(_GBLOCK_FILLER) + 1))[:need]
    return dict(seq=block, off_index=oi, amplicon_identity=round(ident, 1), span=(a, e))


def probe_span(template, assay):
    """Exact base coordinates of the probe within the template -- the literal location the probe
    was taken from in design_assay (find_probe scans template[fend:rstart] on both strands).
    Returns [start, end] (so template[start:end] == probe for a +-strand probe, or == revcomp(probe)
    for a --strand probe), or None. No estimation: it is a substring search for the real oligo."""
    pi = assay.get("probe_info"); probe = assay.get("probe")
    if not pi or not probe:
        return None
    template = template.upper()
    needle = probe if pi.get("strand", "+") == "+" else T.revcomp(probe)
    fx, rx = assay["f_xy"][1], assay["r_xy"][0]
    i = template[fx:rx].find(needle)
    if i >= 0:
        return [fx + i, fx + i + len(needle)]
    j = template.find(needle)
    return [j, j + len(needle)] if j >= 0 else None


def _spread_starts(starts):
    """Return coordinates in an endpoint/midpoint order so a time budget samples the
    whole template before filling local gaps."""
    starts = sorted(set(int(x) for x in starts))
    out, seen = [], set()

    def visit(lo, hi):
        if lo > hi:
            return
        mid = (lo + hi) // 2
        for idx in (lo, hi, mid):
            if idx not in seen:
                seen.add(idx)
                out.append(starts[idx])
        visit(lo + 1, mid - 1)
        visit(mid + 1, hi - 1)

    if starts:
        visit(0, len(starts) - 1)
    return out


def _candidate_rank(a, c):
    """Lower is better. Rank complete assays globally rather than returning the
    first acceptable design window."""
    ft = float(a.get("f_tm", 0.0)); rt = float(a.get("r_tm", 0.0))
    score = abs((ft + rt) / 2.0 - float(c.get("tm_opt", (ft + rt) / 2.0)))
    score += 1.5 * abs(ft - rt)
    amp = int(a.get("amplicon") or 0)
    if amp > 150:
        score += 0.20 * (amp - 150)
    pair_dimer = a.get("pair_dimer")
    if isinstance(pair_dimer, (int, float)):
        score += 2.0 * max(0.0, -5.5 - float(pair_dimer))
    pi = a.get("probe_info")
    if pi:
        target = min(max(float(c.get("probe_offset_min", 0.0)), 9.0),
                     float(c.get("probe_offset_max", 10.0)))
        score += 0.35 * abs(float(pi.get("offset", target)) - target)
        score += 1.5 * max(0.0, -3.0 - float(pi.get("hairpin_dg", 0.0)))
    return score


def design_candidates(template, c, n=5, window=400, step=120, budget=12.0, objective=None):
    """Return diverse, structured-ranked candidates across the complete template.

    The old viewer path kept one greedy assay per search window and ranked those
    preliminary winners.  It now shares the same multi-pair/multi-probe candidate
    search and evidence hierarchy as automatic, manual and rescue design.
    """
    from . import candidate_search as CSEARCH
    from . import ranking as RANK
    template = (template or "").upper()
    n = max(1, int(n))
    rows, ledger = CSEARCH.search(
        template, c, window=window, step=step, budget_s=float(budget),
        pair_limit=max(4, min(8, n + 1)),
        probes_per_pair=(1 if c.get("no_probe") else 2),
        triplets_per_window=max(8, min(18, n * 3)),
        retained_limit=max(18, min(60, n * 8)), max_windows=14)
    ranked = _rank_single_target_assays(rows, template, c, objective=objective,
                                         search_ledger=ledger, annotation_limit=max(12, min(24, n * 4)))
    finalists = RANK.select_finalists(ranked, n=n)
    return [x["assay"] for x in finalists]
