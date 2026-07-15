"""Evidence-structured ranking for complete primer/probe assays.

The ranker never permits a favorable Tm or amplicon length to compensate for a
hard failure.  Candidates are ordered lexicographically from explicit evidence,
with Pareto-front information and a subordinate display score for convenience.
"""
from math import prod
from . import thermo as T
from . import specificity as SP
from .ranking_profiles import get_profile, RANKER_VERSION, PROFILE_VERSION
from . import provenance as PROV

RANKING_SCHEMA_VERSION = "1.2.0"


def _fasta(seqs, prefix):
    return "".join(">%s%d\n%s\n" % (prefix, i + 1, s) for i, s in enumerate(seqs or []))


def _coverage(forward, reverse, probe, seqs, amp_min, amp_max, max_mm=2):
    n = len(seqs or [])
    if not n:
        return dict(n=0, product_subjects=0, signal_subjects=0,
                    product_fraction=None, signal_fraction=None, products=[])
    res = SP.in_silico_pcr_offline(forward, reverse, _fasta(seqs, "s"), probe=probe,
                                   max_mm=max_mm, min_product=max(20, int(amp_min)),
                                   max_product=max(int(amp_max), int(amp_min) + 1))
    if res.get("error"):
        return dict(n=n, error=res["error"], product_subjects=0, signal_subjects=0,
                    product_fraction=0.0, signal_fraction=0.0, products=[])
    products = res.get("products") or []
    psub = {p.get("subject") for p in products}
    ssub = {p.get("subject") for p in products if (not probe or p.get("probe_binds"))}
    return dict(n=n, product_subjects=len(psub), signal_subjects=len(ssub),
                product_fraction=round(len(psub) / n, 4), signal_fraction=round(len(ssub) / n, 4),
                products=products[:50], uncertain_3prime_hits=res.get("uncertain_3prime_hits", 0))


def _condition_robustness(assay, chemistry):
    """Evaluate complete oligo validity across a declared reaction envelope.

    Each scenario re-evaluates primer/probe Tm, hairpins, self-dimers, primer
    cross-dimers and primer-probe interactions at that scenario's annealing
    temperature.  Earlier releases perturbed salt concentrations but checked only
    Tm windows, which overstated robustness for structure-sensitive candidates.
    Immutable snapshots are passed directly to the cached thermodynamic functions;
    process-global reaction conditions are never mutated.
    """
    snap = T._snapshot()
    mv, dv, dntp, dna, anneal = snap
    scenarios = [
        (mv, dv, dntp, dna, anneal, "nominal"),
        (max(0.0, mv * 0.8), max(0.0, dv - 0.5),
         min(dntp + 0.2, max(0.0, dv - 0.55)), max(1e-6, dna * 0.75),
         anneal - 2, "low-stability"),
        (mv * 1.2, dv + 0.5, max(0.0, dntp - 0.2), dna * 1.25,
         anneal + 2, "high-stability"),
    ]
    rows, valid = [], 0
    f, r, p = assay.get("forward"), assay.get("reverse"), assay.get("probe")
    for smv, sdv, sdn, sdc, sa, name in scenarios:
        ss = (smv, sdv, sdn, sdc, sa)
        ft, rt = T._tm_at(f, ss), T._tm_at(r, ss)
        fhp = T._hairpin_full_at(f, sa, ss)[1]
        rhp = T._hairpin_full_at(r, sa, ss)[1]
        fsd = T._self_dimer_full_at(f, sa, ss)[1]
        rsd = T._self_dimer_full_at(r, sa, ss)[1]
        frd = T._hetero_dimer_full_at(f, r, sa, ss)[1]
        failures = []
        if not (chemistry["tm_min"] <= ft <= chemistry["tm_max"]): failures.append("forward_tm")
        if not (chemistry["tm_min"] <= rt <= chemistry["tm_max"]): failures.append("reverse_tm")
        if abs(ft - rt) > chemistry["pair_tm_gap_max"]: failures.append("primer_tm_gap")
        if fhp <= chemistry["hairpin_min"]: failures.append("forward_hairpin")
        if rhp <= chemistry["hairpin_min"]: failures.append("reverse_hairpin")
        if fsd <= chemistry["self_dimer_min"]: failures.append("forward_self_dimer")
        if rsd <= chemistry["self_dimer_min"]: failures.append("reverse_self_dimer")
        if frd <= chemistry["pair_dimer_min"]: failures.append("primer_cross_dimer")
        pt = poff = php = psd = pfd = prd = None
        if p:
            pt = T._tm_at(p, ss)
            poff = pt - max(ft, rt)
            php = T._hairpin_full_at(p, sa, ss)[1]
            psd = T._self_dimer_full_at(p, sa, ss)[1]
            pfd = T._hetero_dimer_full_at(p, f, sa, ss)[1]
            prd = T._hetero_dimer_full_at(p, r, sa, ss)[1]
            if not (chemistry["probe_offset_min"] <= poff <= chemistry["probe_offset_max"]): failures.append("probe_tm_offset")
            if php <= chemistry["probe_hairpin_min"]: failures.append("probe_hairpin")
            if psd <= chemistry["self_dimer_min"]: failures.append("probe_self_dimer")
            if pfd <= chemistry["pair_dimer_min"]: failures.append("probe_forward_dimer")
            if prd <= chemistry["pair_dimer_min"]: failures.append("probe_reverse_dimer")
        ok = not failures
        valid += int(ok)
        rows.append(dict(
            name=name, conditions=dict(mv_conc=smv, dv_conc=sdv,
                                       dntp_conc=sdn, dna_conc=sdc,
                                       anneal_c=sa),
            f_tm=round(ft, 2), r_tm=round(rt, 2),
            probe_tm=(round(pt, 2) if pt is not None else None),
            probe_offset=(round(poff, 2) if poff is not None else None),
            structure_at_anneal=dict(
                forward_hairpin_dg=round(fhp, 3), reverse_hairpin_dg=round(rhp, 3),
                forward_self_dimer_dg=round(fsd, 3), reverse_self_dimer_dg=round(rsd, 3),
                primer_cross_dimer_dg=round(frd, 3),
                probe_hairpin_dg=(round(php, 3) if php is not None else None),
                probe_self_dimer_dg=(round(psd, 3) if psd is not None else None),
                probe_forward_dimer_dg=(round(pfd, 3) if pfd is not None else None),
                probe_reverse_dimer_dg=(round(prd, 3) if prd is not None else None)),
            valid=ok, failure_reasons=failures))
    return dict(valid_fraction=round(valid / len(scenarios), 4),
                valid_scenarios=valid, n_scenarios=len(scenarios), scenarios=rows,
                model_scope="Tm and oligo secondary/cross-structure across three declared reaction-condition scenarios")


def degeneracy_fold(seq):
    vals = T._IUPAC_SETS
    try:
        return int(prod(len(vals.get(b, "ACGT")) for b in (seq or "").upper()))
    except Exception:
        return 1



def evaluate_panel_fit(assay, panel, dimer_threshold=-6.0, amp_tm_gap=2.0):
    """Candidate-versus-panel interactions using the same thermodynamic engine.

    Panel records accept either ``oligos:[{name,seq}]`` or forward/reverse/probe
    fields.  Separate assays remain separate even when display names collide.
    """
    panel = panel or []
    candidate = [("F", assay.get("forward")), ("R", assay.get("reverse"))]
    if assay.get("probe"):
        candidate.append(("P", assay.get("probe")))
    cross = []
    for ci, cs in candidate:
        if not cs:
            continue
        for ai, row in enumerate(panel):
            oligos = row.get("oligos") or [
                {"name": "F", "seq": row.get("forward")},
                {"name": "R", "seq": row.get("reverse")},
                {"name": "P", "seq": row.get("probe")},
            ]
            for po in oligos:
                ps = po.get("seq")
                if not ps:
                    continue
                dg = T.hetero_dimer(cs, ps)
                if dg <= dimer_threshold:
                    end_dg = min(T.end_stability(cs, ps), T.end_stability(ps, cs))
                    cross.append(dict(oligo=ci, assay=row.get("name") or "assay_%d" % (ai + 1),
                                      assay_index=ai, oligo_b=po.get("name") or "?",
                                      dg=round(dg, 2), end_dg=round(end_dg, 2),
                                      three_prime=end_dg <= -5.0))
    cross.sort(key=lambda x: (x["dg"], x["assay_index"], x["oligo_b"]))
    melt = []
    if not assay.get("probe") and assay.get("amplicon_tm") is not None:
        for ai, row in enumerate(panel):
            if row.get("amplicon_tm") is None:
                continue
            is_sybr = row.get("sybr")
            if is_sybr is None:
                is_sybr = not bool(row.get("probe") or any(o.get("name") == "P" for o in row.get("oligos", [])))
            if not is_sybr:
                continue
            delta = abs(float(assay["amplicon_tm"]) - float(row["amplicon_tm"]))
            if delta < amp_tm_gap:
                melt.append(dict(assay=row.get("name") or "assay_%d" % (ai + 1),
                                 assay_index=ai, tm=round(float(row["amplicon_tm"]), 1),
                                 delta=round(delta, 1)))
    return dict(n_panel=len(panel), cross=cross[:20], melt=melt,
                worst_dg=(cross[0]["dg"] if cross else None),
                three_prime=any(x["three_prime"] for x in cross))

def build_evidence(scored, targets, offs, chemistry, objective=None, junction=None, panel_fit=None):
    assay = scored["assay"]
    obj = objective or get_profile(no_probe=not assay.get("probe"))
    amp_min, amp_max = chemistry["amp_min"], chemistry["amp_max"]
    bf, br, bp = assay["forward"], assay["reverse"], assay.get("probe")
    base_tgt = _coverage(bf, br, bp, targets, amp_min, amp_max)
    base_off = _coverage(bf, br, bp, offs or [], amp_min, amp_max)
    df = assay.get("forward_deg") or bf
    dr = assay.get("reverse_deg") or br
    # A degenerate alternative may modify only one component.  Unchanged
    # components must remain in the complete triplet; using ``None`` here used to
    # drop an unchanged hydrolysis probe and mislabel primer-only products as signal.
    dp = (assay.get("probe_deg") or bp) if bp else None
    has_deg = (df != bf or dr != br or dp != bp)
    deg_tgt = _coverage(df, dr, dp, targets, amp_min, amp_max) if has_deg else base_tgt
    deg_off = _coverage(df, dr, dp, offs or [], amp_min, amp_max) if has_deg else base_off

    def _signal(x):
        return x["signal_fraction"] if bp else x["product_fraction"]
    use_deg = False
    if has_deg and degeneracy_fold(df) <= int(obj.get("max_degeneracy_fold", 64)) and degeneracy_fold(dr) <= int(obj.get("max_degeneracy_fold", 64)) and degeneracy_fold(dp or "A") <= int(obj.get("max_degeneracy_fold", 64)):
        coverage_gain = (_signal(deg_tgt) or 0.0) > (_signal(base_tgt) or 0.0) + 1e-12
        no_new_signal = deg_off["signal_subjects"] <= base_off["signal_subjects"]
        no_new_product = deg_off["product_subjects"] <= base_off["product_subjects"]
        exclusive_ok = no_new_signal and (not obj.get("require_no_product_offtargets") or no_new_product)
        # Multi-target objective profiles may promote a documented degenerate pool,
        # but only when it improves coherent recovery without worsening a required
        # exclusivity criterion.  The selected order sequences are explicit in evidence.
        use_deg = bool(coverage_gain and exclusive_ok and len(targets or []) > 1)
    ef, er, ep = (df, dr, dp) if use_deg else (bf, br, bp)
    tgt, off = (deg_tgt, deg_off) if use_deg else (base_tgt, base_off)
    cons = scored.get("conservation") or {}
    f3 = (cons.get("F") or {}).get("worst_3prime")
    r3 = (cons.get("R") or {}).get("worst_3prime")
    pmean = (cons.get("P") or {}).get("mean_ident") if assay.get("probe") else None
    robust = _condition_robustness(assay, chemistry)
    amp = int(assay.get("amplicon") or 0)
    worst_dimer = min(T.self_dimer(assay["forward"]), T.self_dimer(assay["reverse"]),
                      T.hetero_dimer(assay["forward"], assay["reverse"]),
                      T.hetero_dimer(assay.get("probe") or assay["forward"], assay["forward"]),
                      T.hetero_dimer(assay.get("probe") or assay["reverse"], assay["reverse"]))
    fold = max(degeneracy_fold(ef), degeneracy_fold(er), degeneracy_fold(ep or "A"))

    hard = []
    synth = ef + er + (ep or "")
    if any(b not in T._IUPAC_SETS for b in synth):
        hard.append("invalid synthesized-oligo base")
    if "N" in synth:
        hard.append("unresolved N in synthesized oligo")
    fxy, rxy = assay.get("f_xy"), assay.get("r_xy")
    if not (fxy and rxy and fxy[0] < fxy[1] <= rxy[0] < rxy[1]):
        hard.append("incoherent primer geometry")
    if not (amp_min <= amp <= amp_max):
        hard.append("amplicon outside declared limits")
    cov = tgt["signal_fraction"] if assay.get("probe") else tgt["product_fraction"]
    if cov is None or cov + 1e-12 < float(obj["min_target_coverage"]):
        hard.append("target coverage below objective requirement")
    if assay.get("probe") and float(tgt.get("signal_fraction") or 0.0) + 1e-12 < float(obj["min_probe_coverage"]):
        hard.append("probe-bearing target coverage below objective requirement")
    if offs:
        if obj.get("require_no_signal_offtargets") and off["signal_subjects"]:
            hard.append("signal-generating off-target product predicted")
        if obj.get("require_no_product_offtargets") and off["product_subjects"]:
            hard.append("off-target product predicted")
    if obj.get("require_junction") and not (junction and junction.get("level") == "strong"):
        hard.append("required junction-spanning oligo absent")
    if fold > int(obj.get("max_degeneracy_fold", 64)):
        hard.append("degenerate pool exceeds objective limit")

    # The automatic search applies these chemistry rules as hard gates.  Manual,
    # imported, and experimentally annotated assays must be judged by the exact
    # same rules rather than receiving weaker post-hoc analysis.
    def _oligo_gate(seq, role):
        failures = []
        if not seq:
            return failures
        if role == "primer":
            if not (int(chemistry["len_min"]) <= len(seq) <= int(chemistry["len_max"])):
                failures.append("primer length outside chemistry limits")
            gc = T.gc_percent(seq)
            if not (float(chemistry["gc_min"]) <= gc <= float(chemistry["gc_max"])):
                failures.append("primer GC outside chemistry limits")
            tm = T.tm(seq)
            if not (float(chemistry["tm_min"]) <= tm <= float(chemistry["tm_max"])):
                failures.append("primer Tm outside chemistry limits")
            if chemistry.get("no_three_prime_T") and seq[-1:] == "T":
                failures.append("3-prime terminal T prohibited by chemistry profile")
            if T.last5_gc(seq) > int(chemistry["max_3prime_gc"]):
                failures.append("primer 3-prime GC clamp exceeds chemistry limit")
        else:
            if not (int(chemistry["probe_len_min"]) <= len(seq) <= int(chemistry["probe_len_max"])):
                failures.append("probe length outside chemistry limits")
        if T.max_run(seq, "G") >= int(chemistry["max_g_run"]):
            failures.append("G-run exceeds chemistry limit")
        if T.max_run(seq) >= int(chemistry["max_any_run"]):
            failures.append("homopolymer run exceeds chemistry limit")
        ac = float(chemistry.get("anneal_c", T.ANNEAL_C))
        if T.hairpin_full(seq, ac)[1] <= float(chemistry["probe_hairpin_min"] if role == "probe" else chemistry["hairpin_min"]):
            failures.append("%s hairpin persists at annealing temperature" % role)
        if T.self_dimer_full(seq, ac)[1] <= float(chemistry["self_dimer_min"]):
            failures.append("%s self-dimer persists at annealing temperature" % role)
        return failures

    # Gate the actual reference triplet on the same rules used during automatic
    # enumeration.  For a suggested degenerate pool, additionally test whether
    # its concrete Tm intervals have any feasible overlap; do not pretend that an
    # arbitrary IUPAC resolution is the whole pool.
    af, ar, ap = assay["forward"], assay["reverse"], assay.get("probe")
    hard.extend(_oligo_gate(af, "primer"))
    hard.extend(_oligo_gate(ar, "primer"))
    if ap:
        hard.extend(_oligo_gate(ap, "probe"))
    if abs(T.tm(af) - T.tm(ar)) > float(chemistry["pair_tm_gap_max"]):
        hard.append("primer Tm gap exceeds chemistry limit")
    ac = float(chemistry.get("anneal_c", T.ANNEAL_C))
    if T.hetero_dimer_full(af, ar, ac)[1] <= float(chemistry["pair_dimer_min"]):
        hard.append("primer cross-dimer persists at annealing temperature")
    if ap:
        ptm = T.tm(ap)
        poff = ptm - max(T.tm(af), T.tm(ar))
        if not (float(chemistry["probe_offset_min"]) <= poff <= float(chemistry["probe_offset_max"])):
            hard.append("probe Tm offset outside chemistry limits")
        if T.hetero_dimer_full(ap, af, ac)[1] <= float(chemistry["pair_dimer_min"]):
            hard.append("probe-forward interaction persists at annealing temperature")
        if T.hetero_dimer_full(ap, ar, ac)[1] <= float(chemistry["pair_dimer_min"]):
            hard.append("probe-reverse interaction persists at annealing temperature")
    elif not chemistry.get("no_probe"):
        hard.append("probe required by chemistry profile")

    degenerate_tm = None
    if ef != af or er != ar or ep != ap:
        fr, rr = T.tm_range(ef), T.tm_range(er)
        pr = T.tm_range(ep) if ep else None
        pair_min_gap = max(0.0, fr["min"] - rr["max"], rr["min"] - fr["max"])
        if pair_min_gap > float(chemistry["pair_tm_gap_max"]):
            hard.append("all degenerate primer resolutions exceed the pair Tm-gap limit")
        probe_offset_range = None
        if pr:
            lo = pr["min"] - max(fr["max"], rr["max"])
            hi = pr["max"] - max(fr["min"], rr["min"])
            probe_offset_range = [round(lo, 2), round(hi, 2)]
            if hi < float(chemistry["probe_offset_min"]) or lo > float(chemistry["probe_offset_max"]):
                hard.append("all degenerate probe resolutions miss the required Tm offset")
        degenerate_tm = dict(forward=fr, reverse=rr, probe=pr,
                             minimum_possible_pair_gap=round(pair_min_gap, 2),
                             probe_offset_range=probe_offset_range,
                             note="interval feasibility; target-specific oligo combinations may occupy only part of this range")
    hard = sorted(set(hard))

    panel_risk = 0
    if panel_fit:
        panel_risk = 2 if panel_fit.get("three_prime") else (1 if panel_fit.get("cross") or panel_fit.get("melt") else 0)
    triplet_penalty = (abs(float(assay.get("pair_tm_gap") or 0.0)) +
                       max(0.0, -5.5 - worst_dimer) * 2.0)
    if assay.get("probe_info"):
        pi = assay["probe_info"]
        target_offset = min(max(float(chemistry.get("probe_offset_min", 0)), 9.0),
                            float(chemistry.get("probe_offset_max", 10)))
        triplet_penalty += 0.35 * abs(float(pi.get("offset", target_offset)) - target_offset)
    practical_penalty = abs(amp - int(obj.get("prefer_short_amplicon", 110))) / 25.0 + max(0, fold - 1) / 16.0

    e = dict(
        schema_version=RANKING_SCHEMA_VERSION,
        hard_valid=not hard,
        hard_failures=hard,
        objective=obj["key"],
        target=tgt,
        offtarget=off,
        target_coverage=cov or 0.0,
        worst_isolate_3prime=min(x for x in (f3, r3) if x is not None) if any(x is not None for x in (f3, r3)) else 0.0,
        probe_mean_identity=pmean,
        condition_robustness=robust,
        worst_dimer=round(worst_dimer, 3),
        triplet_penalty=round(triplet_penalty, 4),
        practical_penalty=round(practical_penalty, 4),
        degeneracy_fold=fold,
        panel_risk=panel_risk,
        junction=junction,
        effective_oligos=dict(forward=ef, reverse=er, probe=ep, uses_degeneracy=bool(use_deg),
                               base_target_coverage=_signal(base_tgt), degenerate_target_coverage=_signal(deg_tgt),
                               base_offtarget_products=base_off.get("product_subjects", 0),
                               degenerate_offtarget_products=deg_off.get("product_subjects", 0)),
        degenerate_tm=degenerate_tm,
        preliminary_score=scored.get("score_raw"),
        legacy_score=scored.get("score"),
        evaluations=dict(target_epcr=True, offtarget_epcr=bool(offs), conservation=bool(cons),
                         condition_robustness=True, panel=bool(panel_fit), junction=junction is not None),
    )
    return e


def _component(e, name):
    if name == "offtarget":
        return (e["offtarget"]["signal_subjects"], e["offtarget"]["product_subjects"])
    if name == "coverage":
        return (-e["target_coverage"], -e["worst_isolate_3prime"], -(e["probe_mean_identity"] or 0.0))
    if name == "robustness":
        return (-e["condition_robustness"]["valid_fraction"],)
    if name == "multiplex":
        return (e["panel_risk"],)
    if name == "junction":
        level = ((e.get("junction") or {}).get("level") or "none")
        return ({"strong": 0, "size": 1, "none": 2}.get(level, 2),)
    if name == "triplet":
        return (e["triplet_penalty"], -e["worst_dimer"])
    if name == "practical":
        return (e["practical_penalty"], e["degeneracy_fold"])
    return ()


def rank_key(evidence, objective):
    key = (0 if evidence["hard_valid"] else 1, len(evidence["hard_failures"]))
    for part in objective["priority"]:
        key += tuple(_component(evidence, part))
    return key




def structured_rank_vector(evidence, objective):
    """Return the named evidence vector used by :func:`rank_key`.

    This is deliberately redundant with the tuple key: the tuple is convenient for
    deterministic sorting, while the named vector is required for audit, UI
    explanation and reproducible comparisons.
    """
    return dict(
        hard_valid=bool(evidence.get("hard_valid")),
        hard_failure_count=len(evidence.get("hard_failures") or []),
        hard_failures=list(evidence.get("hard_failures") or []),
        priorities=[dict(name=name, value=list(_component(evidence, name)), direction="lower_is_better")
                    for name in objective.get("priority", ())],
    )


def compare_rank_vectors(better, worse, objective):
    """Identify the first authoritative ordering difference between two candidates."""
    bk, wk = rank_key(better, objective), rank_key(worse, objective)
    labels = ["hard_validity", "hard_failure_count"]
    spans = [1, 1]
    for name in objective.get("priority", ()):
        labels.append(name); spans.append(len(_component(better, name)))
    pos = 0
    for label, width in zip(labels, spans):
        bpart = tuple(bk[pos:pos+width]); wpart = tuple(wk[pos:pos+width]); pos += width
        if bpart != wpart:
            return dict(decisive_component=label, candidate_value=list(bpart),
                        competitor_value=list(wpart), direction="lower_is_better")
    return dict(decisive_component="deterministic_tie_breaker",
                candidate_value=None, competitor_value=None,
                direction="coordinate_then_sequence_lexicographic")

def _pareto_vector(e):
    # Hard validity must dominate all soft trade-offs.  An invalid assay can be
    # useful diagnostically, but it must never occupy the same Pareto tier as a
    # hard-valid assay merely because one soft metric is attractive.
    return (0 if e["hard_valid"] else 1, len(e.get("hard_failures") or []),
            -e["target_coverage"], e["offtarget"]["signal_subjects"],
            e["offtarget"]["product_subjects"], -e["condition_robustness"]["valid_fraction"],
            e["triplet_penalty"], e["practical_penalty"], e["panel_risk"])


def pareto_fronts(items):
    remaining = set(range(len(items)))
    fronts = []
    while remaining:
        front = []
        for i in sorted(remaining):
            vi = _pareto_vector(items[i]["evidence"])
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                vj = _pareto_vector(items[j]["evidence"])
                if all(a <= b for a, b in zip(vj, vi)) and any(a < b for a, b in zip(vj, vi)):
                    dominated = True; break
            if not dominated:
                front.append(i)
        fronts.append(front)
        remaining.difference_update(front)
    for fi, front in enumerate(fronts, 1):
        for i in front:
            items[i]["evidence"]["pareto_front"] = fi
    return fronts


def _display_score(e):
    # Convenience only; rank_key remains authoritative.
    s = 100.0
    s -= 35.0 * min(1.0, e["offtarget"]["signal_subjects"])
    s -= 15.0 * max(0, e["offtarget"]["product_subjects"] - e["offtarget"]["signal_subjects"])
    s -= 30.0 * max(0.0, 1.0 - e["target_coverage"])
    s -= 12.0 * max(0.0, 1.0 - e["condition_robustness"]["valid_fraction"])
    s -= min(15.0, e["triplet_penalty"])
    s -= min(8.0, e["practical_penalty"])
    s -= 8.0 * e["panel_risk"]
    if not e["hard_valid"]:
        s = min(s, 49.0)
    return round(max(0.0, min(100.0, s)), 1)


def rank_candidates(scored_candidates, targets, offs, chemistry, objective_name=None,
                    junction_by_identity=None, panel_fit_by_identity=None, panel=None):
    obj = get_profile(objective_name, no_probe=bool(chemistry.get("no_probe")))
    items = []
    for sc in scored_candidates:
        a = sc["assay"]
        ident = (a.get("forward"), a.get("reverse"), a.get("probe"))
        pf = (panel_fit_by_identity or {}).get(ident)
        if pf is None and panel:
            pf = evaluate_panel_fit(a, panel)
        ev = build_evidence(sc, targets, offs, chemistry, obj,
                            junction=(junction_by_identity or {}).get(ident),
                            panel_fit=pf)
        if pf is not None:
            ev["panel_fit"] = pf
        sc["evidence"] = ev
        sc["rank_key"] = rank_key(ev, obj)
        sc["display_score"] = _display_score(ev)
        items.append(sc)
    pareto_fronts(items)
    items.sort(key=lambda x: (x["rank_key"], x["evidence"].get("pareto_front", 99),
                              (x["assay"].get("f_xy") or [0])[0],
                              x["assay"].get("forward", ""), x["assay"].get("reverse", "")))
    for i, it in enumerate(items, 1):
        it["rank"] = i
        it["ranking_evidence_vector"] = structured_rank_vector(it["evidence"], obj)
        nxt = items[i] if i < len(items) else None
        it["rank_trace"] = dict(
            schema_version="oligoforge-ranking-trace/v1",
            rank=i,
            objective=obj.get("key"),
            priority_order=list(obj.get("priority", ())),
            authoritative_rank_key=list(it["rank_key"]),
            evidence_vector=it["ranking_evidence_vector"],
            pareto_front=it["evidence"].get("pareto_front"),
            comparison_to_next=(compare_rank_vectors(it["evidence"], nxt["evidence"], obj) if nxt else None),
            deterministic_tie_breakers=["forward_start", "forward_sequence", "reverse_sequence", "probe_sequence"],
            optimum_status="heuristic_bounded_retained_pool",
        )

    # Consecutive candidates that cannot be defensibly separated receive an
    # explicit rank band.  The deterministic order is preserved for reproducible
    # display/export, but users are not told that a soft tie-break proves a real
    # biological difference.
    from . import ranking_explain as REXPLAIN
    pair_states = []
    for idx, it in enumerate(items):
        nxt = items[idx + 1] if idx + 1 < len(items) else None
        assessment = REXPLAIN.preference_strength(it, nxt)
        it["rank_uncertainty"] = assessment
        pair_states.append(assessment.get("state"))
    group_id = 1
    start = 0
    while start < len(items):
        end = start
        while end < len(items) - 1 and pair_states[end] in {
            "near-equivalent alternatives", "insufficient evidence to distinguish"
        }:
            end += 1
        ranks = list(range(start + 1, end + 2))
        state = ("rank-indistinguishable band" if len(ranks) > 1 else "single deterministic rank")
        group = dict(group_id=group_id, state=state, ranks=ranks,
                     size=len(ranks), schema_version="oligoforge-rank-band/v1")
        for j in range(start, end + 1):
            items[j]["equivalence_group"] = group
        group_id += 1
        start = end + 1
    return items, obj


def select_finalists(items, n=10):
    """Select truthful category winners, then diverse ranked alternatives.

    A category label always belongs to the actual best candidate for that
    category.  Earlier code skipped an already-selected category winner and
    mislabeled the next unused assay as, for example, ``best_specificity``.
    Remaining display slots prefer new primer pairs and genomic regions before
    additional probes on an already represented pair.
    """
    if not items:
        return []
    n = max(1, int(n))
    selectors = [
        ("recommended_balanced", lambda x: x["rank"]),
        ("best_specificity", lambda x: (x["evidence"]["offtarget"]["signal_subjects"], x["evidence"]["offtarget"]["product_subjects"], x["rank"])),
        ("best_inclusivity", lambda x: (-x["evidence"]["target_coverage"], -x["evidence"]["worst_isolate_3prime"], x["rank"])),
        ("most_condition_robust", lambda x: (-x["evidence"]["condition_robustness"]["valid_fraction"], x["rank"])),
        ("best_short_amplicon", lambda x: (x["assay"].get("amplicon", 10**9), x["rank"])),
        ("minimal_degeneracy", lambda x: (x["evidence"]["degeneracy_fold"], x["rank"])),
        ("best_multiplex_fit", lambda x: (x["evidence"]["panel_risk"], x["rank"])),
    ]
    for x in items:
        x["finalist_categories"] = []
    valid = [x for x in items if x["evidence"]["hard_valid"]]
    if not valid:
        return []

    def ident(x):
        a=x["assay"]
        return (a.get("forward"),a.get("reverse"),a.get("probe"))
    def pair_ident(x):
        a=x["assay"]
        return (a.get("forward"),a.get("reverse"))
    def region(x):
        return int((x["assay"].get("f_xy") or [0])[0]) // 250

    chosen=[]; by_ident={}
    for label, fn in selectors:
        winner=min(valid,key=fn)
        key=ident(winner)
        if key not in by_ident:
            chosen.append(winner);by_ident[key]=winner
        if label not in by_ident[key]["finalist_categories"]:
            by_ident[key]["finalist_categories"].append(label)

    # If category winners exceed the requested display budget, authoritative rank
    # decides which are shown; labels remain truthful for the retained winners.
    if len(chosen) >= n:
        return sorted(chosen,key=lambda x:x["rank"])[:n]

    used=set(by_ident); used_pairs={pair_ident(x) for x in chosen}; used_regions={region(x) for x in chosen}
    # First seek alternate regions with new primer pairs.
    for x in valid:
        if len(chosen) >= n: break
        k=ident(x); pk=pair_ident(x); reg=region(x)
        if k not in used and pk not in used_pairs and reg not in used_regions:
            x["finalist_categories"].append("alternate_region")
            chosen.append(x);used.add(k);used_pairs.add(pk);used_regions.add(reg)
    # Then new primer pairs, even within an already represented region.
    for x in valid:
        if len(chosen) >= n: break
        k=ident(x);pk=pair_ident(x)
        if k not in used and pk not in used_pairs:
            x["finalist_categories"].append("alternate_primer_pair")
            chosen.append(x);used.add(k);used_pairs.add(pk);used_regions.add(region(x))
    # Finally allow distinct probe alternatives on represented primer pairs.
    for x in valid:
        if len(chosen) >= n: break
        k=ident(x)
        if k not in used:
            x["finalist_categories"].append("alternate_probe")
            chosen.append(x);used.add(k)
    return sorted(chosen,key=lambda x:x["rank"])[:n]


def manifest(objective, candidate_limits, input_hashes=None, *, external_databases=None,
             warnings=None, fallbacks=None, constraints=None):
    # Imports remain local to avoid turning ranking.py into a dependency hub.
    from .candidate_search import SEARCH_VERSION
    from .candidate_retention import RETENTION_VERSION
    from .manual_design import MANUAL_DESIGN_VERSION
    from .assay_rescue import RESCUE_VERSION
    return PROV.build_manifest(
        ranker_version=RANKER_VERSION,
        ranking_schema=RANKING_SCHEMA_VERSION,
        profile_version=PROFILE_VERSION,
        objective=objective.get("key"),
        candidate_limits=candidate_limits,
        input_hashes=input_hashes or {},
        external_databases=external_databases,
        warnings=warnings,
        fallbacks=fallbacks,
        constraints=constraints,
        search_version=SEARCH_VERSION,
        retention_version=RETENTION_VERSION,
        manual_design_version=MANUAL_DESIGN_VERSION,
        rescue_version=RESCUE_VERSION,
        random_seed=None,
    )
