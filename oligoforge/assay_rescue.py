"""Diagnosis and minimally disruptive redesign of existing qPCR assays."""
import time
from . import manual_design as MD
from . import thermo as T

RESCUE_VERSION = "2.2.0"


def _diagnoses(base, observed=None):
    observed = observed or {}
    candidate = base.get("candidate") or {}
    a = candidate.get("assay") or {}
    e = candidate.get("evidence") or {}
    maps = base.get("mappings") or {}
    out = []
    def add(code, severity, evidence, inference):
        out.append(dict(code=code, severity=severity, computational_evidence=evidence,
                        experimental_inference=inference))
    if not e.get("hard_valid", False):
        add("hard_requirement_failure", "critical", "; ".join(e.get("hard_failures") or []),
            "the assay should not be treated as fit for the declared objective until these requirements are resolved")
    if abs(float(a.get("pair_tm_gap") or 0)) > 2.0:
        add("primer_tm_imbalance", "high", "primer Tm gap %.1f C" % a.get("pair_tm_gap"),
            "may reduce balanced amplification; optimize experimentally")
    if e.get("worst_dimer", 0) <= -6.0:
        add("dimer_risk", "high", "worst modeled dimer ΔG %.2f kcal/mol" % e.get("worst_dimer"),
            "may contribute to primer consumption or nonspecific signal")
    if e.get("target_coverage", 0) < 1.0:
        add("target_coverage", "high", "coherent modeled target coverage %.1f%%" % (100 * e.get("target_coverage", 0)),
            "some declared targets may amplify poorly or not at all")
    off = e.get("offtarget") or {}
    if off.get("signal_subjects"):
        add("false_positive_risk", "critical", "%d off-target sequence(s) predict signal" % off["signal_subjects"],
            "redesign before interpreting the assay as exclusive")
    elif off.get("product_subjects"):
        add("offtarget_product", "high", "%d off-target sequence(s) predict a primer product" % off["product_subjects"],
            "a hydrolysis probe may suppress fluorescence, but product competition and assay-purpose requirements still need review")
    for role in ("forward", "reverse"):
        hits = maps.get(role) or []
        if len(hits) > 1:
            add("multiple_%s_placements" % role, "high", "%d viable template placements" % len(hits),
                "multiple primer placements can create unintended products")
        if hits and min(h.get("mismatches", 0) for h in hits) > 0:
            add("%s_target_mismatch" % role, "high", "best mapped placement has %d mismatch(es)" % min(h.get("mismatches", 0) for h in hits),
                "confirm variant frequency and repair 3-prime-proximal mismatches first")
    if a.get("probe") and len(maps.get("probe") or []) > 1:
        add("multiple_probe_placements", "moderate", "%d probe placements" % len(maps.get("probe") or []),
            "repetitive probe binding can complicate signal interpretation")
    if a.get("amplicon", 0) > 180:
        add("long_amplicon", "moderate", "%d bp modeled product" % a["amplicon"],
            "long products can be less efficient, especially for degraded templates")
    robust = (e.get("condition_robustness") or {}).get("valid_fraction")
    if robust is not None and robust < 1.0:
        add("condition_sensitivity", "moderate", "valid in %.0f%% of tested reaction-condition scenarios" % (100 * robust),
            "small changes in salt, Mg, oligo concentration, or annealing temperature may alter performance")
    if observed.get("efficiency") is not None:
        eff = float(observed["efficiency"])
        if eff < 90 or eff > 110:
            add("observed_efficiency", "high", "reported efficiency %.1f%%" % eff,
                "may reflect inhibition, nonspecific amplification, poor standard preparation, or assay design")
    if observed.get("r2") is not None and float(observed["r2"]) < 0.98:
        add("standard_curve_linearity", "high", "reported R² %.3f" % float(observed["r2"]),
            "review dilution accuracy, inhibition, range, and replicate outliers before assigning causality to sequence design")
    if observed.get("melt_peaks") is not None and int(observed["melt_peaks"]) != 1:
        add("melt_abnormality", "high", "%d reported melt peaks" % int(observed["melt_peaks"]),
            "supports nonspecific products or primer dimers but does not identify them")
    if observed.get("probe_signal_problem"):
        add("probe_signal_problem", "high", "low or abnormal probe fluorescence was reported",
            "probe hydrolysis, fluorophore/quencher selection, probe integrity, position, and structure should be tested")
    if observed.get("multiplex_only_failure") or observed.get("multiplex_failure"):
        add("multiplex_only_failure", "high", "failure reported only in multiplex",
            "cross-dimers, reagent competition, target abundance imbalance, or channel compensation may contribute")
    if not out:
        add("no_decisive_computational_fault", "informational", "no modeled hard defect identified",
            "bench failure may arise from chemistry, template quality, inhibitors, synthesis, or unmodeled context")
    return out



def _one_base_repairs(forward, reverse, probe, template, profile, base, objective, targets, offs):
    """Generate target-matched one-base repairs and keep only evidence improvements.

    A substitution is proposed only when an oligo has exactly one mismatch to a
    concrete template placement in the physically correct orientation.  Every
    proposal is then re-run through the full manual analysis, including supplied
    target/off-target panels; no edit is accepted merely because it matches the
    intended template better.
    """
    old = {"forward": forward, "reverse": reverse, "probe": probe}
    expected = {"forward": "+", "reverse": "-", "probe": None}
    base_row = (base or {}).get("candidate") or {}
    base_key = tuple(base_row.get("rank_key") or ())
    rows, seen = [], set()
    for role in ("forward", "reverse", "probe"):
        seq = old.get(role)
        if not seq:
            continue
        # anchor3=False is deliberate: a terminal mismatch is exactly the case a
        # rescue tool must detect rather than hide from its candidate generator.
        hits = MD.map_oligo(seq, template, role, max_mm=1, anchor3=False)
        for hit in hits:
            if hit.get("mismatches") != 1 or hit.get("uncertain"):
                continue
            if expected[role] and hit.get("strand") != expected[role]:
                continue
            site = str(hit.get("site") or "").upper()
            repl = T.revcomp(site) if hit.get("strand") == "-" else site
            if len(repl) != len(seq) or any(b not in "ACGT" for b in repl):
                continue
            diffs = [i for i, (a, b) in enumerate(zip(seq, repl)) if a != b]
            if len(diffs) != 1:
                continue
            ident = (role, repl, hit.get("start"), hit.get("end"))
            if ident in seen:
                continue
            seen.add(ident)
            edited = dict(old); edited[role] = repl
            try:
                analysis = MD.analyze_assay(
                    edited["forward"], edited["reverse"], template, profile,
                    edited.get("probe"), targets=targets, offs=offs,
                    objective=objective, max_mm=2)
            except Exception:
                continue
            cand = analysis.get("candidate") or {}
            new_key = tuple(cand.get("rank_key") or ())
            if base_key and new_key and not (new_key < base_key):
                continue
            pos = diffs[0]
            change = dict(component=role, old=seq, new=repl,
                          nucleotide_position_1based=pos + 1,
                          old_base=seq[pos], new_base=repl[pos],
                          template_placement=[hit.get("start"), hit.get("end")],
                          three_prime_distance=(len(seq) - 1 - pos))
            cand["changes"] = [change]
            cand["components_changed"] = [role]
            cand["components_retained"] = [x for x in ("forward", "reverse", "probe")
                                            if x != role and old.get(x)]
            rows.append(dict(
                disruption_level="one_base_target_match_repair",
                disruption_order=1,
                rationale="Replace the single mismatching base with the concrete intended-template base, then re-evaluate the complete assay.",
                components_retained=cand["components_retained"],
                components_changed=[role], changes=[change], candidate=cand,
                evidence_level="computational hypothesis with complete reranking",
                new_risks=(cand.get("rank_explanation") or {}).get("weaknesses") or [],
                wet_lab_confirmation=["target-variant frequency", "specificity panel",
                                      "efficiency/linearity", "product identity"]))
    rows.sort(key=lambda x: (tuple((x["candidate"] or {}).get("rank_key") or ()),
                             x["components_changed"], x["changes"][0]["nucleotide_position_1based"]))
    return rows

def rescue(forward, reverse, template, profile, probe=None, objective="balanced",
           observed=None, targets=None, offs=None, max_results=4, max_runtime_s=60.0):
    base = MD.analyze_assay(forward, reverse, template, profile, probe,
                            targets=targets, offs=offs, objective=objective)
    diagnoses = _diagnoses(base, observed)
    started = time.monotonic()
    max_runtime_s = max(5.0, min(float(max_runtime_s), 180.0))
    plans = [
        ("small_positional_shift", 2, {}, 3, "Move unlocked oligos no more than three bases."),
        ("replace_forward_primer", 3, {"reverse": True}, None, "Preserve the reverse primer; redesign the forward primer and probe."),
        ("replace_reverse_primer", 3, {"forward": True}, None, "Preserve the forward primer; redesign the reverse primer and probe."),
        ("replace_probe_only", 4, {"primer_pair": True}, None, "Preserve both primers; search for a stronger probe."),
        ("replace_primer_pair_keep_probe", 5, ({"probe": True} if probe else {}), None,
         "When a probe exists, preserve it while replacing the primer pair; otherwise redesign the pair in the same locus."),
        ("new_amplicon_same_target", 6, {}, None, "Select a new amplicon within the supplied target."),
    ]
    redesigns = _one_base_repairs(forward, reverse, probe, template, profile, base,
                                  objective, targets, offs)
    seen = {tuple((x.get("candidate") or {}).get("assay", {}).get(k)
                  for k in ("forward", "reverse", "probe")) for x in redesigns}
    ledgers = []
    stopped_for_budget = False
    for level, order, locks, shift, rationale in plans:
        if time.monotonic() - started >= max_runtime_s:
            stopped_for_budget = True
            break
        rr = MD.constrained_redesign(forward, reverse, template, profile, probe,
                                     locks=locks, objective=objective,
                                     max_results=max_results, max_shift=shift,
                                     targets=targets, offs=offs,
                                     pair_limit=12, probes_per_pair=2,
                                     full_annotation_limit=12,
                                     base_analysis=base)
        ledgers.append(dict(disruption_level=level, ledger=rr.get("search_ledger")))
        level_rows = []
        for row in rr.get("candidates") or []:
            changed = set(row.get("components_changed") or [])
            if level == "replace_probe_only" and changed != {"probe"}:
                continue
            if level == "small_positional_shift" and not changed:
                continue
            ident = tuple((row.get("assay") or {}).get(k) for k in ("forward", "reverse", "probe"))
            if ident in seen:
                continue
            seen.add(ident)
            level_rows.append(dict(disruption_level=level, disruption_order=order,
                                   rationale=rationale,
                                   components_retained=row.get("components_retained"),
                                   components_changed=row.get("components_changed"),
                                   changes=row.get("changes"), candidate=row,
                                   evidence_level="computational hypothesis",
                                   new_risks=(row.get("rank_explanation") or {}).get("weaknesses") or [],
                                   wet_lab_confirmation=["specificity panel", "efficiency/linearity",
                                                         "replicate precision", "product identity"]))
        if level_rows:
            redesigns.extend(level_rows[:max_results])
        if len(redesigns) >= max_results * 3:
            break
    elapsed = round(time.monotonic() - started, 3)
    redesigns.sort(key=lambda x: (int(x.get("disruption_order", 99)),
                                  tuple(((x.get("candidate") or {}).get("rank_key") or ())),
                                  str(x.get("components_changed"))))
    hard_valid = [x for x in redesigns if ((x.get("candidate") or {}).get("evidence") or {}).get("hard_valid")]
    escalation = []
    if not hard_valid:
        escalation.append(dict(disruption_order=7, recommendation="select_new_target_region",
                               reason="No returned redesign was hard-valid for the declared objective; do not force a weak assay within this region."))
    elif any((d.get("code") in {"false_positive_risk", "target_coverage"}) for d in diagnoses):
        escalation.append(dict(disruption_order=7, recommendation="compare_new_target_region",
                               reason="A new region should be benchmarked against the repaired locus because specificity or inclusivity was a primary failure mode."))
    return dict(version=RESCUE_VERSION, base=base, diagnoses=diagnoses,
                redesigns=redesigns[:max_results * 3], escalation_recommendations=escalation,
                search_ledgers=ledgers,
                ranker_manifest=(base.get("ranker_manifest") if isinstance(base, dict) else None),
                runtime=dict(elapsed_s=elapsed, budget_s=max_runtime_s,
                             stopped_for_budget=stopped_for_budget),
                caveat="Computational diagnoses are hypotheses, not causal proof of the observed bench behavior.")
