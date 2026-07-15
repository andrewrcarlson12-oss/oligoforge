#!/usr/bin/env python3
"""Manual mapping, run-comparison, feedback-summary and uncertainty gates."""
import copy
import json
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oligoforge import manual_design as MD
from oligoforge import run_compare as RC
from oligoforge import ranking as R
from oligoforge import ranking_profiles as RP
from oligoforge import ranking_explain as RX
from oligoforge import experimental_feedback as EF
from oligoforge import ranking_benchmark as RB
from oligoforge import profiles as PROFILES

N = F = 0

def check(name, ok, detail=None):
    global N, F
    N += 1
    if ok:
        print("PASS", name)
    else:
        F += 1
        print("FAIL", name, detail if detail is not None else "")

# Manual mapping must expose a 3'-terminal mismatch rather than hiding it.
primer = "ACGTACGT"
template = "TTTACGTACGATTT"
hits = MD.map_oligo(primer, template, "forward", max_mm=1)
terminal = [x for x in hits if x["strand"] == "+" and x["site"] == "ACGTACGA"]
check("manual map reports 3-prime mismatch placement", len(terminal) == 1, hits)
check("3-prime mismatch is not extension eligible",
      terminal and terminal[0]["three_prime_status"] == "mismatch" and not terminal[0]["extension_eligible"], terminal)
check("explicit amplification-only map still filters terminal mismatch",
      not any(x["site"] == "ACGTACGA" and x["strand"] == "+"
              for x in MD.map_oligo(primer, template, "forward", max_mm=1, anchor3=True)))


def ev(cov=1.0, off=0, robust=1.0, trip=1.0, practical=0.0, valid=True):
    return {
        "hard_valid": valid, "hard_failures": [] if valid else ["invalid"],
        "objective": "balanced", "target_coverage": cov, "worst_isolate_3prime": 1.0,
        "probe_mean_identity": 1.0,
        "offtarget": {"signal_subjects": off, "product_subjects": off},
        "condition_robustness": {"valid_fraction": robust},
        "triplet_penalty": trip, "practical_penalty": practical,
        "degeneracy_fold": 1, "panel_risk": 0,
        "evaluations": {"target_epcr": True, "offtarget_epcr": True,
                        "condition_robustness": True, "panel": False, "junction": False},
    }


def cand(f, rank, evidence):
    return {"rank": rank, "assay": {"forward": f, "reverse": "TGCATGCATGCATGCATGCA",
            "probe": "CGTACGTACGTACGTACGTA", "amplicon": 90, "f_xy": [rank * 10, rank * 10 + 20]},
            "evidence": evidence}

obj = RP.get_profile("balanced")
manifest = R.manifest(obj, {"full_annotation_limit": 20}, {"template_sha256": "abc"})
run_a = {"ranker_manifest": manifest, "candidates": [cand("ACGTACGTACGTACGTACGA", 1, ev()), cand("ACGTACGTACGTACGTACGC", 2, ev(trip=2))]}
run_b = copy.deepcopy(run_a)
cmp = RC.compare_runs(run_a, run_b)
check("identical design runs reproduce", cmp["reproducibility"]["state"] == "reproduced", cmp)
check("identical design runs have perfect rank stability",
      cmp["ranking_stability"]["spearman_shared"] == 1.0 and cmp["ranking_stability"]["n_reversed_pairs"] == 0, cmp)

run_reordered = copy.deepcopy(run_a)
run_reordered["candidates"] = list(reversed(run_reordered["candidates"]))
for i, row in enumerate(run_reordered["candidates"], 1): row["rank"] = i
cmp_bad = RC.compare_runs(run_a, run_reordered)
check("identical manifest plus reordered candidates is critical non-reproducibility",
      cmp_bad["reproducibility"]["state"] == "critical_non_reproducibility", cmp_bad["reproducibility"])

obj2 = RP.get_profile("broad_inclusivity")
manifest2 = R.manifest(obj2, {"full_annotation_limit": 20}, {"template_sha256": "abc"})
run_context = copy.deepcopy(run_reordered); run_context["ranker_manifest"] = manifest2
cmp_context = RC.compare_runs(run_a, run_context)
check("declared objective change explains changed context",
      cmp_context["reproducibility"]["state"] == "context_changed" and
      any(x["field"] == "objective" for x in cmp_context["context_differences"]), cmp_context)

# Rank explanations should identify measurable reversal conditions.
item = {"assay": {}, "evidence": ev(cov=.95, off=0, robust=.33, practical=2.0)}
competitor = {"assay": {}, "evidence": ev(cov=1.0, off=1, robust=1.0, practical=.5)}
scenarios = RX.rank_reversal_scenarios(item, competitor)
check("rank reversal analysis names inclusivity tradeoff",
      any("inclusivity" in x["trigger"] for x in scenarios), scenarios)
check("rank reversal analysis names condition robustness tradeoff",
      any("reaction conditions" in x["trigger"] for x in scenarios), scenarios)

# Feedback summaries remain local and only infer preferences in matched contexts.
conditions = {"mv_conc_mM": 50, "dv_conc_mM": 3}
records = [
    {"assay_id": "A", "target_group": "gene1", "status": "success", "conditions": conditions, "efficiency": 98},
    {"assay_id": "A", "target_group": "gene1", "status": "success", "conditions": conditions, "efficiency": 100},
    {"assay_id": "B", "target_group": "gene1", "status": "failed", "conditions": conditions, "efficiency": 70},
    {"assay_id": "B", "target_group": "gene1", "status": "failed", "conditions": conditions, "efficiency": 72},
]
summary = EF.evidence_summary(records)
byid = {x["assay_id"]: x for x in summary["assays"]}
check("feedback summary identifies consistent local success", byid["A"]["evidence_state"] == "consistent_local_success", byid)
check("feedback summary identifies consistent local failure", byid["B"]["evidence_state"] == "consistent_local_failure", byid)
check("feedback summary creates matched-context pairwise preference",
      len(summary["pairwise_preferences"]) == 1 and summary["pairwise_preferences"][0]["preferred_assay_id"] == "A", summary)
check("feedback summary does not silently alter ranker", "does not silently" in summary["ranker_policy"], summary["ranker_policy"])

ci = RB.wilson_interval(11, 11)
check("benchmark 100 percent result carries finite uncertainty", ci["lower"] < 1.0 and ci["upper"] == 1.0, ci)

# Manual edits must be compared through complete evidence, not only sequence/Tm deltas.
original_analyze = MD.analyze_assay
try:
    def fake_analysis(forward, reverse, template, profile, probe=None, targets=None, offs=None,
                      objective="balanced", max_mm=2):
        edited = str(forward).endswith("C")
        evidence = ev(cov=(1.0 if edited else 0.9), off=(0 if edited else 1),
                      robust=(1.0 if edited else 0.33), trip=(0.5 if edited else 2.0),
                      valid=edited)
        if not edited:
            evidence["hard_failures"] = ["no coherent target product"]
        row = {"rank_key": ([0, 0, 0] if edited else [1, 1, 1]),
               "assay": {"forward": forward, "reverse": reverse, "probe": probe,
                         "amplicon": (90 if edited else 0), "f_tm": 61.0, "r_tm": 61.2,
                         "pair_tm_gap": 0.2, "probe_info": ({"tm": 69.0} if probe else None)},
               "evidence": evidence}
        return {"candidate": row, "predicted_products": ([{"size": 90}] if edited else []),
                "ambiguous_mapping": False, "ranker_manifest": {"manifest_hash": "x"}}
    MD.analyze_assay = fake_analysis
    cmp_edit = MD.compare_edits("ACGTACGTACGTACGTACGA", "TGCATGCATGCATGCATGCA",
                                "ACGTACGTACGTACGTACGC", "TGCATGCATGCATGCATGCA",
                                "ACGT" * 50, {}, baseline_probe=None, edited_probe=None)
    check("manual edit comparison prefers complete-evidence improvement",
          cmp_edit["preference"] == "edited_assay_preferred", cmp_edit)
    check("manual edit comparison reports resolved hard failure",
          "no coherent target product" in cmp_edit["hard_failures_resolved"], cmp_edit)
    check("manual edit comparison records nucleotide operation",
          cmp_edit["sequence_edits"] and cmp_edit["sequence_edits"][0]["operations"], cmp_edit)
finally:
    MD.analyze_assay = original_analyze

HMBS = ("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
        "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
actual_same = MD.compare_edits(
    "GAGCTATACCCCGACCTCTG", "CTTCTCTCCAATCTTGGAAAGCG",
    "GAGCTATACCCCGACCTCTG", "CTTCTCTCCAATCTTGGAAAGCG",
    HMBS, PROFILES.PROFILES["sybr_generic"])
check("real manual comparison recalculates identical assay deterministically",
      actual_same["preference"] == "computationally_indistinguishable" and
      actual_same["sequence_edits"] == [], actual_same.get("preference"))

# API contracts for the new tools.
try:
    from fastapi.testclient import TestClient
    import app as APP
    client = TestClient(APP.app)
    r = client.post("/api/design-runs/compare", json={"left": run_a, "right": run_b, "top_k": 10})
    check("run comparison API succeeds", r.status_code == 200 and r.json()["reproducibility"]["state"] == "reproduced", (r.status_code, r.text))
    r = client.post("/api/experimental-feedback/summary", json={"records": records})
    check("feedback summary API succeeds", r.status_code == 200 and r.json()["n_assays"] == 2, (r.status_code, r.text))
    r = client.post("/api/design-runs/compare", json={"left": {"candidates": []}, "right": {"candidates": []}})
    check("empty run comparison remains explicit rather than crashing", r.status_code == 200 and r.json()["winner"]["left"] is None, (r.status_code, r.text))

    # Batch design must use the same structured ranker and preserve its evidence.
    original_design = APP._AD_design.design_from_sequences
    try:
        def fake_design(targets, profile, offs=None, min_ident=0.6, n_candidates=5,
                        objective="balanced", junctions=None, panel=None,
                        search_budget_s=35.0):
            assay = {"forward": "ACGTACGTACGTACGTACGA", "reverse": "TGCATGCATGCATGCATGCA",
                     "probe": "CGTACGTACGTACGTACGTA", "amplicon": 91,
                     "f_tm": 61.2, "r_tm": 61.5, "gblock": "A" * 125,
                     "probe_info": {"tm": 69.8}}
            return {"candidates": [{"rank": 1, "display_score": 94,
                     "assay": assay, "evidence": {"hard_valid": True, "hard_failures": []},
                     "rank_explanation": {"uncertainty": "strong preference",
                                          "strongest_feature": "clean paired specificity",
                                          "weakest_feature": "single target sequence"}}],
                    "ranker_manifest": {"manifest_hash": "abc"},
                    "candidate_attrition": {"stages": []},
                    "n_candidates_screened": 17,
                    "ranking_statement": "structured retained-pool winner"}
        APP._AD_design.design_from_sequences = fake_design
        r = client.post("/api/batch_design", json={"items": [{"name": "sample <1>",
            "template": "ACGT" * 40, "profile": "idt_taqman", "objective": "balanced"}]})
        data = r.json()
        check("batch design uses authoritative structured ranker",
              r.status_code == 200 and data.get("pipeline") == "authoritative_structured_ranker", data)
        check("batch winner preserves ranking evidence and manifest",
              data["results"][0].get("ranker_manifest", {}).get("manifest_hash") == "abc" and
              data["results"][0].get("alternatives_evaluated") == 17 and
              data["results"][0].get("hard_valid") is True, data)
        r = client.post("/api/batch_design", json={"items": [
            {"name": str(i), "template": "ACGT" * 40} for i in range(9)]})
        check("batch request enforces item budget", r.status_code == 422, (r.status_code, r.text))
    finally:
        APP._AD_design.design_from_sequences = original_design

    original_compare = APP.MDS.compare_edits
    try:
        APP.MDS.compare_edits = lambda *args, **kwargs: {
            "preference": "edited_assay_preferred", "improvements": [{"metric": "target_coverage"}],
            "worsenings": [], "hard_failures_introduced": []}
        r = client.post("/api/manual-design/compare-edit", json={
            "baseline_forward": "ACGTACGTACGTACGTACGA",
            "baseline_reverse": "TGCATGCATGCATGCATGCA",
            "edited_forward": "ACGTACGTACGTACGTACGC",
            "edited_reverse": "TGCATGCATGCATGCATGCA",
            "template": "ACGT" * 50, "profile": "idt_taqman", "objective": "balanced"})
        check("manual edit comparison API succeeds",
              r.status_code == 200 and r.json().get("preference") == "edited_assay_preferred",
              (r.status_code, r.text))
    finally:
        APP.MDS.compare_edits = original_compare
except Exception as exc:
    check("decision-analysis API contracts", False, str(exc))

print("DECISION_ANALYSIS: %d passed / %d failed / %d total" % (N-F, F, N))
sys.exit(1 if F else 0)
