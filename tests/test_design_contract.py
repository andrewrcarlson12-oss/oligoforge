#!/usr/bin/env python3
"""Focused canonical design-contract and search-policy regression gates."""
from copy import deepcopy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oligoforge import autodesign as AD
from oligoforge import design_contract as DC
from oligoforge import profiles as P
from oligoforge import ranking_profiles as RP


passed = failed = 0


def check(name, condition, detail=None):
    global passed, failed
    if condition:
        passed += 1
        print("PASS", name)
    else:
        failed += 1
        print("FAIL", name, detail if detail is not None else "")


TARGET = (
    "GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
    "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA"
)
OFF_TARGET = TARGET[:74] + "A" + TARGET[75:]
PROFILE = P.PROFILES["idt_taqman"]


def synthetic_result():
    return {
        "candidates": [{
            "rank": 1,
            "assay": {
                "forward": "GAGCTATACCCCGACCTCTG",
                "reverse": "CTTCTCTCCAATCTTGGAAAGCG",
                "probe": "TGGTGATGTCAACAACTGGGGAC",
                "amplicon": 93,
                "f_tm": 61.2,
                "r_tm": 61.4,
            },
            "evidence": {
                "hard_valid": True,
                "hard_failures": [],
                "target_coverage": 1.0,
                "offtarget": {"signal_subjects": 0, "product_subjects": 0},
                "condition_robustness": {"valid_fraction": 1.0},
                "panel_fit": None,
                "junction": None,
            },
            "rank_trace": {"priority_order": ["offtarget", "coverage"]},
        }],
        "objective_profile": {"key": "balanced", "label": "Balanced hydrolysis-probe assay"},
        "ranker_manifest": {
            "run_id": "ofrun_contract_test",
            "manifest_sha256": "a" * 64,
            "candidate_limits": {"retained_pool": 96, "full_annotation_pool": 28},
            "scientific_models": {"reaction_condition_snapshot": {
                "mv_conc_mM": 50.0,
                "dv_conc_mM": 3.0,
                "dntp_conc_mM": 0.8,
                "total_oligo_conc_nM": 250.0,
                "anneal_c": 60.0,
            }},
        },
        "search_status": "heuristic_bounded",
    }


kwargs = dict(
    workflow="contract_test",
    profile=PROFILE,
    profile_key="idt_taqman",
    objective="balanced",
    targets=[TARGET],
    off_targets=[OFF_TARGET],
    search_tier="interactive",
    search_budget_seconds=12.0,
    constraints={"declared_test": True},
)

result = synthetic_result()
contract = DC.build_contract(result, **kwargs)
contract_again = DC.build_contract(deepcopy(result), **kwargs)
verification = DC.verify_contract(contract)

check("contract build is deterministic", contract == contract_again, (contract, contract_again))
check("contract self-hash verifies", verification.get("valid") is True, verification)
check("qualified synthetic winner records declared computational limits",
      contract.get("status") == "computationally_qualified_with_declared_limits", contract.get("status"))
check("contract declares wet-lab confirmation",
      contract.get("wet_lab_confirmation_required") is True and
      "laboratory" in contract.get("scope_statement", "").lower(), contract.get("scope_statement"))
serialized = json.dumps(contract, sort_keys=True)
check("contract hashes sequence identity without embedding target sequence",
      TARGET not in serialized and OFF_TARGET not in serialized and
      contract["evidence_scope"]["targets"]["count"] == 1 and
      len(contract["evidence_scope"]["targets"]["sha256"]) == 64,
      contract.get("evidence_scope"))

tampered = deepcopy(contract)
tampered["status"] = "not_computationally_qualified"
check("contract verification detects tampering", not DC.verify_contract(tampered).get("valid"),
      DC.verify_contract(tampered))
tampered_comparison = DC.compare_contracts(contract, tampered)
check("contract comparison refuses to trust a tampered contract",
      tampered_comparison.get("state") == "invalid_contract" and
      tampered_comparison.get("checks", {}).get("right_contract_valid") is False,
      tampered_comparison)

attached = DC.attach_contract(result, **kwargs)
check("contract attachment is additive and does not mutate caller result",
      "design_contract" not in result and DC.verify_contract(attached["design_contract"])["valid"],
      attached.keys())

same = DC.compare_contracts(contract, deepcopy(contract))
check("identical contracts compare as equivalent reproduction",
      same.get("state") == "equivalent_reproduction" and not same.get("differences"), same)

changed_context = DC.build_contract(result, **dict(kwargs, targets=[TARGET + "ACGT"]))
comparison = DC.compare_contracts(contract, changed_context)
check("context changes remain explicit rather than falsely equivalent",
      comparison.get("state") == "scientifically_comparable_with_declared_context_difference" and
      "same_context" in comparison.get("differences", []), comparison)

display_variant = deepcopy(contract)
display_variant["search"]["displayed_candidates"] = 8
display_variant.pop("contract_sha256")
display_variant["contract_sha256"] = DC._sha(display_variant)
display_comparison = DC.compare_contracts(contract, display_variant)
check("presentation count is distinguished from scientific search execution",
      display_comparison.get("checks", {}).get("same_search") is True and
      display_comparison.get("checks", {}).get("same_display_count") is False and
      display_comparison.get("state") == "equivalent_reproduction", display_comparison)

discrimination_result = deepcopy(result)
discrimination_result["objective_profile"] = {
    "key": "discrimination", "label": "Species/strain discrimination",
}
missing_exclusivity = DC.build_contract(
    discrimination_result, workflow="contract_test", profile=PROFILE, profile_key="idt_taqman",
    objective="discrimination", targets=[TARGET], off_targets=[],
)
check("objective prerequisites prevent unsupported qualification claims",
      missing_exclusivity.get("status") == "insufficient_declared_evidence" and
      "off-target corpus" in missing_exclusivity["objective"]["missing_prerequisites"],
      missing_exclusivity.get("objective"))

# The generic UI default must never weaken a probe-less chemistry into the old
# hydrolysis-probe baseline.  This resolver is shared by every ranking entry point.
sybr_key = RP.resolve_objective("balanced", no_probe=True)
sybr = RP.get_profile("balanced", no_probe=True)
check("probe-less balanced default canonically resolves to SYBR", sybr_key == "sybr" and sybr["key"] == "sybr", sybr)
check("canonical SYBR objective rejects predicted off-target products",
      sybr.get("require_no_product_offtargets") is True and sybr.get("min_probe_coverage") == 0.0, sybr)
check("explicit non-default probe-less objective remains explicit",
      RP.resolve_objective("degraded_template", no_probe=True) == "degraded_template")

# Display count is presentation policy, not scientific search depth.  A fake
# candidate search makes the invariant fast and proves n is absent from the cache
# key and from the retained-corpus limit.
original_search = AD.CSEARCH.search
calls = []


def fake_search(reference, profile, **search_kwargs):
    calls.append(deepcopy(search_kwargs))
    return ([{"forward": "A" * 20, "reverse": "C" * 20, "amplicon": 90}],
            {"candidate_limits": {"retained_limit": search_kwargs["retained_limit"]}})


try:
    AD.clear_design_caches()
    AD.CSEARCH.search = fake_search
    small = AD._candidates_with_ledger(TARGET, PROFILE, n=1, budget_s=3.0)
    large = AD._candidates_with_ledger(TARGET, PROFILE, n=10, budget_s=3.0)
finally:
    AD.CSEARCH.search = original_search
    AD.clear_design_caches()

check("display count does not change or repeat the canonical search corpus",
      small == large and len(calls) == 1, {"calls": calls, "small": small, "large": large})
check("candidate search always uses the canonical retained-pool limit",
      bool(calls) and calls[0].get("retained_limit") == AD.CANONICAL_RETAINED_LIMIT == 96, calls)

print("DESIGN_CONTRACT: %d passed / %d failed / %d total" % (passed, failed, passed + failed))
raise SystemExit(1 if failed else 0)
