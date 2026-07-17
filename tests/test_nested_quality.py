#!/usr/bin/env python3
"""Focused nested-PCR structured-quality regression gates."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oligoforge import autodesign as AD
from oligoforge import design as D
from oligoforge import design_contract as DC
from oligoforge import profiles as P
from oligoforge import thermo as T


passed = failed = 0


def check(name, condition, detail=None):
    global passed, failed
    if condition:
        passed += 1
        print("PASS", name)
    else:
        failed += 1
        print("FAIL", name, detail if detail is not None else "")


# Use the long-standing HMBS regression primers as a thermodynamically valid
# outer pair, embedded around an artificial inner product.  Candidate enumeration
# is isolated so this test exercises nested geometry and the real structured
# rank/evidence/manifest engines without doing an exhaustive sliding search.
OUTER_F = "GAGCTATACCCCGACCTCTG"
OUTER_R = "CTTCTCTCCAATCTTGGAAAGCG"
OUTER_R_SITE = T.revcomp(OUTER_R)
INNER_F = "AACCGTTAAGGCTACGTTCA"
REFERENCE = (
    "ATCG" * 25 + OUTER_F + "AATTCCGG" * 15 + INNER_F +
    "CGAT" * 18 + "GCTA" * 25 + OUTER_R_SITE + "TACG" * 40
)
INNER = {"forward": INNER_F, "reverse": "ACGT" * 5, "amplicon": 90}

original_enumerate = D.enumerate_primers


def controlled_enumeration(template, profile, audit=None):
    fwd, rev = [], []
    fi = template.find(OUTER_F)
    if fi >= 0:
        fwd.append((fi, fi + len(OUTER_F), OUTER_F))
    ri = template.find(OUTER_R_SITE)
    if ri >= 0:
        rev.append((ri, ri + len(OUTER_R), OUTER_R))
    return fwd, rev


try:
    D.enumerate_primers = controlled_enumeration
    result = AD.design_nested(
        REFERENCE, P.PROFILES["idt_taqman"], INNER,
        targets=[REFERENCE], objective="balanced", n_candidates=3,
    )
    excluded = AD.design_nested(
        REFERENCE, P.PROFILES["idt_taqman"], INNER,
        targets=[REFERENCE], offs=[REFERENCE], objective="balanced",
    )
finally:
    D.enumerate_primers = original_enumerate


check("nested structured design returns an outer recommendation",
      bool(result and result.get("outer") and result.get("candidates")),
      result and result.get("constraint_note"))
if result and result.get("outer"):
    outer = result["outer"]
    evidence = outer.get("evidence") or {}
    check("probe-less balanced default resolves to the SYBR objective",
          result["objective_profile"].get("key") == "sybr" and
          evidence.get("objective") == "sybr",
          (result.get("objective_profile"), evidence.get("objective")))
    check("outer winner is hard-valid with full target/robustness evidence",
          evidence.get("hard_valid") is True and
          (evidence.get("target") or {}).get("product_subjects") == 1 and
          (evidence.get("evaluations") or {}).get("target_epcr") is True and
          (evidence.get("evaluations") or {}).get("condition_robustness") is True,
          evidence.get("hard_failures"))
    check("nested geometry survives structured candidate conversion",
          outer["f_xy"][1] <= REFERENCE.find(INNER_F) - 8 and
          outer["r_xy"][0] >= REFERENCE.find(INNER_F) + INNER["amplicon"] + 8 and
          outer.get("f_outside", 0) >= 8 and outer.get("r_outside", 0) >= 8 and
          (outer.get("nested_geometry") or {}).get("fully_nested") is True,
          (outer.get("f_xy"), outer.get("r_xy"), outer.get("nested_geometry")))
    check("compatibility outer exposes authoritative structured rank provenance",
          outer.get("rank") == 1 and outer.get("rank_trace") and
          outer.get("rank_explanation") and
          "structured rank" in outer.get("score_semantics", ""),
          outer)

manifest = (result or {}).get("ranker_manifest") or {}
check("nested path emits a versioned ranker manifest with primer-only constraints",
      bool(manifest.get("manifest_sha256")) and
      (manifest.get("constraints") or {}).get("workflow") == "nested_outer_primer_design" and
      (manifest.get("constraints") or {}).get("primer_only_semantics") is True,
      manifest)
contract = (result or {}).get("design_contract") or {}
check("nested path emits a verifiable canonical design contract",
      DC.verify_contract(contract).get("valid") is True and
      contract.get("workflow") == "nested_outer_primer_design",
      DC.verify_contract(contract))

stages = ((result or {}).get("candidate_attrition") or {}).get("stages") or []
check("every nested attrition stage balances entered/retained/rejected",
      bool(stages) and all(
          stage.get("entered") == stage.get("retained") + stage.get("rejected")
          for stage in stages
          if all(key in stage for key in ("entered", "retained", "rejected"))
      ), stages)

# The identical sequence is a known off-target product.  Under canonical SYBR
# semantics it must be a hard failure and must never leak through ``outer`` merely
# because the pair has attractive Tm or geometry.
rejected = ((excluded or {}).get("top_rejected") or [{}])[0]
rejected_evidence = rejected.get("evidence") or {}
check("predicted outer-pair off-target product is a hard SYBR failure",
      (excluded or {}).get("objective_profile", {}).get("key") == "sybr" and
      (excluded or {}).get("outer") is None and
      "off-target product predicted" in rejected_evidence.get("hard_failures", []),
      rejected_evidence.get("hard_failures"))


print("NESTED QUALITY:", passed, "passed /", failed, "failed")
if failed:
    sys.exit(1)
