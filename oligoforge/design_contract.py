"""Canonical, additive quality contract for every OligoForge design workflow.

The contract does not claim laboratory performance.  It proves which shared
search/ranking policy ran, which evidence was actually supplied, whether the
recommended candidate cleared computational hard gates, and whether two results
are scientifically comparable.  Sequence content is represented only by hashes.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any, Iterable, Mapping, Optional

from . import __version__
from . import thermo as T
from .candidate_search import SEARCH_VERSION
from .candidate_retention import RETENTION_VERSION
from .provenance import THERMODYNAMIC_MODEL_VERSION, SPECIFICITY_MODEL_VERSION
from .ranking import RANKING_SCHEMA_VERSION
from .ranking_profiles import RANKER_VERSION, PROFILE_VERSION


CONTRACT_SCHEMA = "oligoforge-design-contract/v1"
CONTRACT_VERSION = "1.0.0"
POLICY_ID = "oligoforge-canonical-design-policy/v1"


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _sha(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _corpus_hash(values: Optional[Iterable[str]]) -> dict[str, Any]:
    rows = [str(x or "").upper() for x in (values or [])]
    return {"count": len(rows), "sha256": hashlib.sha256("\n".join(rows).encode()).hexdigest()}


def _profile_record(profile: Optional[Mapping[str, Any]], key: Optional[str]) -> dict[str, Any]:
    values = dict(profile or {})
    return {
        "key": key or values.get("key") or "custom",
        "name": values.get("name") or key or "custom profile",
        "snapshot_sha256": _sha(values),
        "no_probe": bool(values.get("no_probe")),
        "effective_anneal_c": float(values.get("anneal_c", T._snapshot()[4])),
    }


def _candidate_rows(result: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    if isinstance(result.get("candidates"), list):
        return [x for x in result["candidates"] if isinstance(x, Mapping)]
    if isinstance(result.get("candidate"), Mapping):
        return [result["candidate"]]
    if isinstance(result.get("redesigns"), list):
        return [x.get("candidate") for x in result["redesigns"]
                if isinstance(x, Mapping) and isinstance(x.get("candidate"), Mapping)]
    return []


def _assay(row: Mapping[str, Any]) -> Mapping[str, Any]:
    return row.get("assay") if isinstance(row.get("assay"), Mapping) else row


def _evidence(row: Mapping[str, Any]) -> Mapping[str, Any]:
    assay = _assay(row)
    value = row.get("evidence")
    if not isinstance(value, Mapping):
        value = assay.get("ranking_evidence")
    return value if isinstance(value, Mapping) else {}


def _identity(row: Mapping[str, Any]) -> dict[str, Any]:
    assay = _assay(row)
    return {
        "forward_sha256": hashlib.sha256(str(assay.get("forward") or "").encode()).hexdigest(),
        "reverse_sha256": hashlib.sha256(str(assay.get("reverse") or "").encode()).hexdigest(),
        "probe_sha256": (hashlib.sha256(str(assay.get("probe")).encode()).hexdigest()
                         if assay.get("probe") else None),
        "amplicon": assay.get("amplicon"),
    }


def _objective_record(result: Mapping[str, Any], objective: Optional[str], *,
                      target_count: int, off_target_count: int,
                      panel_count: int, junction_count: int) -> dict[str, Any]:
    declared = result.get("objective_profile")
    key = (declared.get("key") if isinstance(declared, Mapping) else None) or objective or "balanced"
    missing = []
    if key in {"discrimination", "confirmatory"} and off_target_count < 1:
        missing.append("off-target corpus")
    if key == "broad_inclusivity" and target_count < 2:
        missing.append("multiple representative target sequences")
    if key == "multiplex" and panel_count < 1:
        missing.append("comparison panel")
    if key == "transcript_specific" and junction_count < 1:
        missing.append("declared exon-junction coordinates")
    return {
        "key": key,
        "label": declared.get("label") if isinstance(declared, Mapping) else key,
        "prerequisites_met": not missing,
        "missing_prerequisites": missing,
    }


def build_contract(result: Mapping[str, Any], *, workflow: str,
                   profile: Optional[Mapping[str, Any]] = None,
                   profile_key: Optional[str] = None,
                   objective: Optional[str] = None,
                   targets: Optional[Iterable[str]] = None,
                   off_targets: Optional[Iterable[str]] = None,
                   panel_count: int = 0, junction_count: int = 0,
                   search_tier: str = "interactive",
                   search_budget_seconds: Optional[float] = None,
                   constraints: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
    """Build a deterministic contract without embedding sequence content."""
    rows = _candidate_rows(result)
    target_scope = _corpus_hash(targets)
    off_scope = _corpus_hash(off_targets)
    profile_record = _profile_record(profile, profile_key)
    objective_record = _objective_record(
        result, objective, target_count=target_scope["count"],
        off_target_count=off_scope["count"], panel_count=int(panel_count or 0),
        junction_count=int(junction_count or 0),
    )
    hard_states = [(_evidence(row).get("hard_valid")
                    if "hard_valid" in _evidence(row) else None) for row in rows]
    winner_hard_valid = hard_states[0] if hard_states else None
    ranks = [row.get("rank") for row in rows]
    winner_rank_ok = not rows or ranks[0] in (None, 1)
    rank_order_ok = all(rank is None or rank == idx + 1 for idx, rank in enumerate(ranks))

    manifest = result.get("ranker_manifest") if isinstance(result.get("ranker_manifest"), Mapping) else {}
    candidate_limits = dict(manifest.get("candidate_limits") or {})
    if not candidate_limits:
        attrition = result.get("candidate_attrition") or result.get("search_ledger") or {}
        candidate_limits = dict((attrition.get("candidate_limits") or attrition.get("limits") or {})
                                if isinstance(attrition, Mapping) else {})
    if search_budget_seconds is None:
        search_budget_seconds = candidate_limits.get("search_budget_seconds")

    performed = {"thermodynamics": False, "structured_ranking": False,
                 "target_coverage": False, "offtarget_specificity": False,
                 "condition_robustness": False, "multiplex_fit": False,
                 "junction_evidence": False}
    for row in rows:
        assay, evidence = _assay(row), _evidence(row)
        performed["thermodynamics"] |= bool(assay.get("f_tm") is not None and assay.get("r_tm") is not None)
        performed["structured_ranking"] |= bool(evidence or row.get("rank_trace") or assay.get("rank_trace"))
        performed["target_coverage"] |= evidence.get("target_coverage") is not None
        performed["offtarget_specificity"] |= isinstance(evidence.get("offtarget"), Mapping)
        performed["condition_robustness"] |= isinstance(evidence.get("condition_robustness"), Mapping)
        performed["multiplex_fit"] |= bool(evidence.get("panel_fit"))
        performed["junction_evidence"] |= bool(evidence.get("junction"))

    warnings = []
    if off_scope["count"] == 0:
        warnings.append("No off-target corpus was supplied; exclusivity is unresolved beyond the declared evidence.")
    if str(result.get("search_status") or "heuristic_bounded") == "heuristic_bounded":
        warnings.append("Candidate search is deterministic and bounded, not a proof of the global biological optimum.")
    if not objective_record["prerequisites_met"]:
        warnings.append("Objective evidence is incomplete: " + ", ".join(objective_record["missing_prerequisites"]) + ".")
    if rows and any(state is not True for state in hard_states):
        warnings.append("One or more returned candidates lack a recorded passing hard-validity state.")

    if not rows:
        status = "no_candidate"
    elif not objective_record["prerequisites_met"]:
        status = "insufficient_declared_evidence"
    elif winner_hard_valid is True and winner_rank_ok and rank_order_ok:
        status = "computationally_qualified_with_declared_limits"
    elif winner_hard_valid is False:
        status = "not_computationally_qualified"
    else:
        status = "qualification_not_reported"

    conditions = (manifest.get("scientific_models") or {}).get("reaction_condition_snapshot")
    if not isinstance(conditions, Mapping):
        snap = T._snapshot()
        conditions = {"mv_conc_mM": snap[0], "dv_conc_mM": snap[1],
                      "dntp_conc_mM": snap[2], "total_oligo_conc_nM": snap[3],
                      "anneal_c": snap[4]}
    context_identity = {
        "policy_id": POLICY_ID, "profile_sha256": profile_record["snapshot_sha256"],
        "objective": objective_record["key"], "conditions": dict(conditions),
        "targets": target_scope, "off_targets": off_scope,
        "constraints": dict(constraints or {}),
    }
    result_identity = [_identity(row) | {"rank": row.get("rank"),
                                        "hard_valid": _evidence(row).get("hard_valid")}
                       for row in rows]
    contract = {
        "schema_version": CONTRACT_SCHEMA,
        "contract_version": CONTRACT_VERSION,
        "application_version": __version__,
        "workflow": str(workflow),
        "policy_id": POLICY_ID,
        "status": status,
        "scope_statement": ("Computational design qualification for the declared inputs, models, and bounded search; "
                            "laboratory confirmation is required."),
        "wet_lab_confirmation_required": True,
        "engine": {
            "candidate_generation": "joint_multi_pair_multi_probe",
            "ranking": "structured_lexicographic_retained_pool",
            "ranker_version": RANKER_VERSION,
            "ranking_schema": RANKING_SCHEMA_VERSION,
            "objective_profile_version": PROFILE_VERSION,
            "search_version": SEARCH_VERSION,
            "retention_version": RETENTION_VERSION,
            "thermodynamic_model": THERMODYNAMIC_MODEL_VERSION,
            "specificity_model": SPECIFICITY_MODEL_VERSION,
        },
        "profile": profile_record,
        "objective": objective_record,
        "evidence_scope": {
            "targets": target_scope, "off_targets": off_scope,
            "panel_assays": int(panel_count or 0), "junctions": int(junction_count or 0),
            "performed": performed,
        },
        "search": {
            "tier": str(search_tier),
            "status": result.get("search_status") or "heuristic_bounded",
            "budget_seconds": search_budget_seconds,
            "candidate_limits": candidate_limits,
            "displayed_candidates": len(rows),
        },
        "conformance": {
            "winner_hard_valid": winner_hard_valid,
            "winner_rank_is_one": winner_rank_ok,
            "rank_order_consistent": rank_order_ok,
            "hard_valid_candidates": sum(x is True for x in hard_states),
            "candidate_count": len(rows),
        },
        "manifest": {
            "run_id": manifest.get("run_id"),
            "sha256": manifest.get("manifest_sha256"),
        },
        "warnings": warnings,
        "context_fingerprint": _sha(context_identity),
        "result_fingerprint": _sha(result_identity),
        "winner_fingerprint": (_sha(result_identity[0]) if result_identity else None),
    }
    contract["contract_sha256"] = _sha(contract)
    return contract


def attach_contract(result: Mapping[str, Any], **kwargs: Any) -> dict[str, Any]:
    out = deepcopy(dict(result))
    out["design_contract"] = build_contract(out, **kwargs)
    return out


def verify_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    value = deepcopy(dict(contract or {}))
    claimed = value.pop("contract_sha256", None)
    actual = _sha(value)
    required = {"schema_version", "contract_version", "application_version", "workflow",
                "policy_id", "status", "scope_statement", "wet_lab_confirmation_required",
                "engine", "profile", "objective", "evidence_scope", "search", "conformance",
                "context_fingerprint", "result_fingerprint"}
    missing = sorted(required - set(value))
    valid = (not missing and value.get("schema_version") == CONTRACT_SCHEMA and
             isinstance(claimed, str) and claimed == actual)
    return {"valid": valid, "claimed_sha256": claimed, "computed_sha256": actual,
            "missing_fields": missing, "schema_version": value.get("schema_version")}


def compare_contracts(left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    left_verification = verify_contract(left)
    right_verification = verify_contract(right)
    left_search = dict(left.get("search") or {})
    right_search = dict(right.get("search") or {})
    # Display count is presentation policy.  It must remain visible, but it does
    # not imply that the scientific search/annotation execution differed.
    left_display = left_search.pop("displayed_candidates", None)
    right_display = right_search.pop("displayed_candidates", None)
    checks = {
        "left_contract_valid": bool(left_verification.get("valid")),
        "right_contract_valid": bool(right_verification.get("valid")),
        "same_policy": left.get("policy_id") == right.get("policy_id"),
        "same_engine": left.get("engine") == right.get("engine"),
        "same_profile": (left.get("profile") or {}).get("snapshot_sha256") == (right.get("profile") or {}).get("snapshot_sha256"),
        "same_objective": (left.get("objective") or {}).get("key") == (right.get("objective") or {}).get("key"),
        "same_context": left.get("context_fingerprint") == right.get("context_fingerprint"),
        "same_search": left_search == right_search,
        "same_display_count": left_display == right_display,
        "same_winner": left.get("winner_fingerprint") == right.get("winner_fingerprint"),
        "same_result": left.get("result_fingerprint") == right.get("result_fingerprint"),
    }
    if not checks["left_contract_valid"] or not checks["right_contract_valid"]:
        state = "invalid_contract"
    elif checks["same_context"] and checks["same_search"] and checks["same_result"]:
        state = "equivalent_reproduction"
    elif checks["same_context"] and checks["same_search"] and checks["same_winner"]:
        state = "winner_equivalent_with_declared_result_set_difference"
    elif checks["same_context"] and checks["same_search"]:
        state = "configuration_equivalent_result_drift"
    elif checks["same_policy"] and checks["same_engine"] and checks["same_profile"] and checks["same_objective"]:
        state = "scientifically_comparable_with_declared_context_difference"
    else:
        state = "not_directly_comparable"
    return {"schema_version": "oligoforge-design-contract-comparison/v1",
            "state": state, "checks": checks,
            "differences": sorted(key for key, same in checks.items() if not same),
            "verification": {"left": left_verification, "right": right_verification},
            "interpretation": ("Matching contracts prove software/configuration reproducibility only; "
                               "they do not establish wet-lab equivalence.")}
