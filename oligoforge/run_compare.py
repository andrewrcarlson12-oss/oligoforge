"""Deterministic comparison of two OligoForge design/ranking runs.

The comparison distinguishes scientific-context changes from unexplained ordering
changes.  It never treats two runs as directly comparable merely because their
visible assay names match; candidates are identified by normalized oligo triplets.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from . import provenance as PROV

RUN_COMPARE_VERSION = "1.0.0"
MAX_COMPARE_CANDIDATES = 500


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _triplet(row: Mapping[str, Any]) -> Tuple[str, str, str]:
    assay = row.get("assay") if isinstance(row.get("assay"), Mapping) else row
    return tuple(str((assay or {}).get(k) or "").upper() for k in ("forward", "reverse", "probe"))


def candidate_id(row: Mapping[str, Any]) -> str:
    f, r, p = _triplet(row)
    return "ofcand_" + hashlib.sha256((f + "|" + r + "|" + p).encode()).hexdigest()[:20]


def _candidate_source(run: Any) -> List[Mapping[str, Any]]:
    if isinstance(run, list):
        rows = run
    elif isinstance(run, Mapping):
        rows = None
        for key in ("candidates", "ranked_candidates", "finalists", "results", "rows", "designs"):
            if isinstance(run.get(key), list):
                rows = run.get(key)
                break
        if rows is None and isinstance(run.get("candidate"), Mapping):
            rows = [run.get("candidate")]
        if rows is None and any(k in run for k in ("forward", "reverse", "assay")):
            rows = [run]
        rows = rows or []
    else:
        raise ValueError("each design run must be an object or candidate array")
    if len(rows) > MAX_COMPARE_CANDIDATES:
        raise ValueError("design run contains %d candidates; comparison limit %d" %
                         (len(rows), MAX_COMPARE_CANDIDATES))
    out = []
    for index, raw in enumerate(rows):
        if not isinstance(raw, Mapping):
            raise ValueError("candidate %d is not an object" % (index + 1))
        f, r, p = _triplet(raw)
        if not f or not r:
            raise ValueError("candidate %d requires forward and reverse primer sequences" % (index + 1))
        out.append(raw)
    return out


def _finite(value: Any) -> Optional[float]:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _evidence_snapshot(row: Mapping[str, Any]) -> Dict[str, Any]:
    e = row.get("evidence") if isinstance(row.get("evidence"), Mapping) else {}
    a = row.get("assay") if isinstance(row.get("assay"), Mapping) else row
    off = e.get("offtarget") if isinstance(e.get("offtarget"), Mapping) else {}
    robust = e.get("condition_robustness") if isinstance(e.get("condition_robustness"), Mapping) else {}
    return {
        "hard_valid": bool(e.get("hard_valid")),
        "hard_failures": sorted(str(x) for x in (e.get("hard_failures") or [])),
        "target_coverage": _finite(e.get("target_coverage")),
        "worst_isolate_3prime": _finite(e.get("worst_isolate_3prime")),
        "probe_mean_identity": _finite(e.get("probe_mean_identity")),
        "signal_offtargets": int(off.get("signal_subjects", 0) or 0),
        "product_offtargets": int(off.get("product_subjects", 0) or 0),
        "condition_robustness": _finite(robust.get("valid_fraction")),
        "triplet_penalty": _finite(e.get("triplet_penalty")),
        "practical_penalty": _finite(e.get("practical_penalty")),
        "degeneracy_fold": int(e.get("degeneracy_fold", 1) or 1),
        "panel_risk": int(e.get("panel_risk", 0) or 0),
        "amplicon": int((a or {}).get("amplicon", 0) or 0),
        "display_score": _finite(row.get("display_score")),
        "preference_state": ((row.get("rank_explanation") or {}).get("preference_state")
                             if isinstance(row.get("rank_explanation"), Mapping) else None),
    }


def normalize_run(run: Any, label: str = "run") -> Dict[str, Any]:
    rows = _candidate_source(run)
    normalized, seen = [], set()
    duplicate_ids = []
    for index, raw in enumerate(rows):
        cid = candidate_id(raw)
        if cid in seen:
            duplicate_ids.append(cid)
            continue
        seen.add(cid)
        assay = raw.get("assay") if isinstance(raw.get("assay"), Mapping) else raw
        normalized.append({
            "candidate_id": cid,
            "rank": int(raw.get("rank") or raw.get("candidate_rank") or index + 1),
            "forward": str((assay or {}).get("forward") or "").upper(),
            "reverse": str((assay or {}).get("reverse") or "").upper(),
            "probe": str((assay or {}).get("probe") or "").upper() or None,
            "region_start": int(((assay or {}).get("f_xy") or [0])[0] or 0),
            "finalist_categories": sorted(str(x) for x in (raw.get("finalist_categories") or [])),
            "evidence": _evidence_snapshot(raw),
        })
    normalized.sort(key=lambda x: (x["rank"], x["candidate_id"]))
    for index, row in enumerate(normalized):
        row["normalized_rank"] = index + 1

    manifest = None
    if isinstance(run, Mapping):
        manifest = run.get("ranker_manifest") or run.get("manifest")
    manifest = dict(manifest or {}) if isinstance(manifest, Mapping) else {}
    verification = PROV.verify_manifest(manifest) if manifest else {
        "valid": False, "digest_valid": False, "run_id_valid": False,
        "stored_run_id": None, "stored_manifest_sha256": None,
    }
    return {
        "version": RUN_COMPARE_VERSION,
        "label": str(label),
        "manifest": manifest,
        "manifest_verification": verification,
        "candidates": normalized,
        "n_candidates": len(normalized),
        "duplicate_candidate_ids_suppressed": duplicate_ids,
    }


def _context(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    models = manifest.get("scientific_models") if isinstance(manifest.get("scientific_models"), Mapping) else {}
    return {
        "application_version": manifest.get("application_version"),
        "ranker_version": manifest.get("ranker_version"),
        "ranking_schema": manifest.get("ranking_schema"),
        "scoring_profile_version": manifest.get("scoring_profile_version"),
        "search_version": manifest.get("search_version"),
        "retention_version": manifest.get("retention_version"),
        "objective": manifest.get("objective"),
        "scientific_models": models,
        "candidate_limits": manifest.get("candidate_limits") or {},
        "constraints": manifest.get("constraints") or {},
        "input_hashes": manifest.get("input_hashes") or {},
        "external_databases": manifest.get("external_databases") or {},
        "external_database_state": manifest.get("external_database_state"),
        "random_seed": manifest.get("random_seed"),
        "deterministic": manifest.get("deterministic"),
    }


def _context_differences(a: Mapping[str, Any], b: Mapping[str, Any]) -> List[Dict[str, Any]]:
    ca, cb = _context(a), _context(b)
    rows = []
    for key in ca:
        if _canonical(ca.get(key)) != _canonical(cb.get(key)):
            rows.append({"field": key, "left": ca.get(key), "right": cb.get(key)})
    return rows


def _spearman(shared: Iterable[str], left_rank: Mapping[str, int], right_rank: Mapping[str, int]) -> Optional[float]:
    ids = list(shared)
    n = len(ids)
    if n < 2:
        return None
    d2 = sum((left_rank[x] - right_rank[x]) ** 2 for x in ids)
    return round(1.0 - (6.0 * d2) / (n * (n * n - 1)), 6)


def _pairwise_reversals(shared: Iterable[str], left_rank: Mapping[str, int], right_rank: Mapping[str, int]) -> Dict[str, Any]:
    ids = sorted(shared)
    reversed_pairs = []
    total = 0
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            total += 1
            left_order = left_rank[a] < left_rank[b]
            right_order = right_rank[a] < right_rank[b]
            if left_order != right_order:
                reversed_pairs.append([a, b])
    return {
        "n_shared_pairs": total,
        "n_reversed_pairs": len(reversed_pairs),
        "fraction_reversed": round(len(reversed_pairs) / total, 6) if total else None,
        "reversed_pairs": reversed_pairs[:100],
        "truncated": len(reversed_pairs) > 100,
    }


def _numeric_delta(left: Any, right: Any) -> Optional[float]:
    a, b = _finite(left), _finite(right)
    return round(b - a, 6) if a is not None and b is not None else None


def _candidate_changes(left: Mapping[str, Any], right: Mapping[str, Any]) -> Dict[str, Any]:
    le, re = left["evidence"], right["evidence"]
    return {
        "candidate_id": left["candidate_id"],
        "left_rank": left["normalized_rank"],
        "right_rank": right["normalized_rank"],
        "rank_shift": left["normalized_rank"] - right["normalized_rank"],
        "hard_valid_changed": le["hard_valid"] != re["hard_valid"],
        "hard_failures_added": sorted(set(re["hard_failures"]) - set(le["hard_failures"])),
        "hard_failures_removed": sorted(set(le["hard_failures"]) - set(re["hard_failures"])),
        "evidence_deltas": {
            key: _numeric_delta(le.get(key), re.get(key))
            for key in ("target_coverage", "worst_isolate_3prime", "probe_mean_identity",
                        "signal_offtargets", "product_offtargets", "condition_robustness",
                        "triplet_penalty", "practical_penalty", "degeneracy_fold",
                        "panel_risk", "amplicon", "display_score")
        },
        "finalist_categories_added": sorted(set(right["finalist_categories"]) - set(left["finalist_categories"])),
        "finalist_categories_removed": sorted(set(left["finalist_categories"]) - set(right["finalist_categories"])),
    }


def compare_runs(left: Any, right: Any, top_k: int = 10) -> Dict[str, Any]:
    top_k = max(1, min(int(top_k), 100))
    a, b = normalize_run(left, "left"), normalize_run(right, "right")
    amap = {x["candidate_id"]: x for x in a["candidates"]}
    bmap = {x["candidate_id"]: x for x in b["candidates"]}
    shared = set(amap) & set(bmap)
    added = sorted(set(bmap) - set(amap), key=lambda x: bmap[x]["normalized_rank"])
    removed = sorted(set(amap) - set(bmap), key=lambda x: amap[x]["normalized_rank"])
    arank = {x: amap[x]["normalized_rank"] for x in shared}
    brank = {x: bmap[x]["normalized_rank"] for x in shared}
    context_diffs = _context_differences(a["manifest"], b["manifest"])
    both_manifests_valid = bool(a["manifest_verification"].get("valid") and
                                b["manifest_verification"].get("valid"))
    exact_manifest_match = bool(both_manifests_valid and
                                a["manifest"].get("manifest_sha256") == b["manifest"].get("manifest_sha256"))
    same_scientific_context = bool(both_manifests_valid and not context_diffs)
    left_order = [x["candidate_id"] for x in a["candidates"]]
    right_order = [x["candidate_id"] for x in b["candidates"]]
    same_candidate_set = set(left_order) == set(right_order)
    same_order = left_order == right_order

    if exact_manifest_match and (not same_candidate_set or not same_order):
        reproducibility = {
            "state": "critical_non_reproducibility",
            "reason": "identical self-hashed manifests produced a different candidate set or order",
        }
    elif same_scientific_context and same_candidate_set and not same_order:
        reproducibility = {
            "state": "unexplained_order_change",
            "reason": "scientific context matches but shared candidates changed order",
        }
    elif same_scientific_context and same_order:
        reproducibility = {
            "state": "reproduced",
            "reason": "scientific context and deterministic ordering match",
        }
    elif context_diffs:
        reproducibility = {
            "state": "context_changed",
            "reason": "ordering changes must be interpreted in light of the listed scientific-context differences",
        }
    else:
        reproducibility = {
            "state": "insufficient_provenance",
            "reason": "one or both runs lack a complete comparable manifest",
        }

    left_top = left_order[:top_k]; right_top = right_order[:top_k]
    top_union = set(left_top) | set(right_top)
    top_intersection = set(left_top) & set(right_top)
    changes = [_candidate_changes(amap[x], bmap[x]) for x in shared]
    changes.sort(key=lambda x: (-abs(x["rank_shift"]), x["right_rank"], x["candidate_id"]))
    lw = a["candidates"][0] if a["candidates"] else None
    rw = b["candidates"][0] if b["candidates"] else None
    winner_changed = bool(lw and rw and lw["candidate_id"] != rw["candidate_id"])

    return {
        "version": RUN_COMPARE_VERSION,
        "left": {k: a[k] for k in ("label", "n_candidates", "manifest_verification", "duplicate_candidate_ids_suppressed")},
        "right": {k: b[k] for k in ("label", "n_candidates", "manifest_verification", "duplicate_candidate_ids_suppressed")},
        "reproducibility": reproducibility,
        "both_manifests_valid": both_manifests_valid,
        "manifest_exact_match": exact_manifest_match,
        "same_scientific_context": same_scientific_context,
        "context_differences": context_diffs,
        "candidate_set": {
            "n_shared": len(shared), "n_added": len(added), "n_removed": len(removed),
            "added": [bmap[x] for x in added], "removed": [amap[x] for x in removed],
        },
        "ranking_stability": {
            "spearman_shared": _spearman(shared, arank, brank),
            **_pairwise_reversals(shared, arank, brank),
            "top_k": top_k,
            "top_k_overlap": len(top_intersection),
            "top_k_jaccard": round(len(top_intersection) / len(top_union), 6) if top_union else None,
        },
        "winner": {
            "changed": winner_changed,
            "left": lw,
            "right": rw,
            "explanation": (
                "winner changed after a declared scientific-context change" if winner_changed and context_diffs else
                "winner changed without a declared scientific-context change" if winner_changed else
                "same normalized primer/probe triplet remains rank 1"
            ),
        },
        "shared_candidate_changes": changes,
        "interpretation": (
            "Rank movement is descriptive. It does not establish that either run is biologically superior; inspect context differences and candidate evidence."
        ),
    }
