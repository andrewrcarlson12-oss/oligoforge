"""Transparent, versioned records for assay bench outcomes.

Feedback remains evidence attached to an assay.  It never silently changes the
authoritative structured ranker.  This module validates records, deduplicates them,
reports conflicts/missingness, creates leakage-safe target-group splits, and decides
whether a future calibrated reranker is even eligible for development.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import Counter, defaultdict
from typing import Any, Dict, Iterable, List

SCHEMA_VERSION = "1.3.0"
DATASET_SCHEMA_VERSION = "1.0.0"
ALLOWED_STATUS = {"untested", "success", "failed", "partial", "retired"}
ALLOWED_DESIGNATION = {
    "none", "locked", "preferred", "experimental", "previously_ordered",
    "previously_validated", "failed_at_bench", "retired",
}
ALLOWED_SIGNAL = {None, "", "poor", "weak", "acceptable", "good", "strong", "not_measured"}


def _finite(value, field, *, minimum=None, maximum=None):
    if value is None or value == "":
        return None
    try:
        val = float(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be numeric" % field)
    if not math.isfinite(val):
        raise ValueError("%s must be finite" % field)
    if minimum is not None and val < minimum:
        raise ValueError("%s must be >= %s" % (field, minimum))
    if maximum is not None and val > maximum:
        raise ValueError("%s must be <= %s" % (field, maximum))
    return val


def _melt_peaks(value):
    if value in (None, ""):
        return None
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except Exception:
            value = [x.strip() for x in value.split(";") if x.strip()]
    if not isinstance(value, (list, tuple)):
        raise ValueError("melt_peaks must be a list of temperatures")
    rows = [_finite(x, "melt_peaks", minimum=0.0, maximum=120.0) for x in value]
    return rows[:100]


def _int_nonnegative(value, field):
    if value in (None, ""):
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be a non-negative integer" % field)
    if not math.isfinite(f) or f < 0 or int(f) != f:
        raise ValueError("%s must be a non-negative integer" % field)
    return int(f)


def normalize(record):
    r = dict(record or {})
    status = str(r.get("status") or "untested").strip().lower()
    if status not in ALLOWED_STATUS:
        raise ValueError("status must be one of: %s" % ", ".join(sorted(ALLOWED_STATUS)))
    designation = str(r.get("designation") or "none").strip().lower()
    if designation not in ALLOWED_DESIGNATION:
        raise ValueError("designation must be one of: %s" % ", ".join(sorted(ALLOWED_DESIGNATION)))
    signal = str(r.get("probe_signal_quality") or "").strip().lower() or None
    if signal not in ALLOWED_SIGNAL:
        raise ValueError("probe_signal_quality must be poor, weak, acceptable, good, strong, or not_measured")
    out = dict(
        schema_version=SCHEMA_VERSION,
        status=status,
        designation=designation,
        assay_id=(str(r.get("assay_id")).strip() if r.get("assay_id") is not None else None),
        target_group=(str(r.get("target_group")).strip() if r.get("target_group") is not None else None),
        design_run_id=(str(r.get("design_run_id")).strip() if r.get("design_run_id") is not None else None),
        oligoforge_version=r.get("oligoforge_version"),
        ranker_version=r.get("ranker_version"),
        objective=r.get("objective"),
        input_hash=r.get("input_hash"),
        conditions=r.get("conditions"),
        efficiency=_finite(r.get("efficiency"), "efficiency", minimum=0.0, maximum=300.0),
        r2=_finite(r.get("r2"), "r2", minimum=0.0, maximum=1.0),
        replicate_cv=_finite(r.get("replicate_cv"), "replicate_cv", minimum=0.0, maximum=500.0),
        cq_at_declared_input=_finite(r.get("cq_at_declared_input"), "cq_at_declared_input", minimum=0.0, maximum=100.0),
        lod_observation=r.get("lod_observation"),
        nonspecific_products=_int_nonnegative(r.get("nonspecific_products"), "nonspecific_products"),
        product_identity_method=r.get("product_identity_method"),
        melt_peak_count=_int_nonnegative(
            r.get("melt_peak_count") if r.get("melt_peak_count") not in (None, "")
            else (r.get("melt_peaks") if isinstance(r.get("melt_peaks"), (int, float)) else None),
            "melt_peak_count"),
        melt_peaks=(None if isinstance(r.get("melt_peaks"), (int, float)) else _melt_peaks(r.get("melt_peaks"))),
        probe_signal_quality=signal,
        multiplex_failure=bool(r.get("multiplex_failure", False)),
        notes=(str(r.get("notes"))[:10000] if r.get("notes") is not None else None),
    )
    if status in {"success", "failed", "partial"} and not out["assay_id"]:
        raise ValueError("labeled feedback requires assay_id")
    if out["r2"] is not None and out["efficiency"] is None:
        # R2 alone is allowed but is a data-completeness warning, not an error.
        pass
    out["record_sha256"] = hashlib.sha256(
        json.dumps(out, sort_keys=True, default=str, separators=(",", ":")).encode()
    ).hexdigest()
    return out


def parse_records(payload: str, format_hint: str = "auto") -> List[Dict[str, Any]]:
    """Parse a JSON array/object or CSV text into raw feedback records."""
    text = str(payload or "").lstrip("\ufeff").strip()
    if not text:
        return []
    hint = str(format_hint or "auto").lower()
    if hint not in {"auto", "json", "csv"}:
        raise ValueError("format_hint must be auto, json, or csv")
    if hint == "json" or (hint == "auto" and text[:1] in "[{"):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError("invalid feedback JSON: %s" % exc.msg)
        if isinstance(data, dict):
            data = data.get("records", [data])
        if not isinstance(data, list):
            raise ValueError("feedback JSON must contain an array of records")
        return [dict(x or {}) for x in data]
    try:
        return [dict(row) for row in csv.DictReader(io.StringIO(text))]
    except Exception as exc:
        raise ValueError("invalid feedback CSV: %s" % type(exc).__name__)


def dataset_status(records):
    normalized, rejected = [], []
    for idx, record in enumerate(records or []):
        try:
            normalized.append(normalize(record))
        except ValueError as exc:
            rejected.append({"row": idx + 1, "error": str(exc)})
    unique, duplicate_hashes = [], []
    seen = set()
    for row in normalized:
        h = row["record_sha256"]
        if h in seen:
            duplicate_hashes.append(h)
        else:
            seen.add(h); unique.append(row)

    by_identity = defaultdict(list)
    for row in unique:
        key = (row.get("assay_id"), row.get("design_run_id"),
               json.dumps(row.get("conditions"), sort_keys=True, default=str))
        by_identity[key].append(row)
    conflicts = []
    for key, rows in by_identity.items():
        outcomes = sorted({r["status"] for r in rows if r["status"] in {"success", "failed", "partial"}})
        if len(outcomes) > 1:
            conflicts.append({"assay_id": key[0], "design_run_id": key[1],
                              "outcomes": outcomes, "record_hashes": [r["record_sha256"] for r in rows]})

    labeled = [r for r in unique if r["status"] in {"success", "failed", "partial"}]
    required_model_fields = ("target_group", "assay_id", "ranker_version", "objective")
    missingness = {field: sum(not bool(r.get(field)) for r in labeled) for field in required_model_fields}
    quantitative_fields = ("efficiency", "r2", "replicate_cv", "cq_at_declared_input")
    quantitative_completeness = {
        field: (sum(r.get(field) is not None for r in labeled) / len(labeled) if labeled else 0.0)
        for field in quantitative_fields
    }
    return {
        "schema_version": DATASET_SCHEMA_VERSION,
        "n_input": len(records or []),
        "n_normalized": len(normalized),
        "n_unique": len(unique),
        "n_duplicates": len(duplicate_hashes),
        "n_rejected": len(rejected),
        "rejected": rejected,
        "duplicate_record_hashes": duplicate_hashes,
        "conflicts": conflicts,
        "n_conflicts": len(conflicts),
        "missingness_labeled": missingness,
        "quantitative_completeness": quantitative_completeness,
        "outcomes": dict(Counter(r["status"] for r in labeled)),
        "normalized_records": unique,
    }


def target_group_split(records, train=0.70, validation=0.15):
    """Deterministically split whole target groups; no group can leak across sets."""
    ds = dataset_status(records)
    if ds["rejected"]:
        raise ValueError("feedback dataset contains rejected records")
    if train <= 0 or validation < 0 or train + validation >= 1:
        raise ValueError("split fractions must leave a non-zero test set")
    assignments, rows = {}, []
    for row in ds["normalized_records"]:
        group = row.get("target_group")
        if not group:
            raise ValueError("every split record requires target_group")
        if group not in assignments:
            x = int(hashlib.sha256(group.encode()).hexdigest()[:12], 16) / float(16**12)
            assignments[group] = "train" if x < train else ("validation" if x < train + validation else "test")
        copy = dict(row); copy["dataset_split"] = assignments[group]; rows.append(copy)
    sets = defaultdict(set)
    for group, split in assignments.items():
        sets[split].add(group)
    leakage = bool((sets["train"] & sets["validation"]) or (sets["train"] & sets["test"]) or
                   (sets["validation"] & sets["test"]))
    return {
        "schema_version": DATASET_SCHEMA_VERSION,
        "policy": "sha256 target-group split; deterministic and leakage-safe",
        "fractions": {"train": train, "validation": validation, "test": 1-train-validation},
        "n_groups": {k: len(v) for k, v in sets.items()},
        "group_leakage": leakage,
        "assignments": assignments,
        "records": rows,
    }


def calibration_status(records):
    ds = dataset_status(records)
    rows = ds["normalized_records"]
    groups = {r.get("target_group") for r in rows if r.get("target_group")}
    labeled = [r for r in rows if r["status"] in {"success", "failed", "partial"}]
    successes = sum(r["status"] == "success" for r in labeled)
    failures = sum(r["status"] == "failed" for r in labeled)
    decisive = successes + failures
    minority_fraction = min(successes, failures) / decisive if decisive else 0.0
    complete_groups = sum(bool(r.get("target_group")) for r in labeled) == len(labeled)
    versioned = sum(bool(r.get("ranker_version")) for r in labeled) == len(labeled)
    conflict_free = ds["n_conflicts"] == 0
    adequate = (
        len(labeled) >= 100 and len(groups) >= 10 and successes >= 25 and failures >= 25 and
        minority_fraction >= 0.20 and complete_groups and versioned and conflict_free and not ds["rejected"]
    )
    unmet = []
    if len(labeled) < 100: unmet.append("fewer than 100 labeled records")
    if len(groups) < 10: unmet.append("fewer than 10 target groups")
    if successes < 25 or failures < 25: unmet.append("fewer than 25 successes or 25 failures")
    if minority_fraction < 0.20: unmet.append("outcome imbalance exceeds the 80/20 gate")
    if not complete_groups: unmet.append("target_group missing from labeled records")
    if not versioned: unmet.append("ranker_version missing from labeled records")
    if not conflict_free: unmet.append("conflicting outcomes require adjudication")
    if ds["rejected"]: unmet.append("invalid records were rejected")
    return dict(
        schema_version=SCHEMA_VERSION,
        n_records=ds["n_unique"], n_labeled=len(labeled), n_target_groups=len(groups),
        n_success=successes, n_failed=failures, n_partial=sum(r["status"] == "partial" for r in labeled),
        minority_outcome_fraction=round(minority_fraction, 4),
        learned_reranker_allowed=adequate,
        evidence_gate_only=True,
        reason=(
            "minimum evidence gate met; target-group split, leakage-controlled cross-validation, calibration, ablation, and held-out improvement are still required"
            if adequate else "insufficient or internally inconsistent evidence for a defensible learned reranker"
        ),
        unmet_requirements=unmet,
        fallback="transparent structured ranker remains authoritative",
        dataset=ds,
        normalized_records=rows,
    )


def _median(values):
    vals = sorted(float(x) for x in values if x is not None and math.isfinite(float(x)))
    if not vals:
        return None
    mid = len(vals) // 2
    return vals[mid] if len(vals) % 2 else (vals[mid - 1] + vals[mid]) / 2.0


def evidence_summary(records):
    """Aggregate local bench evidence without silently changing the ranker.

    Summaries are grouped by assay ID and by comparable target/condition contexts.
    A success in one matrix or target group is never generalized to another.
    """
    ds = dataset_status(records)
    if ds["rejected"]:
        raise ValueError("feedback dataset contains rejected records")
    rows = ds["normalized_records"]
    by_assay = defaultdict(list)
    by_context = defaultdict(lambda: defaultdict(list))
    for row in rows:
        aid = row.get("assay_id")
        if not aid:
            continue
        by_assay[aid].append(row)
        context = (
            row.get("target_group") or "unassigned",
            json.dumps(row.get("conditions") or {}, sort_keys=True, default=str, separators=(",", ":")),
        )
        by_context[context][aid].append(row)

    assay_rows = []
    for aid, vals in sorted(by_assay.items()):
        labeled = [x for x in vals if x["status"] in {"success", "failed", "partial"}]
        decisive = [x for x in labeled if x["status"] in {"success", "failed"}]
        counts = Counter(x["status"] for x in labeled)
        if len(decisive) < 2:
            state = "insufficient_local_evidence"
        elif counts["success"] and counts["failed"]:
            state = "mixed_local_evidence"
        elif counts["success"]:
            state = "consistent_local_success"
        else:
            state = "consistent_local_failure"
        contexts = {(x.get("target_group"), json.dumps(x.get("conditions") or {}, sort_keys=True, default=str))
                    for x in vals}
        assay_rows.append({
            "assay_id": aid,
            "evidence_state": state,
            "n_records": len(vals),
            "n_labeled": len(labeled),
            "n_contexts": len(contexts),
            "outcomes": dict(counts),
            "success_fraction_decisive": (round(counts["success"] / len(decisive), 4) if decisive else None),
            "median_efficiency": _median(x.get("efficiency") for x in vals),
            "median_r2": _median(x.get("r2") for x in vals),
            "median_replicate_cv": _median(x.get("replicate_cv") for x in vals),
            "n_nonspecific_product_records": sum((x.get("nonspecific_products") or 0) > 0 for x in vals),
            "n_multiplex_failures": sum(bool(x.get("multiplex_failure")) for x in vals),
            "record_hashes": [x["record_sha256"] for x in vals],
            "scope": "local assay evidence; do not generalize across target groups, matrices, conditions, or laboratories",
        })

    pair_votes = defaultdict(lambda: Counter())
    comparison_rows = []
    for (target_group, condition_key), assays in sorted(by_context.items()):
        summaries = {}
        for aid, vals in assays.items():
            decisive = [x for x in vals if x["status"] in {"success", "failed"}]
            if not decisive:
                continue
            counts = Counter(x["status"] for x in decisive)
            # Only create a preference when the local context has a unanimous
            # decisive outcome for each assay. Mixed histories are left unresolved.
            if counts["success"] and counts["failed"]:
                continue
            summaries[aid] = "success" if counts["success"] else "failed"
        aids = sorted(summaries)
        for i, a in enumerate(aids):
            for b in aids[i + 1:]:
                if summaries[a] == summaries[b]:
                    continue
                preferred = a if summaries[a] == "success" else b
                other = b if preferred == a else a
                pair = tuple(sorted((a, b)))
                pair_votes[pair][preferred] += 1
                comparison_rows.append({
                    "target_group": target_group,
                    "conditions_sha256": hashlib.sha256(condition_key.encode()).hexdigest(),
                    "preferred_assay_id": preferred,
                    "nonpreferred_assay_id": other,
                    "basis": "unanimous success versus unanimous failure in the same declared target/condition context",
                })

    pairwise = []
    for pair, votes in sorted(pair_votes.items()):
        total = sum(votes.values())
        winner, nwin = votes.most_common(1)[0]
        conflict = len(votes) > 1
        pairwise.append({
            "assay_ids": list(pair),
            "preferred_assay_id": None if conflict else winner,
            "n_comparable_contexts": total,
            "votes": dict(votes),
            "conflict": conflict,
            "state": "conflicting_local_preferences" if conflict else "local_pairwise_preference",
            "scope": "descriptive local evidence only; not a universal assay-quality label",
        })

    return {
        "schema_version": SCHEMA_VERSION,
        "n_assays": len(assay_rows),
        "n_contexts": len(by_context),
        "assays": assay_rows,
        "pairwise_preferences": pairwise,
        "context_comparisons": comparison_rows,
        "dataset": {k: v for k, v in ds.items() if k != "normalized_records"},
        "ranker_policy": "feedback is displayed and exported but does not silently alter the authoritative structured rank",
    }
