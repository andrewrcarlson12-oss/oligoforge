"""Complete-product sequence-snapshot monitoring for registered assays."""
from __future__ import annotations

from collections import defaultdict

from .. import isolates as ISO
from ..provenance import sha256_value, SPECIFICITY_MODEL_VERSION
from .assaysbom import build_assaysbom
from .snapshots import snapshot_delta, validate_snapshot


DRIFTGUARD_VERSION = "1.0.0"


def _oligos(assay):
    by_role = defaultdict(list)
    for c in assay.get("components", []):
        by_role[c.get("role")].append(c)
    return by_role


def _primary_triplet(assay):
    roles = _oligos(assay)
    return (roles["forward_primer"][0], roles["reverse_primer"][0],
            roles["probe"][0] if roles["probe"] else None)


def _record_eval(assay, record, model):
    f, r, p = _primary_triplet(assay)
    result = ISO.amplify(f["sequence"], r["sequence"], (p or {}).get("sequence") or "",
                         record["sequence"], max_mm=model["max_mm"], clamp_n=model["clamp_n"],
                         min_product=model["min_product"], max_product=model["max_product"],
                         min_probe_ident=model["min_probe_ident"])
    signal = bool(result.get("amplifies") and (not p or result.get("probe_binds")))
    component_concerns = []
    reason = str(result.get("reason") or "")
    if not result.get("amplifies"):
        if "forward" in reason:
            component_concerns.append(f["component_id"])
        if "reverse" in reason:
            component_concerns.append(r["component_id"])
    elif p and not result.get("probe_binds"):
        component_concerns.append(p["component_id"])
    return {
        "record_id": record["record_id"], "sequence_sha256": record["sequence_sha256"],
        "group": str((record.get("metadata") or {}).get("group") or
                     (record.get("metadata") or {}).get("lineage") or "ungrouped"),
        "coherent_product": bool(result.get("amplifies")), "signal_product": signal,
        "product_length": result.get("product"), "product_count": result.get("n_products"),
        "forward_identity": result.get("f_ident"), "reverse_identity": result.get("r_ident"),
        "probe_identity": result.get("probe_ident"), "probe_recognized": result.get("probe_binds"),
        "modeled_reason": reason or None, "affected_components": sorted(set(component_concerns)),
        "terminal_mismatch_concern": "3\u2032" in reason or "mismatch" in reason,
    }


def _scan_assay(assay, snapshot, model):
    rows = [_record_eval(assay, r, model) for r in snapshot.get("unique_records", [])]
    by_group = {}
    for group in sorted({x["group"] for x in rows}):
        gr = [x for x in rows if x["group"] == group]
        by_group[group] = {"n": len(gr), "coherent_products": sum(x["coherent_product"] for x in gr),
                           "signal_products": sum(x["signal_product"] for x in gr)}
    return {"assay_id": assay["assay_id"], "n_records": len(rows),
            "coherent_products": sum(x["coherent_product"] for x in rows),
            "signal_products": sum(x["signal_product"] for x in rows),
            "records": rows, "groups": by_group}


def _new_record_reasons(assay_id, role, baseline_eval, current_eval):
    old = {x["sequence_sha256"]: x for x in baseline_eval["records"]}
    reasons = []
    for cur in current_eval["records"]:
        if cur["sequence_sha256"] in old:
            continue
        if role == "target" and not cur["coherent_product"]:
            reasons.append({"code": "new_target_lost_coherent_product", "assay_id": assay_id,
                            "severity": "action_review_recommended", "record": cur,
                            "predicted_consequence": "possible target dropout under the declared product model"})
        elif role == "target" and cur["coherent_product"] and not cur["signal_product"]:
            reasons.append({"code": "new_target_lost_probe_recognition", "assay_id": assay_id,
                            "severity": "action_review_recommended", "record": cur,
                            "predicted_consequence": "possible probe-signal dropout despite a modeled primer product"})
        elif role == "off_target" and cur["signal_product"]:
            reasons.append({"code": "new_signal_generating_off_target", "assay_id": assay_id,
                            "severity": "action_review_recommended", "record": cur,
                            "predicted_consequence": "possible signal-generating off-target under the declared model"})
    return reasons


def scan_drift(assaysbom, baseline_target, current_target, *, baseline_offtarget=None,
               current_offtarget=None, model=None, scan_complete=True):
    sbom = build_assaysbom(assaysbom)
    for snap in [baseline_target, current_target] + [x for x in (baseline_offtarget, current_offtarget) if x]:
        if not validate_snapshot(snap)["valid"]:
            raise ValueError("all DriftGuard snapshots must have valid immutable hashes")
    if baseline_target.get("role") != "target" or current_target.get("role") != "target":
        raise ValueError("target snapshot roles are required")
    if bool(baseline_offtarget) != bool(current_offtarget):
        raise ValueError("baseline and current off-target snapshots must be supplied together")
    if baseline_offtarget and (baseline_offtarget.get("role") != "off_target" or current_offtarget.get("role") != "off_target"):
        raise ValueError("off-target snapshot roles are required")
    declared = {"max_mm": 5, "clamp_n": 2, "min_product": 40,
                "max_product": 3000, "min_probe_ident": 85.0}
    declared.update(dict(model or {}))
    assay_rows, reasons = [], []
    for assay in sbom["assays"]:
        bt = _scan_assay(assay, baseline_target, declared)
        ct = _scan_assay(assay, current_target, declared)
        reasons.extend(_new_record_reasons(assay["assay_id"], "target", bt, ct))
        row = {"assay_id": assay["assay_id"], "baseline_target": bt, "current_target": ct}
        if baseline_offtarget:
            bo = _scan_assay(assay, baseline_offtarget, declared)
            co = _scan_assay(assay, current_offtarget, declared)
            reasons.extend(_new_record_reasons(assay["assay_id"], "off_target", bo, co))
            row.update({"baseline_off_target": bo, "current_off_target": co})
        assay_rows.append(row)
    target_delta = snapshot_delta(baseline_target, current_target)
    off_delta = snapshot_delta(baseline_offtarget, current_offtarget) if baseline_offtarget else None
    if not scan_complete:
        state = "Scan incomplete"
    elif not current_target.get("unique_records"):
        state = "Evidence insufficient"
    elif any(x["code"] == "new_signal_generating_off_target" for x in reasons):
        state = "Possible signal-generating off-target"
    elif any(x["code"] in {"new_target_lost_coherent_product", "new_target_lost_probe_recognition"} for x in reasons):
        state = "Possible target dropout"
    elif target_delta["counts"]["added"] or (off_delta and off_delta["counts"]["added"]):
        state = "Stable with new variation"
    else:
        state = "Stable"
    body = {
        "schema_version": "oligoforge-drift-scan/v1", "driftguard_version": DRIFTGUARD_VERSION,
        "assaysbom_id": sbom["assaysbom_id"],
        "baseline_target_snapshot_id": baseline_target["snapshot_id"],
        "current_target_snapshot_id": current_target["snapshot_id"],
        "baseline_offtarget_snapshot_id": baseline_offtarget.get("snapshot_id") if baseline_offtarget else None,
        "current_offtarget_snapshot_id": current_offtarget.get("snapshot_id") if current_offtarget else None,
        "state": state, "action_review_recommended": bool(reasons),
        "scan_complete": bool(scan_complete), "declared_model": declared,
        "model_version": SPECIFICITY_MODEL_VERSION, "target_delta": target_delta,
        "offtarget_delta": off_delta, "assay_results": assay_rows, "reason_records": reasons,
        "search_completeness": {"status": "complete" if scan_complete else "incomplete",
                                "bounded_to_supplied_snapshots": True},
        "limitations": ["Sequence-level modeled evidence only; experimental confirmation required.",
                        "No numeric biological risk probability is calculated.",
                        "Snapshot sampling and metadata may be incomplete or unrepresentative."],
    }
    body["scan_sha256"] = sha256_value(body)
    body["scan_id"] = "ofscan_" + body["scan_sha256"][:24]
    return body
