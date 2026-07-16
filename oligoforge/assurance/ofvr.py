"""Versioned, reason-coded OligoForge Molecular Vulnerability Records."""
from __future__ import annotations

import datetime

from ..provenance import sha256_value


OFVR_VERSION = "1.0.0"
OFVR_SCHEMA = "oligoforge-ofvr/v1"


def generate_ofvrs(scan, *, issuance_year=None):
    if not isinstance(scan, dict) or scan.get("schema_version") != "oligoforge-drift-scan/v1":
        raise ValueError("a DriftGuard scan is required")
    year = int(issuance_year or datetime.datetime.now(datetime.timezone.utc).year)
    if year < 2000 or year > 9999:
        raise ValueError("issuance_year must be four digits")
    records, seen = [], set()
    for reason in scan.get("reason_records", []):
        record = reason.get("record") or {}
        fingerprint_body = {"scan_id": scan.get("scan_id"), "assay_id": reason.get("assay_id"),
                            "code": reason.get("code"), "sequence_sha256": record.get("sequence_sha256"),
                            "affected_components": record.get("affected_components") or []}
        fingerprint = sha256_value(fingerprint_body)
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        body = {
            "schema_version": OFVR_SCHEMA, "ofvr_version": OFVR_VERSION,
            "assaysbom_id": scan.get("assaysbom_id"), "assay_id": reason.get("assay_id"),
            "detecting_scan_id": scan.get("scan_id"),
            "baseline_snapshot_ids": [x for x in [scan.get("baseline_target_snapshot_id"),
                                                   scan.get("baseline_offtarget_snapshot_id")] if x],
            "current_snapshot_ids": [x for x in [scan.get("current_target_snapshot_id"),
                                                  scan.get("current_offtarget_snapshot_id")] if x],
            "reason_code": reason.get("code"), "affected_components": record.get("affected_components") or [],
            "affected_records": [{"record_id": record.get("record_id"),
                                  "sequence_sha256": record.get("sequence_sha256"),
                                  "group": record.get("group")}],
            "exact_modeled_observation": {k: record.get(k) for k in
                                           ("coherent_product", "signal_product", "product_length",
                                            "forward_identity", "reverse_identity", "probe_identity",
                                            "probe_recognized", "terminal_mismatch_concern", "modeled_reason")},
            "predicted_consequence": reason.get("predicted_consequence"),
            "evidence_strength": "computationally_supported_within_declared_model",
            "severity_rationale": reason.get("severity"),
            "search_completeness": scan.get("search_completeness"),
            "model_versions": {"driftguard": scan.get("driftguard_version"),
                               "product_model": scan.get("model_version")},
            "limitations": list(scan.get("limitations") or []),
            "review_status": "unreviewed", "linked_repairs": [], "disposition_history": [],
            "fingerprint_sha256": fingerprint,
            "standard_status": "OligoForge-local record; not a recognized external vulnerability standard",
        }
        body["record_sha256"] = sha256_value(body)
        body["ofvr_id"] = "OFVR-%04d-%s" % (year, fingerprint[:12].upper())
        records.append(body)
    return records
