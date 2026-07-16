"""Deterministic machine-readable and human-readable Assurance evidence bundles."""
from __future__ import annotations

import html
import json

from ..provenance import sha256_value


EVIDENCE_PACKAGE_VERSION = "1.0.0"


def build_evidence_package(*, assaysbom, snapshots, deltas=None, drift_scans=None,
                           vulnerabilities=None, validation_plans=None, repairs=None):
    artifacts = {
        "assaysbom": assaysbom, "snapshots": list(snapshots or []), "deltas": list(deltas or []),
        "drift_scans": list(drift_scans or []), "vulnerabilities": list(vulnerabilities or []),
        "validation_plans": list(validation_plans or []), "repairs": list(repairs or []),
    }
    manifest_rows = []
    for category, value in artifacts.items():
        rows = value if isinstance(value, list) else [value]
        for index, row in enumerate(rows):
            if row is not None:
                manifest_rows.append({"category": category, "index": index,
                                      "sha256": sha256_value(row),
                                      "schema_version": row.get("schema_version") if isinstance(row, dict) else None})
    body = {
        "schema_version": "oligoforge-assurance-evidence-package/v1",
        "evidence_package_version": EVIDENCE_PACKAGE_VERSION,
        "artifacts": artifacts, "manifest": manifest_rows,
        "scope_statement": "Reproducible computational evidence package; expert and laboratory confirmation required.",
    }
    body["package_sha256"] = sha256_value(body)
    body["package_id"] = "ofpkg_" + body["package_sha256"][:24]
    return body


def verify_evidence_package(package):
    if not isinstance(package, dict):
        return {"valid": False, "errors": ["package must be an object"]}
    supplied = package.get("package_sha256")
    body = {k: v for k, v in package.items() if k not in {"package_sha256", "package_id"}}
    digest = sha256_value(body)
    rows = []
    artifacts = package.get("artifacts") or {}
    for item in package.get("manifest") or []:
        value = artifacts.get(item.get("category"))
        values = value if isinstance(value, list) else [value]
        try:
            actual = sha256_value(values[int(item.get("index", 0))])
        except (IndexError, TypeError, ValueError):
            actual = None
        rows.append({"item": item, "valid": actual == item.get("sha256"), "calculated_sha256": actual})
    valid = supplied == digest and package.get("package_id") == "ofpkg_" + digest[:24] and all(x["valid"] for x in rows)
    return {"valid": valid, "package_digest_valid": supplied == digest, "artifact_checks": rows}


def evidence_package_html(package):
    checked = verify_evidence_package(package)
    rows = "".join("<tr><td>%s</td><td>%s</td><td><code>%s</code></td></tr>" %
                   (html.escape(str(x["category"])), html.escape(str(x["index"])), html.escape(str(x["sha256"])))
                   for x in package.get("manifest", []))
    state = "verified" if checked["valid"] else "verification failed"
    return """<!doctype html><meta charset=\"utf-8\"><title>OligoForge Assurance evidence</title>
<style>body{font:14px system-ui;max-width:1100px;margin:2rem auto;color:#17212d}table{border-collapse:collapse;width:100%%}th,td{border:1px solid #ccd3dd;padding:.5rem;text-align:left}code{word-break:break-all}.note{background:#f2f5f8;padding:1rem}</style>
<h1>OligoForge Assurance evidence package</h1><p class=\"note\">%s</p><p><b>Package:</b> %s<br><b>Verification:</b> %s</p>
<table><thead><tr><th>Artifact type</th><th>Index</th><th>SHA-256</th></tr></thead><tbody>%s</tbody></table>""" % (
        html.escape(str(package.get("scope_statement") or "")), html.escape(str(package.get("package_id") or "")),
        html.escape(state), rows)


def package_json(package):
    return json.dumps(package, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
