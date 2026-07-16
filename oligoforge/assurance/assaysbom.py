"""Portable, deterministic molecular assay bill of materials (AssaySBOM)."""
from __future__ import annotations

import copy
import html
import json
import re

from .. import __version__ as APP_VERSION
from .. import profiles as PROFILES
from .. import thermo as T
from ..provenance import sha256_value, software_versions


ASSAYSBOM_VERSION = "1.0.0"
ASSAYSBOM_SCHEMA = "oligoforge-assaysbom/v1"
MAX_ASSAYS = 64
MAX_COMPONENTS = 512
_ROLES = {"forward_primer", "reverse_primer", "probe", "alternate_primer",
          "internal_control_primer", "internal_control_probe"}
_STATUSES = {"draft", "active", "under_review", "retired", "invalid"}


def _bare(value):
    return re.sub(r"[^ACGTRYSWKMBDHVN]", "", T.strip_mods(str(value or "")).upper().replace("U", "T"))


def _component(component, assay_id, index):
    if not isinstance(component, dict):
        raise ValueError("component %d in %s must be an object" % (index + 1, assay_id))
    role = str(component.get("role") or "").strip().lower().replace(" ", "_")
    role = {"forward": "forward_primer", "reverse": "reverse_primer"}.get(role, role)
    if role not in _ROLES:
        raise ValueError("component %d in %s has an unsupported role" % (index + 1, assay_id))
    order_sequence = str(component.get("order_sequence") or component.get("sequence") or "").strip()
    bare = _bare(order_sequence)
    if not bare:
        raise ValueError("component %d in %s has no usable nucleotide sequence" % (index + 1, assay_id))
    if len(bare) > 200:
        raise ValueError("component %d in %s exceeds 200 bases" % (index + 1, assay_id))
    component_id = str(component.get("component_id") or ("%s-%s-%d" % (assay_id, role, index + 1))).strip()
    return {
        "component_id": component_id, "role": role, "sequence": bare,
        "order_sequence": order_sequence or bare,
        "modifications_preserved": bool(order_sequence and order_sequence.upper() != bare),
        "locked_legacy_component": bool(component.get("locked_legacy_component") or component.get("locked")),
        "channel": component.get("channel"), "dye": component.get("dye"),
        "quencher": component.get("quencher"), "vendor": component.get("vendor"),
        "scale": component.get("scale"), "purification": component.get("purification"),
        "notes": str(component.get("notes") or "")[:2000],
    }


def _legacy_components(assay):
    rows = []
    for key, role in (("forward", "forward_primer"), ("reverse", "reverse_primer"), ("probe", "probe")):
        if assay.get(key):
            rows.append({"role": role, "sequence": assay[key]})
    return rows


def migrate_assaysbom(value):
    """Migrate the documented pre-v1 shape without guessing scientific fields."""
    if not isinstance(value, dict):
        raise ValueError("AssaySBOM must be a JSON object")
    src = copy.deepcopy(value)
    if src.get("schema_version") == ASSAYSBOM_SCHEMA:
        return src
    assays = src.get("assays")
    if not isinstance(assays, list):
        assays = [src]
    return {
        "schema_version": ASSAYSBOM_SCHEMA,
        "portfolio_name": src.get("portfolio_name") or src.get("name") or "Imported assay portfolio",
        "portfolio_version": str(src.get("portfolio_version") or src.get("version") or "1"),
        "status": src.get("status") or "draft", "review_state": src.get("review_state") or "unreviewed",
        "assays": assays, "interpretation_rules": src.get("interpretation_rules") or [],
        "existing_evidence": src.get("existing_evidence") or [],
        "repair_history": src.get("repair_history") or [],
        "migration": {"from_schema": src.get("schema_version") or "legacy-unversioned",
                      "lossless_for_declared_fields": True},
    }


def build_assaysbom(value):
    src = migrate_assaysbom(value)
    assays_in = src.get("assays")
    if not isinstance(assays_in, list) or not assays_in:
        raise ValueError("AssaySBOM requires at least one assay")
    if len(assays_in) > MAX_ASSAYS:
        raise ValueError("AssaySBOM exceeds %d assays" % MAX_ASSAYS)
    assays, total_components, ids = [], 0, set()
    for i, raw in enumerate(assays_in):
        if not isinstance(raw, dict):
            raise ValueError("assay %d must be an object" % (i + 1))
        assay_id = str(raw.get("assay_id") or raw.get("id") or ("assay-%d" % (i + 1))).strip()
        if not assay_id or assay_id in ids:
            raise ValueError("assay_id values must be non-empty and unique")
        ids.add(assay_id)
        chemistry = str(raw.get("chemistry") or raw.get("profile") or "unspecified")
        profile_key = raw.get("profile_key") or (chemistry if chemistry in PROFILES.PROFILES else None)
        components_in = raw.get("components") if isinstance(raw.get("components"), list) else _legacy_components(raw)
        components = [_component(x, assay_id, j) for j, x in enumerate(components_in)]
        total_components += len(components)
        roles = [x["role"] for x in components]
        if "forward_primer" not in roles or "reverse_primer" not in roles:
            raise ValueError("assay %s requires a forward and reverse primer" % assay_id)
        no_probe = "probe" not in roles
        assay_type = str(raw.get("assay_type") or ("sybr" if no_probe else "probe"))
        if assay_type == "sybr" and not no_probe:
            raise ValueError("assay %s is declared SYBR but contains a probe" % assay_id)
        assays.append({
            "assay_id": assay_id, "name": str(raw.get("name") or assay_id), "assay_type": assay_type,
            "chemistry": chemistry, "profile_key": profile_key,
            "components": components,
            "intended_target_groups": list(raw.get("intended_target_groups") or []),
            "near_neighbor_groups": list(raw.get("near_neighbor_groups") or []),
            "host_background_groups": list(raw.get("host_background_groups") or []),
            "multiplex_group": raw.get("multiplex_group"), "internal_control": bool(raw.get("internal_control")),
            "instrument": raw.get("instrument"), "channel_metadata": dict(raw.get("channel_metadata") or {}),
            "reaction_conditions": dict(raw.get("reaction_conditions") or src.get("reaction_conditions") or {}),
            "interpretation_rules": list(raw.get("interpretation_rules") or []),
            "existing_evidence": list(raw.get("existing_evidence") or []),
            "status": str(raw.get("status") or "draft"), "review_state": str(raw.get("review_state") or "unreviewed"),
        })
    if total_components > MAX_COMPONENTS:
        raise ValueError("AssaySBOM exceeds %d components" % MAX_COMPONENTS)
    status = str(src.get("status") or "draft")
    if status not in _STATUSES:
        raise ValueError("portfolio status is invalid")
    body = {
        "schema_version": ASSAYSBOM_SCHEMA, "assaysbom_version": ASSAYSBOM_VERSION,
        "application_version": APP_VERSION,
        "portfolio_name": str(src.get("portfolio_name") or "Assay portfolio"),
        "portfolio_version": str(src.get("portfolio_version") or "1"),
        "status": status, "review_state": str(src.get("review_state") or "unreviewed"),
        "assays": assays,
        "interpretation_rules": list(src.get("interpretation_rules") or []),
        "existing_evidence": list(src.get("existing_evidence") or []),
        "repair_history": list(src.get("repair_history") or []),
        "software_versions": software_versions(),
        "scope_statement": "Molecular and computational bill of materials; not proof of assay performance.",
        "migration": src.get("migration"),
    }
    body["content_sha256"] = sha256_value(body)
    body["assaysbom_id"] = "ofsbom_" + body["content_sha256"][:24]
    return body


def validate_assaysbom(value):
    try:
        normalized = build_assaysbom(value)
        return {"valid": True, "errors": [], "assaysbom": normalized}
    except (TypeError, ValueError) as exc:
        return {"valid": False, "errors": [{"path": "$", "message": str(exc), "remedy": "Correct the declared assay record and retry."}]}


def assaysbom_html(sbom):
    validated = validate_assaysbom(sbom)
    if not validated["valid"]:
        raise ValueError(validated["errors"][0]["message"])
    doc = validated["assaysbom"]
    rows = []
    for assay in doc["assays"]:
        for c in assay["components"]:
            rows.append("<tr><td>%s</td><td>%s</td><td>%s</td><td><code>%s</code></td><td>%s</td></tr>" %
                        tuple(html.escape(str(x)) for x in (assay["assay_id"], assay["assay_type"], c["role"],
                                                             c["order_sequence"], "yes" if c["locked_legacy_component"] else "no")))
    return """<!doctype html><meta charset=\"utf-8\"><title>OligoForge AssaySBOM</title>
<style>body{font:14px system-ui;max-width:1100px;margin:2rem auto;color:#18212f}table{border-collapse:collapse;width:100%%}th,td{border:1px solid #ccd3dd;padding:.45rem;text-align:left}code{word-break:break-all}.note{padding:1rem;background:#f2f5f8}</style>
<h1>%s</h1><p class=\"note\">%s</p><p><b>ID:</b> %s<br><b>Version:</b> %s<br><b>Status:</b> %s</p>
<table><thead><tr><th>Assay</th><th>Type</th><th>Component</th><th>Order sequence</th><th>Locked</th></tr></thead><tbody>%s</tbody></table>""" % (
        html.escape(doc["portfolio_name"]), html.escape(doc["scope_statement"]),
        html.escape(doc["assaysbom_id"]), html.escape(doc["portfolio_version"]),
        html.escape(doc["status"]), "".join(rows))


def assaysbom_json(sbom):
    return json.dumps(build_assaysbom(sbom), sort_keys=True, indent=2, ensure_ascii=False) + "\n"
