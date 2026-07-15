"""RDML 1.3 assay-definition export for a qPCR panel.

The exporter creates a well-formed RDML ZIP containing target, dye, primer, probe and optional
amplification-efficiency metadata. It does not invent reference-gene roles or reporter dyes:
``target_type`` and ``dye`` should be supplied explicitly when known. Modified order notation is
reduced to the nucleotide backbone because RDML sequence fields do not encode vendor order syntax.
"""
import base64
import datetime
import io
import re
import xml.etree.ElementTree as ET
import zipfile

from . import thermo as T
from . import provenance as PROV

RDML_NS = "http://www.rdml.org"
XSI = "http://www.w3.org/2001/XMLSchema-instance"
SCHEMA = "http://www.rdml.org http://www.rdml.org/files/rdml/RDML_v1_3_REC.xsd"
_DYES = ["FAM", "SYBR", "VIC", "HEX", "JOE", "TET", "CY5", "CY5.5", "CY3", "ROX", "TAMRA",
         "TEX615", "ABY", "NED", "QUASAR670", "MAX", "TYE665"]


def _id_clean(s, fallback):
    s = re.sub(r"\s+", "_", str(s or "").strip())
    s = re.sub(r"[^A-Za-z0-9_.+-]", "", s)
    return s or fallback


def _dye_of(assay):
    explicit = str(assay.get("dye") or assay.get("fluorophore") or "").strip()
    if explicit:
        return _id_clean(explicit, "UNSPECIFIED")
    blob = str(assay.get("chemistry") or "").upper().replace(" ", "")
    for dye in _DYES:
        if dye.replace(".", "") in blob.replace(".", ""):
            return dye
    if not assay.get("probe") and any(x in blob for x in ("SYBR", "EVAGREEN", "INTERCALAT")):
        return "SYBR"
    return "UNSPECIFIED"


def _dye_chemistry(assay):
    blob = str(assay.get("chemistry") or "").upper()
    if assay.get("probe"):
        if "BEACON" in blob:
            return "hybridization probe"
        return "hydrolysis probe"
    if "EVA" in blob or "SATURATING" in blob:
        return "saturating DNA binding dye"
    if "SYBR" in blob or "INTERCALAT" in blob:
        return "non-saturating DNA binding dye"
    return None


def _ttype(assay):
    t = str(assay.get("target_type") or "").strip().lower()
    return t if t in ("ref", "toi") else "toi"


def _sub(parent, tag, text=None):
    e = ET.SubElement(parent, "{%s}%s" % (RDML_NS, tag))
    if text is not None:
        e.text = str(text)
    return e


def _plain_sequence(seq, field):
    if not seq:
        return ""
    bare, _notes, err = T.clean_seq(T.strip_mods(str(seq)))
    if err:
        raise ValueError("%s: %s" % (field, err))
    if not bare:
        raise ValueError("%s is empty after removing modification notation" % field)
    return bare.upper()



def _provenance_description(assay):
    """Compact, schema-safe provenance text for RDML's target description."""
    a = assay or {}
    m = a.get("ranker_manifest") or {}
    if not m:
        return "OligoForge rank provenance: not attached (manual or legacy assay)"
    v = PROV.verify_manifest(m)
    objective = a.get("objective_profile") or a.get("objective") or {}
    if isinstance(objective, dict):
        objective = objective.get("label") or objective.get("key") or "declared objective"
    rank = a.get("candidate_rank")
    return ("OligoForge %s; ranker %s; run %s; manifest_sha256 %s; manifest %s; "
            "candidate rank %s; objective %s" % (
                m.get("application_version", "unrecorded"),
                m.get("ranker_version", "unrecorded"),
                m.get("run_id", "unrecorded"),
                m.get("manifest_sha256", "unrecorded"),
                "verified" if v.get("valid") else "invalid_or_altered",
                rank if rank is not None else "unrecorded",
                objective or "unrecorded"))

def _oligo(seqparent, tag, seq, field):
    if not seq:
        return
    o = _sub(seqparent, tag)
    _sub(o, "sequence", _plain_sequence(seq, field))


def build(panel, meta=None):
    meta = meta or {}
    panel = list(panel or [])
    if len(panel) > 500:
        raise ValueError("RDML export is capped at 500 assays")
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    ET.register_namespace("", RDML_NS)
    ET.register_namespace("xsi", XSI)
    root = ET.Element("{%s}rdml" % RDML_NS,
                      {"version": "1.3", "{%s}schemaLocation" % XSI: SCHEMA})
    _sub(root, "dateMade", now)
    _sub(root, "dateUpdated", now)

    dyes, seen_ids, parsed = {}, set(), []
    for i, assay in enumerate(panel, 1):
        tid = _id_clean(assay.get("name") or assay.get("gene") or "assay_%d" % i, "assay_%d" % i)
        base, suffix = tid, 2
        while tid in seen_ids:
            tid = "%s_%d" % (base, suffix); suffix += 1
        seen_ids.add(tid)
        dye = _dye_of(assay)
        dyes.setdefault(dye, _dye_chemistry(assay))
        # Validate all sequence values before generating a partial document.
        for key in ("forward", "reverse", "probe"):
            if assay.get(key):
                _plain_sequence(assay[key], "%s %s" % (tid, key))
        parsed.append((tid, dye, assay))

    # RDML element order matters: all dye declarations precede all targets.
    for dye, chemistry in dyes.items():
        node = ET.SubElement(root, "{%s}dye" % RDML_NS, {"id": dye})
        if chemistry:
            _sub(node, "dyeChemistry", chemistry)

    for tid, dye, assay in parsed:
        tgt = ET.SubElement(root, "{%s}target" % RDML_NS, {"id": tid})
        desc_parts = [str(x) for x in (assay.get("gene"), assay.get("organism")) if x]
        desc_parts.append(_provenance_description(assay))
        _sub(tgt, "description", " / ".join(desc_parts))
        _sub(tgt, "type", _ttype(assay))
        val = assay.get("validation") or {}
        eff = val.get("efficiency_pct")
        if eff is not None:
            try:
                e = 1.0 + float(eff) / 100.0
                if 1.0 <= e <= 3.0:
                    _sub(tgt, "amplificationEfficiency", round(e, 4))
            except (TypeError, ValueError):
                pass
        amp_tm = assay.get("observed_amplicon_tm", val.get("observed_amplicon_tm"))
        if amp_tm is not None:
            try:
                tm = float(amp_tm)
                if 0 < tm < 120:
                    _sub(tgt, "meltingTemperature", round(tm, 2))
            except (TypeError, ValueError):
                pass
        ET.SubElement(tgt, "{%s}dyeId" % RDML_NS, {"id": dye})
        if assay.get("forward") or assay.get("reverse") or assay.get("probe"):
            seqs = _sub(tgt, "sequences")
            _oligo(seqs, "forwardPrimer", assay.get("forward"), "%s forward" % tid)
            _oligo(seqs, "reversePrimer", assay.get("reverse"), "%s reverse" % tid)
            _oligo(seqs, "probe1", assay.get("probe"), "%s probe" % tid)

    xml_bytes = b'<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(root, encoding="utf-8")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("rdml_data.xml", xml_bytes)
    return dict(xml=xml_bytes.decode("utf-8"), rdml_b64=base64.b64encode(buf.getvalue()).decode("ascii"),
                n_assays=len(panel), n_dyes=len(dyes), dyes=list(dyes), version="RDML 1.3",
                filename=_id_clean(meta.get("title", "oligoforge_panel"), "oligoforge_panel") + ".rdml")
