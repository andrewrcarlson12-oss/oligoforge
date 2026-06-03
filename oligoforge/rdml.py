"""RDML 1.2 export for a qPCR assay panel.

RDML (Real-time PCR Data Markup Language, https://www.rdml.org) is the standard XML interchange
format that qPCR analysis software (LinRegPCR, RDML-Ninja, Bio-Rad CFX, Roche LightCycler, ...)
imports. OligoForge's report is human-readable HTML + a flat CSV; RDML is the machine-readable
sibling, so a designed panel can be loaded straight into instrument/analysis software as predefined
targets instead of being retyped.

This emits an ASSAY-DEFINITION RDML: one <target> per assay carrying its forward/reverse primers,
probe, detection dye, and (when known) amplification efficiency. That is a valid RDML document
(samples/experiments/run data are all optional in the schema) and is exactly the part a user wants
to pre-populate before a run. A proper .rdml file is a ZIP containing the XML, so build() returns
both the raw XML (for inspection) and a base64 .rdml zip (for download).

HONEST SCOPE: the XML is built to the published RDML 1.2 element model and is verified well-formed
and structurally checked by tests/test_rdml.py. Full XSD-schema validation against rdml.org's XSD
is not performed in this build environment (no bundled schema / validator); confirm in your target
software on first import. Targets at this stage carry no Cq/run data -- add those after the run.
"""
import base64
import datetime
import io
import re
import xml.etree.ElementTree as ET
import zipfile

RDML_NS = "http://www.rdml.org"
XSI = "http://www.w3.org/2001/XMLSchema-instance"
SCHEMA = "http://www.rdml.org http://www.rdml.org/files/rdml/RDML_v1_2_REC.xsd"

# common reference-gene symbols -> RDML target type "ref"; everything else is "toi"
_REF_GENES = {"RPL13", "YWHAZ", "GAPDH", "ACTB", "B2M", "HMBS", "SDHA", "TBP", "HPRT1",
              "RPL4", "RPL19", "EEF1A1", "UBC", "PGK1", "18S", "28S", "RPS18", "TUBB"}
# fluorophore tokens we recognise in a dye/chemistry string
_DYES = ["FAM", "SYBR", "VIC", "HEX", "JOE", "TET", "Cy5", "Cy5.5", "Cy3", "ROX", "TAMRA",
         "TEX615", "TEX 615", "ABY", "NED", "Quasar 670", "Quasar670", "MAX", "TYE665", "TYE 665"]


def _id_clean(s, fallback):
    """RDML ids are free strings but must be unique and non-empty; keep them tidy."""
    s = re.sub(r"\s+", "_", (s or "").strip())
    s = re.sub(r"[^A-Za-z0-9_.+-]", "", s)
    return s or fallback


def _dye_of(assay):
    blob = " ".join(str(assay.get(k, "")) for k in ("dye", "chemistry", "fluorophore")).upper()
    for d in _DYES:
        if d.upper() in blob:
            return d.replace(" ", "")
    if assay.get("probe"):
        return "FAM"                 # probe assay, dye unstated -> FAM is the OligoForge default
    return "SYBR"                    # no probe -> intercalating dye


def _ttype(assay):
    t = (assay.get("target_type") or "").strip().lower()
    if t in ("ref", "toi"):
        return t
    g = (assay.get("gene") or "").upper().strip()
    return "ref" if g in _REF_GENES else "toi"


def _sub(parent, tag, text=None):
    e = ET.SubElement(parent, "{%s}%s" % (RDML_NS, tag))
    if text is not None:
        e.text = str(text)
    return e


def _oligo(seqparent, tag, seq):
    if not seq:
        return
    o = _sub(seqparent, tag)
    _sub(o, "sequence", str(seq).upper())


def build(panel, meta=None):
    meta = meta or {}
    panel = panel or []
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    ET.register_namespace("", RDML_NS)
    ET.register_namespace("xsi", XSI)
    root = ET.Element("{%s}rdml" % RDML_NS, {"version": "1.2", "{%s}schemaLocation" % XSI: SCHEMA})
    _sub(root, "dateMade", now)
    _sub(root, "dateUpdated", now)

    # element ORDER matters in RDML: all <dye> before all <target>
    dyes_used, seen_ids = [], set()
    parsed = []
    for i, a in enumerate(panel):
        tid = _id_clean(a.get("name") or a.get("gene") or "assay_%d" % (i + 1), "assay_%d" % (i + 1))
        base, n = tid, 2
        while tid in seen_ids:                       # guarantee unique target ids
            tid = "%s_%d" % (base, n); n += 1
        seen_ids.add(tid)
        dye = _dye_of(a)
        if dye not in dyes_used:
            dyes_used.append(dye)
        parsed.append((tid, dye, a))

    for d in dyes_used:
        ET.SubElement(root, "{%s}dye" % RDML_NS, {"id": d})

    for tid, dye, a in parsed:
        tgt = ET.SubElement(root, "{%s}target" % RDML_NS, {"id": tid})
        desc = " / ".join(x for x in (a.get("gene"), a.get("organism")) if x)
        if desc:
            _sub(tgt, "description", desc)
        _sub(tgt, "type", _ttype(a))                  # ref | toi
        val = a.get("validation") or {}
        eff = val.get("efficiency_pct")
        if eff is not None:
            try:
                _sub(tgt, "amplificationEfficiency", round(1.0 + float(eff) / 100.0, 4))  # RDML: E=2 is 100%
            except (TypeError, ValueError):
                pass
        ET.SubElement(tgt, "{%s}dyeId" % RDML_NS, {"id": dye})
        if a.get("forward") or a.get("reverse") or a.get("probe"):
            seqs = _sub(tgt, "sequences")
            _oligo(seqs, "forwardPrimer", a.get("forward"))
            _oligo(seqs, "reversePrimer", a.get("reverse"))
            _oligo(seqs, "probe1", a.get("probe"))

    xml_bytes = b'<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(root, encoding="utf-8")
    xml_text = xml_bytes.decode("utf-8")

    # .rdml is a zip containing the XML
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("rdml_data.xml", xml_bytes)
    rdml_b64 = base64.b64encode(buf.getvalue()).decode("ascii")

    return dict(xml=xml_text, rdml_b64=rdml_b64, n_assays=len(panel),
                n_dyes=len(dyes_used), dyes=dyes_used, version="RDML 1.2",
                filename=_id_clean(meta.get("title", "oligoforge_panel"), "oligoforge_panel") + ".rdml")
