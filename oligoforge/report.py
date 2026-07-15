"""Self-contained assay-readiness report for a qPCR panel.

The report recomputes transparent oligo QC, distinguishes computational screening from empirical
validation, escapes all HTML fields, and neutralizes spreadsheet-formula injection in CSV exports.
"""
import csv
import datetime
import html as _h
import io
import math

from . import thermo as T
from . import nn as NN
from . import __version__

TOOL_VERSION = "OligoForge v%s" % __version__


def _finite(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _csv_cell(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value
    s = str(value if value is not None else "")
    return "'" + s if s[:1] in ("=", "+", "-", "@") else s


def _oligo_qc(seq):
    if not seq:
        return None
    raw = str(seq).strip()
    bare, _notes, err = T.clean_seq(T.strip_mods(raw))
    if err:
        return dict(seq=raw, error=err)
    lna = NN.params_lna(raw) if "+" in raw else None
    tm = lna.get("tm") if lna else T.tm_acc(bare)
    hp37, hp_an, hptm = T.hairpin_full(bare)
    sd37, sd_an, sdtm = T.self_dimer_full(bare)
    return dict(seq=raw, backbone=bare, length=len(bare), tm=round(tm, 1),
                tm_basis=("McTigue LNA estimate" if lna else "DNA nearest-neighbor estimate"),
                gc=round(T.gc_percent(bare), 1), hairpin=round(hp37, 2), hairpin_anneal=round(hp_an, 2),
                hairpin_tm=round(hptm, 1), self_dimer=round(sd37, 2),
                self_dimer_anneal=round(sd_an, 2), self_dimer_tm=round(sdtm, 1))


def _checklist(a):
    chk = a.get("checks") or {}
    val = a.get("validation") or {}
    empirical_product = val.get("observed_peaks") is not None or val.get("product_identity_method")
    return [
        ("Target identity recorded", bool(a.get("gene") and a.get("organism"))),
        ("Primer sequences recorded", bool(a.get("forward") and a.get("reverse"))),
        ("Detection chemistry recorded", bool(a.get("probe")) or "SYBR" in str(a.get("chemistry") or "").upper()),
        ("Amplicon length/sequence recorded", bool(a.get("amplicon"))),
        ("Computational specificity screen attached", bool(chk.get("specificity"))),
        ("gDNA/exon-junction assessment attached", bool(chk.get("intron"))),
        ("Empirical PCR efficiency attached", val.get("efficiency_pct") is not None),
        ("Empirical linearity attached", val.get("r2") is not None),
        ("Empirical product-identity evidence attached", bool(empirical_product)),
        ("Replicated detection-limit study attached", val.get("validated_lod95") is not None),
    ]


def _orow(label, q):
    if not q:
        return ""
    if q.get("error"):
        return '<tr><td>%s</td><td class="mono">%s</td><td colspan="7">invalid: %s</td></tr>' % (
            _h.escape(label), _h.escape(q.get("seq", "")), _h.escape(q["error"]))
    vals = [label, q["seq"], q["length"], q["tm"], q["gc"], q["hairpin"],
            q["hairpin_anneal"], q["self_dimer"], q["self_dimer_anneal"]]
    esc = [_h.escape(str(v)) for v in vals]
    return ('<tr><td>%s</td><td class="mono">%s</td><td>%s</td><td>%s</td><td>%s</td>'
            '<td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>') % tuple(esc)


def build(panel, meta=None):
    panel = list(panel or [])
    if len(panel) > 500:
        raise ValueError("report export is capped at 500 assays")
    meta = meta or {}
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    c = dict(T.COND)
    salt = ("Mg2+ %.1f mM, monovalent %.0f mM, dNTP %.2f mM, total oligo %.0f nM"
            % (c["dv_conc"], c["mv_conc"], c["dntp_conc"], c["dna_conc"]))

    blocks, csv_rows = [], []
    for a in panel:
        f, r, p = a.get("forward", ""), a.get("reverse", ""), a.get("probe")
        qf, qr, qp = _oligo_qc(f), _oligo_qc(r), _oligo_qc(p)
        chk, val = a.get("checks") or {}, a.get("validation") or {}
        spec = chk.get("specificity")
        spec_txt = ("%s predicted product(s)" % spec.get("n_products")
                    if isinstance(spec, dict) and "n_products" in spec
                    else ("screen attached" if spec else "not run"))
        intr = chk.get("intron")
        intr_txt = (str(intr.get("verdict")) if isinstance(intr, dict) and intr.get("verdict")
                    else ("assessment attached" if intr else "not run"))
        eff, r2 = _finite(val.get("efficiency_pct")), _finite(val.get("r2"))
        lod = val.get("validated_lod95")
        lod_screen = val.get("lowest_fully_detected_standard", val.get("lod_practical", val.get("lod")))
        amp = a.get("amplicon")
        amp_display = len(amp) if isinstance(amp, str) else amp

        oligo_tbl = ("<table class=sub><tr><th></th><th>sequence/order notation (5'→3')</th><th>nt</th>"
                     "<th>Tm °C</th><th>GC%</th><th>hairpin ΔG37</th><th>hairpin ΔG@anneal</th>"
                     "<th>self-dimer ΔG37</th><th>self-dimer ΔG@anneal</th></tr>"
                     + _orow("Forward", qf) + _orow("Reverse", qr) + (_orow("Probe", qp) if qp else "")
                     + "</table>")
        checklist = "".join("<li>%s %s</li>" % ("&#10003;" if ok else "&#9744;", _h.escape(lbl))
                            for lbl, ok in _checklist(a))
        fields = {
            "name": _h.escape(str(a.get("name") or "(unnamed)")),
            "gene": _h.escape(str(a.get("gene") or "")),
            "organism": _h.escape(str(a.get("organism") or "")),
            "amp": _h.escape(str(amp_display if amp_display is not None else "—")),
            "chem": _h.escape(str(a.get("chemistry") or "—")),
            "status": _h.escape(str(a.get("status") or "designed")),
            "spec": _h.escape(spec_txt), "intr": _h.escape(intr_txt),
            "eff": (("%g" % eff) + "%" if eff is not None else "—"),
            "r2": ("%.4f" % r2 if r2 is not None else "—"),
            "lod": (_h.escape(str(lod)) if lod is not None else "—"),
            "lod_screen": (_h.escape(str(lod_screen)) if lod_screen is not None else "—"),
        }
        blocks.append(
            '<div class="assay"><h3>{name} <span class="meta">{gene}{org}</span></h3>'
            '<div class="kv"><b>Amplicon</b> {amp} nt &nbsp; <b>Chemistry</b> {chem} &nbsp; '
            '<b>Status</b> {status}</div>{tbl}'
            '<div class="kv"><b>Computational specificity</b> {spec} &nbsp; <b>gDNA/intron</b> {intr}</div>'
            '<div class="kv"><b>Efficiency</b> {eff} &nbsp; <b>R²</b> {r2} &nbsp; '
            '<b>Validated LOD95</b> {lod} &nbsp; <b>lowest fully detected standard (screen)</b> {lod_screen}</div>'
            '<ul class="chk">{check}</ul></div>'.format(
                org=(" / " + fields["organism"] if fields["organism"] else ""), tbl=oligo_tbl,
                check=checklist, **fields))
        csv_rows.append([a.get("name", ""), a.get("gene", ""), a.get("organism", ""),
                         a.get("chemistry", ""), f, r, p or "", amp_display if amp_display is not None else "",
                         qf.get("tm", "") if qf else "", qr.get("tm", "") if qr else "", qp.get("tm", "") if qp else "",
                         spec_txt, intr_txt, eff if eff is not None else "", r2 if r2 is not None else "",
                         lod if lod is not None else "", lod_screen if lod_screen is not None else ""])

    css = ("<style>body{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1080px;"
           "margin:30px auto;color:#111;padding:0 16px}h1{font-size:22px;margin:0 0 2px}h3{margin:16px 0 6px}"
           ".meta{font-weight:400;color:#666;font-size:14px}.assay{border:1px solid #ddd;border-radius:8px;"
           "padding:14px 16px;margin:14px 0}.kv{margin:4px 0}table.sub{border-collapse:collapse;width:100%;"
           "margin:8px 0}table.sub th,table.sub td{border-bottom:1px solid #eee;text-align:left;padding:4px 8px;"
           "font-size:12px}.mono{font-family:ui-monospace,Menlo,Consolas,monospace}ul.chk{columns:2;font-size:13px;"
           "color:#333;margin:8px 0}.hdr{color:#666;font-size:12px;margin-bottom:14px}</style>")
    hdr = ('<div class="hdr">%s &middot; %s &middot; %s &middot; thermodynamic estimates at %s &middot; '
           'anneal %.0f&deg;C &middot; %d assay(s)<br>This is an assay-readiness summary, not MIQE certification. '
           'Modified-oligo predictions require vendor-specific confirmation; specificity, efficiency, product identity, '
           'and detection limits require empirical validation.</div>') % (
               _h.escape(str(meta.get("title", "OligoForge"))), now, TOOL_VERSION, _h.escape(salt), T.ANNEAL_C, len(panel))
    body = "".join(blocks) if blocks else "<p>(empty panel)</p>"
    doc = ("<!doctype html><html><head><meta charset=utf-8><title>OligoForge qPCR panel report</title>"
           + css + "</head><body><h1>qPCR panel report</h1>" + hdr + body + "</body></html>")

    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(["name", "gene", "organism", "chemistry", "forward", "reverse", "probe", "amplicon_nt",
                "fwd_Tm", "rev_Tm", "probe_Tm", "specificity", "intron", "efficiency_pct", "R2",
                "validated_LOD95", "lowest_fully_detected_standard_screen"])
    w.writerows([[_csv_cell(v) for v in row] for row in csv_rows])
    return dict(html=doc, csv=buf.getvalue(), n_assays=len(panel), tool=TOOL_VERSION, generated=now)
