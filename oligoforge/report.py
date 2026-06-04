"""Self-contained MIQE-style report for a qPCR assay panel.

Builds one HTML document (inline CSS, no external assets) summarizing each assay's
design, recomputed oligo QC, specificity/intron checks, and validation (efficiency,
R^2, dynamic range, LOD), plus a MIQE 2.0 completeness checklist. Also returns a flat
CSV of the key reporting fields.
"""
import csv
import datetime
import html as _h
import io

from . import thermo as T

TOOL_VERSION = "OligoForge v1.21.4"


def _oligo_qc(seq):
    if not seq:
        return None
    s = seq.upper()
    return dict(seq=s, length=len(s), tm=round(T.tm_acc(s), 1), gc=round(T.gc_percent(s), 1),
                hairpin=round(T.hairpin(s)[0], 2), self_dimer=round(T.self_dimer(s), 2))


def _checklist(a):
    chk = a.get("checks") or {}
    val = a.get("validation") or {}
    return [
        ("Target gene & organism", bool(a.get("gene") and a.get("organism"))),
        ("Primer sequences + Tm", bool(a.get("forward") and a.get("reverse"))),
        ("Probe / detection chemistry", bool(a.get("probe")) or "SYBR" in (a.get("chemistry") or "").upper()),
        ("Amplicon length", bool(a.get("amplicon"))),
        ("Specificity (in-silico PCR / BLAST)", bool(chk.get("specificity"))),
        ("Exon/intron location (gDNA)", bool(chk.get("intron"))),
        ("PCR efficiency", val.get("efficiency_pct") is not None),
        ("R^2 / dynamic range / LOD", val.get("r2") is not None),
    ]


def _orow(label, q):
    if not q:
        return ""
    return ("<tr><td>%s</td><td class=mono>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
            % (label, q["seq"], q["length"], q["tm"], q["gc"], q["hairpin"]))


def build(panel, meta=None):
    meta = meta or {}
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    c = T.COND
    salt = ("Mg2+ %.1f mM, monovalent %.0f mM, dNTP %.2f mM, oligo %.0f nM"
            % (c["dv_conc"], c["mv_conc"], c["dntp_conc"], c["dna_conc"]))

    blocks, csv_rows = [], []
    for a in panel:
        f, r, p = a.get("forward", ""), a.get("reverse", ""), a.get("probe")
        qf, qr, qp = _oligo_qc(f), _oligo_qc(r), _oligo_qc(p)
        chk, val = a.get("checks") or {}, a.get("validation") or {}
        spec = chk.get("specificity")
        spec_txt = ("%s predicted product(s)" % spec.get("n_products")
                    if isinstance(spec, dict) and "n_products" in spec
                    else ("checked" if spec else "not run"))
        intr = chk.get("intron")
        intr_txt = (_h.escape(str(intr.get("verdict"))) if isinstance(intr, dict) and intr.get("verdict")
                    else ("checked" if intr else "not run"))
        eff, r2 = val.get("efficiency_pct"), val.get("r2")
        lod = val.get("lod_practical") if val.get("lod_practical") is not None else val.get("lod")

        oligo_tbl = ("<table class=sub><tr><th></th><th>sequence (5'-&gt;3')</th><th>nt</th>"
                     "<th>Tm</th><th>GC%</th><th>hairpin dG</th></tr>"
                     + _orow("Forward", qf) + _orow("Reverse", qr) + (_orow("Probe", qp) if qp else "")
                     + "</table>")
        checklist = "".join("<li>%s %s</li>" % ("&#10003;" if ok else "&#9744;", _h.escape(lbl))
                            for lbl, ok in _checklist(a))
        blocks.append(
            "<div class=assay><h3>%s <span class=meta>%s%s</span></h3>"
            "<div class=kv><b>Amplicon</b> %s nt &nbsp; <b>Chemistry</b> %s &nbsp; <b>Status</b> %s</div>"
            "%s"
            "<div class=kv><b>Specificity</b> %s &nbsp; <b>gDNA / intron</b> %s</div>"
            "<div class=kv><b>Efficiency</b> %s &nbsp; <b>R\u00b2</b> %s &nbsp; <b>LOD</b> %s</div>"
            "<ul class=chk>%s</ul></div>"
            % (_h.escape(a.get("name") or "(unnamed)"),
               _h.escape(a.get("gene") or ""),
               (" / " + _h.escape(a.get("organism")) if a.get("organism") else ""),
               a.get("amplicon") or "\u2014", _h.escape(a.get("chemistry") or "") or "\u2014",
               _h.escape(a.get("status") or "designed"),
               oligo_tbl, spec_txt, intr_txt,
               ("%s%%" % eff if eff is not None else "\u2014"),
               (r2 if r2 is not None else "\u2014"),
               (lod if lod is not None else "\u2014"),
               checklist))
        csv_rows.append([a.get("name", ""), a.get("gene", ""), a.get("organism", ""),
                         a.get("chemistry", ""), f, r, p or "", a.get("amplicon", ""),
                         qf["tm"] if qf else "", qr["tm"] if qr else "", qp["tm"] if qp else "",
                         spec_txt, intr_txt, eff if eff is not None else "",
                         r2 if r2 is not None else "", lod if lod is not None else ""])

    css = ("<style>body{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:920px;"
           "margin:30px auto;color:#111;padding:0 16px}h1{font-size:22px;margin:0 0 2px}h3{margin:16px 0 6px}"
           ".meta{font-weight:400;color:#666;font-size:14px}.assay{border:1px solid #ddd;border-radius:8px;"
           "padding:14px 16px;margin:14px 0}.kv{margin:4px 0}table.sub{border-collapse:collapse;width:100%;"
           "margin:8px 0}table.sub th,table.sub td{border-bottom:1px solid #eee;text-align:left;padding:4px 8px;"
           "font-size:13px}.mono{font-family:ui-monospace,Menlo,Consolas,monospace}ul.chk{columns:2;font-size:13px;"
           "color:#333;margin:8px 0}.hdr{color:#666;font-size:12px;margin-bottom:14px}</style>")
    hdr = ("<div class=hdr>" + _h.escape(meta.get("title", "OligoForge")) + " &middot; " + now + " &middot; "
           + TOOL_VERSION + " &middot; thermodynamics: SantaLucia 1998 + Owczarzy at " + _h.escape(salt)
           + " &middot; structure &Delta;G is reported at 37&deg;C (vendor-comparable) and evaluated at the "
           + ("%.0f" % T.ANNEAL_C) + "&deg;C anneal &mdash; a hairpin/dimer whose melting Tm is below the "
           "anneal temperature is largely absent during priming"
           + " &middot; " + str(len(panel)) + " assay(s)<br>MIQE-style summary for assembly and review; "
           "confirm final oligos in your vendor tool, and cross-check reference-gene choice in "
           "NormFinder / RefFinder.</div>")
    body = "".join(blocks) if blocks else "<p>(empty panel)</p>"
    doc = ("<!doctype html><html><head><meta charset=utf-8><title>OligoForge qPCR panel report</title>"
           + css + "</head><body><h1>qPCR panel report</h1>" + hdr + body + "</body></html>")

    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["name", "gene", "organism", "chemistry", "forward", "reverse", "probe", "amplicon_nt",
                "fwd_Tm", "rev_Tm", "probe_Tm", "specificity", "intron", "efficiency_pct", "R2", "LOD"])
    w.writerows(csv_rows)
    return dict(html=doc, csv=buf.getvalue(), n_assays=len(panel), tool=TOOL_VERSION, generated=now)
