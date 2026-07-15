"""Self-contained assay-readiness report for a qPCR panel.

The report recomputes transparent oligo QC, distinguishes computational screening from empirical
validation, escapes all HTML fields, and neutralizes spreadsheet-formula injection in CSV exports.
"""
import csv
import datetime
import html as _h
import io
import math
import json

from . import thermo as T
from . import nn as NN
from . import __version__
from . import provenance as PROV

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



def _objective_label(value):
    if isinstance(value, dict):
        return value.get("label") or value.get("key") or value.get("name") or "declared objective"
    return str(value or "unrecorded")


def _provenance_fields(assay):
    """Return export-safe rank provenance without inventing missing metadata."""
    a = assay or {}
    manifest = a.get("ranker_manifest") or {}
    explanation = a.get("rank_explanation") or {}
    completeness = explanation.get("evidence_completeness") or {}
    if manifest:
        verification = PROV.verify_manifest(manifest)
        manifest_state = "verified" if verification.get("valid") else "invalid / altered"
    else:
        manifest_state = "not attached"
    return {
        "application_version": manifest.get("application_version") or a.get("oligoforge_version") or "unrecorded",
        "ranker_version": manifest.get("ranker_version") or a.get("ranker_version") or "unrecorded",
        "design_run_id": manifest.get("run_id") or a.get("design_run_id") or "unrecorded",
        "manifest_sha256": manifest.get("manifest_sha256") or a.get("manifest_sha256") or "",
        "manifest_state": manifest_state,
        "objective": _objective_label(a.get("objective_profile") or a.get("objective")),
        "candidate_rank": a.get("candidate_rank") if a.get("candidate_rank") is not None else "unrecorded",
        "preference_state": explanation.get("preference_state") or "unrecorded",
        "evidence_state": completeness.get("state") or "unrecorded",
        "rank_reversal_scenarios": explanation.get("rank_reversal_scenarios") or [],
        "rank_reversal_summary": explanation.get("ranking_may_reverse_if") or "unrecorded",
        "source_workflow": a.get("source_workflow") or "manual_or_legacy",
    }


def _condition_snapshot(assay):
    """Use the assay's recorded design conditions; fall back explicitly for legacy records."""
    m = (assay or {}).get("ranker_manifest") or {}
    raw = ((m.get("scientific_models") or {}).get("reaction_condition_snapshot") or {})
    keys = ("mv_conc_mM", "dv_conc_mM", "dntp_conc_mM", "total_oligo_conc_nM", "anneal_c")
    vals = [_finite(raw.get(k)) for k in keys]
    if all(v is not None for v in vals) and vals[0] >= 0 and vals[1] >= 0 and vals[2] >= 0 and vals[3] > 0:
        return tuple(vals), "recorded rank-manifest conditions"
    return T._snapshot(), "current-session fallback (no complete recorded conditions)"


def _condition_text(snap, source):
    return ("Mg2+ %.1f mM, monovalent %.0f mM, dNTP %.2f mM, total oligo %.0f nM, anneal %.1f °C; %s"
            % (snap[1], snap[0], snap[2], snap[3], snap[4], source))

def _oligo_qc(seq, snap=None):
    if not seq:
        return None
    raw = str(seq).strip()
    bare, _notes, err = T.clean_seq(T.strip_mods(raw))
    if err:
        return dict(seq=raw, error=err)
    snap = snap or T._snapshot()
    cond = {"mv_conc": snap[0], "dv_conc": snap[1], "dntp_conc": snap[2], "dna_conc": snap[3]}
    lna = NN.params_lna(raw, cond=cond, anneal_c=snap[4]) if "+" in raw else None
    tm = lna.get("tm") if lna else T._tm_acc_at(bare, snap)
    hp37, hp_an, hptm = T._hairpin_full_at(bare, snap[4], snap)
    sd37, sd_an, sdtm = T._self_dimer_full_at(bare, snap[4], snap)
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
        ("Versioned design provenance attached and self-verifying",
         bool(a.get("ranker_manifest")) and PROV.verify_manifest(a.get("ranker_manifest") or {}).get("valid", False)),
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
    fallback_snap = T._snapshot()
    fallback_salt = _condition_text(fallback_snap, "current session")

    blocks, csv_rows = [], []
    for a in panel:
        snap, condition_source = _condition_snapshot(a)
        condition_text = _condition_text(snap, condition_source)
        f, r, p = a.get("forward", ""), a.get("reverse", ""), a.get("probe")
        qf, qr, qp = _oligo_qc(f, snap), _oligo_qc(r, snap), _oligo_qc(p, snap)
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
        prov = _provenance_fields(a)

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
            "app_version": _h.escape(str(prov["application_version"])),
            "ranker_version": _h.escape(str(prov["ranker_version"])),
            "run_id": _h.escape(str(prov["design_run_id"])),
            "manifest_state": _h.escape(str(prov["manifest_state"])),
            "manifest_hash": _h.escape(str(prov["manifest_sha256"] or "unrecorded")),
            "objective": _h.escape(str(prov["objective"])),
            "candidate_rank": _h.escape(str(prov["candidate_rank"])),
            "preference": _h.escape(str(prov["preference_state"])),
            "evidence_state": _h.escape(str(prov["evidence_state"])),
            "rank_reversal_summary": _h.escape(str(prov["rank_reversal_summary"])),
            "source_workflow": _h.escape(str(prov["source_workflow"])),
            "condition_text": _h.escape(condition_text),
        }
        blocks.append(
            '<div class="assay"><h3>{name} <span class="meta">{gene}{org}</span></h3>'
            '<div class="kv"><b>Amplicon</b> {amp} nt &nbsp; <b>Chemistry</b> {chem} &nbsp; '
            '<b>Status</b> {status}</div>{tbl}'
            '<div class="kv"><b>Computational specificity</b> {spec} &nbsp; <b>gDNA/intron</b> {intr}</div>'
            '<div class="kv"><b>Efficiency</b> {eff} &nbsp; <b>R²</b> {r2} &nbsp; '
            '<b>Validated LOD95</b> {lod} &nbsp; <b>lowest fully detected standard (screen)</b> {lod_screen}</div>'
            '<div class="prov"><b>Design provenance</b> {app_version} &middot; ranker {ranker_version} &middot; '
            'run <span class="mono">{run_id}</span> &middot; manifest {manifest_state}<br>'
            '<b>Thermodynamic conditions</b> {condition_text}<br>'
            '<b>Rank context</b> rank {candidate_rank} &middot; {objective} &middot; {preference} &middot; {evidence_state} '
            '&middot; workflow {source_workflow}<br><b>Rank-reversal sensitivity</b> {rank_reversal_summary}<br>'
            '<span class="mono small">SHA-256 {manifest_hash}</span></div>'
            '<ul class="chk">{check}</ul></div>'.format(
                org=(" / " + fields["organism"] if fields["organism"] else ""), tbl=oligo_tbl,
                check=checklist, **fields))
        csv_rows.append([a.get("name", ""), a.get("gene", ""), a.get("organism", ""),
                         a.get("chemistry", ""), f, r, p or "", amp_display if amp_display is not None else "",
                         qf.get("tm", "") if qf else "", qr.get("tm", "") if qr else "", qp.get("tm", "") if qp else "",
                         spec_txt, intr_txt, eff if eff is not None else "", r2 if r2 is not None else "",
                         lod if lod is not None else "", lod_screen if lod_screen is not None else "",
                         prov["application_version"], prov["ranker_version"], prov["design_run_id"],
                         prov["manifest_sha256"], prov["manifest_state"], prov["objective"],
                         prov["candidate_rank"], prov["preference_state"], prov["evidence_state"],
                         prov["rank_reversal_summary"],
                         json.dumps(prov["rank_reversal_scenarios"], sort_keys=True, default=str),
                         prov["source_workflow"], condition_text])

    css = ("<style>body{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1080px;"
           "margin:30px auto;color:#111;padding:0 16px}h1{font-size:22px;margin:0 0 2px}h3{margin:16px 0 6px}"
           ".meta{font-weight:400;color:#666;font-size:14px}.assay{border:1px solid #ddd;border-radius:8px;"
           "padding:14px 16px;margin:14px 0}.kv{margin:4px 0}table.sub{border-collapse:collapse;width:100%;"
           "margin:8px 0}table.sub th,table.sub td{border-bottom:1px solid #eee;text-align:left;padding:4px 8px;"
           "font-size:12px}.mono{font-family:ui-monospace,Menlo,Consolas,monospace}ul.chk{columns:2;font-size:13px;"
           "color:#333;margin:8px 0}.prov{background:#f6f8fa;border:1px solid #e2e6ea;border-radius:6px;padding:8px;margin:8px 0;font-size:12px}.small{font-size:10px;word-break:break-all}.hdr{color:#666;font-size:12px;margin-bottom:14px}</style>")
    hdr = ('<div class="hdr">%s &middot; %s &middot; %s &middot; %d assay(s)<br>'
           "Thermodynamic estimates use each assay's recorded rank-manifest conditions when complete; legacy assays fall back to %s.<br>"
           'This is an assay-readiness summary, not MIQE certification. Modified-oligo predictions require vendor-specific confirmation; '
           'specificity, efficiency, product identity, and detection limits require empirical validation.</div>') % (
               _h.escape(str(meta.get("title", "OligoForge"))), now, TOOL_VERSION, len(panel), _h.escape(fallback_salt))
    body = "".join(blocks) if blocks else "<p>(empty panel)</p>"
    doc = ("<!doctype html><html><head><meta charset=utf-8><title>OligoForge qPCR panel report</title>"
           + css + "</head><body><h1>qPCR panel report</h1>" + hdr + body + "</body></html>")

    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(["name", "gene", "organism", "chemistry", "forward", "reverse", "probe", "amplicon_nt",
                "fwd_Tm", "rev_Tm", "probe_Tm", "specificity", "intron", "efficiency_pct", "R2",
                "validated_LOD95", "lowest_fully_detected_standard_screen",
                "oligoforge_version", "ranker_version", "design_run_id", "manifest_sha256",
                "manifest_state", "objective", "candidate_rank", "preference_state",
                "evidence_state", "rank_reversal_summary", "rank_reversal_scenarios_json",
                "source_workflow", "thermodynamic_conditions"])
    w.writerows([[_csv_cell(v) for v in row] for row in csv_rows])
    return dict(html=doc, csv=buf.getvalue(), n_assays=len(panel), tool=TOOL_VERSION, generated=now)
