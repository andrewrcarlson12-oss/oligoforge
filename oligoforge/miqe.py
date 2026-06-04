"""Consolidate an assay's identity, engine-predicted thermodynamics, and (optional) empirical run data
into a single MIQE-aligned, provenance-pinned validation record.

This is the join between the physical-chemistry engine (oligoforge/nn.py — predicted Tm, ΔG, fraction
bound) and the empirical analyses (oligoforge/quant.standard_curve, melt) that the rest of the app
already produces. It does not re-derive any thermodynamics or statistics; it calls the validated
primitives, applies MIQE acceptance criteria (Bustin et al. 2009), and stamps a reproducible provenance
block (sequence checksums, tool version, UTC timestamp, reaction conditions) so a record can be cited.
"""
import hashlib
from datetime import datetime, timezone
from . import thermo as T
from . import nn as NN
from .__init__ import __version__ as _VERSION


def _csum(seq):
    s = T.strip_mods(seq or "").upper()
    return ("sha256:" + hashlib.sha256(s.encode()).hexdigest()[:12]) if s else None


def _predict(seq, role):
    """Engine prediction for one oligo: LNA-aware Tm/fraction-bound if it carries '+' bases, else the
    plain NN readout. Returns a compact dict or None."""
    if not seq:
        return None
    if "+" in seq:
        p = NN.params_lna(seq)
        if p:
            return dict(role=role, tm=p["tm"], dg37=p["dg37"], frac_bound=p["frac_bound"],
                        lna_n=p["lna_n"], basis="McTigue LNA")
    p = NN.params(T.strip_mods(seq))
    if p:
        return dict(role=role, tm=p["tm"], dg37=p["dg37"], frac_bound=p["frac_bound"], basis="SantaLucia NN")
    return None


def _chk(cid, label, status, detail=""):
    return dict(id=cid, label=label, status=status, detail=detail)


def validate_assay(assay, cond=None, anneal_c=None, observed_amplicon_tm=None, observed_peaks=None):
    """Build a MIQE validation record for one assay dict (forward, reverse, probe, amplicon[/_tm],
    name, chemistry, dye, validation={efficiency_pct,r2,slope,lod}). Optional observed_amplicon_tm and
    observed_peaks (e.g. from a melt run) drive the predicted-vs-observed and specificity checks.
    Conditions default to the live reaction conditions. Returns the record (also carries a Markdown
    rendering under 'markdown')."""
    a = assay or {}
    F, Rv, P = a.get("forward", ""), a.get("reverse", ""), a.get("probe", "")
    cond = dict(T.COND, **(cond or {}))
    ac = T.ANNEAL_C if anneal_c is None else float(anneal_c)
    val = a.get("validation") or {}
    if observed_amplicon_tm is None:
        observed_amplicon_tm = val.get("observed_amplicon_tm")

    # predicted thermodynamics (engine)
    pred = {k: v for k, v in (("forward", _predict(F, "forward")), ("reverse", _predict(Rv, "reverse")),
                              ("probe", _predict(P, "probe"))) if v}
    amp = a.get("amplicon")
    pred_amp_tm = a.get("amplicon_tm")
    if pred_amp_tm is None and isinstance(amp, str) and len(amp) >= 20:
        try:
            pred_amp_tm = round(T.amplicon_tm(amp), 1)
        except Exception:
            pred_amp_tm = None

    checks = []
    # --- empirical amplification metrics (from an attached standard curve) ---
    eff = val.get("efficiency_pct")
    if eff is not None:
        st = "pass" if 90.0 <= eff <= 110.0 else ("warn" if 85.0 <= eff <= 115.0 else "fail")
        checks.append(_chk("efficiency", "Amplification efficiency 90-110%", st, "%.1f%%" % eff))
    else:
        checks.append(_chk("efficiency", "Amplification efficiency 90-110%", "na", "no standard curve attached"))
    r2 = val.get("r2")
    if r2 is not None:
        st = "pass" if r2 >= 0.98 else ("warn" if r2 >= 0.95 else "fail")
        checks.append(_chk("linearity", "Calibration linearity R² ≥ 0.98", st, "%.4f" % r2))
    else:
        checks.append(_chk("linearity", "Calibration linearity R² ≥ 0.98", "na", "no standard curve attached"))
    slope = val.get("slope")
    if slope is not None:
        st = "pass" if -3.58 <= slope <= -3.10 else "warn"
        checks.append(_chk("slope", "Calibration slope −3.58 to −3.10", st, "%.3f" % slope))
    lod = val.get("lod")
    if lod is not None:
        checks.append(_chk("lod", "Limit of detection reported", "pass", "lowest fully-detected standard: %s" % lod))
    else:
        checks.append(_chk("lod", "Limit of detection reported", "na", "no LOD on file"))

    # --- specificity (observed melt) ---
    if observed_peaks is not None:
        st = "pass" if observed_peaks == 1 else "fail"
        checks.append(_chk("single_product", "Single specific product (one melt peak)", st, "%d peak(s) observed" % observed_peaks))
    else:
        checks.append(_chk("single_product", "Single specific product (one melt peak)", "na", "no melt curve provided"))

    # --- predicted vs observed amplicon Tm (closes the loop on the Tm engine) ---
    if observed_amplicon_tm is not None and pred_amp_tm is not None:
        d = abs(observed_amplicon_tm - pred_amp_tm)
        st = "pass" if d <= 3.0 else ("warn" if d <= 5.0 else "fail")
        checks.append(_chk("amplicon_tm", "Observed amplicon Tm matches prediction (±3 °C)", st,
                           "observed %.1f vs predicted %.1f °C (Δ %.1f)" % (observed_amplicon_tm, pred_amp_tm, observed_amplicon_tm - pred_amp_tm)))
    elif pred_amp_tm is not None:
        checks.append(_chk("amplicon_tm", "Observed amplicon Tm matches prediction (±3 °C)", "na",
                           "predicted %.1f °C; no observed melt Tm attached" % pred_amp_tm))

    # --- engine binding checks at the annealing temperature ---
    if "probe" in pred:
        fb = pred["probe"]["frac_bound"]
        st = "pass" if fb is not None and fb >= 0.5 else ("warn" if fb is not None and fb >= 0.2 else "fail")
        checks.append(_chk("probe_binding", "Probe predicted bound at anneal (≥50%)", st,
                           ("%d%% bound @%g °C" % (round(fb * 100), ac)) if fb is not None else "n/a"))
    prim = [pred[k]["frac_bound"] for k in ("forward", "reverse") if k in pred and pred[k]["frac_bound"] is not None]
    if prim:
        st = "pass" if all(x >= 0.5 for x in prim) else "warn"
        checks.append(_chk("primer_binding", "Primers predicted bound at anneal (≥50%)", st,
                           "F/R: " + " / ".join("%d%%" % round(x * 100) for x in prim)))

    # --- provenance completeness ---
    prov_ok = bool(F and Rv)
    checks.append(_chk("provenance", "Provenance complete (sequences, version, conditions)",
                       "pass" if prov_ok else "fail", "checksums + tool version + reaction conditions recorded"))

    statuses = [c["status"] for c in checks if c["status"] != "na"]
    n_na = sum(1 for c in checks if c["status"] == "na")
    miqe_status = "fail" if "fail" in statuses else ("review" if "warn" in statuses else "pass")

    record = dict(
        assay=dict(name=a.get("name", "assay"), gene=a.get("gene", ""), organism=a.get("organism", ""),
                   chemistry=a.get("chemistry", ""), dye=a.get("dye", ""),
                   forward=F, reverse=Rv, probe=P,
                   amplicon_bp=(amp if isinstance(amp, int) else (len(amp) if isinstance(amp, str) else None))),
        provenance=dict(oligoforge_version=_VERSION, generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                        conditions=dict(mv_conc_mM=cond.get("mv_conc"), dv_conc_mM=cond.get("dv_conc"),
                                        dntp_conc_mM=cond.get("dntp_conc"), oligo_nM=cond.get("dna_conc"), anneal_c=ac),
                        checksums=dict(forward=_csum(F), reverse=_csum(Rv), probe=_csum(P))),
        predicted=dict(oligos=pred, amplicon_tm=pred_amp_tm),
        empirical=dict(efficiency_pct=eff, r2=r2, slope=slope, lod=lod, observed_amplicon_tm=observed_amplicon_tm,
                       observed_peaks=observed_peaks),
        checks=checks, n_missing=n_na, miqe_status=miqe_status,
    )
    record["markdown"] = to_markdown(record)
    return record


def to_markdown(rec):
    """Render a validation record as a self-contained MIQE Markdown report."""
    a, pv, pr = rec["assay"], rec["provenance"], rec["predicted"]
    sym = {"pass": "PASS", "fail": "FAIL", "warn": "REVIEW", "na": "n/a"}
    L = []
    L.append("# MIQE validation — %s" % a["name"])
    status_word = {"pass": "MEETS MIQE acceptance criteria", "review": "REVIEW (warnings present)", "fail": "FAILS one or more criteria"}[rec["miqe_status"]]
    L.append("**Overall: %s**%s\n" % (status_word, ("  ·  %d item(s) not yet assessed" % rec["n_missing"]) if rec["n_missing"] else ""))
    L.append("## Assay")
    if a["gene"] or a["organism"]:
        L.append("- Target: %s %s" % (a["gene"], ("(%s)" % a["organism"]) if a["organism"] else ""))
    if a["chemistry"]:
        L.append("- Chemistry: %s%s" % (a["chemistry"], ("  ·  dye %s" % a["dye"]) if a["dye"] else ""))
    L.append("- Forward: `%s`" % a["forward"])
    L.append("- Reverse: `%s`" % a["reverse"])
    if a["probe"]:
        L.append("- Probe: `%s`" % a["probe"])
    if a["amplicon_bp"]:
        L.append("- Amplicon: %s bp" % a["amplicon_bp"])
    L.append("\n## Predicted thermodynamics (engine, at reaction conditions)")
    for k in ("forward", "reverse", "probe"):
        o = pr["oligos"].get(k)
        if o:
            fb = ("%d%% bound" % round(o["frac_bound"] * 100)) if o["frac_bound"] is not None else "n/a"
            L.append("- %s: Tm %.1f °C, ΔG°37 %.1f kcal/mol, %s (%s)" % (k.capitalize(), o["tm"], o["dg37"], fb, o["basis"]))
    if pr["amplicon_tm"] is not None:
        L.append("- Amplicon Tm: %.1f °C" % pr["amplicon_tm"])
    L.append("\n## MIQE checklist")
    for c in rec["checks"]:
        L.append("- [%s] %s%s" % (sym[c["status"]], c["label"], ("  —  %s" % c["detail"]) if c["detail"] else ""))
    L.append("\n## Provenance")
    co = pv["conditions"]
    L.append("- OligoForge %s  ·  %s" % (pv["oligoforge_version"], pv["generated_utc"]))
    L.append("- Conditions: %g mM monovalent, %g mM Mg²⁺, %g mM dNTP, %g nM oligo, anneal %g °C" %
             (co["mv_conc_mM"], co["dv_conc_mM"], co["dntp_conc_mM"], co["oligo_nM"], co["anneal_c"]))
    cs = pv["checksums"]
    L.append("- Checksums: F %s · R %s%s" % (cs["forward"], cs["reverse"], (" · P %s" % cs["probe"]) if cs["probe"] else ""))
    L.append("\n*MIQE acceptance per Bustin et al. 2009. Predicted thermodynamics are from nearest-neighbor "
             "models (SantaLucia 1998; McTigue 2004 for LNA) and are guidance, not a substitute for empirical validation.*")
    return "\n".join(L)
