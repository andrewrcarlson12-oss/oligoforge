"""Multiplex planner: detection-channel conflicts + cross-assay dimer gate.

Given several assays, each with a reporter dye and its oligos, flag (a) assays that
share a dye and so cannot be distinguished in one reaction, and (b) any oligo pair
from *different* assays whose heterodimer dG is at or below a threshold. Self-dimers
and within-assay dimers are out of scope here (use QC / the panel matrix for those).
"""
from . import thermo as T


def check(assays, dimer_threshold=-6.0, amp_tm_gap=2.0):
    by_dye = {}
    for a in assays:
        dye = (a.get("dye") or "").strip()
        if dye:
            by_dye.setdefault(dye, []).append(a.get("name") or "(unnamed)")
    conflicts = [dict(dye=d, assays=names) for d, names in by_dye.items() if len(names) > 1]

    flat = []
    for a in assays:
        nm = a.get("name") or "(unnamed)"
        for o in a.get("oligos", []):
            if o.get("seq"):
                flat.append((nm, o.get("name") or "?", o["seq"].upper()))

    cross = []
    for i in range(len(flat)):
        for j in range(i + 1, len(flat)):
            a1, n1, s1 = flat[i]
            a2, n2, s2 = flat[j]
            if a1 == a2:
                continue
            dg = T.hetero_dimer(s1, s2)
            if dg <= dimer_threshold:
                # Annotate which flagged dimers are 3'-ENGAGED (extension-competent, the dangerous
                # kind) vs internal-only. end_dg is the stronger of the two 3'-end stabilities;
                # three_prime flags a 3'-anchored dimer at the conventional ~-5 kcal/mol screen.
                end_dg = min(T.end_stability(s1, s2), T.end_stability(s2, s1))
                cross.append(dict(assay_a=a1, oligo_a=n1, assay_b=a2, oligo_b=n2, dg=round(dg, 2),
                                  end_dg=round(end_dg, 2), three_prime=end_dg <= -5.0))
    cross.sort(key=lambda x: x["dg"])

    # SYBR melt overlap: SYBR assays are told apart by melt-peak Tm, so two SYBR
    # amplicons with near-equal predicted Tm can't be resolved in one reaction.
    def _is_sybr(a):
        if a.get("sybr") is not None:
            return bool(a["sybr"])
        return not any((o.get("name") == "P") for o in a.get("oligos", []))
    sybr = [a for a in assays if _is_sybr(a) and a.get("amplicon_tm") is not None]
    melt = []
    for i in range(len(sybr)):
        for j in range(i + 1, len(sybr)):
            d = abs(float(sybr[i]["amplicon_tm"]) - float(sybr[j]["amplicon_tm"]))
            if d < amp_tm_gap:
                melt.append(dict(assay_a=sybr[i].get("name") or "(unnamed)",
                                 assay_b=sybr[j].get("name") or "(unnamed)",
                                 tm_a=round(float(sybr[i]["amplicon_tm"]), 1),
                                 tm_b=round(float(sybr[j]["amplicon_tm"]), 1),
                                 delta=round(d, 1)))
    melt.sort(key=lambda x: x["delta"])
    return dict(n_assays=len(assays), n_oligos=len(flat), threshold=dimer_threshold,
                channel_conflicts=conflicts, cross_dimers=cross[:50], n_flagged=len(cross),
                n_sybr=len(sybr), amp_tm_gap=amp_tm_gap, melt_overlaps=melt)
