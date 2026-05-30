"""Multiplex planner: detection-channel conflicts + cross-assay dimer gate.

Given several assays, each with a reporter dye and its oligos, flag (a) assays that
share a dye and so cannot be distinguished in one reaction, and (b) any oligo pair
from *different* assays whose heterodimer dG is at or below a threshold. Self-dimers
and within-assay dimers are out of scope here (use QC / the panel matrix for those).
"""
from . import thermo as T


def check(assays, dimer_threshold=-9.0):
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
                cross.append(dict(assay_a=a1, oligo_a=n1, assay_b=a2, oligo_b=n2, dg=round(dg, 2)))
    cross.sort(key=lambda x: x["dg"])
    return dict(n_assays=len(assays), n_oligos=len(flat), threshold=dimer_threshold,
                channel_conflicts=conflicts, cross_dimers=cross[:50], n_flagged=len(cross))
