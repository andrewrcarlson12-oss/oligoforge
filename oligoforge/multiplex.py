"""Multiplex compatibility screening.

The planner evaluates detection-channel collisions, cross-assay oligo dimers and
SYBR melt-peak overlap.  It is deliberately conservative and rejects malformed
oligos rather than letting the thermodynamic layer silently resolve or discard
characters.  Modified oligos are analysed as their DNA backbone and the result
reports that limitation explicitly.
"""
from . import thermo as T


def _clean_oligo(raw, assay_name, oligo_name):
    """Return a validated DNA/IUPAC backbone plus audit metadata.

    Thermodynamic functions accept only a DNA backbone.  IDT modification blocks,
    LNA ``+`` markers and phosphorothioate markers are therefore stripped only
    after the original string has been retained for provenance.  Any other invalid
    symbol is a hard error; it is never silently removed or reinterpreted.
    """
    original = str(raw or "").strip()
    if not original:
        return None, "%s / %s: empty sequence" % (assay_name, oligo_name)
    modified = any(mark in original for mark in ("/", "+", "*"))
    bare = T.strip_mods(original)
    clean, notes, err = T.clean_seq(bare)
    if err:
        return None, "%s / %s: %s" % (assay_name, oligo_name, err)
    if len(clean) < 6:
        return None, "%s / %s: too short for multiplex thermodynamics (need >=6 nt)" % (
            assay_name, oligo_name)
    if len(clean) > 60:
        return None, "%s / %s: %d nt exceeds the 60-nt nearest-neighbour dimer limit" % (
            assay_name, oligo_name, len(clean))
    return dict(seq=clean, original=original, modified=modified, notes=notes), None


def check(assays, dimer_threshold=-6.0, amp_tm_gap=2.0):
    """Screen a collection of assays for multiplex conflicts.

    Assay identity is based on list position (or an explicit ``id``), not display
    name.  This prevents two distinct assays that happen to share a name from being
    mistaken for one assay and skipped during cross-dimer analysis.
    """
    assays = list(assays or [])
    if len(assays) > 100:
        return dict(error="too many assays for one multiplex screen (limit 100)",
                    n_assays=len(assays))
    try:
        dimer_threshold = float(dimer_threshold)
        amp_tm_gap = float(amp_tm_gap)
    except (TypeError, ValueError):
        return dict(error="dimer_threshold and amp_tm_gap must be numeric")
    if not (-100.0 <= dimer_threshold <= 20.0):
        return dict(error="dimer_threshold out of range (-100..20 kcal/mol)")
    if not (0.0 <= amp_tm_gap <= 30.0):
        return dict(error="amp_tm_gap out of range (0..30 C)")

    normalized = []
    invalid = []
    modified_count = 0
    for ai, assay in enumerate(assays):
        if not isinstance(assay, dict):
            invalid.append(dict(assay_index=ai, error="assay must be an object"))
            continue
        name = str(assay.get("name") or "(unnamed)")
        assay_id = str(assay.get("id") or "assay-%d" % (ai + 1))
        dye_raw = str(assay.get("dye") or "").strip()
        dye_key = dye_raw.upper()
        clean_oligos = []
        for oi, oligo in enumerate(assay.get("oligos") or []):
            if not isinstance(oligo, dict):
                invalid.append(dict(assay_index=ai, assay=name, oligo_index=oi,
                                    error="oligo must be an object"))
                continue
            oname = str(oligo.get("name") or "?")
            cleaned, err = _clean_oligo(oligo.get("seq"), name, oname)
            if err:
                invalid.append(dict(assay_index=ai, assay=name, oligo_index=oi,
                                    oligo=oname, error=err))
                continue
            modified_count += int(cleaned["modified"])
            clean_oligos.append(dict(name=oname, **cleaned))
        normalized.append(dict(index=ai, id=assay_id, name=name, dye=dye_raw,
                               dye_key=dye_key, sybr=assay.get("sybr"),
                               amplicon_tm=assay.get("amplicon_tm"), oligos=clean_oligos))

    if sum(len(a["oligos"]) for a in normalized) > 400:
        return dict(error="too many oligos for one multiplex screen (limit 400)",
                    n_assays=len(assays))
    if invalid:
        return dict(error="one or more multiplex oligos are invalid", invalid_oligos=invalid,
                    n_assays=len(assays), n_invalid=len(invalid))

    by_dye = {}
    for assay in normalized:
        if assay["dye_key"]:
            by_dye.setdefault(assay["dye_key"], []).append(assay)
    conflicts = [
        dict(dye=dye, assays=[a["name"] for a in members], assay_ids=[a["id"] for a in members])
        for dye, members in sorted(by_dye.items()) if len(members) > 1
    ]

    flat = []
    for assay in normalized:
        for oligo in assay["oligos"]:
            flat.append(dict(assay_index=assay["index"], assay_id=assay["id"],
                             assay_name=assay["name"], oligo_name=oligo["name"],
                             seq=oligo["seq"], modified=oligo["modified"]))

    cross = []
    for i in range(len(flat)):
        for j in range(i + 1, len(flat)):
            one, two = flat[i], flat[j]
            if one["assay_index"] == two["assay_index"]:
                continue
            dg37, dg_anneal, dimer_tm = T.hetero_dimer_full(one["seq"], two["seq"])
            if dg37 <= dimer_threshold:
                end_dg = min(T.end_stability(one["seq"], two["seq"]),
                             T.end_stability(two["seq"], one["seq"]))
                cross.append(dict(
                    assay_a=one["assay_name"], assay_a_id=one["assay_id"],
                    oligo_a=one["oligo_name"], assay_b=two["assay_name"],
                    assay_b_id=two["assay_id"], oligo_b=two["oligo_name"],
                    dg=round(dg37, 2), dg_anneal=round(dg_anneal, 2),
                    dimer_tm=round(dimer_tm, 1), anneal_c=T.ANNEAL_C,
                    end_dg=round(end_dg, 2), three_prime=end_dg <= -5.0,
                    modified_backbone_only=bool(one["modified"] or two["modified"])))
    cross.sort(key=lambda row: (row["dg"], row["end_dg"]))

    def _is_sybr(assay):
        if assay["sybr"] is not None:
            return bool(assay["sybr"])
        return not any(o["name"].upper() == "P" for o in assay["oligos"])

    sybr = []
    invalid_tm = []
    for assay in normalized:
        if not _is_sybr(assay) or assay["amplicon_tm"] is None:
            continue
        try:
            tm = float(assay["amplicon_tm"])
        except (TypeError, ValueError):
            invalid_tm.append(dict(assay=assay["name"], assay_id=assay["id"],
                                   value=assay["amplicon_tm"]))
            continue
        if tm != tm or tm in (float("inf"), float("-inf")):
            invalid_tm.append(dict(assay=assay["name"], assay_id=assay["id"], value=str(tm)))
            continue
        sybr.append((assay, tm))

    melt = []
    for i in range(len(sybr)):
        for j in range(i + 1, len(sybr)):
            a, ta = sybr[i]; b, tb = sybr[j]
            delta = abs(ta - tb)
            if delta < amp_tm_gap:
                melt.append(dict(assay_a=a["name"], assay_a_id=a["id"],
                                 assay_b=b["name"], assay_b_id=b["id"],
                                 tm_a=round(ta, 1), tm_b=round(tb, 1), delta=round(delta, 1)))
    melt.sort(key=lambda row: row["delta"])

    note = None
    if modified_count:
        note = ("Modified oligos were screened as their unmodified DNA backbones. Dye, quencher, "
                "LNA and linkage effects require vendor-specific confirmation.")
    return dict(n_assays=len(normalized), n_oligos=len(flat), threshold=dimer_threshold,
                anneal_c=T.ANNEAL_C, channel_conflicts=conflicts,
                cross_dimers=cross[:50], n_flagged=len(cross),
                n_sybr=len(sybr), amp_tm_gap=amp_tm_gap, melt_overlaps=melt,
                invalid_amplicon_tms=invalid_tm, modified_oligos=modified_count,
                note=note)
