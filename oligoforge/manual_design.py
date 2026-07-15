"""Authoritative manual assay analysis and constrained redesign."""
from copy import deepcopy
import difflib
import hashlib
from . import thermo as T, specificity as SP, conservation as C, design as D
from . import ranking as RANK, ranking_explain as REXPLAIN
from . import candidate_retention as CRET

MANUAL_DESIGN_VERSION = "1.4.0"


def _clean_oligo(value, role):
    bare = T.strip_mods(value or "")
    seq, notes, err = T.clean_seq(bare)
    if err:
        raise ValueError("%s: %s" % (role, err))
    if len(seq) > T.MAX_OLIGO_LEN:
        raise ValueError("%s is too long for an oligo" % role)
    return seq, notes


def _clean_template(value):
    seq, notes, err = T.clean_seq(value or "")
    if err:
        raise ValueError("template: %s" % err)
    if len(seq) > T.MAX_TEMPLATE_LEN:
        raise ValueError("template exceeds %d nt" % T.MAX_TEMPLATE_LEN)
    return seq, notes


def map_oligo(seq, template, role, max_mm=2, anchor3=None):
    """Map an oligo to every near-match within ``max_mm`` mismatches.

    Manual diagnosis must show failed 3'-anchor placements rather than hiding them.
    Therefore the default is an exhaustive near-match inventory; callers that are
    constructing viable PCR products still filter on ``extension_eligible``.  The
    legacy ``anchor3=True`` option remains available for explicit amplification-only
    scans.
    """
    primer_role = role in {"forward", "reverse", "primer", "F", "R"}
    scan_anchor = bool(anchor3) if anchor3 is not None else False
    hits = SP.scan_primer_sites(seq, [("template", template)], max_mm=max_mm,
                                anchor3=scan_anchor, label=role)
    rows = []
    for h in hits:
        q3 = bool(h["q3"])
        q3_def = bool(h.get("q3_definite", True))
        extension_eligible = (not primer_role) or (q3 and q3_def)
        if not primer_role:
            status = "not_applicable"
        elif q3 and q3_def:
            status = "exact"
        elif q3:
            status = "possible_ambiguous"
        else:
            status = "mismatch"
        rows.append(dict(
            strand=h["strand"], start=h["lo"], end=h["hi"] + 1,
            mismatches=h["mm"], exact_match=(int(h["mm"]) == 0),
            three_prime_match=q3, three_prime_definite=q3_def,
            three_prime_status=status, extension_eligible=extension_eligible,
            uncertain=h.get("uncertain", False),
            ambiguous_subject_bases=int(h.get("ambiguous_subject_bases", 0)),
            site=h.get("site")))
    rows.sort(key=lambda x: (x["mismatches"], not x["extension_eligible"],
                             x["start"], x["end"], x["strand"]))
    return rows


def _assay_record(forward, reverse, probe, template, profile, products, mappings=None):
    """Build a coordinate-aware record from all-site mapping results.

    Manual assays may legitimately contain target mismatches.  Therefore the
    selected product is mapped through the same near-match scanner used by
    specificity rather than being forced through an exact-string lookup.
    """
    mappings = mappings or {}
    prod = products[0] if len(products) == 1 else None
    fsite = rsite = pxy = None
    if prod:
        lo, hi = prod["span"]
        fsites = mappings.get("forward") or map_oligo(forward, template, "F", max_mm=2)
        rsites = mappings.get("reverse") or map_oligo(reverse, template, "R", max_mm=2)
        fsite = next((x for x in fsites if x["strand"] == "+" and x["start"] == lo), None)
        rsite = next((x for x in rsites if x["strand"] == "-" and x["end"] - 1 == hi), None)
    amp = int(prod["size"]) if prod else 0
    pi = None
    if probe:
        ps = mappings.get("probe") or map_oligo(probe, template, "P", max_mm=2, anchor3=False)
        if prod:
            inside = [x for x in ps if prod["span"][0] <= x["start"] and x["end"] - 1 <= prod["span"][1]]
            # Do not silently choose among multiple placements.
            pxy = inside[0] if len(inside) == 1 else None
        pt = T.tm(probe)
        pi = dict(probe=probe, strand=(pxy or {}).get("strand"), tm=pt,
                  offset=pt - max(T.tm(forward), T.tm(reverse)),
                  hairpin_dg=T.hairpin(probe)[0], self_dimer=T.self_dimer(probe),
                  dimer_f=T.hetero_dimer(probe, forward), dimer_r=T.hetero_dimer(probe, reverse))
    return dict(forward=forward, reverse=reverse, probe=probe or None, amplicon=amp,
                amplicon_tm=(round(T.amplicon_tm(template[prod["span"][0]:prod["span"][1] + 1]), 1) if prod else None),
                pair_tm_gap=abs(T.tm(forward) - T.tm(reverse)),
                pair_dimer=min(T.self_dimer(forward), T.self_dimer(reverse), T.hetero_dimer(forward, reverse)),
                f_tm=T.tm(forward), r_tm=T.tm(reverse), probe_info=pi,
                f_xy=([fsite["start"], fsite["end"]] if fsite else None),
                r_xy=([rsite["start"], rsite["end"]] if rsite else None),
                probe_xy=([pxy["start"], pxy["end"]] if pxy else None),
                amplicon_xy=([prod["span"][0], prod["span"][1] + 1] if prod else None),
                gblock=None, quality_flags=[])


def analyze_assay(forward, reverse, template, profile, probe=None, targets=None, offs=None,
                  objective="balanced", max_mm=2):
    """Analyze a manually entered assay using the same ranking evidence as auto-design."""
    tmpl, tnotes = _clean_template(template)
    f, fnotes = _clean_oligo(forward, "forward")
    r, rnotes = _clean_oligo(reverse, "reverse")
    p = None; pnotes = []
    if probe:
        p, pnotes = _clean_oligo(probe, "probe")
    maps = dict(forward=map_oligo(f, tmpl, "F", max_mm=max_mm),
                reverse=map_oligo(r, tmpl, "R", max_mm=max_mm),
                probe=(map_oligo(p, tmpl, "P", max_mm=max_mm, anchor3=False) if p else []))
    fasta = ">template\n%s\n" % tmpl
    epcr = SP.in_silico_pcr_offline(f, r, fasta, probe=p, max_mm=max_mm,
                                    min_product=max(20, int(profile["amp_min"])),
                                    max_product=max(int(profile["amp_max"]), int(profile["amp_min"]) + 1))
    products = epcr.get("products") or []
    assay = _assay_record(f, r, p, tmpl, profile, products, mappings=maps)
    target_set = targets or [tmpl]
    cons = {"F": C.conservation(f, target_set, 0.6), "R": C.conservation(r, target_set, 0.6)}
    if p:
        cons["P"] = C.conservation(p, target_set, 0.6)
    scored = dict(score=0.0, score_raw=0.0, quality_penalty=0.0, amplicon_penalty=0.0,
                  assay=assay, conservation=cons, discrimination=None)
    ranked, obj = RANK.rank_candidates([scored], target_set, offs or [], profile, objective_name=objective)
    row = ranked[0]
    row["rank_explanation"] = REXPLAIN.explain(row)
    manual_failures = []
    if len(products) == 0:
        manual_failures.append("no coherent forward/reverse product maps to the template")
    elif len(products) > 1:
        manual_failures.append("multiple coherent products map to the template")
    if p:
        within = [x for x in maps["probe"] if any(pr["span"][0] <= x["start"] and x["end"] <= pr["span"][1] + 1 for pr in products)]
        if not within:
            manual_failures.append("probe does not map inside a predicted product")
        elif len(within) > 1:
            manual_failures.append("probe has multiple placements inside predicted products")
    if manual_failures:
        row["evidence"]["hard_valid"] = False
        row["evidence"]["hard_failures"] = sorted(set(row["evidence"]["hard_failures"] + manual_failures))
        row["rank_explanation"] = REXPLAIN.explain(row)
    input_hashes = {
        "template_sha256": hashlib.sha256(tmpl.encode()).hexdigest(),
        "forward_sha256": hashlib.sha256(f.encode()).hexdigest(),
        "reverse_sha256": hashlib.sha256(r.encode()).hexdigest(),
        "probe_sha256": (hashlib.sha256(p.encode()).hexdigest() if p else None),
        "target_corpus_sha256": hashlib.sha256("\n".join(target_set).encode()).hexdigest(),
        "offtarget_corpus_sha256": hashlib.sha256("\n".join(offs or []).encode()).hexdigest(),
    }
    run_manifest = RANK.manifest(
        obj, {"fully_annotated_candidates": 1}, input_hashes=input_hashes,
        constraints={"workflow": "manual_analysis", "max_mismatches": int(max_mm),
                     "probe_required": bool(p)})
    return dict(version=MANUAL_DESIGN_VERSION, template=tmpl,
                notes=tnotes + fnotes + rnotes + pnotes, mappings=maps,
                predicted_products=products, epcr=epcr, candidate=row,
                objective_profile=obj, ranker_manifest=run_manifest,
                ambiguous_mapping=(len(products) != 1),
                authoritative_calculations="shared thermodynamics, specificity, conservation, and structured ranking")


def _sequence_edit(old, new, role):
    old = str(old or ""); new = str(new or "")
    matcher = difflib.SequenceMatcher(a=old, b=new, autojunk=False)
    operations = []
    for tag, a0, a1, b0, b1 in matcher.get_opcodes():
        if tag == "equal":
            continue
        operations.append({"operation": tag, "old_span": [a0, a1], "new_span": [b0, b1],
                           "old": old[a0:a1], "new": new[b0:b1]})
    return {"component": role, "changed": old != new, "old": old or None, "new": new or None,
            "old_length": len(old), "new_length": len(new), "operations": operations}


def _comparison_snapshot(analysis):
    row = (analysis or {}).get("candidate") or {}
    assay = row.get("assay") or {}
    evidence = row.get("evidence") or {}
    off = evidence.get("offtarget") or {}
    robust = evidence.get("condition_robustness") or {}
    return {
        "hard_valid": bool(evidence.get("hard_valid")),
        "hard_failures": list(evidence.get("hard_failures") or []),
        "coherent_products": len((analysis or {}).get("predicted_products") or []),
        "ambiguous_mapping": bool((analysis or {}).get("ambiguous_mapping")),
        "target_coverage": evidence.get("target_coverage"),
        "worst_isolate_3prime": evidence.get("worst_isolate_3prime"),
        "probe_mean_identity": evidence.get("probe_mean_identity"),
        "signal_offtargets": int(off.get("signal_subjects", 0)),
        "product_offtargets": int(off.get("product_subjects", 0)),
        "condition_robustness": robust.get("valid_fraction"),
        "triplet_penalty": evidence.get("triplet_penalty"),
        "practical_penalty": evidence.get("practical_penalty"),
        "forward_tm": assay.get("f_tm"),
        "reverse_tm": assay.get("r_tm"),
        "primer_tm_gap": assay.get("pair_tm_gap"),
        "probe_tm": ((assay.get("probe_info") or {}).get("tm") if assay.get("probe") else None),
        "amplicon_bp": assay.get("amplicon"),
        "worst_dimer_dg": evidence.get("worst_dimer"),
        "rank_key": list(row.get("rank_key") or []),
    }


def _metric_change(name, before, after, higher_is_better=None):
    delta = None
    if isinstance(before, (int, float)) and isinstance(after, (int, float)):
        delta = round(float(after) - float(before), 6)
    state = "changed" if before != after else "unchanged"
    if higher_is_better is not None and delta not in (None, 0):
        state = "improved" if ((delta > 0) == bool(higher_is_better)) else "worsened"
    return {"metric": name, "before": before, "after": after, "delta": delta,
            "interpretation": state}


def compare_edits(baseline_forward, baseline_reverse, edited_forward, edited_reverse,
                  template, profile, baseline_probe=None, edited_probe=None,
                  targets=None, offs=None, objective="balanced", max_mm=2):
    """Compare an exact manual edit using complete authoritative evidence.

    The function intentionally recalculates both assays under the same target,
    off-target, chemistry and reaction-condition context.  It does not infer that
    a computational improvement guarantees better efficiency or fluorescence.
    """
    baseline = analyze_assay(baseline_forward, baseline_reverse, template, profile,
                             baseline_probe, targets=targets, offs=offs,
                             objective=objective, max_mm=max_mm)
    edited = analyze_assay(edited_forward, edited_reverse, template, profile,
                           edited_probe, targets=targets, offs=offs,
                           objective=objective, max_mm=max_mm)
    brow, erow = baseline["candidate"], edited["candidate"]
    bkey, ekey = tuple(brow.get("rank_key") or ()), tuple(erow.get("rank_key") or ())
    if ekey < bkey:
        preference = "edited_assay_preferred"
        pref_detail = REXPLAIN.preference_strength(erow, brow)
    elif bkey < ekey:
        preference = "baseline_assay_preferred"
        pref_detail = REXPLAIN.preference_strength(brow, erow)
    else:
        preference = "computationally_indistinguishable"
        pref_detail = REXPLAIN.preference_strength(erow, brow)

    bs, es = _comparison_snapshot(baseline), _comparison_snapshot(edited)
    direction = {
        "target_coverage": True, "worst_isolate_3prime": True,
        "probe_mean_identity": True, "condition_robustness": True,
        "signal_offtargets": False, "product_offtargets": False,
        "triplet_penalty": False, "practical_penalty": False,
        "primer_tm_gap": False, "coherent_products": None,
        "worst_dimer_dg": True,
    }
    metrics = [_metric_change(k, bs.get(k), es.get(k), direction.get(k)) for k in (
        "hard_valid", "coherent_products", "ambiguous_mapping", "target_coverage",
        "worst_isolate_3prime", "probe_mean_identity", "signal_offtargets",
        "product_offtargets", "condition_robustness", "triplet_penalty",
        "practical_penalty", "forward_tm", "reverse_tm", "primer_tm_gap",
        "probe_tm", "amplicon_bp", "worst_dimer_dg")]
    resolved = sorted(set(bs["hard_failures"]) - set(es["hard_failures"]))
    introduced = sorted(set(es["hard_failures"]) - set(bs["hard_failures"]))
    improvements = [m for m in metrics if m["interpretation"] == "improved"]
    worsenings = [m for m in metrics if m["interpretation"] == "worsened"]
    if resolved:
        improvements.insert(0, {"metric": "hard_failures_resolved", "before": bs["hard_failures"],
                                "after": es["hard_failures"], "delta": None,
                                "interpretation": "improved", "items": resolved})
    if introduced:
        worsenings.insert(0, {"metric": "hard_failures_introduced", "before": bs["hard_failures"],
                              "after": es["hard_failures"], "delta": None,
                              "interpretation": "worsened", "items": introduced})
    edits = [_sequence_edit(baseline_forward, edited_forward, "forward"),
             _sequence_edit(baseline_reverse, edited_reverse, "reverse"),
             _sequence_edit(baseline_probe, edited_probe, "probe")]
    return {
        "version": MANUAL_DESIGN_VERSION,
        "preference": preference,
        "preference_detail": pref_detail,
        "sequence_edits": [x for x in edits if x["changed"]],
        "hard_failures_resolved": resolved,
        "hard_failures_introduced": introduced,
        "improvements": improvements,
        "worsenings": worsenings,
        "metric_changes": metrics,
        "baseline": baseline,
        "edited": edited,
        "edited_maps_to_intended_target": (es["coherent_products"] == 1 and not es["ambiguous_mapping"]),
        "interpretation": ("The ordering is conditional on the supplied template/panels, reaction conditions, "
                           "chemistry profile and modeled evidence. Re-test efficiency, specificity, product identity, "
                           "probe signal and multiplex behavior at the bench."),
    }


def _changes(old, new):
    rows = []
    for role, key in (("forward", "forward"), ("reverse", "reverse"), ("probe", "probe")):
        a, b = old.get(key), new.get(key)
        if a != b:
            rows.append(dict(component=role, old=a, new=b))
    return rows


def _mapped_primer_tuples(seq, template, role, max_mm=2):
    """Return design.py-compatible tuples for every viable locked-primer placement."""
    strand = "+" if role == "forward" else "-"
    return [(h["start"], h["end"], seq) for h in map_oligo(seq, template, role, max_mm=max_mm)
            if h["strand"] == strand and h.get("extension_eligible")]


def _pair_pool(template, profile, forward, reverse, locks, max_mm=2, limit=60):
    lock_f = bool(locks.get("forward") or locks.get("primer_pair"))
    lock_r = bool(locks.get("reverse") or locks.get("primer_pair"))
    if not lock_f and not lock_r:
        fwd, rev = D.enumerate_primers(template, profile)
        pairs = D.pair_primers(fwd, rev, profile)
        kept, _ledger = CRET.retain_pairs_diverse(
            pairs, limit=max(1, int(limit)),
            region_size=max(60, int(profile.get("amp_max", 150))),
            amplicon_bin=25, per_near=2)
        return kept

    # Enumerate only the unlocked side.  A locked component must enter the
    # construction pool even when it would not survive preliminary top-N pruning.
    if lock_f and lock_r:
        fwd = _mapped_primer_tuples(forward, template, "forward", max_mm=max_mm)
        rev = _mapped_primer_tuples(reverse, template, "reverse", max_mm=max_mm)
    else:
        all_f, all_r = D.enumerate_primers(template, profile)
        fwd = _mapped_primer_tuples(forward, template, "forward", max_mm=max_mm) if lock_f else all_f
        rev = _mapped_primer_tuples(reverse, template, "reverse", max_mm=max_mm) if lock_r else all_r
    ac = float(profile.get("anneal_c", T.ANNEAL_C))
    rows = []
    both_locked = lock_f and lock_r
    # Pre-limit an unlocked side by Tm fit and coordinate plausibility before
    # costly dimer calculations.  This does not affect the locked component.
    if lock_f and not lock_r:
        rev = sorted(rev, key=lambda x: (abs(T.tm(x[2]) - T.tm(forward)), x[0], x[2]))[:600]
    elif lock_r and not lock_f:
        fwd = sorted(fwd, key=lambda x: (abs(T.tm(x[2]) - T.tm(reverse)), x[0], x[2]))[:600]
    for fs, fe, f in fwd:
        ft = T.tm(f)
        for rs, re, r in rev:
            amp = re - fs
            if not (int(profile["amp_min"]) <= amp <= int(profile["amp_max"])):
                continue
            if rs - fe < int(profile["min_probe_gap"]):
                continue
            rt = T.tm(r)
            gap = abs(ft - rt)
            if not both_locked and gap > float(profile["pair_tm_gap_max"]):
                continue
            dga = T.hetero_dimer_full(f, r, ac)[1]
            if not both_locked and dga <= float(profile["pair_dimer_min"]):
                continue
            worst = min(T.self_dimer(f), T.self_dimer(r), T.hetero_dimer(f, r))
            score = abs((ft + rt) / 2 - float(profile["tm_opt"])) + gap
            score += 2.0 * max(0.0, -5.5 - worst)
            rows.append(dict(score=score, fstart=fs, fend=fe, f=f,
                             rstart=rs, rend=re, r=r, amp=amp, gap=gap,
                             dimer=round(worst, 2)))
    rows.sort(key=lambda x: (x["score"], x["fstart"], x["rstart"], x["f"], x["r"]))
    kept, _ledger = CRET.retain_pairs_diverse(
        rows, limit=max(1, int(limit)),
        region_size=max(60, int(profile.get("amp_max", 150))),
        amplicon_bin=25, per_near=2)
    return kept


def _locked_probe_infos(probe, template, pair, profile, max_mm=2):
    ac = float(profile.get("anneal_c", T.ANNEAL_C))
    out = []
    for hit in map_oligo(probe, template, "probe", max_mm=max_mm, anchor3=False):
        if not (pair["fend"] <= hit["start"] and hit["end"] <= pair["rstart"]):
            continue
        pt = T.tm(probe)
        hd37, hda, htm = T.hairpin_full(probe, ac)
        sd37, sda, sdtm = T.self_dimer_full(probe, ac)
        df37, dfa, dftm = T.hetero_dimer_full(probe, pair["f"], ac)
        dr37, dra, drtm = T.hetero_dimer_full(probe, pair["r"], ac)
        out.append(dict(probe=probe, strand=hit["strand"], tm=pt,
                        offset=pt - max(T.tm(pair["f"]), T.tm(pair["r"])),
                        hairpin_dg=hd37, hairpin_dg_anneal=hda, hairpin_tm=htm,
                        self_dimer=sd37, self_dimer_anneal=sda, self_dimer_tm=sdtm,
                        dimer_f=df37, dimer_f_anneal=dfa, dimer_f_tm=dftm,
                        dimer_r=dr37, dimer_r_anneal=dra, dimer_r_tm=drtm,
                        start=hit["start"], end=hit["end"],
                        mismatches=hit["mismatches"], uncertain=hit["uncertain"],
                        preliminary_penalty=abs((pt - max(T.tm(pair["f"]), T.tm(pair["r"]))) - 9.0)))
    out.sort(key=lambda x: (x["mismatches"], x["preliminary_penalty"], x["start"], x["strand"]))
    return out


def _overlaps(span, regions):
    if not span:
        return False
    return any(not (int(span[1]) <= int(lo) or int(span[0]) >= int(hi)) for lo, hi in regions)


def constrained_redesign(forward, reverse, template, profile, probe=None, locks=None,
                         objective="balanced", max_results=8, max_shift=None,
                         excluded_regions=None, max_mm=2, amp_min=None, amp_max=None,
                         required_region=None, targets=None, offs=None,
                         pair_limit=60, probes_per_pair=3, full_annotation_limit=60,
                         base_analysis=None):
    """Redesign unlocked components with true component locks.

    Locked primers/probes are injected into candidate construction rather than
    filtered after a generic beam.  This prevents a valid locked component from
    disappearing merely because it was not in the preliminary top-N pool.
    """
    tmpl, _ = _clean_template(template)
    target_set = list(targets or [tmpl])
    off_set = list(offs or [])
    f, _ = _clean_oligo(forward, "forward")
    r, _ = _clean_oligo(reverse, "reverse")
    p = _clean_oligo(probe, "probe")[0] if probe else None
    locks = dict(locks or {})
    if locks.get("primer_pair"):
        locks["forward"] = locks["reverse"] = True
    if locks.get("probe") and not p:
        raise ValueError("probe lock requested but no probe was supplied")
    old = dict(forward=f, reverse=r, probe=p)
    excluded_regions = [(int(lo), int(hi)) for lo, hi in (excluded_regions or []) if int(lo) < int(hi)]
    local_profile = dict(profile)
    if amp_min is not None:
        local_profile["amp_min"] = max(20, int(amp_min))
    if amp_max is not None:
        local_profile["amp_max"] = max(int(local_profile["amp_min"]), int(amp_max))

    pair_limit = max(1, min(int(pair_limit), 200))
    probes_per_pair = max(1, min(int(probes_per_pair), 12))
    full_annotation_limit = max(1, min(int(full_annotation_limit), 240))
    pairs = _pair_pool(tmpl, local_profile, f, r, locks, max_mm=max_mm, limit=pair_limit)
    assays = []
    pairs_without_probe = 0
    probe_placements = 0
    for pair in pairs:
        if local_profile.get("no_probe"):
            probes = [None]
        elif locks.get("probe"):
            probes = _locked_probe_infos(p, tmpl, pair, local_profile, max_mm=max_mm)
            probe_placements += len(probes)
        else:
            probes = D.enumerate_probe_candidates(tmpl, pair["fend"], pair["rstart"],
                                                   pair["f"], pair["r"], local_profile, limit=probes_per_pair)
        if not probes:
            pairs_without_probe += 1
            continue
        for pi in probes:
            a = D.assay_from_pair(tmpl, pair, pi)
            a["amplicon_xy"] = [a["f_xy"][0], a["r_xy"][1]]
            if excluded_regions and any(_overlaps(a.get(k), excluded_regions)
                                        for k in ("f_xy", "r_xy", "probe_xy")):
                continue
            if required_region:
                lo, hi = int(required_region[0]), int(required_region[1])
                if not (a["amplicon_xy"][0] <= lo and hi <= a["amplicon_xy"][1]):
                    continue
            assays.append(a)

    if base_analysis is None:
        base_analysis = analyze_assay(f, r, tmpl, local_profile, p, targets=target_set,
                                      offs=off_set, objective=objective, max_mm=max_mm)
    base = base_analysis["candidate"]
    f0 = (base.get("assay") or {}).get("f_xy")
    r0 = (base.get("assay") or {}).get("r_xy")
    p0 = (base.get("assay") or {}).get("probe_xy")
    eligible = []
    shift_rejected = 0
    for a in assays:
        if max_shift is not None:
            lim = int(max_shift)
            if not locks.get("forward") and f0 and abs(a["f_xy"][0] - f0[0]) > lim:
                shift_rejected += 1; continue
            if not locks.get("reverse") and r0 and abs(a["r_xy"][0] - r0[0]) > lim:
                shift_rejected += 1; continue
            if not locks.get("probe") and p0 and a.get("probe_xy") and abs(a["probe_xy"][0] - p0[0]) > lim:
                shift_rejected += 1; continue
        a["candidate_rank"] = round(D._candidate_rank(a, local_profile), 6)
        a["search_window_start"] = int((a.get("f_xy") or [0])[0])
        eligible.append(a)

    # Full specificity/coverage/robustness annotation is the expensive stage.  Preserve
    # target-region and trade-off diversity before imposing a documented interactive
    # budget; do not simply keep a dense block of nearly identical triplets.
    retained, retention_ledger = CRET.retain_diverse(
        eligible, limit=full_annotation_limit,
        region_size=max(40, int(local_profile.get("amp_max", 150))),
        per_region=max(2, min(12, full_annotation_limit // 4 or 2)), per_near=2)
    kept = []
    for a in retained:
        cons = {"F": C.conservation(a["forward"], target_set, 0.6),
                "R": C.conservation(a["reverse"], target_set, 0.6)}
        if a.get("probe"):
            cons["P"] = C.conservation(a["probe"], target_set, 0.6)
        kept.append(dict(score=0.0, score_raw=0.0, quality_penalty=0.0,
                         amplicon_penalty=0.0, assay=a, conservation=cons,
                         discrimination=None))
    ranked, obj = RANK.rank_candidates(kept, target_set, off_set, local_profile,
                                        objective_name=objective) if kept else ([], RANK.get_profile(objective, no_probe=local_profile.get("no_probe", False)))
    finalists = RANK.select_finalists(ranked, n=max_results)
    for row in finalists:
        row["changes"] = _changes(old, row["assay"])
        row["components_retained"] = [k for k in ("forward", "reverse", "probe")
                                      if old.get(k) == row["assay"].get(k)]
        row["components_changed"] = [x["component"] for x in row["changes"]]
        row["rank_explanation"] = REXPLAIN.explain(row, base)
    ledger = dict(stage="locked_component_joint_search", search_status="heuristic_bounded",
                  pairs_considered=len(pairs), triplets_constructed=len(assays),
                  candidates_after_constraints=len(eligible), shift_rejected=shift_rejected,
                  candidates_fully_annotated=len(ranked), pairs_without_probe=pairs_without_probe,
                  locked_probe_placements=probe_placements, locks=locks,
                  diversity_retention=retention_ledger,
                  limits=dict(pair_limit=pair_limit,
                              probes_per_pair=(1 if locks.get("probe") else probes_per_pair),
                              full_annotation_limit=full_annotation_limit))
    redesign_hashes = {
        "template_sha256": hashlib.sha256(tmpl.encode()).hexdigest(),
        "forward_sha256": hashlib.sha256(f.encode()).hexdigest(),
        "reverse_sha256": hashlib.sha256(r.encode()).hexdigest(),
        "probe_sha256": (hashlib.sha256(p.encode()).hexdigest() if p else None),
        "target_corpus_sha256": hashlib.sha256("\n".join(target_set).encode()).hexdigest(),
        "offtarget_corpus_sha256": hashlib.sha256("\n".join(off_set).encode()).hexdigest(),
    }
    run_manifest = RANK.manifest(
        obj, ledger.get("limits") or {}, input_hashes=redesign_hashes,
        constraints={"workflow": "constrained_redesign", "locks": locks,
                     "max_shift": max_shift, "excluded_regions": excluded_regions,
                     "required_region": required_region,
                     "amp_min": local_profile.get("amp_min"), "amp_max": local_profile.get("amp_max")})
    return dict(version=MANUAL_DESIGN_VERSION, locks=locks, base=base,
                candidates=finalists, n_screened=len(ranked), search_ledger=ledger,
                objective_profile=obj, ranker_manifest=run_manifest,
                note=(None if finalists else "No hard-valid candidate satisfied every lock and constraint. Locked components were not replaced."))
