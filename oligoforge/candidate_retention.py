"""Diversity-preserving candidate retention with auditable attrition.

Preliminary ranking is intentionally cheap.  It may decide which candidates receive
expensive annotation, but it must not allow one dense target region, amplicon-length
band, or near-duplicate family to consume the whole budget.  Every bounded retention
operation returns an aggregate ledger and a deterministic per-candidate decision trace.
"""
from collections import defaultdict, deque
import hashlib
import json

RETENTION_VERSION = "1.2.0"


def identity_key(assay):
    return (assay.get("forward"), assay.get("reverse"), assay.get("probe"))


def pair_identity_key(pair):
    return (pair.get("f"), pair.get("r"), int(pair.get("fstart", 0)), int(pair.get("rend", 0)))


def probe_identity_key(probe):
    return (probe.get("probe"), probe.get("strand"), int(probe.get("start", 0)), int(probe.get("end", 0)))


def near_key(assay, pos_bin=12):
    f = assay.get("f_xy") or [0, 0]
    r = assay.get("r_xy") or [0, 0]
    p = assay.get("probe_xy") or [0, 0]
    return (int(f[0]) // pos_bin, int(r[1]) // pos_bin, int(p[0]) // pos_bin if p else -1,
            int(assay.get("amplicon") or 0) // 10)


def _prelim(assay):
    return float(assay.get("candidate_rank", assay.get("preliminary_penalty", 1e12)))


def _pair_prelim(pair):
    # Mirror design.py's preliminary pair ordering without recomputing thermodynamics.
    return float(pair.get("score", 1e12)) + 2.0 * max(0.0, -5.5 - float(pair.get("dimer", 0.0)))


def _probe_prelim(probe):
    return float(probe.get("preliminary_penalty", 1e12))


def _stable_id(kind, identity, occurrence=None):
    raw = json.dumps([kind, list(identity), occurrence], separators=(",", ":"), default=str)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def _round_robin_buckets(rows, bucket_fn, limit, score_fn):
    """Take one candidate per bucket per pass, then repeat.

    This is deterministic and guarantees broad representation before any bucket can
    consume its second, third, ... slot.  ``rows`` should already be deduplicated.
    """
    buckets = defaultdict(list)
    for row in rows:
        buckets[bucket_fn(row)].append(row)
    queues = {}
    for key, vals in buckets.items():
        queues[key] = deque(sorted(vals, key=lambda x: (score_fn(x), str(x))))
    keys = sorted(queues, key=str)
    selected = []
    while len(selected) < limit and any(queues[k] for k in keys):
        for key in keys:
            if queues[key] and len(selected) < limit:
                selected.append(queues[key].popleft())
    return selected, buckets


def retain_pairs_diverse(pairs, limit=24, region_size=120, amplicon_bin=25, per_near=2):
    """Retain primer pairs across target regions and product-length bands.

    The old path simply used ``pairs[:limit]`` after a Tm/dimer preliminary sort.  A
    dense local region could therefore prevent distant or different-length pairs from
    ever receiving probe search.  This function keeps the preliminary rank as a
    within-stratum ordering while distributing the beam across region/length strata.
    """
    limit = max(1, int(limit)); entered = len(pairs)
    exact = {}
    duplicates = []
    for idx, pair in enumerate(pairs):
        key = pair_identity_key(pair)
        if key not in exact or _pair_prelim(pair) < _pair_prelim(exact[key][1]):
            if key in exact:
                duplicates.append(exact[key])
            exact[key] = (idx, pair)
        else:
            duplicates.append((idx, pair))
    unique = [x[1] for x in exact.values()]
    unique.sort(key=lambda x: (_pair_prelim(x), pair_identity_key(x)))

    near_counts = defaultdict(int); near_kept = []; near_rejected = []
    for pair in unique:
        nk = (int(pair.get("fstart", 0)) // 10,
              int(pair.get("rend", 0)) // 10,
              int(pair.get("amp", 0)) // max(1, int(amplicon_bin)))
        if near_counts[nk] >= max(1, int(per_near)):
            near_rejected.append(pair); continue
        near_counts[nk] += 1; near_kept.append(pair)

    def bucket(pair):
        return (int(pair.get("fstart", 0)) // max(1, int(region_size)),
                int(pair.get("amp", 0)) // max(1, int(amplicon_bin)))
    # Guarantee at least one representative from every available region before
    # spending remaining slots on finer region/amplicon strata.
    by_region = defaultdict(list)
    for pair in near_kept:
        by_region[bucket(pair)[0]].append(pair)
    selected=[]; selected_ids=set()
    for reg in sorted(by_region):
        if len(selected) >= limit: break
        pair=min(by_region[reg], key=lambda x:(_pair_prelim(x), pair_identity_key(x)))
        selected.append(pair); selected_ids.add(pair_identity_key(pair))
    remaining=[x for x in near_kept if pair_identity_key(x) not in selected_ids]
    extra, buckets = _round_robin_buckets(remaining, bucket, max(0,limit-len(selected)), _pair_prelim)
    for pair in extra:
        if pair_identity_key(pair) not in selected_ids:
            selected.append(pair); selected_ids.add(pair_identity_key(pair))
    budget_rejected = [x for x in near_kept if pair_identity_key(x) not in selected_ids]

    decisions = []
    for idx, pair in duplicates:
        decisions.append(dict(candidate_id=_stable_id("pair", pair_identity_key(pair), idx),
                              identity=list(pair_identity_key(pair)), decision="rejected",
                              reason="exact_duplicate", preliminary_rank=_pair_prelim(pair)))
    for pair in near_rejected:
        decisions.append(dict(candidate_id=_stable_id("pair", pair_identity_key(pair)),
                              identity=list(pair_identity_key(pair)), decision="rejected",
                              reason="near_duplicate_budget", preliminary_rank=_pair_prelim(pair)))
    for pair in budget_rejected:
        decisions.append(dict(candidate_id=_stable_id("pair", pair_identity_key(pair)),
                              identity=list(pair_identity_key(pair)), decision="rejected",
                              reason="diversity_beam_budget", preliminary_rank=_pair_prelim(pair),
                              stratum=list(bucket(pair))))
    for pair in selected:
        decisions.append(dict(candidate_id=_stable_id("pair", pair_identity_key(pair)),
                              identity=list(pair_identity_key(pair)), decision="retained",
                              reason="diversity_beam", preliminary_rank=_pair_prelim(pair),
                              stratum=list(bucket(pair))))
    decisions.sort(key=lambda x: (x["decision"] != "retained", x["preliminary_rank"], x["candidate_id"]))
    ledger = dict(version=RETENTION_VERSION, stage="pair_diversity_retention",
                  entered=entered, retained=len(selected), rejected=entered-len(selected),
                  reasons=dict(exact_duplicate=len(duplicates), near_duplicate_budget=len(near_rejected),
                               diversity_beam_budget=len(budget_rejected)),
                  hard_gate=False, reversible=True, candidate_limit=limit,
                  region_size=region_size, amplicon_bin=amplicon_bin,
                  strata_available=len({bucket(x) for x in near_kept}),
                  regions_retained=len({bucket(x)[0] for x in selected}),
                  amplicon_bands_retained=len({bucket(x)[1] for x in selected}),
                  candidate_decisions=decisions)
    return selected, ledger


def retain_probes_diverse(probes, limit=4, position_bin=8, offset_bin=1.5):
    """Retain probe alternatives across position, strand and Tm-offset trade-offs."""
    limit=max(1,int(limit)); entered=len(probes)
    exact={}; duplicates=[]
    for idx, probe in enumerate(probes):
        key=probe_identity_key(probe)
        if key not in exact or _probe_prelim(probe) < _probe_prelim(exact[key][1]):
            if key in exact: duplicates.append(exact[key])
            exact[key]=(idx,probe)
        else: duplicates.append((idx,probe))
    unique=[x[1] for x in exact.values()]
    def bucket(p):
        return (str(p.get("strand") or "+"), int(p.get("start",0))//max(1,int(position_bin)),
                int(round(float(p.get("offset",0.0))/max(0.1,float(offset_bin)))))
    # Preserve both probe orientations when available before distributing the
    # remaining beam across position/Tm-offset strata.
    by_strand=defaultdict(list)
    for probe in unique: by_strand[str(probe.get("strand") or "+")].append(probe)
    selected=[]; selected_ids=set()
    # When the beam is smaller than the number of available strands, select the
    # strongest strand representative first.  The previous lexical ``+`` then
    # ``-`` order could discard the globally best probe whenever ``limit == 1``.
    strand_representatives=[]
    for strand, rows in by_strand.items():
        probe=min(rows,key=lambda x:(_probe_prelim(x),probe_identity_key(x)))
        strand_representatives.append((probe,strand))
    strand_representatives.sort(key=lambda x:(_probe_prelim(x[0]),probe_identity_key(x[0]),x[1]))
    for probe, strand in strand_representatives:
        if len(selected)>=limit: break
        selected.append(probe); selected_ids.add(probe_identity_key(probe))
    remaining=[x for x in unique if probe_identity_key(x) not in selected_ids]
    extra,buckets=_round_robin_buckets(remaining,bucket,max(0,limit-len(selected)),_probe_prelim)
    for probe in extra:
        if probe_identity_key(probe) not in selected_ids:
            selected.append(probe); selected_ids.add(probe_identity_key(probe))
    budget=[x for x in unique if probe_identity_key(x) not in selected_ids]
    decisions=[]
    for idx,p in duplicates:
        decisions.append(dict(candidate_id=_stable_id("probe",probe_identity_key(p),idx),
                              identity=list(probe_identity_key(p)),decision="rejected",
                              reason="exact_duplicate",preliminary_rank=_probe_prelim(p)))
    for p in budget:
        decisions.append(dict(candidate_id=_stable_id("probe",probe_identity_key(p)),
                              identity=list(probe_identity_key(p)),decision="rejected",
                              reason="diversity_beam_budget",preliminary_rank=_probe_prelim(p),
                              stratum=list(bucket(p))))
    for p in selected:
        decisions.append(dict(candidate_id=_stable_id("probe",probe_identity_key(p)),
                              identity=list(probe_identity_key(p)),decision="retained",
                              reason="diversity_beam",preliminary_rank=_probe_prelim(p),
                              stratum=list(bucket(p))))
    decisions.sort(key=lambda x:(x["decision"]!="retained",x["preliminary_rank"],x["candidate_id"]))
    ledger=dict(version=RETENTION_VERSION,stage="probe_diversity_retention",
                entered=entered,retained=len(selected),rejected=entered-len(selected),
                reasons=dict(exact_duplicate=len(duplicates),diversity_beam_budget=len(budget)),
                hard_gate=False,reversible=True,candidate_limit=limit,
                position_bin=position_bin,offset_bin=offset_bin,strata_available=len({bucket(x) for x in unique}),
                positions_retained=len({bucket(x)[1] for x in selected}),
                strands_retained=sorted({str(x.get("strand")) for x in selected}),
                candidate_decisions=decisions)
    return selected,ledger


def retain_diverse(candidates, limit=160, region_size=250, per_region=16, per_near=2):
    """Return a bounded, deterministic, regionally diverse triplet pool plus ledger.

    Lower ``candidate_rank`` is better.  Exact duplicates are removed first, then
    near-identical triplets are capped, then each target region receives a quota.
    Remaining slots are filled globally.  The decision trace is intentionally stored
    because a later rank cannot recover a candidate discarded here.
    """
    limit = max(1, int(limit))
    entered = len(candidates)
    groups = defaultdict(list)
    for idx, assay in enumerate(candidates):
        groups[identity_key(assay)].append((idx, assay))
    exact = {}; duplicate_rows=[]
    for key, rows in groups.items():
        rows=sorted(rows,key=lambda x:(_prelim(x[1]),x[0]))
        exact[key]=rows[0][1]
        duplicate_rows.extend(rows[1:])
    ordered = sorted(exact.values(), key=lambda a: (_prelim(a), a.get("search_window_start", 0), identity_key(a)))

    near_counts = defaultdict(int)
    near_kept, near_removed_rows = [], []
    for a in ordered:
        k = near_key(a)
        if near_counts[k] >= per_near:
            near_removed_rows.append(a)
            continue
        near_counts[k] += 1
        near_kept.append(a)

    regions = defaultdict(list)
    for a in near_kept:
        f = a.get("f_xy") or [a.get("search_window_start", 0), 0]
        regions[int(f[0]) // max(1, int(region_size))].append(a)
    for reg in regions:
        regions[reg].sort(key=lambda a: (_prelim(a), identity_key(a)))

    selected, used = [], set()
    # First guarantee one representative per available region (up to the global
    # limit), ordered by each region's best preliminary candidate rather than by
    # coordinate.  The old coordinate-order loop could fill the budget entirely
    # from the first region when ``per_region >= limit``.
    region_order = sorted(regions, key=lambda reg: (_prelim(regions[reg][0]), reg))
    for reg in region_order:
        if len(selected) >= limit:
            break
        a = regions[reg][0]
        selected.append(a); used.add(identity_key(a))

    # Then spend additional slots round-robin across regions.  No region can take
    # its third candidate before every represented region has had a chance at its
    # second, and ``per_region`` remains an explicit cap.
    depth = 1
    cap = max(1, int(per_region))
    while len(selected) < limit and depth < cap:
        added = False
        for reg in region_order:
            if len(selected) >= limit:
                break
            if depth >= len(regions[reg]):
                continue
            a = regions[reg][depth]; k = identity_key(a)
            if k not in used:
                selected.append(a); used.add(k); added = True
        if not added and all(depth >= len(regions[r]) for r in region_order):
            break
        depth += 1

    # If per-region caps leave spare capacity, fill globally by preliminary rank.
    for a in near_kept:
        k = identity_key(a)
        if len(selected) >= limit:
            break
        if k not in used:
            selected.append(a); used.add(k)

    selected.sort(key=lambda a: (_prelim(a), a.get("search_window_start", 0), identity_key(a)))
    selected_ids={identity_key(a) for a in selected}
    budget_rows=[a for a in near_kept if identity_key(a) not in selected_ids]
    decisions=[]
    for idx,a in duplicate_rows:
        decisions.append(dict(candidate_id=_stable_id("triplet",identity_key(a),idx),
                              identity=list(identity_key(a)),decision="rejected",reason="exact_duplicate",
                              preliminary_rank=_prelim(a),region=int((a.get("f_xy") or [0])[0])//max(1,int(region_size))))
    for a in near_removed_rows:
        decisions.append(dict(candidate_id=_stable_id("triplet",identity_key(a)),
                              identity=list(identity_key(a)),decision="rejected",reason="near_duplicate_budget",
                              preliminary_rank=_prelim(a),region=int((a.get("f_xy") or [0])[0])//max(1,int(region_size))))
    for a in budget_rows:
        decisions.append(dict(candidate_id=_stable_id("triplet",identity_key(a)),
                              identity=list(identity_key(a)),decision="rejected",reason="regional_or_global_budget",
                              preliminary_rank=_prelim(a),region=int((a.get("f_xy") or [0])[0])//max(1,int(region_size))))
    for a in selected:
        decisions.append(dict(candidate_id=_stable_id("triplet",identity_key(a)),
                              identity=list(identity_key(a)),decision="retained",reason="full_annotation_pool",
                              preliminary_rank=_prelim(a),region=int((a.get("f_xy") or [0])[0])//max(1,int(region_size))))
    decisions.sort(key=lambda x:(x["decision"]!="retained",x["preliminary_rank"],x["candidate_id"]))
    ledger = {
        "version": RETENTION_VERSION,
        "stage": "diversity_retention",
        "entered": entered,
        "retained": len(selected),
        "rejected": entered - len(selected),
        "reasons": {
            "exact_duplicate": len(duplicate_rows),
            "near_duplicate_budget": len(near_removed_rows),
            "regional_or_global_budget": len(budget_rows),
        },
        "hard_gate": False,
        "reversible": True,
        "region_size": region_size,
        "per_region": per_region,
        "per_near_family": per_near,
        "candidate_limit": limit,
        "regions_available": len(regions),
        "regions_retained": len({int((a.get("f_xy") or [0])[0]) // max(1, int(region_size)) for a in selected}),
        "selection_policy": "one_per_region_then_round_robin_then_global_fill",
        "candidate_decisions": decisions,
    }
    return selected, ledger
