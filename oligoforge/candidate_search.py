"""Target-wide, diversity-preserving search for complete assay triplets."""
import time
from . import design as D
from .candidate_retention import retain_diverse, identity_key

SEARCH_VERSION = "2.1.1"


def _spread_order(starts):
    starts = sorted(set(int(x) for x in starts))
    out, seen = [], set()
    def visit(lo, hi):
        if lo > hi:
            return
        mid = (lo + hi) // 2
        for idx in (lo, hi, mid):
            if idx not in seen:
                seen.add(idx); out.append(starts[idx])
        visit(lo + 1, mid - 1)
        visit(mid + 1, hi - 1)
    if starts:
        visit(0, len(starts) - 1)
    return out


def _starts(length, window, max_windows=18, step=140):
    if length <= window:
        return [0]
    span = length - window
    even = [round(i * span / (max_windows - 1)) for i in range(max_windows)]
    grid = list(range(0, span + 1, max(1, int(step)))) + [span]
    merged = sorted(set(even + grid))
    if len(merged) > max_windows:
        merged = [merged[round(i * (len(merged) - 1) / (max_windows - 1))]
                  for i in range(max_windows)]
    return _spread_order(merged)


def search(reference, profile, *, window=420, step=140, budget_s=35.0,
           pair_limit=12, probes_per_pair=3, triplets_per_window=30,
           retained_limit=180, max_windows=18):
    """Search a reference broadly enough to preserve strong complete triplets.

    The result is heuristic-bounded, not an exact global optimum.  Every limit and
    attrition count is returned in a machine-readable ledger.
    """
    ref = (reference or "").upper()
    window = max(int(window), min(int(profile.get("amp_max", 150)) + 140, 2200))
    starts = _starts(len(ref), window, max_windows=max_windows, step=step)
    t0 = time.monotonic()
    raw, window_rows = [], []
    expired = False
    for wi, start in enumerate(starts):
        # Always evaluate the first three target-spanning work units.  This is
        # a deterministic minimum (5', 3', midpoint in the spread ordering),
        # so cold/warm cache state cannot change the minimal candidate corpus.
        # Larger searches retain the declared soft runtime ceiling.
        if wi >= 3 and time.monotonic() - t0 >= float(budget_s):
            expired = True
            break
        sub = ref[start:start + window]
        try:
            rows, led = D.generate_assay_candidates(
                sub, profile, pair_limit=pair_limit,
                probes_per_pair=probes_per_pair,
                triplet_limit=triplets_per_window)
        except Exception as exc:
            rows, led = [], {"stage": "window_joint_triplet_search", "error": type(exc).__name__}
        for a in rows:
            a["search_window_start"] = start
            if a.get("f_xy"):
                a["f_xy"] = [a["f_xy"][0] + start, a["f_xy"][1] + start]
            if a.get("r_xy"):
                a["r_xy"] = [a["r_xy"][0] + start, a["r_xy"][1] + start]
            if a.get("probe_xy"):
                a["probe_xy"] = [a["probe_xy"][0] + start, a["probe_xy"][1] + start]
            if a.get("gblock_span"):
                a["gblock_span"] = [a["gblock_span"][0] + start, a["gblock_span"][1] + start]
            if a.get("f_xy") and a.get("r_xy"):
                a["amplicon_xy"] = [a["f_xy"][0], a["r_xy"][1]]
            raw.append(a)
        led = dict(led)
        led.update(window_index=wi, window_start=start, window_end=min(len(ref), start + window),
                   retained_from_window=len(rows))
        window_rows.append(led)
    retained, retention_ledger = retain_diverse(raw, limit=retained_limit)
    # Aggregate stages use consistent units.  Earlier ledgers mixed primer-pair
    # counts with triplet counts, which made apparent attrition impossible to audit.
    primer_entered = sum((x.get("primer_attrition") or {}).get("entered", 0) for x in window_rows)
    primer_retained = sum((x.get("primer_attrition") or {}).get("retained", 0) for x in window_rows)
    pair_cartesian = sum((x.get("pair_attrition") or {}).get("entered", 0) for x in window_rows)
    pair_cheap = sum((x.get("pair_attrition") or {}).get("cheap_gate_survivors", 0) for x in window_rows)
    pair_dimer_retained = sum(x.get("pairs_after_hard_gates", 0) for x in window_rows)
    pair_explored = sum(x.get("pairs_fully_explored", 0) for x in window_rows)
    probe_entered = sum(sum(p.get("entered", 0) for p in (x.get("probe_attrition") or [])) for x in window_rows)
    probe_hard_survivors = sum(sum(p.get("hard_gate_survivors", 0) for p in (x.get("probe_attrition") or [])) for x in window_rows)
    probe_retained = sum(sum(p.get("retained", 0) for p in (x.get("probe_attrition") or [])) for x in window_rows)
    triplets_generated = sum(x.get("triplets_generated", 0) for x in window_rows)
    triplets_window_retained = sum(x.get("triplets_retained", 0) for x in window_rows)
    stages = [
        {"stage": "window_sampling", "unit": "windows", "entered": len(starts), "retained": len(window_rows),
         "rejected": len(starts) - len(window_rows), "hard_gate": False, "reversible": True,
         "reason": "runtime budget" if expired else "none"},
        {"stage": "primer_hard_screen", "unit": "oriented_oligo_windows", "entered": primer_entered,
         "retained": primer_retained, "rejected": primer_entered-primer_retained,
         "hard_gate": True, "reversible": False},
        {"stage": "pair_geometry_tm_screen", "unit": "primer_pairs", "entered": pair_cartesian,
         "retained": pair_cheap, "rejected": pair_cartesian-pair_cheap,
         "hard_gate": True, "reversible": False},
        {"stage": "pair_dimer_screen_and_native_cap", "unit": "primer_pairs", "entered": pair_cheap,
         "retained": pair_dimer_retained, "rejected": pair_cheap-pair_dimer_retained,
         "hard_gate": "mixed", "reversible": True,
         "reason": "annealing-temperature dimer failures plus documented native pair cap"},
        {"stage": "pair_diversity_beam", "unit": "primer_pairs", "entered": pair_dimer_retained,
         "retained": pair_explored, "rejected": pair_dimer_retained-pair_explored,
         "hard_gate": False, "reversible": True},
    ]
    if not profile.get("no_probe"):
        stages.extend([
            {"stage": "probe_hard_screen", "unit": "probe_windows", "entered": probe_entered,
             "retained": probe_hard_survivors, "rejected": probe_entered-probe_hard_survivors,
             "hard_gate": True, "reversible": False},
            {"stage": "probe_diversity_beam", "unit": "probe_candidates", "entered": probe_hard_survivors,
             "retained": probe_retained, "rejected": probe_hard_survivors-probe_retained,
             "hard_gate": False, "reversible": True},
        ])
    stages.extend([
        {"stage": "window_triplet_beam", "unit": "complete_assays", "entered": triplets_generated,
         "retained": triplets_window_retained, "rejected": triplets_generated-triplets_window_retained,
         "hard_gate": False, "reversible": True},
        retention_ledger,
    ])
    ledger = {
        "schema_version": "1.1.0",
        "search_version": SEARCH_VERSION,
        "status": "heuristic_bounded",
        "reference_length": len(ref),
        "window_size": window,
        "windows_planned": len(starts),
        "windows_evaluated": len(window_rows),
        "runtime_budget_seconds": float(budget_s),
        "runtime_budget_expired": expired,
        "deterministic_minimum_windows": min(3, len(starts)),
        "candidate_limits": {
            "pair_limit_per_window": pair_limit,
            "probes_per_pair": probes_per_pair,
            "triplets_per_window": triplets_per_window,
            "full_annotation_pool": retained_limit,
        },
        "stages": stages,
        "window_ledgers": window_rows,
        "raw_triplets": len(raw),
        "retained_for_full_annotation": len(retained),
        "unique_regions": len({(a.get("f_xy") or [0])[0] // 250 for a in retained}),
    }
    # Stable fingerprint for traces/tests.
    ledger["retained_identities"] = ["|".join(x or "" for x in identity_key(a)) for a in retained]
    return retained, ledger
