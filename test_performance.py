"""G9 -- runtime/performance benchmark regression (pins bench_performance.json).

Wall-clock ms are host-dependent, so this test does NOT pin milliseconds. It pins the *portable
shapes* -- the claims a reviewer would recompute on any machine and that must not silently drift:

  1. Reads the committed bench_performance.json and pins its structural findings (OF slower than
     Primer3 on every target; the slowdown is length-dominated; the specificity scan is linear;
     the GC extreme is clean; Tm is sub-millisecond; a whole human genome is hours-scale).
  2. RE-COMPUTES portable implementation facts live:
       * specificity scan work units are affine in subject length (R^2 > 0.999),
       * the frozen extreme templates remain present and Primer3 accepts them,
       * per-oligo Tm is sub-millisecond.
     Relative wall-clock ratios remain in the frozen, environment-recorded artifact. They are not
     live pass/fail gates because scheduler noise in a single fast Primer3 call is not a scientific
     invariant.

Offline; needs primer3-py (a hard dependency of the tool) for the direction checks, and degrades
to a skip-pass on those two if it is somehow absent. Run: PYTHONPATH=. python3 tests/test_performance.py
"""
import os
import sys
import json
import time
import random

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.benchmark import bench_performance as BP  # noqa: E402
from oligoforge import thermo as T, design as D, profiles as P, specificity as S  # noqa: E402

BENCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "benchmark")

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + ("" if cond else f"   [{detail}]"))
    if not cond:
        fails.append(name)


# ---- 1. committed artifact: structural findings -------------------------------------------------
s = json.load(open(os.path.join(BENCH, "bench_performance.json")))
check("schema + env recorded", s.get("schema") == "oligoforge/bench_performance/1"
      and bool(s["env"].get("python")) and s["env"].get("cpu_count") is not None,
      str(s.get("schema")))

dv = s["design_vs_primer3"]
check("committed: OligoForge slower than Primer3 on ALL targets",
      dv["have_primer3"] and dv["of_slower_than_primer3_on_all"] is True)
check("committed: slowdown is a large multiple (min ratio > 10, honest big gap)",
      dv["have_primer3"] and dv["of_over_p3_min"] > 10.0, str(dv["of_over_p3_min"]))
check("committed: every corpus design succeeded (no crash across GC 27-62%)",
      dv["all_designs_succeeded"] is True)

dc = dv["driver_correlations"]
check("committed: length is the dominant driver (r(len,OF ms) >= r(GC,ratio))",
      dc["r_len_vs_of_ms"] is not None and dc["r_gc_vs_ratio"] is not None
      and dc["r_len_vs_of_ms"] >= dc["r_gc_vs_ratio"],
      f"r_len_ms={dc['r_len_vs_of_ms']} r_gc_ratio={dc['r_gc_vs_ratio']}")
ex = dv["gc_extreme"]
check("committed: GC-rich extreme slower than AT-rich extreme (clean extreme, not monotone law)",
      ex.get("GCrich_ratio_exceeds_ATrich") is True
      and ex["most_GC_rich"]["ratio"] > ex["most_AT_rich"]["ratio"],
      f"{ex['most_AT_rich']['ratio']} vs {ex['most_GC_rich']['ratio']}")

sc = s["specificity_scan_scaling"]
check("committed: specificity scan linear (recorded R^2 > 0.99)",
      sc["is_linear_r2_over_0_99"] is True and sc["linear_fit"]["r2"] > 0.99,
      str(sc["linear_fit"]["r2"]))
check("committed: whole human 3.2 Gb genome is hours-scale (linear projection > 1 h)",
      sc["linear_fit"]["projected_human_3.2Gb_s"] > 3600.0
      and sc["linear_fit"]["projected_human_3.2Gb_s"]
      > sc["linear_fit"]["projected_E_coli_4.6Mb_s"],
      str(sc["linear_fit"]["projected_human_3.2Gb_s"]))

tm = s["tm_throughput"]
check("committed: Tm is sub-millisecond per oligo, warm faster than cold",
      tm["cold_us_per_oligo"] < 1000.0 and tm["warm_us_per_oligo"] < tm["cold_us_per_oligo"],
      f"cold={tm['cold_us_per_oligo']} warm={tm['warm_us_per_oligo']}")


# ---- 2. live recompute of portable implementation shapes ----------------------------------------
# (a) Count scan work units rather than timing them.  scan_primer_sites calls _match_at exactly
# twice per candidate window (plus and minus orientation), so this pins linear algorithmic work
# without turning host scheduling noise into a release failure.
fwd, rev = "ACGTGACCTGACTGATCAGT", "TGACTGATCAGTCAGGTCACG"
xs, ys = [], []
original_match = S._match_at
counter = {"n": 0}
def _count_match(*args, **kwargs):
    counter["n"] += 1
    return 99, False, True
try:
    S._match_at = _count_match
    for length in (1000, 2000, 4000, 8000):
        counter["n"] = 0
        S.scan_primer_sites(fwd, [("g", "A" * length)], max_mm=2)
        xs.append(length)
        ys.append(counter["n"])
finally:
    S._match_at = original_match
_, _, r2_live = BP._linfit(xs, ys)
check("live: specificity scan work is linear in subject length (R^2 > 0.999)", r2_live > 0.999,
      f"R^2={r2_live:.5f}")

# (b) The extreme templates remain in the frozen corpus. OligoForge's full design success is
# already pinned in the committed rows and the scientific/golden suites; repeating two expensive
# exhaustive designs here would duplicate those gates. Relative timing is artifact-only.
corpus = {t["id"]: t for t in json.load(open(os.path.join(BENCH, "bench_corpus.json")))["targets"]}
at_t = corpus["plas_cytb_ATrich"]
gc_t = corpus["Mtb_rpoB_GCrich"]
check("live: both frozen extreme templates remain available",
      bool(at_t.get("seq")) and bool(gc_t.get("seq")) and at_t["gc"] < gc_t["gc"])

try:
    import primer3
    have_p3 = True
except Exception:
    have_p3 = False

if have_p3:
    def _p3_result(t):
        seq = t["seq"].upper()
        cfg = {"PRIMER_NUM_RETURN": 5, "PRIMER_MIN_SIZE": 18, "PRIMER_OPT_SIZE": 20,
               "PRIMER_MAX_SIZE": 24, "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]]}
        return primer3.design_primers({"SEQUENCE_TEMPLATE": seq}, cfg)
    at_p3 = _p3_result(at_t)
    gc_p3 = _p3_result(gc_t)
    check("live: Primer3 returns result records for both extreme templates",
          isinstance(at_p3, dict) and isinstance(gc_p3, dict))
else:
    check("live: (primer3 absent) comparison smoke skipped -- primer3 is a hard dependency", True,
          "skip-pass")

# (c) Tm sub-millisecond, live
rng = random.Random(7)
seqs = ["".join(rng.choice("ACGT") for _ in range(22)) for _ in range(200)]
t0 = time.perf_counter()
for q in seqs:
    T.tm_acc(q)
us = (time.perf_counter() - t0) / len(seqs) * 1e6
check("live: per-oligo Tm sub-millisecond (< 1000 us, large margin)", us < 1000.0, f"{us:.1f} us")


print(f"\n{'PASS' if not fails else 'FAIL'} test_performance.py "
      f"({'all checks green' if not fails else str(len(fails)) + ' failed: ' + ', '.join(fails)})")
sys.exit(1 if fails else 0)
