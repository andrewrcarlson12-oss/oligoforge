"""G9 -- runtime/performance benchmark regression (pins bench_performance.json).

Wall-clock ms are host-dependent, so this test does NOT pin milliseconds. It pins the *portable
shapes* -- the claims a reviewer would recompute on any machine and that must not silently drift:

  1. Reads the committed bench_performance.json and pins its structural findings (OF slower than
     Primer3 on every target; the slowdown is length-dominated; the specificity scan is linear;
     the GC extreme is clean; Tm is sub-millisecond; a whole human genome is hours-scale).
  2. RE-COMPUTES the hardware-robust facts live and asserts they still hold:
       * the offline specificity scan is linear in subject length (R^2 > 0.99),
       * OligoForge design is SLOWER than Primer3's C core on both GC extremes (the honest
         direction -- guards against a future change spuriously claiming parity/superiority),
       * the most-GC-rich extreme is slower (higher ratio) than the most-AT-rich extreme,
       * per-oligo Tm is sub-millisecond.
     These hold by large margins regardless of host, so the test is deterministic, not flaky.

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


# ---- 2. live recompute of the hardware-robust shapes --------------------------------------------
# (a) specificity scan linearity on small sizes (fast; the shape holds regardless of host speed)
fix = json.load(open(os.path.join(BENCH, "specificity_fixture.json")))["sequences"]
base = "".join(fix.values())
fwd, rev = "ACGTGACCTGACTGATCAGT", "TGACTGATCAGTCAGGTCACG"
xs, ys = [], []
for mult in (1, 2, 4, 8):
    fasta = ">g\n" + base * mult + "\n"
    t0 = time.perf_counter()
    S.in_silico_pcr_offline(fwd, rev, fasta, max_mm=2)
    xs.append(len(base) * mult)
    ys.append((time.perf_counter() - t0) * 1000.0)
_, _, r2_live = BP._linfit(xs, ys)
check("live: specificity scan is linear in subject length (R^2 > 0.99)", r2_live > 0.99,
      f"R^2={r2_live:.5f}")

# (b) design direction vs Primer3 on the two GC extremes
corpus = {t["id"]: t for t in json.load(open(os.path.join(BENCH, "bench_corpus.json")))["targets"]}
at_t = corpus["plas_cytb_ATrich"]
gc_t = corpus["Mtb_rpoB_GCrich"]

def _of_ms(t):
    prof = P.PROFILES[t["profile"]]
    seq = t["seq"].upper()
    return BP._cold_ms(lambda: D.design_assay(seq, prof), t.get("anneal_c", 60)), \
        D.design_assay(seq, prof)

at_ms, at_res = _of_ms(at_t)
gc_ms, gc_res = _of_ms(gc_t)
T.set_conditions(anneal_c=60)
check("live: both extreme designs succeed", at_res is not None and gc_res is not None)

try:
    import primer3
    have_p3 = True
except Exception:
    have_p3 = False

if have_p3:
    def _p3_ms(t):
        seq = t["seq"].upper()
        # warm the C core once, then time a single call (fast; direction not ms is the claim)
        cfg = {"PRIMER_NUM_RETURN": 5, "PRIMER_MIN_SIZE": 18, "PRIMER_OPT_SIZE": 20,
               "PRIMER_MAX_SIZE": 24, "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]]}
        primer3.design_primers({"SEQUENCE_TEMPLATE": seq}, cfg)
        t0 = time.perf_counter()
        primer3.design_primers({"SEQUENCE_TEMPLATE": seq}, cfg)
        return (time.perf_counter() - t0) * 1000.0
    at_p3 = _p3_ms(at_t)
    gc_p3 = _p3_ms(gc_t)
    at_ratio = at_ms / at_p3
    gc_ratio = gc_ms / gc_p3
    check("live: OligoForge slower than Primer3 on AT-rich extreme (ratio > 5, honest direction)",
          at_ratio > 5.0, f"{at_ratio:.1f}x")
    check("live: OligoForge slower than Primer3 on GC-rich extreme (ratio > 5, honest direction)",
          gc_ratio > 5.0, f"{gc_ratio:.1f}x")
    check("live: GC-rich extreme slower than AT-rich extreme (clean extreme reproduces)",
          gc_ratio > at_ratio, f"AT {at_ratio:.1f}x vs GC {gc_ratio:.1f}x")
else:
    check("live: (primer3 absent) direction checks skipped -- primer3 is a hard dependency", True,
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
