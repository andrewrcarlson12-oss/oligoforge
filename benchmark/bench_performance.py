"""Performance / runtime benchmark (G9).

Measures OligoForge's runtime on its real code paths and, honestly, against Primer3's C core.
Writes ``bench_performance.json`` + ``bench_performance.csv`` and (if matplotlib is present) the
figure ``performance_benchmark.png``.

HONESTY (this is the point of the benchmark, per the release brief -- a loss reported plainly is
worth more than an engineered win):

  * OligoForge's exhaustive *Python* candidate enumeration is **1-2 orders of magnitude slower**
    than Primer3's C core on the same template, and the gap **grows with template GC** because
    the anneal-temperature structure gate (differentiator D2) admits *more* GC-rich candidates,
    which then cost more to score. That is the compute price of the D2 design decision. It is not
    hidden: it is the headline of this benchmark.
  * The tool's use case is **interactive single-assay design**, where OligoForge's absolute cost
    (a fraction of a second to a few seconds per assay) is fine. It is **not** a batch/genome-scale
    design engine, and this benchmark says so with numbers.
  * The offline in-silico-PCR specificity scan is **linear** in subject length (~10 us/bp), which
    bounds it: fine for transcriptomes and small genomes, impractical for a whole 3.2 Gb human
    genome -- which is exactly why the tool checks specificity against a *supplied* FASTA rather
    than claiming a genome-wide guarantee (D3).

REPRODUCIBILITY: absolute wall-clock depends on the host (recorded under ``env`` in the JSON).
The *portable* claims -- the ones ``test_performance.py`` re-checks on any machine -- are SHAPES,
not milliseconds: scan linearity (R^2), OF-slower-than-Primer3 (direction), ratio-grows-with-GC
(monotonic direction), and sub-millisecond Tm. Those hold by large margins regardless of hardware.

Offline and self-contained. Run:  PYTHONPATH=. python3 tests/benchmark/bench_performance.py
"""
import os
import sys
import json
import csv
import time
import platform
import statistics
import random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from oligoforge import thermo as T, design as D, profiles as P, specificity as S  # noqa: E402

CORPUS = json.load(open(os.path.join(HERE, "bench_corpus.json")))["targets"]
SPEC_FIXTURE = json.load(open(os.path.join(HERE, "specificity_fixture.json")))["sequences"]

# A fixed, well-formed primer pair used only to drive the specificity SCAN timing (its biological
# meaning is irrelevant here -- we are timing the scan cost per base, not making a specificity call).
_SCAN_FWD = "ACGTGACCTGACTGATCAGT"
_SCAN_REV = "TGACTGATCAGTCAGGTCACG"


def _median_ms(fn, n=3):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts) * 1000.0


def _cold_ms(fn, anneal_c):
    """One cold call: clear the thermo caches first (set_conditions clears them) so we measure the
    real first-time-a-user-designs cost, not a warm-cache repeat."""
    T.set_conditions(anneal_c=anneal_c)   # also clears all snapshot-keyed thermo caches
    t0 = time.perf_counter()
    fn()
    return (time.perf_counter() - t0) * 1000.0


def _linfit(xs, ys):
    """Least-squares slope/intercept + R^2 for y = slope*x + intercept."""
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    ybar = sy / n
    sstot = sum((y - ybar) ** 2 for y in ys)
    ssres = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ssres / sstot if sstot > 0 else 1.0
    return slope, intercept, r2


# --------------------------------------------------------------------------------------------------
# B1 + B2 -- single-assay design latency, and the honest head-to-head vs Primer3's C core
# --------------------------------------------------------------------------------------------------
def bench_design_vs_primer3():
    try:
        import primer3
        have_p3 = True
    except Exception:
        have_p3 = False

    rows = []
    for t in CORPUS:
        prof = P.PROFILES[t["profile"]]
        seq = t["seq"].upper()
        anneal = t.get("anneal_c", 60)
        of_ms = _cold_ms(lambda: D.design_assay(seq, prof), anneal)
        succeeded = D.design_assay(seq, prof) is not None  # warm; correctness, not timing

        p3_ms = None
        ratio = None
        if have_p3:
            def _p3():
                primer3.design_primers(
                    {"SEQUENCE_TEMPLATE": seq},
                    {"PRIMER_NUM_RETURN": 5, "PRIMER_MIN_SIZE": 18, "PRIMER_OPT_SIZE": 20,
                     "PRIMER_MAX_SIZE": 24, "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]]})
            p3_ms = _median_ms(_p3, n=3)
            ratio = of_ms / p3_ms if p3_ms > 0 else None

        rows.append({
            "id": t["id"], "organism": t.get("organism", ""), "len": len(seq),
            "gc": float(t["gc"]), "profile": t["profile"],
            "of_design_ms": round(of_ms, 1),
            "primer3_ms": round(p3_ms, 1) if p3_ms is not None else None,
            "of_over_p3": round(ratio, 1) if ratio is not None else None,
            "design_succeeded": succeeded,
        })
    T.set_conditions(anneal_c=60)  # restore session default
    return rows, have_p3


# --------------------------------------------------------------------------------------------------
# B3 -- Tm throughput (the hottest inner-loop primitive)
# --------------------------------------------------------------------------------------------------
def bench_tm_throughput(n=400):
    rng = random.Random(1998)  # deterministic corpus of oligos (seed => reproducible set)
    seqs = ["".join(rng.choice("ACGT") for _ in range(22)) for _ in range(n)]
    # cold: every seq distinct => no cache reuse (worst case = real first-touch cost)
    t0 = time.perf_counter()
    for s in seqs:
        T.tm_acc(s)
    cold_us = (time.perf_counter() - t0) / n * 1e6
    # warm: same seq repeated => lru_cache hit (best case)
    s0 = seqs[0]
    t0 = time.perf_counter()
    for _ in range(n):
        T.tm_acc(s0)
    warm_us = (time.perf_counter() - t0) / n * 1e6
    return {"n": n, "cold_us_per_oligo": round(cold_us, 2), "warm_us_per_oligo": round(warm_us, 3)}


# --------------------------------------------------------------------------------------------------
# B4 -- offline in-silico-PCR specificity scan: scaling in subject length
# --------------------------------------------------------------------------------------------------
def bench_specificity_scaling():
    base = "".join(SPEC_FIXTURE.values())          # ~15.8 kb of real gene sequence, tiled
    tile = len(base)
    points = []
    for mult in (1, 2, 5, 10, 20, 40):
        fasta = ">g\n" + base * mult + "\n"
        bp = tile * mult
        ms = _median_ms(lambda: S.in_silico_pcr_offline(_SCAN_FWD, _SCAN_REV, fasta, max_mm=2), n=3)
        points.append({"bp": bp, "ms": round(ms, 1)})
    xs = [p["bp"] for p in points]
    ys = [p["ms"] for p in points]
    slope, intercept, r2 = _linfit(xs, ys)         # slope in ms/bp
    us_per_bp = slope * 1000.0
    # projected envelope (LINEAR extrapolation -- labelled as such, not measured)
    def proj_seconds(bp):
        return round((slope * bp + intercept) / 1000.0, 1)
    envelope = {
        "measured_max_bp": max(xs),
        "us_per_bp": round(us_per_bp, 3),
        "r2": round(r2, 5),
        "projected_E_coli_4.6Mb_s": proj_seconds(4_600_000),
        "projected_250Mb_chromosome_s": proj_seconds(250_000_000),
        "projected_human_3.2Gb_s": proj_seconds(3_200_000_000),
    }
    return points, envelope


def _env():
    return {
        "python": platform.python_version(),
        "system": platform.system(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "note": ("Single-host wall-clock; absolute ms are host-dependent. Portable claims are the "
                 "SHAPES (linearity, OF>P3 direction, ratio-vs-GC monotonicity, sub-ms Tm), which "
                 "test_performance.py re-checks live on any machine."),
    }


def _pearson(a, b):
    n = len(a)
    if n < 2:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    va = sum((x - ma) ** 2 for x in a) ** 0.5
    vb = sum((y - mb) ** 2 for y in b) ** 0.5
    return cov / (va * vb) if va > 0 and vb > 0 else None


def run():
    design_rows, have_p3 = bench_design_vs_primer3()
    tm = bench_tm_throughput()
    scan_points, scan_env = bench_specificity_scaling()

    of_times = [r["of_design_ms"] for r in design_rows]

    # HONEST driver analysis. The OF/Primer3 slowdown is driven mainly by template LENGTH (more
    # windows to enumerate) and by candidate COUNT; GC contributes because the anneal gate (D2)
    # admits more GC-rich candidates -- a mechanism already quantified independently in
    # bench_gate_impact.csv (Mtb rpoB GC-rich +49.3% admitted vs AT-rich Plasmodium +0.0%). We do
    # NOT claim the ratio rises monotonically with GC (within a fixed length band the mid-GC points
    # are noisy); we report the correlations plainly and pin only the robust extreme.
    r_gc_ratio = r_len_ratio = r_len_ofms = None
    extreme = {}
    if have_p3:
        gcs = [r["gc"] for r in design_rows]
        lens = [r["len"] for r in design_rows]
        ratios_all = [r["of_over_p3"] for r in design_rows]
        r_gc_ratio = _pearson(gcs, ratios_all)
        r_len_ratio = _pearson(lens, ratios_all)
        r_len_ofms = _pearson(lens, of_times)
        most_at = min(design_rows, key=lambda r: r["gc"])
        most_gc = max(design_rows, key=lambda r: r["gc"])
        extreme = {
            "most_AT_rich": {"id": most_at["id"], "gc": most_at["gc"], "ratio": most_at["of_over_p3"]},
            "most_GC_rich": {"id": most_gc["id"], "gc": most_gc["gc"], "ratio": most_gc["of_over_p3"]},
            "GCrich_ratio_exceeds_ATrich": most_gc["of_over_p3"] > most_at["of_over_p3"],
        }

    summary = {
        "schema": "oligoforge/bench_performance/1",
        "description": ("Runtime benchmark (G9). OligoForge single-assay design latency and an honest "
                        "head-to-head vs Primer3's C core; Tm throughput; and the offline in-silico-PCR "
                        "specificity scan's scaling in subject length. Absolute ms are host-dependent "
                        "(see env); portable claims are the shapes."),
        "env": _env(),
        "design_vs_primer3": {
            "have_primer3": have_p3,
            "rows": design_rows,
            "of_design_ms_median": round(statistics.median(of_times), 1),
            "of_design_ms_min": round(min(of_times), 1),
            "of_design_ms_max": round(max(of_times), 1),
            "of_slower_than_primer3_on_all": (
                all(r["of_over_p3"] is not None and r["of_over_p3"] > 1.0 for r in design_rows)
                if have_p3 else None),
            "of_over_p3_min": (round(min(r["of_over_p3"] for r in design_rows), 1) if have_p3 else None),
            "of_over_p3_max": (round(max(r["of_over_p3"] for r in design_rows), 1) if have_p3 else None),
            "driver_correlations": {
                "r_gc_vs_ratio": round(r_gc_ratio, 3) if r_gc_ratio is not None else None,
                "r_len_vs_ratio": round(r_len_ratio, 3) if r_len_ratio is not None else None,
                "r_len_vs_of_ms": round(r_len_ofms, 3) if r_len_ofms is not None else None,
                "note": ("Slowdown is driven mainly by template length and candidate count; GC "
                         "contributes weakly-to-moderately via the anneal gate admitting more "
                         "GC-rich candidates (mechanism quantified in bench_gate_impact.csv). NOT "
                         "a monotone GC law -- mid-GC points are noisy; only the extreme is robust."),
            },
            "gc_extreme": extreme,
            "all_designs_succeeded": all(r["design_succeeded"] for r in design_rows),
        },
        "tm_throughput": tm,
        "specificity_scan_scaling": {
            "points": scan_points,
            "linear_fit": scan_env,
            "is_linear_r2_over_0_99": scan_env["r2"] > 0.99,
        },
        "honest_reading": (
            "OligoForge's exhaustive Python design is %sx-%sx slower than Primer3's C core on the "
            "same templates (slower on all %d). The slowdown is driven mainly by template length "
            "(r(len,OF ms)=%s) and candidate count; GC contributes weakly-to-moderately "
            "(r(GC,ratio)=%s) via the anneal gate (D2) admitting more GC-rich candidates -- the "
            "extreme is clean (most-AT-rich %s %sx vs most-GC-rich %s %sx) but it is NOT a monotone "
            "GC law. Absolute design cost (%.0f-%.0f ms on a %d-core %s host) is fine for "
            "interactive single-assay use, the tool's actual scope; it is NOT a batch/genome-scale "
            "design engine. The specificity scan is linear (R^2=%.4f, ~%.1f us/bp): fine for "
            "transcriptomes and small genomes, impractical for a whole 3.2 Gb human genome (~%.0f h "
            "projected) -- which is why specificity is checked against a supplied FASTA, not claimed "
            "genome-wide." % (
                summary_ratio_str(design_rows, have_p3, "min"),
                summary_ratio_str(design_rows, have_p3, "max"),
                len(design_rows),
                round(r_len_ofms, 2) if r_len_ofms is not None else "?",
                round(r_gc_ratio, 2) if r_gc_ratio is not None else "?",
                extreme.get("most_AT_rich", {}).get("id", "?"),
                extreme.get("most_AT_rich", {}).get("ratio", "?"),
                extreme.get("most_GC_rich", {}).get("id", "?"),
                extreme.get("most_GC_rich", {}).get("ratio", "?"),
                min(of_times), max(of_times), os.cpu_count(), platform.machine(),
                scan_env["r2"], scan_env["us_per_bp"],
                scan_env["projected_human_3.2Gb_s"] / 3600.0)
        ),
    }
    return summary


def summary_ratio_str(rows, have_p3, which):
    if not have_p3:
        return "?"
    vals = [r["of_over_p3"] for r in rows if r["of_over_p3"] is not None]
    return str(int(min(vals) if which == "min" else max(vals)))


def write_artifacts(summary):
    # JSON
    with open(os.path.join(HERE, "bench_performance.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    # CSV (per-target design head-to-head)
    with open(os.path.join(HERE, "bench_performance.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "organism", "len", "gc", "profile", "of_design_ms", "primer3_ms",
                    "of_over_p3", "design_succeeded"])
        for r in summary["design_vs_primer3"]["rows"]:
            w.writerow([r["id"], r["organism"], r["len"], r["gc"], r["profile"],
                        r["of_design_ms"], r["primer3_ms"], r["of_over_p3"], r["design_succeeded"]])
    # scan CSV
    with open(os.path.join(HERE, "bench_performance_scan.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bp", "ms"])
        for p in summary["specificity_scan_scaling"]["points"]:
            w.writerow([p["bp"], p["ms"]])


def make_figure(summary):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False

    rows = summary["design_vs_primer3"]["rows"]
    have_p3 = summary["design_vs_primer3"]["have_primer3"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    # Panel A: OF vs Primer3 design time per target (sorted by length), log-y
    rs = sorted(rows, key=lambda r: r["len"])
    x = range(len(rs))
    ax[0].plot(list(x), [r["of_design_ms"] for r in rs], "o-", color="#b5651d", label="OligoForge")
    if have_p3:
        ax[0].plot(list(x), [r["primer3_ms"] for r in rs], "s--", color="#1f3a5f", label="Primer3 (C)")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("target (sorted by template length)")
    ax[0].set_ylabel("design time (ms, log)")
    ax[0].set_title("A. Single-assay design latency")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, which="both", alpha=0.25)

    # Panel B: OF/Primer3 ratio vs template GC (weak positive; length co-drives it -- honest scatter)
    if have_p3:
        gcs = [r["gc"] for r in rows]
        ratios = [r["of_over_p3"] for r in rows]
        ax[1].scatter(gcs, ratios, color="#b5651d", zorder=3)
        # trend line
        s, i, r2 = _linfit(gcs, ratios)
        rp = summary["design_vs_primer3"]["driver_correlations"]["r_gc_vs_ratio"]
        xr = [min(gcs), max(gcs)]
        ax[1].plot(xr, [s * v + i for v in xr], "--", color="#555", lw=1,
                   label="r(GC,ratio)=%.2f" % (rp if rp is not None else 0.0))
        ax[1].set_xlabel("template GC (%)")
        ax[1].set_ylabel("OligoForge / Primer3 (x)")
        ax[1].set_title("B. Slowdown vs GC (weak +; length co-drives)")
        ax[1].legend(fontsize=8)
        ax[1].grid(True, alpha=0.25)

    # Panel C: specificity scan scaling with linear fit + R^2
    pts = summary["specificity_scan_scaling"]["points"]
    fit = summary["specificity_scan_scaling"]["linear_fit"]
    xb = [p["bp"] / 1000.0 for p in pts]     # kb
    yb = [p["ms"] for p in pts]
    ax[2].scatter(xb, yb, color="#1f3a5f", zorder=3)
    s, i, r2 = _linfit([p["bp"] for p in pts], yb)
    ax[2].plot([min(x for x in [p["bp"] for p in pts]) / 1000.0,
                max(p["bp"] for p in pts) / 1000.0],
               [s * min(p["bp"] for p in pts) + i, s * max(p["bp"] for p in pts) + i],
               "--", color="#b5651d", lw=1)
    ax[2].set_xlabel("subject length (kb)")
    ax[2].set_ylabel("scan time (ms)")
    ax[2].set_title("C. Specificity scan is linear\n(%.1f us/bp, R^2=%.4f)"
                    % (fit["us_per_bp"], fit["r2"]))
    ax[2].grid(True, alpha=0.25)

    fig.suptitle("OligoForge runtime (G9) — honest: slower than Primer3's C core, linear specificity scan",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(HERE, "performance_benchmark.png"), dpi=130)
    plt.close(fig)
    return True


if __name__ == "__main__":
    s = run()
    write_artifacts(s)
    ok = make_figure(s)
    print("wrote bench_performance.json / .csv / _scan.csv" + (" + performance_benchmark.png" if ok else " (no figure: matplotlib missing)"))
    dv = s["design_vs_primer3"]
    ex = dv["gc_extreme"]
    print("design: median %.0f ms (%.0f-%.0f); OF/P3 %s-%sx; r(GC,ratio)=%s; AT-rich %sx vs GC-rich %sx"
          % (dv["of_design_ms_median"], dv["of_design_ms_min"], dv["of_design_ms_max"],
             dv["of_over_p3_min"], dv["of_over_p3_max"],
             dv["driver_correlations"]["r_gc_vs_ratio"],
             ex["most_AT_rich"]["ratio"], ex["most_GC_rich"]["ratio"]))
    sc = s["specificity_scan_scaling"]["linear_fit"]
    print("scan: %.1f us/bp, R^2=%.4f; proj E.coli 4.6Mb ~%ss, human 3.2Gb ~%ss"
          % (sc["us_per_bp"], sc["r2"], sc["projected_E_coli_4.6Mb_s"], sc["projected_human_3.2Gb_s"]))
    print("tm: cold %.1f us/oligo, warm %.3f us/oligo" %
          (s["tm_throughput"]["cold_us_per_oligo"], s["tm_throughput"]["warm_us_per_oligo"]))
