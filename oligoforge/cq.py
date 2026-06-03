"""Cq (quantification cycle) determination from raw per-cycle amplification data.

Most cyclers export per-well fluorescence by cycle; the Cq itself is usually called
on the instrument. This module makes OligoForge able to call it too, so a raw run can
flow straight into the standard-curve / efficiency tools without leaving the app. Two
independent methods are provided:

  threshold     -- baseline-subtract, place a threshold at sd_mult x baseline-noise SD
                   (or an explicit value), and interpolate the fractional crossing cycle
                   on a log10(fluorescence) scale (the exponential phase is linear there).
  SDM           -- threshold-free second-derivative maximum (the cycle where d2F/dn2 peaks),
                   the method LightCycler-style instruments use; refined by parabolic
                   interpolation of the discrete second derivative.

Pure Python (no numpy/scipy) to match the rest of the local-first engine. Reference:
Bustin SA et al. 2009, Clin Chem 55:611-622 (MIQE).
"""
import math
import statistics


def _coerce(cycles, fluor):
    f = [float(v) for v in fluor]
    if cycles is None:
        c = [float(i + 1) for i in range(len(f))]      # default 1..N
    else:
        c = [float(v) for v in cycles]
    if len(c) != len(f):
        raise ValueError("cycles and fluor must have equal length")
    if len(f) < 6:
        raise ValueError("need >= 6 cycles of fluorescence")
    return c, f


def _baseline_subtract(c, f, baseline):
    lo, hi = baseline
    idx = [i for i, x in enumerate(c) if lo <= x <= hi]
    if len(idx) >= 2:
        base = statistics.mean(f[i] for i in idx)
    else:
        k = max(1, len(f) // 5)
        base = statistics.mean(f[:k])
        idx = list(range(k))
    return [v - base for v in f], float(base), idx


def cq_threshold(cycles, fluor, threshold=None, baseline=(3, 15), sd_mult=10.0):
    """(Cq, threshold_used). Cq is NaN if the curve never crosses the threshold."""
    c, f = _coerce(cycles, fluor)
    fb, _base, bidx = _baseline_subtract(c, f, baseline)
    base_vals = [fb[i] for i in bidx]
    if len(base_vals) > 2:
        noise = statistics.stdev(base_vals)
    else:
        pos = [v for v in fb if v > 0]
        noise = statistics.stdev(pos) if len(pos) > 2 else 0.0
    thr = float(threshold) if threshold is not None else max(sd_mult * noise, 1e-12)

    cross = next((k for k in range(len(fb)) if fb[k] >= thr), None)
    if cross is None or cross == 0:
        return float("nan"), thr
    f1, f2 = fb[cross - 1], fb[cross]
    c1, c2 = c[cross - 1], c[cross]
    if f1 <= 0:                                  # linear bridge across the baseline
        cq = c1 + (thr - f1) / (f2 - f1) * (c2 - c1)
    else:                                        # log-linear within the exponential phase
        cq = c1 + (math.log10(thr) - math.log10(f1)) / (math.log10(f2) - math.log10(f1)) * (c2 - c1)
    return float(cq), thr


def _smooth3(y, passes=1):
    for _ in range(passes):
        out = y[:]
        for i in range(1, len(y) - 1):
            out[i] = (y[i - 1] + y[i] + y[i + 1]) / 3.0
        y = out
    return y


def cq_second_derivative(cycles, fluor, baseline=(3, 15), smooth=True):
    """Threshold-free Cq = cycle at the maximum of the second derivative (SDM),
    refined by parabolic interpolation. Light 3-point smoothing optional. NaN if
    the curve is flat (no amplification)."""
    c, f = _coerce(cycles, fluor)
    fb, _base, _bidx = _baseline_subtract(c, f, baseline)
    y = _smooth3(fb, passes=2) if (smooth and len(fb) >= 7) else fb
    if max(y) - min(y) <= 0:
        return float("nan")
    d1 = [0.0] * len(y)
    for i in range(1, len(y) - 1):
        d1[i] = (y[i + 1] - y[i - 1]) / 2.0
    d2 = [0.0] * len(y)
    for i in range(1, len(y) - 1):
        d2[i] = y[i + 1] - 2.0 * y[i] + y[i - 1]
    k = max(range(len(d2)), key=lambda i: d2[i])
    if 1 <= k < len(d2) - 1:
        y0, y1, y2 = d2[k - 1], d2[k], d2[k + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        # keep the refinement within the neighbouring cycle spacing
        if -1.0 <= delta <= 1.0:
            return float(c[k] + delta * (c[k] - c[k - 1]))
    return float(c[k])


def analyze(fluor, cycles=None, threshold=None, baseline=(3, 15), sd_mult=10.0):
    """Both Cq methods plus diagnostics, for one amplification trace.

    Returns threshold/SDM Cq, the threshold used, baseline window, and a simple
    amplified flag (final signal rises clearly above baseline noise)."""
    try:
        c, f = _coerce(cycles, fluor)
    except ValueError as e:
        return dict(error=str(e))
    cq_t, thr = cq_threshold(c, f, threshold=threshold, baseline=baseline, sd_mult=sd_mult)
    cq_s = cq_second_derivative(c, f, baseline=baseline)
    fb, base, bidx = _baseline_subtract(c, f, baseline)
    base_vals = [fb[i] for i in bidx]
    noise = statistics.stdev(base_vals) if len(base_vals) > 2 else 0.0
    amplified = (max(fb) > 10.0 * noise) if noise > 0 else (max(fb) > 0)
    out = dict(n_cycles=len(c),
               cq_threshold=round(cq_t, 2) if cq_t == cq_t else None,
               cq_sdm=round(cq_s, 2) if cq_s == cq_s else None,
               threshold=round(thr, 4),
               baseline=[baseline[0], baseline[1]], baseline_fluor=round(base, 4),
               amplified=bool(amplified))
    if not out["amplified"]:
        out["note"] = "no clear amplification above baseline noise; Cq may be unreliable"
    elif out["cq_threshold"] is not None and out["cq_sdm"] is not None:
        d = abs(out["cq_threshold"] - out["cq_sdm"])
        if d > 1.5:
            out["note"] = ("threshold and SDM Cq differ by %.1f cycles; check baseline window "
                           "and curve shape" % d)
    return out
