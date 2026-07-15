"""Cq determination from raw amplification traces.

The caller provides threshold-crossing and second-derivative estimates, but only
reports them when the trace passes a conservative amplification-shape gate.  The
baseline is linearly detrended, crossings must persist, and malformed/non-monotonic
cycle axes are rejected.  These are screening calculations; instrument-specific
algorithms and run-level QC remain authoritative for validated workflows.
"""
import math
import statistics


def _coerce(cycles, fluor):
    try:
        f = [float(v) for v in fluor]
        c = [float(i + 1) for i in range(len(f))] if cycles is None else [float(v) for v in cycles]
    except (TypeError, ValueError):
        raise ValueError("cycles and fluorescence must be numeric")
    if len(c) != len(f):
        raise ValueError("cycles and fluor must have equal length")
    if len(f) < 6:
        raise ValueError("need >= 6 cycles of fluorescence")
    if any(not math.isfinite(v) for v in c + f):
        raise ValueError("cycles and fluorescence must be finite")
    if any(c[i] <= c[i - 1] for i in range(1, len(c))):
        raise ValueError("cycle values must be strictly increasing")
    return c, f


def _baseline_model(c, f, baseline):
    lo, hi = float(baseline[0]), float(baseline[1])
    idx = [i for i, x in enumerate(c) if lo <= x <= hi]
    fallback = False
    if len(idx) < 3:
        k = min(len(f), max(3, len(f) // 5))
        idx = list(range(k))
        fallback = True
    xs = [c[i] for i in idx]
    ys = [f[i] for i in idx]
    xb = statistics.mean(xs); yb = statistics.mean(ys)
    sxx = sum((x - xb) ** 2 for x in xs)
    slope = (sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / sxx) if sxx > 0 else 0.0
    intercept = yb - slope * xb
    residual = [v - (intercept + slope * x) for x, v in zip(c, f)]
    return residual, dict(indices=idx, intercept=intercept, slope=slope,
                          mean=yb, fallback=fallback,
                          used=[c[idx[0]], c[idx[-1]]])


def _baseline_subtract(c, f, baseline):
    """Compatibility wrapper returning detrended values, baseline mean and indices."""
    fb, info = _baseline_model(c, f, baseline)
    return fb, float(info["mean"]), info["indices"]


def _noise(fb, idx):
    vals = [fb[i] for i in idx]
    return statistics.stdev(vals) if len(vals) > 2 else 0.0


def _sustained_crossing(fb, threshold, start, sustain=3):
    sustain = max(2, int(sustain))
    for k in range(max(1, start), len(fb) - sustain + 1):
        if all(fb[j] >= threshold for j in range(k, k + sustain)):
            return k
    return None


def _interpolate_crossing(c, fb, k, threshold):
    f1, f2 = fb[k - 1], fb[k]
    c1, c2 = c[k - 1], c[k]
    if f2 == f1:
        return c[k]
    if f1 > 0 and f2 > 0 and threshold > 0:
        den = math.log10(f2) - math.log10(f1)
        if den != 0:
            return c1 + (math.log10(threshold) - math.log10(f1)) / den * (c2 - c1)
    return c1 + (threshold - f1) / (f2 - f1) * (c2 - c1)


def cq_threshold(cycles, fluor, threshold=None, baseline=(3, 15), sd_mult=10.0,
                 sustain=3):
    """Return ``(Cq, threshold_used)``; Cq is NaN without a sustained crossing."""
    c, f = _coerce(cycles, fluor)
    fb, info = _baseline_model(c, f, baseline)
    noise = _noise(fb, info["indices"])
    try:
        thr = float(threshold) if threshold is not None else max(float(sd_mult) * noise, 1e-12)
    except (TypeError, ValueError):
        raise ValueError("threshold and sd_mult must be numeric")
    if not math.isfinite(thr) or thr <= 0:
        raise ValueError("threshold must be finite and greater than zero")
    cross = _sustained_crossing(fb, thr, max(info["indices"]) + 1, sustain=sustain)
    if cross is None:
        return float("nan"), thr
    return float(_interpolate_crossing(c, fb, cross, thr)), thr


def _smooth3(y, passes=1):
    y = list(y)
    for _ in range(max(0, int(passes))):
        out = y[:]
        for i in range(1, len(y) - 1):
            out[i] = (y[i - 1] + y[i] + y[i + 1]) / 3.0
        y = out
    return y


def cq_second_derivative(cycles, fluor, baseline=(3, 15), smooth=True):
    """Threshold-free second-derivative maximum with uneven-cycle support."""
    c, f = _coerce(cycles, fluor)
    fb, _info = _baseline_model(c, f, baseline)
    y = _smooth3(fb, passes=2) if (smooth and len(fb) >= 7) else fb
    if max(y) - min(y) <= 1e-15:
        return float("nan")
    d2 = [float("-inf")] * len(y)
    for i in range(1, len(y) - 1):
        dxl = c[i] - c[i - 1]
        dxr = c[i + 1] - c[i]
        left = (y[i] - y[i - 1]) / dxl
        right = (y[i + 1] - y[i]) / dxr
        d2[i] = 2.0 * (right - left) / (dxl + dxr)
    valid = range(1, len(y) - 1)
    k = max(valid, key=lambda i: d2[i])
    if not math.isfinite(d2[k]) or d2[k] <= 0:
        return float("nan")
    # Local parabola in index space; then scale by neighbouring cycle spacing.
    if 1 <= k < len(d2) - 1 and all(math.isfinite(d2[j]) for j in (k - 1, k, k + 1)):
        y0, y1, y2 = d2[k - 1], d2[k], d2[k + 1]
        denom = y0 - 2.0 * y1 + y2
        delta = 0.5 * (y0 - y2) / denom if denom else 0.0
        if -1.0 <= delta <= 1.0:
            step = (c[k + 1] - c[k - 1]) / 2.0
            return float(c[k] + delta * step)
    return float(c[k])


def _amplification_gate(c, fb, info, threshold):
    idx = info["indices"]
    noise = _noise(fb, idx)
    start = max(idx) + 1
    cross = _sustained_crossing(fb, threshold, start, sustain=3)
    tail_n = min(5, max(3, len(fb) // 10))
    tail = statistics.mean(fb[-tail_n:])
    amplitude = max(fb[start:] or fb) - min(fb[start:] or fb)
    post = fb[start:]
    rises = sum(1 for i in range(1, len(post)) if post[i] > post[i - 1])
    rise_fraction = rises / max(1, len(post) - 1)
    signal_floor = max(10.0 * noise, 3.0 * threshold, 1e-12)
    amplified = bool(cross is not None and amplitude >= signal_floor and
                     tail >= max(1.5 * threshold, 5.0 * noise) and rise_fraction >= 0.45)
    return amplified, dict(noise=noise, tail=tail, amplitude=amplitude,
                           rise_fraction=rise_fraction, crossing_index=cross)


def analyze(fluor, cycles=None, threshold=None, baseline=(3, 15), sd_mult=10.0):
    """Return two Cq estimates plus trace-shape diagnostics for one well."""
    try:
        c, f = _coerce(cycles, fluor)
        fb, info = _baseline_model(c, f, baseline)
        cq_t, thr = cq_threshold(c, f, threshold=threshold, baseline=baseline,
                                 sd_mult=sd_mult, sustain=3)
        cq_s = cq_second_derivative(c, f, baseline=baseline)
    except ValueError as e:
        return dict(error=str(e))

    amplified, gate = _amplification_gate(c, fb, info, thr)
    if not amplified:
        cq_t = cq_s = float("nan")

    out = dict(
        n_cycles=len(c),
        cq_threshold=round(cq_t, 2) if math.isfinite(cq_t) else None,
        cq_sdm=round(cq_s, 2) if math.isfinite(cq_s) else None,
        threshold=round(thr, 6),
        baseline=[baseline[0], baseline[1]],
        baseline_used=[round(info["used"][0], 3), round(info["used"][1], 3)],
        baseline_fluor=round(info["mean"], 6),
        baseline_slope=round(info["slope"], 8),
        baseline_fallback=bool(info["fallback"]),
        baseline_noise=round(gate["noise"], 6),
        amplitude=round(gate["amplitude"], 6),
        rise_fraction=round(gate["rise_fraction"], 3),
        amplified=amplified,
    )
    notes = []
    if info["fallback"]:
        notes.append("requested baseline contained fewer than 3 cycles; used the first fifth of the trace")
    if not amplified:
        notes.append("trace failed the sustained-crossing/plateau amplification gate; Cq values suppressed")
    elif out["cq_threshold"] is not None and out["cq_sdm"] is not None:
        delta = abs(out["cq_threshold"] - out["cq_sdm"])
        if delta > 1.5:
            notes.append("threshold and SDM Cq differ by %.1f cycles; inspect baseline and curve shape" % delta)
    if notes:
        out["note"] = " · ".join(notes)
    return out
