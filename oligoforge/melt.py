"""Melt-curve (dissociation) analysis for SYBR qPCR.

MIQE requires a melt curve for SYBR assays to confirm a single specific product. As dsDNA
denatures with rising temperature the intercalating dye is released and fluorescence falls; the
melt temperature is the peak of -dF/dT. One sharp peak means one product. Extra peaks usually mean
primer-dimers (which melt lower, ~70-80 C) or off-target amplicons. This takes raw fluorescence vs
temperature and returns the melt peak temperature(s), so an observed run can be checked against the
amplicon melt Tm OligoForge predicts at design time (amplicon_tm). Pure Python (no numpy/scipy).
"""
import statistics


def _coerce(temps, fluor):
    f = [float(v) for v in fluor]
    t = [float(v) for v in temps]
    if len(t) != len(f):
        raise ValueError("temps and fluor must have equal length")
    if len(f) < 8:
        raise ValueError("need >= 8 temperature points")
    if t[0] > t[-1]:                       # normalize to ascending temperature
        t = t[::-1]; f = f[::-1]
    return t, f


def _smooth(y, passes=2):
    for _ in range(passes):
        out = y[:]
        for i in range(1, len(y) - 1):
            out[i] = (y[i - 1] + y[i] + y[i + 1]) / 3.0
        y = out
    return y


def melt_peaks(temps, fluor, min_prominence_frac=0.10, smooth=True):
    """Melt peaks as maxima of -dF/dT. Returns (peaks, temps, neg_deriv) where peaks is a list of
    (Tm, height) sorted by height, keeping only peaks at/above min_prominence_frac x the tallest
    (filters noise wiggles) and merging peaks within ~0.8 C."""
    t, f = _coerce(temps, fluor)
    y = _smooth(f) if (smooth and len(f) >= 7) else f
    d = [0.0] * len(t)
    for i in range(1, len(t) - 1):
        dt = t[i + 1] - t[i - 1]
        d[i] = -(y[i + 1] - y[i - 1]) / dt if dt != 0 else 0.0
    dmax = max(d) if d else 0.0
    if dmax <= 0:
        return [], t, d
    thr = min_prominence_frac * dmax
    peaks = []
    for i in range(2, len(d) - 2):
        if d[i] >= thr and d[i] >= d[i - 1] and d[i] >= d[i + 1] and d[i] > d[i - 2] and d[i] > d[i + 2]:
            y0, y1, y2 = d[i - 1], d[i], d[i + 1]
            denom = (y0 - 2 * y1 + y2)
            delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
            tm = t[i] + (delta * (t[i + 1] - t[i - 1]) / 2.0 if -1.0 <= delta <= 1.0 else 0.0)
            peaks.append((round(tm, 2), round(d[i], 4)))
    peaks.sort(key=lambda p: -p[1])
    merged = []
    for tm, h in peaks:
        if all(abs(tm - m[0]) > 0.8 for m in merged):
            merged.append((tm, h))
    return merged, t, d


def analyze(fluor, temps, min_prominence_frac=0.10, predicted_tm=None):
    """Melt-curve summary for one trace: peak Tm(s), single-product flag, and (optional)
    comparison to the amplicon Tm predicted at design time."""
    try:
        peaks, t, _d = melt_peaks(temps, fluor, min_prominence_frac=min_prominence_frac)
    except ValueError as e:
        return dict(error=str(e))
    out = dict(n_peaks=len(peaks),
               peaks=[dict(tm=tm, height=h) for tm, h in peaks],
               dominant_tm=peaks[0][0] if peaks else None,
               single_product=(len(peaks) == 1),
               temp_range=[round(min(t), 1), round(max(t), 1)])
    if not peaks:
        out["note"] = ("no clear melt peak; check this is a dissociation curve "
                       "(fluorescence falling as temperature rises)")
    elif len(peaks) > 1:
        out["note"] = ("%d peaks: a primer-dimer or off-target product is melting alongside the "
                       "amplicon -- SYBR signal is not specific" % len(peaks))
    if predicted_tm is not None and peaks:
        out["predicted_tm"] = round(float(predicted_tm), 1)
        out["dominant_minus_predicted"] = round(peaks[0][0] - float(predicted_tm), 1)
    return out
