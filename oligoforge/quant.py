"""Standard-curve quantification: copy number from concentration, and dilution series."""
AVOGADRO = 6.02214076e23
DS_MW_PER_BP = 650.0          # g/mol per bp, dsDNA approximation
# Average molecular weight per base/bp by molecule type. Using the dsDNA constant
# for an ssDNA or RNA standard (e.g. an RNA-virus transcript) over-states the mass
# per molecule and throws the copy number off by ~2x.
MW_PER_UNIT = {"dsDNA": 650.0, "ssDNA": 330.0, "RNA": 340.0}

# Two-sided t critical values t_{0.975, df} for the standard-curve slope CI (no SciPy dep).
_T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
         8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
         15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
         22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056, 27: 2.052, 28: 2.048,
         29: 2.045, 30: 2.042}


def _mw(molecule_type):
    return MW_PER_UNIT.get(molecule_type or "dsDNA", DS_MW_PER_BP)


def copies_per_ul(ng_per_ul, length_bp, molecule_type="dsDNA"):
    grams = ng_per_ul * 1e-9
    moles = grams / (length_bp * _mw(molecule_type))
    return moles * AVOGADRO


def ng_for_copies(copies_per_ul, length_bp, molecule_type="dsDNA"):
    moles = copies_per_ul / AVOGADRO
    grams = moles * length_bp * _mw(molecule_type)
    return grams / 1e-9


def dilution_series(start_copies_per_ul, factor=10, points=6):
    return [start_copies_per_ul / (factor ** i) for i in range(points)]


def _lod95_logistic(points, target=0.95, iters=60):
    """MIQE-style limit of detection: the input quantity at which detection
    probability reaches `target` (0.95), from a logistic fit of detect/non-detect
    against log10(quantity). Returns the quantity or None when the data can't
    support it (no detect/non-detect transition, too few points, non-convergence,
    or an estimate that would be pure extrapolation beyond the tested range).
    Pure-Python Newton-Raphson on standardized x for numerical stability."""
    import math, statistics
    data = [(math.log10(q), 1.0 if c is not None else 0.0) for q, c in points if q and q > 0]
    if len(data) < 4:
        return None
    ys = [y for _, y in data]
    if all(y == 1.0 for y in ys) or all(y == 0.0 for y in ys):
        return None                      # no transition: LOD is unidentifiable from these data
    xs = [x for x, _ in data]
    mx = statistics.mean(xs); sx = statistics.pstdev(xs) or 1.0
    xz = [(x - mx) / sx for x in xs]
    b0 = b1 = 0.0
    for _ in range(iters):
        g0 = g1 = h00 = h01 = h11 = 0.0
        for x, y in zip(xz, ys):
            z = b0 + b1 * x
            p = 1.0 / (1.0 + math.exp(-z)) if -700.0 < z < 700.0 else (1.0 if z > 0 else 0.0)
            w = p * (1.0 - p)
            g0 += (y - p); g1 += (y - p) * x
            h00 += w; h01 += w * x; h11 += w * x * x
        det = h00 * h11 - h01 * h01
        if abs(det) < 1e-12:
            break
        db0 = (h11 * g0 - h01 * g1) / det
        db1 = (-h01 * g0 + h00 * g1) / det
        b0 += db0; b1 += db1
        if abs(db0) + abs(db1) < 1e-9:
            break
    if b1 <= 0:                          # detection must increase with quantity
        return None
    zt = math.log(target / (1.0 - target))
    x95 = ((zt - b0) / b1) * sx + mx      # back to log10(quantity)
    if not (min(xs) - 1.0 <= x95 <= max(xs) + 1.0):
        return None                       # would be >1 log of extrapolation; don't claim it
    return 10.0 ** x95


def standard_curve(points):
    """Fit Cq = slope*log10(quantity) + intercept over a standard dilution series.

    points: list of (quantity, Cq). Cq may be None for a non-detect. Repeat a
    quantity for replicates. Reports efficiency E = 10^(-1/slope) - 1, R2, dynamic
    range, per-level replicate stats, and a practical LOD (lowest standard detected
    in every replicate). A formal MIQE LOD needs ~60 replicates near the limit."""
    import math, statistics
    det = [(q, c) for q, c in points if c is not None and q and q > 0]
    if len(det) < 2:
        return dict(error="need >=2 detected points with positive quantity")
    xs = [math.log10(q) for q, _ in det]
    ys = [c for _, c in det]
    n = len(xs)
    xbar = sum(xs) / n; ybar = sum(ys) / n
    Sxx = sum((x - xbar) ** 2 for x in xs)
    Syy = sum((y - ybar) ** 2 for y in ys)
    Sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    if Sxx == 0:
        return dict(error="all standards at one quantity; need a range")
    slope = Sxy / Sxx
    intercept = ybar - slope * xbar
    r2 = (Sxy * Sxy) / (Sxx * Syy) if Syy else 0.0
    eff = 10 ** (-1 / slope) - 1 if slope else 0.0

    # Slope standard error, 95% CI, and propagated efficiency CI (MIQE reporting).
    # SS_resid = Syy - slope*Sxy; MSE = SS_resid/(n-2); SE_slope = sqrt(MSE/Sxx).
    # Efficiency is monotonic increasing in slope, so the slope-CI endpoints map
    # directly to the efficiency CI. Needs n>=3 (>=1 residual degree of freedom).
    slope_se = None; slope_ci = None; eff_ci_pct = None
    if n >= 3 and Sxx > 0:
        ss_resid = max(Syy - slope * Sxy, 0.0)
        mse = ss_resid / (n - 2)
        slope_se = (mse / Sxx) ** 0.5
        tcrit = _T975.get(n - 2, 1.96)
        s_lo = slope - tcrit * slope_se      # more negative slope bound
        s_hi = slope + tcrit * slope_se      # less negative slope bound
        slope_ci = [round(s_lo, 4), round(s_hi, 4)]
        if s_lo < 0 and s_hi < 0:
            e_lo = 10 ** (-1 / s_lo) - 1
            e_hi = 10 ** (-1 / s_hi) - 1
            eff_ci_pct = [round(e_lo * 100, 1), round(e_hi * 100, 1)]

    groups = {}
    for q, c in points:
        groups.setdefault(q, []).append(c)
    levels = []
    for q in sorted(groups, reverse=True):
        cs = [c for c in groups[q] if c is not None]
        nrep = len(groups[q]); ndet = len(cs)
        levels.append(dict(quantity=q, n=nrep, detected=ndet,
                           mean_cq=round(statistics.mean(cs), 2) if cs else None,
                           sd_cq=round(statistics.pstdev(cs), 3) if len(cs) > 1 else None,
                           detection_rate=round(100 * ndet / nrep, 0) if nrep else 0))
    full = [l["quantity"] for l in levels if l["n"] and l["detected"] == l["n"]]
    eff_pct = round(eff * 100, 1)
    lod95 = _lod95_logistic(points)
    return dict(slope=round(slope, 4), intercept=round(intercept, 3),
                efficiency=round(eff, 4), efficiency_pct=eff_pct, r2=round(r2, 5),
                slope_se=round(slope_se, 4) if slope_se is not None else None,
                slope_ci95=slope_ci, efficiency_ci_pct=eff_ci_pct,
                amp_factor=round(10 ** (-1 / slope), 4) if slope else None,
                n_points=n, dynamic_range=[min(q for q, _ in det), max(q for q, _ in det)],
                lod_practical=min(full) if full else None,
                lod95=round(lod95, 4) if lod95 is not None else None,
                efficiency_ok=90.0 <= eff_pct <= 110.0, r2_ok=r2 >= 0.98,
                slope_ok=-3.58 <= slope <= -3.10, levels=levels,
                notes="MIQE-aligned acceptance: efficiency 90-110% (slope -3.58 to -3.10), R2>=0.98. "
                      "slope_ci95 / efficiency_ci_pct are 95% CIs (t, n-2 df); they need n>=3 points. "
                      "lod95 is the 95%-detection LOD from a logistic fit of detect/non-detect vs "
                      "log10(quantity) (needs a detection transition across levels; null otherwise). "
                      "lod_practical is the lowest standard detected in all replicates (Forootan 2017).")
