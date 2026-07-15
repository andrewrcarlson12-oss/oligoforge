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
    """Fit Cq against log10(quantity) using one mean Cq per dilution level.

    Technical replicates characterize precision and detection probability; they are
    not independent x-axis levels.  Regressing every replicate separately allows an
    unbalanced replicate count to overweight one concentration and produces an
    artificially large residual degree of freedom.  This implementation therefore
    aggregates detected Cq values by quantity for slope/R2/CI, while retaining all
    replicate calls for level statistics and the exploratory detection model.
    """
    import math, statistics

    cleaned = []
    for q, c in points:
        try:
            qv = float(q)
            cv = None if c is None else float(c)
        except (TypeError, ValueError):
            return dict(error="quantity and Cq values must be numeric (or Cq blank for non-detect)")
        if not math.isfinite(qv) or qv <= 0:
            return dict(error="every standard quantity must be finite and greater than zero")
        if cv is not None and (not math.isfinite(cv) or cv <= 0 or cv > 80):
            return dict(error="detected Cq values must be finite and in (0, 80]")
        cleaned.append((qv, cv))

    groups = {}
    for q, c in cleaned:
        groups.setdefault(q, []).append(c)

    levels = []
    regression = []
    for q in sorted(groups, reverse=True):
        all_calls = groups[q]
        detected = [c for c in all_calls if c is not None]
        nrep, ndet = len(all_calls), len(detected)
        mean_cq = statistics.mean(detected) if detected else None
        levels.append(dict(
            quantity=q, n=nrep, detected=ndet,
            mean_cq=round(mean_cq, 2) if mean_cq is not None else None,
            sd_cq=round(statistics.stdev(detected), 3) if len(detected) > 1 else None,
            detection_rate=round(100.0 * ndet / nrep, 1) if nrep else 0.0,
        ))
        if mean_cq is not None:
            regression.append((math.log10(q), mean_cq))

    if len(regression) < 2:
        return dict(error="need detected amplification at two or more distinct standard quantities",
                    levels=levels)

    xs = [x for x, _ in regression]
    ys = [y for _, y in regression]
    n_levels = len(xs)
    xbar = statistics.mean(xs); ybar = statistics.mean(ys)
    Sxx = sum((x - xbar) ** 2 for x in xs)
    Syy = sum((y - ybar) ** 2 for y in ys)
    Sxy = sum((x - xbar) * (y - ybar) for x, y in regression)
    if Sxx <= 0:
        return dict(error="all detected standards are at one quantity; need a dilution range",
                    levels=levels)

    slope = Sxy / Sxx
    intercept = ybar - slope * xbar
    ss_resid = sum((y - (slope * x + intercept)) ** 2 for x, y in regression)
    r2 = 1.0 - ss_resid / Syy if Syy > 0 else (1.0 if ss_resid <= 1e-15 else 0.0)
    r2 = min(1.0, max(0.0, r2))

    def _efficiency(s):
        if not math.isfinite(s) or s >= -0.25:
            return None
        exponent = -1.0 / s
        if exponent > 4.0:                 # prevents overflow / physically meaningless output
            return None
        return 10.0 ** exponent - 1.0

    eff = _efficiency(slope)
    slope_se = None; slope_ci = None; eff_ci_pct = None
    if n_levels >= 3:
        mse = ss_resid / (n_levels - 2)
        slope_se = math.sqrt(max(mse / Sxx, 0.0))
        tcrit = _T975.get(n_levels - 2, 1.96)
        s_lo = slope - tcrit * slope_se
        s_hi = slope + tcrit * slope_se
        slope_ci = [round(s_lo, 4), round(s_hi, 4)]
        e_lo, e_hi = _efficiency(s_lo), _efficiency(s_hi)
        if e_lo is not None and e_hi is not None:
            eff_ci_pct = [round(min(e_lo, e_hi) * 100.0, 1),
                          round(max(e_lo, e_hi) * 100.0, 1)]

    fully_detected = [l["quantity"] for l in levels if l["n"] and l["detected"] == l["n"]]
    lod_screen = _lod95_logistic(cleaned)
    eff_pct = round(eff * 100.0, 1) if eff is not None else None
    amp_factor = round(eff + 1.0, 4) if eff is not None else None
    detected_observations = sum(l["detected"] for l in levels)

    return dict(
        slope=round(slope, 4), intercept=round(intercept, 3),
        efficiency=round(eff, 4) if eff is not None else None,
        efficiency_pct=eff_pct, r2=round(r2, 5),
        slope_se=round(slope_se, 4) if slope_se is not None else None,
        slope_ci95=slope_ci, efficiency_ci_pct=eff_ci_pct,
        amp_factor=amp_factor,
        n_points=n_levels, n_levels=n_levels,
        n_observations=len(cleaned), n_detected_observations=detected_observations,
        dynamic_range=[min(10.0 ** x for x in xs), max(10.0 ** x for x in xs)],
        lod_practical=min(fully_detected) if fully_detected else None,
        lowest_fully_detected_standard=min(fully_detected) if fully_detected else None,
        lod95=round(lod_screen, 4) if lod_screen is not None else None,
        lod95_status=("exploratory logistic estimate" if lod_screen is not None else "not estimable"),
        efficiency_ok=(eff_pct is not None and 90.0 <= eff_pct <= 110.0),
        r2_ok=r2 >= 0.98,
        slope_ok=-3.58 <= slope <= -3.10,
        levels=levels,
        notes=("Regression uses the mean detected Cq at each distinct quantity; technical replicates "
               "contribute to SD and detection rate, not independent slope degrees of freedom. "
               "The lowest fully detected standard is a screening descriptor, not a validated LOD. "
               "lod95 is an exploratory logistic estimate and requires dedicated replicate-rich "
               "experiments near the detection boundary before publication or clinical claims."),
    )
