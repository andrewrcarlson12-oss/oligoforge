"""Standard-curve quantification: copy number from concentration, and dilution series."""
AVOGADRO = 6.02214076e23
DS_MW_PER_BP = 650.0          # g/mol per bp, dsDNA approximation
# Average molecular weight per base/bp by molecule type. Using the dsDNA constant
# for an ssDNA or RNA standard (e.g. an RNA-virus transcript) over-states the mass
# per molecule and throws the copy number off by ~2x.
MW_PER_UNIT = {"dsDNA": 650.0, "ssDNA": 330.0, "RNA": 340.0}


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
    return dict(slope=round(slope, 4), intercept=round(intercept, 3),
                efficiency=round(eff, 4), efficiency_pct=eff_pct, r2=round(r2, 5),
                amp_factor=round(10 ** (-1 / slope), 4) if slope else None,
                n_points=n, dynamic_range=[min(q for q, _ in det), max(q for q, _ in det)],
                lod_practical=min(full) if full else None,
                efficiency_ok=90.0 <= eff_pct <= 110.0, r2_ok=r2 >= 0.98,
                slope_ok=-3.58 <= slope <= -3.10, levels=levels,
                notes="MIQE-aligned acceptance: efficiency 90-110% (slope -3.58 to -3.10), R2>=0.98. "
                      "A formal LOD needs ~60 replicates near the limit (Forootan 2017); lod_practical "
                      "is the lowest standard detected in all replicates here.")
