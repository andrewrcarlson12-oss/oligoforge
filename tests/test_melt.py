"""Melt-curve regression. Offline. Run from repo root:  python tests/test_melt.py

Validates SYBR melt analysis (oligoforge.melt): one product -> one peak at its Tm; a dimer +
amplicon -> two peaks at the right temperatures; a flat trace -> no peak. Property-based.
"""
import sys, math
sys.path.insert(0, ".")
from oligoforge import melt as M

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        fails.append(name)


def melt_curve(tms_heights, k=0.8, t0=65.0, t1=96.0, step=0.5):
    """Falling fluorescence vs temperature = sum of decreasing sigmoids, one per product."""
    temps, fl = [], []
    T = t0
    while T <= t1 + 1e-9:
        v = 0.0
        for tm, amp in tms_heights:
            v += amp / (1.0 + math.exp(k * (T - tm)))
        temps.append(round(T, 2)); fl.append(v)
        T += step
    return temps, fl


# one product at 85 C
t1, f1 = melt_curve([(85.0, 1000.0)])
r1 = M.analyze(f1, t1)
print("  single:", {k: r1[k] for k in ("n_peaks", "dominant_tm", "single_product")})
check("single product -> one peak", r1["n_peaks"] == 1, r1["n_peaks"])
check("single product Tm near 85 C", abs(r1["dominant_tm"] - 85.0) <= 1.0, r1["dominant_tm"])
check("single_product flag true", r1["single_product"] is True, r1)

# dimer at 76 C + amplicon at 85 C
t2, f2 = melt_curve([(76.0, 450.0), (85.0, 650.0)])
r2 = M.analyze(f2, t2)
print("  double:", {k: r2[k] for k in ("n_peaks", "peaks", "single_product")})
check("two products -> two peaks", r2["n_peaks"] == 2, r2["n_peaks"])
check("two-peak not flagged single", r2["single_product"] is False, r2)
_tms = sorted(p["tm"] for p in r2["peaks"])
check("peaks at ~76 and ~85", abs(_tms[0] - 76.0) <= 1.5 and abs(_tms[1] - 85.0) <= 1.5, _tms)

# predicted-Tm comparison
r3 = M.analyze(f1, t1, predicted_tm=84.0)
check("predicted-Tm delta computed", abs(r3["dominant_minus_predicted"] - (r3["dominant_tm"] - 84.0)) < 1e-6, r3)

# flat trace -> no peak
flat = M.analyze([500.0] * 40, [65.0 + 0.5 * i for i in range(40)])
check("flat trace -> no peak", flat["n_peaks"] == 0, flat)

# descending input is normalized (temps high->low still works)
r4 = M.analyze(f1[::-1], t1[::-1])
check("descending temps normalized", r4["n_peaks"] == 1 and abs(r4["dominant_tm"] - 85.0) <= 1.0, r4)

check("short input returns error", "error" in M.analyze([1, 2, 3], [80, 81, 82]))

if fails:
    print("MELT TEST FAILURES:", fails); sys.exit(1)
print("ALL MELT ASSERTS PASS")
