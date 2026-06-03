"""Cq-calling regression. Offline. Run from repo root:  python tests/test_cq.py

Validates raw-fluorescence -> Cq (oligoforge.cq) on synthetic amplification curves with
PROPERTY-based asserts, not a hard-coded magic Cq: (1) a curve that inflects later (less
template) gives a LARGER Cq, for BOTH methods; (2) the threshold and SDM methods agree within
a few cycles on a clean curve; (3) Cq is noise-robust; (4) a flat trace is reported as
non-amplified. Deterministic checks use an explicit threshold so the crossing-cycle ordering
is tested directly; the auto (noise x sd_mult) threshold path is covered by the amplified flag.
"""
import sys, math, random
sys.path.insert(0, ".")
from oligoforge import cq as CQ

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        fails.append(name)


def curve(n0, k=0.6, amp=1000.0, base=50.0, n=45, noise=0.0, seed=7):
    """Logistic amplification trace inflecting at cycle n0, flat baseline before it."""
    random.seed(seed)
    out = []
    for i in range(1, n + 1):
        v = base + amp / (1.0 + math.exp(-k * (i - n0)))
        if noise:
            v += random.uniform(-noise, noise)
        out.append(v)
    return out


THR = 100.0   # explicit threshold (~10% of plateau, mid-exponential) for deterministic ordering
early = CQ.analyze(curve(24), threshold=THR)
late = CQ.analyze(curve(32), threshold=THR)
print("  early(n0=24):", {k: early[k] for k in ("cq_threshold", "cq_sdm", "amplified")})
print("  late (n0=32):", {k: late[k] for k in ("cq_threshold", "cq_sdm", "amplified")})

check("early curve amplified (auto-threshold path)", CQ.analyze(curve(24)).get("amplified") is True)
check("late curve amplified (auto-threshold path)", CQ.analyze(curve(32)).get("amplified") is True)
check("threshold Cq larger for later curve",
      late["cq_threshold"] > early["cq_threshold"] + 3, (early["cq_threshold"], late["cq_threshold"]))
check("SDM Cq larger for later curve",
      late["cq_sdm"] > early["cq_sdm"] + 3, (early["cq_sdm"], late["cq_sdm"]))
check("8-cycle inflection shift moves threshold Cq ~8 cycles",
      6.0 <= (late["cq_threshold"] - early["cq_threshold"]) <= 10.0,
      late["cq_threshold"] - early["cq_threshold"])
check("threshold and SDM agree within 3 cycles",
      abs(early["cq_threshold"] - early["cq_sdm"]) <= 3.0, (early["cq_threshold"], early["cq_sdm"]))

noisy = CQ.analyze(curve(26, noise=8.0), threshold=THR)
clean = CQ.analyze(curve(26), threshold=THR)
check("noisy threshold Cq within 1.5 cycles of clean",
      abs(noisy["cq_threshold"] - clean["cq_threshold"]) <= 1.5,
      (noisy["cq_threshold"], clean["cq_threshold"]))

flat = CQ.analyze([100.0 + (i % 2) * 0.5 for i in range(40)])
check("flat trace not amplified", flat.get("amplified") is False, flat)
check("flat trace Cq null", flat.get("cq_threshold") is None, flat)

expl = CQ.analyze(curve(25), cycles=list(range(1, 46)), threshold=200.0)
check("explicit threshold honored", abs(expl["threshold"] - 200.0) < 1e-6, expl["threshold"])
check("explicit-threshold Cq is finite", expl["cq_threshold"] is not None, expl)

check("short input returns error", "error" in CQ.analyze([1, 2, 3]))

if fails:
    print("CQ TEST FAILURES:", fails); sys.exit(1)
print("ALL CQ ASSERTS PASS")
