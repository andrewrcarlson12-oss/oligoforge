"""Probe-site accessibility endpoint (/api/accessibility) — folds the amplicon and grades how
base-paired each binding site is at the annealing temperature. Skips cleanly when ViennaRNA is
not installed. Run: OLIGOFORGE_EMAIL=you@x python3 tests/test_accessibility.py  (exit 0 = pass)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import structure as STR, thermo as T

if not STR.available():
    print("ViennaRNA not installed -> accessibility is optional; skipping (OK)")
    sys.exit(0)

from fastapi.testclient import TestClient
import app as A
c = TestClient(A.app)
fails = []

def check(name, cond):
    if not cond:
        fails.append(name)
    print(("PASS " if cond else "FAIL ") + name)

def access(amp, F, R, probe=None):
    return c.post("/api/accessibility", json={"amplicon": amp, "forward": F, "reverse": R,
                                              "probe": probe, "anneal_c": 60}).json()

# 1) strong GC hairpin spanning the probe site -> inaccessible (risk)
F = "AAACAAACAAACAAACAA"; rcR = "TTTGTTTGTTTGTTTGTT"; R = T.revcomp(rcR)
hp = "GCGCGCGCGC" + "TTTTT" + "GCGCGCGCGC"
d1 = access(F + hp + rcR, F, R, hp)
check("hairpin folds", d1.get("folded") is True)
check("hairpin probe is paired at anneal temp", d1.get("probe", {}).get("anneal_paired", 0) >= 0.5)
check("hairpin verdict == risk", d1.get("verdict") == "risk")
check("hairpin uses DNA params", d1.get("dna_params") is True)

# 2) purine-only amplicon -> no Watson-Crick pairs possible -> accessible (ok)
amp = "AGAGAGAAGGAAGAGGAAGAGAGAAGGAAGAGGAAGAGAGAAGGAAGAGG"
F2 = amp[:18]; R2 = T.revcomp(amp[-18:]); probe2 = amp[18:38]
d2 = access(amp, F2, R2, probe2)
check("unfoldable probe is unpaired", d2.get("probe", {}).get("anneal_paired") == 0.0)
check("unfoldable verdict == ok", d2.get("verdict") == "ok")

# 3) SYBR (no probe) -> verdict comes from the primer 3' ends, no probe block
d3 = access(F + hp + rcR, F, R, None)
check("SYBR has no probe block", "probe" not in d3)
check("SYBR reports primer 3' sites", bool(d3.get("f3")) and bool(d3.get("r3")))

# 4) too-short amplicon -> guarded, not folded
d4 = c.post("/api/accessibility", json={"amplicon": "ACGT"}).json()
check("short amplicon guarded", d4.get("available") is True and d4.get("folded") is False)

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL ACCESSIBILITY ASSERTS PASS")
