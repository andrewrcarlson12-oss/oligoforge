"""Offline in-silico PCR (autodesign.epcr_offline) — deterministic, no network.
Run: OLIGOFORGE_EMAIL=you@x python3 tests/test_epcr_offline.py   (exit 0 = pass)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import autodesign as AD, thermo as T

F = "ACGTACGTACGTACGTACGT"          # 20-mer forward
R = "TTGGCCAATTGGCCAATTGG"          # 20-mer reverse
FILL = "GATCGATCGATCAGCTAGCTAGGCATCGATCGATCGTAGCTAGCATCGATCGATCG"   # 56 nt, probe lives here
PROBE = FILL[10:30]

# off1 carries a real amplicon: F ... fill ... revcomp(R)  -> convergent product
off1 = F + FILL + T.revcomp(R)
# off2 is unrelated -> no primer placement, no product
off2 = "AAAAACCCCCGGGGGTTTTT" * 6
# off3 has F + fill + revcomp(R) but the base read by R's 3' end is broken -> not extension-competent.
# R's 3' terminus pairs with the LEFT edge of the revcomp(R) region (revcomp(R)[0]).
_rcR = list(T.revcomp(R)); _rcR[0] = {"A": "C", "C": "A", "G": "T", "T": "G"}[_rcR[0]]
off3 = F + FILL + "".join(_rcR)

fails = []

def check(name, cond):
    if not cond:
        fails.append(name)
    print(("PASS " if cond else "FAIL ") + name)

p1 = AD.epcr_offline(F, R, [off1], probe=PROBE)
check("off1 predicts exactly one product", len(p1) == 1)
check("off1 product size == full amplicon", p1 and p1[0]["size"] == len(off1))
check("off1 probe binds inside product", p1 and p1[0]["probe_binds"] is True)
check("off1 product references subject 0", p1 and p1[0]["subject"] == 0)

p2 = AD.epcr_offline(F, R, [off2], probe=PROBE)
check("off2 (unrelated) predicts no product", len(p2) == 0)

# off3: reverse primer 3' terminal base mismatches -> require_3prime drops the product
p3 = AD.epcr_offline(F, R, [off3], probe=PROBE)
check("off3 (broken 3' end) predicts no product", len(p3) == 0)

# mixed set: only off1 should yield a product, and it should point at the right subject index
pm = AD.epcr_offline(F, R, [off2, off1, off2], probe=PROBE)
check("mixed set finds the single real product", len(pm) == 1 and pm[0]["subject"] == 1)

# probe=None -> probe_binds is None (no false claim)
pn = AD.epcr_offline(F, R, [off1], probe=None)
check("probe=None -> probe_binds None", pn and pn[0]["probe_binds"] is None)

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL EPCR_OFFLINE ASSERTS PASS")
