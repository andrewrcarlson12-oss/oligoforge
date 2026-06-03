"""Locked-panel integrity gate.  Offline.  Run from repo root:  python tests/test_locked_panel.py

WHY THIS EXISTS
The five FSJ assays that go on the real IDT order (Wilson/Burtt $3,000 grant) are seeded
verbatim by `seedFSJ` in static/index.html, and a second copy of the host primers lives in
the Python regression's LOCKED dict.  Nothing previously enforced that those two copies agree,
and the UI harnesses only checked that the assay NAMES render -- not that the bases are correct.
So a one-character edit to a sequence in one place (and not the other) could ship a wrong-base
order with the whole gate still green.  This test is the single authority: the canonical panel
lives HERE, and the test fails if static/index.html's seedFSJ OR the regression's LOCKED dict
drifts from it.  It also pins each oligo's Tm/GC against the current primer3 so a silent
numerical drift (a primer3 bump, an accidental COND edit) is caught, not just a sequence edit.

If you intentionally change an ordered sequence, update CANON below and re-confirm in IDT
OligoAnalyzer -- that is the deliberate, reviewed path.  An accidental drift fails loudly.
"""
import os, re, sys
sys.path.insert(0, ".")
from oligoforge import thermo as T

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INDEX = os.path.join(ROOT, "static", "index.html")
REGRESSION = os.path.join(HERE, "test_regression.py")

TM_TOL = 0.2   # C; tight enough to catch a real drift, loose enough to ignore float noise
GC_TOL = 0.1   # %

# ---- THE CANONICAL ORDERED PANEL (single source of truth) ----
# name: (forward, reverse, probe, amplicon_bp, chemistry_keyword)
CANON = {
    "IFNG":  ("AGTCATTCTGATGTCGCTGATG", "ACCTGTCAGTGTTTTCAAGCA",
              "TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT", 136, "ZEN"),
    "IL4":   ("AACTTGCTCAGCCTGGTTTG", "ATTCTTTAGTGAGGTGGTGCTG",
              "CTTGTGCCCTGCTCTGGTCCC", 84, "ZEN"),
    "RPL13": ("TCGCTGGCATCAACAAGAAG", "TCGGGAAGAGGATGAGCTTG",
              "AACAAGTCCACCGAGTCCCTGCA", 138, "ZEN"),
    "YWHAZ": ("CCGTTACTTGGCTGAGGTTG", "GATGGGATGTGTTGGTTGCA",
              "CCACTATCCCTTTCTTGTCATCTCCAGCAG", 121, "ZEN"),
    "Plasmodium": ("TACCTGGACTWGTTTCATGG", "AAAGGATTTGTGCTACCTTG",
                   "CTTACAAGATATCCACCACA", 157, "LNA"),
}

# Golden Tm/GC pinned against the current primer3 + qPCR COND.
# (name, role, gc, tm_or_None, tm_min_or_None, tm_max_or_None) -- degenerate oligos pin a range.
GOLDEN = [
    ("IFNG", "F", 45.5, 60.7, None, None), ("IFNG", "R", 42.9, 60.5, None, None),
    ("IFNG", "P", 45.5, 69.7, None, None),
    ("IL4", "F", 50.0, 61.5, None, None), ("IL4", "R", 45.5, 60.7, None, None),
    ("IL4", "P", 66.7, 68.5, None, None),
    ("RPL13", "F", 50.0, 61.2, None, None), ("RPL13", "R", 55.0, 62.0, None, None),
    ("RPL13", "P", 56.5, 68.8, None, None),
    ("YWHAZ", "F", 55.0, 61.5, None, None), ("YWHAZ", "R", 50.0, 61.4, None, None),
    ("YWHAZ", "P", 50.0, 68.5, None, None),
    ("Plasmodium", "F", 45.0, None, 56.4, 57.5),   # degenerate W -> Tm range
    ("Plasmodium", "R", 40.0, 56.0, None, None), ("Plasmodium", "P", 40.0, 54.4, None, None),
]

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)


def _oligos_in(text):
    """Every IUPAC oligo string (>=15 nt) literal in a chunk of source, as a multiset."""
    return sorted(re.findall(r'["\']([ACGTWRYSKMBDHVN]{15,40})["\']', text))


# 1) the canonical multiset of all 15 ordered oligos
canon_oligos = sorted(s for v in CANON.values() for s in v[:3])
check("canonical panel holds 15 oligos", len(canon_oligos) == 15, len(canon_oligos))

# 2) static/index.html seedFSJ matches the canonical panel EXACTLY (byte-for-byte, no add/drop/edit)
html = open(INDEX, encoding="utf-8").read()
m = re.search(r"function seedFSJ\s*\(\)\s*\{(.*?)\n\s*\}", html, re.S)
check("seedFSJ() found in index.html", bool(m))
seed_oligos = _oligos_in(m.group(1)) if m else []
check("seedFSJ oligos == canonical (no drift in the ordered sequences)",
      seed_oligos == canon_oligos, "seed=%s" % seed_oligos)
# amplicons present in the seed body
if m:
    body = m.group(1)
    for name, (_f, _r, _p, amp, _c) in CANON.items():
        check("seedFSJ %s amplicon %d present" % (name, amp), str(amp) in body, name)

# 3) the Python regression's keyed panel sequences agree with canonical, name+role for name+role.
# We match only keys named like the ordered panel (IFNG_F, IL4_R, Plas*_P, ...), so synthetic
# dimer/matrix fixtures and the HMBS validation anchor in the regression are correctly ignored.
reg = open(REGRESSION, encoding="utf-8").read()
ALIAS = {"PLAS": "Plasmodium", "PCYTB": "Plasmodium", "PLASMODIUM": "Plasmodium"}
keyed = re.findall(r'["\'](IFNG|IL4|RPL13|YWHAZ|PLAS\w*|PCYTB\w*|PLASMODIUM\w*)_([FRP])["\']\s*:\s*["\']([ACGTWRYSKMBDHVN]{15,40})["\']',
                   reg, re.I)
n_checked = 0
for raw_name, role, seq in keyed:
    base = raw_name.upper()
    name = ALIAS.get(base, ALIAS.get(base.rstrip("_FRP"), raw_name if raw_name in CANON else None))
    if name is None and raw_name in CANON:
        name = raw_name
    if name in CANON:
        want = CANON[name][{"F": 0, "R": 1, "P": 2}[role.upper()]]
        check("regression %s_%s == canonical" % (name, role.upper()), seq == want, "reg=%s want=%s" % (seq, want))
        n_checked += 1
check("regression defines at least the host primer set (>=8 keyed seqs)", n_checked >= 8, n_checked)

# 4) pin Tm/GC of every ordered oligo against the current primer3
ROLE = {"F": 0, "R": 1, "P": 2}
for name, role, gc, tm, tmlo, tmhi in GOLDEN:
    seq = CANON[name][ROLE[role]]
    g = T.gc_percent(seq)
    check("%s %s GC %.1f" % (name, role, gc), abs(g - gc) <= GC_TOL, "%.1f" % g)
    if tm is not None:
        t = T.tm(seq)
        check("%s %s Tm %.1f +-%.2f" % (name, role, tm, TM_TOL), abs(t - tm) <= TM_TOL, "%.2f" % t)
    else:
        tr = T.tm_range(seq)
        check("%s %s degenerate Tm range [%.1f,%.1f]" % (name, role, tmlo, tmhi),
              tr["degenerate"] and abs(tr["min"] - tmlo) <= TM_TOL and abs(tr["max"] - tmhi) <= TM_TOL,
              "%s" % tr)

print("")
if _fails:
    print("LOCKED-PANEL GATE FAILED:", ", ".join(_fails)); sys.exit(1)
print("ALL LOCKED-PANEL ASSERTS PASS")
