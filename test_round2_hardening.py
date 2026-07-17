"""Regression gates for the v1.31.1 audit repairs.

Standalone by project convention: run with ``python tests/test_round2_hardening.py``.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoforge import specificity as SP
from oligoforge import multiplex as MX
from oligoforge import thermo as T
from oligoforge import autodesign as AD

fails = []

def check(name, condition, detail=""):
    print(("PASS" if condition else "FAIL"), "-", name,
          "" if condition else "[%s]" % detail)
    if not condition:
        fails.append(name)

# Ambiguous bases in a subject are possible matches, not definite mismatches.
F = "ACGTTGCACTGA"
R = "TCCGATGACCTA"
RC_R = SP._rc_iupac(R)
mid = "GATCACGATTCGATGCTAGCATCGATCGTACGATCA"
reference = F + mid + RC_R
ambiguous = F[:-1] + "N" + mid + RC_R
mismatch = F[:-1] + ("C" if F[-1] != "C" else "A") + mid + RC_R
ra = SP.in_silico_pcr_offline(F, R, ">ambiguous\n%s\n" % ambiguous, max_mm=0)
rm = SP.in_silico_pcr_offline(F, R, ">mismatch\n%s\n" % mismatch, max_mm=0)
check("ambiguous 3-prime subject base is retained as possible off-target",
      ra.get("n_products") == 1, ra)
check("ambiguous anchor is reported as uncertain",
      ra.get("uncertain_3prime_hits", 0) >= 1 and ra.get("uncertainty_note"), ra)
check("definite 3-prime mismatch blocks product",
      rm.get("n_products") == 0, rm)

bad_fasta = SP.in_silico_pcr_offline(F, R, ">bad\n%sX%s\n" % (F, RC_R))
check("non-IUPAC FASTA symbols are rejected", "invalid nucleotide" in bad_fasta.get("error", ""), bad_fasta)

# Distinct assays with the same display name must still be cross-screened.
a = "ACGATCAGTTGCATCAGGTA"
b = T.revcomp(a)
mx = MX.check([
    {"id": "a1", "name": "duplicate", "dye": "fam", "oligos": [{"name": "F", "seq": a}]},
    {"id": "a2", "name": "duplicate", "dye": "FAM", "oligos": [{"name": "R", "seq": b}]},
])
check("same-name assays are treated as distinct", mx.get("n_flagged", 0) >= 1, mx)
check("dye comparison is case-insensitive", len(mx.get("channel_conflicts", [])) == 1, mx)
row = (mx.get("cross_dimers") or [{}])[0]
check("multiplex dimer reports annealing context",
      all(k in row for k in ("dg", "dg_anneal", "dimer_tm", "end_dg", "three_prime")), row)

mx_bad = MX.check([
    {"name": "bad", "oligos": [{"name": "F", "seq": "ACGTACGTXYZ!!"}]},
    {"name": "ok", "oligos": [{"name": "R", "seq": b}]},
])
check("malformed multiplex oligo is rejected rather than sanitized",
      bool(mx_bad.get("error")) and mx_bad.get("n_invalid") == 1, mx_bad)

# The current collector lives in candidate_search and must inspect the full target before
# retention, even when only one finalist was requested. Pin the target-spanning behavior
# without monkeypatching the retired first-hit _design_one helper.
orig = AD.CSEARCH.D.generate_assay_candidates
calls = []
def fake_generate(window, profile, **kwargs):
    calls.append(window[0])
    base = window[0]
    row = {"forward": base * 12, "reverse": ("A" if base != "A" else "C") * 12,
           "probe": None, "amplicon": 80, "amplicon_tm": 80.0,
           "f_tm": 60.0, "r_tm": 60.0, "probe_info": None,
           "candidate_rank": 0.0, "f_xy": [0, 12], "r_xy": [68, 80],
           "amplicon_xy": [0, 80]}
    return [row], {"pairs_after_hard_gates": 1, "pairs_truncated": 0,
                   "triplets_truncated": 0}
try:
    AD.CSEARCH.D.generate_assay_candidates = fake_generate
    fake_ref = "A" * 350 + "C" * 350 + "G" * 350
    pool = AD._candidates(fake_ref, {"amp_max": 150}, n=1, window=350, step=350, budget_s=5)
finally:
    AD.CSEARCH.D.generate_assay_candidates = orig
starts = {a.get("search_window_start") for a in pool}
check("candidate collection does not stop after the first success", len(pool) > 1, starts)
check("candidate collection reaches the downstream target region", max(starts or {0}) >= 700, starts)

if fails:
    print("ROUND-2 HARDENING FAILURES:", ", ".join(fails))
    sys.exit(1)
print("ALL ROUND-2 HARDENING ASSERTS PASS")
