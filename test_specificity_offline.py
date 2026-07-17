"""
Validation of the OFFLINE in-silico-PCR specificity engine (specificity.in_silico_pcr_offline).

Fully deterministic and offline: RefSeq sequences are committed in
tests/benchmark/specificity_fixture.json (fetched once from NCBI, provenance recorded there).
Controlled positive/negative set with a confusion matrix.

  Negative controls (specific primers must clear): 8 geNorm assays vs an 8-gene mini-transcriptome
    -> each must hit ONLY its own gene, exactly one on-size product.
  Positive controls (off-target MUST be flagged):
    - a near-identical processed pseudogene (>=2 subjects)
    - a paralog with 2 internal mismatches, 3' intact (mismatch tolerance)
    - an off-size mispriming product (n_off_size >= 1)
  Discrimination negatives (must NOT amplify): a 3'-terminal-mismatch variant, and a random sequence.

Target: sensitivity = 100% (all off-targets caught) and specificity = 100% (all specific assays clear).
"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import specificity as SP

fails = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("" if ok else f"  [{detail}]"))
    if not ok: fails.append(name)

HERE = os.path.dirname(os.path.abspath(__file__))
fix = json.load(open(os.path.join(HERE, "benchmark", "specificity_fixture.json")))
seqs = fix["sequences"]

PRIMERS = {
    "ACTB": ("CTGGAACGGTGAAGGTGACA", "AAGGGACTTCCTGTAACAATGCA"),
    "B2M": ("TGCTGTCTCCATGTTTGATGTATCT", "TCTCTGCTCCCCACCTCTAAGT"),
    "GAPDH": ("TGCACCACCAACTGCTTAGC", "GGCATGGACTGTGGTCATGAG"),
    "HMBS": ("GGCAATGCGGCTGCAA", "GGGTACCCACGCGAATCAC"),
    "HPRT1": ("TGACACTGGCAAAACAATGCA", "GGTCCTTTTCACCAGCAAGCT"),
    "RPL13A": ("CCTGGAGGAGAAGAGGAAAGAGA", "TTGAGGACCTCTGTGTATTTGTCAA"),
    "SDHA": ("TGGGAACAAGAGGGCATCTG", "CCACCACTGCATCAAATTCATG"),
    "YWHAZ": ("ACTTTTGGTACATTGTGGCTTCAA", "CCGCCAGGACAAACCAGTAT"),
}

TP = FN = TN = FP = 0
mini = "".join(">%s\n%s\n" % (g, seqs[g]) for g in PRIMERS)

# ---- negative controls: each specific assay must hit only its own gene ----
for g, (F, R) in PRIMERS.items():
    r = SP.in_silico_pcr_offline(F, R, mini, max_mm=2)
    subs = sorted({p["subject"] for p in r["products"]})
    clean = (subs == [g] and r["n_products"] == 1)
    check("negative control clears: %s hits only itself" % g, clean, "subjects=%s" % subs)
    if clean: TN += 1
    else: FP += 1

# ---- positive control 1: near-identical processed pseudogene ----
gapdh = seqs["GAPDH"]; F, R = PRIMERS["GAPDH"]
r = SP.in_silico_pcr_offline(F, R, ">GAPDH\n%s\n>GAPDHP\n%s\n" % (gapdh, gapdh), max_mm=2)
caught = len({p["subject"] for p in r["products"]}) >= 2
check("positive: near-identical pseudogene flags >=2 subjects", caught, "n=%d" % r["n_products"])
TP += 1 if caught else 0; FN += 0 if caught else 1

# ---- positive control 2: paralog with 2 internal mismatches (3' intact) ----
actb = seqs["ACTB"]; F, R = PRIMERS["ACTB"]
r0 = SP.in_silico_pcr_offline(F, R, ">ACTB\n%s\n" % actb)
sp0, sp1 = r0["products"][0]["span"]; amp = list(actb[sp0:sp1 + 1])
amp[5] = "A" if amp[5] != "A" else "C"; amp[8] = "A" if amp[8] != "A" else "C"
decoy = "GGGATCCTTT" * 2 + "".join(amp) + "AAACCCGGGT" * 2
r = SP.in_silico_pcr_offline(F, R, ">ACTB\n%s\n>ACTB_paralog\n%s\n" % (actb, decoy), max_mm=2)
caught = "ACTB_paralog" in {p["subject"] for p in r["products"]}
check("positive: 2-mismatch paralog caught (mismatch tolerance)", caught, str(sorted({p["subject"] for p in r["products"]})))
TP += 1 if caught else 0; FN += 0 if caught else 1

# ---- positive control 3: off-size mispriming ----
hmbs = seqs["HMBS"]; F, R = PRIMERS["HMBS"]
rh = SP.in_silico_pcr_offline(F, R, ">HMBS\n%s\n" % hmbs)
h0, h1 = rh["products"][0]["span"]
misprime = hmbs[h0:h1 + 1] + "A" * 250 + SP._rc_iupac(R)
r = SP.in_silico_pcr_offline(F, R, ">HMBS_misprime\n%s\n" % misprime, max_mm=2)
caught = r["n_off_size"] >= 1
check("positive: off-size mispriming flagged", caught, "n_off_size=%d" % r["n_off_size"])
TP += 1 if caught else 0; FN += 0 if caught else 1

# ---- discrimination negative: 3'-terminal mismatch must NOT amplify ----
amp2 = list(actb[sp0:sp1 + 1]); k = len(PRIMERS["ACTB"][0]) - 1
amp2[k] = "A" if amp2[k] != "A" else "C"
dv = "GGGATCCTTT" * 2 + "".join(amp2) + "AAACCCGGGT" * 2
r = SP.in_silico_pcr_offline(PRIMERS["ACTB"][0], PRIMERS["ACTB"][1], ">ACTB\n%s\n>ACTB_3pvar\n%s\n" % (actb, dv), max_mm=2)
ok = "ACTB_3pvar" not in {p["subject"] for p in r["products"]}
check("discrimination: 3'-terminal-mismatch variant rejected", ok, str(sorted({p["subject"] for p in r["products"]})))
TN += 1 if ok else 0; FP += 0 if ok else 1

# ---- discrimination negative (MINUS strand): a 3'-terminal mismatch on the REVERSE primer's
# binding site must also block amplification. The reverse primer anneals to the plus strand
# (a '-'-orientation hit), so its 3' terminus maps to the START of its footprint; this exercises
# the strand-aware 3' anchor. Corrupt the base at the reverse primer's 3' end in a decoy copy.
Fp, Rp = PRIMERS["ACTB"]
rc_R = SP._rc_iupac(Rp)                          # reverse primer's footprint on the plus strand
amp3 = list(actb[sp0:sp1 + 1])
# the reverse footprint sits at the 3' end of the amplicon; its plus-strand form is rc_R.
# The reverse primer's 3' terminus aligns to the FIRST base of rc_R within the amplicon.
r_start = len(amp3) - len(rc_R)
amp3[r_start] = "A" if amp3[r_start] != "A" else "C"   # corrupt reverse primer's 3' terminus
dv2 = "GGGATCCTTT" * 2 + "".join(amp3) + "AAACCCGGGT" * 2
r = SP.in_silico_pcr_offline(Fp, Rp, ">ACTB\n%s\n>ACTB_rev3pvar\n%s\n" % (actb, dv2), max_mm=2)
ok = "ACTB_rev3pvar" not in {p["subject"] for p in r["products"]}
check("discrimination: reverse-primer 3'-mismatch (minus strand) rejected", ok, str(sorted({p["subject"] for p in r["products"]})))
TN += 1 if ok else 0; FP += 0 if ok else 1

# ---- discrimination negative: random sequence gives 0 products ----
random.seed(1)
rand = "".join(random.choice("ACGT") for _ in range(2000))
r = SP.in_silico_pcr_offline("CTGGAACGGTGAAGGTGACA", "AAGGGACTTCCTGTAACAATGCA", ">rand\n%s\n" % rand, max_mm=2)
ok = r["n_products"] == 0
check("discrimination: unrelated random sequence gives 0 products", ok, "n=%d" % r["n_products"])
TN += 1 if ok else 0; FP += 0 if ok else 1

# ---- confusion matrix ----
sens = TP / (TP + FN) if (TP + FN) else 0.0
spec = TN / (TN + FP) if (TN + FP) else 0.0
print("CONFUSION: TP=%d FN=%d TN=%d FP=%d  sensitivity=%.0f%% specificity=%.0f%%" % (TP, FN, TN, FP, 100*sens, 100*spec))
check("sensitivity == 100%% (all 3 off-targets flagged)", TP == 3 and FN == 0, "TP=%d FN=%d" % (TP, FN))
check("specificity == 100%% (all 11 specific controls clear)", TN == 11 and FP == 0, "TN=%d FP=%d" % (TN, FP))

if fails:
    print("SPECIFICITY VALIDATION FAILURES:", fails); sys.exit(1)
print("ALL SPECIFICITY VALIDATION ASSERTS PASS")
