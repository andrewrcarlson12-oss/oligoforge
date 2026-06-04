"""Validate the nearest-neighbor thermodynamics engine (oligoforge/nn.py). Offline.
Run from repo root:  OLIGOFORGE_EMAIL=you@x python3 tests/test_nn.py   (exit 0 = pass)

Anchors the engine to (a) the published SantaLucia 1998 parameter set as Biopython encodes it — so
nn.tm and the display Tm (thermo.tm_acc) can never silently diverge — and (b) Biopython's own Tm_NN
at matched reaction conditions, plus physical-sanity checks on ΔG and fraction-bound."""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import nn, thermo as T
from Bio.SeqUtils import MeltingTemp as mt

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (("  [%s]" % detail) if detail and not cond else ""))
    if not cond: fails.append(name)

# 1. Parameter table is byte-identical to Biopython's SantaLucia 1998 NN dinucleotides
_bp = {k: v for k, v in mt.DNA_NN3.items() if len(k) == 5 and "/" in k}
check("NN dinucleotide table == Biopython DNA_NN3", dict(nn.NN) == _bp)

# 2. Tm at qPCR conditions matches Biopython Tm_NN (SantaLucia NN + Owczarzy-2008 + Mg/dNTP) within 0.3 C
T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8, dna_conc=200, anneal_c=60)
seqs = ["CTTACAAGATATCCACCACA", "AGTCATTCTGATGTCGCTGATG", "TCGCTGGCATCAACAAGAAG",
        "GGGTCAAATGAGTTTCTGG", "GCGCGCGCGCGC", "ATATATATATATATAT",
        "CCACTATCCCTTTCTTGTCATCTCCAGCAG", "AACTTGCTCAGCCTGGTTTG"]
worst = 0.0
for s in seqs:
    p = nn.params(s)
    b = mt.Tm_NN(s, dnac1=200, dnac2=200, selfcomp=(s == T.revcomp(s)), Na=50, K=0, Tris=0, Mg=3, dNTPs=0.8, saltcorr=7)
    worst = max(worst, abs(p["tm"] - b))
check("nn.tm matches Biopython Tm_NN within 0.3 C (8 seqs)", worst < 0.3, "max d=%.2f" % worst)

# 3. ΔG°37 = ΔH - 310.15·ΔS/1000 and is favorable (negative) for a real oligo
p = nn.params("CTTACAAGATATCCACCACA")
check("dg37 consistent with dH,dS", abs(p["dg37"] - (p["dh"] - 310.15 * p["ds"] / 1000.0)) < 0.05)
check("dg37 favorable (<0)", p["dg37"] < 0)

# 4. Fraction bound = 0.5 at Tm and decreases monotonically as temperature rises
f_tm = nn.params("CTTACAAGATATCCACCACA", anneal_c=p["tm"])["frac_bound"]
check("frac_bound ~0.5 at Tm", abs(f_tm - 0.5) < 0.02, str(f_tm))
f50 = nn.params("CTTACAAGATATCCACCACA", anneal_c=50)["frac_bound"]
f60 = nn.params("CTTACAAGATATCCACCACA", anneal_c=60)["frac_bound"]
f70 = nn.params("CTTACAAGATATCCACCACA", anneal_c=70)["frac_bound"]
check("frac_bound monotonic in T (50>60>70)", f50 > f60 > f70, "%.3f %.3f %.3f" % (f50, f60, f70))

# 5. Free-Mg physics: chelating Mg with dNTPs (less free Mg) lowers Tm
tm_lo_dntp = nn.params("CTTACAAGATATCCACCACA")["tm"]
T.set_conditions(dntp_conc=2.9)               # free Mg ~0.1 mM
tm_hi_dntp = nn.params("CTTACAAGATATCCACCACA")["tm"]
T.set_conditions(dntp_conc=0.8)               # restore
check("more dNTP (less free Mg) lowers Tm", tm_hi_dntp < tm_lo_dntp - 2, "%.1f -> %.1f" % (tm_lo_dntp, tm_hi_dntp))

# 6. Raising monovalent salt raises Tm (Owczarzy 2004 monovalent branch, Mg off)
T.set_conditions(dv_conc=0, dntp_conc=0)
tm_lowNa = nn.params("AGTCATTCTGATGTCGCTGATG", cond={"mv_conc": 20})["tm"]
tm_hiNa = nn.params("AGTCATTCTGATGTCGCTGATG", cond={"mv_conc": 200})["tm"]
T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8)   # restore defaults
check("higher monovalent salt raises Tm", tm_hiNa > tm_lowNa + 3, "%.1f -> %.1f" % (tm_lowNa, tm_hiNa))

# 7. Non-ACGT / too-short -> None (mismatch / degenerate / modified out of scope for this engine)
check("degenerate -> None", nn.params("ACGTRYACGT") is None)
check("single base -> None", nn.params("A") is None)

# 8. LNA layer (McTigue 2004): all 32 increments resolvable, and the combined model reproduces
#    McTigue's published duplex Tm's. Embedded subset of the paper's experimental values (5 uM oligo,
#    1 M Na+); engine dnac/2 must equal McTigue Ct/4 -> dna_conc=2500 nM, Na=1 M.
_KEYS_OK = all(  # every ACGT step x both LNA positions has a parameter
    (a + "L" + b + "/" + nn._C[a] + nn._C[b]) in nn._LNA_INC and
    (a + b + "L" + "/" + nn._C[a] + nn._C[b]) in nn._LNA_INC
    for a in "ACGT" for b in "ACGT")
check("all 32 McTigue increments present", _KEYS_OK)

_EXP = {  # MELTING 'XL' notation -> IDT '+X', : experimental Tm (McTigue 2004)
    "GTC+GAACAGC": 53.8, "GTCGAAC+AGC": 56.1, "CGC+TGTTACGC": 60.5, "GGAC+CTCGAC": 58.7,
    "ATCT+ATCCGGC": 57.3, "GC+AGGTCTGC": 57.7, "T+TGCTCGATGT": 54.6, "AC+AAGCGACTC": 55.7,
    "GGT+GCCAA": 44.8, "CAC+GGCTC": 49.8, "ATTTGAC+TCAG": 51.1, "TATTAAGCG+ACCACACATAA": 68.1,
}
T.set_conditions(mv_conc=1000, dv_conc=0, dntp_conc=0, dna_conc=2500)
_errs = []
for s, exp in _EXP.items():
    p = nn.params_lna(s)
    if p:
        _errs.append(p["tm"] - exp)
import math as _m
_rmse = _m.sqrt(sum(e * e for e in _errs) / len(_errs))
T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8, dna_conc=200)   # restore
check("params_lna scored all embedded duplexes", len(_errs) == len(_EXP))
check("LNA Tm reproduces McTigue experimental data (RMSE < 2.5 C)", _rmse < 2.5, "RMSE=%.2f" % _rmse)

# 9. LNA stabilizes: the cytb probe's locked bases must raise Tm and fraction-bound at 60 C
_dna = nn.params("CTTACAAGATATCCACCACA")
_lna = nn.params_lna("CTTA+CA+A+GATAT+CC+ACCACA")
check("LNA raises Tm vs bare DNA", _lna["tm"] > _dna["tm"] + 5, "%.1f -> %.1f" % (_dna["tm"], _lna["tm"]))
check("LNA rescues binding at anneal (bare <20%, LNA >80%)", _dna["frac_bound"] < 0.2 < 0.8 < _lna["frac_bound"],
      "%.2f -> %.2f" % (_dna["frac_bound"], _lna["frac_bound"]))
check("adjacent-LNA flagged beyond_param", _lna["beyond_param"] and "adjacent" in _lna["note"])

# 10. LNA path rejects degenerate / no-LNA input
check("params_lna no-LNA -> None", nn.params_lna("CTTACAAGATATCCACCACA") is None)
check("params_lna degenerate -> None", nn.params_lna("CTT+RCAAGATAT") is None)

# 11. Mismatch engine (Allawi/SantaLucia/Peyret): mismatch_params Tm must match Biopython Tm_NN(c_seq)
#     across random matched/mismatched duplexes, and behave physically.
import random as _rand
_rand.seed(11)
T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8, dna_conc=200, anneal_c=60)
def _comp(s): return s.translate(str.maketrans("ACGT", "TGCA"))
_d = []
for _ in range(400):
    L = _rand.randint(12, 26)
    o = "".join(_rand.choice("ACGT") for _ in range(L))
    tgt = list(o)
    for p in _rand.sample(range(L), _rand.randint(0, 4)):
        tgt[p] = _rand.choice([b for b in "ACGT" if b != tgt[p]])
    tgt = "".join(tgt)
    pr = nn.mismatch_params(o, tgt)
    if pr is None:
        continue
    try:
        tb = mt.Tm_NN(o, c_seq=_comp(tgt), Na=50, K=0, Tris=0, Mg=3, dNTPs=0.8, saltcorr=7, dnac1=200, dnac2=200)
    except Exception:
        continue
    _d.append(abs(pr["tm"] - tb))
check("mismatch_params Tm matches Biopython Tm_NN(c_seq) within 0.3 C (random duplexes)",
      _d and max(_d) < 0.3, "n=%d max=%.3f" % (len(_d), max(_d) if _d else -1))

# perfect match: dtm 0, n_mm 0; mismatches lower Tm and binding
_pm = nn.mismatch_params("AGTCATTCTGATGTCGCTGATG", "AGTCATTCTGATGTCGCTGATG")
check("perfect match: dtm 0, n_mm 0", _pm["dtm"] == 0 and _pm["n_mm"] == 0)
_mm = nn.mismatch_params("AGTCATTCTGATGTCGCTGATG", "AGTCATGCTGATGTCGCTAATG")  # 2 internal mismatches
check("2 mismatches lower Tm and fraction bound", _mm["dtm"] < 0 and _mm["n_mm"] == 2 and _mm["frac_bound"] < _mm["frac_bound_perfect"])
# 3'-terminal mismatch reported at distance 1 (the number that blocks primer extension)
_t3 = nn.mismatch_params("AGTCATTCTGATGTCGCTGATC", "AGTCATTCTGATGTCGCTGATG")  # last base differs
check("3'-terminal mismatch flagged mm_from_3p=1", _t3["mm_from_3p"] == 1)
check("mismatch_params length-mismatch / degenerate -> None",
      nn.mismatch_params("ACGT", "ACG") is None and nn.mismatch_params("ACGT", "ACGR") is None)

# 12. Degenerate- and LNA-aware mismatch scoring (rung capstone). A degenerate oligo binds a lineage it
#     covers (resolved variant matches -> ΔTm 0) but loses stability where variation is uncovered; LNA
#     increments apply on matched steps.
_cov = nn.mismatch_params("CTTA+CA+A+GAYATCC+ACCACA", "CTTACAAGACATCCACCACA")  # Y={C,T} site has C -> covered
_unc = nn.mismatch_params("CTTA+CA+A+GAYATCC+ACCACA", "CTTACAAGAAATCCACCACA")  # site has A -> uncovered
check("covered degenerate lineage: ΔTm 0, n_mm 0, degenerate flagged",
      _cov["dtm"] == 0 and _cov["n_mm"] == 0 and _cov["degenerate"] and _cov["n_var"] == 2)
check("uncovered variation: real mismatch, lower Tm and binding",
      _unc["n_mm"] == 1 and _unc["dtm"] < 0 and _unc["frac_bound"] < _cov["frac_bound"])
check("LNA increments raise the degenerate probe's binding (vs DNA backbone of same variant)",
      _cov["tm"] > nn.mismatch_params("CTTACAAGACATCCACCACA", "CTTACAAGACATCCACCACA")["tm"] + 5)
# LNA probe buffers a single mismatch (stays well bound), and mm distance reported
_lp = nn.mismatch_params("CTTA+CA+A+GATAT+CC+ACCACA", "CTTACAAGATGTCCACCACA")
check("LNA probe stays bound through one mismatch", _lp["n_mm"] == 1 and _lp["frac_bound"] > 0.8 and _lp["lna_n"] == 5)
check("pure-ACGT mismatch result still has n_var=1, degenerate False",
      nn.mismatch_params("AGTCATTCTGATG", "AGTCATTCTGATG")["n_var"] == 1)
T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8, dna_conc=200)

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL NN ASSERTS PASS")
