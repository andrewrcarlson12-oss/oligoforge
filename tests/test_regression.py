"""Hard regression asserts. Offline (no network). Run from repo root:  python tests/test_regression.py
Pins the locked science so a code change that shifts a number fails loudly."""
import sys; sys.path.insert(0, ".")
from oligoforge import thermo as T, design as D, profiles as P, conservation as C, specificity as SP

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond: fails.append(name)

# 1. Locked host-primer Tms stay in the design window (60-62.5 C at qPCR salt)
host = {
 "IFNG_F":"AGTCATTCTGATGTCGCTGATG","IFNG_R":"ACCTGTCAGTGTTTTCAAGCA",
 "IL4_F":"AACTTGCTCAGCCTGGTTTG","IL4_R":"ATTCTTTAGTGAGGTGGTGCTG",
 "RPL13_F":"TCGCTGGCATCAACAAGAAG","RPL13_R":"TCGGGAAGAGGATGAGCTTG",
 "YWHAZ_F":"CCGTTACTTGGCTGAGGTTG","YWHAZ_R":"GATGGGATGTGTTGGTTGCA",
}
for k,s in host.items():
    tm=T.tm(s); check(f"Tm {k} in [60,62.5]", 60.0<=tm<=62.5, f"{tm:.2f}")

# 2. Probe Tms sit >=5 C above their primers
for probe,primtm in [("TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT",62),("CTTGTGCCCTGCTCTGGTCCC",62),
                     ("AACAAGTCCACCGAGTCCCTGCA",62),("CCACTATCCCTTTCTTGTCATCTCCAGCAG",62)]:
    check(f"probe Tm >= primer+5 ({probe[:10]}..)", T.tm(probe) >= primtm+5, f"{T.tm(probe):.1f}")

# 3. design_assay reproduces the hand-built HMBS assay exactly (engine anchor)
HMBS=("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
      "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
a=D.design_assay(HMBS, P.PROFILES["idt_taqman"])
check("HMBS forward", a and a["forward"]=="GAGCTATACCCCGACCTCTG", a and a["forward"])
check("HMBS reverse", a and a["reverse"]=="CTTCTCTCCAATCTTGGAAAGCG", a and a["reverse"])
check("HMBS amplicon 93", a and a["amplicon"]==93, a and a["amplicon"])

# 4. in-silico PCR combiner: 157 bp target found, far/same-strand excluded
hits=[dict(primer="F",subject="t",lo=100,hi=120,strand="+",q3=True),
      dict(primer="R",subject="t",lo=240,hi=256,strand="-",q3=True),
      dict(primer="F",subject="far",lo=100,hi=120,strand="+",q3=True),
      dict(primer="R",subject="far",lo=9000,hi=9016,strand="-",q3=True)]
prod=SP.epcr(hits,max_product=3000)
check("ePCR finds 157 bp", any(p["size"]==157 for p in prod), [p["size"] for p in prod])
check("ePCR excludes 8917 bp", not any(p["subject"]=="far" for p in prod))

# 5. conservation: perfect-match target set -> 100%; off-target mismatch detected
cons=C.conservation("ACGTACGT", ["TTACGTACGTAA","CCACGTACGTGG"], min_ident=0.5)
check("conservation perfect=100%", cons["mean_ident"]==100.0, cons["mean_ident"])
disc=C.discrimination("ACGTACGT", ["GGACGTTCGTAA"])
check("discrimination counts mismatch", disc["n"]==1 and disc["min_mismatch"]>=1, disc)

# 6. LNA Tm: backbone + range, range above backbone
l=T.tm_lna("CTTACAAGATATCCACCACA", n_lna=4)
check("LNA range above backbone", l["est_tm_high"]>l["dna_backbone_tm"]>0, l)

# 7. degenerate-base safety (Plas_F carries a W) — must not crash
check("degenerate Tm computes", isinstance(T.tm("TACCTGGACTWGTTTCATGG"), float))
check("degenerate revcomp safe", set(T.revcomp("ACGTW")) <= set("ACGTWN"))
# 8. SYBR (no-probe) profile designs primers with probe=None
_sb="ATGGGTTATGTATTACCTTGGACTAGTTTCATGGTTTACAAGATATCCACCACATTTGGGTCACTTACAAGATATCCAAGCTTGGATCCGAATTC"
_sa=D.design_assay(_sb, P.PROFILES["sybr_generic"])
check("SYBR designs primers, no probe", bool(_sa) and bool(_sa["forward"]) and _sa["probe"] is None, _sa)

# 9. clean_seq: cleaning preserves clean input, strips junk, preserves degenerate, rejects bad
from oligoforge import thermo as _T
check("clean_seq no-op on clean", _T.clean_seq("GAGCTATACCCCGACCTCTG")[0]=="GAGCTATACCCCGACCTCTG")
check("clean_seq strips FASTA+space+case", _T.clean_seq(">x\n gagct atac ")[0]=="GAGCTATAC")
check("clean_seq RNA->DNA", _T.clean_seq("GAGCUAUACC")[0]=="GAGCTATACC")
check("clean_seq keeps degenerate W", _T.clean_seq("TACCTGGACTWGTTTCATGG")[0]=="TACCTGGACTWGTTTCATGG")
check("clean_seq rejects bad char", _T.clean_seq("GAGCTXZ")[2] is not None)

# 10. refgenes: unstable gene ranks last; geNorm M matches pairwise formula
import math, statistics as _st
from oligoforge import refgenes as _RG
_rg=_RG.analyze({"REF_A":[20,20,20,20],"REF_B":[22,22,22,22],"BAD":[18,24,19,26]})
_order=[g["gene"] for g in _rg["ranking"]]; _M={g["gene"]:g["M"] for g in _rg["ranking"]}
_qbad=[2.0**(18-c) for c in [18,24,19,26]]; _V=_st.stdev([-math.log2(x) for x in _qbad])
check("refgenes ranks unstable last", _order[-1]=="BAD", _order)
check("refgenes M matches pairwise formula", abs(_M["BAD"]-round(_V,4))<1e-4, (_M["BAD"],_V))
check("refgenes stable pair M~0", _M["REF_A"]<1e-6 and _M["REF_B"]<1e-6, _M)

# 11. MIQE report builds HTML + CSV with checklist marks
from oligoforge import report as _RPT
_rep=_RPT.build([{"name":"X","gene":"IFNG","organism":"A. coerulescens","forward":"AGTCATTCTGATGTCGCTGATG","reverse":"ACCTGTCAGTGTTTTCAAGCA","probe":"TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT","amplicon":136,"validation":{"efficiency_pct":99.2,"r2":0.998}}])
check("report builds html+csv", "IFNG" in _rep["html"] and "99.2%" in _rep["html"] and _rep["csv"].startswith("name,gene"), _rep.get("n_assays"))

# 12. multiplex flags a shared dye and a cross-assay dimer
from oligoforge import multiplex as _MX
_mx=_MX.check([{"name":"A","dye":"FAM","oligos":[{"name":"F","seq":"GGGGGGGGGGCCCCCCCCCC"}]},
               {"name":"B","dye":"FAM","oligos":[{"name":"R","seq":"GGGGGGGGGGCCCCCCCCCC"}]}])
check("multiplex flags shared dye", len(_mx["channel_conflicts"])==1, _mx["channel_conflicts"])
check("multiplex flags cross-dimer", _mx["n_flagged"]>=1, _mx["n_flagged"])

# 13. degenerate Tm reported as a real range; non-degenerate collapses to one
from oligoforge import thermo as _TH
_tr=_TH.tm_range("TACCTGGACTWGTTTCATGG")
check("tm_range degenerate spread", _tr["degenerate"] and _tr["n"]==2 and _tr["max"]>_tr["min"], _tr)
check("tm_range non-degenerate single", _TH.tm_range("AGTCATTCTGATGTCGCTGATG")["n"]==1)
# 14. set_conditions retargets Tm then restores
_b=_TH.tm("AGTCATTCTGATGTCGCTGATG"); _TH.set_conditions(dv_conc=8.0); _hi=_TH.tm("AGTCATTCTGATGTCGCTGATG"); _TH.set_conditions(dv_conc=3.0)
check("set_conditions changes Tm", round(_hi,2)!=round(_b,2), (round(_b,2),round(_hi,2)))

# 15. project save/list/load/delete roundtrip on a non-jay panel (organism-agnostic)
try:
    from fastapi.testclient import TestClient as _TC
    import app as _APP2
    _c2 = _TC(_APP2.app)
    _pj = [{"name": "x", "gene": "IL1B", "organism": "Lithobates pipiens",
            "forward": "ACGTACGTACGTACGTAC", "reverse": "TGCATGCATGCATGCATG", "amplicon": 110}]
    check("project save", _c2.post("/api/project/save", json={"name": "__rtest", "assays": _pj}).json().get("n") == 1)
    check("project load roundtrip", _c2.post("/api/project/load", json={"name": "__rtest"}).json()["assays"][0]["organism"] == "Lithobates pipiens")
    check("project delete", _c2.post("/api/project/delete", json={"name": "__rtest"}).json().get("deleted") == "__rtest")
except Exception as _e:
    check("project roundtrip", False, str(_e))

print()
if fails: print("REGRESSION FAILURES:", fails); sys.exit(1)
print("ALL REGRESSION ASSERTS PASS")
