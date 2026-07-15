"""Hard regression asserts. Offline (no network). Run from repo root:  python tests/test_regression.py
Pins the locked science so a code change that shifts a number fails loudly."""
import sys; sys.path.insert(0, ".")
from oligoforge import thermo as T, design as D, profiles as P, conservation as C, specificity as SP, multiplex as MX

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

# 15b. QC accepts a pasted IDT order string / LNA oligo (strips mod notation, scores DNA backbone)
try:
    from fastapi.testclient import TestClient as _TC3
    import app as _APP3
    _c3 = _TC3(_APP3.app)
    _q = _c3.post("/api/qc", json={"seq": "/56-FAM/CTTA+CA+A+GATAT+CC+ACCACA/3IABkFQ/", "role": "probe"}).json()
    check("QC parses IDT order string (no error)", not _q.get("error"))
    check("QC recovers bare 20-mer from order string", _q.get("seq") == "CTTACAAGATATCCACCACA")
    check("QC warns LNA Tm is DNA-backbone", bool(_q.get("lna_note")))
    _q2 = _c3.post("/api/qc", json={"seq": "ACGTACGTACGTACGTACGT", "role": "primer"}).json()
    check("QC plain oligo unaffected (no lna_note)", not _q2.get("lna_note"))
except Exception as _e:
    check("QC order-string handling", False, str(_e))

# 16. NCBI api_key wiring: per-request, env fallback, and clearing
try:
    import os as _os
    _os.environ.pop("OLIGOFORGE_NCBI_KEY", None)
    import app as _APP
    from oligoforge import ncbi as _NCBI
    _APP._set_email("e@x.com", "KEY123")
    check("api_key set from request", _NCBI.Entrez.api_key == "KEY123")
    _APP._set_email("e@x.com", None)
    check("api_key cleared when none given", _NCBI.Entrez.api_key is None)
    _os.environ["OLIGOFORGE_NCBI_KEY"] = "ENVKEY"
    _APP._set_email("e@x.com", None)
    check("api_key falls back to env", _NCBI.Entrez.api_key == "ENVKEY")
    _APP._set_email("e@x.com", "REQWINS")
    check("request key beats env", _NCBI.Entrez.api_key == "REQWINS")
    _os.environ.pop("OLIGOFORGE_NCBI_KEY", None)
    _APP._set_email("e@x.com", None)
except Exception as _e:
    check("api_key wiring", False, str(_e))

print()

# --- v1.5.4: gapped cross-check + multiplex melt overlap ---
_oli="ACGTACGTACGTACGT"
_di=C.discrimination(_oli, [_oli[:8]+"T"+_oli[8:]])   # 1-bp insertion, otherwise identical
check("gapped cross-check flags indel-masked off-target", _di.get("indel_masked",0)>=1, _di.get("indel_masked"))
_ds2=C.discrimination(_oli, ["ACGTACGTAAGTAAGT"])      # substitutions only, no indel
check("no false indel flag on substitutions", _ds2.get("indel_masked",0)==0, _ds2.get("indel_masked"))
_mxm=MX.check([{"name":"A","sybr":True,"amplicon_tm":78.0,"oligos":[{"name":"F","seq":"ACGTACGTACGTACGT"},{"name":"R","seq":"TTTTGGGGCCCCAAAA"}]},
               {"name":"B","sybr":True,"amplicon_tm":78.4,"oligos":[{"name":"F","seq":"GGGGCCCCAAAATTTT"},{"name":"R","seq":"CCCCAAAATTTTGGGG"}]}])
check("multiplex flags SYBR melt overlap", len(_mxm.get("melt_overlaps",[]))==1, _mxm.get("melt_overlaps"))
_mxs=MX.check([{"name":"A","sybr":True,"amplicon_tm":78.0,"oligos":[{"name":"F","seq":"ACGTACGTACGTACGT"}]},
               {"name":"B","sybr":True,"amplicon_tm":84.0,"oligos":[{"name":"F","seq":"GGGGCCCCAAAATTTT"}]}])
check("multiplex passes well-separated SYBR amplicons", len(_mxs.get("melt_overlaps",[]))==0, _mxs.get("melt_overlaps"))


# --- v1.16/v1.17: anneal-temp thermo, logistic LOD, raw Cq, melt peaks, ensemble structure ---
import math as _mm
_seq_sd = "GGGGCCCCAAAATTTTGGGGCCCC"
_sdf = _TH.self_dimer_full(_seq_sd)
check("self_dimer_full dG@37 == legacy self_dimer", abs(_sdf[0] - _TH.self_dimer(_seq_sd)) < 1e-9, _sdf)
check("self_dimer_full anneal dG less negative than 37C", _sdf[1] > _sdf[0], _sdf)
_hpf = _TH.hairpin_full("TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT")
check("hairpin_full dG@37 == legacy hairpin", abs(_hpf[0] - _TH.hairpin("TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT")[0]) < 1e-9, _hpf)
check("end_stability returns float", isinstance(_TH.end_stability("AGTCATTCTGATGTCGCTGATG", "ACCTGTCAGTGTTTTCAAGCA"), float))

from oligoforge import quant as _Q
_ptsT = []
for _q, _nd in [(1e5, 5), (1e4, 5), (1e3, 5), (1e2, 4), (1e1, 1), (1, 0)]:
    for _i in range(5):
        _ptsT.append((_q, (40 - 3.32 * _mm.log10(_q)) if _i < _nd else None))
check("lod95 finite on a detection transition", _Q.standard_curve(_ptsT)["lod95"] is not None)
check("lod95 None when every level fully detected",
      _Q.standard_curve([(q, 40 - 3.32 * _mm.log10(q)) for q in [1e5, 1e4, 1e3, 1e2, 1e1] for _ in range(3)])["lod95"] is None)

from oligoforge import cq as _CQ
_cq = _CQ.analyze([50 + 1000 / (1 + _mm.exp(-0.6 * (i - 25))) for i in range(1, 46)], threshold=100.0)
check("cq.analyze finite Cq + amplified on a sigmoid", _cq["cq_threshold"] is not None and _cq["amplified"] is True, _cq)

from oligoforge import melt as _ML
_mt = [65 + 0.5 * i for i in range(62)]
_mr = _ML.analyze([1000 / (1 + _mm.exp(0.8 * (t - 85))) for t in _mt], _mt)
check("melt.analyze single peak near 85C", _mr["n_peaks"] == 1 and abs(_mr["dominant_tm"] - 85) <= 1.0, _mr)

from oligoforge import structure as _STR
if _STR.available():
    _amp_s = "GAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCCATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGG"
    _fe = _STR.fold_ensemble(_amp_s, anneal_c=60.0); _f0 = _STR.fold(_amp_s)
    check("fold_ensemble MFE paired identical to fold()", _fe["paired"] == _f0["paired"])
    check("fold_ensemble anneal mfe >= 37C mfe", _fe["mfe_anneal"] >= _fe["mfe"], (_fe["mfe"], _fe["mfe_anneal"]))
    check("fold_ensemble paired_prob length matches amplicon", len(_fe["paired_prob"]) == len(_amp_s))
else:
    check("structure ensemble (ViennaRNA absent -- skipped)", True)

# ---- v1.27.0: structure gates evaluated at the annealing temperature ----
# The reject gates in design.py now judge hairpins/dimers at ANNEAL_C (54-60 C), not primer3's
# 37 C default. Pin (a) the hand-validated golden designs are UNCHANGED by the switch, and
# (b) the switch actually admits the melted-structure candidates the 37 C gate wrongly discarded.
import os as _os2, json as _json2
_HMBS = ("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
         "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
T.set_conditions(anneal_c=60)
_a = D.design_assay(_HMBS, P.PROFILES["idt_taqman"])
check("gate@anneal: HMBS golden forward unchanged", _a["forward"] == "GAGCTATACCCCGACCTCTG", _a["forward"])
check("gate@anneal: HMBS golden reverse unchanged", _a["reverse"] == "CTTCTCTCCAATCTTGGAAAGCG", _a["reverse"])
check("gate@anneal: HMBS golden amplicon 93 unchanged", _a["amplicon"] == 93, _a["amplicon"])
_fw, _rv = D.enumerate_primers(_HMBS, P.PROFILES["idt_taqman"])
check("gate@anneal: admits more candidates than the 37C gate (HMBS >150 fwd)", len(_fw) > 150, len(_fw))
# a primer whose hairpin is real at 37 C but fully melted at 60 C must now be ACCEPTED
_melt = "GAGCTATACCCCGACCTCTG"   # golden fwd: passes both; a positive control that acceptance still works
check("gate@anneal: golden fwd still passes _ok_primer", D._ok_primer(_melt, P.PROFILES["idt_taqman"]))

# ---- v1.27.0: reported Tm is divalent-aware (Owczarzy-2008 + von Ahsen free Mg) ----
# The salt-model audit: the DISPLAYED Tm must respond to Mg2+ and to dNTP chelation; the SELECTION
# Tm (primer3) is near-Mg-insensitive by design. Pin both so a future salt-model edit is caught.
_S = "AGTCATTCTGATGTCGCTGATG"
T.set_conditions(mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, anneal_c=60)
_acc_mg3 = T.tm_acc(_S); _sel_mg3 = T.tm(_S)
T.set_conditions(dv_conc=6.0)
_acc_mg6 = T.tm_acc(_S); _sel_mg6 = T.tm(_S)
T.set_conditions(dv_conc=3.0, dntp_conc=0.0)
_acc_dntp0 = T.tm_acc(_S)
T.set_conditions(mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, anneal_c=60)   # restore
check("salt: reported Tm rises with Mg2+ (3->6 mM, >=0.8 C)", _acc_mg6 - _acc_mg3 >= 0.8, (_acc_mg3, _acc_mg6))
check("salt: dNTP chelation lowers reported Tm (0.8 mM dNTP drops it >=0.2 C)",
      _acc_dntp0 - _acc_mg3 >= 0.2, (_acc_mg3, _acc_dntp0))
check("salt: selection Tm is near-Mg-insensitive by design (3->6 mM, <0.4 C)",
      abs(_sel_mg6 - _sel_mg3) < 0.4, (_sel_mg3, _sel_mg6))
check("salt: reported Tm reads above selection Tm (Owczarzy-2008 vs primer3, ~1-2 C)",
      1.0 <= _acc_mg3 - _sel_mg3 <= 2.5, (_sel_mg3, _acc_mg3))

# ---- intron_check returns a verdict on the no-record path without requiring live NCBI ----
_orig_exons = SP.exon_junctions_mrna
try:
    SP.exon_junctions_mrna = lambda *a, **k: (None, "offline fixture: no exon structure", None)
    _ir = SP.intron_check("NOSUCHGENE___", "Nonexistent organism", 1, 100)
finally:
    SP.exon_junctions_mrna = _orig_exons
check("intron: 'verdict' present on the degraded path (no KeyError)", "verdict" in _ir, sorted(_ir))
check("intron: degraded verdict is a non-empty string", isinstance(_ir.get("verdict"), str) and _ir["verdict"], _ir.get("verdict"))

# ---- v1.27.0: Plasmodium-vs-Haemoproteus conservation/discrimination reproduces (offline fixtures) ----
_FX = _os2.path.join(_os2.path.dirname(_os2.path.abspath(__file__)), "fixtures")
try:
    _plas = _json2.load(open(_os2.path.join(_FX, "plasmodium_cytb.json")))["sequences"]
    _haem = _json2.load(open(_os2.path.join(_FX, "haemoproteus_cytb.json")))["sequences"]
except Exception:
    _plas = _haem = None
if _plas and _haem:
    _cons = C.conservation("CTTACAAGATATCCACCACA", _plas)         # locked Plas probe
    _mi = _cons["mean_ident"] if isinstance(_cons, dict) else None
    check("cons: locked Plas probe >=98% conserved across Plasmodium", _mi is not None and _mi >= 98.0, _mi)
    _disc = C.discrimination("CTTACAAGATATCCACCACA", _haem)        # vs Haemoproteus off-target
    check("cons: Plas probe >=3 mismatch to closest Haemoproteus", (_disc.get("min_mismatch") or 0) >= 3, _disc.get("min_mismatch"))
    check("cons: Plas probe max off-target identity <=90%", (_disc.get("max_ident") or 100) <= 90.0, _disc.get("max_ident"))
else:
    check("cons: Plasmodium/Haemoproteus fixtures present", False, "fixtures missing")

# ---- v1.27.1: secondary-structure REJECTION is gated at each ASSAY'S anneal temp, not a session global ----
# The host TaqMan/SYBR/GC panel anneals at 60 C; the AT-rich parasite mtDNA assays at ~54 C. Before this
# fix every profile was judged at the module global (default 60 C), so the 54 C parasite assays were gated
# 6 C too hot and under-rejected structure that survives at 54 C. Guard: profiles carry the right temp, the
# gate honors it per-assay, and it no longer depends on the mutable global (thread-safe under the endpoint pool).
for _k, _p in P.PROFILES.items():
    _want = 54.0 if _k.startswith("parasite") else 60.0
    check(f"anneal_c on profile '{_k}' == {_want:.0f} C", _p.get("anneal_c") == _want, _p.get("anneal_c"))

# A self-dimer that persists at 54 C but is melted by 60 C: REJECTED by the parasite profile, ADMITTED by the
# host profile. Same -6.0 self-dimer floor on both, so the ONLY thing that can differ is the anneal temperature.
_FLIP = "ACAGGAAGCTTCATGGCCTAT"
_par, _host = P.PROFILES["parasite_mtdna"], P.PROFILES["idt_taqman"]
_sd54 = T.self_dimer_full(_FLIP, _par["anneal_c"])[1]
_sd60 = T.self_dimer_full(_FLIP, _host["anneal_c"])[1]
check("gate: flip-oligo self-dimer rejected at parasite 54 C", _sd54 <= _par["self_dimer_min"], f"{_sd54:.2f}")
check("gate: same oligo admitted at host 60 C (per-assay, not global)", _sd60 > _host["self_dimer_min"], f"{_sd60:.2f}")

# Global-independence: perturbing the session anneal global must NOT move a profile-temp evaluation.
_before = T.self_dimer_full(_FLIP, 54.0)[1]
T.set_conditions(anneal_c=37.0)
_after = T.self_dimer_full(_FLIP, 54.0)[1]
T.set_conditions(anneal_c=60.0)   # restore
check("gate: 54 C evaluation is independent of the session global", abs(_before - _after) < 1e-9, (_before, _after))

# The parasite autodesign winner is preserved now that the gate genuinely runs at 54 C (design gated through
# the profile's anneal_c, not the default global). Byte-identical to the pinned discrimination winner.
if _plas and _haem:
    from oligoforge import autodesign as _AD
    _d54 = _AD.design_from_sequences(_plas, P.PROFILES["parasite_mtdna"], offs=_haem, n_candidates=5)
    _w = (_d54.get("candidates") or [{}])[0].get("assay", {})
    check("autodesign: parasite disc winner preserved at 54 C gate",
          _w.get("forward") == "TTTCCATTTATAGCCTTATGTATTG" and _w.get("amplicon") == 96,
          (_w.get("forward"), _w.get("amplicon")))

if fails: print("REGRESSION FAILURES:", fails); sys.exit(1)
print("ALL REGRESSION ASSERTS PASS")
