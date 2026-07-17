"""Offline invariant and versioned-winner gate for discrimination ranking.

The expected winners are not treated as biologically universal optima. They are a
versioned regression record for the current fully annotated evidence hierarchy. A
winner may change only with a documented validation update, not merely to satisfy a
legacy exact-sequence test.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from oligoforge import autodesign as AD, profiles as P

HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixtures")
fails=[]
def check(name,cond,detail=""):
    print(("  PASS " if cond else "  FAIL ")+name+("" if cond else "   [%s]" % (detail,)),flush=True)
    if not cond:fails.append(name)

plas=json.load(open(os.path.join(FX,"plasmodium_cytb.json")))["sequences"]
haem=json.load(open(os.path.join(FX,"haemoproteus_cytb.json")))["sequences"]
expected=json.load(open(os.path.join(FX,"autodesign_expected_v1.37.json")))
prof=P.PROFILES["parasite_mtdna"]
check("fixtures load",len(plas)==11 and len(haem)==12,(len(plas),len(haem)))

d=AD.design_from_sequences(plas,prof,offs=haem,n_candidates=5,objective="discrimination")
check("discrimination design returns",not d.get("error") and d.get("candidates"),d.get("error"))
if d.get("candidates"):
    check("every finalist received full annotation",all(c.get("evidence",{}).get("evaluations",{}).get("target_epcr") and c.get("evidence",{}).get("evaluations",{}).get("offtarget_epcr") for c in d["candidates"]))
    check("all returned finalists are hard-valid",all(c["evidence"]["hard_valid"] for c in d["candidates"]),[c["evidence"]["hard_failures"] for c in d["candidates"]])
    check("rank 1 has no Haemoproteus product",d["candidates"][0]["evidence"]["offtarget"]["product_subjects"]==0,d["candidates"][0]["evidence"]["offtarget"])
    check("rank 1 preserves declared target coverage",d["candidates"][0]["evidence"]["target_coverage"]>=0.90,d["candidates"][0]["evidence"]["target_coverage"])
    check("run manifest is versioned",d["ranker_manifest"]["ranker_version"]==expected["ranker_version"],d["ranker_manifest"])
    check("search corpus policy is versioned",d["ranker_manifest"]["search_version"]==expected["search_version"],d["ranker_manifest"])
    check("attrition ledger balances final stage",d["candidate_attrition"]["stages"][-1]["entered"]==d["candidate_attrition"]["stages"][-1]["retained"]+d["candidate_attrition"]["stages"][-1]["rejected"],d["candidate_attrition"]["stages"][-1])
    _aug=next(x for x in d["candidate_attrition"]["stages"] if x.get("stage")=="objective_probe_augmentation")
    check("objective probe augmentation is fully accounted",
          _aug["entered"]==_aug["retained"]+_aug["rejected"] and
          sum(x.get("decision") in {"retained","rejected"} for x in _aug.get("candidate_decisions",[]))==_aug["entered"],_aug)
    _ann=next(x for x in d["candidate_attrition"]["stages"] if x.get("stage")=="full_annotation_diversity_retention")
    check("full annotation beam is fully accounted",
          _ann["entered"]==_ann["retained"]+_ann["rejected"] and len(_ann.get("candidate_decisions",[]))==_ann["entered"],_ann)
    ids=[[c["assay"]["forward"],c["assay"]["reverse"],c["assay"].get("probe")] for c in d["candidates"]]
    check("versioned discrimination ordering is deterministic",ids==expected["discrimination_ids"],(ids,expected["discrimination_ids"]))

if fails:
    print("AUTODESIGN DISCRIMINATION GATE FAILED:",", ".join(fails),flush=True);os._exit(1)
print("ALL AUTODESIGN DISCRIMINATION INVARIANTS PASS",flush=True)
# CPython 3.13 + primer3-py may stall in native interpreter finalization after an
# exhaustive search. This file is already an isolated test process, so terminate
# immediately after flushing the verified result.
os._exit(0)
