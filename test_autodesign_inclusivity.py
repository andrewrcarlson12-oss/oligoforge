"""Independent broad-inclusivity ranking regression for the frozen isolate corpus."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from oligoforge import autodesign as AD, profiles as P

HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixtures")
fails=[]
def check(name,cond,detail=""):
    print(("  PASS " if cond else "  FAIL ")+name+("" if cond else "   [%s]" % (detail,)),flush=True)
    if not cond:fails.append(name)

plas=json.load(open(os.path.join(FX,"plasmodium_cytb.json")))["sequences"]
expected=json.load(open(os.path.join(FX,"autodesign_expected_v1.33.json")))
g=AD.design_from_sequences(plas,P.PROFILES["parasite_mtdna"],n_candidates=5,objective="broad_inclusivity")
check("genus design returns",not g.get("error") and g.get("candidates"),g.get("error"))
if g.get("candidates"):
    e=g["candidates"][0]["evidence"]
    check("genus rank 1 is hard-valid",e["hard_valid"],e["hard_failures"])
    check("genus rank 1 coherent coverage >=95%",e["target_coverage"]>=.95,e["target_coverage"])
    check("effective order pool is explicit",set(e["effective_oligos"])>={"forward","reverse","probe","uses_degeneracy"},e["effective_oligos"])
    check("finalists are not exact duplicates",len({(c["assay"]["forward"],c["assay"]["reverse"],c["assay"].get("probe")) for c in g["candidates"]})==len(g["candidates"]))
    _aug=next(x for x in g["candidate_attrition"]["stages"] if x.get("stage")=="objective_probe_augmentation")
    check("inclusivity probe augmentation is fully accounted",
          _aug["entered"]==_aug["retained"]+_aug["rejected"] and
          sum(x.get("decision") in {"retained","rejected"} for x in _aug.get("candidate_decisions",[]))==_aug["entered"],_aug)
    ids=[[c["assay"]["forward"],c["assay"]["reverse"],c["assay"].get("probe")] for c in g["candidates"]]
    check("versioned inclusivity ordering is deterministic",ids==expected["broad_inclusivity_ids"],(ids,expected["broad_inclusivity_ids"]))
if fails:
    print("AUTODESIGN INCLUSIVITY GATE FAILED:",", ".join(fails),flush=True);os._exit(1)
print("ALL AUTODESIGN INCLUSIVITY INVARIANTS PASS",flush=True)
os._exit(0)
