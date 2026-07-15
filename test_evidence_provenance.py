#!/usr/bin/env python3
"""Evidence uncertainty, provenance, feedback, and benchmark-integrity gates."""
import csv, io, json, sys
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oligoforge import __version__ as APP_VERSION
from oligoforge import ranking as R
from oligoforge import ranking_explain as RX
from oligoforge import ranking_profiles as RP
from oligoforge import experimental_feedback as EF
from oligoforge import ranking_benchmark as RB
from oligoforge import provenance as PROV
from oligoforge import report as REPORT
from oligoforge import rdml as RDML
from oligoforge import thermo as T

passed=failed=0

def check(name, ok, detail=None):
    global passed, failed
    if ok:
        passed+=1; print("PASS",name)
    else:
        failed+=1; print("FAIL",name,detail if detail is not None else "")


def ev(*, valid=True, cov=1.0, robust=1.0, trip=1.0, practical=0.0,
       off_eval=True, signal=0, product=0):
    return dict(hard_valid=valid, hard_failures=[] if valid else ["bad geometry"],
                objective="balanced", target_coverage=cov, worst_isolate_3prime=1.0,
                probe_mean_identity=1.0,
                offtarget=dict(signal_subjects=signal, product_subjects=product),
                condition_robustness=dict(valid_fraction=robust), triplet_penalty=trip,
                worst_dimer=-4.0, practical_penalty=practical, panel_risk=0,
                degeneracy_fold=1, junction=None,
                evaluations=dict(target_epcr=True, offtarget_epcr=off_eval,
                                 condition_robustness=True, panel=False, junction=False))

# Evidence-completeness and explicit indifference.
a={"evidence":ev(off_eval=False,trip=1.0),"assay":{}}
b={"evidence":ev(off_eval=False,trip=1.2),"assay":{}}
p=RX.preference_strength(a,b)
check("missing exclusivity evidence prevents false moderate preference",
      p["state"]=="insufficient evidence to distinguish" and "offtarget_epcr" in p["missing_evidence"],p)

c={"evidence":ev(off_eval=True,trip=1.0),"assay":{}}
d={"evidence":ev(off_eval=True,trip=1.2),"assay":{}}
p=RX.preference_strength(c,d)
check("complete close candidates are near-equivalent",p["state"]=="near-equivalent alternatives",p)

p=RX.preference_strength(c,{"evidence":ev(valid=False),"assay":{}})
check("hard validity remains a strong preference",p["state"]=="strong preference",p)

# Rank bands retain deterministic order while admitting modeled indifference.
prof={"amp_min":70,"amp_max":160,"len_min":18,"len_max":28,"gc_min":20,"gc_max":80,
      "tm_min":40,"tm_max":75,"tm_opt":60,"pair_tm_gap_max":8,"no_three_prime_T":False,
      "max_3prime_gc":5,"max_g_run":5,"max_any_run":6,"hairpin_min":-20,"self_dimer_min":-20,
      "pair_dimer_min":-20,"probe_len_min":10,"probe_len_max":40,"probe_offset_min":-20,
      "probe_offset_max":30,"probe_hairpin_min":-20,"anneal_c":60,"no_probe":False}
# Use prebuilt evidence to avoid native calculations: reproduce the rank-band loop logic through rank candidates is
# tested elsewhere; here verify explanation carries an injected group.
c["equivalence_group"]={"group_id":1,"ranks":[1,2],"size":2}
x=RX.explain(c,d)
check("explanation carries equivalence group",x["equivalence_group"]["ranks"]==[1,2],x)
check("explanation reports computational scope",x["evidence_completeness"]["scope"].startswith("computational"),x)

# Deterministic provenance is self-hashing and complete.
obj=RP.get_profile("balanced")
m1=R.manifest(obj,{"full_annotation_pool":25},{"corpus":"abc"},
              external_databases={"ncbi":{"accessions":["X1"]}},
              constraints={"workflow":"test"})
m2=R.manifest(obj,{"full_annotation_pool":25},{"corpus":"abc"},
              external_databases={"ncbi":{"accessions":["X1"]}},
              constraints={"workflow":"test"})
check("provenance manifest deterministic",m1==m2,(m1,m2))
check("manifest has stable run id and self hash",m1["run_id"].startswith("ofrun_") and len(m1["manifest_sha256"])==64,m1)
check("manifest records native scientific versions",bool(m1["software_versions"].get("primer3_py")) and bool(m1["software_versions"].get("viennarna")),m1)
check("manifest records models and conditions",m1["scientific_models"]["reaction_condition_snapshot"]["anneal_c"] is not None,m1)
check("manifest records external database state",m1["external_database_state"]=="declared",m1)
check("manifest verifier accepts untampered manifest",PROV.verify_manifest(m1)["valid"],PROV.verify_manifest(m1))
tampered=dict(m1); tampered["objective"]="tampered"
check("manifest verifier rejects tampering",not PROV.verify_manifest(tampered)["valid"],PROV.verify_manifest(tampered))


# Workbench exports preserve and verify the ranker's chain of custody.
export_assay={
    "name":"=unsafe name", "gene":"G", "organism":"O",
    "forward":"ACGTACGTACGTACGTACGT", "reverse":"TGCATGCATGCATGCATGCA",
    "probe":"CGTACGTACGTACGTACGTA", "amplicon":95,
    "ranker_manifest":m1, "objective_profile":obj, "candidate_rank":1,
    "rank_explanation":{
        "preference_state":"near-equivalent alternatives",
        "evidence_completeness":{"state":"complete for declared computational objective"}},
    "source_workflow":"automatic_design",
}
rep=REPORT.build([export_assay],{"title":"provenance test"})
check("HTML report preserves verified run provenance",
      m1["run_id"] in rep["html"] and m1["manifest_sha256"] in rep["html"] and "manifest verified" in rep["html"],rep["html"][-1500:])
rows=list(csv.DictReader(io.StringIO(rep["csv"])))
check("CSV report preserves ranker manifest fields",
      len(rows)==1 and rows[0]["design_run_id"]==m1["run_id"] and rows[0]["manifest_state"]=="verified" and rows[0]["ranker_version"]==RP.RANKER_VERSION,rows)
check("CSV report still neutralizes formula injection",rows[0]["name"].startswith("'="),rows[0]["name"])
original_snap=T._snapshot()
try:
    T.set_conditions(mv_conc=80,dv_conc=6,dntp_conc=1.2,dna_conc=500,anneal_c=64)
    rep_conditions=REPORT.build([export_assay])
    condition_row=list(csv.DictReader(io.StringIO(rep_conditions["csv"])))[0]
    ms=m1["scientific_models"]["reaction_condition_snapshot"]
    recorded_snap=(ms["mv_conc_mM"],ms["dv_conc_mM"],ms["dntp_conc_mM"],ms["total_oligo_conc_nM"],ms["anneal_c"])
    expected_tm=round(T._tm_acc_at(export_assay["forward"],recorded_snap),1)
    check("report recomputes QC under recorded assay conditions, not mutable session state",
          abs(float(condition_row["fwd_Tm"])-expected_tm)<1e-9 and "recorded rank-manifest conditions" in condition_row["thermodynamic_conditions"],condition_row)
finally:
    T.set_conditions(mv_conc=original_snap[0],dv_conc=original_snap[1],dntp_conc=original_snap[2],dna_conc=original_snap[3],anneal_c=original_snap[4])
rd=RDML.build([export_assay],{"title":"provenance_test"})
root=ET.fromstring(rd["xml"]); ns="{http://www.rdml.org}"
desc=root.find(ns+"target/"+ns+"description")
check("RDML description preserves verified rank provenance",
      desc is not None and m1["run_id"] in desc.text and m1["manifest_sha256"] in desc.text and "manifest verified" in desc.text,desc.text if desc is not None else None)
tampered_export=dict(export_assay); tampered_export["ranker_manifest"]=dict(m1,objective="tampered")
rep_bad=REPORT.build([tampered_export])
check("report surfaces altered manifests instead of trusting them","invalid / altered" in rep_bad["html"],rep_bad["html"][-1000:])

# Feedback validation, dedupe, conflict audit, and leakage-safe group splits.
base={"assay_id":"A","target_group":"G1","design_run_id":"R1","ranker_version":RP.RANKER_VERSION,
      "objective":"balanced","status":"success","efficiency":95,"r2":0.99}
ds=EF.dataset_status([base,base])
check("feedback exact duplicates removed",ds["n_unique"]==1 and ds["n_duplicates"]==1,ds)
conf=dict(base); conf["status"]="failed"
ds=EF.dataset_status([base,conf])
check("feedback conflicts surfaced",ds["n_conflicts"]==1,ds)
try:
    EF.normalize(dict(base,r2=1.2)); bad=False
except ValueError: bad=True
check("feedback rejects impossible R squared",bad)
try:
    EF.normalize(dict(base,efficiency=-1)); bad=False
except ValueError: bad=True
check("feedback rejects negative efficiency",bad)
parsed=EF.parse_records('assay_id,target_group,status,ranker_version,objective\nA,G1,success,2.2.0,balanced\n','csv')
check("feedback CSV import parses records",len(parsed)==1 and parsed[0]["assay_id"]=="A",parsed)
rows=[]
for i in range(20):
    rows.append(dict(base,assay_id="A%d"%i,target_group="G%d"%(i//2),design_run_id="R%d"%i,
                     status="success" if i%2==0 else "failed"))
split=EF.target_group_split(rows)
by={}
for row in split["records"]:
    by.setdefault(row["target_group"],set()).add(row["dataset_split"])
check("target-group split prevents leakage",not split["group_leakage"] and all(len(v)==1 for v in by.values()),split)
status=EF.calibration_status([base])
check("small feedback remains non-authoritative",not status["learned_reranker_allowed"] and status["unmet_requirements"],status)

# Benchmark schema validator catches leakage and malformed expectations.
corpus={"cases":[
    {"id":"a","class":"experimental","target_group":"G","expected":["x"],"candidates":[{"id":"x"}]},
    {"id":"b","class":"experimental","target_group":"G","expected":["y"],"candidates":[{"id":"y"}]},
],"splits":{"development":["a"],"held_out":["b"]}}
v=RB.validate_corpus(corpus)
check("benchmark validator detects target-group leakage",not v["valid"] and v["target_group_leakage"],v)
corpus["cases"][1]["target_group"]="H"
v=RB.validate_corpus(corpus)
check("benchmark validator accepts isolated groups",v["valid"] and not v["target_group_leakage"],v)
corpus["cases"][1]["expected"]=["missing"]
v=RB.validate_corpus(corpus)
check("benchmark validator rejects unknown accepted candidate",not v["valid"],v)

# Public API surfaces use the same validators and leakage-safe splitter.
from fastapi.testclient import TestClient
from app import app
client=TestClient(app)
payload=json.dumps([base,base])
r=client.post("/api/experimental-feedback/import",json={"payload":payload,"format_hint":"json"})
body=r.json()
check("feedback import API normalizes and deduplicates",r.status_code==200 and body.get("n_unique")==1 and body.get("n_duplicates")==1,body)
r=client.post("/api/experimental-feedback/import",json={"payload":json.dumps([dict(base,r2=2.0)]),"format_hint":"json"})
body=r.json()
check("feedback import API reports rejected scientific values",r.status_code==200 and body.get("n_rejected")==1 and "r2" in body.get("rejected",[{}])[0].get("error",""),body)
r=client.post("/api/experimental-feedback/split",json={"records":rows,"train_fraction":0.70,"validation_fraction":0.15})
body=r.json()
check("feedback split API preserves target-group isolation",r.status_code==200 and not body.get("group_leakage") and len(body.get("records") or [])==len(rows),body)
r=client.post("/api/experimental-feedback/split",json={"records":rows,"train_fraction":0.90,"validation_fraction":0.15})
check("feedback split API rejects impossible fractions",r.status_code==422 and "split fractions" in r.json().get("error",""),r.json())

# Frozen biological trace carries current provenance and a verified ordering.
trace=json.loads((Path(__file__).resolve().parent/"benchmark"/"plasmodium_ranking_trace.json").read_text())
tm=trace.get("ranker_manifest") or {}
check("biological trace uses current schema",trace.get("schema")=="oligoforge-biological-ranking-trace/v3",trace.get("schema"))
check("biological trace ordering matches frozen winner fixture",trace.get("expected_ordering_match") is True,trace.get("expected_ordering_match"))
check("biological trace carries current ranker provenance",tm.get("application_version")==APP_VERSION and tm.get("ranker_version")==RP.RANKER_VERSION and tm.get("scoring_profile_version")==RP.PROFILE_VERSION,tm)
check("biological trace manifest is self-hashing",PROV.verify_manifest(tm)["valid"],PROV.verify_manifest(tm))

print("EVIDENCE_PROVENANCE: %d passed / %d failed / %d total"%(passed,failed,passed+failed))
raise SystemExit(1 if failed else 0)
