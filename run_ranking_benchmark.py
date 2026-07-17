#!/usr/bin/env python3
"""Frozen ranking-truth benchmark with deterministic metrics and ablations.

The adversarial expected preferences were authored before running this script.
They establish whether the ranker obeys declared constraints; they do not prove
that rank 1 will be the best wet-lab assay.  Existing biological/published corpus
checks remain separate because a documented viable assay is not a universal optimum.

The regression gate writes into temporary output directories and reuses the
committed figure bytes.  This keeps ordinary testing non-mutating and byte
deterministic without making Matplotlib a runtime dependency.  The command-line
default refreshes the committed evidence; figure regeneration is a separate
maintenance action with fixed metadata and a stable SVG hash salt.
"""
import argparse, csv, hashlib, json, math, shutil, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oligoforge import ranking as R
from oligoforge import ranking_profiles as RP
from oligoforge import ranking_benchmark as RB

ROOT = Path(__file__).resolve().parent
CORPUS = ROOT / "benchmark" / "ranking_truth_corpus.json"
OUT_JSON = ROOT / "benchmark" / "ranking_truth_results.json"
OUT_CSV = ROOT / "benchmark" / "ranking_truth_results.csv"
OUT_PNG = ROOT / "benchmark" / "ranking_truth_topk.png"
OUT_SVG = ROOT / "benchmark" / "ranking_truth_topk.svg"
OUT_MANIFEST = ROOT / "benchmark" / "ranking_truth_manifest.json"
OUTPUT_NAMES = (
    "ranking_truth_results.json",
    "ranking_truth_results.csv",
    "ranking_truth_topk.png",
    "ranking_truth_topk.svg",
    "ranking_truth_manifest.json",
)


def _clone(value):
    return json.loads(json.dumps(value))


def ordered(case, objective=None, ablate=None, perturbation=None):
    profile = RP.get_profile(objective or case["objective"])
    rows = []
    for c in case["candidates"]:
        e = _clone(c["evidence"])
        if ablate == "specificity":
            e["offtarget"] = {"signal_subjects": 0, "product_subjects": 0}
            if not e["hard_valid"] and any("off-target" in x for x in e["hard_failures"]):
                e["hard_failures"] = [x for x in e["hard_failures"] if "off-target" not in x]
                e["hard_valid"] = not e["hard_failures"]
        elif ablate == "coverage":
            e.update(target_coverage=1.0, worst_isolate_3prime=1.0, probe_mean_identity=1.0)
            if not e["hard_valid"] and any("coverage" in x for x in e["hard_failures"]):
                e["hard_failures"] = [x for x in e["hard_failures"] if "coverage" not in x]
                e["hard_valid"] = not e["hard_failures"]
        elif ablate == "three_prime":
            e["worst_isolate_3prime"] = 1.0
        elif ablate == "probe":
            e["probe_mean_identity"] = 1.0
        elif ablate == "robustness":
            e["condition_robustness"] = {"valid_fraction": 1.0}
        elif ablate == "triplet":
            e["triplet_penalty"] = 0.0; e["worst_dimer"] = -3.0
        elif ablate == "practical":
            e["practical_penalty"] = 0.0; e["degeneracy_fold"] = 1
        elif ablate == "multiplex":
            e["panel_risk"] = 0.0
        elif ablate == "junction":
            e["junction"] = {"level":"strong"}
            if not e["hard_valid"] and any("junction" in x for x in e["hard_failures"]):
                e["hard_failures"] = [x for x in e["hard_failures"] if "junction" not in x]
                e["hard_valid"] = not e["hard_failures"]
        if perturbation:
            perturbation(case["id"], c["id"], e)
        rows.append((R.rank_key(e, profile), c["id"], c["region"], e))
    rows.sort(key=lambda x: (x[0], x[1]))
    return rows


def _mean(values):
    vals=[float(x) for x in values if x is not None and math.isfinite(float(x))]
    return (sum(vals)/len(vals)) if vals else None


def metric(cases, method="new", ablate=None):
    n = len(cases); top = {1:0,3:0,5:0,10:0}; clean = 0; diversity = []
    pair_good=pair_total=0; records=[]; winner_evidence=[]
    for case in cases:
        expected=set(case["expected"])
        if method == "legacy":
            rows=sorted(case["candidates"], key=lambda x:(-x["legacy_score"], x["id"]))
            ids=[x["id"] for x in rows]; top1e=rows[0]["evidence"]; regions=[x["region"] for x in rows[:3]]
            rankpos={x["id"]:i for i,x in enumerate(rows)}
        else:
            rows=ordered(case, ablate=ablate)
            ids=[x[1] for x in rows]; top1e=rows[0][3]; regions=[x[2] for x in rows[:3]]
            rankpos={x[1]:i for i,x in enumerate(rows)}
        for k in top: top[k] += int(bool(expected & set(ids[:k])))
        clean += int(top1e["hard_valid"] and top1e["offtarget"]["signal_subjects"] == 0)
        diversity.append(len(set(regions)))
        for good in expected:
            if good not in rankpos: continue
            for bad in set(ids)-expected:
                pair_total += 1; pair_good += int(rankpos[good] < rankpos[bad])
        winner_evidence.append(top1e)
        records.append(dict(case=case["id"], expected=sorted(expected), top1=ids[0], top3=ids[:3],
                            correct=ids[0] in expected,
                            top1_hard_valid=bool(top1e.get("hard_valid")),
                            top1_signal_offtargets=(top1e.get("offtarget") or {}).get("signal_subjects",0),
                            top1_product_offtargets=(top1e.get("offtarget") or {}).get("product_subjects",0),
                            top1_target_coverage=top1e.get("target_coverage"),
                            top1_worst_3prime=top1e.get("worst_isolate_3prime"),
                            top1_robustness=(top1e.get("condition_robustness") or {}).get("valid_fraction")))
    return dict(n=n, top1=top[1]/n, top3=top[3]/n, top5=top[5]/n, top10=top[10]/n,
                top1_wilson_95=RB.wilson_interval(top[1], n),
                top3_wilson_95=RB.wilson_interval(top[3], n),
                clean_top1=clean/n, clean_top1_wilson_95=RB.wilson_interval(clean, n),
                pairwise_preference_accuracy=(pair_good/pair_total if pair_total else None),
                pairwise_preference_wilson_95=RB.wilson_interval(pair_good, pair_total),
                top1_signal_offtarget_rate=_mean([int((e.get("offtarget") or {}).get("signal_subjects",0)>0) for e in winner_evidence]),
                top1_product_offtarget_rate=_mean([int((e.get("offtarget") or {}).get("product_subjects",0)>0) for e in winner_evidence]),
                serious_3prime_mismatch_rate=_mean([int(float(e.get("worst_isolate_3prime",1.0))<0.8) for e in winner_evidence]),
                mean_top1_target_coverage=_mean([e.get("target_coverage") for e in winner_evidence]),
                mean_top1_probe_coverage=_mean([e.get("probe_mean_identity") for e in winner_evidence]),
                mean_top1_condition_robustness=_mean([(e.get("condition_robustness") or {}).get("valid_fraction") for e in winner_evidence]),
                mean_top3_regions=sum(diversity)/n, records=records)


def _noise(case_id, candidate_id, run, field):
    raw=(case_id+"|"+candidate_id+"|"+str(run)+"|"+field).encode()
    u=int(hashlib.sha256(raw).hexdigest()[:12],16)/float(16**12-1)
    return 2.0*u-1.0


def rank_stability(cases, runs=24):
    rows=[]; stable=0; total=0
    for case in cases:
        baseline=ordered(case)[0][1]; same=0
        for run in range(runs):
            def perturb(cid, cand, e, run=run):
                e["target_coverage"]=min(1.0,max(0.0,float(e.get("target_coverage",0))+0.003*_noise(cid,cand,run,"coverage")))
                e["worst_isolate_3prime"]=min(1.0,max(0.0,float(e.get("worst_isolate_3prime",0))+0.005*_noise(cid,cand,run,"3prime")))
                if e.get("probe_mean_identity") is not None:
                    e["probe_mean_identity"]=min(1.0,max(0.0,float(e["probe_mean_identity"])+0.003*_noise(cid,cand,run,"probe")))
                e["condition_robustness"]["valid_fraction"]=min(1.0,max(0.0,float(e["condition_robustness"].get("valid_fraction",0))+0.03*_noise(cid,cand,run,"robust")))
                e["triplet_penalty"]=max(0.0,float(e.get("triplet_penalty",0))+0.25*_noise(cid,cand,run,"triplet"))
                e["practical_penalty"]=max(0.0,float(e.get("practical_penalty",0))+0.15*_noise(cid,cand,run,"practical"))
            same += int(ordered(case, perturbation=perturb)[0][1] == baseline)
        frac=same/runs; stable += same; total += runs
        rows.append(dict(case=case["id"], baseline_top1=baseline, stable_fraction=frac))
    return dict(runs_per_case=runs, overall_stable_fraction=stable/total if total else None, cases=rows,
                note="small deterministic soft-metric perturbations; hard validity/off-target counts are not perturbed")


def _failure_categories(records):
    out={}
    for r in records:
        if r["correct"]: continue
        if not r["top1_hard_valid"]: cat="hard_invalid_winner"
        elif r["top1_signal_offtargets"]: cat="signal_offtarget_winner"
        elif (r["top1_target_coverage"] or 0)<1: cat="coverage_tradeoff"
        else: cat="soft_or_tie_break_ordering"
        out[cat]=out.get(cat,0)+1
    return out


def _output_paths(output_dir):
    output_dir = Path(output_dir).resolve()
    return {
        "json": output_dir / OUTPUT_NAMES[0],
        "csv": output_dir / OUTPUT_NAMES[1],
        "png": output_dir / OUTPUT_NAMES[2],
        "svg": output_dir / OUTPUT_NAMES[3],
        "manifest": output_dir / OUTPUT_NAMES[4],
    }


def _write_figures(result, png_path, svg_path):
    """Regenerate figures with deterministic metadata under one renderer stack."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels=["Top-1","Top-3","Clean top-1","Pairwise"]
    legacy=[result["all"]["legacy"]["top1"],result["all"]["legacy"]["top3"],result["all"]["legacy"]["clean_top1"],result["all"]["legacy"]["pairwise_preference_accuracy"]]
    new=[result["all"]["new"]["top1"],result["all"]["new"]["top3"],result["all"]["new"]["clean_top1"],result["all"]["new"]["pairwise_preference_accuracy"]]
    with matplotlib.rc_context({"svg.hashsalt": "oligoforge-ranking-truth-v2"}):
        x=range(len(labels)); width=.36
        fig,ax=plt.subplots(figsize=(7.5,4.5)); ax.bar([i-width/2 for i in x],legacy,width,label="Legacy preliminary rank")
        ax.bar([i+width/2 for i in x],new,width,label="Structured ranker")
        ax.set_xticks(list(x),labels); ax.set_ylim(0,1.05); ax.set_ylabel("Fraction / accuracy")
        ax.set_title("OligoForge frozen adversarial ranking benchmark"); ax.legend(); fig.tight_layout()
        fig.savefig(png_path, dpi=180, metadata={"Software":"OligoForge", "Creation Time":None})
        fig.savefig(svg_path, metadata={"Creator":"OligoForge", "Date":None})
        plt.close(fig)


def _stage_figures(paths, regenerate_figures):
    if regenerate_figures:
        return
    for key, canonical in (("png", OUT_PNG), ("svg", OUT_SVG)):
        destination = paths[key]
        if not canonical.is_file():
            raise FileNotFoundError("frozen benchmark figure is missing: %s" % canonical)
        if destination != canonical.resolve():
            shutil.copyfile(canonical, destination)


def main(output_dir=None, regenerate_figures=False):
    output_dir = Path(output_dir or OUT_JSON.parent).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = _output_paths(output_dir)
    corpus=json.load(open(CORPUS))
    corpus_validation = RB.validate_corpus(corpus)
    if not corpus_validation["valid"]:
        print(json.dumps(corpus_validation, indent=2))
        return 1
    byid={x["id"]:x for x in corpus["cases"]}
    t0=time.perf_counter()
    splits={name:[byid[x] for x in ids] for name,ids in corpus["splits"].items()}
    all_cases=corpus["cases"]
    result={"schema":"oligoforge-ranking-results/v2", "corpus":CORPUS.name,
            "ranker_version":RP.RANKER_VERSION, "profile_version":RP.PROFILE_VERSION,
            "warning":"Synthetic/adversarial benchmark; not evidence of universal wet-lab optimality.",
            "corpus_validation": corpus_validation,
            "comparator":"frozen legacy preliminary scalar supplied in each fixture",
            "all":{"legacy":metric(all_cases,"legacy"), "new":metric(all_cases,"new")},
            "splits":{}, "ablations":{}, "rank_stability":rank_stability(all_cases)}
    for name,cases in splits.items():
        result["splits"][name]={"legacy":metric(cases,"legacy"), "new":metric(cases,"new")}
    for feature in ("specificity","coverage","three_prime","probe","robustness","triplet","practical","multiplex","junction"):
        result["ablations"][feature]=metric(all_cases,"new",ablate=feature)
    result["all"]["legacy"]["failure_categories"]=_failure_categories(result["all"]["legacy"]["records"])
    result["all"]["new"]["failure_categories"]=_failure_categories(result["all"]["new"]["records"])
    oc=byid["objective_switch"]
    alt=ordered(oc, objective=oc["alternate_objective"])
    result["objective_switch"]={"objective":oc["alternate_objective"], "top1":alt[0][1],
                                "expected":oc["alternate_expected"], "correct":alt[0][1] in oc["alternate_expected"]}
    elapsed = round(time.perf_counter()-t0, 4)
    paths["json"].write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    with paths["csv"].open("w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh); w.writerow(["case","split","legacy_top1","new_top1","expected","legacy_correct","new_correct",
                                      "new_target_coverage","new_signal_offtargets","new_robustness"])
        split_of={cid:s for s,ids in corpus["splits"].items() for cid in ids}
        old={x["case"]:x for x in result["all"]["legacy"]["records"]}
        new={x["case"]:x for x in result["all"]["new"]["records"]}
        for case in all_cases:
            nr=new[case["id"]]
            w.writerow([case["id"],split_of[case["id"]],old[case["id"]]["top1"],nr["top1"],
                        "|".join(case["expected"]),old[case["id"]]["correct"],nr["correct"],
                        nr["top1_target_coverage"],nr["top1_signal_offtargets"],nr["top1_robustness"]])
    if regenerate_figures:
        _write_figures(result, paths["png"], paths["svg"])
    else:
        _stage_figures(paths, regenerate_figures=False)
    files = [CORPUS, paths["json"], paths["csv"], paths["png"], paths["svg"]]
    trace = ROOT / "benchmark" / "plasmodium_ranking_trace.json.gz"
    if trace.exists(): files.append(trace)
    manifest = {
        "schema":"oligoforge-ranking-benchmark-manifest/v2",
        "frozen":corpus.get("frozen"), "ranker_version":RP.RANKER_VERSION,
        "profile_version":RP.PROFILE_VERSION, "random_seed":None, "deterministic":True,
        "dataset_note":"Synthetic/adversarial preference fixtures plus a biological pipeline trace; not wet-lab optimality labels.",
        "corpus_validation": corpus_validation,
        "splits":corpus.get("splits"), "files":[]}
    for path in files:
        raw=path.read_bytes(); manifest["files"].append({"path":path.name,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()})
    paths["manifest"].write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8", newline="\n")
    new_all=result["all"]["new"]
    held=result["splits"].get("held_out",{}).get("new",{})
    final=result["splits"].get("final_test",{}).get("new",{})
    ok=(new_all["top1"] >= .9 and new_all["clean_top1"] == 1.0 and
        (held.get("top1") or 0)>=.75 and (final.get("top1") or 0)>=.5 and result["objective_switch"]["correct"])
    print(json.dumps({"legacy_top1":result["all"]["legacy"]["top1"],"new_top1":new_all["top1"],
                      "new_pairwise":new_all["pairwise_preference_accuracy"],
                      "new_clean_top1":new_all["clean_top1"],"held_out_top1":held.get("top1"),
                      "final_test_top1":final.get("top1"),"rank_stability":result["rank_stability"]["overall_stable_fraction"],
                      "objective_switch":result["objective_switch"]["correct"],"runtime_seconds_observed":elapsed},indent=2))
    return 0 if ok else 1


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT_JSON.parent)
    parser.add_argument(
        "--regenerate-figures",
        action="store_true",
        help="regenerate PNG/SVG with Matplotlib instead of reusing frozen figure bytes",
    )
    return parser

if __name__ == "__main__":
    args = _parser().parse_args()
    raise SystemExit(main(args.output_dir, args.regenerate_figures))
