"""Regression gates for the v1.33 ranking-truth and manual-rescue release."""
import json, os, sys, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from oligoforge import candidate_retention as CR
from oligoforge import ranking as R
from oligoforge import ranking_profiles as RP
from oligoforge import design as D
from oligoforge import manual_design as MD
from oligoforge import assay_rescue as AR
from oligoforge import experimental_feedback as EF
from oligoforge import profiles as P
from oligoforge import thermo as T

N = F = 0
def check(name, ok, detail=None):
    global N, F
    N += 1
    if ok:
        print("PASS", name)
    else:
        F += 1
        print("FAIL", name, detail)


def ev(valid=True, cov=1.0, signal=0, products=0, robust=1.0, triplet=0.0,
       practical=0.0, panel=0.0, fold=1, f3=1.0):
    return dict(hard_valid=valid, hard_failures=[] if valid else ["hard failure"],
                target_coverage=cov, worst_isolate_3prime=f3,
                offtarget=dict(signal_subjects=signal, product_subjects=products),
                condition_robustness=dict(valid_fraction=robust), triplet_penalty=triplet,
                practical_penalty=practical, panel_risk=panel, degeneracy_fold=fold,
                probe_mean_identity=1.0, worst_dimer=-4.0, junction=None)

# Attrition and diversity retention.
def assay(i, reg=0, score=None, f=None, r=None, p=None):
    start = reg * 250 + i
    return dict(forward=f or ("ACGTACGTACGTACG%02d" % i),
                reverse=r or ("TGCATGCATGCATGC%02d" % i), probe=p or ("GATTACAGATTACA%02d" % i),
                f_xy=[start, start + 18], r_xy=[start + 80, start + 100],
                probe_xy=[start + 35, start + 55], amplicon=100,
                candidate_rank=float(i if score is None else score), search_window_start=start)
rows = [assay(i, reg=0) for i in range(12)] + [assay(20+i, reg=1) for i in range(4)]
rows += [dict(rows[0]), dict(rows[0])]
kept, led = CR.retain_diverse(rows, limit=8, per_region=3, per_near=1)
check("attrition ledger balances", led["entered"] == led["retained"] + led["rejected"], led)
check("exact duplicates counted", led["reasons"]["exact_duplicate"] >= 2, led)
check("regional diversity preserved", led["regions_retained"] >= 2, led)
check("candidate budget enforced", 2 <= len(kept) <= 8, len(kept))
check("triplet decision trace accounts for every input", len(led.get("candidate_decisions", [])) == len(rows), led)
# A dense early region must not consume a tight global budget before later regions
# receive one representative.
_dense = [assay(i, reg=0, score=i/100.0) for i in range(20)] + [assay(80, reg=1, score=5), assay(90, reg=2, score=6)]
_dkeep, _dled = CR.retain_diverse(_dense, limit=3, region_size=250, per_region=16, per_near=3)
check("tight triplet beam gives each available region a representative",
      _dled["regions_retained"] == 3, (_dled, [x["f_xy"] for x in _dkeep]))

# Pair and probe beams preserve trade-off diversity before expensive annotation.
pairs = []
for reg in range(3):
    for i in range(5):
        pairs.append(dict(f="F%d%d"%(reg,i), r="R%d%d"%(reg,i), fstart=reg*180+i,
                          fend=reg*180+i+20, rstart=reg*180+70+i, rend=reg*180+95+i,
                          amp=95+(i%3)*25, score=float(i)+reg*.1, gap=1.0, dimer=-4.0))
pkeep, pled = CR.retain_pairs_diverse(pairs, limit=6, region_size=150, amplicon_bin=25)
check("pair beam preserves target regions", pled["regions_retained"] >= 3, pled)
check("pair beam decision trace balances", len(pled["candidate_decisions"]) == len(pairs) and pled["entered"] == pled["retained"]+pled["rejected"], pled)
probes = [dict(probe="ACGTACGTACGT%02d"%i, strand=("+" if i%2==0 else "-"),
               start=20+i*5, end=32+i*5, offset=6+(i%4), preliminary_penalty=float(i%3)) for i in range(10)]
prkeep, prled = CR.retain_probes_diverse(probes, limit=4, position_bin=8)
check("probe beam preserves positional alternatives", prled["positions_retained"] >= 3, prled)
check("probe beam preserves both strands when available", len(prled["strands_retained"]) == 2, prled)
# A one-slot beam must choose the best probe overall, not the first strand label.
_one_slot = [dict(probe="AAAAAAAAAAAAAAAAAAAA", strand="+", start=10, end=30, offset=8,
                  preliminary_penalty=4.0),
             dict(probe="CCCCCCCCCCCCCCCCCCCC", strand="-", start=12, end=32, offset=8,
                  preliminary_penalty=0.1)]
_one_keep, _one_led = CR.retain_probes_diverse(_one_slot, limit=1)
check("one-slot probe beam selects strongest overall candidate",
      _one_keep[0]["probe"] == "CCCCCCCCCCCCCCCCCCCC", (_one_keep, _one_led))

# Hard-invalid candidates never share the first Pareto front with a dominating valid candidate.
items = [dict(evidence=ev(False, cov=1.0), assay=assay(40), rank=2),
         dict(evidence=ev(True, cov=1.0), assay=assay(41), rank=1)]
R.pareto_fronts(items)
check("hard validity dominates Pareto tier", items[1]["evidence"]["pareto_front"] == 1 and items[0]["evidence"]["pareto_front"] > 1, items)

# Lexicographic hard constraints and objective-specific reversals.
bal = RP.get_profile("balanced")
check("hard failure outranks no soft advantage", R.rank_key(ev(True, triplet=50), bal) < R.rank_key(ev(False, triplet=0), bal))
check("specificity reverses preliminary preference", R.rank_key(ev(True, signal=0, triplet=8), bal) < R.rank_key(ev(False, signal=1, triplet=0), bal))
inc = RP.get_profile("broad_inclusivity")
check("coverage is primary in inclusivity profile", R.rank_key(ev(True, cov=1.0, triplet=8), inc) < R.rank_key(ev(True, cov=.96, triplet=0), inc))

# Finalist selection is deterministic/idempotent and never fills with invalid assays.
valid_rows = []
for i in range(4):
    valid_rows.append(dict(evidence=ev(True, cov=1-.01*i, triplet=i), assay=assay(50+i, reg=i), rank=i+1))
valid_rows.append(dict(evidence=ev(False), assay=assay(60), rank=5))
a = R.select_finalists(valid_rows, 5)
b = R.select_finalists(valid_rows, 5)
check("invalid finalists excluded", all(x["evidence"]["hard_valid"] for x in a), a)
check("finalist categories do not accumulate", [x["finalist_categories"] for x in a] == [x["finalist_categories"] for x in b])
# Category labels must remain attached to the actual category winner instead of
# being reassigned to the next unused assay merely to fill the display list.
_top=min(valid_rows[:-1],key=lambda x:x["rank"])
check("category labels describe the true winner",
      {"recommended_balanced","best_specificity","best_inclusivity","most_condition_robust","minimal_degeneracy","best_multiplex_fit"}.issubset(set(_top["finalist_categories"])),
      _top["finalist_categories"])
check("finalists prefer distinct primer pairs when available",
      len({(x["assay"]["forward"],x["assay"]["reverse"]) for x in a}) == len(a),
      [(x["assay"]["forward"],x["assay"]["reverse"],x["finalist_categories"]) for x in a])

# Joint triplet search retains several probe alternatives for at least one pair.
HMBS=("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCC"
      "ATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
prof = P.PROFILES["idt_taqman"]
triplets, tled = D.generate_assay_candidates(HMBS, prof, pair_limit=8, probes_per_pair=4, triplet_limit=32)
pairs = {}
for x in triplets:
    pairs.setdefault((x["forward"], x["reverse"]), set()).add(x.get("probe"))
check("joint search constructs multiple probes per primer pair", any(len(v) > 1 for v in pairs.values()), {k:len(v) for k,v in pairs.items()})
check("joint-search ledger declares bounded status", tled.get("search_status") == "heuristic_bounded", tled)

# Manual mapping reports all placements, including reverse strand.
rep = "ACGTACGTACGTACGTACGTGGGG" * 2
hits = MD.map_oligo("ACGTACGTACGTACGTACGT", rep, "forward", max_mm=0)
check("manual mapping preserves multiple exact placements", len(hits) >= 2, hits)
rhits = MD.map_oligo(T.revcomp("ACGTACGTACGTACGTACGT"), rep, "reverse", max_mm=0)
check("manual mapping includes reverse-strand placements", any(h["strand"] == "-" for h in rhits), rhits)

# Manual and automatic analysis use the same authoritative Tm functions.
base = D.design_assay(HMBS, prof)
man = MD.analyze_assay(base["forward"], base["reverse"], HMBS, prof, base["probe"])
ma = man["candidate"]["assay"]
check("manual-auto forward Tm identity", abs(ma["f_tm"] - T.tm(base["forward"])) < 1e-9, ma)
check("manual-auto reverse Tm identity", abs(ma["r_tm"] - T.tm(base["reverse"])) < 1e-9, ma)
_rob = man["candidate"]["evidence"]["condition_robustness"]
check("condition robustness re-evaluates structure in every scenario",
      len(_rob.get("scenarios", [])) == 3 and all("structure_at_anneal" in x and "failure_reasons" in x for x in _rob["scenarios"]), _rob)

# A primer-only degenerate alternative must retain the unchanged probe.
_iupac = {frozenset("AG"):"R",frozenset("CT"):"Y",frozenset("CG"):"S",frozenset("AT"):"W",
          frozenset("GT"):"K",frozenset("AC"):"M"}
mut = list(HMBS); fdeg = list(base["forward"])
for j in range(3):
    pos = base["f_xy"][0] + j
    old = mut[pos]; new = next(x for x in "ACGT" if x != old)
    mut[pos] = new; fdeg[j] = _iupac.get(frozenset((old,new)), "N")
deg_assay = dict(base); deg_assay["forward_deg"] = "".join(fdeg)
deg_sc = dict(score=0, score_raw=0, assay=deg_assay, conservation={})
deg_ev = R.build_evidence(deg_sc, [HMBS,"".join(mut)], [], prof, RP.get_profile("broad_inclusivity"))
check("primer-only degeneracy retains unchanged probe", deg_ev["effective_oligos"]["probe"] == base["probe"], deg_ev["effective_oligos"])

# True locked-component redesigns.
rr = MD.constrained_redesign(base["forward"], base["reverse"], HMBS, prof, base["probe"],
                             locks={"primer_pair": True}, max_results=2)
check("locked primer pair is never replaced", bool(rr["candidates"]) and all(x["assay"]["forward"] == base["forward"] and x["assay"]["reverse"] == base["reverse"] for x in rr["candidates"]), rr.get("note"))
rrp = MD.constrained_redesign(base["forward"], base["reverse"], HMBS, prof, base["probe"],
                              locks={"probe": True}, max_results=2)
check("locked probe is never replaced", bool(rrp["candidates"]) and all(x["assay"]["probe"] == base["probe"] for x in rrp["candidates"]), rrp.get("note"))

# Rescue diagnosis separates model evidence from causal claims.
synthetic = {"candidate": {"assay": {"pair_tm_gap": 4.5, "probe": "A", "amplicon": 220},
                           "evidence": ev(False, cov=.8, signal=1, robust=.33)},
             "mappings": {"forward": [{"mismatches": 1}], "reverse": [{"mismatches": 0}]}}
diag = AR._diagnoses(synthetic, {"efficiency": 80, "r2": .95})
check("rescue identifies multiple defensible hypotheses", len(diag) >= 5, diag)
check("rescue labels experimental inference separately", all("computational_evidence" in x and "experimental_inference" in x for x in diag))
rescued = AR.rescue(base["forward"], base["reverse"], HMBS, prof, base["probe"],
                    observed={"efficiency": 80}, max_results=1, max_runtime_s=20)
check("rescue completes inside declared interactive budget",
      rescued["runtime"]["elapsed_s"] <= rescued["runtime"]["budget_s"], rescued["runtime"])
check("rescue records bounded search provenance",
      bool(rescued.get("search_ledgers")) and all((x.get("ledger") or {}).get("limits", {}).get("full_annotation_limit") == 12
                                                  for x in rescued["search_ledgers"]),
      rescued.get("search_ledgers"))
# A one-base mismatch to a concrete intended-template site is repaired only after
# the edited complete assay outranks the original under the same evidence model.
_badf = list(base["forward"]); _badf[5] = next(x for x in "ACGT" if x != _badf[5]); _badf = "".join(_badf)
_badbase = MD.analyze_assay(_badf, base["reverse"], HMBS, prof, base["probe"])
_repairs = AR._one_base_repairs(_badf, base["reverse"], base["probe"], HMBS, prof,
                                _badbase, "balanced", None, None)
check("one-base rescue restores the concrete intended-template primer",
      any(x["candidate"]["assay"]["forward"] == base["forward"] and
          x["components_changed"] == ["forward"] for x in _repairs), _repairs)
check("rescue disruption order is nondecreasing",
      [x.get("disruption_order", 99) for x in rescued.get("redesigns", [])] ==
      sorted(x.get("disruption_order", 99) for x in rescued.get("redesigns", [])),
      [x.get("disruption_order") for x in rescued.get("redesigns", [])])

# Feedback schema and no-ML evidence gate.
fb = EF.calibration_status([{"assay_id":"A", "target_group":"G", "status":"failed", "efficiency":80}])
check("small feedback set cannot activate learned reranker", not fb["learned_reranker_allowed"], fb)
rec = EF.normalize({"assay_id":"A", "target_group":"G", "status":"failed", "notes":"retained failure"})
check("feedback record is reproducibly hashed", len(rec["record_sha256"]) == 64, rec)

# Manifest stability and version provenance.
m1 = R.manifest(bal, {"pool": 40}, {"template":"abc"})
m2 = R.manifest(bal, {"pool": 40}, {"template":"abc"})
check("rank manifest deterministic", m1 == m2 and m1["ranker_version"] == RP.RANKER_VERSION, m1)
ranked_trace, _ = R.rank_candidates([
    dict(score=0, score_raw=0, assay=assay(70), conservation={}),
    dict(score=0, score_raw=0, assay=assay(71), conservation={})],
    [], [], P.PROFILES["idt_taqman"], objective_name="balanced")
check("ranking trace reconstructs named hierarchy", bool(ranked_trace[0].get("rank_trace")) and
      ranked_trace[0]["rank_trace"]["priority_order"] == list(bal["priority"]), ranked_trace[0].get("rank_trace"))

# API contracts: invalid ranking objectives fail visibly instead of returning a success envelope.
try:
    from fastapi.testclient import TestClient
    import app as APP
    client = TestClient(APP.app)
    bad = client.post("/api/manual-design/analyze", json={"forward":base["forward"],"reverse":base["reverse"],
        "probe":base["probe"],"template":HMBS,"objective":"not_a_real_objective"})
    check("manual API rejects unknown objective", bad.status_code == 422, (bad.status_code,bad.text))
    bad2 = client.post("/api/manual-design/redesign", json={"forward":base["forward"],"reverse":base["reverse"],
        "probe":base["probe"],"template":HMBS,"objective":"not_a_real_objective"})
    check("redesign API rejects unknown objective", bad2.status_code == 422, (bad2.status_code,bad2.text))
except Exception as exc:
    check("ranking-truth API contract", False, str(exc))

print("RANKING_TRUTH: %d passed / %d failed / %d total" % (N-F, F, N))
sys.exit(1 if F else 0)
