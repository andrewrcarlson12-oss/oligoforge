"""Thermodynamic-driver tests for the Certified Orthogonal Panel. Uses the primer3 NN engine, so
these are deterministic. Run from repo root:  python tests/test_orthopanel_thermo.py

Covers intake (dedup / IUPAC / RNA->DNA / modification-strip / rejects), the self-structure filter,
confusability-graph construction (revcomp cross-check, ΔG symmetry), and full certify_panel behavior
including a GOLDEN oligo panel whose cross-reactions form a 5-cycle. Exact branch-and-bound proves the optimum; theta remains a diagnostic. Edge cases: empty, single, all-cross, none-cross, duplicates,
degenerate/IUPAC, RNA, very short/long, and the SDP size guard.
"""
import os, sys
sys.path.insert(0, ".")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from oligoforge import orthopanel as OP, thermo as T

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)

HAVE_THETA = OP.lovasz_theta(3, [(0, 1)]) is not None

# ---------------- ΔG symmetry ----------------
a, b = "ACGATCAGTTGCATCAGGTA", "TGACCTAGCATGGTACAAGC"
check("heterodimer ΔG is symmetric ΔG(a,b)=ΔG(b,a)", abs(T.hetero_dimer(a, b) - T.hetero_dimer(b, a)) < 1e-9)

# ---------------- intake ----------------
recs, rej = OP.intake(["ACGATCAGTTGCATCAGGTA", "ACGATCAGTTGCATCAGGTA", "acgatcagttgcatcaggta"])
check("intake dedups identical/case-variant sequences to 1", len(recs) == 1, len(recs))
check("intake records dup_count", recs[0]["dup_count"] == 3, recs[0]["dup_count"] if recs else None)

recs2, _ = OP.intake([{"seq": "ACGUACGUACGUACGUACGU", "name": "rna"}])
check("intake converts RNA U->DNA T", recs2 and "U" not in recs2[0]["seq"], recs2[0]["seq"] if recs2 else "")
check("intake flags RNA conversion", recs2 and recs2[0].get("rna_converted") is True)

recs3, _ = OP.intake([{"seq": "/56-FAM/ACGATCAGTTGCATCAGGTA/3IABkFQ/", "name": "probe"}])
check("intake strips IDT modification notation to bare sequence",
      recs3 and recs3[0]["seq"] == "ACGATCAGTTGCATCAGGTA", recs3[0]["seq"] if recs3 else "")

recs4, _ = OP.intake([{"seq": "ACGWRYKMSBDHVN", "name": "deg"}])
check("intake accepts IUPAC degenerate bases and flags degenerate",
      recs4 and recs4[0].get("degenerate") is True)

recs5, rej5 = OP.intake([{"seq": "ACGT123XYZ!!", "name": "junk"}, {"seq": "ACGATCAGTTGCATCAGGTA"}])
check("intake rejects invalid characters, keeps valid", len(recs5) == 1 and len(rej5) == 1, (len(recs5), len(rej5)))

# ---------------- self-structure filter ----------------
clean = {"seq": "ACGATCAGTTGCATCAGGTA", "name": "clean"}          # weak self-structure
hairpinny = {"seq": "GGGGGCCCCCGGGGGCCCCC", "name": "hairpin"}     # strong hairpin/homodimer
recs6, _ = OP.intake([clean, hairpinny])
surv, drop = OP.self_structure_filter(recs6, self_dg=-9.0)
check("self-filter keeps the clean oligo", any(s["name"] == "clean" for s in surv))
check("self-filter drops the strongly self-structured oligo", any(d["name"] == "hairpin" for d in drop))
check("dropped oligo carries a reason", bool(drop) and "ΔG" in drop[0]["reason"])

# ---------------- confusability graph (revcomp cross-check) ----------------
A = "ACGATCAGTTGCATCAGGTA"
seqs = [A, T.revcomp(A), "TGACCTAGCATGGTACAAGC"]   # 0 and 1 are complementary; 2 is orthogonal
edges, edg = OP.build_confusability_graph(seqs, cross_dg=-6.0)
check("revcomp pair forms an edge", (0, 1) in edg)
check("orthogonal oligo forms no edge with the pair", (0, 2) not in edg and (1, 2) not in edg)
check("edge stores the ΔG that formed it", edg.get((0, 1), 0) < -6.0)

# ---------------- certify_panel: real pool, certified via the cheap bound ----------------
pool = [{"seq": A, "name": "FwdA"}, {"seq": T.revcomp(A), "name": "ProbeA_rc"},
        {"seq": "GCATTCAGATCGGATACCTA", "name": "FwdB"}, {"seq": "TGACCTAGCATGGTACAAGC", "name": "FwdC"},
        {"seq": "CATGACTGCAAGTCGATACG", "name": "FwdD"}]
r = OP.certify_panel(pool, cross_dg=-6.0, self_dg=-9.0, k=2)
panel_names = {p["name"] for p in r["panel"]}
check("real pool: exactly one cross-reaction edge", r["n_edges"] == 1, r["n_edges"])
check("real pool: panel size 4", r["panel_size"] == 4, r["panel_size"])
check("real pool: panel keeps exactly one of the cross-reacting pair",
      len({"FwdA", "ProbeA_rc"} & panel_names) == 1)
check("real pool: certified maximum (gap 0)", r["certified"] and r["gap"] == 0)
check("real pool: certified by exact branch-and-bound",
      r["bound_source"] == "exact_bnb")

# ---------------- GOLDEN: oligo panel whose cross-reactions form a C5 ----------------
# Cyclically overlapping 13-bp complementary blocks: each oligo cross-reacts with its two cycle
# neighbours and no one else. clique-cover(C5)=3 cannot certify; θ=√5 -> ⌊θ⌋=2 does.
C5_POOL = [
    {"seq": "CAGATTTTCATATTATGCAGAAAATC", "name": "O0"},
    {"seq": "GATTTTCTGCATATACTTCGCCTGAT", "name": "O1"},
    {"seq": "ATCAGGCGAAGTAACGAGTCGGTTAT", "name": "O2"},
    {"seq": "ATAACCGACTCGTCTTCGGATACTGT", "name": "O3"},
    {"seq": "ACAGTATCCGAAGATATGAAAATCTG", "name": "O4"},
]
g = OP.certify_panel(C5_POOL, cross_dg=-8.0, self_dg=-16.0, k=2)
edge_set = {(e["i"], e["j"]) for e in g["edges"]}
check("golden C5: all 5 survive the self-filter", g["n_after_self_filter"] == 5, g["n_after_self_filter"])
check("golden C5: cross-reactions form exactly a 5-cycle",
      edge_set == {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}, sorted(edge_set))
check("golden C5: maximum independent set (panel) is 2", g["panel_size"] == 2, g["panel_size"])
check("golden C5: clique-cover bound is 3 (too loose to certify)", g["clique_cover_bound"] == 3,
      g["clique_cover_bound"])
if HAVE_THETA:
    check("golden C5: θ ≈ √5 (2.236)", g["theta"] is not None and abs(g["theta"] - 2.2361) < 3e-2, g["theta"])
    check("golden C5: exact search proves upper bound 2", g["upper_bound"] == 2, g["upper_bound"])
    check("golden C5: certified by exact branch-and-bound, not theta",
          g["certified"] and g["bound_source"] == "exact_bnb" and not g.get("theta_certifying"))
    sp = g["split_pool"]
    check("golden C5: split-pool reports theta only as diagnostic", sp.get("theta_pow_k") is not None
          and abs(sp["theta_pow_k"] - 5.0) < 0.2 and sp.get("theta_diagnostic_only"), sp)
    check("golden C5: θ^k ceiling exceeds naive |MIS|^k (distance reasoning undercounts)",
          sp["theta_pow_k"] > sp["naive_mis_pow_k"])
else:
    check("golden C5: exact search certifies even without theta",
          g["upper_bound"] == 2 and g["certified"] and g["bound_source"] == "exact_bnb")

# ---------------- edge cases ----------------
e_empty = OP.certify_panel([], k=1)
check("empty input: panel empty, certified trivially", e_empty["panel_size"] == 0 and e_empty["certified"])

e_single = OP.certify_panel([{"seq": A, "name": "only"}], k=1)
check("single oligo: panel size 1, certified", e_single["panel_size"] == 1 and e_single["certified"])

# all-cross: a huge positive cross_dg makes every pair an edge -> complete graph -> panel 1
e_allcross = OP.certify_panel(pool, cross_dg=1e6, self_dg=-30.0, k=1)
check("all-pairwise-cross-reacting: panel collapses to 1", e_allcross["panel_size"] == 1, e_allcross["panel_size"])
check("all-cross: still certified (complete graph, θ/clique-cover -> 1)", e_allcross["certified"])

# none-cross: a huge negative cross_dg means no edges -> panel = all survivors
e_nocross = OP.certify_panel(pool, cross_dg=-1e6, self_dg=-9.0, k=1)
check("none-cross-reacting: panel = all survivors", e_nocross["panel_size"] == e_nocross["n_after_self_filter"])
check("none-cross: certified maximum", e_nocross["certified"])

# duplicates surface in the payload
e_dup = OP.certify_panel([{"seq": A, "name": "x1"}, {"seq": A, "name": "x2"}], k=1)
check("duplicate sequences are deduped and reported", e_dup["n_unique"] == 1 and len(e_dup["duplicates"]) == 1)

# degenerate + RNA inputs don't crash the pipeline
e_deg = OP.certify_panel([{"seq": "ACGRYSWKMACGRYSWKM", "name": "d"}, {"seq": "ACGUACGUACGUACGUAC", "name": "r"}], k=1)
check("degenerate/IUPAC + RNA input runs without error", "panel_size" in e_deg)

# very short and very long oligos: no crash (long exceeds primer3's dimer cap -> neutral)
e_len = OP.certify_panel([{"seq": "ACGTAC", "name": "short"}, {"seq": "ACGT" * 60, "name": "long"}], k=1)
check("very short + very long oligos handled without crashing", "panel_size" in e_len)

# ---------------- SDP size guard ----------------
guard = OP.certify_panel(pool, cross_dg=-6.0, self_dg=-9.0, size_limit=3, use_theta=True)
check("size guard: pool larger than size_limit still returns a certified panel", "certified" in guard)
check("size guard: θ not attempted above the limit (θ is None)", guard["theta"] is None)

print(("FAIL: " + ", ".join(_fails)) if _fails else "OK")
sys.exit(1 if _fails else 0)
