"""Track-C benchmark regression: pins the invariants the manual-vs-OligoForge comparison proved.
Offline, deterministic. Run:  PYTHONPATH=. python3 tests/test_benchmark.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.benchmark import bench_score as B
from oligoforge import thermo as T, design as D, profiles as P

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + ("" if cond else f"   [{detail}]"))
    if not cond: fails.append(name)

rows = B.run()
by = {r["id"]: r for r in rows}

# 1. Displayed Tm (Biopython path) agrees with the INDEPENDENT from-scratch NN engine to <0.3 C
#    on every designed oligo across the whole GC range -- the core Tm-accuracy guarantee.
worst = 0.0
for r in rows:
    for k in ("of_fwd_tm_vs_OA", "of_rev_tm_vs_OA", "ref_fwd_tm_vs_OA"):
        v = r.get(k)
        if v is not None: worst = max(worst, v)
check("Tm: displayed vs independent NN engine agree <0.3 C on all oligos", worst < 0.3, f"worst={worst}")

# 2. OligoForge designs a valid assay for every target across GC 27-62%
for tid in ("plas_cytb_ATrich","HMBS_host_balanced","SDHA_intron_spanning","Mtb_rpoB_GCrich","human_ACTB_published"):
    check(f"design produced for {tid}", bool(by[tid].get("of_forward")), by[tid].get("of_error"))

# 3. Published PrimerBank ACTB pair LOCATES and its Tm is computed correctly, even though the
#    stricter modern SYBR profile declines it (documented class-(b) 'both valid' divergence).
actb = by["human_ACTB_published"]
check("ACTB published pair locates in NM_001101.5", actb.get("ref_locates") is True)
check("ACTB published fwd Tm matches independent engine <0.3 C", (actb.get("ref_fwd_tm_vs_OA") or 9) < 0.3)
check("ACTB divergence is a documented 3'-rule decline (not a Tm error)",
      actb.get("ref_fwd_admitted") is False and "3'" in (actb.get("ref_fwd_reject_reason") or ""),
      actb.get("ref_fwd_reject_reason"))

# 4. SDHA design spans the exon-exon junction at mRNA 1707 (gDNA won't co-amplify)
import json
corp = json.load(open(os.path.join(os.path.dirname(B.__file__), "bench_corpus.json")))
sdha = [t for t in corp["targets"] if t["id"]=="SDHA_intron_spanning"][0]["seq"]
T.set_conditions(anneal_c=60)
a = D.design_assay(sdha, P.PROFILES["idt_taqman"])
fpos = sdha.find(a["forward"]); rend = sdha.find(T.revcomp(a["reverse"])) + len(a["reverse"])
check("SDHA amplicon straddles junction 1707", fpos < 1707 < rend, f"[{fpos},{rend}]")
T.set_conditions(anneal_c=60)

# 5. Determinism: two independent design passes give identical picks
a1 = D.design_assay(sdha, P.PROFILES["idt_taqman"]); a2 = D.design_assay(sdha, P.PROFILES["idt_taqman"])
check("design deterministic run-to-run", (a1["forward"],a1["reverse"],a1["amplicon"])==(a2["forward"],a2["reverse"],a2["amplicon"]))

# 6. Structure dG is threaded to primer3 at the anneal temp EXACTLY (no wrapper units/rounding bug).
#    hairpin_full(seq, 60)[1] must equal a direct primer3 calc_hairpin(temp_c=60).dg/1000.
import primer3
T.set_conditions(mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, anneal_c=60)
worst_hp = 0.0
for tid in ("HMBS_host_balanced", "human_GAPDH", "Mtb_rpoB_GCrich"):
    seq = by[tid]["of_forward"]
    of = T.hairpin_full(seq, 60)[1]
    p3 = primer3.calc_hairpin(seq, mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, temp_c=60).dg / 1000.0
    worst_hp = max(worst_hp, abs(of - p3))
check("structure dG threaded to primer3 at anneal temp exactly (<0.01)", worst_hp < 0.01, f"worst={worst_hp:.4f}")

# 7. Probe science: every TaqMan/GC probe is >=3 C above its primers and does NOT start with 5'-G.
bad_probe = []
for r in rows:
    if r.get("of_probe") and r["profile"] in ("idt_taqman", "gc_rich"):
        pr = r["of_probe"]
        ptm = T.tm_acc(pr); ftm = T.tm_acc(r["of_forward"]); rtm = T.tm_acc(r["of_reverse"])
        if pr[0] == "G" or (ptm - max(ftm, rtm)) < 3.0:
            bad_probe.append((r["id"], pr[0], round(ptm - max(ftm, rtm), 1)))
check("all TaqMan/GC probes: 5'-non-G and Tm >=3C above primers", not bad_probe, str(bad_probe[:3]))

# 8. Corpus is broad enough to be meaningful: >=15 targets spanning >=25 GC-percentage-points.
gcs = [float(r["gc"]) for r in rows]
check("corpus >=15 targets", len(rows) >= 15, len(rows))
check("corpus spans >=25 GC points", (max(gcs) - min(gcs)) >= 25, f"{min(gcs)}-{max(gcs)}")
T.set_conditions(anneal_c=60)

if fails:
    print("BENCHMARK FAILURES:", fails); sys.exit(1)
print("ALL BENCHMARK ASSERTS PASS")
