"""
Head-to-head validation: OligoForge vs Primer3 on the PUBLISHED-primer corpus.

Pins the differentiator numbers reported in bench_headtohead_report.md so they are
regression-protected:
  - displayed Tm agrees with an INDEPENDENT NN implementation to <0.1 C (internal validation)
  - displayed Tm runs warmer than Primer3 (divalent-aware salt), a real physical gap
  - hairpin dG at 37 C is byte-identical to Primer3 (same backend)
  - structure is less stable at the true annealing temperature
  - the anneal-temp gate admits materially more candidates on a GC-rich template than a 37 C gate,
    and ~0 more on an AT-rich template (the effect scales with GC as physics requires)

Offline except for the corpus JSON (no network). Standalone script (sys.exit), matching the suite.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import thermo as T, nn as NN, profiles as P, design as D
import primer3

fails = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("" if ok else f"  [{detail}]"))
    if not ok: fails.append(name)

HERE = os.path.dirname(os.path.abspath(__file__))
pub = json.load(open(os.path.join(HERE, "benchmark", "bench_corpus_published.json")))["assays"]
tc  = json.load(open(os.path.join(HERE, "benchmark", "bench_corpus.json")))["targets"]

# qPCR master-mix salt for a fair comparison
T.set_conditions(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=200.0, anneal_c=60)

def nondeg(seq): return all(b in "ACGT" for b in seq.upper())

# ---- 1. displayed Tm vs independent NN engine ----
worst_nn = 0.0
for a in pub:
    for seq in [a["forward"], a["reverse"]] + ([a["probe"]] if a["probe"] else []):
        if not nondeg(seq): continue
        d = abs(T.tm_acc(seq) - NN.params(seq)["tm"])
        worst_nn = max(worst_nn, d)
check("displayed Tm agrees with independent NN engine to <0.1 C", worst_nn < 0.1, f"worst={worst_nn:.3f}")

# ---- 2. displayed Tm runs warmer than Primer3 at matched qPCR salt (divalent-aware) ----
diffs = []
for a in pub:
    for seq in [a["forward"], a["reverse"]] + ([a["probe"]] if a["probe"] else []):
        if not nondeg(seq): continue
        p3 = primer3.calc_tm(seq, mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=200.0)
        diffs.append(T.tm_acc(seq) - p3)
mean_gap = sum(diffs) / len(diffs)
check("displayed Tm runs warmer than Primer3 at qPCR salt (divalent-aware)", 0.2 < mean_gap < 1.5, f"mean_gap={mean_gap:.2f} C over {len(diffs)} oligos")

# ---- 3. hairpin dG at 37 C is identical to Primer3 (same backend) ----
worst_37 = 0.0
for a in pub:
    for seq in [a["forward"], a["reverse"]] + ([a["probe"]] if a["probe"] else []):
        if not nondeg(seq): continue
        of = T.hairpin_full(seq, 37)[0]
        p3 = primer3.calc_hairpin(seq, mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, temp_c=37).dg / 1000.0
        worst_37 = max(worst_37, abs(of - p3))
check("hairpin dG at 37 C identical to Primer3 backend (<0.01)", worst_37 < 0.01, f"worst={worst_37:.4f}")

# ---- 4. structure is (on average) less stable at the true annealing temp ----
shifts = []
for a in pub:
    ann = a["anneal_c"]
    for seq in [a["forward"], a["reverse"]] + ([a["probe"]] if a["probe"] else []):
        if not nondeg(seq): continue
        dg37, dg_ann, _ = T.hairpin_full(seq, ann)
        shifts.append(dg_ann - dg37)
mean_shift = sum(shifts) / len(shifts)
check("structure less stable at true Ta (mean dG shift > 0)", mean_shift > 0.0, f"mean_shift={mean_shift:.3f}")

# ---- 5. anneal-temp gate admits more candidates than a 37 C gate, scaling with GC ----
def admits(template, prof, anneal):
    f, r = D.enumerate_primers(template, {**prof, "anneal_c": anneal})
    return len(f) + len(r)

def seq_of(tid): return [t for t in tc if t["id"] == tid][0]["seq"]

gc_seq = seq_of("Mtb_rpoB_GCrich"); prof_gc = P.PROFILES["idt_taqman"]
at_seq = seq_of("plas_cytb_ATrich"); prof_at = P.PROFILES["parasite_mtdna"]
gc_gain = admits(gc_seq, prof_gc, prof_gc.get("anneal_c", 60)) - admits(gc_seq, prof_gc, 37)
at_gain = admits(at_seq, prof_at, prof_at.get("anneal_c", 54)) - admits(at_seq, prof_at, 37)
check("anneal-temp gate recovers many candidates on a GC-rich template", gc_gain > 300, f"gc_gain={gc_gain}")
check("anneal-temp gate effect ~0 on an AT-rich template (scales with GC)", at_gain == 0, f"at_gain={at_gain}")

T.set_conditions(anneal_c=60)
if fails:
    print("HEAD-TO-HEAD FAILURES:", fails); sys.exit(1)
print("ALL HEAD-TO-HEAD ASSERTS PASS")
