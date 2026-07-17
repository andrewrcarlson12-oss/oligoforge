"""
G6b + G3 — Primer-BLAST specificity head-to-head on a real genome (pins bench_specificity_*).

The head-to-head oracle is a faithful offline Primer-BLAST-equivalent built on blastn 2.17.0 (the
same engine Primer-BLAST wraps). blastn is not part of this offline suite, so this test:

  1. RE-RUNS the OligoForge offline in-silico-PCR deterministically against the committed genome
     FASTA fixtures (human transcriptome with real paralogs; apicomplexan 18S) and pins the subject
     calls per pair.
  2. PINS the committed head-to-head result (bench_specificity_realgenome.json), where OligoForge
     and the Primer-BLAST-equivalent made IDENTICAL calls on all 11 human-transcriptome pairs plus
     the pan-Plasmodium case (100% concordance), both scoring sensitivity 100% / specificity 98.5%
     vs a biologically-grounded ground truth and SHARING the same 2 conservative false positives.
  3. Pins the real project result: pan-Plasmodium 18S primers amplify all 3 Plasmodium spp. and do
     NOT cross-react with Haemoproteus / Leucocytozoon / Toxoplasma (genus-specific).

Fully offline (reads committed FASTA + JSON; no blastn, no network). Standalone script.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import specificity as SP

fails = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("" if ok else f"  [{detail}]"))
    if not ok: fails.append(name)

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "benchmark")
FIX = os.path.join(BENCH, "genome_fixtures")
human = open(os.path.join(FIX, "human_transcriptome.fasta")).read()
api = open(os.path.join(FIX, "apicomplexa_18S.fasta")).read()
res = json.load(open(os.path.join(BENCH, "bench_specificity_realgenome.json")))
recs = json.load(open(os.path.join(BENCH, "bench_specificity_records.json")))

# ---- 1. deterministic OligoForge in-silico-PCR on the specific geNorm pairs ----
genorm = {
    "ACTB":   ("CTGGAACGGTGAAGGTGACA","AAGGGACTTCCTGTAACAATGCA"),
    "GAPDH":  ("TGCACCACCAACTGCTTAGC","GGCATGGACTGTGGTCATGAG"),
    "YWHAZ":  ("ACTTTTGGTACATTGTGGCTTCAA","CCGCCAGGACAAACCAGTAT"),
}
for gene, (f, r) in genorm.items():
    of = SP.in_silico_pcr_offline(f, r, human, max_mm=2)
    subj = sorted(set(p["subject"] for p in of.get("products", [])))
    check(f"{gene} amplifies ONLY {gene} (specific, despite paralogs in DB)",
          subj == [gene], f"got {subj}")

# ---- 2. deterministic cross-reactivity detection on a designed paralog pair ----
# ACTB/ACTG1 shared-CDS pair must amplify BOTH (real cross-reaction OligoForge must catch)
F, R = "CAACGGCTCCGGCATGTGCA", "CCTGGGGCGCCCCACGATGG"
of = SP.in_silico_pcr_offline(F, R, human, max_mm=2)
subj = sorted(set(p["subject"] for p in of.get("products", [])))
check("designed ACTB/ACTG1 pair CATCHES the cross-reaction (both subjects)",
      subj == ["ACTB", "ACTG1"], f"got {subj}")

# ---- 3. pin committed head-to-head concordance + confusion ----
ht = res["human_transcriptome"]
check("OligoForge vs Primer-BLAST: 100% concordance on human pairs",
      ht["tool_concordance"] == "11/11", ht["tool_concordance"])
check("OligoForge sensitivity 100% vs ground truth", ht["oligoforge"]["sensitivity"] == 1.0)
check("OligoForge specificity 98.5% vs ground truth (2 shared conservative FP)",
      abs(ht["oligoforge"]["specificity"] - 0.985) < 0.005, str(ht["oligoforge"]["specificity"]))
check("OligoForge and Primer-BLAST share the SAME accuracy (identical confusion)",
      ht["oligoforge"] == ht["primerblast"])
check("the 2 false positives are the borderline paralogs (ACTG1, YWHAB)",
      sorted(x[1] for x in ht["shared_false_positives"]) == ["ACTG1", "YWHAB"])

# ---- 4. pin the real pan-Plasmodium project result (deterministic OligoForge run) ----
corpus = json.load(open(os.path.join(BENCH, "bench_corpus_published.json")))["assays"]
plas = next(a for a in corpus if "Plasmodium" in a["id"])
kF, kR, kP = "GCTCTTTCTTGATTTCTTGGATG", plas["reverse"], plas.get("probe", "")
of = SP.in_silico_pcr_offline(kF, kR, api, probe=kP or None, max_mm=2)
amp = sorted(set(p["subject"] for p in of.get("products", [])))
plasmodium = [s for s in amp if s.startswith(("Pfalciparum","Pvivax","Pberghei"))]
nontarget = [s for s in amp if s.startswith(("Hcatharti","Leucocytozoon","Toxoplasma"))]
check("pan-Plasmodium primers amplify all 3 Plasmodium spp.", len(plasmodium) == 3, str(amp))
check("pan-Plasmodium primers do NOT cross-react with related genera (Haemoproteus/Leuco/Toxo)",
      len(nontarget) == 0, f"cross-reacts: {nontarget}")

print(f"\n{'ALL PRIMER-BLAST SPECIFICITY ASSERTS PASS' if not fails else str(len(fails))+' FAILED'}")
sys.exit(len(fails))
