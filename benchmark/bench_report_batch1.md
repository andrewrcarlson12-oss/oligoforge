# Track-C benchmark — batch 1 (5 targets)

**Method.** Real RefSeq/GenBank templates (accessions recorded), committed as offline fixtures so
scoring is deterministic and network-free. Each target is designed by OligoForge at its *real*
annealing temperature and scored on a fixed rubric. OligoForge's displayed Tm (`tm_acc`, Biopython
Owczarzy-2008 path) is cross-checked against an **independent** from-scratch NN engine
(`oligoforge.nn.params`, SantaLucia-1998 + Owczarzy-2004/2008) — a genuine two-implementation
agreement test, not the same code twice. Where a target carries a published or golden reference
assay, the harness tests whether OligoForge *admits* those exact primers and, if not, records the
precise gate that declines them.

## Corpus (batch 1)

| id | organism | accession | GC% | profile | Ta | reference |
|----|----------|-----------|-----|---------|----|-----------|
| plas_cytb_ATrich | Plasmodium cytb | fixture (haemosporidian) | 26.8 | parasite_mtdna | 54 | hand (AT-rich anchor) |
| HMBS_host_balanced | A. coerulescens HMBS | XM_068994916.1 | 52.8 | idt_taqman | 60 | golden (host anchor) |
| SDHA_intron_spanning | A. coerulescens SDHA | XM_069004197.1 | 50.0 | idt_taqman | 60 | hand (junction 1707) |
| Mtb_rpoB_GCrich | M. tuberculosis rpoB | NC_000962.3:759807-761000 | 61.9 | gc_rich | 60 | hand (GC-rich stress) |
| human_ACTB_published | H. sapiens ACTB | NM_001101.5 | 55.2 | sybr_generic | 60 | **published (PrimerBank, Spandidos et al. 2010)** |

The corpus spans GC 27→62 %, three chemistries (LNA/SYBR parasite, TaqMan, generic SYBR), and both
intron-spanning and single-exon templates.

## Results

**1. Tm accuracy — PASS across the board.** OligoForge's displayed Tm agrees with the independent
NN engine to **0.008–0.06 °C** on every designed oligo and on both published reference primers,
across the full 27–62 % GC range. The two Tm implementations are separate code; this is real
cross-validation, and it says the thermodynamics are correct end to end.

**2. OligoForge designs a valid assay for all five targets**, including the GC-rich Mtb rpoB
(61.9 % GC) and the AT-rich Plasmodium cytb (26.8 % GC) at their correct anneal temperatures.

**3. SDHA intron-spanning — correct.** Unprompted, OligoForge's SDHA amplicon (mRNA 1636–1766)
**straddles the exon–exon junction at 1707**, so genomic DNA will not co-amplify — the right call
for an intron-aware host-gene assay.

**4. Published PrimerBank ACTB pair — one divergence, classified (b) both valid.** The published
2003 ACTB primers locate exactly in NM_001101.5 (250 bp, as published), and OligoForge computes
their Tm correctly (0.008 °C from the independent engine). OligoForge nonetheless **declines** them:
- forward `CATGTACGTTGCTATCCAGGC` — 3′ GC-clamp too strong (4 of last 5 are G/C, limit 3);
- reverse `CTCCTTAATGTCACGCACGAT` — ends in 3′-T (disallowed by the SYBR profile).

Both are legitimate, current SYBR design guidelines (weak 3′ clamp to limit mispriming; avoid 3′-T).
The published pair predates these conventions. This is **not an engine bug**: the Tm is right, the
primers locate, and OligoForge is applying a defensible stricter rule. No fix warranted; recorded as
a class-(b) divergence.

## Verdict for batch 1

No engine defect surfaced. The Tm math cross-validates to <0.06 °C, designs are produced across the
full GC range, intron-spanning works, and the single reference-primer divergence is a documented
stylistic-rule difference, not a correctness error. Harness is deterministic (byte-identical
scorecard run-to-run) and committed for reuse as the corpus scales.

## Deliverables
- `tests/benchmark/bench_corpus.json` — corpus with recorded accessions/provenance
- `tests/benchmark/bench_score.py` — offline deterministic scorer
- `tests/benchmark/bench_scorecard.csv` — full per-target scorecard
- `tests/test_benchmark.py` — 11 regression asserts pinning these invariants (suite now 25/25)
