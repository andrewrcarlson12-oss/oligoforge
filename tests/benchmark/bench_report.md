# Track-C benchmark — full corpus report

## Corpus
**18 real templates**, all fetched fresh from NCBI RefSeq/GenBank with accessions recorded,
then committed as offline fixtures (`bench_corpus.json`) so scoring is deterministic and network-free.

- **GC coverage:** 26.8% – 61.9%
- **Organism classes:** apicomplexan mtDNA/rRNA (Plasmodium, Toxoplasma), vertebrate housekeeping
  mRNAs (8 human ref-genes), bacterial (M. tuberculosis rpoB, E. coli 16S), viral (SARS-CoV-2 N),
  plant (Arabidopsis ACT2).
- **Structure/chemistry:** intron-spanning (SDHA, junction 1707) vs single-exon; TaqMan, generic SYBR,
  GC-rich, AT-rich parasite/LNA.
- **Reference types:** 1 golden (HMBS host anchor), 1 published (PrimerBank ACTB, independently
  confirmed + verified to locate in NM_001101.5), 16 hand/profile-rubric.

## Method
Each target is designed by OligoForge at its real annealing temperature. The rubric scores Tm-in-profile,
GC, amplicon size, 3′ clamp, run limits, and hairpin/self-dimer ΔG **at the anneal temperature**. The
displayed Tm (`tm_acc`, Biopython Owczarzy-2008 path) is cross-checked against an **independent
from-scratch NN engine** (`oligoforge.nn.params`) — genuine two-implementation agreement, not the same
code twice. Reference primers are tested for admission, with the exact declining gate recorded.

## Headline results

**1. Tm accuracy — PASS, decisively.** Across **37 designed + reference oligos**, OligoForge's
displayed Tm agrees with the independent NN engine to **max 0.060 °C, mean 0.032 °C,
p95 0.056 °C** — over the full 27–62 % GC range and every organism class.
Two separate thermodynamic implementations converging this tightly is strong evidence the Tm math is correct.

**2. Design production — 18/18.** OligoForge produced a valid assay for every template, AT-rich to GC-rich.

**3. Structure ΔG threading — exact.** OligoForge's `hairpin_full(seq, 60 °C)` reproduces a direct
primer3 `calc_hairpin(temp_c=60)` to **0.0000 kcal/mol**. The v1.27.0/1.27.1 anneal-temperature gate
passes the temperature through correctly; there is no wrapper rounding/units bug.

**4. Probe science — correct on all 16 TaqMan/GC targets.** Every probe sits **+6 to +9 °C above its
primers** (the TaqMan ΔTm rule) and **none starts with 5′-G** (fluorophore-quenching rule). Enforced
automatically.

**5. Intron awareness — correct.** The SDHA amplicon straddles the exon–exon junction at mRNA 1707
unprompted (gDNA won't co-amplify).

**6. Degenerate handling — correct.** For a degenerate primer, QC surfaces the full Tm span
(e.g. `W` → 58.0–59.2 °C over 2 resolutions; 2×`Y` → 55.1–58.5 °C over 4) with an explicit note. The
single scalar `tm()` is the internal ranking value; the user-facing number is the range. No under-reporting.

**7. LNA — reasonable.** McTigue-informed increments (~+3.5 °C/LNA) with honest ±2 °C bands; the
independent `nn.params_lna` reproduces the locked Plas probe's LNA-effective Tm (66.9 °C) matching the
QC endpoint.

## Divergences (classified)

Exactly **one** reference-primer divergence across the corpus:

- **human_ACTB_published (class b — both valid).** The published PrimerBank ACTB pair locates exactly
  in NM_001101.5 (250 bp) and OligoForge computes its Tm correctly (0.008 °C from the independent
  engine), but OligoForge **declines** it: forward has a 3′ GC-clamp of 4/5 (profile limit 3), reverse
  ends in 3′-T (disallowed). Both are defensible modern SYBR rules (weak 3′ clamp, avoid 3′-T) that the
  older published pair predates. **Not an engine bug** — the Tm is right and the primers locate;
  OligoForge is applying a stricter, legitimate guideline. No fix.

No class-(a) divergences (OligoForge scientifically wrong) were found.

## Verdict
Across 18 real templates spanning the full GC range and six organism classes, **no engine defect
surfaced.** The Tm math cross-validates to <0.06 °C against an independent implementation; structure
ΔG is exact against primer3 at the anneal temperature; probe, intron, degenerate, and LNA handling are
all scientifically correct. The single divergence is a documented stylistic-rule difference, not an
error. Harness is deterministic and committed for reuse.

## Scorecard

| target | organism | GC% | profile | Ta | OF amplicon | fwd Tm |Δ|vs NN | rev Tm |Δ| | design? |
|--------|----------|-----|---------|----|-----------|--------------|-----------|---------|
| plas_cytb_ATrich | — | 26.8 | parasite_mtdna | 54 | 122 | 0.06 | 0.014 | yes |
| HMBS_host_balanced | — | 52.8 | idt_taqman | 60 | 135 | 0.028 | 0.044 | yes |
| SDHA_intron_spanning | — | 50.0 | idt_taqman | 60 | 130 | 0.027 | 0.029 | yes |
| Mtb_rpoB_GCrich | — | 61.9 | gc_rich | 60 | 148 | 0.034 | 0.031 | yes |
| human_ACTB_published | — | 55.2 | sybr_generic | 60 | 109 | 0.032 | 0.039 | yes |
| human_GAPDH | — | 56.1 | idt_taqman | 60 | 108 | 0.056 | 0.028 | yes |
| human_B2M | — | 39.1 | idt_taqman | 60 | 106 | 0.029 | 0.016 | yes |
| human_HPRT1 | — | 40.4 | idt_taqman | 60 | 96 | 0.031 | 0.028 | yes |
| human_TBP | — | 46.6 | idt_taqman | 60 | 105 | 0.055 | 0.042 | yes |
| human_GUSB | — | 55.5 | idt_taqman | 60 | 103 | 0.035 | 0.022 | yes |
| human_RPLP0 | — | 53.0 | idt_taqman | 60 | 103 | 0.026 | 0.02 | yes |
| human_PPIA | — | 46.0 | idt_taqman | 60 | 84 | 0.029 | 0.023 | yes |
| Tgondii_B1 | — | 48.1 | idt_taqman | 60 | 126 | 0.024 | 0.036 | yes |
| Pfal_18S | — | 35.6 | idt_taqman | 60 | 139 | 0.031 | 0.025 | yes |
| Athaliana_ACT2 | — | 42.5 | idt_taqman | 60 | 116 | 0.017 | 0.05 | yes |
| SARSCoV2_Ngene | — | 47.2 | idt_taqman | 60 | 98 | 0.045 | 0.047 | yes |
| human_GC_GNAS | — | 56.2 | idt_taqman | 60 | 107 | 0.037 | 0.026 | yes |
| Ecoli_16S | — | 54.4 | idt_taqman | 60 | 116 | 0.05 | 0.003 | yes |

## Deliverables
- `tests/benchmark/bench_corpus.json` — 18-target corpus, accessions + provenance
- `tests/benchmark/bench_score.py` — offline deterministic scorer (independent-Tm cross-check + reject-reason)
- `tests/benchmark/bench_scorecard.csv` — full per-target scorecard
- `tests/benchmark/bench_report.md` — this report
- `tests/test_benchmark.py` — regression asserts (suite 25/25)
