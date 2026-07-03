# OligoForge — Paper-Readiness Assessment

**Version:** 1.28.0 · **Date:** 2026-07-02 · **Scope:** honest evaluation of what is
publishable now, what a reviewer will push back on, and what is missing.

> This document is deliberately critical. Every quantitative claim below is traceable to a
> committed artifact (benchmark CSV/JSON + test). Where OligoForge does **not** beat the
> incumbent (Primer3), that is stated plainly — an honest benchmark that shows a loss on some
> axis is more useful than one engineered to win.

---

## 1. What OligoForge is (and is not)

OligoForge is a **local-first, qPCR-specialized primer/probe design and QC application**: a
FastAPI backend + a single-file vanilla-JS UI wrapping a 23-module engine. It is built on
Florida Scrub-Jay RefSeq data for a real 5-target Th1/Th2 immune-polarization panel.

It is **not** a new thermodynamic method. Primer3 (`primer3-py`) is the secondary-structure
engine underneath, and the nearest-neighbor Tm sum is SantaLucia-1998, the same core Primer3
uses. OligoForge's contribution is a **qPCR-correct application layer** on top of that engine:
a divalent-aware display Tm, an annealing-temperature structure gate, offline in-silico PCR,
LNA/degenerate awareness, and MIQE-compliant reporting. The publication case rests on that
layer being *correct* and *useful*, not on a novel algorithm.

---

## 2. Named differentiators (real, quantified)

Each is backed by the head-to-head benchmark (`bench_headtohead_*.{csv,json}`, `test_headtohead.py`)
over an 11-assay published-primer corpus (25 non-degenerate oligos), or by the validation suites.

### D1 — qPCR-specific salt model (divalent-aware displayed Tm)
The **displayed/reported** Tm uses Owczarzy-2008 magnesium correction with free-Mg²⁺ computed as
[Mg²⁺]−[dNTP] (von Ahsen 2001), the quantity that actually sets duplex stability in a master mix.
- vs an **independent** NN reimplementation (`nn.params`, a separate code path from the display
  `tm_acc`): **mean 0.030 °C, max 0.060 °C** — genuine two-implementation agreement, not the
  same code twice.
- vs Primer3 at **matched qPCR salt**: **+0.44 °C mean** (Primer3's SantaLucia salt path is
  near-Mg-insensitive; the divalent term is the difference).
- vs Primer3 at **its defaults** (dv=1.5, dntp=0.6): **+3.97 °C mean** — this is the *naive-user*
  gap, i.e. what a user gets if they never change Primer3's salt defaults. It is **not** a
  method-vs-method advantage and must not be presented as one.

**Honest magnitude:** the fair, matched-salt advantage is **modest (~0.4 °C)**. The headline
+4 °C only exists against unconfigured Primer3.

### D2 — annealing-temperature structure gate
Dimer/hairpin ΔG reject-gates are evaluated at the **real annealing temperature** (60 °C host /
54 °C parasite), not Primer3's default 37 °C. Structure is less stable at Ta, so the 37 °C gate
is over-conservative and discards viable primers.
- Structure ΔG at 37 °C matches Primer3 **exactly (0.0000 kcal/mol)** — same backend, so the
  anneal-temperature evaluation is the *only* difference (no reimplementation drift).
- ΔG(Ta) − ΔG(37) = **+0.25 kcal/mol mean, up to +1.60** across published oligos.
- **Candidate-admission impact** scales with template GC exactly as physics requires:
  Mtb rpoB (61.9 % GC) **+49.3 %**, HMBS (52.8 %) **+32.3 %**, GNAS (56.2 %) **+30.8 %**,
  Plasmodium cytb (26.8 % AT-rich) **+0.0 %**.

**Honest caveat:** this is measured at *candidate-screening* time — more candidates survive the
gate. It does **not** by itself prove the *final selected* primers are better; a larger admitted
pool can also contain marginal candidates. The claim is "physically-correct gating recovers
viable GC-rich candidates a 37 °C gate wrongly rejects," not "produces better assays." 0/25
published oligos crossed the structure floor at either temperature (they are already clean), so
the value is at design time, not at QC of already-good primers.

### D3 — offline in-silico-PCR specificity (new in 1.28.0)
Deterministic primer-site scan of a supplied FASTA: both strands, IUPAC-aware, ≤ `max_mm`
internal mismatches with a **required exact 3′-terminal anchor** (the physical requirement for
polymerase extension), verified on both `+` and `−` orientations. No network, no BLAST install.
Emits hits in the same shape the existing `epcr()`/`assay_verdict()` consume, so amplicon +
probe-binding logic is shared.
- Controlled validation (`test_specificity_offline.py`, 9 checks): **sensitivity 100 %,
  specificity 100 % on 14 controls** — 8 specific assays clear, a near-identical pseudogene / a
  2-mismatch paralog / an off-size misprime are all flagged, and a 3′-terminal-mismatch variant
  is correctly rejected on both strands.

### D4 — LNA support (McTigue 2004 nearest-neighbor)
Full 32-term McTigue increment set on the SantaLucia core.
- vs McTigue's experimental duplex Tm's: **RMSE 1.86 °C, MAE 1.62 °C** (n = 12). **8/12 within
  ±2 °C, 12/12 within ±3 °C**; worst single duplex 3.0 °C (an outlier *beyond* ±2 °C, flagged in
  the figure). Small negative bias (−0.52 °C).

### D5 — degenerate-primer awareness end to end
`tm_range` surfaces the true per-resolution Tm spread (validated exact vs brute-force enumeration);
genus autodesign collapses aligned variants into IUPAC codes with a min-count/min-frequency guard
against sequencing-noise false degeneracy; the real Plasmodium-cytb genus workflow reproduces its
golden (F `TTTCTACATTTACAAGGTAGCA`, amp 86, n_degenerate 4).

### D6 — engineering quality
Local-first (no cloud, no data egress — a genuine advantage for clinical/field use), single-file
no-build UI, MIQE-structured reporting, per-request-safe conditions (immutable-snapshot cache
key; 0 torn/stale reads across 20 000 read cycles), input-validation hardening (length caps,
non-finite salt rejection), and a **38-module test suite** (30 Python + 8 Node) with byte-locked
golden designs.

---

## 3. Figures & tables inventory (paper-ready)

| # | Artifact | Type | Shows |
|---|----------|------|-------|
| F1 | `tm_comparison.png` | figure (3-panel) | NN agreement (y=x), Bland-Altman vs Primer3, gate impact by GC |
| F2 | `structure_comparison.png` | figure (2-panel) | exact 37 °C agreement; candidate-dimer flip at true Ta |
| F3 | `specificity_validation.png` | figure (2-panel) | confusion matrix (100/100); per-control-mode resolution |
| F4 | `lna_degenerate_validation.png` | figure (3-panel) | LNA vs McTigue; per-duplex error; degenerate Tm spread |
| T1 | `bench_headtohead_tm.csv` | table | per-oligo Tm: OligoForge / independent NN / Primer3 (qPCR + default salt) |
| T2 | `bench_headtohead_structure.csv` | table | per-oligo hairpin ΔG at 37 °C and Ta vs Primer3 |
| T3 | `bench_gate_impact.csv` | table | candidate admission at 37 °C vs Ta, by template GC |
| T4 | `bench_corpus_published.json` | data | 11 published assays, dual-verified (locate + DOI in OpenAlex) |
| T5 | `corpus_provenance.csv` | table | citation + accession provenance for every corpus entry |
| T6 | `specificity_validation.json` | data | 14-control confusion matrix |
| T7 | `lna_validation.json` | data | McTigue per-duplex errors, RMSE/MAE/bias |

Supporting: `bench_headtohead_report.md`, `bench_report.md`. Tests pinning every finding:
`test_headtohead.py` (6 checks), `test_specificity_offline.py` (9), `test_lna_degenerate.py` (13),
`test_benchmark.py` (design-side), `test_nn.py` (LNA + NN core).

---

## 4. Honest gaps — what a reviewer will push back on

**G1 — No wet-lab validation.** Everything is in-silico. The 5-target panel is being ordered but
results are not in. For BMC Bioinformatics / Bioinformatics Advances this is the single biggest
limitation: the qPCR-salt and anneal-gate advantages are *thermodynamically motivated and
self-consistent*, but not shown to improve real assay performance (Cq, efficiency, specificity on
a bench). **Mitigation for now:** frame claims as "thermodynamically correct" not "empirically
better"; add wet-lab data before a methods-journal submission.

**G2 — OligoForge is a layer on Primer3, not a new engine.** The 0.03 °C NN agreement partly
reflects a shared SantaLucia backend. A reviewer will (correctly) note the novelty is the qPCR
application layer, not the thermodynamics. **Mitigation:** claim the layer, cite Primer3 as the
engine, do not overstate.

**G3 — Specificity validation is a constructed control set, not a genome-wide benchmark.**
100 %/100 % is on **14 designed controls** (positive controls are engineered pseudogene /
paralog / off-size constructs). It demonstrates the *mechanism* is correct; it is **not** a claim
of genome-wide off-target accuracy, and it has not been compared against Primer-BLAST or a BLAST
ground truth on a large real off-target set. **Mitigation:** benchmark against Primer-BLAST calls
on a real genome; report ROC on a larger natural set.

**G4 — The salt advantage is modest at matched salt.** The fair comparison is ~0.4 °C; only the
naive-default comparison reaches ~4 °C. Do not lead with the +4 °C number.

**G5 — Small, human-skewed corpus.** 11 published assays, mostly human housekeeping genes +
SARS-CoV-2 + Plasmodium. Dropped entries (UBC repeat junction, RdRp pan-sarbecovirus degenerate)
are documented, but the corpus is not a broad cross-species or cross-chemistry benchmark.

**G6 — Only Primer3 as a comparator.** No head-to-head vs Primer-BLAST, PrimerQuest, Beacon
Designer, or IDT OligoAnalyzer (the tool users actually cross-check against). A reviewer will ask
why. **Mitigation:** at minimum add OligoAnalyzer Tm spot-checks and a Primer-BLAST specificity
comparison.

**G7 — LNA sample is small and the two LNA layers differ in rigor.** 12 duplexes; the
base-averaged `tm_lna` (thermo.py) is cruder than the NN `params_lna` (nn.py). MELTING/rmelting
(the reference implementation) could not be run here (needs Java) — validation is against the
embedded McTigue experimental subset, not a live MELTING cross-check. **Mitigation:** run
`rmelting` where Java is available; expand the LNA panel.

**G8 — "More candidates admitted" ≠ "better assay" (see D2 caveat).** The gate result is a
design-space claim, not an assay-quality claim.

**G9 — Performance not benchmarked.** Exhaustive candidate enumeration in Python is slower than
Primer3's C core; no runtime numbers are reported. Fine for interactive single-assay design;
untested at batch/genome scale.

---

## 5. Venue fit (honest)

- **JOSS** — **strong fit now.** JOSS evaluates software quality, tests, docs, and utility, not
  novel biology. OligoForge has a real user, a 38-module test suite with byte-locked goldens,
  reproducible benchmarks, and a clear scope. This is the most defensible near-term target. Needs:
  a short paper, install/usage docs, contribution guidelines.
- **Bioinformatics Advances / BMC Bioinformatics (Software article)** — **plausible with work.**
  Needs G1 (wet-lab or at least a real-assay retrospective), G3 + G6 (comparison vs Primer-BLAST
  and one qPCR-specific tool), and framing that foregrounds the qPCR application layer. The
  differentiators are real but individually modest; the case is "correct, integrated, local-first
  qPCR workflow," not "new algorithm."
- **A methods/thermodynamics journal** — **not without wet-lab data and a novel method.** The
  thermodynamics is applied, not new.

---

## 6. Bottom line

OligoForge is **publishable as a software tool (JOSS) now**, and **can reach a bioinformatics
software-article bar** with a Primer-BLAST/OligoAnalyzer comparison and, ideally, wet-lab data
from the panel already being ordered. The engine is scientifically sound and its qPCR-specific
choices (divalent salt, anneal-temperature gating, offline specificity, LNA/degenerate handling)
are correct and quantified. The claims that survive scrutiny are **correctness and integration**,
not **algorithmic novelty** or **empirical superiority** — and the benchmarks are built to show
exactly where the modest, real advantages are and where they are not.
