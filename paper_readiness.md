# OligoForge — Paper-Readiness Assessment

**Version:** 1.30.0 · **Date:** 2026-07-03 · **Scope:** honest evaluation of what is
publishable now, what a reviewer will push back on, and what is missing.

> This document is deliberately critical. Every quantitative claim below is traceable to a
> committed artifact (benchmark CSV/JSON + test). Where OligoForge does **not** beat the
> incumbent (Primer3, Primer-BLAST, or OligoAnalyzer/MELTING), that is stated plainly — an
> honest benchmark that shows a *tie or a loss* on some axis is more useful than one engineered
> to win.

> **v1.30.0 update:** G9 (performance/runtime never benchmarked) is now **closed** — with a
> result that is deliberately unflattering. On the same templates OligoForge's exhaustive Python
> design is **65×–315× slower than Primer3's C core**, and the specificity scan, while linear
> (R²≈1.0, ~10 µs/bp), extrapolates to **~9 h for a whole 3.2 Gb human genome**. These are
> reported as limitations, not spun: OligoForge is an *interactive single-assay* tool, not a
> batch/genome-scale engine. See §2 (D9) and the revised §4/§5. The remaining code-closable gap
> toward the software-article bar is **G5 (corpus breadth)**; the principal remaining requirement
> overall is still **G1 (wet-lab validation)**, which is being closed at the bench.

> **v1.29.0 update:** three of the gaps flagged in the v1.28.0 assessment (G3 real-genome
> specificity, G6 comparator breadth, G7 LNA independent cross-check) have been **closed with
> publication-grade benchmarks**. The results are honest: on specificity OligoForge **matches**
> Primer-BLAST rather than beating it (100 % concordance, identical accuracy and identical
> failure mode); on Tm it agrees with an independent implementation of OligoAnalyzer's documented
> algorithm to a sub-degree *mean* with a wider spread; on LNA it ties MELTING against real
> experimental data. See §2 (D3, D7, D8) and the revised §4 gap ledger.

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
- **Real-genome benchmark (new in 1.29.0, `test_specificity_primerblast.py`, 9 checks):** run
  against a 13-gene human transcriptome subset containing **natural paralog families** (ACTB/ACTG1,
  YWHAZ/YWHAB/YWHAE, TUBB/TUBB4B) and a 6-taxon apicomplexan 18S panel with related-genus
  off-targets, with a faithful offline **Primer-BLAST-equivalent** (blastn 2.17.0, the same engine
  Primer-BLAST wraps) as the head-to-head comparator. See D8.

### D4 — LNA support (McTigue 2004 nearest-neighbor), independently cross-checked in 1.29.0
Full 32-term McTigue increment set on the SantaLucia core.
- vs McTigue's experimental duplex Tm's: **RMSE 1.86 °C, MAE 1.62 °C** (n = 12). **8/12 within
  ±2 °C, 9/12 within ±2.5 °C, 12/12 within ±3 °C**; worst single duplex 3.0 °C (an outlier
  *beyond* ±2 °C, flagged in the figure). Small negative bias (−0.52 °C).
- **Independent MELTING 5 cross-check (new in 1.29.0, `test_lna_hardening.py`, 13 checks):** the
  MELTING Java engine (via `rmelting`, `method.locked=mct04`) — an independent, peer-reviewed
  implementation of the *same* McTigue parameters — was run on the 12 experimental duplexes and on
  an **expanded 96-oligo panel** (length 8–22, GC 30–72 %, all four locked bases). At matched
  strand concentration the two engines' **LNA increments agree to RMSE 0.075 °C (max 0.31 °C)** —
  proving OligoForge's implementation is faithful — and against experiment OligoForge (RMSE 1.86)
  and MELTING (RMSE 2.10) are an **honest tie**, neither dominating. See D8-note.

### D5 — degenerate-primer awareness end to end
`tm_range` surfaces the true per-resolution Tm spread (validated exact vs brute-force enumeration);
genus autodesign collapses aligned variants into IUPAC codes with a min-count/min-frequency guard
against sequencing-noise false degeneracy; the real Plasmodium-cytb genus workflow reproduces its
golden (F `TTTCTACATTTACAAGGTAGCA`, amp 86, n_degenerate 4).

### D6 — engineering quality
Local-first (no cloud, no data egress — a genuine advantage for clinical/field use), single-file
no-build UI, MIQE-structured reporting, per-request-safe conditions (immutable-snapshot cache
key; 0 torn/stale reads across 20 000 read cycles), input-validation hardening (length caps,
non-finite salt rejection), and a **42-module test suite** (34 Python + 8 Node) with byte-locked
golden designs.

### D7 — OligoAnalyzer Tm head-to-head (new in 1.29.0)
IDT OligoAnalyzer is the tool most qPCR users cross-check Tm against. It is a hosted service with
no offline API (its REST endpoint needs an IDT account + OAuth), so the comparison is against its
**documented algorithm** — SantaLucia-1998 NN + Owczarzy monovalent/divalent salt, per IDT's own
OligoAnalyzer *Definitions* page — as independently implemented by **MELTING 5**. Across the 25
non-degenerate corpus oligos at matched qPCR salt (`bench_oligoanalyzer_tm.csv`,
`test_oligoanalyzer_tm.py`, 9 checks):
- **Central tendency is sub-degree:** mean |Δ| 0.73 °C, median 0.82 °C, r = 0.96.
- **But the spread is wider and reported honestly:** **7/25 oligos exceed ±1 °C**, worst = 1.74 °C,
  95 % limits of agreement −0.68 … +1.80 °C. A small **positive** mean bias (+0.56 °C) is
  consistent with OligoForge using *free* Mg²⁺ (von Ahsen) while a divalent model applied to total
  Mg runs slightly cooler.
- Framed explicitly as agreement with OligoAnalyzer's **documented algorithm**, not with the live
  hosted tool (which could not be reached). This is the expected level of agreement between two
  correct-but-not-identical NN implementations, **not** a claim of identity.

### D8 — Primer-BLAST specificity head-to-head on a real genome (new in 1.29.0)
The most-cited reviewer objection to any new specificity checker is "why not just use
Primer-BLAST?" OligoForge's offline in-silico-PCR was benchmarked directly against a faithful
offline **Primer-BLAST-equivalent** (blastn 2.17.0 + Primer-BLAST's documented 3′-anchor / total-
mismatch rules) on 11 human-transcriptome primer pairs (8 specific + 3 designed paralog cross-
reactions) plus the real pan-Plasmodium 18S case (`bench_specificity_realgenome.json`,
`test_specificity_primerblast.py`):
- **OligoForge and Primer-BLAST make IDENTICAL calls on every pair — 100 % concordance (11/11 +
  the Plasmodium case).** This is the headline: OligoForge does **not** beat Primer-BLAST, it
  *matches* it — offline, reproducibly, with no NCBI round-trip.
- Both score **sensitivity 100 %, specificity 98.5 %** vs a biologically-grounded ground truth
  (3′-contiguity priming rule), and they **share the same 2 false positives** (ACTG1, YWHAB) —
  borderline paralog cross-reactions where each flags a 3′-anchored match with an internal
  mismatch ~7 bp from the terminus. Both err toward caution, the correct behavior for assay QC.
- **Real project result:** the pan-Plasmodium 18S primers amplify all three *Plasmodium* spp. and
  do **not** cross-react with *Haemoproteus* / *Leucocytozoon* / *Toxoplasma* — genus-specific as
  designed, directly validating the Florida Scrub-Jay panel's parasite assay. (Honest nuance: the
  probe has an *isolated* affinity for *Toxoplasma* 18S, but *Toxoplasma* forms no amplicon, so it
  produces no qPCR signal — which is why amplicon-context probe evaluation beats isolated probe
  BLAST.)

**D8-note (what these three add, honestly):** D7 and D8 answer the two "why not the standard tool"
objections directly, and the answer is *not* "OligoForge is better." It is "OligoForge is offline,
reproducible, and agrees with the standard tools where they are right — including agreeing with
Primer-BLAST's own conservative bias, and tying MELTING against real LNA data." That is a
defensible, honest publication claim; a superiority claim would not be.

### D9 — runtime/performance, benchmarked honestly (new in 1.30.0)
A software article needs runtime numbers, and OligoForge's are **unflattering by design** — reported
as such rather than hidden (`bench_performance.{json,csv}`, `test_performance.py`, 15 checks). Timed
on the 18-target `bench_corpus` (real templates, GC 27–62 %) against Primer3's C core on the same
templates:
- **OligoForge design is 65×–315× slower than Primer3** — slower on **all 18** targets, no exceptions.
  This is the compute price of the tool's design: exhaustive candidate enumeration in *Python* vs
  Primer3's *C* core. The slowdown is driven **mainly by template length** (r(len, OF ms) ≈ 0.49) and
  candidate count; GC contributes **weakly-to-moderately** (r(GC, ratio) ≈ 0.38) via the anneal gate
  (D2) admitting more GC-rich candidates — the mechanism already quantified in `bench_gate_impact.csv`.
  The **extreme is clean** (most-AT-rich Plasmodium ~106–114× vs most-GC-rich Mtb rpoB ~305–315×) but
  it is explicitly **not** claimed as a monotone GC law (mid-GC points are noisy; a lucky-run monotone
  claim would be dishonest).
- **Absolute cost:** ~1.3–6 s per assay on a constrained 1-core host (sub-second to ~2.5 s on a normal
  workstation). Fine for **interactive single-assay design — the tool's actual scope.** OligoForge is
  **not** a batch or genome-scale design engine, and this benchmark says so with numbers rather than
  leaving it to a reviewer to discover.
- **Tm throughput:** ~40 µs/oligo cold, ~0.5 µs warm (lru-cached). Not a bottleneck.
- **Specificity scan scales linearly** in subject length: **R² ≈ 1.0, ~10 µs/bp** across 16 kb → 630 kb.
  This bounds it honestly: a ~5 Mb bacterial genome ≈ 45 s; a whole **3.2 Gb human genome ≈ 9 h**
  (linear projection) — impractical, which is precisely why specificity is checked against a *supplied*
  FASTA (transcriptome / target genome), not claimed as a genome-wide guarantee (consistent with D3).

**D9-note:** this is the "OligoForge loses to a comparator on some axis" result. It costs 1–2 orders of
magnitude of speed relative to Primer3, and does not scale to whole large genomes. Both are stated as
limitations. They do **not** undermine the tool's claims (correctness, parity, integration, offline
reproducibility), because those claims were never about speed or genome-scale throughput — but a
software article must report them, and now does.

---

## 3. Figures & tables inventory (paper-ready)

| # | Artifact | Type | Shows |
|---|----------|------|-------|
| F1 | `tm_comparison.png` | figure (3-panel) | NN agreement (y=x), Bland-Altman vs Primer3, gate impact by GC |
| F2 | `structure_comparison.png` | figure (2-panel) | exact 37 °C agreement; candidate-dimer flip at true Ta |
| F3 | `specificity_validation.png` | figure (2-panel) | confusion matrix (100/100); per-control-mode resolution |
| F4 | `lna_degenerate_validation.png` | figure (3-panel) | LNA vs McTigue; per-duplex error; degenerate Tm spread |
| **F5** | `oligoanalyzer_tm_comparison.png` | figure (2-panel) | **NEW** — vs MELTING (OligoAnalyzer algorithm): y=x + Bland-Altman with honest spread |
| **F6** | `primerblast_headtohead.png` | figure (2-panel) | **NEW** — per-pair concordance grid (OF = Primer-BLAST); identical sens/spec bars |
| **F7** | `realgenome_specificity.png` | figure (2-panel) | **NEW** — confusion matrix; pan-Plasmodium genus-specificity + probe nuance |
| **F8** | `lna_validation_hardened.png` | figure (3-panel) | **NEW** — OF & MELTING vs McTigue experiment; 96-oligo increment agreement; error dists |
| **F9** | `performance_benchmark.png` | figure (3-panel) | **NEW (1.30.0)** — design latency OF vs Primer3 (log-y); slowdown vs GC (weak +); linear specificity scan |
| T1 | `bench_headtohead_tm.csv` | table | per-oligo Tm: OligoForge / independent NN / Primer3 (qPCR + default salt) |
| T2 | `bench_headtohead_structure.csv` | table | per-oligo hairpin ΔG at 37 °C and Ta vs Primer3 |
| T3 | `bench_gate_impact.csv` | table | candidate admission at 37 °C vs Ta, by template GC |
| T4 | `bench_corpus_published.json` | data | 11 published assays, dual-verified (locate + DOI in OpenAlex) |
| T5 | `corpus_provenance.csv` | table | citation + accession provenance for every corpus entry |
| T6 | `specificity_validation.json` | data | 14-control confusion matrix |
| T7 | `lna_validation.json` | data | McTigue per-duplex errors, RMSE/MAE/bias |
| **T8** | `bench_oligoanalyzer_tm.csv` | table | **NEW** — per-oligo Tm: OligoForge tm_acc / nn / MELTING (owcmix08, ahs01) |
| **T9** | `bench_specificity_primerblast.csv` | table | **NEW** — per-pair subject calls: ground truth / OligoForge / Primer-BLAST |
| **T10** | `bench_specificity_realgenome.json` | data | **NEW** — confusion, concordance, shared-FP explanation, Plasmodium case |
| **T11** | `lna_validation_v2.json` + `lna_expanded_panel.csv` | data | **NEW** — MELTING cross-check (12 + 96 oligos), increment agreement |
| **T12** | `bench_performance.csv` | table | **NEW (1.30.0)** — per-target design latency: OligoForge / Primer3 / ratio, by length & GC |
| **T13** | `bench_performance.json` (+ `bench_performance_scan.csv`) | data | **NEW (1.30.0)** — driver correlations, GC extreme, Tm throughput, scan linear fit + projected envelope |

Supporting: `bench_headtohead_report.md`, `bench_report.md`. Tests pinning every finding:
`test_headtohead.py` (6 checks), `test_specificity_offline.py` (9), `test_lna_degenerate.py` (13),
`test_oligoanalyzer_tm.py` (9), `test_specificity_primerblast.py` (9), `test_lna_hardening.py` (13),
`test_performance.py` (15), `test_benchmark.py` (design-side), `test_nn.py` (LNA + NN core). (Check
counts are grep-verifiable `check(` call-sites, not runtime PASS lines, which inflate through loop
bodies.)

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

**G3 — Real-genome specificity benchmark. ✅ CLOSED in 1.29.0.** The v1.28.0 validation was 14
designed controls. It is now supplemented by a **real transcriptome benchmark** with natural
paralog off-targets (13 human genes incl. ACTB/ACTG1, YWHAZ/YWHAB/YWHAE, TUBB/TUBB4B) and a
related-genus apicomplexan 18S panel, scored against a biologically-grounded ground truth: **11
pairs, sensitivity 100 %, specificity 98.5 %** (see D8). The 2 residual false positives are
borderline paralog cross-reactions that **Primer-BLAST also flags** — a shared conservative bias,
not an OligoForge defect. *Remaining honest limit:* the natural off-target set is a curated
paralog-rich subset, not a whole-genome scan; a genome-scale ROC is still future work, but the
mechanism is now validated on real off-targets, not just constructs.

**G4 — The salt advantage is modest at matched salt.** The fair comparison is ~0.4 °C; only the
naive-default comparison reaches ~4 °C. Do not lead with the +4 °C number.

**G5 — Small, human-skewed corpus.** 11 published assays, mostly human housekeeping genes +
SARS-CoV-2 + Plasmodium. Dropped entries (UBC repeat junction, RdRp pan-sarbecovirus degenerate)
are documented, but the corpus is not a broad cross-species or cross-chemistry benchmark.

**G6 — Comparator breadth. ✅ CLOSED in 1.29.0.** v1.28.0 compared only against Primer3. Added:
(a) a **Primer-BLAST specificity head-to-head** (D8) — 100 % concordance on a real transcriptome;
and (b) an **OligoAnalyzer Tm head-to-head** (D7) via MELTING's implementation of OligoAnalyzer's
documented algorithm — sub-degree mean agreement with an honestly-reported wider spread. Both are
the tools qPCR users actually cross-check against. *Remaining honest limit:* PrimerQuest / Beacon
Designer (closed commercial tools) are still not compared, and OligoAnalyzer is compared via its
documented algorithm, not the live hosted service (which requires an IDT account). The two
highest-value comparators a reviewer asks for are now covered.

**G7 — LNA independent cross-check + panel size. ✅ CLOSED in 1.29.0.** The v1.28.0 blocker
("MELTING needs Java, unavailable") is resolved: `rmelting` + the MELTING 5 Java engine now run,
providing the **independent reference cross-check** (D4). The NN `params_lna` layer agrees with
MELTING's implementation of the same McTigue parameters to **RMSE 0.075 °C** across an **expanded
96-oligo panel** (up from 12), and ties MELTING against real experimental data (1.86 vs 2.10 RMSE).
*Remaining honest limits:* (i) the *experimental* ground truth is still the 12 McTigue duplexes —
McTigue's full 100-duplex set and Owczarzy 2011 are closed-access, so the 96-oligo expansion is
validated against MELTING (an independent predictor), not against 96 new *measurements*; (ii) the
cruder base-averaged `tm_lna` (thermo.py) still exists alongside the rigorous NN `params_lna` —
the NN path is the one validated here and used for reporting.

**G8 — "More candidates admitted" ≠ "better assay" (see D2 caveat).** The gate result is a
design-space claim, not an assay-quality claim.

**G9 — Performance/runtime — ✅ CLOSED (v1.30.0, honestly — OligoForge loses to Primer3).**
`bench_performance.py` times design, Tm, and specificity-scan on 18 real `bench_corpus` templates
(GC 27–62 %) against Primer3's C core. Findings, reported as-measured: (i) OligoForge's exhaustive
Python design is **65×–315× slower than Primer3, slower on all 18 templates**; the slowdown is
driven **mainly by template length** (r(len,OF ms)=0.49) and candidate count, with **GC a weak-to-
moderate secondary factor** (r(GC,ratio)=0.38) via the anneal gate (D2) admitting more GC-rich
candidates — the GC *extreme* is clean (AT-rich 114× vs GC-rich 315×) but this is **not** a monotone
GC law (mid-GC noisy). (ii) Absolute design cost 1.3–6.0 s/assay on a 1-core x86_64 host —
fine for **interactive single-assay use, the tool's actual scope; not a batch/genome-scale engine**.
(iii) Tm throughput ~41 µs/oligo cold, ~0.5 µs warm. (iv) The specificity scan is **linear**
(R²=0.9999, ~10 µs/bp), projecting to ~46 s for a 4.6 Mb bacterial genome but **~9 h for the whole
3.2 Gb human genome** — impractical, which is exactly why specificity is checked against a supplied
FASTA (D3), not claimed genome-wide. *Remaining honest limit:* the tool is correct and interactive-
fast but is **one-to-two orders of magnitude slower than the C incumbent and does not scale to a
whole mammalian genome** — a real ceiling, reported rather than hidden. (`test_performance.py`,
15 checks; artifacts `bench_performance.json/.csv`, `performance_benchmark.png`.)

---

## 5. Venue fit (honest)

- **JOSS** — **strong fit now.** JOSS evaluates software quality, tests, docs, and utility, not
  novel biology. OligoForge has a real user, a 42-module test suite (34 Python + 8 Node) with
  byte-locked goldens, reproducible benchmarks, and a clear scope. This is the most defensible
  near-term target. Needs: a short paper, install/usage docs, contribution guidelines.
- **Bioinformatics Advances / BMC Bioinformatics (Software article)** — **materially closer after
  1.29.0 and 1.30.0.** The comparator gaps a reviewer flags first (G3 real-genome specificity, G6
  Primer-BLAST + OligoAnalyzer head-to-heads, G7 independent LNA cross-check) are **closed with
  honest benchmarks**, and **runtime is now benchmarked (G9)** — reported honestly, including that
  OligoForge is 65–315× slower than Primer3 and does not scale to a whole human genome. What remains
  is **G1 (wet-lab / real-assay validation)** — still the single biggest limitation — plus
  **broadening the corpus (G5)**. The submittable case is now concrete: a local-first, offline,
  reproducible qPCR design+QC tool whose numbers **match the standard tools** (Primer-BLAST,
  OligoAnalyzer/MELTING, Primer3) where those are right — including matching Primer-BLAST's own
  conservative specificity bias — with its performance ceiling measured and stated rather than
  hidden. The honest framing is **parity + integration + reproducibility**, not speed, algorithmic
  novelty, or superiority — and the benchmarks now *demonstrate* that parity rather than asserting
  it. With wet-lab data from the panel being ordered, this becomes a solid submission.
- **A methods/thermodynamics journal** — **not without wet-lab data and a novel method.** The
  thermodynamics is applied, not new; the 1.29.0 benchmarks confirm parity with reference engines,
  which strengthens a *software* article but is not itself a methods contribution.

---

## 6. Bottom line

OligoForge is **publishable as a software tool (JOSS) now**. After 1.29.0 and 1.30.0 it is
**materially closer to a bioinformatics software-article bar**: the Primer-BLAST and OligoAnalyzer
comparisons a reviewer asks for, the independent LNA cross-check, and now the runtime benchmark are
done — leaving **wet-lab validation (G1)** as the principal remaining requirement (with corpus
breadth, G5, the next code increment) before a methods-journal software submission. The engine is
scientifically sound and its qPCR-specific choices (divalent salt, anneal-temperature gating,
offline specificity, LNA/degenerate handling) are correct and quantified. The claims that survive
scrutiny are **correctness, parity with the standard tools, and integration** — not **speed**,
**algorithmic novelty**, or **empirical superiority**. The benchmarks were built to test exactly
that: OligoForge matches Primer-BLAST (100 % concordance, same failure mode), agrees with
OligoAnalyzer's documented algorithm to a sub-degree mean, ties MELTING on LNA against real data —
and, measured head-to-head, is **65–315× slower than Primer3 and does not scale to a whole human
genome (~9 h projected)**. That performance ceiling is now reported, not hidden. An honest tie with
the incumbents on *correctness* — offline, reproducible, locally run, and slower — is the claim, and
it is now demonstrated rather than asserted.
