# Changelog

## v1.29.0 — gap-closing benchmarks: Primer-BLAST + OligoAnalyzer head-to-heads, real-genome specificity, independent LNA cross-check

A **benchmark-hardening** release that closes three of the gaps the v1.28.0 `paper_readiness.md`
flagged as blocking a bioinformatics software-article submission (G3 real-genome specificity, G6
comparator breadth, G7 independent LNA cross-check). **No feature was added and no validated design
changed** — the locked panel, the HMBS host regression anchor, and the two parasite autodesign
winners remain **byte-identical** (`test_locked_panel.py`, `test_regression.py`,
`test_autodesign_golden.py`). Suite grew from **30 → 33 Python scripts** (+ 8 Node UI harnesses),
all green. `paper_readiness.md` is updated (v4): G3/G6/G7 marked closed, venue verdict revised.

### Why this release

v1.28.0 established differentiation vs Primer3 and made OligoForge JOSS-ready. The next bar — a
bioinformatics software article — needs the comparisons a reviewer asks for first: *"why not just
use Primer-BLAST / OligoAnalyzer?"* and *"is the LNA layer independently verified?"* This release
answers all three with publication-grade benchmarks. The governing principle was **honesty over
salesmanship**: a benchmark that shows OligoForge *matches* or *ties* the incumbent is more useful
than one engineered to win. All three results are ties or matches — reported as such.

### G6b + G3 — Primer-BLAST specificity head-to-head on a real genome *(test_specificity_primerblast.py, 9 checks)*

- Built a faithful **offline Primer-BLAST-equivalent** oracle on **blastn 2.17.0** (the same engine
  Primer-BLAST wraps), applying Primer-BLAST's documented 3′-anchor + total-mismatch priming rules
  to filter blastn alignments to primer-relevant sites (raw `blastn-short` word-hits are *not* used
  as priming calls — they were filtered to length + 3′-anchor + terminal-mismatch, correcting an
  early over-count).
- Benchmark set: a **13-gene human transcriptome** subset with **natural paralog families**
  (ACTB/ACTG1, YWHAZ/YWHAB/YWHAE, TUBB/TUBB4B) + a 6-taxon apicomplexan 18S panel. 11 primer pairs
  (8 geNorm specific + 3 designed paralog cross-reactions) + the real pan-Plasmodium case.
- **Result (honest — a match, not a win):** OligoForge's offline in-silico-PCR and the
  Primer-BLAST-equivalent make **IDENTICAL subject calls on all 11 pairs plus the Plasmodium case —
  100 % concordance**. Both score **sensitivity 100 %, specificity 98.5 %** vs a biologically-
  grounded ground truth (3′-contiguity priming rule), and they **share the same 2 false positives**
  (ACTG1, YWHAB) — borderline paralog cross-reactions where each flags a 3′-anchored match with an
  internal mismatch ~7 bp from the terminus. Both err toward caution (correct for assay QC).
  OligoForge does not beat Primer-BLAST; it matches it — offline, reproducibly, with no NCBI call.
- **Real project result:** the pan-Plasmodium 18S primers amplify all three *Plasmodium* spp. and
  do **not** cross-react with *Haemoproteus* / *Leucocytozoon* / *Toxoplasma* — genus-specific as
  designed. Honest nuance recorded: the probe has an *isolated* affinity for *Toxoplasma* 18S, but
  *Toxoplasma* forms no amplicon → no qPCR signal (amplicon-context probe eval beats isolated BLAST).
- Artifacts: `bench_specificity_realgenome.json`, `bench_specificity_primerblast.csv`,
  `bench_specificity_records.json`, figures `primerblast_headtohead.png` + `realgenome_specificity.png`.
  Genome FASTA fixtures committed under `tests/benchmark/genome_fixtures/` (blastn DBs rebuilt at
  test time; the pinning test is pure-Python and needs no blastn).

### G6a — OligoAnalyzer Tm head-to-head *(test_oligoanalyzer_tm.py, 9 checks)*

- IDT OligoAnalyzer is a hosted tool with no offline API (its REST endpoint needs an IDT account +
  OAuth, and the domain is unreachable here), so the comparison is against its **documented
  algorithm** — SantaLucia-1998 NN + Owczarzy salt, per IDT's OligoAnalyzer *Definitions* page — as
  independently implemented by **MELTING 5**. Framed explicitly as documented-algorithm agreement,
  **not** identity with the live tool.
- **Result (honest — sub-degree mean, wider spread):** across 25 corpus oligos at matched qPCR salt,
  mean |Δ| **0.73 °C**, median **0.82 °C**, r = 0.96 — but the spread is reported honestly:
  **7/25 exceed ±1 °C**, worst 1.74 °C, 95 % limits of agreement **−0.68 … +1.80 °C**. A small
  **positive** mean bias (+0.56 °C) is consistent with OligoForge using *free* Mg²⁺ (von Ahsen)
  while a divalent model on total Mg runs slightly cooler. The panel-b figure title states "mean
  bias sub-degree; spread wider" rather than overstating uniformity.
- Artifacts: `bench_oligoanalyzer_tm.csv`, `bench_oligoanalyzer_summary.json`,
  figure `oligoanalyzer_tm_comparison.png`.

### G7 — LNA validation hardening with an independent MELTING 5 cross-check *(test_lna_hardening.py, 13 checks)*

- The v1.28.0 blocker ("MELTING needs Java, unavailable") is **resolved**: `rmelting` 1.26.0 + the
  MELTING 5 Java engine now run, giving the independent reference cross-check.
- **Primary validation unchanged and still honest:** vs McTigue 2004 **experimental** data (12
  duplexes) — RMSE **1.86 °C**, MAE **1.62 °C**, bias **−0.52 °C**; **8/12 within ±2 °C, 9/12 within
  ±2.5 °C, 12/12 within ±3 °C** (worst duplex 3.0 °C, an outlier *beyond* ±2, not hidden).
- **Independent cross-check (a tie, not a win):** at matched strand concentration (anchored by
  making MELTING's plain-DNA-core Tm equal OligoForge's — an early 2 µM setting produced a spurious
  −2.9 °C bias that was corrected), OligoForge (RMSE 1.86) and MELTING's `mct04` implementation of
  the same parameters (RMSE 2.10) are both close to experiment; **neither dominates**. The **LNA
  increment** (locked − core Tm, which cancels the backbone salt convention) agrees to **mean |Δ|
  0.03 °C, max 0.08 °C** on the 12 duplexes.
- **Panel expanded 12 → 96 oligos** (length 8–22, GC 30–72 %, all four locked bases, internal
  positions), cross-checked against MELTING: LNA-increment agreement **RMSE 0.075 °C, max 0.31 °C**,
  uniform across locked-base identity and length — proving the McTigue increment layer is faithfully
  implemented. (McTigue's full 100-duplex set and Owczarzy 2011 are closed-access, so the 96-oligo
  expansion is validated against MELTING as an independent *predictor*, not 96 new *measurements* —
  stated plainly in `paper_readiness.md`.)
- Artifacts: `lna_validation_v2.json`, `lna_expanded_panel.csv`, `lna_mctigue12_crosscheck.csv`,
  figure `lna_validation_hardened.png`.

### Tests & docs

- **+3 Python test scripts**, all deterministic and fully offline (they re-compute the OligoForge
  side and pin the committed cross-tool reference values; no Java, blastn, or network at CI time):
  `test_oligoanalyzer_tm.py` (9 checks), `test_specificity_primerblast.py` (9),
  `test_lna_hardening.py` (13). Check counts are **grep-verifiable `check(` call-sites**, not
  runtime PASS lines (which inflate through loop bodies).
- Each new test pins the *honest* number a reviewer would recompute — including the OligoAnalyzer
  7/25-beyond-±1 spread, the shared Primer-BLAST false positives, and the LNA 8/12-within-±2 count —
  so none can silently drift.
- `paper_readiness.md` v4: D7 (OligoAnalyzer) + D8 (Primer-BLAST) added; D3/D4/D6 updated;
  G3/G6/G7 marked ✅ CLOSED each with a "remaining honest limit"; venue verdict revised to
  "materially closer" with wet-lab (G1) as the principal remaining requirement.
- **CI: 33 Python + 8 Node = 41/41 green.** Goldens byte-identical.

## v1.28.0 — publication track: published-primer benchmark vs Primer3, offline in-silico-PCR specificity, LNA/degenerate validation

A **validation and differentiation** release aimed at making OligoForge publishable (JOSS now;
a bioinformatics software article with further work). No validated design changed: the locked
panel, the HMBS host regression anchor, and the two parasite autodesign winners remain
**byte-identical** (pinned in `test_locked_panel.py`, `test_regression.py`,
`test_autodesign_golden.py`). Suite grew from 27 to **30 Python scripts + 8 Node UI harnesses**,
all green. A companion `paper_readiness.md` gives the honest publishability assessment, including
where OligoForge does **not** beat Primer3.

### Why this release

v1.27.2 proved the engine self-consistent. The open question for publication was **differentiation**:
what, specifically and quantifiably, does OligoForge do that Primer3 does not — and where does that
advantage stop? This release answers that with head-to-head benchmarks against a published-primer
corpus, adds the one capability the tool was missing (offline specificity), and validates the two
under-tested paths (LNA, degenerate).

### Science — published-primer corpus + head-to-head vs Primer3 *(tests/benchmark/, test_headtohead.py)*

- **Published-primer corpus** (`bench_corpus_published.json`): **11 assays** from peer-reviewed
  sources, each **dual-verified** — every primer *locates* in its cited RefSeq/GenBank accession
  (forward + reverse-complement, correct orientation and amplicon) **and** the citation DOI
  resolves in OpenAlex. Sequences were taken from the source document / OA full text, never from
  memory (a first from-memory batch was caught with wrong amplicons and discarded). Sources:
  geNorm (Vandesompele 2002, 8 housekeeping genes), Corman 2020 (SARS-CoV-2 E + N), Kamau 2011
  (pan-Plasmodium 18S). Dropped entries documented (UBC repeat junction; RdRp pan-sarbecovirus
  degenerate reverse). Provenance in `corpus_provenance.csv`.
- **Tm head-to-head** (`bench_headtohead_tm.csv`): on 25 non-degenerate published oligos —
  - OligoForge display Tm vs an **independent** NN reimplementation (`nn.params`, a separate code
    path from `tm_acc`): **mean 0.030 °C, max 0.060 °C** — genuine two-implementation agreement.
  - vs Primer3 at **matched qPCR salt**: **+0.44 °C mean** (the divalent-aware term; Primer3's
    SantaLucia salt path is near-Mg-insensitive).
  - vs Primer3 at **its defaults**: **+3.97 °C mean** — reported explicitly as the *naive-user*
    gap (unconfigured Primer3), **not** a method-vs-method advantage.
- **Structure head-to-head** (`bench_headtohead_structure.csv`): hairpin ΔG at 37 °C matches
  Primer3 **exactly (0.0000 kcal/mol)** — same backend, so the annealing-temperature evaluation
  is the *only* difference. ΔG(Ta)−ΔG(37) = **+0.25 kcal/mol mean, up to +1.60**.
- **Design-time differentiator** (`bench_gate_impact.csv`): candidate primers admitted at the true
  annealing temperature vs a 37 °C gate scale with template GC exactly as physics requires — Mtb
  rpoB (61.9 % GC) **+49.3 %**, HMBS (52.8 %) **+32.3 %**, GNAS (56.2 %) **+30.8 %**, Plasmodium
  cytb (26.8 % AT-rich) **+0.0 %**. Reported as a candidate-admission effect, not an assay-quality
  claim (0/25 published oligos cross the structure floor at either temperature).
- Figures: `tm_comparison.png` (3-panel), `structure_comparison.png` (2-panel).

### Feature — offline in-silico-PCR specificity engine *(oligoforge/specificity.py, test_specificity_offline.py)*

The biggest capability gap: `in_silico_pcr` defaulted to remote BLAST with no deterministic offline
path. Added a self-contained engine that scans a supplied FASTA (genome/transcriptome):

- `parse_fasta()`, `scan_primer_sites()`, `in_silico_pcr_offline()`. Both strands, IUPAC-aware,
  ≤ `max_mm` internal mismatches with a **required exact 3′-terminal anchor** — the physical
  requirement for polymerase extension. The 3′ anchor is **strand-aware**: for a minus-strand hit
  the primer is reverse-complemented, so its 3′ terminus maps to the *start* of the footprint;
  this was fixed after an initial version checked the wrong end on `−` hits (a real bug, now
  regression-guarded on both orientations).
- Emits hits in the exact `{subject, lo, hi, strand, q3}` shape the existing `epcr()` /
  `assay_verdict()` consume, so all amplicon + probe-binding logic is shared and unchanged.
- Threaded `mode="offline"` + `fasta=`/`max_mm=` through `in_silico_pcr`, `assay_specificity`,
  and both API endpoints (`/api/epcr`, `/api/assay_specificity`). No network, no BLAST install.
- **Validation**: **sensitivity 100 %, specificity 100 % on 14 controls** — 8 specific assays
  clear (each hits only its own gene in an 8-gene mini-transcriptome), a near-identical
  pseudogene / a 2-mismatch paralog / an off-size misprime are all flagged, and 3′-terminal-mismatch
  variants are correctly rejected on **both** strands. Sequences committed as
  `specificity_fixture.json` (deterministic, offline). Figure: `specificity_validation.png`.
- Honest framing preserved: verdicts read "a BLAST screen, not wet-lab proof," and the capability
  is specificity against a *supplied* genome, not a genome-wide off-target guarantee.

### Science — LNA + degenerate path validation *(test_lna_degenerate.py)*

- **LNA Tm vs McTigue et al. 2004** (the primary source the 32-term NN increment set is transcribed
  from): **RMSE 1.86 °C, MAE 1.62 °C** on 12 experimental duplexes at the paper's conditions.
  **8/12 within ±2 °C, 12/12 within ±3 °C**; worst single duplex 3.0 °C (an outlier *beyond*
  ±2 °C, flagged in the figure). One assert was corrected from "all per-LNA increments > 0" to a
  mean-in-band check — McTigue's model is genuinely context-dependent and has negative increments
  (a weak-A/T-context LNA can be mildly destabilizing). MELTING/rmelting (the reference
  implementation) needs Java, unavailable here, so validation is against the embedded experimental
  subset.
- **Degenerate path stressed**: `tm_range` expands every IUPAC resolution **exactly** (verified
  vs brute-force enumeration on the display Tm scale) and caps astronomically-degenerate input
  (4¹⁵ ≈ 10⁹ resolutions) without blow-up; `_resolve` maps every IUPAC code to ACGT;
  autodesign's `_degenerate` collapses a real recurrent minor allele to an IUPAC code but ignores
  a singleton (sequencing-noise guard); the real Plasmodium-cytb genus workflow reproduces its
  golden. Figure: `lna_degenerate_validation.png`.

### Tests

- New: `test_headtohead.py` (6 checks), `test_specificity_offline.py` (9), `test_lna_degenerate.py`
  (13). All findings above are pinned. Goldens byte-identical; full suite **30 Python + 8 Node**
  green via `run_tests.py`.

### Honest limitations (see `paper_readiness.md` for the full list)

No wet-lab validation yet; OligoForge is a qPCR-correct **application layer on Primer3**, not a new
thermodynamic method; specificity is validated on a constructed control set, not a genome-wide
benchmark; the matched-salt Tm advantage is **modest (~0.4 °C)** — only the naive-default
comparison reaches ~4 °C; Primer3 is the sole comparator so far. The publishable claims are
**correctness and integration**, not algorithmic novelty or empirical superiority.

## v1.27.2 — engine validation benchmark, concurrency hardening, input validation, UI/print polish

A hardening release. No validated design changed: the locked panel, the HMBS host
regression anchor, and the two parasite autodesign winners are all byte-identical
(pinned in `test_locked_panel.py`, `test_regression.py`, `test_autodesign_golden.py`).
Suite grew from 24 to **27 Python scripts + 8 Node UI harnesses — 35/35 green**, run
by a new `run_tests.py` CI entry point.

### Science — engine validation benchmark (Track C)

**A standing 18-target benchmark now cross-validates the engine against real data and
an independent thermodynamics implementation.** *(tests/benchmark/, test_benchmark.py)*

- **What it is.** 18 real templates pulled fresh from NCBI RefSeq/GenBank (accessions +
  provenance recorded, then committed as an offline fixture so scoring is deterministic
  and network-free): apicomplexan mtDNA/rRNA (Plasmodium, Toxoplasma), 8 human
  housekeeping mRNAs, bacterial (M. tuberculosis rpoB, E. coli 16S), viral (SARS-CoV-2
  N), and plant (Arabidopsis ACT2), spanning **26.8–61.9 % GC**. One golden anchor
  (HMBS), one published pair (PrimerBank human ACTB), 16 hand/profile-rubric.
- **What it checks.** Each target is designed at its real annealing temperature; the
  displayed Tm (`tm_acc`, Biopython Owczarzy-2008 path) is cross-checked against the
  repo's **independent from-scratch NN engine** (`oligoforge.nn.params`) — genuine
  two-implementation agreement, not the same code twice.
- **Result.** Tm agrees to **max 0.060 °C / mean 0.032 °C across 37 oligos**; 18/18 valid
  designs; the anneal-temperature structure ΔG reproduces a direct primer3
  `calc_hairpin(temp_c=…)` to **0.0000 kcal/mol**; every TaqMan/GC probe sits +6–9 °C
  above its primers and avoids 5′-G; SDHA straddles its exon junction unprompted;
  degenerate Tm spans and LNA increments are reported correctly. **No engine defect
  surfaced.** The single reference divergence — OligoForge declines the published ACTB
  pair — is classification (b), *both valid*: the Tm is right (0.008 °C) and the primers
  locate (250 bp), but OligoForge applies stricter modern SYBR 3′ rules (weak GC-clamp,
  no 3′-T) that the older pair predates. Recorded, engine unchanged (no cosmetic overfit).

### Concurrency — race-free reaction conditions

**Tm and structure caches are now keyed on an immutable conditions snapshot, and
`set_conditions` swaps the salt atomically.** *(thermo.py, test_concurrency.py)*

- **The problem.** `COND` (salt) and `ANNEAL_C` are process-global and were mutated by
  `/api/conditions`; the Tm/structure `lru_cache`s were cleared on every change. Under
  the sync-endpoint threadpool (48 of 52 routes), a `/api/conditions` call concurrent
  with a design/QC request exposed two hazards: a **torn read** (a Tm computed while
  `COND` was half-updated) and a **stale cache** (a value computed under the old
  conditions served under the new ones).
- **The fix.** `set_conditions` builds a whole new `COND` dict and rebinds it in one
  atomic assignment under a writer-only lock (readers never lock). Every cached Tm/
  structure function is split into a snapshot-keyed worker + a thin public wrapper that
  folds an immutable `(mv, dv, dntp, dna, anneal)` snapshot into the cache key, so a
  conditions change yields a *new* key and can never corrupt an in-flight compute.
  Correct-by-construction, no read-path locks. **All 118 spot-checked thermo values are
  byte-identical** to the pre-refactor implementation across two conditions (same
  conditions → same snapshot → same primer3 call → same number), and a stress test shows
  **0 torn/stale reads across 20 000 read cycles** (5 threads × 4 000) under three threads flipping the
  master mix. This is a single-user local app, so shared master-mix state is intended —
  the fix guarantees it is never *torn or stale*, not per-request isolation.

### Robustness — input validation + uniform error contract

*(app.py, thermo.py, test_fuzz.py)*

- **Length caps.** A giant accidental paste no longer runs unbounded thermodynamics on a
  sync worker: `/api/qc` caps oligos at 300 nt (a primer/probe is ≤60 nt) and
  `/api/design` caps templates at 50 000 nt, each with a clean explanatory error. (A
  400 kB QC paste previously ran ~0.6 s of NN math per request.)
- **Non-finite salt.** `Infinity` / `NaN` / `-Infinity` in a `/api/conditions` body are
  rejected with `{"error": "… must be a finite number"}` (HTTP 200), not a 500.
- **Free-Mg warning.** `dNTP ≥ Mg²⁺` was silently accepted even though it drives free
  Mg²⁺ (`[Mg]−[dNTP]`, von Ahsen) to ~0 and collapses the divalent salt correction.
  Still allowed (a user may model it), but now returns a non-fatal `warning` so the Tm
  isn't misread.
- **Error contract pinned.** 13 hostile sequences × 4 endpoints + 7 bad-condition
  payloads → **zero 5xx**; every bad input yields HTTP 200 + `{error}` or a clean 4xx.
- **Documented `design_assay` contract.** `design_assay` returns `None` on failure
  (empty/too-short/all-N template); the docstring now states this and warns against
  "fixing" it to a truthy error-dict (which would silently break the existing `if not a`
  guards in `batch_design`, `autodesign`, and `design_candidates`). No public endpoint
  500s on such input; behavior unchanged.

### UI — accessibility + print, within the no-build constraint

*(static/index.html — CSS/markup only; all 8 Node UI harnesses still pass)*

- **Print styles.** The cockpit had no `@media print` rules, so PDF-ing a panel or report
  for the lab notebook printed the dark chrome, nav rail, and buttons. Print now drops
  the chrome and renders clean black-on-white content with page-break-avoid on cards/rows.
- **Accessibility.** Added a skip-to-content link, an `aria-label` on the section nav, an
  `sr-only` utility, and `id`/`tabindex` on `<main>`. (Focus-visible states,
  reduced-motion handling, and ARIA roles were already present.)
- **Micro-polish.** A subtle, theme-aware card hover-lift (respects the existing global
  reduced-motion rule). No JS function or handler-referenced element ID was touched.

### Tests / CI

- New: `test_benchmark.py` (15 asserts), `test_concurrency.py` (7), `test_fuzz.py` (12).
- New `run_tests.py` runs the full Python suite **and** the Node UI harnesses (when
  `node` is on PATH), exiting non-zero on any failure. 27 Python + 8 Node = 35/35 green.

## v1.27.1 — per-assay annealing temperature for the structure gate

Follow-up to the v1.27.0 anneal-temperature fix. That release moved
secondary-structure *rejection* off primer3's 37 °C default and onto "the
annealing temperature" — but "the annealing temperature" was a single
module-global (`ANNEAL_C`, default 60 °C) that nothing tied to the assay being
designed. This release makes the gate use **each assay's own** anneal temp.

### Scientific correctness

**Structure-rejection gates now run at the profile's annealing temperature, not
a session global.** *(profiles.py, design.py, thermo.py — `_ok_primer`,
`pair_primers`, `find_probe`, `hairpin_full` / `self_dimer_full` /
`hetero_dimer_full`)*

- **The problem.** The panel has two anneal temps — 60 °C for the host TaqMan/
  SYBR genes, ~54 °C for the AT-rich parasite mtDNA assays — but the parasite
  profiles carried "~54 °C" only in free-text `notes`; no profile carried the
  value, and `design_assay` never read one. Every design path (interactive,
  batch, autodesign) gated structure at the global 60 °C. Judging a 54 °C assay
  at 60 °C melts ~6 °C more structure than reality, so the gate was **too
  permissive** for the parasite assays: a self-dimer/hairpin that still exists at
  54 °C but is melted by 60 °C was scored as gone and admitted. This is the exact
  mirror of the 37 °C **over**-rejection fixed in 1.27.0 — same class of error,
  opposite direction, confined to the 54 °C assays. Demonstrated: on AT-rich
  oligos the self-dimer gate verdict flips between 54 °C and 60 °C for ~0.25 %
  of candidates (e.g. `ACAGGAAGCTTCATGGCCTAT`: ΔG −6.62 at 54 °C → reject,
  −4.96 at 60 °C → admit, against the −6.0 floor).
- **The fix.** Each profile carries `anneal_c` (60 °C for host/SYBR/GC, 54 °C for
  the three `parasite_*` profiles). The design gates read `c["anneal_c"]` and
  pass it explicitly into the anneal-aware structure functions, which now take
  `anneal_c` as an argument that is **part of their cache key**. The gate is
  therefore evaluated at the assay's real Ta and is independent of the mutable
  session global — so it is correct even when a host and a parasite design run
  concurrently under the endpoint threadpool (the previous global-mutation route
  would have been a data race). Thresholds and the 37 °C-basis ranking/tie-breaks
  are unchanged, preserving the "reject on physics at the real Ta, rank on the
  validated basis" separation from 1.27.0.
- **What is preserved (verified).** Every non-parasite profile is set to 60 °C,
  identical to the prior global default, so all host/SYBR/GC designs — including
  the HMBS ref-gene golden (F `GAGCTATACCCCGACCTCTG` / R
  `CTTCTCTCCAATCTTGGAAAGCG` / amplicon 93) — are byte-for-byte unchanged. The
  parasite autodesign winners were checked at 54 °C vs 60 °C on the real
  Plasmodium/Haemoproteus cytb fixtures and are **identical**: discrimination
  winner `TTTCCATTTATAGCCTTATGTATTG` (amplicon 96) with unchanged top-3 ranking
  and discrimination call, and genus winner `TTTCTACATTTACAAGGTAGCA`
  (amplicon 86). No golden needed re-capturing. The locked panel is QC-only (its
  oligos are never re-designed), so nothing that feeds a real IDT order can move.

### Tests

- Added 5 offline regression asserts (tests/test_regression.py): every profile
  carries the correct `anneal_c` (parasite 54, else 60); a self-dimer that
  survives 54 °C is rejected under the parasite profile yet admitted under the
  host profile (proving per-assay, not global, gating); the 54 °C evaluation is
  independent of the session global; and the parasite discrimination winner is
  preserved with the gate genuinely running at 54 °C. The autodesign golden now
  actually exercises the 54 °C gate its own comment already described.
- **Suite status: 24/24 test scripts pass; 4/4 Node UI harnesses pass.**

### Version

- Bumped 1.27.0 → **1.27.1** (patch: robustness/correctness fix; every validated
  and golden design is preserved unchanged). Updated in `app.py`,
  `oligoforge/__init__.py`, `launcher.py`, and `static/index.html`.

## v1.27.0 — diagnostic + scientific-correctness pass

This release is the result of a full fix→verify audit: get it running, run the
suite, exercise the real workflows, and check the **science** (thermodynamics,
salt model, and the temperature at which secondary-structure gates are judged),
not just the plumbing. Two of the changes are scientific-correctness fixes with
reasoning recorded below; one is a robustness bug fix; the rest is test
hardening and an audit that deliberately changed *nothing* because the model was
already correct.

### Scientific correctness

**1. Secondary-structure gates are now evaluated at the annealing temperature,
not primer3's 37 °C default.** *(design.py — `_ok_primer`, `pair_primers`,
`find_probe`)*

- **The problem.** Priming happens at the assay's annealing temperature
  (60 °C for the host-gene TaqMan panel, 54 °C for the parasite LNA assay), but
  the hairpin / self-dimer / cross-dimer *rejection* gates were reading ΔG at
  primer3's default 37 °C. A hairpin or dimer that is fully melted at 54–60 °C
  does not exist during priming, yet the 37 °C gate counted it against the
  candidate. The result was systematic **over-rejection** of good primers, and
  the same absolute ΔG floor was applied at a 54 °C parasite anneal as at a
  60 °C host anneal — which is not physically consistent.
- **The fix.** The reject gates now evaluate ΔG at `ANNEAL_C` using the
  anneal-temperature-aware functions the codebase already provided
  (`hairpin_full` / `self_dimer_full` / `hetero_dimer_full`, which return
  `(ΔG@37, ΔG@anneal, structure_Tm)`). Threshold numbers are unchanged; only the
  temperature at which they are judged changed. Effect on the HMBS host template:
  admitted forward-primer candidates rose from 91 to 262 (reverse 93 → 249) —
  the melted-structure primers the old gate wrongly discarded.
- **What was deliberately *not* changed.** The candidate **ranking / tie-break**
  metrics (the dimer-floor penalty in `pair_primers`, the weakest-hairpin
  tie-break in `find_probe`) still use the 37 °C basis they were calibrated on.
  Moving *ranking* to the anneal temperature shifted the pinned HMBS regression
  design; keeping ranking fixed while moving only *rejection* preserves every
  pinned/golden design exactly. This separation — "reject on physics at the real
  anneal temp, rank on the validated basis" — is the crux of the fix.
- **Preserved golden designs (verified byte-for-byte):** HMBS ref-gene TaqMan
  (F `GAGCTATACCCCGACCTCTG` / R `CTTCTCTCCAATCTTGGAAAGCG` / amplicon 93), the
  SYBR design, and both parasite autodesign winners (discrimination winner
  `TTTCCATTTATAGCCTTATGTATTG`, amplicon 96; genus winner
  `TTTCTACATTTACAAGGTAGCA`, amplicon 86) with their IUPAC-degenerate forms and
  degeneracy counts. The autodesign top-3 *ranking* pin was re-captured (winner
  unchanged; runners-up are new, legitimately-admitted candidates) per that
  test's own re-capture instruction.

**2. Salt-correction model at qPCR Mg²⁺ — audited, verified correct, left
unchanged.** *(thermo.py — audit stamp added)*

- The reported/displayed Tm everywhere the user reads or exports a number (QC,
  pair, viewer, report, MIQE — via `tm_acc` / `nn.params`) already uses the
  **divalent-aware Owczarzy-2008** salt correction with free
  Mg²⁺ = [Mg²⁺] − [dNTP] (von Ahsen 1:1 chelation), the quantity that actually
  sets duplex stability in a PCR master mix. This was **verified to within
  0.03 °C** against a from-scratch, independent Owczarzy-2008 reimplementation
  across all 14 locked-panel oligos, and it responds correctly to Mg²⁺
  (3 → 6 mM raises a 22-mer ≈ 1.3 °C) and to dNTP chelation (0.8 mM lowers it
  ≈ 0.5 °C).
- primer3's **selection** Tm is near-insensitive to Mg²⁺ (its Owczarzy-2004
  monovalent-equivalent path; Mg 1.5 → 6 mM moves it < 0.3 °C, dNTP chelation
  ignored). That is acceptable **because it is used only to rank and gate
  candidates against a fixed acceptance window — never shown as an accurate
  number.**
- **Ruling:** changing the displayed model would change nothing (it is already
  Owczarzy-2008); changing the selection model would shift the hand-validated
  locked panel and autodesign goldens for zero accuracy gain. The correct call
  is to leave the numbers as-is and document the separation, which is now
  recorded as an audit block in `thermo.py`. The Mg²⁺/dNTP sensitivity of the
  reported Tm is pinned by new regression asserts so a future edit can't silently
  break it.

### Bug fixes

**3. `intron_check` now returns a `verdict` on every path (was: `KeyError`).**
*(specificity.py, tests/test_intron.py)*

- `intron_check` returned a `verdict` key only on the success path. On any
  graceful-degradation path — exon structure unresolved, mRNA un-fetchable,
  primer not locatable in the isoform, or coordinates missing — it returned a
  dict without `verdict`, so any caller reading `r["verdict"]` (including the
  shipped `test_intron.py`) crashed with `KeyError` instead of degrading
  cleanly. Added a `verdict` key to all four degraded return paths with a clear
  "could not determine junction spanning" message. Rewrote `test_intron.py` to be
  network-tolerant and to assert the verdict-on-every-path contract as a
  regression guard.

### Tests

- Added 14 offline regression asserts (tests/test_regression.py): gate-at-anneal
  preserves the HMBS golden and admits > 150 candidates; reported Tm rises with
  Mg²⁺ and falls with dNTP chelation while the selection Tm stays Mg-insensitive;
  reported Tm reads ~1–2 °C above the selection Tm; `intron_check` verdict on
  every path; Plasmodium probe ≥ 98 % conserved across Plasmodium and ≥ 3
  mismatch / ≤ 90 % identity vs the closest Haemoproteus off-target (offline
  fixtures).
- Re-captured the autodesign top-3 discrimination ranking golden
  (tests/test_autodesign_golden.py) for the new anneal-temperature gate; winner
  and its discrimination call unchanged.
- **Suite status: 24/24 test scripts pass** (was 23/24 — `test_intron.py` was the
  failure, now fixed).

### Version

- Bumped 1.26.2 → **1.27.0** (minor: design output for new templates changes;
  all validated designs preserved). Updated in `app.py`, `oligoforge/__init__.py`,
  `launcher.py`, and `static/index.html`.
