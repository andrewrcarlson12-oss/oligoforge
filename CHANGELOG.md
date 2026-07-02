# Changelog

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
