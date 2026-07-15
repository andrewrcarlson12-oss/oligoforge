# OligoForge 1.34.0 Release Summary

## Release objective

This release addresses the central product weakness: technically acceptable candidates could be returned without evidence that the highest-ranked complete assay was the strongest candidate preserved under the declared objective.

OligoForge now reports the best-supported candidate among a bounded, fully evaluated retained pool, together with assumptions, attrition, hard failures, competing alternatives, rank rationale, uncertainty, and evaluations not performed. Version 1.34 retains those protections and adds exhaustive manual near-match visibility, exact edit-versus-baseline analysis, reproducibility-aware design-run comparison, rank-reversal scenarios, local feedback evidence summaries, finite benchmark uncertainty intervals, and an authoritative evidence-preserving Batch Design path. It does not claim a universal wet-lab optimum.

## Ranking defects repaired

1. Early-window and dense-region candidate monopoly.
2. Early primer and pair truncation before complete-triplet evidence.
3. Greedy first-probe attachment.
4. One-slot probe strand-order bias.
5. Discrimination-specialist probe monopoly.
6. Target/off-target-superior probes lost outside the local thermodynamic beam.
7. Exact and near duplicates consuming retention budgets.
8. Plain top-N truncation before expensive annotation.
9. Specificity and conservation displayed after they could no longer change rank.
10. Weighted-score compensation for hard failures.
11. Degenerate primers evaluated with the unchanged probe accidentally omitted.
12. Robustness claims based only on Tm windows.
13. Finalist categories attached to non-category winners.
14. Automatic, direct, batch, viewer, manual, and rescue workflow divergence.
15. Repeated-target coordinate ambiguity.
16. Locked components disappearing from generic candidate beams.
17. Multiplex interaction underweighting.
18. False decimal-score precision and missing comparison rationale.
19. Shared mutable reaction-condition state.
20. Opaque runtime/native-call failure behavior from excessive full-annotation batches.
21. Deterministic tie-breakers displayed as stronger evidence than they contained.
22. Missing off-target, panel, or junction evidence not always suppressing comparison confidence.
23. Incomplete scientific/library/database provenance in run manifests.
24. Duplicate, conflicting, impossible, or target-leaking experimental-feedback datasets.
25. Benchmark split integrity assumed rather than machine-validated.
26. Rank manifests and evidence could be dropped when a candidate entered the Workbench or a report/export.
27. Reports could recompute Tm and structures under mutable session conditions instead of the assay's recorded design conditions.
28. Exact manually analyzed assays could not be transferred directly into the Workbench without first generating a redesign or rescue result.
29. Manual mapping hid near placements with terminal 3′ mismatches instead of showing why extension was unsafe.
30. Exact manual sequence edits lacked authoritative before/after evidence comparison.
31. Two saved runs could not be audited for winner changes, rank shifts, candidate loss, or context drift.
32. Rank explanations did not state concrete conditions under which the ordering might reverse.
33. Feedback records could be stored but lacked bounded, context-local descriptive evidence summaries.
34. Batch Design bypassed structured ranking and discarded attrition, explanations, and run provenance.

See `RANKING_AUDIT.md` for reproductions, mechanisms, and repair details.

## Evidence provenance and feedback integrity

- Consecutive candidates can share an explicit equivalence group while retaining deterministic export order.
- Required evidence is objective-aware; missing specificity, panel, or junction evaluation produces an insufficient-evidence state.
- Run manifests are deterministic, self-hashing, and independently verifiable, and record application/ranker/schema/search/retention/manual/rescue versions, native scientific-library versions, scientific model identifiers, reaction conditions, constraints, candidate limits, input hashes, external-database state, warnings, and fallbacks.
- Selecting a candidate no longer strips that provenance: automatic, direct, manual, and rescue candidates retain the manifest, objective, candidate rank, evidence vector, explanation, trace, attrition, and search state in the Workbench.
- HTML and CSV reports verify and expose the complete run ID and manifest digest; RDML 1.3 target descriptions carry compact verified provenance. Missing or altered manifests are labeled rather than trusted.
- Report thermodynamics use each assay's recorded manifest condition snapshot. Legacy assays are explicitly labeled as using the current-session fallback rather than silently mixing conditions.
- JSON and CSV feedback imports are normalized through one schema with bounded values, exact deduplication, conflict detection, completeness reporting, and versioned record hashes.
- Dataset splitting assigns complete target groups—not individual assay rows—to deterministic train, validation, or test sets.
- A learned reranker remains disabled unless the minimum evidence gate is met and subsequent leakage-controlled validation demonstrates held-out improvement.
- Benchmark manifests now fail on duplicate or missing IDs, incomplete/overlapping splits, target-group leakage, or expected candidates absent from the case.

## Manual Design Studio implemented

- Direct forward-primer, reverse-primer, optional probe, and template entry.
- Shared authoritative sequence normalization and condition snapshot.
- Every exact and allowed near placement, strand, orientation, coordinates, mismatch count, 3′ status, and uncertainty flag.
- Every coherent forward/reverse product; no silent placement selection when several exist.
- Hydrolysis-probe and SYBR analysis through the same ranker used by automatic design.
- Direct transfer of the exact analyzed assay into the Workbench with its objective, evidence, rank trace, condition snapshot, source workflow, and verified run manifest intact.
- Component locks for forward, reverse, primer pair, and probe.
- Product-length, local-shift, and excluded-region constraints.
- Locked-primer, probe-only, primer-only, pair, and local redesign.
- Immediate complete recalculation after sequence changes.
- Exhaustive display of allowed near placements, including non-extension-eligible terminal 3′ mismatches.
- Exact baseline-versus-edited comparison through the complete mapping, specificity, conservation, robustness, and ranking path.
- Side-by-side base-versus-redesign changes, retained/changed components, new risks, and rank explanations.
- Manual design designations and versioned experimental-feedback records.

## Assay Diagnosis and Rescue implemented

- Existing assay analysis with optional observed efficiency, standard-curve, melt, amplification, nonspecific-product, probe-signal, replicate, and multiplex information.
- Computational diagnoses separated from experimental inference.
- Evidence-checked one-base intended-template repair.
- Small positional shifts.
- Forward-primer replacement while retaining reverse primer.
- Reverse-primer replacement while retaining forward primer.
- Probe-only replacement while retaining both primers.
- Primer-pair replacement while retaining an existing probe when applicable.
- New amplicon within the same target.
- Explicit escalation to a new target region when no defensible local repair survives.
- Disruption order, components retained/changed, exact sequence changes, predicted improvement, new risks, evidence level, and wet-lab confirmation requirements.

## Benchmark outcome

Frozen synthetic/adversarial benchmark:

- Legacy Top-1 expected preference: 18.2%.
- Structured ranker 2.2.0 Top-1 expected preference: 100%.
- Structured pairwise preference accuracy: 100%.
- Held-out synthetic Top-1: 100%.
- Final synthetic-test Top-1: 100%.
- Structured Top-1 supplied signal/product off-target rate: 0%.
- Small-soft-perturbation Top-1 stability: 79.2%, with instability concentrated in intentionally near-equivalent or knife-edge fixtures.

This is measurable improvement in frozen software-decision fixtures. It is not evidence of improved held-out biological performance.

## Files added

- `RANKING_AUDIT.md`
- `RANKING_VALIDATION_REPORT.md`
- `RELEASE_SUMMARY.md`
- `RELEASE_TEST_REPORT.md`
- `run_tests.sh`
- `tests/run_ranking_benchmark.py`
- `tests/generate_biological_ranking_trace.py`
- `tests/test_ranking_benchmark.py`
- `tests/test_ranking_truth.py`
- `tests/test_autodesign_inclusivity.py`
- `tests/test_evidence_provenance.py`
- `oligoforge/provenance.py`
- `oligoforge/run_compare.py`
- `tests/test_decision_analysis.py`
- `tests/ui_manual_studio.js`
- `tests/fixtures/autodesign_expected_v1.33.json`
- Frozen benchmark corpus, manifest, JSON/CSV results, trace, PNG, and SVG under `tests/benchmark/`.

## Substantially changed files

- Scientific design/ranking: `oligoforge/autodesign.py`, `design.py`, `candidate_search.py`, `candidate_retention.py`, `ranking.py`, `ranking_profiles.py`, `ranking_explain.py`, `provenance.py`, `thermo.py`.
- Manual/rescue/feedback/decision audit: `manual_design.py`, `assay_rescue.py`, `experimental_feedback.py`, `run_compare.py`.
- Product/API/UI: `app.py`, `static/index.html`, `launcher.py`, `oligoforge/__init__.py`.
- Tests and release tooling: `run_tests.py`, `package.json`, `package-lock.json`, multiple scientific regression programs, and release documentation.

No scientific module was intentionally removed. The frozen automatic-design fixture was renamed from `autodesign_expected_v1.32.json` to `autodesign_expected_v1.33.json`; both biological regression programs now reference the version-matched file.

## Deployment

### Local

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
npm ci
python run_tests.py
uvicorn app:app --host 127.0.0.1 --port 8000
```

### Render

Use `render.yaml` or the included Dockerfile. Set `OLIGOFORGE_HOSTED=1`, provide an NCBI contact email through environment variables, and store all secrets in Render configuration. Hosted server-side project storage, reaction-setting mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment.

## Known limitations

- Heuristic-bounded search can miss an unretained triplet.
- No target-grouped held-out wet-lab preference dataset is included.
- A near-equivalent or insufficient-evidence label is conditional on the modeled evidence and supplied corpora, not a biological equivalence proof.
- The feedback eligibility gate is not a validated learned model.
- Live NCBI and remote BLAST availability is external to the release.
- Modified-probe and degenerate-pool predictions require vendor and experimental confirmation.
- The three-scenario robustness screen does not exhaust all reagent and matrix effects.
- Live Render deployment could not be confirmed without repository/service access.

## Recommended wet-lab comparison

Preregister a target-grouped comparison of rank 1 against diverse alternatives from ranks 2–10 and an external comparator under matched sequences, chemistry, concentrations, and annealing conditions. Measure amplification success, efficiency, linearity, fixed-input Cq, replicate precision, adequate LOD/LOQ, product identity, inclusivity, exclusivity, probe signal, inhibition sensitivity, synthesis failures, and multiplex performance. Preserve failed assays and redesign history.

## Release conclusion

The new ranker **did improve held-out synthetic/adversarial selection performance** and repaired confirmed candidate-loss mechanisms. The release does **not** establish improved held-out biological selection performance or universal assay optimality.

## 1.34.0 decision-analysis additions

- Every allowed manual near-match is visible, including failed terminal 3′ placements.
- Exact manual edits receive complete baseline-versus-edited evidence comparison.
- Complete design runs can be compared for reproducibility and rank stability.
- Rank explanations identify concrete conditions that could reverse an ordering.
- Experimental outcomes receive local, context-bounded summaries without activating an unvalidated learned model.
- Benchmark proportions include finite uncertainty intervals.
- Batch Design now uses the authoritative structured ranker and preserves evidence/provenance; its smaller search budget is explicit and auditable.

These changes improve the ability to understand, override, repair, and reproduce rankings. They do not add wet-lab labels to the benchmark and do not convert computational preference into a validated biological claim.

