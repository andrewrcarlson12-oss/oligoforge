# OligoForge v1.37.0 Ranking Validation Report

Version 1.37.0 retains ranker 2.2.0 and versions the display-independent candidate corpus as search 2.2.0. No rank weights, lexicographic priorities, hard constraints, objective profiles, retention priorities, or benchmark labels were changed. New conformance tests verify that supported workflows resolve effective objectives consistently, preserve the canonical search pool independent of display count, never reorder after rank annotation, and emit verifiable quality contracts. Nested outer-pair recommendations now enter the existing structured primer-only ranker rather than a separate scalar path. The changed frozen biological finalists document the wider, consistently annotated corpus; they are not new wet-lab labels.

## 1. Existing ranker description

The legacy ranker used an early scalar assembled from conservation, discrimination, sequence quality, product length, and selected dimer terms. Candidate loss occurred before complete target/off-target, robustness, junction, and panel annotation. Rank 1 therefore was not demonstrably the strongest fully evaluated triplet.

## 2. Confirmed failure mechanisms

The forensic audit confirmed early-window and early-beam loss, greedy probe selection, one-slot strand bias, specialist probe monopoly, corpus-superior probes outside the local beam, duplicate and regional budget consumption, post-hoc specificity/conservation, weighted compensation, incomplete robustness, degenerate-triplet corruption, manual/automatic divergence, ambiguous mapping, disappearing locks, panel underweighting, category mislabeling, false score precision, shared reaction state, and opaque native-call runtime. Reproductions and repairs are in `RANKING_AUDIT.md`.

## 3. New ranker description

Ranker 2.2.0 / ranking schema 1.2.0 uses:

1. true hard requirements;
2. objective-profile-specific lexicographic priorities;
3. complete target and supplied off-target product evidence;
4. three-scenario Tm, hairpin, self-dimer, cross-dimer, and primer–probe robustness screening;
5. coherent coverage and effective degenerate-triplet evidence;
6. Pareto-front annotation;
7. deterministic tie-breakers;
8. truthful finalist categories and region/pair diversity;
9. subordinate rounded display scores;
10. adjacent-candidate explanations, preference-strength states, and versioned manifests.

## 4. Objective profiles

Implemented profiles: balanced, single-target detection, broad inclusivity, discrimination, confirmatory exclusivity, screening, transcript-specific, degraded-template, multiplex, and SYBR. Profiles may change legitimate requirements and evidence order; they cannot make a true hard failure compensable.

## 5. Benchmark composition

The frozen ranking-truth corpus contains 11 synthetic/adversarial cases covering hard-gate compensation, greedy probe choice, specificity reversal, coverage reversal, condition robustness, transcript junctions, multiplex interaction, confirmatory off-target products, regional diversity, near-equivalent alternatives, and objective switching. Biological regression fixtures independently exercise AT-rich, GC-rich, transcript, paralog, published, and multi-isolate workflows.

## 6. Dataset split

- Development: 2 cases.
- Tuning: 3 cases.
- Held-out synthetic validation: 4 cases.
- Final untouched synthetic test: 2 cases.

The frozen split is recorded in `tests/benchmark/ranking_truth_corpus.json`. These are software-decision fixtures, not wet-lab labels.

## 7. Comparator settings

The legacy comparator uses the frozen preliminary scalar stored with each fixture. The new comparator uses ranker 2.2.0 and profile version `2026-07-ranking-truth-3`. Ordering is deterministic and uses no random seed.

## 8. Top-k metrics

| Metric | Legacy | Structured 2.2 |
|---|---:|---:|
| Top-1 expected preference | 18.2% (2/11) | 100% (11/11) |
| Top-3 recovery | 100% | 100% |
| Top-5 recovery | 100% | 100% |
| Top-10 recovery | 100% | 100% |
| Pairwise preference accuracy | 9.1% | 100% |
| Hard-valid, signal-clean Top-1 | 45.5% | 100% |
| Supplied signal-off-target rate at Top-1 | 9.1% | 0% |
| Supplied product-off-target rate at Top-1 | 27.3% | 0% |

Held-out synthetic Top-1 and final synthetic-test Top-1 were both 100%. This large gain is expected because the corpus targets known ranking defects; it must not be translated into a wet-lab success probability.

## 9. Pairwise preference results

The structured ranker recovered all 11 frozen expected preferences; the legacy order recovered two Top-1 choices and 1/11 pairwise preferences. Every winner change is recorded in the benchmark CSV.

## 10. Off-target results

No structured-ranker Top-1 candidate in the adversarial corpus had a supplied signal-generating or disqualifying product. Removing specificity evidence reduced Top-1 recovery to 72.7%.

## 11. Coverage results

The broad-inclusivity fixture promoted coherent complete target/probe coverage over a narrower candidate with a better preliminary scalar. Removing coverage evidence reduced Top-1 recovery to 90.9%.

## 12. Robustness results

The robust candidate defeated a nominally attractive knife-edge candidate. Removing robustness reduced Top-1 recovery to 90.9%. Under small deterministic soft-metric perturbations, overall Top-1 stability was 79.2%; instability was concentrated in intentionally near-equivalent or knife-edge cases. The UI therefore reports near-equivalence or insufficient evidence rather than biological certainty.

## 13. Runtime results

The synthetic evidence-order benchmark uses precomputed evidence vectors and is measured only as a local smoke signal, not frozen evidence. The 1.37 Plasmodium/Haemoproteus trace fully annotated 28 complete assays and returned five finalists in 172.753 seconds in the recorded release environment. Runtime depends on host load, sequence length, chemistry, corpus size, native libraries, and beam limits; it is not a throughput qualification.

## 14. Ablation results

| Feature removed | Top-1 recovery |
|---|---:|
| None | 100% |
| Specificity | 72.7% |
| Coverage | 90.9% |
| Robustness | 90.9% |
| Multiplex interaction | 90.9% |
| Junction relationship | 90.9% |

These fixtures establish that the named evidence changes the intended adversarial decisions. They do not estimate real-world biological effect size.

## 15. Legacy-versus-new ranking changes

Nine of 11 adversarial Top-1 choices changed. Changes were attributable to hard validity, paired specificity, coherent coverage, condition robustness, panel interaction, junction relationship, or regional evidence—not an unexplained score adjustment. Versioned biological winner fixtures were also updated only after the newly retained candidate received complete target/off-target and robustness annotation.

## 16. Cases where the new ranker lost

None in the frozen synthetic corpus. This is not proof that it will not lose on biological data. The repository does not yet contain a sufficiently large, leakage-controlled, target-grouped wet-lab preference dataset.

## 17. Indistinguishable cases

The near-equivalent fixture remains explicitly near-equivalent. Deterministic tie-breaking provides stable presentation without claiming a meaningful biological advantage. Perturbation testing also identifies cases where rank can reverse under plausible soft-metric uncertainty.


## Evidence uncertainty and provenance validation

The release adds a separate uncertainty layer without changing the hard/lexicographic authority of the ranker. Regression fixtures verify that:

- a hard-valid assay remains a strong preference over an invalid assay;
- candidates with small complete-evidence differences can be labeled near-equivalent;
- missing required off-target or panel evidence forces an insufficient-evidence state;
- shared equivalence groups preserve deterministic export order without implying biological separation;
- identical scientific inputs produce identical run IDs and manifest hashes;
- manifests record application, ranker, schema, native-library, scientific-model, reaction-condition, candidate-limit, constraint, input-hash, database-state, warning, and fallback provenance.

Experimental-feedback validation rejects impossible quantitative values, removes exact duplicates, surfaces conflicting outcomes, reports field completeness, and assigns whole target groups to deterministic train/validation/test sets. The minimum evidence gate does not activate a learned reranker; leakage-controlled model development, calibration, ablation, and held-out improvement remain separately required.

The benchmark corpus validator rejects duplicate/missing case identifiers, incomplete or overlapping splits, target-group leakage, and expected winners absent from the candidate set. This validates benchmark structure, not biological labels.

## 18. Limitations

- No held-out wet-lab dataset proves that rank 1 outperforms ranks 2–20.
- Published assays are viable or characterized references, not universal optima.
- External specificity depends on database coverage, versions, and query settings.
- Robustness spans three declared condition scenarios, not all chemistry interactions.
- Search is heuristic-bounded and can miss an unretained global optimum.
- Modified probes and degenerate pools require vendor review and bench validation.
- Efficiency, fluorescence, inhibition, synthesis behavior, and biological matrix effects are incompletely predictable.

## 19. Wet-lab validation requirements

Use target-group-separated assays and compare rank 1 with diverse alternatives from ranks 2–10 under matched conditions. Measure amplification success, efficiency, linearity, Cq at fixed input, replicate precision, LOD/LOQ with adequate replication, product identity, inclusivity, exclusivity, probe signal, inhibition sensitivity, and multiplex performance. Preserve failed assays and redesign history.

## 20. Historical 1.34.0 release recommendation

The 1.34.0 decision was to release a computational ranking-truth and evidence-provenance pre-release. Ranker 2.2.0 measurably improved held-out **synthetic/adversarial selection**, repaired confirmed candidate-loss mechanisms, and exposed its assumptions and attrition. It did not demonstrate improved held-out biological selection performance.

The 1.35.0 release preserved that ranking decision while adding staged orchestration and lifecycle engines. Version 1.36.0 added first-class browser workspaces over those engines. Version 1.37.0 repairs cross-workflow conformance, versions the canonical search tier as 2.2.0, and exposes machine-verifiable contracts but provides no new rank labels or biological-accuracy claim. The regenerated 1.37 biological trace is a current software regression fixture, not a public historical replay.

## Machine-readable outputs and figures

- `tests/benchmark/ranking_truth_corpus.json`
- `tests/benchmark/ranking_truth_manifest.json`
- `tests/benchmark/ranking_truth_results.json`
- `tests/benchmark/ranking_truth_results.csv`
- `tests/benchmark/ranking_truth_topk.png`
- `tests/benchmark/ranking_truth_topk.svg`
- `tests/benchmark/plasmodium_ranking_trace.json.gz`

## Provenance persistence validation

The release now tests the full chain from ranked candidate to Workbench to export. Automatic, direct exact-manual, constrained-manual, and rescue candidates retain their run manifest, objective, rank, evidence vector, trace and explanation. An exact manually entered assay can be selected directly after authoritative analysis; redesign is not required merely to preserve it. HTML/CSV reports independently verify the self-hash and label altered or missing manifests; their Tm and structure calculations use each assay's recorded condition snapshot rather than mutable session settings. RDML 1.3 descriptions preserve the run ID, ranker version, manifest SHA-256, verification state, candidate rank and objective. This prevents a reproducible rank from becoming an untraceable sequence pair after selection.

## 21. Decision-analysis validation added in v1.34.0

The release adds independent regressions for exhaustive manual near-match reporting, 3′ extension eligibility, run reproducibility classification, context-aware rank changes, concrete rank-reversal scenarios, feedback-local evidence summaries, Wilson benchmark intervals, edit-delta reporting, and authoritative batch-path provenance.

The ranking corpus remains small: 11 frozen synthetic/adversarial preference fixtures. Proportion estimates therefore include Wilson 95% intervals. An observed 11/11 Top-1 recovery has a finite lower confidence bound and must not be described as certainty about unseen biological targets. Published and biological fixtures remain viability/regression references, not labels of universal optimality.

Version 1.34.0 does **not** claim improved held-out biological selection performance. It improves decision traceability, manual repair analysis, path consistency, and honest uncertainty while preserving the ranker 2.2.0 ordering rules.
