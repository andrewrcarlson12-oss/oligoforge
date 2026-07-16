# Changelog

## 1.35.0 — staged design jobs, Validation Studio, and assay assurance

- Replaced the browser's monolithic automatic-design request with a bounded, single-worker job backend that returns a capability identifier promptly and records actual retrieval, design, enrichment, and optional specificity stages.
- Added idempotent submission, bounded queueing, required-stage deadlines, optional-BLAST deadlines, cancellation, terminal TTL cleanup, sanitized snapshots, and BLAST-only retry without repeating the completed primary design.
- Preserved the legacy synchronous API for compatibility while refactoring its scientific work into the same reusable stage functions used by the job worker.
- Added browser job submission, real stage progress, cancellation, resume after navigation/reload, non-secret form preservation, visible deployment limits, and primary-result retention when optional BLAST is unavailable.
- Replaced raw Manual Design ranking structures in ordinary workflows with structured mapping, hard-requirement, target/off-target, thermodynamic, robustness, rank, uncertainty, and provenance presentations; machine evidence remains downloadable in Advanced evidence.
- Added Validation Studio candidate/case normalization, duplicate suppression, complete-product disagreement evaluation, deterministic bounded case selection, candidate-interleaved 96/384-well layouts, control handling, edge warnings, injection-safe CSV, result import, and conservative interpretation.
- Added the first offline OligoForge Assurance workflow: versioned AssaySBOM, bounded FASTA/FASTA.GZ snapshots, accepted/rejected record ledgers, exact deduplication, immutable hashes, deterministic snapshot deltas, complete-product DriftGuard, OFVR records, and self-verifying evidence packages.
- Added JSON schemas, offline Assurance CLI, methods/limits/threat-model documentation, API and deployment guidance, current-source evidence mapping, and data attribution.
- Made Python and Node test discovery recursive and added a release gate that rejects any source directory containing 100 or more direct child files.
- Bumped candidate search to 2.1.1 and made its minimum candidate corpus deterministic: up to the first three spread-ordered target windows (5′, 3′, midpoint) are evaluated before the soft runtime cutoff. The ledger records this minimum; search remains heuristic-bounded.
- Replaced a cache/order-sensitive live wall-clock performance ratio with deterministic scan-work linearity while retaining the frozen environment-recorded benchmark and scientific-library acceptance checks.
- The structured ranker remains version 2.2.0; this release does not change authoritative candidate ordering or claim new held-out biological accuracy.
- The shipped job backend is explicitly process-local and non-durable, with stage-boundary cancellation and one scientific worker; jobs are lost on restart and horizontal scaling requires a shared authenticated queue/store.
- Aegis multi-edit search, Repair Compiler orchestration, FutureProof/`assurance_futureproof`, enterprise Assurance persistence, a complete Assurance browser workspace, and public historical sequence replay are not included.

## 1.34.0 — decision analysis, edit comparison, and path unification

- Manual mapping now reports every allowed near-match, including primer placements with a terminal 3′ mismatch; such placements are explicitly marked non-extension-eligible rather than hidden.
- Added exact baseline-versus-edited assay comparison using the same mapping, product, thermodynamic, conservation, specificity, robustness, and structured-rank calculations as automatic design.
- Edit comparison reports sequence operations, resolved and introduced hard failures, improvements, worsenings, unique target mapping, and conditional rank preference.
- Added concrete evidence-linked rank-reversal scenarios instead of generic uncertainty boilerplate.
- Added design-run comparison for manifest validity, scientific-context changes, candidate additions/losses, winner changes, shared-candidate rank shifts, Spearman stability, and pair-order reversals.
- Added local experimental-evidence summaries and matched-context pairwise preferences without silently modifying the authoritative ranker.
- Added Wilson 95% intervals to ranking benchmark proportions so the small adversarial corpus is not presented with false certainty.
- Unified Batch Design with the authoritative retained-pool structured ranker. Batch winners now retain hard-validity evidence, explanations, attrition, and a self-hashing manifest.
- Added an eight-template request cap and declared 12-second per-template batch enumeration budget; the search limit is recorded in candidate provenance.
- Removed the remaining unused first-pair/greedy-probe helper from the automatic-design module.
- Added decision-analysis regressions and browser integration coverage for edit comparison and authoritative batch design.

## 1.33.0 — evidence provenance, rank uncertainty, and feedback integrity

- Added explicit near-equivalent rank bands and an insufficient-evidence state when required specificity, panel, or junction evidence is missing.
- Added deterministic, verifiable self-hashing run manifests with application, ranker, schema, scientific-model, native-library, condition, input-hash, constraint, candidate-limit, database-state, warning, and fallback provenance.
- Propagated manifests through automatic, manual, constrained-redesign, assay-rescue, Workbench selection, HTML/CSV reports, and RDML 1.3 exports.
- Added validated JSON/CSV experimental-feedback import, quantitative range checks, exact deduplication, conflicting-outcome detection, completeness reporting, and deterministic target-group-isolated train/validation/test splitting.
- Kept the structured computational ranker authoritative; a learned reranker remains disabled unless minimum evidence and leakage-control gates are met.
- Added benchmark-corpus validation for unique identifiers, complete/disjoint splits, target-group leakage, and valid expected candidates.
- Added browser controls for feedback dataset audit, split, and export.
- Direct exact manual analyses can now be added to the Workbench without requiring a redesign, while preserving the objective, complete evidence, condition snapshot, workflow, and self-verifying run manifest.
- Added 36 evidence/provenance/API invariants and expanded browser integration coverage, including Workbench chain-of-custody tests.

## 1.32.0 — ranking-truth, manual design and assay rescue

- Replaced first-hit/window-biased ranking with target-spanning, bounded joint triplet search.
- Retain multiple primer pairs and multiple probes per pair, including target/off-target-aware probe augmentation.
- Added exact and near-duplicate suppression, regional quotas and machine-readable attrition ledgers.
- Added objective profiles, hard constraints, lexicographic evidence ranking, Pareto fronts and deterministic manifests.
- Full target/off-target products, coherent isolate coverage, three-scenario Tm/structure/interaction robustness, junction and panel evidence can change final rank.
- Added rank explanations, closest-competitor rationale, uncertainty categories, truthful category winners, and pair/region-diverse finalists.
- Added authoritative manual assay mapping and analysis with every exact/near placement.
- Added locked forward, reverse, pair and probe redesign; local shift, excluded-region and amplicon constraints.
- Added assay diagnosis and minimally disruptive rescue plans.
- Added versioned experimental-feedback records without enabling an unvalidated learned reranker.
- Added frozen ranking benchmark, ablations, trace files and validation documentation.
- Preserved the v1.31.1 scientific-correctness, hosted-security, RDML, export and quantitative-analysis repairs.

## 1.31.1 — scientific-correctness and release-hardening audit

See repository history and `RELEASE_AUDIT.md` for the prior release record.
