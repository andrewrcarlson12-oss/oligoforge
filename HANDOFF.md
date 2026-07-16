# OligoForge 1.36.0 engineering handoff

Release identity: application 1.36.0, authoritative ranker 2.2.0, search 2.1.1, ranking schema 1.2.0, and objective profile `2026-07-ranking-truth-3`.

Authoritative search/ranking modules are `candidate_search.py`, `candidate_retention.py`, `ranking_profiles.py`, `ranking.py`, `ranking_explain.py`, `provenance.py`, and `ranking_benchmark.py`.

Manual/evidence modules are `manual_design.py`, `assay_rescue.py`, and `experimental_feedback.py`.

Staged automatic-design orchestration is in `jobs.py`; scientific stage implementations remain in `autodesign.py`. Validation planning is in `validation_studio.py`. Assurance records are under `oligoforge/assurance/`, with the offline entry point in `assurance_cli.py`.

Do not reintroduce a UI-only score or a parallel Tm, specificity, sequence-normalization, coordinate, or ranking implementation. Automatic, batch, viewer, manual, rescue, API, and report paths must continue to use the same scientific modules.

Any authoritative rank-order change requires:

1. a new ranker/profile or scientific-model version;
2. updated frozen benchmark and biological traces;
3. a documented winner-change rationale;
4. deterministic regression tests;
5. confirmation that hard constraints and evidence-completeness semantics were not weakened.

A deterministic display order is not proof that close candidates are biologically distinguishable. Preserve equivalence groups and explicit missing-evidence states.

Experimental feedback is assay-specific evidence. Do not activate a learned reranker merely because records exist; require target-group isolation, outcome balance, conflict adjudication, leakage-controlled validation, calibration, ablation, and held-out improvement.

Current status: computational pre-release. Synthetic/adversarial ranking behavior improved; held-out wet-lab rank-1 superiority is not established.

## 1.36.0 operational invariants

- The in-memory job manager has one scientific worker, matching the process-wide Primer3/Entrez constraints. Do not increase worker count without removing shared native/global state and adding isolation tests.
- The shipped queue is process-local and non-durable. Restart loses jobs and idempotency records; terminal records expire; multiple processes have independent queues and capability namespaces. Use one process, or replace this backend with an authenticated shared queue/store before horizontal scaling.
- Job identifiers are capability tokens. Do not add an unauthenticated list endpoint or place target sequences, credentials, local paths, or native exception details in job snapshots.
- A required-stage timeout/cancel must never expose a partial finalist as complete. Optional BLAST may fail while retaining the completed primary result with an explicit warning.
- Cancellation and deadlines are stage-boundary signals. A native or network call already in progress may drain before the single worker advances; do not claim preemptive termination.
- Snapshot and AssaySBOM identifiers verify canonical content. Migration or schema changes require new versions and deterministic-hash regressions.
- DriftGuard uses reconstructed coherent products from `isolates.amplify`; do not replace this with independent oligo-percent-identity thresholds or a cosmetic biological-risk score.
- Validation Studio selection is a bounded deterministic greedy coverage method, not a globally optimal experiment claim.
- The Validation Studio and Assurance pages are first-party clients of server authorities. Keep primary evidence human-readable, downloads machine-complete, workflow state explicitly session-local, snapshot baseline/retrieval fields intact, and every visible control programmatically named.

## Search determinism invariant

Search 2.1.1 evaluates up to the first three spread-ordered target windows (5′, 3′, midpoint) before applying its soft runtime budget and records `deterministic_minimum_windows`. The direct `design_assay` path intentionally declares three windows; automatic and other broader paths pass larger explicit budgets. Do not move the wall-clock cutoff ahead of this invariant corpus. Any change requires cold-cache/warm-cache determinism regressions and regenerated traces.

This repair does not make search exhaustive and does not change ranker 2.2.0. Keep search-version changes separate from ranker/profile changes in manifests and release claims.

## Workbench and exports

Ranked candidates must retain `ranker_manifest`, `objective_profile`, `candidate_rank`, `ranking_evidence`, `rank_trace`, and `rank_explanation` when entering the Workbench. HTML/CSV reports verify the manifest, recompute under its recorded reaction-condition snapshot, and RDML descriptions carry compact provenance. Do not remove these fields during UI refactors or panel import/export.

Ordinary Manual Design panels now render evidence cards/tables; machine records remain available through Advanced copy/download actions. Preserve the stale-state guards and bind Workbench transfer to the analyzed record, not the current editable form. The UI display cap for placements and engine cap for recorded products must remain visible when reached.

## Assurance and Validation Studio boundaries

- AssaySBOM, snapshot, delta, DriftGuard, OFVR, and package records are deterministic local evidence formats. Hash verification establishes integrity, not authenticity, representativeness, biological validity, or regulatory compliance.
- Snapshot ingestion must preserve accepted/rejected and raw/unique/group counts separately; do not silently discard duplicates or rejected-record evidence.
- Validation Studio interpretations are scoped to the declared experiment and must not feed rank order automatically.
- Keep formulas and user text escaped in CSV/HTML outputs and preserve package tamper verification.

## Explicit omissions

Aegis multi-edit mutation search, Repair Compiler orchestration, the `assurance_futureproof` objective/FutureProof design, enterprise portfolio persistence, and identity-backed review/approval are not implemented. Existing constrained redesign and Assay Rescue must not be renamed or documented as those workflows.

No public historical sequence replay was performed for 1.36.0. The frozen biological trace is a regression fixture only; do not add `DEMO_HISTORICAL_REPLAY.md` or a retrospective performance claim without a genuine source-versioned public replay and documented methods.

## 1.34.0 maintenance notes

- Keep `manual_design.compare_edits` on the same `analyze_assay` path; do not create a lightweight Tm-only edit checker.
- Keep manual mapping exhaustive by default and filter with `extension_eligible` only when constructing PCR products.
- Batch Design intentionally uses the same structured ranker with a smaller, manifested search-time budget. Do not restore `design_assay()` as an endpoint shortcut.
- `run_compare` candidate identity is the normalized F/R/P triplet, not the display name.
- Experimental-feedback summaries are descriptive and context-local; they must not silently feed ranker weights.
