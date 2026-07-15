# OligoForge 1.34.0 engineering handoff

Authoritative search/ranking modules are `candidate_search.py`, `candidate_retention.py`, `ranking_profiles.py`, `ranking.py`, `ranking_explain.py`, `provenance.py`, and `ranking_benchmark.py`.

Manual/evidence modules are `manual_design.py`, `assay_rescue.py`, and `experimental_feedback.py`.

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

## Workbench and exports

Ranked candidates must retain `ranker_manifest`, `objective_profile`, `candidate_rank`, `ranking_evidence`, `rank_trace`, and `rank_explanation` when entering the Workbench. HTML/CSV reports verify the manifest, recompute under its recorded reaction-condition snapshot, and RDML descriptions carry compact provenance. Do not remove these fields during UI refactors or panel import/export.

## 1.34.0 maintenance notes

- Keep `manual_design.compare_edits` on the same `analyze_assay` path; do not create a lightweight Tm-only edit checker.
- Keep manual mapping exhaustive by default and filter with `extension_eligible` only when constructing PCR products.
- Batch Design intentionally uses the same structured ranker with a smaller, manifested search-time budget. Do not restore `design_assay()` as an endpoint shortcut.
- `run_compare` candidate identity is the normalized F/R/P triplet, not the display name.
- Experimental-feedback summaries are descriptive and context-local; they must not silently feed ranker weights.

