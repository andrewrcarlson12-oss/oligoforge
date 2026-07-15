# OligoForge 1.34.0 Release Test Report

## Release candidate

- Application version: 1.34.0
- Ranker version: 2.2.0
- Ranking schema: 1.2.0
- Objective-profile version: `2026-07-ranking-truth-3`
- Provenance schema: 1.0.0
- Experimental-feedback schema: 1.3.0
- Manual-design schema: 1.4.0
- Design-run comparison schema: 1.0.0
- Search status: heuristic-bounded

## Automated gates

| Gate | Result |
|---|---|
| Python scientific/software programs | 44 passed / 0 failed |
| Browser/UI harnesses | 9 passed / 0 failed |
| Ranking-truth invariants | 41 passed / 0 failed |
| Evidence/provenance/API invariants | 36 passed / 0 failed |
| Decision-analysis/API invariants | 25 passed / 0 failed |
| Manual Design Studio browser checks | 38 passed / 0 failed |
| Frozen ranking benchmark | Passed |
| Held-out synthetic Top-1 | 100% |
| Final synthetic-test Top-1 | 100% |
| Structured Top-1 Wilson 95% interval | 74.1%–100% (11 fixtures) |
| Python byte compilation | Passed |
| npm high-severity audit | Passed; 0 known vulnerabilities |
| Clean pinned Python environment import and `pip check` | Passed; 63 application routes and no broken requirements |

## Exercised workflows

The automated suites exercise automatic hydrolysis-probe design, SYBR design, multi-isolate inclusivity, discrimination against near neighbors, transcript-junction design, authoritative batch design, exact and near manual mapping, visible terminal 3′ mismatch placements, reverse-strand mapping, manual probe and SYBR analysis, exact edit-versus-baseline comparison, locked-primer and probe-only redesign, local/excluded-region constraints, assay diagnosis and rescue, candidate attrition, rank bands, rank-reversal scenarios, evidence-completeness states, deterministic run manifests, design-run comparison, JSON/CSV feedback import, descriptive feedback evidence summaries, duplicate/conflict detection, target-group-isolated feedback splits, benchmark leakage checks, direct exact-manual transfer into the Workbench, Workbench provenance persistence, manifest verification in HTML/CSV reports, provenance-bearing RDML 1.3, order exports, Cq, melt, standard curves, expression, multiplex calculations, hosted isolation/security, and browser event handlers.

## Determinism and uncertainty

Frozen biological winner fixtures and the ranking corpus rerun deterministically under the pinned environment. Deterministic order is separated from evidence strength: close candidates can share an equivalence group, and missing objective-required evidence produces an insufficient-evidence state. Rank-stability perturbation testing reports 79.2% Top-1 stability, with reversals concentrated in intentionally near-equivalent or knife-edge fixtures.

The synthetic/adversarial benchmark contains 11 fixtures. Its 100% Top-1 estimate has a Wilson 95% lower bound of 74.1%; the result must not be interpreted as proof of biological superiority or as an adequately powered wet-lab validation dataset.

## Network and deployment boundary

Deterministic release gates run offline. Live NCBI and remote BLAST checks are optional integration tests. The local application imports with 63 routes. The actual Render service cannot be confirmed without the repository/service URL and deployment access.

## Scientific boundary

Passing these gates proves software invariants and computational decision behavior. It does not biologically validate a primer/probe set, prove that rank 1 is the best wet-lab assay, or validate a learned reranker.
