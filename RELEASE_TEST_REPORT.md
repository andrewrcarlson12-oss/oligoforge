# OligoForge 1.35.0 release test report

## Release candidate identity

- Application version: 1.35.0
- Ranker version: 2.2.0 (unchanged)
- Search version: 2.1.1
- Ranking schema: 1.2.0
- Objective-profile version: `2026-07-ranking-truth-3`
- Provenance schema: 1.0.0
- Experimental-feedback schema: 1.3.0
- Manual-design schema: 1.4.0
- Validation plan schema: `oligoforge-validation-plan/v1`
- Search status: heuristic-bounded
- Application route declarations: 72

## Final automated gates

All results below were measured on 2026-07-16. Earlier-release counts were not carried forward.

| Gate | Result |
|---|---|
| Recursive Python suite (`run_tests.sh`) | Passed; 49 / 49 programs |
| Node/browser harnesses | Passed; 9 / 9 harnesses |
| Async automatic-design job integration | Passed; 10 scenarios / 55 assertions |
| Validation Studio focused regression | Passed |
| Assurance focused regression | Passed |
| New 1.35 API integration regression | Passed |
| Manual Design Studio targeted browser checks | Passed; 75 checks |
| Real-DOM async job integration | Passed; 8 checks with 111 event handlers registered |
| Frozen ranking benchmark | Passed after search 2.1.1 determinism repair |
| Evidence/provenance regression and biological trace verification | Passed; regenerated trace identifies application 1.35.0 and search 2.1.1 |
| Performance gate | Passed; deterministic scan-work linearity and frozen environment-recorded benchmark retained |
| Python byte compilation | Passed; `app.py`, `launcher.py`, `oligoforge/`, `tests/`, and `tools/` |
| `pip check` in working environment | Passed; no broken requirements |
| Clean pinned-environment install/import and `pip check` | Passed in a new Python 3.12 virtual environment |
| npm high-severity audit | Passed; 0 vulnerabilities at every severity across 63 dependencies |
| Per-directory direct-child count (<100) | Passed; 18 directories / 217 files scanned, maximum 60 in `tests/` |
| Live Uvicorn health/limits/representative API smoke | Passed; health, staged-design limits, and AssaySBOM returned HTTP 200 |
| Staged release-manifest verification | Passed; every listed path, byte count, and SHA-256 verified after extraction |
| ZIP path/symlink safety and SHA-256 verification | Passed; one release root, no unsafe paths or symlinks; archive digest recorded in the delivery handoff |
| Re-extracted archive full suite | Passed; 49 / 49 Python programs and 9 / 9 Node harnesses |

## 1.35.0 workflows exercised

The focused and recursive suites cover:

- prompt job submission, real stage transitions, bounded queueing, idempotent replay/conflict, capability lookup, queued/running cancellation, primary deadlines, terminal expiry, sanitized snapshots, optional-BLAST degradation, and BLAST-only retry;
- browser submit/poll/cancel/resume behavior, preserved non-secret inputs, visible limits, and retention of primary results when optional specificity is unavailable;
- Manual Design mapping, hard requirements, products, target/off-target evidence, thermodynamics, interactions, robustness, ranking rationale, reversal uncertainty, provenance, stale-state guards, Workbench transfer, reports, and downloads;
- Validation Studio normalization, duplicate suppression, complete-product candidate disagreement, deterministic bounded case selection, candidate-interleaved 96/384-well layouts, controls, edge warnings, CSV injection neutralization, import, and conservative interpretation;
- AssaySBOM normalization and hashing, order notation and locks, bounded FASTA/FASTA.GZ snapshots, accepted/rejected ledgers, exact deduplication, metadata groups, deterministic deltas, complete-product DriftGuard, OFVR deduplication, package verification/tamper detection, HTML escaping, schemas, APIs, and offline CLI;
- legacy automatic, batch, viewer, manual, constrained-redesign, rescue, isolate, specificity, multiplex, quantitative, RDML, report, hosted-hardening, and ranking/provenance regressions.

## Determinism repair

The search loop previously could honor its soft wall-clock cutoff before an invariant candidate corpus existed. Cache warmth could therefore change the number of windows evaluated. Search version 2.1.1 always evaluates up to the first three spread-ordered target windows (5′, 3′, midpoint) before applying the soft runtime cutoff. The ledger records `deterministic_minimum_windows`. The direct `design_assay` path declares a three-window corpus, while broader workflows continue to pass larger explicit search budgets.

This repair stabilizes the bounded search input; it does not make the search exhaustive. The authoritative ranker remains 2.2.0, and its ordering semantics did not change.

## Performance-test interpretation

The live release gate measures deterministic scientific work units rather than a cross-target wall-clock ratio that varied with host scheduling and cache state. The environment-recorded benchmark remains frozen, Primer3 acceptance of the frozen templates is checked, and the Tm smoke threshold remains enforced. These are software-performance regressions, not a clinical throughput qualification.

## Network and deployment boundary

Deterministic release gates run offline. Live NCBI retrieval and remote BLAST are optional external integrations and are not required for archive acceptance. The job backend is process-local and non-durable; restart loss, stage-boundary cancellation, single-worker throughput, and horizontally scaled routing limitations are expected behavior and are disclosed in `DEPLOYMENT.md`.

No live Render deployment result is claimed without an authorized repository and service target.

## Runtime and dependency audit

The verified runtime was Python 3.12.13 with OligoForge 1.35.0, Primer3-py 2.3.0, Biopython 1.87, ViennaRNA 2.7.2, FastAPI 0.139.0, Uvicorn 0.51.0, Pydantic 2.13.4, and HTTPX 0.28.1. Browser checks used Node 24.14.0 and npm 11.9.0. The clean environment installed directly from `requirements.txt`, imported the application and scientific libraries, and passed `pip check`.

The npm lockfile audit reported zero info, low, moderate, high, or critical findings. A heuristic repository scan found no private-key markers, credential-like filenames, strong token signatures, or quoted credential assignments. Dedicated secret scanners and a Python advisory-database scanner were not available, so these checks do not constitute a complete supply-chain or secret audit; Python transitive dependencies are not fully locked.

## Scientific boundary

Passing these gates demonstrates specified software behavior and deterministic computational evidence handling. It does not biologically validate any primer/probe set, establish that rank 1 is the best wet-lab assay, validate a learned reranker, prove sequence-corpus representativeness, or establish analytical, clinical, future-variant, or regulatory performance.
