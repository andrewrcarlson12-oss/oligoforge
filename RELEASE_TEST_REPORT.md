# OligoForge 1.37.0 release test report

## Release candidate identity

- Application version: 1.37.0
- Canonical design policy: `oligoforge-canonical-design-policy/v1`
- Design contract: `oligoforge-design-contract/v1`
- API problem contract: `oligoforge-problem/v1`
- Ranker version: 2.2.0
- Search version: 2.2.0
- Ranking schema: 1.2.0
- Objective-profile version: `2026-07-ranking-truth-3`
- Provenance schema: 1.0.0
- Search status: heuristic-bounded
- Runtime route objects with methods: 79
- OpenAPI surface: 72 paths / 75 HTTP operations

## Final automated gates

All results below were measured for this source tree on 2026-07-16 (America/New_York). Earlier-release counts were not carried forward.

| Gate | Result |
|---|---|
| Recursive Python suite (`run_tests.sh`) | Passed; 54 / 54 programs |
| Node/browser harnesses | Passed; 10 / 10 harnesses |
| Canonical design-contract regression | Passed; 17 checks |
| Structured error-contract regression | Passed; 13 checks |
| Nested outer-primer quality regression | Passed; 9 checks |
| Frozen discrimination / inclusivity order | Passed; 12 / 12 and 8 / 8 checks |
| Evidence/provenance regression | Passed; 37 checks, including current compressed biological trace |
| Frozen ranking-truth regression | Passed; 41 / 41 checks |
| Ranking benchmark reproducibility | Passed; two fresh outputs were byte-identical and the test left committed evidence unchanged |
| Async automatic-design job integration | Passed; 10 scenarios / 55 assertions |
| Manual Design browser regression | Passed; 75 checks |
| Viewer / Workbench browser regressions | Passed; 34 / 34 and 13 / 13 checks |
| Real-DOM async integration | Passed; 8 checks and 131 registered handlers audited |
| Real-DOM lifecycle integration | Passed; 20 checks |
| Release identity / source packaging | Passed; 23 checks |
| Release-manifest regression | Passed; 11 checks |
| Python byte compilation | Passed for application, package, tests, and tools |
| `pip check` | Passed; no broken requirements |
| npm high-severity audit | Passed; 0 vulnerabilities at every severity across 63 dependencies |
| Per-directory direct-child limit (<100) | Passed; maximum 66 direct files in `tests/` |
| Live Uvicorn smoke | Passed; health, diagnostics, limits, profiles, signed-contract verification, and OpenAPI returned HTTP 200 with request IDs |
| Committed release-manifest verification | Passed; every selected path, byte count, and SHA-256 matched |
| Deterministic source build | Passed; two independent builds produced byte-identical ZIPs and sidecars |
| ZIP safety and SHA-256 verification | Passed; one versioned root, no unsafe paths or symlinks, and a verified adjacent digest |
| Re-extracted archive full suite | Passed; 54 / 54 Python programs and 10 / 10 Node harnesses |

## 1.37.0 workflows exercised

The focused and recursive suites cover:

- one canonical search/ranking contract across direct, batch, automatic, sequence-viewer, manual, constrained-redesign, rescue, and nested workflows;
- display-independent retained and full-annotation depth, consistent SYBR objective resolution, stable post-annotation rank order, contract verification/comparison, and preserved Workbench provenance;
- nested outer-primer diversity retention, primer-only hard gates, target/off-target complete-product evidence, robustness, rank trace, explanation, attrition ledger, manifest, and qualification contract;
- structured validation, capacity, conflict, capability, stage, and unexpected problems with safe request identifiers, field detail, retry guidance, and recovery actions;
- browser recovery cards, escaped diagnostics, visible startup and cancellation failures, bounded transient polling retries, continuous readiness monitoring, and honest persistence/reset outcomes;
- prompt job submission, stage transitions, bounded queueing, idempotent replay/conflict, queued/running cancellation, deadlines, terminal expiry, sanitized snapshots, optional-BLAST degradation, and BLAST-only retry;
- Manual Design evidence, Validation Studio planning/plate/interpretation, Assurance snapshots/deltas/DriftGuard/OFVR/package verification, and lifecycle accessibility;
- deterministic non-mutating benchmark evidence, deterministic compressed biological trace, source inventory, secret/cache exclusions, checksums, CI tag/version enforcement, and atomic release assets.

## Search and ranking interpretation

Search 2.2.0 preserves the deterministic three-window minimum from 2.1.1 and fixes the shared automatic-design tier at up to 96 retained candidates, 30 discrimination specialists, 20 objective-aware primer pairs, and 28 fully annotated candidates independent of display count. The authoritative ranker remains 2.2.0; no rank weights, hard gates, objective priorities, or synthetic benchmark labels changed.

The wider, consistently annotated corpus changed the frozen biological finalist set. Every new discrimination finalist is hard-valid on the declared fixture, with complete target coverage, no supplied Haemoproteus product, and full condition-envelope validity. The new rank 1 has higher modeled probe identity than the former winner under the unchanged discrimination key. This is a versioned computational regression result, not a wet-lab superiority claim.

The complete biological evidence trace is stored as deterministic gzip, reducing it from about 7.2 MiB to about 0.47 MiB without removing records. The benchmark gate regenerates into temporary directories, compares all JSON/CSV/PNG/SVG/manifest bytes, and proves ordinary testing does not mutate committed release evidence.

## Runtime and dependency audit

The verified runtime was CPython 3.12.13 with Primer3-py 2.3.0, Biopython 1.87, ViennaRNA 2.7.2, FastAPI 0.139.0, Uvicorn 0.51.0, Pydantic 2.13.4, HTTPX 0.28.1, and NumPy 2.5.1. Browser checks used Node 24.14.0 and npm 11.9.0. The working pinned environment imported the application/scientific libraries and passed `pip check`; no fresh package-download claim is made for this release run.

The npm lockfile audit reported zero info, low, moderate, high, or critical findings. Source selection excludes common environment files, registry credentials, private-key containers, coverage output, caches, build output, dependencies, and runtime user-data directories. The resulting archive was independently checked for unsafe paths and symlinks. These checks do not replace a full third-party supply-chain or secret audit, and Python transitive dependencies are not fully hash-locked.

## Network and deployment boundary

Deterministic release gates run offline. Live NCBI retrieval and remote BLAST are optional external integrations and are not required for archive acceptance. The shipped job backend is process-local, single-worker, and non-durable; restart loss, stage-boundary cancellation, and horizontally scaled routing limitations remain disclosed in `DEPLOYMENT.md`.

The local live smoke verified request-correlated health, diagnostics, limits, profiles, contract verification, and OpenAPI behavior. Local BLAST was correctly reported unavailable; remote NCBI contact configuration was reported without exposing credentials. No Render deployment or merge to an auto-deploying branch is claimed.

## Scientific boundary

Passing these gates demonstrates specified software behavior, internal consistency, reproducible computational evidence handling, and release integrity. It does not biologically validate any primer/probe set, establish universal or wet-lab rank-1 optimality, prove sequence-corpus representativeness, or establish analytical, clinical, future-variant, or regulatory performance.
