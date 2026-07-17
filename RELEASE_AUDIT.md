# OligoForge 1.37.0 release audit

## Audited scope

Version 1.37.0 is the design-consistency and operational-reliability release. The audit covers:

- first-class Validation Studio and Assurance destinations in the main browser navigation;
- readable candidate-disagreement, case-selection, full-plate, control, interpretation, AssaySBOM, snapshot, delta, DriftGuard, OFVR, and package evidence;
- exact end-to-end API wiring, snapshot baseline/retrieval provenance, active-plan packaging, downloads, workflow status, responsive layout, and accessibility semantics;
- preservation of the 1.35 staged-design, Manual Design, lifecycle engines, authoritative ranker, provenance chain, hosted restrictions, recursive test discovery, and release directory limits.
- one verifiable quality contract across supported design surfaces, consistent probe-less objective resolution, display-independent evidence depth, stable rank ordering, and structured nested outer-primer ranking;
- one structured API problem format with request correlation, bounded browser retries, visible diagnostics, honest persistence outcomes, and continuous readiness monitoring; and
- deterministic source packaging, checksum verification, release-identity gates, and least-privilege CI/release automation.

The ranker remains version 2.2.0. The canonical, display-independent retained and annotation corpus is versioned as search 2.2.0. Version 1.37.0 does not change rank weights, objective priorities, hard constraints, or benchmark labels; it corrects inconsistent routing into those authorities and makes each execution contract explicit.

## Implemented boundary

- Job capability identifiers are unguessable access tokens; there is no public job-list endpoint.
- The queue is bounded, in memory, and served by one scientific worker. Jobs and idempotency records are lost on process restart and expire after the terminal TTL.
- Required-stage failure, timeout, or cancellation cannot be presented as a completed primary design. Optional BLAST failure, timeout, or cancellation retains the completed primary result with an explicit warning and supports BLAST-only retry.
- Cancellation and deadlines are observed at stage boundaries. A native or network call already in progress can continue draining before the sole worker accepts another job.
- Validation Studio uses deterministic bounded greedy case selection. It does not claim a globally optimal experimental design.
- DriftGuard reconstructs coherent complete products through the existing isolate engine. Its reason-coded states are computational observations, not biological-risk probabilities.
- Snapshot and evidence-package hashes provide integrity and reproducibility, not authenticity, representativeness, regulatory compliance, or wet-lab validity.
- Browser lifecycle state is session-local. The new workspaces are not an authenticated registry, durable approval system, or second scientific implementation.

## Verification status

| Gate | Status at documentation freeze |
|---|---|
| Focused automatic-design job integration scenarios | Passed |
| Focused Validation Studio and Assurance unit/API regressions | Passed |
| Ranking benchmark determinism and 1.37 search 2.2.0 biological fixture | Passed |
| Evidence/provenance regression after trace regeneration | Passed |
| Browser/UI harness set | Passed; 10 harnesses, including 20 lifecycle workflow/accessibility checks |
| Manual Design targeted browser checks | Passed; 75 checks |
| Real-DOM async job integration check | Passed; submit/idempotency/resume-state/result path exercised |
| Recursive source-tree Python suite | Passed; see final measured counts in `RELEASE_TEST_REPORT.md` |
| Clean pinned-environment install/import and `pip check` | Passed |
| Live Uvicorn health and representative API smoke test | Passed; six endpoints returned HTTP 200 (health, diagnostics, limits, profiles, contract verification, and OpenAPI) |
| npm high-severity audit | Passed; 0 vulnerabilities across 63 dependencies |
| Directory direct-child count gate | Passed |
| Staged-copy manifest, archive safety, re-extraction, and extracted-tree tests | Passed; all manifest hashes verified and the full extracted suite passed |

Exact counts, dependency versions, and extracted-tree results are recorded in `RELEASE_TEST_REPORT.md`, `RELEASE_TEST_RESULTS.json`, and `RELEASE_MANIFEST.json`. The final archive digest is recorded in the delivery handoff because an archive cannot contain its own stable digest.

## Honest omissions

The release does not implement Aegis multi-edit mutation search, Repair Compiler orchestration, the `assurance_futureproof` objective/FutureProof design, enterprise or multi-tenant assurance persistence, identity-backed approval, or durable lifecycle-session recovery. Existing constrained redesign and Assay Rescue remain distinct workflows and are not represented as those features.

No public historical sequence replay was performed, so no `DEMO_HISTORICAL_REPLAY.md` or retrospective performance claim is included. The frozen Plasmodium/Haemoproteus trace remains a software regression fixture, not a new public historical replay or prospective validation.

## Release decision

This is a computational pre-release for primer/probe design, decision support, experiment planning, and evidence packaging. The measured software and packaging gates pass. This does not establish analytical sensitivity, analytical specificity, LOD/LOQ, clinical performance, future-variant resilience, regulatory compliance, or rank-1 wet-lab superiority.

## Deployment boundary

The supplied Dockerfile, Render blueprint, and deployment guidance are documented and locally testable. No live Render service was modified or verified because no service URL, credentials, or deployment authorization were supplied. The application has no built-in authentication, tenant isolation, shared durable queue, or regulatory electronic-record controls; public or sensitive-data deployments must add them externally.
