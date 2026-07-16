# OligoForge 1.35.0 release summary

## Outcome

OligoForge 1.35.0 converts automatic design from a single long browser request into a bounded staged job, makes Manual Design evidence readable without discarding its machine record, adds a reproducible Validation Studio, and delivers an offline assay-assurance vertical slice from AssaySBOM through immutable sequence evidence, DriftGuard, OFVR, and a self-verifying evidence package.

The authoritative structured ranker remains version 2.2.0. This release does not change its weights, hard constraints, lexicographic priorities, objective profiles, or deterministic tie-breaking. Search version 2.1.1 repairs cache-sensitive candidate-corpus variation by evaluating a deterministic minimum of up to three spread-ordered target windows before honoring the soft runtime cutoff. Search remains heuristic-bounded.

## Delivered in 1.35.0

### Staged automatic design

- `POST /api/autodesign/jobs` returns a capability identifier promptly and queues sequence retrieval, candidate design/ranking, enrichment, and optional specificity as actual stages.
- The process-local backend provides a bounded queue, one scientific worker, idempotent submission, per-job capability access, required-stage deadlines, optional-BLAST deadlines, terminal expiry, sanitized snapshots, and no public list endpoint.
- A required-stage failure, timeout, or cancellation never exposes a partial finalist as complete.
- Optional BLAST failure, timeout, or cancellation retains the completed primary design with a warning. `POST /api/autodesign/jobs/{job_id}/retry-blast` retries specificity without repeating retrieval, design, or enrichment.
- The browser submits, polls, cancels, and resumes by capability identifier; preserves non-secret form inputs; shows deployment limits and actual stage state; and keeps the primary result visible when optional specificity is unavailable.
- The synchronous `/api/autodesign` route remains for compatibility and composes the same refactored scientific stage functions.

This is not a durable or distributed queue. Jobs and idempotency records are lost on restart, capability namespaces are per process, and cancellation/deadline observation is stage-boundary. A native or network call already in progress can drain before the sole worker advances.

### Manual Design evidence presentation

Ordinary analysis, redesign, rescue, edit comparison, feedback, and run-comparison views now present structured scientific evidence rather than primary raw JSON:

- mapping status and all allowed placements;
- observed/required hard requirements, rationale, remedy, and source;
- coherent target and off-target products;
- thermodynamics, structures, interactions, and robustness;
- evaluation state, rank rationale, closest-competitor evidence, reversal conditions, uncertainty, and provenance.

Machine-readable records remain available through Advanced copy/download actions. Stale-state guards prevent an edited form from silently reusing an earlier result, and Workbench transfer binds the analyzed candidate and its provenance. Placement display is capped at 100 rows and engine product records at 50; those caps are disclosed. Modified order notation is analyzed as the supported bare backbone and requires vendor/experimental confirmation.

### Validation Studio

- Normalizes complete assays with chemistry and suppresses exact molecular-and-chemistry duplicates.
- Reuses coherent complete-product reconstruction to identify candidate disagreements over supplied target and off-target cases.
- Selects a bounded case set deterministically using new candidate-pair distinctions, modeled-state diversity, group diversity, and stable tie-breaks.
- Produces candidate-interleaved 96- or 384-well layouts with 1–12 replicates, declared controls, deterministic seed, optional edge exclusion/warnings, and spreadsheet-injection-safe CSV.
- Imports completed results and labels predictions supported, contradicted, mixed, missing, invalid, or inconclusive under the declared experiment.

The selector is deterministic greedy coverage, not a proof of global experimental optimality. Interpretation is local to the supplied cases, controls, conditions, observations, and acceptance criteria; it does not alter the ranker.

### Assay Assurance vertical slice

- A versioned AssaySBOM records normalized forward/reverse/probe components, chemistry, order notation, and locks with a deterministic content identifier.
- Bounded offline FASTA/FASTA.GZ ingestion records source metadata, accepted and rejected records, exact duplicates, unique haplotypes, and group metrics separately.
- Immutable target/off-target sequence snapshots have deterministic hashes; deterministic delta records report added, removed, and unchanged unique sequences.
- DriftGuard reconstructs coherent complete products through `isolates.amplify` and emits structured, reason-coded nonnumeric states for target dropout and signal-capable off-target observations.
- OFVR creates deterministic OligoForge-local vulnerability records without representing them as an external standard.
- Evidence packages combine AssaySBOMs, snapshots, deltas, scans, OFVRs, and Validation Studio plans with per-artifact and package SHA-256 verification and escaped HTML rendering.
- JSON schemas, an offline CLI, methods, limits, threat model, source/evidence mapping, attribution, and API/deployment guidance are included.

Hashes establish content integrity, not authorship, database completeness, population representativeness, clinical validity, regulatory acceptability, or wet-lab performance.

### Release engineering

- Python and Node test discovery is recursive.
- The test runner prunes caches and dependency/build directories and applies a per-program timeout.
- A release gate rejects any source directory with 100 or more direct children.
- The performance gate uses deterministic scientific work units for its live linearity assertion while retaining the frozen environment-recorded benchmark and scientific-library acceptance checks.
- Current API, deployment, licensing/attribution, prior-art/requirements, and regulatory evidence-mapping documents cover the implemented routes and operational boundaries.

## Ranking and scientific behavior

Ranker 2.2.0 remains authoritative. Validation Studio observations, Assurance scans, OFVRs, and imported experimental feedback do not retrain, reweight, or silently change candidate order. A learned reranker remains disabled unless target-group isolation, conflict adjudication, leakage-controlled validation, calibration, ablation, and held-out improvement are demonstrated.

Search version 2.1.1 always evaluates up to the first three windows in the spread ordering (5′, 3′, midpoint) before applying the soft time budget. The ledger records `deterministic_minimum_windows`. Direct `design_assay` declares this three-window corpus; broader workflows retain their explicit larger window and retention budgets. This establishes deterministic bounded inputs, not an exhaustive or globally optimal assay search.

The frozen synthetic/adversarial ranking benchmark and Plasmodium/Haemoproteus trace remain regression evidence for software behavior. The trace was regenerated with application 1.35.0 and search 2.1.1 after the determinism repair. These fixtures do not establish improved held-out biological selection performance.

## Main files added

- Job orchestration: `oligoforge/jobs.py`, `tests/integration/test_autodesign_jobs.py`.
- Validation Studio: `oligoforge/validation_studio.py`, `VALIDATION_STUDIO_METHODS.md`, `schemas/validation_plan.schema.json`, `tests/test_validation_studio.py`.
- Assurance: `oligoforge/assurance/`, `oligoforge/assurance_cli.py`, Assurance guides/methods/limits/threat-model documents, five Assurance schemas, `templates/assaysbom_components.csv`, and `tests/test_assurance.py`.
- API/release guidance: `API.md`, `DEPLOYMENT.md`, `DATA_LICENSING_AND_ATTRIBUTION.md`, `ASSURANCE_PRIOR_ART_AND_REQUIREMENTS.md`, and `REGULATORY_EVIDENCE_MAPPING.md`.
- Release checks: `tests/integration/test_new_api.py`, `tests/test_release_engineering.py`, and `tools/check_directory_file_counts.py`.

## Main files changed

- Scientific/orchestration: `oligoforge/autodesign.py`, `oligoforge/candidate_search.py`, and `oligoforge/design.py`.
- API/UI: `app.py`, `static/index.html`, `launcher.py`, and `oligoforge/__init__.py`.
- Test/release tooling: `run_tests.py`, `run_tests.sh`, `package.json`, `package-lock.json`, `tests/test_performance.py`, UI harnesses, benchmark trace, and release documentation.

No production scientific module was intentionally removed. The final clean-archive file inventory and digest are generated during packaging rather than inferred in this summary.

## Deployment

Local use remains the safest default. The included Dockerfile and `render.yaml` support a one-process service. Set `OLIGOFORGE_HOSTED=1` for hosted restrictions, configure NCBI contact credentials through environment variables, and keep secrets outside source and job snapshots. Hosted project persistence, shared reaction-condition mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment.

Run one application process for the shipped job backend. Multiple workers or replicas create independent queues and capability namespaces; use sticky routing for a temporary constrained deployment or replace the backend with an authenticated shared durable queue and store. The application does not provide authentication, authorization, tenant isolation, TLS termination, a shared database, or regulatory electronic-record controls.

No live Render service was deployed or verified because no authorized repository/service target was provided.

## Explicitly not delivered

- Aegis multi-edit mutation search.
- Repair Compiler orchestration.
- The `assurance_futureproof` objective or FutureProof design.
- Enterprise portfolio persistence, multi-tenant access control, or a complete Assurance browser workspace.
- A public historical sequence replay or retrospective performance demonstration.

Existing constrained redesign and Assay Rescue remain available, but they are not relabeled as Aegis or Repair Compiler. No `DEMO_HISTORICAL_REPLAY.md` is included because no genuine public replay was run.

## Scientific and product limitations

- Search and retention are bounded heuristics and can miss an unretained triplet.
- A supplied sequence snapshot can be incomplete, biased, duplicated, or unrepresentative even when its hash is valid.
- Complete-product modeling does not determine amplification efficiency, fluorescence, inhibition, matrix effects, LOD/LOQ, analytical sensitivity/specificity, clinical performance, or multiplex competition.
- Polymerase-specific mismatch tolerance, modified oligos, degenerate pools, and vendor notation require appropriate experimental confirmation.
- Live NCBI and BLAST availability, content, licensing, and rate limits are external.
- No target-grouped held-out wet-lab preference dataset or validated learned reranker is included.
- Rank equivalence and uncertainty labels are conditional on modeled and supplied evidence; they are not biological-equivalence claims.
- No regulatory compliance or intended-use claim follows from the evidence mappings or package structure.

## Recommended wet-lab comparison

Preregister a target-grouped comparison of rank 1 against region/pair-diverse alternatives and an external comparator under matched templates, chemistry, concentrations, cycling, and analysis rules. Measure amplification success, efficiency, linearity, fixed-input Cq, replicate precision, product identity, inclusivity, exclusivity, probe signal, inhibition sensitivity, adequate LOD/LOQ, synthesis failures, and multiplex performance. Preserve failed assays, invalid runs, sequence-corpus provenance, and redesign history.

## Historical continuity

Version 1.34.0 added decision analysis, exhaustive allowed manual near-placement visibility, exact edit comparison, run comparison, rank-reversal conditions, context-local feedback summaries, benchmark uncertainty intervals, and authoritative Batch Design. Version 1.35.0 preserves those behaviors while adding staged jobs, evidence-oriented UI presentation, Validation Studio, and Assurance. See `CHANGELOG.md`, `RANKING_AUDIT.md`, and `RANKING_VALIDATION_REPORT.md` for the historical ranking-truth record.
