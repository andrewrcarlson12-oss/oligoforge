# OligoForge 1.37.0 release summary

## Outcome

OligoForge 1.37.0 makes scientific behavior and operational failure handling consistent across the product. Supported primer/probe workflows now carry a canonical, machine-verifiable design contract; probe-less chemistry resolves to the same objective everywhere; display count no longer changes the retained evidence corpus; and junction annotations cannot reorder an already-ranked result. Structured problems, request IDs, recovery actions, and non-sensitive diagnostics replace opaque failure notifications with supportable evidence.

The authoritative structured ranker remains version 2.2.0. This release corrects workflow wiring and evidence consistency while preserving ranker weights, hard constraints, lexicographic priorities, and deterministic tie-breaking. It does not claim new held-out biological accuracy, wet-lab performance, or universal optimality.

## Delivered in 1.37.0

### One design-quality contract

- Each supported design result declares workflow, requested and resolved objective, chemistry identity/hash, reaction context, evidence scope, engine versions, bounded search limits, qualification state, conformance checks, and deterministic fingerprints.
- The contract deliberately hashes sequence corpora rather than embedding raw target or off-target sequences.
- Verification and comparison endpoints distinguish an intact contract from changed scientific context, evidence scope, or workflow.
- Probe-less chemistries consistently resolve generic balanced requests to the SYBR objective; unknown profiles are rejected instead of silently substituted.
- Candidate enumeration and full-annotation depth use canonical retained pools rather than changing with the requested display count.
- Junction preference remains an annotation on the authoritative ranked order, preventing stale rank numbers and downstream winner drift.
- Nested outer primers are now retained as a diverse bounded pool and evaluated by the same primer-only structured ranker under SYBR semantics; raw Tm/dimer scalar selection can no longer promote an invalid outer pair.

### Failures that can be diagnosed and recovered

- API failures use a versioned problem schema with a stable code, category, stage, retryability, recovery actions, field errors, and request identifier while retaining a top-level compatibility error string.
- Validation, capability denial, queue conflicts, lost jobs, retry exhaustion, and unexpected server failures share the same response contract.
- The browser preserves actionable error details, safely escapes server text, avoids leaking credentials to unrelated endpoints, and exposes request IDs that can be copied into support reports.
- A system-diagnostics endpoint reports non-sensitive readiness, capability, queue, condition, version, and deployment-limit state.
- State-changing settings, reset, panel, job, and cancellation actions report their true outcome instead of showing success after a rejected request.
- Automatic-job polling retries only bounded transient failures while preserving the job/stage snapshot, and the readiness lamp continues probing so a formerly healthy service cannot remain falsely green.

### Release integrity

- The desktop launcher now derives its release identity from `oligoforge.__version__`, removing a duplicate version source.
- Package metadata, browser lifecycle fixtures, README, changelog, and release summary carry the 1.37.0 identity.
- A deterministic source builder normalizes member order, timestamps, file modes, and storage and writes a standard SHA-256 sidecar.
- A path-safe verifier rejects malformed, escaping, missing, or modified checksum targets.
- Focused tests reproduce the archive twice, verify byte identity, and prove that tampering is detected.
- GitHub Actions uses least-privilege permissions, explicit timeouts, dependency integrity checks, tag/version enforcement, and source-artifact verification before release attachment.

## Inherited 1.36.0 lifecycle workspaces

### Visible Validation Studio

- A dedicated navigation destination accepts two or more pasted or Workbench assays and supplied target/near-neighbor FASTA cases.
- The form declares objective, plate format, replicates, controls, randomization seed, model bounds, Cq precision, and efficiency acceptance criteria.
- The primary result explains why cases were selected, which candidate pairs they distinguish, the modeled state for each candidate/case, and the exact plan hash.
- A full responsive 96- or 384-well map shows interleaved candidate assignments, controls, unused wells, and edge-well warnings.
- The workspace downloads the machine plan and fillable CSV, imports or accepts completed CSV, and reports control validity, supported/contradicted predictions, missing observations, and remaining uncertainty without changing the ranker.

### Visible Assurance lifecycle

- A dedicated navigation destination registers the displayed assay or entire Workbench as an AssaySBOM and shows normalized components, order sequences, locks, review state, identifiers, and content hashes.
- Users can freeze baseline and follow-up target and optional off-target FASTA evidence. Follow-ups preserve baseline linkage, and browser-originated records explicitly declare offline retrieval with `network_used: false`.
- Snapshot identity, raw/unique/duplicate/rejected counts, exact deltas, and download points are visible before scanning.
- DriftGuard runs the real complete-product API, presents state and assay support in readable tables, exposes reason records by affected assay/record/group/component, and issues local OFVRs.
- One action assembles and verifies the AssaySBOM, snapshots, deltas, scan, OFVRs, and any active Validation Studio plan. JSON and escaped review HTML are downloadable.

### Product and accessibility finish

- Both workspaces expose six-stage progress, scope boundaries, session-local storage disclosure, human-readable primary evidence, and explicit machine-download surfaces.
- Every visible lifecycle input has a programmatic accessible name; dynamic statuses use polite live announcements; each workflow exposes one current step.
- The API now carries snapshot `baseline_snapshot_id` and structured `retrieval` provenance from browser through immutable record.
- A real-DOM lifecycle harness exercises the complete visible workflow in 20 checks, including 96-well rendering, API payloads, provenance linkage, deltas, OFVRs, package verification, readable evidence, and accessibility contracts.

## Inherited 1.35.0 foundation

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

Search version 2.2.0 preserves the deterministic minimum of up to the first three windows in spread order (5′, 3′, midpoint) before applying the soft time budget and fixes the shared automatic-design tier at up to 96 retained candidates, 30 discrimination specialists, 20 objective-aware primer pairs, and 28 fully annotated candidates independent of display count. The ledger records those limits and `deterministic_minimum_windows`. Direct `design_assay` declares the three-window corpus; broader workflows retain explicit larger window budgets. This establishes deterministic bounded inputs, not an exhaustive or globally optimal assay search.

The frozen synthetic/adversarial ranking benchmark and Plasmodium/Haemoproteus trace remain regression evidence for software behavior. The trace is regenerated for application 1.37.0 and search 2.2.0 against a newly versioned expected ordering after the canonical corpus expansion. These fixtures document deterministic software behavior; they do not establish improved held-out biological selection performance.

## Foundation files retained from 1.35.0

- Job orchestration: `oligoforge/jobs.py`, `tests/integration/test_autodesign_jobs.py`.
- Validation Studio: `oligoforge/validation_studio.py`, `VALIDATION_STUDIO_METHODS.md`, `schemas/validation_plan.schema.json`, `tests/test_validation_studio.py`.
- Assurance: `oligoforge/assurance/`, `oligoforge/assurance_cli.py`, Assurance guides/methods/limits/threat-model documents, five Assurance schemas, `templates/assaysbom_components.csv`, and `tests/test_assurance.py`.
- API/release guidance: `API.md`, `DEPLOYMENT.md`, `DATA_LICENSING_AND_ATTRIBUTION.md`, `ASSURANCE_PRIOR_ART_AND_REQUIREMENTS.md`, and `REGULATORY_EVIDENCE_MAPPING.md`.
- Release checks: `tests/integration/test_new_api.py`, `tests/test_release_engineering.py`, and `tools/check_directory_file_counts.py`.

## Main files changed for 1.37.0

- Scientific consistency: `oligoforge/design_contract.py`, `oligoforge/autodesign.py`, and `oligoforge/ranking_profiles.py`.
- API and browser reliability: `oligoforge/api_errors.py`, `app.py`, and `static/index.html`.
- Release integrity: `oligoforge/__init__.py`, `launcher.py`, package metadata, GitHub Actions, deterministic source/checksum tools, and focused release regressions.

## Main files inherited from 1.36.0

- API/UI: `app.py`, `static/index.html`, `launcher.py`, and `oligoforge/__init__.py`.
- Verification: new `tests/ui_lifecycle.js`, expanded `tests/integration/test_new_api.py`, and regenerated versioned biological trace.
- Release identity and guidance: `package.json`, `package-lock.json`, README, changelog, specification, API, deployment, handoff, ranking, validation, requirements, regulatory, and release documents.

No production scientific module was intentionally removed. The committed clean-source inventory is regenerated and verified immediately before packaging; the finished archive digest is generated during packaging rather than inferred in this summary.

## Deployment

Local use remains the safest default. The included Dockerfile and `render.yaml` support a one-process service. Set `OLIGOFORGE_HOSTED=1` for hosted restrictions, configure NCBI contact credentials through environment variables, and keep secrets outside source and job snapshots. Hosted project persistence, shared reaction-condition mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment.

Run one application process for the shipped job backend. Multiple workers or replicas create independent queues and capability namespaces; use sticky routing for a temporary constrained deployment or replace the backend with an authenticated shared durable queue and store. The application does not provide authentication, authorization, tenant isolation, TLS termination, a shared database, or regulatory electronic-record controls.

No live service is deployed or verified merely by building this source release. Repository publication and any auto-deploying service remain review-gated operational actions.

## Explicitly not delivered

- Aegis multi-edit mutation search.
- Repair Compiler orchestration.
- The `assurance_futureproof` objective or FutureProof design.
- Enterprise portfolio persistence, multi-tenant access control, identity-backed approval, or durable lifecycle-session recovery.
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

Version 1.34.0 added decision analysis, exhaustive allowed manual near-placement visibility, exact edit comparison, run comparison, rank-reversal conditions, context-local feedback summaries, benchmark uncertainty intervals, and authoritative Batch Design. Version 1.35.0 added staged jobs, evidence-oriented Manual Design, Validation Studio, and Assurance engines. Version 1.36.0 made both lifecycle engines visible and operational in the browser. Version 1.37.0 adds one declared design-quality contract and one supportable failure contract across those surfaces. See `CHANGELOG.md`, `RANKING_AUDIT.md`, and `RANKING_VALIDATION_REPORT.md` for the historical ranking-truth record.
