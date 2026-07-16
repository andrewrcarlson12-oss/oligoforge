# OligoForge computational model specification

## Product boundary

OligoForge 1.36.0 is a computational workbench for primer/probe design, assay comparison, experiment planning, and local evidence packaging. It reports modeled evidence and reproducible provenance through both APIs and visible lifecycle workspaces. It does not establish wet-lab amplification, analytical or clinical performance, future-variant resilience, regulatory acceptability, or a globally optimal design.

Release identity:

- application: 1.36.0;
- authoritative ranker: 2.2.0;
- candidate search: 2.1.1;
- ranking schema: 1.2.0;
- objective profile: `2026-07-ranking-truth-3`;
- provenance schema: 1.0.0.

`API.md` specifies the HTTP surface. This document specifies the scientific and orchestration invariants behind it.

## Shared scientific invariants

1. Oligo sequences are normalized and validated through shared engine functions; the browser must not maintain a second scientific implementation.
2. Thermodynamics and structures use a request/run reaction-condition snapshot. Reports prefer the assay's recorded design conditions and label legacy current-session fallbacks.
3. A PCR observation is based on a coherent forward/reverse product. Probe signal additionally requires probe recognition within that product; independent oligo similarity is not a substitute.
4. Automatic, batch, viewer, manual, constrained-redesign, rescue, Workbench, and report paths preserve the authoritative evidence and provenance fields appropriate to their workflow.
5. Missing objective-required evidence suppresses comparison confidence. Deterministic export order is not evidence that close candidates are biologically distinguishable.

## Candidate search and ranking

The search constructs complete forward-primer/reverse-primer/probe triplets over spread-ordered target windows, applies hard screens and bounded native-library limits, suppresses exact/near duplicates, retains regional and pair diversity, completes expensive annotation on a bounded retained pool, and ranks through objective-specific structured evidence.

Search is heuristic-bounded. Search 2.1.1 evaluates up to the first three windows in spread order—5′, 3′, midpoint—before applying the soft runtime budget. The ledger records planned/evaluated windows, the deterministic minimum, runtime expiry, stage attrition, candidate limits, and retained identities. Direct `design_assay` declares three windows; broader workflows pass explicit larger limits.

The ranker remains 2.2.0. Hard constraints cannot be compensated by soft advantages. Objective profiles determine lexicographic evidence priority, followed by deterministic tie-breaks. Pareto, equivalence, missing-evidence, uncertainty, closest-competitor, and reversal information explain order without changing it. Any authoritative ordering change requires a version change and updated frozen evidence.

## Staged automatic-design jobs

### State model

The ordered stages are:

1. `resolve_fetch` (required);
2. `design` (required);
3. `enrich` (required);
4. `blast` (optional).

Public job states are queued/running and the terminal states `succeeded`, `succeeded_with_warnings`, `failed`, `timed_out`, or `cancelled`. A terminal snapshot contains timestamps, per-stage records, warnings, a sanitized structured error, and a result only when allowed by the stage outcome.

### Queue and access model

- One process owns one bounded queue and one scientific worker.
- Submission accepts `Idempotency-Key`; reuse with identical input returns the existing job, while reuse for different input is a conflict.
- Random capability identifiers are the only job lookup mechanism. There is no public list endpoint.
- Credentials, local database paths, fetched corpora, target templates, raw sequences, and native exception details are private job state and are removed from public snapshots.
- Jobs, capabilities, results, and idempotency records are non-durable and expire after a terminal TTL or process restart.

### Failure semantics

A required-stage failure, deadline, or cancellation terminates the job without presenting a partial primary result. Optional BLAST failure/deadline/cancellation returns the completed primary result as `succeeded_with_warnings`. BLAST retry creates a new capability-addressed job whose required stages are explicitly skipped and whose input is the retained primary result.

Cancellation and deadlines are observed at stage boundaries. Python does not preempt a thread inside Primer3, NCBI, or BLAST; such a call may drain after the public job becomes terminal, and the sole worker does not start another scientific stage until it returns.

This backend is neither durable nor horizontally scalable. Multiple processes create independent queues and capability namespaces.

## Manual Design and Assay Rescue

Manual Design resolves every allowed exact/near placement with strand, orientation, coordinates, mismatch count, terminal-3′ status, extension eligibility, and ambiguity. Only extension-eligible coherent primer placements form products. Hydrolysis-probe and SYBR analyses use the shared product, thermodynamic, target/off-target, robustness, constraint, and ranking paths.

The ordinary UI presents evidence as mapping, requirement, product, thermodynamic/interaction, robustness, evaluation, rank/uncertainty, and provenance cards/tables. Machine evidence remains available in Advanced exports. Editing an input invalidates stale analysis; Workbench transfer uses the analyzed record, not an unanalysed form state.

Constrained redesign supports component locks, product/local-shift/excluded-region constraints, and bounded replacement strategies. Assay Rescue separates computational diagnosis from experimental inference and orders surviving changes by disruption. Neither workflow is Aegis or Repair Compiler.

## Validation Studio model

Validation Studio accepts 2–16 distinct complete candidate assays and at most 96 supplied target/off-target cases, subject to per-record and total-base limits. Candidate identity includes normalized molecular components and chemistry.

Each candidate/case combination is evaluated with the shared coherent-product isolate model. Cases are selected deterministically by uncovered candidate-pair distinctions, number of distinct modeled states, target/off-target group diversity, and stable identifier tie-break. This greedy bounded cover is not globally optimal.

Plate generation supports 96/384 wells, 1–12 replicates, candidate interleaving, declared positive/extraction/no-template controls, a recorded deterministic seed, and optional edge exclusion or warning. CSV text is neutralized against spreadsheet formula execution.

Result interpretation is conservative. Invalid controls invalidate the comparison; otherwise observations are compared with modeled states and reported as supported, contradicted, mixed, missing, or inconclusive. Any preference is scoped to the declared experiment and never silently changes ranker 2.2.0.

The browser workspace is an input, orchestration, and evidence-presentation layer over these server authorities. It exposes candidate and case identity, selection rationale, the complete plate, plan hash, fillable CSV, control status, support/contradiction counts, and unresolved uncertainty. It does not calculate an independent browser-side scientific result.

## Assay Assurance model

### AssaySBOM

A versioned AssaySBOM normalizes the molecular components and chemistry of an assay while preserving order notation and declared component locks. Its identifier is derived from canonical content. It is a computational bill of materials, not proof of performance.

### Sequence snapshots and deltas

Offline bounded FASTA/FASTA.GZ ingestion creates immutable target or off-target snapshots. Accepted records, rejected records, exact duplicates, unique haplotypes, metadata groups, source declarations, baseline linkage, retrieval provenance, and their metrics remain distinct. Snapshot hashes detect content change; they do not establish source authenticity or corpus representativeness.

A deterministic delta compares compatible baseline and follow-up snapshots and reports added, removed, and unchanged unique sequences. It does not infer prevalence or evolutionary direction.

### DriftGuard and OFVR

DriftGuard evaluates AssaySBOM components against complete products reconstructed from supplied target/off-target snapshots. It reports structured states and reason codes, including modeled target dropout and modeled signal-capable off-target observations. It does not emit a probability-like biological-risk score.

OFVR is a deterministic, deduplicated OligoForge-local Molecular Vulnerability Record derived from declared scan evidence. It is not an external vulnerability standard.

### Evidence package

An Assurance evidence package combines declared AssaySBOMs, snapshots, deltas, scans, OFVRs, and Validation Studio plans. Per-artifact and package SHA-256 digests support independent tamper checks; HTML output escapes supplied text. Hash validity demonstrates internal integrity only, not authorship, approval, regulatory electronic-record compliance, or biological validity.

The browser Assurance workspace calls these same APIs in the explicit order register → baseline → follow-up → scan → OFVR → package. It presents primary evidence in readable tables and reason records while retaining machine artifacts as downloads. Its state is browser-session-local; clearing or reloading the session can discard unexported lifecycle state, and no enterprise registry or approval ledger is implied.

## Hosted and security boundary

Hosted mode disables shared server project storage, shared reaction-condition mutation, and local BLAST paths by default. Request-size and scientific input caps limit resource exposure, and API errors avoid reflecting secrets, local paths, or stack details.

Capability identifiers are bearer secrets, not user authentication. The application has no built-in authorization, tenant isolation, shared durable store, TLS termination, audit-signature identity, or regulatory electronic-record controls. Operators must provide those controls before a public or sensitive-data deployment.

## Features outside 1.36.0

The following are not implemented: Aegis multi-edit mutation search, Repair Compiler orchestration, the `assurance_futureproof` objective/FutureProof design, enterprise portfolio persistence, identity-backed review/approval, and public historical sequence replay. Frozen biological fixtures are software regressions, not retrospective or prospective validation studies.

## Orthogonal-panel graph model

`oligoforge/orthopanel.py` selects a mutually non-confusable subset of candidate oligos under a thresholded, pairwise thermodynamic graph model. An edge means that the modeled cross-dimer is more stable than the configured threshold. A valid panel is an independent set of that graph.

This module does not establish wet-lab multiplex compatibility. Concentration, polymerase, templates, cycling, higher-order competition, modifications, and matrix effects are outside the pairwise graph model.

### Formal outputs

The implementation reports three quantities separately:

1. **Constructive lower bound:** the size of the returned independent set. This is always a real feasible panel under the graph model.
2. **Rigorous upper bound:** either the completed exact branch-and-bound result or the size of a valid clique cover. A clique cover is an upper bound because an independent set can contain at most one vertex from each clique.
3. **Numerical Lovász-theta diagnostic:** optional floating-point SDP output. The mathematical theta value is an upper bound, but the ordinary numerical solver result is not rounded into a formal certificate in this implementation.

`certified=true` means the constructive lower bound equals a rigorous upper bound. It refers only to the graph model.

### Pipeline

1. Parse and strictly normalize candidate oligos.
2. Merge exact duplicate sequences while retaining names and counts.
3. Exclude candidates whose modeled hairpin or homodimer crosses the self-structure threshold.
4. Build all pairwise cross-dimer edges under a request-local thermodynamic-condition snapshot.
5. Solve maximum independent set by exact branch-and-bound when the configured size and expansion limits allow it; otherwise produce a greedy feasible panel.
6. Compute a greedy clique cover as a rigorous upper bound when exact search does not complete.
7. Optionally compute numerical Lovász theta for diagnostic context only.

### Split-pool output

`|panel|^k` is a constructive combinatorial count assuming independent reuse of the selected panel across `k` rounds. `theta^k`, when present, is a numerical graph diagnostic. Neither quantity is an experimentally demonstrated barcode capacity.

### Computational scope

- Graph construction is quadratic in the number of surviving candidates.
- Exact branch-and-bound is exponential in the worst case and may fall back to a greedy lower bound.
- Inputs are intended to be curated primer/probe-sized DNA oligos, not genome-scale libraries.
- Degenerate and modified oligos are reduced to the thermodynamic models explicitly supported by the application; model limitations remain visible in the output.

### Verification

Run `python tests/test_orthopanel.py` and `python tests/test_orthopanel_thermo.py`. The tests cover known graph families, exact/heuristic paths, rigorous bounds, request-local reaction conditions, input normalization, and edge construction.
