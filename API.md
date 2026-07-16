# OligoForge 1.35.0 API

This document describes the HTTP surface implemented for the 1.35.0 release. The executable authority is `app.py`; request and response details generated from that code are available at `/openapi.json`, `/docs`, and `/redoc` while the service is running. The JSON schemas under `schemas/` are authorities for Assurance artifacts, not substitutes for the HTTP OpenAPI document.

OligoForge is research and engineering software. An API response is not analytical, clinical, or regulatory evidence by itself.

## Service contract

- Local Docker default: `http://127.0.0.1:8111` when port 8111 is published from the container.
- Content type: JSON for API requests and responses unless an endpoint returns an encoded export inside JSON.
- Authentication: none is implemented. Do not expose the service publicly without an authentication and authorization layer.
- Request limit: `OLIGOFORGE_MAX_REQUEST_BYTES`, default 5 MiB. An oversized declared body receives HTTP 413; an invalid `Content-Length` receives HTTP 400.
- Validation: Pydantic request errors receive HTTP 422 with field locations and messages but not rejected values. Unhandled errors receive a generic HTTP 500 response; details remain in server logs.
- Legacy domain errors: several scientific endpoints return `{"error": "..."}` with HTTP 200. Clients must inspect the response body and must not treat status 200 alone as scientific success.
- Security headers: responses set content-type sniffing, framing, referrer, permissions, and content-security-policy headers. TLS, identity, authorization, rate limiting, audit logging, backups, and tenant isolation remain deployment responsibilities.
- Hosted mode: `OLIGOFORGE_HOSTED=1` disables server project/panel storage, shared reaction-condition mutation, and client-supplied local BLAST paths by default.

## Endpoint catalog

The request-model column names the model in `app.py`. Required fields and all defaults remain machine-readable in `/openapi.json`.

### Service, profiles, and settings

| Method | Path | Request model | Purpose |
|---|---|---|---|
| GET | `/` | — | Serve the browser workbench. |
| GET | `/healthz` | — | Readiness, build/version, Primer3, data-directory, and NCBI configuration summary. |
| GET | `/api/profiles` | — | List assay design profiles. |
| GET | `/api/ranking-profiles` | — | List ranking objectives/profiles. |
| GET | `/api/conditions` | — | Read process-wide reaction conditions. |
| POST | `/api/conditions` | `CondReq` | Change shared reaction conditions when deployment policy permits. |

### Oligo and assay design

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/qc` | `OligoReq` | Single-oligo thermodynamic and sequence QC. |
| POST | `/api/pair` | `PairReq` | Primer-pair QC. |
| POST | `/api/matrix` | `MatrixReq` | Pairwise oligo interaction matrix. |
| POST | `/api/design` | `DesignReq` | Design against supplied template/FASTA input. |
| POST | `/api/accessibility` | `AccessReq` | Accessibility analysis. |
| POST | `/api/batch_design` | `BatchReq` | Bounded batch design, at most eight templates per request. |
| POST | `/api/viewer_design` | `ViewerDesignReq` | Design candidates on a user-selected sequence with coordinates. |

### Retrieval and sequence context

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/fetch` | `FetchReq` | Fetch an accession or resolve a gene/organism query through NCBI. |
| POST | `/api/fetch_nuc` | `FetchNucReq` | Search and fetch nucleotide records. |
| POST | `/api/gene_lookup` | `GeneLookupReq` | Query NCBI Gene context. |
| POST | `/api/intron` | `IntronReq` | Resolve intron/exon context for a proposed amplicon. |
| POST | `/api/blast` | `BlastReq` | Run remote or permitted local BLAST for one oligo. |
| POST | `/api/scan_markers` | `MarkerReq` | Scan marker candidates using NCBI-derived sequence context. |
| POST | `/api/suggest_genes` | `MarkerReq` | Suggest marker genes dynamically. |
| POST | `/api/refgenes` | `RefGenesReq` | Analyze a supplied reference-gene table. |
| POST | `/api/isolate_genomes` | `IsolateGenomesReq` | Search genome accessions, capped at 200 returned records. |

NCBI-facing request models may accept `email` and `ncbi_key`. Prefer server-side `OLIGOFORGE_EMAIL` and secret `OLIGOFORGE_NCBI_KEY` instead of transmitting credentials in request bodies.

### Specificity, conservation, and multiplex analysis

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/pair_specificity` | `PairSpecReq` | Primer-specific BLAST results. |
| POST | `/api/conservation` | `ConsReq` | Match oligos across supplied target/off-target sequences. |
| POST | `/api/epcr` | `EpcrReq` | In-silico PCR using remote/local search or supplied FASTA. |
| POST | `/api/assay_specificity` | `AssaySpecReq` | Complete-product and optional probe-recognition analysis. |
| POST | `/api/isolate_check` | `IsolateCheckReq` | Check one assay against up to six accessions per request. |
| POST | `/api/multiplex` | `MultiplexReq` | Multiplex compatibility analysis. |
| POST | `/api/orthogonal-panel` | `OrthoPanelReq` | Select/certify bounded oligo panels under declared thermodynamic conditions. |

### Quantification, reporting, and export

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/copies` | `CopiesReq` | Copy-number and dilution-series calculation. |
| POST | `/api/lna_tm` | `LnaReq` | LNA-aware melting-temperature calculation. |
| POST | `/api/lna_suggest` | `LnaReq` | Suggest bounded LNA placements. |
| POST | `/api/standard_curve` | `StdCurveReq` | Fit a standard curve. |
| POST | `/api/cq` | `CqReq` | Cq-related calculations. |
| POST | `/api/expression` | `ExpressionReq` | Relative-expression analysis. |
| POST | `/api/melt` | `MeltReq` | Melt-curve analysis. |
| POST | `/api/validate` | `ValidateReq` | Assay-readiness checks against supplied observations. |
| POST | `/api/order_csv` | `OrderReq` | Produce oligo CSV and gBlock FASTA strings. |
| POST | `/api/report` | `ReportReq` | Produce a structured panel report. |
| POST | `/api/rdml` | `RdmlReq` | Produce an RDML 1.3 assay-definition export. |

### Local project and panel state

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/project/save` | `ProjectSaveReq` | Save a project when server storage is enabled. |
| GET, POST | `/api/project/list` | — | List server projects. |
| POST | `/api/project/load` | `ProjectNameReq` | Load a named project. |
| POST | `/api/project/delete` | `ProjectNameReq` | Delete a named project. |
| POST | `/api/panel/save` | `PanelSaveReq` | Save a named oligo panel. |
| GET | `/api/panel/list` | — | List server panels. |
| POST | `/api/panel/load` | `PanelLoadReq` | Load a named panel. |
| POST | `/api/factory_reset` | — | Delete JSON projects and panels in the two managed server directories. |

These routes are single-instance file storage, not authenticated multi-tenant records. Hosted mode disables them unless an operator explicitly overrides the default.

### Automatic design jobs

| Method | Path | Request model | Purpose |
|---|---|---|---|
| GET | `/api/autodesign/limits` | — | Report queue, timeout, retention, privacy, and stage limits. |
| POST | `/api/autodesign/jobs` | `AutoDesignReq` | Submit an in-memory design job; returns HTTP 202. |
| GET | `/api/autodesign/jobs/{job_id}` | — | Read a public job snapshot. |
| DELETE | `/api/autodesign/jobs/{job_id}` | — | Request cancellation. |
| POST | `/api/autodesign/jobs/{job_id}/retry-blast` | `BlastRetryReq` | Retry specificity after a completed primary design; returns HTTP 202. |
| POST | `/api/autodesign` | `AutoDesignReq` | Synchronous compatibility endpoint. |

`AutoDesignReq` requires `target_query`; notable fields are `profile`, `off_query`, `n_fetch` (1–30), `min_ident` (0–1), optional BLAST settings, `prefer_junction`, `nested`, and `objective`. The asynchronous endpoint accepts `Idempotency-Key`. Reusing a key with different inputs returns HTTP 409; a full queue returns HTTP 429 with `Retry-After`; missing, expired, or restart-lost jobs return HTTP 404.

The queue is deliberately single-process and non-durable. Defaults are eight queued jobs, 1,800-second terminal retention, a 240-second primary timeout, and a 360-second BLAST timeout. Jobs and their capability-style identifiers disappear on restart. Cancellation is observed at stage boundaries; an already-running native or network call may drain before the worker starts another scientific stage.

### Manual design and decision analysis

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/manual-design/analyze` | `ManualDesignReq` | Analyze a supplied primer/probe design. |
| POST | `/api/manual-design/redesign` | `RedesignReq` | Generate bounded redesign alternatives under locks/regions. |
| POST | `/api/assay-rescue` | `RescueReq` | Rank redesign candidates from supplied observed issues. |
| POST | `/api/manual-design/compare-edit` | `ManualEditCompareReq` | Compare baseline and edited designs. |
| POST | `/api/experimental-feedback/status` | `FeedbackReq` | Summarize feedback readiness. |
| POST | `/api/experimental-feedback/import` | `FeedbackImportReq` | Parse supplied feedback content. |
| POST | `/api/experimental-feedback/split` | `FeedbackSplitReq` | Deterministically split supplied feedback records. |
| POST | `/api/experimental-feedback/summary` | `FeedbackReq` | Summarize supplied feedback. |
| POST | `/api/design-runs/compare` | `RunCompareReq` | Compare two supplied design-run records. |

The `/api/assay-rescue` endpoint is a bounded redesign/ranking helper. It is not an implemented Assurance Repair subsystem and does not autonomously modify an operational assay.

### Validation Studio

| Method | Path | Request model | Purpose |
|---|---|---|---|
| POST | `/api/validation-studio/plan` | `ValidationPlanReq` | Create a deterministic, bounded candidate-comparison plan and plate CSV. |
| POST | `/api/validation-studio/interpret` | `ValidationResultsReq` | Parse and interpret a supplied results CSV against a supplied plan. |

The plan model includes `candidates`, `cases`, objective, reaction conditions, plate format, replicates, controls, acceptance criteria, model, maximum cases, seed, edge-well policy, and existing evidence. Validation Studio plans experiments and interprets supplied data; it does not conduct experiments or create analytical or clinical validity.

### Assurance

| Method | Path | Required or notable request fields | Output |
|---|---|---|---|
| POST | `/api/assurance/assaysbom` | `assay` | Deterministic `oligoforge-assaysbom/v1` record. |
| POST | `/api/assurance/snapshots` | `fasta`; optional `name`, `role`, `source` object, and delimited `metadata` text | Immutable `oligoforge-sequence-snapshot/v1` record. |
| POST | `/api/assurance/snapshots/delta` | `baseline`, `followup` | Deterministic set delta after hash verification. |
| POST | `/api/assurance/drift-scan` | AssaySBOM and baseline/current target snapshots; optional off-target snapshots and model | Bounded `oligoforge-drift-scan/v1` record. |
| POST | `/api/assurance/ofvr` | `drift_scan`; optional `issuance_year` | OligoForge-local vulnerability records. |
| POST | `/api/assurance/package` | `assaysbom`; optional artifact lists | Evidence package, hash verification result, and escaped HTML rendering. |

Assurance routes are local computations; they do not fetch new sequence data. The file-oriented CLI is the reference workflow for FASTA.GZ and metadata-table input:

```bash
python -m oligoforge.assurance_cli --help
```

Assurance semantics are deliberately bounded:

- AssaySBOM is a molecular/computational inventory, not proof of performance.
- Snapshots describe only supplied records and do not establish database completeness or population representativeness.
- DriftGuard compares supplied snapshots under a declared complete-product model. It does not emit a numeric probability of clinical failure or future evolution.
- OFVR identifiers are local records, unreviewed by default, and are not a recognized external vulnerability standard.
- Evidence package hashes detect modification; they are not digital signatures and do not establish author identity.
- The package request may carry caller-supplied opaque `repairs` artifacts. OligoForge 1.35.0 does not implement an Assurance Repair generator or workflow.
- Aegis, Repair, and FutureProof subsystems are not implemented.

## Environment variables

| Variable | Default | Effect |
|---|---:|---|
| `OLIGOFORGE_HOSTED` | `0` | Enables hosted hardening behavior. |
| `OLIGOFORGE_ALLOW_SERVER_STORAGE` | on locally; off hosted | Enables project/panel file storage. |
| `OLIGOFORGE_ALLOW_SHARED_CONDITIONS` | on locally; off hosted | Enables process-wide condition changes. |
| `OLIGOFORGE_DATA_PATH` | application/temp default | Data, projects/panels, and NCBI cache location. |
| `OLIGOFORGE_EMAIL` | unset | NCBI Entrez identity. |
| `OLIGOFORGE_NCBI_KEY` | unset | NCBI API key. |
| `OLIGOFORGE_NCBI_TIMEOUT` | 30 s | NCBI socket timeout, clamped to 3–300 s. |
| `OLIGOFORGE_NCBI_CACHE` | `1` | Enables NCBI response caching. |
| `OLIGOFORGE_NCBI_CACHE_TTL` | 604,800 s | Cache retention. |
| `OLIGOFORGE_NCBI_RETRIES` | `3` | Retry count, bounded to 1–5. |
| `OLIGOFORGE_JOB_QUEUE` | `8` | Automatic-design queue capacity. |
| `OLIGOFORGE_JOB_TTL_SECONDS` | `1800` | Terminal job retention. |
| `OLIGOFORGE_DESIGN_TIMEOUT_SECONDS` | `240` | Primary design timeout. |
| `OLIGOFORGE_BLAST_TIMEOUT_SECONDS` | `360` | Optional specificity-stage timeout. |
| `OLIGOFORGE_MAX_REQUEST_BYTES` | 5 MiB | Declared request-body ceiling. |
| `OLIGOFORGE_LOG_LEVEL` | `INFO` | Application log level. |
| `PORT` | `8111` in Docker | Uvicorn listen port. |

## Implementation sources

Repository sources reviewed for this contract: `app.py`, `oligoforge/jobs.py`, `oligoforge/assurance/`, `oligoforge/assurance_cli.py`, `oligoforge/validation_studio.py`, `schemas/*.schema.json`, `Dockerfile`, `render.yaml`, and `SECURITY.md`.

Official external sources, all accessed 2026-07-15:

- FastAPI, automatic OpenAPI documentation: <https://fastapi.tiangolo.com/features/#automatic-docs>
- NCBI, E-utilities usage guidance: <https://www.ncbi.nlm.nih.gov/books/NBK25497/>
- NCBI, website and data usage policies: <https://www.ncbi.nlm.nih.gov/home/about/policies/>
