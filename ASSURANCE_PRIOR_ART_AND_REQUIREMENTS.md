# Assurance prior art and requirements

## Purpose and boundary

This document defines the implemented requirements baseline for OligoForge Assurance 1.35.0 and compares its artifacts with relevant public standards and reporting practices. It is not a patent prior-art search, freedom-to-operate opinion, standards certification, regulatory determination, or claim of novelty.

Only the following Assurance capabilities are implemented:

1. deterministic AssaySBOM generation;
2. bounded offline sequence snapshots and snapshot deltas;
3. DriftGuard comparison of supplied baseline and follow-up snapshots under a declared complete-product model;
4. OligoForge-local OFVR generation from reason records;
5. Validation Studio planning and interpretation of caller-supplied results; and
6. deterministic evidence-package assembly and integrity verification.

Aegis multi-edit search, an Assurance Repair system, and FutureProof design are not implemented. A caller may place opaque records in an evidence package's `repairs` list, but OligoForge does not generate, validate, approve, or apply a repair.

## Implemented requirements baseline

“Implemented” below means present in the reviewed code and schemas. It does not mean analytically or clinically validated.

| ID | Requirement | Implemented behavior | Verification source |
|---|---|---|---|
| ASR-001 | Deterministic serialization and identity | Artifacts use canonical value hashing and SHA-256-derived identifiers; fixed normalized inputs, declared versions, and parameters produce unchanged artifact content. OFVR callers must fix `issuance_year` when reproducibility across calendar years is required. | `oligoforge/provenance.py`, `tests/test_assurance.py` |
| ASR-002 | Explicit schema and artifact versions | AssaySBOM, snapshot, drift scan, OFVR, validation plan, and evidence package declare schema/version fields; JSON schemas are published under `schemas/`. | `schemas/*.schema.json`, `oligoforge/assurance/` |
| ASR-003 | Bounded assay inventory | An AssaySBOM accepts 1–64 assays, at most 512 components total, and at most 200 bare nucleotide bases per component. Every assay requires forward and reverse primers. | `oligoforge/assurance/assaysbom.py` |
| ASR-004 | Preserve molecular ordering information | Each component stores a normalized bare sequence separately from its `order_sequence`; modification notation is preserved rather than silently interpreted as bases. | `oligoforge/assurance/assaysbom.py` |
| ASR-005 | Record declared operational context | Portfolio/assay status, review state, chemistry/profile, target/background groups, reaction conditions, interpretation rules, existing evidence, and software versions have explicit fields. Most are caller declarations, not independently verified facts. | `schemas/assaysbom.schema.json`, `oligoforge/assurance/assaysbom.py` |
| ASR-006 | Bounded snapshot ingestion | Snapshot input is limited to 20 MB compressed, 100 MB uncompressed, 10,000 records, 500,000 nt per record, and 10,000,000 total bases. The nucleotide alphabet is validated. | `oligoforge/assurance/snapshots.py` |
| ASR-007 | Transparent corpus disposition | Accepted, rejected, unique, and exact-duplicate records remain separately visible; rejection reasons and raw/unique/group counts are retained. No sampling representativeness is inferred. | `schemas/sequence_snapshot.schema.json`, `oligoforge/assurance/snapshots.py` |
| ASR-008 | Immutable snapshot identity | A snapshot records input/content hashes, role, source/retrieval mappings, and a content-addressed ID. Hash validation detects content changes but does not authenticate an author or source. | `oligoforge/assurance/snapshots.py` |
| ASR-009 | Reproducible baseline/follow-up delta | A delta is calculated only from hash-valid snapshots of the same role and reports added, removed, and unchanged exact sequence hashes. | `oligoforge/assurance/snapshots.py` |
| ASR-010 | Complete-product modeled comparison | DriftGuard evaluates forward and reverse primer coherence and optional probe recognition using the existing isolate amplification model and declared mismatch/product/probe parameters. | `oligoforge/assurance/driftguard.py`, `oligoforge/isolates.py` |
| ASR-011 | Reason-coded, bounded states | DriftGuard emits record-level reasons and one of: Scan incomplete, Evidence insufficient, Possible signal-generating off-target, Possible target dropout, Stable with new variation, or Stable. The state is bounded to supplied snapshots. | `schemas/drift_scan.schema.json`, `oligoforge/assurance/driftguard.py` |
| ASR-012 | No false risk precision | DriftGuard explicitly states that it calculates no numeric biological-risk probability and that experimental confirmation is required. | `oligoforge/assurance/driftguard.py`, `ASSURANCE_VALIDATION_LIMITS.md` |
| ASR-013 | Local vulnerability records | OFVR records are generated from DriftGuard reason records, include model/snapshot links and limitations, start as `unreviewed`, and state that they are not a recognized external vulnerability standard. An explicit `issuance_year` makes IDs reproducible; otherwise the current UTC year is used. | `schemas/ofvr.schema.json`, `oligoforge/assurance/ofvr.py` |
| ASR-014 | Tamper-evident evidence packaging | Evidence packages contain artifact hashes, a manifest, a package hash/ID, verification results, and escaped HTML rendering. Hashes detect modification but are not digital signatures. | `schemas/assurance_evidence_package.schema.json`, `oligoforge/assurance/evidence_package.py` |
| ASR-015 | Bounded experiment planning | Validation Studio lays out declared candidates, cases, replicates, controls, acceptance criteria, model, and plate policy; it interprets a caller-supplied results CSV against the plan. | `schemas/validation_plan.schema.json`, `oligoforge/validation_studio.py` |
| ASR-016 | Offline file workflow | The Assurance CLI has explicit file arguments and performs no network retrieval for `build-assaysbom`, `build-snapshot`, `delta`, `drift-scan`, `ofvr`, or `package`. | `oligoforge/assurance_cli.py`, `ASSURANCE_CLI.md` |
| ASR-017 | Honest scope statements | Artifacts and reports distinguish computational evidence from analytical, clinical, population, regulatory, and future-performance evidence. | `ASSURANCE_VALIDATION_LIMITS.md`, `ASSURANCE_THREAT_MODEL.md`, artifact builders |

## Prior-art and standards comparison

### Software bills of materials

NTIA's minimum elements, SPDX, CycloneDX, and the IMDRF medical-device SBOM principles address software component transparency and supply-chain/cybersecurity uses. OligoForge borrows the useful inventory, versioning, supplier/context, dependency/provenance, and machine-readable-manifest ideas for a molecular assay portfolio.

AssaySBOM is nevertheless a separate OligoForge schema. It inventories primers, probes, declared chemistry, conditions, target groups, and evidence links. It is not an NTIA-complete software SBOM, an SPDX document, a CycloneDX BOM, or an IMDRF medical-device software SBOM. It does not replace a software SBOM for OligoForge or for a regulated device that embeds OligoForge.

### Provenance models

W3C PROV provides a general model for entities, activities, agents, derivations, and responsibility. OligoForge artifacts use content identifiers, source/retrieval mappings, baseline/follow-up relationships, model versions, and manifest links that can support later provenance mapping. OligoForge does not emit PROV-N, PROV-XML, PROV-O, or another W3C PROV serialization, and no PROV conformance is claimed. Source and retrieval fields remain caller-supplied and are not automatically corroborated.

### Quality, risk, and software lifecycle practices

ISO 13485, ISO 14971, IEC 62304, FDA's Quality Management System Regulation, FDA device-software guidance, FDA computer-software-assurance guidance, NIST SSDF, and FDA cybersecurity guidance establish or describe broader quality, risk, software-lifecycle, validation, and cybersecurity activities. Assurance artifacts can support configuration identification, traceability, change review, risk inputs, and verification records within a controlled process.

The implementation does not supply a quality management system, design controls, risk-management file, software safety classification, lifecycle plan, requirements approval, independent review, electronic signatures, access-controlled audit trail, complaint handling, CAPA, supplier controls, release authorization, postmarket surveillance program, or regulatory submission. A content hash is not a compliant record system by itself.

The ISO and IEC public pages cited below describe the standards but do not contain the full copyrighted requirements. Clause-level assessment must use authorized copies of the applicable editions and organization-specific procedures.

### MIQE and RDML

The MIQE paper describes minimum information needed to evaluate quantitative real-time PCR experiments, and RDML provides a machine-readable data format for qPCR information. OligoForge can export RDML 1.3 assay definitions and can record selected design and experimental context. That export is not evidence that every MIQE item is present, correctly measured, reviewed, or reported. OligoForge does not automatically establish sample provenance, extraction quality, inhibition, assay efficiency, dynamic range, limit of detection/quantitation, precision, clinical performance, or inter-laboratory reproducibility.

### Sequence surveillance and change assessment

DriftGuard applies an established software pattern—versioned baseline, immutable input snapshot, deterministic set delta, bounded model execution, reason records, and review state—to molecular sequence inputs. It performs no scheduled retrieval, epidemiological surveillance, prevalence estimation, phylogenetics, forecasting, or future-evolution prediction. “Stable” means no modeled concern in the supplied snapshots and declared model; it is not a guarantee of future assay performance.

## Explicit non-requirements for 1.35.0

The following are out of scope and must not be represented as present:

- Aegis multi-edit search or adversarial evolutionary search;
- automated Repair generation, laboratory bridging, approval, deployment, or rollback;
- FutureProof assay generation or a guarantee against future variants;
- automated database retrieval, scheduling, continuous monitoring, alert delivery, or corpus licensing clearance;
- numerical probability of dropout, false positivity, clinical harm, or future mutation;
- external registration or review of OFVR identifiers;
- conformance to SPDX, CycloneDX, W3C PROV, MIQE, FDA, IVDR, ISO, IEC, NIST, GxP, HIPAA, SOC 2, or ISO 27001;
- analytical sensitivity/specificity, precision, LOD/LOQ, matrix, inhibition, stability, clinical sensitivity/specificity, or clinical utility evidence;
- authenticated multi-tenant evidence storage, qualified electronic signatures, trusted timestamps, or immutable audit trails; and
- legal clearance of sequences, publications, software dependencies, patents, or third-party data.

## Change-control requirements for future work

A future feature must not reuse an existing schema version if it changes scientific meaning. It must define inputs, bounds, deterministic behavior, failure states, provenance, model version, limitations, and migration behavior; add schema and negative tests; and update the regulatory mapping without converting a planning relationship into a conformance claim.

Any future retrieval or monitoring function must additionally define source terms, authentication, rate limits, query/retrieval timestamps, record-version handling, completeness status, retries, cache/retention, privacy/egress, and reproducible corpus freezing. Any future repair function must separate proposal, experimental verification, expert approval, operational release, and rollback, with no automatic promotion from in-silico ranking to operational use.

## Official and primary sources

All sources were accessed 2026-07-15.

### Inventory, provenance, and secure development

- US Department of Commerce/NTIA, *The Minimum Elements for a Software Bill of Materials*: <https://www.ntia.gov/report/2021/minimum-elements-software-bill-materials-sbom>
- SPDX 3.0.1 specification: <https://spdx.github.io/spdx-spec/v3.0.1/>
- CycloneDX specification overview, current published version 1.7: <https://cyclonedx.org/specification/overview/>
- IMDRF, *Principles and Practices for Software Bill of Materials for Medical Device Cybersecurity* (N73): <https://www.imdrf.org/documents/principles-and-practices-software-bill-materials-medical-device-cybersecurity>
- W3C, PROV overview: <https://www.w3.org/TR/prov-overview/>
- NIST SP 800-218, Secure Software Development Framework: <https://csrc.nist.gov/pubs/sp/800/218/final>

### Quality, risk, software, and laboratory context

- FDA, Quality Management System Regulation, effective 2026-02-02: <https://www.fda.gov/medical-devices/postmarket-requirements-devices/quality-management-system-regulation-qmsr>
- FDA, *Content of Premarket Submissions for Device Software Functions* (June 2023): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/content-premarket-submissions-device-software-functions>
- FDA, *Computer Software Assurance for Production and Quality Management System Software* (February 2026): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/computer-software-assurance-production-and-quality-management-system-software>
- FDA, *Cybersecurity in Medical Devices: Quality System Considerations and Content of Premarket Submissions* (February 2026): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/cybersecurity-medical-devices-quality-management-system-considerations-and-content-premarket>
- ISO 13485:2016 official page, reviewed and confirmed 2025: <https://www.iso.org/standard/59752.html>
- ISO 14971:2019 official page, reviewed and confirmed 2025: <https://www.iso.org/standard/72704.html>
- IEC 62304:2006+A1:2015 consolidated version: <https://webstore.iec.ch/en/publication/22794>
- ISO 15189:2022 official page: <https://www.iso.org/standard/76677.html>

### qPCR reporting

- Bustin et al., *The MIQE Guidelines: Minimum Information for Publication of Quantitative Real-Time PCR Experiments* (primary paper): <https://pubmed.ncbi.nlm.nih.gov/19246619/>
- RDML Consortium, MIQE/RDML information: <https://rdml.org/miqe.html>
