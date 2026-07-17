# Regulatory evidence mapping

## Status and intended use of this mapping

This mapping is a planning aid for organizations evaluating OligoForge 1.37.0 in an assay-development or quality workflow. It is not legal advice, a regulatory submission, a declaration of conformity, a quality-system procedure, or a determination that OligoForge is a medical device or is suitable for a particular regulated use.

OligoForge has not been shown in this repository to be FDA cleared, approved, authorized, listed, or registered; CE marked under the EU IVDR; certified to ISO 13485, ISO 14971, IEC 62304, ISO 15189, or another standard; or validated as a 21 CFR Part 11 electronic-record/signature system. Product classification and obligations depend on intended purpose, claims, jurisdiction, deployment, users, and how software output is used.

The mapping uses three contribution levels:

- **Supporting** — an implemented artifact may be placed in an organization's controlled evidence set after review.
- **Partial** — the artifact addresses only a narrow input or traceability need and substantial external evidence is required.
- **Not provided** — the release does not implement the required evidence or control.

## Artifact-to-evidence crosswalk

| OligoForge artifact | Potential controlled use | Contribution | What it does not establish |
|---|---|---|---|
| AssaySBOM (`oligoforge-assaysbom/v1`) | Configuration identification for primers, probes, declared chemistry/conditions, intended groups, status, review state, evidence links, and software versions | Supporting for an assay design record after independent review | Intended-use approval, manufacturing specifications, supplier qualification, analytical/clinical performance, or a software SBOM |
| Sequence snapshot (`oligoforge-sequence-snapshot/v1`) | Frozen input corpus, accession/source/retrieval fields, accepted/rejected ledger, exact deduplication, and content identity | Supporting for computational-input traceability if source metadata is completed and checked | Database completeness, population prevalence, independence, geographic/temporal representativeness, consent, or source licensing clearance |
| Snapshot delta | Exact added/removed/unchanged sequence comparison between two hash-valid snapshots | Supporting change input | Biological significance, prevalence, causality, or an approved change decision |
| DriftGuard scan (`oligoforge-drift-scan/v1`) | Declared-model evaluation and reason-coded review trigger for new supplied sequences | Partial risk/change-assessment input | Probability of failure or harm, analytical sensitivity/specificity, clinical performance, surveillance completeness, or future stability |
| OFVR (`oligoforge-ofvr/v1`) | Local issue record linking modeled observation, affected components, snapshots, model, and limitations | Partial issue-triage input | External vulnerability registration, severity consensus, regulatory reportability, complaint/CAPA disposition, or validated repair |
| Validation Studio plan (`oligoforge-validation-plan/v1`) | Deterministic candidate/case/control/replicate/plate plan and acceptance-criteria container | Supporting protocol draft | Executed testing, specimen suitability, instrument qualification, raw-data integrity, statistical adequacy, approvals, or analytical/clinical validity |
| Validation Studio interpretation | Calculation over a caller-supplied CSV and plan | Partial analysis record | Source-data authenticity, protocol adherence, investigation of deviations, independent statistical validation, or release authorization |
| Assurance evidence package | Artifact manifest, SHA-256 integrity checks, machine-readable bundle, and escaped HTML rendering | Supporting compilation/integrity record | Digital signature, trusted timestamp, author identity, electronic approval, access control, audit trail, archival validation, or regulatory acceptance |
| `SECURITY.md` and `ASSURANCE_THREAT_MODEL.md` | Security assumptions, implemented hardening, protected assets, and residual-risk inputs | Partial cybersecurity/risk input | Full security architecture, software SBOM, penetration testing, coordinated vulnerability process, patch SLA, monitoring, incident response, or secure-product lifecycle evidence |
| `API.md`, schemas, tests, and release records | Interface description, data contracts, verification inputs, and configuration traceability | Partial software documentation input | Complete software requirements, architecture, safety classification, hazard analysis, trace matrix, unresolved-anomaly assessment, independent V&V, or submission-ready documentation |

The evidence package accepts caller-supplied opaque `repairs` records. That field is only a packaging slot. OligoForge 1.37.0 does not implement an Assurance Repair subsystem. Aegis and FutureProof are also not implemented.

## United States mapping

### FDA Quality Management System Regulation

FDA's Quality Management System Regulation (QMSR) became effective 2026-02-02 and incorporates ISO 13485:2016 by reference, with FDA-specific provisions. Applicability and record requirements must be assessed within the manufacturer's actual quality system.

| Regulatory area | OligoForge contribution | Required external evidence or control |
|---|---|---|
| Design/development configuration and change inputs | AssaySBOM, versioned schemas, immutable snapshots, deltas, DriftGuard reason records, validation plans, and package manifests can support traceability. | Approved design/development procedures, user needs, intended use, design inputs/outputs, reviews, verification/validation, transfer, change authorization, competence, record retention, and management controls. |
| Risk-based decisions | DriftGuard/OFVR can identify bounded computational observations for risk-file review. | ISO 14971-aligned risk-management plan/file, hazard analysis, probability/severity rationale, risk controls, verification of controls, benefit-risk evaluation, production/postproduction feedback, and accountable approval. |
| Production and quality-system software assurance | Tests, deterministic artifacts, schemas, documented limits, and package verification can be inputs to a risk-based assurance record. | Intended-use definition for the software, process-risk analysis, assurance plan, objective evidence in the deployed configuration, deviations, approvals, supplier/cloud qualification, change control, periodic review, and validated backup/recovery where needed. |
| Records | JSON artifacts and SHA-256 checks can show that content changed. | Controlled record ownership, contemporaneous review/approval, retention, retrieval, backups, access control, auditability, and any FDA inspection requirements. Hashes alone are insufficient. |

Primary sources: FDA QMSR page and FAQ listed below.

### 21 CFR Part 11 electronic records and signatures

Part 11 applicability is an organization- and record-specific determination. Section 11.10 describes controls for closed systems, including validation, accurate and complete copies, record protection and retrieval, access limitation, secure time-stamped audit trails, authority/device checks, training/accountability, documentation controls, and signature-related controls.

OligoForge supplies deterministic export and hash verification only. It has no built-in authentication, role-based authorization, electronic signatures, signature/record binding, secure computer-generated time-stamped audit trail, trusted time source, training control, or validated long-term archive. Therefore:

- do not describe an Assurance package as “Part 11 compliant”;
- if OligoForge output becomes a regulated electronic record, place it in a validated record system with approved procedures and preserve the original artifact/hash relationship; and
- document whether a paper record, predicate-rule electronic record, or signed record in another system is the official record.

### IVD labeling under 21 CFR 809.10

AssaySBOM fields and OligoForge reports may provide controlled inputs for identity, reagent/component, instruction, warning, and performance-document drafting. They do not produce compliant IVD labeling. The responsible organization must establish the applicable intended use, indications, specimen, method, reagent/instrument details, storage/stability, warnings/limitations, procedural instructions, expected values, performance characteristics, and labeling approvals. Computational model output must not be substituted for required performance characteristics.

### Device software documentation

FDA's June 2023 guidance describes recommended documentation for premarket submissions for device software functions. Depending on applicability and documentation level, an organization may need software description, system/context and architecture information, requirements/specifications, risk management, cybersecurity, verification/validation, revision history, and unresolved-anomaly information.

`API.md`, schemas, deterministic artifact formats, model/version fields, tests, and release records can support portions of interface, configuration, and verification documentation. The repository does not present a complete submission package, device description, architecture set, software requirements specification, safety classification, hazard trace, V&V report for the intended device use, unresolved-anomaly assessment, or regulatory rationale.

### Cybersecurity and software bill of materials

FDA's February 2026 cybersecurity guidance and IMDRF SBOM principles concern security risk management and software-component transparency across the product lifecycle. `SECURITY.md`, the Assurance threat model, hosted-mode restrictions, request bounds, and sanitized errors are useful inputs. They are not a full cybersecurity management plan or evidence set.

AssaySBOM is a molecular/computational assay inventory. It must not be submitted or labeled as the software SBOM for OligoForge or a medical device. A separate software SBOM must cover the exact released software and transitive/native dependencies, with vulnerability monitoring and remediation processes appropriate to the product.

## European Union IVDR mapping

The consolidated Regulation (EU) 2017/746 text linked below was current through 2025-01-10 on EUR-Lex when accessed. Transition provisions, delegated/implementing acts, MDCG guidance, national requirements, and later amendments must be assessed separately.

| IVDR area | OligoForge contribution | Required external evidence or control |
|---|---|---|
| Article 10 manufacturer obligations and QMS | Versioned assay/configuration artifacts, change inputs, validation-plan drafts, and evidence packaging may support controlled records. | Legal manufacturer responsibilities, conformity route, QMS, regulatory strategy, PRRC where applicable, supplier/manufacturing controls, risk management, performance evaluation, registration, vigilance, CAPA, and lifecycle procedures. |
| Annex I general safety and performance requirements | Assay component inventory and bounded target/off-target model results may inform design/risk analysis. | Applicable-requirements checklist, objective evidence for every claimed requirement, risk controls, chemical/biological/physical safety, usability, labeling, performance, stability, and review/approval. |
| Annex II device description, design/manufacturing, and technical documentation | AssaySBOM and API/schemas can support part of device/design configuration; package manifests can organize referenced records. | Intended purpose, classification/rationale, variants/accessories, prior generations, design/manufacturing information, GSPR evidence, benefit-risk/risk file, complete product verification/validation, labeling, and submission structure. |
| Article 56 and Annex XIII performance evaluation | Validation Studio can draft a comparison experiment and hold declared criteria; snapshots/drift scans may identify computational challenge cases. | Performance-evaluation plan/report and continuous updates; scientific validity; executed analytical-performance studies; clinical-performance evidence; state of the art; statistical methods; specimen and population justification; literature/data appraisal; and qualified approval. |
| Annex II analytical and clinical performance documentation | Modeled complete-product/probe observations may be exploratory computational evidence. | Analytical sensitivity/specificity, trueness/bias, precision, interference/cross-reactivity, measuring range, cut-off, specimen handling, stability, clinical sensitivity/specificity, predictive values, likelihood ratios, expected values, clinical utility/association, and all applicable study records. |
| Annex III post-market surveillance | A reviewed snapshot/delta/DriftGuard/OFVR set can be one input to signal assessment when an organization supplies a controlled corpus. | Approved PMS plan/system, data sources and schedules, complaint/vigilance integration, trend/reporting rules, signal validation, benefit-risk and performance-evaluation updates, CAPA, reporting, and effectiveness checks. OligoForge performs no automated surveillance or retrieval. |

“Stable” in DriftGuard is limited to supplied snapshots and the declared model. It cannot support an IVDR claim that performance will remain stable across future variants or populations.

## Consensus standards and reporting practices

| Source | Possible use of OligoForge artifacts | Gap that remains |
|---|---|---|
| ISO 13485:2016 | Configuration, change, verification, and record inputs | An implemented and audited QMS, organization-specific procedures, responsibilities, review, and objective evidence |
| ISO 14971:2019 | Hazard/change inputs and links to affected components | Complete lifecycle risk management, estimates/evaluations, controls, residual/overall risk, benefit-risk, and production/postproduction review |
| IEC 62304:2006+A1:2015 | Interface/schema/version/test inputs | Software safety classification, lifecycle plans/processes, requirements/architecture, unit/integration/system V&V, problem resolution, configuration/change management, and maintenance evidence |
| ISO 15189:2022 | Method-selection and verification/validation planning inputs for a medical laboratory | Laboratory accreditation scope, competence, equipment/reagent control, examination procedures, QC/EQA, uncertainty where applicable, result reporting, and authorized validation/verification |
| MIQE | Some assay identity, oligo, conditions, standard-curve, and report fields | Complete MIQE information and the underlying measured, reviewed experimental evidence |
| RDML | RDML 1.3 assay-definition export | Completeness/correctness of all required experimental metadata and no independent conformance certification |
| NIST SP 800-218 | Repository tests, dependency pins, threat model, and release controls as practice inputs | Organization-wide secure-development practices, provenance, review, build integrity, vulnerability response, and measured implementation |

Clause-level conformity cannot be assessed from public ISO/IEC landing pages. Use licensed copies of the standards, the exact applicable editions, and qualified regulatory/quality review.

## Claims language

| Supported wording after record-specific review | Wording not supported by this release |
|---|---|
| “The attached snapshot is a hash-verified record of the supplied FASTA corpus.” | “The snapshot is complete or representative of the circulating population.” |
| “The scan found no modeled concern in the supplied snapshots under the declared parameters.” | “The assay is future-proof,” “will not drop out,” or “is clinically safe/effective.” |
| “The OFVR is an unreviewed OligoForge-local reason record.” | “The issue has an externally recognized vulnerability identifier or regulatory severity.” |
| “The evidence package verifies its manifest and content hashes.” | “The package is signed, Part 11 compliant, regulator approved, or submission ready.” |
| “Validation Studio generated a bounded experiment plan.” | “OligoForge analytically or clinically validated the assay.” |
| “AssaySBOM records the declared molecular configuration.” | “AssaySBOM is an SPDX/CycloneDX/IMDRF software SBOM or proves device compliance.” |
| “A ranked redesign candidate was generated by the assay-rescue helper.” | “An Assurance Repair was validated, approved, or deployed.” |

## Evidence and controls that must come from outside OligoForge

Before regulated reliance, the responsible organization should determine and document at least:

1. intended purpose, claims, users, patient/specimen population, jurisdiction, product classification, and regulatory route;
2. controlled quality, software-lifecycle, risk-management, supplier, change, record, training, complaint, CAPA, and release procedures;
3. the validated deployed configuration, including OS/container, exact dependencies, external services, databases, models, settings, infrastructure, access, backup/recovery, and cybersecurity controls;
4. independent verification of source data, accessions, licenses, corpus selection, transformations, model applicability, and search completeness;
5. executed analytical-performance studies with prespecified protocols, raw data, deviations, statistics, acceptance decisions, and authorized approvals;
6. scientific-validity and clinical-performance evidence appropriate to the claimed purpose;
7. manufacturing/lot, reagent, instrument, specimen, stability, transport, interference, cross-reactivity, and quality-control evidence;
8. human-factors/usability, labeling, warnings, limitations, and operator-training evidence where applicable;
9. a software SBOM, security risk assessment, vulnerability monitoring/remediation, penetration/security testing, incident response, and supported-lifecycle plan;
10. validated electronic-record/signature controls if Part 11 or analogous requirements apply; and
11. a postmarket/performance follow-up and surveillance system that validates signals and feeds risk, performance evaluation, labeling, and CAPA.

## Official and primary sources

All sources were accessed 2026-07-15.

### United States

- FDA, Quality Management System Regulation (QMSR): <https://www.fda.gov/medical-devices/postmarket-requirements-devices/quality-management-system-regulation-qmsr>
- FDA, QMSR Frequently Asked Questions: <https://www.fda.gov/medical-devices/quality-management-system-regulation-qmsr/quality-management-system-regulation-frequently-asked-questions>
- eCFR, 21 CFR Part 11, Electronic Records; Electronic Signatures: <https://www.ecfr.gov/current/title-21/chapter-I/subchapter-A/part-11>
- eCFR, 21 CFR 11.10, Controls for closed systems: <https://www.ecfr.gov/current/title-21/chapter-I/subchapter-A/part-11/subpart-B/section-11.10>
- FDA, Part 11 Scope and Application guidance: <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/part-11-electronic-records-electronic-signatures-scope-and-application>
- eCFR, 21 CFR 809.10, Labeling for in vitro diagnostic products: <https://www.ecfr.gov/current/title-21/chapter-I/subchapter-H/part-809/subpart-B/section-809.10>
- FDA, *Content of Premarket Submissions for Device Software Functions* (June 2023): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/content-premarket-submissions-device-software-functions>
- FDA, *Computer Software Assurance for Production and Quality Management System Software* (February 2026): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/computer-software-assurance-production-and-quality-management-system-software>
- FDA, *Cybersecurity in Medical Devices: Quality System Considerations and Content of Premarket Submissions* (February 2026): <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/cybersecurity-medical-devices-quality-management-system-considerations-and-content-premarket>

### European Union

- EUR-Lex, consolidated Regulation (EU) 2017/746 on in vitro diagnostic medical devices, version dated 2025-01-10: <https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX%3A02017R0746-20250110>

### Standards and technical practice

- ISO 13485:2016 official page, reviewed and confirmed 2025: <https://www.iso.org/standard/59752.html>
- ISO 14971:2019 official page, reviewed and confirmed 2025: <https://www.iso.org/standard/72704.html>
- IEC 62304:2006+A1:2015 consolidated version: <https://webstore.iec.ch/en/publication/22794>
- ISO 15189:2022 official page: <https://www.iso.org/standard/76677.html>
- NIST SP 800-218, Secure Software Development Framework: <https://csrc.nist.gov/pubs/sp/800/218/final>
- IMDRF, *Principles and Practices for Software Bill of Materials for Medical Device Cybersecurity* (N73): <https://www.imdrf.org/documents/principles-and-practices-software-bill-materials-medical-device-cybersecurity>
- Bustin et al., MIQE primary paper: <https://pubmed.ncbi.nlm.nih.gov/19246619/>
- RDML Consortium, MIQE/RDML information: <https://rdml.org/miqe.html>
