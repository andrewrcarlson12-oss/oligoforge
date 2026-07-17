# Data licensing and attribution

This document records the data and direct-dependency obligations relevant to the implemented OligoForge 1.37.0 release. It is an engineering inventory, not legal advice. A deployer or distributor must review the terms that apply to its jurisdiction, data sources, users, and distribution format.

The repository's MIT `LICENSE` covers OligoForge-authored source code. It does not relicense third-party libraries, sequence records, publications, user inputs, local BLAST databases, names, logos, or generated content derived from them.

## Implemented data flows

| Source or input | Implemented use | Rights and attribution handling |
|---|---|---|
| User-supplied FASTA, FASTA.GZ, JSON, CSV/TSV, assay definitions, and experimental records | Design, QC, in-silico analysis, Validation Studio, and Assurance artifacts | The user/deployer must have authority to process and share the data. Preserve source, retrieval date, accession/version, and any study or consent restrictions. Do not place controlled or identifiable data in a public service. |
| NCBI Nucleotide/GenBank/RefSeq, Gene, Taxonomy, Entrez E-utilities, and remote BLAST | Accession/gene retrieval, marker discovery, automatic design, specificity, and isolate analysis | Follow NCBI policies and database-specific notices. Record stable accession.version values, database, retrieval date, query, and transformations. Identify OligoForge as an independent client; do not imply NCBI endorsement. |
| Operator-supplied local BLAST database | Local specificity and in-silico PCR | The operator is responsible for the database license, updates, access control, attribution, and provenance. OligoForge does not inspect or grant database rights. |
| Committed benchmark/test fixtures | Offline deterministic tests and comparisons | Fixture manifests record NCBI accessions and, where applicable, publication DOI/source. Preserve those manifests with redistributed fixtures and re-check source/publication rights before repurposing them outside testing. |
| Generated AssaySBOM, snapshot, delta, DriftGuard, OFVR, validation-plan, and evidence-package records | Reproducible local evidence artifacts | These may embed sequences, metadata, citations, and source identifiers. Their hashes do not change the ownership or confidentiality of embedded material. Review and redact metadata before sharing. |

No Ensembl, UniProt, GISAID, commercial sequence database, or automated epidemiological feed integration is implemented in this release. Their terms are therefore not asserted here.

## NCBI and GenBank handling

NCBI states that information created by or for the US government on its sites is generally public domain, while some content may be supplied or copyrighted by third parties. NCBI requests acknowledgment as the source, and individual database records may carry additional notices. GenBank states that NCBI places no restrictions on use or distribution of GenBank data, but submitters may claim patent, copyright, or other intellectual-property rights and NCBI cannot grant permission for unrestricted use. Those qualifications must travel with a release or evidence package; “NCBI data” is not a blanket rights warranty.

For an NCBI-derived sequence or corpus, retain at least:

- database and record accession.version;
- record title/organism and sequence role (target, off-target, or background);
- retrieval timestamp and the query or selection rule;
- source URL or API endpoint and any record-level restrictions;
- raw-input hash, accepted/rejected disposition, deduplication method, and subsequent transformations;
- citation or submitter attribution when supplied.

Suggested acknowledgment: “Sequence records were obtained from the NCBI [database name], National Library of Medicine; accession versions and retrieval dates are recorded in the accompanying provenance manifest.” This wording does not imply NCBI endorsement.

NCBI's E-utilities guidance asks software using E-utilities to display a disclaimer/copyright notice, identify its tool and contact email, and stay within rate limits. The published defaults are no more than three requests per second without an API key and ten requests per second with a key unless NCBI grants a higher rate. Set `OLIGOFORGE_EMAIL`; keep `OLIGOFORGE_NCBI_KEY` in deployment secrets. Respect retry/backoff and cache policies, and do not share keys in project files, logs, screenshots, or evidence packages.

Recommended in-product notice for an NCBI-enabled deployment:

> This product uses NCBI services but is not endorsed or certified by NCBI. Users remain responsible for verifying records, versions, database notices, and permitted use.

## Publications and assay sequences

Facts such as an accession or DOI are not a license to reproduce a paper, table, figure, or complete supplementary dataset. Published primer/probe sequences and assay descriptions may be covered by the article's license, database terms, patent rights, or other restrictions. Store the DOI, bibliographic source, location within the source, and the article/data license. Quote only what is necessary and permitted. A publication reporting an assay does not establish that OligoForge may market, manufacture, or clinically use it.

The benchmark corpus includes accessions and literature provenance in `tests/benchmark/corpus_provenance.csv`, `tests/benchmark/bench_corpus_published.json`, and `tests/benchmark/genome_fixtures/provenance.json`. Those files are provenance aids, not a legal clearance opinion.

## Direct software dependencies

The versions below are pinned in `requirements.txt` or declared in `package.json`. License metadata is recorded from the named package pages as accessed 2026-07-15. A distributor must retain the actual license and notice files from the installed artifacts; package-index metadata alone is not a complete notices bundle.

| Dependency | Release | Recorded license metadata | Distribution action |
|---|---:|---|---|
| primer3-py | 2.3.0 | GPLv2 | Preserve GPL notices and corresponding-source obligations applicable to the distributed combination. Obtain legal review before distributing a bundled executable/container under an incompatible licensing theory. |
| Biopython | 1.87 | `LicenseRef-Biopython-License-Agreement` | Include the license included with Biopython; do not replace it with a guessed SPDX identifier. |
| ViennaRNA | 2.7.2 | “Disclaimer and Copyright”; package page directs users to `COPYING` | Include the exact `COPYING`/copyright material from the installed release and review component-specific notices. |
| FastAPI | 0.139.0 | MIT | Preserve the MIT copyright and permission notice. |
| Uvicorn | 0.51.0 | BSD-3-Clause | Preserve copyright, conditions, and disclaimer. |
| HTTPX | 0.28.1 | BSD-3-Clause | Preserve copyright, conditions, and disclaimer. |
| jsdom (development/test only) | 24.1.x resolved as 24.1.3 in `package-lock.json` | MIT | Preserve the MIT notice when redistributing the test dependency or a bundle containing it. It is not a Python runtime dependency. |

Transitive dependencies also require review. Generate a notices/SBOM inventory from the exact built environment and verify it against the container or executable that will be shipped. Do not infer the installed graph solely from top-level pins.

## Confidentiality and controlled data

- Assurance commands are offline, but their outputs can contain complete sequences and free-text metadata.
- The web API has no built-in authentication. Use a local single-user deployment for sensitive work unless a controlled deployment supplies identity, authorization, TLS, tenant isolation, audit, retention, deletion, and incident-response controls.
- Do not submit protected human genomic data, patient identifiers, confidential assay designs, export-controlled data, or restricted pathogen data unless the deployment and intended use have been formally approved.
- Remote BLAST sends sequences to NCBI. Automatic design sends only the winning primer pair when the caller explicitly requests remote BLAST. A local BLAST workflow avoids upload but shifts database licensing and security to the operator.
- Content hashes provide integrity checking, not anonymization. Short oligos and hashed full sequences may still be sensitive or linkable.

## Release and evidence-package checklist

1. Freeze exact dependency versions and preserve their license/notice files.
2. Inventory transitive dependencies and any native code bundled in wheels or executables.
3. For each external corpus, record owner, license/terms URL, accession/version, query, retrieval date, and transformation.
4. Preserve rejected-record and deduplication ledgers; do not describe a filtered corpus as the source database.
5. Verify publication/table reuse rights separately from sequence-database rights.
6. Remove credentials, local paths, personal data, and unnecessary metadata before export.
7. Document whether any remote service received sequences.
8. Re-check time-varying policies at release and at each corpus refresh.

## Official and primary sources

All sources were accessed 2026-07-15.

### Data and service terms

- NCBI, Website and Data Usage Policies: <https://www.ncbi.nlm.nih.gov/home/about/policies/>
- NCBI, GenBank Overview and Data Usage: <https://www.ncbi.nlm.nih.gov/genbank/>
- NCBI, Entrez Programming Utilities Help, Usage Guidelines and Requirements: <https://www.ncbi.nlm.nih.gov/books/NBK25497/>

### Direct dependency metadata

- primer3-py 2.3.0: <https://pypi.org/project/primer3-py/2.3.0/>
- Biopython 1.87: <https://pypi.org/project/biopython/1.87/>
- ViennaRNA 2.7.2: <https://pypi.org/project/ViennaRNA/2.7.2/>
- FastAPI 0.139.0: <https://pypi.org/project/fastapi/0.139.0/>
- Uvicorn 0.51.0: <https://pypi.org/project/uvicorn/0.51.0/>
- HTTPX 0.28.1: <https://pypi.org/project/httpx/0.28.1/>
- jsdom package metadata: <https://registry.npmjs.org/jsdom>
