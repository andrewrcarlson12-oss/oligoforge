# OligoForge 1.31.1 release audit

Audit date: 2026-07-15

## Automated release gates

- Python scientific/software test programs: **39 passed, 0 failed** (`python run_tests.py --python`).
- Browser/UI harnesses: **8 passed, 0 failed** (`python run_tests.py --node`).
- Inline browser handlers audited by the integration harness: **90**.
- Python byte-compilation: passed (`python -m compileall`).
- npm dependency audit: **0 known vulnerabilities** after lockfile update.
- Archive integrity: verify the distributed ZIP with `unzip -t` and the SHA-256 supplied alongside it.

## What the offline gate covers

Design and ranking, exact coordinates, thermodynamics, LNA/degenerate handling, exon-junction design,
conservation, ambiguity-aware specificity, offline PCR, Primer-BLAST comparison fixtures, multiplex
screening, orthogonal graph bounds, Cq/melt/standard-curve analysis, reference genes, MIQE-aligned
readiness records, RDML/order/report export, hosted isolation, request validation, concurrency,
performance fixtures, fuzzing, and browser workflows.

## Deliberately separate integration work

Live NCBI retrieval and remote BLAST are not part of the deterministic offline release gate. Enable
those checks with `OLIGOFORGE_LIVE_NCBI=1` and a real `OLIGOFORGE_EMAIL`; local BLAST additionally
requires a configured database. These external services can fail independently of OligoForge.

## Scientific boundary

Passing these gates demonstrates internal consistency against the included models and fixtures. It
does not validate a primer/probe assay experimentally. Read `VALIDATION_LIMITS.md` before using a
candidate in a publication, diagnostic workflow, or production laboratory method.
