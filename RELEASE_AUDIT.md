# OligoForge 1.34.0 Release Audit

## Scope

Ranking truth, complete-triplet retention, objective profiles, structured ranking, attrition provenance, explicit rank uncertainty, deterministic self-hashing manifests, Manual Design Studio, exact edit comparison, design-run comparison, constrained redesign, Assay Rescue, experimental-feedback dataset audit and context-local summaries, authoritative Batch Design, benchmark-integrity checks, browser integration, deployment hardening, and release documentation.

## Final gate results

- Python scientific/software programs: 44 passed, 0 failed.
- Browser/UI harnesses: 9 passed, 0 failed.
- Ranking-truth assertions: 41 passed, 0 failed.
- Evidence/provenance/API assertions: 36 passed, 0 failed.
- Decision-analysis/API assertions: 25 passed, 0 failed.
- Manual Design Studio browser assertions: 38 passed, 0 failed.
- Workbench/report provenance persistence: automatic, batch, direct, manual and rescue candidates retain rank manifests; HTML/CSV and RDML exports surface verified or altered state.
- Frozen synthetic/adversarial benchmark: passed; 11/11 Top-1 with Wilson 95% interval 74.1%–100%.
- Frozen Plasmodium/Haemoproteus biological regressions: passed.
- Clean pinned Python environment import and dependency check: passed; 63 routes.
- Python byte compilation: passed.
- npm high-severity audit: 0 known vulnerabilities.
- Feedback records preserve design run ID, application/ranker versions, manifest hash, and reaction-condition snapshot.
- Benchmark corpus validation rejects split overlap, target-group leakage, invalid expectations, and malformed identifiers.
- Final staged-copy, file-manifest, and archive-integrity checks are completed during packaging and recorded in `RELEASE_MANIFEST.json`.

## Release decision

The release is suitable as a computational ranking-truth, decision-analysis, and evidence-provenance pre-release. It improves the frozen synthetic/adversarial selection benchmark, fixes confirmed candidate-loss and false-confidence pathways, makes missing evidence and close ranks explicit, and provides auditable manual edits, run comparisons, and batch results. It is not a biologically validated assay platform and does not establish rank-1 wet-lab superiority.

## Deployment limitation

The actual Render service was not modified or verified because no GitHub repository URL, Render service URL, or deployment credentials were available in the supplied archive.
