# OligoForge 1.35.0

OligoForge is a computational primer, probe, and qPCR assay-development workbench. Version 1.35.0 is the **assay assurance release**: automatic design now runs as a bounded, cancellable staged job; Manual Design presents evidence as human-readable scientific decisions; Validation Studio creates reproducible candidate-comparison experiments; and the first offline Assurance workflow produces versioned AssaySBOMs, immutable sequence snapshots, deterministic deltas, complete-product DriftGuard scans, OligoForge Molecular Vulnerability Records, and self-verifying evidence packages.

## Core behavior

- Target-spanning, bounded joint search over complete forward-primer/reverse-primer/probe triplets.
- Diversity-preserving retention across regions, primer pairs, probe positions, product lengths, and trade-off profiles.
- True hard scientific constraints before objective-specific evidence and soft preferences.
- Coherent target coverage, paired all-site specificity, condition-envelope thermodynamics, junction evidence, and multiplex interactions can change final rank.
- Candidate attrition ledger, structured rank trace, Pareto fronts, truthful finalist categories, and concrete rank-reversal scenarios.
- Explicit near-equivalent/insufficient-evidence states rather than false numerical certainty.
- Deterministic self-hashing manifests with application/ranker/model versions, reaction conditions, inputs, constraints, candidate limits, database state, warnings, and fallbacks.
- Automatic, batch, sequence-viewer, manual-design, constrained-redesign, and rescue workflows use the authoritative scientific modules and structured ranker.
- Manual mapping reports every exact and allowed near-match, including non-extension-eligible 3′ mismatches, instead of silently hiding failed placements.
- Manual edits can be compared against the last analyzed baseline using complete recalculated evidence: resolved/new hard failures, target/off-target coverage, robustness, interactions, mapping, and rank preference.
- Complete design runs can be compared for candidate loss, winner changes, rank reversals, context changes, and unexplained non-reproducibility.
- Experimental feedback is validated, deduplicated, conflict-checked, locally summarized, and split by target group without silently changing the ranker. A learned reranker remains disabled without sufficient leakage-controlled evidence.
- Long automatic designs return a capability-style job ID promptly, expose real stage states, support cancellation and idempotent retry, and retain the completed primary design when optional BLAST is unavailable.
- Validation Studio selects a bounded, diverse set of supplied cases on which leading candidates disagree, interleaves candidates on deterministic 96/384-well layouts, exports injection-safe CSV, and interprets completed observations conservatively.
- Assurance freezes assay and sequence evidence with deterministic hashes and reports complete-product sequence concerns without converting them into an unvalidated probability or clinical-performance claim.

## Run locally

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
uvicorn app:app --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000`.

## Test

```bash
npm ci
python run_tests.py
python tools/check_directory_file_counts.py
python tests/run_ranking_benchmark.py
npm audit --audit-level=high
pip check
```

The runners discover Python and Node tests recursively below `tests/` while pruning cache and
dependency directories. The directory file-count gate rejects source directories with 100 or more
direct child files, excluding configured dependency, cache, and build-output directories.
Deterministic release tests are offline. Live NCBI and remote BLAST checks are optional integration tests.

## Render

Use `render.yaml` or the included Dockerfile. Configure `OLIGOFORGE_HOSTED=1`, provide an NCBI contact email through environment variables, and store secrets only in Render configuration. Hosted server-side project storage, reaction-setting mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment. Automatic-design jobs are held only in the application process, expire after a configured TTL, and are lost on restart; see `DEPLOYMENT.md` and `API.md`.

## Scientific boundary

OligoForge ranks the **best-supported computational candidate among the retained and fully evaluated pool under declared assumptions**. Search is heuristic-bounded, and software cannot establish a universal wet-lab optimum. Validation Studio and Assurance remain sequence/model evidence: they do not establish efficiency, LOD/LOQ, analytical or clinical sensitivity/specificity, matrix behavior, or regulatory acceptability. Read `RANKING_AUDIT.md`, `RANKING_VALIDATION_REPORT.md`, `ASSURANCE_VALIDATION_LIMITS.md`, `VALIDATION_LIMITS.md`, and `RELEASE_SUMMARY.md` before publication or large-scale ordering.
