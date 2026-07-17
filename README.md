# OligoForge 1.37.0

OligoForge is a computational primer, probe, and qPCR assay-development workbench. Version 1.37.0 is the **design consistency and reliability release**: supported design workflows now declare the same versioned scientific contract, probe-less chemistry resolves consistently, candidate evidence depth no longer depends on display count, and failures carry actionable recovery guidance plus a request identifier. The lifecycle workspaces from 1.36 remain available for reproducible validation plans and sequence-evidence assurance.

## Core behavior

- Target-spanning, bounded joint search over complete forward-primer/reverse-primer/probe triplets.
- Diversity-preserving retention across regions, primer pairs, probe positions, product lengths, and trade-off profiles.
- True hard scientific constraints before objective-specific evidence and soft preferences.
- Coherent target coverage, paired all-site specificity, condition-envelope thermodynamics, junction evidence, and multiplex interactions can change final rank.
- Candidate attrition ledger, structured rank trace, Pareto fronts, truthful finalist categories, and concrete rank-reversal scenarios.
- Explicit near-equivalent/insufficient-evidence states rather than false numerical certainty.
- Deterministic self-hashing manifests with application/ranker/model versions, reaction conditions, inputs, constraints, candidate limits, database state, warnings, and fallbacks.
- A canonical design contract reports the resolved objective, chemistry-profile hash, evidence scope, search limits, engine versions, conformance checks, and qualification state across supported design surfaces.
- Structured API problems and non-sensitive system diagnostics replace opaque red failures with stable codes, retry guidance, recovery actions, and request IDs suitable for support handoff.
- Automatic, batch, sequence-viewer, manual-design, constrained-redesign, and rescue workflows use the authoritative scientific modules and structured ranker.
- Manual mapping reports every exact and allowed near-match, including non-extension-eligible 3′ mismatches, instead of silently hiding failed placements.
- Manual edits can be compared against the last analyzed baseline using complete recalculated evidence: resolved/new hard failures, target/off-target coverage, robustness, interactions, mapping, and rank preference.
- Complete design runs can be compared for candidate loss, winner changes, rank reversals, context changes, and unexplained non-reproducibility.
- Experimental feedback is validated, deduplicated, conflict-checked, locally summarized, and split by target group without silently changing the ranker. A learned reranker remains disabled without sufficient leakage-controlled evidence.
- Long automatic designs return a capability-style job ID promptly, expose real stage states, support cancellation and idempotent retry, and retain the completed primary design when optional BLAST is unavailable.
- The visible Validation Studio accepts Workbench or pasted candidates and FASTA cases, explains why cases were selected, renders the full interleaved 96/384-well plate, exports a fillable CSV, imports completed observations, and reports support, contradiction, control failure, and unresolved uncertainty in plain language.
- The visible Assurance workspace registers the displayed assay or entire Workbench, links baseline and follow-up target/off-target snapshots with explicit offline provenance, shows exact deltas, runs complete-product DriftGuard, issues reason-coded OFVRs, and packages all evidence—including an active validation plan—for verified download.
- Lifecycle progress, form labels, status announcements, responsive plate views, browser-session boundaries, hashes, scope statements, and download points are explicit in the UI. Primary views remain human-readable; machine records remain downloadable.

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
python tests/test_ranking_benchmark.py
python tools/build_release_manifest.py --verify
python tools/build_source_release.py --output-dir dist
python tools/verify_release_checksums.py dist/oligoforge-1.37.0-source.zip.sha256
npm audit --audit-level=high
pip check
```

The runners discover Python and Node tests recursively below `tests/` while pruning cache and
dependency directories. The directory file-count gate rejects source directories with 100 or more
direct child files, excluding configured dependency, cache, and build-output directories.
`test_ranking_benchmark.py` verifies fresh benchmark output in temporary directories and leaves the
source tree unchanged. Run `tests/run_ranking_benchmark.py` only when intentionally regenerating the
committed evidence; add `--regenerate-figures` only for an explicit frozen-figure refresh.
Deterministic release tests are offline. The source builder normalizes archive order, timestamps,
permissions, and storage so identical source bytes produce an identical ZIP; the adjacent SHA-256
sidecar is verified before publication. Live NCBI and remote BLAST checks are optional integration tests.

## Render

Use `render.yaml` or the included Dockerfile. Configure `OLIGOFORGE_HOSTED=1`, provide an NCBI contact email through environment variables, and store secrets only in Render configuration. Hosted server-side project storage, reaction-setting mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment. Automatic-design jobs are held only in the application process, expire after a configured TTL, and are lost on restart; see `DEPLOYMENT.md` and `API.md`.

## Scientific boundary

OligoForge ranks the **best-supported computational candidate among the retained and fully evaluated pool under declared assumptions**. A `computationally_qualified_with_declared_limits` contract status is conditional on the recorded sequence corpus, evidence, chemistry, constraints, and bounded search; it is not a wet-lab or clinical claim. Search is heuristic-bounded, and software cannot establish a universal wet-lab optimum. Validation Studio and Assurance remain sequence/model evidence: they do not establish efficiency, LOD/LOQ, analytical or clinical sensitivity/specificity, matrix behavior, or regulatory acceptability. The browser lifecycle state is session-local and is not an enterprise assay registry. Read `RANKING_AUDIT.md`, `RANKING_VALIDATION_REPORT.md`, `ASSURANCE_VALIDATION_LIMITS.md`, `VALIDATION_LIMITS.md`, and `RELEASE_SUMMARY.md` before publication or large-scale ordering.
