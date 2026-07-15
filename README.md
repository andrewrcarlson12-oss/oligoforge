# OligoForge 1.34.0

OligoForge is a computational primer, probe, and qPCR assay-development workbench. Version 1.34.0 is the **decision-analysis release**: it extends the ranking-truth architecture with edit-specific evidence comparison, complete near-match placement reporting, design-run reproducibility comparison, local experimental-evidence summaries, honest benchmark uncertainty, and an authoritative batch-design path.

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
python tests/run_ranking_benchmark.py
npm audit --audit-level=high
pip check
```

Deterministic release tests are offline. Live NCBI and remote BLAST checks are optional integration tests.

## Render

Use `render.yaml` or the included Dockerfile. Configure `OLIGOFORGE_HOSTED=1`, provide an NCBI contact email through environment variables, and store secrets only in Render configuration. Hosted server-side project storage, reaction-setting mutation, and local BLAST paths remain disabled unless deliberately enabled in a private authenticated deployment.

## Scientific boundary

OligoForge ranks the **best-supported computational candidate among the retained and fully evaluated pool under declared assumptions**. Search is heuristic-bounded, and software cannot establish a universal wet-lab optimum. Read `RANKING_AUDIT.md`, `RANKING_VALIDATION_REPORT.md`, `VALIDATION_LIMITS.md`, and `RELEASE_SUMMARY.md` before publication or large-scale ordering.
