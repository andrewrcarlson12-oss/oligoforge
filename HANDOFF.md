# OligoForge 1.31.1 engineering handoff

## Current state

This tree is the repaired continuation of the user-supplied `oligoforge_v1.31.0.zip`. It is a computational pre-release. No wet-lab performance claim is implied.

## Highest-impact changes

- Full-target candidate collection and global ranking replaced first-success/early-window acceptance.
- Exact selected-site coordinates are propagated through design and display.
- Ambiguous template bases are not silently converted into arbitrary oligo bases.
- Offline PCR performs ambiguity-aware, all-site, exact-3′-anchor screening and reports uncertain matches conservatively.
- Multiplex analysis strictly validates oligos, distinguishes assays by identity rather than display name, and reports annealing-context dimer metrics.
- Raw Cq, standard-curve, reference-gene, MIQE-readiness, report, order, and RDML calculations were corrected to avoid false statistical or validation confidence.
- Orthogonal-panel exact proofs, rigorous bounds, and numerical theta diagnostics are separated.
- Hosted deployment isolates shared state, blocks local database paths, sanitizes errors, limits requests, and adds browser security headers.

## Release gate

Run:

```bash
python run_tests.py
```

Network-dependent tools should be manually smoke-tested with a real NCBI email and, where relevant, a local BLAST database. Render/GitHub deployment was not performed from this archive.

## Required empirical work

Read `VALIDATION_LIMITS.md`. Final assay candidates require paired specificity review, target/off-target panels, reaction optimization, efficiency/linearity experiments, product-identity evidence, and dedicated LOD/LOQ studies.
