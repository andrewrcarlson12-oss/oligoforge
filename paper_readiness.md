# Paper-readiness assessment — OligoForge 1.31.1

## Verdict

The software is suitable for continued internal development, computational benchmarking, and generation of candidate assays. It is **not yet sufficient by itself for a methods-paper claim of assay validation**.

## What can be reported now

- Software architecture and deterministic algorithms.
- Thermodynamic models and stated concentration conventions.
- Offline regression fixtures and benchmark design.
- Input validation, hosted-security controls, and reproducible exports.
- Computational comparisons against independent primary implementations where those comparisons are actually run and reported.

## What still needs empirical evidence

- A preregistered or frozen benchmark panel of targets and near-neighbor off-targets.
- Comparison against Primer3/Primer-BLAST and selected vendor tools using the same input conditions.
- Wet-lab success rate for top-ranked versus lower-ranked or external-tool assays.
- Efficiency, linearity, product identity, inclusivity, exclusivity, matrix effects, reproducibility, LOD, and LOQ.
- Failure-case analysis rather than only successful examples.

See `VALIDATION_LIMITS.md` for the minimum evidence expected before assay-level publication claims.
