# OligoForge Assurance validation report

## Implemented vertical slice

This release implements one complete offline sequence-evidence lifecycle:

1. Normalize and validate a versioned AssaySBOM while preserving order notation and locked components.
2. Build bounded immutable target and off-target sequence snapshots from FASTA/FASTA.GZ plus optional CSV/TSV metadata.
3. Preserve accepted/rejected dispositions, exact duplicates, unique haplotypes and metadata-group counts separately.
4. Calculate a deterministic baseline-to-follow-up delta.
5. Reconstruct complete assay products for every unique supplied sequence through the existing isolate engine.
6. Issue reason-coded DriftGuard states without a cosmetic probability score.
7. Create deterministic OligoForge-local Molecular Vulnerability Records.
8. Package the AssaySBOM, snapshots, deltas, scans, records and Validation Studio plans with per-artifact and package SHA-256 verification.

The frozen regression uses non-pathogen synthetic sequences to exercise an unchanged record, a newly supplied target sequence with a terminal-primer concern, a newly supplied signal-capable off-target, exact-sequence deduplication, group metadata, tamper detection, HTML escaping and the offline CLI.

## Validated properties

- Identical normalized inputs create identical content-addressed identifiers.
- Snapshot hashes fail verification after mutation.
- Incremental delta results equal exact full-set differences for added, removed and unchanged unique sequences.
- Exact duplicate records remain visible while unique-haplotype metrics are not double counted.
- DriftGuard examines complete forward/reverse/probe products, not isolated oligo similarity.
- New modeled target dropout and new modeled signal-generating off-target observations create distinct reason codes and deduplicated OFVRs.
- Evidence-package tampering invalidates package verification.
- Offline tests do not require NCBI or BLAST.

No wet-lab biological-accuracy claim follows from these software regressions.
