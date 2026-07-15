# Validation limits

OligoForge 1.31.1 produces computational candidates and assay-readiness records. The following evidence remains necessary before describing an assay as validated or publication-ready.

## Sequence and specificity

1. Confirm the intended reference build, transcript/isoform, strand, and amplicon coordinates.
2. Re-run final primer pairs through NCBI Primer-BLAST or an equivalent paired-primer search against the relevant database and near-neighbor taxa.
3. Test representative target isolates and the biologically plausible off-target panel.
4. Verify product identity by sequencing, size plus a specific probe, or another orthogonal method. A single melt peak alone is insufficient.
5. For RT-qPCR, establish how genomic DNA is excluded: exon-junction placement, DNase treatment, no-RT controls, or a validated intron-spanning design.

## Reaction optimization

1. Confirm final oligo Tm and modifications with the selected vendor and actual master-mix conditions.
2. Run an annealing-temperature gradient.
3. Evaluate primer/probe concentrations and multiplex interactions experimentally.
4. Include no-template, no-RT where applicable, positive, extraction, and inhibition controls.

## Quantitative validation

1. Use independent dilution levels spanning the intended dynamic range.
2. Treat technical replicates as precision measurements, not independent dilution levels.
3. Report slope, efficiency, linearity, residual behavior, replicate dispersion, and excluded data with reasons.
4. Estimate LOD/LOQ using a dedicated replicate-rich experiment around the transition region. OligoForge’s lowest fully detected standard and logistic estimate are screening descriptors only.
5. Validate performance in the actual sample matrix and extraction workflow.

## Modified probes and degenerate oligos

- LNA estimates use explicit `+N` positions and a published model. Sequence context, placement, dye/quencher chemistry, and vendor parameterization matter.
- MGB and other proprietary modifications are not predicted as exact physical Tm values.
- Degenerate pools must be reviewed for synthesis complexity, component concentration, target coverage, and off-target coverage.

## Expression studies

- Average technical replicates before reference-gene analysis.
- Supply assay-specific amplification efficiencies when available.
- Use biological replicates for inferential statistics.
- geNorm-style M/V results should be cross-checked with biological knowledge and an independent method where publication claims depend on reference stability.

## Regulatory and clinical use

OligoForge is not a medical device, diagnostic system, laboratory-developed-test validation package, or regulatory submission system. Do not use its output alone for patient care, release testing, or a regulated performance claim.
