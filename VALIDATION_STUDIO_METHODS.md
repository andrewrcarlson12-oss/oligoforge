# Validation Studio methods

Validation Studio converts a disagreement among two or more computationally evaluated assays into a bounded comparison experiment. It is an experiment-planning and evidence-interpretation tool, not an assay-performance predictor.

## Selection method

1. Normalize each complete forward/reverse/probe assay together with its chemistry. Exact molecular-and-chemistry duplicates are suppressed; candidate identifiers must remain unique.
2. Evaluate every supplied target and off-target case using OligoForge's existing isolate-level coherent-product reconstruction. A probe assay distinguishes a product from a signal-capable product.
3. Record candidate-specific modeled states, terminal-primer concerns, product sizes, primer identities and probe recognition.
4. Select cases deterministically by new candidate-pair distinctions, number of distinct modeled states, target/off-target group diversity and a stable identifier tie-break.

This bounded greedy selection is not claimed to be globally optimal. Every selected case records why it was retained and which candidate pairs it distinguishes.

## Plate layouts

The planner supports 96- and 384-well plates, 1–12 replicates, candidate interleaving, no-template controls, declared positive controls and extraction controls. Candidate order rotates across case/replicate blocks; blocks are deterministically randomized with a recorded seed. Edge wells are either excluded or explicitly flagged. CSV exports neutralize spreadsheet-formula prefixes.

## Result interpretation

Completed CSV files may report amplified/not amplified, Cq, efficiency, melt abnormalities, unexpected products, probe signal, missing observations and notes. An amplified no-template control or failed positive control invalidates the comparison. Otherwise, predictions are labeled supported, contradicted, mixed or missing. Any candidate preference is limited to the declared experiment and does not modify the deterministic ranker.

## Scientific limits

- Modeled products do not establish amplification efficiency, LOD/LOQ, fluorescence, inhibition, matrix effects, clinical sensitivity/specificity or multiplex competition.
- Natural, synthetic and computational cases are labeled separately.
- Candidate comparison is conditional on the supplied cases, controls, reaction conditions and user-declared acceptance criteria.
- Laboratory confirmation remains required before consequential assay decisions.
