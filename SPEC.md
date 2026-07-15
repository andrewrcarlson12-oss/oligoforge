# Orthogonal-panel model specification

## Purpose

`oligoforge/orthopanel.py` selects a mutually non-confusable subset of candidate oligos under a
thresholded, pairwise thermodynamic graph model. An edge means that the modeled cross-dimer is more
stable than the configured threshold. A valid panel is an independent set of that graph.

This module does **not** establish wet-lab multiplex compatibility. Concentration, polymerase,
templates, cycling, higher-order competition, modifications, and matrix effects are outside the
pairwise graph model.

## Formal outputs

The implementation reports three quantities separately:

1. **Constructive lower bound:** the size of the returned independent set. This is always a real
   feasible panel under the graph model.
2. **Rigorous upper bound:** either the completed exact branch-and-bound result or the size of a
   valid clique cover. A clique cover is an upper bound because an independent set can contain at
   most one vertex from each clique.
3. **Numerical Lovász-theta diagnostic:** optional floating-point SDP output. The mathematical theta
   value is an upper bound, but the ordinary numerical solver result is not rounded into a formal
   certificate in this implementation.

`certified=true` means the constructive lower bound equals a rigorous upper bound. It refers only to
the graph model.

## Pipeline

1. Parse and strictly normalize candidate oligos.
2. Merge exact duplicate sequences while retaining names and counts.
3. Exclude candidates whose modeled hairpin or homodimer crosses the self-structure threshold.
4. Build all pairwise cross-dimer edges under a request-local thermodynamic-condition snapshot.
5. Solve maximum independent set by exact branch-and-bound when the configured size and expansion
   limits allow it; otherwise produce a greedy feasible panel.
6. Compute a greedy clique cover as a rigorous upper bound when exact search does not complete.
7. Optionally compute numerical Lovász theta for diagnostic context only.

## Split-pool output

`|panel|^k` is a constructive combinatorial count assuming independent reuse of the selected panel
across `k` rounds. `theta^k`, when present, is a numerical graph diagnostic. Neither quantity is an
experimentally demonstrated barcode capacity.

## Computational scope

- Graph construction is quadratic in the number of surviving candidates.
- Exact branch-and-bound is exponential in the worst case and may fall back to a greedy lower bound.
- Inputs are intended to be curated primer/probe-sized DNA oligos, not genome-scale libraries.
- Degenerate and modified oligos are reduced to the thermodynamic models explicitly supported by
  the application; model limitations remain visible in the output.

## Verification

Run `python tests/test_orthopanel.py` and `python tests/test_orthopanel_thermo.py`. The tests cover
known graph families, exact/heuristic paths, rigorous bounds, request-local reaction conditions,
input normalization, and edge construction.
