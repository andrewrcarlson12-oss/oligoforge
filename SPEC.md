# SPEC — Certified Orthogonal Panel (`orthopanel`)

Selects a maximum set of mutually non-cross-hybridizing oligos from a candidate pool under a
thermodynamic (nearest-neighbor ΔG) confusability model, and returns an **optimality certificate**:
how close the returned panel is to the provably largest possible panel *under that model*. No other
primer/probe tool certifies optimality — they pick panels by greedy threshold filtering and stop.

## What the certificate does and does not mean (read this first)

The certificate is exact with respect to the **model**, not the wet lab. It proves the returned panel
is the maximum independent set of the ΔG-threshold confusability graph. That graph is a pairwise
proxy; real multiplex cross-priming is not pairwise (it depends on concentrations, competing
templates, polymerase, cycling). So "CERTIFIED MAXIMUM" means *no larger panel exists under this
pairwise thermodynamic model*, not *this panel is experimentally optimal*. The UI and payload state
this explicitly. This honesty is the point, not a disclaimer to bury.

## Where things live (matches existing conventions)

- **Backend module:** `oligoforge/orthopanel.py` — pure functions. Two layers, deliberately split so
  the graph/certificate math is unit-testable with no thermodynamics:
  - *math core* (no `thermo` import): `max_independent_set`, `greedy_clique_cover`,
    `clique_cover_bound`, `lovasz_theta` (optional), `integer_upper_bound`, `strong_product_capacity`.
  - *thermodynamic driver* (uses `thermo as T`): `intake`, `self_structure_filter`,
    `build_confusability_graph`, and the orchestrator `certify_panel`.
- **Endpoint:** `POST /api/orthogonal-panel` in `app.py`, pydantic `OrthoPanelReq`, same
  `try/except -> JSONResponse(status_code=200)` shape as `/api/autodesign`.
- **Frontend:** one `<section id="orthopanel">` in `static/index.html`, registered in the `TABS`
  array under the **Check** group (it is an analysis, like Panel matrix / in-silico PCR). Vanilla JS
  handler `doOrthoPanel()` + a dependency-free `<canvas>` force-directed view of the confusability
  graph. Matches the existing `card / hint / row / ph / err / pill` styles and the `api()` helper.
- **Tests:** `tests/test_orthopanel.py` (math core: θ known values, invariants, edge cases, golden)
  and `tests/test_orthopanel_thermo.py` (driver: symmetry, filter, graph build). Standalone scripts
  with the repo's `check(name, cond, detail)` convention; auto-discovered by `run_tests.py`.

## Research (verified against primary sources, 2026-07)

- **seqwalk** (Gowri et al., *Nat. Comput. Sci.* 4:423–428, 2024; doi 10.1038/s43588-024-00646-z):
  minimizes sequence symmetry (SSM) via a de Bruijn graph, scales to >10^6 25-nt barcodes in <20 s.
  The paper states outright that SSM "does not explicitly capture thermodynamic properties of
  sequences, and cannot guarantee low off-target binding energies." **This module occupies exactly
  that gap, at curated-panel scale.** We do not reimplement seqwalk and do not compete at 10^6 scale.
- **Prior art on thermodynamic multiplex design exists** and we are not the first there: SADDLE (Xie
  et al., *Nat. Commun.* 13:1881, 2022) does simulated-annealing multiplex primer design on dimer
  likelihood; Prider (2022) does set-cover primer selection. What is **not** in the literature (no
  prior art found) is a Lovász-θ / certified-maximum-independent-set treatment of thermodynamic oligo
  orthogonality. The novel contribution is the *certificate*, not "graphs for panel design."
- **Lovász θ facts, confirmed:** sandwich α(G) ≤ θ(G) ≤ clique-cover(G); primal SDP
  θ = max ⟨J,X⟩ s.t. X⪰0, tr(X)=1, X_ij=0 ∀{i,j}∈E; empty→n, complete→1, C5→√5≈2.236, Petersen→4;
  strong-product multiplicativity θ(G⊠H)=θ(G)·θ(H), so θ(G)^k for k identical rounds. Independent
  implementations (Sage `lovasz_theta`) return Petersen 4.0 and C5 2.236068 — the C5=√5 unit test
  genuinely catches a wrong SDP formulation.

## Judgment calls (and why)

1. **Certificate order: cheap first, θ only if it must.** The sandwich theorem gives a *free* upper
   bound — α(G) ≤ (number of cliques in a clique cover of G) — computable in pure Python by greedy
   clique cover, **zero new dependencies**. On curated panels this bound frequently already equals the
   exact MIS, giving gap 0 and a certified maximum with no SDP at all. θ (tighter) is attempted only
   when the free bound does not close the gap. This preserves OligoForge's minimal-dependency,
   local-first identity: the SDP is the exception, not the default.
2. **θ solver is an optional, import-guarded dependency.** `cvxpy` + `SCS` (both OSS; SCS is
   Apache-2.0, cvxpy Apache-2.0) are imported lazily inside `lovasz_theta`. If absent, the module
   runs and certifies via the clique-cover bound; the payload flags θ as unavailable. MOSEK is a
   faster optional backend if a license is present; not required. **The app never hard-requires
   cvxpy.**
3. **Exact MIS in pure Python, not an ILP-solver dependency.** The spec suggested PuLP/OR-Tools.
   Instead we use a branch-and-bound maximum-independent-set with greedy coloring bounds — exact for
   the small–moderate n where certification is meaningful (sparse panels to ~150 vertices), and
   again **no new dependency**. Above a configurable limit it falls back to greedy min-degree MIS,
   flagged as a lower bound only. (If lower==upper even from greedy, the panel is still provably
   maximum, since a greedy MIS is a valid independent set.)
4. **Thermodynamic engine: primer3-py** (`thermo.calc_heterodimer/homodimer/hairpin`, SantaLucia NN),
   already the engine the whole app uses — no new engine, results consistent with the rest of
   OligoForge. NUPACK (license-restricted) and ViennaRNA (already present, RNA-oriented) were
   considered; no reason to switch for DNA dimer ΔG.
5. **ΔG conventions.** More-negative ΔG = more stable. A candidate is dropped if its hairpin or
   homodimer ΔG is below (more stable than) `self_dg` (default −9 kcal/mol). An edge {i,j} is added if
   the most stable of heterodimer(i,j) and heterodimer(i, revcomp(j)) is below `cross_dg` (default
   −6 kcal/mol, matching the app's existing dimer PASS/WARN/FAIL cutoffs). Conditions (temperature,
   [Na+], [Mg2+], [dNTP]) are applied via `thermo.set_conditions` for the computation and restored
   afterward (snapshot/restore in a `finally`), the same pattern autodesign uses.
6. **Split-pool bound is labeled a theoretical ceiling, not a guarantee.** θ^k is reported as the
   certified collision-free k-round barcode capacity *under the assumption* that two k-round barcodes
   collide iff they cross-hybridize in ≥1 round (the strong-product confusability relation). That
   modeling assumption is stated next to the number. We also report the naive product-of-independent-
   sets count |MIS|^k to show where pairwise-distance reasoning overcounts relative to the θ ceiling.
7. **Positioning.** "Optimality certificate for curated thermodynamic panels," complementary to
   seqwalk (massive SSM libraries), not a competitor at 10^6 scale — because the SDP is infeasible at
   that scale and the O(n²) graph build is too. The certificate is most useful at the n where it is
   computable; that is the honest scope.

## Assumptions

- Inputs are short oligos (primers/probes, ~15–40 nt); very long inputs exceed primer3's dimer
  length cap and return ΔG 0 (handled: they simply never form edges, and this is surfaced).
- DNA with IUPAC degenerate bases supported (`thermo._resolve` expands to a concrete representative
  for the NN calc, as elsewhere in the app); RNA input is uppercased/handled as sequence, dimer ΔG
  via the DNA NN engine (flagged, since primer3 is a DNA model).
- Curated pools: tens to low-hundreds of candidates. The size guard trips above a configurable limit.

## Definition of done (tracked)

Backend module + endpoint; frontend tab with certificate banner + graph canvas + split-pool number;
all tests green (new + existing 43); a runnable demo (`python -m oligoforge.orthopanel` prints a
certificate on a sample panel); this SPEC; a README section. Version bumped. Goldens byte-identical.
