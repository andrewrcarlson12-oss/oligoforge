# OligoForge 1.31.1

OligoForge is a local-first qPCR primer/probe design, QC, specificity-screening, and assay-readiness application. A FastAPI backend runs the scientific calculations; a single-page browser interface provides design, review, panel, analysis, and export workflows.

OligoForge is a **computational design and documentation tool**. It does not certify an assay, establish clinical performance, or replace empirical validation.

## Core capabilities

- Whole-target primer/probe candidate search with global ranking rather than first-window acceptance.
- Exact primer, probe, amplicon, and search-window coordinates, including repeated-sequence targets.
- Primer3 and published nearest-neighbor thermodynamic calculations at configurable reaction conditions.
- Hydrolysis-probe, modified-probe, SYBR, low-Tm, degenerate, and exon-junction-aware design paths.
- Ambiguity-aware target conservation, off-target discrimination, offline in-silico PCR, remote NCBI BLAST, and optional local BLAST.
- Multi-isolate inclusivity/exclusivity screening and coherent full-assay coverage calculations.
- Multiplex dye, amplicon-melt, cross-dimer, 3′-engagement, and annealing-temperature checks.
- Orthogonal-panel graph analysis with exact branch-and-bound model proofs when available, rigorous clique-cover bounds otherwise, and optional non-certifying Lovász-theta diagnostics.
- Raw fluorescence Cq screening, standard-curve analysis, melt-peak analysis, relative-expression analysis, and geNorm-style reference-gene screening.
- Strict IDT-style order CSV, unambiguous synthetic-fragment FASTA, escaped HTML/CSV assay-readiness reports, and RDML 1.3 assay-definition export.

## Important scientific limits

- Selection Tm and displayed Tm use transparent published models but may differ from vendor calculators because concentration conventions, parameter sets, and modification handling differ.
- LNA calculations require explicit `+N` positions and remain model estimates. MGB and proprietary modified-probe behavior require the selected vendor’s calculation.
- A predicted single amplicon, a single melt peak, or a favorable graph model is not proof of biological specificity.
- The reported lowest fully detected standard and exploratory logistic detection estimate are not validated LOD95 values.
- Reference-gene output is geNorm-style pairwise M/V plus a Cq-SD screen. It is not full BestKeeper or NormFinder.
- MIQE-aligned outputs are readiness records, not certification of MIQE compliance.

See [VALIDATION_LIMITS.md](VALIDATION_LIMITS.md) before using results in a publication, regulated workflow, diagnostic claim, or ordering decision.

## Run from source

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
pip install -r requirements.txt
export OLIGOFORGE_EMAIL="you@institution.edu"
uvicorn app:app --host 127.0.0.1 --port 8111
```

Open `http://127.0.0.1:8111`.

Run the complete source-level test gate:

```bash
python run_tests.py
```

Optional local BLAST requires NCBI BLAST+ and a local nucleotide database. Public hosted deployments intentionally block arbitrary server-side BLAST database paths.

## Render / multi-user deployment

`render.yaml` enables hosted mode. In hosted mode, server-side project/panel storage and process-wide reaction-condition mutation are disabled unless the operator explicitly opts in behind appropriate authentication and isolation.

Relevant environment variables:

```text
OLIGOFORGE_HOSTED=1
OLIGOFORGE_ALLOW_SERVER_STORAGE=0
OLIGOFORGE_ALLOW_SHARED_CONDITIONS=0
OLIGOFORGE_EMAIL=<contact email for NCBI>
OLIGOFORGE_NCBI_KEY=<optional NCBI API key>
OLIGOFORGE_MAX_REQUEST_BYTES=5242880
OLIGOFORGE_NCBI_TIMEOUT=30
OLIGOFORGE_NCBI_RETRIES=2
```

Read [SECURITY.md](SECURITY.md) before exposing the service publicly.

## Main modules

```text
app.py                         FastAPI API and hosted-mode controls
oligoforge/design.py           candidate enumeration, pairing, probing, ranking
oligoforge/autodesign.py       full target/off-target workflow
oligoforge/thermo.py           reaction-aware thermodynamics
oligoforge/specificity.py      BLAST, offline PCR, intron/exon checks
oligoforge/conservation.py     target conservation and discrimination
oligoforge/multiplex.py        channel, dimer, and melt compatibility
oligoforge/orthopanel.py       graph-based orthogonal-panel analysis
oligoforge/cq.py               raw-curve Cq screening
oligoforge/quant.py            copies, dilution series, standard curves
oligoforge/refgenes.py         geNorm-style reference-gene screening
oligoforge/report.py           escaped assay-readiness HTML/CSV export
oligoforge/orders.py           strict order CSV and gBlock FASTA export
oligoforge/rdml.py             RDML 1.3 assay-definition export
```

## Reproducibility and release discipline

Dependencies are pinned in `requirements.txt`. The release version is synchronized across the Python package, API, launcher, and UI. `run_tests.py` runs every standalone Python regression script and every Node UI harness, with per-test timeouts so a stalled network or numerical call cannot hang CI indefinitely.

## License

See [LICENSE](LICENSE).
