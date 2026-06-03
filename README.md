# OligoForge

One tool for qPCR primer/probe design and QC, run locally. Real primer3
thermodynamics, multi-vendor design-rule linting, NCBI fetch, BLAST specificity,
and intron/exon-junction checking — the whole pipeline in one browser cockpit
backed by a small Python server.

## Why local
A browser-only page can't BLAST, can't reach NCBI reliably, and can't map exon
structure. Running the same UI on a local FastAPI backend removes all three
limits at once — and the thermodynamics are the actual primer3 numbers, not an
approximation.

## What it does
Works on **any organism** — you supply the sequence, accession, NCBI query, and target/off-target sets; nothing is hard-coded to a species. Built and validated on non-model wildlife (a draft-genome songbird and its blood parasites), where model-organism web tools have no annotations.

- **Oligo QC** — Tm (SantaLucia 1998 + Owczarzy salt, at qPCR conditions), GC, hairpin, self/hetero-dimer, homopolymer runs, 3' clamp
- **Reaction-condition aware** — set Mg2+, monovalent, dNTP, and oligo concentration in the header to match your master mix; every Tm in the tool recomputes for your bench, not a generic default
- **Degenerate-aware Tm** — a primer carrying an IUPAC base (e.g. W = A/T) reports the Tm *range* across its resolutions instead of one misleading number
- **Vendor lint** — IDT PrimeTime / Affinity Plus, Thermo MGB, Bio-Rad, generic SYBR; pass/fail against each chemistry
- **Primer-pair / amplicon** QC (Tm match, F×R dimer, amplicon size)
- **Panel cross-dimer matrix** — every oligo against every other
- **Auto-design** from a template: enumerate -> gate -> pair -> probe -> gBlock, then one click into the Workbench or the Order tab (carrying its gBlock standard). Design and auto-design results draw a to-scale **amplicon map** (forward / probe / reverse in genomic context) so you can eyeball product size and probe placement before ordering
- **NCBI fetch** by accession or gene+organism, with isoform-common region
- **Intron / exon-junction check** — flags gDNA-contamination risk
- **BLAST specificity** — remote (NCBI) or local (blastn)
- **In-silico PCR** — predicted amplicons from both primers (catches off-target products a per-primer BLAST misses)
- **Conservation / discrimination** — per-base oligo conservation across a target set, and mismatch scoring vs an off-target set (e.g. Plasmodium vs Haemoproteus)
- **LNA-aware Tm** — honest effective-Tm range for LNA probes
- **Standard-curve copies & IDT order export** — copy number / dilution series; IDT bulk CSV + gBlock FASTA
- **Reference-gene stability** — geNorm M + pairwise V (how many to use) + BestKeeper Cq SD/CV, from a Cq table; picks your >=2 MIQE reference genes
- **Workbench** — the hub. Add a designed assay once, then from its card run QC, in-silico PCR (opens the tab with results), the intron/gDNA check, conservation (oligos + a seeded target query prefilled), a standard curve, generate its gBlock standard, and export an order — every result saves onto that specific assay. Edit name/gene/organism/dye inline; when the in-silico intron check can't resolve (common for non-model organisms NCBI hasn't annotated), "mark handled" records gDNA exclusion as user-asserted so the report stays honest. Save and reload named **projects** (one per study system / lab member), export/import the whole panel as JSON, and load a one-click example panel to start from
- **MIQE report** — one button emits a self-contained HTML report (per-oligo QC recomputed, specificity, validation, MIQE 2.0 completeness checklist) plus a CSV, ready for a supplement
- **RDML export** — one button writes the panel as an RDML 1.2 file (the standard machine-readable qPCR interchange format), so designed assays load straight into instrument/analysis software (LinRegPCR, RDML-Ninja, Bio-Rad CFX, Roche LightCycler) as predefined targets carrying their primers, probe, detection dye, and amplification efficiency — the machine-readable sibling of the MIQE report
- **Help tab** — a short, accurate in-app guide: the workflow, probe/SYBR design targets, suggested cycling, exactly which run ticks each MIQE checklist box, a worked example, and honest limits
- **Multiplex planner** — flags detection-channel conflicts (assays sharing a dye) and cross-assay primer/probe dimers below threshold
- **Paste-safe input** — strips FASTA headers, whitespace, and numbering from pasted sequences and converts RNA->DNA; flags (never silently drops) invalid characters

## Get it — no Python needed
Grab the file for your OS from the project's **Releases** page and double-click it.
A small window opens ("OligoForge is running") and your browser opens the cockpit.
Close that window to stop. That is the whole install.

- **Windows** (`OligoForge-windows-x64.exe`): SmartScreen may flag an "unrecognized app"
  because the build is not code-signed — click **More info -> Run anyway**.
- **macOS**: take `OligoForge-macos-arm64` for Apple Silicon (M1/M2/M3/M4) or
  `OligoForge-macos-intel` for older Intel Macs. On first run, right-click -> **Open** ->
  **Open** to clear Gatekeeper (unsigned). On recent macOS you may instead need
  **System Settings -> Privacy & Security -> Open Anyway**. If it reports the app is
  "damaged" with no Open option, clear the quarantine flag in Terminal:
  `xattr -cr OligoForge-macos-arm64` (until the binaries are code-signed).
- Saved panels live in an `OligoForge/panels` folder under your user directory
  (`%LOCALAPPDATA%` on Windows, your home folder otherwise).

Set your NCBI email once in the header field in the UI; it is remembered after that.

## Run from source (developers)
Python 3.10+.
```
pip install -r requirements.txt
export OLIGOFORGE_EMAIL="you@university.edu"   # Entrez requires an email
uvicorn app:app --reload --port 8111
```
Open http://localhost:8111. Optional, for fast offline specificity: install NCBI
BLAST+ (`blastn`, `makeblastdb`).

## Build the executables
One-off, locally:
```
pip install -r requirements.txt -r requirements-build.txt
pyinstaller oligoforge.spec          # -> dist/OligoForge(.exe)
python ci_smoke.py                   # sanity-check the frozen binary
```
Releases are automated: push a version tag and GitHub Actions builds and attaches
Windows / macOS / Linux binaries to the release.
```
git tag v0.1.0 && git push origin v0.1.0
```
For a native app window instead of the browser, `pip install pywebview` before
building; the launcher picks it up automatically when present.

## BLAST modes
- **Remote** (default): NCBI over the network, no install, queued (30 s to a few minutes). Optional organism filter narrows the search.
- **Local**: point at a blastn database. Build one from a genome FASTA:
  ```
  makeblastdb -in GCF_xxxxx_genomic.fna -dbtype nucl -out fsj_genome
  ```
  then choose mode = local and put the DB path (`fsj_genome`) in the db field.

## Editing
- qPCR salt conditions -> `oligoforge/thermo.py` (`COND`)
- Vendor rules / add a chemistry -> `oligoforge/profiles.py`
- Design + scoring logic -> `oligoforge/design.py`

## Accuracy
Tm tracks IDT OligoAnalyzer to ~1-2 C absolute, exact for relative comparisons.
Hairpin/dimer dG are primer3 estimates. Confirm final sequences in your vendor's
tool before ordering.

## Validated against
Reproduces the hand-built Florida Scrub-Jay HMBS and SDHA assays exactly
(primers, probe, amplicon, QC) and the intron check confirmed both amplicons are
exon-junction-spanning. See `tests/`.

## Layout
```
oligoforge/
  app.py                FastAPI server
  launcher.py           desktop entrypoint (boots server, opens window/browser)
  oligoforge/
    thermo.py           thermodynamics (primer3)
    design.py           enumerate / pair / probe / gBlock
    profiles.py         vendor chemistry profiles + linter
    ncbi.py             Entrez fetch + isoform-common
    specificity.py      intron check + BLAST + in-silico PCR
    conservation.py     per-base conservation / discrimination
    quant.py            copies, dilutions, standard curve
    orders.py           IDT CSV + gBlock FASTA
    autodesign.py       Target -> assay (Auto chemistry + auto-BLAST)
    refgenes.py         reference-gene stability (geNorm + BestKeeper)
  static/index.html     browser cockpit
  tests/                validation scripts (run from repo root)
  oligoforge.spec       PyInstaller build
  ci_smoke.py           frozen-binary smoke test
  .github/workflows/    release build (Windows / macOS / Linux)
  requirements.txt      runtime deps
  requirements-build.txt  build deps (pyinstaller)
  CLAUDE.md             context for Claude Code
```
