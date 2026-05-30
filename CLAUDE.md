# OligoForge — context for Claude Code

This file is for a Claude Code session picking up this project. It carries the
design context so you don't restart from zero.

## What this is
A local qPCR primer/probe design + QC cockpit. FastAPI backend (`app.py`) +
vanilla HTML/JS front-end (`static/index.html`). The engine is in the
`oligoforge/` package. Everything was built and validated against real Florida
Scrub-Jay (*Aphelocoma coerulescens*) RefSeq data; see `tests/`.

Run: `pip install -r requirements.txt` then `uvicorn app:app --reload --port 8111`,
open http://127.0.0.1:8111. Set `OLIGOFORGE_EMAIL` for NCBI Entrez.

## Architecture
- `oligoforge/thermo.py` — primer3 (SantaLucia + Owczarzy) at qPCR salt (`COND`). Tm, hairpin, dimers, runs. Plus `tm_lna` (LNA-aware Tm range).
- `oligoforge/design.py` — `design_assay(template, profile)`: enumerate -> gate -> pair -> probe-search -> gBlock.
- `oligoforge/profiles.py` — five chemistry profiles + `lint_oligo` / `lint_pair`. Add vendors here.
- `oligoforge/ncbi.py` — Entrez: `fetch_accessions`, `fetch_isoforms`, `common_region`, `gene_id`.
- `oligoforge/specificity.py` — `intron_check` + `blast_remote`/`blast_local`, plus `epcr` + `in_silico_pcr` (in-silico PCR).
- `oligoforge/conservation.py` — per-base oligo conservation across a target set + discrimination vs an off-target set (oligo-anchored, no MSA).
- `app.py` — endpoints incl. /api/{qc,pair,matrix,design,fetch,intron,blast,conservation,epcr,lna_tm,fetch_nuc,order_csv,copies,batch_design,panel/*}.

## Scientific context (the owner's project)
Senior honors thesis: helminth-mediated Th1/Th2 immune polarization in the
Florida Scrub-Jay. A 5-target qPCR panel is being ordered from IDT, budget-tight.
Targets: 4 host genes by RT-qPCR expression + 1 parasite detection target.
Reference-gene strategy: screen 4 candidates by SYBR (cheap primers), pick the
best 2 by geNorm/NormFinder/BestKeeper, then buy ZEN probes + gBlocks only for
the two winners. Avoid GAPDH/ACTB/B2M (immune-regulated).

### The locked panel (verified by hand and reproduced by this tool)
Host genes — IDT PrimeTime ZEN (5'FAM/ZEN/3'IBFQ), primers 25 nmol std desalt, 60 C 2-step:
- IFNG : F AGTCATTCTGATGTCGCTGATG | R ACCTGTCAGTGTTTTCAAGCA | P TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT | amp 136
- IL4  : F AACTTGCTCAGCCTGGTTTG | R ATTCTTTAGTGAGGTGGTGCTG | P CTTGTGCCCTGCTCTGGTCCC | amp 84
- RPL13: F TCGCTGGCATCAACAAGAAG | R TCGGGAAGAGGATGAGCTTG | P AACAAGTCCACCGAGTCCCTGCA | amp 138
- YWHAZ: F CCGTTACTTGGCTGAGGTTG | R GATGGGATGTGTTGGTTGCA | P CCACTATCCCTTTCTTGTCATCTCCAGCAG | amp 121
Parasite — IDT Affinity Plus (5'FAM/3'IBFQ, LNA), 54 C singleplex:
- Plas : F TACCTGGACTWGTTTCATGG (W=A/T) | R AAAGGATTTGTGCTACCTTG | P CTTACAAGATATCCACCACA | amp 157

### Reference-gene candidates (screen these 4 by SYBR, keep 2)
- RPL13, YWHAZ (above)
- HMBS  (XM_068994916): F GAGCTATACCCCGACCTCTG | R CTTCTCTCCAATCTTGGAAAGCG | P ATCTTGTCCCCAGTTGTTGACATGGCC | amp 93 | intron-spanning (junctions 253,303)
- SDHA  (XM_069004197): F ACTCCAAGAAGGCTGTGAGAAA | R ACAGAGCATCAGGTTCTGCA | P CAAGGGTCTCCACCAAGTCAGTGTTCC | amp 130 | intron-spanning (junction 1707)
gBlock standards:
- HMBS (173 bp): GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATTGTGGCCATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCTTCACCAAAGAGTTGGAAAATGCACTTGAAAGAA
- SDHA (210 bp): CAATGCAAAACCATGCTGCTGTATTTCGTACTGGTTCTGTACTCCAAGAAGGCTGTGAGAAACTTAGCCAAATTTATGGTGATCTGGCTCATCTAAAGACTTTTGACAGAGGTATTGTGTGGAACACTGACTTGGTGGAGACCCTTGAACTGCAGAACCTGATGCTCTGTGCTCTACAAACCATTTATGCTGCGGAGGCTCGCAAAGAGT

## Built in stage 5 (endpoints live, all with UI)
- **Standard-curve copies** — `/api/copies`: copies/uL from ng/uL + length, with dilution series (`quant.py`). UI: Quant tab.
- **IDT order export** — `/api/order_csv`: bulk CSV (Name,Sequence,Scale,Purification) + gBlock FASTA; ZEN probes wrap /56-FAM/.../3IABkFQ/ (`orders.py`). UI: Order tab. Order ZEN probes via PrimeTime for the double-quench + BOGO.
- **Save / load panels** — `/api/panel/{save,list,load}`: named oligo sets as JSON on disk (no localStorage). UI: Save/Load on the matrix.
- **Batch design** — `/api/batch_design`: design_assay over many templates at once. UI: Batch card under Design.
- **Pair specificity** — `/api/pair_specificity`: BLAST forward + reverse, remote or local. UI: "check specificity" button in design output.

## Built in stage 6 (accuracy + specificity)
- **Conservation / discrimination** — `/api/conservation` + `conservation.py`. Per-base conservation across a target set and mismatch scoring vs an off-target set. Validated on real Plasmodium vs Haemoproteus cytb (probe 99% conserved across Plasmodium, >=3 mismatches to closest Haemoproteus). UI: Conservation tab, with NCBI fetch.
- **In-silico PCR** — `/api/epcr` + `specificity.in_silico_pcr`. BLAST both primers, report predicted products (same subject, opposite strands, convergent, in size range). UI: in-silico PCR tab. Combiner unit-tested.
- **LNA-aware Tm** — `thermo.tm_lna` + `/api/lna_tm`. DNA-backbone Tm plus an honest effective-Tm range for LNA oligos, not a misleading point value.
- **NCBI nucleotide fetch** — `/api/fetch_nuc` pulls target/off-target sets into the conservation tab.
- **Regression suite** — `tests/test_regression.py` pins the locked panel, HMBS reproduction, ePCR, conservation, LNA. Offline.

## Genuine next steps (still open)
1. **Multiplex checker** — dye/channel assignment + cross-assay dimer gate for pooling. Panel is singleplex; every probe is 5'FAM, so nothing co-runs as-is.
2. **Efficiency / LOD from Cq** — take the standard-curve run's Cq values, fit, report E = 10^(-1/slope)-1, R2, LOD (MIQE-required).
3. **SYBR melt-curve Tm** — predict amplicon Tm, flag secondary products for the reference-gene screen.
4. **More vendors** — Thermo / Sigma order formats beside IDT.
5. **Empirical calibration** — feed observed efficiency/Tm back to compare against predicted over time.

## Conventions
Keep thermodynamics flowing through primer3 so numbers stay consistent with the
hand-verified panel. Profiles are guideline defaults — keep them editable, don't
bake assumptions into the engine. The owner verifies finals in OligoAnalyzer +
the Complexity Tool before ordering, so this tool is a fast first pass and linter,
not the final authority.
