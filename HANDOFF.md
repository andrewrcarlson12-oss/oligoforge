# OligoForge — handoff brief (paste this as your first message in a new chat)

You are continuing work on **OligoForge**, an existing project. Everything you need is below and in the
attached `oligoforge.zip`. Do not start over. Read the code before changing it.

---

## 0. START HERE
1. I have uploaded `oligoforge.zip` with this message.
2. Extract it and verify, running these first:
```
cd /home/claude && unzip -o /mnt/user-data/uploads/oligoforge.zip -d /home/claude
cd /home/claude/oligoforge && pip install --break-system-packages -q primer3-py biopython fastapi "uvicorn[standard]" httpx
OLIGOFORGE_EMAIL=arcarl27@colby.edu python3 tests/test_regression.py        # expect: ALL REGRESSION ASSERTS PASS
node tests/ui_handlers.js && node tests/ui_workbench.js && node tests/ui_conditions.js && node tests/ui_projects.js
```
(The zip's top folder is `oligoforge/`, so it lands at `/home/claude/oligoforge`.)

---

## 1. WHO I AM / HOW TO WORK WITH ME
- Andrew. Senior bio major (cell/molecular + biochem), Colby '27, pre-med. Senior RA in Dr. Anna Forsman's
  Wild Symbioses Lab. Florida Scrub-Jay (Aphelocoma coerulescens) eco-immunology; honors thesis on
  helminth-driven Th1/Th2 polarization.
- Tone: terse, direct, no fluff. Avoid AI-isms (delve, crucial, leverage, comprehensive, robust, nuanced,
  foster, moreover, furthermore, additionally as openers). Active voice. Sparse em dashes.
- Do all the work IN THE CHAT — build, test, fix here. Don't hand me homework.
- Two standing bars: keep the science/data perfect (zero mistakes), and make my life easier.
- My pattern: I ask "what's worth adding / can we make it better?" then say "go" or "do it." When I say go,
  build it AND verify it before telling me it's done.
- I'm on Windows 11, Python 3.14.5, Firefox. I use `localhost`, not `127.0.0.1` (Firefox balks at the IP).
- GitHub: github.com/andrewrcarlson12-oss. Repo: **andrewrcarlson12-oss/oligoforge** (public). NCBI email: arcarl27@colby.edu.
- This OligoForge work is technical/conversational. Do NOT apply the human-writing-style skill to it — that
  skill is only for writing prose for me (manuscripts, essays, emails, statements).

## 2. WHAT OLIGOFORGE IS
Local-first qPCR primer/probe design + QC web app for non-model organisms. FastAPI backend + a single-file
browser cockpit (`static/index.html`) + engine package `oligoforge/`. It designs/validates primer+probe
assays and runs the full chain: specificity (in-silico PCR / BLAST), exon/intron placement, conservation
across isoforms, reference-gene stability, multiplex compatibility, standard curves, and a MIQE-style
HTML+CSV report. Built on real NCBI data. Differentiator vs Primer-BLAST / QuantPrime / Beacon Designer:
the assembled offline pipeline for non-model organisms + parasite discrimination + no data hand-off.

**Current version: v1.2.2.** Three places hold the version and MUST stay in sync on every bump:
footer in `static/index.html`, `TOOL_VERSION` in `oligoforge/report.py`, `APP_VERSION` in `launcher.py`.

## 3. HOW I RUN IT (Windows, confirmed working)
Files extract to a NESTED folder: `C:\Users\andre\OneDrive\Desktop\oligoforge\oligoforge` (zip top dir is
`oligoforge/`, so Windows "extract" makes oligoforge\oligoforge; app.py is in the nested folder).
- Source: `cd` into the nested folder -> `python -m pip install -r requirements.txt` (first time) ->
  `python -m uvicorn app:app --port 8111` -> open http://localhost:8111.
- Double-click: `OligoForge.bat` (runs launcher.py, opens localhost; close the black window to stop).
- Standalone .exe (built + confirmed): `python -m pip install -r requirements-build.txt` then
  `python -m PyInstaller oligoforge.spec --noconfirm` -> `dist\OligoForge.exe`. The
  "Permission denied / retrying set_exe_build_timestamp" warning is harmless (OneDrive touching the file).
- Update from a new zip: stop, delete the old Desktop\oligoforge folder, extract the new zip there, restart.
  Ctrl+Shift+R in Firefox to drop cached CSS.

## 4. REPO LAYOUT (/home/claude/oligoforge)
- `app.py` — FastAPI, 30 `/api` routes, frozen-aware paths (sys._MEIPASS; %LOCALAPPDATA%/OligoForge when frozen),
  PANELS_DIR + PROJECTS_DIR. Endpoints return `{error}` JSON 200 on failure.
- `oligoforge/` — engine: thermo.py, design.py, profiles.py (7 chemistry profiles), ncbi.py, quant.py,
  orders.py, conservation.py, specificity.py, autodesign.py, refgenes.py, report.py, multiplex.py.
- `static/index.html` — browser cockpit: 3 style blocks (`<style>`, `#polish`, `#polish2`), data-driven tabs,
  global `api(path,body)` (always POSTs JSON; now also injects email+key centrally). Workbench is a CARD layout.
- `launcher.py` — finds a free port (8111+), serves uvicorn (binds 127.0.0.1), opens http://localhost:PORT,
  optional pywebview window, startup GitHub update-check. APP_VERSION + GITHUB_REPO live here.
- `oligoforge.spec`, `Dockerfile`, `render.yaml`, `.github/workflows/build.yml`, `ci_smoke.py`, `.gitignore`,
  `requirements.txt`, `requirements-build.txt`, `OligoForge.bat`, `README.md`, `tests/`.

## 5. ENGINE + NCBI SPECIFICS
- `thermo.py`: primer3 wrappers at qPCR salt COND=dict(mv_conc=50,dv_conc=3,dntp_conc=0.8,dna_conc=200)
  (SantaLucia 1998 + Owczarzy). IUPAC-resolving, so degenerate bases (e.g. W) never crash. `set_conditions()`
  mutates COND; `tm_range()` enumerates the degenerate Tm spread.
- `specificity.py`: `intron_check` auto-locates amplicon coords from the forward/reverse primers via IUPAC
  regex when coords are absent. in_silico_pcr / remote BLAST present.
- `report.py` (TOOL_VERSION="OligoForge v1.1.2"): `build(panel, meta)` -> {html, csv, n_assays, tool, generated}.
  Self-contained MIQE-style HTML — inline CSS is built by STRING CONCATENATION on purpose, because the CSS
  contains a literal "100%" and %-formatting would break. Do not switch it to %-format.
- `multiplex.py`: `check(assays, dimer_threshold=-9.0)` -> channel_conflicts (shared dye) + cross_dimers + n_flagged.
- `refgenes.py`: geNorm + BestKeeper. `ncbi.py`: socket.setdefaulttimeout(45) so fetches fail fast.
- **NCBI auth:** `app._set_email(email=None, key=None)` sets Entrez.email from arg or env OLIGOFORGE_EMAIL,
  and Entrez.api_key from arg or env OLIGOFORGE_NCBI_KEY (request key beats env; None -> 3 req/s, key -> 10 req/s).
  Called at startup `_set_email(None)` and per-request `_set_email(r.email, r.ncbi_key)`. All 7 NCBI request
  models carry `email` + `ncbi_key`. The UI `api()` injects email (if absent) + key into every body; non-NCBI
  endpoints ignore the extra field. Header has masked `#ncbi_key` (localStorage of_ncbi_key) next to `#email`.

## 6. THE LOCKED qPCR PANEL — DO NOT ALTER THESE SEQUENCES
Seeded verbatim by the "Load example (FSJ)" button (`seedFSJ`). For the IDT order (Wilson/Burtt $3,000 grant).
Host (IDT PrimeTime ZEN, 5'FAM/ZEN/3'IBFQ, 60 C):
- IFNG:  F AGTCATTCTGATGTCGCTGATG   R ACCTGTCAGTGTTTTCAAGCA      P TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT  amp 136
- IL4:   F AACTTGCTCAGCCTGGTTTG     R ATTCTTTAGTGAGGTGGTGCTG     P CTTGTGCCCTGCTCTGGTCCC              amp 84
- RPL13: F TCGCTGGCATCAACAAGAAG     R TCGGGAAGAGGATGAGCTTG       P AACAAGTCCACCGAGTCCCTGCA            amp 138
- YWHAZ: F CCGTTACTTGGCTGAGGTTG     R GATGGGATGTGTTGGTTGCA       P CCACTATCCCTTTCTTGTCATCTCCAGCAG     amp 121
Plasmodium (IDT Affinity Plus LNA, low-Tm TaqMan, 54 C, FAM, genus-specific vs Haemoproteus):
- F TACCTGGACTWGTTTCATGG (W=A/T at pos 11)   R AAAGGATTTGTGCTACCTTG   P CTTACAAGATATCCACCACA   amp 157
Ref-gene candidates HMBS + SDHA (intron-spanning, gBlocks) noted, not yet ordered. Bench gates still open
(annotation confirm, SYBR ref screen, sequence-confirm positives, OligoAnalyzer final check before ordering).
Locked metrics: S/B 15-20, intra-assay CV 4.4%, inter-assay CV 9.9%, avidity 45.7-81.8%, ACU median ~7.20,
seroprevalence 75-83% (threshold-dependent). Plasmodium W primer Tm 56.4-57.5 C across its 2 resolutions.

## 7. DEPLOYMENT STATE
- **GitHub: DONE.** Code is pushed to andrewrcarlson12-oss/oligoforge (public). Future updates via
  GitHub Desktop: edit -> Commit -> Push.
- **Hosted URL (the "any computer, always latest" path): NOT yet deployed** (needs my Render account/click).
  `Dockerfile` (python:3.12-slim, `uvicorn app:app --host 0.0.0.0 --port ${PORT:-8111}`) + `render.yaml`
  (Render web/docker/free, autoDeploy:true, healthCheckPath:/, OLIGOFORGE_EMAIL=arcarl27@colby.edu,
  OLIGOFORGE_NCBI_KEY as a `sync: false` dashboard secret). Deploy path: render.com -> New+ -> Blueprint ->
  pick the repo -> Apply -> get URL; every push auto-redeploys. Free tier sleeps ~15 min idle (cold start
  ~30-60 s). Server-side projects are ephemeral on free tier; the browser Workbench (localStorage) persists
  per-user. Dockerfile/render.yaml are standard config, NOT build-tested in the sandbox.
- **Installers: NOT yet built on CI.** `.github/workflows/build.yml` builds Windows / macOS-arm64 /
  macOS-Intel / Linux on a tag push (py3.12), runs ci_smoke.py, uploads uniquely-named assets to a Release.
  Trigger by creating + pushing a tag (e.g. v1.1.2). The Windows .exe has been built locally already.
- **launcher.py update-check:** on startup queries GitHub releases/latest and prints
  "Update available: vX — download <url>" if newer than APP_VERSION; silent on failure/no releases.

## 8. WHAT'S DONE
Workbench as wrapping assay CARDS (was a wide table that ran off-page) + global overflow containment so wide
results scroll inside their panel; dark-industrial visual polish + tidier tabs; user-settable reaction
conditions + degenerate Tm range; Load-example(FSJ) seed; copy buttons; 45 s NCBI timeout; server-side named
Projects (save/list/load/delete); footer version scheme; localhost launcher + OligoForge.bat + verified
Windows .exe; Dockerfile + render.yaml (hosted auto-deploy, email baked, key as secret); CI builds incl.
Intel-Mac with unique asset names; startup update-check; intron amplicon auto-location from primers; MIQE
HTML+CSV report + multiplex checker; NCBI api_key support (env + masked UI field + per-request, key beats env);
.gitignore; UI test harnesses moved into tests/.

**v1.1.3 (this revision):** fixed a savePanel name collision — the Panel-matrix "Save" button was silently calling the workbench's localStorage saver (duplicate `function savePanel`), so named oligo sets never reached /api/panel/save; the matrix saver is now `saveOligoPanel` and Save works again. Added workbench Export/Import as a JSON file (client-side backup + portability, since hosted projects are ephemeral on Render's free tier): `exportPanel()` downloads `oligoforge_panel.json`, `importPanel()` reads one back. useTmpl now calls `gotab("design")` instead of a hardcoded nav index. Assay-id generation consolidated into one `uid()` helper (was three slightly different inline generators). The MIQE report's two downloads are staggered so the browser doesn't suppress the CSV. Removed a duplicate `os.makedirs(PANELS_DIR)` and a duplicate `IntronReq.mrna_acc` field. Aligned the two stray version markers (`oligoforge/__init__.py` `__version__`, FastAPI `version=` in app.py) to the same number as the three tracked spots. Engine left untouched (no offline test covers it the way the API layer does).

**v1.2.0 (this revision — UI/UX redesign, engine still untouched):** full cosmetic reskin to a "precision lab-instrument" look — one cohesive stylesheet replacing the three stacked `<style>` blocks, IBM Plex Sans/Mono web fonts (with offline fallbacks), deeper slate palette, real depth/shadows, refined cards/tables/inputs/buttons. Navigation regrouped + reordered into Create / Check / Quantify / Finish with group labels, led by the create-from-target flow (the `auto` tab is now TABS[0] and the default `on` section; `#qc` no longer carries `class="on"`). The `auto` tab gained a guided hero (title + plain-language description + 4 step chips). New `toast()` + `flash()` helpers: most informational `alert()`s are now toasts (a `<div id="toast" class="toast-wrap">` lives just before `</body>`); destructive `confirm()` gates are unchanged. `fld()` now flashes every field it autofills, so Pair/in-silico-PCR/Intron prefills visibly highlight. Slow remote ops (doAuto, doEpcr, doPairSpec) show an animated `.ofbar` progress bar. gBlock made one-click: the Design result has copy + "→ Order as standard" buttons (`gblockToOrder()`), `designAdd`/`autoAdd` now carry `gblock` onto the assay, and the Workbench Order button prefills the assay's gBlock standard when present. TABS entries are now `[id,label,group]`; the nav-build inserts a `.navgroup` span per group and tags the first tab `lead`. All harness contracts preserved (one `<script>`, one `</head>`, every parsed function/id intact, top-level localStorage still try/caught).

**v1.2.1 (this revision):** raised the global NCBI socket timeout from 45 s to 120 s and made it env-configurable via `OLIGOFORGE_NCBI_TIMEOUT` (`oligoforge/ncbi.py`) — the 45 s cap was too aggressive for a real efetch of several/large records. The **Workbench is now the hub**: each assay can run **QC primers** (Tm / gap / dimers via /api/pair, `wbPair`), **Specificity (in-silico PCR)** (existing `checkAssay`), and **Intron / junction** (/api/intron, `wbIntron`) in place, and each result is SAVED onto that specific assay (`a.results.pair`, `a.checks.specificity`, `a.results.intron`) so it persists in localStorage and renders in a per-assay results block (`assayResultsHtml`); each result links back to its full tab (`toPair`/`toEpcr`/`toIntron` are unchanged and still harness-tested). Every assay shows a suggested cycling protocol (`suggestProtocol`) — TaqMan two-step 95/60 for probe assays, SYBR three-step with Ta ≈ lowest primer Tm − 4 °C + melt for no-probe — which closes the one MIQE-2.0 design gap (annealing-temperature guidance). New CSS `.acard-res/.resrow/.reslab/.resval`. Engine untouched apart from the `ncbi.py` timeout constant; all harness contracts still hold.

**v1.2.2 (this revision — correctness + robustness pass, no UI redesign):** five fixes found by a full module audit. (1) `orders.py` `_wrap` now emits an LNA / Affinity Plus probe (`probe_lna`) with `/56-FAM/{seq}/3IABkFQ/`, the same 5'FAM/3'IBFQ ends as a ZEN probe — previously an LNA probe row in the order CSV was a bare, unlabeled (non-functional) oligo. `/api/order_csv` also now accepts IDT `+N` LNA notation for `probe_lna` (it validates the DNA backbone and keeps the +N positions in the ordered sequence) instead of rejecting the `+` via `clean_seq`. (2) `/api/copies` now rejects a dilution `factor <= 1` with a clean error and clamps `points` to 1–40; a direct call with `factor=0` previously hit a ZeroDivisionError → 500, and a huge `points` could blow memory (the UI never triggered either). (3) the Workbench intron check (`wbIntron`) now writes a definitive verdict to `a.checks.intron` — the field `report.py` reads — so the MIQE report's "Exon/intron location (gDNA)" checklist line finally reflects a run (it was previously unsatisfiable because nothing wrote `a.checks.intron`); error / could-not-locate states stay in `a.results.intron`, and `assayResultsHtml` reads `a.checks.intron` first then falls back. (4) untrusted NCBI/BLAST free text (fetch descriptions, BLAST subject titles in `doFetch`/`doBlast`/`doPairSpec`) is now `esc()`-escaped before going into `innerHTML`. (5) removed a pointless `from .autodesign import _reference` self-import in `autodesign.py`. Audited and found clean (no change): `thermo`, `design`, `quant` (standard-curve + copies math is fully div-by-zero-guarded), `conservation`, `profiles`, `multiplex`, `specificity` (epcr convergence + intron junction-span geometry verified correct), `refgenes`, `report`. Note for the future: the public Render instance still runs every visitor's BLAST/fetch/autodesign under the baked-in NCBI email + key and writes panels/projects with no per-user quota — fine for personal use, but would need auth/rate-limiting before wider exposure.

## 9. NOT DONE / NEXT STEPS
- Deploy on Render (my account/click) -> get the live URL. After it's live, confirm NCBI fetch + BLAST fire server-side.
- Tag + push a release (v1.1.2) to auto-build the Win/Mac/Linux installers.
- Optional, on request: persistent disk for hosted projects (paid); true RDML XML export (needs schema
  validation that can't be done in-sandbox; currently MIQE HTML+CSV only); full silent self-updating binary
  (fragile/platform-specific — deliberately not shipped, would risk bricking on update); empirical Tm
  calibration; SNP/VCF masking under primer 3' ends; pywebview windowless native app; free NCBI key already wired.
- Bench gates before the IDT order (listed in section 6).

## 10. HARD RULES / GOTCHAS (do not break these)
- Keep TOOL_VERSION (report.py) + APP_VERSION (launcher.py) + footer (static/index.html) IN SYNC every bump.
- Preserve every class/id name the JS depends on when restyling — the harnesses parse them.
- NEVER use a truncate-before-read one-liner like `open(p,'w').write(open(p).read()...)`; it empties the file
  (this emptied report.py once — regression caught it). Read into a variable first, then write.
- After ANY edit: re-run `tests/test_regression.py` + all four `tests/ui_*.js` harnesses + a clean-unzip boot,
  and REBUILD the zip (it is built from /home/claude/oligoforge) BEFORE presenting.
- Deliverables go to: `/mnt/user-data/outputs/oligoforge.zip` and `/mnt/user-data/outputs/OligoForge-README.md`.
  Build the zip with:
  `cd /home/claude && zip -rqX /mnt/user-data/outputs/oligoforge.zip oligoforge -x "*/build/*" "*/dist/*" "*/__pycache__/*" "*.pyc" "*/panels/*" "*/projects/*"`

## 11. TESTING (run all from the repo root)
- `OLIGOFORGE_EMAIL=x@y.com python3 tests/test_regression.py` -> "ALL REGRESSION ASSERTS PASS" (offline;
  pins locked Tms, clean_seq, refgenes, report build, multiplex, tm_range, set_conditions, project roundtrip,
  NCBI api_key wiring).
- `node tests/ui_handlers.js` (32 tab-handler renders), `node tests/ui_workbench.js` (12, card flow),
  `node tests/ui_conditions.js` (5, conditions + FSJ seed), `node tests/ui_projects.js` (4, projects).
  They read static/index.html relative to CWD, so run from the repo root.
