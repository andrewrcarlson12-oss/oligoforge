"""OligoForge backend — FastAPI server exposing the engine to the browser cockpit.

Run:  uvicorn app:app --reload --port 8111
Then open http://127.0.0.1:8111
Set your NCBI email once (env var or the field in the UI):  export OLIGOFORGE_EMAIL=you@uni.edu
"""
import os, sys
from typing import Dict, List, Optional
from fastapi import FastAPI
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

from oligoforge import thermo as T, design as D, profiles as P, ncbi, specificity as SP, isolates as ISO

app = FastAPI(title="OligoForge", version="1.21.2")
HERE = os.path.dirname(os.path.abspath(__file__))
# When frozen by PyInstaller: read-only resources (static/) live under sys._MEIPASS,
# and user data (saved panels) must go somewhere writable, not the temp unpack dir.
RES_DIR = getattr(sys, "_MEIPASS", HERE)
if os.environ.get("OLIGOFORGE_DATA_PATH"):
    DATA_DIR = os.environ["OLIGOFORGE_DATA_PATH"]   # mount a persistent volume here in Docker/Render
elif getattr(sys, "frozen", False):
    _base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    DATA_DIR = os.path.join(_base, "OligoForge")
else:
    DATA_DIR = HERE


# ---------- build identity (so a deployed instance can prove which commit it is) ----------
import subprocess, datetime


def _build_commit():
    """Short commit of the running build. Render/most hosts inject this as an env var;
    fall back to git if a checkout is present, else 'unknown'. Surfaced in /healthz and
    the page footer so a stale deploy or a browser-cached page is immediately obvious."""
    for _k in ("RENDER_GIT_COMMIT", "GIT_COMMIT", "SOURCE_COMMIT", "COMMIT_SHA"):
        _c = os.environ.get(_k)
        if _c:
            return _c[:7]
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       cwd=HERE, stderr=subprocess.DEVNULL, text=True).strip() or "unknown"
    except Exception:
        return "unknown"


BUILD_COMMIT = _build_commit()
BOOT_TIME = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _set_email(email=None, key=None):
    e = email or os.environ.get("OLIGOFORGE_EMAIL")
    if e:
        ncbi.Entrez.email = SP.Entrez.email = e
    k = key or os.environ.get("OLIGOFORGE_NCBI_KEY")
    ncbi.Entrez.api_key = SP.Entrez.api_key = (k or None)


_set_email(None)


# ---------- logging (visibility on a hosted instance; level via OLIGOFORGE_LOG_LEVEL) ----------
import logging, time as _time
logging.basicConfig(level=os.environ.get("OLIGOFORGE_LOG_LEVEL", "INFO").upper(),
                    format="%(asctime)s %(levelname)s oligoforge %(message)s")
log = logging.getLogger("oligoforge")


@app.middleware("http")
async def _log_requests(request, call_next):
    t0 = _time.perf_counter()
    try:
        resp = await call_next(request)
    except Exception:
        log.exception("unhandled error %s %s", request.method, request.url.path)
        raise
    log.info("%s %s -> %s %.0fms", request.method, request.url.path,
             resp.status_code, (_time.perf_counter() - t0) * 1000.0)
    return resp


# ---------- request models ----------
class OligoReq(BaseModel):
    seq: str; role: str = "primer"; profile: str = "idt_taqman"

class PairReq(BaseModel):
    forward: str; reverse: str; amplicon: Optional[int] = None; profile: str = "idt_taqman"

class MatrixReq(BaseModel):
    oligos: Dict[str, str]

class DesignReq(BaseModel):
    template: str; profile: str = "idt_taqman"

class FetchReq(BaseModel):
    accession: Optional[str] = None; gene: Optional[str] = None
    organism: Optional[str] = None; isoform_common: bool = False; email: Optional[str] = None; ncbi_key: Optional[str] = None

class IntronReq(BaseModel):
    gene: str; organism: str
    amp_start: Optional[int] = None
    amp_end: Optional[int] = None
    forward: Optional[str] = None
    reverse: Optional[str] = None
    mrna_acc: Optional[str] = None; email: Optional[str] = None; ncbi_key: Optional[str] = None

class BlastReq(BaseModel):
    seq: str; organism: Optional[str] = None; mode: str = "remote"
    db: str = "nt"; db_path: Optional[str] = None; top: int = 40
    email: Optional[str] = None; ncbi_key: Optional[str] = None


# ---------- routes ----------
@app.get("/")
def index():
    # no-store: the cockpit HTML must never be served stale, or a redeploy looks like it
    # "didn't take" (footer/version frozen). The page reads its real version from /healthz.
    return FileResponse(os.path.join(RES_DIR, "static", "index.html"),
                        headers={"Cache-Control": "no-store, no-cache, must-revalidate",
                                 "Pragma": "no-cache", "Expires": "0"})

@app.get("/healthz")
def healthz():
    """Cheap readiness probe (for Render/uptime checks): version, engine, NCBI auth/cache, data dir."""
    try:
        primer3_ok = isinstance(T.tm("ACGTACGTACGTACGTACGT"), float)
    except Exception:
        primer3_ok = False
    data_ok = os.path.isdir(DATA_DIR) and os.access(DATA_DIR, os.W_OK)
    email = str(getattr(ncbi.Entrez, "email", "") or "")
    return dict(ok=bool(primer3_ok and data_ok), version=app.version, commit=BUILD_COMMIT, booted=BOOT_TIME,
                primer3=primer3_ok,
                ncbi_email_set=bool(email and "example" not in email),
                ncbi_api_key=bool(getattr(ncbi.Entrez, "api_key", None)),
                ncbi_cache=getattr(ncbi, "_CACHE_ON", None),
                ncbi_cache_ttl_s=getattr(ncbi, "_CACHE_TTL", None),
                data_dir_writable=data_ok,
                routes=len([r for r in app.routes if getattr(r, "methods", None)]))

@app.get("/api/profiles")
def profiles():
    return {k: {"name": v["name"], "no_probe": v.get("no_probe", False),
                "notes": v.get("notes", "")} for k, v in P.PROFILES.items()}

@app.post("/api/qc")
def qc(r: OligoReq):
    s, notes, err = T.clean_seq(r.seq)
    if err:
        return JSONResponse({"error": err}, status_code=200)
    hdg37, hdg_an, htm = T.hairpin_full(s)
    sd37, sd_an, sdtm = T.self_dimer_full(s)
    out = dict(seq=s, length=len(s), gc=round(T.gc_percent(s), 1), tm=round(T.tm_acc(s), 1),
               hairpin_dg=round(hdg37, 2), hairpin_tm=round(htm, 0),
               self_dimer=round(sd37, 2), max_run=T.max_run(s),
               last5_gc=T.last5_gc(s), revcomp=T.revcomp(s),
               anneal_c=T.ANNEAL_C, hairpin_dg_anneal=round(hdg_an, 2),
               self_dimer_dg_anneal=round(sd_an, 2), self_dimer_tm=round(sdtm, 0))
    if T.has_degenerate(s):
        _tr = T.tm_range(s)
        if _tr["degenerate"] and not _tr["capped"] and _tr["min"] != _tr["max"]:
            out["tm_min"] = _tr["min"]; out["tm_max"] = _tr["max"]; out["tm_n"] = _tr["n"]
            notes.append("degenerate: Tm spans %s-%s C across %d resolutions" % (_tr["min"], _tr["max"], _tr["n"]))
        else:
            notes.append("contains degenerate bases — Tm/structure resolved to a representative base")
    if r.profile in P.PROFILES:
        out["lint"] = [dict(rule=n, status=st, detail=d)
                       for n, st, d in P.lint_oligo(s, r.role, P.PROFILES[r.profile])]
        if r.profile in ("idt_affinity", "parasite_mtdna") and r.role == "probe":
            out["lna_note"] = ("LNA probe: the Tm above is the DNA-backbone value. Each LNA base adds "
                               "~2-8 C; set LNA positions and confirm in IDT OligoAnalyzer.")
    if notes:
        out["note"] = " · ".join(notes)
    return out

@app.post("/api/pair")
def pair(r: PairReq):
    f, nf, ef = T.clean_seq(r.forward)
    rev, nr, er = T.clean_seq(r.reverse)
    if ef:
        return JSONResponse({"error": "forward primer: " + ef}, status_code=200)
    if er:
        return JSONResponse({"error": "reverse primer: " + er}, status_code=200)
    hx37, hx_an, hxtm = T.hetero_dimer_full(f, rev)
    out = dict(f_tm=round(T.tm_acc(f), 1), r_tm=round(T.tm_acc(rev), 1),
               pair_gap=round(abs(T.tm_acc(f) - T.tm_acc(rev)), 1),
               fxr=round(hx37, 2),
               f_self=round(T.self_dimer(f), 2), r_self=round(T.self_dimer(rev), 2),
               anneal_c=T.ANNEAL_C, fxr_anneal=round(hx_an, 2), fxr_tm=round(hxtm, 0),
               fxr_end_dg=round(min(T.end_stability(f, rev), T.end_stability(rev, f)), 2))
    for _seq, _k in ((f, "f"), (rev, "r")):
        _tr = T.tm_range(_seq)
        if _tr["degenerate"] and not _tr["capped"] and _tr["min"] != _tr["max"]:
            out[_k + "_tm_min"] = _tr["min"]; out[_k + "_tm_max"] = _tr["max"]
    if r.amplicon and r.profile in P.PROFILES:
        out["lint"] = [dict(rule=n, status=st, detail=d)
                       for n, st, d in P.lint_pair(f, rev, r.amplicon, P.PROFILES[r.profile])]
    if nf + nr:
        out["note"] = " · ".join(nf + nr)
    return out

class CondReq(BaseModel):
    mv_conc: Optional[float] = None
    dv_conc: Optional[float] = None
    dntp_conc: Optional[float] = None
    dna_conc: Optional[float] = None
    anneal_c: Optional[float] = None

@app.get("/api/conditions")
def get_conditions():
    return dict(T.COND, anneal_c=T.ANNEAL_C)

@app.post("/api/conditions")
def post_conditions(r: CondReq):
    return T.set_conditions(mv_conc=r.mv_conc, dv_conc=r.dv_conc, dntp_conc=r.dntp_conc, dna_conc=r.dna_conc, anneal_c=r.anneal_c)

@app.post("/api/matrix")
def matrix(r: MatrixReq):
    clean = {}
    for name, seq in r.oligos.items():
        cs, _, err = T.clean_seq(seq)
        if err:
            return JSONResponse({"error": "%s: %s" % (name, err)}, status_code=200)
        clean[name] = cs
    names = list(clean.keys())
    cells = []
    for i, a in enumerate(names):
        for b in names[i:]:
            dg = T.self_dimer(clean[a]) if a == b else T.hetero_dimer(clean[a], clean[b])
            cells.append(dict(a=a, b=b, dg=round(dg, 2)))
    return dict(names=names, cells=cells)

from oligoforge import autodesign as _AD_design
@app.post("/api/design")
def design(r: DesignReq):
    tmpl, notes, err = T.clean_seq(r.template)
    if err:
        return JSONResponse({"error": "template: " + err}, status_code=200)
    notes = list(notes or [])
    gc = 100.0 * sum(c in "GC" for c in tmpl) / max(1, len(tmpl))
    if (r.profile or "").lower() == "auto":
        order, _gc = _AD_design._auto_order(tmpl)
        if _gc < 40.0 and "parasite_lna" not in order:      # offer the LNA-core option for AT-rich targets
            order = order[:1] + ["parasite_lna"] + order[1:]
        a = None; used = None
        for pk in order:
            try:
                a = D.design_assay(tmpl, P.PROFILES[pk])
            except Exception:
                a = None
            if a:
                used = pk; break
        if not a:
            return JSONResponse({"error": "no clean assay found under any Auto chemistry (tried: %s; template GC %.0f%%). "
                                          "Try a longer/cleaner region." % (", ".join(order), gc)}, status_code=200)
        prof = P.PROFILES[used]
        notes.append("Auto-selected chemistry: %s (template GC %.0f%%)." % (prof["name"], gc))
    else:
        prof = P.PROFILES.get(r.profile, P.PROFILES["idt_taqman"])
        try:
            a = D.design_assay(tmpl, prof)
        except Exception as e:
            return JSONResponse({"error": f"design failed: {e}"}, status_code=200)
        if not a:
            if gc < 40:
                hint = " This template is AT-rich (GC %.0f%%) \u2014 switch to Auto, 'low-Tm TaqMan', or 'AT-rich + LNA probe'." % gc
            elif gc > 62:
                hint = " This template is GC-rich (GC %.0f%%) \u2014 switch to Auto or the 'GC-rich' profile." % gc
            else:
                hint = " Try Auto, a longer region, or a different chemistry."
            return JSONResponse({"error": "no clean assay found in this template under %s.%s" % (prof["name"], hint)},
                                status_code=200)
    pi = a.get("probe_info")
    return dict(forward=a["forward"], reverse=a["reverse"], probe=a["probe"],
                amplicon=a["amplicon"], amplicon_tm=a.get("amplicon_tm"),
                pair_tm_gap=round(a["pair_tm_gap"], 1),
                f_tm=round(a["f_tm"], 1), r_tm=round(a["r_tm"], 1),
                probe_tm=round(pi["tm"], 1) if pi else None,
                probe_offset=round(pi["offset"], 1) if pi else None,
                probe_hairpin=round(pi["hairpin_dg"], 2) if pi else None,
                probe_dimer_f=round(pi["dimer_f"], 2) if pi else None,
                probe_dimer_r=round(pi["dimer_r"], 2) if pi else None,
                gblock=a["gblock"], f_xy=a["f_xy"], r_xy=a["r_xy"], profile=prof["name"],
                note=(" \u00b7 ".join(notes) if notes else None))

@app.post("/api/fetch")
def fetch(r: FetchReq):
    _set_email(r.email, r.ncbi_key)
    try:
        if r.accession:
            recs = ncbi.fetch_accessions([x.strip() for x in r.accession.split(",") if x.strip()])
            return dict(records=[dict(id=x.id, desc=x.description, length=len(x.seq), seq=str(x.seq)) for x in recs])
        if r.gene and r.organism:
            recs, q = ncbi.fetch_isoforms(r.gene, r.organism)
            res = dict(query=q, records=[dict(id=x.id, desc=x.description, length=len(x.seq), seq=str(x.seq)) for x in recs])
            if r.isoform_common and len(recs) > 1:
                res["common_region"] = ncbi.common_region([x.seq for x in recs])
            return res
        return JSONResponse({"error": "provide an accession, or gene + organism"}, status_code=200)
    except Exception as e:
        return JSONResponse({"error": f"NCBI fetch failed: {e}"}, status_code=200)

@app.post("/api/intron")
def intron(r: IntronReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return SP.intron_check(r.gene, r.organism, r.amp_start, r.amp_end, r.mrna_acc, r.forward, r.reverse)
    except Exception as e:
        return JSONResponse({"error": f"intron check failed: {e}"}, status_code=200)

@app.post("/api/blast")
def blast(r: BlastReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return SP.blast_summary(r.seq, mode=r.mode, db=r.db, db_path=r.db_path,
                                organism=r.organism, top=r.top)
    except Exception as e:
        return JSONResponse({"error": f"BLAST failed: {e}"}, status_code=200)


# ============ stage-5 enhancements ============
import json, re
from oligoforge import quant as Q, orders as O

PANELS_DIR = os.path.join(DATA_DIR, "panels")
os.makedirs(PANELS_DIR, exist_ok=True)
PROJECTS_DIR = os.path.join(DATA_DIR, "projects")
os.makedirs(PROJECTS_DIR, exist_ok=True)

class ProjectSaveReq(BaseModel):
    name: str
    assays: List[dict]

class ProjectNameReq(BaseModel):
    name: str

@app.post("/api/project/save")
def project_save(r: ProjectSaveReq):
    import datetime
    if not r.name.strip():
        return JSONResponse({"error": "project name required"}, status_code=200)
    fn = os.path.join(PROJECTS_DIR, _safe(r.name) + ".json")
    json.dump({"name": r.name, "assays": r.assays,
               "saved": datetime.datetime.now().strftime("%Y-%m-%d %H:%M")},
              open(fn, "w"), indent=1)
    return dict(saved=r.name, n=len(r.assays))

def _project_list():
    out = []
    if os.path.isdir(PROJECTS_DIR):
        for f in sorted(os.listdir(PROJECTS_DIR)):
            if f.endswith(".json"):
                try:
                    d = json.load(open(os.path.join(PROJECTS_DIR, f)))
                    out.append(dict(name=d.get("name", f[:-5]),
                                    n=len(d.get("assays", [])), saved=d.get("saved", "")))
                except Exception:
                    pass
    return dict(projects=out)


@app.get("/api/project/list")
def project_list_get():
    return _project_list()


@app.post("/api/project/list")
def project_list_post():
    return _project_list()

@app.post("/api/project/load")
def project_load(r: ProjectNameReq):
    fn = os.path.join(PROJECTS_DIR, _safe(r.name) + ".json")
    if not os.path.exists(fn):
        return JSONResponse({"error": "no project named " + r.name}, status_code=200)
    return json.load(open(fn))

@app.post("/api/project/delete")
def project_delete(r: ProjectNameReq):
    fn = os.path.join(PROJECTS_DIR, _safe(r.name) + ".json")
    if os.path.exists(fn):
        os.remove(fn)
        return dict(deleted=r.name)
    return JSONResponse({"error": "no project named " + r.name}, status_code=200)
def _safe(n): return re.sub(r"[^A-Za-z0-9_.-]", "_", n or "")[:80]


class CopiesReq(BaseModel):
    ng_per_ul: float; length_bp: int; factor: float = 10.0; points: int = 6; molecule_type: str = "dsDNA"

class BatchItem(BaseModel):
    name: str; template: str; profile: str = "idt_taqman"
class BatchReq(BaseModel):
    items: List[BatchItem]

class PanelSaveReq(BaseModel):
    name: str; oligos: Dict[str, str]
class PanelLoadReq(BaseModel):
    name: str

class OrderEntry(BaseModel):
    name: str; seq: str; kind: str = "primer"
class GBlockEntry(BaseModel):
    name: str; seq: str
class OrderReq(BaseModel):
    oligos: List[OrderEntry] = []; gblocks: List[GBlockEntry] = []

class PairSpecReq(BaseModel):
    forward: str; reverse: str; organism: Optional[str] = None
    mode: str = "remote"; db: str = "nt"; db_path: Optional[str] = None; email: Optional[str] = None; ncbi_key: Optional[str] = None


@app.post("/api/copies")
def copies(r: CopiesReq):
    if r.length_bp < 1 or r.ng_per_ul <= 0:
        return JSONResponse({"error": "enter a positive concentration and a length in bp"}, status_code=200)
    if r.factor <= 1:
        return JSONResponse({"error": "dilution factor must be greater than 1 (e.g. 10 for a 10-fold series)"}, status_code=200)
    pts = max(1, min(int(r.points), 40))
    mt = r.molecule_type if r.molecule_type in Q.MW_PER_UNIT else "dsDNA"
    c = Q.copies_per_ul(r.ng_per_ul, r.length_bp, mt)
    series = Q.dilution_series(c, r.factor, pts)
    return dict(copies_per_ul=c, ng_per_ul=r.ng_per_ul, length_bp=r.length_bp, molecule_type=mt,
                series=[dict(point=i, copies_per_ul=v) for i, v in enumerate(series)])

@app.post("/api/batch_design")
def batch_design(r: BatchReq):
    out = []
    for it in r.items:
        prof = P.PROFILES.get(it.profile, P.PROFILES["idt_taqman"])
        try:
            a = D.design_assay(it.template.upper().strip(), prof)
        except Exception:
            a = None
        if a:
            pi = a.get("probe_info")
            out.append(dict(name=it.name, ok=True, forward=a["forward"], reverse=a["reverse"],
                            probe=a["probe"], amplicon=a["amplicon"], f_tm=round(a["f_tm"], 1),
                            r_tm=round(a["r_tm"], 1),
                            probe_tm=round(pi["tm"], 1) if pi else None, gblock=a["gblock"]))
        else:
            out.append(dict(name=it.name, ok=False, error="no clean assay under that profile"))
    return dict(results=out)

@app.post("/api/panel/save")
def panel_save(r: PanelSaveReq):
    fn = os.path.join(PANELS_DIR, _safe(r.name) + ".json")
    json.dump({"name": r.name, "oligos": r.oligos}, open(fn, "w"), indent=1)
    return dict(saved=r.name)

@app.get("/api/panel/list")
def panel_list():
    return dict(panels=sorted(f[:-5] for f in os.listdir(PANELS_DIR) if f.endswith(".json")))


@app.post("/api/factory_reset")
def factory_reset():
    """Delete every saved panel and project on the server -- the server half of a factory reset.
    Only *.json inside the two managed directories are touched (names come from os.listdir, so no
    path traversal). The browser clears its own localStorage/session separately."""
    removed = {"panels": 0, "projects": 0}
    for key, d in (("panels", PANELS_DIR), ("projects", PROJECTS_DIR)):
        try:
            names = os.listdir(d)
        except OSError:
            names = []
        for fn in names:
            if fn.endswith(".json"):
                try:
                    os.remove(os.path.join(d, fn)); removed[key] += 1
                except OSError:
                    pass
    return dict(ok=True, panels=removed["panels"], projects=removed["projects"])

@app.post("/api/panel/load")
def panel_load(r: PanelLoadReq):
    fn = os.path.join(PANELS_DIR, _safe(r.name) + ".json")
    if not os.path.exists(fn):
        return JSONResponse({"error": "panel not found"}, status_code=200)
    return json.load(open(fn))

@app.post("/api/order_csv")
def order_csv(r: OrderReq):
    oligos = []
    for e in r.oligos:
        if e.kind == "probe_lna":
            # LNA / Affinity Plus probes use IDT '+N' notation; validate the DNA
            # backbone but keep the +N positions in the ordered sequence.
            raw = "".join(e.seq.upper().split())
            backbone, _, _ = T.strip_lna(raw)
            _, _, err = T.clean_seq(backbone)
            seq_for_order = raw
        else:
            seq_for_order, _, err = T.clean_seq(e.seq)
        if err:
            return JSONResponse({"error": "%s: %s" % (e.name or "oligo", err)}, status_code=200)
        d = e.dict(); d["seq"] = seq_for_order; oligos.append(d)
    gbs = []
    for g in r.gblocks:
        cs, _, err = T.clean_seq(g.seq)
        if err:
            return JSONResponse({"error": "%s: %s" % (g.name or "gBlock", err)}, status_code=200)
        d = g.dict(); d["seq"] = cs; gbs.append(d)
    return dict(oligo_csv=O.oligo_csv(oligos) if oligos else "",
                gblock_fasta=O.gblock_fasta(gbs) if gbs else "")

@app.post("/api/pair_specificity")
def pair_specificity(r: PairSpecReq):
    _set_email(r.email, r.ncbi_key)
    def one(seq):
        if r.mode == "local":
            return SP.blast_local(seq, r.db_path or "")
        return SP.blast_remote(seq, organism=r.organism, db=r.db)
    try:
        return dict(forward=one(r.forward.upper().strip()), reverse=one(r.reverse.upper().strip()))
    except Exception as e:
        return JSONResponse({"error": f"specificity check failed: {e}"}, status_code=200)


# ============ stage-6: conservation, in-silico PCR, LNA Tm ============
from oligoforge import conservation as CONS, thermo as T

class ConsReq(BaseModel):
    oligos: Dict[str, str]; targets: List[str]
    off_targets: List[str] = []; min_ident: float = 0.6
class EpcrReq(BaseModel):
    forward: str; reverse: str; mode: str = "remote"; db: str = "nt"
    db_path: Optional[str] = None; organism: Optional[str] = None
    min_product: int = 40; max_product: int = 3000; email: Optional[str] = None; ncbi_key: Optional[str] = None
class LnaReq(BaseModel):
    seq: str; n_lna: Optional[int] = None; snp_pos: Optional[int] = None

@app.post("/api/conservation")
def api_conservation(r: ConsReq):
    tg = [t for t in r.targets if t.strip()]
    off = [o for o in r.off_targets if o.strip()] or None
    return CONS.analyze(r.oligos, tg, off, r.min_ident)

@app.post("/api/epcr")
def api_epcr(r: EpcrReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return SP.in_silico_pcr(r.forward, r.reverse, r.mode, r.db, r.db_path,
                                r.organism, r.min_product, r.max_product)
    except Exception as e:
        return JSONResponse({"error": f"in-silico PCR failed: {e}"}, status_code=200)

class AssaySpecReq(BaseModel):
    forward: str; reverse: str; probe: Optional[str] = None
    mode: str = "remote"; db: str = "nt"; db_path: Optional[str] = None
    organism: Optional[str] = None
    min_product: int = 40; max_product: int = 3000
    email: Optional[str] = None; ncbi_key: Optional[str] = None
@app.post("/api/assay_specificity")
def api_assay_specificity(r: AssaySpecReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return SP.assay_specificity(r.forward, r.reverse, r.probe, r.mode, r.db, r.db_path,
                                    r.organism, r.min_product, r.max_product)
    except Exception as e:
        return JSONResponse({"error": f"assay specificity failed: {e}"}, status_code=200)

@app.post("/api/lna_tm")
def api_lna_tm(r: LnaReq):
    return T.tm_lna(r.seq, r.n_lna)

@app.post("/api/lna_suggest")
def api_lna_suggest(r: LnaReq):
    return T.suggest_lna(r.seq, snp_pos=r.snp_pos, max_lna=r.n_lna)

class GeneLookupReq(BaseModel):
    gene: str; organism: Optional[str] = None; email: Optional[str] = None; ncbi_key: Optional[str] = None

@app.post("/api/gene_lookup")
def api_gene_lookup(r: GeneLookupReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return ncbi.gene_lookup(r.gene, r.organism)
    except Exception as e:
        return JSONResponse({"error": f"gene lookup failed: {e}"}, status_code=200)

from oligoforge import refgenes as RG
class RefGenesReq(BaseModel):
    text: str
from oligoforge import report as RPT, multiplex as MX, refmarkers as RM, markerscan as MS
class ReportReq(BaseModel):
    panel: List[dict]; meta: Optional[dict] = None
@app.post("/api/report")
def api_report(r: ReportReq):
    try:
        return RPT.build(r.panel, r.meta)
    except Exception as e:
        return JSONResponse({"error": "report failed: %s" % e}, status_code=200)

from oligoforge import rdml as RDML
class RdmlReq(BaseModel):
    panel: List[dict]; meta: Optional[dict] = None
@app.post("/api/rdml")
def api_rdml(r: RdmlReq):
    """Export the panel as RDML 1.2 (machine-readable assay definitions for qPCR software)."""
    try:
        return RDML.build(r.panel, r.meta)
    except Exception as e:
        return JSONResponse({"error": "RDML export failed: %s" % e}, status_code=200)

class MultiplexReq(BaseModel):
    assays: List[dict]; dimer_threshold: float = -6.0
class MarkerReq(BaseModel):
    organism: str; exclude: Optional[str] = None; intent: Optional[str] = "any"
    email: Optional[str] = None; ncbi_key: Optional[str] = None
@app.post("/api/scan_markers")
def api_scan_markers(r: MarkerReq):
    _set_email(r.email, r.ncbi_key)
    org = (r.organism or "").strip()
    if not org:
        return JSONResponse({"error": 'enter an organism (a genus or species, e.g. "Plasmodium" or "Aphelocoma coerulescens")'}, status_code=200)
    try:
        base_info = RM.suggest(org, r.exclude, r.intent)
        sc = MS.scan(base_info["resolved"], base_info["markers"], r.exclude)
        by = {x["gene"]: x for x in sc["results"]}
        for m in base_info["markers"]:
            d = by.get(m["gene"])
            if d:
                m.update(count=d["count"], off_count=d.get("off_count"), n_seqs=d["n"],
                         median_len=d["median_len"], min_len=d["min_len"], max_len=d["max_len"], scanned=True)
        base_info["scanned"] = True; base_info["n_scanned"] = len(sc["results"])
        return base_info
    except Exception as e:
        return JSONResponse({"error": "scan_markers failed: %s" % e}, status_code=200)

@app.post("/api/suggest_genes")
def api_suggest_genes(r: MarkerReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return RM.suggest(r.organism, r.exclude, r.intent)
    except Exception as e:
        return JSONResponse({"error": "suggest_genes failed: %s" % e}, status_code=200)

@app.post("/api/multiplex")
def api_multiplex(r: MultiplexReq):
    try:
        return MX.check(r.assays, r.dimer_threshold)
    except Exception as e:
        return JSONResponse({"error": "multiplex failed: %s" % e}, status_code=200)

@app.post("/api/refgenes")
def api_refgenes(r: RefGenesReq):
    try:
        data = RG.parse_table(r.text)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=200)
    if len(data) < 2:
        return JSONResponse({"error": "need at least 2 candidate reference genes (rows)"}, status_code=200)
    try:
        return RG.analyze(data)
    except Exception as e:
        return JSONResponse({"error": "stability analysis failed: %s" % e}, status_code=200)

from oligoforge import ncbi as NC
class FetchNucReq(BaseModel):
    query: str; n: int = 10; email: Optional[str] = None; ncbi_key: Optional[str] = None
@app.post("/api/fetch_nuc")
def api_fetch_nuc(r: FetchNucReq):
    _set_email(r.email, r.ncbi_key)
    try:
        recs = NC.search_fetch_fasta(r.query, min(r.n, 30))
        return dict(n=len(recs), fasta="\n".join(f">{d}\n{s}" for d, s in recs),
                    sequences=[s for _, s in recs])
    except Exception as e:
        return JSONResponse({"error": f"fetch failed: {e}"}, status_code=200)

from oligoforge import quant as QT
class CqPoint(BaseModel):
    quantity: float; cq: Optional[float] = None
class StdCurveReq(BaseModel):
    points: List[CqPoint]
@app.post("/api/standard_curve")
def api_standard_curve(r: StdCurveReq):
    return QT.standard_curve([(p.quantity, p.cq) for p in r.points])

from oligoforge import cq as CQ
class CqReq(BaseModel):
    fluor: List[float]
    cycles: Optional[List[float]] = None
    threshold: Optional[float] = None
    sd_mult: float = 10.0
    baseline_start: int = 3
    baseline_end: int = 15
@app.post("/api/cq")
def api_cq(r: CqReq):
    return CQ.analyze(r.fluor, cycles=r.cycles, threshold=r.threshold,
                      sd_mult=r.sd_mult, baseline=(r.baseline_start, r.baseline_end))

from oligoforge import expression as EXP
class ExpressionReq(BaseModel):
    csv: Optional[str] = None
    samples: Optional[List[dict]] = None
    reference_genes: List[str]
    control_group: str
    efficiencies: Optional[Dict[str, float]] = None
@app.post("/api/expression")
def api_expression(r: ExpressionReq):
    try:
        samples = r.samples or EXP.parse_table(r.csv or "")
        return EXP.analyze(samples, r.reference_genes, r.control_group, r.efficiencies)
    except Exception as e:
        return JSONResponse({"error": "expression analysis failed: %s" % e}, status_code=200)

from oligoforge import melt as MELT
class MeltReq(BaseModel):
    fluor: List[float]
    temps: List[float]
    min_prominence_frac: float = 0.10
    predicted_tm: Optional[float] = None
@app.post("/api/melt")
def api_melt(r: MeltReq):
    return MELT.analyze(r.fluor, r.temps, min_prominence_frac=r.min_prominence_frac,
                        predicted_tm=r.predicted_tm)

from oligoforge import autodesign as AD
class AutoDesignReq(BaseModel):
    target_query: str; profile: str = "auto"
    off_query: Optional[str] = None; n_fetch: int = 20
    min_ident: float = 0.6
    run_blast: bool = False; blast_mode: str = "remote"
    blast_db: str = "nt"; blast_db_path: Optional[str] = None; organism: Optional[str] = None
    email: Optional[str] = None; ncbi_key: Optional[str] = None; prefer_junction: bool = False; nested: bool = False
@app.post("/api/autodesign")
def api_autodesign(r: AutoDesignReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return AD.design_from_query(r.target_query, r.profile, r.off_query, min(r.n_fetch, 30),
                                    r.min_ident, r.run_blast, r.blast_mode, r.blast_db,
                                    r.blast_db_path, r.organism, prefer_junction=r.prefer_junction,
                                    nested=r.nested)
    except Exception as e:
        return JSONResponse({"error": f"autodesign failed: {e}"}, status_code=200)


# ============ isolate panel validation (inclusivity / exclusivity in-silico PCR) ============
class IsolateGenomesReq(BaseModel):
    query: str; retmax: int = 40
    email: Optional[str] = None; ncbi_key: Optional[str] = None

@app.post("/api/isolate_genomes")
def api_isolate_genomes(r: IsolateGenomesReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return {"genomes": ncbi.search_genomes(r.query, r.retmax)}
    except Exception as e:
        return JSONResponse({"error": f"genome search failed: {e}"}, status_code=200)

class IsolateCheckReq(BaseModel):
    forward: str; reverse: str; probe: str = ""
    accessions: List[str]; role: str = "target"
    max_mm: int = 5; min_product: int = 40; max_product: int = 3000; min_probe_ident: float = 85.0
    email: Optional[str] = None; ncbi_key: Optional[str] = None

@app.post("/api/isolate_check")
def api_isolate_check(r: IsolateCheckReq):
    """Run in-silico PCR of one assay against each accession (one genome at a time, freed after).
    The frontend chunks the panel into small batches so no single request fetches many genomes."""
    _set_email(r.email, r.ncbi_key)
    accs = [a.strip() for a in (r.accessions or []) if a.strip()][:6]    # hard cap; frontend paginates
    maxp = max(int(r.min_product), min(int(r.max_product), 50000))       # qPCR amplicons are small; bound the
    mm = max(0, min(int(r.max_mm), 10))                                  # probe scan and primer search so a
    pid = max(0.0, min(float(r.min_probe_ident), 100.0))                 # pathological request cannot stall a worker
    out = []
    for acc in accs:
        try:
            title, seq = ncbi.fetch_one(acc)
            if not seq:
                out.append(dict(acc=acc, title=title, slen=0, role=r.role, amplifies=None,
                                error="no sequence (accession may be an assembly master / empty record)"))
                continue
            res = ISO.amplify(r.forward, r.reverse, r.probe, seq, max_mm=mm,
                              min_product=r.min_product, max_product=maxp,
                              min_probe_ident=pid)
            res.update(acc=acc, title=title, slen=len(seq), role=r.role, error=None)
            out.append(res)
        except Exception as e:
            out.append(dict(acc=acc, title="", slen=0, role=r.role, amplifies=None, error=str(e)))
    return {"results": out}


# ---- SnapGene-style viewer: strict-rules design on a user-picked sequence, with base coordinates ----
class ViewerDesignReq(BaseModel):
    sequence: str
    tm_min: float = 59.0
    tm_max: float = 64.5
    gc_min: float = 35.0
    gc_max: float = 65.0
    amp_min: int = 70
    amp_max: int = 150
    len_min: int = 18
    len_max: int = 24
    probe: bool = True
    probe_offset_min: float = 5.0
    probe_offset_max: float = 10.5
    n: int = 5


@app.post("/api/viewer_design")
def viewer_design(r: ViewerDesignReq):
    seq, notes, err = T.clean_seq(r.sequence)
    if err:
        return JSONResponse({"error": "sequence: " + err}, status_code=200)
    if len(seq) < 50:
        return JSONResponse({"error": "sequence too short to design on (need at least 50 nt)"}, status_code=200)
    # strict config from the TaqMan base, with the user's multiplex rules clamped to sane ranges
    c = dict(P.PROFILES["idt_taqman"])
    tmlo, tmhi = sorted((float(r.tm_min), float(r.tm_max)))
    gclo, gchi = sorted((float(r.gc_min), float(r.gc_max)))
    amlo, amhi = sorted((int(r.amp_min), int(r.amp_max)))
    lnlo, lnhi = sorted((int(r.len_min), int(r.len_max)))
    c.update(tm_min=tmlo, tm_max=tmhi, tm_opt=(tmlo + tmhi) / 2.0,
             gc_min=max(0.0, gclo), gc_max=min(100.0, gchi),
             amp_min=max(40, amlo), amp_max=min(2000, amhi),
             len_min=max(12, min(lnlo, 40)), len_max=max(12, min(lnhi, 40)))
    if r.probe:
        polo, pohi = sorted((float(r.probe_offset_min), float(r.probe_offset_max)))
        c.update(no_probe=False, probe_offset_min=polo, probe_offset_max=pohi)
    else:
        c["no_probe"] = True
    n = max(1, min(int(r.n), 8))
    try:
        cands = D.design_candidates(seq, c, n=n)
    except Exception as e:
        return JSONResponse({"error": f"design failed: {e}"}, status_code=200)

    def _gc(s):
        return round(T.gc_percent(s), 1) if s else None

    def _inw(tm):
        return bool(tmlo <= tm <= tmhi)

    out = []
    for a in cands:
        pi = a.get("probe_info")
        out.append(dict(
            forward=a["forward"], reverse=a["reverse"], probe=a.get("probe"),
            amplicon=a["amplicon"], amplicon_tm=a.get("amplicon_tm"),
            f_xy=a["f_xy"], r_xy=a["r_xy"], probe_xy=a.get("probe_xy"), amplicon_xy=a["amplicon_xy"],
            f_tm=round(a["f_tm"], 1), r_tm=round(a["r_tm"], 1), pair_tm_gap=round(a["pair_tm_gap"], 1),
            f_gc=_gc(a["forward"]), r_gc=_gc(a["reverse"]), probe_gc=_gc(a.get("probe")),
            probe_tm=round(pi["tm"], 1) if pi else None,
            probe_offset=round(pi["offset"], 1) if pi else None,
            probe_strand=(pi.get("strand") if pi else None),
            f_in=_inw(a["f_tm"]), r_in=_inw(a["r_tm"]),
            gblock=a.get("gblock")))
    return dict(seq_len=len(seq), tm_window=[tmlo, tmhi], gc_window=[gclo, gchi],
                candidates=out, note=(" · ".join(notes) if notes else None))
