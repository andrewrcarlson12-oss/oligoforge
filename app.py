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

from oligoforge import thermo as T, design as D, profiles as P, ncbi, specificity as SP

app = FastAPI(title="OligoForge", version="1.11.5")
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


def _set_email(email=None, key=None):
    e = email or os.environ.get("OLIGOFORGE_EMAIL")
    if e:
        ncbi.Entrez.email = SP.Entrez.email = e
    k = key or os.environ.get("OLIGOFORGE_NCBI_KEY")
    ncbi.Entrez.api_key = SP.Entrez.api_key = (k or None)


_set_email(None)


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
    return FileResponse(os.path.join(RES_DIR, "static", "index.html"))

@app.get("/api/profiles")
def profiles():
    return {k: {"name": v["name"], "no_probe": v.get("no_probe", False),
                "notes": v.get("notes", "")} for k, v in P.PROFILES.items()}

@app.post("/api/qc")
def qc(r: OligoReq):
    s, notes, err = T.clean_seq(r.seq)
    if err:
        return JSONResponse({"error": err}, status_code=200)
    hdg, htm = T.hairpin(s)
    out = dict(seq=s, length=len(s), gc=round(T.gc_percent(s), 1), tm=round(T.tm(s), 1),
               hairpin_dg=round(hdg, 2), hairpin_tm=round(htm, 0),
               self_dimer=round(T.self_dimer(s), 2), max_run=T.max_run(s),
               last5_gc=T.last5_gc(s), revcomp=T.revcomp(s))
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
    out = dict(f_tm=round(T.tm(f), 1), r_tm=round(T.tm(rev), 1),
               pair_gap=round(abs(T.tm(f) - T.tm(rev)), 1),
               fxr=round(T.hetero_dimer(f, rev), 2),
               f_self=round(T.self_dimer(f), 2), r_self=round(T.self_dimer(rev), 2))
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

@app.get("/api/conditions")
def get_conditions():
    return dict(T.COND)

@app.post("/api/conditions")
def post_conditions(r: CondReq):
    return T.set_conditions(mv_conc=r.mv_conc, dv_conc=r.dv_conc, dntp_conc=r.dntp_conc, dna_conc=r.dna_conc)

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

@app.post("/api/design")
def design(r: DesignReq):
    tmpl, notes, err = T.clean_seq(r.template)
    if err:
        return JSONResponse({"error": "template: " + err}, status_code=200)
    prof = P.PROFILES.get(r.profile, P.PROFILES["idt_taqman"])
    try:
        a = D.design_assay(tmpl, prof)
    except Exception as e:
        return JSONResponse({"error": f"design failed: {e}"}, status_code=200)
    if not a:
        return JSONResponse({"error": "no clean assay found in this template under that profile"}, status_code=200)
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
                note=(" · ".join(notes) if notes else None))

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

@app.api_route("/api/project/list", methods=["GET", "POST"])
def project_list():
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

class MultiplexReq(BaseModel):
    assays: List[dict]; dimer_threshold: float = -9.0
class MarkerReq(BaseModel):
    organism: str; exclude: Optional[str] = None; intent: Optional[str] = "any"
    email: Optional[str] = None; ncbi_key: Optional[str] = None
@app.post("/api/scan_markers")
def api_scan_markers(r: MarkerReq):
    _set_email(r.email, r.ncbi_key)
    try:
        base_info = RM.suggest(r.organism, r.exclude, r.intent)
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
