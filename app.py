"""OligoForge backend — FastAPI server exposing the engine to the browser cockpit.

Run:  uvicorn app:app --reload --port 8111
Then open http://127.0.0.1:8111
Set your NCBI email once (env var or the field in the UI):  export OLIGOFORGE_EMAIL=you@uni.edu
"""
import os, sys
from typing import Dict, List, Optional
from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.exceptions import RequestValidationError
from pydantic import BaseModel

from oligoforge import thermo as T, design as D, profiles as P, ncbi, specificity as SP, isolates as ISO, multiplex as MX, structure as STR, nn as NN
from oligoforge import ranking_profiles as RPROF, manual_design as MDS, assay_rescue as ARES, experimental_feedback as EFB, run_compare as RCOMP

app = FastAPI(title="OligoForge", version="1.34.0")
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


def _env_flag(name, default=False):
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


# The local desktop/server workflow may use process-wide reaction settings and JSON
# storage.  A public Render instance is multi-user, so those shared mutable features
# are disabled unless the operator deliberately opts in behind authentication.
HOSTED_MODE = _env_flag("OLIGOFORGE_HOSTED", False)
ALLOW_SERVER_STORAGE = _env_flag("OLIGOFORGE_ALLOW_SERVER_STORAGE", not HOSTED_MODE)
ALLOW_SHARED_CONDITIONS = _env_flag("OLIGOFORGE_ALLOW_SHARED_CONDITIONS", not HOSTED_MODE)


def _shared_feature_disabled(feature):
    return JSONResponse(
        {"error": "%s is disabled on this multi-user hosted deployment; use browser-local storage "
                  "or run a private instance" % feature},
        status_code=403,
    )


def _api_failure(area, exc, status_code=200):
    """Log full diagnostics server-side without reflecting secrets, paths or stack details publicly."""
    log.exception("%s failed", area)
    detail = "%s failed" % area
    if not HOSTED_MODE:
        detail += ": %s" % exc
    return JSONResponse({"error": detail}, status_code=status_code)


def _local_blast_denied(mode=None, db_path=None):
    if HOSTED_MODE and (str(mode or "").lower() == "local" or bool(db_path)):
        return JSONResponse({"error": "local BLAST database paths are disabled on hosted deployments"},
                            status_code=403)
    return None


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


MAX_REQUEST_BYTES = int(os.environ.get("OLIGOFORGE_MAX_REQUEST_BYTES", 5 * 1024 * 1024))


def _security_headers(resp):
    resp.headers.setdefault("X-Content-Type-Options", "nosniff")
    resp.headers.setdefault("X-Frame-Options", "DENY")
    resp.headers.setdefault("Referrer-Policy", "no-referrer")
    resp.headers.setdefault("Permissions-Policy", "camera=(), microphone=(), geolocation=()")
    resp.headers.setdefault("Content-Security-Policy",
        "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
        "font-src 'self' https://fonts.gstatic.com; img-src 'self' data:; connect-src 'self'; object-src 'none'; "
        "base-uri 'none'; frame-ancestors 'none'")
    return resp


@app.exception_handler(RequestValidationError)
async def _validation_error(_request: Request, exc: RequestValidationError):
    # Do not echo rejected payload values (which may contain huge sequences, API keys or local paths).
    details = [{"loc": list(e.get("loc", ())), "msg": e.get("msg", "invalid value"),
                "type": e.get("type", "validation_error")} for e in exc.errors()]
    return _security_headers(JSONResponse({"error": "request validation failed", "details": details},
                                         status_code=422))


@app.middleware("http")
async def _log_requests(request, call_next):
    t0 = _time.perf_counter()
    raw_len = request.headers.get("content-length")
    try:
        if raw_len is not None and int(raw_len) > MAX_REQUEST_BYTES:
            return _security_headers(JSONResponse({"error": "request body too large"}, status_code=413))
    except ValueError:
        return _security_headers(JSONResponse({"error": "invalid Content-Length header"}, status_code=400))
    try:
        resp = await call_next(request)
    except Exception:
        log.exception("unhandled error %s %s", request.method, request.url.path)
        resp = JSONResponse({"error": "internal server error"}, status_code=500)
    log.info("%s %s -> %s %.0fms", request.method, request.url.path,
             resp.status_code, (_time.perf_counter() - t0) * 1000.0)
    return _security_headers(resp)


# ---------- request models ----------
class OligoReq(BaseModel):
    seq: str; role: str = "primer"; profile: str = "idt_taqman"

class PairReq(BaseModel):
    forward: str; reverse: str; amplicon: Optional[int] = None; profile: str = "idt_taqman"

class MatrixReq(BaseModel):
    oligos: Dict[str, str]

class DesignReq(BaseModel):
    template: str; profile: str = "idt_taqman"; off_targets: Optional[str] = None
    panel: Optional[List[dict]] = None
    objective: str = "balanced"

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
                data_dir_writable=data_ok, hosted_mode=HOSTED_MODE,
                server_storage_enabled=ALLOW_SERVER_STORAGE,
                shared_conditions_enabled=ALLOW_SHARED_CONDITIONS,
                routes=len([r for r in app.routes if getattr(r, "methods", None)]))

@app.get("/api/profiles")
def profiles():
    return {k: {"name": v["name"], "no_probe": v.get("no_probe", False),
                "notes": v.get("notes", "")} for k, v in P.PROFILES.items()}

@app.get("/api/ranking-profiles")
def ranking_profiles():
    return RPROF.public_profiles()


def _require_objective(name):
    key = str(name or "balanced").strip().lower()
    if key not in RPROF.OBJECTIVE_PROFILES:
        raise ValueError("unknown assay objective: %s" % key)
    return key


def _require_profile(key):
    if key not in P.PROFILES:
        raise ValueError("unknown chemistry profile: %s" % key)
    return P.PROFILES[key]


def _bounded_sequence_list(rows, label, max_records=50):
    rows = list(rows or [])
    if len(rows) > max_records:
        raise ValueError("%s contains %d records; limit %d" % (label, len(rows), max_records))
    out = []
    for i, row in enumerate(rows):
        if len(row or "") > T.MAX_TEMPLATE_LEN:
            raise ValueError("%s record %d exceeds %d nt" % (label, i + 1, T.MAX_TEMPLATE_LEN))
        c, _n, err = T.clean_seq(row or "")
        if err or not c:
            raise ValueError("%s record %d is invalid: %s" % (label, i + 1, err or "empty"))
        out.append(c)
    return out

@app.post("/api/qc")
def qc(r: OligoReq):
    raw = r.seq or ""
    if len(raw) > T.MAX_OLIGO_LEN:
        return JSONResponse({"error": "oligo too long (%d nt; limit %d). A primer/probe is <=60 nt — "
                                      "this looks like a template. Use the Design tab for a full sequence."
                                      % (len(raw), T.MAX_OLIGO_LEN)}, status_code=200)
    mod_bits = []
    if "/" in raw: mod_bits.append("5'/3' or internal modification codes")
    if "+" in raw: mod_bits.append("LNA (+) bases")
    if "*" in raw: mod_bits.append("phosphorothioate (*) linkages")
    s, notes, err = T.clean_seq(T.strip_mods(raw))   # accept a pasted IDT order string / LNA oligo
    if err:
        return JSONResponse({"error": err}, status_code=200)
    if mod_bits:
        notes.append("modification notation stripped to score the DNA backbone ("
                     + ", ".join(mod_bits) + "); Tm / \u0394G shown are for the unmodified sequence")
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
        if r.profile in ("idt_affinity", "parasite_lna") and r.role == "probe" and "+" not in raw:
            out["lna_note"] = ("The displayed Tm is for the unmodified DNA backbone. LNA effects are "
                               "sequence- and position-dependent; enter explicit +N positions for the "
                               "McTigue-model estimate and confirm the final order with the vendor.")
    if "+" in raw:
        _nnl = NN.params_lna(raw)   # context-specific McTigue 2004 LNA estimate
        if _nnl:
            out["nn_lna"] = _nnl
            out["lna_note"] = ("LNA estimate calculated from the explicit +N positions using the "
                               "McTigue nearest-neighbour model. It is a computational estimate, not "
                               "a vendor certificate or empirical hybridization measurement.")
        else:
            out["lna_note"] = ("LNA notation was detected but could not be fully parameterized. Review "
                               "the +N syntax and confirm the sequence with the vendor calculator.")
    _nn = NN.params(s)   # transparent NN thermodynamics at the reaction conditions (None if degenerate)
    if _nn:
        out["nn"] = _nn
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
    if not ALLOW_SHARED_CONDITIONS:
        return _shared_feature_disabled("process-wide reaction-condition changes")
    return T.set_conditions(mv_conc=r.mv_conc, dv_conc=r.dv_conc, dntp_conc=r.dntp_conc,
                            dna_conc=r.dna_conc, anneal_c=r.anneal_c)

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
def _parse_fasta_targets(text):
    """Accept either raw sequence or multi-record FASTA in one box. Returns (cleaned_seqs>60nt, notes).
    Each record is run through the same clean_seq the engine uses, so a pasted alignment / RNA is corrected."""
    if not text or not text.strip():
        return [], []
    has_header = ">" in text
    recs = []
    if has_header:
        cur = []
        for line in text.splitlines():
            if line.startswith(">"):
                if cur:
                    recs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip())
        if cur:
            recs.append("".join(cur))
    else:
        recs = [text]
    out = []
    for raw in recs:
        c, _n, err = T.clean_seq(raw)
        if c and not err and len(c) > 60:
            out.append(c)
    notes = []
    if has_header and len(out) < len(recs):
        notes.append("%d of %d FASTA record(s) usable (others too short / unreadable)." % (len(out), len(recs)))
    return out, notes


def _panel_fit(cand_oligos, amplicon_tm, is_sybr, panel, dimer_threshold=-6.0, amp_tm_gap=2.0):
    """Cross-check ONE Design candidate against the existing workbench panel: heterodimers between the
    candidate's oligos and every panel oligo (same thresholds + 3'-engagement test as the Multiplex tab),
    plus SYBR amplicon-Tm proximity. Panel-internal conflicts are the Multiplex tab's job; this reports
    only what dropping the NEW assay into the panel would add."""
    cross = []
    for cn, cs in cand_oligos:
        cs = (cs or "").upper()
        if not cs:
            continue
        for pa in panel:
            anm = pa.get("name") or "assay"
            for po in pa.get("oligos", []):
                ps = (po.get("seq") or "").upper()
                if not ps:
                    continue
                dg = T.hetero_dimer(cs, ps)
                if dg <= dimer_threshold:
                    end_dg = min(T.end_stability(cs, ps), T.end_stability(ps, cs))
                    cross.append(dict(oligo=cn, assay=anm, oligo_b=po.get("name") or "?",
                                      dg=round(dg, 2), end_dg=round(end_dg, 2), three_prime=end_dg <= -5.0))
    cross.sort(key=lambda x: x["dg"])
    melt = []
    if is_sybr and amplicon_tm is not None:
        for pa in panel:
            if pa.get("amplicon_tm") is None:
                continue
            psybr = pa.get("sybr")
            if psybr is None:
                psybr = not any((o.get("name") == "P") for o in pa.get("oligos", []))
            if not psybr:
                continue
            d = abs(float(amplicon_tm) - float(pa["amplicon_tm"]))
            if d < amp_tm_gap:
                melt.append(dict(assay=pa.get("name") or "assay", tm=round(float(pa["amplicon_tm"]), 1),
                                 delta=round(d, 1)))
        melt.sort(key=lambda x: x["delta"])
    return dict(n_panel=len(panel), cross=cross[:12], melt=melt,
                worst_dg=(cross[0]["dg"] if cross else None), three_prime=any(x["three_prime"] for x in cross))


def _parse_junction_template(text):
    """A single transcript whose exon boundaries are marked with '|'. Splits on '|', cleans each segment
    with the same normalizer used everywhere else, and records the boundary positions in the CLEANED
    sequence (clean_seq is character-wise, so piecewise cleaning concatenates to the same string and the
    positions line up with the oligo spans). Returns (cleaned_seq, [junction_positions], notes)."""
    segs = (text or "").split("|")
    cleaned = []
    for s in segs:
        c, _n, err = T.clean_seq(s)
        cleaned.append(c if (c and not err) else "")
    ref = "".join(cleaned)
    pos, run = [], 0
    for c in cleaned[:-1]:
        run += len(c)
        if 0 < run < len(ref):
            pos.append(run)
    pos = sorted(set(pos))
    notes = []
    if "|" in (text or "") and not pos:
        notes.append("exon marker '|' found but it isn't inside the sequence \u2014 ignoring it.")
    elif pos:
        notes.append("%d exon junction%s marked \u2014 candidates ranked for gDNA exclusion." %
                     (len(pos), "s" if len(pos) != 1 else ""))
    return ref, pos, notes


@app.post("/api/design")
def design(r: DesignReq):
    try:
        r.objective = _require_objective(r.objective)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    junctions = []
    if len(r.template or "") > T.MAX_TEMPLATE_LEN:
        return JSONResponse({"error": "template too long (%d nt; limit %d). Paste the transcript / "
                                      "amplicon region, not a whole chromosome." % (len(r.template or ""), T.MAX_TEMPLATE_LEN)},
                            status_code=200)
    if "|" in (r.template or ""):                            # single transcript with marked exon boundaries
        ref_single, junctions, tnotes = _parse_junction_template(r.template)
        targets = [ref_single] if (ref_single and len(ref_single) >= 60) else []
    else:
        targets, tnotes = _parse_fasta_targets(r.template)
    if not targets:
        return JSONResponse({"error": "no usable template \u2014 paste at least ~60 nt of the transcript / region "
                                      "(raw sequence, multi-record FASTA, or one sequence with '|' at exon boundaries)."},
                            status_code=200)
    offs, _onotes = _parse_fasta_targets(r.off_targets) if r.off_targets else ([], [])
    panel = [a for a in (r.panel or []) if a and a.get("oligos")]   # workbench assays to check multiplex fit against
    notes = list(tnotes or [])
    multi = len(targets) >= 2 and not junctions
    ref = _AD_design._reference(targets)                     # the sequence all spans / highlights map onto
    gc = 100.0 * sum(c in "GC" for c in ref) / max(1, len(ref))
    NCAND = 10
    # Same ranked-candidate engine as Target->assay; here the target set + off-target set are pasted, so multiple
    # targets unlock conservation + IUPAC-degenerate genus primers and an off-target set unlocks discrimination +
    # an offline in-silico-PCR verdict. One pasted sequence with no off-targets behaves exactly as before.
    if (r.profile or "").lower() == "auto":
        order, _gc = _AD_design._auto_order(ref)
        if _gc < 40.0 and "parasite_lna" not in order:
            order = order[:1] + ["parasite_lna"] + order[1:]
        res = None; used = None
        for pk in order:
            try:
                res = _AD_design.design_from_sequences(targets, P.PROFILES[pk], offs=offs or None,
                                                        n_candidates=NCAND, objective=r.objective,
                                                        junctions=junctions, panel=panel)
            except Exception:
                res = None
            if res and res.get("candidates"):
                used = pk; break
        if not (res and res.get("candidates")):
            return JSONResponse({"error": "no clean assay found under any Auto chemistry (tried: %s; reference GC %.0f%%). "
                                          "Try a longer / cleaner region." % (", ".join(order), gc)}, status_code=200)
        prof = P.PROFILES[used]
        notes.append("Auto-selected chemistry: %s (reference GC %.0f%%)." % (prof["name"], gc))
    else:
        pk = r.profile if r.profile in P.PROFILES else "idt_taqman"
        prof = P.PROFILES[pk]
        try:
            res = _AD_design.design_from_sequences(targets, prof, offs=offs or None,
                                                    n_candidates=NCAND, objective=r.objective,
                                                    junctions=junctions, panel=panel)
        except Exception as e:
            return _api_failure("design", e)
        if not (res and res.get("candidates")):
            if gc < 40:
                hint = " The reference is AT-rich (GC %.0f%%) \u2014 switch to Auto, 'low-Tm TaqMan', or 'AT-rich + LNA probe'." % gc
            elif gc > 62:
                hint = " The reference is GC-rich (GC %.0f%%) \u2014 switch to Auto or the 'GC-rich' profile." % gc
            else:
                hint = " Try Auto, a longer region, or a different chemistry."
            return JSONResponse({"error": "no clean assay found under %s.%s" % (prof["name"], hint)}, status_code=200)

    def _loc(sub):                                          # fallback for legacy candidates lacking coordinates
        i = ref.find((sub or "").upper())
        return [i, i + len(sub)] if (sub and i >= 0) else None

    def _span(assay, key, fallback):
        """Use the exact design-window coordinate when available.

        Falling back to sequence search is retained only for old/imported candidate
        records.  It is ambiguous when a target contains repeated motifs.
        """
        sp = assay.get(key)
        if isinstance(sp, (list, tuple)) and len(sp) == 2 and 0 <= sp[0] < sp[1] <= len(ref):
            return [int(sp[0]), int(sp[1])]
        return _loc(fallback)

    off_ctrl = None                                         # worst-case off-target amplicon as a bench negative control
    if offs:
        try:
            a0 = res["candidates"][0]["assay"]
            f0 = _span(a0, "f_xy", a0["forward"])
            r0 = _span(a0, "r_xy", T.revcomp(a0["reverse"]))
            if f0 and r0:
                g = D.build_offtarget_gblock(ref[f0[0]:r0[1]], offs)   # amplicon SEQUENCE, not its length
                if g:
                    off_ctrl = dict(seq=g["seq"], identity=g.get("amplicon_identity"), off_index=g.get("off_index"))
        except Exception:
            off_ctrl = None

    cands = []
    for sc in res["candidates"]:
        a = sc["assay"]; pi = a.get("probe_info")
        fspan = _span(a, "f_xy", a["forward"])
        rspan = _span(a, "r_xy", T.revcomp(a["reverse"]))  # reverse primer binds the antisense strand
        probe_rc = bool((a.get("probe_info") or {}).get("strand") == "-")
        pspan = _span(a, "probe_xy", (T.revcomp(a["probe"]) if probe_rc else a.get("probe"))) if a.get("probe") else None
        amp_span = [fspan[0], rspan[1]] if (fspan and rspan) else None

        cons = sc.get("conservation") or {}                 # per-oligo per-position match across targets
        cons_out = None
        if multi and cons:
            def _pp(k):
                d = cons.get(k)
                return [pp["pct_match"] for pp in d["per_pos"]] if (d and d.get("per_pos")) else None
            cons_out = dict(F=_pp("F"), R=_pp("R"), P=(_pp("P") if a.get("probe") else None),
                            mean=dict(F=(cons.get("F") or {}).get("mean_ident"),
                                      R=(cons.get("R") or {}).get("mean_ident"),
                                      P=(cons.get("P") or {}).get("mean_ident")))

        disc = sc.get("discrimination") or None             # per-oligo off-target mismatch (3'-block = good)
        disc_out = None
        if disc:
            def _dd(k):
                d = disc.get(k)
                if not d or not d.get("n"):
                    return None
                return dict(n=d.get("n"), max_ident=d.get("max_ident"), min_3p_mm=d.get("min_3prime_mismatch"))
            disc_out = dict(F=_dd("F"), R=_dd("R"), P=(_dd("P") if a.get("probe") else None))

        spec = None                                         # offline in-silico PCR over the pasted off-targets
        if offs:
            try:
                prods = _AD_design.epcr_offline(a["forward"], a["reverse"], offs, probe=a.get("probe"))
                spec = dict(n_products=len(prods),
                            products=[dict(off=p["subject"] + 1, size=p["size"], probe_binds=p["probe_binds"])
                                      for p in prods[:8]])
            except Exception:
                spec = None

        deg = None                                          # IUPAC-degenerate genus variants (>=2 targets)
        if a.get("forward_deg") or a.get("reverse_deg") or a.get("probe_deg"):
            deg = dict(forward=a.get("forward_deg"), reverse=a.get("reverse_deg"), probe=a.get("probe_deg"),
                       n=a.get("n_degenerate"), targets=a.get("deg_targets"))

        pf = None                                           # multiplex fit vs the current workbench panel
        if panel:
            try:
                cand_oligos = [("F", a["forward"]), ("R", a["reverse"])] + ([("P", a["probe"])] if a.get("probe") else [])
                pf = _panel_fit(cand_oligos, a.get("amplicon_tm"), not a.get("probe"), panel)
            except Exception:
                pf = None

        jx = None                                           # exon-junction / gDNA-exclusion relationship
        if junctions:
            def _crosses(sp):
                return [j for j in junctions if sp and sp[0] < j < sp[1]]
            pj = _crosses(pspan); fj = _crosses(fspan); rj = _crosses(rspan)
            aj = [j for j in junctions if amp_span and amp_span[0] < j < amp_span[1]]
            level = "strong" if (pj or fj or rj) else ("size" if aj else "none")
            jx = dict(level=level, probe=bool(pj), forward=bool(fj), reverse=bool(rj),
                      amplicon=bool(aj), junctions=sorted(set(pj + fj + rj + aj)), n_junctions=len(junctions))

        cands.append(dict(
            forward=a["forward"], reverse=a["reverse"], probe=a.get("probe"),
            amplicon=a["amplicon"], amplicon_tm=a.get("amplicon_tm"),
            pair_tm_gap=round(a["pair_tm_gap"], 1), f_tm=round(a["f_tm"], 1), r_tm=round(a["r_tm"], 1),
            probe_tm=round(pi["tm"], 1) if pi else None,
            probe_offset=round(pi["offset"], 1) if pi else None,
            probe_hairpin=round(pi["hairpin_dg"], 2) if pi else None,
            probe_dimer_f=round(pi["dimer_f"], 2) if pi else None,
            probe_dimer_r=round(pi["dimer_r"], 2) if pi else None,
            gblock=a["gblock"], score=sc.get("display_score", sc.get("score")),
            legacy_score=sc.get("legacy_score", sc.get("score")), quality_flags=a.get("quality_flags") or [],
            f_span=fspan, r_span=rspan, probe_span=pspan, probe_rc=probe_rc, amp_span=amp_span,
            score_breakdown=dict(raw=sc.get("score_raw"), quality=sc.get("quality_penalty"),
                                 amplicon=sc.get("amplicon_penalty")),
            conservation=cons_out, discrimination=disc_out, specificity=spec, degenerate=deg,
            panel_fit=(sc.get("evidence") or {}).get("panel_fit") or pf,
            junction=jx or (sc.get("evidence") or {}).get("junction"),
            rank=sc.get("rank"), pareto_front=(sc.get("evidence") or {}).get("pareto_front"),
            hard_valid=(sc.get("evidence") or {}).get("hard_valid"),
            hard_failures=(sc.get("evidence") or {}).get("hard_failures") or [],
            evidence=sc.get("evidence"), finalist_categories=sc.get("finalist_categories") or [],
            rank_explanation=sc.get("rank_explanation")))
    return dict(template=ref, profile=prof["name"], n=len(cands), candidates=cands,
                n_targets=len(targets), n_offs=len(offs), n_panel=len(panel),
                junctions=junctions, n_junctions=len(junctions), off_control=off_ctrl,
                objective_profile=res.get("objective_profile"), candidate_attrition=res.get("candidate_attrition"),
                ranker_manifest=res.get("ranker_manifest"), search_status=res.get("search_status"),
                ranking_statement=res.get("ranking_statement"),
                note=(" · ".join(notes) if notes else None), constraint_note=res.get("constraint_note"))


class AccessReq(BaseModel):
    amplicon: str
    forward: Optional[str] = None; reverse: Optional[str] = None; probe: Optional[str] = None
    anneal_c: Optional[float] = None


@app.post("/api/accessibility")
def api_accessibility(r: AccessReq):
    """Probe-site / primer-3' accessibility for ONE candidate: fold the amplicon (only) and report how
    base-paired each binding site is at the annealing temperature, where the probe must actually hybridize.
    On-demand (the Design call stays fast); folds the amplicon, never the template. ViennaRNA-guarded."""
    if not STR.available():
        return dict(available=False)
    amp = ("".join(c for c in (r.amplicon or "").upper() if c in "ACGTUN")).replace("U", "T")
    if len(amp) < 8:
        return dict(available=True, folded=False, note="amplicon too short to fold")
    anneal = r.anneal_c if r.anneal_c is not None else T.ANNEAL_C
    fe = STR.fold_ensemble(amp, anneal_c=anneal)
    if not fe:
        return dict(available=True, folded=False, note="amplicon outside the foldable range (8\u20131000 nt)")
    pa, pp = fe.get("paired_anneal"), fe.get("paired_prob")

    def locate(sub):
        if not sub:
            return None
        s = sub.upper().replace("U", "T")
        i = amp.find(s)
        if i < 0:
            i = amp.find(T.revcomp(s))
        return i if i >= 0 else None

    def site(start, end):
        start = max(0, start); end = min(len(amp), end)
        if end <= start:
            return None
        return dict(anneal_paired=(STR.site_paired_fraction(pa, start, end) if pa is not None else None),
                    ens_paired=STR.site_paired_prob(pp, start, end))

    out = dict(available=True, folded=True, anneal_c=anneal, dna_params=fe["dna_params"],
               mfe=fe["mfe"], mfe_anneal=fe.get("mfe_anneal"), n=fe["n"])
    if r.probe:
        pi = locate(r.probe)
        if pi is not None:
            out["probe"] = site(pi, pi + len(r.probe))
    if r.forward:
        L = len(r.forward)
        out["f3"] = site(L - 5, L)                          # forward 3' end sits at the amplicon's 5' start
    if r.reverse:
        L = len(r.reverse)
        out["r3"] = site(len(amp) - L, len(amp) - L + 5)    # reverse 3' end = left edge of revcomp(R) at the amplicon's 3' end

    def grade(frac):
        if frac is None:
            return None
        return "ok" if frac <= 0.0 else ("warn" if frac <= 0.4 else "risk")
    # verdict: probe accessibility if there's a probe, else the worse of the primer 3' ends
    if out.get("probe"):
        out["verdict"] = grade(out["probe"].get("anneal_paired"))
    else:
        fr = [s.get("anneal_paired") for s in (out.get("f3"), out.get("r3")) if s and s.get("anneal_paired") is not None]
        out["verdict"] = grade(max(fr)) if fr else None
    return out

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
        return _api_failure("NCBI fetch", e)

@app.post("/api/intron")
def intron(r: IntronReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return SP.intron_check(r.gene, r.organism, r.amp_start, r.amp_end, r.mrna_acc, r.forward, r.reverse)
    except Exception as e:
        return _api_failure("intron check", e)

@app.post("/api/blast")
def blast(r: BlastReq):
    _set_email(r.email, r.ncbi_key)
    denied = _local_blast_denied(r.mode, r.db_path)
    if denied:
        return denied
    try:
        return SP.blast_summary(r.seq, mode=r.mode, db=r.db, db_path=r.db_path,
                                organism=r.organism, top=r.top)
    except Exception as e:
        return _api_failure("BLAST", e)


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
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project storage")
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
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project storage")
    return _project_list()


@app.post("/api/project/list")
def project_list_post():
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project storage")
    return _project_list()

@app.post("/api/project/load")
def project_load(r: ProjectNameReq):
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project storage")
    fn = os.path.join(PROJECTS_DIR, _safe(r.name) + ".json")
    if not os.path.exists(fn):
        return JSONResponse({"error": "no project named " + r.name}, status_code=200)
    return json.load(open(fn))

@app.post("/api/project/delete")
def project_delete(r: ProjectNameReq):
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project storage")
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
    objective: str = "balanced"
    off_targets: Optional[str] = None
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
    if not r.items:
        return JSONResponse({"error": "batch contains no design items"}, status_code=422)
    if len(r.items) > 8:
        return JSONResponse({"error": "batch is limited to 8 templates per request"}, status_code=422)
    out = []
    for it in r.items:
        item_name = str(it.name or "unnamed")[:120]
        try:
            objective = _require_objective(it.objective)
        except ValueError as exc:
            out.append(dict(name=item_name, ok=False, error=str(exc)))
            continue
        if len(it.template or "") > T.MAX_TEMPLATE_LEN:
            out.append(dict(name=item_name, ok=False,
                            error="template exceeds %d nt" % T.MAX_TEMPLATE_LEN))
            continue
        targets, _notes = _parse_fasta_targets(it.template)
        if not targets:
            out.append(dict(name=item_name, ok=False,
                            error="no usable target sequence (need at least ~60 nt)"))
            continue
        offs, _off_notes = _parse_fasta_targets(it.off_targets) if it.off_targets else ([], [])
        ref = _AD_design._reference(targets)
        try:
            if str(it.profile or "").lower() == "auto":
                order, _gc = _AD_design._auto_order(ref)
                result = None; used_key = None
                for profile_key in order:
                    candidate_result = _AD_design.design_from_sequences(
                        targets, P.PROFILES[profile_key], offs=offs or None,
                        n_candidates=3, objective=objective, search_budget_s=12.0)
                    if candidate_result and candidate_result.get("candidates"):
                        result = candidate_result; used_key = profile_key; break
                prof = P.PROFILES[used_key] if used_key else None
            else:
                profile_key = it.profile if it.profile in P.PROFILES else "idt_taqman"
                prof = P.PROFILES[profile_key]
                result = _AD_design.design_from_sequences(
                    targets, prof, offs=offs or None,
                    n_candidates=3, objective=objective, search_budget_s=12.0)
        except Exception as exc:
            log.exception("batch design failed for %s", item_name)
            detail = "batch design failed" if HOSTED_MODE else str(exc)
            out.append(dict(name=item_name, ok=False, error=detail))
            continue
        if not (result and result.get("candidates")):
            out.append(dict(name=item_name, ok=False,
                            error=(result or {}).get("error") or "no hard-valid assay under that profile"))
            continue
        ranked = result["candidates"][0]
        a = ranked["assay"]
        pi = a.get("probe_info")
        evidence = ranked.get("evidence") or {}
        explanation = ranked.get("rank_explanation") or {}
        out.append(dict(
            name=item_name, ok=True, forward=a["forward"], reverse=a["reverse"],
            probe=a.get("probe"), amplicon=a["amplicon"],
            f_tm=round(a["f_tm"], 1), r_tm=round(a["r_tm"], 1),
            probe_tm=round(pi["tm"], 1) if pi else None, gblock=a.get("gblock"),
            profile=(prof or {}).get("name"), objective=objective,
            rank=ranked.get("rank"), display_score=ranked.get("display_score"),
            hard_valid=evidence.get("hard_valid"), hard_failures=evidence.get("hard_failures") or [],
            uncertainty=(explanation.get("uncertainty") or explanation.get("preference_strength")),
            strongest_feature=explanation.get("strongest_feature"),
            weakest_feature=explanation.get("weakest_feature"),
            ranking_statement=result.get("ranking_statement"),
            ranker_manifest=result.get("ranker_manifest"),
            candidate_attrition=result.get("candidate_attrition"),
            alternatives_evaluated=result.get("n_candidates_screened"),
        ))
    return dict(results=out, pipeline="authoritative_structured_ranker",
                policy=("Every successful batch winner passed the same retained-pool annotation and "
                        "structured ranking path used by interactive automatic design. Batch search uses "
                        "a declared 12-second per-template enumeration budget; the manifest records that limit."))

@app.post("/api/panel/save")
def panel_save(r: PanelSaveReq):
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server panel storage")
    fn = os.path.join(PANELS_DIR, _safe(r.name) + ".json")
    json.dump({"name": r.name, "oligos": r.oligos}, open(fn, "w"), indent=1)
    return dict(saved=r.name)

@app.get("/api/panel/list")
def panel_list():
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server panel storage")
    return dict(panels=sorted(f[:-5] for f in os.listdir(PANELS_DIR) if f.endswith(".json")))


@app.post("/api/factory_reset")
def factory_reset():
    """Delete every saved panel and project on the server -- the server half of a factory reset.
    Only *.json inside the two managed directories are touched (names come from os.listdir, so no
    path traversal). The browser clears its own localStorage/session separately."""
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server project/panel storage")
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
    if not ALLOW_SERVER_STORAGE:
        return _shared_feature_disabled("server panel storage")
    fn = os.path.join(PANELS_DIR, _safe(r.name) + ".json")
    if not os.path.exists(fn):
        return JSONResponse({"error": "panel not found"}, status_code=200)
    return json.load(open(fn))

@app.post("/api/order_csv")
def order_csv(r: OrderReq):
    try:
        return dict(oligo_csv=O.oligo_csv([e.dict() for e in r.oligos]) if r.oligos else "",
                    gblock_fasta=O.gblock_fasta([g.dict() for g in r.gblocks]) if r.gblocks else "")
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=200)

@app.post("/api/pair_specificity")
def pair_specificity(r: PairSpecReq):
    _set_email(r.email, r.ncbi_key)
    denied = _local_blast_denied(r.mode, r.db_path)
    if denied:
        return denied
    def one(seq):
        if r.mode == "local":
            return SP.blast_local(seq, r.db_path or "")
        return SP.blast_remote(seq, organism=r.organism, db=r.db)
    try:
        return dict(forward=one(r.forward.upper().strip()), reverse=one(r.reverse.upper().strip()))
    except Exception as e:
        return _api_failure("specificity check", e)


# ============ stage-6: conservation, in-silico PCR, LNA Tm ============
from oligoforge import conservation as CONS, thermo as T

class ConsReq(BaseModel):
    oligos: Dict[str, str]; targets: List[str]
    off_targets: List[str] = []; min_ident: float = 0.6
class EpcrReq(BaseModel):
    forward: str; reverse: str; mode: str = "remote"; db: str = "nt"
    db_path: Optional[str] = None; organism: Optional[str] = None
    min_product: int = 40; max_product: int = 3000; email: Optional[str] = None; ncbi_key: Optional[str] = None
    fasta: Optional[str] = None; max_mm: int = 2
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
    denied = _local_blast_denied(r.mode, r.db_path)
    if denied:
        return denied
    try:
        return SP.in_silico_pcr(r.forward, r.reverse, r.mode, r.db, r.db_path,
                                r.organism, r.min_product, r.max_product,
                                fasta=r.fasta, max_mm=r.max_mm)
    except Exception as e:
        return _api_failure("in-silico PCR", e)

class AssaySpecReq(BaseModel):
    forward: str; reverse: str; probe: Optional[str] = None
    mode: str = "remote"; db: str = "nt"; db_path: Optional[str] = None
    organism: Optional[str] = None
    min_product: int = 40; max_product: int = 3000
    email: Optional[str] = None; ncbi_key: Optional[str] = None
    fasta: Optional[str] = None; max_mm: int = 2
@app.post("/api/assay_specificity")
def api_assay_specificity(r: AssaySpecReq):
    _set_email(r.email, r.ncbi_key)
    denied = _local_blast_denied(r.mode, r.db_path)
    if denied:
        return denied
    try:
        return SP.assay_specificity(r.forward, r.reverse, r.probe, r.mode, r.db, r.db_path,
                                    r.organism, r.min_product, r.max_product,
                                    fasta=r.fasta, max_mm=r.max_mm)
    except Exception as e:
        return _api_failure("assay specificity", e)

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
        return _api_failure("gene lookup", e)

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
        return _api_failure("report", e)

from oligoforge import rdml as RDML
class RdmlReq(BaseModel):
    panel: List[dict]; meta: Optional[dict] = None
@app.post("/api/rdml")
def api_rdml(r: RdmlReq):
    """Export the panel as RDML 1.3 (machine-readable assay definitions for qPCR software)."""
    try:
        return RDML.build(r.panel, r.meta)
    except Exception as e:
        return _api_failure("RDML export", e)

class MultiplexReq(BaseModel):
    assays: List[dict]
    dimer_threshold: float = -6.0
    amp_tm_gap: float = 2.0
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
        return _api_failure("marker scan", e)

@app.post("/api/suggest_genes")
def api_suggest_genes(r: MarkerReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return RM.suggest_dynamic(r.organism, r.exclude, r.intent)
    except Exception as e:
        return _api_failure("gene suggestion", e)

@app.post("/api/multiplex")
def api_multiplex(r: MultiplexReq):
    try:
        return MX.check(r.assays, r.dimer_threshold, r.amp_tm_gap)
    except Exception as e:
        return _api_failure("multiplex analysis", e)

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
        return _api_failure("reference-gene stability analysis", e)

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
        return _api_failure("sequence fetch", e)

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
        return _api_failure("expression analysis", e)

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

from oligoforge import miqe as MIQE
class ValidateReq(BaseModel):
    assay: dict
    observed_amplicon_tm: Optional[float] = None
    observed_peaks: Optional[int] = None
@app.post("/api/validate")
def api_validate(r: ValidateReq):
    try:
        return MIQE.validate_assay(r.assay, observed_amplicon_tm=r.observed_amplicon_tm,
                                   observed_peaks=r.observed_peaks)
    except Exception as e:
        return _api_failure("assay-readiness validation", e)

from oligoforge import autodesign as AD
class AutoDesignReq(BaseModel):
    target_query: str; profile: str = "auto"
    off_query: Optional[str] = None; n_fetch: int = 20
    min_ident: float = 0.6
    run_blast: bool = False; blast_mode: str = "remote"
    blast_db: str = "nt"; blast_db_path: Optional[str] = None; organism: Optional[str] = None
    email: Optional[str] = None; ncbi_key: Optional[str] = None; prefer_junction: bool = False; nested: bool = False
    objective: str = "balanced"
@app.post("/api/autodesign")
def api_autodesign(r: AutoDesignReq):
    try:
        r.objective = _require_objective(r.objective)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    _set_email(r.email, r.ncbi_key)
    denied = _local_blast_denied(r.blast_mode, r.blast_db_path) if r.run_blast else None
    if denied:
        return denied
    try:
        return AD.design_from_query(r.target_query, r.profile, r.off_query, min(r.n_fetch, 30),
                                    r.min_ident, r.run_blast, r.blast_mode, r.blast_db,
                                    r.blast_db_path, r.organism, prefer_junction=r.prefer_junction,
                                    nested=r.nested, objective=r.objective)
    except Exception as e:
        return _api_failure("automatic design", e)


# ============ ranking-truth manual design and assay rescue ============
class ManualDesignReq(BaseModel):
    forward: str; reverse: str; template: str
    probe: Optional[str] = None; profile: str = "idt_taqman"; objective: str = "balanced"
    targets: Optional[List[str]] = None; off_targets: Optional[List[str]] = None; max_mm: int = 2

class RedesignReq(ManualDesignReq):
    locks: Optional[Dict[str, bool]] = None; max_results: int = 8
    max_shift: Optional[int] = None; excluded_regions: Optional[List[List[int]]] = None
    amp_min: Optional[int] = None; amp_max: Optional[int] = None
    required_region: Optional[List[int]] = None

class RescueReq(ManualDesignReq):
    observed: Optional[dict] = None; max_results: int = 4

class ManualEditCompareReq(BaseModel):
    baseline_forward: str; baseline_reverse: str
    edited_forward: str; edited_reverse: str
    template: str
    baseline_probe: Optional[str] = None; edited_probe: Optional[str] = None
    profile: str = "idt_taqman"; objective: str = "balanced"
    targets: Optional[List[str]] = None; off_targets: Optional[List[str]] = None
    max_mm: int = 2

class FeedbackReq(BaseModel):
    records: List[dict]

class FeedbackImportReq(BaseModel):
    payload: str
    format_hint: str = "auto"

class FeedbackSplitReq(FeedbackReq):
    train_fraction: float = 0.70
    validation_fraction: float = 0.15


class RunCompareReq(BaseModel):
    left: dict
    right: dict
    top_k: int = 10


def _profile_or_default(key):
    return _require_profile(key)


def _manual_inputs(r):
    objective = _require_objective(r.objective)
    if len(r.template or "") > T.MAX_TEMPLATE_LEN:
        raise ValueError("template exceeds %d nt" % T.MAX_TEMPLATE_LEN)
    targets = _bounded_sequence_list(r.targets, "target set") if r.targets is not None else None
    offs = _bounded_sequence_list(r.off_targets, "off-target set") if r.off_targets is not None else None
    return _profile_or_default(r.profile), objective, targets, offs

@app.post("/api/manual-design/analyze")
def api_manual_design(r: ManualDesignReq):
    try:
        profile, objective, targets, offs = _manual_inputs(r)
        return MDS.analyze_assay(r.forward, r.reverse, r.template, profile,
                                 r.probe, targets=targets, offs=offs,
                                 objective=objective, max_mm=max(0, min(int(r.max_mm), 4)))
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("manual assay analysis", e)

@app.post("/api/manual-design/redesign")
def api_manual_redesign(r: RedesignReq):
    try:
        profile, objective, targets, offs = _manual_inputs(r)
        excluded = [(int(x[0]), int(x[1])) for x in (r.excluded_regions or []) if len(x) == 2]
        required = r.required_region if (r.required_region and len(r.required_region) == 2) else None
        return MDS.constrained_redesign(r.forward, r.reverse, r.template, profile,
                                        r.probe, locks=r.locks, objective=objective,
                                        max_results=max(1, min(int(r.max_results), 20)),
                                        max_shift=r.max_shift, excluded_regions=excluded,
                                        max_mm=max(0, min(int(r.max_mm), 4)),
                                        amp_min=r.amp_min, amp_max=r.amp_max,
                                        required_region=required, targets=targets, offs=offs)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("constrained redesign", e)

@app.post("/api/assay-rescue")
def api_assay_rescue(r: RescueReq):
    try:
        profile, objective, targets, offs = _manual_inputs(r)
        return ARES.rescue(r.forward, r.reverse, r.template, profile,
                           r.probe, objective=objective, observed=r.observed,
                           targets=targets, offs=offs,
                           max_results=max(1, min(int(r.max_results), 10)))
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("assay rescue", e)


@app.post("/api/manual-design/compare-edit")
def api_manual_compare_edit(r: ManualEditCompareReq):
    try:
        objective = _require_objective(r.objective)
        profile = _profile_or_default(r.profile)
        if len(r.template or "") > T.MAX_TEMPLATE_LEN:
            raise ValueError("template exceeds %d nt" % T.MAX_TEMPLATE_LEN)
        targets = _bounded_sequence_list(r.targets, "target set") if r.targets is not None else None
        offs = _bounded_sequence_list(r.off_targets, "off-target set") if r.off_targets is not None else None
        return MDS.compare_edits(
            r.baseline_forward, r.baseline_reverse,
            r.edited_forward, r.edited_reverse, r.template, profile,
            baseline_probe=r.baseline_probe, edited_probe=r.edited_probe,
            targets=targets, offs=offs, objective=objective,
            max_mm=max(0, min(int(r.max_mm), 4)))
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("manual assay edit comparison", e)

@app.post("/api/experimental-feedback/status")
def api_feedback_status(r: FeedbackReq):
    try:
        return EFB.calibration_status(r.records)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("experimental feedback", e)


@app.post("/api/experimental-feedback/import")
def api_feedback_import(r: FeedbackImportReq):
    try:
        records = EFB.parse_records(r.payload, r.format_hint)
        return EFB.dataset_status(records)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("experimental feedback import", e)

@app.post("/api/experimental-feedback/split")
def api_feedback_split(r: FeedbackSplitReq):
    try:
        return EFB.target_group_split(r.records, r.train_fraction, r.validation_fraction)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("experimental feedback split", e)


@app.post("/api/experimental-feedback/summary")
def api_feedback_summary(r: FeedbackReq):
    try:
        return EFB.evidence_summary(r.records)
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("experimental feedback summary", e)


@app.post("/api/design-runs/compare")
def api_run_compare(r: RunCompareReq):
    try:
        return RCOMP.compare_runs(r.left, r.right, top_k=max(1, min(int(r.top_k), 100)))
    except ValueError as e:
        return JSONResponse({"error": str(e)}, status_code=422)
    except Exception as e:
        return _api_failure("design run comparison", e)

# ============ orthogonal panel analysis (exact proofs where available; diagnostics otherwise) ============
from oligoforge import orthopanel as OP

def _parse_ortho_candidates(text):
    """Parse a pasted block into candidate dicts. Supports FASTA (>name / sequence lines), one
    'name<comma|tab|space>sequence' per line, or a bare sequence per line."""
    if not text:
        return []
    lines = text.splitlines()
    if any(l.lstrip().startswith(">") for l in lines):
        cands, name, seq = [], None, []
        for raw in lines:
            l = raw.strip()
            if not l:
                continue
            if l.startswith(">"):
                if seq:
                    cands.append({"name": name, "seq": "".join(seq)})
                name, seq = (l[1:].strip() or None), []
            else:
                seq.append(l)
        if seq:
            cands.append({"name": name, "seq": "".join(seq)})
        return cands
    cands = []
    for raw in lines:
        l = raw.strip()
        if not l:
            continue
        parts = None
        for sep in (",", "\t"):
            if sep in l:
                parts = [p.strip() for p in l.split(sep, 1)]
                break
        if parts is None:
            toks = l.split()
            if len(toks) == 2 and len(toks[1]) >= len(toks[0]):
                parts = toks
        if parts and len(parts) == 2 and parts[1]:
            cands.append({"name": parts[0] or None, "seq": parts[1]})
        else:
            cands.append({"seq": l})
    return cands

class OrthoPanelReq(BaseModel):
    candidates: Optional[List[Dict]] = None; text: Optional[str] = None
    cross_dg: float = -6.0; self_dg: float = -9.0; k: int = 1
    size_limit: int = 600; use_theta: bool = True
    mv_conc: Optional[float] = None; dv_conc: Optional[float] = None
    dntp_conc: Optional[float] = None; dna_conc: Optional[float] = None; anneal_c: Optional[float] = None

@app.post("/api/orthogonal-panel")
def api_orthogonal_panel(r: OrthoPanelReq):
    try:
        cands = list(r.candidates or []) + _parse_ortho_candidates(r.text)
        if not cands:
            return JSONResponse({"error": "no candidate oligos provided"}, status_code=200)
        if len(cands) > 1000:
            return JSONResponse({"error": "too many candidates (%d); this tool caps at 1000"
                                          % len(cands)}, status_code=200)
        k = max(1, min(int(r.k or 1), 12))
        vals = dict(T.COND)
        limits = {"mv_conc": (0.0, 2000.0), "dv_conc": (0.0, 200.0),
                  "dntp_conc": (0.0, 100.0), "dna_conc": (1e-6, 1e6)}
        for name, raw in (("mv_conc", r.mv_conc), ("dv_conc", r.dv_conc),
                          ("dntp_conc", r.dntp_conc), ("dna_conc", r.dna_conc)):
            if raw is None:
                continue
            value = float(raw)
            lo, hi = limits[name]
            if not (lo <= value <= hi):
                return JSONResponse({"error": "%s out of range (%g..%g)" % (name, lo, hi)}, status_code=200)
            vals[name] = value
        anneal = T.ANNEAL_C if r.anneal_c is None else float(r.anneal_c)
        if not (30.0 <= anneal <= 85.0):
            return JSONResponse({"error": "anneal_c out of range (30..85)"}, status_code=200)
        if vals["mv_conc"] + vals["dv_conc"] <= 0:
            return JSONResponse({"error": "need some salt: monovalent + divalent must be > 0 mM"}, status_code=200)
        snap = (vals["mv_conc"], vals["dv_conc"], vals["dntp_conc"], vals["dna_conc"], anneal)
        out = OP.certify_panel(cands, cross_dg=float(r.cross_dg), self_dg=float(r.self_dg),
                               k=k, size_limit=int(r.size_limit or 600), use_theta=bool(r.use_theta),
                               thermo_snapshot=snap)
        out["conditions"] = dict(vals, anneal_c=anneal)
        return out
    except Exception as e:
        return _api_failure("orthogonal-panel analysis", e)


# ============ isolate panel validation (inclusivity / exclusivity in-silico PCR) ============
class IsolateGenomesReq(BaseModel):
    query: str; retmax: int = 60
    email: Optional[str] = None; ncbi_key: Optional[str] = None

@app.post("/api/isolate_genomes")
def api_isolate_genomes(r: IsolateGenomesReq):
    _set_email(r.email, r.ncbi_key)
    try:
        return {"genomes": ncbi.search_genomes(r.query, min(max(int(r.retmax or 60), 5), 200))}
    except Exception as e:
        return _api_failure("genome search", e)

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
            log.exception("isolate check failed for accession %s", acc)
            detail = "accession analysis failed" if HOSTED_MODE else str(e)
            out.append(dict(acc=acc, title="", slen=0, role=r.role, amplifies=None, error=detail))
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
        return _api_failure("design", e)

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
