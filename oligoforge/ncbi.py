"""NCBI retrieval via Entrez (biopython): fetch by accession, search a gene in
an organism, pull all mRNA isoforms, and find the region common to every isoform
(so an assay reads total transcript, not one variant)."""
import time, difflib, socket, os, re
from Bio import Entrez, SeqIO

# Bound every NCBI HTTP call so a genuinely stalled connection eventually fails
# instead of hanging forever — but generously, because a remote blastn against
# nt is queued server-side at NCBI and a single poll can legitimately take
# several minutes overall, but no single HTTP operation should monopolize a worker.
# The default is 30 seconds per network operation; retry logic provides bounded recovery.
# Override OLIGOFORGE_NCBI_TIMEOUT for a known slow private deployment.
try:
    NCBI_TIMEOUT = float(os.environ.get("OLIGOFORGE_NCBI_TIMEOUT", "30"))
except ValueError:
    NCBI_TIMEOUT = 30.0
NCBI_TIMEOUT = min(max(NCBI_TIMEOUT, 3.0), 300.0)
socket.setdefaulttimeout(NCBI_TIMEOUT)
from io import StringIO

Entrez.email = "set-me@example.com"   # set to your address; NCBI requires it
Entrez.api_key = None                  # optional: set for 10 req/s instead of 3


# ---- transient-failure retry + on-disk FASTA cache (reliability layer) ----
# NCBI is the single remote dependency. A transient 502 / timeout used to fail a whole design, and
# an identical re-run always re-hit the network. _net() retries the network-OPENING call with
# backoff on transient errors only (a 400/404 is a logical error and is NOT retried); efetch FASTA
# text (immutable per accession) is cached on disk so repeats are instant and survive a brief NCBI
# blip. Env knobs: OLIGOFORGE_NCBI_CACHE (1/0, default on), OLIGOFORGE_NCBI_CACHE_TTL (s, default 7d),
# OLIGOFORGE_NCBI_RETRIES (default 3). Cache dir = OLIGOFORGE_DATA_PATH/ncbi_cache or a temp dir.
import hashlib, json as _json, tempfile, urllib.error, http.client


def _int_env(name, default):
    try:
        return int(float(os.environ.get(name, default)))
    except (TypeError, ValueError):
        return int(default)


_CACHE_ON = os.environ.get("OLIGOFORGE_NCBI_CACHE", "1") not in ("0", "false", "False", "")
_CACHE_TTL = _int_env("OLIGOFORGE_NCBI_CACHE_TTL", 7 * 24 * 3600)
_RETRIES = min(5, max(1, _int_env("OLIGOFORGE_NCBI_RETRIES", 3)))
_CACHE_DIR = os.path.join(os.environ.get("OLIGOFORGE_DATA_PATH") or tempfile.gettempdir(),
                          "oligoforge_ncbi_cache")
_TRANSIENT = (socket.timeout, urllib.error.URLError, http.client.HTTPException,
              ConnectionError, TimeoutError)
_RETRY_HTTP = {400, 429, 500, 502, 503, 504}   # NCBI intermittently 400s a well-formed query under load (clears on retry)
# NCBI E-utilities periodically return an <ERROR> element instead of results when their
# backend is briefly overloaded (more common on broad queries over huge taxa like Plasmodium).
# Biopython raises this as a RuntimeError at PARSE time -- AFTER the HTTP call returns -- so it
# escapes the per-call retry in _net. These substrings mark such transient server-side errors,
# which clear on a retry; matched case-insensitively against str(e).
_TRANSIENT_MSGS = ("search backend failed", "error reading from backend", "backend failed",
                   "unable to obtain", "temporarily unavailable", "service unavailable",
                   "server is temporarily unable", "timed out", "try again")
# esearch+read is retried more aggressively than _net's default: each attempt is cheap and the
# backend error is intermittent (~half the time on the worst queries), so a few extra tries take
# the residual failure rate to near zero.
_SEARCH_RETRIES = max(_RETRIES, 3)


def _is_transient(e):
    if isinstance(e, urllib.error.HTTPError):
        return getattr(e, "code", None) in _RETRY_HTTP
    if isinstance(e, _TRANSIENT):
        return True
    msg = str(e).lower()
    return any(s in msg for s in _TRANSIENT_MSGS)


def _net(fn, *a, **k):
    """Call an Entrez network function, retrying transient failures with backoff; re-raise others."""
    last = None
    for attempt in range(_RETRIES):
        try:
            return fn(*a, **k)
        except Exception as e:
            last = e
            if not _is_transient(e) or attempt == _RETRIES - 1:
                raise
            time.sleep(min(8.0, 0.5 * (3 ** attempt)))   # 0.5s, 1.5s, 4.5s ...
    raise last


def _retry_read(call, *a, **k):
    """Run an Entrez query AND Entrez.read() inside one retry loop, so a transient backend
    error surfaced at parse time (e.g. 'Search Backend failed') is retried -- _net only wraps
    the HTTP call and returns before the read happens. Backoff: 0.4, 0.8, 1.6, 3.2, 6s ...
    Re-raises a non-transient error immediately; on exhausting retries against a transient
    server error, raises a clean, user-facing message instead of the raw NCBI string."""
    last = None
    for attempt in range(_SEARCH_RETRIES):
        try:
            h = call(*a, **k)
            r = Entrez.read(h)
            try:
                h.close()
            except Exception:
                pass
            return r
        except Exception as e:
            last = e
            if not _is_transient(e) or attempt == _SEARCH_RETRIES - 1:
                break
            time.sleep(min(6.0, 0.4 * (2 ** attempt)))
    msg = str(last).lower()
    if any(s in msg for s in _TRANSIENT_MSGS):
        raise RuntimeError("NCBI search is temporarily unavailable (retried %d times). "
                           "Please try again in a moment." % _SEARCH_RETRIES) from last
    raise last


def _esearch_read(db, term, retmax, **extra):
    """esearch + read with parse-time retry. Returns the parsed esearch result (has IdList, Count)."""
    return _retry_read(Entrez.esearch, db=db, term=term, retmax=retmax, **extra)


def _esummary_read(db, ids, **extra):
    """esummary + read with parse-time retry. `ids` may be a list or comma-joined string."""
    if isinstance(ids, (list, tuple)):
        ids = ",".join(str(x) for x in ids)
    return _retry_read(Entrez.esummary, db=db, id=ids, **extra)


def _cache_key(**k):
    norm = dict(k)
    if isinstance(norm.get("id"), (list, tuple)):
        norm["id"] = ",".join(str(x) for x in norm["id"])
    return hashlib.sha256(repr(sorted(norm.items())).encode("utf-8")).hexdigest()


def _cache_get(key):
    if not _CACHE_ON:
        return None
    try:
        with open(os.path.join(_CACHE_DIR, key + ".json"), encoding="utf-8") as fh:
            rec = _json.load(fh)
        if (time.time() - rec.get("ts", 0)) <= _CACHE_TTL:
            return rec.get("text")
    except Exception:
        return None
    return None


def _cache_put(key, text):
    if not _CACHE_ON:
        return
    try:
        os.makedirs(_CACHE_DIR, exist_ok=True)
        tmp = os.path.join(_CACHE_DIR, key + ".tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            _json.dump({"ts": time.time(), "text": text}, fh)
        os.replace(tmp, os.path.join(_CACHE_DIR, key + ".json"))
    except Exception:
        pass


def _efetch_text(**k):
    """efetch -> decoded text, with retry + on-disk cache (FASTA/GB is immutable by accession)."""
    key = _cache_key(op="efetch", **k)
    hit = _cache_get(key)
    if hit is not None:
        return hit
    h = _net(Entrez.efetch, **k)
    try:
        raw = h.read()
    finally:
        h.close()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", "replace")
    _cache_put(key, raw)
    return raw


def _records(ids, rettype="fasta"):
    raw = _efetch_text(db="nucleotide", id=list(ids), rettype=rettype, retmode="text")
    if rettype == "gb":
        return list(SeqIO.parse(StringIO(raw), "genbank"))
    cut = raw.find(">")                     # drop any NCBI notice/comment before the first FASTA record
    return list(SeqIO.parse(StringIO(raw[cut:] if cut >= 0 else ""), "fasta"))


def fetch_accessions(accs, rettype="fasta"):
    return _records(accs, rettype)


def search_mrna(gene, organism, retmax=20):
    """Return mRNA accessions for a gene in an organism (RefSeq predicted/curated)."""
    for q in (f'{organism}[Organism] AND {gene}[Gene] AND biomol_mrna[PROP]',
              f'{organism}[Organism] AND {gene}[All Fields] AND biomol_mrna[PROP]',
              f'{organism}[Organism] AND {gene} AND biomol_mrna[PROP]'):
        r = _esearch_read("nucleotide", q, retmax)
        if r["IdList"]:
            return r["IdList"], q
        time.sleep(0.11 if Entrez.api_key else 0.34)
    return [], None


def fetch_isoforms(gene, organism):
    ids, q = search_mrna(gene, organism)
    if not ids:
        return [], None
    recs = _records(ids[:8], "fasta")
    recs = [r for r in recs if len(r.seq) <= 50000]   # belt+braces: drop any genomic contig
    return recs, q


def common_region(seqs):
    """Longest contiguous block present in EVERY sequence (isoform-common)."""
    if not seqs:
        return ""
    seqs = [str(s).upper() for s in seqs]
    block = seqs[0]
    for other in seqs[1:]:
        sm = difflib.SequenceMatcher(None, block, other, autojunk=False)
        m = sm.find_longest_match(0, len(block), 0, len(other))
        block = block[m.a:m.a + m.size]
        if not block:
            break
    return block


def gene_id(gene, organism):
    """Entrez Gene ID for a gene in an organism. Tries progressively looser field
    tags so a symbol that isn't indexed as an exact [Gene Name] still resolves."""
    for term in (f'{gene}[Gene Name] AND {organism}[Organism]',
                 f'{gene}[Gene] AND {organism}[Organism]',
                 f'{gene}[All Fields] AND {organism}[Organism]'):
        r = _esearch_read("gene", term, 20)
        if r["IdList"]:
            return r["IdList"][0]
        time.sleep(0.11 if Entrez.api_key else 0.34)
    return None


_GENE_DECOY = ("receptor", "antagonist", "binding", "associated", "induced", "inducible",
               "like", "regulator", "activator", "inhibitor", "antisense", "converting",
               "pseudogene", "interacting", "downstream")


def _gene_match_score(doc, gene_q):
    """Higher = closer to the gene the user named. An exact description/symbol match wins;
    identity-changing modifiers the user did NOT type (receptor, antagonist, -like ...) are
    penalised, so 'interferon gamma' resolves to IFNG, not IFNGR1/IFNGR2."""
    desc = " ".join(str(doc.get("Description") or "").lower().split())
    name = str(doc.get("Name") or "").lower().strip()
    q = " ".join((gene_q or "").lower().split())
    qn = q.replace("-", "").replace(" ", "")
    s = 0.0
    if name == q or desc == q:
        s += 100
    if name.replace("-", "") == qn or desc.replace("-", "").replace(" ", "") == qn:
        s += 80
    qwords = q.split()
    if qwords and all(w in desc.split() for w in qwords):
        s += 25
    for d in _GENE_DECOY:
        if d in desc and d not in q:
            s -= 40
    s -= max(0, len(desc.split()) - len(qwords)) * 3
    return s


def resolve_gene(gene, organism=None):
    """Resolve a free-text gene name to NCBI's OFFICIAL gene symbol via the Gene database,
    so a design targets the right locus (IFNG, not interferon-gamma RECEPTOR 1 / IFNGR1).
    Among all candidates it picks the best description match (see _gene_match_score). Returns
    dict(found, symbol, gene_id, organism, description, in_requested_organism, clean). When the
    organism's only matches are decoys (e.g. just the receptor), it looks the gene up without the
    organism so the caller can say where the real gene exists instead of designing on a paralog.
    """
    gene = (gene or "").strip()
    if not gene:
        return dict(found=False)
    org = (organism or "").strip()
    tiers = []
    if org:
        tiers += [(f'"{gene}"[Gene Name] AND "{org}"[Organism]', True),
                  (f'{gene}[Gene Name] AND {org}[Organism]', True),
                  (f'{gene}[Gene Full Name] AND {org}[Organism]', True),
                  (f'{gene}[All Fields] AND {org}[Organism] AND alive[prop]', True)]
    tiers += [(f'"{gene}"[Gene Name]', False), (f'{gene}[Gene Name]', False),
              (f'{gene}[Gene Full Name]', False),
              (f'{gene}[All Fields] AND alive[prop]', False)]

    def best_of(ids):
        try:
            s = _esummary_read("gene", ids[:50])
            docs = list(s["DocumentSummarySet"]["DocumentSummary"])
        except Exception:
            return None, 0.0
        if not docs:
            return None, 0.0
        docs.sort(key=lambda d: _gene_match_score(d, gene), reverse=True)
        return docs[0], _gene_match_score(docs[0], gene)

    for term, in_org in tiers:
        try:
            ids = _esearch_read("gene", term, 50)["IdList"]
        except Exception:
            ids = []
        if ids:
            doc, score = best_of(ids)
            if doc is None:
                continue
            clean = score >= 25                      # all query words present, no unrequested modifier
            if in_org and not clean:
                continue                             # only decoys here -> look the gene up elsewhere
            try:
                uid = doc.attributes.get("uid", ids[0])
            except Exception:
                uid = ids[0]
            return dict(found=True, symbol=str(doc.get("Name") or gene), gene_id=uid,
                        organism=str((doc.get("Organism") or {}).get("ScientificName") or org or ""),
                        description=str(doc.get("Description") or ""),
                        in_requested_organism=in_org, clean=clean)
        time.sleep(0.11 if Entrez.api_key else 0.34)
    return dict(found=False)


def gene_lookup(gene, organism=None, max_candidates=8):
    """User-facing gene resolution for the specific-gene finder. Returns the best-match official
    gene PLUS every other gene whose name matched the text (so IFNG vs IFNGR1/IFNGR2 is laid out
    for the user), and how many nucleotide mRNA records exist for the best match. Honest about
    whether the gene is annotated in the requested organism or only found elsewhere."""
    gene = (gene or "").strip(); org = (organism or "").strip()
    if not gene:
        return dict(found=False, error="enter a gene name or symbol")
    ids, scoped = [], False
    if org:
        for term in (f'{gene}[Gene Name] AND {org}[Organism]',
                     f'{gene}[All Fields] AND {org}[Organism] AND alive[prop]'):
            try:
                ids = _esearch_read("gene", term, max_candidates)["IdList"]
            except Exception:
                ids = []
            if ids:
                scoped = True; break
            time.sleep(0.11 if Entrez.api_key else 0.34)
    if not ids:
        for term in (f'{gene}[Gene Name]', f'{gene}[All Fields] AND alive[prop]'):
            try:
                ids = _esearch_read("gene", term, max_candidates)["IdList"]
            except Exception:
                ids = []
            if ids:
                break
            time.sleep(0.11 if Entrez.api_key else 0.34)
    if not ids:
        return dict(found=False, gene=gene, organism=org,
                    error="NCBI Gene has no record matching \u201c%s\u201d%s" % (gene, (" in " + org) if org else ""))
    try:
        docs = list(_esummary_read("gene", ids[:max_candidates])["DocumentSummarySet"]["DocumentSummary"])
    except Exception:
        docs = []
    cand = []
    for d in docs:
        try:
            uid = d.attributes.get("uid")
        except Exception:
            uid = None
        cand.append(dict(symbol=str(d.get("Name") or ""), description=str(d.get("Description") or ""),
                         organism=str((d.get("Organism") or {}).get("ScientificName") or ""),
                         aliases=str(d.get("OtherAliases") or ""), gene_id=uid,
                         score=_gene_match_score(d, gene)))
    cand.sort(key=lambda c: c["score"], reverse=True)
    best = cand[0] if cand else None
    # Ortholog fallback: the typed name resolved only in another organism (e.g. chicken "INFG"),
    # but the requested organism may carry the same gene under its own symbol — find it by the
    # resolved description so "INFG" + Aphelocoma still lands on Aphelocoma IFNG.
    if org and not scoped and best and best.get("description"):
        desc = best["description"]
        for term in (f'{desc}[All Fields] AND {org}[Organism] AND alive[prop]',
                     f'"{desc}"[Gene Name] AND {org}[Organism]'):
            try:
                oids = _esearch_read("gene", term, max_candidates)["IdList"]
            except Exception:
                oids = []
            if oids:
                try:
                    odocs = list(_esummary_read("gene", oids[:max_candidates])["DocumentSummarySet"]["DocumentSummary"])
                except Exception:
                    odocs = []
                ocand = []
                for d in odocs:
                    try:
                        uid = d.attributes.get("uid")
                    except Exception:
                        uid = None
                    ocand.append(dict(symbol=str(d.get("Name") or ""), description=str(d.get("Description") or ""),
                                      organism=str((d.get("Organism") or {}).get("ScientificName") or ""),
                                      aliases=str(d.get("OtherAliases") or ""), gene_id=uid,
                                      score=_gene_match_score(d, desc)))
                ocand.sort(key=lambda c: c["score"], reverse=True)
                if ocand and ocand[0]["score"] >= 25:
                    best = ocand[0]; scoped = True
                    cand = ocand + [c for c in cand if c["symbol"] != ocand[0]["symbol"]]
                    break
            time.sleep(0.11 if Entrez.api_key else 0.34)
    n_rec, acc = 0, []
    if best and best["symbol"]:
        rorg = org if scoped else best["organism"]       # count where the gene actually is
        nq = (f'{best["symbol"]}[Gene Name] AND {rorg}[Organism]') if rorg else f'{best["symbol"]}[Gene Name]'
        try:
            n_rec = count_hits(nq)
            acc = [a.split()[0] for a, _ in search_fetch_fasta(nq, 3)]
        except Exception:
            pass
    in_org = bool(org) and scoped and bool(best) and best["score"] >= 25
    return dict(found=True, gene=gene, organism=org, scoped=scoped, in_requested_organism=in_org,
                best=best, candidates=cand[:max_candidates], records=n_rec, accessions=acc)


def gene_id_from_accession(acc):
    """Entrez Gene ID linked to a nucleotide accession, via elink (organism-independent).
    This lets the intron check work for a non-model organism that has a RefSeq mRNA but no
    Gene record searchable under the typed organism name. None on failure."""
    if not acc:
        return None
    try:
        recs = _retry_read(Entrez.elink, dbfrom="nucleotide", db="gene", id=acc.strip())
        for r in recs:
            for ls in r.get("LinkSetDb", []):
                links = ls.get("Link", [])
                if links:
                    return links[0]["Id"]
    except Exception:
        pass
    return None


# ---- gene / marker synonym canonicalization ----------------------------------------------------
# One gene is written many ways: "cytochrome b" / "cyt b" / "cyt-b" / "cytb" / "cob" / "MT-CYB". As raw
# NCBI text these match DIFFERENT record sets, so a design or an isolate panel silently changed with the
# user's wording. Every free-text gene query is routed through this table first, so the retrieved set is
# the union of all spellings and no longer depends on which one was typed. Each group is (canonical label,
# [equivalent forms, lowercase]); matching is word-boundary, longest-form-first, case-insensitive, and an
# explicitly field-tagged token ("cytb[Gene]") is left untouched.
_SYNONYM_GROUPS = [
    ("cytochrome b", ["cytochrome b", "cytochrome-b", "cyt b", "cyt-b", "cytb", "cob", "mt-cyb", "mtcyb", "cyb"]),
    ("cytochrome c oxidase subunit I", ["cytochrome c oxidase subunit 1", "cytochrome c oxidase subunit i",
        "cytochrome oxidase subunit 1", "cytochrome oxidase subunit i", "cox1", "coxi", "coi", "co1", "mt-co1", "mtco1"]),
    ("cytochrome c oxidase subunit II", ["cytochrome c oxidase subunit 2", "cytochrome c oxidase subunit ii",
        "cytochrome oxidase subunit 2", "cox2", "coxii", "coii", "co2", "mt-co2"]),
    ("cytochrome c oxidase subunit III", ["cytochrome c oxidase subunit 3", "cytochrome c oxidase subunit iii",
        "cytochrome oxidase subunit 3", "cox3", "coxiii", "coiii", "co3", "mt-co3"]),
    ("NADH dehydrogenase subunit 1", ["nadh dehydrogenase subunit 1", "nad1", "nd1", "mt-nd1"]),
    ("NADH dehydrogenase subunit 2", ["nadh dehydrogenase subunit 2", "nad2", "nd2", "mt-nd2"]),
    ("NADH dehydrogenase subunit 4", ["nadh dehydrogenase subunit 4", "nad4", "nd4", "mt-nd4"]),
    ("NADH dehydrogenase subunit 5", ["nadh dehydrogenase subunit 5", "nad5", "nd5", "mt-nd5"]),
    ("cytochrome c oxidase subunit", []),  # (placeholder, intentionally empty: see specific subunits above)
    ("16S ribosomal RNA", ["16s ribosomal rna", "16s rrna", "16s rdna", "16s"]),
    ("18S ribosomal RNA", ["18s ribosomal rna", "18s rrna", "18s rdna", "18s"]),
    ("28S ribosomal RNA", ["28s ribosomal rna", "28s rrna", "28s rdna", "28s"]),
    ("12S ribosomal RNA", ["12s ribosomal rna", "12s rrna", "12s rdna", "12s"]),
    ("internal transcribed spacer", ["internal transcribed spacer", "its1", "its2", "its region", "its"]),
]
_SYNONYM_GROUPS = [g for g in _SYNONYM_GROUPS if g[1]]   # drop empty placeholders


def _or_terms(forms, field="[All Fields]"):
    """An OR group of synonym forms, each carrying `field`; multi-word / hyphenated forms are quoted."""
    parts = []
    for f in forms:
        t = ('"%s"' % f) if any(ch not in "abcdefghijklmnopqrstuvwxyz0123456789" for ch in f) else f
        parts.append(t + field)
    return "(" + " OR ".join(parts) + ")"


_SYN_FORMS = sorted(((m, gi) for gi, (_lab, ms) in enumerate(_SYNONYM_GROUPS) for m in ms),
                    key=lambda x: -len(x[0]))


def canonicalize_query(query, field="[All Fields]"):
    """Expand recognised gene/marker names to the OR of all their synonyms so an NCBI fetch returns the
    same records regardless of spelling. Returns dict(expanded, taxon, groups): `expanded` is the query
    with each gene phrase replaced by its synonym OR-group (each term carries `field`); `taxon` is the
    query with those phrases stripped (the organism remainder); `groups` is [(canonical_label,
    or_expansion), ...]. A query with no known gene name comes back unchanged (groups empty)."""
    q = query or ""
    low = q.lower()
    consumed = [False] * len(q)
    hits, spans = {}, []
    for form, gi in _SYN_FORMS:
        if gi in hits:
            continue                                   # group already matched once -> same gene, skip
        pat = r"(?<![A-Za-z0-9])" + re.escape(form) + r"(?![A-Za-z0-9])"
        for mt in re.finditer(pat, low):
            s, e = mt.start(), mt.end()
            if any(consumed[s:e]):
                continue
            if e < len(q) and q[e] == "[":             # user field-tagged this token -> respect it
                continue
            for i in range(s, e):
                consumed[i] = True
            hits[gi] = (_SYNONYM_GROUPS[gi][0], _or_terms(_SYNONYM_GROUPS[gi][1], field))
            spans.append((s, e, gi))
            break
    if not spans:
        return dict(expanded=q, taxon=q.strip(), groups=[])
    expanded = q
    for s, e, gi in sorted(spans, key=lambda x: -x[0]):
        expanded = expanded[:s] + hits[gi][1] + expanded[e:]
    taxon = re.sub(r"\s+", " ", "".join(c for i, c in enumerate(q) if not consumed[i])).strip(" ,;")
    groups = [hits[gi] for (_, _, gi) in sorted(spans, key=lambda x: x[0])]
    return dict(expanded=expanded, taxon=taxon, groups=groups)


def search_fetch_fasta(query, n=10):
    """esearch + efetch a nucleotide query -> [(description, sequence), ...].

    Bounds the work three ways so a broad query (a gene name like cox3 that also matches
    genome/chromosome records) can't run a small instance out of memory: (1) an esearch
    sequence-length filter excludes genome-scale records up front, (2) the efetch response
    is STREAM-parsed rather than read whole into memory, and (3) a per-record size skip and
    a total-bases ceiling break early. Before this, a stray Plasmodium genome hit pulled
    ~375 MB through `h.read()` and OOM-killed the worker (HTTP 502)."""
    from Bio import SeqIO
    n = max(1, min(int(n), 40))
    term = f"({canonicalize_query(query)['expanded']}) AND 50:120000[SLEN]"   # canonicalize gene synonyms (cytb == cytochrome b == cob ...) so spelling can't change the record set; drop genome-scale records
    ids = _esearch_read("nucleotide", term, n)["IdList"]
    if not ids:
        return []
    raw = _efetch_text(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")  # SLEN filter caps each record at 120 kb
    cut = raw.find(">")                            # NCBI sometimes prepends a notice/comment that the strict
    raw = raw[cut:] if cut >= 0 else ""            # FASTA parser rejects; drop anything before the first record
    recs, total = [], 0
    for r in SeqIO.parse(StringIO(raw), "fasta"):
        s = str(r.seq)
        if len(s) > 120000:                       # defensive: skip any oversized record that slips through
            continue
        recs.append((r.description, s)); total += len(s)
        if len(recs) >= n or total > 4_000_000:
            break
    return recs


def taxonomy_lineage(name):
    """Resolve a taxon name to (scientific_name, rank, lineage_string) via NCBI
    Taxonomy, or None if not found / offline. Used by the marker recommender so
    it can classify any taxon, not just a hard-coded genus list."""
    try:
        r = _esearch_read("taxonomy", name, 20)
        ids = r.get("IdList") or []
        if not ids:
            return None
        recs = _retry_read(Entrez.efetch, db="taxonomy", id=ids[0], retmode="xml")
        if not recs:
            return None
        rec = recs[0]
        return (rec.get("ScientificName", "") or name,
                rec.get("Rank", "") or "",
                rec.get("Lineage", "") or "")
    except Exception:
        return None


def count_hits(query):
    """Number of nucleotide records matching a query (availability signal). 0 on failure."""
    try:
        r = _esearch_read("nucleotide", query, 0)
        return int(r.get("Count", 0))
    except Exception:
        return None


def genes_for_organism(organism, n=20):
    """Real gene symbols catalogued for an organism in NCBI Gene -- live, specific to exactly what
    was typed (best-effort). Returns [{symbol, name}] in NCBI relevance order; empty for poorly
    annotated taxa or on failure. Lets gene suggestions reflect the actual organism instead of a
    fixed per-clade list."""
    org = (organism or "").strip()
    if not org:
        return []
    try:
        ids = _esearch_read("gene", "%s[Organism]" % org, n, sort="relevance").get("IdList", [])
    except Exception:
        return []
    if not ids:
        return []
    try:
        recs = _esummary_read("gene", ids)
    except Exception:
        return []
    rows = recs
    if isinstance(recs, dict):                       # gene esummary nests under DocumentSummarySet
        rows = (recs.get("DocumentSummarySet") or {}).get("DocumentSummary") or []
    out, seen = [], set()
    for d in rows:
        try:
            sym = str(d.get("Name") or d.get("NomenclatureSymbol") or "").strip()
            desc = str(d.get("Description") or d.get("NomenclatureName")
                       or d.get("OtherDesignations") or "").strip()
        except Exception:
            continue
        if not sym or sym.lower() in seen:
            continue
        seen.add(sym.lower())
        if desc and "|" in desc:                     # OtherDesignations is pipe-delimited; take the first
            desc = desc.split("|", 1)[0].strip()
        out.append({"symbol": sym, "name": desc or sym})
    return out


_MARKER_HINT = {"cytb", "cob", "cytochrome", "coi", "co1", "cox1", "coii", "cox2", "coiii", "cox3",
                "nad1", "nad5", "nad4", "nad2", "atp6", "atp8", "18s", "28s", "16s", "23s", "12s",
                "its", "its1", "its2", "rrna", "ribosomal", "rbcl", "matk", "trnl", "trnh", "psba",
                "d-loop", "dhfr", "msp1", "msp2", "ama1", "csp"}


def search_genomes(query, retmax=40, gene=None):
    """List isolate-level records for a taxon, for an inclusivity/exclusivity panel.

    Two modes, picked automatically:
      * WHOLE-GENOME (default) -- complete genomes / chromosomes >= 100 kb (bacterial / viral isolate
        panels, e.g. Salmonella enterica strains).
      * MARKER -- when the query (or `gene`) names a barcode/marker (cytb, COI, 18S, ITS, ...), pulls many
        short single-gene records across the taxon (e.g. every Plasmodium cytb on file) instead of the
        almost-nonexistent assembled genomes. This is what an apicomplexan / haemosporidian lineage panel
        actually needs.

    Metadata only via esummary -- NO sequence fetched here. Returns [{acc, title, slen}]; the picker shows
    each title so the user confirms exactly what goes into the run."""
    q = (query or "").strip()
    if not q:
        return []
    # Canonicalize the gene/marker so "cytb", "cyt b" and "cytochrome b" pull the SAME records (the OR of
    # every spelling) — and so the full name triggers MARKER mode just like the abbreviation does, instead
    # of silently dropping to whole-genome mode and returning a different (wrong) set for e.g. haemosporidians.
    marker_term, taxon = "", q
    cg = canonicalize_query((gene or "").strip() or q)
    if cg["groups"]:
        marker_term = cg["groups"][0][1]                          # OR group, already [All Fields]-tagged
        taxon = q if (gene or "").strip() else (cg["taxon"] or q)
    elif (gene or "").strip():
        marker_term, taxon = gene.strip() + "[All Fields]", q     # caller named a marker we don't expand
    else:
        toks = q.split()                                          # fall back to single-token hints (matK, rpoB, D-loop, ...)
        mk = [t for t in toks if t.lower().strip(",.;") in _MARKER_HINT]
        if mk:
            marker_term = mk[0] + "[All Fields]"
            taxon = " ".join(t for t in toks if t.lower().strip(",.;") not in _MARKER_HINT).strip() or q
    if marker_term:
        tiers = [
            f'({taxon}[Organism]) AND {marker_term} AND 150:30000[SLEN]',
            f'({taxon}) AND {marker_term} AND 150:30000[SLEN]',
            f'({taxon}[Organism]) AND 150:30000[SLEN] NOT genome[Title]',
        ]
    else:
        tiers = [
            f'({q}[Organism]) AND complete genome[Title] AND 500000:25000000[SLEN] NOT plasmid[Title]',
            f'({q}[Organism]) AND (chromosome[Title] OR complete sequence[Title]) AND 500000:25000000[SLEN] NOT plasmid[Title]',
            f'({q}) AND 100000:25000000[SLEN] NOT plasmid[Title]',
        ]
    ids = []
    for term in tiers:
        ids = _esearch_read("nucleotide", term, int(retmax))["IdList"]
        if ids:
            break
    if not ids:
        return []
    recs = _esummary_read("nucleotide", ids)
    out = []
    for d in recs:
        out.append(dict(acc=str(d.get("AccessionVersion") or d.get("Caption") or ""),
                        title=str(d.get("Title") or ""), slen=int(d.get("Slen") or 0)))
    out.sort(key=lambda x: x["title"])
    return out


def fetch_one(acc):
    """Fetch a single accession as (title, sequence). Streamed for one record; a leading NCBI
    notice/comment is stripped. Used by isolate validation, which scans one genome then frees it
    (peak memory = one record, so a 40-isolate panel never holds 40 genomes at once)."""
    raw = _efetch_text(db="nucleotide", id=str(acc).strip(), rettype="fasta", retmode="text")
    cut = raw.find(">")
    recs = list(SeqIO.parse(StringIO(raw[cut:] if cut >= 0 else ""), "fasta"))
    if not recs:
        return (str(acc), "")
    return (recs[0].description, str(recs[0].seq).upper())
