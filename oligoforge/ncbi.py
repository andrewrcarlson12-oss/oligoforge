"""NCBI retrieval via Entrez (biopython): fetch by accession, search a gene in
an organism, pull all mRNA isoforms, and find the region common to every isoform
(so an assay reads total transcript, not one variant)."""
import time, difflib, socket, os
from Bio import Entrez, SeqIO

# Bound every NCBI HTTP call so a genuinely stalled connection eventually fails
# instead of hanging forever — but generously, because a remote blastn against
# nt is queued server-side at NCBI and a single poll can legitimately take
# several minutes (and a large efetch can too). 600 s is the default so BLAST
# from the Target->assay tab completes instead of being killed mid-poll.
# Override with the env var OLIGOFORGE_NCBI_TIMEOUT (seconds); raise it higher
# if your network is slow. Async server sockets are non-blocking and ignore this.
try:
    NCBI_TIMEOUT = float(os.environ.get("OLIGOFORGE_NCBI_TIMEOUT", "600"))
except ValueError:
    NCBI_TIMEOUT = 600.0
socket.setdefaulttimeout(NCBI_TIMEOUT)
from io import StringIO

Entrez.email = "set-me@example.com"   # set to your address; NCBI requires it
Entrez.api_key = None                  # optional: set for 10 req/s instead of 3


def _records(ids, rettype="fasta"):
    h = Entrez.efetch(db="nucleotide", id=list(ids), rettype=rettype, retmode="text")
    try:
        raw = h.read()
    finally:
        h.close()
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
        h = Entrez.esearch(db="nucleotide", term=q, retmax=retmax)
        r = Entrez.read(h); h.close()
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
        h = Entrez.esearch(db="gene", term=term)
        r = Entrez.read(h); h.close()
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
            h = Entrez.esummary(db="gene", id=",".join(ids[:50])); s = Entrez.read(h); h.close()
            docs = list(s["DocumentSummarySet"]["DocumentSummary"])
        except Exception:
            return None, 0.0
        if not docs:
            return None, 0.0
        docs.sort(key=lambda d: _gene_match_score(d, gene), reverse=True)
        return docs[0], _gene_match_score(docs[0], gene)

    for term, in_org in tiers:
        try:
            h = Entrez.esearch(db="gene", term=term, retmax=50)
            ids = Entrez.read(h)["IdList"]; h.close()
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
                h = Entrez.esearch(db="gene", term=term, retmax=max_candidates)
                ids = Entrez.read(h)["IdList"]; h.close()
            except Exception:
                ids = []
            if ids:
                scoped = True; break
            time.sleep(0.11 if Entrez.api_key else 0.34)
    if not ids:
        for term in (f'{gene}[Gene Name]', f'{gene}[All Fields] AND alive[prop]'):
            try:
                h = Entrez.esearch(db="gene", term=term, retmax=max_candidates)
                ids = Entrez.read(h)["IdList"]; h.close()
            except Exception:
                ids = []
            if ids:
                break
            time.sleep(0.11 if Entrez.api_key else 0.34)
    if not ids:
        return dict(found=False, gene=gene, organism=org,
                    error="NCBI Gene has no record matching \u201c%s\u201d%s" % (gene, (" in " + org) if org else ""))
    try:
        h = Entrez.esummary(db="gene", id=",".join(ids[:max_candidates]))
        docs = list(Entrez.read(h)["DocumentSummarySet"]["DocumentSummary"]); h.close()
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
                h = Entrez.esearch(db="gene", term=term, retmax=max_candidates)
                oids = Entrez.read(h)["IdList"]; h.close()
            except Exception:
                oids = []
            if oids:
                try:
                    h = Entrez.esummary(db="gene", id=",".join(oids[:max_candidates]))
                    odocs = list(Entrez.read(h)["DocumentSummarySet"]["DocumentSummary"]); h.close()
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
        h = Entrez.elink(dbfrom="nucleotide", db="gene", id=acc.strip())
        recs = Entrez.read(h); h.close()
        for r in recs:
            for ls in r.get("LinkSetDb", []):
                links = ls.get("Link", [])
                if links:
                    return links[0]["Id"]
    except Exception:
        pass
    return None


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
    term = f"({query}) AND 50:120000[SLEN]"      # keep genes / mitogenomes, drop chromosome & genome records
    h = Entrez.esearch(db="nucleotide", term=term, retmax=n)
    ids = Entrez.read(h)["IdList"]; h.close()
    if not ids:
        return []
    h = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
    try:
        raw = h.read()                            # SLEN filter caps each record at 120 kb, so this is bounded
    finally:
        h.close()
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
        h = Entrez.esearch(db="taxonomy", term=name)
        r = Entrez.read(h); h.close()
        ids = r.get("IdList") or []
        if not ids:
            return None
        h = Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml")
        recs = Entrez.read(h); h.close()
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
        h = Entrez.esearch(db="nucleotide", term=query, retmax=0)
        r = Entrez.read(h); h.close()
        return int(r.get("Count", 0))
    except Exception:
        return None


def search_genomes(query, retmax=40):
    """List isolate-level nucleotide records (complete genomes / chromosomes) for a taxon.
    Metadata only via esummary -- NO sequence is fetched here, so it stays light even for a
    broad taxon. Returns [{acc, title, slen}]. The picker shows title (organism+strain) so the
    user confirms each isolate before it goes into a validation run."""
    q = (query or "").strip()
    if not q:
        return []
    tiers = [
        f'({q}[Organism]) AND complete genome[Title] AND 500000:25000000[SLEN] NOT plasmid[Title]',
        f'({q}[Organism]) AND (chromosome[Title] OR complete sequence[Title]) AND 500000:25000000[SLEN] NOT plasmid[Title]',
        f'({q}) AND 100000:25000000[SLEN] NOT plasmid[Title]',
    ]
    ids = []
    for term in tiers:
        h = Entrez.esearch(db="nucleotide", term=term, retmax=int(retmax)); ids = Entrez.read(h)["IdList"]; h.close()
        if ids:
            break
    if not ids:
        return []
    h = Entrez.esummary(db="nucleotide", id=",".join(ids)); recs = Entrez.read(h); h.close()
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
    h = Entrez.efetch(db="nucleotide", id=str(acc).strip(), rettype="fasta", retmode="text")
    try:
        raw = h.read()
    finally:
        h.close()
    cut = raw.find(">")
    recs = list(SeqIO.parse(StringIO(raw[cut:] if cut >= 0 else ""), "fasta"))
    if not recs:
        return (str(acc), "")
    return (recs[0].description, str(recs[0].seq).upper())
