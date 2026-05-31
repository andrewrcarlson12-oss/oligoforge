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

Entrez.email = "set-me@example.com"   # set to your address; NCBI requires it
Entrez.api_key = None                  # optional: set for 10 req/s instead of 3


def _records(ids, rettype="fasta"):
    h = Entrez.efetch(db="nucleotide", id=list(ids), rettype=rettype, retmode="text")
    recs = list(SeqIO.parse(h, "genbank" if rettype == "gb" else "fasta"))
    h.close()
    return recs


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


def search_fetch_fasta(query, n=10):
    """esearch + efetch a nucleotide query -> [(description, sequence), ...].
    Caps the request count and skips oversized records (genomic contigs /
    chromosomes) with a total-bases ceiling, so a stray large hit can't run the
    server out of memory on a small instance."""
    from io import StringIO
    from Bio import SeqIO
    n = max(1, min(int(n), 40))
    h = Entrez.esearch(db="nucleotide", term=query, retmax=n)
    ids = Entrez.read(h)["IdList"]; h.close()
    if not ids:
        return []
    h = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
    recs, total = [], 0
    for r in SeqIO.parse(StringIO(h.read()), "fasta"):
        s = str(r.seq)
        if len(s) > 120000:                 # skip genomic contigs / chromosomes
            continue
        recs.append((r.description, s)); total += len(s)
        if total > 4_000_000:               # ~4 Mb total cap
            break
    h.close()
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
