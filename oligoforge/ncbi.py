"""NCBI retrieval via Entrez (biopython): fetch by accession, search a gene in
an organism, pull all mRNA isoforms, and find the region common to every isoform
(so an assay reads total transcript, not one variant)."""
import time, difflib, socket, os
from Bio import Entrez, SeqIO

# Bound every NCBI HTTP call so a genuinely stalled connection eventually fails
# instead of hanging forever — but generously, because an efetch of several or
# large records can legitimately take a while. Override with the env var
# OLIGOFORGE_NCBI_TIMEOUT (seconds). Async server sockets are non-blocking and
# ignore this.
try:
    NCBI_TIMEOUT = float(os.environ.get("OLIGOFORGE_NCBI_TIMEOUT", "120"))
except ValueError:
    NCBI_TIMEOUT = 120.0
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
        time.sleep(0.34)
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
    h = Entrez.esearch(db="gene", term=f'{gene}[Gene Name] AND {organism}[Organism]')
    r = Entrez.read(h); h.close()
    return r["IdList"][0] if r["IdList"] else None


def search_fetch_fasta(query, n=10):
    """esearch + efetch a nucleotide query -> [(description, sequence), ...]."""
    from io import StringIO
    from Bio import SeqIO
    h = Entrez.esearch(db="nucleotide", term=query, retmax=n)
    ids = Entrez.read(h)["IdList"]; h.close()
    if not ids:
        return []
    h = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
    recs = [(r.description, str(r.seq)) for r in SeqIO.parse(StringIO(h.read()), "fasta")]
    h.close()
    return recs
