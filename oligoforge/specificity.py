"""Specificity + intron/exon-junction checks — the parts a browser sandbox cannot do.

intron_check : reads RefSeq exon structure (NCBI gene_table) and reports whether an
               amplicon spans an exon-exon junction. Within one exon => gDNA will
               co-amplify (bad for an expression assay). Spanning a junction => safe.
blast_remote : NCBI BLAST over the network (no install; slower, queued).
blast_local  : subprocess to a local blastn DB (fast/offline; needs BLAST+ + genome).
"""
import re, shutil, subprocess
from Bio import Entrez
from .ncbi import fetch_accessions
try:
    from Bio.Blast import NCBIWWW, NCBIXML
except Exception:                       # pragma: no cover
    NCBIWWW = NCBIXML = None


_RC_MAP = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "R": "Y", "Y": "R", "S": "S",
           "W": "W", "K": "M", "M": "K", "B": "V", "D": "H", "H": "D", "V": "B", "N": "N"}
_LOC_RE = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "[AG]", "Y": "[CT]",
           "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
           "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def _rc_iupac(seq):
    return "".join(_RC_MAP.get(b, "N") for b in reversed(seq.upper()))


def _locate(primer, seq):
    """0-based start of primer within seq (IUPAC-aware exact match), or None."""
    if not primer or not seq:
        return None
    pat = "".join(_LOC_RE.get(b, b) for b in primer.upper().replace("U", "T"))
    try:
        m = re.search(pat, seq.upper().replace("U", "T"))
    except re.error:
        return None
    return m.start() if m else None


def exon_junctions_mrna(gene, organism, mrna_acc=None):
    """Exon-junction positions in mRNA coordinates from NCBI gene_table.
    Sections are anchored to their transcript accession so the right isoform
    is parsed. Junctions = cumulative spliced-exon lengths."""
    from .ncbi import gene_id, gene_id_from_accession
    gid = gene_id_from_accession(mrna_acc) if mrna_acc else None
    if not gid:
        gid = gene_id(gene, organism)
    if not gid:
        where = (" in " + organism) if organism else ""
        return None, ("NCBI has no gene/exon annotation for %s%s (common for non-model "
                      "species). Give a RefSeq mRNA accession from this or a closely related "
                      "annotated species, or mark the gDNA control as handled." %
                      (gene or "this gene", where)), None
    h = Entrez.efetch(db="gene", id=gid, rettype="gene_table", retmode="text")
    txt = h.read(); h.close()

    # locate each "Exon table for  mRNA  <ACC>" and the span until the next one
    heads = list(re.finditer(r"Exon table for\s+mRNA\s+(\S+)", txt))
    if not heads:
        return None, "no exon tables found", None
    sections = []
    for i, m in enumerate(heads):
        s = m.end()
        e = heads[i + 1].start() if i + 1 < len(heads) else len(txt)
        sections.append((m.group(1), txt[s:e]))

    chosen = None
    if mrna_acc:
        key = mrna_acc.split(".")[0]
        for acc, body in sections:
            if acc.split(".")[0] == key:
                chosen = (acc, body); break
    if chosen is None:
        chosen = sections[0]
    acc, body = chosen

    exon_lens = []
    for line in body.splitlines():
        if not re.match(r"^\s*\d+-\d+", line):
            continue
        toks = [t for t in re.split(r"\t+|\s{2,}", line.strip()) if t]
        for t in toks:                       # exon length = first bare integer
            if re.fullmatch(r"\d+", t):
                exon_lens.append(int(t)); break
    if not exon_lens:
        return None, f"no exon lengths parsed for {acc}", acc

    junctions, cum = [], 0
    for L in exon_lens[:-1]:
        cum += L
        junctions.append(cum)
    return junctions, f"{acc}: {len(exon_lens)} exons, spliced length {sum(exon_lens)}", acc


def intron_check(gene, organism, amp_start=None, amp_end=None, mrna_acc=None,
                 forward=None, reverse=None):
    """Report whether an amplicon spans an exon-exon junction.

    amp_start/amp_end are 1-based mRNA coordinates. If they are not supplied but the
    forward and reverse primers are, the amplicon is located on the mRNA automatically
    (forward and the reverse-complement of reverse), so the coordinates fill themselves.
    """
    junctions, info, acc = exon_junctions_mrna(gene, organism, mrna_acc)
    if junctions is None:
        return dict(ok=None, info=info, junctions=None, spanned=[])

    located = False
    if (amp_start is None or amp_end is None) and forward and reverse:
        try:
            recs = fetch_accessions([acc], "fasta")
            mrna = str(recs[0].seq) if recs else ""
        except Exception as e:
            mrna = ""
        if not mrna:
            return dict(ok=None, junctions=junctions, spanned=[],
                        info=f"{info}; couldn't fetch mRNA {acc} to locate the amplicon — enter amp start/end manually")
        fs = _locate(forward, mrna)
        rrc = _rc_iupac(reverse)
        rs = _locate(rrc, mrna)
        if fs is None or rs is None:
            miss = " and ".join(n for n, ok in (("forward", fs is not None), ("reverse", rs is not None)) if not ok)
            return dict(ok=None, junctions=junctions, spanned=[], amp_located=False,
                        info=f"{info}; could not locate the {miss} primer in {acc} (different isoform?) — enter amp start/end manually")
        amp_start = fs + 1
        amp_end = rs + len(rrc)
        located = True

    if amp_start is None or amp_end is None:
        return dict(ok=None, junctions=junctions, spanned=[],
                    info=f"{info}; enter amplicon start/end (mRNA coordinates), or paste both primers to auto-locate")

    spanned = [j for j in junctions if amp_start <= j < amp_end]
    return dict(ok=bool(spanned), info=info, junctions=junctions, spanned=spanned,
                amp_start=amp_start, amp_end=amp_end, amplicon=amp_end - amp_start + 1,
                amp_located=located,
                verdict=("spans an exon-exon junction (gDNA will NOT co-amplify)"
                         if spanned else
                         "lies WITHIN a single exon (gDNA can co-amplify -> DNase-treat RNA "
                         "or redesign across a junction)"))


def blast_remote(seq, organism=None, hitlist=10, program="blastn", db="nt"):
    if NCBIWWW is None:
        return dict(error="biopython Blast unavailable")
    entrez_q = f"{organism}[Organism]" if organism else None
    h = NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist,
                       entrez_query=entrez_q, megablast=True)
    rec = NCBIXML.read(h); h.close()
    hits = []
    for al in rec.alignments:
        hsp = al.hsps[0]
        hits.append(dict(title=al.title[:80], length=al.length,
                         identity=hsp.identities, align_len=hsp.align_length,
                         evalue=hsp.expect))
    return dict(query_len=len(seq), n_hits=len(hits), hits=hits)


def blast_local(seq, db_path, hitlist=10):
    if not shutil.which("blastn"):
        return dict(error="blastn not on PATH (install BLAST+)")
    p = subprocess.run(
        ["blastn", "-db", db_path, "-outfmt", "6 sacc pident length evalue stitle",
         "-max_target_seqs", str(hitlist)],
        input=f">q\n{seq}\n", capture_output=True, text=True)
    rows = [l.split("\t") for l in p.stdout.strip().splitlines() if l]
    return dict(query_len=len(seq), n_hits=len(rows), hits=rows, stderr=p.stderr[:200])


# ---- in-silico PCR: do the two primers actually bracket a product? ----
def _blast_remote_coords(seq, organism=None, hitlist=20, db="nt"):
    if NCBIWWW is None:
        return dict(error="biopython Blast unavailable")
    eq = f"{organism}[Organism]" if organism else None
    h = NCBIWWW.qblast("blastn", db, seq, hitlist_size=hitlist, entrez_query=eq, megablast=False)
    rec = NCBIXML.read(h); h.close()
    hits = []
    for al in rec.alignments:
        for hsp in al.hsps:
            ss, se = hsp.sbjct_start, hsp.sbjct_end
            lo, hi = (ss, se) if ss <= se else (se, ss)
            hits.append(dict(subject=al.accession, lo=lo, hi=hi,
                             strand="+" if ss <= se else "-",
                             q3=hsp.query_end >= len(seq) - 1,
                             ident=hsp.identities, alen=hsp.align_length))
    return dict(n_hits=len(hits), hits=hits)


def _blast_local_coords(seq, db_path, hitlist=20):
    if not shutil.which("blastn"):
        return dict(error="blastn not on PATH (install BLAST+)")
    fmt = "6 sacc pident length evalue sstart send qstart qend sstrand"
    p = subprocess.run(["blastn", "-db", db_path, "-outfmt", fmt,
                        "-max_target_seqs", str(hitlist)],
                       input=f">q\n{seq}\n", capture_output=True, text=True)
    hits = []
    for l in p.stdout.strip().splitlines():
        c = l.split("\t")
        if len(c) < 9:
            continue
        ss, se = int(c[4]), int(c[5]); lo, hi = sorted((ss, se))
        hits.append(dict(subject=c[0], lo=lo, hi=hi,
                         strand="+" if c[8].lower().startswith("plus") else "-",
                         q3=int(c[7]) >= len(seq) - 1, ident=c[1], alen=int(c[2])))
    return dict(n_hits=len(hits), hits=hits, stderr=p.stderr[:200])


def epcr(hits, min_product=40, max_product=3000, require_3prime=True):
    """Predicted amplicons from a pooled hit list: two hits on one subject, opposite
    strands, convergent (plus-strand hit upstream of minus-strand hit), product size
    in range. hits: [{primer, subject, lo, hi, strand('+'/'-'), q3}]. The intended
    target shows up as a product of the expected size; anything else is a flag."""
    by = {}
    for h in hits:
        by.setdefault(h["subject"], []).append(h)
    products = []
    for subj, hs in by.items():
        plus = [h for h in hs if h["strand"] == "+" and (not require_3prime or h.get("q3", True))]
        minus = [h for h in hs if h["strand"] == "-" and (not require_3prime or h.get("q3", True))]
        for P in plus:
            for M in minus:
                if P["lo"] <= M["hi"] and P["hi"] <= M["lo"] + 5:
                    size = M["hi"] - P["lo"] + 1
                    if min_product <= size <= max_product:
                        products.append(dict(subject=subj, size=size,
                                             left=P["primer"], right=M["primer"],
                                             span=[P["lo"], M["hi"]]))
    products.sort(key=lambda x: x["size"])
    return products


def in_silico_pcr(forward, reverse, mode="remote", db="nt", db_path=None,
                  organism=None, min_product=40, max_product=3000):
    def run(seq):
        return (_blast_local_coords(seq, db_path or "") if mode == "local"
                else _blast_remote_coords(seq, organism=organism, db=db))
    fr, rr = run(forward.upper().strip()), run(reverse.upper().strip())
    if fr.get("error") or rr.get("error"):
        return dict(error=fr.get("error") or rr.get("error"))
    hits = ([dict(primer="F", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in fr["hits"]] +
            [dict(primer="R", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in rr["hits"]])
    prod = epcr(hits, min_product, max_product)
    return dict(forward_hits=fr["n_hits"], reverse_hits=rr["n_hits"],
                n_products=len(prod), products=prod[:25])
