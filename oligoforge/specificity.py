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
    _UNKNOWN = "could not determine junction spanning (exon structure unavailable)"
    junctions, info, acc = exon_junctions_mrna(gene, organism, mrna_acc)
    if junctions is None:
        return dict(ok=None, info=info, junctions=None, spanned=[], verdict=_UNKNOWN)

    located = False
    if (amp_start is None or amp_end is None) and forward and reverse:
        try:
            recs = fetch_accessions([acc], "fasta")
            mrna = str(recs[0].seq) if recs else ""
        except Exception as e:
            mrna = ""
        if not mrna:
            return dict(ok=None, junctions=junctions, spanned=[], verdict=_UNKNOWN,
                        info=f"{info}; couldn't fetch mRNA {acc} to locate the amplicon — enter amp start/end manually")
        fs = _locate(forward, mrna)
        rrc = _rc_iupac(reverse)
        rs = _locate(rrc, mrna)
        if fs is None or rs is None:
            miss = " and ".join(n for n, ok in (("forward", fs is not None), ("reverse", rs is not None)) if not ok)
            return dict(ok=None, junctions=junctions, spanned=[], amp_located=False, verdict=_UNKNOWN,
                        info=f"{info}; could not locate the {miss} primer in {acc} (different isoform?) — enter amp start/end manually")
        amp_start = fs + 1
        amp_end = rs + len(rrc)
        located = True

    if amp_start is None or amp_end is None:
        return dict(ok=None, junctions=junctions, spanned=[], verdict=_UNKNOWN,
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
    seq = (seq or "").strip()
    entrez_q = f"{organism}[Organism]" if organism else None
    h = NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist,
                       entrez_query=entrez_q, **_short_blast_kw(seq))
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
    seq = (seq or "").strip()
    cmd = ["blastn", "-db", db_path, "-outfmt", "6 sacc pident length evalue stitle",
           "-max_target_seqs", str(hitlist)]
    if len(seq) <= 30:                       # short-query params so primers return hits
        cmd += ["-word_size", "7", "-evalue", "1000", "-reward", "1", "-penalty", "-3",
                "-gapopen", "5", "-gapextend", "2", "-dust", "no"]
    p = subprocess.run(cmd, input=f">q\n{seq}\n", capture_output=True, text=True)
    rows = [l.split("\t") for l in p.stdout.strip().splitlines() if l]
    return dict(query_len=len(seq), n_hits=len(rows), hits=rows, stderr=p.stderr[:200])


# ---- in-silico PCR: do the two primers actually bracket a product? ----
def _short_blast_kw(seq):
    """Parameters so a ~20 nt primer actually returns hits against nt (NCBI's
    'blastn-short' behaviour): small word size, no low-complexity masking, large
    E-value. Without these, default blastn silently drops short queries and reports
    zero hits — which is why primers came back empty while a long sequence worked."""
    if len(seq) <= 30:
        return dict(word_size=7, expect=1000.0, nucl_reward=1, nucl_penalty=-3,
                    gapcosts="5 2", filter="F", megablast=False)
    return dict(megablast=True)


def _hits_from_record(rec, qlen):
    hits = []
    for al in rec.alignments:
        for hsp in al.hsps:
            ss, se = hsp.sbjct_start, hsp.sbjct_end
            lo, hi = (ss, se) if ss <= se else (se, ss)
            hits.append(dict(subject=al.accession, lo=lo, hi=hi,
                             strand="+" if ss <= se else "-",
                             q3=hsp.query_end >= qlen - 1,
                             ident=hsp.identities, alen=hsp.align_length))
    return hits


def _blast_remote_pair(fwd, rev, organism=None, db="nt", hitlist=50):
    """BLAST BOTH primers in a single NCBI job (one network round-trip instead of two)
    using short-query parameters. Returns ({'F':[...], 'R':[...]}, error)."""
    if NCBIWWW is None:
        return None, "biopython Blast unavailable"
    eq = f"{organism}[Organism]" if organism else None
    query = ">F\n%s\n>R\n%s\n" % (fwd, rev)
    kw = _short_blast_kw(fwd if len(fwd) >= len(rev) else rev)
    try:
        h = NCBIWWW.qblast("blastn", db, query, hitlist_size=hitlist, entrez_query=eq, **kw)
        recs = list(NCBIXML.parse(h)); h.close()
    except Exception as e:
        return None, f"NCBI BLAST error: {e}"
    out, lengths = {"F": [], "R": []}, {"F": len(fwd), "R": len(rev)}
    for i, rec in enumerate(recs[:2]):
        lab = ("F", "R")[i]
        out[lab] = _hits_from_record(rec, getattr(rec, "query_letters", 0) or lengths[lab])
    return out, None


def _blast_local_pair(fwd, rev, db_path, hitlist=50):
    if not shutil.which("blastn"):
        return None, "blastn not on PATH (install BLAST+)"
    fmt = "6 qseqid sacc pident length evalue sstart send qstart qend sstrand"
    q = ">F\n%s\n>R\n%s\n" % (fwd, rev)
    p = subprocess.run(["blastn", "-db", db_path, "-outfmt", fmt, "-word_size", "7",
                        "-evalue", "1000", "-reward", "1", "-penalty", "-3", "-gapopen", "5",
                        "-gapextend", "2", "-dust", "no", "-max_target_seqs", str(hitlist)],
                       input=q, capture_output=True, text=True)
    out, ql = {"F": [], "R": []}, {"F": len(fwd), "R": len(rev)}
    for l in p.stdout.strip().splitlines():
        c = l.split("\t")
        if len(c) < 10:
            continue
        lab = c[0] if c[0] in ("F", "R") else "F"
        ss, se = int(c[5]), int(c[6]); lo, hi = sorted((ss, se))
        out[lab].append(dict(subject=c[1], lo=lo, hi=hi,
                             strand="+" if c[9].lower().startswith("plus") else "-",
                             q3=int(c[8]) >= ql[lab] - 1, ident=c[2], alen=int(c[3])))
    return out, None


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
    fwd, rev = forward.upper().strip(), reverse.upper().strip()
    pair, err = (_blast_local_pair(fwd, rev, db_path or "") if mode == "local"
                 else _blast_remote_pair(fwd, rev, organism=organism, db=db))
    if err:
        return dict(error=err)
    F, R = pair.get("F", []), pair.get("R", [])
    hits = ([dict(primer="F", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in F] +
            [dict(primer="R", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in R])
    prod = epcr(hits, min_product, max_product)
    sizes = [p["size"] for p in prod]
    summ = {}
    if sizes:
        ss = sorted(sizes)
        med = ss[len(ss) // 2]
        tol = max(10, int(round(0.10 * med)))
        summ = dict(size_min=ss[0], size_max=ss[-1], size_modal=med, size_tol=tol,
                    n_on_size=sum(1 for s in sizes if abs(s - med) <= tol),
                    n_subjects=len({p["subject"] for p in prod}))
    return dict(forward_hits=len(F), reverse_hits=len(R), organism=(organism or None),
                n_products=len(prod), products=prod[:25], **summ)


def blast_summary(seq, mode="remote", db="nt", db_path=None, organism=None, top=40):
    """BLAST one sequence and return a ranked hit list (top N) with accession, title,
    % identity, query coverage, E-value and bit score — a readable BLAST report. Short
    sequences (primers/probes) automatically use short-query parameters."""
    seq = "".join(c for c in (seq or "").upper().replace("U", "T") if c.isalpha())
    if len(seq) < 10:
        return dict(error="sequence too short to BLAST (need >= 10 nt)")
    top = max(1, min(int(top), 100))
    if mode == "local":
        if not shutil.which("blastn"):
            return dict(error="blastn not on PATH (install BLAST+)")
        fmt = "6 sacc pident length evalue bitscore qstart qend sstrand stitle"
        cmd = ["blastn", "-db", db_path or "", "-outfmt", fmt, "-max_target_seqs", str(top)]
        if len(seq) <= 30:
            cmd += ["-word_size", "7", "-evalue", "1000", "-reward", "1", "-penalty", "-3",
                    "-gapopen", "5", "-gapextend", "2", "-dust", "no"]
        p = subprocess.run(cmd, input=f">q\n{seq}\n", capture_output=True, text=True)
        hits = []
        for l in p.stdout.strip().splitlines():
            c = l.split("\t")
            if len(c) < 9:
                continue
            hits.append(dict(accession=c[0], pident=float(c[1]), length=int(c[2]), evalue=c[3],
                             bitscore=c[4], strand=c[7],
                             qcov=round(100 * (abs(int(c[6]) - int(c[5])) + 1) / len(seq), 1),
                             title=(c[9] if len(c) > 9 else c[0])[:95]))
        return dict(query_len=len(seq), n_hits=len(hits), hits=hits[:top])
    if NCBIWWW is None:
        return dict(error="biopython Blast unavailable")
    eq = f"{organism}[Organism]" if organism else None
    try:
        h = NCBIWWW.qblast("blastn", db, seq, hitlist_size=top, entrez_query=eq, **_short_blast_kw(seq))
        rec = NCBIXML.read(h); h.close()
    except Exception as e:
        return dict(error=f"NCBI BLAST error: {e}")
    qlen = getattr(rec, "query_letters", 0) or len(seq)
    hits = []
    for al in rec.alignments:
        hsp = al.hsps[0]
        hits.append(dict(accession=al.accession, title=al.hit_def[:95], length=al.length,
                         pident=round(100 * hsp.identities / max(1, hsp.align_length), 1),
                         qcov=round(100 * hsp.align_length / max(1, qlen), 1),
                         evalue=("%.0e" % hsp.expect if hsp.expect else "0"),
                         bitscore=round(hsp.bits), strand="plus" if hsp.sbjct_start <= hsp.sbjct_end else "minus"))
    return dict(query_len=qlen, n_hits=len(hits), hits=hits[:top])


# ---- probe-aware assay specificity (the whole TaqMan assay vs NCBI, not just the primers) ----
def _blast_remote_set(seqs, organism=None, db="nt", hitlist=50):
    """BLAST several oligos (ordered [(label, seq), ...]) in ONE NCBI job, short-query params.
    Returns ({label: [hits]}, error). Records come back in query order."""
    if NCBIWWW is None:
        return None, "biopython Blast unavailable"
    seqs = [(l, (s or "").strip()) for l, s in seqs if (s or "").strip()]
    if not seqs:
        return None, "no oligos to BLAST"
    eq = f"{organism}[Organism]" if organism else None
    query = "".join(">%s\n%s\n" % (l, s) for l, s in seqs)
    longest = max(seqs, key=lambda x: len(x[1]))[1]
    try:
        h = NCBIWWW.qblast("blastn", db, query, hitlist_size=hitlist,
                           entrez_query=eq, **_short_blast_kw(longest))
        recs = list(NCBIXML.parse(h)); h.close()
    except Exception as e:
        return None, f"NCBI BLAST error: {e}"
    out = {l: [] for l, _ in seqs}
    for i, rec in enumerate(recs[:len(seqs)]):
        lab, sq = seqs[i]
        out[lab] = _hits_from_record(rec, getattr(rec, "query_letters", 0) or len(sq))
    return out, None


def assay_verdict(products, probe_hits, size_tol_frac=0.10):
    """Pure logic (no network): annotate predicted products with on-size / probe-binding and
    summarise the off-target picture for a TaqMan assay.

    products   : from epcr() -> [{subject, size, left, right, span:[lo,hi]}, ...]
    probe_hits : [{subject, lo, hi, ...}] for the probe (may be empty). A probe "binds" a product
                 when a probe hit lies wholly inside that product's span on the same subject --
                 i.e. the product would actually generate fluorescence, the real false-positive risk.

    The intended amplicon is taken to be the modal product size; off-size products flag possible
    mispriming. A genus assay legitimately hits many subjects at the on-size length, so subject
    count alone is never treated as failure -- only off-size products and probe cross-binding are.
    """
    products = [dict(p) for p in products]
    sizes = sorted(p["size"] for p in products)
    modal = sizes[len(sizes) // 2] if sizes else None
    tol = max(10, int(round(size_tol_frac * modal))) if modal else 0
    phits = probe_hits or []
    n_off = n_probe = 0
    for p in products:
        p["on_size"] = (modal is not None and abs(p["size"] - modal) <= tol)
        if not p["on_size"]:
            n_off += 1
        sp = p.get("span") or [0, 0]
        binds = any(h.get("subject") == p["subject"] and sp[0] <= h.get("lo", 1e18) and h.get("hi", -1) <= sp[1]
                    for h in phits)
        p["probe_binds"] = bool(binds)
        if binds and not p["on_size"]:
            n_probe += 1
    n_subj = len({p["subject"] for p in products})
    if n_probe:
        verdict = ("probe binds inside %d off-size predicted product(s) -- inspect those subjects for "
                   "cross-reactivity" % n_probe)
    elif n_off:
        verdict = ("%d off-size product(s) predicted -- possible mispriming; the intended amplicon is the "
                   "~%d bp product" % (n_off, modal))
    elif products:
        verdict = "no off-size or probe-cross-reactive products in the searched set (a BLAST screen, not wet-lab proof)"
    else:
        verdict = "no predicted products in the searched set"
    return dict(products=products, modal_size=modal, size_tol=tol, n_products=len(products),
                n_subjects=n_subj, n_off_size=n_off, n_probe_binding=n_probe, verdict=verdict)


def assay_specificity(forward, reverse, probe=None, mode="remote", db="nt", db_path=None,
                      organism=None, min_product=40, max_product=3000):
    """Full-assay specificity: BLAST F + R (+ probe) against NCBI (mode='remote', the default) or a
    local DB, predict amplicons (epcr), and check whether the probe binds inside any of them."""
    fwd, rev = forward.upper().strip(), reverse.upper().strip()
    pr = (probe or "").upper().strip()
    if mode == "local":
        pair, err = _blast_local_pair(fwd, rev, db_path or "")
        phits = []
        if not err and pr:
            pp, perr = _blast_local_pair(pr, pr, db_path or "")
            phits = (pp or {}).get("F", []) if not perr else []
    else:
        seqs = [("F", fwd), ("R", rev)] + ([("P", pr)] if pr else [])
        sets, err = _blast_remote_set(seqs, organism=organism, db=db)
        pair = {"F": (sets or {}).get("F", []), "R": (sets or {}).get("R", [])} if not err else None
        phits = (sets or {}).get("P", []) if (not err and pr) else []
    if err:
        return dict(error=err)
    hits = ([dict(primer="F", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in pair.get("F", [])] +
            [dict(primer="R", **{k: h[k] for k in ("subject", "lo", "hi", "strand", "q3")}) for h in pair.get("R", [])])
    products = epcr(hits, min_product, max_product)
    v = assay_verdict(products, phits)
    v.update(forward_hits=len(pair.get("F", [])), reverse_hits=len(pair.get("R", [])),
             probe_hit_count=len(phits), probe_used=bool(pr), organism=(organism or None))
    v["products"] = v["products"][:25]
    return v
