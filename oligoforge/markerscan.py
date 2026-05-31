"""Optional data-driven check on the gene-finder's curated markers.

For each candidate marker, report what GenBank actually holds for the taxon:
  count                          nucleotide records matching the query (availability)
  n                              records sampled for length stats
  median_len / min_len / max_len typical record size (bp)

That is deliberately all it measures. A cheap cross-record "conservation" or
target-vs-relative "separation" score was prototyped (k-mer containment and local-
alignment identity, with length-banding) and dropped: GenBank records for a single marker
cover different sub-regions and lengths, so any quick cross-record similarity conflates
"different region sampled" with "divergent" and misranks markers -- e.g. 18S came out as
both the best and the worst genus separator depending on the method, neither true, while
cytb compared honestly only because the whole field sequences the same MalAvi region.
Availability and length are the signals that hold up. Rigorous per-position conservation
and discrimination is the Conservation tab (it aligns a chosen region); the gene-finder
routes there via the marker chip rather than fabricating a number here.
"""
from . import ncbi


def scan(base, markers, exclude=None, k_seqs=6, max_markers=6):
    out = []
    for m in markers[:max_markers]:
        q = ("%s %s" % (base, m["qwords"])).strip()
        count = ncbi.count_hits(q)
        try:
            seqs = [s for _, s in ncbi.search_fetch_fasta(q, k_seqs)]
        except Exception:
            seqs = []
        lens = sorted(len(s) for s in seqs)
        out.append(dict(
            gene=m["gene"], query=q, count=count, n=len(seqs),
            median_len=(lens[len(lens) // 2] if lens else None),
            min_len=(lens[0] if lens else None),
            max_len=(lens[-1] if lens else None),
        ))
    return dict(base=base, exclude=(exclude or "").strip(), results=out)
