"""Reference-gene stability from a Cq matrix.

geNorm (Vandesompele 2002): gene-stability measure M (mean pairwise variation of a
gene against all others; lower = more stable) and pairwise variation V(n/n+1) for
choosing how many reference genes to use (threshold 0.15). Plus per-gene Cq SD and
CV (BestKeeper-style). The model-based NormFinder estimate is intentionally not
reimplemented; cross-check the ranking in NormFinder / RefFinder. Relative
quantities are derived from Cq assuming efficiency E = 2:  Q = 2^(minCq - Cq).
"""
import math
import statistics as stats


def parse_table(text):
    """Parse 'GENE cq cq cq ...' lines (space or comma separated) -> {gene: [cq,...]}.
    Raises ValueError on malformed input."""
    data, n = {}, None
    for ln in text.strip().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.replace(",", " ").split()
        if len(parts) < 3:
            raise ValueError("row '%s' needs a gene name and at least 2 Cq values" % ln[:30])
        gene = parts[0]
        try:
            cqs = [float(x) for x in parts[1:]]
        except ValueError:
            raise ValueError("non-numeric Cq in row '%s'" % gene)
        if any(c <= 0 or c > 50 for c in cqs):
            raise ValueError("row '%s': Cq values should be in (0, 50]" % gene)
        if n is None:
            n = len(cqs)
        elif len(cqs) != n:
            raise ValueError("row '%s' has %d Cq values; expected %d (one per sample)" % (gene, len(cqs), n))
        if gene in data:
            raise ValueError("duplicate gene '%s'" % gene)
        data[gene] = cqs
    if not data:
        raise ValueError("no data rows found")
    return data


def _rel_quant(cqs):
    m = min(cqs)
    return [2.0 ** (m - c) for c in cqs]


def _pairwise_V(qj, qk):
    A = [math.log2(qj[s] / qk[s]) for s in range(len(qj))]
    return stats.stdev(A) if len(A) > 1 else 0.0


def _M_for(genes, q):
    M = {}
    for j in genes:
        vs = [_pairwise_V(q[j], q[k]) for k in genes if k != j]
        M[j] = sum(vs) / len(vs) if vs else 0.0
    return M


def analyze(cq_by_gene):
    genes = list(cq_by_gene)
    q = {g: _rel_quant(cq_by_gene[g]) for g in genes}
    n_samples = len(next(iter(cq_by_gene.values())))

    remaining = list(genes)
    elimination = []                       # least-stable first
    while len(remaining) > 2:
        M = _M_for(remaining, q)
        worst = max(M, key=lambda g: (M[g], g))
        elimination.append((worst, round(M[worst], 4)))
        remaining.remove(worst)
    finalM = _M_for(remaining, q)
    for g in sorted(remaining, key=lambda g: (finalM[g], g), reverse=True):
        elimination.append((g, round(finalM[g], 4)))
    ranking = list(reversed(elimination))  # most stable first
    order = [g for g, _ in ranking]
    Mby = {g: m for g, m in ranking}

    def NF(k):
        chosen = order[:k]
        return [math.exp(sum(math.log(q[g][s]) for g in chosen) / k) for s in range(n_samples)]

    V = []
    for k in range(2, len(order)):
        a, b = NF(k), NF(k + 1)
        ratios = [math.log2(a[s] / b[s]) for s in range(n_samples)]
        V.append(dict(step="V%d/%d" % (k, k + 1),
                      v=round(stats.stdev(ratios) if n_samples > 1 else 0.0, 4)))

    rec = next((int(item["step"][1:].split("/")[0]) for item in V if item["v"] < 0.15), len(order))

    per_gene = []
    for g in order:
        cqs = cq_by_gene[g]
        mean = stats.mean(cqs)
        sd = stats.stdev(cqs) if n_samples > 1 else 0.0
        per_gene.append(dict(gene=g, M=Mby[g], cq_mean=round(mean, 2), cq_sd=round(sd, 3),
                             cq_cv=round(100 * sd / mean, 2) if mean else 0.0,
                             bestkeeper="unstable (SD>1 Cq)" if sd > 1.0 else "ok"))

    return dict(n_genes=len(genes), n_samples=n_samples, ranking=per_gene, pairwise_V=V,
                recommended_n=rec,
                note=("M = geNorm gene-stability (lower is more stable; <1.5 is the usual cutoff). "
                      "V(n/n+1) < 0.15 means n reference genes are enough. BestKeeper flag uses Cq SD > 1. "
                      "Cross-check the ranking in NormFinder / RefFinder."))
