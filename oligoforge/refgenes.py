"""geNorm-style reference-gene stability screening from a Cq matrix.

This module implements pairwise geNorm M and V(n/n+1), plus a simple per-gene
Cq-SD diagnostic.  It does *not* implement the full BestKeeper or NormFinder
models.  Columns must represent independent biological samples; technical
replicates should be averaged before input.
"""
import math
import statistics as stats


def parse_table(text):
    """Parse ``GENE cq cq ...`` rows into ``{gene: [Cq, ...]}``."""
    data, n = {}, None
    for ln in (text or "").strip().splitlines():
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
        if any(not math.isfinite(c) or c <= 0 or c > 50 for c in cqs):
            raise ValueError("row '%s': Cq values should be finite and in (0, 50]" % gene)
        if n is None:
            n = len(cqs)
        elif len(cqs) != n:
            raise ValueError("row '%s' has %d Cq values; expected %d (one per sample)" %
                             (gene, len(cqs), n))
        if gene in data:
            raise ValueError("duplicate gene '%s'" % gene)
        data[gene] = cqs
    if len(data) < 2:
        raise ValueError("need at least two candidate reference genes")
    return data


def _rel_quant(cqs, efficiency=1.0):
    """Relative quantity using amplification factor ``1 + efficiency``.

    ``efficiency=1.0`` is 100% efficiency (factor 2).  Values may be supplied per
    assay when validated; otherwise the 100% assumption is explicitly reported.
    """
    factor = 1.0 + float(efficiency)
    if not math.isfinite(factor) or factor <= 1.0 or factor > 3.0:
        raise ValueError("efficiency must be a fraction in (0, 2], e.g. 0.95 for 95%")
    m = min(cqs)
    return [factor ** (m - c) for c in cqs]


def _pairwise_V(qj, qk):
    ratios = [math.log2(qj[s] / qk[s]) for s in range(len(qj))]
    return stats.stdev(ratios) if len(ratios) > 1 else 0.0


def _M_for(genes, q):
    out = {}
    for gene in genes:
        values = [_pairwise_V(q[gene], q[other]) for other in genes if other != gene]
        out[gene] = sum(values) / len(values) if values else 0.0
    return out


def analyze(cq_by_gene, efficiencies=None):
    """Run geNorm-style ranking.

    ``efficiencies`` optionally maps gene -> amplification efficiency fraction
    (1.0 = 100%).  Missing entries use 1.0 and are listed in the result.
    """
    genes = list(cq_by_gene)
    if len(genes) < 2:
        return dict(error="need at least two candidate reference genes")
    n_samples = len(next(iter(cq_by_gene.values())))
    if n_samples < 2:
        return dict(error="need at least two biological samples")
    if any(len(cq_by_gene[g]) != n_samples for g in genes):
        return dict(error="every gene must have one Cq value for each biological sample")

    efficiencies = dict(efficiencies or {})
    assumed = []
    q = {}
    try:
        for gene in genes:
            if gene not in efficiencies:
                assumed.append(gene)
            q[gene] = _rel_quant(cq_by_gene[gene], efficiencies.get(gene, 1.0))
    except (TypeError, ValueError) as exc:
        return dict(error=str(exc))

    remaining = list(genes)
    elimination = []
    while len(remaining) > 2:
        M = _M_for(remaining, q)
        worst = max(remaining, key=lambda gene: (M[gene], gene))
        elimination.append((worst, M[worst]))
        remaining.remove(worst)

    # geNorm cannot distinguish the final two genes: their pairwise variation and
    # therefore M are identical.  Preserve that tie rather than inventing a winner.
    final_m = _M_for(remaining, q)
    final_pair = sorted(remaining)
    ranking_names = final_pair + [g for g, _ in reversed(elimination)]
    Mby = {g: final_m[g] for g in final_pair}
    Mby.update({g: m for g, m in elimination})

    def normalization_factor(k):
        chosen = ranking_names[:k]
        return [math.exp(sum(math.log(q[g][s]) for g in chosen) / k)
                for s in range(n_samples)]

    pairwise = []
    for k in range(2, len(ranking_names)):
        a, b = normalization_factor(k), normalization_factor(k + 1)
        ratios = [math.log2(a[s] / b[s]) for s in range(n_samples)]
        pairwise.append(dict(step="V%d/%d" % (k, k + 1),
                             n=k,
                             v=round(stats.stdev(ratios), 4)))

    crossing = next((item for item in pairwise if item["v"] < 0.15), None)
    recommended_n = crossing["n"] if crossing else (2 if len(genes) == 2 else None)
    if recommended_n is None:
        recommendation = ("No V(n/n+1) value crossed the commonly used 0.15 guideline. "
                          "Do not automatically use every candidate; review biological context, "
                          "sample size, efficiency and an independent stability method.")
    else:
        recommendation = ("The first V(n/n+1) below 0.15 occurs at n=%d. Treat 0.15 as a "
                          "context-dependent guideline, not a universal pass/fail threshold." %
                          recommended_n)

    per_gene = []
    for gene in ranking_names:
        cqs = cq_by_gene[gene]
        mean = stats.mean(cqs)
        sd = stats.stdev(cqs)
        per_gene.append(dict(
            gene=gene,
            M=round(Mby[gene], 4),
            tied_final_pair=gene in final_pair,
            cq_mean=round(mean, 2), cq_sd=round(sd, 3),
            cq_cv=round(100.0 * sd / mean, 2) if mean else 0.0,
            cq_sd_screen=("review (SD > 1 Cq)" if sd > 1.0 else "within 1 Cq"),
            # Retain the old key for frontend/backward compatibility, but label it
            # as a screen rather than claiming a full BestKeeper calculation.
            bestkeeper=("unstable (SD>1 Cq)" if sd > 1.0 else "ok"),
            efficiency=efficiencies.get(gene, 1.0),
            efficiency_assumed=gene in assumed,
        ))

    return dict(
        n_genes=len(genes), n_samples=n_samples,
        ranking=per_gene, final_pair=final_pair,
        pairwise_V=pairwise, recommended_n=recommended_n,
        recommendation=recommendation,
        assumed_efficiency_genes=assumed,
        note=("M and V are geNorm-style pairwise stability statistics. The final two genes are tied by "
              "definition. Cq SD is only a simple dispersion screen, not full BestKeeper. Columns must "
              "be independent biological samples; average technical replicates first. Missing assay "
              "efficiencies were assumed to be 100%. Cross-check with biological knowledge and an "
              "independent method such as NormFinder when making publication claims."),
    )
