"""Certified Orthogonal Panel — pick a maximum mutually non-cross-hybridizing oligo set from a
candidate pool under a nearest-neighbor ΔG confusability model, and CERTIFY how close the chosen
panel is to the true maximum *under that model*.

Two layers, kept separate so the graph/certificate math is testable with zero thermodynamics:

  math core (no `thermo`):  max_independent_set, greedy_clique_cover, clique_cover_bound,
                            lovasz_theta (optional cvxpy), integer_upper_bound,
                            strong_product_capacity
  thermo driver (uses T):   intake, self_structure_filter, build_confusability_graph, certify_panel

WHAT THE CERTIFICATE MEANS
θ(G) and the clique-cover number are UPPER bounds on the independence number α(G); a maximum
independent set is a LOWER bound. When lower == upper, the panel is the provable maximum *of the
ΔG-threshold graph*. That graph is a pairwise proxy for a non-pairwise phenomenon (real multiplex
cross-priming depends on concentrations, competing templates, polymerase, cycling), so a certified
maximum is optimal for the model, not guaranteed optimal at the bench. That distinction is reported,
not hidden. See SPEC.md.

CERTIFICATE STRATEGY (cheap first)
The sandwich theorem gives a FREE upper bound — α(G) ≤ number of cliques in a clique cover of G —
computable in pure Python. On curated panels it often already equals the exact MIS (gap 0, certified,
no SDP). The Lovász-θ SDP (tighter) is attempted only when the free bound doesn't close the gap, and
only if cvxpy is installed; otherwise the module still certifies via the clique-cover bound. cvxpy is
never a hard requirement.
"""
import math

# ------------------------------------------------------------------------------------------------
# MATH CORE — graphs as (n, edges) where edges is a list of (i, j) with i < j. No thermodynamics.
# ------------------------------------------------------------------------------------------------

def _adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for e in edges:
        i, j = e[0], e[1]
        if i != j:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def _greedy_mis(n, adj):
    """Min-degree greedy independent set. Returns a valid independent set (a genuine lower bound)."""
    remaining = set(range(n))
    chosen = []
    while remaining:
        v = min(remaining, key=lambda x: len(adj[x] & remaining))
        chosen.append(v)
        remaining.discard(v)
        remaining -= adj[v]
    return sorted(chosen)


def max_independent_set(n, edges, exact_limit=150, budget=1_500_000):
    """Maximum independent set of G.

    Exact branch-and-bound for n <= exact_limit (with an expansion budget that guarantees
    termination); greedy min-degree otherwise, or if the exact search exhausts its budget on a
    pathologically dense graph. Returns (vertex_list, method, exact: bool). A greedy result is a
    valid independent set and thus a valid LOWER bound; `exact` is False so callers don't overclaim.
    """
    adj = _adjacency(n, edges)
    if n == 0:
        return [], "trivial", True
    if n > exact_limit:
        return _greedy_mis(n, adj), "greedy_min_degree", False

    best = []
    calls = [0]
    aborted = [False]

    def expand(P, cur):
        if aborted[0]:
            return
        calls[0] += 1
        if calls[0] > budget:
            aborted[0] = True
            return
        if not P:
            if len(cur) > len(best[0] if best else []):
                best[:] = [list(cur)]
            return
        if len(cur) + len(P) <= (len(best[0]) if best else 0):
            return
        # branch on a min-degree-in-P vertex; include-first tends to reach good solutions fast
        v = min(P, key=lambda u: len(adj[u] & P))
        nbrs_in_P = adj[v] & P
        # include v
        expand(P - nbrs_in_P - {v}, cur + [v])
        # exclude v — only meaningful if v has neighbors in P (an isolated vertex is always kept)
        if nbrs_in_P:
            expand(P - {v}, cur)

    best = [[]]
    expand(set(range(n)), [])
    if aborted[0]:
        return _greedy_mis(n, adj), "greedy_after_budget", False
    return sorted(best[0]), "exact_bnb", True


def greedy_clique_cover(n, edges):
    """Partition the vertices of G into cliques of G (greedy, seeding from high-degree vertices and
    extending maximally). Returns a list of cliques (each a sorted vertex list). The number of
    cliques is an UPPER bound on α(G): each clique holds at most one independent-set vertex
    (sandwich theorem, α(G) ≤ clique-cover number). Pure Python, no dependencies."""
    adj = _adjacency(n, edges)
    remaining = set(range(n))
    cover = []
    while remaining:
        seed = max(remaining, key=lambda v: len(adj[v] & remaining))
        clique = [seed]
        cand = adj[seed] & remaining
        while cand:
            u = max(cand, key=lambda v: len(adj[v] & cand))
            clique.append(u)
            cand = cand & adj[u]
        cover.append(sorted(clique))
        remaining -= set(clique)
    return cover


def clique_cover_bound(n, edges):
    """Free (pure-Python) upper bound on α(G): the number of cliques in a greedy clique cover."""
    if n == 0:
        return 0
    return len(greedy_clique_cover(n, edges))


def lovasz_theta(n, edges):
    """Lovász θ(G) via the primal SDP  θ = max ⟨J,X⟩  s.t. X⪰0, tr(X)=1, X_ij=0 ∀{i,j}∈E.
    An UPPER bound on α(G), tighter than the clique-cover bound. Requires cvxpy (+ a conic solver);
    returns None if cvxpy is not installed or the solve fails, so the caller degrades to the
    clique-cover certificate. ⟨J,X⟩ is the sum of all entries of X."""
    try:
        import cvxpy as cp  # optional dependency — imported lazily on purpose
    except Exception:
        return None
    if n == 0:
        return 0.0
    if n == 1:
        return 1.0
    try:
        X = cp.Variable((n, n), symmetric=True)
        cons = [X >> 0, cp.trace(X) == 1]
        for e in edges:
            i, j = e[0], e[1]
            if i != j:
                cons.append(X[i, j] == 0)
        prob = cp.Problem(cp.Maximize(cp.sum(X)), cons)
        try:
            prob.solve(solver=cp.SCS, verbose=False)
        except Exception:
            prob.solve(verbose=False)
        if X.value is None:
            return None
        val = float(X.value.sum())
        # θ ∈ [1, n]; clamp tiny solver excursions outside the valid range
        return min(max(val, 1.0), float(n))
    except Exception:
        return None


def integer_upper_bound(theta, snap_tol=2e-2, eps=1e-6):
    """Integer upper bound on α(G) from a (noisy) θ. CRITICAL: do not naive-floor — a solver can
    return 3.9999 for a true 4, and floor()=3 would understate the bound and could violate the
    α ≤ ⌊θ⌋ invariant. If θ is within snap_tol of an integer, snap to it; otherwise floor(θ+eps).
    Snapping never understates: θ ≥ α, and round(θ) for θ within snap_tol (<0.5) of an integer is
    ≥ α."""
    if theta is None:
        return None
    r = round(theta)
    if abs(theta - r) <= snap_tol:
        return int(r)
    return int(math.floor(theta + eps))


def strong_product_capacity(value, k):
    """θ(G⊠H)=θ(G)·θ(H); for k identical rounds, θ(G)^k — the certified collision-free k-round
    combinatorial-barcode capacity UNDER the strong-product confusability model (two k-round
    barcodes collide iff they cross-hybridize in ≥1 round). Also used with |MIS| to compute the
    naive product-of-independent-sets count for comparison. Returns None for value None."""
    if value is None:
        return None
    return float(value) ** k


# ------------------------------------------------------------------------------------------------
# THERMODYNAMIC DRIVER — turns oligos into a confusability graph and runs the certificate.
# ------------------------------------------------------------------------------------------------

from . import thermo as T


def intake(candidates):
    """Normalize a candidate pool. Accepts strings or dicts {seq|sequence, role?, target?, name?}.
    Strips IDT modification/LNA notation, cleans (IUPAC preserved, RNA->DNA, FASTA header dropped),
    dedups identical sequences (recording dup_count), and flags degenerate/RNA. Returns
    (records, rejects) where each record is {seq, role, target, name, degenerate, dup_count} and
    each reject is {input, name, reason}."""
    seen = {}
    records = []
    rejects = []
    for idx, c in enumerate(candidates or []):
        if isinstance(c, str):
            raw, role, target, name = c, None, None, None
        elif isinstance(c, dict):
            raw = c.get("seq") or c.get("sequence")
            role = c.get("role")
            target = c.get("target")
            name = c.get("name")
        else:
            continue
        nm = name or ("oligo%d" % (idx + 1))
        if not raw or not str(raw).strip():
            continue
        bare = T.strip_mods(str(raw))
        s, notes, err = T.clean_seq(bare)
        if err:
            rejects.append({"input": str(raw)[:60], "name": nm, "reason": err})
            continue
        if s in seen:
            seen[s]["dup_count"] += 1
            if nm not in seen[s]["names"]:
                seen[s]["names"].append(nm)
            continue
        rec = {"seq": s, "role": role, "target": target, "name": nm,
               "degenerate": T.has_degenerate(s), "dup_count": 1, "names": [nm],
               "rna_converted": any("RNA" in x for x in notes)}
        seen[s] = rec
        records.append(rec)
    return records, rejects


def self_structure_filter(records, self_dg=-9.0):
    """Drop candidates whose worst self-structure (hairpin or homodimer) ΔG is below (more stable
    than) self_dg, under the current thermo conditions. Records why each was dropped. Returns
    (survivors, dropped). Each survivor gains hairpin_dg / homodimer_dg."""
    survivors = []
    dropped = []
    for r in records:
        s = r["seq"]
        hp = T.hairpin(s)[0]          # hairpin() returns (dG, Tm); we gate on dG
        hd = T.self_dimer(s)
        worst = min(hp, hd)
        if worst < self_dg:
            which = "hairpin" if hp <= hd else "homodimer"
            dropped.append({**r, "hairpin_dg": round(hp, 2), "homodimer_dg": round(hd, 2),
                            "reason": "%s ΔG %.1f < %.1f kcal/mol (too stable a self-structure)"
                                      % (which, worst, self_dg)})
        else:
            survivors.append({**r, "hairpin_dg": round(hp, 2), "homodimer_dg": round(hd, 2)})
    return survivors, dropped


def build_confusability_graph(seqs, cross_dg=-6.0):
    """Confusability graph over surviving candidate sequences. Edge {i,j} iff the most stable
    heterodimer between i and j — the more negative of ΔG(i,j) and ΔG(i, revcomp(j)), covering both
    the direct duplex and i binding j's complement — is below (more stable than) cross_dg. Returns
    (edges, edge_dg) with edge_dg[(i,j)] = the ΔG that formed the edge. O(n²) NN calls; serial (the
    primer3 engine is lock-guarded, so threads would not parallelize it anyway) — fine at curated
    scale."""
    n = len(seqs)
    rc = [T.revcomp(s) for s in seqs]
    edges = []
    edge_dg = {}
    for i in range(n):
        for j in range(i + 1, n):
            dg = min(T.hetero_dimer(seqs[i], seqs[j]), T.hetero_dimer(seqs[i], rc[j]))
            if dg < cross_dg:
                edges.append((i, j))
                edge_dg[(i, j)] = round(dg, 2)
    return edges, edge_dg


def certify_panel(candidates, cross_dg=-6.0, self_dg=-9.0, k=1,
                  size_limit=600, use_theta=True):
    """Full pipeline: intake -> self-structure filter -> confusability graph -> maximum independent
    set (lower bound) -> clique-cover bound (free upper) -> optional θ (tighter upper) -> gap /
    certified flag -> split-pool θ^k ceiling. Pure w.r.t. thermo conditions except that callers set
    conditions before calling. Returns a JSON-able dict."""
    records, rejects = intake(candidates)
    survivors, dropped = self_structure_filter(records, self_dg=self_dg)
    seqs = [r["seq"] for r in survivors]
    n = len(seqs)

    result = {
        "params": {"cross_dg": cross_dg, "self_dg": self_dg, "k": k, "size_limit": size_limit},
        "n_input": len(candidates or []),
        "n_unique": len(records),
        "n_after_self_filter": n,
        "rejects": rejects,
        "dropped_self_structure": [{"name": d["name"], "seq": d["seq"], "reason": d["reason"]}
                                   for d in dropped],
        "duplicates": [{"seq": r["seq"], "names": r["names"], "count": r["dup_count"]}
                       for r in records if r["dup_count"] > 1],
    }

    if n == 0:
        result.update({"panel": [], "panel_size": 0, "edges": [], "theta": None,
                       "upper_bound": 0, "gap": 0, "certified": True,
                       "note": "no candidates survived the self-structure filter"})
        return result
    if n == 1:
        result.update({"panel": _panel_out(survivors, [0]), "panel_size": 1, "edges": [],
                       "theta": None, "upper_bound": 1, "gap": 0, "certified": True,
                       "bound_source": "trivial",
                       "split_pool": _split_pool(1, None, k)})
        return result

    edges, edge_dg = build_confusability_graph(seqs, cross_dg=cross_dg)

    mis, mis_method, mis_exact = max_independent_set(n, edges)
    lower = len(mis)

    # free upper bound first (pure Python, always available)
    cc_bound = clique_cover_bound(n, edges)
    upper = cc_bound
    bound_source = "clique_cover"
    theta_raw = None
    theta_int = None

    # tighten with θ only if the free bound hasn't already closed the gap, θ is wanted,
    # cvxpy is available, and the graph is within the SDP size guard
    if use_theta and (cc_bound - lower) > 0 and n <= size_limit:
        theta_raw = lovasz_theta(n, edges)
        theta_int = integer_upper_bound(theta_raw)
        if theta_int is not None and theta_int < upper:
            upper = theta_int
            bound_source = "lovasz_theta"

    theta_unavailable = use_theta and (cc_bound - lower) > 0 and (
        n > size_limit or (theta_raw is None))

    gap = upper - lower
    certified = (gap == 0)  # a valid MIS whose size equals a valid upper bound IS maximum

    result.update({
        "panel": _panel_out(survivors, mis),
        "panel_size": lower,
        "mis_method": mis_method,
        "mis_exact": mis_exact,
        "edges": [{"i": i, "j": j, "dg": edge_dg[(i, j)],
                   "a": survivors[i]["name"], "b": survivors[j]["name"]} for (i, j) in edges],
        "n_edges": len(edges),
        "clique_cover_bound": cc_bound,
        "theta": (round(theta_raw, 4) if theta_raw is not None else None),
        "theta_int": theta_int,
        "upper_bound": upper,
        "bound_source": bound_source,
        "theta_unavailable": theta_unavailable,
        "gap": gap,
        "certified": certified,
        "split_pool": _split_pool(lower, (theta_raw if bound_source == "lovasz_theta" else None), k),
        "nodes": [{"id": i, "name": survivors[i]["name"], "seq": survivors[i]["seq"],
                   "role": survivors[i].get("role"), "target": survivors[i].get("target"),
                   "in_panel": i in set(mis)} for i in range(n)],
    })
    if not certified:
        result["note"] = ("gap %d: the panel is at least %d and at most %d oligos; not proven "
                          "maximum under this model" % (gap, lower, upper))
    return result


def _panel_out(survivors, idxs):
    out = []
    for i in idxs:
        r = survivors[i]
        out.append({"name": r["name"], "seq": r["seq"], "role": r.get("role"),
                    "target": r.get("target"), "hairpin_dg": r.get("hairpin_dg"),
                    "homodimer_dg": r.get("homodimer_dg")})
    return out


def _split_pool(mis_size, theta_raw, k):
    """Split-pool capacity: certified θ^k ceiling (when θ is available) vs the naive
    product-of-independent-sets count |MIS|^k, to show where distance reasoning overcounts."""
    naive = strong_product_capacity(mis_size, k)
    certified_ceiling = strong_product_capacity(theta_raw, k) if theta_raw is not None else None
    out = {"k": k, "naive_mis_pow_k": (round(naive, 4) if naive is not None else None)}
    if certified_ceiling is not None:
        out["theta_pow_k"] = round(certified_ceiling, 4)
        out["certified_max_barcodes"] = int(math.floor(certified_ceiling + 1e-6))
        out["assumption"] = ("θ^k is the certified collision-free k-round barcode capacity under the "
                             "strong-product model: two barcodes collide iff they cross-hybridize in "
                             "at least one round. A modeling assumption, not a wet-lab guarantee.")
    else:
        out["note"] = ("θ unavailable (cvxpy not installed or graph too large); split-pool ceiling "
                       "not certified. Naive |MIS|^k shown for reference only.")
    return out


# ------------------------------------------------------------------------------------------------
# DEMO — `python -m oligoforge.orthopanel` prints a certificate on a small sample panel.
# ------------------------------------------------------------------------------------------------

def _demo():
    print("=== Certified Orthogonal Panel — demo ===\n")

    # Part 1 — a real oligo pool. Five candidates, all passing the self-structure filter, where one
    # is the reverse complement of another (they cross-hybridize and cannot both be in the panel).
    A = "ACGATCAGTTGCATCAGGTA"
    pool = [
        {"seq": A, "name": "FwdA"},
        {"seq": T.revcomp(A), "name": "ProbeA_rc"},   # revcomp of FwdA -> edge
        {"seq": "GCATTCAGATCGGATACCTA", "name": "FwdB"},
        {"seq": "TGACCTAGCATGGTACAAGC", "name": "FwdC"},
        {"seq": "CATGACTGCAAGTCGATACG", "name": "FwdD"},
    ]
    res = certify_panel(pool, cross_dg=-6.0, self_dg=-9.0, k=2)
    print("[1] real oligo pool")
    print("    %d unique, %d after self-filter, %d cross-reaction edge(s)"
          % (res["n_unique"], res["n_after_self_filter"], res.get("n_edges", 0)))
    for e in res.get("edges", []):
        print("    edge %s–%s  ΔG %.1f kcal/mol" % (e["a"], e["b"], e["dg"]))
    print("    panel (size %d): %s" % (res["panel_size"], ", ".join(p["name"] for p in res["panel"])))
    verdict = "CERTIFIED MAXIMUM" if res["certified"] else ("gap %d" % res["gap"])
    print("    certificate: %d ≤ α ≤ %d  (bound: %s)  ->  %s"
          % (res["panel_size"], res["upper_bound"], res.get("bound_source"), verdict))
    print("    (θ not needed here — the free clique-cover bound already closed the gap)\n")

    # Part 2 — a confusability graph that is a 5-cycle (each oligo cross-reacts with two neighbours).
    # This is the case the free bound CANNOT certify: clique-cover of C5 = 3, but α(C5) = 2. Only the
    # Lovász-θ SDP (θ = √5 ≈ 2.236, ⌊θ⌋ = 2) closes the gap and proves the panel of 2 is maximum.
    C5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    mis, _, _ = max_independent_set(5, C5)
    cc = clique_cover_bound(5, C5)
    th = lovasz_theta(5, C5)
    print("[2] a 5-way cross-reaction cycle (C5) — where θ earns its keep")
    print("    max independent set (panel): %d oligos" % len(mis))
    print("    free clique-cover bound: %d  (too loose — cannot certify)" % cc)
    if th is not None:
        print("    Lovász θ: %.4f  ->  ⌊θ⌋ = %d  (CERTIFIES the panel of %d is maximum)"
              % (th, integer_upper_bound(th), len(mis)))
        for kk in (2, 3):
            cap = strong_product_capacity(th, kk)
            print("    split-pool k=%d: certified ≤ %d collision-free barcodes (θ^%d=%.3f) vs "
                  "naive |MIS|^%d=%d" % (kk, math.floor(cap + 1e-6), kk, cap, kk, len(mis) ** kk))
    else:
        print("    Lovász θ: unavailable (cvxpy not installed) — install cvxpy for the tight "
              "certificate; the free bound above still applies")


if __name__ == "__main__":
    _demo()
