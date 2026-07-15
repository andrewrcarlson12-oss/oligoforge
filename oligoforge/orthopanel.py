"""Orthogonal-panel selection under a pairwise thermodynamic confusability model.

The thermodynamic layer builds a graph whose edges represent candidate cross-hybridization. A
valid panel is an independent set. Exact branch-and-bound results are reported as proven maxima;
otherwise a greedy independent set is a constructive lower bound and a clique cover supplies a
rigorous upper bound. An optional Lovász-theta SDP is retained only as a numerical diagnostic: an
approximate solver value is never rounded into a formal certificate.

All guarantees are guarantees about the thresholded pairwise model, not wet-lab multiplex behavior.
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
    """Numerically estimate Lovász θ(G) from the primal SDP.

    The mathematical optimum is an upper bound on α(G), but ordinary floating-point conic-solver
    output is treated here as a diagnostic only. Formal model bounds come from completed exact
    branch-and-bound or a valid clique cover. Returns None when cvxpy/solver support is unavailable.
    """
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
    """Return a display-only integerization of a numerical theta estimate.

    This helper is retained for backwards compatibility and tests. The returned number is NOT used
    as a rigorous upper bound or certificate because ordinary SDP solver output is approximate. If
    theta is within snap_tol of an integer, snap to it; otherwise floor(θ+eps).
    Because the input is a numerical solver estimate, this display value must not be used as a
    formal graph bound."""
    if theta is None:
        return None
    r = round(theta)
    if abs(theta - r) <= snap_tol:
        return int(r)
    return int(math.floor(theta + eps))


def strong_product_capacity(value, k):
    """Raise a graph quantity to k for split-pool diagnostics.

    With an exact mathematical theta value, multiplicativity applies to the strong graph product;
    this implementation receives approximate numerical values and therefore reports theta^k only
    as a diagnostic. The same helper computes the constructive |panel|^k count. Returns None for
    value None.
    """
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


def self_structure_filter(records, self_dg=-9.0, thermo_snapshot=None):
    """Drop candidates whose worst self-structure (hairpin or homodimer) ΔG is below (more stable
    than) self_dg, under the current thermo conditions. Records why each was dropped. Returns
    (survivors, dropped). Each survivor gains hairpin_dg / homodimer_dg."""
    snap = thermo_snapshot or T._snapshot()
    survivors = []
    dropped = []
    for r in records:
        s = r["seq"]
        hp = T._hairpin_at(s, snap)[0]
        hd = T._self_dimer_at(s, snap)
        worst = min(hp, hd)
        if worst < self_dg:
            which = "hairpin" if hp <= hd else "homodimer"
            dropped.append({**r, "hairpin_dg": round(hp, 2), "homodimer_dg": round(hd, 2),
                            "reason": "%s ΔG %.1f < %.1f kcal/mol (too stable a self-structure)"
                                      % (which, worst, self_dg)})
        else:
            survivors.append({**r, "hairpin_dg": round(hp, 2), "homodimer_dg": round(hd, 2)})
    return survivors, dropped


def build_confusability_graph(seqs, cross_dg=-6.0, thermo_snapshot=None):
    """Confusability graph over surviving candidate sequences. Edge {i,j} iff the most stable
    heterodimer between i and j — the more negative of ΔG(i,j) and ΔG(i, revcomp(j)), covering both
    the direct duplex and i binding j's complement — is below (more stable than) cross_dg. Returns
    (edges, edge_dg) with edge_dg[(i,j)] = the ΔG that formed the edge. O(n²) NN calls; serial (the
    primer3 engine is lock-guarded, so threads would not parallelize it anyway) — fine at curated
    scale."""
    snap = thermo_snapshot or T._snapshot()
    n = len(seqs)
    rc = [T.revcomp(s) for s in seqs]
    edges = []
    edge_dg = {}
    for i in range(n):
        for j in range(i + 1, n):
            dg = min(T._hetero_dimer_at(seqs[i], seqs[j], snap),
                     T._hetero_dimer_at(seqs[i], rc[j], snap))
            if dg < cross_dg:
                edges.append((i, j))
                edge_dg[(i, j)] = round(dg, 2)
    return edges, edge_dg


def certify_panel(candidates, cross_dg=-6.0, self_dg=-9.0, k=1,
                  size_limit=600, use_theta=True, thermo_snapshot=None):
    """Select a mutually non-confusable panel and report rigorous model bounds.

    Exact branch-and-bound proves optimality when it completes. Otherwise the selected independent
    set is a lower bound and a clique cover is a rigorous upper bound. Lovász theta is calculated,
    when requested and feasible, as a diagnostic only and never changes the certificate.
    """
    snap = thermo_snapshot or T._snapshot()
    records, rejects = intake(candidates)
    survivors, dropped = self_structure_filter(records, self_dg=self_dg, thermo_snapshot=snap)
    seqs = [r["seq"] for r in survivors]
    n = len(seqs)

    result = {
        "params": {"cross_dg": cross_dg, "self_dg": self_dg, "k": k, "size_limit": size_limit},
        "n_input": len(candidates or []), "n_unique": len(records), "n_after_self_filter": n,
        "rejects": rejects,
        "dropped_self_structure": [{"name": d["name"], "seq": d["seq"], "reason": d["reason"]}
                                   for d in dropped],
        "duplicates": [{"seq": r["seq"], "names": r["names"], "count": r["dup_count"]}
                       for r in records if r["dup_count"] > 1],
    }
    if n == 0:
        result.update({"panel": [], "panel_size": 0, "edges": [], "theta": None,
                       "upper_bound": 0, "gap": 0, "certified": True,
                       "bound_source": "trivial", "note": "no candidates survived the self-structure filter",
                       "split_pool": _split_pool(0, None, k)})
        return result
    if n == 1:
        result.update({"panel": _panel_out(survivors, [0]), "panel_size": 1, "edges": [],
                       "theta": None, "upper_bound": 1, "gap": 0, "certified": True,
                       "bound_source": "trivial", "split_pool": _split_pool(1, None, k)})
        return result

    edges, edge_dg = build_confusability_graph(seqs, cross_dg=cross_dg, thermo_snapshot=snap)
    mis, mis_method, mis_exact = max_independent_set(n, edges)
    lower = len(mis)
    cc_bound = clique_cover_bound(n, edges)

    # A completed exact search is its own proof. Otherwise only the clique cover is used as a
    # rigorous upper bound; approximate SDP values remain diagnostics.
    if mis_exact:
        upper, bound_source, certified = lower, "exact_bnb", True
    else:
        upper, bound_source = cc_bound, "clique_cover"
        certified = (upper == lower)
    gap = upper - lower

    theta_raw = None
    theta_display_int = None
    if use_theta and n <= size_limit:
        theta_raw = lovasz_theta(n, edges)
        theta_display_int = integer_upper_bound(theta_raw)
    theta_unavailable = bool(use_theta and (n > size_limit or theta_raw is None))

    result.update({
        "panel": _panel_out(survivors, mis), "panel_size": lower,
        "mis_method": mis_method, "mis_exact": mis_exact,
        "edges": [{"i": i, "j": j, "dg": edge_dg[(i, j)],
                   "a": survivors[i]["name"], "b": survivors[j]["name"]} for (i, j) in edges],
        "n_edges": len(edges), "clique_cover_bound": cc_bound,
        "theta": (round(theta_raw, 4) if theta_raw is not None else None),
        "theta_display_int": theta_display_int,
        "theta_certifying": False,
        "upper_bound": upper, "bound_source": bound_source,
        "theta_unavailable": theta_unavailable, "gap": gap, "certified": certified,
        "split_pool": _split_pool(lower, theta_raw, k),
        "nodes": [{"id": i, "name": survivors[i]["name"], "seq": survivors[i]["seq"],
                   "role": survivors[i].get("role"), "target": survivors[i].get("target"),
                   "in_panel": i in set(mis)} for i in range(n)],
    })
    if not certified:
        result["note"] = ("model optimum not proven: constructed panel has %d oligos and the "
                          "rigorous clique-cover upper bound is %d" % (lower, upper))
    elif bound_source == "exact_bnb":
        result["note"] = "exact branch-and-bound proved this panel is maximum for the graph model"
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
    """Report constructive and diagnostic k-round quantities without claiming achieved capacity."""
    naive = strong_product_capacity(mis_size, k)
    out = {"k": k, "constructed_panel_pow_k": (round(naive, 4) if naive is not None else None),
           "naive_mis_pow_k": (round(naive, 4) if naive is not None else None)}
    if theta_raw is not None:
        theta_power = strong_product_capacity(theta_raw, k)
        out["theta_pow_k"] = round(theta_power, 4)
        out["theta_diagnostic_only"] = True
        out["note"] = ("theta^k is shown as a numerical graph diagnostic. It is not an achieved "
                       "barcode count and is not used as a formal certificate in this implementation.")
    else:
        out["note"] = ("No theta diagnostic was computed. The constructed |panel|^k count assumes "
                       "independent use of the selected panel in each round; it is not a wet-lab capacity claim.")
    return out


# ------------------------------------------------------------------------------------------------
# DEMO — `python -m oligoforge.orthopanel` prints a certificate on a small sample panel.
# ------------------------------------------------------------------------------------------------

def _demo():
    print("=== Orthogonal Panel — model audit demo ===\n")

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
    verdict = "MODEL MAXIMUM PROVEN" if res["certified"] else ("gap %d" % res["gap"])
    print("    certificate: %d ≤ α ≤ %d  (bound: %s)  ->  %s"
          % (res["panel_size"], res["upper_bound"], res.get("bound_source"), verdict))
    print("    (θ not needed here — the free clique-cover bound already closed the gap)\n")

    # Part 2 — a confusability graph that is a 5-cycle (each oligo cross-reacts with two neighbours).
    # This is the case the free bound CANNOT certify: clique-cover of C5 = 3, but α(C5) = 2. The exact branch-and-bound search proves the panel of 2 is maximum; theta is shown only as a diagnostic.
    C5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    mis, _, _ = max_independent_set(5, C5)
    cc = clique_cover_bound(5, C5)
    th = lovasz_theta(5, C5)
    print("[2] a 5-way cross-reaction cycle (C5)")
    print("    exact branch-and-bound panel: %d oligos (model maximum proven)" % len(mis))
    print("    greedy clique-cover upper bound: %d  (valid but loose)" % cc)
    if th is not None:
        print("    numerical Lovász θ diagnostic: %.4f  (not used as a formal certificate)" % th)
        for kk in (2, 3):
            cap = strong_product_capacity(th, kk)
            print("    split-pool k=%d: θ^%d diagnostic %.3f; constructed |MIS|^%d=%d"
                  % (kk, kk, cap, kk, len(mis) ** kk))
    else:
        print("    numerical Lovász θ diagnostic unavailable (optional cvxpy dependency)")


if __name__ == "__main__":
    _demo()
