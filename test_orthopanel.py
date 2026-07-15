"""Math-core tests for the Certified Orthogonal Panel (orthopanel). Offline, no thermodynamics.
Run from repo root:  python tests/test_orthopanel.py

Covers the graph/certificate math independently of any oligo input:
  * Lovász θ against known values (empty->n, complete->1, C5->√5, Petersen->4). The C5=√5 test
    genuinely catches a wrong SDP formulation. θ tests SKIP cleanly if cvxpy is not installed
    (θ is an optional dependency), so this file passes on a machine without cvxpy.
  * The load-bearing INVARIANT: |MIS| <= ⌊θ+eps⌋ and |MIS| <= clique-cover bound, on many random
    graphs. A violation means a bug in the MIS, the SDP, or the bounds.
  * integer_upper_bound rounding: a value tuned to sit just below an integer must NOT floor down.
  * MIS exactness vs brute force on small graphs; the sandwich α <= θ <= clique-cover.
  * Edge cases: empty graph, single vertex, complete, edgeless, and the greedy large-n fallback flag.
"""
import math, os, sys, random, itertools
sys.path.insert(0, ".")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from oligoforge import orthopanel as OP

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)

HAVE_THETA = OP.lovasz_theta(3, [(0, 1)]) is not None
if not HAVE_THETA:
    print("  NOTE cvxpy not installed — θ tests will be skipped (θ is optional; clique-cover certifies)")


def brute_alpha(n, edges):
    """Exact independence number by brute force (small n only), for cross-checking the B&B."""
    adj = OP._adjacency(n, edges)
    best = 0
    for r in range(n, -1, -1):
        for comb in itertools.combinations(range(n), r):
            cs = set(comb)
            if all(not (adj[v] & cs) for v in comb):
                return r
    return best


# ---------------- θ known values ----------------
if HAVE_THETA:
    check("θ(empty n=5) = 5", abs(OP.lovasz_theta(5, []) - 5) < 3e-2, OP.lovasz_theta(5, []))
    K4 = list(itertools.combinations(range(4), 2))
    check("θ(K4 complete) = 1", abs(OP.lovasz_theta(4, K4) - 1) < 3e-2, OP.lovasz_theta(4, K4))
    C5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    check("θ(C5) = √5 ≈ 2.236 (catches a wrong SDP)",
          abs(OP.lovasz_theta(5, C5) - math.sqrt(5)) < 3e-2, OP.lovasz_theta(5, C5))
    petersen = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
                (5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]
    check("θ(Petersen) = 4", abs(OP.lovasz_theta(10, petersen) - 4) < 5e-2, OP.lovasz_theta(10, petersen))
    check("⌊θ(C5)⌋ = 2", OP.integer_upper_bound(OP.lovasz_theta(5, C5)) == 2)
    check("⌊θ(Petersen)⌋ = 4 (snaps 4.000x, no off-by-one)",
          OP.integer_upper_bound(OP.lovasz_theta(10, petersen)) == 4)

# ---------------- integer_upper_bound rounding (the +eps / snap guard) ----------------
check("integer_upper_bound(3.9999) = 4 (not floored to 3)", OP.integer_upper_bound(3.9999) == 4)
check("integer_upper_bound(4.0001) = 4", OP.integer_upper_bound(4.0001) == 4)
check("integer_upper_bound(2.236) = 2 (genuine non-integer floors)", OP.integer_upper_bound(2.236) == 2)
check("integer_upper_bound(4.60) = 4", OP.integer_upper_bound(4.60) == 4)
check("integer_upper_bound(None) = None", OP.integer_upper_bound(None) is None)

# ---------------- sandwich + invariant on random graphs ----------------
random.seed(1234)
inv_ok = True
sand_ok = True
mis_exact_ok = True
for trial in range(40):
    n = random.randint(2, 11)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.35]
    mis, method, exact = OP.max_independent_set(n, edges)
    lo = len(mis)
    # returned set is genuinely independent
    adj = OP._adjacency(n, edges)
    if any(adj[u] & (set(mis) - {u}) for u in mis):
        inv_ok = False
    cc = OP.clique_cover_bound(n, edges)
    if lo > cc:
        inv_ok = False
    # MIS is exact and equals brute-force α on these small graphs
    if exact and lo != brute_alpha(n, edges):
        mis_exact_ok = False
    if HAVE_THETA:
        th = OP.lovasz_theta(n, edges)
        ub = OP.integer_upper_bound(th)
        if lo > ub:
            inv_ok = False
        # sandwich α <= θ <= clique-cover (allow small solver slack on θ)
        if not (lo <= th + 1e-2 and th <= cc + 1e-2):
            sand_ok = False
check("invariant |MIS| <= ⌊θ⌋ and <= clique-cover over 40 random graphs", inv_ok)
check("MIS is exact = brute-force α on small graphs", mis_exact_ok)
if HAVE_THETA:
    check("sandwich α <= θ <= clique-cover holds", sand_ok)

# ---------------- edge cases ----------------
mis0, m0, e0 = OP.max_independent_set(0, [])
check("empty graph: MIS = [] and exact", mis0 == [] and e0)
check("clique_cover_bound(empty) = 0", OP.clique_cover_bound(0, []) == 0)
mis1, _, e1 = OP.max_independent_set(1, [])
check("single vertex: MIS size 1", len(mis1) == 1 and e1)
Kn = list(itertools.combinations(range(6), 2))
misK, _, _ = OP.max_independent_set(6, Kn)
check("complete K6: MIS size 1", len(misK) == 1)
check("clique_cover_bound(K6) = 1", OP.clique_cover_bound(6, Kn) == 1)
misE, _, _ = OP.max_independent_set(6, [])
check("edgeless n=6: MIS size 6", len(misE) == 6)
check("clique_cover_bound(edgeless 6) = 6", OP.clique_cover_bound(6, []) == 6)

# large-n greedy fallback is flagged non-exact
big_n = 200
big_edges = [(i, j) for i in range(big_n) for j in range(i + 1, big_n) if random.random() < 0.02]
_, big_method, big_exact = OP.max_independent_set(big_n, big_edges, exact_limit=150)
check("n>exact_limit falls back to greedy, flagged non-exact", (not big_exact) and "greedy" in big_method)

# ---------------- strong-product capacity ----------------
check("strong_product_capacity(√5, 2) = 5", abs(OP.strong_product_capacity(math.sqrt(5), 2) - 5) < 1e-9)
check("strong_product_capacity(2, 3) = 8", OP.strong_product_capacity(2, 3) == 8)
check("strong_product_capacity(None, k) = None", OP.strong_product_capacity(None, 3) is None)

print(("FAIL: " + ", ".join(_fails)) if _fails else "OK")
sys.exit(1 if _fails else 0)
