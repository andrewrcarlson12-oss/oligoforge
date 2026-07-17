"""Regression gate for the autodesign query-path annotation step. Offline, no NCBI.
Run from repo root:  python tests/test_autodesign_annotate.py

WHY THIS EXISTS
`design_from_query` (the "Autodesign from a target" path in the UI) calls the private
`_annotate(out, ref, prefer_junction)` helper to attach amplicon secondary-structure (and,
optionally, exon-junction) information to each candidate. That helper folded the amplicon with
`profile.get("anneal_c", ...)` -- but `profile` was never a parameter or local of `_annotate`,
so the call raised `NameError: name 'profile' is not defined` for EVERY candidate whenever
ViennaRNA was present (i.e. on any real install), surfacing to the user as
"autodesign failed: name 'profile' is not defined".

The rest of the suite never caught it because the golden/regression tests exercise
`design_assay` and `design_from_sequences` DIRECTLY; none of them route through `_annotate`,
which lives only on the fetch->design->annotate query path. This test closes that gap by
building a real `out` from an offline design and calling `_annotate` itself.

The portable assertions are (1) `_annotate` does not raise, and (2) it runs the per-candidate
loop (amp_span populated). When ViennaRNA is available it must ALSO attach the `structure`
block -- proving the previously-crashing fold line now executes. No absolute structure numbers
are pinned, so the test is not flaky across machines/ViennaRNA versions.
"""
import json, os, sys
sys.path.insert(0, ".")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from oligoforge import autodesign as AD, profiles as P, structure as STR

FX = os.path.join(HERE, "fixtures")

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)

# ---- build a real, offline design (same engine the query path uses) --------------------------
plas = json.load(open(os.path.join(FX, "plasmodium_cytb.json")))["sequences"]
prof = P.PROFILES["parasite_mtdna"]
out = AD.design_from_sequences(plas, prof, n_candidates=5)
check("offline design yields candidates", not out.get("error") and out.get("candidates"), out.get("error"))

# the query path sets profile_used before calling _annotate; mirror that
out["profile_used"] = "parasite_mtdna"
ref = AD._reference(plas)
check("reference builds", bool(ref) and len(ref) > 100, len(ref) if ref else 0)

# ---- the actual regression: _annotate must not raise NameError -------------------------------
raised = None
try:
    AD._annotate(out, ref, prefer_junction=False)
except Exception as e:
    raised = e
check("_annotate runs without raising (was NameError: 'profile')", raised is None,
      f"{type(raised).__name__}: {raised}" if raised else "")

cands = out.get("candidates") or []
check("_annotate populated amp_span on candidates (loop ran)",
      bool(cands) and any(c.get("amp_span") for c in cands),
      "no amp_span on any candidate")

# ---- when ViennaRNA is present, the previously-crashing fold line must now execute ------------
if STR.available():
    has_struct = any(c.get("structure") for c in cands)
    check("structure block attached (fold line executed under ViennaRNA)", has_struct,
          "no candidate got a structure block")
    # the anneal temperature actually used should be the profile's, or the global default
    _prof_anneal = prof.get("anneal_c")
    if has_struct and _prof_anneal is not None:
        s = next(c["structure"] for c in cands if c.get("structure"))
        check("fold used the profile anneal temperature", s.get("anneal_c") == _prof_anneal,
              f"{s.get('anneal_c')} != profile {_prof_anneal}")
else:
    print("  SKIP structure assertions (ViennaRNA not available on this host)")

# ---- prefer_junction=True must also be crash-safe offline (NCBI unreachable -> handled) -------
out2 = AD.design_from_sequences(plas, prof, n_candidates=3)
out2["profile_used"] = "parasite_mtdna"
out2["source_accession"] = None  # no accession -> junction branch is skipped, must not raise
raised2 = None
try:
    AD._annotate(out2, ref, prefer_junction=True)
except Exception as e:
    raised2 = e
check("_annotate(prefer_junction=True) crash-safe with no accession", raised2 is None,
      f"{type(raised2).__name__}: {raised2}" if raised2 else "")

print(("FAIL: " + ", ".join(_fails)) if _fails else "OK")
sys.exit(1 if _fails else 0)
