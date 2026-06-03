"""Golden-output gate for the auto-design discrimination + degeneracy engine. Offline.
Run from repo root:  python tests/test_autodesign_golden.py

WHY THIS EXISTS
`design_assay` (the single-template path) is pinned byte-for-byte by test_regression's HMBS
anchor. But the scientific brain -- `design_from_sequences` with its multi-candidate ranking,
3'-block discrimination scoring, and IUPAC-degenerate genus design -- was, until now, only ever
checked by live NCBI runs recorded in the handoff. That path was revised across ~8 consecutive
releases (v1.11.11 -> v1.11.20). A scorer change could silently flip the winning design, the
ranking, the degenerate codes, or the discrimination mechanism, and no offline test would catch
it. This freezes the *decisions* of two canonical designs against SAVED FASTA fixtures (no
network, fully deterministic since primer3 is deterministic): the Plasmodium-vs-Haemoproteus
cytb discrimination design, and the Plasmodium genus (degenerate) design. Tolerant on float
scores; exact on the sequences, ranking, degeneracy, and the 3'-block call -- the things that
matter scientifically.

If a deliberate engine change moves these, re-capture the goldens and review the diff. An
accidental regression fails loudly.
"""
import json, os, sys
sys.path.insert(0, ".")
from oligoforge import autodesign as AD, profiles as P

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
FX = os.path.join(HERE, "fixtures")

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)

plas = json.load(open(os.path.join(FX, "plasmodium_cytb.json")))["sequences"]
haem = json.load(open(os.path.join(FX, "haemoproteus_cytb.json")))["sequences"]
prof = P.PROFILES["parasite_mtdna"]
check("fixtures load (11 Plasmodium / 12 Haemoproteus cytb)", len(plas) == 11 and len(haem) == 12,
      "%d / %d" % (len(plas), len(haem)))

# ---------------- 1) DISCRIMINATION design: Plasmodium target vs Haemoproteus off-target ----------------
d = AD.design_from_sequences(plas, prof, offs=haem, n_candidates=5)
check("disc design returns candidates", not d.get("error") and d.get("candidates"), d.get("error"))
if d.get("candidates"):
    check("disc n_targets/n_offs", d["n_targets"] == 11 and d["n_offs"] == 12, (d["n_targets"], d["n_offs"]))
    check("disc enumerates 10 candidates", d["n_candidates"] == 10, d["n_candidates"])
    w = d["candidates"][0]["assay"]
    check("disc winner F", w["forward"] == "TTTCCATTTATAGCCTTATGTATTG", w["forward"])
    check("disc winner R", w["reverse"] == "TTTTAAAGCTGTATCATACCCT", w["reverse"])
    check("disc winner P", w.get("probe") == "ACATTTACAAGGTAGCACAAATCCTTT", w.get("probe"))
    check("disc winner amplicon 96", w.get("amplicon") == 96, w.get("amplicon"))
    # IUPAC degenerate genus forms across the 11 targets
    check("disc winner F_deg", w.get("forward_deg") == "TTTCCWTTTATAGCYTTATGTATTG", w.get("forward_deg"))
    check("disc winner P_deg", w.get("probe_deg") == "ACATTTACAAGGTAGCACWAATCCTTT", w.get("probe_deg"))
    check("disc winner n_degenerate 3", w.get("n_degenerate") == 3, w.get("n_degenerate"))
    # ranking (the scorer's order)
    rank = [c["assay"]["forward"] for c in d["candidates"][:3]]
    check("disc top-3 ranking stable",
          rank == ["TTTCCATTTATAGCCTTATGTATTG", "TTTCTACATTTACAAGGTAGCA", "ATATCAATAGTTACTGCTTTTATGG"], rank)
    # THE discrimination decision: the forward 3'-blocks Haemoproteus (>=1 terminal mismatch),
    # and the off-target is genuinely separated (median identity well below 100%).
    disc = d["candidates"][0].get("discrimination") or {}
    fdisc = disc.get("F") or {}
    check("disc forward 3'-blocks the off-target (>=1 terminal mismatch)",
          (fdisc.get("min_3prime_mismatch") or 0) >= 1, fdisc.get("min_3prime_mismatch"))
    check("disc forward off-target median identity in [80,90]",
          80.0 <= (fdisc.get("median_ident") or 0) <= 90.0, fdisc.get("median_ident"))

# ---------------- 2) GENUS design (no off-target): degenerate multi-template ----------------
g = AD.design_from_sequences(plas, prof, n_candidates=5)
check("genus design returns candidates", not g.get("error") and g.get("candidates"), g.get("error"))
if g.get("candidates"):
    check("genus n_offs == 0", g["n_offs"] == 0, g["n_offs"])
    gw = g["candidates"][0]["assay"]
    check("genus winner F", gw["forward"] == "TTTCTACATTTACAAGGTAGCA", gw["forward"])
    check("genus winner R", gw["reverse"] == "TCAAGACTTAATAGATTTGGATAGA", gw["reverse"])
    check("genus winner P", gw.get("probe") == "CAAATCCTTTAGGGTATGATACAGCTT", gw.get("probe"))
    check("genus winner amplicon 86", gw.get("amplicon") == 86, gw.get("amplicon"))
    check("genus winner F_deg", gw.get("forward_deg") == "TTYYTACATTTACAAGGTAGCA", gw.get("forward_deg"))
    check("genus winner R_deg", gw.get("reverse_deg") == "TCAAGACTTAAWAGATTTGGATAGA", gw.get("reverse_deg"))
    check("genus winner n_degenerate 4", gw.get("n_degenerate") == 4, gw.get("n_degenerate"))
    check("genus deg_targets 11", gw.get("deg_targets") == 11, gw.get("deg_targets"))

print("")
if _fails:
    print("AUTODESIGN GOLDEN GATE FAILED:", ", ".join(_fails)); sys.exit(1)
print("ALL AUTODESIGN GOLDEN ASSERTS PASS")
