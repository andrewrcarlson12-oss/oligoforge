"""Offline tests for the NCBI retry classifier and the dynamic GenBank-driven gene suggester.

No network: the NCBI access points (count_hits / genes_for_organism) are injected or monkeypatched,
so this exercises the pure merge/rank/fallback logic deterministically. Script-style: exit 0 = pass.
"""
import sys, urllib.error
sys.path.insert(0, ".")
from oligoforge import ncbi as N, refmarkers as RM

fail = 0
def check(name, cond):
    global fail
    print(("  ok  " if cond else "  FAIL ") + name)
    if not cond:
        fail += 1

# --- 1) retry classifier recognises NCBI's parse-time backend error and transient HTTP codes ---
check("'Search Backend failed' is transient",
      N._is_transient(RuntimeError("Search Backend failed: \t\t\n")) is True)
check("'error reading from backend' is transient",
      N._is_transient(RuntimeError("Error reading from backend")) is True)
check("HTTP 400 is transient (NCBI flakes under load)",
      N._is_transient(urllib.error.HTTPError("u", 400, "Bad Request", {}, None)) is True)
check("HTTP 502 is transient",
      N._is_transient(urllib.error.HTTPError("u", 502, "x", {}, None)) is True)
check("HTTP 404 is NOT transient",
      N._is_transient(urllib.error.HTTPError("u", 404, "x", {}, None)) is False)
check("a plain unrelated error is NOT transient",
      N._is_transient(ValueError("nope")) is False)

# --- 2) merge_and_rank: dedup synonyms, drop real 0, keep failed (None), sort desc, off-target ---
base = [
    {"gene": "cytb", "name": "cytochrome b", "qwords": "cytochrome b", "conservation": "variable",
     "chemistry": "parasite_mtdna", "good_for": "x", "copies": "mtDNA", "level": "genotyping", "region": ""},
    {"gene": "18S rRNA", "name": "18S", "qwords": "18S ribosomal RNA", "conservation": "conserved",
     "chemistry": "idt_taqman", "good_for": "y", "copies": "multicopy", "level": "detection", "region": ""},
    {"gene": "rbcL", "name": "rbcL", "qwords": "rbcL", "conservation": "moderate",
     "chemistry": "sybr_generic", "good_for": "z", "copies": "plastid", "level": "id", "region": ""},
]
disc = [{"symbol": "COB", "name": "cytochrome b (dup via synonym -> dropped)"},
        {"symbol": "msp1", "name": "merozoite surface protein 1"}]
counts = {
    "Plasmodium[Organism] AND cytochrome b AND 50:120000[SLEN]": 14000,
    "Plasmodium[Organism] AND 18S ribosomal RNA AND 50:120000[SLEN]": 900,
    "Plasmodium[Organism] AND rbcL AND 50:120000[SLEN]": 0,         # absent -> dropped
    "Plasmodium[Organism] AND msp1 AND 50:120000[SLEN]": None,      # call failed -> kept as '-'
}
off = {"Haemoproteus[Organism] AND cytochrome b AND 50:120000[SLEN]": 1200}
res = RM.merge_and_rank(base, disc, "Plasmodium", lambda t: counts.get(t),
                        exclude="Haemoproteus", off_count_fn=lambda t: off.get(t))
genes = [m["gene"] for m in res]
check("zero-count marker dropped (rbcL)", "rbcL" not in genes)
check("synonym duplicate dropped (COB == cytb)", "COB" not in genes)
check("failed-count (None) marker kept (msp1)", "msp1" in genes)
check("sorted by records desc -> cytb first", genes[0] == "cytb")
check("exact order is cytb, 18S, msp1", genes == ["cytb", "18S rRNA", "msp1"])
cytb = [m for m in res if m["gene"] == "cytb"][0]
check("cytb carries live record count", cytb.get("count") == 14000)
check("cytb carries off-target count", cytb.get("off_count") == 1200)
check("discovered gene tagged source=ncbi_gene",
      [m for m in res if m["gene"] == "msp1"][0].get("source") == "ncbi_gene")

# --- 3) suggest_dynamic falls back to the curated list when GenBank is unreachable ---
_orig_count, _orig_genes = N.count_hits, N.genes_for_organism
try:
    N.count_hits = lambda q: None              # every count call fails
    N.genes_for_organism = lambda o, n=18: []  # no discovery
    d = RM.suggest_dynamic("Plasmodium", None, "any")
    check("fallback: does not claim scanned when no live counts", not d.get("scanned"))
    check("fallback: still returns curated markers", bool(d.get("markers")))
    check("fallback: note is honest about not reaching GenBank", "curated" in (d.get("note") or "").lower())
finally:
    N.count_hits, N.genes_for_organism = _orig_count, _orig_genes

# --- 4) suggest_dynamic uses live counts when available (still injected, no network) ---
try:
    N.genes_for_organism = lambda o, n=18: [{"symbol": "msp1", "name": "merozoite surface protein 1"}]
    N.count_hits = lambda q: 500 if "msp1" in q or "cytochrome" in q else 0
    d = RM.suggest_dynamic("Plasmodium", None, "any")
    check("scanned=True when counts available", d.get("scanned") is True)
    check("dynamic flag set", d.get("dynamic") is True)
    check("discovered gene present in markers", any(m.get("gene") == "msp1" for m in d.get("markers", [])))
finally:
    N.count_hits, N.genes_for_organism = _orig_count, _orig_genes

print("RESULT:", "ALL PASS" if fail == 0 else ("%d FAILED" % fail))
sys.exit(1 if fail else 0)
