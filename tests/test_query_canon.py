"""Gene/marker query canonicalization (oligoforge/ncbi.canonicalize_query). OFFLINE — pure string logic,
no network. Guards the reproducibility fix: a gene must retrieve the SAME records regardless of spelling
("cytb" vs "cytochrome b"), without over-matching lookalike symbols or clobbering explicit field tags.
Run from repo root:  OLIGOFORGE_EMAIL=you@x python3 tests/test_query_canon.py   (exit 0 = pass)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import ncbi as N

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (("  [%s]" % detail) if detail and not cond else ""))
    if not cond:
        fails.append(name)

def exp(s):
    return N.canonicalize_query(s)["expanded"]

# 1. Every spelling of cytochrome b -> one identical expanded term (the reported bug)
cytb = ["plasmodium cytb", "plasmodium cytochrome b", "plasmodium cyt b", "plasmodium cyt-b",
        "plasmodium CYTB", "plasmodium Cytochrome B", "plasmodium cob", "plasmodium MT-CYB"]
outs = {exp(s) for s in cytb}
check("all cytb spellings converge to one term", len(outs) == 1, "%d distinct" % len(outs))
check("expanded term is an OR group with [All Fields]",
      outs and outs.copy().pop().count(" OR ") >= 6 and "[All Fields]" in outs.copy().pop())

# 2. Offtarget phrasing converges too (target AND offtarget both went through this)
check("haemoproteus cytb == haemoproteus cytochrome b",
      exp("haemoproteus cytb") == exp("haemoproteus cytochrome b"))

# 3. COI / barcode family converges across abbreviation and full name
coi = {exp("Aphelocoma COI"), exp("Aphelocoma cox1"), exp("Aphelocoma CO1"),
       exp("Aphelocoma cytochrome c oxidase subunit 1")}
check("COI spellings converge", len(coi) == 1, "%d distinct" % len(coi))
# COII must NOT collide with COI
check("COI and COII are different genes", exp("x COI") != exp("x COII"))

# 4. rRNA / ITS converge
check("16S spellings converge", exp("Bacterium 16S") == exp("Bacterium 16S rRNA") == exp("Bacterium 16s ribosomal RNA"))
check("ITS spellings converge", exp("Fungus ITS") == exp("Fungus internal transcribed spacer"))

# 5. No over-matching: unknown symbols and lookalikes are returned unchanged
for s in ["plasmodium AMA1", "plasmodium MSP1", "human COX10", "bacterium 216s region", "Homo sapiens IFNG"]:
    check("unchanged (no false expansion): %s" % s, exp(s) == s)

# 6. Explicit field tags are respected (the precise autodesign [Gene Name] path must survive)
for s in ["MT-CYB[Gene Name] AND Plasmodium[Organism]", "CYTB[Gene] AND foo", "cox1[Title]"]:
    check("field-tagged token untouched: %s" % s, exp(s) == s)

# 7. Determinism
check("deterministic", N.canonicalize_query("plasmodium cytb") == N.canonicalize_query("plasmodium cytb"))

# 8. taxon extraction strips the gene phrase (drives isolate marker mode)
r = N.canonicalize_query("plasmodium cytochrome b")
check("taxon is the organism remainder", r["taxon"].strip().lower() == "plasmodium" and bool(r["groups"]))

# 9. search_genomes now enters MARKER mode on the full name just like the abbreviation
def mode(q):
    cg = N.canonicalize_query((q,) and q)
    return "MARKER" if cg["groups"] else "GENOME"
check("isolate panel: 'cytb' and 'cytochrome b' both -> MARKER mode",
      mode("plasmodium cytb") == mode("plasmodium cytochrome b") == "MARKER")

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL QUERY-CANONICALIZATION ASSERTS PASS")
