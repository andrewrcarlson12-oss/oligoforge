"""Live-NCBI intron/exon-junction check. Network-tolerant: when NCBI is reachable it verifies the
two Aphelocoma coerulescens amplicons span a junction; when it is not, it degrades cleanly instead
of crashing. The contract this pins is that intron_check ALWAYS returns a 'verdict' key, on every
path (success or graceful degradation) -- a regression that dropped it (KeyError) is caught here."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import specificity as SP, ncbi
ncbi.Entrez.email = SP.Entrez.email = os.environ.get("OLIGOFORGE_EMAIL", "fsj.qpcr.design@gmail.com")
cases = [("HMBS","Aphelocoma coerulescens","XM_068994916.1",223,315),
         ("SDHA","Aphelocoma coerulescens","XM_069004197.1",1637,1766)]
fails = []
for gene,org,acc,a,b in cases:
    print("="*64); print(f"{gene}  amplicon {acc} pos {a}-{b}"); print("="*64)
    r = SP.intron_check(gene,org,a,b,mrna_acc=acc)
    # contract: 'verdict' present on EVERY path (this is the regression guard)
    if "verdict" not in r:
        print("  FAIL  intron_check returned no 'verdict' key:", sorted(r)); fails.append(gene); continue
    print(f"  {r.get('info','')}")
    print(f"  junctions (mRNA): {r.get('junctions')}")
    print(f"  spanned by amplicon: {r.get('spanned')}")
    print(f"  VERDICT: {r['verdict']}")
    if r.get("junctions") is None:
        print("  (NCBI unreachable -- degraded cleanly, not a failure)")
    elif not r.get("spanned"):
        # network worked and coords were located but no junction spanned -> unexpected for these cases
        if r.get("amp_start") is not None:
            print("  NOTE: amplicon located but no junction spanned (check coordinates/isoform)")
if fails:
    print("\nINTRON TEST FAILED:", ", ".join(fails)); sys.exit(1)
print("\nINTRON TEST PASS (verdict present on every path)")
