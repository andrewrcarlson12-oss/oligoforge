import sys; sys.path.insert(0,".")
from oligoforge import specificity as SP, ncbi
ncbi.Entrez.email = SP.Entrez.email = "fsj.qpcr.design@gmail.com"
cases = [("HMBS","Aphelocoma coerulescens","XM_068994916.1",223,315),
         ("SDHA","Aphelocoma coerulescens","XM_069004197.1",1637,1766)]
for gene,org,acc,a,b in cases:
    print("="*64); print(f"{gene}  amplicon {acc} pos {a}-{b}"); print("="*64)
    r = SP.intron_check(gene,org,a,b,mrna_acc=acc)
    print(f"  {r['info']}")
    print(f"  junctions (mRNA): {r['junctions']}")
    print(f"  spanned by amplicon: {r['spanned']}")
    print(f"  VERDICT: {r['verdict']}")
