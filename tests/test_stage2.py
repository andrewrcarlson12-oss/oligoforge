import sys, time; sys.path.insert(0, ".")
from oligoforge import thermo as T, ncbi, specificity as SP
from oligoforge import profiles as P
ncbi.Entrez.email = SP.Entrez.email = "fsj.qpcr.design@gmail.com"

print("="*72, "\n1) MULTI-VENDOR LINT (HMBS forward primer + probe)\n", "="*72)
F="GAGCTATACCCCGACCTCTG"; Pr="ATCTTGTCCCCAGTTGTTGACATGGCC"
for prof in ("idt_taqman","thermo_taqman","sybr_generic"):
    print(f"\n  [{P.PROFILES[prof]['name']}]  forward primer:")
    for n,s,d in P.lint_oligo(F,"forward",P.PROFILES[prof]):
        print(f"     {s:4} {n:32} {d}")
print("\n  [IDT PrimeTime] probe:")
for n,s,d in P.lint_oligo(Pr,"probe",P.PROFILES["idt_taqman"]):
    print(f"     {s:4} {n:32} {d}")

print("\n"+"="*72, "\n2) NCBI FETCH + ISOFORM-COMMON REGION (HMBS)\n", "="*72)
recs, q = ncbi.fetch_isoforms("HMBS","Aphelocoma coerulescens")
print(f"  query: {q}\n  isoforms: {[r.id for r in recs]}  lengths {[len(r.seq) for r in recs]}")
common = ncbi.common_region([r.seq for r in recs])
print(f"  common-to-all block: {len(common)} bp")

print("\n"+"="*72, "\n3) gene_table raw head (to confirm exon-parse format)\n", "="*72)
gid = ncbi.gene_id("HMBS","Aphelocoma coerulescens")
print(f"  gene id: {gid}")
h=ncbi.Entrez.efetch(db="gene", id=gid, rettype="gene_table", retmode="text"); raw=h.read(); h.close()
print("  ---- first 1200 chars ----")
print(raw[:1200])
