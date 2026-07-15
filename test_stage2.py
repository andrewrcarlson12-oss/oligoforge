"""Stage-2 smoke checks.

The default release gate is deterministic and offline. Set
OLIGOFORGE_LIVE_NCBI=1 to additionally exercise live NCBI retrieval.
"""
import os
import sys
sys.path.insert(0, ".")

from oligoforge import ncbi, profiles as P

F = "GAGCTATACCCCGACCTCTG"
PR = "ATCTTGTCCCCAGTTGTTGACATGGCC"

print("=" * 72, "\n1) MULTI-VENDOR LINT\n", "=" * 72)
for prof in ("idt_taqman", "thermo_taqman", "sybr_generic"):
    findings = P.lint_oligo(F, "forward", P.PROFILES[prof])
    assert isinstance(findings, list)
    print(f"  {prof}: {len(findings)} finding(s)")
probe_findings = P.lint_oligo(PR, "probe", P.PROFILES["idt_taqman"])
assert isinstance(probe_findings, list)
print(f"  IDT probe: {len(probe_findings)} finding(s)")

print("\n" + "=" * 72, "\n2) ISOFORM-COMMON REGION (OFFLINE)\n", "=" * 72)
seqs = [
    "TTTT" + "ACGT" * 30 + "AAAA",
    "GGGG" + "ACGT" * 30 + "CCCC",
    "NNNN" + "ACGT" * 30 + "TTTT",
]
common = ncbi.common_region(seqs)
assert common and common in seqs[0] and common in seqs[1] and common in seqs[2]
print(f"  common block: {len(common)} bp")

if os.environ.get("OLIGOFORGE_LIVE_NCBI") == "1":
    ncbi.Entrez.email = os.environ.get("OLIGOFORGE_EMAIL", "ci@example.com")
    print("\n" + "=" * 72, "\n3) LIVE NCBI OPTIONAL CHECK\n", "=" * 72)
    recs, query = ncbi.fetch_isoforms("HMBS", "Aphelocoma coerulescens")
    assert recs
    print(f"  query: {query}\n  isoforms: {len(recs)}")
    gid = ncbi.gene_id("HMBS", "Aphelocoma coerulescens")
    assert gid
    print(f"  gene id: {gid}")
else:
    print("\nSKIP live NCBI checks (set OLIGOFORGE_LIVE_NCBI=1 to enable).")

print("ALL STAGE-2 ASSERTS PASS")
