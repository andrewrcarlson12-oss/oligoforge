import sys; sys.path.insert(0, ".")
from oligoforge import thermo as T, design as D
from oligoforge.profiles import IDT_TAQMAN as C

# Real FSJ HMBS sequence (the verified gBlock region, XM_068994916 vicinity).
HMBS = ("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATT"
        "GTGGCCATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCT"
        "TCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")

print("=== engine sanity: known oligo values (should match prior manual QC) ===")
for label, seq in [("HMBS_F","GAGCTATACCCCGACCTCTG"),
                   ("HMBS_R","CTTCTCTCCAATCTTGGAAAGCG"),
                   ("HMBS_P","ATCTTGTCCCCAGTTGTTGACATGGCC")]:
    print(f"  {label:7s} Tm={T.tm(seq):.1f} GC={T.gc_percent(seq):.0f} "
          f"hairpin={T.hairpin(seq)[0]:.2f} selfD={T.self_dimer(seq):.2f}")

print("\n=== full auto-design on the real HMBS template ===")
a = D.design_assay(HMBS, C)
if not a:
    print("  NO ASSAY FOUND"); sys.exit(1)
pi = a["probe_info"]
print(f"  F  {a['forward']}   Tm={a['f_tm']:.1f}")
print(f"  R  {a['reverse']}   Tm={a['r_tm']:.1f}")
print(f"  P  {a['probe']}   Tm={pi['tm']:.1f}  offset=+{pi['offset']:.1f}  "
      f"hairpin={pi['hairpin_dg']:.2f}  5'={a['probe'][:4]}")
print(f"  amplicon={a['amplicon']}bp  pair_Tm_gap={a['pair_tm_gap']:.1f}  "
      f"FxR={T.hetero_dimer(a['forward'],a['reverse']):.2f}")
print(f"  PxF={pi['dimer_f']:.2f}  PxR={pi['dimer_r']:.2f}")
print(f"  gBlock ({len(a['gblock'])}bp): {a['gblock']}")

# does it recover the hand-designed assay?
match = (a['forward']=="GAGCTATACCCCGACCTCTG" and a['reverse']=="CTTCTCTCCAATCTTGGAAAGCG")
print(f"\n  recovers hand-designed primers exactly: {match}")
print("  (probe may differ — engine optimizes hairpin/offset over the interior)")
