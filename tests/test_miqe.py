"""Validate the MIQE consolidation module (oligoforge/miqe.py). Offline.
Run from repo root:  OLIGOFORGE_EMAIL=you@x python3 tests/test_miqe.py   (exit 0 = pass)

Checks that the record correctly fuses engine-predicted thermodynamics (incl. the LNA probe), applies
MIQE acceptance criteria to attached empirical data, runs the predicted-vs-observed amplicon-Tm check,
and stamps a deterministic provenance block."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import miqe as MIQE, thermo as T

fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (("  [%s]" % detail) if detail and not cond else ""))
    if not cond:
        fails.append(name)

def status(rec, cid):
    return next(c["status"] for c in rec["checks"] if c["id"] == cid)

T.set_conditions(mv_conc=50, dv_conc=3, dntp_conc=0.8, dna_conc=200, anneal_c=60)
IFNG = dict(name="IFNG", gene="IFNG", organism="Aphelocoma coerulescens",
            chemistry="IDT PrimeTime (ZEN 5'FAM)", dye="FAM", amplicon=136, amplicon_tm=82.0,
            forward="AGTCATTCTGATGTCGCTGATG", reverse="ACCTGTCAGTGTTTTCAAGCA",
            probe="TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT")

# 1. No empirical data: predicted/binding/provenance scored; empirical checks are n/a; overall not fail
r = MIQE.validate_assay(IFNG)
check("no-empirical: overall pass (no fails)", r["miqe_status"] == "pass")
check("no-empirical: efficiency n/a", status(r, "efficiency") == "na")
check("no-empirical: provenance pass", status(r, "provenance") == "pass")
check("predicted thermodynamics present for F/R/probe",
      all(k in r["predicted"]["oligos"] for k in ("forward", "reverse", "probe")))
check("probe binding scored from engine", status(r, "probe_binding") in ("pass", "warn", "fail"))

# 2. Full empirical, good assay: every criterion passes incl. predicted-vs-observed amplicon Tm
full = dict(IFNG); full["validation"] = {"efficiency_pct": 98.5, "r2": 0.997, "slope": -3.36, "lod": 10}
rf = MIQE.validate_assay(full, observed_amplicon_tm=81.5, observed_peaks=1)
check("full-empirical: overall pass", rf["miqe_status"] == "pass" and rf["n_missing"] == 0)
check("efficiency pass", status(rf, "efficiency") == "pass")
check("linearity pass", status(rf, "linearity") == "pass")
check("single product pass", status(rf, "single_product") == "pass")
check("amplicon Tm predicted-vs-observed pass (|81.5-82.0|<3)", status(rf, "amplicon_tm") == "pass")

# 3. Out-of-spec efficiency fails; far-off observed Tm fails the predicted-vs-observed check
bad = dict(IFNG); bad["validation"] = {"efficiency_pct": 70.0, "r2": 0.90, "slope": -4.5}
rb = MIQE.validate_assay(bad, observed_amplicon_tm=70.0)
check("low efficiency -> fail", status(rb, "efficiency") == "fail")
check("R2 0.90 -> fail", status(rb, "linearity") == "fail")
check("observed Tm far from predicted -> fail", status(rb, "amplicon_tm") == "fail")
check("overall fail", rb["miqe_status"] == "fail")

# 4. LNA probe scored via McTigue (not skipped), and predicted bound
genus = dict(name="Plasmodium", gene="cytb", organism="Plasmodium", chemistry="LNA TaqMan", dye="FAM",
             amplicon=157, amplicon_tm=78.0, forward="TACCTGGACTWGTTTCATGG", reverse="AAAGGATTTGTGCTACCTTG",
             probe="/56-FAM/CTTA+CA+A+GATAT+CC+ACCACA/3IABkFQ/")
rg = MIQE.validate_assay(genus)
pp = rg["predicted"]["oligos"].get("probe")
check("LNA probe scored via McTigue", pp is not None and pp["basis"] == "McTigue LNA" and pp.get("lna_n") == 5)
check("LNA probe predicted bound (pass)", status(rg, "probe_binding") == "pass")

# 5. Provenance is deterministic and complete
r2 = MIQE.validate_assay(IFNG)
check("checksums deterministic", r2["provenance"]["checksums"] == r["provenance"]["checksums"])
check("checksum format", str(r["provenance"]["checksums"]["forward"]).startswith("sha256:"))
check("version + utc stamped", bool(r["provenance"]["oligoforge_version"]) and r["provenance"]["generated_utc"].endswith("Z"))

# 6. Markdown report renders the key sections
md = rf["markdown"]
check("markdown has MIQE checklist + provenance + predicted",
      "MIQE checklist" in md and "Provenance" in md and "Predicted thermodynamics" in md)

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL MIQE ASSERTS PASS")
