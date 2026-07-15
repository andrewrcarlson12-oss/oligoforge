"""
G7 — LNA validation hardening with an independent MELTING 5 cross-check (pins lna_validation_v2).

MELTING 5 (via rmelting + Java) is an independent, peer-reviewed implementation of the McTigue 2004
LNA parameters. That engine is not part of this offline suite, so this test:

  1. RE-COMPUTES the OligoForge LNA Tm deterministically on the 12 McTigue experimental duplexes and
     pins the PRIMARY validation vs real experimental data (RMSE 1.86, MAE 1.62, bias -0.52, and the
     honest within-band count 8/12 within +/-2 C, 12/12 within +/-3 C).
  2. RE-COMPUTES OligoForge LNA increments on the committed expanded 96-oligo panel and pins that
     they match MELTING's independent implementation of the SAME McTigue parameters to a tiny RMSE
     (< 0.15 C) -- i.e. the LNA layer is correctly implemented, verified against a second engine.
  3. Pins the honest cross-check conclusion: OligoForge (RMSE 1.86) and MELTING (RMSE 2.10) are a
     TIE vs experiment; neither dominates.

Fully offline (reads committed JSON/CSV; no Java, no network). Standalone script.
"""
import sys, os, json, csv, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import thermo as T, nn as NN

fails = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("" if ok else f"  [{detail}]"))
    if not ok: fails.append(name)

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "benchmark")
val = json.load(open(os.path.join(BENCH, "lna_validation_v2.json")))
mct12 = list(csv.DictReader(open(os.path.join(BENCH, "lna_mctigue12_crosscheck.csv"))))
panel = list(csv.DictReader(open(os.path.join(BENCH, "lna_expanded_panel.csv"))))

# McTigue conditions: 5 uM oligo, 1 M Na; engine dna_conc=2500 nM so dnac/2 == Ct/4.
T.set_conditions(mv_conc=1000, dv_conc=0, dntp_conc=0, dna_conc=2500)

# ---- 1. deterministic re-computation vs McTigue EXPERIMENTAL data (primary validation) ----
errs = []
for row in mct12:
    plus = row["oligo"]
    pred = NN.params_lna(plus)["tm"]
    errs.append(pred - float(row["exp_tm"]))
    # committed pred must match a fresh recompute
    check_dev = abs(pred - float(row["oligoforge_tm"]))
    if check_dev >= 0.1:
        check(f"recompute matches committed for {plus}", False, f"dev {check_dev:.2f}")
rmse = math.sqrt(sum(e*e for e in errs)/len(errs))
mae = sum(abs(e) for e in errs)/len(errs)
bias = sum(errs)/len(errs)
within2 = sum(1 for e in errs if abs(e) <= 2.0)
within3 = sum(1 for e in errs if abs(e) <= 3.0)
check("12 McTigue experimental duplexes present", len(mct12) == 12, f"{len(mct12)}")
check("LNA Tm vs experiment: RMSE ~1.86 C", abs(rmse - 1.86) < 0.1, f"RMSE={rmse:.2f}")
check("LNA Tm vs experiment: MAE ~1.62 C", abs(mae - 1.62) < 0.1, f"MAE={mae:.2f}")
check("LNA Tm vs experiment: bias ~-0.52 C", abs(bias - (-0.52)) < 0.1, f"bias={bias:+.2f}")
check("HONEST within-band: 8/12 within +/-2 C (NOT 10/12)", within2 == 8, f"{within2}/12")
check("HONEST within-band: 12/12 within +/-3 C", within3 == 12, f"{within3}/12")

# ---- 2. deterministic LNA increment re-computation vs MELTING (expanded 96 panel) ----
inc_dev = []
for row in panel:
    plus = row["oligo_plus"]
    core = plus.replace("+", "")
    of_inc = NN.params_lna(plus)["tm"] - NN.params(core)["tm"]
    # committed OligoForge increment must match fresh recompute
    inc_dev.append(abs(of_inc - float(row["oligoforge_incr"])))
check("expanded panel has 96 oligos", len(panel) == 96, f"{len(panel)}")
check("recomputed LNA increments match committed values (<0.05 C)", max(inc_dev) < 0.05,
      f"max dev {max(inc_dev):.3f}")

# increment agreement vs MELTING (from committed CSV: incr_diff column)
diffs = [float(r["incr_diff"]) for r in panel]
inc_rmse = math.sqrt(sum(d*d for d in diffs)/len(diffs))
inc_max = max(abs(d) for d in diffs)
check("LNA increment matches MELTING to RMSE < 0.15 C (implementation verified)",
      inc_rmse < 0.15, f"RMSE={inc_rmse:.3f}")
check("LNA increment max|delta| vs MELTING < 0.5 C across all 96", inc_max < 0.5, f"max={inc_max:.3f}")

# ---- 3. pin the honest tie conclusion from the committed validation JSON ----
of_exp = val["melting_crosscheck_12"]["oligoforge_vs_experiment"]
mt_exp = val["melting_crosscheck_12"]["melting_mct04_vs_experiment"]
check("committed: OligoForge RMSE 1.86 vs experiment", abs(of_exp["rmse"] - 1.86) < 0.05)
check("committed: MELTING RMSE ~2.10 vs experiment (honest tie, neither dominates)",
      2.0 < mt_exp["rmse"] < 2.2, f"MELTING RMSE={mt_exp['rmse']}")

T.set_conditions(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=200.0)  # restore
print(f"\n{'ALL LNA HARDENING ASSERTS PASS' if not fails else str(len(fails))+' FAILED'}")
sys.exit(len(fails))
