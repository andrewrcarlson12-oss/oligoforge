"""Tm accuracy calibration against IDT OligoAnalyzer. Offline. Run from repo root:
    python tests/test_tm_calibration.py

WHY THIS EXISTS
The README claims OligoForge's Tm "tracks IDT OligoAnalyzer to ~1-2 C absolute." Until this
file that claim was never tested -- and the only way to test it is with real OligoAnalyzer
numbers, which only a human can pull from IDT's tool. The drift-catching half (that OligoForge's
own Tm doesn't silently shift between builds) is already pinned by test_locked_panel.py against
primer3. THIS file is the *external accuracy* half: it compares OligoForge's Tm to the
OligoAnalyzer Tm for a fixed reference set, asserts every oligo is within OA_TOL_C, and prints
the mean-absolute-error and worst deviation so the accuracy claim is backed by a measured number.

HOW TO ACTIVATE (one-time, ~10 min):
  1. Open https://www.idtdna.com/calc/analyzer
  2. Set the same conditions OligoForge uses (Target type DNA; Oligo conc 0.2 uM = 200 nM;
     Na+ 50 mM; Mg++ 3 mM; dNTPs 0.8 mM) so the comparison is apples-to-apples.
  3. For each sequence below, read OligoAnalyzer's Tm and put it in place of None.
  4. Re-run. The test then asserts agreement and reports MAE / max deviation.
Leave any entry as None to skip it. With the table empty the test SKIPS (passes) so it never
blocks the gate before you've supplied data.
"""
import sys
sys.path.insert(0, ".")
from oligoforge import thermo as T

OA_TOL_C = 2.0   # the README's stated bound; tighten to 1.5 once you've seen the real spread

# (sequence, oligoanalyzer_tm_C_or_None). Seeded with the locked-panel oligos; add any others.
# Replace each None with the value IDT OligoAnalyzer reports under the conditions above.
CALIB = [
    ("AGTCATTCTGATGTCGCTGATG",            None),   # IFNG F
    ("ACCTGTCAGTGTTTTCAAGCA",             None),   # IFNG R
    ("AACTTGCTCAGCCTGGTTTG",              None),   # IL4 F
    ("ATTCTTTAGTGAGGTGGTGCTG",            None),   # IL4 R
    ("TCGCTGGCATCAACAAGAAG",              None),   # RPL13 F
    ("TCGGGAAGAGGATGAGCTTG",              None),   # RPL13 R
    ("CCGTTACTTGGCTGAGGTTG",              None),   # YWHAZ F
    ("GATGGGATGTGTTGGTTGCA",              None),   # YWHAZ R
    ("AAAGGATTTGTGCTACCTTG",              None),   # Plasmodium R
    ("CTTACAAGATATCCACCACA",              None),   # Plasmodium P
    # ("...add more oligos with OligoAnalyzer Tms here...", None),
]

pairs = [(s, oa) for s, oa in CALIB if oa is not None]
if not pairs:
    print("  SKIP  Tm calibration not activated -- no OligoAnalyzer values supplied yet.")
    print("        See the docstring: paste IDT OligoAnalyzer Tms into CALIB to turn this on.")
    print("\nTm CALIBRATION SKIPPED (no data) -- not a failure.")
    sys.exit(0)

_fails = []
devs = {"owczarzy": [], "santalucia": []}
print("  conditions: DNA, 200 nM oligo, 50 mM Na+, 3 mM Mg++, 0.8 mM dNTP (must match OligoAnalyzer)")
print("  comparing both monovalent salt models against IDT (shipping = %s)\n" % T.SALT_METHOD)
for s, oa in pairs:
    t_ow = T._calc_tm(s, "owczarzy")
    t_sl = T._calc_tm(s, "santalucia")
    devs["owczarzy"].append(abs(t_ow - oa))
    devs["santalucia"].append(abs(t_sl - oa))
    ship = t_ow if T.SALT_METHOD == "owczarzy" else t_sl
    dev = abs(ship - oa)
    ok = dev <= OA_TOL_C
    print("  %s OA=%.1f  owczarzy=%.1f santalucia=%.1f  ship|dev|=%.2f %s" %
          (s[:18].ljust(18), oa, t_ow, t_sl, dev, "OK" if ok else "OVER"))
    if not ok:
        _fails.append(s)

mae_ow = sum(devs["owczarzy"]) / len(pairs)
mae_sl = sum(devs["santalucia"]) / len(pairs)
better = "santalucia" if mae_sl < mae_ow else "owczarzy"
print("\n  MAE vs IDT:  owczarzy=%.2f C   santalucia=%.2f C   (lower is better -> %s)" % (mae_ow, mae_sl, better))
if better != T.SALT_METHOD:
    print("  NOTE: '%s' tracks IDT better here. To switch, set thermo.SALT_METHOD='%s' AND update "
          "test_locked_panel.GOLDEN (the ordered Tms move with it)." % (better, better))
ship_devs = devs[T.SALT_METHOD]
print("  shipping (%s): n=%d  MAE=%.2f C  worst=%.2f C  (claim: within %.1f C)" %
      (T.SALT_METHOD, len(pairs), sum(ship_devs) / len(pairs), max(ship_devs), OA_TOL_C))
if _fails:
    print("Tm CALIBRATION FAILED for %d oligo(s) beyond %.1f C:" % (len(_fails), OA_TOL_C), ", ".join(_fails))
    sys.exit(1)
print("Tm CALIBRATION PASS -- OligoForge tracks OligoAnalyzer within %.1f C across %d oligos." % (OA_TOL_C, len(pairs)))
