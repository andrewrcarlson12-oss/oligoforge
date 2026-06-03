"""Regression pins for the isolate-validation engine (oligoforge.isolates.amplify).

Correctness (amplifies / 3'-block / orientation / probe / sanitization) plus the performance
cap that keeps a low-complexity template from blowing up the F x R product search. No network.
"""
import time, random
import sys; sys.path.insert(0, ".")
import oligoforge.isolates as ISO

F = "ACGTACGTACGTACGTACGT"; R = "TTTTGGGGCCCCAAAATTTT"
random.seed(1); bg = "".join(random.choice("ACGT") for _ in range(400))
tmpl = bg[:100] + F + bg[120:200] + ISO._rc_iupac(R) + bg[220:]

# clean amplicon: convergent F / rc(R), expected size, perfect identity, probe binds
c = ISO.amplify(F, R, "GTACGTACGTACGTACGTAC", tmpl)
assert c["amplifies"] and 110 <= c["product"] <= 130 and c["f_ident"] == 100.0 and c["r_ident"] == 100.0, c
# a 3'-terminal mismatch must stop extension -> no product
assert ISO.amplify(F[:-1] + ("C" if F[-1] != "C" else "A"), R, "", tmpl)["amplifies"] is False, "3' mismatch must block"
# divergent orientation (primers point away) -> no product
assert ISO.amplify(ISO._rc_iupac(R), ISO._rc_iupac(F), "", tmpl)["amplifies"] is False, "divergent must give no product"
# empty / N-only inputs are handled, never raise
assert ISO.amplify("", "", "", "")["amplifies"] is False
assert ISO.amplify(F, R, "", "N" * 3000)["amplifies"] is False
# probe on the minus strand still scores as binding
assert ISO.amplify(F, R, ISO._rc_iupac("GTACGTACGTACGTACGTAC"), tmpl)["probe_binds"] is True
# sanitization: regex metacharacters in a primer must neither crash nor corrupt matching
ISO.amplify("ACGT.()|+ACGTACGT", "TTTT$^GGGGCCCCAAAA", "", tmpl)

# performance: low-complexity templates give thousands of seed hits; the site cap must keep this bounded
homo = ("A" * 50 + "C" * 50) * 2000
t = time.time(); ISO.amplify("A" * 20, "T" * 20, "", homo); dt = time.time() - t
assert dt < 2.0, "pathological homopolymer not bounded: %.2fs" % dt
rep = ("ACGTACGTACGTACGTACGT" + "N" * 100) * 5000
t = time.time(); ISO.amplify("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", "", rep); dt = time.time() - t
assert dt < 2.0, "pathological repeat not bounded: %.2fs" % dt

print("ALL ISOLATE ASSERTS PASS")
