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

# --- LNA / IDT-modified oligos must be scored by their BARE bases, not the '+' / mod notation
# (the '+' previously became 'N', frame-shifting the probe into a spurious uniform "weak" identity) ---
import oligoforge.thermo as T
CORE = "GTACGTACGTACGTACGTAC"                                  # 20-mer probe, embedded once below
assert T.strip_mods("G+TAC*GT AC") == "GTACGTAC"
assert T.strip_mods("/56-FAM/GTAC+GT+ACGTACGTAC+GTAC/3IABkFQ/") == CORE
random.seed(7)
_lf = "".join(random.choice("ACGT") for _ in range(82)); _lr = "".join(random.choice("ACGT") for _ in range(82))
tlna = F + _lf + CORE + _lr + ISO._rc_iupac(R)                 # F ... CORE ... rc(R): probe present at 100%
_bare = ISO.amplify(F, R, CORE, tlna)
_plus = ISO.amplify(F, R, "GTAC+GT+ACGTACGTAC+GTAC", tlna)
_ordr = ISO.amplify(F, R, "/56-FAM/GTAC+GT+ACGTACGTAC+GTAC/3IABkFQ/", tlna)
assert _bare["amplifies"] and _bare["probe_ident"] == 100.0 and _bare["probe_binds"], _bare
assert _plus["probe_ident"] == _bare["probe_ident"] and _plus["probe_binds"], ("LNA '+' probe mis-scored", _plus)
assert _ordr["probe_ident"] == _bare["probe_ident"] and _ordr["probe_binds"], ("IDT-order probe mis-scored", _ordr)
# the fix doesn't blanket-pass: a genuinely mismatched probe is still flagged weak
_bad = ISO.amplify(F, R, CORE[:6] + "TTTTTT" + CORE[12:], tlna)
assert _bad["probe_ident"] < 100.0 and not _bad["probe_binds"], ("real mismatch should stay weak", _bad)

# --- 'amplifies' uses a 3'-terminal extension gate (last ~6 nt exact + <=max_mm total), not a long
#     exact seed; a no-product result now reports REAL best homology + a reason (absent vs 3' mismatch
#     vs size window), never a misleading 0% for a primer whose region is present but won't prime ---
import random as _rnd
_rnd.seed(101); _r=lambda n:"".join(_rnd.choice("ACGT") for _ in range(n))
gF=_r(20); gR=_r(20); _m={"A":"C","C":"A","G":"T","T":"G"}
_tw=lambda Fv: _r(15)+Fv+_r(60)+ISO._rc_iupac(gR)+_r(15)
_mk=lambda d: gF[:len(gF)-1-d]+_m[gF[len(gF)-1-d]]+gF[len(gF)-d:]    # single mismatch d nt from the 3' end
_far=ISO.amplify(gF,gR,"",_tw(_mk(6)))                              # 6 nt from 3' -> tolerated
assert _far["amplifies"] and _far["f_ident"]==95.0, ("mismatch 6 nt from 3' should amplify", _far)
_near=ISO.amplify(gF,gR,"",_tw(_mk(3)))                             # 3 nt from 3' -> blocks priming
assert (not _near["amplifies"]) and _near["f_ident"]==95.0 and "mismatch" in (_near.get("reason") or ""), ("3' mismatch", _near)
_abs=ISO.amplify(gF,gR,"",_r(120))                                  # region absent
assert (not _abs["amplifies"]) and _abs["f_ident"]==0.0 and "absent" in (_abs.get("reason") or ""), ("absent", _abs)
_big=ISO.amplify(gF,gR,"",_r(15)+gF+_r(4000)+ISO._rc_iupac(gR)+_r(15),max_product=2000)
assert (not _big["amplifies"]) and "size window" in (_big.get("reason") or ""), ("size window", _big)

print("ALL ISOLATE ASSERTS PASS")
