"""Exon-junction / gDNA-exclusion design: '|' markers in the template -> per-candidate junction
classification + gDNA-exclusion re-ranking. Offline + deterministic.
Run: OLIGOFORGE_EMAIL=you@x python3 tests/test_junction.py   (exit 0 = pass)."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import app as A
from fastapi.testclient import TestClient
c = TestClient(A.app)
fails = []

def check(name, cond):
    if not cond:
        fails.append(name)
    print(("PASS " if cond else "FAIL ") + name)

# --- parser unit checks (clean_seq is char-wise, so positions line up with the cleaned ref) ---
ref, pos, _ = A._parse_junction_template("AAAACCCC|GGGGTTTT")
check("parser: ref strips marker", ref == "AAAACCCCGGGGTTTT")
check("parser: one boundary at 8", pos == [8])
ref2, pos2, _ = A._parse_junction_template("AA AC CC C|GG GG TT TT")          # whitespace tolerated, alignment holds
check("parser: whitespace tolerated", ref2 == "AAACCCCGGGGTTTT" and pos2 == [7])
_, pos3, _ = A._parse_junction_template("|ABCDEF")                            # leading marker -> no boundary at 0
check("parser: leading marker ignored", pos3 == [])
_, pos4, _ = A._parse_junction_template("ABCDEF|")                            # trailing marker -> no boundary at len
check("parser: trailing marker ignored", pos4 == [])
_, pos5, _ = A._parse_junction_template("AAA|CCC|GGG")
check("parser: two boundaries", pos5 == [3, 6])

# --- endpoint: classification follows where the marker falls relative to the chosen oligos ---
random.seed(21)
tmpl = "".join(random.choice("ACGT") for _ in range(520))
base = c.post("/api/design", json={"template": tmpl, "profile": "idt_taqman"}).json()
check("no marker -> junction None", base["n_junctions"] == 0 and base["candidates"][0]["junction"] is None)
c0 = base["candidates"][0]
f, r, p, amp = c0["f_span"], c0["r_span"], c0["probe_span"], c0["amp_span"]

def design_marked(at):
    t = tmpl[:at] + "|" + tmpl[at:]
    d = c.post("/api/design", json={"template": t, "profile": "idt_taqman"}).json()
    same = [x for x in d["candidates"] if x["forward"] == c0["forward"]]
    return d, (same[0] if same else None)

# (a) inside the probe -> strong, and a gDNA-excluding candidate is ranked first
da, ca = design_marked((p[0] + p[1]) // 2)
check("probe-spanning -> strong", ca and ca["junction"]["level"] == "strong" and ca["junction"]["probe"])
check("gDNA-excluding candidate ranked #1", da["candidates"][0]["junction"]["level"] == "strong")
check("n_junctions surfaced", da["n_junctions"] == 1 and da["junctions"] == [(p[0] + p[1]) // 2])

# (b) inter-oligo gap (inside amplicon, no oligo crosses) -> size
gap = (f[1] + p[0]) // 2 if (p[0] - f[1]) >= 2 else (p[1] + r[0]) // 2
_, cb = design_marked(gap)
check("inter-oligo gap -> size", cb and cb["junction"]["level"] == "size"
      and cb["junction"]["amplicon"] and not (cb["junction"]["probe"] or cb["junction"]["forward"] or cb["junction"]["reverse"]))

# (c) outside the amplicon -> none
outside = max(1, f[0] // 2) if f[0] > 3 else min(len(tmpl) - 1, r[1] + 10)
_, cc = design_marked(outside)
check("outside amplicon -> none", cc and cc["junction"]["level"] == "none" and not cc["junction"]["amplicon"])

if fails:
    print("\nFAILED:", ", ".join(fails)); sys.exit(1)
print("\nALL JUNCTION ASSERTS PASS")
