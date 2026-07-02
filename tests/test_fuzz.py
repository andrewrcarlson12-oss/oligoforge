"""Input-validation + error-contract fuzz (v1.27.2).

Every sequence/parameter endpoint must degrade cleanly on hostile input: a clean {"error": ...}
(HTTP 200) or a validation 4xx -- NEVER a 500 or an unhandled exception. Also pins the three
concrete gaps closed this release: oligo/template length caps, non-finite salt rejection, and the
dNTP>=Mg free-Mg warning. In-process TestClient, offline, deterministic.
Run standalone: `python tests/test_fuzz.py`.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("OLIGOFORGE_EMAIL", "test@example.com")
from fastapi.testclient import TestClient
import app as A
from oligoforge import thermo as T

cl = TestClient(A.app)
cl.post("/api/conditions", json={"mv_conc": 50, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200, "anneal_c": 60})
fails = []
def check(name, ok, detail=""):
    print(("  PASS " if ok else "  FAIL ") + name + ("" if ok else "  -> " + str(detail)))
    if not ok:
        fails.append(name)

# 1. no endpoint 500s on hostile sequence input
BAD = ["", " ", "\n\t ", "N" * 400, "ACGTX", "12345", ">only header", "U" * 50, "ACGT" * 100000,
       "RYSWKM", "A", "-" * 30, "ACGT\x00GT"]
fivexx = []
for s in BAD:
    for ep, payload in (("/api/qc", {"seq": s, "role": "primer", "profile": "idt_taqman"}),
                        ("/api/design", {"template": s, "profile": "idt_taqman"}),
                        ("/api/pair", {"forward": s, "reverse": s, "profile": "idt_taqman"}),
                        ("/api/matrix", {"oligos": {"a": s, "b": s}})):
        try:
            r = cl.post(ep, json=payload)
            if r.status_code >= 500:
                fivexx.append((ep, repr(s)[:24], r.status_code))
        except Exception as e:
            fivexx.append((ep, repr(s)[:24], "EXC " + type(e).__name__))
check("no 5xx across %d bad seqs x 4 endpoints" % len(BAD), not fivexx, fivexx[:3])

# 2. oligo length cap
big = cl.post("/api/qc", json={"seq": "ACGT" * 1000, "role": "primer", "profile": "idt_taqman"}).json()
check("QC rejects over-long oligo with clean error", bool(big.get("error")) and "too long" in big["error"], big.get("error"))
ok = cl.post("/api/qc", json={"seq": "AGTCATTCTGATGTCGCTGATG", "role": "primer", "profile": "idt_taqman"}).json()
check("QC still accepts a normal 22-mer", ok.get("tm") and not ok.get("error"), ok.get("error"))

# 3. template length cap
bigt = cl.post("/api/design", json={"template": "ACGT" * 20000, "profile": "idt_taqman"}).json()
check("design rejects over-long template with clean error", bool(bigt.get("error")) and "too long" in bigt["error"], bigt.get("error"))

# 4. non-finite salt rejected (raw JSON tokens, the way a real client can send them), no 500
for label, body in (("Infinity", '{"mv_conc": Infinity, "dv_conc":3, "dntp_conc":0.8, "dna_conc":200, "anneal_c":60}'),
                    ("NaN", '{"dv_conc": NaN, "mv_conc":50, "dntp_conc":0.8, "dna_conc":200, "anneal_c":60}')):
    r = cl.post("/api/conditions", content=body, headers={"content-type": "application/json"})
    j = r.json() if r.status_code == 200 else {}
    check("non-finite salt (%s) -> clean error, no 500" % label,
          r.status_code == 200 and "finite" in (j.get("error") or ""), (r.status_code, r.text[:60]))

# 5. dNTP >= Mg -> non-fatal warning (free Mg ~0), and normal conditions produce none
warn = cl.post("/api/conditions", json={"mv_conc": 50, "dv_conc": 1.0, "dntp_conc": 5.0, "dna_conc": 200, "anneal_c": 60}).json()
check("dNTP>=Mg returns a free-Mg warning (not silent)", bool(warn.get("warning")) and "free Mg2+" in warn["warning"], warn)
normal = cl.post("/api/conditions", json={"mv_conc": 50, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200, "anneal_c": 60}).json()
check("normal conditions produce no spurious warning", not normal.get("warning"), normal.get("warning"))

# 6. out-of-range salt still rejected cleanly (regression on the pre-existing validation)
for c, why in (({"mv_conc": -5}, "negative mv"), ({"dv_conc": 99999}, "huge Mg"),
               ({"anneal_c": 200}, "anneal too high"), ({"mv_conc": 0, "dv_conc": 0}, "no salt")):
    base = {"mv_conc": 50, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200, "anneal_c": 60}; base.update(c)
    j = cl.post("/api/conditions", json=base).json()
    check("out-of-range rejected: %s" % why, bool(j.get("error")), j)
cl.post("/api/conditions", json={"mv_conc": 50, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200, "anneal_c": 60})

if fails:
    print("FUZZ FAILURES:", fails); sys.exit(1)
print("ALL FUZZ / VALIDATION ASSERTS PASS")
