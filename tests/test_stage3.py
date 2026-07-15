"""Stage-3 API smoke checks.

Offline by default. Set OLIGOFORGE_LIVE_NCBI=1 to include live NCBI routes.
"""
import os
import sys
sys.path.insert(0, ".")
os.environ.setdefault("OLIGOFORGE_EMAIL", "ci@example.com")

from fastapi.testclient import TestClient
import app as A

c = TestClient(A.app)

def require_ok(title, response):
    print(f"  {title}: HTTP {response.status_code}")
    assert response.status_code == 200, response.text[:500]
    return response.json()

require_ok("profiles", c.get("/api/profiles"))
require_ok("qc", c.post("/api/qc", json={
    "seq": "ATCTTGTCCCCAGTTGTTGACATGGCC", "role": "probe", "profile": "idt_taqman"}))
require_ok("pair", c.post("/api/pair", json={
    "forward": "GAGCTATACCCCGACCTCTG", "reverse": "CTTCTCTCCAATCTTGGAAAGCG",
    "amplicon": 93, "profile": "idt_taqman"}))
require_ok("matrix", c.post("/api/matrix", json={"oligos": {
    "HMBS_F": "GAGCTATACCCCGACCTCTG",
    "HMBS_R": "CTTCTCTCCAATCTTGGAAAGCG",
    "HMBS_P": "ATCTTGTCCCCAGTTGTTGACATGGCC"}}))

template = ("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATT"
            "GTGGCCATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCT"
            "TCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
design = require_ok("design", c.post("/api/design", json={
    "template": template, "profile": "idt_taqman"}))
assert design.get("candidates") and design["candidates"][0].get("forward") and design["candidates"][0].get("reverse")

if os.environ.get("OLIGOFORGE_LIVE_NCBI") == "1":
    require_ok("fetch", c.post("/api/fetch", json={
        "gene": "HMBS", "organism": "Aphelocoma coerulescens", "isoform_common": True}))
    require_ok("intron", c.post("/api/intron", json={
        "gene": "HMBS", "organism": "Aphelocoma coerulescens",
        "mrna_acc": "XM_068994916.1", "amp_start": 223, "amp_end": 315}))
else:
    print("  SKIP live NCBI routes (set OLIGOFORGE_LIVE_NCBI=1 to enable).")

print("ALL STAGE-3 ASSERTS PASS")
