import sys, os; sys.path.insert(0,".")
os.environ["OLIGOFORGE_EMAIL"]="fsj.qpcr.design@gmail.com"
from fastapi.testclient import TestClient
import app as A
c = TestClient(A.app)

def show(title, r):
    print(f"\n--- {title} -> HTTP {r.status_code} ---")
    j=r.json(); 
    import json; print(json.dumps(j, indent=1)[:700])

print("="*60,"\nLOCAL ENDPOINTS (no network)\n","="*60)
show("GET /api/profiles", c.get("/api/profiles"))
show("POST /api/qc (HMBS probe, IDT)", c.post("/api/qc", json={"seq":"ATCTTGTCCCCAGTTGTTGACATGGCC","role":"probe","profile":"idt_taqman"}))
show("POST /api/pair (HMBS F/R, amp 93)", c.post("/api/pair", json={"forward":"GAGCTATACCCCGACCTCTG","reverse":"CTTCTCTCCAATCTTGGAAAGCG","amplicon":93,"profile":"idt_taqman"}))
show("POST /api/matrix (3 HMBS oligos)", c.post("/api/matrix", json={"oligos":{"HMBS_F":"GAGCTATACCCCGACCTCTG","HMBS_R":"CTTCTCTCCAATCTTGGAAAGCG","HMBS_P":"ATCTTGTCCCCAGTTGTTGACATGGCC"}}))

HMBS_TMPL=("GGCCCGGATTCAGACTGATAGTGTAGTTATGATGCTCCGTGAGCTATACCCCGACCTCTGCTTTGAGATT"
           "GTGGCCATGTCAACAACTGGGGACAAGATCTTGGATACAGCGCTTTCCAAGATTGGAGAGAAGAGTCTCT"
           "TCACCAAAGAGTTGGAAAATGCACTTGAAAGAA")
show("POST /api/design (HMBS template)", c.post("/api/design", json={"template":HMBS_TMPL,"profile":"idt_taqman"}))

print("\n"+"="*60,"\nNETWORK ENDPOINTS (NCBI)\n","="*60)
show("POST /api/fetch (HMBS isoforms + common)", c.post("/api/fetch", json={"gene":"HMBS","organism":"Aphelocoma coerulescens","isoform_common":True}))
show("POST /api/intron (HMBS 223-315)", c.post("/api/intron", json={"gene":"HMBS","organism":"Aphelocoma coerulescens","mrna_acc":"XM_068994916.1","amp_start":223,"amp_end":315}))
print("\nALL ENDPOINTS RESPONDED.")
