"""Deterministic intron-check contract; optional live NCBI smoke test.

The release gate is offline. Set ``OLIGOFORGE_LIVE_NCBI=1`` to additionally exercise the two
published accessions against NCBI.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import specificity as SP, ncbi

cases = [("HMBS", "Aphelocoma coerulescens", "XM_068994916.1", 223, 315),
         ("SDHA", "Aphelocoma coerulescens", "XM_069004197.1", 1637, 1766)]
fails=[]

def check_case(gene, org, acc, a, b):
    r=SP.intron_check(gene, org, a, b, mrna_acc=acc)
    ok=isinstance(r.get('verdict'), str) and bool(r['verdict'])
    print(('  PASS ' if ok else '  FAIL ') + gene + ' verdict on degraded path')
    if not ok: fails.append(gene)

if os.environ.get('OLIGOFORGE_LIVE_NCBI') == '1':
    ncbi.Entrez.email = SP.Entrez.email = os.environ.get('OLIGOFORGE_EMAIL', 'ci@example.com')
    for case in cases:
        check_case(*case)
else:
    original=SP.exon_junctions_mrna
    try:
        SP.exon_junctions_mrna=lambda *a, **k: (None, 'offline fixture: no exon structure', None)
        for case in cases:
            check_case(*case)
    finally:
        SP.exon_junctions_mrna=original

print(('INTRON TEST FAILED: '+', '.join(fails)) if fails else 'INTRON TEST PASS')
raise SystemExit(1 if fails else 0)
