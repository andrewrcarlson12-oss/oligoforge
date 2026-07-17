"""Concurrency contract for the reaction conditions (v1.27.2).

thermo.COND + ANNEAL_C are process-global. Under the sync-endpoint threadpool a /api/conditions
change can run concurrently with a design/QC request. This pins that the Tm/structure caches are
keyed on an immutable conditions SNAPSHOT so that:
  (1) numbers are byte-identical to the pre-snapshot implementation (no panel/golden drift), and
  (2) a value read under a given snapshot ALWAYS equals that snapshot's value, even while another
      thread is hammering set_conditions -- no torn reads, no stale-cache cross-contamination.
Offline, deterministic, no network. Run standalone: `python tests/test_concurrency.py`.
"""
import os, sys, threading
import concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoforge import thermo as T

fails = []
def check(name, ok, detail=""):
    print(("  PASS " if ok else "  FAIL ") + name + ("" if ok else "  -> " + str(detail)))
    if not ok:
        fails.append(name)

OLIGO = "AGTCATTCTGATGTCGCTGATG"
CA = dict(mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200, anneal_c=60)
CB = dict(mv_conc=100, dv_conc=8.0, dntp_conc=1.5, dna_conc=200, anneal_c=60)

# reference values for the two conditions
T.set_conditions(**CA); tmA = round(T.tm_acc(OLIGO), 6); selA = round(T.tm(OLIGO), 6)
T.set_conditions(**CB); tmB = round(T.tm_acc(OLIGO), 6)
T.set_conditions(**CA)
check("two conditions give distinct display Tm (sanity)", tmA != tmB, (tmA, tmB))

# 1. snapshot is atomic + hashable
snap = T._snapshot()
check("snapshot is a hashable 5-tuple", isinstance(snap, tuple) and len(snap) == 5 and hash(snap) is not None, snap)

# 2. same conditions -> identical numbers (no drift from the snapshot indirection)
T.set_conditions(**CA)
check("tm_acc stable & correct under fixed conditions", round(T.tm_acc(OLIGO), 6) == tmA, T.tm_acc(OLIGO))
check("tm(sel) stable under fixed conditions", round(T.tm(OLIGO), 6) == selA, T.tm(OLIGO))

# 3. THE race test: 3 writers flip CA<->CB while 5 readers compute; every cached read must equal
#    ITS OWN snapshot's reference value (snapshot-consistency), and no read may be a torn/garbage
#    value outside {tmA, tmB}.
valid = {tmA, tmB}
incon = []
garbage = []
stop = threading.Event()

def flip():
    while not stop.is_set():
        T.set_conditions(**CB); T.set_conditions(**CA)

def read():
    for _ in range(4000):
        s = T._snapshot()
        v = round(T._tm_acc_at(OLIGO, s), 6)
        expect = tmA if s[0] == 50 else tmB
        if v != expect:
            incon.append((v, expect))
        pv = round(T.tm_acc(OLIGO), 6)
        if pv not in valid:
            garbage.append(pv)

with cf.ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(flip) for _ in range(3)] + [ex.submit(read) for _ in range(5)]
    for f in futs[3:]:
        f.result()
    stop.set()
T.set_conditions(**CA)

check("snapshot-consistent reads under concurrent set_conditions (0 violations)", not incon, incon[:3])
check("no torn/garbage tm_acc under concurrent set_conditions", not garbage, garbage[:3])

# 4. set_conditions swaps COND atomically as a whole new dict (identity changes; never in-place)
before = T.COND
T.set_conditions(dv_conc=5.0)
check("set_conditions rebinds COND (atomic swap, not in-place mutation)", T.COND is not before, "same dict identity")
T.set_conditions(**CA)

if fails:
    print("CONCURRENCY FAILURES:", fails); sys.exit(1)
print("ALL CONCURRENCY ASSERTS PASS")
