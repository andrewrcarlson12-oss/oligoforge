"""NCBI reliability-layer gate (cache + retry). Offline -- all network calls are monkeypatched.
Run from repo root:  python tests/test_ncbi_cache.py
"""
import os, sys, tempfile, urllib.error
sys.path.insert(0, ".")
from oligoforge import ncbi

_fails = []
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        _fails.append(name)


class _FakeHandle:
    def __init__(self, text): self._t = text
    def read(self): return self._t
    def close(self): pass


# isolate the cache in a fresh temp dir, caching ON
ncbi._CACHE_DIR = tempfile.mkdtemp(prefix="of_cache_test_")
ncbi._CACHE_ON = True
ncbi._CACHE_TTL = 3600

# ---- 1) cache hit: a second identical efetch does NOT touch the network ----
calls = {"efetch": 0}
def fake_efetch(**k):
    calls["efetch"] += 1
    return _FakeHandle(">acc\nACGTACGTACGT\n")
ncbi.Entrez.efetch = fake_efetch

r1 = ncbi._efetch_text(db="nucleotide", id="ZZ_TEST_1", rettype="fasta", retmode="text")
r2 = ncbi._efetch_text(db="nucleotide", id="ZZ_TEST_1", rettype="fasta", retmode="text")
check("efetch text returned", ">acc" in r1 and r1 == r2, r1[:12])
check("second identical efetch served from cache (network hit once)", calls["efetch"] == 1, calls["efetch"])

# a DIFFERENT id is a cache miss -> one more network call
ncbi._efetch_text(db="nucleotide", id="ZZ_TEST_2", rettype="fasta", retmode="text")
check("different id is a cache miss (network hit twice total)", calls["efetch"] == 2, calls["efetch"])

# caching OFF -> always hits the network
ncbi._CACHE_ON = False
ncbi._efetch_text(db="nucleotide", id="ZZ_TEST_1", rettype="fasta", retmode="text")
check("cache disabled re-fetches even for a known id", calls["efetch"] == 3, calls["efetch"])
ncbi._CACHE_ON = True

# ---- 2) retry: a transient URLError is retried, then succeeds ----
state = {"n": 0}
def flaky_esearch(**k):
    state["n"] += 1
    if state["n"] < 3:
        raise urllib.error.URLError("temporary failure in name resolution")
    return _FakeHandle("OK")
# speed: no real backoff sleep
import time as _t
_real_sleep = _t.sleep
ncbi.time.sleep = lambda *_a, **_k: None
try:
    h = ncbi._net(flaky_esearch, db="nucleotide", term="x")
    check("transient failure retried to success", h.read() == "OK" and state["n"] == 3, state["n"])
except Exception as e:
    check("transient failure retried to success", False, "raised %r" % e)

# ---- 3) a logical HTTP 404 is NOT retried (fails fast) ----
state2 = {"n": 0}
def hard_404(**k):
    state2["n"] += 1
    raise urllib.error.HTTPError("http://x", 404, "Not Found", {}, None)
raised = False
try:
    ncbi._net(hard_404, db="nucleotide", term="x")
except urllib.error.HTTPError:
    raised = True
check("HTTP 404 raises immediately (no retry)", raised and state2["n"] == 1, state2["n"])

# ---- 4) a retryable HTTP 503 IS retried ----
state3 = {"n": 0}
def flaky_503(**k):
    state3["n"] += 1
    if state3["n"] < 2:
        raise urllib.error.HTTPError("http://x", 503, "Service Unavailable", {}, None)
    return _FakeHandle("OK503")
try:
    h = ncbi._net(flaky_503, db="nucleotide", term="x")
    check("HTTP 503 retried to success", h.read() == "OK503" and state3["n"] == 2, state3["n"])
except Exception as e:
    check("HTTP 503 retried to success", False, "raised %r" % e)

ncbi.time.sleep = _real_sleep
print("")
if _fails:
    print("NCBI CACHE/RETRY GATE FAILED:", ", ".join(_fails)); sys.exit(1)
print("ALL NCBI CACHE/RETRY ASSERTS PASS")
