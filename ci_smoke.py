"""CI smoke test for the frozen OligoForge binary.

Launches dist/OligoForge(.exe), waits for the local server, checks that the
cockpit is served and that primer3 thermodynamics work inside the bundle, then
shuts it down. Exits non-zero on any failure so the build fails loudly.
"""
import contextlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
import urllib.request

PORT = 8111
exe = os.path.join("dist", "OligoForge" + (".exe" if platform.system() == "Windows" else ""))
if not os.path.exists(exe):
    sys.exit("frozen binary not found: " + exe)


def listening():
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        return s.connect_ex(("127.0.0.1", PORT)) == 0


proc = subprocess.Popen([exe], env=dict(os.environ, OLIGOFORGE_EMAIL="ci@example.org"))
try:
    for _ in range(160):                       # up to ~40 s for first-run unpack
        if listening():
            break
        time.sleep(0.25)
    else:
        raise SystemExit("server never started listening on %d" % PORT)

    home = urllib.request.urlopen("http://127.0.0.1:%d/" % PORT, timeout=15)
    body = home.read()
    assert home.status == 200 and b"OligoForge" in body, "home page not served"

    req = urllib.request.Request(
        "http://127.0.0.1:%d/api/qc" % PORT,
        data=json.dumps({"seq": "GAGCTATACCCCGACCTCTG", "role": "primer",
                         "profile": "idt_taqman"}).encode(),
        headers={"Content-Type": "application/json"},
    )
    d = json.loads(urllib.request.urlopen(req, timeout=15).read())
    assert abs(d["tm"] - 61.8) < 1.5, "primer3 Tm wrong in frozen binary: %r" % d
    assert d.get("hairpin_dg") is not None, "primer3 hairpin missing in frozen binary"
    print("SMOKE OK on %s: served + primer3 thermo working (tm=%.2f)" % (platform.system(), d["tm"]))
finally:
    proc.terminate()
    try:
        proc.wait(timeout=10)
    except Exception:
        proc.kill()
