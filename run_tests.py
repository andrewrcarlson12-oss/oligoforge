#!/usr/bin/env python3
"""Full source-level test runner / CI entry point for OligoForge.

Runs BOTH suites and exits non-zero if anything fails, so CI (or a pre-release check) fails loudly:
  * every tests/test_*.py  -- standalone assert scripts (engine, thermo, design, concurrency, fuzz,
    benchmark, goldens, locked panel). Each is a script that sys.exit(1)s on failure, so we run it
    as a subprocess and check the return code.
  * every tests/ui_*.js    -- Node UI harnesses (handler logic + a real-DOM jsdom integration audit),
    IFF a `node` binary is on PATH. Skipped-with-warning if node is absent (the Python suite is the
    hard gate; the JS harnesses are an additional guard).

Usage:  python run_tests.py            # full suite
        python run_tests.py --python   # Python suite only
        python run_tests.py --node     # Node harnesses only
Offline and deterministic: sets OLIGOFORGE_EMAIL so no test prompts for a contact address, and the
network-touching tests (test_intron, live epcr) are written to tolerate NCBI being unreachable.
"""
import os
import sys
import glob
import shutil
import subprocess
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ)
ENV.setdefault("OLIGOFORGE_EMAIL", "ci@example.com")
ENV["PYTHONPATH"] = ROOT + os.pathsep + ENV.get("PYTHONPATH", "")

GREEN, RED, DIM, RESET = "\033[32m", "\033[31m", "\033[2m", "\033[0m"


PY_TEST_TIMEOUT = int(os.environ.get("OLIGOFORGE_TEST_TIMEOUT", "300"))
NODE_TEST_TIMEOUT = int(os.environ.get("OLIGOFORGE_NODE_TEST_TIMEOUT", "120"))


def _run(cmd, timeout):
    try:
        p = subprocess.run(cmd, cwd=ROOT, env=ENV, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired as exc:
        out = (exc.stdout or "") + (exc.stderr or "")
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        return 124, out + "\nTIMEOUT after %d seconds" % timeout


def run_python():
    tests = sorted(glob.glob(os.path.join(ROOT, "tests", "test_*.py")))
    npass = nfail = 0
    failed = []
    for t in tests:
        name = os.path.basename(t)
        rc, out = _run([sys.executable, t], PY_TEST_TIMEOUT)
        if rc == 0:
            npass += 1
            print("  %sPASS%s %s" % (GREEN, RESET, name))
        else:
            nfail += 1
            failed.append(name)
            print("  %sFAIL%s %s" % (RED, RESET, name))
            print(DIM + "\n".join("      " + ln for ln in out.strip().splitlines()[-8:]) + RESET)
    print("PYTHON: %d passed / %d failed / %d total" % (npass, nfail, len(tests)))
    return nfail, failed


def run_node():
    node = shutil.which("node")
    harnesses = sorted(glob.glob(os.path.join(ROOT, "tests", "ui_*.js")))
    if not node:
        print("  %s(node not on PATH — skipping %d UI harnesses; Python suite is the hard gate)%s"
              % (DIM, len(harnesses), RESET))
        return 0, []
    npass = nfail = 0
    failed = []
    for h in harnesses:
        name = os.path.basename(h)
        rc, out = _run([node, h], NODE_TEST_TIMEOUT)
        if rc == 0:
            npass += 1
            print("  %sPASS%s %s" % (GREEN, RESET, name))
        else:
            nfail += 1
            failed.append(name)
            print("  %sFAIL%s %s" % (RED, RESET, name))
            print(DIM + "\n".join("      " + ln for ln in out.strip().splitlines()[-8:]) + RESET)
    print("NODE:   %d passed / %d failed / %d total" % (npass, nfail, len(harnesses)))
    return nfail, failed


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else ""
    total_fail = 0
    allfailed = []
    if only in ("", "--python"):
        print("== Python suite ==")
        f, names = run_python(); total_fail += f; allfailed += names
    if only in ("", "--node"):
        print("== Node UI harnesses ==")
        f, names = run_node(); total_fail += f; allfailed += names
    print()
    if total_fail:
        print("%sCI FAILED%s: %d suite(s) failing: %s" % (RED, RESET, total_fail, ", ".join(allfailed)))
        sys.exit(1)
    print("%sCI PASS%s: all suites green" % (GREEN, RESET))
