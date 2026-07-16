#!/usr/bin/env python3
"""Full source-level test runner / CI entry point for OligoForge.

Runs BOTH suites and exits non-zero if anything fails, so CI (or a pre-release check) fails loudly:
  * every test_*.py recursively under tests/ -- standalone assert scripts (engine, thermo, design,
    concurrency, fuzz, benchmark, goldens, locked panel). Each is a script that sys.exit(1)s on
    failure, so we run it as a subprocess and check the return code.
  * every ui_*.js recursively under tests/ -- Node UI harnesses (handler logic + a real-DOM jsdom
    integration audit), IFF a `node` binary is on PATH. Skipped-with-warning if node is absent (the
    Python suite is the hard gate; the JS harnesses are an additional guard).

Test discovery prunes dependency and test-cache directories so generated cache contents cannot
silently become release tests.

Usage:  python run_tests.py            # full suite
        python run_tests.py --python   # Python suite only
        python run_tests.py --node     # Node harnesses only
Offline and deterministic: sets OLIGOFORGE_EMAIL so no test prompts for a contact address, and the
network-touching tests (test_intron, live epcr) are written to tolerate NCBI being unreachable.
"""
import fnmatch
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ)
ENV.setdefault("OLIGOFORGE_EMAIL", "ci@example.com")
ENV["PYTHONPATH"] = ROOT + os.pathsep + ENV.get("PYTHONPATH", "")

GREEN, RED, DIM, RESET = "\033[32m", "\033[31m", "\033[2m", "\033[0m"


PY_TEST_TIMEOUT = int(os.environ.get("OLIGOFORGE_TEST_TIMEOUT", "600"))
NODE_TEST_TIMEOUT = int(os.environ.get("OLIGOFORGE_NODE_TEST_TIMEOUT", "120"))

TEST_DISCOVERY_EXCLUDED_DIRS = frozenset({
    "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".hypothesis",
    ".cache", "cache", ".tox", ".nox", ".nyc_output", "node_modules",
})
TEST_DISCOVERY_EXCLUDED_DIR_PATTERNS = ("*.egg-info",)


def _is_excluded_test_directory(name):
    return (name in TEST_DISCOVERY_EXCLUDED_DIRS or
            any(fnmatch.fnmatchcase(name, pattern)
                for pattern in TEST_DISCOVERY_EXCLUDED_DIR_PATTERNS))


def discover_tests(filename_pattern, root=ROOT):
    """Return deterministic absolute paths matching *filename_pattern* below tests/."""
    tests_root = os.path.join(root, "tests")
    matches = []
    for current, dirs, files in os.walk(tests_root):
        dirs[:] = sorted(d for d in dirs if not _is_excluded_test_directory(d))
        matches.extend(
            os.path.join(current, name)
            for name in sorted(files)
            if fnmatch.fnmatchcase(name, filename_pattern)
        )
    return sorted(matches)


def _test_name(path):
    return os.path.relpath(path, ROOT).replace(os.sep, "/")


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
    tests = discover_tests("test_*.py")
    npass = nfail = 0
    failed = []
    for t in tests:
        name = _test_name(t)
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
    harnesses = discover_tests("ui_*.js")
    if not node:
        print("  %s(node not on PATH — skipping %d UI harnesses; Python suite is the hard gate)%s"
              % (DIM, len(harnesses), RESET))
        return 0, []
    npass = nfail = 0
    failed = []
    for h in harnesses:
        name = _test_name(h)
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
    # On POSIX, replace this process with the shell runner so exhaustive Primer3
    # tests are direct children of the shell. CPython 3.13 can deadlock primer3-py
    # when the same workload runs in a nested Python subprocess. Windows retains
    # the pure-Python fallback below.
    shell_runner = os.path.join(ROOT, "run_tests.sh")
    if os.name != "nt" and os.path.exists(shell_runner) and os.environ.get("OLIGOFORGE_PY_RUNNER") != "1":
        os.execv("/bin/bash", ["bash", shell_runner] + sys.argv[1:])
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
