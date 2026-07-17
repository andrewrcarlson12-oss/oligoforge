"""Focused tests for recursive suite discovery and the directory file-count release gate."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load_module("oligoforge_test_runner", ROOT / "run_tests.py")
FILE_GATE = load_module("oligoforge_file_count_gate", ROOT / "tools" / "check_directory_file_counts.py")


def touch(path, content=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def check(label, condition, detail=None):
    if not condition:
        raise AssertionError("%s: %r" % (label, detail))
    print("PASS", label)


with tempfile.TemporaryDirectory() as tmp:
    repo = Path(tmp)
    touch(repo / "tests" / "test_root.py")
    touch(repo / "tests" / "unit" / "test_nested.py")
    touch(repo / "tests" / "unit" / "ui_nested.js")
    touch(repo / "tests" / "ui_root.js")
    for cache_name in RUNNER.TEST_DISCOVERY_EXCLUDED_DIRS:
        touch(repo / "tests" / cache_name / "test_cached.py")
        touch(repo / "tests" / cache_name / "ui_cached.js")
    touch(repo / "tests" / "fixture.egg-info" / "test_cached.py")
    touch(repo / "tests" / "fixture.egg-info" / "ui_cached.js")

    python_tests = [Path(path).relative_to(repo).as_posix()
                    for path in RUNNER.discover_tests("test_*.py", str(repo))]
    node_tests = [Path(path).relative_to(repo).as_posix()
                  for path in RUNNER.discover_tests("ui_*.js", str(repo))]
    check("Python discovery is recursive and prunes caches",
          python_tests == ["tests/test_root.py", "tests/unit/test_nested.py"], python_tests)
    check("Node discovery is recursive and prunes caches",
          node_tests == ["tests/ui_root.js", "tests/unit/ui_nested.js"], node_tests)


with tempfile.TemporaryDirectory() as tmp:
    repo = Path(tmp)
    runner = repo / "run_tests.sh"
    shutil.copyfile(ROOT / "run_tests.sh", runner)
    touch(repo / "tests" / "test_root.py")
    touch(repo / "tests" / "unit" / "test_nested.py")
    touch(repo / "tests" / "__pycache__" / "test_cached.py", "raise SystemExit(9)\n")
    touch(repo / "tests" / "ui_root.js")
    touch(repo / "tests" / "browser" / "ui_nested.js")
    touch(repo / "tests" / ".cache" / "ui_cached.js")
    env = dict(os.environ)
    env["PYTHON"] = sys.executable
    result = subprocess.run(["bash", str(runner), "--python"], cwd=repo, env=env,
                            capture_output=True, text=True, timeout=30)
    check("shell runner recursively discovers Python tests", result.returncode == 0, result.stdout + result.stderr)
    check("shell runner excludes cached Python tests",
          "PYTHON: 2 passed / 0 failed / 2 total" in result.stdout, result.stdout)

    fake_bin = repo / "bin"
    fake_node = fake_bin / "node"
    touch(fake_node, "#!/bin/sh\nexit 0\n")
    fake_node.chmod(0o755)
    env["PATH"] = str(fake_bin) + os.pathsep + env.get("PATH", "")
    result = subprocess.run(["bash", str(runner), "--node"], cwd=repo, env=env,
                            capture_output=True, text=True, timeout=30)
    check("shell runner recursively discovers Node harnesses", result.returncode == 0,
          result.stdout + result.stderr)
    check("shell runner excludes cached Node harnesses",
          "NODE:   2 passed / 0 failed / 2 total" in result.stdout, result.stdout)


with tempfile.TemporaryDirectory() as tmp:
    repo = Path(tmp)
    below = repo / "below"
    at_limit = repo / "source" / "at_limit"
    excluded_cache = repo / ".pytest_cache" / "bulk"
    excluded_dependency = repo / "node_modules" / "bulk"
    excluded_packaging = repo / "package.egg-info" / "bulk"
    excluded_build = repo / "build" / "bulk"
    for number in range(FILE_GATE.DEFAULT_LIMIT - 1):
        touch(below / ("file_%03d" % number))
    for directory in (at_limit, excluded_cache, excluded_dependency, excluded_packaging,
                      excluded_build):
        for number in range(FILE_GATE.DEFAULT_LIMIT):
            touch(directory / ("file_%03d" % number))

    overfull = FILE_GATE.find_overfull_directories(repo)
    relative = [(path.relative_to(repo).as_posix(), count) for path, count in overfull]
    check("file-count gate uses an inclusive 100-file threshold",
          relative == [("source/at_limit", FILE_GATE.DEFAULT_LIMIT)], relative)
    result = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "check_directory_file_counts.py"), str(repo)],
        capture_output=True, text=True, timeout=30)
    check("file-count gate command fails and identifies an overfull directory",
          result.returncode == 1 and "source/at_limit: 100 direct child files" in result.stdout,
          result.stdout + result.stderr)
current_overfull = FILE_GATE.find_overfull_directories(ROOT)
check("current repository passes the file-count release gate",
      current_overfull == [], current_overfull)

workflow = (ROOT / ".github" / "workflows" / "build.yml").read_text(encoding="utf-8")
check("CI verifies the committed source inventory before testing and packaging",
      workflow.count("python tools/build_release_manifest.py --verify") >= 2,
      workflow.count("python tools/build_release_manifest.py --verify"))

print("release engineering tests passed")
