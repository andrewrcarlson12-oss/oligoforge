"""Standalone regression checks for the deterministic release inventory."""

import importlib.util
import json
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def touch(path, content="fixture\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")


def check(label, condition, detail=None):
    if not condition:
        raise AssertionError("%s: %r" % (label, detail))
    print("PASS", label)


manifest_tool = load_module(
    "oligoforge_release_manifest",
    ROOT / "tools" / "build_release_manifest.py",
)


with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "source"
    touch(source / "oligoforge" / "__init__.py", '__version__ = "1.37.0"\n')
    touch(source / "README.md", "release fixture\n")
    touch(source / "oligoforge" / "engine.py", "ENGINE = 1\n")
    touch(source / "tests" / "test_fixture.py", "pass\n")

    # Baseline build-source exclusions.
    touch(source / "dist" / "oligoforge.zip")
    touch(source / "build" / "generated.bin")
    touch(source / "__pycache__" / "cached.pyc")
    touch(source / "node_modules" / "dependency.js")
    touch(source / "panels" / "saved-panel.json")
    touch(source / "projects" / "saved-project.json")
    touch(source / ".env", "TOKEN=do-not-inventory\n")

    # Shared source-package safety exclusions.
    touch(source / ".env.production", "TOKEN=also-private\n")
    touch(source / "credentials.json", "{}\n")
    touch(source / "private.key", "private\n")
    touch(source / "htmlcov" / "index.html")
    touch(source / ".cache" / "generated.json")
    touch(source / "coverage.xml")

    manifest_path = manifest_tool.write_release_manifest(source)
    first_bytes = manifest_path.read_bytes()
    document = json.loads(first_bytes)
    paths = [record["path"] for record in document["files"]]

    check("manifest identifies the current release series",
          document["schema"] == "oligoforge-release-inventory/v1" and
          document["release"] == {
              "application": "OligoForge", "version": "1.37.0", "series": "1.37"
          }, document["release"])
    check("archive digest is delegated to an external sidecar",
          document["inventory_scope"]["archive_digest"] == "external_sidecar" and
          document["inventory_scope"]["archive_digest_sidecar"] ==
          "oligoforge-1.37.0-source.zip.sha256", document["inventory_scope"])
    check("manifest is excluded from its own inventory",
          "RELEASE_MANIFEST.json" not in paths, paths)
    check("source records are stable and ordered", paths == sorted(paths), paths)
    check("expected source files are inventoried",
          paths == ["README.md", "oligoforge/__init__.py", "oligoforge/engine.py",
                    "tests/test_fixture.py"], paths)
    check("summary fields match file records",
          document["file_count"] == len(document["files"]) and
          document["content_bytes"] == sum(item["bytes"] for item in document["files"]))

    manifest_tool.write_release_manifest(source)
    check("repeated generation is byte deterministic", manifest_path.read_bytes() == first_bytes)
    valid, issues = manifest_tool.verify_release_manifest(source)
    check("freshly generated manifest verifies", valid and not issues, issues)

    touch(source / "oligoforge" / "engine.py", "ENGINE = 2\n")
    valid, issues = manifest_tool.verify_release_manifest(source)
    check("verification detects source drift",
          not valid and any("changed files: oligoforge/engine.py" in issue for issue in issues),
          issues)

    try:
        manifest_tool.write_release_manifest(source, source / "generated" / "manifest.json")
    except ValueError:
        unsafe_output_rejected = True
    else:
        unsafe_output_rejected = False
    check("alternate in-tree output cannot create self-reference", unsafe_output_rejected)


with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "wrong-series"
    touch(source / "oligoforge" / "__init__.py", '__version__ = "1.38.0"\n')
    try:
        manifest_tool.build_release_manifest(source)
    except ValueError:
        wrong_series_rejected = True
    else:
        wrong_series_rejected = False
    check("1.37 schema rejects a different release series", wrong_series_rejected)


print("release manifest tests passed")
