"""Release-identity and reproducible-source-package regression checks."""

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(label, condition, detail=None):
    if not condition:
        raise AssertionError("%s: %r" % (label, detail))
    print("PASS", label)


from oligoforge import __version__
import launcher

builder = load_module("oligoforge_source_release", ROOT / "tools" / "build_source_release.py")
verifier = load_module("oligoforge_checksum_verifier", ROOT / "tools" / "verify_release_checksums.py")
writer = load_module("oligoforge_checksum_writer", ROOT / "tools" / "write_release_checksum.py")

package = json.loads((ROOT / "package.json").read_text(encoding="utf-8"))
lock = json.loads((ROOT / "package-lock.json").read_text(encoding="utf-8"))

check("Python package has the intended release identity", __version__ == "1.37.0", __version__)
check("desktop launcher derives its identity from the Python package",
      launcher.APP_VERSION == "v" + __version__, launcher.APP_VERSION)
check("package.json identity matches Python", package["version"] == __version__, package["version"])
check("package-lock root identity matches Python", lock["version"] == __version__, lock["version"])
check("package-lock package identity matches Python",
      lock["packages"][""]["version"] == __version__, lock["packages"][""]["version"])
check("README publishes the current identity",
      (ROOT / "README.md").read_text(encoding="utf-8").startswith("# OligoForge " + __version__))
check("release summary publishes the current identity",
      (ROOT / "RELEASE_SUMMARY.md").read_text(encoding="utf-8").startswith(
          "# OligoForge " + __version__ + " release summary"))
check("changelog leads with the current identity",
      ("## " + __version__) in (ROOT / "CHANGELOG.md").read_text(encoding="utf-8").splitlines()[2])
check("lifecycle browser fixture reports the current health version",
      'version:"%s"' % __version__ in
      (ROOT / "tests" / "ui_lifecycle.js").read_text(encoding="utf-8"))
ui = (ROOT / "static" / "index.html").read_text(encoding="utf-8")
check("browser header and footer publish the current identity",
      ("v" + __version__) in ui and "design consistency and reliability release" in ui)
try:
    builder.build_source_release(ROOT, ROOT)
except ValueError:
    root_output_rejected = True
else:
    root_output_rejected = False
check("source builder rejects the repository root as its output directory", root_output_rejected)

with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "source"
    (source / "oligoforge").mkdir(parents=True)
    (source / "oligoforge" / "__init__.py").write_text(
        '__version__ = "1.37.0"\n', encoding="utf-8"
    )
    (source / "README.md").write_text("safe source\n", encoding="utf-8")
    unsafe_files = {
        ".env.production": "TOKEN=private\n",
        ".ENV.staging": "TOKEN=also-private\n",
        "credentials.json": "{}\n",
        "Credentials.JSON": "{}\n",
        "private.key": "private\n",
        ".npmrc": "//registry/:_authToken=private\n",
        "coverage.xml": "<coverage/>\n",
        ".cache/generated.json": "{}\n",
        "htmlcov/index.html": "generated\n",
    }
    for relative, content in unsafe_files.items():
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    archive_path, _, _ = builder.build_source_release(source, Path(tmp) / "out")
    with zipfile.ZipFile(archive_path) as archive:
        archived = {entry.filename.split("/", 1)[1] for entry in archive.infolist()}
    check("source builder excludes common credentials and generated local state",
          archived == {"README.md", "oligoforge/__init__.py"}, archived)

with tempfile.TemporaryDirectory() as tmp:
    first_dir = Path(tmp) / "first"
    second_dir = Path(tmp) / "second"
    first_archive, first_manifest, first_digest = builder.build_source_release(ROOT, first_dir)
    second_archive, second_manifest, second_digest = builder.build_source_release(ROOT, second_dir)
    check("source archive name is versioned", first_archive.name == "oligoforge-1.37.0-source.zip")
    check("identical source produces byte-identical archives", first_digest == second_digest)
    check("identical source produces byte-identical checksum records",
          first_manifest.read_text(encoding="utf-8") == second_manifest.read_text(encoding="utf-8"))
    with zipfile.ZipFile(first_archive) as archive:
        members = archive.infolist()
        check("source archive members use a stable path order",
              [entry.filename for entry in members] == sorted(entry.filename for entry in members))
        check("source archive has one versioned root",
              all(entry.filename.startswith("oligoforge-1.37.0/") for entry in members))
        check("source archive excludes generated and secret paths",
              not any("/__pycache__/" in entry.filename or "/.git/" in entry.filename or
                      "/node_modules/" in entry.filename or entry.filename.endswith("/.env")
                      for entry in members))
        check("source archive timestamps are normalized",
              {entry.date_time for entry in members} == {(1980, 1, 1, 0, 0, 0)})
    records = verifier.verify_manifest(first_manifest)
    check("generated checksum verifies", len(records) == 1 and records[0][1], records)
    first_archive.write_bytes(first_archive.read_bytes() + b"tamper")
    records = verifier.verify_manifest(first_manifest)
    check("checksum verification detects artifact modification",
          len(records) == 1 and not records[0][1], records)
    unsafe_manifest = Path(tmp) / "unsafe.sha256"
    unsafe_manifest.write_text("%s  ../escape.zip\n" % ("0" * 64), encoding="utf-8")
    try:
        verifier.verify_manifest(unsafe_manifest)
    except ValueError:
        unsafe_rejected = True
    else:
        unsafe_rejected = False
    check("checksum verification rejects escaping paths", unsafe_rejected)
    generic_artifact = Path(tmp) / "OligoForge-test-binary"
    generic_artifact.write_bytes(b"frozen-binary-fixture")
    generic_manifest = writer.write_sidecar(generic_artifact)
    records = verifier.verify_manifest(generic_manifest)
    check("generic release artifacts receive verifiable checksum sidecars",
          len(records) == 1 and records[0][1], records)

print("release identity tests passed")
