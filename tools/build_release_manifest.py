#!/usr/bin/env python3
"""Generate or verify OligoForge's deterministic source-tree inventory.

``RELEASE_MANIFEST.json`` inventories source files and their SHA-256 digests.  It
does not inventory itself: a manifest cannot contain its own stable digest.
Likewise, the finished source archive is authenticated by its adjacent
``.sha256`` sidecar, which is created only after the archive has been built.

The source walk and baseline exclusions are shared with
``tools/build_source_release.py`` so the inventory and release builder do not
silently develop different path-selection rules.

Usage::

    python tools/build_release_manifest.py
    python tools/build_release_manifest.py --verify
    python tools/build_release_manifest.py --source-root /path/to/source
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_NAME = "RELEASE_MANIFEST.json"
SCHEMA = "oligoforge-release-inventory/v1"
FORMAT_VERSION = 1
SUPPORTED_RELEASE_SERIES = "1.37"

def _load_source_builder():
    """Load the sibling builder without requiring ``tools`` to be a package."""
    path = Path(__file__).with_name("build_source_release.py")
    spec = importlib.util.spec_from_file_location("oligoforge_source_release", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load source release builder: %s" % path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOURCE_BUILDER = _load_source_builder()


def _manifest_excluded(relative: Path) -> bool:
    """Exclude only self-reference; source safety is centralized in the ZIP builder."""
    return relative.as_posix() == MANIFEST_NAME


def manifest_source_files(source_root: Path = ROOT) -> list[Path]:
    """Return the canonical, stable file list for the source-tree manifest."""
    files = SOURCE_BUILDER.source_files(source_root)
    return [relative for relative in files if not _manifest_excluded(relative)]


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _inventory_digest(records: list[dict[str, Any]]) -> str:
    """Hash canonical file records, not the manifest containing that hash."""
    canonical = json.dumps(
        records,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return _sha256_bytes(canonical)


def build_release_manifest(source_root: Path = ROOT) -> dict[str, Any]:
    """Build the in-memory deterministic manifest for *source_root*."""
    source_root = source_root.resolve()
    version = SOURCE_BUILDER._version_from_source(source_root)
    series = ".".join(version.split(".")[:2])
    if series != SUPPORTED_RELEASE_SERIES:
        raise ValueError(
            "release manifest schema targets series %s, found application version %s"
            % (SUPPORTED_RELEASE_SERIES, version)
        )

    records: list[dict[str, Any]] = []
    for relative in manifest_source_files(source_root):
        data = (source_root / relative).read_bytes()
        records.append({
            "path": relative.as_posix(),
            "bytes": len(data),
            "sha256": _sha256_bytes(data),
        })

    archive_name = "oligoforge-%s-source.zip" % version
    return {
        "schema": SCHEMA,
        "manifest_format_version": FORMAT_VERSION,
        "release": {
            "application": "OligoForge",
            "version": version,
            "series": series,
        },
        "inventory_scope": {
            "kind": "source_tree",
            "path_format": "UTF-8 POSIX paths relative to the source root",
            "ordering": "ascending path",
            "digest_algorithm": "SHA-256",
            "source_selection": (
                "tools/build_source_release.py:source_files; manifest self-reference excluded"
            ),
            "self_reference": "%s is intentionally excluded" % MANIFEST_NAME,
            "archive_digest": "external_sidecar",
            "archive_digest_sidecar": archive_name + ".sha256",
            "archive_digest_note": (
                "The archive digest is computed after packaging and is not stored in this "
                "source-tree inventory; verify the adjacent external sidecar."
            ),
        },
        "file_count": len(records),
        "content_bytes": sum(record["bytes"] for record in records),
        "inventory_sha256": _inventory_digest(records),
        "files": records,
    }


def render_release_manifest(manifest: dict[str, Any]) -> str:
    """Render one canonical UTF-8 JSON representation with a final newline."""
    return json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"


def write_release_manifest(
    source_root: Path = ROOT,
    output_path: Path | None = None,
) -> Path:
    """Generate and write the canonical manifest, returning its resolved path."""
    source_root = source_root.resolve()
    output_path = (output_path or source_root / MANIFEST_NAME).resolve()
    try:
        relative_output = output_path.relative_to(source_root)
    except ValueError:
        relative_output = None
    if relative_output is not None and relative_output.as_posix() != MANIFEST_NAME:
        raise ValueError(
            "an in-tree manifest must be written as %s so it can be excluded from itself"
            % MANIFEST_NAME
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        render_release_manifest(build_release_manifest(source_root)),
        encoding="utf-8",
        newline="\n",
    )
    return output_path


def verify_release_manifest(
    source_root: Path = ROOT,
    manifest_path: Path | None = None,
) -> tuple[bool, list[str]]:
    """Verify canonical bytes and current file content; return success and issues."""
    source_root = source_root.resolve()
    manifest_path = (manifest_path or source_root / MANIFEST_NAME).resolve()
    try:
        actual = manifest_path.read_bytes()
    except OSError as exc:
        return False, ["cannot read manifest: %s" % exc]

    try:
        parsed = json.loads(actual.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        return False, ["manifest is not valid UTF-8 JSON: %s" % exc]

    try:
        expected_document = build_release_manifest(source_root)
    except (OSError, ValueError) as exc:
        return False, ["cannot inventory source tree: %s" % exc]
    expected = render_release_manifest(expected_document).encode("utf-8")
    if actual == expected:
        return True, []

    issues: list[str] = []
    if not isinstance(parsed, dict):
        issues.append("manifest root must be a JSON object")
    else:
        if parsed.get("schema") != SCHEMA:
            issues.append("schema is not %s" % SCHEMA)
        if parsed.get("release") != expected_document["release"]:
            issues.append("release identity does not match source")
        actual_records = {
            record.get("path"): record
            for record in parsed.get("files", [])
            if isinstance(record, dict) and isinstance(record.get("path"), str)
        }
        expected_records = {record["path"]: record for record in expected_document["files"]}
        missing = sorted(set(expected_records) - set(actual_records))
        extra = sorted(set(actual_records) - set(expected_records))
        changed = sorted(
            path for path in set(actual_records) & set(expected_records)
            if actual_records[path] != expected_records[path]
        )
        if missing:
            issues.append("missing files: %s" % ", ".join(missing[:10]))
        if extra:
            issues.append("unexpected files: %s" % ", ".join(extra[:10]))
        if changed:
            issues.append("changed files: %s" % ", ".join(changed[:10]))
    if not issues:
        issues.append("manifest JSON is not in canonical form or summary metadata is stale")
    return False, issues


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path,
                        help="manifest path (default: SOURCE_ROOT/RELEASE_MANIFEST.json)")
    parser.add_argument("--verify", action="store_true",
                        help="verify the existing manifest instead of rewriting it")
    return parser


def main() -> int:
    args = _parser().parse_args()
    manifest_path = args.output or args.source_root / MANIFEST_NAME
    if args.verify:
        valid, issues = verify_release_manifest(args.source_root, manifest_path)
        if valid:
            print("OK  ", manifest_path)
            return 0
        for issue in issues:
            print("FAIL", issue)
        return 1
    try:
        written = write_release_manifest(args.source_root, manifest_path)
    except (OSError, ValueError) as exc:
        print("FAIL", exc)
        return 1
    document = json.loads(written.read_text(encoding="utf-8"))
    print("wrote", written)
    print("files", document["file_count"])
    print("inventory_sha256", document["inventory_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
