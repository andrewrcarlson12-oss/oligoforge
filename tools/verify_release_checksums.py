#!/usr/bin/env python3
"""Verify standard SHA-256 sidecar manifests without trusting their paths.

Usage:
    python tools/verify_release_checksums.py dist/*.sha256
"""

from __future__ import annotations

import argparse
import hashlib
import hmac
from pathlib import Path
import re


CHECKSUM_LINE = re.compile(r"^([0-9a-fA-F]{64}) [ *](.+)$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_target(manifest: Path, name: str) -> Path:
    candidate_name = Path(name)
    if candidate_name.is_absolute() or ".." in candidate_name.parts:
        raise ValueError("unsafe checksum path: %s" % name)
    root = manifest.resolve().parent
    candidate = (root / candidate_name).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError("checksum path escapes its manifest directory: %s" % name) from exc
    return candidate


def verify_manifest(manifest: Path) -> list[tuple[Path, bool, str, str]]:
    """Return ``(path, valid, expected, actual)`` records for one manifest."""
    manifest = manifest.resolve()
    if not manifest.is_file():
        raise FileNotFoundError("checksum manifest not found: %s" % manifest)
    records: list[tuple[Path, bool, str, str]] = []
    for line_number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        match = CHECKSUM_LINE.fullmatch(raw)
        if match is None:
            raise ValueError("%s:%d: invalid SHA-256 line" % (manifest, line_number))
        expected, name = match.groups()
        target = _safe_target(manifest, name)
        if not target.is_file():
            records.append((target, False, expected.lower(), "missing"))
            continue
        actual = sha256_file(target)
        records.append((target, hmac.compare_digest(expected.lower(), actual),
                        expected.lower(), actual))
    if not records:
        raise ValueError("checksum manifest has no records: %s" % manifest)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifests", nargs="+", type=Path)
    args = parser.parse_args()
    failed = False
    for manifest in args.manifests:
        try:
            records = verify_manifest(manifest)
        except (OSError, ValueError) as exc:
            print("FAIL", exc)
            failed = True
            continue
        for path, valid, expected, actual in records:
            if valid:
                print("OK  ", path)
            else:
                print("FAIL", path, "expected=%s actual=%s" % (expected, actual))
                failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
