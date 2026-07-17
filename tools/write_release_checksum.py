#!/usr/bin/env python3
"""Write standard SHA-256 sidecars for release artifacts.

Usage:
    python tools/write_release_checksum.py OligoForge-linux-x64
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_sidecar(artifact: Path) -> Path:
    artifact = artifact.resolve()
    if not artifact.is_file():
        raise FileNotFoundError("release artifact not found: %s" % artifact)
    sidecar = artifact.with_name(artifact.name + ".sha256")
    sidecar.write_text("%s  %s\n" % (sha256_file(artifact), artifact.name),
                       encoding="utf-8", newline="\n")
    return sidecar


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args()
    try:
        for artifact in args.artifacts:
            print("wrote", write_sidecar(artifact))
    except OSError as exc:
        print("FAIL", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
