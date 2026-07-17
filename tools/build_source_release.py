#!/usr/bin/env python3
"""Build a byte-reproducible OligoForge source ZIP and SHA-256 sidecar.

The archive uses a stable path order, a fixed timestamp, normalized Unix file
modes, and stored (uncompressed) entries. Given identical source bytes, it is
therefore identical across operating systems and Python/zlib versions.

Usage:
    python tools/build_source_release.py --output-dir dist
    SOURCE_DATE_EPOCH=1704067200 python tools/build_source_release.py
"""

from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import re
import stat
import zipfile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DATE_EPOCH = 315532800  # 1980-01-01T00:00:00Z, ZIP's lower bound.
MAX_SOURCE_DATE_EPOCH = 4354819198  # 2107-12-31T23:59:58Z, ZIP's upper bound.
EXCLUDED_DIRECTORY_NAMES = frozenset({
    ".cache", ".git", ".hypothesis", ".mypy_cache", ".nox", ".nyc_output",
    ".pytest_cache", ".ruff_cache", ".tox", ".venv", ".vscode", "__pycache__",
    "build", "coverage", "dist", "htmlcov", "node_modules", "panels", "projects",
    "venv",
})
EXCLUDED_FILE_NAMES = frozenset({
    ".coverage", ".ds_store", ".env", ".npmrc", ".pypirc", "thumbs.db",
    "coverage.xml", "credentials.json", "desktop.ini", "junit.xml", "secrets.json",
    "service-account.json",
})
EXCLUDED_SUFFIXES = (".egg-info", ".pyc", ".pyo")
SECRET_SUFFIXES = (".key", ".p12", ".pfx")
ENV_FILE = re.compile(r"^\.env(?:\..+)?$")


def _version_from_source(source_root: Path) -> str:
    """Read ``__version__`` without importing application dependencies."""
    init_path = source_root / "oligoforge" / "__init__.py"
    tree = ast.parse(init_path.read_text(encoding="utf-8"), filename=str(init_path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "__version__"
               for target in node.targets):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                return node.value.value
    raise ValueError("oligoforge.__version__ is missing or is not a string literal")


def _is_excluded(relative: Path) -> bool:
    if any(part.lower() in EXCLUDED_DIRECTORY_NAMES or part.lower().endswith(".egg-info")
           for part in relative.parts[:-1]):
        return True
    name = relative.name
    lowered = name.lower()
    return (
        lowered in EXCLUDED_FILE_NAMES
        or lowered.endswith(EXCLUDED_SUFFIXES)
        or lowered.endswith(SECRET_SUFFIXES)
        or ENV_FILE.fullmatch(lowered) is not None
    )


def source_files(source_root: Path, output_dir: Path | None = None) -> list[Path]:
    """Return a stable list of distributable regular files below *source_root*."""
    source_root = source_root.resolve()
    resolved_output = output_dir.resolve() if output_dir is not None else None
    files: list[Path] = []
    for current, directories, names in os.walk(source_root, topdown=True):
        current_path = Path(current)
        for name in directories:
            directory = current_path / name
            if directory.is_symlink():
                raise ValueError("source archives do not follow symbolic links: %s" %
                                 directory.relative_to(source_root))
        directories[:] = sorted(
            name for name in directories
            if name.lower() not in EXCLUDED_DIRECTORY_NAMES and
            not name.lower().endswith(".egg-info")
        )
        for name in sorted(names):
            path = current_path / name
            relative = path.relative_to(source_root)
            if _is_excluded(relative):
                continue
            if resolved_output is not None:
                try:
                    path.resolve().relative_to(resolved_output)
                    continue
                except ValueError:
                    pass
            if path.is_symlink():
                raise ValueError("source archives do not follow symbolic links: %s" % relative)
            if path.is_file():
                files.append(relative)
    return sorted(files, key=lambda item: item.as_posix())


def _zip_timestamp(source_date_epoch: int) -> tuple[int, int, int, int, int, int]:
    normalized_epoch = min(
        max(source_date_epoch, DEFAULT_SOURCE_DATE_EPOCH), MAX_SOURCE_DATE_EPOCH
    )
    when = datetime.fromtimestamp(normalized_epoch, timezone.utc)
    # ZIP stores seconds at two-second precision. Normalizing avoids implementation variance.
    return when.year, when.month, when.day, when.hour, when.minute, when.second // 2 * 2


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_source_release(
    source_root: Path = ROOT,
    output_dir: Path | None = None,
    source_date_epoch: int = DEFAULT_SOURCE_DATE_EPOCH,
) -> tuple[Path, Path, str]:
    """Create the deterministic archive, checksum sidecar, and return their paths/hash."""
    source_root = source_root.resolve()
    output_dir = (output_dir or source_root / "dist").resolve()
    if output_dir == source_root:
        raise ValueError("output directory must not be the source root")
    output_dir.mkdir(parents=True, exist_ok=True)
    version = _version_from_source(source_root)
    archive_name = "oligoforge-%s-source.zip" % version
    archive_path = output_dir / archive_name
    checksum_path = output_dir / (archive_name + ".sha256")
    archive_root = "oligoforge-%s" % version
    timestamp = _zip_timestamp(source_date_epoch)

    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_STORED,
                         allowZip64=True, strict_timestamps=True) as archive:
        archive.comment = b""
        for relative in source_files(source_root, output_dir):
            member = "%s/%s" % (archive_root, relative.as_posix())
            info = zipfile.ZipInfo(member, date_time=timestamp)
            info.create_system = 3
            info.compress_type = zipfile.ZIP_STORED
            mode = 0o755 if relative.suffix == ".sh" else 0o644
            info.external_attr = (stat.S_IFREG | mode) << 16
            info.flag_bits = 0x800  # UTF-8 member names.
            archive.writestr(info, (source_root / relative).read_bytes())

    digest = sha256_file(archive_path)
    checksum_path.write_text("%s  %s\n" % (digest, archive_path.name),
                             encoding="utf-8", newline="\n")
    return archive_path, checksum_path, digest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=ROOT,
                        help="repository root (default: inferred from this script)")
    parser.add_argument("--output-dir", type=Path,
                        help="artifact directory (default: SOURCE_ROOT/dist)")
    parser.add_argument(
        "--source-date-epoch",
        type=int,
        default=int(os.environ.get("SOURCE_DATE_EPOCH", DEFAULT_SOURCE_DATE_EPOCH)),
        help="normalized ZIP timestamp (default: SOURCE_DATE_EPOCH or 1980-01-01)",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    output_dir = args.output_dir or args.source_root / "dist"
    try:
        archive, checksum, digest = build_source_release(
            source_root=args.source_root,
            output_dir=output_dir,
            source_date_epoch=args.source_date_epoch,
        )
    except (OSError, ValueError) as exc:
        print("FAIL", exc)
        return 1
    print("built", archive)
    print("sha256", digest)
    print("checksum", checksum)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
