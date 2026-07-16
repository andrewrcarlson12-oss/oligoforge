#!/usr/bin/env python3
"""Fail a release when a source directory has too many direct child files."""

import argparse
import os
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LIMIT = 100

# These directory names contain dependencies, caches, or generated build output. Entire subtrees
# are pruned. Additional project-specific names can be supplied with --exclude.
DEPENDENCY_DIRECTORY_NAMES = frozenset({
    ".venv", "venv", "node_modules", "site-packages", ".eggs",
})
CACHE_DIRECTORY_NAMES = frozenset({
    "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".hypothesis",
    ".cache", "cache", ".tox", ".nox", ".nyc_output", ".parcel-cache",
})
BUILD_DIRECTORY_NAMES = frozenset({
    "build", "dist", "_build", "htmlcov", ".next", ".nuxt", "target",
})
EXCLUDED_DIRECTORY_NAMES = frozenset({".git", ".hg", ".svn"}).union(
    DEPENDENCY_DIRECTORY_NAMES, CACHE_DIRECTORY_NAMES, BUILD_DIRECTORY_NAMES)
EXCLUDED_DIRECTORY_SUFFIXES = (".egg-info",)


def _is_excluded_directory(name, excluded_names):
    return (name in excluded_names or
            name.endswith(EXCLUDED_DIRECTORY_SUFFIXES))


def find_overfull_directories(root, limit=DEFAULT_LIMIT, excluded_names=EXCLUDED_DIRECTORY_NAMES):
    """Return ``(path, direct_file_count)`` pairs, sorted by repository-relative path."""
    if limit < 1:
        raise ValueError("limit must be at least 1")

    root = Path(root).resolve()
    excluded_names = frozenset(excluded_names)
    overfull = []
    for current, dirs, files in os.walk(root):
        dirs[:] = sorted(name for name in dirs if not _is_excluded_directory(name, excluded_names))
        if len(files) >= limit:
            overfull.append((Path(current), len(files)))
    return sorted(overfull, key=lambda item: item[0].relative_to(root).as_posix())


def _positive_integer(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", type=Path, default=REPOSITORY_ROOT,
                        help="tree to scan (default: repository root)")
    parser.add_argument("--limit", type=_positive_integer, default=DEFAULT_LIMIT,
                        help="fail at this many direct child files (default: 100)")
    parser.add_argument("--exclude", action="append", default=[], metavar="DIR_NAME",
                        help="additional directory basename to prune (repeatable)")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    if not root.is_dir():
        parser.error("root is not a directory: %s" % root)
    excluded = EXCLUDED_DIRECTORY_NAMES.union(args.exclude)
    overfull = find_overfull_directories(root, args.limit, excluded)
    if overfull:
        print("DIRECTORY FILE-COUNT GATE FAILED: %d director%s at or above the %d-file limit"
              % (len(overfull), "y" if len(overfull) == 1 else "ies", args.limit))
        for path, count in overfull:
            relative = path.relative_to(root).as_posix() or "."
            print("  %s: %d direct child files" % (relative, count))
        return 1

    print("DIRECTORY FILE-COUNT GATE PASS: no scanned directory has %d or more direct child files"
          % args.limit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
