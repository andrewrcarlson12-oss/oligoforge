# OligoForge tests — run from the repository root

Run the complete deterministic suite with:

    python run_tests.py

`run_tests.py` and `run_tests.sh` recursively discover every `test_*.py` and `ui_*.js` below
`tests/`. Python, pytest, coverage, and Node cache/dependency directories are pruned during
discovery. UI harnesses read `static/index.html` relative to the current directory, and the runners
therefore execute every test from the repository root.

Run one family, one focused regression, or the source-layout release gate with:

    python run_tests.py --python
    python run_tests.py --node
    python tests/test_release_engineering.py
    python tools/check_directory_file_counts.py

The file-count gate fails when any scanned directory contains 100 or more direct child files.
Dependency, cache, and generated build directories listed in the script are pruned; repeat
`--exclude DIR_NAME` only for an additional project-specific generated directory.
