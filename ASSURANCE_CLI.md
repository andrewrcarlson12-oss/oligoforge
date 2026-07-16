# Assurance CLI

`python -m oligoforge.assurance_cli --help` exposes offline commands:

- `build-assaysbom`
- `build-snapshot`
- `delta`
- `drift-scan`
- `ofvr`
- `package`

Every input and output path is explicit. The CLI performs no network request. JSON is emitted with deterministic key ordering and a trailing newline.

Minimal offline lifecycle:

```bash
python -m oligoforge.assurance_cli build-assaysbom assay.json assay.sbom.json
python -m oligoforge.assurance_cli build-snapshot baseline.fasta baseline.json --role target --metadata baseline.csv
python -m oligoforge.assurance_cli build-snapshot followup.fasta followup.json --role target --metadata followup.csv
python -m oligoforge.assurance_cli delta baseline.json followup.json delta.json
python -m oligoforge.assurance_cli drift-scan assay.sbom.json baseline.json followup.json scan.json
python -m oligoforge.assurance_cli ofvr scan.json vulnerabilities.json --year 2026
python -m oligoforge.assurance_cli package assay.sbom.json evidence-package.json --snapshots baseline.json followup.json --deltas delta.json --scans scan.json --vulnerabilities vulnerabilities.json
```

The commands overwrite the explicitly named output file. Run them in a version-controlled or otherwise governed evidence directory when audit history matters.
