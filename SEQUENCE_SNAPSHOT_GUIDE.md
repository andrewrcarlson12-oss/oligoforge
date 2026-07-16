# Sequence snapshot guide

Build snapshots from local FASTA or bounded FASTA.GZ. Optional CSV/TSV metadata joins by `record_id`, `id` or `accession`. The record disposition ledger distinguishes unique accepted sequences, accepted exact duplicates and rejected records. Metrics report raw record, unique sequence and metadata-group counts separately.

```bash
python -m oligoforge.assurance_cli build-snapshot target.fasta target.snapshot.json --role target --metadata target.csv
python -m oligoforge.assurance_cli delta baseline.snapshot.json followup.snapshot.json delta.json
```

Snapshots are immutable by convention and content-addressed. Editing any field invalidates the stored hash.
