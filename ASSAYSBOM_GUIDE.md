# AssaySBOM guide

An AssaySBOM is a portable versioned record of assays, components, chemistry, reaction conditions, channels, intended groups, near neighbors, interpretation rules, evidence, review status and repair history. Modified order notation is retained separately from the bare sequence used by computational models. The deterministic `assaysbom_id` changes whenever normalized content changes.

Legacy objects containing `forward`, `reverse` and optional `probe` fields migrate to the v1 component list. Migration never invents validation evidence.

CLI example:

```bash
python -m oligoforge.assurance_cli build-assaysbom assay.json assay.sbom.json
```

The JSON schema is `schemas/assaysbom.schema.json`. `templates/assaysbom_components.csv` is a spreadsheet-friendly intake template for collecting component records before converting them to the documented JSON structure; the v1 CLI deliberately accepts validated JSON rather than guessing how partially populated CSV rows should be grouped.

AssaySBOM proves content identity and chain of custody only. It is not proof of amplification, fluorescence, analytical sensitivity/specificity, clinical performance, or regulatory acceptability.
