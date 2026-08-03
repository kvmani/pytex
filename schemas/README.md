# Schemas

This directory stores versioned machine-readable schemas used by PyTex manifests, benchmark fixtures, and workflow interchange artifacts.

Current stable schemas:

- `ebsd_import_manifest.schema.json`
- `experiment_manifest.schema.json`
- `benchmark_manifest.schema.json`
- `validation_manifest.schema.json`
- `workflow_result_manifest.schema.json`
- `transformation_manifest.schema.json`
- `composite_saed_manifest.schema.json`
- `measured_saed_pattern.schema.json`

Stable schemas must remain versioned here so adapters and workflows do not invent incompatible ad hoc locations.
