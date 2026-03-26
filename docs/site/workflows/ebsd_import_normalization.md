# EBSD: Import Normalization And Manifests

PyTex now includes both the stable manifest contract for EBSD normalization and the first real object-backed bridge entry points for vendor ecosystems.

## Scope

- `EBSDImportManifest` as the stable manifest for normalized EBSD inputs
- `NormalizedEBSDDataset` as the bundle tying a canonical `CrystalMap` to its import manifest
- `normalize_kikuchipy_payload()` and `normalize_pyebsdindex_payload()` as the first vendor-named normalization entry points
- `normalize_kikuchipy_dataset()` and `normalize_pyebsdindex_result()` as object-backed bridge entry points
- JSON manifest write/read helpers for stable interchange
- versioned schema storage under `schemas/`
- explicit validation for required convention, frame, and provenance fields

## Why This Matters

EBSD import is scientifically fragile because the raw arrays are never enough on their own. A normalized dataset must retain:

- source system identity
- source file identity
- phase and point-group declarations
- crystal and specimen frame declarations
- orientation convention and angle unit
- canonical-convention normalization intent

PyTex now treats that context as a stable machine-readable contract rather than notebook folklore.

## Example

```python
from pytex import EBSDImportManifest, NormalizedEBSDDataset

manifest = EBSDImportManifest(
    source_system="kikuchipy",
    source_file="demo.h5",
    phase_name="fcc_demo",
    point_group="m-3m",
    crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
    specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
)

dataset = NormalizedEBSDDataset(crystal_map=crystal_map, manifest=manifest)
```

For vendor-facing normalization, PyTex now also exposes named adapter entry points that accept structured payloads and return normalized datasets:

- `normalize_kikuchipy_payload(...)`
- `normalize_pyebsdindex_payload(...)`
- `normalize_kikuchipy_dataset(...)`
- `normalize_pyebsdindex_result(...)`

The object-backed bridges are intentionally explicit and light-touch. They do not introduce hard import-time dependencies on KikuchiPy or PyEBSDIndex; instead, they extract the required semantic fields from vendor-like objects and convert them into the stable PyTex manifest-plus-`CrystalMap` bundle.

## Manifest IO

The manifest is now not just an in-memory dataclass. It is also a stable interchange artifact:

- `manifest.write_json(path)` writes the schema-backed import manifest
- `read_ebsd_import_manifest(path)` restores and validates it
- `dataset.write_manifest_json(path)` writes the manifest attached to a normalized dataset

## Current Limits

- the bridge layer currently extracts orientation/map semantics only; detector-pattern semantics are still outside this contract
- multi-phase import normalization is still ahead
- dependency-pinned integration tests against live vendor package versions are still ahead

## Related Material

- `docs/standards/data_contracts_and_manifests.md`
- [../../../schemas/ebsd_import_manifest.schema.json](../../../schemas/ebsd_import_manifest.schema.json)
