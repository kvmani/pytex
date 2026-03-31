# EBSD: Import Normalization And Manifests

PyTex now treats EBSD normalization as a stable semantic contract rather than an ad hoc adapter
step.

## Stable Surfaces

- `EBSDImportManifest`
- `NormalizedEBSDDataset`
- `CrystalMapPhase`
- `normalize_kikuchipy_payload(...)`
- `normalize_kikuchipy_dataset(...)`
- `normalize_pyebsdindex_payload(...)`
- `normalize_pyebsdindex_result(...)`
- `normalize_ebsd(...)`

## What The Normalized Surface Preserves

- source system and source file identity
- crystal, specimen, and map frame declarations
- orientation convention and angle unit
- phase declarations, including multiphase tables when present
- per-point phase identifiers for multiphase crystal maps
- reproducible manifest IO through stable JSON schemas

## Single-Phase Example

```python
from pytex import EBSDImportManifest

manifest = EBSDImportManifest(
    source_system="kikuchipy",
    source_file="indexed.h5",
    phase_name="fcc_demo",
    point_group="m-3m",
    crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
    specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
)
```

## Multiphase Example

When the incoming dataset carries an ORIX or KikuchiPy phase table plus per-point `phase_id`
values, PyTex preserves that structure explicitly.

```python
dataset = normalize_ebsd(
    ebsd_signal,
    frames=(crystal_frame, specimen_frame, map_frame),
    phases={
        1: ferrite_phase,
        2: austenite_phase,
    },
)

summary = dataset.crystal_map.phase_summary()
alpha_map = dataset.crystal_map.select_phase("ferrite")
```

The normalized crystal map remains one object, but texture extraction from a multiphase map
requires an explicit phase selector so PyTex never mixes phase populations silently.

## Experiment And Transformation Manifests

Normalized EBSD datasets can now feed the broader manifest family directly:

- `dataset.write_manifest_json(...)`
- `dataset.crystal_map.to_experiment_manifest(...)`
- phase-resolved acquisition context carried into `ExperimentManifest`

## Current Limits

- detector-pattern semantics are still outside the normalized EBSD dataset contract
- live dependency-pinned integration tests against specific KikuchiPy releases remain ahead
- per-point confidence weighting is not yet part of graph construction or segmentation

## Related Material

- {doc}`ebsd_kam`
- {doc}`ebsd_grains`
- {doc}`ebsd_to_texture_outputs`
- {doc}`orix_kikuchipy_interop`
- [../../../schemas/ebsd_import_manifest.schema.json](../../../schemas/ebsd_import_manifest.schema.json)
