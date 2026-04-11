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

## What Is Verified Today

The current validation surface checks three different things separately:

- manifest schema stability and JSON round-trips
- payload normalization into canonical PyTex frames, phases, and crystal-map semantics
- downstream EBSD-to-texture behavior such as phase-aware ODF, PF, and IPF extraction

The highest-signal executable coverage currently lives in:

- `tests/unit/test_plotting_and_adapters.py`
- `tests/integration/test_cli_and_optional_adapters.py`
- {doc}`../validation/index`

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
- ORIX and KikuchiPy bridge coverage is strongest at the canonicalization boundary, not yet as a
  broad package-version compatibility matrix

## Interpretation Rule

Treat `normalize_ebsd(...)` as the point where PyTex fixes scientific meaning, not as a guarantee
that every source-system nuance has been preserved or validated. If a workflow depends on
vendor-specific metadata, detector calibration state, or package-version-specific behavior, record
that dependency explicitly and add validation before treating it as stable.

## Related Material

- {doc}`ebsd_kam`
- {doc}`ebsd_grains`
- {doc}`ebsd_to_texture_outputs`
- {doc}`orix_kikuchipy_interop`
- [../../../schemas/ebsd_import_manifest.schema.json](../../../schemas/ebsd_import_manifest.schema.json)
