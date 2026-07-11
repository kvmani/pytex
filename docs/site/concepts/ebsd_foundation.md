# EBSD Foundation

This page drills into the EBSD layer as the map-and-neighborhood specialization built on the canonical orientation model.

## What The EBSD Layer Owns

- `CrystalMap` as the canonical orientation-map container
- map-frame metadata and regular-grid semantics
- neighborhood topology on 2D grids
- KAM, segmentation, GROD, boundaries, and cleanup
- grain-graph aggregation and reproducible normalization contracts
- import-manifest normalization for vendor or third-party inputs

## EBSD Flow

:::{figure} ../../figures/ebsd_foundation_flow.svg
:alt: EBSD foundation flow showing adapter normalization, import manifests, CrystalMap, neighborhoods, grains, boundaries, and texture outputs.
:class: architecture-poster-figure
:::

## Why This Layer Matters

EBSD is where PyTex must combine scientific semantics with measurement topology.

- the map frame must stay distinct from the specimen frame unless a workflow explicitly links them
- neighborhood logic must be deterministic
- grain cleanup must remain reproducible
- normalization must preserve source-system meaning rather than flattening it away

## Current State

- `CrystalMap` exists as the canonical map container
- KAM, segmentation, GROD, boundary extraction, cleanup, and grain graphs are implemented
- stable import manifests exist
- richer vendor detector/pattern metadata normalization remains ahead

## Related Material

- {doc}`../workflows/ebsd_grains`
- {doc}`../workflows/ebsd_kam`
- {doc}`../architecture/ebsd_foundation`
- {doc}`../architecture/multimodal_characterization_foundation`

## References

### Normative

- {doc}`../architecture/ebsd_foundation`
- {doc}`../standards/reference_canon`

### Informative

- {doc}`../workflows/ebsd_import_normalization`
