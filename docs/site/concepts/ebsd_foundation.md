# EBSD Foundation

This page drills into the EBSD layer as the map-and-neighborhood specialization built on the canonical orientation model.

## What The EBSD Layer Owns

- `CrystalMap` as the canonical orientation-map container
- map-frame metadata plus rectangular, staggered-hexagonal, and coordinate-graph semantics
- four/eight-neighbour rectangular and six-neighbour hexagonal topology on 2D scans
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
- direct `.ang` import preserves square or EDAX/TSL `HexGrid` topology; `.ctf` import preserves
  rectangular topology when indexed rows are complete; EDAX OIM HDF5 `.oh5`/`.h5` import (one
  container under two extensions) preserves both, and carries every per-point scalar channel the
  file holds rather than only those a text row has room for
- KAM, segmentation, GROD, boundary extraction, cleanup, and grain graphs are implemented on the
  shared topology graph
- stable import manifests exist
- hex curvature/GND stencils, the remaining HDF5-family readers (Oxford H5OINA, Bruker), and
  richer detector/pattern metadata remain ahead

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
