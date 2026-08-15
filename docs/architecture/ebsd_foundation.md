# EBSD Foundation

This document records the initial Phase 3 implementation posture for EBSD-facing workflows.

## Implemented

- `CrystalMap` as the canonical container for coordinates, orientations, and map-frame metadata
- rectangular-grid validation through `grid_shape`
- ragged hexagonal-grid validation through `grid_kind` and immutable `row_lengths`
- deterministic neighbor-pair generation for rectangular 4/8-connectivity and hexagonal
  6-connectivity, including cumulative graph order
- kernel-average misorientation (KAM) and thresholded grain segmentation on both topologies
- grain reference orientations via within-grain representative selection
- grain reference orientation deviation (GROD) maps
- grain-boundary extraction from segmentation-adjacent pixel pairs
- small-grain cleanup through adjacency-based merging
- grain-graph aggregation from boundary connectivity
- stable EBSD import-manifest and normalization contract surfaces
- manifest JSON IO for stable interchange
- object-backed KikuchiPy/PyEBSDIndex bridge entry points

## Deliberate Current Limits

- no denoising workflows yet
- no dependency-pinned live-package integration tests for KikuchiPy or PyEBSDIndex yet
- no vendor-specific detector/pattern metadata normalization contract yet
- no hexagonal finite-difference curvature/GND stencil or cell-boundary perimeter model yet

## Why This Is The Right First Step

The EBSD layer proves that workflows can reuse the same orientation, frame, and symmetry semantics
already established in the core model while keeping acquisition topology explicit. Rectangular and
hexagonal KAM are useful boundary tests because they require neighborhood, misorientation, and scan
metadata without allowing one grid convention to masquerade as another.

## References

### Normative

- [Canonical Data Model](canonical_data_model.md)
- [Multimodal Characterization Foundation](multimodal_characterization_foundation.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- <a href="../site/theory/ebsd_kam_parameterization.md">EBSD KAM Parameterization</a>
- <a href="../site/theory/ebsd_grain_segmentation_and_grod.md">EBSD Grain Segmentation And GROD</a>
