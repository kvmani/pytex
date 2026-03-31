# EBSD Foundation

This document records the initial Phase 3 implementation posture for EBSD-facing workflows.

## Implemented

- `CrystalMap` as the canonical container for coordinates, orientations, and map-frame metadata
- regular-grid validation through `grid_shape`
- deterministic neighbor-pair generation for 4- and 8-connectivity
- kernel-average misorientation (KAM) for regular 2D grids
- thresholded grain segmentation on regular 2D grids
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

## Why This Is The Right First Step

Phase 3 should start by proving that EBSD workflows can reuse the same orientation, frame, and symmetry semantics already established in the core model. KAM on regular grids is a useful first boundary test because it requires neighbor topology, misorientation semantics, and map-shape metadata without forcing premature commitment to full grain-analysis infrastructure.

## References

### Normative

- [Canonical Data Model](../site/architecture/canonical_data_model.md)
- [Multimodal Characterization Foundation](../site/architecture/multimodal_characterization_foundation.md)
- [Reference Canon](../site/standards/reference_canon.md)

### Informative

- <a href="../tex/algorithms/ebsd_kam_parameterization.tex">EBSD KAM Parameterization</a>
- <a href="../tex/algorithms/ebsd_grain_segmentation_and_grod.tex">EBSD Grain Segmentation And GROD</a>
