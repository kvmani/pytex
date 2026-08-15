# EBSD KAM Parameterization

PyTex implements kernel-average misorientation (KAM) on rectangular, staggered hexagonal, and
coordinate-graph 2D maps with the MTEX-relevant controls needed for parity hardening.

## Supported Parameters

- neighborhood order on square grids
- cumulative graph-distance order on six-neighbour hexagonal grids
- optional misorientation threshold
- mean or maximum aggregation
- optional restriction to within-grain neighbors through a segmentation mask

## Neighborhood Policy

For four-connected square grids, order is cumulative Manhattan distance; for eight-connected grids,
it is cumulative Chebyshev distance. Hexagonal order is cumulative shortest-path distance on the
six-neighbour staggered-row graph. These definitions make `order=n` a topology statement rather
than a Euclidean-radius approximation.

## Current Limits

- Unstructured point clouds still use a coordinate-radius approximation.
- Property/confidence filtering is an explicit map-selection step rather than an implicit weight.
- Cross-phase edges and optionally cross-grain edges are excluded, but interphase orientation
  relationships are not reinterpreted as same-phase KAM.
