# EBSD KAM Parameterization

PyTex currently implements kernel-average misorientation (KAM) on regular 2D grids with the MTEX-relevant controls needed for parity hardening.

## Supported Parameters

- neighborhood order on square grids
- optional misorientation threshold
- mean or maximum aggregation
- optional restriction to within-grain neighbors through a segmentation mask

## Neighborhood Policy

For square grids, PyTex currently interprets KAM order through cumulative Manhattan neighborhoods. This matches the parity-fixture assumptions used in the current validation suite.

## Current Limits

- Only regular 2D grids are supported.
- Hexagonal-grid neighborhood semantics are not implemented yet.
- Confidence- or phase-aware KAM filtering is not implemented yet.
