# EBSD Grain Boundaries And Cleanup Foundations

This note records the initial PyTex definition of grain-boundary extraction and small-grain cleanup on regular EBSD grids.

## Boundary Extraction

PyTex currently defines grain-boundary segments from segmentation-adjacent neighbor pairs whose grain labels differ. Each segment records endpoint indices, grain ids, midpoint, and the misorientation across that local boundary.

## Grain Graph Aggregation

PyTex can aggregate the boundary network into a grain graph by grouping boundary segments with the same unordered grain pair. Each grain-graph edge records total shared boundary length, mean boundary misorientation, and the fraction of contributing segments classified as high-angle.

## Cleanup Rule

Small-grain cleanup currently merges the smallest grain below a user-defined size threshold into an adjacent grain chosen by:

1. highest shared boundary count,
2. lowest mean boundary misorientation as a tie-breaker,
3. lowest grain id as a final deterministic tie-breaker.

## Current Limits

- No graph simplification or boundary clustering is implemented yet.
- No curvature, length, or boundary-plane attributes are computed yet.
- Cleanup does not yet use confidence or phase metadata.
