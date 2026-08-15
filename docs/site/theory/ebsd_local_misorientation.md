# EBSD Local Misorientation Foundations

This note records the PyTex definition of local misorientation on rectangular, staggered
hexagonal, and coordinate-graph EBSD maps.

## Scope

The current implementation supports:

- regular two-dimensional grids with explicit grid shape metadata,
- four- and eight-neighbor connectivity,
- hexagonal grids with explicit row lengths and six-neighbor connectivity,
- cumulative graph-distance neighborhood order, thresholding, and phase/grain masking,
- kernel-average misorientation (KAM) computed as the arithmetic mean of neighbor misorientation angles.

## Definition

For a measurement site $i$ with valid neighbor set $\mathcal{N}(i)$, PyTex currently defines

$$
\mathrm{KAM}(i) = \frac{1}{|\mathcal{N}(i)|} \sum_{j \in \mathcal{N}(i)} \omega(g_i, g_j)
$$

where $\omega(g_i, g_j)$ is the misorientation angle between orientations $g_i$ and $g_j$, optionally reduced by the active crystal symmetry.

On a hexagonal scan, an interior point's first shell contains six sites. PyTex retains the ragged
logical rows rather than padding them into a false rectangle; KAM therefore returns one value per
measurement. Higher orders mean shortest-path distance on this graph.

## Current Limits

- Irregular point clouds use a coordinate radius rather than a reconstructed acquisition lattice.
- Confidence/property thresholds are explicit preprocessing masks, not implicit KAM weights.
- Hexagonal curvature/GND finite differences remain outside this local-angle surface.

## Normative And Informative References

- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*.
- MTEX documentation and public test categories for EBSD neighborhood metrics, as a parity baseline rather than a semantic source of truth.
