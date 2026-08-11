# EBSD Local Misorientation Foundations

This note records the initial PyTex definition of local misorientation on regular EBSD grids.

## Scope

The current implementation supports:

- regular two-dimensional grids with explicit grid shape metadata,
- four- and eight-neighbor connectivity,
- kernel-average misorientation (KAM) computed as the arithmetic mean of neighbor misorientation angles.

## Definition

For a measurement site $i$ with valid neighbor set $\mathcal{N}(i)$, PyTex currently defines

$$
\mathrm{KAM}(i) = \frac{1}{|\mathcal{N}(i)|} \sum_{j \in \mathcal{N}(i)} \omega(g_i, g_j)
$$

where $\omega(g_i, g_j)$ is the misorientation angle between orientations $g_i$ and $g_j$, optionally reduced by the active crystal symmetry.

## Current Limits

- No thresholded KAM truncation is applied yet.
- No confidence-index masking or phase filtering is implemented yet.
- No grain-aware neighborhood filtering is implemented yet.

## Normative And Informative References

- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*.
- MTEX documentation and public test categories for EBSD neighborhood metrics, as a parity baseline rather than a semantic source of truth.
