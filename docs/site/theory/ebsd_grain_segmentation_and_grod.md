# EBSD Grain Segmentation And GROD Foundations

This note records the initial PyTex definition of thresholded grain segmentation and grain reference orientation deviation (GROD) on regular EBSD grids.

## Segmentation Rule

Two neighboring pixels are assigned to the same connected component when their misorientation angle is less than or equal to a user-specified threshold. Connected components are then promoted to grains.

## Representative Grain Orientation

PyTex currently selects a representative measured orientation from each grain by minimizing the summed within-grain misorientation to the other member orientations. This is a pragmatic medoid-style choice that avoids introducing an under-specified orientation averaging convention at this stage.

Two members of one grain that carry the same orientation to within the resolution of the arithmetic — under $10^{-6}$ rad, a ten-thousandth of a degree — have equal totals, and a cluster made entirely of such members has no medoid at all. In that case the lowest member index is taken, as a definition rather than as a consequence of summation order, so the grain reference orientation of a uniform grain is reproducible across machines and BLAS builds.

### Evaluating the medoid without searching the symmetry group

The disorientation is $\min_{S_l, S_r \in G} \omega(S_l M S_r)$, so a literal medoid costs the full operator search on every pair. It does not have to. The rotation angle is a bi-invariant metric on $SO(3)$, so for any $S_l, S_r$

$$
\omega(S_l M S_r) \ \ge\ \omega(S_l S_r) - \omega(M).
$$

The product $S_l S_r$ ranges over $G$. Where it is the identity the conjugate has exactly $\omega(M)$, so those branches never improve on the one in hand; where it is not, the angle is at least $\theta_\mathrm{min} - \omega(M)$, with $\theta_\mathrm{min}$ the smallest non-identity rotation angle in the group ($90^\circ$ for cubic, $60^\circ$ for hexagonal, $180^\circ$ for orthorhombic). Hence an angle measured on a single branch and found below $\theta_\mathrm{min}/2$ **is** the disorientation.

PyTex therefore brings every member of a grain onto the symmetry branch nearest the grain's mean, after which the pair angle is $2\arccos|q_i \cdot q_j|$ and the whole grain is one dense Gram matrix. The certificate is then an $O(n)$ test: if every member lies within $\theta_\mathrm{min}/4$ of the mean, the triangle inequality places every pair below the threshold. A grain too spread to certify falls back to the full group search, so the definition above is what is computed in every case.

## GROD Definition

For pixel $i$ in grain $G$, PyTex currently defines

$$
\mathrm{GROD}(i) = \omega(g_i, g_G^\mathrm{ref})
$$

where $g_G^\mathrm{ref}$ is the representative grain orientation and $\omega$ is optionally symmetry-reduced.

## Current Limits

- No grain-cleaning or minimum-size post-processing is applied yet.
- No orientation-mean estimator is used yet.
- No boundary classification or graph abstraction is implemented yet.
