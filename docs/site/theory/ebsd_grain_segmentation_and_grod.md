# EBSD Grain Segmentation And GROD Foundations

This note records the initial PyTex definition of thresholded grain segmentation and grain reference orientation deviation (GROD) on regular EBSD grids.

## Segmentation Rule

Two neighboring pixels are assigned to the same connected component when their misorientation angle is less than or equal to a user-specified threshold. Connected components are then promoted to grains.

## Representative Grain Orientation

PyTex currently selects a representative measured orientation from each grain by minimizing the summed within-grain misorientation to the other member orientations. This is a pragmatic medoid-style choice that avoids introducing an under-specified orientation averaging convention at this stage.

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
