# Orientation Space, Symmetry Reduction, and Disorientation in PyTex

## Purpose

PyTex represents orientations by unit quaternions and rotation matrices while keeping crystal and specimen reference frames explicit. The central Phase 2 question is how these rotations are reduced under symmetry.

## Orientation Equivalence

If $g$ is an orientation mapping crystal coordinates into specimen coordinates, then symmetry-equivalent orientations take the form

$$
g' = s g c
$$

where $s$ belongs to specimen symmetry and $c$ belongs to crystal symmetry.

## Disorientation

Given two orientations, PyTex computes a base misorientation and then searches the symmetry orbit for the minimum-angle representative. This minimum-angle representative is treated as the current disorientation used in kernel evaluation and neighborhood calculations.

## Exact Orientation Representative

For stable orientation-space comparison, PyTex now exposes an exact orbit-reduction rule in the quaternion hemisphere. For each symmetry-equivalent orientation candidate, it computes a canonical quaternion with non-negative scalar part and then selects the representative with maximal scalar part, equivalently minimal rotation angle to the identity.

This gives an exact minimum-angle representative even when no workflow-specific reference orientation is provided, with a deterministic tie-break for symmetry-boundary cases.

## Current Implementation Boundary

PyTex presently implements exact orbit reduction for the supported proper point groups already present in the codebase and minimum-angle disorientation. What remains ahead is a proof-oriented class-by-class boundary catalog and broader external parity on those exact boundaries.

## Normative References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, 1969. DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- M. De Graef, *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, 2003. DOI: <https://doi.org/10.1017/CBO9780511615092>.
