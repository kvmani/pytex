# Diffraction Foundation

This page drills into the diffraction layer as the detector, beam, and reciprocal-space specialization built on the canonical core.

## What The Diffraction Layer Owns

- diffraction geometry for specimen, detector, and laboratory frames
- scattering-vector and reciprocal-space semantics
- detector-space `2θ` and azimuth
- Bragg ring prediction
- kinematic spot simulation
- reflection-family grouping and indexing candidates

## Diffraction Flow

:::{figure} ../../figures/diffraction_foundation_flow.svg
:alt: Diffraction foundation flow showing canonical core objects, diffraction geometry, reciprocal-space semantics, XRD, SAED, indexing, and validation evidence.
:class: architecture-poster-figure
:::

## Why This Layer Matters

Diffraction is where PyTex has to keep geometry honest:

- detector coordinates are not specimen coordinates
- reciprocal-space quantities are not detector-plane offsets
- kinematic spot generation must preserve frame ownership
- local indexing should remain interpretable rather than becoming a black box

## Current State

- detector and beam geometry are implemented
- reciprocal-space primitives are implemented
- kinematic spot generation, family grouping, and local refinement scaffolding exist
- calibrated detector distortion, dynamical intensity, and fuller external-baseline validation remain ahead

## Related Material

- {doc}`../workflows/diffraction_geometry`
- {doc}`../workflows/diffraction_spots`
- {doc}`../architecture/diffraction_foundation`
- {doc}`../architecture/multimodal_characterization_foundation`

## References

### Normative

- {doc}`../architecture/diffraction_foundation`
- {doc}`../standards/reference_canon`

### Informative

- {doc}`../workflows/diffraction_geometry`
