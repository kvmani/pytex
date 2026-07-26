# Core Foundation

This page drills into the canonical scientific core of PyTex.

The core is not a convenience layer around arrays. It is the semantic contract that every other subsystem depends on.

## Core Responsibilities

- define frame domains and transforms explicitly
- keep symmetry operators reusable and cacheable
- represent direct and reciprocal lattice semantics without ambiguity
- keep orientation semantics attached to the correct crystal and specimen frames
- preserve provenance so boundary conversions remain inspectable

## Core Object Map

:::{figure} ../../figures/core_foundation_map.svg
:alt: Core foundation map showing conventions, frames, symmetry, lattice, orientation, and provenance around the canonical core.
:class: architecture-poster-figure
:::

## Design Principles

- `ReferenceFrame` and `FrameTransform` keep the coordinate system explicit.
- `SymmetrySpec` centralizes operators and fundamental-sector logic.
- `Lattice` and `Phase` keep structure and symmetry attached to the correct crystal domain.
- `Rotation`, `Orientation`, and `Misorientation` keep orientation meaning stable across workflows.
- `ProvenanceRecord` keeps import and normalization context attached to scientific objects.

## What This Enables

- texture workflows can speak the same orientation language as EBSD workflows
- diffraction workflows can reuse the same lattice and phase semantics
- adapters can normalize vendor objects into one canonical model
- later multimodal and transformation layers can extend rather than replace the core

## Related Material

- {doc}`core_model`
- {doc}`../architecture/canonical_data_model`
- {doc}`../architecture/overview`
- {doc}`../standards/reference_canon`

## References

### Normative

- {doc}`../architecture/canonical_data_model`
- {doc}`../standards/reference_canon`

### Informative

- [../../figures/reference_frames_vectors.svg](../../figures/reference_frames_vectors.svg)
- [../../figures/reference_frame_catalog.svg](../../figures/reference_frame_catalog.svg)
- [../../figures/sample_frame_rd_td_nd.svg](../../figures/sample_frame_rd_td_nd.svg)
