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

```{mermaid}
flowchart TB
    conventions["conventions.py<br/>Canonical convention set and frame domains"]
    frames["frames.py<br/>ReferenceFrame and FrameTransform"]
    symmetry["symmetry.py<br/>SymmetrySpec and symmetry reduction"]
    lattice["lattice.py<br/>Lattice, Basis, UnitCell, Phase,<br/>directions, planes, reciprocal vectors"]
    orientation["orientation.py<br/>Rotation, Orientation, Misorientation, OrientationSet"]
    provenance["provenance.py<br/>ProvenanceRecord"]

    conventions --> frames
    conventions --> lattice
    frames --> lattice
    frames --> orientation
    symmetry --> orientation
    lattice --> orientation
    provenance --> frames
    provenance --> lattice
    provenance --> orientation

    classDef core fill:#f0f4f8,stroke:#334e68,color:#102a43,stroke-width:1.5px;
    class conventions,frames,symmetry,lattice,orientation,provenance core;
```

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
- `docs/architecture/canonical_data_model.md`
- `docs/architecture/overview.md`
- `docs/standards/reference_canon.md`

## References

### Normative

- `docs/architecture/canonical_data_model.md`
- `docs/standards/reference_canon.md`

### Informative

- [../../figures/reference_frames_vectors.svg](../../figures/reference_frames_vectors.svg)
