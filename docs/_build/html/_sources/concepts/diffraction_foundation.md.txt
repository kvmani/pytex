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

```{mermaid}
flowchart LR
    core["Canonical core<br/>Frames, symmetry, lattice, phase, orientation"]
    geometry["DiffractionGeometry<br/>Detector, specimen, laboratory"]
    reciprocal["Reciprocal-space semantics<br/>Vectors, planes, zone axes"]
    spots["KinematicSimulation<br/>Spots, families, acceptance masks"]
    indexing["Indexing surface<br/>Candidate ranking and refinement"]
    docs["Workflow and theory notes<br/>Geometry, spots, Bragg rings"]

    core --> geometry
    core --> reciprocal
    geometry --> spots
    reciprocal --> spots
    spots --> indexing
    docs --> geometry
    docs --> spots
    docs --> indexing

    classDef core fill:#f0f4f8,stroke:#334e68,color:#102a43,stroke-width:1.5px;
    classDef domain fill:#d9f2e6,stroke:#2d6a4f,color:#1b4332,stroke-width:1.5px;
    classDef evidence fill:#edf6f9,stroke:#2a9d8f,color:#1d3557,stroke-width:1.5px;
    class core core;
    class geometry,reciprocal,spots,indexing domain;
    class docs evidence;
```

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
- `docs/architecture/diffraction_foundation.md`
- `docs/architecture/multimodal_characterization_foundation.md`

## References

### Normative

- `docs/architecture/diffraction_foundation.md`
- `docs/standards/reference_canon.md`

### Informative

- {doc}`../workflows/diffraction_geometry`
