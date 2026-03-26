# Implementation Roadmap

This document is the authoritative phase-by-phase guide for building PyTex from foundation work through the first stable scientific releases.

Current implementation posture:

- Phases 0 and 1 are in place as repository foundations.
- Phase 2 has a working implementation for rotations, disorientation, exact orbit-based symmetry reduction in orientation space, PF/IPF containers, class-specific IPF sector reduction, and discrete kernel ODF support.
- Full external-baseline parity for exact class-by-class fundamental-region boundaries and harmonic PF-to-ODF inversion remains ahead.
- Phase 3 now has a substantial regular-grid implementation through KAM, grain segmentation, GROD, grain-boundary, cleanup, grain-graph aggregation, import-manifest normalization contracts, manifest IO, and object-backed vendor bridge entry points.
- Phase 4 now has a real foundation through detector-space geometry, scattering vectors, reciprocal-space primitives, kinematic spot simulation, reflection families, detector masks, proxy intensity ranking, detector-space indexing association primitives, local candidate refinement, and family-level indexing reports.

## Phase 0: Documentation And Standards Foundation

Goals:

- finalize mission, scope, and product principles
- establish notation and convention rules
- establish the hybrid documentation architecture: Sphinx for the primary docs surface, LaTeX for scientific notes, SVG for canonical figures
- create validation doctrine and MTEX parity matrix

Deliverables:

- `docs/standards/documentation_architecture.md`
- `docs/site/README.md`
- updated foundational repo contracts reflecting the hybrid docs policy

## Phase 1: Canonical Data Model

Goals:

- implement stable primitives for frames, symmetry, lattice, orientation, provenance, and domain containers
- define canonical internal conventions
- enforce explicit semantics at public API boundaries

Deliverables:

- `src/pytex/core/`
- core tests
- canonical data-model docs

## Phase 2: Orientation And Texture Foundation

Goals:

- rotation and misorientation operations
- PF and IPF container semantics
- ODF container and reconstruction foundations
- kernels, harmonics, and inversion roadmap

## Phase 3: EBSD Post-Processing Foundation

Goals:

- crystal-map workflows
- grain and neighborhood semantics
- KAM/GROD and map-cleaning roadmap
- adapter boundaries for KikuchiPy and PyEBSDIndex

## Phase 4: Diffraction Foundation

Goals:

- diffraction geometry model
- stereographic and reciprocal-space foundations
- kinematic XRD and SAED workflow scaffolding
- teaching-grade geometry figures

## Phase 5: Experimental Incubator

Goals:

- unstable research methods with explicit quality gates
- ML-assisted indexing
- parent reconstruction
- distortion correction
- dynamical TEM and related frontier work
