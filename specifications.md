# Specifications

## 1. Purpose

This document defines the initial architecture, product constraints, data-model expectations, documentation posture, and validation rules for PyTex.

## 2. Executive Summary

PyTex should be implemented as a modular Python package with:

- a canonical core data model under `src/pytex/core/`
- texture-domain models under `src/pytex/texture/`
- EBSD-domain models under `src/pytex/ebsd/`
- diffraction-domain models under `src/pytex/diffraction/`
- optional adapters under `src/pytex/adapters/`
- unstable research methods under `src/pytex/experimental/`

The first implementation stage focuses on the data model and governance required to build serious scientific functionality safely.

## 3. Stable Public Primitive Set

Stable public primitives must include:

- `ReferenceFrame`
- `FrameTransform`
- `SymmetrySpec`
- `Lattice`
- `UnitCell`
- `Basis`
- `Phase`
- `MillerIndex`
- `CrystalDirection`
- `CrystalPlane`
- `Rotation`
- `Orientation`
- `Misorientation`
- `OrientationSet`
- `CrystalMap`
- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `DiffractionGeometry`
- `DiffractionPattern`
- `ProvenanceRecord`

Additional supporting primitives may be introduced where needed, but stable APIs must compose from these concepts.

## 4. Canonical Convention Policy

PyTex must define one canonical internal convention set and normalize imported data into it once.

The initial canonical policy is:

- right-handed Cartesian reference frames
- Bunge-style `phi1`, `Phi`, `phi2` orientation convention for Euler-angle labeling
- unit quaternions stored in `w, x, y, z` order
- reciprocal basis normalized such that `a*_i dot a_j = delta_ij`
- explicit distinction between crystal, specimen, map, detector, laboratory, and reciprocal frames
- preservation of source-system provenance and original convention metadata for round-trip traceability

## 5. Data-Model Requirements

- Metadata objects should be immutable or effectively immutable.
- Vectorized containers should use contiguous NumPy arrays.
- Derived quantities may be lazy where eager recomputation would be wasteful.
- Public APIs must prefer named domain types over ambiguous positional arrays.
- Reference-frame transforms must be reusable objects, not repeated ad hoc matrix math.
- Symmetry operators must be cacheable and reusable across computations.
- Stable domain objects must reject inconsistent frames, phases, or symmetries at construction time.
- Machine-readable manifests and schemas must be introduced before adapter-heavy workflows become stable.

## 6. Documentation Requirements

- `docs/site/` is the planned Sphinx root for the primary browsable and searchable documentation surface.
- `docs/tex/` is the canonical source for major scientific notes.
- `docs/figures/` contains canonical SVG figure sources.
- Sphinx/MyST pages are the default home for concepts, tutorials, workflows, and curated API guidance.
- Every major stable feature requires:
  - concept or usage page in the Sphinx layer
  - theory note
  - algorithm/implementation note
  - validation/limitations note
- Root and web-facing docs must link to canonical LaTeX artifacts, concise summaries, and relevant figures.
- Major scientific documents must contain explicit citations and separate normative references from informative ones.
- Conventions for hexagonal and trigonal systems must be fixed centrally and reused across code, figures, examples, and tests.
- Development principles and data-contract rules must be documented explicitly and treated as stable repo policy.
- The documentation architecture must explicitly support HTML browsing/search and PDF-grade scientific notes without treating either layer as optional.

## 7. Validation Requirements

- `docs/testing/mtex_parity_matrix.md` is the authoritative parity ledger.
- Every relevant MTEX public test category must map to either:
  - a PyTex test or benchmark
  - an explicit non-applicability note
- PyTex must add tests that MTEX does not cover, especially:
  - provenance and manifest integrity
  - vendor reference-frame conversion safety
  - adapter interoperability
  - figure regression and documentation asset integrity
  - notation normalization and integer-index round-trips for special crystal systems such as HCP

## 8. Delivery Phases

1. Documentation and standards foundation
2. Canonical data model
3. Orientation and texture foundation
4. EBSD post-processing and adapters
5. Diffraction foundation
6. Experimental incubator

The detailed phase guide lives in `docs/roadmap/implementation_roadmap.md`.
