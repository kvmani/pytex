# Specifications

## 1. Purpose

This document defines the stable architectural posture, public primitive set, documentation requirements, and validation doctrine for PyTex.

## 2. Executive Summary

PyTex is a modular Python package with:

- a canonical crystallographic core under `src/pytex/core/`
- orientation and texture-domain models under `src/pytex/texture/`
- EBSD-domain models under `src/pytex/ebsd/`
- diffraction-domain models under `src/pytex/diffraction/`
- optional adapters under `src/pytex/adapters/`
- unstable research methods under `src/pytex/experimental/`

The present repository is beyond a pure scaffold. Core, texture, EBSD, diffraction, and manifest foundations exist today, but the project still requires broader foundational doctrine before it can credibly claim world-class multimodal scope.

## 3. Stable Public Primitive Set

Stable APIs must compose from named scientific primitives rather than naked arrays where scientific meaning would be ambiguous.

### Core Crystallography

- `ReferenceFrame`
- `FrameTransform`
- `ConventionSet`
- `VectorSet`
- `SymmetrySpec`
- `SpaceGroupSpec`
- `Lattice`
- `Basis`
- `UnitCell`
- `Phase`
- `MillerIndex`
- `CrystalDirection`
- `CrystalPlane`
- `ReciprocalLatticeVector`
- `ZoneAxis`
- `ProvenanceRecord`

### Orientation And Texture

- `EulerSet`
- `QuaternionSet`
- `Rotation`
- `RotationSet`
- `Orientation`
- `Misorientation`
- `OrientationSet`
- `PoleFigure`
- `InversePoleFigure`
- `ODF`

### Mapping And Measurement

- `CrystalMap`
- `AcquisitionGeometry`
- `CalibrationRecord`
- `MeasurementQuality`

### Diffraction

- `DiffractionGeometry`
- `DiffractionPattern`
- `ReflectionFamily`
- `ScatteringSetup`

### Transformation

- `OrientationRelationship`
- `TransformationVariant`
- `PhaseTransformationRecord`

### Stable Manifest Surface

- `EBSDImportManifest`
- `ExperimentManifest`
- `BenchmarkManifest`
- `ValidationManifest`
- `WorkflowResultManifest`

Some of these objects are already implemented; others are mandatory architectural targets that must be defined on paper before corresponding stable APIs are added.

`SymmetrySpec` and `SpaceGroupSpec` serve different roles. `SymmetrySpec` is the orientation and reduction surface. `SpaceGroupSpec` is the structure-definition surface used by phases and CIF-backed construction.

## 4. Canonical Convention Policy

PyTex must define one canonical internal convention set and normalize imported data into it once.

The current canonical policy is:

- right-handed Cartesian reference frames
- Bunge-style `phi1`, `Phi`, `phi2` orientation labeling for the canonical Euler convenience path
- unit quaternions stored in `w, x, y, z` order
- reciprocal basis normalized such that `a*_i dot a_j = delta_ij`
- explicit distinction between crystal, specimen, map, detector, laboratory, and reciprocal frames
- preservation of source-system provenance, original convention metadata, and boundary mappings for traceability

The broader frame-chain doctrine now lives in `docs/standards/notation_and_conventions.md` and the multimodal interpretation doctrine now lives in `docs/architecture/multimodal_characterization_foundation.md`.

## 5. Data-Model Requirements

- Metadata objects should be immutable or effectively immutable.
- Vectorized containers should use contiguous NumPy arrays, but stable public vectorized workflows should prefer semantic batch types over naked arrays.
- Derived quantities may be lazy where eager recomputation would be wasteful.
- Public APIs must prefer named domain types over ambiguous positional arrays.
- Stable batch APIs must preserve frame, convention, symmetry, and provenance metadata when scientific meaning depends on them.
- Reference-frame transforms must be reusable objects, not repeated ad hoc matrix math.
- Symmetry operators must be cacheable and reusable across computations.
- Stable domain objects must reject inconsistent frames, phases, symmetries, and calibration metadata at construction time.
- Stable manifests and schemas must be introduced before workflow interchange is treated as stable.
- Space-group, acquisition, calibration, uncertainty, and transformation semantics must not be invented separately inside downstream subsystems.

## 6. Documentation Requirements

- `docs/site/` is the live Sphinx root for the primary browsable and searchable documentation surface.
- `docs/tex/` is the canonical source for major scientific notes.
- `docs/figures/` contains canonical SVG figure sources.
- Sphinx/MyST pages are the default home for concepts, tutorials, workflows, and curated API guidance.
- Core documentation must explain both scalar and batched semantics when a primitive supports vectorized operations.
- Every major stable feature requires:
  - a concept or workflow page in the Sphinx layer
  - a theory or architecture note
  - an algorithm or implementation note when algorithms exist
  - a validation and limitations note
- Root and web-facing docs must link to canonical LaTeX artifacts, concise summaries, and relevant figures.
- Major scientific documents must contain explicit citations and separate normative references from informative ones.
- Conventions for hexagonal and trigonal systems must be fixed centrally and reused across code, figures, examples, and tests.
- Development principles, reference canon rules, and data-contract rules are stable repo policy.
- The documentation architecture must support HTML browsing and PDF-grade scientific exposition without treating either layer as optional.

## 7. Validation Requirements

- `docs/testing/mtex_parity_matrix.md` is the authoritative parity ledger for texture and EBSD categories aligned to MTEX.
- `docs/testing/diffraction_validation_matrix.md` is the authoritative validation ledger for diffraction-facing workflows.
- `docs/testing/structure_validation_matrix.md` is the authoritative validation ledger for structure-import workflows.
- Every relevant external validation category must map to either:
  - a PyTex test or benchmark
  - an explicit foundational status
  - an explicit non-applicability note
- PyTex must add tests that external tools do not cover, especially:
  - provenance and manifest integrity
  - vendor reference-frame conversion safety
  - adapter interoperability
  - figure and documentation asset integrity
  - notation normalization and integer-index round-trips for special crystal systems such as HCP
  - structure-import and space-group semantics
  - calibration, uncertainty, and workflow reproducibility contracts

## 8. Delivery Posture

The detailed implementation ladder lives in `docs/roadmap/implementation_roadmap.md`. Repository claims, README status text, and CI expectations must remain synchronized with that roadmap.

## References

### Normative

- `docs/standards/reference_canon.md`
- `docs/standards/notation_and_conventions.md`
- `docs/standards/data_contracts_and_manifests.md`

### Informative

- `docs/roadmap/implementation_roadmap.md`
