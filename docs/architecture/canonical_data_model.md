# Canonical Data Model

The canonical data model is the scientific backbone of PyTex.

## Why It Exists

Texture and diffraction workflows are unusually sensitive to hidden assumptions:

- which frame a vector belongs to
- whether a basis is direct or reciprocal
- what symmetry group is active
- which Euler convention is implied
- whether detector coordinates and specimen coordinates have already been reconciled

Allowing these assumptions to travel as unnamed arrays creates avoidable scientific risk. PyTex therefore defines stable primitives with explicit semantics.

## Core Design Rules

- Metadata objects are immutable or effectively immutable.
- Vectorized data is stored in contiguous NumPy arrays.
- Frames, transforms, and symmetry are reusable objects.
- Provenance is attached where scientific traceability matters.
- Stable APIs must prefer domain types over generic arrays.
- Import normalization must preserve both canonical internal semantics and original source semantics.
- Special-system notation, such as hexagonal four-index forms, must be normalized at the boundary and documented centrally.
- Domain objects should fail fast on inconsistent frame, phase, and symmetry combinations.
- Stable workflow interchange should use versioned manifests and schemas rather than ad hoc dicts.
- Major stable objects and results that cross workflow or tool boundaries should expose versioned JSON contracts that preserve enough semantics for reconstruction.

## Canonical Primitives

- `ReferenceFrame`: named right-handed coordinate frame with explicit domain
- `FrameTransform`: reusable rigid transform between frames
- `VectorSet`: batch of vectors sharing one explicit reference frame
- `SymmetrySpec`: point-group or specimen-symmetry description plus operators
- `SpaceGroupSpec`: structure-facing space-group identity attached to a crystal frame
- `Basis`: direct or reciprocal basis matrix plus frame information
- `Lattice`: lattice parameters and derived bases
- `UnitCell`: lattice plus atomic basis
- `Phase`: named material phase with lattice, symmetry, crystal frame, and optional space-group or formula metadata
- `AcquisitionGeometry`: shared multimodal experiment-frame and transform context
- `CalibrationRecord`: explicit calibration state and uncertainty metadata
- `MeasurementQuality`: stable quality, validity, and masking surface
- `EulerSet`: batch of convention-aware Euler triples
- `QuaternionSet`: batch of canonical unit quaternions
- `ScatteringSetup`: shared scattering-experiment context
- `OrientationRelationship`, `TransformationVariant`, `PhaseTransformationRecord`: phase-transformation foundation primitives
- `Rotation`, `RotationSet`, `Orientation`, `Misorientation`, `OrientationSet`: orientation-domain primitives
- `CrystalMap`, `PoleFigure`, `InversePoleFigure`, `ODF`, `HarmonicODF`, `ODFInversionReport`, `HarmonicODFReconstructionReport`: higher-level texture and EBSD models
- `DiffractionGeometry`, `DiffractionPattern`: diffraction-facing models
- `ProvenanceRecord`: source-system and transform traceability
- `ExperimentManifest`, `BenchmarkManifest`, `ValidationManifest`, `WorkflowResultManifest`: stable workflow-interchange and evidence artifacts

## Public API Rule

If a function would be ambiguous without knowledge of frame, symmetry, or basis meaning, it should not accept a naked array as its stable public input.

If a vectorized function would be ambiguous without shared frame, convention, symmetry, or provenance meaning, it should expose or prefer a semantic batch primitive as its stable public input.

If a major stable result must be archived, benchmarked, exchanged with another tool, or consumed by a research agent, it should expose a canonical JSON contract rather than relying on private Python object layout.

PyTex keeps structure-facing and orientation-facing symmetry distinct. `SymmetrySpec` is the point-group-facing surface for direction and orientation reduction. `SpaceGroupSpec` is the structure-facing surface for phase definition and CIF-backed construction.

## Import-Normalization Rule

Every importer must preserve:

- source system and file provenance
- original frame and notation metadata
- canonicalized PyTex representation
- the transform or mapping used to reconcile them

## Structure Import Rule

Crystallographic structure import, including CIF-backed construction, should terminate in canonical `Lattice`, `UnitCell`, and `Phase` objects rather than exposing source-library structure objects as the stable PyTex surface.

## References

### Normative

- [Notation And Conventions](../standards/notation_and_conventions.md)
- [Data Contracts And Manifests](../standards/data_contracts_and_manifests.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- `multimodal_characterization_foundation.md`
