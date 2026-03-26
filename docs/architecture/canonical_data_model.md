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

## Canonical Primitives

- `ReferenceFrame`: named right-handed coordinate frame with explicit domain
- `FrameTransform`: reusable rigid transform between frames
- `SymmetrySpec`: point-group or specimen-symmetry description plus operators
- `Basis`: direct or reciprocal basis matrix plus frame information
- `Lattice`: lattice parameters and derived bases
- `UnitCell`: lattice plus atomic basis
- `Phase`: named material phase with lattice, symmetry, crystal frame, and optional space-group or formula metadata
- `Rotation`, `Orientation`, `Misorientation`, `OrientationSet`: orientation-domain primitives
- `CrystalMap`, `PoleFigure`, `InversePoleFigure`, `ODF`: higher-level texture and EBSD models
- `DiffractionGeometry`, `DiffractionPattern`: diffraction-facing models
- `ProvenanceRecord`: source-system and transform traceability

## Public API Rule

If a function would be ambiguous without knowledge of frame, symmetry, or basis meaning, it should not accept a naked array as its stable public input.

## Import-Normalization Rule

Every importer must preserve:

- source system and file provenance
- original frame and notation metadata
- canonicalized PyTex representation
- the transform or mapping used to reconcile them

## Structure Import Rule

Crystallographic structure import, including CIF-backed construction, should terminate in canonical `Lattice`, `UnitCell`, and `Phase` objects rather than exposing source-library structure objects as the stable PyTex surface.
