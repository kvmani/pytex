# Notation And Conventions

## Canonical Internal Defaults

- handedness: right-handed Cartesian frames
- quaternion storage order: `w, x, y, z`
- Euler-angle labeling: Bunge `phi1`, `Phi`, `phi2`
- reciprocal basis normalization: `a*_i dot a_j = delta_ij`
- frame domains: crystal, specimen, map, detector, laboratory, reciprocal
- units at stable boundaries: angstrom for direct lengths, inverse angstrom for
  reciprocal lengths, degrees for angles (radians only where a name ends in `_rad`)
- point-group and space-group naming: Hermann-Mauguin symbols

## Crystallographic Notation

PyTex writes crystallographic quantities the way the international literature writes them. These
rules are implemented once in `pytex.core.notation` and must not be re-implemented inline: a
repository that formats indices in three places will eventually print three different spellings of
the same plane.

### The reciprocal star marks the basis, not the indices

This is the rule most easily got wrong, so it is stated explicitly.

| Quantity | Written | Starred |
| --- | --- | --- |
| direct basis vectors | **a**, **b**, **c** | no |
| **reciprocal basis vectors** | **a\***, **b\***, **c\*** | **yes** |
| axis labels of a reciprocal-domain frame | `a*`, `b*`, `c*` | **yes** |
| reciprocal lattice parameters | `a*`, `b*`, `c*`, `alpha*`, `beta*`, `gamma*` | **yes** |
| Miller plane indices | `(hkl)` | **no** |
| lattice direction indices | `[uvw]` | **no** |
| reciprocal-lattice vector | **g**\_hkl = h**a\*** + k**b\*** + l**c\*** | bold `g`, starred basis |

Miller indices are *already* reciprocal-basis components by definition, so starring them would
name a different quantity. Any surface quoting reciprocal-basis components must instead **say** the
basis is reciprocal — which is what the starred axis labels and `BasisKind.RECIPROCAL` accomplish.

Use `format_reciprocal_axis_label` / `format_reciprocal_axis_labels` rather than appending `"*"` by
hand; the helpers are idempotent, so a label cannot acquire `a**` by passing through two layers.

### Specific quantities versus symmetry families

| | specific | symmetry family |
| --- | --- | --- |
| lattice plane | `(hkl)` | `{hkl}` |
| lattice direction | `[uvw]` | `<uvw>` |
| Miller-Bravais plane | `(hkil)` | `{hkil}` |
| Miller-Bravais direction | `[uvtw]` | `<uvtw>` |

The distinction is scientific, not cosmetic. A pole figure, a powder reflection, a slip-system
family, and any multiplicity-bearing quantity denote the whole symmetry-related orbit and take the
family brackets. Naming a single member where the science means the orbit misstates the quantity.
Where an object can be either — a `PoleFigure` can be built with or without family expansion — the
object records which it is and the notation follows the record rather than an assumption.

### Negative indices and separators

- A negative index is written with an **overbar**, not a leading minus, in publication-facing
  output. The `"mathtext"` style emits the overbar; the `"plain"` style falls back to `-1`
  because a terminal has no overbar.
- Components are concatenated only when that is unambiguous. `(110)` is the classical form, but
  `[1-10]` could be read as `[1, -1, 0]` or `[1, -10]`, and `(1210)` as `(1, 2, 1, 0)` or
  `(12, 1, 0)`. A separator is therefore inserted whenever a component is negative in plain style
  or any component has more than one digit.

### Zone law

A direction `[uvw]` lies in a plane `(hkl)` when `h u + k v + l w = 0`. This is the statement
relating the two index families and the basis of zone-axis reasoning in the diffraction module.

### Four-index hexagonal forms

Miller-Bravais forms carry the redundancy constraints `h + k + i = 0` and `u + v + t = 0`. The
conversion rules and the direct/reciprocal basal-plane rotation are centralized in
`hexagonal_and_trigonal_conventions.md`.

## Canonical Frame Chain

PyTex uses one repository-wide frame-chain doctrine:

`crystal -> specimen -> map -> detector -> laboratory -> reciprocal`

The arrows do not imply that every workflow instantiates every frame. They state the admissible vocabulary and the order in which PyTex expects those domains to be related when a workflow spans them.

- `crystal -> specimen`
  Normative by PyTex adoption of the orientation-as-crystal-to-specimen mapping used throughout the core model.
- `specimen -> map`
  Normative by PyTex architectural rule: map coordinates are not assumed identical to specimen coordinates unless a workflow declares that relationship.
- `specimen -> detector -> laboratory`
  Normative by PyTex diffraction-geometry contracts: detector and laboratory semantics must remain explicit and separate from specimen semantics.
- `crystal -> reciprocal`
  Normative from IUCr-style crystallographic basis duality and the PyTex reciprocal normalization rule.
- `crystal or reciprocal -> visualization`
  Normative by PyTex plotting doctrine: visualization views may re-express already defined scientific geometry, but they do not define new crystallographic meaning or replace crystal, detector, or reciprocal frames.

No subsystem may invent additional stable frame domains or silently collapse these distinctions.

## Explicitness Rules

- Every reference frame must be named and domain-typed.
- Every transform must declare source and target frames.
- Every stable batch primitive must declare the shared frame or convention metadata required to interpret the batch.
- Every crystallographic vector or plane representation must state whether it is in direct or reciprocal basis.
- Every diffraction plot must make clear whether plotted coordinates live in detector, reciprocal, or angular coordinates.
- Every 3D crystal visualization must state whether camera alignment is specified by a visualization view only or by a crystallographic direction such as a `CrystalDirection`.
- Every imported dataset must preserve source-system provenance and original convention notes.
- Hexagonal and trigonal notation rules are centralized in `hexagonal_and_trigonal_conventions.md`.
- Repository-wide symbol names and glossary terms are centralized in `terminology_and_symbol_registry.md`.

## Frame Geometry And The Standard Catalog

- Every `ReferenceFrame` carries the components of its three labelled axes in the **canonical
  right-handed Cartesian reference** `X, Y, Z`. The default is the identity triad, meaning the
  frame's axes coincide with the canonical Cartesian axes.
- Frame axis vectors are dimensionless axis *orientation*. Physical axis lengths belong to
  `Basis`, which carries a `BasisKind` and a unit; the two must not be conflated.
- A declared handedness must agree with the sign of the axis-vector determinant; a contradiction
  is a construction-time error, not a downstream surprise.
- Non-orthonormal frame axes are permitted (an oblique crystal frame is legitimate) and reported
  through `is_orthonormal` rather than rejected.
- The standard frames are built once in `pytex.core.frame_catalog` and reused: `cartesian`,
  `specimen`, `sample` (`RD, TD, ND`), `crystal`, `map`, `detector`, `laboratory`, plus
  `reciprocal_frame_for(...)`. Subsystems use the catalog rather than constructing frames inline,
  so frame identity is stable across module boundaries.
- A `FrameTransform` maps **components in its source frame to components in its target frame**
  (`v_target = R v_source + t`). Directions, plane normals, and poles are translation-invariant and
  must use the direction-only application path.
- Workflows spanning more than one declared relationship should register their transforms in a
  `FrameGraph` rather than composing chains by hand; the graph composes the shortest declared
  chain.

See [Reference Frame Foundation](../architecture/reference_frame_foundation.md) for the full model.

## Batch Primitive Rule

PyTex treats vectorized semantics as first-class scientific meaning.

- vector batches that share one frame should use a frame-aware primitive such as `VectorSet`
- Euler batches should keep convention metadata through `EulerSet`
- quaternion batches should keep canonical normalization and storage semantics through `QuaternionSet`
- rotation batches should use `RotationSet`
- orientation batches should use `OrientationSet`

Raw arrays may still be accepted at boundaries for compatibility, but the stable semantic model prefers named batch primitives whenever the batch meaning would otherwise be implicit.

## Literature And Tooling Alignment

PyTex aims to remain compatible in meaning, not merely in syntax, with:

- MTEX terminology and validation categories
- ORIX orientation and symmetry workflows
- EBSD vendor frame-conversion conventions
- diffsims and diffraction geometry expectations

## References

### Normative

- Hahn, Th. (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000100>.
- Hall, S. R. and McMahon, B. (eds.), *International Tables for Crystallography, Volume G: Definition and Exchange of Crystallographic Data*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000107>.
- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- `reference_canon.md`

### Informative

- De Graef, M., *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, DOI: <https://doi.org/10.1017/CBO9780511615092>.
- Nolze et al., *Journal of Applied Crystallography* (2023), DOI: <https://doi.org/10.1107/S1600576723009275>.
- MTEX documentation: <https://mtex-toolbox.github.io/>
