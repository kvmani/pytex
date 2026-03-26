# Notation And Conventions

## Canonical Internal Defaults

- handedness: right-handed Cartesian frames
- quaternion storage order: `w, x, y, z`
- Euler-angle labeling: Bunge `phi1`, `Phi`, `phi2`
- reciprocal basis normalization: `a*_i dot a_j = delta_ij`
- frame domains: crystal, specimen, map, detector, laboratory, reciprocal

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

No subsystem may invent additional stable frame domains or silently collapse these distinctions.

## Explicitness Rules

- Every reference frame must be named and domain-typed.
- Every transform must declare source and target frames.
- Every crystallographic vector or plane representation must state whether it is in direct or reciprocal basis.
- Every imported dataset must preserve source-system provenance and original convention notes.
- Hexagonal and trigonal notation rules are centralized in `hexagonal_and_trigonal_conventions.md`.

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
