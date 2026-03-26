# Notation And Conventions

This document records the initial canonical convention set for PyTex.

## Canonical Internal Defaults

- handedness: right-handed Cartesian frames
- quaternion storage order: `w, x, y, z`
- Euler-angle labeling: Bunge `phi1`, `Phi`, `phi2`
- reciprocal basis normalization: `a*_i dot a_j = delta_ij`
- frame domains: crystal, specimen, map, detector, laboratory, reciprocal

## Explicitness Rules

- Every reference frame must be named and domain-typed.
- Every transform must declare source and target frames.
- Every crystallographic vector or plane representation must state whether it is in direct or reciprocal basis.
- Every imported dataset must preserve source-system provenance and original convention notes.
- Hexagonal and trigonal notation rules are centralized in `hexagonal_and_trigonal_conventions.md`.

## Literature And Tooling Alignment

PyTex will remain explicit about where it aligns with:

- MTEX terminology and validation categories
- ORIX orientation and symmetry workflows
- EBSD vendor frame-conversion conventions
- diffsims and diffraction geometry expectations

When PyTex intentionally differs from another tool's default, the difference must be documented.

## Normative References

- Hahn, Th. (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000100>.
- Hall, S. R. and McMahon, B. (eds.), *International Tables for Crystallography, Volume G: Definition and Exchange of Crystallographic Data*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000107>.
- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- De Graef, M., *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, DOI: <https://doi.org/10.1017/CBO9780511615092>.
