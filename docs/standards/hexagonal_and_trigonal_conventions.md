# Hexagonal And Trigonal Conventions

This document fixes the initial PyTex policy for hexagonal and trigonal indexing so algorithms and importers do not drift.

## Canonical Rule

- PyTex stores internal direct-space vectors as three-component coordinates in the chosen crystal direct basis.
- PyTex stores internal reciprocal-space plane normals as three-component coordinates in the reciprocal basis.
- Four-index notation is accepted and documented at the interface boundary, but normalized into canonical internal forms immediately.

## HCP Crystal Frame

PyTex fixes the conventional HCP crystal frame as:

- `a1` along the positive crystal `x` direction
- `a2` in the basal plane at `120°` from `a1`
- `c` along the positive crystal `z` direction
- right-handed Cartesian embedding for the crystal frame

This is the frame assumed by the canonical HCP figure and by the conversion helpers for basal-plane notation.

## Direction Notation

For hexagonal and trigonal lattices, PyTex recognizes:

- three-index direct-basis directions `[u v w]`
- four-index Weber or Miller-Bravais directions `[U V T W]`, with `U + V + T = 0`

PyTex adopts the De Graef / ORIX conversion convention:

```text
U = (2u - v) / 3
V = (2v - u) / 3
T = -(u + v) / 3
W = w
```

At the human-facing notation layer, integer Miller-Bravais indices are obtained by clearing denominators and reducing to the smallest common integer set.

## Plane Notation

For hexagonal reciprocal-space planes, PyTex recognizes:

- three-index reciprocal-basis planes `(h k l)`
- four-index Miller-Bravais planes `(h k i l)`, with `i = -(h + k)`

Internally, PyTex stores the reciprocal-space form as three components and preserves the original four-index notation in provenance or importer metadata where relevant.

## Why This Is Centralized

- The direct and reciprocal hexagonal bases are rotated by `30°` relative to one another in the basal plane.
- Four-index notation restores manifest sixfold symmetry in the basal plane, but it is redundant.
- Internal storage and human-facing notation should therefore not be conflated.

## Normative References

- De Graef, M., *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, 2003, DOI: <https://doi.org/10.1017/CBO9780511615092>.
- Hahn, Th. (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000100>.

## Informative References

- ORIX documentation, “Crystal directions”, for an explicit software expression of the De Graef convention: <https://orix.readthedocs.io/en/stable/tutorials/crystal_directions.html>.
- Oak Ridge single-crystal diffraction notes on the hexagonal lattice, for geometric explanation of the `30°` direct/reciprocal rotation: <https://single-crystal.ornl.gov/more/hexagonal/index.html>.
