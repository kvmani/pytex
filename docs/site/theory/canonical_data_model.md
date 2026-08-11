# Canonical Data Model for PyTex

## Motivation

Texture and diffraction workflows are sensitive to convention drift. A vector without an attached frame, a plane without a declared reciprocal basis, or an orientation without explicit crystal and specimen frames is scientifically under-specified.

## Design

PyTex therefore defines stable primitives for:

- reference frames and transforms,
- symmetry specifications,
- lattice and basis representations,
- rotations, orientations, and misorientations,
- texture-domain containers,
- EBSD-domain containers,
- diffraction-domain containers,
- provenance records.

## Canonical Convention

The initial canonical convention set is:

\begin{align*}
\text{handedness} &= \text{right-handed}, \\
  q &= (w, x, y, z), \\
  \text{Euler labels} &= (\phi_1, \Phi, \phi_2), \\
  \mathbf{a}_i^\ast \cdot \mathbf{a}_j &= \delta_{ij}.
\end{align*}

## Implementation Rule

Public APIs should not expose raw arrays where domain types are needed to recover scientific meaning.
