# Fundamental Region Reduction

PyTex currently distinguishes between two reduction problems:

- reducing crystal directions into class-specific inverse-pole-figure sectors
- reducing orientations into a symmetry-equivalent representative in orientation space

## Direction Reduction

Let $\mathbf{d}_c$ be a unit crystal direction and let $\mathcal{C}$ denote the crystal symmetry group. PyTex forms the symmetry orbit

$$
\mathcal{O}(\mathbf{d}_c) = \{ c \mathbf{d}_c \; | \; c \in \mathcal{C} \}
$$

Direction reduction is then implemented by deterministic selection of the representative lying inside the supported direction-space fundamental sector

$$
\mathbf{d}_{\mathrm{FR}} \in \mathcal{F}_{\mathrm{dir}}
$$

This is the reduction used for inverse pole figures and IPF color keys.

## Orientation Reduction

Let $g$ denote an orientation mapping crystal coordinates into specimen coordinates. With crystal symmetry $\mathcal{C}$ and optional specimen symmetry $\mathcal{S}$, PyTex forms the orbit

$$
\mathcal{O}(g) = \{ s g c \; | \; s \in \mathcal{S},\; c \in \mathcal{C} \}
$$

PyTex currently supports two selection modes.

**Exact orbit reduction in the quaternion hemisphere.**

If no reference orientation is supplied, each candidate is canonicalized into the unit-quaternion hemisphere and PyTex selects the symmetry-equivalent representative with maximal scalar part:

$$
q_0 = \max_{q' \in \mathcal{O}(q)} q'_0
$$

Because the rotation angle is

$$
\omega = 2 \arccos(q_0)
$$

on that hemisphere, this is exactly the minimum-angle representative relative to the identity. PyTex then applies a deterministic lexicographic tie-break for orbit points lying on symmetry bisectors.

**Reference-aware projection.**

If a reference orientation $g_{\mathrm{ref}}$ is supplied, PyTex selects the candidate minimizing the unsymmetrized orientation distance to that reference.

## Current Limits

- PyTex now implements the exact orbit-reduction rule for the supported proper point groups already generated in the codebase, and the implementation is regression-tested class by class across that supported set.
- What remains ahead is a class-by-class closed-form boundary catalog and broader external parity fixtures for those exact boundaries.
- Current parity hardening therefore treats orientation projection as scientifically strong, but not yet exhaustively benchmarked against every external boundary reference.
