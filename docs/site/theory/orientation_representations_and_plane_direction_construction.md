# Orientation Representations and Plane–Direction Construction in PyTex

## Scope

This note records two related PyTex foundations:

- vectorized conversion among rotation matrices, unit quaternions, axis–angle pairs, and Rodrigues / Rodrigues–Frank coordinates,
- construction of crystal-to-specimen orientations from a crystal plane and a crystal direction.

The implementation is intentionally explicit about mapping direction, singularities, and default specimen references so the batch APIs remain scientifically auditable.

The *equal-volume* representations — homochoric and cubochoric — are not covered here; they, the invariant-measure argument that motivates them, and the inverse construction that recovers $(hkl)[uvw]$ indices from an orientation are in 	exttt{docs/tex/theory/orientation_representations.tex}.

## Canonical Rotation Representation Rules

PyTex stores unit quaternions in $(w, x, y, z)$ order and interprets them as active crystal-to-specimen rotations. For an axis–angle pair $(\hat{\mathbf{n}}, \omega)$ the quaternion is

$$
q = \left(\cos\frac{\omega}{2},\ \hat{\mathbf{n}}\sin\frac{\omega}{2}\right)
$$

For batch export to axis–angle form, PyTex canonicalizes the quaternion sign so $w \ge 0$ before extracting the angle. This yields

$$
\omega \in [0, \pi]
$$

with a default axis of $\hat{\mathbf{e}}_z$ for the identity rotation. That default axis is a storage convention only; the identity rotation has no unique physical axis.

## Rodrigues and Rodrigues–Frank Coordinates

PyTex uses the Rowenhorst et al.\ definition of Rodrigues coordinates:

$$
\boldsymbol{\rho} = \hat{\mathbf{n}}\tan\frac{\omega}{2}
$$

The Rodrigues–Frank form stores the same information as a four-component vector

$$
\left(\hat{\mathbf{n}}, \tan\frac{\omega}{2}\right)
$$

This makes the two-fold singularity explicit. When $\omega = \pi$, the Rodrigues vector magnitude tends to infinity and the Rodrigues–Frank scale is represented as $+\infty$ in PyTex. That representation is deliberate: it preserves the mathematical singularity instead of hiding it behind an arbitrary clip value.

## Plane–Direction Orientation Construction

Given a crystal plane normal $\hat{\mathbf{n}}_c$ and a crystal direction $\mathbf{d}_c$, PyTex first projects the direction into the plane:

$$
\hat{\mathbf{x}}_c =
\frac{\mathbf{d}_c - (\mathbf{d}_c \cdot \hat{\mathbf{n}}_c)\hat{\mathbf{n}}_c}
{\left\lVert \mathbf{d}_c - (\mathbf{d}_c \cdot \hat{\mathbf{n}}_c)\hat{\mathbf{n}}_c \right\rVert}
$$

The corresponding right-handed crystal basis is

$$
\hat{\mathbf{y}}_c = \hat{\mathbf{n}}_c \times \hat{\mathbf{x}}_c, \qquad
\hat{\mathbf{z}}_c = \hat{\mathbf{n}}_c
$$

PyTex constructs an analogous specimen basis from the requested specimen plane normal and specimen in-plane direction. The current default is:

- crystal plane normal $\rightarrow$ specimen $Z$,
- in-plane crystal direction $\rightarrow$ specimen $X$.

The rotation matrix is then

$$
\mathbf{R}_{c \rightarrow s} =
\mathbf{B}_s \mathbf{B}_c^{\mathsf{T}}
$$

where $\mathbf{B}_c = [\hat{\mathbf{x}}_c\ \hat{\mathbf{y}}_c\ \hat{\mathbf{z}}_c]$ and
$\mathbf{B}_s = [\hat{\mathbf{x}}_s\ \hat{\mathbf{y}}_s\ \hat{\mathbf{z}}_s]$.

This guarantees

$$
\mathbf{R}_{c \rightarrow s}\hat{\mathbf{n}}_c = \hat{\mathbf{n}}_s,
\qquad
\mathbf{R}_{c \rightarrow s}\hat{\mathbf{x}}_c = \hat{\mathbf{x}}_s
$$

and preserves the PyTex orientation convention that orientations map crystal vectors into specimen vectors.

## Failure Modes

PyTex rejects the construction if:

- the supplied plane and direction do not share a phase,
- the projected crystal direction becomes zero,
- the projected specimen direction becomes zero,
- the resulting triads do not define proper right-handed bases.

These checks are performed at construction time so invalid frame or phase semantics cannot survive into later EBSD, PF, IPF, or ODF workflows.

## Normative References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, 1969. DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- D. Rowenhorst et al., “Consistent Representations of and Conversions Between 3D Rotations”, *Modelling and Simulation in Materials Science and Engineering*, 23(8), 2015. DOI: <https://doi.org/10.1088/0965-0393/23/8/083501>.
