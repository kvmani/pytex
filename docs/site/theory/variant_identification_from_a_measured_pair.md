# Identifying The Variant Of A Measured Child Grain

This note derives the calculation implemented in
`pytex.core.variant_identification.identify_variant` and shown step by step by the workbench
operation *Identify the variant of a measured child grain* (Variants panel). Given the measured
orientations of one parent grain and of one child grain next to it, and an orientation
relationship, it answers: **which variant of that relationship is the child?**

Symbols follow {doc}`../standards/terminology_and_symbol_registry`, section *Variant
identification from one measured pair*. The fitting and naming of an unknown relationship from
many pairs is the subject of {doc}`orientation_relationship_determination`; here the
relationship is given and one pair is explained.

## Purpose And When To Use It

Use it for a handful of hand-picked grain pairs, when the reasoning has to be followed or
documented: teaching, checking a parent reconstruction, or writing the methods section of a
paper. Every intermediate quantity is kept and displayed. For whole maps, the same assignment is
made without the intermediates by `pytex.core.transformation.or_deviation` and
`pytex.core.parent_reconstruction.select_variants`; the tests hold the three to the same answer.

**Inputs.** Two Euler triples in one convention (Bunge or Matthies), the two phases, the
relationship (catalogued or stated as a plane pair and a direction pair), and an assignment
tolerance $\tau$.

**Outputs.** The variant index and the parallelism it realizes, the misfit $\theta^{*}$, the
runner-up and the margin, a verdict, and every matrix of the derivation.

## Step 1–2: Orientation Matrices

PyTex orientations are crystal-to-specimen matrices. For the Bunge (ZXZ) convention

$$
\mathbf{P} = \mathbf{R}_{z}(\varphi_{1})\,\mathbf{R}_{x}(\Phi)\,\mathbf{R}_{z}(\varphi_{2}),
\qquad
\mathbf{R}_{z}(\omega) =
\begin{bmatrix} \cos\omega & -\sin\omega & 0 \\ \sin\omega & \cos\omega & 0 \\ 0 & 0 & 1 \end{bmatrix},
\quad
\mathbf{R}_{x}(\omega) =
\begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos\omega & -\sin\omega \\ 0 & \sin\omega & \cos\omega \end{bmatrix}
$$ (eq-vi-bunge)

and likewise for $\mathbf{C}$. The Matthies (ZYZ) convention replaces the middle factor by
$\mathbf{R}_{y}$. The columns of $\mathbf{P}$ are the parent crystal axes written in specimen
coordinates (Bunge 1982, chapter 2; see {doc}`euler_convention_handling`).

## Step 3: The Relationship As A Rotation

A relationship stated as $(hkl)_{p} \parallel (hkl)_{c}$ with
$[uvw]_{p} \parallel [uvw]_{c}$ fixes one orthonormal triad per phase. With $\hat{\mathbf{n}}$
the unit plane normal and $\hat{\mathbf{d}}$ the unit direction, both in the crystal's Cartesian
frame,

$$
\hat{\mathbf{x}} =
\frac{\hat{\mathbf{d}} - (\hat{\mathbf{d}}\cdot\hat{\mathbf{n}})\,\hat{\mathbf{n}}}
{\lVert \hat{\mathbf{d}} - (\hat{\mathbf{d}}\cdot\hat{\mathbf{n}})\,\hat{\mathbf{n}} \rVert},
\qquad
\hat{\mathbf{y}} = \hat{\mathbf{n}} \times \hat{\mathbf{x}},
\qquad
\mathbf{F} = [\hat{\mathbf{x}}\;\hat{\mathbf{y}}\;\hat{\mathbf{n}}],
\qquad
\mathbf{R} = \mathbf{F}_{c}\,\mathbf{F}_{p}^{\mathsf{T}}.
$$ (eq-vi-triad)

$\mathbf{R}$ maps parent-crystal Cartesian components to child-crystal ones; this is the
construction `OrientationRelationship.from_parallel_plane_direction` uses, and the workbench
shows it only when it reproduces the relationship's stored rotation (it does for every
catalogued plane-and-direction relationship). Projecting $\hat{\mathbf{d}}$ into the plane makes
a slightly inconsistent literature statement still define a proper rotation.

## Step 4: The Measured Rotation

A parent-crystal vector $\mathbf{u}$ reaches the specimen as $\mathbf{P}\mathbf{u}$ and is read
back in the child's axes as $\mathbf{C}^{\mathsf{T}}\mathbf{P}\mathbf{u}$, so

$$
\mathbf{V}_{\mathrm{m}} = \mathbf{C}^{\mathsf{T}}\,\mathbf{P}.
$$ (eq-vi-measured)

$\mathbf{V}_{\mathrm{m}}$ is not yet comparable with $\mathbf{R}$: each measured crystal could
have been described in any of its symmetry-equivalent settings, and each choice changes it.

## Step 5: The Variants

Every parent operator $\mathbf{S}_{p} \in G_{p}$ gives an equivalent description of the parent
and carries the relationship onto a child orientation:

$$
\mathbf{V}_{k} = \mathbf{R}\,\mathbf{S}_{p,k}^{\mathsf{T}}, \qquad \mathbf{S}_{p,k} \in G_{p}.
$$ (eq-vi-variants)

Candidates that differ only by a child operator,
$\mathbf{V}_{j} \in G_{c}\,\mathbf{V}_{k}$, describe the same child crystal and are kept once.
This reproduces the literature counts: 24 for Kurdjumov-Sachs and Greninger-Troiano, 12 for
Nishiyama-Wassermann, Pitsch and Burgers, 3 for Bain (`OrientationRelationship.generate_variants`).

## Step 6: The Misfit Of Each Variant

Variant $k$ predicts the child orientation $\mathbf{C}_{k} = \mathbf{P}\,\mathbf{V}_{k}^{\mathsf{T}}$.
The measured child is the grain $\mathbf{C}\,\mathbf{S}_{c}$ for any $\mathbf{S}_{c} \in G_{c}$,
so the rotation from prediction to measurement is

$$
\mathbf{C}_{k}^{\mathsf{T}}\,\mathbf{C}\,\mathbf{S}_{c}
= \mathbf{V}_{k}\,\mathbf{P}^{\mathsf{T}}\mathbf{C}\,\mathbf{S}_{c}
= \mathbf{V}_{k}\,\mathbf{V}_{\mathrm{m}}^{\mathsf{T}}\,\mathbf{S}_{c},
$$ (eq-vi-residual)

and the misfit of variant $k$ is its smallest rotation angle,

$$
\theta_{k} = \min_{\mathbf{S}_{c} \in G_{c}}
\arccos \frac{\operatorname{tr}\!\left(\mathbf{V}_{k}\,\mathbf{V}_{\mathrm{m}}^{\mathsf{T}}\,\mathbf{S}_{c}\right) - 1}{2}.
$$ (eq-vi-misfit)

Only the child group is minimized over. The parent group is already spent: it is what generated
the variants in {eq}`eq-vi-variants`, and minimizing over it again would make every variant
equivalent and the question meaningless.

**Numerical form.** Near the identity $\operatorname{tr}\Delta\mathbf{g} = 3 - \theta^{2}$ to
leading order, so the $\arccos$ form recovers $\theta$ from a quantity carrying it squared and
loses half its digits. PyTex evaluates

$$
\theta = \operatorname{atan2}\!\left(\tfrac{1}{2}\lVert\operatorname{vec}(\Delta\mathbf{g} -
\Delta\mathbf{g}^{\mathsf{T}})\rVert,\; \tfrac{1}{2}(\operatorname{tr}\Delta\mathbf{g} - 1)\right),
$$ (eq-vi-atan2)

which is the same angle and keeps its digits at zero (`pytex.core._angles`). The workbench prints
both forms for the identified variant.

## Step 7: The Identified Variant

$k^{*} = \arg\min_{k}\theta_{k}$, $\theta^{*} = \theta_{k^{*}}$, and $\mathbf{S}_{c}^{*}$ is the
child operator achieving it. The residual

$$
\Delta\mathbf{g} = \mathbf{V}_{k^{*}}\,\mathbf{V}_{\mathrm{m}}^{\mathsf{T}}\,\mathbf{S}_{c}^{*}
$$ (eq-vi-winner)

is the rotation, in child crystal axes, by which the measured child departs from the prediction.
$\mathbf{S}_{c}^{*}$ is the operator that undoes the setting in which the measured child happened
to be indexed: a child built as $\mathbf{P}\mathbf{V}_{k}^{\mathsf{T}}\mathbf{S}$ is matched with
$\mathbf{S}_{c}^{*} = \mathbf{S}^{\mathsf{T}}$ and $\Delta\mathbf{g} = \mathbf{I}$.

## Step 8: The Parallelism On The Measured Grains

The variant holds the parent member $\mathbf{S}_{p,k^{*}}\hat{\mathbf{n}}_{p}$ parallel to the
child member $\hat{\mathbf{n}}_{c}$; in the measured child's own indices the latter is
$\mathbf{S}_{c}^{*}\hat{\mathbf{n}}_{c}$. Carried into the specimen frame,

$$
\delta = \angle\!\left(\mathbf{P}\,\mathbf{S}_{p,k^{*}}\hat{\mathbf{n}}_{p},\;
\mathbf{C}\,\mathbf{S}_{c}^{*}\hat{\mathbf{n}}_{c}\right) \le \theta^{*},
$$ (eq-vi-parallelism)

and the same for the direction pair. The bound holds because a rotation by $\theta^{*}$ moves no
unit vector by more than $\theta^{*}$; equality holds for a vector normal to the residual axis.
This is the check that needs no symmetry bookkeeping: two measured crystals and two planes.

## Step 9: The Decision

The variant is **assigned** when $\theta^{*} \le \tau$ and every other variant is outside $\tau$;
**ambiguous** when another variant is inside $\tau$ too; and **no variant fits** when
$\theta^{*} > \tau$. The default $\tau = 5^{\circ}$ sits above EBSD orientation noise (about
$0.5^{\circ}$) and above the natural scatter of a transformation about its nominal relationship
(Greninger-Troiano lies $2.40^{\circ}$ from Kurdjumov-Sachs), and below the smallest
Kurdjumov-Sachs intervariant disorientation, $10.53^{\circ}$ (Morito et al. 2003), so two KS
variants can never both be inside it.

## Variant Numbering Is A Convention

Variant indices are PyTex's generation order of {eq}`eq-vi-variants`, not a published table (only
Kurdjumov-Sachs V1 is pinned to Morito et al.). They are also **relative to the parent setting**:
replacing $\mathbf{P}$ by the equivalent $\mathbf{P}\mathbf{S}$ changes
$\mathbf{V}_{\mathrm{m}}$ to $\mathbf{V}_{\mathrm{m}}\mathbf{S}$, and the matching variant
becomes the one generated by $\mathbf{S}\,\mathbf{S}_{p,k^{*}}$. The same child grain is then
reported under a different number while its parallelism {eq}`eq-vi-parallelism`, measured in the
specimen frame, is unchanged. The workbench therefore states every variant by the parallelism it
realizes beside its number, and that statement is the identity to quote. The test
`test_the_number_follows_the_parent_setting_but_the_parallelism_does_not` pins this.

## Failure Modes

- **Wrong convention or setting.** Euler angles exported in another convention, or a hexagonal
  phase indexed with a different Cartesian setting, give a large $\theta^{*}$ and the verdict
  *no variant fits*.
- **Wrong relationship.** A Kurdjumov-Sachs child tested against Nishiyama-Wassermann misses by
  $5.26^{\circ}$, the separation of the two relationships; the warning says so rather than naming
  the nearest variant as an answer.
- **Loose tolerance.** A tolerance above the smallest intervariant angle admits two variants; the
  verdict is then *ambiguous*.

## Validation

`tests/unit/test_variant_identification.py` builds children as
$\mathbf{P}\,\mathbf{V}_{k}^{\mathsf{T}}\mathbf{S}_{c}$ for chosen $k$ and $\mathbf{S}_{c}$, turns
them by known angles, and checks the identified index, the misfit, the matching operator, the
parallelism bound and agreement with `or_deviation`.
`tests/unit/test_app_variant_identification.py` reads the matrices back out of the typeset
equations and compares them with the library's. The worked example
`or-variant-identified-from-a-measured-pair` is in {doc}`../examples/index`.

## References

Normative:

- H.-J. Bunge, *Texture Analysis in Materials Science*, Butterworths (1982), chapter 2.
- S. Morito, H. Tanaka, R. Konishi, T. Furuhara, T. Maki, *Acta Materialia* 51 (2003) 1789.

Informative:

- N. Gey, M. Humbert, *Acta Materialia* 50 (2002) 277.
- C. Cayron, *Journal of Applied Crystallography* 40 (2007) 1183.
- W. G. Burgers, *Physica* 1 (1934) 561.
