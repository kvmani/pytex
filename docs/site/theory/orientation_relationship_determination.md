# Determining An Orientation Relationship From Measured Orientations

The index-correspondence note fixes what an orientation relationship *does* once its
rotation $\mathbf{R}$ is known. This note fixes how $\mathbf{R}$ is *recovered* from
measured parent and child orientations, and how it is then named and stated. The algorithm is
implemented in `pytex.core.transformation.characterize_orientation_relationship`.

## The Measured Quantity

PyTex stores orientations as crystal-to-specimen matrices in the Bunge convention. For a parent
grain $\mathbf{P}_{i}$ and a child grain $\mathbf{C}_{i}$ formed from it, the canonical
composition is $\mathbf{C} = \mathbf{P}\,\mathbf{V}^{\mathsf{T}}$, so the parent-to-child
rotation that pair exhibits is

$$
\mathbf{V}_{i} = \mathbf{C}_{i}^{\mathsf{T}}\,\mathbf{P}_{i}.
$$ (eq-or-measured-pair)

Equation {eq}`eq-or-measured-pair` has one implementation in the library; no downstream surface
re-derives the placement of the transpose.

## Two Ambiguities, One Group Structure

The $\mathbf{V}_{i}$ are not directly comparable, for two distinct reasons that happen to be
resolved by the same construction.

First, each orientation is defined only up to its own crystal symmetry: replacing
$\mathbf{P}_{i}$ by $\mathbf{P}_{i}\mathbf{S}_{p}$ or $\mathbf{C}_{i}$ by
$\mathbf{C}_{i}\mathbf{S}_{c}$ describes the same grain.

Second, distinct child grains generally form through distinct *variants*. A variant is the
relationship composed with a parent symmetry operation, $\mathbf{V} = \mathbf{R}\,\mathbf{S}_{p}$,
so two grains obeying one relationship can show $\mathbf{V}_{i}$ tens of degrees apart.

Both ambiguities are absorbed by the double coset

$$
G_{c}\,\mathbf{V}_{i}\,G_{p}
= \left\{ \mathbf{S}_{c}\,\mathbf{V}_{i}\,\mathbf{S}_{p} \;:\;
\mathbf{S}_{c} \in G_{c},\ \mathbf{S}_{p} \in G_{p} \right\}
$$ (eq-or-double-coset)

where $G_{p}$ and $G_{c}$ are the parent and child point groups as sets of proper rotations. The
variant operation lies inside $G_{p}$, so {eq}`eq-or-double-coset` contains every description
of every variant of one relationship. This is the fact the whole algorithm rests on. Its size is
$\lvert G_{p}\rvert \lvert G_{c}\rvert$ before deduplication: $576$ for a cubic-to-cubic pair,
$288$ for cubic-to-hexagonal.

## Seeding Without A Nominal Relationship

Earlier fitting surfaces required the caller to supply a nominal relationship as the starting
estimate, which presupposes the answer. Here the estimate is taken from the data: one measured
pair is reduced to its minimum-angle representative,

$$
\mathbf{R}_{0}
= \arg\max_{\mathbf{S}_{c},\,\mathbf{S}_{p}}
\operatorname{tr}\!\left( \mathbf{S}_{c}\,\mathbf{V}_{0}\,\mathbf{S}_{p} \right)
$$ (eq-or-seed)

maximum trace being minimum rotation angle since $\operatorname{tr}\mathbf{R} = 1 + 2\cos\theta$.
Equation {eq}`eq-or-seed` returns the disorientation description of the relationship that pair
shows; every other pair then has an equivalent description near it.

**Why exactly one pair is reduced.** Reducing all pairs independently and averaging the
representatives is incorrect. The maximiser in {eq}`eq-or-seed` is not unique when
$\mathbf{R}$ is itself symmetric under conjugation by part of the group, and different pairs then
select different tied representatives whose mean is a rotation none of them exhibits. The Bain
correspondence is the concrete failure case: $45^{\circ}$ about $\langle 100 \rangle$ with three
variants, whose independently-reduced representatives average to approximately $26.9^{\circ}$,
which is subsequently identified as Kurdjumov–Sachs. Seeding from a single pair and resolving
the remainder against it breaks the ties consistently.

## Symmetry-Aware Rotation Averaging

Two steps alternate to convergence. The alignment step replaces each measurement by the
description nearest the current estimate,

$$
\tilde{\mathbf{V}}_{i}
= \arg\max_{\mathbf{S}_{c},\,\mathbf{S}_{p}}
\operatorname{tr}\!\left( \mathbf{S}_{c}\,\mathbf{V}_{i}\,\mathbf{S}_{p}\,
\mathbf{R}^{\mathsf{T}} \right)
$$ (eq-or-align)

and the averaging step replaces the estimate by the quaternion eigen-mean of the aligned set. With
$\mathbf{q}_{i}$ the unit quaternion of $\tilde{\mathbf{V}}_{i}$, form the scatter matrix

$$
\mathbf{M} = \sum_{i} \mathbf{q}_{i}\,\mathbf{q}_{i}^{\mathsf{T}}
$$ (eq-or-scatter)

and take the eigenvector of largest eigenvalue as the mean quaternion. This is Markley's
attitude average: the maximum-likelihood estimate for small isotropic noise, and, unlike averaging
rotation matrices, it needs no re-orthogonalization because a unit quaternion is a rotation by
construction. Note $\mathbf{q}$ and $-\mathbf{q}$ describe the same rotation and
$\mathbf{M}$ is invariant under that sign, so no sign convention is required.

Convergence is declared when the assignment set in {eq}`eq-or-align` repeats, the mean being a
deterministic function of it, or when the step angle falls below the tolerance. Testing the
assignments rather than the step is what makes the criterion robust to the
$\sim 10^{-6}$ degree matrix-to-quaternion round-trip floor.

Each iteration costs $\mathcal{O}\!\left(n \lvert G_{p}\rvert \lvert G_{c}\rvert\right)$ and is
evaluated as a single contraction over all pairs and both groups; convergence is typically
attained in two to four iterations.

## Naming: Symmetry-Reduced Catalog Distance

The fitted rotation is compared with each candidate under both groups,

$$
d\!\left(\mathbf{R}, \mathbf{R}_{\mathrm{cand}}\right)
= \min_{\mathbf{S}_{c},\,\mathbf{S}_{p}}
\angle\!\left( \mathbf{S}_{c}\,\mathbf{R}\,\mathbf{S}_{p},\ \mathbf{R}_{\mathrm{cand}} \right)
$$ (eq-or-catalog-distance)

the catalog itself being selected from the two crystal systems through one dispatch table. The
separations {eq}`eq-or-catalog-distance` must resolve are fixed crystallography: within the
face-centred to body-centred cubic family the closest pair is Kurdjumov–Sachs and
Greninger–Troiano at $2.404^{\circ}$, followed by Kurdjumov–Sachs to Nishiyama–Wassermann at
$5.264^{\circ}$. That smallest gap sets the usable orientation-noise budget. For the
cubic-to-hexagonal family the two catalog members, Burgers and Shoji–Nishiyama, are separated by
$42.848^{\circ}$, so identification there is far less demanding.

## Stating The Relationship

A rotation has three degrees of freedom. One plane parallelism removes two and one in-plane
direction parallelism removes the third, which is why the classical statement of an orientation
relationship takes exactly that form and why it is complete.

Recovery searches canonical-sign primitive parent triples up to a bound $N$, maps each into the
child basis by the index correspondences of the companion note, and retains a clause when the
angle between the image and a candidate child triple satisfies

$$
\left| \cos \angle\!\left(
\hat{\mathbf{g}}_{\mathrm{image}},\ \hat{\mathbf{g}}_{\mathrm{child}} \right) \right|
\ \ge\ \cos \varepsilon .
$$

The absolute value is correct because the canonical-sign filter has already collapsed each
antiparallel pair to one representative.

**Non-uniqueness of the statement.** A rotation generally satisfies several exact
low-index parallelisms simultaneously. For Kurdjumov–Sachs both $(111)\parallel(011)$ and
$(10\bar{1})\parallel(11\bar{1})$ hold exactly, and index magnitude alone cannot choose between
them. Which clause the literature quotes is determined by the two *structures* — their
close-packed planes and directions — and not by the rotation, which carries no information about
atomic positions. The search therefore accepts a preference: the relationship's own recorded
defining families, or those of the matched catalog member. Fit quality outranks preference in the
ordering, so a nominated family cannot promote a visibly worse clause above an exact one.

## Conclusiveness

An identification is reported as conclusive only when the winning candidate both fits within the
stated tolerance and leads the runner-up by more than the measurement scatter and its own misfit,

$$
\text{margin} \;>\; \max\!\left( \overline{\rho},\ d_{\mathrm{best}} \right)
$$

with $\overline{\rho}$ the mean per-pair residual. Those are precisely the two quantities that
could otherwise account for the lead. On planted Kurdjumov–Sachs data with added Gaussian
orientation scatter the verdict remains conclusive to $2^{\circ}$ of scatter and correctly
becomes inconclusive at $5^{\circ}$, which is comparable to the $2.404^{\circ}$ catalog spacing.

## Assumptions And Limits

- Pairs must be row-matched and share a specimen frame; grain-mean orientations are assumed, and the method does not itself perform grain segmentation.
- The cubic-to-cubic catalog assumes an fcc-to-bcc transformation, because point-group symmetry cannot distinguish an fcc phase from a bcc one. An explicit catalog must be supplied when that assumption fails.
- The parallelism search is bounded; a relationship defined by higher-index parallelisms reports no statement rather than an invented one.
- Validation is synthetic. Measured-EBSD fixtures and a MTEX `calcParent2Child` parity comparison remain outstanding, and no PyTex document claims that parity.

## Normative references

International Tables for Crystallography, Vol. A (point groups and basis conventions).

## Informative references

Markley, F. L., Cheng, Y., Crassidis, J. L., Oshman, Y., J. Guid. Control Dyn. 30 (2007) 1193.
Kurdjumov, G., Sachs, G., Z. Phys. 64 (1930) 325.
Burgers, W. G., Physica 1 (1934) 561.
Morito, S., Tanaka, H., Konishi, R., Furuhara, T., Maki, T., Acta Mater. 51 (2003) 1789.
