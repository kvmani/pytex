# Determining An Orientation Relationship From Measured Orientations

This note derives the mathematical formulation and optimization framework implemented in
`pytex.core.transformation.characterize_orientation_relationship` to determine orientation
relationships from experimental parent and child orientation data. It covers double-coset symmetry
quotienting, weighted quaternion eigen-mean refinement,
automated classification against classical cataloged relationships, and recovery of rational
plane and direction parallelisms.

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
pair with positive weight is reduced to its minimum-angle representative,

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
\mathbf{Q}_{\mathrm{OR}} = \sum_{i} w_i\,\mathbf{q}_{i}\,\mathbf{q}_{i}^{\mathsf{T}}
$$ (eq-or-scatter)

and take the eigenvector of largest eigenvalue as the mean quaternion. The non-negative scalar
weights sum to one; equal weights are the default. This is Markley's attitude average,
minimizing the weighted squared Frobenius distance between aligned rotation matrices. It is
a chordal objective, rather than the squared geodesic-angle objective of a Karcher mean.
For small isotropic errors, inverse-variance weights motivate a local statistical interpretation;
arbitrary reliability weights do not establish that interpretation. Unlike averaging
rotation matrices, it needs no re-orthogonalization because a unit quaternion is a rotation by
construction. Note $\mathbf{q}$ and $-\mathbf{q}$ describe the same rotation and
$\mathbf{Q}_{\mathrm{OR}}$ is invariant under that sign, so no sign convention is required.

Convergence is declared when the assignment set in {eq}`eq-or-align` repeats, the mean being a
deterministic function of it, or when the step angle falls below the tolerance. Testing the
assignments also detects a fixed point directly. Step and residual angles use the shared
skew/trace `atan2` recovery, avoiding the loss of precision of `arccos` near zero.

Each iteration costs $\mathcal{O}\!\left(n \lvert G_{p}\rvert \lvert G_{c}\rvert\right)$ and is
evaluated in bounded blocks over pairs and both groups. Only one aligned rotation per pair
survives each block, so the symmetry-expanded temporary is independent of dataset size.
Final residuals are realigned to the final estimate, including when the iteration limit is reached.

## Measurement Weights And Excluded Pairs

Use `pair_weights` on `fit_orientation_relationship`, `characterize_orientation_relationship`,
or `orientation_relationship_from_euler` when there is an independent reason to give pairs
different influence. A zero weight excludes a pair from estimation and seed selection; it
does not remove the pair from `residuals_deg`, the all-pair mean, or the maximum residual.
The report owns a normalized, read-only copy of the weights.

The weighted mean residual is the sum of each residual multiplied by its normalized weight.
The effective pair count is the reciprocal of the sum of squared normalized weights. It
measures concentration of evidence: equal weights give the supplied pair count, while one
nonzero weight gives one. Repeated observations from one grain are still correlated; this
quantity is neither an independence test nor a confidence interval. Inspect both weighted
and all-pair residuals, and retain the reason for exclusions with the experimental record.

For two rotations about one axis with angles zero and a known angle, the fitted angle is
the argument of the weighted sum of their unit complex numbers. This independent circular
identity tests the quaternion implementation without comparing it against an earlier run.
The [weighted worked examples](../examples/generated/weighted-or-fitting.md) compute this identity
and demonstrate an excluded pair retained in the diagnostics.

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

## From A Statement To An Object, And What That Costs

A recovered statement is prose. `ORCharacterizationReport.as_rational_relationship` turns it into a
genuine `OrientationRelationship`, built from the integer plane pair and the integer direction pair
by `from_parallel_plane_direction` — the same constructor a catalog relationship is built with, so
the result behaves like one in every downstream call.

That construction is an **idealization**, and the returned `RationalizedORResult` prices it. The
integer statement names a nearby *exact* relationship $\mathbf{R}_{\mathrm{ideal}}$; the measurement
is $\mathbf{R}_{\mathrm{fit}}$; and the cost is the symmetry-reduced angle between them,

$$
\Delta\omega \;=\; \min_{S_c \in G_c,\ S_p \in G_p}
\angle\!\left( S_c\, \mathbf{R}_{\mathrm{ideal}}\, S_p,\ \mathbf{R}_{\mathrm{fit}} \right),
$$

reported beside the per-clause deviations. **A rational OR handed back without that number reads as
a measurement of the relationship it was rounded to.** The number is meaningful against one
reference: the scatter of the data it came from. An idealization costing less than the scatter is
one the measurement cannot distinguish; one costing several times the scatter is a claim the
measurement contradicts, however tidy the integers look.

The trade is real and visible. Six exact Greninger–Troiano pairs, held to $|index| \le 2$,
rationalize to the *Kurdjumov–Sachs* statement — $\{111\}\parallel\{110\}$ with
$\langle 110 \rangle \parallel \langle 111 \rangle$ — at a cost of $2.404^{\circ}$, which is exactly
the catalog spacing between the two. Allowing index 3 finds $[321]\parallel[223]$ at $0.44^{\circ}$,
and index 4 finds a statement at $0.02^{\circ}$. Tidier integers cost more, and the caller chooses.

**The zone law is a constraint, not a refinement.** `from_parallel_plane_direction` removes the
component of the direction along the plane normal, so a plane and a direction that is not in it
would build a relationship the two printed labels do not describe. The direction clause is
therefore selected from those satisfying $\mathbf{h}\cdot\mathbf{u} = 0$ against the chosen plane,
and the result carries the resulting angular departure so a reader can check it is zero. When no
clause in the chosen plane exists within the bounds, the method refuses with a message
distinguishing that case from having found no plane clause at all: the two ask for different
responses.

## Conclusiveness

An identification is reported as conclusive only after convergence, when both the weighted
mean scatter and the best catalog distance are within the naming tolerance. With more than
one catalog candidate the winner must also lead the runner-up by more than the weighted
scatter and its own misfit,

$$
\text{margin} \;>\; \max\!\left( \overline{\rho},\ d_{\mathrm{best}} \right)
$$

with $\overline{\rho}$ the weighted mean per-pair residual. Those are precisely the two quantities that
could otherwise account for the lead. On planted Kurdjumov–Sachs data with added Gaussian
orientation scatter the verdict remains conclusive to $2^{\circ}$ of scatter and correctly
becomes inconclusive at $5^{\circ}$, which is comparable to the $2.404^{\circ}$ catalog spacing.

## Assumptions And Limits

- Pairs must be row-matched and share a specimen frame; grain-mean orientations are assumed, and the method does not itself perform grain segmentation.
- The cubic-to-cubic catalog assumes an fcc-to-bcc transformation, because point-group symmetry cannot distinguish an fcc phase from a bcc one. An explicit catalog must be supplied when that assumption fails.
- The parallelism search is bounded; a relationship defined by higher-index parallelisms reports no statement rather than an invented one.
- A rationalized relationship is an idealization and is always returned with its cost; it is named with a `_rationalized` suffix so that a later report cannot mistake it for the measurement.
- Validation is synthetic. Measured-EBSD fixtures and a MTEX `calcParent2Child` parity comparison remain outstanding, and no PyTex document claims that parity.

## Normative references

International Tables for Crystallography, Vol. A (point groups and basis conventions).

## Informative references

Markley, F. L., Cheng, Y., Crassidis, J. L., Oshman, Y., *Averaging Quaternions*,
J. Guid. Control Dyn. 30 (2007) 1193–1197, [doi:10.2514/1.28949](https://doi.org/10.2514/1.28949).
[Author manuscript at NASA](https://ntrs.nasa.gov/citations/20070017872).
Kurdjumov, G., Sachs, G., Z. Phys. 64 (1930) 325.
Burgers, W. G., Physica 1 (1934) 561.
Morito, S., Tanaka, H., Konishi, R., Furuhara, T., Maki, T., Acta Mater. 51 (2003) 1789.
The symbols follow the {doc}`central registry </standards/terminology_and_symbol_registry>`.
