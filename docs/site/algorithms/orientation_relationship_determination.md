# Determining An Orientation Relationship From Measured Orientations

**Surface:** `pytex.core.transformation.characterize_orientation_relationship`,
`orientation_relationship_from_euler`, `describe_orientation_relationship`.

Given parent and child grain orientations measured by EBSD, this recovers the
operative parent-to-child rotation, names it against the standard catalog, states
it as parallel planes and directions, and says whether the identification can be
trusted. No nominal relationship has to be supplied — the starting estimate comes
from the data.

```{figure} ../../figures/or_determination_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 turns paired parent and child orientations
  into per-pair rotations. Lane 2 seeds the fit from the double-coset reduction
  of one pair and refines it by alternating symmetry alignment with the
  quaternion eigen-mean. Lane 3 ranks the catalog, extracts the parallelism
  statement, and emits a report with a conclusiveness verdict.
:width: 100%

The algorithm, with the constraint governing each stage.
```

## 1. What is being computed

An orientation relationship is a rotation $\mathbf{R}$ carrying a vector
expressed in the parent crystal's Cartesian frame into the child crystal's
frame. A measurement gives orientations, not $\mathbf{R}$: for pair $i$ the
parent's crystal-to-specimen matrix $\mathbf{P}_i$ and the child's
$\mathbf{C}_i$. In the canonical convention
$\mathbf{C} = \mathbf{P}\mathbf{V}^{\mathsf{T}}$, so the rotation that pair shows
is

$$\mathbf{V}_i = \mathbf{C}_i^{\mathsf{T}} \mathbf{P}_i .$$

This expression has exactly one definition in the library
(`_measured_parent_to_child`); no call site re-derives where the transpose goes.

The $\mathbf{V}_i$ are not directly comparable. Each orientation is defined only
up to its own crystal symmetry, and different grains formed through different
**variants** of the same relationship, which differ by a parent symmetry
operation. Two pairs obeying one relationship can therefore produce
$\mathbf{V}_i$ tens of degrees apart. Recovering the relationship means undoing
both ambiguities.

## 2. Symmetry, and what it does to the problem

Write $G_p$ and $G_c$ for the parent and child point groups as sets of proper
rotations. Every description of one measurement is

$$\mathbf{S}_c\,\mathbf{V}_i\,\mathbf{S}_p, \qquad
\mathbf{S}_c \in G_c,\; \mathbf{S}_p \in G_p ,$$

the **double coset** $G_c \mathbf{V}_i G_p$. Its size is the number of
descriptions the algorithm must choose among for each pair:

| system | example | $\lvert G_p \rvert$ | $\lvert G_c \rvert$ | descriptions per pair |
| --- | --- | --- | --- | --- |
| cubic → cubic | Kurdjumov-Sachs, fcc → bcc | 24 | 24 | 576 |
| cubic → hexagonal | Burgers, bcc → hcp | 24 | 12 | 288 |

The variant operation lives *inside* this coset: $\mathbf{V}_i = \mathbf{R}\,
\mathbf{S}_{p,i}$, and $\mathbf{S}_{p,i} \in G_p$. That is the fact the whole
algorithm rests on — absorbing the coset absorbs the variant, so pairs formed
through different variants become comparable.

## 3. The starting estimate, without a nominal relationship

Reduce **one** pair to its minimum-angle representative in its double coset:

$$\mathbf{R}_0 = \arg\max_{\mathbf{S}_c,\, \mathbf{S}_p}
\operatorname{tr}\!\left( \mathbf{S}_c\,\mathbf{V}_0\,\mathbf{S}_p \right),$$

maximum trace being minimum rotation angle, since
$\operatorname{tr}\mathbf{R} = 1 + 2\cos\theta$. This is the *disorientation*
description of the relationship that pair shows. Every other pair has an
equivalent description close to it, which the next step finds.

:::{admonition} Constraint: reduce one pair, not all of them
:class: warning

Reducing every pair independently and averaging the results looks more robust
and is wrong. The maximum-trace element is **not unique** when the
relationship's own rotation is symmetric, so different pairs land on different
tied representatives and their mean is a rotation none of them shows.

Bain is the concrete failure: $45^{\circ}$ about $\langle 100 \rangle$ with three
variants averages to a meaningless $26.9^{\circ}$, which then reads as
Kurdjumov-Sachs. Seeding from one pair and resolving the rest against it breaks
the ties consistently, and `test_bain_survives_the_double_coset_tie` fails if
that regresses.
:::

## 4. Refinement: align, average, iterate

Two steps alternate until the assignment stops changing.

**Align.** Replace each measurement by the description nearest the current
estimate $\mathbf{R}$:

$$\tilde{\mathbf{V}}_i = \arg\max_{\mathbf{S}_c,\, \mathbf{S}_p}
\operatorname{tr}\!\left( \mathbf{S}_c\,\mathbf{V}_i\,\mathbf{S}_p\,
\mathbf{R}^{\mathsf{T}} \right).$$

**Average.** Take the quaternion eigen-mean (Markley) of the aligned set: with
$\mathbf{q}_i$ the unit quaternion of $\tilde{\mathbf{V}}_i$, the mean rotation
is the eigenvector of largest eigenvalue of

$$\mathbf{M} = \sum_i \mathbf{q}_i \mathbf{q}_i^{\mathsf{T}} .$$

This is the maximum-likelihood rotation average for small isotropic noise, and
unlike averaging matrices it needs no re-orthogonalization.

Convergence is declared when the alignment **assignments** repeat — the mean is
then a deterministic function of them — or when the step falls below the angular
tolerance. Testing the assignments rather than the step is what makes it robust
to the $\sim 10^{-6}$ degree matrix-to-quaternion round-trip floor.

Cost is $\mathcal{O}\!\left(n \lvert G_p \rvert \lvert G_c \rvert\right)$ per
iteration, evaluated as a single `einsum` over all pairs and both groups at once,
and convergence is typically two to four iterations.

### Worked behaviour: the fit averages noise

Twelve child grains planted through mixed Kurdjumov-Sachs variants of one
parent, with Gaussian orientation scatter added per grain:

| added scatter | catalog winner | deviation of the fit | pair scatter | conclusive |
| --- | --- | --- | --- | --- |
| 0.0° | Kurdjumov-Sachs | 0.000° | 0.000° | yes |
| 0.5° | Kurdjumov-Sachs | 0.166° | 0.368° | yes |
| 2.0° | Kurdjumov-Sachs | 0.745° | 1.716° | yes |
| 5.0° | Kurdjumov-Sachs | 0.469° | 4.011° | **no** |

The fitted rotation sits consistently *closer* to the true relationship than the
individual pairs do — that gap is the averaging. At 5° the verdict correctly
degrades to inconclusive, because the scatter is comparable to the 2.40°
separating Kurdjumov-Sachs from Greninger-Troiano.

## 5. Naming it

Each catalog member is compared with the fit under both symmetry groups:

$$d\!\left(\mathbf{R}, \mathbf{R}_{\text{cand}}\right) =
\min_{\mathbf{S}_c,\, \mathbf{S}_p} \angle\!\left(
\mathbf{S}_c\,\mathbf{R}\,\mathbf{S}_p,\; \mathbf{R}_{\text{cand}} \right).$$

The catalog itself is chosen from the two crystal systems through one dispatch
table (`default_relationship_catalog`), and the separations it must resolve are
fixed crystallography, not tuning:

| | Bain | KS | NW | GT | Pitsch |
| --- | --- | --- | --- | --- | --- |
| **Bain** | 0.000 | 11.065 | 9.736 | 10.146 | 9.736 |
| **Kurdjumov-Sachs** | 11.065 | 0.000 | 5.264 | **2.404** | 5.264 |
| **Nishiyama-Wassermann** | 9.736 | 5.264 | 0.000 | 2.861 | 7.444 |
| **Greninger-Troiano** | 10.146 | 2.404 | 2.861 | 0.000 | 5.787 |
| **Pitsch** | 9.736 | 5.264 | 7.444 | 5.787 | 0.000 |

*Symmetry-reduced separations in degrees, fcc → bcc.* The smallest gap, 2.404°
between Kurdjumov-Sachs and Greninger-Troiano, is what sets the usable noise
budget. For the hexagonal product the cubic→hexagonal catalog holds Burgers and
Shoji-Nishiyama, separated by 42.848°, so that identification is far easier.

:::{admonition} Constraint: only the candidates offered can win
:class: warning

Cubic-to-cubic resolves to the fcc→bcc family because point-group symmetry cannot
distinguish an fcc phase from a bcc one. Supply an explicit catalog when that
assumption is wrong. A rotation matching nothing is reported as matching
nothing — `matches_catalog` is `False` and `describe()` says so.
:::

## 6. Stating it crystallographically

A rotation matrix is unreadable; "$(111)_{\gamma}$ parallel to
$(011)_{\alpha}$" is the working fact. `describe_orientation_relationship`
recovers that statement from the rotation alone.

Every canonical-sign primitive parent triple up to `max_index` is carried into
the child basis at once — planes on the reciprocal basis, directions on the
direct basis, so the zone law is preserved — and compared with every candidate
child triple in a single cosine matrix. A clause is kept when

$$\left| \cos \angle\!\left( \hat{\mathbf{g}}_{\text{image}},\,
\hat{\mathbf{g}}_{\text{child}} \right) \right| \ge
\cos(\texttt{tolerance\_deg}).$$

Comparison uses $\lvert\cos\rvert$ because the canonical-sign filter has already
collapsed each antiparallel pair to one representative.

A rotation has three degrees of freedom, so **one plane clause fixes two and one
in-plane direction clause fixes the third** — which is exactly the classical form
of an orientation relationship. Everything else follows.

| system | relationship | recovered statement |
| --- | --- | --- |
| cubic | Kurdjumov-Sachs | $(111) \parallel (011)$, $[10\bar{1}] \parallel [11\bar{1}]$ |
| cubic | Nishiyama-Wassermann | $(111) \parallel (011)$, $[1\bar{1}0] \parallel [100]$ |
| hexagonal | Burgers | $(011)_{\beta} \parallel (0001)_{\alpha}$, $[\bar{1}11]_{\beta} \parallel [\bar{1}2\bar{1}0]_{\alpha}$ |

all at zero deviation, because these are the defining parallelisms. Hexagonal
phases are labelled in four-index Miller-Bravais form throughout; a three-index
hexagonal label hides the symmetry of the family and is not how the hcp
literature states a relationship.

:::{admonition} Constraint: several statements are true at once
:class: important

A rotation typically satisfies *several* exact low-index parallelisms
simultaneously. For Kurdjumov-Sachs both $(111) \parallel (011)$ and
$(10\bar{1}) \parallel (11\bar{1})$ are exact, and index magnitude alone
tie-breaks between them arbitrarily.

Which one the literature quotes is a fact about the two **structures** — their
close-packed planes and directions — not about the rotation, which does not know
them. The search therefore takes a preference: by default the relationship's own
recorded defining families, and for a fitted relationship those of the catalog
member it matched. Fit quality outranks preference in the sort, so a nominated
family can never promote a visibly worse clause above an exact one.
:::

## 7. The verdict, and how to read it

`is_conclusive` requires the winner both to fit within `catalog_tolerance_deg`
**and** to lead the runner-up by more than the measurement scatter and its own
misfit:

$$\text{margin} > \max\!\left( \overline{\text{residual}},\;
d_{\text{best}} \right).$$

Those are the two quantities that could otherwise explain the lead away. The
report states the verdict in words, not only as a flag.

**Failure modes, all reported rather than hidden.**

- *Scatter comparable to the catalog spacing.* Verdict becomes inconclusive.
  Measured above: fine to 2°, inconclusive at 5°.
- *A single pair.* Scatter is zero by construction and says nothing about
  measurement quality; `describe()` states this explicitly.
- *No catalog for the phase pair.* The fit is reported without a name rather
  than forced onto an inapplicable list.
- *No low-index parallelism within tolerance.* Reported as such — informative in
  itself, since it means the relationship is not of the classical
  parallel-planes type at that index bound.

## 8. Limits

- Validation is **synthetic**: planted variants of a known relationship,
  recovered with the relationship withheld. Measured-EBSD fixtures remain
  outstanding, and no claim of MTEX parity is made.
- The parallelism search is bounded by `max_index` (default 3 on both sides). A
  relationship defined by higher-index parallelisms reports none rather than
  inventing one.
- Pairs must be row-matched and share a specimen frame. Grain-mean orientations
  are assumed; the method does not itself segment grains.

## Verification

| Claim | Where it is checked |
| --- | --- |
| Defining parallelisms recovered at zero deviation (KS, NW, Burgers) | `tests/unit/test_or_characterization.py` |
| Published catalog separations reproduced (KS-NW 5.26°, KS-GT 2.40°) | same |
| All five fcc→bcc members and Burgers recovered from their own variants | same |
| Bain survives the double-coset tie | `test_bain_survives_the_double_coset_tie` |
| A random rotation matches no catalog member | same file |
| Live numerical demonstration | [worked examples](../examples/generated/transformation.md) |

## See also

- {doc}`../concepts/orientation_relationships` — the concept-level tour
- {doc}`variant_correspondence` — what the relationship does to an arbitrary
  plane or direction
- {doc}`../tutorials/notebooks/23_transformation_crystallography_end_to_end` —
  this algorithm inside the end-to-end Burgers workflow

## References

### Normative

- {doc}`../architecture/orientation_relationship_analysis_foundation`
- {doc}`../standards/notation_and_conventions`

### Informative

- Markley, Cheng, Crassidis and Oshman, *Journal of Guidance, Control and
  Dynamics* 30 (2007) 1193 — quaternion averaging.
- Kurdjumov and Sachs, *Zeitschrift für Physik* 64 (1930) 325.
- Burgers, *Physica* 1 (1934) 561.
- Morito, Tanaka, Konishi, Furuhara and Maki, *Acta Materialia* 51 (2003) 1789 —
  variant numbering and packet structure.
