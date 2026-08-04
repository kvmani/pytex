# Variant-Resolved Plane And Direction Correspondence

**Surface:** `pytex.core.transformation.variant_correspondence_table`, built on
`map_plane_to_child` / `map_direction_to_child` and their parent-inverses.

Given an orientation relationship and any parent plane $(hkl)$ or direction
$[uvw]$, this tabulates what that object becomes in **every** product variant,
grouped so the answer is readable rather than 24 rows of indices.

```{figure} ../../figures/variant_correspondence_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 takes the nominated object and the variant
  set. Lane 2 maps each object through each variant on the correct basis,
  rationalizes the irrational image to the nearest primitive integer triple by
  true angle, and groups variants whose images are symmetry-equivalent. Lane 3
  emits the table, the exactly-parallel subset, and the exports.
:width: 100%

The algorithm, with the constraint governing each stage.
```

## 1. The two index maps, and why they differ

A relationship is a rotation $\mathbf{R}$ on Cartesian vectors, but crystal
objects are given as *indices* — components in a lattice basis. Two different
bases are involved, and using the wrong one is the classic error.

Let $\mathbf{A}$ be the direct structure matrix whose columns are the lattice
vectors in Cartesian coordinates, and $\mathbf{A}^{*}$ the reciprocal one. A
**direction** $[uvw]$ has Cartesian image $\mathbf{A}\mathbf{u}$; a **plane**
$(hkl)$ has Cartesian *normal* $\mathbf{A}^{*}\mathbf{h}$, because Miller indices
are already reciprocal-basis components. So the two index maps are

$$\mathbf{M} = \mathbf{A}_c^{-1}\,\mathbf{R}\,\mathbf{A}_p , \qquad
\mathbf{M}^{*} = \left(\mathbf{A}_c^{*}\right)^{-1}\mathbf{R}\,\mathbf{A}_p^{*} ,$$

and they satisfy $\mathbf{M}^{*} = \mathbf{M}^{-\mathsf{T}}$, which is exactly
what preserves the zone law: if $\mathbf{h}\cdot\mathbf{u} = 0$ in the parent
then $\mathbf{h}_c\cdot\mathbf{u}_c = 0$ in the child. A plane containing a
direction still contains its image.

For a cubic-to-cubic relationship $\mathbf{A}^{*} \propto \mathbf{A}^{-\mathsf{T}}$
with a scalar factor, and the two maps coincide up to that scale — which is why
the error is invisible in cubic tests and appears the moment a hexagonal phase is
involved. The library routes every mapping through one pair of helpers so the
choice cannot be made per call site.

## 2. Variants

A variant is the relationship composed with a parent symmetry operation,
$\mathbf{V}_k = \mathbf{R}\,\mathbf{S}_{p,k}$, reduced by the child-symmetry
orbit so that two operations giving crystallographically identical children count
once. That reduction is what produces the literature counts:

| relationship | system | variants |
| --- | --- | --- |
| Bain | cubic → cubic | 3 |
| Nishiyama-Wassermann, Pitsch | cubic → cubic | 12 |
| Kurdjumov-Sachs, Greninger-Troiano | cubic → cubic | 24 |
| **Burgers** | cubic → hexagonal | **12** |
| Shoji-Nishiyama | cubic → hexagonal | 4 |

Burgers has 12 because six $\{110\}_{\beta}$ planes can become the basal plane
and each admits two $\langle 111 \rangle_{\beta}$ directions as
$\langle 11\bar{2}0 \rangle_{\alpha}$.

## 3. Rationalizing the image

The image of a low-index parent object is, in general, **irrational** — its
components in the child basis are not integers. The exact components are kept,
and a nearest integer triple is reported alongside with the angle between them:

$$\mathbf{t}^{\star} = \arg\max_{\mathbf{t}}
\frac{\left( \mathbf{B}\mathbf{t} \right) \cdot \mathbf{g}_{\text{exact}}}
{\lVert \mathbf{B}\mathbf{t} \rVert\; \lVert \mathbf{g}_{\text{exact}} \rVert},
\qquad
\rho = \operatorname{atan2}\!\left(
\lVert \hat{\mathbf{t}} \times \hat{\mathbf{g}} \rVert,\;
\hat{\mathbf{t}} \cdot \hat{\mathbf{g}} \right),$$

over primitive integer triples with entries bounded by `max_index`, with
$\mathbf{B}$ the reciprocal basis for planes and the direct basis for directions.

Three details matter:

- Candidates are compared through their **Cartesian images**, so the residual
  $\rho$ is a true angle in the relevant space and is meaningful for any lattice
  metric — not an index-space distance, which would be metric-blind.
- `atan2` rather than `arccos`, because `arccos` loses precision exactly where
  these residuals live, near zero.
- The match is **sign-sensitive**: the triple whose image points *along* the
  exact image wins. This matters because a plane's normal direction is a
  diffraction vector $\mathbf{g}$, and $\mathbf{g}$ and $-\mathbf{g}$ are
  different reflections. A consequence is that the hexagonal basal image may be
  reported as $(0001)$ or $(000\bar{1})$ depending on the variant; both name the
  same plane.

The default bound is 17, chosen to cover the full standard catalog including the
Greninger-Troiano $\langle 5\,12\,17 \rangle$ direction family.

:::{admonition} Constraint: what the index bound does and does not change
:class: warning

Raising `max_index` **never worsens** a residual — the larger candidate set
contains the smaller — and **never changes which correspondences are exact**. The
set of exactly-parallel variants is identical at bound 3 and bound 17.

What it does change is how the *irrational* images are labelled, and therefore
how many index families they fall into. "Four distinct images" is a statement
partly about the bookkeeping, not purely about the crystallography, and
`describe()` says so rather than letting the count be over-read.
:::

## 4. Grouping: what makes 24 rows readable

Each rationalized image is reduced to its symmetry-canonical family
representative under the *image phase* point group; variants sharing a
representative share an `equivalence_group`. This turns a wall of indices into
the actual answer.

### Cubic: the packet structure of lath martensite

Kurdjumov-Sachs, the austenite $(111)$ across all 24 variants:

| group | variants | image | residual |
| --- | --- | --- | --- |
| 0 | 6 | $(011)$ | **0.0000°** |
| 1 | 6 | $(5\bar{1}4)$ | 0.3646° |
| 2 | 6 | $(\overline{17}\,\bar{1}\,10)$ | 0.6868° |
| 3 | 6 | $(1\,11\,\bar{5})$ | 1.0313° |

Four distinct answers, six variants each. The parent symmetry acts transitively,
so the split is even. The six exact ones are the variants whose close-packed
plane *is* this $(111)$ — one packet in Morito's sense — and the test asserts
they are the same six that `variant_close_packed_groups` returns, comparing two
independent computations rather than a stored constant.

A nominated *direction* is more selective: $[1\bar{1}0]$ lies in two of the four
$\{111\}$ planes and each contributes two variants, so it is exactly parallel in
**4 of 24**.

### Hexagonal: six packets of two

Burgers, the $\beta$ $(011)$ across all 12 variants:

| group | variants | image | residual |
| --- | --- | --- | --- |
| 2 | 2 | $(0001)$ | **0.0000°** |
| 3 | 2 | $(\overline{10}\,1\,9\,0)$ | 0.0557° |
| 1 | 4 | $(\overline{16}\,16\,0\,17)$ | 0.0917° |
| 0 | 4 | $(3\,13\,\overline{16}\,16)$ | 0.7171° |

Exactly two variants make this $\{110\}$ the basal plane — six packets of two,
the hexagonal counterpart of the cubic four-by-six. The groups are *not* evenly
sized here, because the hexagonal child group is smaller (12 proper operators
against the cubic 24) and reduces fewer images together.

### The reverse map is not selective

Mapping the child $(0001)$ back with `sense="child_to_parent"` lands on a
$\{110\}_{\beta}$ plane in **all 12 variants**, one equivalence group. Every
variant's basal plane came from *some* $\{110\}$, so the reverse question has one
answer where the forward question has several. Worth knowing before reading a
forward table as if it were symmetric.

## 5. Cost

Each row is one call to the shared index-map helper, so the table is
$\mathcal{O}(n_{\text{objects}} \times n_{\text{variants}})$ rationalizations,
each a vectorized comparison against the primitive-triple candidate set. For 24
variants and the default bound this is milliseconds; the surface is not on any
hot path, and correctness through the shared helper was preferred to a bespoke
vectorization.

## 6. Limits

- Objects must be homogeneous in kind and phase; call once for planes and once
  for directions. Mixing them is rejected rather than guessed.
- Direction indices must be integer-valued, because the family enumeration is an
  integer orbit.
- The table describes a *geometric* correspondence. It says nothing about which
  variants actually form — that is variant selection, driven by stress and
  interface energy, and is a separate question.

## Verification

| Claim | Where it is checked |
| --- | --- |
| KS $(111)$ gives 4 groups of 6, with 6 exact $\{011\}$ | `tests/unit/test_variant_correspondence_table.py` |
| Those 6 are exactly the packet `variant_close_packed_groups` returns | same |
| $[1\bar{1}0]$ exact in 4 of 24 | same |
| Burgers $(011)_{\beta}$ basal in 2 of 12 | same |
| Reverse map exact in all 12, one group | same |
| Raising `max_index` never worsens a residual; exact set unchanged | same |
| Forward-then-reverse round trip returns the source | same |
| Live numerical demonstration | [worked examples](../examples/generated/transformation.md) |

## See also

- {doc}`orientation_relationship_determination` — where the relationship itself
  comes from
- {doc}`composite_saed_assembly` — the same variant machinery producing a
  diffraction observable
- {doc}`../concepts/orientation_relationships`

## References

### Normative

- {doc}`../architecture/orientation_relationship_analysis_foundation`
- {doc}`../standards/hexagonal_and_trigonal_conventions`

### Informative

- Morito, Tanaka, Konishi, Furuhara and Maki, *Acta Materialia* 51 (2003) 1789 —
  packet structure.
- Burgers, *Physica* 1 (1934) 561.
- International Tables for Crystallography, Vol. A — reciprocal bases and the
  zone law.
