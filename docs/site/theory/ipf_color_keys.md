# Inverse-Pole-Figure Colour Keys

An inverse-pole-figure map is the most reproduced image in texture analysis, and the least
specified. Every textbook shows the coloured standard triangle; almost none states what function
carries a direction to a colour. That function is not physics. It is a chain of four conventions,
and each one changes every pixel of the map.

This note fixes that chain for `pytex.plotting.ipf.IPFColorKey`: what the colour *is*, the
closed form it reduces to for cubic symmetry, and the consequences of the choices — including two
that constrain what an IPF map is able to show at all.

## What The Colouring Has To Achieve

Colouring a direction $\hat{\mathbf{d}}$ in the crystal frame is not free. Four requirements fix
almost everything:

1. **Symmetry invariance.** Two directions related by a crystal symmetry operation are the same
   physical direction and must receive the same colour. Otherwise the map would show contrast
   where there is no crystallography, and an arbitrary choice of equivalent index would change
   the picture.
2. **Corner anchoring.** The corners of the fundamental sector — for cubic symmetry $[001]$,
   $[101]$, $[111]$ — must take the three primaries, so a reader can name a colour without
   consulting the legend.
3. **Continuity.** Nearby directions must take nearby colours, or grain interiors would speckle.
4. **Determinism.** The same direction must give the same colour in every run, with no dependence
   on which symmetry operator happened to be enumerated first.

Requirements 1 and 2 are met by construction below. Requirement 3 holds everywhere except on the
sector boundary, where it holds only because the boundary directions are themselves
symmetry-equivalent across the edge.

## Stage 1: Fold Into The Fundamental Sector

The symmetry group $G$ partitions the sphere into orbits. Reduction picks the one representative
of each orbit lying in the fundamental sector $\mathcal{T}$ — the standard stereographic triangle,
generalised to any point group:

$$
\hat{\mathbf{d}}_{\mathcal{T}} = \operatorname{reduce}_{G}(\hat{\mathbf{d}}),
\qquad
\hat{\mathbf{d}}_{\mathcal{T}} \in \mathcal{T},
\qquad
\hat{\mathbf{d}}_{\mathcal{T}} = \mathbf{S}_{i}\hat{\mathbf{d}}
\ \ \text{for some } \mathbf{S}_{i} \in G .
$$ (eq-ipf-reduce)

This single step delivers requirement 1: everything downstream is a function of
$\hat{\mathbf{d}}_{\mathcal{T}}$ alone, so symmetric equivalents cannot disagree. In PyTex it is
`SymmetrySpec.reduce_vectors_to_fundamental_sector`, and `antipodal=True` additionally identifies
$\pm\hat{\mathbf{d}}$, which is correct for a direction whose sense is not observable.

**Closed form for $m\bar{3}m$.** For the cubic Laue class with antipodal identification, the
24 proper rotations and the inversion act on a direction by permuting and signing its components.
The orbit therefore contains every signed permutation, and the sector representative is obtained
by taking absolute values and sorting:

$$
\hat{\mathbf{d}}_{\mathcal{T}}
= \bigl(\,|d|_{(2)},\ |d|_{(1)},\ |d|_{(3)}\,\bigr),
\qquad
|d|_{(1)} \le |d|_{(2)} \le |d|_{(3)},
$$ (eq-ipf-cubic-reduce)

so the largest component becomes $z$, the middle becomes $x$, and the smallest becomes $y$. The
standard triangle is then exactly the ordered cone

$$
\mathcal{T}_{m\bar{3}m} = \{\, \hat{\mathbf{d}} : d_{z} \ge d_{x} \ge d_{y} \ge 0 \,\},
$$ (eq-ipf-triangle)

whose three bounding planes are precisely the sector edge normals PyTex stores:
$d_{y} \ge 0$, $d_{z} \ge d_{x}$, and $d_{x} \ge d_{y}$. No general group-theoretic machinery is
needed for the cubic case — a sort suffices — and the vectorized reduction agrees with
{eq}`eq-ipf-cubic-reduce` to $1.5\times10^{-15}$ over 500 random directions.

## Stage 2: Barycentric Position In The Sector

Colour is assigned by *where in the triangle* the direction falls. Let $\mathbf{K}$ be the matrix
whose columns are the three sector corner directions. The barycentric weights $\boldsymbol{\beta}$ solve

$$
\mathbf{K}\,\boldsymbol{\beta} = \hat{\mathbf{d}}_{\mathcal{T}},
\qquad
\mathbf{K} = \bigl[\, \hat{\mathbf{s}}_{1} \ \ \hat{\mathbf{s}}_{2} \ \ \hat{\mathbf{s}}_{3} \,\bigr].
$$ (eq-ipf-bary)

This is a $3\times3$ linear solve per direction, batched by `numpy.linalg.solve`. Because the three
corners of a proper sector are linearly independent, $\mathbf{K}$ is invertible and the weights are
unique.

**Closed form for $m\bar{3}m$.** The cubic corners are
$\hat{\mathbf{s}}_{1} = [001]$, $\hat{\mathbf{s}}_{2} = [101]/\sqrt{2}$,
$\hat{\mathbf{s}}_{3} = [111]/\sqrt{3}$, so $\mathbf{K}$ is triangular enough to invert by
inspection. Reading {eq}`eq-ipf-bary` row by row from the bottom gives

$$
\beta_{1} = d_{z} - d_{x},
\qquad
\beta_{2} = \sqrt{2}\,\left(d_{x} - d_{y}\right),
\qquad
\beta_{3} = \sqrt{3}\,d_{y},
$$ (eq-ipf-cubic-bary)

with all components those of $\hat{\mathbf{d}}_{\mathcal{T}}$. Two things fall out immediately.
First, the weights are non-negative exactly on {eq}`eq-ipf-triangle`, so the sector membership test
and the colour computation are the same three inequalities. Second, the barycentric solve is
avoidable entirely for cubic symmetry. PyTex keeps the general solve because it must serve all
32 point groups, but {eq}`eq-ipf-cubic-bary` is what it computes there, and the two agree to
$2.2\times10^{-16}$.

## Stage 3: Clip, Then Normalise

Negative weights are clipped to zero and the result is scaled to sum to one:

$$
\tilde{\beta}_{i} = \frac{\max(\beta_{i}, 0)}{\sum_{j}\max(\beta_{j}, 0)} .
$$ (eq-ipf-clip)

The clip is not cosmetic. A direction reduced onto a sector edge can land a hair outside it in
floating point, giving a small negative weight and, without the clip, a colour channel outside
$[0,1]$. Clipping projects such a direction onto the boundary, which is where it belongs. The
denominator vanishing is the genuine degeneracy — a direction with no component inside the sector
cone — and is raised rather than silently divided.

## Stage 4: Weights To Colour

The weights are mapped through the corner colours $\mathbf{C}_{\mathrm{rgb}}$, whose rows are the
RGB triples assigned to the three corners:

$$
\mathbf{c}_{0} = \tilde{\boldsymbol{\beta}}^{\mathsf{T}} \mathbf{C}_{\mathrm{rgb}} .
$$ (eq-ipf-mix)

With the default $\mathbf{C}_{\mathrm{rgb}} = \mathbf{I}_{3}$ — red, green, blue on corners 1, 2, 3
— this is the identity, and **the RGB triple simply is the barycentric coordinate**. That is the
whole of the "standard" IPF colouring, and it is worth saying plainly because it is usually left
implicit: the colour of a direction is its position in the triangle, expressed in a basis of
corners.

## Stage 5 And 6: Saturation, And Why Every Colour Is Fully Saturated

Two operations remain, and neither is derivable from anything:

$$
c_{1,i} = c_{0,i}^{\,1/\gamma_{s}},
\qquad
\mathbf{c}_{\mathrm{rgb}} = \frac{\mathbf{c}_{1}}{\max_{j} c_{1,j}} .
$$ (eq-ipf-gamma)

The exponent $1/\gamma_{s}$ with the default $\gamma_{s} = 0.5$ **squares** each channel. Squaring
values in $[0,1]$ suppresses the small ones relative to the large, which pushes each colour toward
its dominant primary and widens the visual separation between orientations near a corner. It is a
contrast control, tuned by eye, and nothing in crystallography selects the value.

The second operation divides by the largest channel, so the brightest channel is always exactly 1.

**These two steps collapse.** Because $x \mapsto x^{1/\gamma_{s}}$ is monotonic on $[0,\infty)$, it
commutes with taking the maximum, and {eq}`eq-ipf-gamma` composed with {eq}`eq-ipf-clip` and
{eq}`eq-ipf-mix` reduces — for the default identity colour basis — to a single expression in the
*unnormalised* weights:

$$
c_{\mathrm{rgb},i}
= \left( \frac{\max(\beta_{i},0)}{\max_{j}\max(\beta_{j},0)} \right)^{1/\gamma_{s}} .
$$ (eq-ipf-closed)

Verified against the implementation to $1.0\times10^{-15}$ over 2000 random directions. Three
consequences follow, and each is a real constraint on what an IPF map can express:

- **The sum normalisation in {eq}`eq-ipf-clip` is inert.** Only ratios of weights survive into
  {eq}`eq-ipf-closed`, and the division by the sum cancels against the division by the maximum.
  It runs for its degeneracy check, not for its arithmetic. This holds for any
  $\mathbf{C}_{\mathrm{rgb}}$, since a linear map commutes with the scaling too.
- **Every IPF colour lies on the saturated surface of the RGB cube**, with
  $\max_{i} c_{\mathrm{rgb},i} = 1$ exactly (verified: the minimum over 2000 random directions of
  the maximum channel is $1.0$). The colouring uses a two-dimensional surface of the
  three-dimensional cube, which is correct — the sector is two-dimensional — but it means
  **brightness carries no information**. Any attempt to modulate an IPF map by a second scalar,
  such as band contrast or confidence index, necessarily breaks the key: the colour is no longer
  invertible to a direction.
- **Colour is invariant to the length of $\hat{\mathbf{d}}$**, so normalisation of the input
  matters only through the reduction.

## Agreement With Known Answers

Every value below is an analytic identity derived from {eq}`eq-ipf-cubic-bary` and
{eq}`eq-ipf-closed`, not a recorded program output.

| Direction | $\boldsymbol{\beta}$ (exact) | Colour (exact) | Computed error |
| --- | --- | --- | ---: |
| $[001]$ | $(1, 0, 0)$ | $(1, 0, 0)$ red | $0$ |
| $[101]$ | $(0, 1, 0)$ | $(0, 1, 0)$ green | $0$ |
| $[111]$ | $(0, 0, 1)$ | $(0, 0, 1)$ blue | $0$ |
| $[112]$ | $\tfrac{1}{\sqrt{6}}(1, 0, \sqrt{3})$ | $(\tfrac{1}{3}, 0, 1)$ | $1.1\times10^{-15}$ |
| $[113]$ | $\tfrac{1}{\sqrt{11}}(2, 0, \sqrt{3})$ | $(1, 0, \tfrac{3}{4})$ | $9.5\times10^{-30}$ |

The $[112]$ and $[113]$ rows are worth following by hand, because they exercise the whole chain.
For $[113]$, normalising gives $\hat{\mathbf{d}} = (1,1,3)/\sqrt{11}$, already ordered as
$d_z \ge d_x \ge d_y$. Then {eq}`eq-ipf-cubic-bary` gives
$\beta_{1} = 2/\sqrt{11}$, $\beta_{2} = 0$, and $\beta_{3} = \sqrt{3}/\sqrt{11}$. The largest is $\beta_{1}$, so
{eq}`eq-ipf-closed` gives $\bigl(1,\ 0,\ (\sqrt{3}/2)\bigr)^{2} = (1, 0, 3/4)$ — a red with
three-quarters blue and no green, which is the magenta-leaning red seen along the $[001]$–$[111]$
edge. $[112]$ lies further along the same edge, and its $\beta_{3}$ now dominates:
$(1/\sqrt{3}, 0, 1)^{2} = (1/3, 0, 1)$, a blue with one-third red.

Symmetry invariance is exact rather than approximate: colouring all 24 equivalents of a general
direction gives colours agreeing to $7.1\times10^{-15}$, which is the accumulated error of the
rotation arithmetic and not a tolerance in the colouring.

## What This Colouring Is Not

- **Not perceptually uniform.** Equal angular steps on the sphere do not give equal perceptual
  colour steps. A gradient that looks abrupt need not correspond to a large misorientation, and
  the eye's poor discrimination in the green–cyan region compresses part of the triangle. Read
  quantitative misorientation from a KAM or GROD map, never from IPF colour.
- **Not comparable across symmetries.** The sector, and therefore the corners and the whole
  mapping, is a property of the point group. Two phases of different symmetry coloured side by
  side share no colour meaning, which is why `IPFColorKey` binds a `SymmetrySpec` and refuses
  orientations carrying a different one.
- **Not a single "the" IPF colouring.** $\gamma_{s}$, the corner-to-primary assignment, and the
  antipodal convention are all free. Maps from two tools agree on the corners and disagree in
  between. Comparisons should quote the key, not just the picture.
- **Not defined by a triangle for every group.** Low-symmetry groups have sectors that are
  hemispheres or wedges with fewer than three corners. PyTex then anchors the colour basis on the
  reference octant $[001], [100], [010]$ rather than on sector corners. The result is continuous
  and symmetry-invariant, but its corners are not crystallographic features of the sector, and it
  should be treated as a display convention rather than a standard.

## Assumptions And Limits

- The specimen direction is fixed per map. An "IPF map" is always *IPF-$\hat{\mathbf{z}}$* or
  IPF-RD and so on; the colour answers "which crystal direction points along this specimen axis",
  and the axis must be stated for the map to mean anything.
- Colour is assigned per orientation, with no spatial smoothing.
- The legend mesh is sampled on a polar/azimuth raster and masked to the sector, so its outline is
  resolution-limited even though the colours themselves are exact.

## References

### Normative

- Th. Hahn (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*,
  IUCr / Springer. DOI: <https://doi.org/10.1107/97809553602060000100>. Point groups and the
  asymmetric unit that defines the fundamental sector.

### Informative

- G. Nolze and R. Hielscher, *Orientations — perfectly colored*, Journal of Applied
  Crystallography **49** (2016) 1786–1802. DOI:
  <https://doi.org/10.1107/S1600576716012942>. The reference treatment of orientation colouring,
  including why a fully general colouring cannot be simultaneously continuous, symmetry-invariant,
  and perceptually uniform.
- T. B. Britton et al., *Tutorial: Crystal orientations and EBSD — or which way is up?*,
  Materials Characterization **117** (2016) 113–126. DOI:
  <https://doi.org/10.1016/j.matchar.2016.04.008>.

## See Also

- {doc}`fundamental_region_reduction` — the reduction {eq}`eq-ipf-reduce` relies on.
- {doc}`discrete_odf_and_pole_figures` — the pole-figure side of the same projection machinery.
- {doc}`/workflows/ipf_colors` — the runtime surface and worked usage.
