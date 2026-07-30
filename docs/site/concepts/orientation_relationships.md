# Orientation Relationships And Index Correspondence

An orientation relationship (OR) fixes how a product ("child") crystal is
oriented relative to its parent during a phase transformation. PyTex treats OR
analysis as a flagship capability: the goal is that indices, orientations,
variants, and (eventually) map-scale reconstructions all move through one
explicit semantic model. This page covers the stable surface implemented
today.

## The three objects an OR carries

PyTex keeps three mathematically distinct objects explicit (see the
{doc}`OR analysis foundation <../architecture/orientation_relationship_analysis_foundation>`
for the full doctrine):

1. **The rigid rotation** $\mathbf{R}$ between parent and child Cartesian
   crystal frames (`parent_to_child_rotation`). It maps unit vectors and
   composes with orientations.
2. **The index correspondence** — the linear maps that carry Miller indices:
   $\mathbf{u}_{c} = \mathbf{M}\,\mathbf{u}_{p}$ for directions with
   $\mathbf{M} = \mathbf{A}_{c}^{-1}\mathbf{R}\,\mathbf{A}_{p}$ (direct
   structure matrices), and $\mathbf{h}_{c} = \mathbf{M}^{*}\,\mathbf{h}_{p}$
   for planes with $\mathbf{M}^{*} = \mathbf{M}^{-\mathsf{T}}$ (reciprocal
   bases). The inverse-transpose relation preserves the zone law
   $\mathbf{h} \cdot \mathbf{u}$, so a direction lying in a plane stays in
   the mapped plane. Correspondence matrices are generally not rotations and
   generally irrational.
3. **The transformation deformation** — `deformation_gradient()` builds the
   nearest-integer lattice correspondence and returns the parent-frame
   gradient with its polar decomposition: for fcc-bcc steel parameters the
   textbook Bain principal strains, and for KS/NW the literature rigid-body
   rotations relative to Bain ($11.06^{\circ}$/$9.74^{\circ}$) as the
   residual polar rotation.

## Composition convention

PyTex orientations map **crystal to specimen** (the normative convention in
the notation standard). The child orientation produced by a variant is
therefore

$$
g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}},
$$

so that corresponding parent and child crystal directions
($\mathbf{d}_{c} = \mathbf{V}\,\mathbf{d}_{p}$) point along the same
specimen direction — the physical meaning of an OR parallelism. Every
prediction, deviation, fitting, and reconstruction surface uses this
composition, crystal-symmetry equivalents act by right multiplication
($g' = g\,S$), and a regression test pins the convention against the
specimen-space parallelism identity.

## Constructing named relationships

`OrientationRelationship` ships literature constructors — Bain,
Kurdjumov-Sachs, Nishiyama-Wassermann, Greninger-Troiano, Pitsch (fcc-bcc),
Burgers (bcc-hcp), and Shoji-Nishiyama (fcc-hcp, the epsilon-martensite
relationship with 4 variants) — plus the generic
`from_parallel_plane_direction(...)` builder for user-defined ORs. Catalogs
(`standard_fcc_bcc_relationships`, `standard_bcc_hcp_relationships`,
`standard_fcc_hcp_relationships`) bundle them per phase pair.

## Mapping indices across the relationship

The correspondence surface answers the canonical OR questions:

- `map_plane_to_child(plane)` / `map_direction_to_child(direction)` — which
  child $(hkl)$ / $[uvw]$ corresponds to this parent plane or direction?
- `map_plane_to_parent(...)` / `map_direction_to_parent(...)` — the inverse,
  for reading product-phase measurements against parent stereography.
- `correspondence_direct()` / `correspondence_reciprocal()` — the raw
  $\mathbf{M}$ and $\mathbf{M}^{*}$ matrices, per relationship or per
  variant.
- `map_plane_across_variants(...)` / `map_direction_across_variants(...)` —
  the variant-resolved tables.

Every mapping returns the **exact** (generally irrational) target components,
the nearest primitive-integer **rationalization** (bounded index search, sign
sensitive, default bound 17 to cover the Greninger-Troiano
$\langle 5\,12\,17 \rangle$ family), and the **angular residual** between
them, so "nearly parallel" is always a quantified statement rather than an
implicit rounding.

Hexagonal phases participate with full index meaning: Burgers maps
$(110)_{\beta}$ to the basal plane $(0001)_{\alpha}$ and
$[\bar{1}11]_{\beta}$ to $[11\bar{2}0]_{\alpha}$ (stored as the
three-index $[110]$).

## Variant pole figures

`variant_pole_figure(parent_orientation, relationship, child_plane)` predicts
where every variant's child plane family lands on the specimen sphere, and
`plot_variant_pole_figure(...)` renders the color-per-variant stereographic
overlay — the standard way to read measured product-phase pole figures for
operative variants and variant selection. The prediction is pinned by the
packet-plane coincidence: each KS variant's $\{011\}$ pole set contains the
specimen-frame normal of its packet's parent $\{111\}$ member.

## Variant packets

`variant_close_packed_groups(relationship, parent_plane)` labels each variant
by the parent family member it carries into exact parallelism — the packet
classification of martensite crystallography. Kurdjumov-Sachs with the
$\{111\}$ family yields the four packets of six variants of lath martensite
(Morito et al.); Burgers with $\{110\}$ yields six groups of two.

## The OR as a misorientation, and deviation from it

`misorientation()` returns the symmetry-reduced (disorientation)
representative of the relationship — the way ORs are measured and reported
from EBSD boundary data. For Kurdjumov-Sachs this is the published
$42.85^{\circ}$ rotation about $\langle 0.968\;0.178\;0.178 \rangle$;
Nishiyama-Wassermann gives $45.99^{\circ}$ and Bain $45^{\circ}$ about
$\langle 100 \rangle$.

`or_deviation(parents, children, relationship)` quantifies how well measured
parent/child orientation pairs obey a nominal OR: for each pair it reports
the smallest child-symmetry-reduced angle to any variant prediction, plus the
winning variant index. Exact synthetic data returns zeros; children generated
with Greninger-Troiano deviate by the documented $2.40^{\circ}$ from
Kurdjumov-Sachs and $2.86^{\circ}$ from Nishiyama-Wassermann — the report's
aggregate statistics are the entry point for OR fitting.

## Which boundaries can share a parent

Reconstructing parent grains, or identifying the operative relationship from a
fully transformed microstructure, both reduce to one question: *could these two
neighbouring product grains have come from the same parent?*

The answer is exact. A child formed from parent $\mathbf{P}$ through variant
$\mathbf{V}_i$ is $\mathbf{C}_i = \mathbf{P}\mathbf{V}_i^{\mathsf{T}}$, so the
crystal-frame boundary misorientation of two same-parent children is

$$\mathbf{C}_i^{\mathsf{T}}\mathbf{C}_j = \mathbf{V}_i \mathbf{V}_j^{\mathsf{T}}.$$

Writing $\mathbf{V}_i = R\,S_{p,i}$ with $R$ the parent-to-child rotation and
$S_p$ the parent point group, the whole family collapses to $R\,G_p\,R^{\mathsf{T}}$
— the parent group *conjugated* by the OR rotation. Each child orientation is
itself defined only up to its own crystal symmetry, so the observable set is the
double coset

$$G_c \left(R\,G_p\,R^{\mathsf{T}}\right) G_c,$$

returned by `intervariant_boundary_fingerprint(relationship)`.
`boundary_fingerprint_distances_deg(...)` scores measured boundaries against it,
and a distance near zero means the two grains are consistent with a shared
parent. Both are used by `reconstruct_parent_grains` and
`identify_orientation_relationship`, so the two workflows share one definition.

**Match the rotation, not just the angle.** It is tempting to compare only the
misorientation angle against the intervariant table from
`intervariant_misorientation_angles_deg`. That discards the axis, and for a
cubic-cubic relationship the angle spectrum is dense enough over its range that
a few-degree window admits most *unrelated* boundaries: measured on uniformly
random rotations, an angle-only test at $3^{\circ}$ accepts about $53\%$ of
them for Kurdjumov-Sachs, against about $7\%$ for the full-rotation test — and
at $1^{\circ}$, $29\%$ against $0.3\%$. At map scale that difference is the
difference between resolving distinct parent grains and silently merging them.
The stricter test costs no sensitivity: true same-parent boundaries score zero
to within the $10^{-6}$-degree round-trip floor.

The set is also a useful object in its own right. The Kurdjumov-Sachs
fingerprint contains the $\Sigma 3$ twin relation — $60^{\circ}$ about
$\langle 111 \rangle$, Morito's V1–V20 intervariant pair — which is why
twin-related martensite variants are so common; this identity is pinned as a
worked example.

## Fitting the operative relationship

`fit_orientation_relationship(parents, children, nominal)` goes one step
further and refines the relationship itself: each measured parent-to-child
map is aligned to the current estimate through both symmetry groups, the
aligned rotations are averaged with the quaternion eigen-mean, and the steps
iterate to convergence. Starting from a Kurdjumov-Sachs nominal on
Greninger-Troiano data, the fit recovers GT exactly and reports the
$2.40^{\circ}$ distance from the assumed nominal; on noisy data the
residual statistics quantify the fit. The returned
`OrientationRelationshipFitReport` carries the fitted
`OrientationRelationship` and a `describe()` summary.

## Parallelism finders

`find_parallel_planes(relationship, parent_plane)` enumerates the parent
plane's symmetry family and reports, per variant, every child plane within an
angular tolerance of exact parallelism — under Kurdjumov-Sachs each of the 24
variants pairs exactly one $\{111\}$ member with a $\{011\}$ child plane
at zero deviation (its close-packed plane), and
`find_parallel_directions(...)` does the same for the
$\langle 110 \rangle \parallel \langle 111 \rangle$ directions. The result
is a typed `ParallelismReport` whose `describe()` prints the parallelism
table.

## Explainable reports

Every OR surface explains itself in convention-explicit prose:
`OrientationRelationship.describe()` states the phases and point groups, the
defining parallelisms, the misorientation representative, and the variant
count; correspondence results, deviation reports, parallelism reports,
variant selections, and parent reconstructions all carry `describe()`
methods whose statements are validated by tests (the explainable-results
doctrine of the development guide).

## Variants

`generate_variants()` enumerates the crystallographically distinct child
orientations (Bain 3; NW, Pitsch, Burgers 12; KS, GT 24), and
`intervariant_misorientations(...)` gives the axis/angle table between
variants. Mapping the parent close-packed plane across all 24 KS variants
shows the expected physics: exactly the six variants sharing that
close-packed plane return $\{011\}$ with zero residual; the other eighteen
land on irrational images.

## Verified numerical examples

The executable worked examples in the
[transformation gallery](../examples/generated/transformation.md) compute the
Kurdjumov-Sachs plane correspondence and the Bain direction correspondence
live from the code and check them against their defining parallelisms.

## The notebook teaching track

Three executed tutorial notebooks walk this whole surface with rendered
outputs, diagrams, and the underlying equations:

- {doc}`../tutorials/notebooks/18_orientation_relationships_fundamentals` —
  conventions, KS variants, intervariant spectrum, packets, variant pole
  figures, deviation and fitting;
- {doc}`../tutorials/notebooks/19_lattice_correspondence_and_transformation_strain`
  — index correspondence, the Bain strain rendered in 3D, polar rotations,
  the Burgers basal-on-(110) overlay;
- {doc}`../tutorials/notebooks/20_or_catalogs_identification_and_reconstruction`
  — every standard catalog with computed separations, then identification,
  rotation refinement, and parent-grain reconstruction end to end.

The scientific diagrams they embed (doctrine pipeline, variant-generation
algorithm, Bain correspondence cells, the F7/F8 pipeline, and the Burgers
family geometry) live under `docs/site/_static/or/` and are reusable in
reports and slides.
