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

1. **The rigid rotation** \(\mathbf{R}\) between parent and child Cartesian
   crystal frames (`parent_to_child_rotation`). It maps unit vectors and
   composes with orientations.
2. **The index correspondence** — the linear maps that carry Miller indices:
   \(\mathbf{u}_{c} = \mathbf{M}\,\mathbf{u}_{p}\) for directions with
   \(\mathbf{M} = \mathbf{A}_{c}^{-1}\mathbf{R}\,\mathbf{A}_{p}\) (direct
   structure matrices), and \(\mathbf{h}_{c} = \mathbf{M}^{*}\,\mathbf{h}_{p}\)
   for planes with \(\mathbf{M}^{*} = \mathbf{M}^{-\mathsf{T}}\) (reciprocal
   bases). The inverse-transpose relation preserves the zone law
   \(\mathbf{h} \cdot \mathbf{u}\), so a direction lying in a plane stays in
   the mapped plane. Correspondence matrices are generally not rotations and
   generally irrational.
3. **The transformation deformation** (Bain strain and friends) — reserved for
   a future release; the API name space is planned so nothing needs renaming.

## Constructing named relationships

`OrientationRelationship` ships literature constructors — Bain,
Kurdjumov-Sachs, Nishiyama-Wassermann, Greninger-Troiano, Pitsch (fcc-bcc)
and Burgers (bcc-hcp) — plus the generic
`from_parallel_plane_direction(...)` builder for user-defined ORs. Catalogs
(`standard_fcc_bcc_relationships`, `standard_bcc_hcp_relationships`) bundle
them per phase pair.

## Mapping indices across the relationship

The correspondence surface answers the canonical OR questions:

- `map_plane_to_child(plane)` / `map_direction_to_child(direction)` — which
  child \((hkl)\) / \([uvw]\) corresponds to this parent plane or direction?
- `map_plane_to_parent(...)` / `map_direction_to_parent(...)` — the inverse,
  for reading product-phase measurements against parent stereography.
- `correspondence_direct()` / `correspondence_reciprocal()` — the raw
  \(\mathbf{M}\) and \(\mathbf{M}^{*}\) matrices, per relationship or per
  variant.
- `map_plane_across_variants(...)` / `map_direction_across_variants(...)` —
  the variant-resolved tables.

Every mapping returns the **exact** (generally irrational) target components,
the nearest primitive-integer **rationalization** (bounded index search, sign
sensitive, default bound 17 to cover the Greninger-Troiano
\(\langle 5\,12\,17 \rangle\) family), and the **angular residual** between
them, so "nearly parallel" is always a quantified statement rather than an
implicit rounding.

Hexagonal phases participate with full index meaning: Burgers maps
\((110)_{\beta}\) to the basal plane \((0001)_{\alpha}\) and
\([\bar{1}11]_{\beta}\) to \([11\bar{2}0]_{\alpha}\) (stored as the
three-index \([110]\)).

## The OR as a misorientation, and deviation from it

`misorientation()` returns the symmetry-reduced (disorientation)
representative of the relationship — the way ORs are measured and reported
from EBSD boundary data. For Kurdjumov-Sachs this is the published
\(42.85^{\circ}\) rotation about \(\langle 0.968\;0.178\;0.178 \rangle\);
Nishiyama-Wassermann gives \(45.99^{\circ}\) and Bain \(45^{\circ}\) about
\(\langle 100 \rangle\).

`or_deviation(parents, children, relationship)` quantifies how well measured
parent/child orientation pairs obey a nominal OR: for each pair it reports
the smallest child-symmetry-reduced angle to any variant prediction, plus the
winning variant index. Exact synthetic data returns zeros; children generated
with Greninger-Troiano deviate by the documented \(2.40^{\circ}\) from
Kurdjumov-Sachs and \(2.86^{\circ}\) from Nishiyama-Wassermann — the report's
aggregate statistics are the entry point for OR fitting.

## Parallelism finders

`find_parallel_planes(relationship, parent_plane)` enumerates the parent
plane's symmetry family and reports, per variant, every child plane within an
angular tolerance of exact parallelism — under Kurdjumov-Sachs each of the 24
variants pairs exactly one \(\{111\}\) member with a \(\{011\}\) child plane
at zero deviation (its close-packed plane), and
`find_parallel_directions(...)` does the same for the
\(\langle 110 \rangle \parallel \langle 111 \rangle\) directions. The result
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
close-packed plane return \(\{011\}\) with zero residual; the other eighteen
land on irrational images.

## Verified numerical examples

The executable worked examples in the
[transformation gallery](../examples/generated/transformation.md) compute the
Kurdjumov-Sachs plane correspondence and the Bain direction correspondence
live from the code and check them against their defining parallelisms.
