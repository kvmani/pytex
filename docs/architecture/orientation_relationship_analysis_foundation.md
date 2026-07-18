# Orientation Relationship Analysis Foundation

This document defines the doctrine and feature program for orientation-relationship (OR)
analysis in PyTex. OR analysis is the designated flagship of the library: the capability set
that lets a researcher move between parent and product phases — indices, orientations, variants,
maps — with explicit crystallographic semantics at every step.

It extends, and is normative over, the transformation portions of the
[Phase Transformation Foundation](phase_transformation_foundation.md). Sequencing and priority
live in the
[Critical Review And Development Guide](../roadmap/critical_review_and_development_guide.md).

## 1. What Exists Today (updated 2026-07-17, post Cycles A-B)

Implemented and tested in `pytex.core.transformation`, `pytex.core.parent_reconstruction`, and
`pytex.experimental`:

- `OrientationRelationship` with named constructors: Bain, Kurdjumov-Sachs (KS),
  Nishiyama-Wassermann (NW), Greninger-Troiano (GT), Pitsch (fcc↔bcc), Burgers (bcc↔hcp), Shoji-Nishiyama (fcc↔hcp, 4 variants), and
  Pitsch-Schrader (hcp↔bcc, 3 variants, 5.26 deg from inverse Burgers), and
  Potter (hcp↔bcc, 12 variants, exact {10-11}||{110} pyramidal parallelism a
  c/a-dependent ~1-2 deg from inverse Burgers), and Bagaryatsky and Isaichev
  (ferrite↔cementite in the Pnma setting, 12 and 24 variants, separated by a
  ~3.6 deg rotation about the cementite a-axis), plus the generic
  `from_parallel_plane_direction(...)` correspondence constructor and standard
  catalogs.
- Variant machinery: symmetry-reduced generation with literature-correct counts, intervariant
  misorientation tables (Morito-validated), variant-selection scoring, and packet
  classification (`variant_close_packed_groups`: KS 4x6 packets, Burgers 6x2).
- **F1-F3 (index correspondence):** `correspondence_direct`/`correspondence_reciprocal`,
  `map_plane_to_child`/`map_direction_to_child` and parent-inverses, across-variant tables,
  with rationalized indices and atan2 angular residuals; hexagonal Miller-Bravais preserved.
- **F4 (parallelism finders):** `find_parallel_planes`/`find_parallel_directions` over
  symmetry families with typed reports.
- **F5 (misorientation + deviation):** `misorientation()` representatives (KS 42.85 deg
  <0.968 0.178 0.178> pinned) and `or_deviation` with best-variant assignment.
- **F6 (fitting):** `fit_orientation_relationship` — symmetry-aligned quaternion eigen-mean
  with iterative realignment; recovers GT exactly from a KS nominal.
- **F8 (reconstruction, experimental):** `reconstruct_parent_grains` /
  `reconstruct_parent_grains_from_graph` — intervariant-fingerprint edges, union-find
  clustering, averaged parent refinement, EBSD grain-graph wiring; validated on the
  24-variant lath-martensite structure fixture.
- **F10 (variant pole figures):** `variant_pole_figure` + `plot_variant_pole_figure`, pinned
  by the packet-plane coincidence.
- **F12 (deformation gradients):** `deformation_gradient()` — nearest-integer lattice
  correspondence, parent-frame gradient with polar decomposition; Bain stretches
  (1.127, 1.127, 0.797) and the literature rigid rotations (KS 11.06 deg, NW 9.74 deg from
  Bain) pinned.
- Explainable `describe()` prose on every report; the canonical composition
  `g_child = g_parent ∘ V^T` is regression-pinned (see the development-guide changelog for
  the convention correction).

## 2. Doctrine

### 2.1 Rotation, correspondence, and deformation are three different objects

An OR analysis surface must keep three mathematically distinct objects explicit and never
conflate them:

1. **The rigid rotation** `R` between parent and child *Cartesian crystal frames* — what
   `parent_to_child_rotation` stores today. It maps unit vectors and orientations.
2. **The index correspondence** — the linear map between *lattice bases* that carries Miller
   indices. For direction indices `[uvw]`: `u_c = A_c⁻¹ R A_p u_p` where `A_p`, `A_c` are the
   direct structure matrices (crystal basis → Cartesian). For plane indices `(hkl)` the map goes
   through the reciprocal bases (equivalently, the inverse-transpose of the direction map).
   Correspondence is generally **not** a rotation matrix and is generally irrational; nearness to
   rational indices is a physical statement about the OR, not a given.
3. **The transformation deformation** — the full deformation gradient (e.g. the Bain strain plus
   rigid rotation for martensite), which requires both lattices' parameters. It is the object
   habit-plane (invariant-plane-strain / PTMC) analysis operates on.

The stable API names these three explicitly (`rotation`, `correspondence`, `deformation`) so a
user can never accidentally push Miller indices through a bare rotation.

### 2.2 ORs are also misorientations

Every OR must be expressible as a symmetry-reduced misorientation (disorientation axis/angle in
stated frames), because that is how ORs are measured, compared, and reported in the literature
(e.g. KS ≈ 42.85° about <0.968 0.178 0.178>). The misorientation representation is the bridge to
EBSD boundary data and to OR fitting.

### 2.3 Variant identity is a convention, not an accident

Variant numbering must be stable, documented, and — where a de facto literature standard exists
(Morito V1–V24 for KS in steels; the standard 12-variant α ordering for Burgers in Ti) —
conforming or explicitly declared non-conforming. Variant tables published from PyTex must be
reproducible from the doctrine alone. Variant grouping structure (Bain groups and
close-packed-plane groups for KS; the six β↔α pair classes for Burgers) is part of variant
identity and must be first-class queryable metadata, not a user-side recomputation.

### 2.4 Every OR result explains itself

Per the explainable-results doctrine in the development guide: OR objects and OR reports carry
`describe()` surfaces that state the defining parallelisms (`{111}_γ ∥ {011}_α′`,
`<-101>_γ ∥ <-1-11>_α′`), the misorientation representation, the variant count and grouping, and
any deviations — with registry terminology and canon citations.

## 3. Feature Program

Phases follow the development-guide cycles. Each feature lands with theory note, worked example,
validation row, and `describe()` support.

### Phase 1 — Correspondence and index mapping (Cycle A)

**F1. Transformation-matrix evaluation surface.** Expose, on `OrientationRelationship` and
`TransformationVariant`:

- `rotation` (exists), returned in stated parent→child Cartesian sense;
- `correspondence_direct()` → the direction-index map (3×3, generally non-orthogonal), built
  from the two phases' structure matrices;
- `correspondence_reciprocal()` → the plane-index map (inverse-transpose relation);
- `misorientation()` → the symmetry-reduced axis/angle representation (§2.2);
- future `deformation()` (Phase 4) reserving the name now.

**F2. Parent→child index mapping, per variant.** Given a parent `CrystalPlane`,
`CrystalDirection`, or `MillerIndex` (Miller-Bravais accepted for hexagonal phases):

- `map_plane_to_child(plane, variant=...)` / `map_direction_to_child(...)` returning the child
  object with **exact** (irrational) components, the **rationalized** nearest low-integer indices
  under a configurable index bound, and the **angular residual** between exact and rationalized
  results. Batch forms over `MillerPlaneSet`/`MillerDirectionSet` and over all variants at once
  (returning a variant-indexed result table).

**F3. Child→parent inverse mapping.** The same surface through the inverse correspondence, so
users can ask "which parent plane became this product plane in variant Vk" — required for trace
analysis and for interpreting product-phase measurements against parent stereography.

**F4. Parallelism finders.** Given an OR and a parent plane/direction *family* (symmetry orbit):
enumerate, per variant, the child planes/directions within a tolerance of parallelism, with
deviations — the general machine behind statements like "{111}_γ ∥ {011}_α′ to within 0.0°" and
behind discovering near-parallelisms of non-defining families (e.g. which {hkl}_α′ lie near
{100}_γ under KS). Output is a typed report with a `describe()` that prints the parallelism
table.

**F5. OR-deviation metric.** `or_deviation(parent_orientations, child_orientations, relation)`
→ per-pair and aggregate angular departure of measured pairs from the nominal OR (minimum over
variants, symmetry-reduced). Zero on exact synthetic data; this is roadmap item I10 and the
entry point to fitting.

### Phase 2 — Fitting and reconstruction (Cycle B)

**F6. OR fitting from measured pairs.** Estimate the operative OR from parent/child orientation
pairs by symmetry-aware rotation averaging (quaternion mean over the per-pair variant-aligned
misorientations, with outlier rejection); report confidence via residual statistics. Parity
floor: MTEX `calcParent2Child`.

**F7. OR determination without the parent.** From child-child boundary misorientations in a
fully transformed microstructure, refine the OR via the MODF peak structure (the intervariant
fingerprint). Validated against lath-martensite literature data.

**F8. Map-scale parent-grain reconstruction.** The flagship: variant-graph voting on the grain
boundary network (`ebsd/models.py` grains + `merge_by_csl`-style union-find), iterative
parent-grain growth, per-grain fit and ambiguity maps, twin-aware voting reusing `ebsd/csl.py`.
Validated on martensite→austenite and α→β Ti fixtures. Promotes the experimental scoring layer
into a stable, evidence-backed surface.

**F9. Variant analysis at map scale.** Variant maps (child grains colored by variant), variant
pair/block statistics, Bain/CP-group maps, variant-selection strength against the uniform
baseline, and (later) Patel-Cohen-type interaction scoring against an applied stress.

### Phase 3 — Visualization and interchange (with Phase 2)

**F10. Variant pole figures and stereographic overlays.** Predicted variant poles over measured
child pole figures; parent-frame IPF coloring of child data; OR parallelism stereograms for
teaching.

**F11. JSON contracts.** Canonical serialization for `OrientationRelationship` (defining
correspondences + rotation + provenance), variant sets, and OR reports, per the data-contracts
standard — reconstruction-grade, not lossy summaries.

### Phase 4 — Deformation and interface crystallography (Cycle C+)

**F12. Transformation deformation gradients.** Lattice-parameter-aware Bain/Burgers distortions;
principal strains; volume change.

**F13. Habit planes and PTMC.** Invariant-plane-strain analysis (lattice-invariant shear
options, habit-plane normals, orientation-relationship prediction), populating the existing
`habit_plane_pairs` slot with computed, provenance-carrying results; invariant-line analysis for
bcc/fcc precipitation.

**F14. OR catalog breadth.** Shoji-Nishiyama (fcc↔hcp), Pitsch-Schrader and Potter (hcp↔bcc),
Bagaryatsky/Isaichev (ferrite↔cementite), and user-defined literature ORs with citation
metadata; the catalog constructors gain a registry keyed by name + phase-family with provenance.

## 4. Validation Program

- Analytic: exact synthetic round-trips (index maps must invert; variant orbits closed under
  parent symmetry; correspondence ∘ inverse = identity to 1e-12).
- Literature: KS/NW/GT/Pitsch/Burgers defining parallelisms and misorientation representatives;
  Morito intervariant and variant-grouping tables; Burgers 12-variant β↔α tables; GT's position
  between KS and NW (2.40°/2.86°).
- Parity: MTEX `calcParent2Child` / parent-grain-reconstruction outputs on shared fixtures;
  rows land in `docs/testing/mtex_parity_matrix.md` and the phase-transformation validation
  matrix.
- Property-based: random lattices and ORs — mapping a plane and its symmetric equivalents must
  produce the same child family; rationalization residuals must be invariant to index scaling.

## 5. Current Limits (honest statement, updated 2026-07-18)

F1-F12 and F14 are implemented and validated as listed in §1. F7 is implemented at both stages
as experimental surfaces: `pytex.experimental.identify_orientation_relationship` ranks candidate
ORs against child-child boundary misorientations via the double-coset intervariant fingerprint,
and `pytex.experimental.refine_orientation_relationship_from_boundaries` refines the winning
rotation by alternating nearest-coset-element assignment with least-squares rotation updates —
both with no parent orientations required. Still **not available**, and no PyTex document may
imply otherwise: habit-plane / PTMC analysis (F13). The
F14 catalog breadth is complete (Shoji-Nishiyama, Pitsch-Schrader, Potter, Bagaryatsky, and
Isaichev all landed). Parent-grain reconstruction (F8) is
experimental: synthetic and literature-structure validation exists, but external measured-data
fixtures and MTEX parity are still required before stabilization. JSON contracts for the newer
report objects (F11) are partial.

## References

### Normative

- [Critical Review And Development Guide](../roadmap/critical_review_and_development_guide.md)
- [Phase Transformation Foundation](phase_transformation_foundation.md)
- [Canonical Data Model](canonical_data_model.md)
- [Notation And Conventions](../standards/notation_and_conventions.md)

### Informative

- Morito, S. et al., "The morphology and crystallography of lath martensite in Fe-C alloys"
  (intervariant tables, V1–V24 convention)
- Burgers, W. G., "On the process of transition of the cubic-body-centered modification into
  the hexagonal-close-packed modification of zirconium"
- Bhadeshia, H. K. D. H., *Geometry of Crystals, Polycrystals, and Phase Transformations*
- MTEX documentation: parent-grain reconstruction and `calcParent2Child`
