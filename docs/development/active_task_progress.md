# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Objective

Execute the **Critical Review And Development Guide**
(`docs/roadmap/critical_review_and_development_guide.md`) cycle by cycle, one phase per
commit+push, as a long-horizon mission. Cycle A (findings 1, 2, 6, 16, 17, 18, 23, 24) is fully
executed; now working Cycle B (findings 3, 4, 7, 19, 20). The OR feature definitions (F1–F14)
live in `docs/architecture/orientation_relationship_analysis_foundation.md`.

## Current Status

- Started: 2026-07-16
- Branch: `main` (tracking `origin/main`, push after each phase)
- Baseline commit: `d1a1561`; Phase 0 pushed as `84a9767`; Phase 1 as `557d5e1`; Phase 2 as
  `cf8391e`; Phase 3 as `cee9ee7`; Phase 4 as `2680679`; Phase 5 as `c3e34bf` (Cycle A done);
  Phase 6 as `b10d3b7`
- Phase 7 pushed as `d40fb14`; Phase 8 as `801af7e`; Phase 9 as `a5e3d73` (Cycle B done);
  Phase 10 as `de58dfe`; Phase 11 as `8ce6049`; Phase 12 as `3266d21` (convention bug fix)
- Phase 13 pushed as `cf91ce9`
- Phase 14 pushed as `ee72419`
- Phase 15 pushed as `551791e`
- Phase 16 pushed as `20dd2dd`
- Phase 17 pushed as `e621fe4`
- Phase 18 pushed as `8247554`
- Phase 19 pushed as `05baa47`
- Phase 20 pushed as `2faa0e7`
- Phase 21 pushed as `d34fef0`
- Phase 22 pushed as `a120111`
- Phase 23 pushed as `4567817`
- Phase 24 pushed as `c6e68e3`
- Phase 25 pushed as `32b225b`
- Phase 26 pushed as `8dd8a81`
- Phase 27 pushed as `d43bc9e` (F14 COMPLETE)
- Phase 28 pushed as `1e47c9b` (F7 COMPLETE, both stages)
- NEW GOAL (2026-07-18): OR scientific documentation program — plan and ledger in
  `docs/roadmap/working_notes_or_documentation_program.md` (phases 29-32: SVG diagram
  assets, three executed OR tutorial notebooks, cross-linking).
- Phase 29 pushed as `b125d40`; Phase 30 as `5a308f5`; Phase 31 as `89ae8ad`
- Phase 32 pushed as `f2efbfc` — OR documentation program v1 COMPLETE (see its
  working notes for outcomes and follow-ons).
- NEW GOAL (2026-07-20): Composite OR diffraction pattern program — vectorized
  kinematic composite SAED engine (parent + OR variants, any parent zone axis)
  with publication-grade configurable plotting and smart annotations. Plan and
  ledger in `docs/roadmap/working_notes_composite_saed_program.md` (phases
  CD0-CD6). Kinematic only; regression tests required for every critical part.
  Phases CD0 `f609c0a`, CD1 `e549bce`, CD2 `d292f1b`, CD3 `202b6da`, CD4
  `88c202e`, CD5 `1a93130`; CD6 committing now. Composite SAED program v1
  COMPLETE (see its working notes for outcomes and follow-ons).
- Next phase candidates (in roadmap order): F13 habit planes / PTMC (large; needs its own
  working notes), reconstruction stabilization (measured-data fixture + MTEX
  `calcParent2Child` parity via scripts/mtex_generators), finding 8 remainder (bump/fibre
  kernels + SO3FunHarmonic), finding 9 (PF->ODF ghost correction), finding 14 (Kikuchi).

## Phase Plan (each phase = verified commit + push)

- **Phase 0 — Foundational docs overhaul (this commit).** New development guide + OR foundation
  doc; mission/AGENTS refresh; docs index completion; stale phase-transformation claims fixed.
- **Phase 1 — Consistency fixes (finding 6).** One shared phase-semantic-match helper (today
  duplicated with different strictness in `core/transformation.py` vs
  `experimental/phase_transformation.py`); type `OrientationRelationship.parallel_directions` as
  `CrystalDirection` pairs instead of naked ndarrays (keep constructor compatibility).
- **Phase 2 — Index correspondence (F1–F3, finding 1).** `correspondence_direct()` /
  `correspondence_reciprocal()` on `OrientationRelationship`/`TransformationVariant`;
  `map_direction_to_child` / `map_plane_to_child` and child→parent inverses with exact
  components, rationalized nearest-low-integer indices (configurable bound), and angular
  residual; all-variant table form. Worked example (Bain/KS known correspondences), LaTeX theory
  note, validation-matrix row.
- **Phase 3 — OR as misorientation + deviation (F5, finding 2).** `misorientation()` (symmetry-
  reduced axis/angle; KS ≈ 42.85° about <0.968 0.178 0.178> as pinned test) and
  `or_deviation(parents, children, relationship)` (min-over-variant symmetry-reduced residuals;
  zero on synthetic data).
- **Phase 4 — Parallelism finders (F4) + `describe()` doctrine start (finding 16).**
  Family-orbit parallelism report; `describe()` on `OrientationRelationship`,
  `VariantSelectionReport`, `ParentReconstructionReport`, intervariant results; tests assert
  conventions + key numbers appear.
- **Phase 5 — Hygiene gates (findings 17, 18).** Autouse matplotlib figure-close fixture;
  `filterwarnings` policy; coverage report + ratchet in CI.

## Verification Gates (every phase, before commit)

- `python -m pytest` (804+ passing, no new warnings)
- `python scripts/check_repo_integrity.py`
- `python -m ruff check .`
- `python -m mypy src`
- Sphinx build when docs/site content changes: `python -m sphinx -b html docs/site docs/_build/html`

## Key facts established (verified against live code, 2026-07-16)

- OR machinery lives in `src/pytex/core/transformation.py` (699 lines): named ORs (Bain, KS, NW,
  GT, Pitsch, Burgers), `generate_variants` (child-symmetry orbit reduction, literature counts),
  `intervariant_misorientation*`, `PhaseTransformationRecord`. Catalogs + variant selection in
  `core/parent_reconstruction.py`; experimental scoring in `experimental/phase_transformation.py`.
- `map_parent_vector_to_child` is Cartesian-only; **no** Miller-index correspondence exists.
- `Lattice.direct_basis().matrix` / `reciprocal_basis().matrix` provide structure matrices;
  `CrystalPlane.normal` goes through the reciprocal basis; Miller-Bravais conversion exists on
  `MillerIndex`/`CrystalDirection`/`CrystalPlane`.
- `_phase_semantically_matches`: transformation.py version compares `symmetry.point_group`;
  experimental version compares full `symmetry` equality — must unify (single helper, decide
  strictness: use full `SymmetrySpec` equality? No — transformation.py's laxer point-group form
  exists so records constructed from phases with equal groups but distinct operator caches still
  match; verify call sites before choosing).
- Test suite: 804 passed, ~27 s, 117 warnings (matplotlib figure leaks among them). CI: ubuntu,
  py3.11, two lanes, no coverage.

## Completed

- Repo-wide critical review; wrote `docs/roadmap/critical_review_and_development_guide.md`
  (25 findings, priorities, explainable-results doctrine, quality bars).
- Wrote `docs/architecture/orientation_relationship_analysis_foundation.md` (doctrine: rotation
  vs correspondence vs deformation; features F1–F14; validation program).
- Refreshed `mission.md` (OR flagship, explainability principle 10, new success criteria) and
  `AGENTS.md` (new primary references, explainability/warnings/convention-pinning rules).
- Completed `docs/README.md` index; fixed stale variant-doctrine claims in
  `docs/architecture/phase_transformation_foundation.md`.
- Verified: integrity check passed, 804 tests passed, ruff clean.

## Phase 1 outcomes (2026-07-16)

- **Found and fixed a latent crash:** `SymmetrySpec.__eq__` (dataclass-generated) raised
  `ValueError` on distinct-but-equal instances because of the ndarray `operators` field; it was
  also unhashable. Now `eq=False` with explicit `__eq__` (name, point_group, specimen_symmetry,
  reference_frame, `np.array_equal` operators; provenance excluded) and `__hash__`.
- Unified phase identity into public `pytex.core.lattice.phases_semantically_match` (None-safe,
  normalized point-group comparison; exported from `pytex.core` and top level). Both duplicated
  `_phase_semantically_matches` copies deleted; all `plane.phase != phase` checks in
  transformation.py now use the helper (avoids latent `Phase.__eq__` array ambiguity too).
- `OrientationRelationship.parallel_directions` now stores typed
  `CrystalDirection` pairs (index meaning preserved, e.g. KS <-101>/<-1-11>, Burgers
  Miller-Bravais-derived [110]); raw Cartesian 3-vectors still accepted and wrapped via new
  `CrystalDirection.from_cartesian` (public inverse of `unit_vector`). Phase membership
  validated. JSON contract now emits typed crystal-direction payloads and still reads legacy
  float-triple payloads.
- Updated consumers: `plotting/scene3d.py`, scene3d test, visualization worked example
  (gallery regenerated). New tests: SymmetrySpec equality/hash; phases_semantically_match;
  typed/legacy/mismatch/Burgers parallel-direction cases. 810 passed; ruff/mypy/integrity green.
- Note: `Phase.__eq__` (and other array-field dataclasses) still have the ambiguous-truth
  hazard when comparing distinct-but-equal instances — flagged for a later dedicated pass.

## Phase 2 outcomes (2026-07-16)

- Index-correspondence surface landed in `core/transformation.py`:
  `correspondence_direct()` (`M = A_c^-1 R A_p`), `correspondence_reciprocal()`
  (`M* = M^-T`, zone-law preserving), `map_direction_to_child/parent`,
  `map_plane_to_child/parent` (each optionally per `variant=`), module functions
  `map_direction_across_variants` / `map_plane_across_variants`; result types
  `DirectionCorrespondence` / `PlaneCorrespondence` with exact components, primitive-integer
  rationalization (bounded search, default `DEFAULT_RATIONALIZATION_MAX_INDEX = 17` for GT
  <5 12 17>), and atan2 angular residuals (arccos floors at ~8.5e-7 deg — use atan2 pattern
  for any near-zero angle work).
- Verified physics: KS (111)→(011) and [-101]→[-1-11] exact; Bain [110]→[100]; Burgers
  (110)→(0001) and [-111]→[11-20]; inverse-transpose + zone-law identities; round trips;
  across 24 KS variants exactly the 6 CP-group variants map (111)γ to {011} (residual 0),
  the other 18 land on irrational images — pinned in tests (9 new tests).
- Docs: registry symbols \(\mathbf{M}\), \(\mathbf{M}^{*}\); theory note
  `docs/tex/algorithms/orientation_relationship_index_correspondence.tex`; concept page
  `docs/site/concepts/orientation_relationships.md` (+ toctree); site include-stubs for the OR
  foundation and the development guide (fixed all Sphinx xref warnings — root-level files are
  referenced as backticked paths, not links, in site-rendered docs); worked-example group
  `transformation` (KS plane + Bain direction identities, gallery regenerated);
  validation-matrix row (implemented).
- Gates: 819 passed; ruff/mypy/integrity green; Sphinx build zero warnings.

## Phase 3 outcomes (2026-07-16)

- `OrientationRelationship.misorientation()` → `Misorientation` (built on the existing
  deterministic `Misorientation.disorientation()` fundamental-zone representative; child
  symmetry left, parent right). Verified against literature: KS 42.848 deg <0.968 0.178 0.178>,
  NW 45.99 deg <0.976 0.201 0.083>, GT 44.23 deg, Bain 45 deg <100>.
- `or_deviation(parents, children, relationship)` + `ORDeviationReport` (per-pair min-over-
  variant child-symmetry-reduced deviations, best-variant indices, mean/median/max). Verified:
  zero + planted-variant recovery on exact GT synthetic children; GT children deviate 2.404 deg
  from KS and 2.861 deg from NW (the documented separations); Bain ≈ 10.15 deg.
- Composition convention confirmed: predicted child = `V @ P` (matches
  `predicted_child_orientations`; `Rotation.compose` is quaternion left-multiplication).
- Note: `OrientationSet[slice]` does NOT return a sub-set (returns Orientation with bad shape)
  — construct sliced sets manually from `.quaternions[...]`; potential later API improvement.
- Docs: worked example `or-ks-misorientation-representation` (42.85/<0.968 0.178 0.178>,
  Verlinden et al. citation); concept-page section; two validation-matrix rows;
  specifications.md Transformation primitive list expanded to the full current surface.
- Gates: 824 passed; ruff/mypy/integrity/Sphinx green.

## Phase 4 outcomes (2026-07-16)

- Parallelism finders: `find_parallel_planes` / `find_parallel_directions` +
  `ParallelismMatch`/`ParallelismReport`, built on `_integer_index_orbit` (symmetry orbit of
  integer indices via Cartesian operator action + integer recovery, antipodal-collapsed).
  Verified: KS pairs exactly one {111} member with a {011} child per variant (24 matches,
  0 deviation), same for <110>||<111> directions.
- `describe()` doctrine landed on: `OrientationRelationship` (phases/point groups, defining
  parallelisms via notation formatters, misorientation representative, variant count),
  `DirectionCorrespondence`, `PlaneCorrespondence`, `ORDeviationReport`, `ParallelismReport`,
  `VariantSelectionReport`, `ParentReconstructionReport`. Substring-validated in tests.
- mypy pattern: notation formatters require `Sequence[int]` — wrap ndarrays with the
  `_index_tuple` helper in transformation.py.
- Docs: concept-page sections (parallelism finders, explainable reports); two validation-matrix
  rows; specifications.md Transformation list extended.
- Gates: 828 passed; ruff/mypy/integrity/Sphinx green.

## Phase 5 outcomes (2026-07-16)

- `tests/conftest.py` autouse fixture closes all matplotlib figures per test (kills the
  "More than 20 figures" leak class).
- `filterwarnings = ["error", ...]` in pyproject with three narrow pymatgen exemptions
  (CIF chatter + dict-interface deprecation).
- **spglib gotcha:** spglib's dict-interface DeprecationWarning is force-enabled inside its own
  filter context — caller/pytest filters cannot silence it. Solved at the adapter boundary:
  `_spglib_dict_shim_silenced()` context in `core/lattice.py` (records, drops only that
  message, re-emits everything else) wrapped around `get_point_group_symbol` and the
  space-group accessor calls.
- Suite now runs **828 passed, zero warnings**. CI base lane gained the coverage ratchet
  (`--cov=pytex --cov-fail-under=87`; measured 88.01%). Ratchet-up policy: raise the floor
  when measured coverage rises, never lower it.
- Development-guide changelog updated: Cycle A findings 1, 2, 6, 16, 17, 18, 23, 24 all closed.

## Phase 6 outcomes (2026-07-17, Cycle B start)

- **Finding 7 closed (vectorization):** `intervariant_misorientations` now computes all pair
  relatives and the full symmetry-product tensor in single einsums (`triu_indices` pairs;
  representative-selection axis order preserved, so results are bit-identical to the historical
  per-pair enumeration); `PhaseTransformationRecord.predicted_child_orientations` and
  `select_variants` now compose variant × parent in matrix space (`einsum` +
  `OrientationSet.from_matrices`) instead of per-element quaternion loops.
- One test updated with justification: predicted-quaternion comparison is now sign-insensitive
  (q and -q are the same rotation; matrix→quaternion conversion picks the canonical branch).
  The pinned behavior is the rotation, not the incidental sign.
- **Finding 19 closed (CI matrix):** base lane now runs ubuntu+macos × Python 3.11–3.13
  (fail-fast off; docs build gated to ubuntu/3.11); classifiers extended to 3.12/3.13.
- Gates: 828 passed, zero warnings; ruff/mypy/integrity green.

## Phase 7 outcomes (2026-07-17)

- **Finding 3 / F6 landed: `fit_orientation_relationship(parents, children, nominal)`** in
  `core/transformation.py` + `OrientationRelationshipFitReport` (fitted OR named
  `<nominal>_fitted`, per-pair residuals, iterations/convergence, symmetry-reduced
  `deviation_from_nominal_deg`, `describe()`).
- Algorithm: precompute all symmetry-equivalent descriptions `S_c (C P^T) S_p` per pair once;
  iterate [align each pair to current estimate by max trace → quaternion eigen-mean (Markley,
  scatter-matrix eigenvector — sign-free) → update] until step < tol.
- Verified: exact GT pairs + GT nominal → zero everything, 1 iteration; exact GT pairs + **KS
  nominal → recovers GT exactly** (distance to true GT = 0, reported nominal distance 2.404
  deg); seeded 0.5-deg noise → fit within 0.15 deg of truth, residuals ≈ noise. 4 new tests.
- Docs: validation-matrix row (MTEX `calcParent2Child` parity comparison + worked example
  explicitly queued); concept-page section; specifications list. Gates: 832 passed, zero
  warnings; ruff/mypy/integrity/Sphinx green.

## Phase 8 outcomes (2026-07-17)

- **Finding 20 started (property-based testing):** `tests/unit/test_property_based.py` with
  Hypothesis (added to the dev extra, >=6.100): rotation composition-inverse law, proper
  orthogonality of quaternion-derived matrices, Miller-Bravais direction round trips (ray
  preserved under GCD reduction; plane round trips exact), correspondence-matrix inversion on
  arbitrary integer planes through any KS variant, zone-law preservation for arbitrary
  plane/zone-direction pairs, and phase-match reflexivity/symmetry. 50 examples per property,
  deadline disabled for CI stability.
- Gates: 839 passed, zero warnings; ruff/mypy/integrity green.

## Phase 9 outcomes (2026-07-17)

- **Finding 4 / F8 v1 landed (experimental):**
  `pytex.experimental.reconstruct_parent_grains(children, adjacency, relationship)` +
  `ParentGrainReconstructionResult` in `experimental/parent_grain_reconstruction.py`.
  Pipeline: (1) edge test — child-symmetry-reduced boundary disorientation matches the
  intervariant fingerprint (0 plus all distinct pair angles from
  `intervariant_misorientation_angles_deg`) within tolerance; (2) union-find clustering over
  linked edges; (3) per cluster, candidate parents `V_k^T C_first` scored against every member
  (min-over-variants disorientation), best candidate wins; per-grain residuals + `describe()`
  with explicit singleton-ambiguity and experimental-status caveats.
- Verified: 3 planted KS parents × 5 children — exact partition recovery, both cross-parent
  edges rejected, parents to 0 deg; 0.3-deg noise keeps the partition (parents within 1.5 deg —
  v1 estimate carries first-member noise; averaging refinement is the queued v2 improvement).
- Validation-matrix row moved planned → foundational with the remaining-work statement
  (literature fixtures, EBSD grain-graph wiring, candidate averaging).
- Gates: 843 passed, zero warnings; ruff/mypy/integrity green.

## Phase 10 outcomes (2026-07-17)

- Parent-reconstruction refinement landed: `_refine_cluster_parent` aligns every member's
  candidate descriptions `S_p V_k^T C_i` to the seed (max trace) and quaternion-eigen-means
  them. Noise-case parent errors drop from ~1.5 deg (first-member inheritance) to
  0.08–0.20 deg (sigma/sqrt(n) behavior); test assertions tightened to <0.5 deg parents,
  <1.0 deg mean residual. Exact case remains exact.
- Development-guide changelog records Cycle B as executed.
- Gates: 843 passed, zero warnings; ruff/mypy/integrity green.

## Phase 11 outcomes (2026-07-17)

- **EBSD wiring landed:** `reconstruct_parent_grains_from_graph(graph, relationship)` takes a
  `GrainSegmentation.grain_graph()` directly — grain-mean orientations become children, graph
  edges become adjacency; result rows follow `graph.node_grain_ids` and the new optional
  `grain_ids` field records the mapping. Phase check relaxed for phase-less map data (accepts
  matching point-group symmetry).
- End-to-end map test: 2x4 pixel `CrystalMap` -> segment_grains -> grain_graph ->
  reconstruction recovers the two planted parents with the cross-parent boundary rejected.
- **Fixture-design learnings (important):** (a) a cube-oriented parent makes some KS variant
  pairs symmetry-degenerate as child orientations — use general parent Eulers in fixtures;
  (b) cross-parent boundaries can coincidentally sit within tolerance of the intervariant
  fingerprint (real reconstruction ambiguity) — fixtures must pick parents whose cross
  boundary is verified far from the fingerprint (used (20,30,40) and (65,20,50): 5.0 deg away).
- **Follow-up to investigate:** `OrientationSet.misorientation_angles_to(symmetry_aware=True)`
  returned 39.92 deg for a same-parent variant pair whose boundary disorientation (C_i C_j^T,
  child ops both sides — the intervariant/table convention) is 57.21 deg. Two different
  relative-rotation conventions coexist; determine whether both are intended semantics
  (orientation-space vs boundary misorientation) and document, or reconcile.

## Phase 12 outcomes (2026-07-17) — REAL BUG FOUND AND FIXED

- The Phase 11 follow-up investigation resolved decisively: PyTex's normative convention is
  **orientation = crystal→specimen** (notation standard; `map_crystal_vector = R @ v`), so
  `misorientation_angles_to` (relative `R_i^T R_j`) was correct, and the transformation
  stack's `child = V @ P` composition was **wrong** (correct: `g_child = g_parent ∘ V^T`,
  which makes corresponding crystal directions coincide in specimen space). Diagnostic on
  canonically built children showed wrong variant selections with 14–16 deg residuals.
- Fixed compositions/relatives in: `predicted_child_orientations`, `select_variants`,
  `or_deviation` (also crystal-frame relative `C^T pred`), `fit_orientation_relationship`
  (measured `V = C^T P`), experimental scoring, parent-grain reconstruction (edge relative
  `C_i^T C_j`; candidates `P = C V`; refinement equivalents right-multiplied `P S_p`).
- Boundary/intervariant tables are convention-invariant (`C_i^T C_j = V_a V_b^T`), so all
  Morito/literature pins stand unchanged.
- Fit convergence made robust: fixed-point criterion = stable alignment assignments (the
  quaternion↔matrix round trip has a ~1.2e-6 deg noise floor that the step-angle test alone
  can't cross).
- Tests: all synthetic builders now construct children canonically; parent-equivalence
  comparisons use RIGHT multiplication (`P S_p`); new regression test pins the
  specimen-space parallelism identity + planted-variant recovery. 845 passed, zero warnings.
- Docs: concept-page "Composition convention" section; development-guide changelog records
  the finding with the lesson (synthetic tests must build inputs through the canonical
  convention, not the code-under-test's own composition).

## Phase 13 outcomes (2026-07-17)

- **F9 packet classification landed:** `variant_close_packed_groups(relationship,
  parent_plane)` (stable surface, `core/transformation.py`) labels each variant by the parent
  family member it carries into exact parallelism. Validated against the Morito hierarchy:
  KS + {111} → four packets of six variants; Burgers + {110} → six groups of two.
- **Lath-martensite literature-structure fixture:** one austenite parent with all 24 KS
  variants as children — reconstruction gathers all 24 into one parent recovered exactly,
  variant selection recovers every planted index 1..24, and packet labels come out 4x6.
  End-to-end validation of reconstruction + selection + packet classification on one fixture
  (external measured-data fixtures remain the queued step before stabilization).
- Two validation-matrix rows; concept-page "Variant packets" section; specifications entry.
- Gates: 847 passed, zero warnings; ruff/mypy/integrity/Sphinx green.

## Phase 14 outcomes (2026-07-17)

- `OrientationSet.__getitem__` now supports slices (returns a metadata-preserving
  `OrientationSet`; typed via `@overload`, int still returns `Orientation`). Closes the
  API wart logged in Phase 3. Gates: 848 passed, zero warnings; mypy/ruff/integrity green.

## Phase 15 outcomes (2026-07-17)

- **F10 landed:** `variant_pole_figure(parent, relationship, child_plane)` (stable core
  surface) computes specimen-frame poles of the child plane family per variant under the
  canonical composition `C = P V^T`, returned as a typed `VariantPoleFigure` with
  `describe()`; `plot_variant_pole_figure(...)` (plotting/runtime) renders the
  color-per-variant stereographic overlay on `plot_stereographic_vectors`.
- Pinned physics: every KS variant's predicted {011} pole set contains the specimen-frame
  normal of its packet's parent {111} member (all 24 variants exact) — ties F10 to the
  packet classification.
- **OR-fitting worked example** `or-fit-recovers-gt-from-ks-nominal` added to the
  transformation gallery (zero residual + 2.4037 deg KS-GT separation, computed live);
  closes the queued docs gap. Gallery regenerated.
- Gates: 851 passed, zero warnings; ruff/mypy/integrity/Sphinx green.

## Next Actions (Cycle C+, per the development guide §3 and world-class roadmap)

## Phase 28 outcomes (2026-07-18) — F7 complete (both stages)

- **F7 second stage (experimental):**
  `pytex.experimental.refine_orientation_relationship_from_boundaries` +
  `ORRefinementReport` in new `experimental/or_refinement.py`. Algorithm:
  alternate (1) per-edge nearest-element assignment over the FULL
  (child-op, parent-op, child-op) coset triple enumeration at the current
  rotation (edge-blocked einsum traces) with (2) scipy LM least-squares on a
  3-parameter left-rotation update, residuals = 2 sin(theta/2) chordal form
  (smooth at zero — exact data is a regular point, unlike arccos);
  convergence = stable assignments + step below tolerance. Rotation
  identifiable only up to coset symmetry: all reported distances
  symmetry-reduced.
- Verified: exact GT boundaries + KS nominal -> recovers the true GT
  rotation (distance ~1.2e-6 deg, the known matrix<->quaternion round-trip
  noise floor; tolerance 1e-5) with the 2.404 deg KS-GT update reported;
  0.3-deg-noise KS boundaries + NW nominal -> within 0.5 deg of true KS.
- Stale-claims sync: or_identification docstring + describe() now point to
  the refinement surface instead of denying it exists (test substring
  updated); foundation §5 marks F7 complete at both stages; validation row.
  Remaining OR program: F13 (habit planes/PTMC), reconstruction
  stabilization (measured-data fixtures, MTEX parity).

## Phase 27 outcomes (2026-07-18) — F14 complete

- **F14 final slice:** `from_bagaryatsky_correspondence` and
  `from_isaichev_correspondence` — the first orthorhombic-child constructors
  (cementite, Pnma/Lipson-Petch setting b > a > c; guard dict gained
  "222": "orthorhombic"). Bagaryatsky pins all three Bhadeshia (MST 34, 2018)
  axis parallelisms exactly ([1-1-1]->[100], [211]->[010], (0-11)->(001));
  Isaichev pins (101)->(031) + shared [1-1-1]->[100]. Symmetry-reduced
  separation 3.586 deg about EXACTLY the cementite a-axis (eigenvector pin;
  literature ~3.8 deg is axial-ratio-dependent). Variants: 12 (Bagaryatsky:
  parent 180@[0-11] pairs with child 180@c) and 24 (Isaichev: irrational
  (031) breaks that stabilizer) — internally derived counts, recorded as
  such. New `standard_ferrite_cementite_relationships` catalog, exported
  from `pytex.core` and top level. The F14 OR catalog program (S-N, PS,
  Potter, Bagaryatsky, Isaichev) is now COMPLETE; foundation doc §5 updated.
  Remaining OR program: F7 second stage (boundary-based rotation
  refinement), F13 (habit planes / PTMC), reconstruction stabilization.

## Phase 26 outcomes (2026-07-18)

- **F14 slice:** `from_potter_correspondence` — exact pyramidal parallelism
  {10-11}_hcp || {110}_bcc with the Burgers close-packed direction pairing
  <2-1-10>_hcp || <111>_bcc (Potter 1973, V-N precipitates in alpha-vanadium;
  622-parent/432-child guards); appended to `standard_hcp_bcc_relationships`
  (names now PS, inverse Burgers, Potter). Pinned: exact pyramidal + direction
  parallelisms; basal-image residual == symmetry-reduced separation from
  inverse Burgers (1.370 deg at c/a = 4.68/2.95 — the literature "~2 deg" is
  c/a-dependent, recorded as such); 12 variants (internally derived orbit
  count). Remaining F14: Bagaryatsky/Isaichev (need orthorhombic cementite
  support — assess before implementing).

## Phase 25 outcomes (2026-07-17)

- **F14 slice:** `from_pitsch_schrader_correspondence` — the first hexagonal-PARENT
  constructor ((0001)||{110}, <11-20>||<001>; 622-parent/432-child guards) +
  `standard_hcp_bcc_relationships` catalog (PS + inverse Burgers). Pinned: exact basal
  defining parallelism; **5.26 deg separation from inverse Burgers** (literature value, the
  hexagonal analogue of KS-Pitsch); 3 variants (internally derived orbit count, noted as
  such in the ledger). Remaining F14: Potter, Bagaryatsky/Isaichev. 872 passed; gates green.

## Phase 24 outcomes (2026-07-17)

- **F7 first stage (experimental):** `pytex.experimental.identify_orientation_relationship`
  + `ORIdentificationReport`. Key derivation: same-parent boundary misorientations populate
  the double coset `G_c (R S_p R^T) G_c` (parent symmetry conjugated by the OR rotation,
  child symmetry both sides); the fingerprint set is generated once per candidate
  (quaternion-key deduped), and per-edge distances need only elementwise traces
  (`einsum("eij,kij->ek")`) — no matrix products at scoring time.
- Verified: KS- and GT-generated microstructures identify their generating OR at 0.000 deg
  mean with >1 deg (KS: 3.67 deg) margins over all other fcc-bcc candidates; 0.3-deg noise
  preserves the ranking. describe() states the honest limit (identification only; rotation
  refinement from boundaries not implemented).
- Docs synced: foundation §5 (F7 first stage), validation row (foundational), CHANGELOG.
  871 passed; all gates green.

## Phase 23 outcomes (2026-07-17)

- **F14 slice:** `from_shoji_nishiyama_correspondence` (fcc->hcp epsilon-martensite;
  {111}||{0001}, <-110>||<11-20>; cubic/hex guards) + `standard_fcc_hcp_relationships`
  catalog. Pinned: 4 variants (one per {111} parent plane, the literature count), exact
  defining-parallelism mapping, one-variant-per-packet structure. Remaining F14: Pitsch-
  Schrader, Potter, Bagaryatsky/Isaichev. Docs synced (validation row, foundation, concept
  page, spec, CHANGELOG). 867 passed; all gates green.

## Phase 22 outcomes (2026-07-17)

- **F12 landed:** `OrientationRelationship.deformation_gradient()` +
  `DeformationGradientReport` — nearest-integer lattice correspondence (rint of the exact
  index correspondence; ray rationalization is WRONG here because magnitudes carry the
  strain), parent-frame gradient `F = R^T A_c M_int A_p^-1`, right-stretch via eigh, polar
  decomposition.
- Pinned physics: Bain principal stretches (sqrt(2) a_c/a_p, x2, a_c/a_p) and volume ratio
  2(a_c/a_p)^3 exact with zero polar rotation; KS/NW/GT all share the identical Bain
  stretches with polar rotations 11.06 / 9.74 / 10.15 deg — the literature rigid-body
  rotations relative to Bain, falling out of the decomposition.
- Docs synced: OR foundation §1/§5 (F12 moved to implemented; remaining gaps now
  F7/F13/F14 + reconstruction stabilization), concept page doctrine item 3, validation-matrix
  row, spec list, CHANGELOG. Gates: 866 passed, zero warnings; ruff/mypy/integrity/Sphinx
  green.

## Phase 21 outcomes (2026-07-17)

- New worked-example gallery group `texture` with
  `texture-gaussian-kernel-normalization-and-halfwidth` (A_0 = 1 and
  psi(halfwidth) = psi(0)/2 computed live — analytic identities), closing the quality-bar
  worked-example obligation for the new kernel surface. Gallery regenerated; 864 passed;
  all gates green.

## Session hand-off summary (2026-07-17, end of session)

Phases 16-21 this session: hexagonal property suites, transformation performance benchmark
lane (finding 21), CHANGELOG + release policy (finding 22), Gaussian/Abel-Poisson kernels
(finding 8 first slice), foundational-docs truth sync, kernel worked example. Remaining
program, in order: finding 8 remainder (bump/fibre kernels, SO3FunHarmonic), finding 9
(PF->ODF zero-range/ghost correction), finding 14 (Kikuchi geometry), F7/F12/F13/F14 of the
OR foundation, reconstruction stabilization (measured-data fixture + MTEX parity), findings
5, 10-13, 15 per the development guide. Resume from this note.

## Phase 20 outcomes (2026-07-17) — docs truth sync

- Stale foundational claims corrected per the stale-claims-are-defects rule:
  OR foundation doc §1 rewritten to the post-Cycles-A/B implemented surface and §5 limits
  updated (only F7/F12/F13/F14 + reconstruction stabilization remain);
  phase-transformation foundation "Current Limits" now defers to the OR foundation;
  repo-audit scorecard row for phase transformation moved to strong/foundationally-ready
  (dated); specifications executive summary updated to the current posture.
- All gates green (863 passed; integrity; Sphinx clean).

## Phase 19 outcomes (2026-07-17)

- **Finding 8 (kernel breadth) first slice landed:** `GaussianSO3Kernel` (spectral
  Gauss-Weierstrass, `A_l = (2l+1) exp(-l(l+1) eps)`) and `AbelPoissonKernel`
  (`A_l = (2l+1) kappa^(2l)`) in `texture/kernels.py`; halfwidth solved by bisection on
  psi(halfwidth)=psi(0)/2; series evaluation via SO(3) characters with the omega->0 limit;
  `KernelSpec` accepts names "gaussian"/"abel_poisson" and `as_so3_kernel()` routes all
  three. Tests: normalization (A_0=1 exact), halfwidth property, quadrature round trip of
  the spectrum, bandwidth monotonicity, Gaussian-vs-Abel tail comparison, spec routing.
  Remaining finding-8 scope: bump/fibre kernels + SO3FunHarmonic integration (with
  finding 9). Parity-ledger and CHANGELOG entries updated. 863 passed; all gates green.

## Phase 18 outcomes (2026-07-17)

- **Finding 22 landed:** root `CHANGELOG.md` (Keep-a-Changelog; Unreleased section captures
  the whole recent program at feature level, with the convention bug honestly under Fixed);
  "Release And Changelog Policy" section added to the API-stability standard (version bump +
  changelog cut + tag + green CI; scientific behavior changes must be categorized honestly);
  changelog linked from the docs index. 855 passed; all gates green.

## Phase 17 outcomes (2026-07-17)

- **Finding 21 (benchmark lane) landed for the transformation stack:**
  `scripts/benchmark_transformation_performance.py` — pinned cases (seed 20260717): KS
  intervariant table (276 pairs), or_deviation + fitting on 5000 pairs, reconstruction on
  400 grains; `--quick` smoke sizes covered by `tests/unit/test_benchmark_lane.py`. Reference
  timings this machine: deviation 1.41 s / fitting 0.18 s (5000 pairs), reconstruction 3.4 s
  (400 grains). Results JSON is local/git-ignored; `benchmarks/` stays reserved for
  schema-validated manifests (the manifest gate globs every JSON there — learned the hard
  way). No CI timing gate by design (noisy runners); the evidence is the runnable lane.

## Phase 16 outcomes (2026-07-17)

- Hexagonal/Burgers property suites (finding 20 breadth): inverse-transpose and zone-law
  invariants hold per variant with the hexagonal child metric; 6/mmm symmetry orbits recover
  primitive antipodal-canonical integer indices for arbitrary (hkl); beta->alpha->beta plane
  round trips invert exactly on the exact components and recover primitive indices when the
  forward image is rational. 854 passed, zero warnings.

1. Queued ledger follow-ups: MTEX `calcParent2Child` parity fixture (needs MATLAB-generated
   fixture data — see scripts/mtex_generators); external measured-data reconstruction
   fixture.
2. Larger Cycle C programs remaining (per the roadmap sequencing): texture kernel breadth
   (finding 8), PF→ODF ghost correction (finding 9), Kikuchi geometry (finding 14). Findings
   21 (benchmark lane) and 22 (release engineering) are now closed.
