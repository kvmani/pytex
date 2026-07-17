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
- Phase 7 pushed as `d40fb14`
- Phase: 8 complete (committing) — Cycle B in progress

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

## Next Actions (Cycle B remainder, per the development guide §3)

1. Map-scale parent-grain reconstruction v1 (F8, finding 4) — variant-graph voting on the
   grain-boundary network (`ebsd/models.py` grains + union-find as in `merge_by_csl`);
   synthetic fixture first (planted parent grains -> variants -> reconstruct), literature
   fixture second. Stage under `pytex.experimental` until validation breadth exists.
   Building blocks now all exist: variants, or_deviation, fit_orientation_relationship,
   vectorized select_variants. This is the last major Cycle B item — an L-size,
   multi-session program; open a fresh session with this note as the entry point.
2. Queued follow-ups recorded in ledgers: MTEX calcParent2Child parity fixture; OR-fitting
   worked example; OrientationSet slicing API (returns malformed Orientation today);
   broaden property suites to hexagonal phases and OrientationSet reductions.
