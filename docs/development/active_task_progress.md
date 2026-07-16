# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Objective

Execute **Cycle A of the Critical Review And Development Guide**
(`docs/roadmap/critical_review_and_development_guide.md`), one phase per commit+push, as a
long-horizon mission. Findings addressed: 1, 2, 6, 16, 17, 18 (23/24 already fixed with the
guide itself). The OR feature definitions (F1–F5) live in
`docs/architecture/orientation_relationship_analysis_foundation.md`.

## Current Status

- Started: 2026-07-16
- Branch: `main` (tracking `origin/main`, push after each phase)
- Baseline commit: `d1a1561`; Phase 0 pushed as `84a9767`
- Phase: 1 complete (committing), next Phase 2

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

## Next Actions

1. Commit Phase 1, push.
2. Phase 2 (F1–F3): `correspondence_direct()`/`correspondence_reciprocal()` on
   `OrientationRelationship` + `TransformationVariant`; `map_direction_to_child` /
   `map_plane_to_child` + child→parent inverses returning exact components, rationalized
   indices (bounded search), angular residual; all-variant table. Worked example (Bain/KS),
   LaTeX theory note `docs/tex/algorithms/orientation_relationship_correspondence.tex`,
   validation-matrix row, parity note. Key math: directions map via
   `A_c^{-1} R A_p` (direct bases), plane normals via reciprocal bases (inverse-transpose).
