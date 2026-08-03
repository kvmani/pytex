# Working Notes: Transformation Crystallography And Diffraction Program (TX)

Running ledger for the **TX** program. The normative specification is
[Transformation Crystallography And Composite Diffraction Program](../architecture/transformation_crystallography_and_diffraction_program.md).

This file exists so an interrupted session can resume without reconstructing context from chat
history (AGENTS.md, "Durable progress and resumability"). Keep it current *before* every long
verification run and every commit.

## Objective

Deliver five user-facing answers on top of existing PyTex primitives:

- **(a) TX1** — measured parent/child Euler angles ⇒ *what is the OR?*
- **(b) TX2** — an OR + an arbitrary parent plane/direction ⇒ the parallel planes/directions in
  every product variant.
- **(c) TX3** — a parent zone axis ⇒ a robust composite kinematic SAED, exportable as graphics
  **and** reflection tables.
- **(d) TX4** — a *product-variant* zone axis ⇒ the same composite, with matrix and sibling
  variants generated around it.
- **(e) TX5** — a measured SAED pattern ⇒ solved, from spots picked interactively or listed in a
  YAML file.

Closed out by **TX6**: a Burgers β↔α notebook demonstrating (a)–(e) end to end.

## Ground rules for this program

- Everything on `main`; no feature branches (user instruction).
- Commit after each phase's gates pass, so no work is ever at risk.
- Reuse before invention: no TX surface may re-derive a rotation convention, a symmetry reduction,
  a rationalization, or a detector basis that already exists in the core.
- Every new report object gets `describe()` and a JSON contract in lockstep.
- Every new numerical surface gets an executable worked example with **independent** provenance.

## Phase status

| Phase | Scope | Status | Commit |
| --- | --- | --- | --- |
| TX0 | Specification + ledger | DONE | `8d75f4c5` |
| TX1 | OR characterization from measured orientations | DONE | `8d623366` |
| TX2 | Variant correspondence tables | DONE | `971c6431` |
| TX3 | Composite SAED robustness + export layer | DONE | `0bb2256f` |
| TX4 | Child-zone-anchored composite patterns | DONE | `d34d2f78` |
| TX5a | Measured-pattern YAML + calibration + solver core | DONE | `d3c03cfc` |
| TX5b | Variant assignment + interactive picker | DONE | `d3c03cfc` |
| TX6 | Burgers notebook + ledger closure | DONE | (this commit) |

## Baseline established at TX0 (verified against live code, 2026-08-03)

What already exists, so later phases do not rebuild it:

- `pytex.core.transformation`: `OrientationRelationship` with eleven correspondence constructors
  (Bain, NW, KS, GT, Pitsch, Shoji–Nishiyama, Burgers, Pitsch–Schrader, Potter, Bagaryatsky,
  Isaichev), `generate_variants()`, `map_{plane,direction}_to_{child,parent}`,
  `map_{plane,direction}_across_variants`, `find_parallel_{planes,directions}`,
  `fit_orientation_relationship`, `or_deviation`, `intervariant_boundary_fingerprint`,
  `boundary_fingerprint_distances_deg`, `deformation_gradient`, `variant_pole_figure`.
- `pytex.core.parent_reconstruction`: `OrientationRelationshipCatalog` and the
  `standard_{fcc_bcc,bcc_hcp,fcc_hcp,hcp_bcc,ferrite_cementite}_relationships` builders.
- `pytex.experimental`: `identify_orientation_relationship` (child–child boundaries),
  `refine_orientation_relationship_from_boundaries`, `reconstruct_parent_grains`.
- `pytex.diffraction.kinematic`: `simulate_zone_axis_spots`, `zone_basis_from_axis`, `SpotTable`,
  `KinematicSimulationConfig`, `centering_allowed_mask`, `electron_structure_factors`.
- `pytex.diffraction.composite`: `simulate_composite_saed`, `CompositeSAEDPattern`,
  `VariantZonePattern`, `rationalize_zone_axis`, `find_spot_coincidences`, `sweep_parent_zone_axes`.
- `pytex.diffraction.models`: `DiffractionGeometry`, `KinematicSimulation`, `index_saed_pattern`,
  `estimate_zone_axis` — a **detector-geometry-driven** indexing path, distinct from TX5's
  calibrated-spot-list path; TX5 must state the difference and reuse what it can.
- `pytex.plotting.composite_saed`: renderer, `CompositeSAEDPlotConfig`, annotation engine.

Identified gaps are enumerated per feature in the specification §3.2, §4.2, §5.2, §6, §7.

## Known pre-existing issues on this machine (inherited, not caused by TX)

Carried forward from the reconstruction-stabilization ledger; do **not** fix inside a TX
scientific-behavior commit:

1. ~~Six phase-fixture SHA-256 mismatches~~ — **FIXED in TX0.** The cause was confirmed, not
   guessed: `fixtures/phases/fe_bcc/phase.cif` held 19 CRLF pairs on disk and hashed to
   `e512334e…`, while the LF form hashes to `8afe4f95…` — exactly the digest pinned in
   `fixtures/phases/catalog.json`. Git's Windows default `core.autocrlf=true` was rewriting
   checksum-pinned artifacts on checkout. A `.gitattributes` marking
   `fixtures/phases/**`, `fixtures/mtex_parity/**` and `*.ipynb` as `-text` disables the
   conversion; `scripts/check_repo_integrity.py` now passes. **The full test suite no longer
   needs the two deselects.**
2. ~20 ruff findings from newer rule versions (RUF022/RUF059/RUF043) in untouched files.
3. Two mypy `to_hex` arg-type errors in `plotting/crystal3d.py` (matplotlib stub drift).
4. No MATLAB/MTEX on this machine, so no new MTEX parity claim can be executed here.

## Verification command set

```
python -m pytest
python scripts/check_repo_integrity.py
python -m ruff check .
python -m mypy src
python -m sphinx -b html docs/site docs/_build/html
python scripts/generate_worked_examples.py
```

## Ledger

### TX0 (2026-08-03) — specification and ledger

- Read the governing documents and audited the live code surface for all five asks; the baseline
  above is what that audit found, not an assumption.
- Wrote the normative specification with per-feature API signatures, algorithms, and validation
  plans, including the honest-limits requirements (`is_conclusive` semantics for TX1 and TX5, the
  `z` vs `-z` SAED ambiguity, kinematic-only intensities).
- Recorded the phase order and its dependencies: TX4 needs TX3's geometry, TX5b needs TX4.
- **Fixed the inherited fixture-hash gate failure** (see item 1 above) by adding `.gitattributes`.
  This was worth doing first: the integrity check is a gate for every TX phase, and leaving it red
  would mean every later phase runs with two deselects and cannot tell a new breakage from the old
  one.
- Artifact hygiene: `.gitignore` now covers `docs/site/_build/`, the regenerated
  `fixtures/mtex_parity/results/pytex/*/` directories, and the stray root `package.json`, all of
  which were sitting untracked in the worktree.

### TX1 (2026-08-03) — OR determination from measured orientations

Shipped, gates green (full suite, integrity, ruff, mypy, Sphinx zero-warning).

- **New core surface**: `characterize_orientation_relationship`,
  `orientation_relationship_from_euler`, `describe_orientation_relationship`,
  `ORCharacterizationReport`, `ORParallelismStatement`, `DEFAULT_OR_TOLERANCE_DEG`, and
  `default_relationship_catalog` (in `parent_reconstruction`, next to the other catalog
  builders). All exported from `pytex.core` and the top level.
- **Refactor, not duplication**: `fit_orientation_relationship`'s align/average loop is now
  `_fit_from_seed`, shared by both entry points, plus `_symmetry_operator_pair` and
  `_measured_parent_to_child` so the `V = C^T P` convention has exactly one definition.
  The existing fitting tests pass unchanged, which is the evidence the refactor is behavior
  preserving.
- **A real defect found and fixed while building the seedless start.** The first design
  reduced *every* pair to its minimum-angle double-coset representative and averaged them.
  That is wrong: the maximum-trace element is not unique when the relationship's own
  rotation is symmetric, so different pairs land on different tied representatives.
  Measured: planted Bain (45 deg / <100>, three variants) averaged to **26.9 deg** and was
  reported as Kurdjumov-Sachs. Fix — reduce **one** pair and let the alignment step resolve
  the rest against it. `test_bain_survives_the_double_coset_tie` fails if this regresses.
- **Statement extraction needed a preference, and the reason is scientific.** A rotation
  satisfies several exact low-index parallelisms simultaneously; for KS both
  `(111)||(011)` and `(10-1)||(11-1)` are exact, and index magnitude alone tie-breaks
  arbitrarily (it picked the latter). Which one the literature quotes depends on the two
  structures' close-packed planes and directions, which a rotation does not know. So the
  search takes a preference: by default the relationship's own recorded defining families,
  and for a fitted relationship those of the matched catalog member. Fit quality outranks
  preference in the sort, so a nominated family can never promote a visibly worse clause
  above an exact one.
- **`is_hexagonal_phase` moved to `pytex.core.hexagonal`** (one-shared-helper rule) so core
  notation can label hexagonal statements in four-index Miller-Bravais form;
  `pytex.diffraction.composite` re-exports it and its tests pass unchanged.
- **26 tests** in `tests/unit/test_or_characterization.py`, every expected value from a
  definitional parallelism, a published separation (KS-NW 5.26 deg, KS-GT 2.40 deg), or an
  analytic identity. **Two worked examples** in the transformation gallery, regenerated.
- **Measured noise envelope** (planted KS, 12 pairs across mixed variants): conclusive at
  0.5 and 2.0 deg scatter; correctly *inconclusive* at 5.0 deg. The failure mode is an
  admitted "cannot separate the candidates", not a confident wrong name.
- **Docs synced**: CHANGELOG (Added + Fixed), concept page (new section with the seeding
  subtlety and the preference rationale), OR foundation §1/§3/§5 (new F6b, honest limits),
  phase-transformation foundation, symbol registry (two new entries), validation matrix
  (five new rows), and a site stub + toctree entry for the program spec.

### Deferred from TX1

- `ORCharacterizationReport` has `to_json_dict()` (a one-way report payload) but is **not**
  registered in `pytex.contracts`' round-trip serializer registry. Report objects generally
  are not; revisit as a batch when TX3's manifests land.
- The cubic-cubic catalog assumes an fcc->bcc transformation because point-group symmetry
  cannot distinguish fcc from bcc. Documented in `default_relationship_catalog` and in the
  foundation limits; a structure-aware dispatch (using the space group / centring) would
  remove the assumption.

### TX2 (2026-08-04) — variant correspondence tables

Shipped, gates green.

- **New core surface**: `variant_correspondence_table`, `VariantCorrespondenceTable`,
  `VariantCorrespondenceRow`. Takes one object or a list (homogeneous kind, one phase),
  either `sense`, an optional variant subset and a rationalization bound; returns rows with
  the exact image, integer indices, residual, labels, and an equivalence-group id.
- **Delegates, does not reimplement**: every row comes from `map_plane_to_child` /
  `map_direction_to_child` and their parent-inverses, so the index-map semantics have one
  definition. The table adds grouping, labels, `describe()` and the exports.
- **The grouping is the value.** KS `(111)` -> 4 distinct images across 24 variants, 6 each;
  the 6 exact ones are `{011}` at zero residual. The test asserts those 6 *are* the packet
  `variant_close_packed_groups` returns, comparing two independent computations rather than
  a stored constant.
- **Asymmetry worth knowing** (now documented and tested): the reverse map is not selective
  — child `(011)` maps back onto `{111}` in **all 24** variants, one equivalence group,
  because each variant's close-packed image came from some `{111}` member.
- **Rationalization behavior pinned**: raising `max_index` never worsens a residual, and the
  set of exactly-parallel variants is identical at bounds 3 and 17. Only the labeling of the
  *irrational* images depends on the bound; `describe()` states this, because a reader
  otherwise reads "6 distinct images" as physics when it is partly a bookkeeping choice.
- **Rationalization is sign-sensitive** (pre-existing, by design — it matches the exact
  image's direction, which matters for a diffraction vector `g`), so the Burgers basal image
  appears as `(0001)` or `(000-1)` depending on the variant. The test accepts both and says
  why.
- **24 tests** in `tests/unit/test_variant_correspondence_table.py`, **one worked example**
  (the packet identity, expected `[24, 4, 6, 6, 6]` with tolerance 0).
- **Docs synced**: CHANGELOG, concept page section, OR foundation §1 (new F2b entry),
  validation matrix (three new rows).

### TX3 (2026-08-04) — composite SAED export, and a real absence bug

Shipped, gates green.

- **The spec was wrong about three of its five "gaps", and has been corrected in place.**
  Auditing the live engine before writing code showed that deterministic spot ordering
  (`np.lexsort((l, k, h, radius, -intensity))`), config/basis guard rails, and the intensity
  normalization contract were all already implemented and documented. A "shared" cross-phase
  normalization option was **rejected**, not deferred: kinematic theory defines no intensity
  ratio between two phases, so the option would manufacture a number the theory does not
  support. Spec §5.2 now records the correction and the reasoning.
- **The one genuine robustness gap was real and had already bitten the repository.**
  `ReflectionCondition.from_phase` falls back to primitive when a phase declares no space
  group, so an undeclared body-centred phase is simulated as primitive and lists forbidden
  reflections silently. The shared Burgers *worked-example* setup had exactly this defect —
  it was simulating beta-titanium without body-centring absences. Both phases now declare
  their space groups, and a worked example pins that no `h + k + l` odd beta reflection
  survives. New `phase_centering_is_declared(phase)` and
  `CompositeSAEDPattern.centering_audit()`; `describe()`, the reflection table and the
  manifest all state the centering applied and whether it was declared or assumed.
- **New module `pytex.diffraction.export`**: `ReflectionTable` / `ReflectionTableRow` /
  `composite_reflection_table`, `composite_saed_manifest`, `CompositeSAEDExport` /
  `export_composite_saed`, plus the public column contract `REFLECTION_TABLE_COLUMNS`.
  Every table value is read from the engine's `SpotTable`, so table and figure cannot drift.
- **New schema** `schemas/composite_saed_manifest.schema.json` (+ README entry and
  `composite_saed_manifest_schema_path()` in the manifests adapter). The manifest is
  validated against it in tests.
- **A test caught a documentation error I introduced.** The first draft of
  `ReflectionTable.describe()` claimed "detector radius = camera constant x |g|". It is the
  camera constant times the *in-plane* part of `g`; the difference is the out-of-plane
  component the excitation error records. Both the prose and the test now say so, and the
  test asserts the projected identity to 1e-12 plus the inequality against the full `|g|`.
- **23 tests** in `tests/unit/test_composite_saed_export.py` (geometry identities, Friedel
  symmetry, thresholding, the centering audit including the undeclared-phase case, file
  export, no leaked figures, CSV/Markdown/JSON contracts, schema validation), **one worked
  example** pinning four identities at once.
- **Docs synced**: CHANGELOG, workflow page (new "Exporting the pattern" section with the
  centering trap), diffraction validation matrix (three new rows), schema README.

### TX4 (2026-08-04) — child-zone anchoring, and a sort-stability defect it exposed

Shipped, gates green.

- **`simulate_composite_saed_from_child_zone(relationship, child_zone_axis,
  anchor_variant_index=k, ...)`** maps the requested child zone back through `R_k^T` and
  delegates to `simulate_composite_saed`, so there is exactly one detector-geometry
  definition. `align_child_g` is given in the child's own indices and mapped internally.
- **`simulate_composite_saed` generalized** to accept an irrational `CrystalDirection` as
  well as a `ZoneAxis`, plus an `align_g_cartesian` escape hatch and an
  `anchor_variant_index` bookkeeping field. `CompositeSAEDPattern` gained
  `anchor_variant_index`, `nearest_parent_zone_axis`, `parent_zone_axis_label()` and
  `anchor_description()`; every label site (coincidence report, reflection table, plotting
  title and legend) now routes through the accessor instead of reading `.indices` directly.
- **The consistency identity found a real defect.** Anchoring on variant `k`'s image of a
  parent zone should reproduce the parent-anchored pattern exactly. The first run showed
  spots displaced by up to **486 mm** — with matching counts. The positions were right; the
  *row order* was not. The sort keys are intensity then radius then `hkl`, and
  symmetry-equivalent reflections have mathematically equal intensity and radius that differ
  by ~1e-14 depending on how the basis was built, so noise decided the order before the exact
  `hkl` tie-break was reached. Both continuous keys are now quantized before `lexsort`
  (1 pm of radius, 1e-12 of full-scale intensity). The identity then holds to **1e-13 mm**
  with identical `hkl` ordering. This is a behavior change for exported table row order and
  is recorded under `Fixed` in the CHANGELOG.
- **18 tests** in `tests/unit/test_composite_saed_child_anchor.py`: the identity for four
  anchor variants, sort-order stability under a 1e-15 perturbation and across repeated runs,
  the Burgers basal-zone geometry, `align_child_g` placing a child reflection on `+u`, the
  anchor recorded in `describe()`/manifest, and validation. **One worked example** pinning
  the identity at tolerance 1e-9 mm.
- **Schema updated**: `parent_zone_axis` is now `number` rather than `integer` (it is the
  exact parent direction, generally irrational), plus required `anchor_variant_index` and
  `parent_zone_axis_nearest`.
- **Docs synced**: CHANGELOG (Added + Fixed), workflow page (new "Anchoring on a product
  zone" section including the spot-order note), diffraction validation matrix (two rows).

### TX5 (2026-08-04) — solving a measured pattern (TX5a and TX5b together)

Shipped, gates green. TX5a and TX5b landed in one commit because the variant
assignment and the picker are both thin layers over the solver core, and splitting them
would have meant a commit whose tests could not exercise the file contract end to end.

- **New module `pytex.diffraction.solving`**: `PatternCalibration`, `MeasuredSpot`,
  `MeasuredSAEDPattern` (YAML in/out), `SolvedSpot`, `PatternSolution`,
  `PatternSolutionReport`, `solve_saed_pattern`, `solve_saed_pattern_file`,
  `assign_transformation_variant`. **New module `pytex.plotting.saed_picker`**:
  `SpotPickerState` (the logic) and `SAEDSpotPicker` (the Matplotlib adapter).
- **New schema** `schemas/measured_saed_pattern.schema.json` +
  `measured_saed_pattern_schema_path()`.
- **Two design decisions worth recording.**
  - *Intensities are never used for indexing.* A kinematic intensity model is not
    reliable enough to index against and a printed pattern rarely carries calibrated
    intensities; geometry alone decides. Intensity is carried through for plotting only.
  - *The picking logic is separated from the GUI.* `SpotPickerState` has no Matplotlib
    dependency and is fully tested headlessly. An interactive tool that cannot be tested
    is a liability, and the tests would otherwise have been skipped everywhere.
- **Two problems found by the first end-to-end run, both fixed.**
  - Symmetry-equivalent descriptions were being reported as *competing* solutions, so an
    unambiguous cubic solve looked contested (`is_conclusive` False with five
    100%-matched entries). Deduplication now compares orientations under the crystal
    point group, and the survivor is rewritten into the conventional description —
    fewest negative indices, then lowest — so `[001]` is reported rather than the
    equally valid `[0-10]` the seed search found first.
  - Spots whose indices exceed the solver's `max_index` are never offered a match, which
    reads as a solver failure. `describe()` now names `max_index` as the first thing to
    raise, before the tolerances.
- **A limitation found and pinned rather than hidden**: a variant seen from a *parent*
  zone axis is off its own zone (its child zone axis is irrational), so its
  excitation-selected spots do not all lie in one ZOLZ and cannot all be indexed. The
  test asserts the partial match — a full match there would mean the solver was
  inventing reflections. The realistic case (tilt the product on zone) uses TX4's
  child-anchored entry point and indexes fully.
- **41 tests** in `tests/unit/test_saed_solving.py`: calibration in all three unit
  systems, the YAML contract and schema, simulate-then-solve closure on five fcc zones
  plus hexagonal alpha-Ti, the cubic sqrt(2)/45-degree identities read back out of the
  solver's own output, fcc-vs-bcc discrimination, an unsolvable pattern returning
  nothing, noise robustness at 0.5 mm and degradation at 25 mm, variant assignment, and
  the picker state machine. **One worked example** pinning the round-trip closure.
- **Docs synced**: new workflow page `saed_pattern_solving.md` (+ toctree), CHANGELOG,
  diffraction validation matrix (five rows), schema README.

### TX6 (2026-08-04) — the end-to-end notebook, and program closure

- **`docs/site/tutorials/notebooks/23_transformation_crystallography_end_to_end.ipynb`**,
  committed executed (21 code cells, 2 figures, zero error outputs), registered in the
  notebook toctree and in the teaching-track summary. It walks Burgers beta->alpha
  through all five asks in one pass, with every printed number computed live.
- **The notebook is also an integration test that a unit test could not be.** It found no
  defects, which is itself the result worth recording: the five surfaces compose without
  glue code. The consistency identity reads 7e-14 mm, the planted variant 3 is assigned
  back at 0.0000 deg, the YAML round trip is exact to 0.0e+00 1/angstrom, and the
  wrong-phase solve returns *no solution at all* rather than a plausible-looking one.
- Each section states the crystallography before computing it, and the notebook closes
  with what was **not** shown — kinematic only, no HOLZ or double diffraction, the
  zone-axis assumption in the solver, synthetic validation for the OR determination, and
  no MTEX parity claim anywhere.

## Program outcome

All five asks are delivered, each with tests, an executable worked example, concept or
workflow documentation, a validation-matrix row, and a CHANGELOG entry:

| Ask | Surface | Tests | Worked examples |
| --- | --- | --- | --- |
| (a) OR from measured Euler angles | `characterize_orientation_relationship`, `orientation_relationship_from_euler`, `describe_orientation_relationship` | 26 | 2 |
| (b) Parallel planes/directions across variants | `variant_correspondence_table` | 24 | 1 |
| (c) Composite SAED + exports | `composite_reflection_table`, `export_composite_saed`, `centering_audit` | 23 | 1 |
| (d) Child-zone-anchored composites | `simulate_composite_saed_from_child_zone` | 18 | 1 |
| (e) Solving a measured pattern | `solve_saed_pattern`, `MeasuredSAEDPattern`, `SAEDSpotPicker`, `assign_transformation_variant` | 41 | 1 |

**Four real defects were found and fixed along the way**, three of them pre-existing:

1. Checksum-pinned fixtures failed on every Windows clone (`core.autocrlf` rewriting
   hash-pinned artifacts). Fixed in TX0; the suite now runs green with no deselects.
2. The seedless OR fit averaged tied double-coset representatives, turning planted Bain
   into a meaningless 26.9 deg that read as Kurdjumov-Sachs. Found in TX1.
3. The shared Burgers worked-example setup declared no space groups, so it had been
   simulating beta-titanium without body-centring absences and listing forbidden
   reflections. Found in TX3 by the centring audit built in the same phase.
4. Kinematic spot ordering was decided by floating-point noise at symmetry-equivalent
   ties, so the same pattern reached two ways came out permuted. Found in TX4 by the
   consistency identity.

**The specification was also corrected twice** where it had claimed gaps that did not
exist (spec §5.2) — auditing the live code before writing it turned out to matter more
than the plan did.

### Open follow-ons (not blockers, deliberately out of scope)

- **Measured-EBSD fixtures** for the OR determination. Validation is synthetic and every
  document says so.
- **JSON round-trip contracts** for the new report objects. They have one-way
  `to_json_dict()` payloads; the `pytex.contracts` registry covers reconstructible
  objects only, and report objects generally are not registered there.
- **Canonical SVG figures** for the OR-statement geometry and the solving flow, and
  **LaTeX theory notes** for the parallelism extraction and the ratio/angle algorithm.
  The prose and the executable examples carry the content today; the figures and notes
  would make it publication-facing.
- **Structure-aware catalog dispatch.** Cubic-to-cubic assumes the fcc->bcc class because
  point-group symmetry cannot distinguish fcc from bcc; using the space group would
  remove the assumption.
- **HOLZ, double diffraction, and dynamical intensities** remain out of scope for both
  the simulator and the solver, as stated from TX0.
