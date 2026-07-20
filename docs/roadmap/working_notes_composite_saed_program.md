# Working Notes: Composite OR Diffraction Pattern Program (2026-07-20)

Purpose: resumable progress ledger for the composite TEM zone-axis diffraction
pattern program. The goal is a sophisticated, fully vectorized, highly
configurable engine that simulates and renders **composite selected-area
electron diffraction (SAED) patterns** — parent phase plus any subset of
orientation-relationship (OR) variants — for **any** parent zone axis, with
publication-grade plotting defaults. Scope is **kinematic only**: no dynamical
(multi-beam/Bloch-wave) effects in this program. This supersedes and greatly
extends the naive composite-pattern capability of the author's earlier
`pycrystallography` package (configurable spot size/shape/symbols, variant
subsets, in-plane rotation, and non-overlapping spot annotation were its
signature features; all are re-designed here on the PyTex canonical model).

If interrupted, resume by reading this file top to bottom; each phase lists
status, deliverables, and exactly what remains. Master context:
`docs/development/active_task_progress.md`.

## Scientific conventions (pinned for the whole program)

- Orientations are crystal→specimen rotations; variant rotations `V_i` map
  parent-crystal-frame Cartesian vectors to child-crystal-frame Cartesian
  vectors of the same physical direction
  (`TransformationVariant.map_parent_vector_to_child`).
- Beam travels antiparallel to the zone axis; the zone axis unit vector
  `z_p` (parent crystal frame, Cartesian) points toward the electron gun.
- The same physical beam direction in child variant *i*'s crystal frame is
  `z_c = V_i z_p`. Child zone axes are in general **irrational**; exact
  Cartesian vectors drive the geometry, and nearest-rational `[uvw]` indices
  (with angular deviation reported) are used only for labeling.
- Shared detector basis: an orthonormal in-plane pair `(u, v)` fixed in the
  **parent** crystal frame (deterministic construction + optional
  `g`-alignment + in-plane rotation). Child reciprocal vectors are pulled back
  to the parent frame (`g_p = V_i^T g_c`) before projection, so all phases and
  variants land on one physically consistent detector.
- Reflection selection is by **excitation error**: `|s_g| <= s_max` with
  `s_g = g_z - g^2/(2k)` for beam wavevector magnitude `k = 1/λ` and `g_z`
  the zone-axis component of `g` (small-angle, kinematic). This treats
  rational parent zones and irrational child zones uniformly and is honest
  about Ewald-sphere curvature. Relativistic electron wavelength
  `λ(V) = h / sqrt(2 m₀ e V (1 + eV/(2 m₀ c²)))`.
- Kinematic intensity: `I ∝ |F_hkl|²` from the electron structure-factor
  proxy already used in `pytex.diffraction.saed` (atomic-number scattering,
  Debye-Waller from `b_iso`), with optional relrod damping
  `1/(1 + (s_g/σ_s)²)` and lattice-centering systematic absences
  (`ReflectionCondition`). Intensities normalized to max = 1 per pattern.
- Detector coordinates: `r_mm = L λ g_⊥ = (camera constant) · g_⊥`, the
  standard SAED small-angle approximation already used by
  `generate_saed_pattern`.

## Program Design (phases; each = verified commit, pushed when remote allows)

| Phase | Deliverable | Status | Commit |
| --- | --- | --- | --- |
| CD0 | This working-notes file + ledger update | complete | f609c0a |
| CD1 | `diffraction/kinematic.py`: vectorized zone-axis engine + `SpotTable` | complete | e549bce |
| CD2 | `diffraction/composite.py`: composite OR pattern assembly | complete | d292f1b |
| CD3 | `plotting/composite_saed.py`: config model + layered renderer | complete | 202b6da |
| CD4 | Annotation engine with coincident-label merging + crowding avoidance | complete | 88c202e |
| CD5 | Spot-coincidence analysis report + zone-axis sweep utilities | complete | 1a93130 |
| CD6 | Worked examples, exports, docs index, CHANGELOG, final verification | complete | (this commit) |

### CD1 — Vectorized kinematic zone-axis engine (`src/pytex/diffraction/kinematic.py`)

- `electron_wavelength_angstrom(beam_energy_kev)` — relativistic, pinned test
  vs published 200 kV value (≈0.02508 Å).
- `KinematicSimulationConfig` (frozen): `beam_energy_kev`,
  `camera_constant_mm_angstrom`, `max_index`, `g_max_inv_angstrom`,
  `max_excitation_error_inv_angstrom`, `intensity_model`
  (`"electron_atomic_number" | "unit"`), `relrod_sigma_inv_angstrom | None`,
  `apply_centering_absences`, `min_relative_intensity`.
- `SpotTable` (frozen, struct-of-arrays, read-only ndarrays): `hkl (N,3)`,
  `g_crystal (N,3)`, `detector_mm (N,2)`, `g_detector_inv_angstrom (N,2)`,
  `d_spacing_angstrom (N,)`, `intensity (N,)` (max-normalized),
  `structure_factor_amplitude (N,)`, `excitation_error_inv_angstrom (N,)`;
  sorted by (-intensity, radius); `describe()`.
- `simulate_zone_axis_spots(phase, zone_axis_cartesian, *, config, basis)` —
  the core vectorized routine: full hkl-cube enumeration, vectorized centering
  mask, vectorized structure factors (broadcast over hkl × sites), excitation
  filter, projection to a supplied 3×3 zone basis. **No Python loop over
  reflections** (site/species loops allowed: few elements).
- `zone_basis_from_axis(zone_cartesian, *, align_g_cartesian=None,
  in_plane_rotation_deg=0.0)` — deterministic detector basis; optional
  alignment of a chosen reciprocal vector along +u; right-handed, pinned.
- Tests (`tests/unit/test_kinematic_engine.py`): wavelength pinned values;
  parity with legacy `generate_saed_pattern` for Ni [011] (same hkl set, same
  detector geometry within tolerance); FCC/BCC forbidden reflections absent;
  d-spacing pinned (Ni 111 ≈ 2.0345 Å); excitation errors satisfy
  `s = -g²/(2k)` for exact ZOLZ; basis orthonormality/right-handedness;
  in-plane rotation equivariance; determinism.

### CD2 — Composite assembly (`src/pytex/diffraction/composite.py`)

- `VariantZonePattern` (frozen): `variant` (TransformationVariant),
  `zone_axis_child_cartesian`, `nearest_zone_axis` (rationalized `[uvw]` +
  `deviation_deg`), `spots: SpotTable` (in shared detector frame).
- `CompositeSAEDPattern` (frozen): `relationship`, `parent_zone_axis`,
  `parent_spots: SpotTable | None`, `variant_patterns: tuple[...]`,
  `zone_basis_parent (3,3)`, `config`, `provenance`; helpers
  `variant_indices`, `select_variants(...)`, `all_detector_coordinates()`;
  `describe()` with convention-explicit prose (frame conventions, beam sense,
  selection rule, per-variant nearest zone axes).
- `simulate_composite_saed(relationship, parent_zone_axis, *,
  variant_indices=None, include_parent=True, config=None,
  align_parent_g=None, in_plane_rotation_deg=0.0, child_config=None)`.
- Rationalization of irrational child zones via bounded integer search
  (adapt `_rationalize_components` logic from `core/transformation.py`).
- Tests (`tests/unit/test_composite_saed.py`): KS fcc→bcc with parent
  [0 1 -1]: a variant whose child zone is exactly [1 1 -1]_bcc exists
  (KS parallelism `<-101>_p || <-1-11>_c` up to sign/orbit); NW parent [011]
  → child [001] pinned; Bain [001]_p → [001]_c with 45° in-plane relation
  pinned via spot coordinates; variant subsetting; shared-basis invariant;
  describe() content.

### CD3 — Plot configuration + renderer (`src/pytex/plotting/composite_saed.py`)

- `SpotStyle` (frozen): `marker`, `color`, `filled`, `size_scale`,
  `size_mode` (`"intensity_area" | "intensity_radius" | "constant"`),
  `min_size_pt`, `alpha`, `edge_color`, `edge_width`, `zorder`.
- `CompositeSAEDPlotConfig`: `parent_style`, `variant_styles`
  (explicit dict or palette cycling; colorblind-safe default palette),
  `variant_indices` subset, `show_transmitted_beam`, `axes_units`
  (`"mm" | "inv_angstrom"`), `show_legend`, `legend_labels`, `title`,
  `background`, `annotation: SpotAnnotationConfig` (CD4), theme integration
  with `plotting/styles.py`.
- `render_composite_saed(pattern, *, config=None, ax=None) -> Figure` —
  layered scatter (parent topmost by default), equal aspect, legend with
  phase + variant + nearest-zone labels, publication defaults.
- Structural tests (`tests/unit/test_composite_saed_plotting.py`): scatter
  collection counts/colors/marker sizes track config; variant subset honored;
  legend entries; axes units switch rescales coordinates; figures closed.

### CD4 — Annotation engine (in `plotting/composite_saed.py`)

- `SpotAnnotationConfig`: `enabled`, `max_labels`, `min_intensity`,
  `format` (`"plain" | "overline"` mathtext), `merge_coincident` (one label
  listing all coincident reflections, phase-tagged),
  `coincidence_tolerance_mm`, `offset_pt`, `font_size`, `leader_lines`,
  `avoid_overlap` (greedy candidate-offset placement with
  `scipy.spatial.cKDTree` collision checks), `label_color_follows_spot`.
- Deterministic placement; labels never overlap each other or spot markers
  beyond tolerance; skipped labels reported on a returned annotation result.
- Tests: coincident KS spots produce merged multi-phase labels; label boxes
  pairwise disjoint (structural assertion via matplotlib bbox extents);
  overline formatting pinned; determinism.

### CD5 — Coincidence analysis + utilities (`diffraction/composite.py` additions)

- `SpotCoincidenceReport` (frozen, `describe()`): pairs of
  (parent hkl, variant index, child hkl, separation_mm) within tolerance,
  grouped clusters, counts per variant — the quantitative statement of which
  reflections superimpose for a given OR/zone (key for OR verification in
  TEM practice).
- `find_spot_coincidences(pattern, *, tolerance_mm)` vectorized via cKDTree.
- `sweep_parent_zone_axes(relationship, zone_axes, ...)` convenience iterator.
- Tests: KS composite has pinned coincidence counts within tight tolerance;
  brute-force cross-check; describe() text.

### CD6 — Integration

- Exports: `pytex.diffraction.__init__`, top-level `pytex.__init__` lazy map,
  `docs/README.md` index entry, CHANGELOG entry, this ledger finalized.
- Worked examples (`worked_examples/examples/`): electron wavelength at
  200 kV (cited standard value); KS child-zone mapping angular identity.
- Regenerate gallery; full gates.

## Verification gates (every phase, before commit)

- `python -m pytest` (full suite green, no new warnings)
- `python -m ruff check .`
- `python -m mypy src`
- `python scripts/check_repo_integrity.py`

## Key facts established (verified against live code, 2026-07-20)

- `core/transformation.py`: `OrientationRelationship` (named constructors:
  Bain, NW, KS, GT, Pitsch, Burgers), `generate_variants` reproduces
  literature counts (Bain 3; NW/Pitsch/Burgers 12; KS/GT 24);
  `TransformationVariant.map_parent_vector_to_child` = child-frame
  re-expression (`R v`); `_rationalize_components` exists for nearest-integer
  index recovery (module-private; adapt, do not import privately across
  modules without promotion).
- `diffraction/saed.py`: loop-based single-phase kinematic SAED
  (`generate_saed_pattern`) with `_choose_zone_basis` (deterministic u,v from
  cross products), integer zone-law filter, electron structure-factor proxy
  `_structure_factor_electron` (Z-scattering + B_iso damping), intensity
  `|F|²/(1+g²)`. Kept intact for backward compatibility; CD1 engine is the
  new vectorized surface and should match it on rational exact zones
  (allowing for the deliberate intensity-model and selection-rule upgrades).
- `diffraction/physics.py`: `ReflectionCondition` centering absences
  (scalar `is_allowed`; CD1 adds a vectorized mask), `ScatteringFactorTable`.
- `plotting/styles.py` `resolve_style(theme=...)` provides themed defaults
  (`common`, `saed` sections); composite plotting should integrate but carry
  its own typed config layer.
- Repo doctrine: frozen slots dataclasses with `__post_init__` validation,
  read-only ndarrays, `describe()` on report objects, vectorized hot paths,
  no naked-array public APIs where frame meaning is ambiguous, tests with
  implementation, figures closed in tests, no byte-level SVG baselines for
  runtime plots.

## Completed

- (CD0, f609c0a) Program designed; conventions pinned above; ledger wired.
- (CD1) `src/pytex/diffraction/kinematic.py`: relativistic
  `electron_wavelength_angstrom` (pinned De Graef Table 2.2 values: 0.037014 /
  0.025079 / 0.019687 angstrom at 100/200/300 kV), `KinematicSimulationConfig`
  (validated frozen config incl. excitation-error half-width, optional relrod
  damping, centering-absence toggle), read-only struct-of-arrays `SpotTable`
  with `describe()`, `zone_basis_from_axis` (legacy-parity deterministic
  construction + g-alignment + CCW in-plane rotation, pinned), vectorized
  `centering_allowed_mask` + `electron_structure_factors` (broadcast over
  hkl; parity with the scalar legacy formula), and `simulate_zone_axis_spots`
  (no per-reflection Python loops; excitation-error selection handles
  irrational zones). 59 regression tests in
  `tests/unit/test_kinematic_engine.py`: wavelength pins, legacy-parity on
  Ni [011] (identical hkl set + detector coordinates to 1e-9 mm), FCC/BCC
  forbidden-reflection absences, centering-mask parity vs scalar reference on
  all seven centerings, Ni d(111)=2.03451 A pin, ZOLZ identity
  s_g = -g^2 lambda/2, HOLZ exclusion, sorting/normalization/determinism,
  shared-basis rotation round-trip, read-only arrays, config validation.
  Gates: 937 passed, ruff clean, mypy clean, integrity passed.

- (CD2) `src/pytex/diffraction/composite.py`: `rationalize_zone_axis` +
  `RationalizedZoneAxis` (bounded primitive-triple search, true angular
  deviation, sign-sensitive), `VariantZonePattern` (exact irrational child
  zone + rational label + shared-frame spots), `CompositeSAEDPattern`
  (variant lookup/subsetting, `iter_spot_tables`, stacked coordinates,
  convention-explicit `describe()`), `simulate_composite_saed` (parent-
  anchored shared basis via child-frame basis rotation `V_i B_p`, variant
  subsetting, separate child config, parent-g alignment, whole-composite
  in-plane rotation). 22 regression tests in
  `tests/unit/test_composite_saed.py`: KS [0 1 -1]_p → exact <111>_c variant
  (24 variants), NW [1 -1 0]_p → exact <100>_c (12 variants), Bain [001]
  collinear (220)_p/(200)_c and pinned 45 deg parent-child (200) split,
  shared-basis identity `B_c = V_i B_p`, composite-wide rotation
  equivariance, g-alignment, child-config override, subset/order/error
  paths, describe() content. Gates: 959 passed, ruff clean, mypy clean,
  integrity passed.

- (CD3) `src/pytex/plotting/composite_saed.py`: `SpotStyle` (marker, color,
  filled/hollow, three intensity->size modes with floor, alpha, edge,
  z-order; validated), colorblind-aware `VARIANT_COLOR_PALETTE` (Tol muted
  + 2) x `VARIANT_MARKER_CYCLE` (8 markers -> 96 distinct combos),
  `CompositeSAEDPlotConfig` (per-variant style overrides, render-time
  variant subset, mm | inv_angstrom axes units, transmitted-beam marker,
  legend cap/outside placement, title/background/figsize/dpi/padding;
  validated), `render_composite_saed` (children first, hollow parent on
  top, equal aspect, machine-readable per-collection gids
  `pytex-composite:*`). Also lowered the default zone-label rationalization
  bound 12 -> 6 (nearest low-index zone labels; the KS variant 5.26-deg-off
  [1 0 0] label now surfaces the classic KS-NW separation). 28 structural
  tests in `tests/unit/test_composite_saed_plotting.py` (gid census, offset
  parity with spot tables, subset/beam/parent toggles, explicit style
  override + palette cycling colors, hollow parent facecolors, size-mode
  arithmetic, unit rescaling + axis labels, legend content/cap/disable,
  title logic, config validation). Visual smoke check rendered KS
  [0 1 -1] composite: exact [1 1 -1] variant overlays coincident parent
  spots as expected. Gates: 987 passed, ruff clean, mypy clean, integrity
  passed.

- (CD4) Annotation engine in `plotting/composite_saed.py`: `format_hkl`
  (plain + crystallographic overline mathtext, compact/spaced, pinned),
  `SpotAnnotationConfig` (budget, intensity floor, coincidence tolerance,
  offset, fonts, leader lines, overlap avoidance; validated),
  cKDTree+union-find coincidence clustering into phase-tagged multi-line
  merged labels ('p' / 'Vn'), deterministic greedy two-ring compass
  placement using measured display-space text extents (label-label and
  label-over-spot collision checks; drop-not-overlap policy), optional
  leader lines for outer-ring placements, `AnnotationResult` report with
  `describe()`, `render_composite_saed(..., return_annotations=True)`.
  Key fix discovered by tests: `fig.tight_layout()` must run *before*
  placement or the measured extents go stale. 25 regression tests in
  `tests/unit/test_composite_saed_annotations.py` (format pins incl.
  overlines, merged parent/variant labels on the KS composite,
  pairwise-disjoint rendered label boxes, budget/floor (floor test uses
  strongly Debye-Waller-damped phases since the Z-only proxy is flat),
  disabled/deterministic/plain-format paths, leader-line engagement,
  result consistency validation). Visual check: KS [0 1 -1] + V2 figure
  shows merged '(-111) p / (011) V2' coincidence labels, color-coded
  single-phase labels, no collisions. Gates: 1012 passed, ruff clean,
  mypy clean, integrity passed.

- (CD5) Coincidence analysis + sweep utilities in
  `diffraction/composite.py`: `SpotCoincidence` (validated pair with both
  hkl, detector coordinates, separation, label), `SpotCoincidenceReport`
  (per-variant totals, exact-count queries, tolerance/consistency
  invariants enforced at construction, convention-explicit `describe()`),
  `find_spot_coincidences` (per-variant cKDTree query_ball_point; sorted by
  separation), `sweep_parent_zone_axes` lazy iterator. 16 regression tests
  in `tests/unit/test_composite_coincidences.py`: pinned totals for the KS
  [0 1 -1] composite (24 pairs at 2.5 mm, 0 at 1 mm, 48 at 5 mm), each
  exactly-oriented variant contributes exactly its two antipodal
  close-packed {111}_p||{011}_c pairs at the analytically derived
  (sqrt(2)/2.87 - sqrt(3)/3.6)*180 = 2.0938 mm separation (simulation
  matches analytic to 1e-13), full brute-force O(N^2) parity on the pair
  set, tolerance monotonicity, sort order, error paths, report validation,
  lazy ordered sweep. Gates: 1028 passed, ruff clean, mypy clean,
  integrity passed.

- (CD6) Integration and documentation. Exports: `pytex.diffraction`
  re-exports the full CD1-CD5 public surface (`simulate_zone_axis_spots`,
  `SpotTable`, `KinematicSimulationConfig`, `zone_basis_from_axis`,
  `electron_wavelength_angstrom`, `centering_allowed_mask`,
  `electron_structure_factors`, `simulate_composite_saed`,
  `CompositeSAEDPattern`, `VariantZonePattern`, `rationalize_zone_axis`,
  `RationalizedZoneAxis`, `find_spot_coincidences`, `SpotCoincidence`,
  `SpotCoincidenceReport`, `sweep_parent_zone_axes`); `pytex.plotting`
  re-exports `render_composite_saed`, `CompositeSAEDPlotConfig`,
  `SpotStyle`, `SpotAnnotationConfig`, `AnnotationResult`, `format_hkl`.
  (Top-level `pytex` intentionally not touched — composite machinery is
  accessed via the subpackages, matching the rest of the diffraction
  surface.) Two worked examples added
  (`worked_examples/examples/composite_diffraction.py`): 200 kV relativistic
  wavelength vs De Graef Table 2.2 (0.02508 A) and the KS exact child-zone
  identity ([0 1 -1]_p -> nearest child zone deviation 0 deg); gallery
  regenerated (`docs/site/examples/generated/composite-diffraction.md` +
  index), `tests/unit/test_worked_examples.py` green. New workflow page
  `docs/site/workflows/composite_or_diffraction.md` wired into the Sphinx
  toctree (fixture-based ni_fcc->fe_bcc KS example verified to run end to
  end). CHANGELOG Unreleased/Added entry. Gates: 1030 passed, ruff clean,
  mypy clean, integrity passed, documentation + reference policy tests
  green, Sphinx build succeeded with no warnings.

### CD7 — Burgers beta->alpha as the canonical hexagonal case (follow-on)

Requested after v1: treat Burgers (bcc beta -> hcp alpha; Ti/Zr/Hf) as a
canonical case alongside KS across tests, demos and illustrations. Burgers
is the hexagonal counterpart that KS cannot exercise: a non-cubic child
lattice, 622 child symmetry, and four-index notation.

- **Miller-Bravais labeling** (the real gap for hexagonal illustration
  correctness): new `is_hexagonal_phase(phase)` (via
  `PointGroup.crystal_system`, trigonal deliberately excluded);
  `RationalizedZoneAxis` gained `indices_bravais` with `U+V+T=0` validation
  and `label()` preferring it; `rationalize_zone_axis` populates it for
  hexagonal phases; `SpotCoincidence` gained `parent_bravais`/`child_bravais`
  flags wired from the relationship phases; `format_hkl(..., bravais=True)`
  expands `(hkl)` to `(h k i l)`. Cubic phases are untouched (three-index).
- **Shared fixture** `make_bcc_hcp_phases()` (beta-Ti a = 3.3065 A;
  alpha-Ti a = 2.9508 A, c = 4.6855 A, P6_3/mmc two-atom motif), reused
  across the composite test modules like `make_fcc_bcc_phases`.
- **Pinned Burgers science** (both defining parallelisms give exactly
  rational views): parent `<110>` -> child `[0001]` basal zone at 0 deg
  ({110}_bcc || (0001)_hcp); parent `<111>` -> child `<11-20>` at 0 deg
  (<-111>_bcc || <11-20>_hcp); 12 variants; the basal view contains only
  hk0 reflections with max d = a*sqrt(3)/2 = 2.5555 A ({10-10}); the basal
  pattern is invariant under 60 deg rotation (six-fold check); and the
  practical TEM signature — {110}_bcc superimposed on (0002)_hcp at
  (sqrt(2)/a_bcc - 2/c_hcp) * 180 = **0.15450 mm**, matched by simulation to
  14 significant figures (a hand-arithmetic slip on this constant was caught
  by the test, which is exactly why it is pinned analytically).
- **Two new worked examples** (Burgers exact basal zone; the {110}/(0002)
  near-coincidence) and 15 new regression tests across the composite,
  coincidence and annotation suites.
- Gates: 1047 passed, ruff clean, mypy clean, integrity passed.

### CD8 — Executed teaching notebook (follow-on)

`docs/site/tutorials/notebooks/21_composite_or_diffraction_patterns.ipynb`
(28 cells, 6 embedded figures), generated by
`scripts/generate_tutorial_notebooks.py` and committed **executed** via
`scripts/execute_notebooks.py --only 21` per repo policy (the site builds
with `nb_execution_mode = "off"`). Structure:

1. Theory — relativistic wavelength and the Ewald-sphere radius; excitation
   error `s_g = g_z - lambda g^2 / 2` with an original two-panel diagram
   (Ewald sphere vs the ZOLZ plane, and the analytic `s_g` curve with the
   simulated spots landing on it); why excitation-error selection rather
   than the integer zone law is required for irrational OR-mapped child
   zones; the shared parent-anchored detector construction.
2. KS (cubic canonical case) — 24 variants, four exactly-oriented on
   <111>, a sorted child-zone deviation chart whose maximum is the 5.26 deg
   KS-NW separation, the rendered composite with merged coincidence labels,
   the coincidence report, and a coincidence-count-vs-tolerance step plot
   annotated with the analytic close-packed separation.
3. Burgers (hexagonal canonical case) — both defining parallelisms shown as
   exactly rational views ([110]->[0001], [111]->[11-20]), the six-fold
   basal pattern (hk0-only, max d = a*sqrt(3)/2, 60 deg invariance to
   1e-13 mm), four-index Miller-Bravais labels in the rendered figure, and
   the {110}_beta/(0002)_alpha fingerprint at 0.1545 mm with simulated and
   analytic values printed side by side.
4. Configuration showcase (styling, per-variant overrides, reciprocal-space
   axes, in-plane rotation) and `describe()` explainability.

Wired into the notebook atlas toctree and the OR teaching-track list, and
cross-linked from the composite workflow page. Gates: 1047 passed, ruff
clean, mypy clean, integrity passed, notebook policy tests green, Sphinx
build succeeded (6 figures + math verified in the generated HTML).

### CD8a — Notebook-output repair and a documented footgun

Found while committing CD8: `scripts/generate_tutorial_notebooks.py` rewrites
**every** notebook with empty outputs, so running it discards stored outputs
of notebooks that were not the target. Regenerating for notebook 21 wiped
notebook 20's executed outputs (2 figures), and that loss was briefly
committed in `407b92f`. Repaired by re-executing 18, 19 and 20; notebooks 18
and 19 turned out to have been **already** unexecuted before this work
(a pre-existing gap, despite earlier ledger entries describing them as
committed executed), so the whole OR teaching track now renders for the
first time: nb18 4 images, nb19 5, nb20 3, nb21 6 in the built HTML.
A warning documenting the regeneration hazard was added to
`scripts/execute_notebooks.py`.

### CD8b — Notebook generator retired; notebooks are hand-authored

Decision (2026-07-20, requested by the repository owner): remove
`scripts/generate_tutorial_notebooks.py` entirely. The `.ipynb` files under
`docs/site/tutorials/notebooks/` are now the source of truth and are edited
directly. Rationale: the generator constrained how notebooks could be
written, added a rewrite-then-execute round trip to every edit, and provided
no benefit that the notebooks themselves did not already provide (the test
suite always read the `.ipynb` files, never the generator).

An audit taken while removing it justified the decision beyond the
ergonomics: **notebooks 01-17 had never been executed at all** — zero
outputs, zero figures — so 17 of the 21 tutorial pages had been rendering on
the Sphinx site as bare code listings. The generator's rewrite-everything
behaviour is what kept resetting them, and because nothing asserted that
notebooks carry outputs, the gap was invisible. Removing the generator,
executing every notebook, and adding an executed-notebook regression test
(`tests/unit/test_notebooks.py`) replaces illusory automation with an actual
guarantee.

Cleaned references: `pyproject.toml` (ruff per-file-ignore),
`docs/site/tutorials/installation_and_build.md` (build instructions rewritten
around editing and executing notebooks), `scripts/execute_notebooks.py`
(docstring), and both roadmap ledgers.

Outcome:

- **All 21 notebooks now execute and are committed executed** — 54 stored
  figures, 59 images rendered across the built site (previously 13), plus
  printed outputs on the text-oriented notebooks. Notebooks 01-17 render
  real results for the first time.
- **One genuine API-drift bug surfaced and was fixed** in notebook 07: its
  `CrystalMap` put orientations in the specimen frame and the grid in the map
  frame, then called `to_experiment_manifest()`, which now (correctly) refuses
  to guess the relationship and demands an `AcquisitionGeometry` carrying an
  explicit `specimen_to_map` transform. The notebook now constructs that
  transform, which also makes the frame contract visible to the reader. This
  bug had been latent precisely because the notebook was never executed.
- **Two guard tests added** (`tests/unit/test_notebooks.py`):
  `test_every_notebook_is_committed_executed` (every non-empty code cell has
  an `execution_count`) and `test_no_notebook_contains_error_output`. Both
  were verified to actually fail when a notebook is sabotaged, and they name
  the offending notebook and cell. This is the durable replacement for the
  generator's illusory guarantee.

## Program v1 outcome and follow-ons

All six phases (CD0-CD6) landed as verified commits. The library now
simulates and renders composite kinematic SAED patterns for any OR / any
parent zone axis with: a vectorized excitation-error-based engine, shared
parent-anchored detector geometry, exact irrational child zones with
nearest-rational labeling, per-variant styling, merged collision-free
annotations, and a quantitative coincidence report — all with `describe()`
explainability, regression tests on every critical part, worked examples,
and a workflow page. Natural follow-ons (explicitly out of scope here):

- Double-diffraction spot prediction (kinematically forbidden spots excited
  via g1+g2 paths) as an optional overlay.
- HOLZ ring / first-order-Laue-zone support (needs relrod + curvature care).
- Interactive backend (plotly / GUI) reusing `CompositeSAEDPlotConfig`, per
  the eventual desktop-GUI goal.
- Dynamical (Bloch-wave / multi-beam) intensities — deliberately excluded
  from this program.

## Next actions

- Program complete; no open actions.
