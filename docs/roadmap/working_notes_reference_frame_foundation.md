# Working Notes — Reference Frame Foundation Program

Durable plan and phase ledger for the reference-frame foundation. This file is the resumption
point: if the working session is interrupted, read this file first, then continue at the first
phase whose status is not `DONE`.

## Objective

Make reference frames a first-class, shared foundation of PyTex, so that every module handles,
creates, manipulates, transforms, represents, and visualizes frames through **one** model:

1. **Handling and creation** — a canonical frame type carrying geometry (axis vectors), semantic
   identity (name, domain, handedness), display labels, and provenance; plus a catalog of the
   standard frames every user expects to reach for by name (Cartesian `X/Y/Z`, sample `RD/TD/ND`,
   crystal `a/b/c`, EBSD map, detector, laboratory, reciprocal).
2. **Manipulation and transformation** — an algebra of frame-to-frame transforms with explicit
   source/target typing, ergonomic constructors (rotation, Euler, axis-angle, axis correspondence,
   vector alignment), inversion, composition, and a **frame graph** that resolves the transform
   between any two registered frames automatically.
3. **Representation** — `describe()` prose for every stable object, JSON contract round-tripping,
   and consistent `repr`.
4. **Visualization and embedding** — one triad model that renders three ways: as scene primitives
   for 3D crystal/world scenes, as a small **embeddable gizmo** for existing 2D axes (SAED
   diffractograms, pole figures, IPF maps, crystal viewer), and as **standalone SVG** for the
   documentation system (no matplotlib required).
5. **Consistency** — every module that previously built frames ad hoc now uses the catalog, so
   frame identity is stable across module boundaries.

## Governing constraints (from `AGENTS.md` and the standards)

- No subsystem may define its own private frame model; the shared core model must be expanded
  instead (`AGENTS.md`, Non-Negotiable Rules).
- The stable frame-domain vocabulary is fixed: `crystal, specimen, map, detector, laboratory,
  reciprocal`. **No new domains may be invented**
  (`docs/standards/notation_and_conventions.md`).
- Frame equality is load-bearing: `VectorSet`, `FrameTransform`, `Orientation`, and `SymmetrySpec`
  compare frames for consistency. New fields therefore must be **hashable, comparable, and
  default to values that preserve equality with every existing construction site**.
- Stable report/result objects carry `describe()`; JSON contracts and `describe()` stay in
  lockstep.
- Canonical SVG must follow `docs/standards/visualization_style_guide.md` (tokens, Arial,
  `<title>`/`<desc>`).
- mypy runs in `strict` mode; ruff `E,F,I,B,UP,N,RUF` at line length 100.
- Tests treat warnings as errors; matplotlib figures opened by tests must be closed.

## Key design decisions

### D1 — `ReferenceFrame` gains geometry without breaking identity

`ReferenceFrame` keeps its role as the *semantic identity* of a frame, and additionally carries the
components of its three labelled axes **in the canonical right-handed Cartesian reference**
`X, Y, Z` (`frame_catalog.CARTESIAN_FRAME`). Quoting every frame's geometry against one shared
reference is what makes `FrameTransform.between_frames` well defined.

- Stored as `axis_vectors: tuple[tuple[float, float, float], ...]` (**not** an `ndarray`) so the
  dataclass stays comparable and hashable; `basis_matrix` exposes a read-only `(3, 3)` array whose
  **columns** are the axis vectors.
- Default is the identity triad, so every frame built by existing code compares exactly equal to
  the catalog frame with the same name/domain/axes. This is what allows migration to the catalog
  without touching downstream equality checks.
- Validation: axes must be linearly independent, and `sign(det)` must agree with `handedness`.
  Construction-time invariants are preferred over downstream recovery.
- Non-orthonormal axis vectors are allowed (an oblique crystal frame is legitimate);
  `is_orthonormal` reports the distinction, and visualization normalizes for legibility.

### D2 — Frame *geometry* vs lattice `Basis`

`Lattice.direct_basis()` already returns a `Basis` with physical units (angstrom) and a
`BasisKind`. `ReferenceFrame.axis_vectors` is dimensionless axis *orientation* only. The two are
complementary and must not be merged: a `Basis` says how long the crystal axes are, a frame says
which directions its labels point.

### D3 — Frame graph, not a frame tree

Workflows relate frames pairwise (`crystal -> specimen`, `specimen -> map`,
`specimen -> detector`), and different workflows declare different subsets. A `FrameGraph` holding
transforms as undirected edges (each edge usable in both directions via `inverse()`) with
breadth-first path resolution matches how the canonical frame chain is actually used, and returns
the *shortest* chain so composed numerical error stays minimal.

### D4 — Three renderers, one triad model

`FrameTriad` geometry (origin, three axis vectors, labels, colors) is computed once in
`pytex.plotting.frames` and consumed by:

- `AxisTriad3D` primitives → existing 3D scene renderers;
- a 2D projected gizmo drawn into any existing matplotlib `Axes` (inset), for embedding;
- a pure-Python SVG emitter with no matplotlib import, for the documentation system.

The SVG path must stay import-light: `pytex.plotting.frames` may not require matplotlib at import
time (matching the existing lazy `_require_matplotlib` pattern).

### D5 — Axis colors

Reuse the existing Okabe-Ito-derived triad palette already fixed in
`pytex.plotting.primitives._DEFAULT_TRIAD_COLORS` (`#1d4ed8`, `#059669`, `#dc2626`) so a triad
looks identical whether it is drawn in a 3D scene, an inset gizmo, or a documentation SVG.

## Phase ledger

Each phase is a verified unit of work. Status values: `TODO`, `IN PROGRESS`, `DONE`.

| Phase | Scope | Status |
| --- | --- | --- |
| RF1 | Working notes + plan (this file); active-task pointer updated | DONE |
| RF2 | `ReferenceFrame` geometry, invariants, helpers, `describe()` | DONE |
| RF3 | `FrameTransform` constructors + `FrameGraph` resolution | DONE |
| RF4 | `pytex.core.frame_catalog` standard frames | DONE |
| RF5 | `pytex.plotting.frames` (primitives, gizmo, SVG) | DONE |
| RF6 | Migrate all ad-hoc frame construction to the catalog | DONE |
| RF7 | JSON contract round-trip for the new fields | DONE |
| RF8 | Unit tests for frames, catalog, graph, visualization | DONE |
| RF9 | Docs, canonical SVG figure, worked examples, CHANGELOG | DONE |
| RF10 | Full verification: ruff, mypy strict, pytest, integrity | DONE |
| RF11 | Wire the gizmo into the SAED, composite-SAED, and crystal renderers | DONE |

## Ad-hoc frame construction sites to migrate (RF6 checklist)

- [x] `src/pytex/adapters/scan_files.py::default_ebsd_frames`
- [x] `src/pytex/plotting/_plotting_validation_cases.py::_make_crystal_frame/_make_specimen_frame`
- [x] `src/pytex/diffraction/saed.py` detector frame
- [x] `src/pytex/core/lattice.py::Lattice.reciprocal_basis` reciprocal frame
- [x] `src/pytex/cli.py::_cmd_core_demo`
- [x] `src/pytex/contracts.py` deserialization (new fields, backward compatible)
- [x] `src/pytex/plotting/primitives.py::reference_frame_triad` (honour `axis_vectors`)

Equality note: catalog defaults are pinned to the exact field values these sites already used
(`crystal`/`a,b,c`, `specimen`/`x,y,z`, `map`/`x,y,z`), so migration is identity-preserving. This
is asserted directly in `tests/unit/test_frame_catalog.py`.

## Verification record

- Baseline before the program (commit `e453e85`): full `pytest` suite green, **1058 tests
  collected, 1058 passed**. Any later count must be at least this.
- RF10 final run (all green):
  - `python -m ruff check .` - All checks passed
  - `python -m mypy src` - no issues in 86 source files (strict mode)
  - `python scripts/check_repo_integrity.py` - passed, including the four new required paths
  - `python -m sphinx -b html docs/site docs/_build/html` - build succeeded; the 5 remaining
    warnings are pre-existing `docs/standards/reference_canon.md` cross-references, unrelated to
    this program
  - `python -m pytest` - **1191 passed** (baseline 1058; +133 from the three new frame test
    modules plus the RF11 renderer-integration tests), no warnings
  - `python -m pytest --cov=pytex --cov-fail-under=87` - **89.31%** total, above the CI gate;
    new modules measured at `core/frame_catalog.py` 100%, `core/frames.py` 99%,
    `plotting/frames.py` 99%
  - `python scripts/generate_worked_examples.py` - gallery regenerated, all five new examples
    within tolerance
  - `python scripts/execute_notebooks.py --only 01` - notebook 01 re-executed with stored outputs

## Planned public surface

Recorded here so a resuming session knows the intended API shape. Items are added to the
"Outcomes" section only once implemented and tested.

- `ReferenceFrame`: `axis_vectors`, `axis_descriptions`, `basis_matrix`, `axis_vector()`,
  `axis_index()`, `is_orthonormal`, `is_right_handed`, `with_axis_vectors()`, `renamed()`,
  `rotated()`, `describe()`.
- `FrameTransform`: `from_rotation`, `from_bunge_euler`, `from_axis_angle`,
  `from_axis_correspondence`, `between_frames`, `as_rotation()`, `rotation_angle_deg`,
  `rotation_axis`, `apply_to_frame()`, `describe()`.
- `FrameGraph`: `add_frame`, `add_transform`, `has_frame`, `frame`, `path`, `transform_between`,
  `convert`, `describe()`.
- `pytex.core.frame_catalog`: `CARTESIAN_FRAME`, `SPECIMEN_FRAME`, `SAMPLE_RD_TD_ND_FRAME`,
  `CRYSTAL_FRAME`, `MAP_FRAME`, `DETECTOR_FRAME`, `LABORATORY_FRAME`, `STANDARD_FRAMES`,
  `get_standard_frame`, `list_standard_frames`, and builders `cartesian_frame`,
  `specimen_frame`, `sample_frame`, `crystal_frame`, `map_frame`, `detector_frame`,
  `laboratory_frame`, `reciprocal_frame_for`, `rolling_frame_graph`.
- `pytex.plotting.frames`: `FrameTriad`, `frame_triad`, `frame_triad_primitives`,
  `plot_reference_frame`, `plot_frame_relationship`, `add_frame_indicator`,
  `reference_frame_svg`, `frame_catalog_svg`.

## Outcomes

Everything in "Planned public surface" landed, with two corrections found during RF5 verification
and recorded here because they changed the design:

**Correction 1 — `from_axis_correspondence` is about components, not canonical geometry.** The
first implementation built the rotation from the two frames' canonical-Cartesian axis vectors
(`R = [+/-b_j] [a_i]^-1`). That is a type error: it feeds canonical components into a map whose
domain is *source-frame* components, and it only happens to be right when both frames have identity
geometry. In its own coordinates a frame's axis `i` is `e_i`, so the declaration "source axis `i` is
target axis `j`" fixes `R e_i = +/- e_j` — a signed permutation matrix, independent of where either
frame's axes point. Fixed, and pinned by
`test_axis_correspondence_is_about_components_not_where_frames_point`.

**Correction 2 — `apply_to_frame` was dropped in favour of `source_axes_in_target`.** The planned
`apply_to_frame` returned a `ReferenceFrame` whose axis vectors were expressed in the *target*
frame, silently violating the D1 convention that axis vectors are canonical-Cartesian components.
It also mapped `a_i` (canonical) rather than `e_i` (source components) through `R`, the same type
error as correction 1. It is replaced by `source_axes_in_target()`, which returns the `(3, 3)`
matrix and documents that those components are target-frame components. `FrameTriad` gained an
optional `basis` override so `plot_frame_relationship` can draw both triads in target-frame
coordinates without fabricating a misleading frame.

Two smaller visualization fixes, both caught by rendering and inspecting the output:

- SVG arrowheads needed `markerUnits="userSpaceOnUse"`; the SVG default scales markers by stroke
  width, so a 3.4-unit stroke produced 30-unit arrowheads. Pinned by
  `test_reference_frame_svg_draws_one_line_per_axis_with_arrowheads`.
- Axis labels are now placed a fixed distance beyond the projected tip with a minimum radius, so an
  axis pointing nearly at the viewer (which projects to a near-zero-length arrow) still gets a label
  clear of the origin. Shared by the gizmo and the SVG path; pinned by
  `test_frame_indicator_labels_stay_clear_of_the_origin`.

Delivered files:

- `src/pytex/core/frames.py` (rewritten), `src/pytex/core/frame_catalog.py` (new),
  `src/pytex/plotting/frames.py` (new).
- RF11 wiring (all opt-in, so no existing figure changes):
  `plot_saed_pattern(show_frame_indicator=True)` draws the detector `u`/`v` axes;
  `CompositeSAEDPlotConfig(show_frame_indicator=True)` draws the **parent crystal** axes as they
  land on the detector (via `zone_basis_parent.T`, the genuinely informative frame for a composite
  pattern — the detector axes there are trivially the page axes);
  `plot_crystal_structure_3d(show_frame_indicator=True)` draws the phase's `a`/`b`/`c` axes from
  the lattice basis at the figure's own view angles. `add_frame_indicator` gained a `basis`
  override to make these possible without violating the canonical axis-vector convention.
- `tests/unit/test_frames.py` (rewritten), `tests/unit/test_frame_catalog.py` (new),
  `tests/unit/test_frame_visualization.py` (new) — 126 tests.
- `docs/architecture/reference_frame_foundation.md` (new) + Sphinx stub and toctree entry;
  `docs/site/concepts/reference_frames_and_conventions.md` extended; API guide, notation standard,
  terminology registry, `docs/README.md` index, and `CHANGELOG.md` updated.
- `scripts/generate_reference_frame_figures.py` (new) producing
  `docs/figures/reference_frame_catalog.svg` and `docs/figures/sample_frame_rd_td_nd.svg`.
- `worked_examples/examples/reference_frames.py` (new, 5 examples) + regenerated gallery.
- `scripts/check_repo_integrity.py` extended with the new required paths.

## Follow-ons (not in this program)

- Attach a declared `FrameGraph` to `CrystalMap`/scan imports so vendor frame conventions
  (`specimen -> map` rotations) are recorded per dataset rather than assumed identical.
- Extend `add_frame_indicator` to hexagonal four-axis (`a1/a2/a3/c`) gizmos once the four-index
  display policy in `hexagonal_and_trigonal_conventions.md` is extended to gizmo labels.
- Offer the gizmo on the pole-figure, IPF-map, and stereographic renderers too. The capability
  works there already (`add_frame_indicator` is tested on polar axes); only the opt-in keyword is
  missing.
- De-conflict overlapping labels in `plot_frame_relationship` when a source axis lands exactly on
  a target axis (the coincident label is currently drawn under the other, which is geometrically
  honest but slightly hard to read).
- Consider promoting `docs/figures/reference_frames_vectors.svg` (hand-authored) to a generated
  asset now that `frame_catalog_svg` exists.
