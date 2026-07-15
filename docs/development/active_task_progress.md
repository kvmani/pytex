# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Objective

Build a **thorough, robust core visualization-primitives layer** so that every basic
crystallographic entity — a bare 3D vector, a Miller direction, a Miller plane, a unit cell, a
reference-frame triad, a point/lattice cloud — is a first-class, *composable* scene primitive, in
**both** real-space 3D and on stereographic / IPF projections. The layer must let complex composite
figures (canonical example: **two crystals in a specific orientation relationship**) be assembled
with minimal code by reusing these primitives, without rewriting the proven single-crystal
`CrystalScene` / `plot_crystal_structure_3d` renderer.

## Current Status

- Started: 2026-07-15
- Branch: `main` (tracking `origin/main`)
- Baseline commit: `cd8834b`
- Phase: implementation

## Key facts established (verified against live code)

- Existing 3D real-space stack: `pytex.plotting.crystal3d` — `CrystalScene` (atoms/bonds/cells/
  planes/directions/polyhedra glyphs, all in a single Phase crystal Cartesian frame),
  `build_crystal_scene`, `plot_crystal_structure_3d` (Blinn-Phong lit, ONE depth-sorted
  `Poly3DCollection`). This is single-crystal, single-frame; there is no placement transform and no
  composite/world scene — that is the core gap.
- Existing generic 2D/3D scene specs: `pytex.plotting._render` (`FigureSpec2D`/`FigureSpec3D` +
  layers, `render_figure_spec`).
- Existing stereographic/IPF: `pytex.plotting.spherical` (`plot_crystal_directions/planes`,
  `plot_wulff_net`, `plot_symmetry_elements`), `pytex.plotting.ipf` (`IPFColorKey`, `plot_ipf_key`).
- Core data model: `CrystalDirection.unit_vector`, `CrystalPlane.normal`, `Lattice.direct_basis()`,
  `Orientation.map_crystal_vector` (crystal->specimen = `g @ v`), `Rotation.as_matrix`.
- **OR placement convention (verified numerically):** with world = parent crystal frame, place the
  child crystal by `relationship.parent_to_child_rotation.inverse().as_matrix()`. This makes both
  parallel directions AND parallel plane normals coincide exactly (dot = 1.0) — see KS fcc->bcc
  check. This is the correctness anchor for composite OR figures.

## Design (renderer-independent primitives + world composition)

New module `src/pytex/plotting/primitives.py`:
- `Transform3D` (linear 3x3 + translation; `identity/from_rotation/from_orientation/from_matrix`,
  `apply_points/apply_vector/compose/inverse`).
- Immutable primitives in a shared world Cartesian (angstrom): `Arrow3D`, `PolyLine3D`,
  `PlanePatch3D`, `PointCloud3D`, `Label3D`, `AxisTriad3D`.
- `PrimitiveScene3D` container: `transformed`, `merge`, `bounds`.
- Builders from crystallographic objects: `vector_arrow`, `direction_arrow`, `plane_normal_arrow`,
  `crystal_plane_patch`, `reference_frame_triad`, `unit_cell_polylines`, `lattice_point_cloud`.
- `render_primitive_scene_3d` (standalone or onto an existing 3D ax).

Refactor `crystal3d.py` (behavior-preserving): extract `_accumulate_crystal_mesh` and
`_draw_crystal_non_mesh`; add `CrystalScene.transformed(transform)` (rigid transforms only).

New module `src/pytex/plotting/scene3d.py`:
- `PlacedCrystal` (CrystalScene + Transform3D + label), `WorldScene3D` (placed crystals + primitives),
  `WorldScene3D.from_orientation_relationship(...)`, `render_world_scene_3d` (ONE global
  depth-sorted mesh across all crystals + primitive rendering + optional world/crystal triads).

2D bridge (spherical.py): `plot_stereographic_vectors` — arbitrary world-frame unit vectors (poles
and/or great-circle traces) on a Wulff net, so the SAME primitives project to a stereogram; enables
a multi-crystal (parent+child) stereographic overlay analog of the OR figure.

## Verification Gates (must all pass before done)

- `python -m pytest -q`
- `python scripts/check_repo_integrity.py`
- `python -m ruff check .`
- `python -m mypy src`
- `python -m sphinx -b html docs/site docs/_build/html`

## Completed

- Explored plotting subsystem, core model, tests, style themes; verified OR placement convention.
- Implemented `src/pytex/plotting/primitives.py`: `Transform3D`; `Arrow3D`, `PolyLine3D`,
  `PlanePatch3D`, `PointCloud3D`, `Label3D`, `AxisTriad3D`; `PrimitiveScene3D`; builders
  (`vector_arrow`, `direction_arrow`, `plane_normal_arrow`, `crystal_plane_patch`,
  `reference_frame_triad`, `unit_cell_polylines`, `lattice_point_cloud`); `render_primitive_scene_3d`
  + shared `_draw_primitive_scene`.
- Refactored `crystal3d.py` behavior-preservingly: extracted `_draw_crystal_frame`,
  `_accumulate_crystal_mesh`, `_draw_crystal_planes_and_directions`; added `CrystalScene.transformed`
  (rigid-only).
- Implemented `src/pytex/plotting/scene3d.py`: `PlacedCrystal`, `WorldScene3D`,
  `WorldScene3D.from_orientation_relationship`, `render_world_scene_3d` (ONE global depth-sorted mesh
  across crystals).
- Added `plot_stereographic_vectors` + `build_vector_stereogram_figure_spec` to `spherical.py`;
  runtime wrapper.
- Wired all exports through `plotting/__init__.py`, `runtime.py`, and top-level `pytex/__init__.py`.
- Tests: `tests/unit/test_plotting_primitives.py` (15) and `tests/unit/test_scene3d_composition.py`
  (17), incl. OR parallelism correctness and single-mesh depth-sort checks.
- Worked examples: `worked_examples/examples/visualization_composition.py` (2 examples:
  Transform3D↔map_crystal_vector identity = 0; KS parallel-direction cosine = 1); regenerated gallery
  (`docs/site/examples/generated/visualization.md`).
- Docs: `docs/site/concepts/visualization_primitives.md` + toctree wiring; registered symbol
  \(\mathbf{T}\) in the terminology registry.

## Phase 2: VESTA parity (2026-07-15, same session)

User directive: match or exceed VESTA quality/functionality for the crystal viewer. Delivered:

- **Render styles** (`render_style=` on `build_crystal_scene` / `plot_crystal_structure_3d`):
  `ball_and_stick` (default), `space_filling` (Slater atomic radii via new
  `_chemistry.atomic_radius_angstrom` / `display_radius_angstrom`, bonds suppressed), `stick`
  (uniform cylinders), `wireframe` (bond lines only, `atom_render_mode="none"`), `polyhedral`
  (auto species selection). Presets merge UNDER user style_overrides (user wins).
- **Occupancy pie-spheres** (VESTA signature): sites sharing a position render as azimuthal sectors
  of one shared-radius sphere; occupancy < 1 leaves a vacancy sector (`vacancy_color`). Automatic
  from `AtomicSite.occupancy`; glyph fields `occupancy` / `sector_start` / `vacancy_fraction`.
- **Atom labels** (`atom_label_mode="species"|"site"`), **site vectors** (`site_vectors={label: v}`,
  moments/displacements on every periodic copy), **depth cueing** (`depth_cue_strength` theme key,
  static per view, both single-crystal and world renderers).
- **Measurement**: `CrystalScene.bond_lengths_angstrom()` + `bond_length_summary()` (species-pair
  stats) — exceeds VESTA's click readout (scriptable/testable).
- **Quality fixes found by visual inspection**: `CrystalScene.bounds()` now includes atom radii
  (space-filling spheres no longer clipped); wireframe hides atom bodies (scatter discs removed).
- **Governance**: `docs/testing/vesta_parity_matrix.md` (authoritative ledger, registered in
  `check_repo_integrity.py` REQUIRED_PATHS + AGENTS.md refs + site validation toctree via include
  stub `docs/site/validation/vesta_parity_matrix.md`); concept page section; worked example
  `viz-scene-bond-length-halite-identity` (a/2 = 2.0 exact); gallery regenerated.
- Known deferred (in matrix as planned): thermal ellipsoids (needs anisotropic ADPs in
  `AtomicSite`), dashed/H-bond styles, angle/dihedral helpers, ionic radii. Spawned separate task
  for the dead pymatgen path in `covalent_radius_angstrom` (changes bond chemistry; needs own
  validation pass).

## Verification Results (all gates green)

- `python -m pytest`: 804 passed (+23 VESTA tests in tests/unit/test_vesta_parity_features.py).
- `python -m ruff check .`: all checks passed.
- `python -m mypy src`: success, 78 source files.
- `python scripts/check_repo_integrity.py`: passed.
- `python -m sphinx -b html docs/site docs/_build/html`: build succeeded.
- Visual proof rendered and inspected: five-style gallery (clipping + wireframe fixed and
  re-verified), occupancy/moments/labels/depth-cue figure, primitives showcase, KS two-crystal OR
  scene, stereographic overlay.

## Next Actions

1. Optional follow-ups (tracked in vesta_parity_matrix.md as planned): thermal ellipsoids,
   dashed/hydrogen-bond styles, angle/dihedral measurement, ionic-radius space filling; promote a
   canonical SVG of the OR schematic into `docs/figures/`.
2. Commit and push when the user requests it (not yet requested).
