# Notebook Improvement Progress

Durable handoff note for the tutorial-notebook overhaul under
`docs/site/tutorials/notebooks/`.

## Objective

Rewrite the tutorial notebooks so each is a genuinely pedagogical, scientifically
rigorous, executable demonstration of the library: theory + conventions + citations,
focused imports (no kitchen-sink boilerplate), live-computed results verified against
analytic/literature values, and at least one figure. Notebooks must execute error-free
(the Sphinx build runs them) and pass `tests/unit/test_notebooks.py`.

## Motivation / baseline problem

Notebooks 01–16 were thin, machine-generated stubs carrying identical dead boilerplate
(an ~80-name kitchen-sink `from pytex import (...)` plus five unused helper functions
`make_context`, `describe_phase_fixture`, `load_zr_hcp_phase`, `load_diamond_phase`,
`publication_crystal_style`). Notebooks 17–21 (OR track) were already hand-authored and
richer — used as the quality template.

## Build tooling

- **Superseded 2026-08-08.** Notebooks are authored and edited by hand as `.ipynb` files,
  with no generator and no execution step. They are committed **without outputs**: the Sphinx
  site sets `nb_execution_mode = "cache"` and executes them at build time, and
  `nb_execution_raise_on_error = True` fails the build on a notebook that no longer runs.
  `scripts/execute_notebooks.py` existed only to bake outputs into the committed file and has
  been removed.
- Check: `tests/unit/test_notebooks.py` rejects any committed output, execution count, or
  run-specific metadata, and smoke-executes the priority notebooks' code cells.

## Status (rewritten + executed error-free, with verified assertions + figure)

- [x] 01 reference frames — frames/domains, FrameTransform apply/inverse/compose, VectorSet guard
- [x] 02 rotations/orientations — parameterizations, geodesic distance, misorientation/disorientation, OrientationSet batch
- [x] 03 symmetry — point/Laue/proper groups, operators, orbits/multiplicity, fundamental sector
- [x] 04 phases/CIF — metric tensor, cell volume, reciprocal d-spacings, atomic basis vs fixture site counts, CIF round-trip
- [x] 05 acquisition/manifests — AcquisitionGeometry, ExperimentManifest JSON round-trip + schema validation
- [x] 06 texture/ODF/PF inversion — kernel ODF, evaluate/volume-fraction, PF/IPF, invert_pole_figures recovers weights
- [x] 07 EBSD grid — CrystalMap, KAM edge detection, grain segmentation, GROD/GOS, IPF+KAM maps
- [x] 08 diffraction geometry — |g|=1/d, electron wavelength/Ewald, Bragg 2theta, kinematic spots + zone law
- [x] 09 phase-transformation foundations — K-S/N-W ORs, describe(), literature misorientation angles, variants, intervariant spectrum, PhaseTransformationRecord

- [x] 10 plotting semantic primitives — Wulff net, directions/planes/symmetry elements, PF/IPF
- [x] 11 powder XRD (priority) — Bragg + FCC selection rules + multiplicity/structure factor + anode shift
- [x] 12 SAED (priority) — zone law + camera constant R=C|g| + forbidden-spot extinction
- [x] 13 crystal visualization (priority) — CrystalScene, plane/direction overlays, Miller-Bravais, zone-axis view
- [x] 14 YAML style customization — theme→file→override precedence, two-theme render
- [x] 15 structure->diffraction pipeline (priority) — one phase → scene+XRD+SAED + workflow/validation manifests
- [x] 16 EBSD -> texture outputs — plane-direction seeding, texture_report ODF/PF/IPF, IPF+KAM maps
- [x] 17 miller vectorized — fully rewritten: d-spacings, cubic angles, families, Miller-Bravais, projection
- [x] 18-21 OR track — already excellent hand-authored; verified green; added house-style cross-links

## STATUS: COMPLETE (all 21 notebooks)

- All 21 notebooks execute error-free and are committed executed (0 errors, 0 unexecuted code cells).
- `tests/unit/test_notebooks.py` passes in full (7 tests), including the priority smoke-execute lane.
- Every rewritten notebook carries: "where this sits" framing, learning goals, theory + citations,
  live-computed results asserted against analytic/literature values, at least one figure (except 05
  manifests and 17 batch-arithmetic, which are inherently non-visual), and forward/back cross-links.
- Fixed: `plt.show()` was removed everywhere — under pytest's warnings-as-errors it raised on the
  Agg backend; the inline backend auto-renders figures at cell end without it.
- Updated `docs/site/tutorials/notebooks.md` (removed stale "notebook generator" reference).

## Constraints learned (API gotchas)

- `Misorientation.disorientation()` and `.angle_deg`/`.angle_rad` — disorientation is a METHOD.
- `misorientation_to(other, reduce_by_symmetry=False)` to show raw vs reduced.
- `Rotation.from_bunge_euler(phi1,Phi,phi2, degrees=True)` (degrees default True).
- `CrystalPlane.d_spacing_angstrom` and `.reciprocal_lattice_vector` are PROPERTIES.
- `Lattice.direct_basis().matrix` (returns Basis object).
- `CrystalMapPhase.name` must equal `phase.name` (e.g. "nickel-fcc"), not an arbitrary label.
- Segmentation `grain_sizes()`/`grain_orientation_spread_deg()` return dicts; `grod_map_deg()` is a method.
- pymatgen emits CIF warnings on the minimal fixtures; suppress with warnings.filterwarnings.
- Priority notebooks (04,11,12,13,15) must contain `get_phase_fixture`/`list_phase_fixtures`.
- notebooks index `docs/site/tutorials/notebooks.md` must mention each notebook stem.
