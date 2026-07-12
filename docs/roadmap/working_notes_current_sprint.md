# Working Notes: EBSD Immediate-Horizon Sprint (2026-07-12)

Purpose: resumable progress ledger for the current long-horizon implementation
sprint toward MTEX feature parity. If the session is interrupted, resume by
reading this file top to bottom; each task lists status, files, and exactly
what remains.

Previous sprint (2026-07-11, commit `80c05a5`) landed orientation statistics,
named texture components, fibres, and the de la Vallee Poussin kernel. Its
notes are preserved in git history (`f7c58f5`).

Master plan: `docs/roadmap/mtex_parity_and_ebsd_feature_roadmap.md`, Section 3.2
(EBSD features, immediate horizon). This sprint works that list end to end.

Baseline at sprint start: `558 passed, 26 deselected` (deselecting the
pre-existing phase-fixture/repo-integrity hash mismatches). Gates:
`python -m ruff check .`, `python -m mypy src`, and
`python -m pytest -q --deselect tests/unit/test_phase_fixtures.py --deselect tests/unit/test_repo_integrity.py`.

## Sprint Task List

| # | Task | Status | Commit |
| --- | --- | --- | --- |
| A | Per-point property channels on `CrystalMap` + reader wiring | pending | — |
| B | `plot_property_map`, `plot_phase_map`, IPF-XYZ triptych | pending | — |
| C | Grain scalar metrics: mean orientation, GOS, GAM, eq. diameter | pending | — |
| D | Cleanup filters: wild-spike removal, property thresholding | pending | — |
| E | Wire dVP kernel into discrete ODF + docs/parity sync | pending | — |

Workflow per task: implement -> unit tests -> `ruff check .` + `mypy src` +
targeted pytest -> commit -> push if remote available -> update this table.

## Task Details And Resume Instructions

### Task A: CrystalMap property channels

- File: `src/pytex/ebsd/models.py`, class `CrystalMap`.
- Add field `properties: Mapping[str, np.ndarray] | None = None` (last field so
  existing positional construction is unaffected). Normalize in `__post_init__`
  to a read-only `MappingProxyType` of contiguous float64 arrays, each of shape
  `(n_points,)`; reject wrong lengths and non-string keys.
- Accessors: `property_names` (tuple), `get_property(name)` (1-D array),
  `property_map(name)` (reshaped to `grid_shape`, requires regular 2-D grid),
  `with_properties(mapping)` (returns a new CrystalMap merging/replacing
  channels).
- Propagate `properties` through `select_phase` (mask each channel) and through
  any new subset constructor added in Task D.
- Reader wiring: `read_ang`/`read_ctf` in `adapters/scan_files.py` already parse
  property arrays into `EBSDScanFileResult.properties`. Attach them to the
  produced `crystal_map` (rebuild the dataset's crystal_map via
  `with_properties`). Keep `EBSDScanFileResult.properties` too for back-compat.
- Tests: `tests/unit/test_crystal_map_properties.py`; also extend the scan-file
  reader test to assert channels reach `result.crystal_map`.

### Task B: Property/phase map plotting

- File: `src/pytex/plotting/ebsd.py`.
- `plot_property_map(crystal_map, name, *, cmap, ...)`: imshow the reshaped
  channel with a labelled colorbar; scatter fallback for graph mode; boundary
  overlay support mirroring `plot_kam_map`.
- `plot_phase_map(crystal_map, ...)`: discrete colors per phase id from
  `phase_id_array`, legend of phase names, ListedColormap.
- `plot_ipf_xyz_maps(crystal_map, ...)`: 1x3 panel of IPF-X/Y/Z maps reusing
  `plot_ipf_map` on shared axes.
- Register cases in `plotting/_plotting_validation_cases.py` so the structural
  validation test covers them; export new functions from `plotting/__init__`.

### Task C: Grain scalar metrics

- File: `src/pytex/ebsd/models.py`, class `GrainSegmentation`.
- `grain_mean_orientation(grain)` / `grain_mean_orientations()`: symmetry-aware
  `OrientationSet.mean_orientation()` over each grain's member orientations.
- `grain_orientation_spread_deg()` (GOS): per grain, mean of
  `spread_angles_deg` of members about the grain mean; broadcast to a per-point
  `gos_map_deg()`.
- `grain_average_misorientation_deg()` (GAM): per grain, mean intragranular KAM
  (reuse `kernel_average_misorientation_deg(segmentation=self)` restricted to
  in-grain pairs, averaged per grain); `gam_map_deg()` broadcast.
- `grain_equivalent_diameters()`: from member pixel count and step area
  (`d = 2 sqrt(area / pi)`, area = size * dx * dy); requires step_sizes.
- Return dicts keyed by grain_id; map broadcasts return grid arrays.
- Tests: single crystal -> GOS/GAM ~0; a two-grain synthetic map -> correct
  per-grain means and diameters.

### Task D: Cleanup filters

- File: `src/pytex/ebsd/models.py`, class `CrystalMap`.
- `remove_wild_spikes(*, threshold_deg, symmetry_aware=True, connectivity=8)`:
  for each point whose min disorientation to every neighbor exceeds
  threshold_deg (an isolated spike), replace its orientation with the neighbor
  that is most representative (neighbor whose orientation is closest to the
  neighbor mean). Returns a new CrystalMap on the same grid, properties carried.
- `property_threshold_mask(name, *, minimum=None, maximum=None)` -> bool array.
- `select_points(mask)` -> subset CrystalMap (graph mode, properties masked),
  used for CI/IQ threshold filtering; shares logic with `select_phase`.
- Tests: a planted spike is corrected; threshold mask/selection counts correct.

### Task E: dVP kernel wiring + docs

- Wire `DeLaValleePoussinKernel` (in `texture/kernels.py`) into the discrete ODF
  evaluation path so `KernelSpec`/`ODF` can use it as the smoothing kernel.
  Inspect `texture/models.py` `ODF`/`KernelSpec` first; keep the existing
  default kernel behavior unless dVP is explicitly requested.
- Update `docs/testing/mtex_parity_matrix.md` EBSD + ODF rows conservatively
  (implemented, not yet MTEX-pinned), mark this table's rows done with commit
  hashes, and note landed features in the roadmap if claims change.

## Known Environment Notes

- Pre-existing failures unrelated to this work: phase-fixture hash mismatches
  (`tests/unit/test_phase_fixtures.py`, `tests/unit/test_repo_integrity.py`);
  a task chip exists to fix them. Always deselect them when gating.
- `CrystalMap` is a frozen, slotted dataclass; construct derived maps with a
  full constructor call (as `select_phase` does) — remember to thread the new
  `properties` field through every such call.
- No pinned MTEX result files exist yet in `fixtures/mtex_parity/results/mtex/`,
  so new metrics are validated against analytic/self-consistent expectations,
  not external pins, this sprint.
