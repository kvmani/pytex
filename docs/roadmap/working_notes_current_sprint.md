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
pre-existing phase-fixture/repo-integrity hash mismatches). Current state after
tasks A-L below: `607 passed, 26 deselected`, ruff + mypy green. Gates:
`python -m ruff check .`, `python -m mypy src`, and
`python -m pytest -q --deselect tests/unit/test_phase_fixtures.py --deselect tests/unit/test_repo_integrity.py`.

This became a single long-horizon session spanning three groups of work:
immediate-horizon EBSD (Tasks A-E), medium-horizon grains/boundaries/MDF
(Tasks F, G, I, J), and a new physical-properties layer (Tasks H, K, L). Each
task landed as its own commit with tests, docs, and parity-matrix updates.
New modules this session: `ebsd/csl.py`, `core/misorientation_distribution.py`,
and the `pytex/properties/` package (`slip.py`, `tensors.py`).

## Sprint Task List

| # | Task | Status | Commit |
| --- | --- | --- | --- |
| A | Per-point property channels on `CrystalMap` + reader wiring | done, tests green | committed |
| B | `plot_property_map`, `plot_phase_map`, IPF-XYZ triptych | done, tests green | committed |
| C | Grain scalar metrics: mean orientation, GOS, GAM, eq. diameter | done, tests green | committed |
| D | Cleanup filters: wild-spike removal, property thresholding | done, tests green | committed |
| E | Wire dVP kernel into discrete ODF + docs/parity sync | done, tests green | this commit |

Sprint outcome: all five tasks landed. Full suite (`577 passed, 26 deselected`),
`ruff check .`, and `mypy src` green; only the pre-existing phase-fixture /
repo-integrity hash mismatches remain deselected (tracked separately).

## Continuation Sprint: Medium Horizon (Section 4), same session

| # | Task | Status | Commit |
| --- | --- | --- | --- |
| F | Moment-based grain shape descriptors (`FittedEllipse`, aspect ratio, bbox) | done, tests green | committed |
| G | CSL/twin boundary classification (`ebsd/csl.py` + network wiring) | done, tests green | committed |
| H | Slip systems + Schmid-factor maps (new `properties/` package) | done, tests green | committed |
| I | Twin/CSL grain merging into parent grains (`merge_by_csl`/`twin_merge`) | done, tests green | committed |
| J | Misorientation distribution (MDF) + random baseline (`core/misorientation_distribution.py`) | done, tests green | this commit |

Note (J): landing the MDF surfaced and fixed a latent batch-axis reshape bug
in the vectorized symmetry reduction that also affected `ebsd/csl.py`
(`classify_misorientations`). The bug was masked everywhere because every prior
CSL test used a single row (S=1); a multi-row regression test now guards it in
`test_csl_classification.py::test_classify_batch_preserves_per_row_assignment`.
The reduction einsums now place the sample axis first (`saik`/`sabik`).

Full suite after this continuation: `585 passed, 26 deselected`; ruff + mypy
green. `FittedEllipse` lives in `ebsd/models.py`; the CSL registry, Brandon
criterion, and matrix classifier live in the new `ebsd/csl.py`, with
`GrainBoundaryNetwork.classify_csl` / `csl_fraction` / `select_csl` as the
consumer surface. CSL matching uses the symmetry-reduced deviation
$\min_{S_a, S_b} \angle\!\left(S_a \mathbf{M} S_b \mathbf{C}^{\mathsf{T}}\right)$ (cubic-only for now).

| K | Elastic tensor layer (`properties/tensors.py`: stiffness/compliance, E(n)) | done, tests green | committed |
| L | Polycrystal elastic homogenization (`homogenize_elastic`, Voigt/Reuss/Hill) | done, tests green | committed |
| M | Grain perimeter / area / shape factor (staircase, regular grid) | done, tests green | committed |
| N | Texture index + entropy on `HarmonicODF` (SO(3) quadrature) | done, tests green | committed |
| O | scipy adopted as core dep + vectorization principle (9a) | done | committed |
| P | Full-constraint Taylor factor (`properties/taylor.py`, scipy LP) | done, tests green | committed |
| Q | Vectorize `OrientationSet.misorientation_angles_to` (principle 9a) | done, tests green | this commit |

Note (Q): the core pairwise-disorientation method was an O(n*m) Python double
loop constructing `Orientation` objects per pair -- the hottest scalar path,
used by ODF evaluation and orientation spread. Replaced with a fully vectorized
`einsum` symmetry reduction (`_reduced_pair_disorientation_angles`), processed
in memory-bounded blocks. Bit-identical to the old result (max diff ~1e-14);
an 80x80 matrix drops from tens of seconds to ~0.14 s. Equivalence regression
test in `test_misorientation_angles_vectorized.py`.

Note (O/P): scipy is now a first-class core dependency (`pyproject.toml`), with
principle 9a in `docs/standards/development_principles.md` mandating highly
vectorized implementations. The Taylor model vectorizes all Schmid-tensor setup
via `einsum` and solves the per-orientation minimum-slip LP with SciPy HiGHS
(`scipy.optimize.linprog`). mypy needs `scipy.*` in the ignore-missing-imports
overrides. Validated: cube orientation = sqrt(6) exactly, random fcc mean ~3.06.

Note (K): slotted frozen dataclass subclasses cannot use zero-arg `super()`
in `__post_init__` (dataclass(slots=True) rebuilds the class, breaking the
`__class__` cell); call the base `__post_init__` explicitly instead. Young's
modulus uses the tensor contraction `1/E(n) = n_i n_j n_k n_l S_ijkl`.

Remaining medium-horizon candidates, in rough priority order:
polyline grain-boundary geometry (true perimeter, convexity/paris shape factor,
`smooth()`) as the rest of Grains 2.0; hex-grid `CrystalMap` support (offset-row
indexing + honeycomb neighbors); twin-merging into parent grains built on the
new CSL layer; `MisorientationDistribution` (MDF) with MacKenzie baseline and
axis/angle marginals; and the `SO3FunHarmonic` quadrature-from-orientations
implementation that lets the dVP kernel drive full harmonic ODF estimation (the
piece deliberately staged in the immediate sprint). `SlipSystem` +
Schmid-factor maps are the highest-demand properties entry point once symmetry
completion (already done) is leveraged.

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
  ($d = 2\sqrt{A/\pi}$, with $A = \text{size} \times \mathrm{d}x \times \mathrm{d}y$); requires step_sizes.
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
