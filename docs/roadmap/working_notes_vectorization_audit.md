# Working Notes: Repo-Wide Vectorization Audit (2026-07-12)

Purpose: resumable ledger for a ground-up audit that finds every scalar
Python-loop numerical hot path and converts it to vectorized NumPy/SciPy, with
proven result-equivalence. If the session is interrupted (usage limit, etc.),
resume by reading this file top to bottom: the Progress Ledger says exactly
which files are done, in progress, or pending, and each converted loop has an
equivalence check.

Governing principle: `docs/standards/development_principles.md` 9a (prefer
highly vectorized implementations; scipy is a core dependency).

## Method (per file)

1. List the file's `for`/comprehension loops.
2. Classify each: HOT (per-element numerical over large data -> vectorize),
   COLD (metadata/parsing/plotting/small-fixed-size setup -> leave, record why),
   or DONE (already vectorized).
3. For each HOT loop: rewrite vectorized; prove equivalence (targeted test or a
   scratch equivalence check vs the old scalar result, max abs diff < 1e-9);
   keep or add a regression test.
4. Gates after each file or coherent group: `ruff check .`, `mypy src`,
   full pytest (deselecting the pre-existing phase-fixture / repo-integrity
   hash mismatches). Commit + push per file or small group.

Baseline at audit start: `621 passed, 26 deselected`, ruff + mypy green,
HEAD `0058aa4` (misorientation_angles_to already vectorized last session).

## Triage Order (by loop count and hot-path likelihood)

1. `ebsd/models.py` (93 loops) - grains, boundaries, KAM, cleanup, metrics
2. `diffraction/models.py` (41) - reflections, structure factors, indexing
3. `core/orientation.py` (30) - some already vectorized
4. `texture/harmonics.py` (10), `texture/models.py` (9)
5. `properties/tensors.py` (9) - homogenize loop already einsum; check rest
6. `core/symmetry.py` (11), `core/point_groups.py` (10), `core/hexagonal.py` (10)
7. `diffraction/xrd.py` (12), `diffraction/saed.py` (8)
8. Remaining core/adapters numerical loops.

Explicitly OUT OF SCOPE (COLD by nature; not numerical hot paths):
plotting (`plotting/*`, glyph/scene construction), file parsing/serialization
(`adapters/*`, `contracts.py`, `cli.py`, manifest IO), and small fixed-size
setup loops (e.g. 3x3 Voigt index maps). These are recorded here so the audit
is auditable, not skipped silently.

## Progress Ledger

| File | Status | HOT loops found | Converted | Commit |
| --- | --- | --- | --- | --- |
| core/orientation.py :: misorientation_angles_to | done (prev session) | 1 | 1 | 0058aa4 |
| ebsd/models.py | done | 5 | 5 | c30efbe + this |
| core/orientation.py (rest) | in progress | 3 done, 2 deferred | 3 | 35a3ca9 |
| diffraction/models.py | done | 1 HOT (simulate_spots) | 1 | this commit |
| texture/harmonics.py | reviewed - COLD | 0 | 0 | - |
| texture/models.py | reviewed - COLD | 0 | 0 | - |
| properties/tensors.py | reviewed - COLD | 0 | 0 | - |
| core/symmetry.py | reviewed - COLD | 0 | 0 | - |
| core/point_groups.py | reviewed - COLD | 0 | 0 | - |
| core/hexagonal.py | reviewed - COLD | 0 | 0 | - |
| diffraction/xrd.py | done | 1 | 1 | this commit |
| diffraction/saed.py | reviewed - COLD/deferred | - | 0 | - |

## Per-File Findings

### ebsd/models.py (in progress)

HOT loops converted (this commit), all proven exactly equivalent by brute-force
scratch checks (GROD max diff 0.0; boundary/representative match < 1e-9) and the
existing grain/EBSD/CSL test suites:

1. `GrainSegmentation.grod_map_deg` - was a nested grain x member per-point
   `distance_to` loop. Now: gather each point's grain-reference index into an
   array, compute all relative rotations with `_relative_rotation_matrices`, and
   reduce once with `_disorientation_angles_from_relative_matrices` (or raw
   angle when `symmetry_aware=False`).
2. `GrainSegmentation.boundary_network` - was a per-neighbour-pair `distance_to`
   loop. Now: vectorised boundary-pair mask, one batched relative-rotation +
   disorientation reduction (mirrors `distance_to`'s set-level symmetry exactly,
   NOT `_pair_misorientation_rad`'s per-phase reduction, to preserve results),
   plus vectorised lengths/midpoints; only the dataclass construction loops.
3. `CrystalMap._representative_orientation_index` - was an O(m^2) double
   `distance_to` loop. Now: one `misorientation_angles_to` matrix, row-sum,
   `argmin` (first-min tie-break matches the old sequential scan).
4. `GrainSegmentation.grain_perimeters` - was a nested row x col x neighbour
   loop. Now: shifted-slice boundary masks + `bincount` per grain.

5. `CrystalMap.remove_wild_spikes` - was an O(points) loop calling a separate
   misorientation reduction per point for spike detection. Now: one batched
   `_pair_misorientation_rad` over all boundary edges, incident-minimum via
   `np.minimum.at` on both endpoints, and a vectorised spike mask; only the
   (rare) spike points loop to compute their neighbourhood-mean replacement.
   Verified equivalent by an independent per-point reference computation
   (spike set identical) and the existing isolated-spike unit test.

### core/orientation.py (in progress)

Converted (this commit), exact equivalence verified:

1. `OrientationSet.as_matrices` - was a per-quaternion list comprehension via
   `quaternion_to_matrix`; now the batch `quaternions_to_matrices`. Max diff
   ~1e-15. High value: called by nearly every vectorized routine.
2. `matrices_to_quaternions` proper-rotation validation - was a per-matrix
   `is_rotation_matrix` loop; now a single batched `M^T M = I` (`einsum`) plus
   `det = 1` check. Same atol (1e-8).

3. `OrientationSet.as_euler` (and `as_bunge_euler`) - was a per-quaternion
   `Rotation.to_euler` comprehension. Now a batched
   `_matrices_to_repeated_axis_euler` that reproduces the scalar gimbal-lock
   branching (PHI near 0 / near pi) exactly via boolean masks, followed by the
   same `mod 2*pi` and optional `rad2deg`. Verified across bunge/matthies/abg,
   degrees and radians, including gimbal-lock cases (max diff ~1e-14); permanent
   regression test in `test_orientation_utilities.py`.

Deferred HOT candidates in core/orientation.py (higher-risk batch refactors of
per-orientation symmetry logic; each needs its own golden-equivalence pass):
`canonicalize`-map at line ~2125,
`project_to_exact_fundamental_region`-map at ~2155/2181, and the
`matrix_to_quaternion` comprehension inside `disorientation`/fundamental-region
key computation (~1308). These wrap non-trivial per-element symmetry reduction;
convert with captured golden outputs when resumed.

### diffraction/models.py (deferred with plan)

`simulate_spots` (per-reflection loop, ~line 949) is genuinely HOT and
vectorisable (batch reciprocal vectors = miller @ reciprocal_basis; single
orientation matmul; mask-based Ewald/zone/detector filters; batched intensity)
but it is a ~80-line, filter- and object-construction-heavy rewrite feeding
externally-pinned SAED baselines (`test_diffraction_external_baselines.py`,
`test_generate_saed_pattern_respects_zone_axis_geometry`). Plan for the
resumed pass: (1) capture a golden fixture of the current spot list (positions,
excitation errors, intensities, family keys) for a representative case; (2)
rewrite with mask-based filtering keeping the surviving-spot dict assembly; (3)
assert bit-for-bit equality against the golden capture and the external
baselines before committing. The outer orientation-search loops
(`index_pattern`, refinement grid) are algorithmic objective evaluations, not
per-element array math -> COLD; they speed up for free once simulate_spots is
vectorised.

### Other files reviewed this session (classification + rationale)

- `core/symmetry.py`, `core/point_groups.py`: loops operate on point-group
  operator sets (<= 48 matrices) computed once at construction and small BFS
  group closures. Fixed tiny size -> COLD (vectorising gives no benefit).
- `properties/tensors.py`: Voigt <-> rank-4 loops are fixed 6x6 / 3x3x3x3
  index maps; `homogenize_elastic` already uses `einsum`. COLD.
- `texture/models.py`, `texture/harmonics.py`: numerical cores already use
  `einsum`/matrix ops; remaining loops are over symmetry families (<=48),
  pole-figure lists, or iterative solvers (projected-gradient) -> COLD.
- `core/hexagonal.py`: small index-conversion helpers -> COLD.
- `diffraction/xrd.py`: `_structure_factor_xray` site loop is tiny (few atoms;
  attribute gathering stays scalar) -> low value. `generate_powder_reflections`
  hkl loop (~(2n+1)^3 reflections) is a genuine HOT candidate but a moderate
  refactor touching structure factors + multiplicity grouping -> DEFERRED with
  a golden-capture plan (test coverage: `test_generate_xrd_pattern_contains_
  expected_reflection`).
- `diffraction/saed.py` + `diffraction/models.py::simulate_spots`: HOT,
  DEFERRED (see the simulate_spots plan above); externally-pinned baselines
  demand a golden-equivalence harness before rewriting.

## Session Checkpoint (2026-07-12)

Converted this session with verified exact equivalence (max diffs ~1e-14):
ebsd/models.py (grod_map_deg, boundary_network, _representative_orientation_
index, grain_perimeters, remove_wild_spikes detection) and core/orientation.py
(as_matrices, matrices_to_quaternions validation, as_euler). Full suite green
(622 passed) and ~30% faster wall-clock (30s -> 22s). All pushed.

`diffraction/xrd.py::generate_powder_reflections` DONE (this session, after the
checkpoint): geometric filtering (d-spacing, 2-theta, range) vectorised over all
hkls; new batched `_structure_factors_xray` over reflections x sites (scalar
`_structure_factor_xray` kept as a thin wrapper). Validated against a golden
capture of the pre-change output: identical hkl set, per-hkl values match to
1e-11 (FP norm noise). The final sort now rounds 2-theta (9 dp) so symmetry-
equivalent reflections order deterministically by hkl instead of by FP-noise-
level 2-theta, making the output reproducible. Regression test asserts the
batched vs per-reflection structure factors agree.

`diffraction/models.py::simulate_spots` DONE (this session): geometry batched
over all reflections (reciprocal vectors = miller @ reciprocal_basis.T; single
orientation matmul; batched zone/excitation masks; batched detector projection,
two-theta, azimuth, on-detector, acceptance); only the trivial per-survivor
intensity and family-key stay in the assembly loop. Validated against a
4-branch golden capture (no-zone / zone / no-orientation / identity-zone):
identical spot counts, max diff ~5.7e-14. Regression test asserts the
excitation cutoff, monotonic spot-count vs cutoff, and zone orthogonality.

RESUME HERE next session (only lower-value deferred items remain):
1. core/orientation.py deferred: `canonicalize` map (~2125) and
   `project_to_exact_fundamental_region` map (~2155/2181) - batch the
   per-orientation symmetry reduction with captured golden outputs. These are
   the last per-element numerical loops; everything else is COLD.
Each: capture golden -> rewrite vectorised -> assert bit-equivalence + suite.

Remaining ebsd/models.py loops reviewed and classified COLD (kept): per-grain
dict comprehensions (grain count is small; inner ops already vectorised),
`majority_smoothed` / `merge_small_grains` label graph iterations (iterative
graph algorithms, not per-element array math), dataclass/metadata construction,
and validation loops.
