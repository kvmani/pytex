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
| ebsd/models.py | in progress | 4 found so far | 4 | this commit |
| diffraction/models.py | pending | - | - | - |
| core/orientation.py (rest) | pending | - | - | - |
| texture/harmonics.py | pending | - | - | - |
| texture/models.py | pending | - | - | - |
| properties/tensors.py | pending | - | - | - |
| core/symmetry.py | pending | - | - | - |
| core/point_groups.py | pending | - | - | - |
| core/hexagonal.py | pending | - | - | - |
| diffraction/xrd.py | pending | - | - | - |
| diffraction/saed.py | pending | - | - | - |

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

Remaining ebsd/models.py loops reviewed and classified COLD (kept): per-grain
dict comprehensions (grain count is small; inner ops already vectorised),
`majority_smoothed` / `merge_small_grains` label graph iterations (iterative
graph algorithms, not per-element array math), dataclass/metadata construction,
and validation loops. `remove_wild_spikes` still has a per-point adjacency loop
building neighbour means -- flagged for a later pass (candidate HOT).
