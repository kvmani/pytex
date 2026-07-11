# Working Notes: Immediate-Roadmap Sprint (2026-07-11)

Purpose: resumable progress ledger for the current fast implementation sprint.
If the session is interrupted, resume by reading this file top to bottom; each
task lists status, files, and exactly what remains.

Session context: continues the session that landed, in order:

1. `aacb50a` MTEX parity roadmap + S2 spherical vector layer (`core/sphere.py`)
2. `cec93f7` full 32-point-group model (`core/point_groups.py`), specimen symmetry
3. `b6d5e61` exact Laue fundamental-sector geometry + IPF key legends
4. `b4520ba` direct `.ang` / `.ctf` EBSD scan readers (`adapters/scan_files.py`)

Master plan: `docs/roadmap/mtex_parity_and_ebsd_feature_roadmap.md` (immediate
horizon, Section 3).

## Sprint Task List (chosen for ~1-2 hours total)

| # | Task | Status | Commit |
| --- | --- | --- | --- |
| 1 | Working-notes file (this file) | done | pending |
| 2 | Symmetry-aware `OrientationSet.mean_orientation()` + spread | todo | - |
| 3 | Named texture components registry + volume fractions (`texture/components.py`) | todo | - |
| 4 | `Fibre` class with named bcc/fcc fibres + fibre distances (`texture/fibres.py`) | todo | - |
| 5 | de la Vallee Poussin kernel with halfwidth/kappa/Chebyshev coefficients (`texture/kernels.py`) | todo | - |
| 6 | Parity-matrix + roadmap status updates, final notes update | todo | - |

Workflow per task: implement -> unit tests -> `ruff check .` + `mypy src` +
targeted pytest -> commit -> push -> update this table.

## Task Details And Resume Instructions

### Task 2: Orientation mean and spread

- File: `src/pytex/core/orientation.py`, class `OrientationSet`.
- Method `mean_orientation()`: symmetry-aware quaternion mean.
  Algorithm: take first orientation's quaternion as reference; for every
  orientation choose the symmetry-equivalent quaternion (right-multiply by
  crystal symmetry operator quaternions, sign-canonicalized) with maximal
  |dot| to the reference; mean = principal eigenvector of the 4x4 outer-product
  sum (Markley method); iterate selection against the new mean until stable
  (max ~10 iterations). Return `Orientation` with the set's frames/symmetry/
  phase.
- Method `spread_angles_deg(reference=None)`: misorientation angles to the
  mean (or given reference) reusing `misorientation_angles_to`.
- Tests in `tests/unit/test_orientation_statistics.py`: cluster around a known
  orientation with symmetry-scattered representatives -> mean recovers it;
  spread near zero for identical orientations; works without symmetry too.

### Task 3: Named texture components

- New file: `src/pytex/texture/components.py`.
- `TextureComponent` frozen dataclass: name, bunge_euler_deg (3-tuple),
  miller_label (e.g. "{011}<100>"), notes. Method
  `orientation(crystal_frame, specimen_frame, symmetry=None, phase=None)`.
- Registry `STANDARD_FCC_ROLLING_COMPONENTS`: cube (0,0,0), rotated_cube
  (45,0,0), goss (0,45,0), brass (35,45,0), copper (90,35,45), s (59,37,63).
  BCC starter registry `STANDARD_BCC_ROLLING_COMPONENTS`: rotated_cube,
  rotated_goss (90,90,45), inverse_brass? keep minimal: rotated_cube plus the
  named fibres cover bcc; document that.
- `component_volume_fractions(orientations, components, tolerance_deg=15,
  weights=None)` -> dict name->fraction using symmetry-aware
  `misorientation_angles_to`.
- Export from `pytex.texture` and top-level `pytex` as appropriate.

### Task 4: Fibre class

- New file: `src/pytex/texture/fibres.py`.
- `Fibre` frozen dataclass: `crystal_direction` (unit vector in crystal
  frame or MillerDirection), `specimen_direction` (str like "ND" or vector),
  name. Methods:
  - `orientations(count, crystal_frame, specimen_frame, symmetry)` -> sample
    the fibre: base rotation aligning crystal dir to specimen dir composed
    with rotations about the specimen dir on a regular grid.
  - `angles_to_deg(orientation_set)`: per-orientation fibre distance = min
    over the symmetry family of the crystal direction of the angle between
    the mapped direction and the specimen direction (antipodal-aware).
  - `volume_fraction(orientation_set, tolerance_deg=10, weights=None)`.
- Named constructors (bcc rolling): `Fibre.alpha_bcc(...)` <110> || RD,
  `Fibre.gamma_bcc(...)` <111> || ND, `Fibre.eta(...)` <100> || RD.
- Tests: sampled fibre orientations have zero fibre distance; perturbed
  orientations measure expected angle; gamma-fibre volume fraction of a
  gamma-sampled set is 1.

### Task 5: de la Vallee Poussin kernel

- New file: `src/pytex/texture/kernels.py`.
- `DeLaValleePoussinKernel` frozen dataclass constructed from
  `halfwidth_deg`; derived `kappa = ln(0.5) / (2 ln cos(halfwidth/2))`.
  Methods: `evaluate(omega_rad)` = C * cos^(2 kappa)(omega/2) with
  normalization C = B(3/2, 1/2) / B(3/2, kappa + 1/2) (scipy-free: use
  math.gamma); `chebyshev_coefficients(bandwidth)` by deterministic
  quadrature A_l = integral psi(omega) chi_l(omega) dmu with
  chi_l = sin((2l+1) omega/2)/sin(omega/2), dmu = (2/pi) sin^2(omega/2) domega;
  `bandwidth(threshold=1e-3)` = first l where |A_l| < threshold.
- Tests: A_0 == 1 (normalization), halfwidth relation psi(hw)/psi(0) == 0.5,
  coefficients decay monotonically-ish, smaller halfwidth -> larger bandwidth.
- Wiring into ODF evaluation/`KernelSpec` is deliberately deferred to the
  SO3Fun refactor (medium term); note this in the parity matrix if edited.

### Task 6: Docs sync

- Update `docs/testing/mtex_parity_matrix.md` ODF row (kernel + components)
  if claims change; keep claims conservative (implemented, not MTEX-parity).
- Mark this file's table rows done with commit hashes; keep for the next
  session; it also seeds the next sprint (property channels on CrystalMap,
  IPF-X/Y/Z map call, phase map plotting are the following candidates).

## Known Environment Notes

- Pre-existing failures unrelated to this work: phase-fixture hash mismatches
  (`tests/unit/test_phase_fixtures.py`, `tests/unit/test_repo_integrity.py`);
  a task chip exists to fix them. Deselect with:
  `python -m pytest -q --deselect tests/unit/test_phase_fixtures.py --deselect tests/unit/test_repo_integrity.py`
- Gates: `python -m ruff check .`, `python -m mypy src`, pytest as above; all
  green as of `b4520ba`.
