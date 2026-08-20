# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history. Governed by the cardinal rule in `AGENTS.md`: ledger plus commit-and-push to `main`
after every substantial increment.

The previous task (Crystal Viewer Orientation Workspace, complete 2026-08-19) is archived at
`docs/development/archive/crystal_viewer_orientation_workspace_2026_08.md`.

## EDAX OIM HDF5 (`.oh5` / `.h5`) EBSD Scan Reader — IN PROGRESS (2026-08-20)

**Objective.** Read EDAX OIM HDF5 scans (`.oh5` and `.h5` — the same container, two extensions)
into the same `NormalizedEBSDDataset` that `read_ang` and `read_ctf` produce, and make the
workbench open them alongside `.ang` and `.ctf`.

### Format facts established by inspecting a real OIM 8.6 file

Inspected `kikuchiBandAnalyzer/testData/DA.oh5` and its byte-identical `DA.h5`, together with the
`DA.ang` that OIM exported from the same scan. Layout:

- Top level: `Manufacturer` (`EDAX`), `Version` (`OIM Analysis 8.6...`), and one group per scan.
- `<scan>/EBSD/Header`: `Grid Type` (`SqrGrid`/`HexGrid`), `nColumns`, `nRows`, `Step X`,
  `Step Y`, `Sample Tilt`, `Working Distance`, and a `Phase` subgroup whose members are named
  `1`, `2`, ... Each phase carries `MaterialName`, `Formula`, `LGsymID`, `Point Group`,
  `Laue Group`, and lattice constants.
- `<scan>/EBSD/Data`: one flat length-`nRows*nColumns` array per channel — `Phi1`, `Phi`, `Phi2`
  (Bunge **radians**, same as `.ang`), `X Position`, `Y Position`, `Phase`, `CI`, `IQ`, `Fit`,
  `SEM Signal`, plus vendor and user channels (`PRIAS ...`, `Valid`, ...) and a `Pattern` stack.
- `LGsymID` uses exactly the TSL symmetry codes the `.ang` `# Symmetry` line uses, so
  `_TSL_SYMMETRY_CODES` in `scan_files.py` is reused rather than duplicated (`43` -> `m-3m`).
- The `Phase` channel matched the `.ang` `Phase index` column value for value in `DA`: `0` for the
  points indexed as the single declared phase and `-1` for the unindexed ones. So OIM writes a
  **zero-based** index into the header phase list with `-1` as the unindexed sentinel. Older
  one-based writers (`0` = unindexed) exist, so the reader infers the base and records which it
  chose in the manifest metadata.

### Decisions

- The reader lives beside its siblings in `pytex/adapters/scan_files.py`, returns the same
  `EBSDScanFileResult`, and reuses `_normalize_vendor_payload` with `source_system="edax_oh5"`.
- `h5py` is an **optional** dependency (extra `hdf5`, also in `dev` so the base lane exercises the
  reader). Import is lazy inside the reader, per the anti-goal against import-time coupling to
  optional scientific stacks.
- Point-keeping matches `read_ang` exactly so the same scan imports the same way from either
  format: a single-phase file keeps every row; a multiphase file keeps only rows referencing a
  declared phase.
- Channel names match `read_ang`'s (`image_quality`, `confidence_index`, `detector_signal`,
  `fit`) so a workflow does not care which format it was fed. Every other per-point scalar in the
  file is exposed too, under a normalized name — that is the reason to prefer the HDF5 export.
- A `read_scan(path)` dispatcher picks the reader by extension, so extension dispatch exists once
  instead of in every caller.

### Plan

1. **Library reader.** `read_oh5` + `read_scan` in `pytex/adapters/scan_files.py`, exports,
   `pyproject.toml` extra and mypy override, unit tests. — **DONE**
2. **Workbench.** Binary (base64) upload transport, `.oh5`/`.h5` in the panel's accepted
   suffixes, service dispatch through `read_scan`, tests. — NOT STARTED
3. **Docs.** Workbench workflow page, adapters docs, CHANGELOG. — NOT STARTED

### Verification of increment 1

Checked `read_oh5` against `kikuchiBandAnalyzer/testData/DA.oh5`, its byte-identical `DA.h5`, and
the `DA.ang` OIM exported from the same scan (these files are outside this repository and are not
committed; the unit tests build their own fixtures with `h5py`). Over the 36 points the `.ang`
export carries:

| quantity | max difference, `.ang` vs `.oh5` |
| --- | --- |
| orientation matrix entries | 7.8e-6 |
| map coordinates | 3.8e-7 |
| `confidence_index` | 4.3e-4 |
| `fit` | 4.9e-4 |
| `image_quality` | 4.9e-2 (the `.ang` writes it to one decimal) |
| `detector_signal` | 0 |

The residuals are the text export's rounding, so the two formats agree. `.h5` and `.oh5` gave
bit-identical results. The `.oh5` additionally carries six points the `.ang` export omitted
entirely and the `valid`, `prias_*` channels the text format has no room for, and it does claim
the `(7, 6)` grid the `.ang` cannot (its 36 rows do not fill the 7x6 the header declares).

`ruff check`, `ruff format --check`, `mypy`, and `pytest tests/unit/test_scan_files.py
tests/unit/test_scan_files_oh5.py tests/unit/test_manifests.py tests/unit/test_repo_integrity.py
tests/unit/test_hex_grid_ebsd.py` all green.

### Worktree state

Increment 1 committed. Clean at `bed1b74` when this task started.
