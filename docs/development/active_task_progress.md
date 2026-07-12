# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Objective

Identify the next high-value PyTex development work from the authoritative roadmap, implement a
substantial coherent slice with tests and documentation, verify the repository, then commit and
push the result.

## Current Status

- Started: 2026-07-12
- Branch: `main` (tracking `origin/main`)
- Baseline commit: `cdc5fd53` (`Add interim semantic primitives and canonical visual docs`)
- Phase: rebased and verified; final lint follow-up commit and push pending
- Pre-existing untracked files to preserve and exclude unless intentionally brought into scope:
  `docs/presentations/`, `package.json`, and `references/britton_up_down_ebsd.pdf`

## Completed

- Read the root `AGENTS.md` instructions and identified the mandated architecture, testing,
  governance, notation, documentation, figure, terminology, and contract references.
- Inspected the worktree, recent commits, package configuration, roadmap, and foundation audit.
- Confirmed that the immediate roadmap prioritizes validation and interoperability hardening over
  broad new stable API surface.
- Added the durable-progress and resumability principle to `AGENTS.md`.
- Created this task handoff record.
- Selected diffraction validation hardening after discovering that cubic-only external tests hid
  duplicate powder-family emission.
- Updated `generate_powder_reflections()` to emit one deterministic representative per symmetry
  orbit while retaining the orbit size as multiplicity.
- Added a pinned `pymatgen` Cu K-alpha powder baseline for non-cubic `zr_hcp` and synchronized the
  benchmark, workflow-result, validation, roadmap, theory, and audit documentation surfaces.
- Updated structural plotting validation to expect the five unique labeled Ni reflection families
  plus the continuous spectrum line instead of accepting duplicated family sticks.

## Decisions

- Preserve all pre-existing untracked user files and do not include them in the task commit without
  a clear repository reason.
- Select work from the Phase 1 validation/interoperability program, favoring a bounded slice with
  deterministic evidence, synchronized manifests, documentation, and tests.
- Do not broaden MTEX, external-tool, or physical-model equivalence claims beyond available
  evidence.
- Treat peak position, multiplicity, and family uniqueness as hard claims; keep intensity agreement
  explicitly non-normative.

## Verification Log

- Initial focused run exposed two expected adjustment needs: the family representative policy had
  to preserve conventional positive-index representatives, and the HCP high-angle peak required a
  `0.2 deg` tolerance rather than the cubic-only `0.15 deg` tolerance.
- Focused verification passed after both adjustments.
- `python -m pytest -q tests/unit/test_diffraction_external_baselines.py tests/unit/test_plotting_structural_validation.py`: passed.
- `python -m pytest -q`: passed (full scientific environment; warnings are existing optional-stack
  CIF/spglib deprecations and fixture notices).
- `python scripts/check_repo_integrity.py`: passed.
- `python -m ruff check .`: passed.
- `python -m mypy src`: passed for 58 source files.
- `python -m sphinx -b html docs/site docs/_build/html`: passed.
- Parsed all benchmark and diffraction fixture JSON files successfully (23 files).
- `git diff --check`: passed; only Git's expected LF-to-CRLF worktree notices were emitted.
- Initial task commit was rebased onto 37 newer `origin/main` commits as `aff1dd82`. The only
  conflict was in the newly vectorized XRD enumeration; the resolution preserves vectorized
  geometry and structure-factor evaluation while deduplicating symmetry families first.
- Post-rebase full test suite, strict mypy (71 source files), repository integrity, and Sphinx HTML
  build passed.
- Cleared 12 mechanical Ruff regressions introduced by the incoming commits (sorted export lists,
  redundant integer casts, and unused test tuple members); final lint verification is pending.

## Next Actions

1. Run final Ruff, mypy, integrity, and diff checks after the mechanical lint cleanup.
2. Commit the lint cleanup and this updated handoff separately.
3. Push `main` to `origin` and record the resulting commit identifiers here if another task resumes
   from this handoff.

## Resume Command Checklist

```powershell
git status --short --branch
Get-Content docs/development/active_task_progress.md
git diff -- AGENTS.md docs/development/active_task_progress.md
```
