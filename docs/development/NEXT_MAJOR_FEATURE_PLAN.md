# Next Major Development Feature: Reconstruction Stabilization (2026-07-30)

**Status:** IDENTIFIED AND READY TO IMPLEMENT  
**Tracked in:** Task #1 (long-horizon), `docs/development/active_task_progress.md` (detailed phases)  
**Baseline Commit:** `50116880` (Burgers beta→alpha zirconium tutorial)

---

## Executive Summary

The next major feasible development feature is **Reconstruction Stabilization Phase**, which completes the OR-analysis flagship work by moving parent-grain reconstruction from experimental to production-ready status. This involves:

1. Creating MTEX fixture-generation scripts
2. Integrating measured-data validation into the test suite
3. Establishing MTEX `calcParent2Child` parity baseline
4. Completing the validation-matrix documentation

**Expected Duration:** 2–3 weeks  
**Complexity:** Medium (well-scoped, clear objectives, minimal external dependencies)  
**Impact:** Unblocks parent reconstruction for production EBSD workflows

---

## Why This Feature?

### Roadmap Alignment

The project follows a **prioritized roadmap** split across:

1. **Feature Opportunities Roadmap** (`references/feature_opportunities.md`): Prioritizes breadth across the library
2. **Critical Review & Development Guide** (`docs/roadmap/critical_review_and_development_guide.md`): Defines the governance cycles for OR analysis (designated flagship)
3. **OR Analysis Foundation** (`docs/architecture/orientation_relationship_analysis_foundation.md`): Detailed feature specification F1–F14

### Strategic Position

- **Orientation-relationship (OR) analysis is the designated flagship** per mission.md and AGENTS.md
- **F1–F12 complete** (index correspondence, misorientation, fitting, parallelism, variants, deformation gradients, packets, pole figures)
- **F4 (parent-grain reconstruction) exists but is experimental** — requires measured-data validation to stabilize
- **Natural continuation** after the Figure and Burgers program (which exhaustively documented Burgers OR in both directions)
- **Unblocks downstream work** on F13 (habit planes), reconstruction-based workflows, and broader Cycle C features

### Feasibility

✓ **Well-scoped:** Clear phases with defined success criteria  
✓ **No blockers:** All required infrastructure exists in code  
✓ **Documented:** Plan fully written in active_task_progress.md  
✓ **Fallback ready:** Synthetic-data expansion path if MATLAB unavailable  
✓ **Good sequencing:** Completes flagship before pivoting to texture/diffraction breadth

---

## Plan Overview

### Phase 1: MTEX Fixture Generation Script (1 week)
**Objective:** Create reference data by running MTEX's `calcParent2Child` algorithm

**Deliverables:**
- `scripts/mtex_generators/generate_martensite_fixture.m` (MATLAB script)
- `scripts/mtex_generators/generate_bcc_hcp_fixture.m` (MATLAB script)
- `scripts/mtex_generators/README.md` (documentation)
- Fixture JSON files (austenite parent + KS/NW/GT variants)

**Verification:** Scripts run in MATLAB R2021+, produce valid JSON fixtures with documented parent counts

---

### Phase 2: Measured Data Integration (1 week)
**Objective:** Load MTEX fixtures and validate reconstruction accuracy

**Deliverables:**
- Extended `tests/unit/test_parent_reconstruction.py` with MTEX data
- Pinned MTEX result directories in `fixtures/`
- Baseline accuracy established (e.g., synthetic = exact, MTEX = ±0.5 deg)

**Verification:** MTEX fixtures load, parents recover within tolerance, all tests pass

---

### Phase 3: Validation Ledger & Documentation (1 week)
**Objective:** Complete the validation matrix and update architecture docs

**Deliverables:**
- Complete row in `docs/testing/mtex_parity_matrix.md`
- Updated `docs/architecture/orientation_relationship_analysis_foundation.md` (F4 → stable)
- CHANGELOG entry
- Cross-links from concepts and worked examples

**Verification:** Validation matrix complete, links work, foundation doc updated

---

## Quality Gates (Every Commit)

All changes must pass:

```bash
python -m pytest                          # No warnings, new tests pass
python scripts/check_repo_integrity.py    # Manifest validity
python -m ruff check .                    # Style + imports
python -m mypy src                        # Type checking strict
python -m sphinx -b html docs/site docs/_build/html  # Docs build clean
```

---

## Fallback Plan (If MATLAB Unavailable)

If MTEX access is not available, the alternative is:

### Phase A: Expand Synthetic-Data Validation (2 weeks)
1. Create large synthetic fixtures with 100+ planted parents
2. Test multiple noise levels (0.5, 1.0, 2.0 deg)
3. Document reconstruction robustness vs grain size/noise
4. Pin accuracy baselines for future comparison

### Phase B: Conservative Ledger Update (1 week)
1. Mark reconstruction as "foundationally ready" with synthetic-only validation
2. Flag MTEX parity as a queued follow-on (future work)
3. Document the validation regime and limitations honestly

**This path is complete and publication-ready, just without external baseline comparison.**

---

## What Comes Next (After This Feature)

In roadmap order:

1. **F13 (Habit Planes and PTMC)** — Large feature; deformation-gradient and invariant-line analysis. Requires dedicated working notes.

2. **Finding 8 remainder (Texture Kernels)** — Add bump/fibre kernels, integrate SO3FunHarmonic for full harmonic ODF estimation.

3. **Finding 9 (PF→ODF Improvements)** — Ghost correction and zero-range handling in harmonic inversion.

4. **Finding 14 (Kikuchi Geometry)** — Gnomonic projection and Kikuchi-band overlays for EBSD pattern indexing.

5. **Human-auditable test documentation (Meta)** — Sphinx pages documenting important tested methods with formulas, sources, examples, and verified outputs. Applies across all subsystems.

---

## How to Resume

If work is interrupted:

1. **Read this file** — Current status and full plan
2. **Check `docs/development/active_task_progress.md`** — Detailed phase ledger and next actions
3. **Check Task #1** — Status and any blockers
4. **Check latest commit** — Understand current worktree state
5. **Run gates** — Verify all quality checks pass before resuming

Each phase is independently landable as a commit+push, so interruption costs at most the current phase.

---

## Key Files

| File | Purpose |
|------|---------|
| `docs/development/active_task_progress.md` | **Read first to resume** — detailed phase plan, ledger, decisions |
| `docs/architecture/orientation_relationship_analysis_foundation.md` | OR feature F4 definition and current status |
| `docs/roadmap/critical_review_and_development_guide.md` | Governance; Cycle C scope and quality bars |
| `mission.md` | Project vision; OR is designated flagship |
| `AGENTS.md` | Working process; resumability principle |
| `references/feature_opportunities.md` | Broader roadmap; context for why OR comes before other features |

---

## Decision Log

**Why OR reconstruction over other candidates?**

- ✓ Completes the flagship work (F1-F12 already done)
- ✓ Unblocks production EBSD workflows
- ✓ Establishes MTEX validation culture for future subsystems
- ✗ Broader features (texture kernels, Kikuchi, PF→ODF) can start after
- ✗ Broader roadmap item #1 (test documentation) is meta; applies after features land

**Why this scoping (not "everything to do" or "one tiny thing")?**

- ✓ Long-horizon task = multiple phases, multi-week effort
- ✓ Well-defined success criteria and acceptance gates
- ✓ Resumable at phase boundaries
- ✓ Aligns with project's explainable-results doctrine and validation culture

---

## Success Criteria

When this feature is complete:

- [ ] MTEX fixtures generated and fixtures/ directory populated
- [ ] Tests load MTEX data and assert parent recovery accuracy
- [ ] Validation-matrix row complete with citations
- [ ] Foundation doc updated; F4 marked as stable
- [ ] CHANGELOG entry recorded
- [ ] All quality gates pass (ruff, mypy, pytest, Sphinx)
- [ ] No coverage decrease
- [ ] No new warnings

At that point, parent-grain reconstruction is production-ready and the OR flagship is 95% complete (only F13 habit planes and F14 cementite ORs remain).
