# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Current Goal (2026-07-30)

**Reconstruction Stabilization Phase** — Make parent-grain reconstruction production-ready by
validating against measured EBSD data and establishing MTEX parity.

This task continues the OR-analysis flagship program after the completion of F1-F12 (index correspondence, misorientation, fitting, parallelism, variant enumeration, deformation gradients, packet classification, and variant pole figures). Phase 4 of Cycle C.

**Scope:**
1. Acquire or generate measured EBSD dataset with known parent/child relationship (e.g., martensite microstructure)
2. Establish MTEX `calcParent2Child` parity baseline
3. Validate reconstruction algorithm against both synthetic and measured data
4. Complete validation-matrix entry and documentation

**Success Criteria:**
- Reconstruction recovers known parents from EBSD maps with <1.5 deg accuracy on measured data
- MTEX parity established via pinned comparison fixtures
- Validation ledger complete and linked from documentation
- All tests pass, ruff/mypy/integrity/Sphinx green

---

## MAJOR FINDING (2026-07-30) — re-scoped Phase 1

**MATLAB/MTEX is not available on this machine** (`which matlab` → not found), so the
measured-data parity path cannot be executed here. While auditing reconstruction to plan the
fallback, a **real scientific defect** was found and quantified, and it is the actual blocker to
stabilization:

**`reconstruct_parent_grains` links adjacency edges using the misorientation ANGLE only.**
`_reference_angles_deg` reduces the intervariant table to its ~11 distinct *angles*; an edge links
if the boundary disorientation angle is within tolerance of any of them. The **axis is discarded.**

Measured (20 000 random unrelated boundaries, cubic-cubic):

| tolerance | angle-only false-accept | axis-aware false-accept | reduction |
| --- | --- | --- | --- |
| 1.0 deg | 28.57 % | 0.26 % | 112x |
| 2.0 deg | 42.98 % | 2.04 % | 21x |
| **3.0 deg (the default)** | **52.83 %** | **7.11 %** | **7.4x** |
| 5.0 deg | 62.62 % | 26.31 % | 2.4x |

At the default tolerance the current test accepts **more than half of all completely unrelated
boundaries** as same-parent. Map-scale impact (12 planted parents x 6 children, KS, tol 3.0,
5 trials): angle-only wrongly links 5-8 of the 11 cross-parent edges and recovers only **4-7 of 12
parents**; the axis-aware test wrongly links 0-2 and recovers **10-12 of 12**. Neither test misses
a single true same-parent edge, so the stricter test costs **nothing** in sensitivity.

The correct fingerprint already exists in the repo: `experimental/or_identification._fingerprint_set`
builds the deduplicated double coset $G_c (R S_p R^{\mathsf{T}}) G_c$ — exactly the set of
misorientations two children of one parent can show. True same-parent edges score 1.2e-6 deg
against it (the known quaternion/matrix round-trip floor); cross-parent edges sit at >= 2.36 deg.
Reconstruction simply does not use it.

**Performance caveat found:** the naive `einsum("eij,kij->ek")` used by the identification module
allocates E x 10664 floats — **4.3 GB at 50 000 edges**. Reformulated as a blocked (9,K) GEMM it is
**2.5x faster and 200x smaller** (1.6 s / 22 MB at block 512-1024), numerically identical to 4e-13.

**Revised Phase 1** is therefore: promote the fingerprint to a shared public core surface with the
blocked distance kernel, fix the reconstruction edge test, and pin the improvement in tests. This is
a **scientific behavior change** (reconstruction groupings change) and must be recorded under
`Fixed` in the CHANGELOG. The MTEX campaign scaffolding follows as Phase 2, and the synthetic
robustness study (documented fallback) as Phase 3.

## Phase 1 outcomes (2026-07-30) — DONE, edge-test defect fixed

- **New shared public core surface** in `core/transformation.py`:
  `intervariant_boundary_fingerprint(relationship)` (deduplicated
  $G_c (R G_p R^{\mathsf{T}}) G_c$) and
  `boundary_fingerprint_distances_deg(relative_matrices, fingerprint)`. Exported from
  `pytex.core` and top level. The private `_fingerprint_set` in
  `experimental/or_identification.py` is **deleted**; identification and reconstruction now
  share one definition (one-shared-helper rule).
- **Distance kernel is a blocked `(512, 9) @ (9, k)` GEMM.** The previous
  `einsum("eij,kij->ek")` in the identification module allocated E x K floats — 4.3 GB at
  50 000 edges. Now bounded near 22 MB, 2.5x faster, identical to 4e-13.
- **Reconstruction edge test now matches the full rotation.** `_reference_angles_deg` deleted.
- **Tests added (all green):**
  - `test_transformation.py`: 5 fingerprint tests — identity membership + symmetry closure
    (asserted as a *set-distance* statement, not by comparing rounded quaternion keys, which
    are brittle at rounding boundaries), exactness on all 276 KS variant-pair boundaries,
    the angle-only vs axis-aware discrimination rates, blocked-vs-direct equivalence at a
    non-multiple edge count, shape validation.
  - `test_parent_grain_reconstruction.py`: 6-parent stress fixture (**seed 187**, chosen
    because its cross-parent boundaries sit >10 deg from the fingerprint — the test *asserts*
    this separability rather than assuming it, per the Phase 11 fixture-design rule), plus a
    multi-seed sweep that skips seeds whose ground truth is genuinely ambiguous.
  - **Revert-verified:** temporarily restoring the angle-only rule makes both new
    reconstruction tests fail (6 parents collapse to 5). The regression bites.
- **Worked example** `or-ks-same-parent-boundary-fingerprint` (transformation group, gallery
  regenerated): the Sigma3 twin (60 deg / <111>) is an admissible KS same-parent boundary —
  Morito's published V1-V20 pair — and all 276 variant-pair boundaries score zero. Both are
  identities with independent provenance, not copied program output.
- **Docs synced:** CHANGELOG (`Fixed`, explicitly flagged as a scientific behavior change with
  the measured rates; plus `Added` for the new surface), symbol registry ($G_p$, $G_c$, and the
  double coset), phase-transformation validation matrix (new fingerprint row + rewritten
  reconstruction row), OR foundation doc F8 entry and the limits paragraph, and a new concept
  page section "Which boundaries can share a parent".
- **Environment fix (unblocked the gate):** spglib 2.7 emits an `OLD_ERROR_HANDLING`
  DeprecationWarning from inside its own filter context, which the warnings-as-errors policy
  turned into 9 failures in `test_transformation.py`. Added a fourth narrow commented
  exemption in `pyproject.toml`. Also installed the missing `hypothesis` dev extra.

### Verification (Phase 1)

- `python -m pytest --deselect tests/unit/test_phase_fixtures.py --deselect
  tests/unit/test_repo_integrity.py` -> **1287 passed, 26 deselected, zero warnings**.
- Sphinx `docs/site` -> **build succeeded**, no warnings.
- ruff/mypy clean on every file this phase touched.

### Known pre-existing failures on this machine (NOT caused by this work; verified at HEAD)

1. 6 phase-fixture SHA-256 mismatches (`test_phase_fixtures.py` hash pinning +
   `test_repo_integrity.py`). Documented in `working_notes_current_sprint.md` as pre-existing;
   the standing instruction is to deselect them when gating. Likely CRLF/LF checkout drift.
2. 20 ruff findings from newer rule versions (RUF022 `__all__` sorting, RUF059 unused unpacked
   variables, RUF043 unescaped `match=` patterns) in files this work did not touch.
3. 2 mypy `to_hex` arg-type errors in `plotting/crystal3d.py` (matplotlib stub drift).

None are in scope here; fixing them inside a scientific-behavior-change commit would obscure
the diff. They are logged as follow-ups.

## Current Status

- **Branch:** `main` (tracking `origin/main`)
- **Baseline for this goal:** commit `50116880` (Burgers beta→alpha zirconium tutorial)
- **Phase 1 committed** (see git log); Phase 2 (MTEX campaign scaffolding) is next.
- **Figure and Burgers program (FX1-FX9):** COMPLETE (pushed as latest commit)
  - All text layout issues fixed
  - Both Burgers directions exhaustively covered
  - Zirconium case with documented transformation strain

## Next Actions

### Phase 1: MTEX Fixture Generation Script
**Objective:** Create `scripts/mtex_generators/` directory with MATLAB scripts to generate reference fixture data.

**Tasks:**
1. Create martensite case: austenite parent with measured KS/NW/GT variants
2. Use MTEX `calcParent2Child` to generate pinned parent/child orientation pairs
3. Export as JSON fixture (parent_orientations, child_orientations, relationship_name)
4. Document variant generation and parent recovery in the script

**Files to Create:**
- `scripts/mtex_generators/generate_martensite_fixture.m` (MATLAB)
- `scripts/mtex_generators/generate_bcc_hcp_fixture.m` (MATLAB)
- `scripts/mtex_generators/README.md` (documentation)

**Verification:** Script runs in MATLAB R2021+, produces valid JSON with documented parent count

---

### Phase 2: Measured Data Integration
**Objective:** Load MTEX fixtures and integrate into test suite.

**Tasks:**
1. Add fixture loading to `tests/unit/test_parent_reconstruction.py`
2. Establish baseline parent-recovery accuracy on MTEX data
3. Document expected residuals and source (MTEX version, parameters)
4. Pin recovered parents in tests (assert within tolerance)

**Files to Modify:**
- `tests/unit/test_parent_reconstruction.py`
- `fixtures/` (add MTEX result directories)
- `docs/testing/mtex_parity_matrix.md` (add row)

**Verification:** MTEX fixtures load, parents recover, tests pass

---

### Phase 3: Validation Ledger Update
**Objective:** Complete the validation-matrix row for parent-grain reconstruction.

**Tasks:**
1. Document MTEX comparison (version, parameters, setup)
2. Record baseline accuracy (synthetic exact, MTEX data ±X deg)
3. Link reconstruction to OR foundation doc (F4 completion)
4. Update CHANGELOG entry

**Files to Modify:**
- `docs/testing/mtex_parity_matrix.md`
- `docs/architecture/orientation_relationship_analysis_foundation.md` (update F4 status)
- `CHANGELOG.md`

**Verification:** Validation matrix complete, links working, foundation doc updated

---

## Alternate Approach (if MTEX access unavailable)

If MATLAB/MTEX is not accessible in the current environment, proceed with synthetic-data expansion:

1. **Phase A1:** Create large synthetic fixtures with known phase transformations
   - 100+ planted parents with known variants
   - Various noise levels (0.5 deg, 1.0 deg, 2.0 deg)
   - Include cross-parent boundary ambiguities
   
2. **Phase A2:** Establish robustness baselines
   - Document reconstruction accuracy vs noise
   - Identify failure modes (grain size too small, noise too high)
   - Compare against literature noise-robustness studies

3. **Phase A3:** Update validation ledger conservatively
   - Mark as "foundationally ready" with synthetic-only validation
   - Flag MTEX parity as queued follow-on

---

## Quality Gates (every phase, before commit)

- `python -m pytest` (pass with no warnings)
- `python scripts/check_repo_integrity.py`
- `python -m ruff check .`
- `python -m mypy src`
- Sphinx build: `python -m sphinx -b html docs/site docs/_build/html`

---

## Phase Tracking

| Phase | Scope | Status | Outcome |
|-------|-------|--------|---------|
| 1 | MTEX fixture generation script | TODO | Fixture JSON files produced |
| 2 | Measured data integration | TODO | Tests pass with MTEX data |
| 3 | Validation ledger update | TODO | Documentation complete |

---

## Parallel Opportunities

After reconstruction stabilization completes, the next major features in roadmap order are:

1. **F13 (Habit Planes and PTMC)** — Deformation-gradient and invariant-line analysis for transformation-induced surface relief mapping. Large scope; requires its own working-notes document.

2. **Finding 8 remainder (Texture kernels)** — Bump and fibre kernels for ODF; SO3FunHarmonic integration for harmonic reconstruction.

3. **Finding 9 (PF→ODF Improvements)** — Ghost correction and zero-range handling in harmonic inversion.

4. **Finding 14 (Kikuchi Geometry)** — Gnomonic projection and Kikuchi-band overlays for EBSD pattern indexing.

5. **Human-auditable automated test documentation** (roadmap priority #1) — Create Sphinx pages documenting important tested methods with formulas, sources, examples, expected outputs, and verified code outputs. Meta-feature applicable to all subsystems.

---

## Key Decisions Logged

- **OR Flagship Status:** F1-F12 complete, F7/F12/F13/F14 remaining. Reconstruction stabilization is the gating item preventing F4 move from experimental to stable.

- **Reconstruction Scope:** MTEX parity is the validation floor, not the ceiling. Synthetic-data robustness studies will supplement measured-data validation.

- **Sequencing:** OR analysis remains the flagship per mission.md and AGENTS.md. Continue OR work until F13 and reconstruction stabilization land; then pivot to broader texture/diffraction/EBSD breadth.

---

## References

- `docs/roadmap/critical_review_and_development_guide.md` — Governance; Cycle C+ includes findings 5, 8-15, 21, 22
- `docs/architecture/orientation_relationship_analysis_foundation.md` — OR features F1-F14; reconstruction is F4
- `references/feature_opportunities.md` — Broader roadmap; human-auditable test documentation is priority #1 (meta)
- `mission.md` — OR analysis is designated flagship capability
- `AGENTS.md` — Working process includes resumability through progress notes
