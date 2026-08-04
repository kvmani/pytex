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

## Phase 2 outcomes (2026-07-30) — DONE, envelope measured + parity campaign defined

- **`scripts/study_reconstruction_robustness.py`** (+ `tests/unit/test_reconstruction_robustness_study.py`
  smoke lane, `--quick`). Sweeps noise x tolerance x grain count against planted ground truth;
  48 cells x 25 seeds. Trials whose ground truth is *genuinely* ambiguous at the tolerance under
  test are counted separately and excluded, so a physical coincidence is never scored as an
  algorithmic failure. Output JSON is git-ignored (matches the `performance_results.json` rule).
- **Headline: false-link rate is exactly 0.0 in all 48 cells** (~700 judged trials). The fixed
  edge test never merged two separable parents anywhere in the sweep.
- **Failure mode inverted:** the only remaining break is a *missed* link (splitting a parent).
  Governed by tolerance vs noise — at 2x noise the partition collapses (0-48% exact, 13-17%
  missed), at 4x it is essentially always exact. **Rule recorded: tolerance >= 4x scatter.**
  Counter-intuitive corollary: *more* children per parent makes exact recovery *harder* at
  marginal tolerance (more intra-parent edges = more chances to miss one).
- **Parent error tracks sigma/sqrt(n)** (0.5 deg noise: 0.316 / 0.200 / 0.135 deg at n = 2 / 5 /
  12, vs 0.354 / 0.224 / 0.144 predicted) — the eigen-mean refinement averages noise correctly.
- **Cost of loose tolerance quantified:** at 5 deg, ~20 of 25 random microstructures contain a
  genuinely ambiguous cross-parent boundary (vs 0-1 of 25 at 1 deg). That bounds the window.
- **Validation note** `docs/testing/reconstruction_robustness_study.md` (+ site stub + toctree +
  docs/README index + validation-matrix link). States its own limitations explicitly: synthetic,
  cubic-cubic KS only, chain-plus-contact adjacency rather than realistic map topology.
- **MTEX parity campaign `or_transformation_v1`** — campaign JSON, MATLAB handler
  (`mtex_parity_transformation.m`), dispatch entries, and PyTex-side generation in
  `generate_pytex_parity_campaign.py` (new `transformation` operation family). PyTex results
  generated and correct: KS 42.848 deg / <0.968 0.178 0.178>, NW 45.988 deg, GT recovered exactly
  from a KS nominal with the 2.4037 deg separation and zero residuals.
- **HONESTY GATE — read before trusting the MATLAB side:** `mtex_parity_transformation.m` has
  **never been executed** (no MATLAB on this machine). It encodes the *intended* comparison only.
  A mismatch on first run is most likely a script bug, NOT a PyTex/MTEX disagreement. Likely
  correction points are listed in `scripts/mtex_generators/README.md`: the
  `orientation.GreningerTrojano` spelling, whether `variants(p2c)` orders as PyTex's
  `generate_variants()`, and the `calcParent2Child` argument form. **No document claims MTEX
  parity for this stack**, and the ledger + foundation doc say so.

### Verification (Phase 2)

- `python -m pytest --deselect .../test_phase_fixtures.py --deselect .../test_repo_integrity.py`
  -> **1288 passed, 26 deselected, zero warnings**.
- Sphinx -> **build succeeded**, zero warnings (needed a site include-stub + toctree entry for the
  new study doc; a `docs/testing/` file linked from a site-rendered page must have one).
- ruff clean on all new files; integrity reports no campaign-schema errors.

## Phase 3 outcomes (2026-07-30) — map-scale sweep found the real limit

Closing the topology and non-cubic gaps that Phase 2 flagged as limitations changed the
conclusion, so Phase 2's headline should not be read on its own.

- **New `_plant_map` / `_map_cell` mode** in the study script: parents tile a square grid, each
  holding a 3x3 patch of child grains, four-connected adjacency over all of them. Every shared
  parent boundary therefore contributes *several* edges instead of one. Burgers (bcc->hcp, hcp
  child phase) added alongside KS.
- **The failure is asymmetric and the sparse sweep hid it.** One chance link anywhere along a
  shared boundary merges two parents irreversibly, so many edges per boundary = many chances.
  At 100 parents / 900 grains / ~1740 edges:

  | relationship | tol 1.0 | tol 2.0 | tol 3.0 (default) |
  | --- | --- | --- | --- |
  | Kurdjumov-Sachs | 97.5 / 100 | 88.5 / 100 | **70.0 / 100** |
  | Burgers | 100 / 100 | 98.5 / 100 | 97.0 / 100 |

- **The false-link rate stayed exactly zero throughout.** Every merge came from a *genuinely*
  ambiguous boundary — unrelated parents that really do share a fingerprint-consistent
  misorientation. Not an edge-test defect; an intrinsic limit of any binary edge test on
  orientations.
- **Therefore the tolerance is two-sided.** Phase 2's `>= 4x scatter` is a *lower* bound only.
  The map-scale upper bound is far more restrictive and, *as of this phase*, the binding one for
  real graphs: for cubic-cubic KS on a dense graph the window is roughly `4*sigma <= tol <= ~1 deg`,
  i.e. the method needs scatter below ~0.25 deg there. Phase 2's stale "default 3.0 is appropriate"
  sentence was corrected in the study doc.
  **-> Superseded by Phase 4:** the consistency vote removed most of this upper bound (69.7 -> 99.7
  of 100 parents at the default tolerance). Read Phase 4 before acting on the window above.
- **Mechanism identified:** the window is relationship-dependent because *fingerprint size* drives
  ambiguity. Burgers has 12 variants vs KS's 24, so 684 vs 10 584 distinct elements and
  <!-- Corrected 2026-08-04: this line originally read "~2 800 vs ~10 700". Both were wrong.
       The Burgers estimate gave the hexagonal child 24 proper operators instead of 12; the
       KS figure came from a quaternion deduplication that double-counted 81 elements. -->
  correspondingly less of orientation space. More variants = intrinsically harder reconstruction.
- **Runtime safeguard shipped (not just docs):** `ParentGrainReconstructionResult.chance_link_probability`
  — the probability two unrelated grains link at this tolerance, estimated from 4096 seeded
  uniform-random misorientations. Depends only on relationship + tolerance, **never on the data**,
  so it works on real maps with no ground truth. `describe()` warns on the *expected count* over
  the tested edges (not the rate) because a small rate over a dense graph is still unreliable.
  Verified against the independent measurement: 7.30% estimated vs 7.11% measured at 3 deg.
- **Default left unchanged at 3.0** deliberately — it suits the sparse/noisy regimes of the first
  sweep, and silently changing a default is worse than reporting the diagnostic. Documented.
- 2 new tests (diagnostic monotonicity/determinism/warning thresholds; graph-wrapper pass-through).

## Current Status

- **Branch:** `main` (tracking `origin/main`)
- **Baseline for this goal:** commit `50116880` (Burgers beta→alpha zirconium tutorial)
- **Phase 1 `7dd77d7b`**, **Phase 2 `2c7c9191`**, Phase 3 committing now.
- **NOT pushed** — every prior phase in this ledger was pushed to `origin/main`, but pushing is
  outward-facing and was left for the user to confirm.

## Phase 4 outcomes (2026-07-30) — the map-scale ceiling is broken

Phase 3 identified confidence-weighted clustering as the main remaining blocker and filed it as
future work. It is now implemented, because it was the one blocker actually within reach here.

- **`_vote_partition_cluster` + `_variant_descriptions`** in the reconstruction module.
  Connectivity proposes, single-parent consistency disposes: every cluster member proposes the
  parent it implies, each proposal is scored by how many members it explains
  (`C_j^T P` near the set `{S_c V_k}`), the best-supported proposal claims its supporters, and
  unexplained members repeat the vote.
- **Key simplification that makes it cheap:** the candidates `C_i V_k` for different `k` differ by
  a parent symmetry operation, and the support set is invariant under that — so they are the same
  hypothesis and **one proposal per member suffices**. Cost is `r^2 x |W|` per round, reusing the
  blocked `boundary_fingerprint_distances_deg` kernel. Verified empirically.
- **Result (identical map-scale sweep, 100 parents / 900 grains / ~1740 edges):**

  | tolerance | KS before | KS after | pure grains after |
  | --- | --- | --- | --- |
  | 1.0 deg | 98.0 | **100.0 (exact partition)** | 100% |
  | 2.0 deg | 88.7 | 99.7 | 97.8% |
  | 3.0 deg (default) | **69.7** | **99.7** | 95.2% |

  Burgers likewise 100.0/exact at <= 1.0 deg, 99.7 at 3.0. Holds with 0.25 deg added scatter.
- **No regression in the sparse sweep** — its cells are already single-parent clusters, which the
  vote returns whole. The three cells at 80-83% exact are the same pre-existing 12-children
  marginal-tolerance splitting cases, unchanged. Max false-link rate still 0.0000% everywhere.
- **Revert-verified regression test** built on a *searched worst case* (seed 1244): two parents
  49.18 deg apart whose children share a boundary only 0.36 deg from the fingerprint. The test
  asserts the edge test correctly *cannot* reject it (all 7 edges link) and that consistency
  separates the parents anyway. Disabling the vote makes it fail (1 parent instead of 2).
- Study script gained `pure_grain_fraction` and `partition_exact` for map cells — cluster *count*
  alone can hide compensating merge/split errors, so purity is reported alongside.
- Task chip `task_cf109b43` (F8 v2) dismissed as implemented.

## Phase 5 outcomes (2026-07-31) — a caveat I wrote was wrong; corrected

Probed the one assumption the new vote rests on that had never been measured: whether a small
parent has enough votes to defend itself against a large neighbour. Phase 4's doc claimed small
parents "remain the weak case". **That claim was wrong and is now corrected.**

- **First probe was invalid** and caught before it was believed: it measured absorption over all
  seeds, but in most of them the joining edge is *rejected* by the fingerprint test, so union-find
  separates the parents trivially and the vote never runs. Mean cluster count was a flat 2.00,
  which is what exposed it.
- **Conditioned on a genuinely ambiguous junction** (the only case where the question arises),
  absorption of the small parent is **5-11% and essentially flat in parent size** — 1 grain: 5.2%
  at 2 deg / 6.0% at 3 deg; 5 grains: 6.0% / 10.7%. If anything the *larger* small parent fares
  slightly worse, the opposite of the worry.
- **Mechanism:** vote counts only decide which proposal is considered first. Claiming a member is a
  per-grain consistency test (`C_j^T P` near the variant-description set), so a grain from a
  genuinely different parent is never claimed however many votes the neighbour has. Absorption
  needs the *additional* coincidence that the absorbing hypothesis also explains those grains,
  which is independent of either parent's size.
- Doc section rewritten with the measured table and the mechanism; test added pinning that a
  one-grain parent survives >80% of ambiguous junctions (measured ~94%, loose bound so the test
  pins the property, not a sampling artefact).

## Remaining work before reconstruction can leave `experimental`

1. **Run the MTEX side** of `or_transformation_v1` on a machine with MATLAB + MTEX 6.0, fix any
   script errors, record fixes in the generator README, then regenerate and
   `python scripts/compare_parity_results.py ...`. Only then may a parity claim be made.
3. **Measured-data fixture** (martensite→austenite or alpha→beta Ti) — the study is synthetic and
   says so; real EBSD noise is neither Gaussian nor independent per grain.
4. **Irregular grain geometry** — the map sweep uses square blocks, which captures edge *density*
   but not real grain shapes, size distribution, or boundary lengths.
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
