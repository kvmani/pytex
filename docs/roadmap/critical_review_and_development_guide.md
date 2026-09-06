# Critical Review And Development Guide (reconciled 2026-08-15)

This document is the governing development guide for PyTex going forward. It records a
repo-wide critical review of the implemented code, algorithms, tests, and documentation. The
original July 2026 findings remain below as the historical baseline; the verified state and
priority program were reconciled against the live tree on 2026-08-15. Every subsequent change —
human- or agent-authored — must follow the current program and the standing doctrine.

Precedence: this guide sits directly below `mission.md` and `specifications.md` and above all
roadmap documents. Where an older roadmap or working note conflicts with this guide, this guide
wins until the conflict is reconciled in writing.

Companion documents introduced with this guide:

- [Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md)
  — the feature and doctrine program for the OR-analysis flagship.
- `AGENTS.md` and `mission.md` (repository root) — refreshed alongside this review; the
  explainable-results doctrine below is now part of both.

## 1. State Of The Repository (verified 2026-08-15)

Facts checked directly against the working tree, not inherited from prior notes:

- The tree contains 138 Python source files, about 82.5k Python source lines, 558 symbols in
  `pytex.__all__`, 31 executable tutorial notebooks, and 2,059 test functions before
  parametrization. The complete base suite at the current head records two intentional skips and
  a runtime of about 12 minutes on the Windows development machine.
- Strict mypy and Ruff lint pass; `py.typed` ships; pytest treats new runtime warnings as errors;
  CI enforces an 87% whole-package coverage floor.
- The base CI lane covers Python 3.11–3.13 on Ubuntu and macOS. Windows coverage and a Sphinx
  warning-count ratchet are part of the 2026-08-15 governance increment; the full-scientific
  no-skip lane remains pinned to Ubuntu/Python 3.11.
- Worked-example framework, parity ledgers (MTEX, VESTA, diffraction, structure, plotting,
  phase-transformation), terminology registry, and citation policy all exist and are enforced in
  tests. This governance layer is the repository's strongest asset and must not erode.
- All original immediate-horizon items I1–I15 in the
  [World-Class Feature Roadmap](world_class_feature_roadmap.md) have landed, including the
  OR-deviation metric. Subsequent programs added OR fitting and experimental map-scale parent
  reconstruction, pole-figure arithmetic, GND density, TEM tilt navigation, kinematic and
  dynamical CBED foundations, stereographic Kikuchi maps, the Kearns parameter, and a shared
  desktop/intranet workbench.

Verdict: the foundations are genuinely strong — explicit semantics, real validation culture,
unusually disciplined documentation, and broad scientific capability. The present risks are
growth risks: governing prose and ledgers drifting behind implementation, stable code depending on
experimental internals, a 558-symbol flat stable namespace, browser behavior validated more by
manual driving than automation, documentation warnings hidden by a permissive build, and external
validation lagging the flagship OR implementation.

### 1.1 Disposition of the July findings

Section 2 is retained because its defect descriptions explain why the later architecture exists;
it is **not** a current missing-feature list. Current disposition:

| July findings | 2026-08-15 disposition |
| --- | --- |
| 1–4, 6–7: OR correspondence, deviation, fitting, reconstruction, consistency, vectorization | Implemented; map-scale reconstruction remains experimental pending measured-data and MTEX evidence. One new boundary defect remains: stable `core.parent_reconstruction` imports an experimental scoring primitive. |
| 5: habit plane / PTMC | Open, long-term; no current parity claim. |
| 8: kernel breadth | Implemented for the planned kernel family. |
| 9–10: harmonic completion, ghost correction, statistics | Partial; harmonic inversion, random-standard defocusing calibration, constrained named-component ODF fitting, and positivity / zero-range ghost correction (`correct_ghosts`, 2026-09-06) exist, but component-shape refinement and uncertainty breadth remain open. |
| 11: EBSD scale and grids | Partial; direct square/hex `.ang`, `.ctf`, `.oh5`/`.h5`, and topology graphs exist, but hex curvature/GND, the remaining HDF5 families, and out-of-core backing remain open. |
| 12: GND density | Implemented with theory, tests, and worked evidence. |
| 13: quantitative XRD | Partial; scattering, profiles, and preferred orientation exist, but measured-pattern I/O, background/multiphase fitting, and refinement remain open. |
| 14–15: TEM geometry and dynamical diffraction | Substantially advanced: Kikuchi, double diffraction, HOLZ, CBED, dynamical foundations, and the plane-parallel finite-thickness kinematic shape factor exist; ring SAED, specimen-thickness distributions/bending, and full dynamical/Rietveld breadth remain open. |
| 16: explainable results | Broadly adopted, but several stable texture/indexing reports and JSON contracts still need closure. |
| 17–20: coverage, warnings, CI, property tests | Coverage, runtime warning hygiene, OS/Python matrix (except Windows), and property suites implemented. Documentation warnings require the new ratchet. |
| 21–22: performance and release engineering | Performance has a runnable quick benchmark but no CI regression lane. Release engineering is closed: `CHANGELOG.md`, the stability and release policy in `docs/standards/api_stability_and_deprecation.md`, and the first cut version `0.1.1`. |
| 23–25: documentation index, stale foundations, OR theory | Original gaps closed; continuing roadmap/ledger drift and Sphinx reference warnings are the active documentation defects. |

## 2. Critical Findings

Ordered by subsystem. Each finding states the defect or gap, why it matters, and the required
direction. Items marked **[P1]** are the next-cycle obligations; **[P2]** medium term; **[P3]**
long term.

### 2.1 Transformation and orientation-relationship analysis

This is the designated flagship, and it is also where the gap between "primitives exist" and
"analyses are possible" is widest. Full program:
[Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md).

1. **[P1] `OrientationRelationship` is rotation-only; there is no lattice correspondence.**
   `map_parent_vector_to_child` rotates Cartesian unit vectors. There is no API that answers the
   canonical OR questions: *given the parent plane (hkl) or direction [uvw], what are the
   corresponding child indices in each variant — and vice versa?* That requires mapping through
   the structure matrices (direct basis for directions, reciprocal basis for plane normals), a
   rationalization step to nearest low-integer indices, and an explicit angular residual. This is
   the single highest-value missing feature in the library.
2. **[P1] No misorientation-space representation or deviation metric for ORs.** Roadmap item I10
   is unlanded: there is no `or_deviation(...)` that reports the mean/max angular departure of
   measured parent-child pairs from a nominal OR, and no way to express an OR *as* a
   symmetry-reduced misorientation (axis/angle, e.g. KS = 42.85° about <0.968 0.178 0.178>).
   Both are prerequisites for OR fitting and for honest literature comparison.
3. **[P2] No OR fitting.** Given measured parent/child orientation pairs (or child-child boundary
   misorientations when the parent is consumed), PyTex cannot estimate the operative OR. The
   experimental `score_parent_orientations` scores *given* candidates only. Symmetry-aware
   quaternion averaging / MODF-peak refinement is the required path (MTEX `calcParent2Child`
   parity floor).
4. **[P2] No map-scale parent-grain reconstruction.** The grain-boundary network, union-find
   merging (`merge_by_csl`), and variant machinery all exist separately; the voting/growth
   algorithm bridging them is the flagship deliverable of the next two cycles.
5. **[P3] Habit-plane slot is dead weight.** `TransformationVariant.habit_plane_pairs` carries
   data but nothing computes it. Deformation-gradient/Bain-strain and invariant-line/PTMC
   analysis are the long-term completion of the subsystem; until then the docstring must say the
   slot is descriptive, not computed.
6. **[P1] Internal inconsistencies to fix now:**
   - `_phase_semantically_matches` is duplicated in `core/transformation.py` and
     `experimental/phase_transformation.py` with *different* strictness (proper-point-group
     equality vs full `SymmetrySpec` equality). One shared helper, one definition.
   - `OrientationRelationship.parallel_directions` stores naked `np.ndarray` pairs while
     `parallel_planes` stores `CrystalPlane` pairs — the asymmetry violates the repository's own
     no-naked-arrays rule and loses phase/frame meaning on the direction side.
   - Variant indices are enumeration-order artifacts. Published variant tables (Morito V1–V24 for
     KS) are a de facto standard; the variant doctrine must state the numbering convention and
     pin it in tests, or explicitly document non-conformance.
7. **[P2] Vectorization debt on the hot path:** `intervariant_misorientations` runs a Python
   double loop over variant pairs; `PhaseTransformationRecord.predicted_child_orientations` and
   `select_variants` build per-element Python lists of `Rotation` objects. Correct today at 24
   variants; wrong shape for map-scale reconstruction where these become inner loops.

### 2.2 Texture

8. **[P2] Kernel breadth.** One kernel (de la Vallée Poussin). Gaussian-on-SO(3), Abel-Poisson,
   von Mises-Fisher, and fibre kernels are required for MTEX-comparable ODF work;
   `KernelSpec.as_so3_kernel()` must generalize.
9. **[P2] Harmonic layer is partial.** `HarmonicODF` lacks quadrature from scattered
   orientations, convolution, and Funk/pole-figure projection at full strength; PF→ODF inversion
   lacks zero-range and ghost correction. These gate any quantitative texture-strength claims.
10. **[P3] Statistics.** No Bingham components, no bootstrap/confidence machinery on ODFs or
    volume fractions. Required before uncertainty language enters public docs.

### 2.3 EBSD

11. **[P2] In-RAM gap.** Hex-grid `CrystalMap`, `.ang`/`.ctf`/`.oh5` import, KAM, and
    segmentation now exist, but hex curvature/GND stencils, the remaining h5ebsd-family readers
    (Oxford H5OINA, Bruker), and chunked backing do not. Chunked backing above all still bounds
    the size of real datasets PyTex can ingest: the EDAX HDF5 reader loads a scan whole.
12. **[P2] No GND density** (curvature-tensor route) despite KAM/GROD existing — the natural
    next step users expect after local-misorientation maps.

### 2.4 Diffraction

13. **[P2] Kinematic-intensity realism stops at the single pattern.** No background models,
    multi-phase summation, or preferred-orientation corrections (March-Dollase and — the unique
    PyTex opportunity — ODF-weighted intensities coupling the texture core into powder XRD).
14. **[P2] TEM: no Kikuchi/gnomonic geometry surface, no ring (Debye-Scherrer) SAED, no
    double-diffraction flagging.** The Kikuchi geometry layer is the highest-priority TEM item.
15. **[P3] Dynamical diffraction and Rietveld** remain the long-term differentiators; the
    medium-term profile/background/Pawley stack is their prerequisite and must land first.

### 2.5 Explainability and verbal outputs (cross-cutting, new doctrine)

16. **[P1] The mission promises interpretable results; the code barely delivers prose.** Exactly
    two `summary()` methods exist (both returning bare dicts in `ebsd/`). Report objects
    (`TextureReport`, `VariantSelectionReport`, `ParentReconstructionReport`, …) carry numbers
    but cannot explain themselves. Adopt the **explainable-results doctrine** (§4) — a uniform
    `describe()` surface producing cited, convention-explicit scientific prose from every stable
    report object, validated in tests like any other output.

### 2.6 Testing and engineering

17. **[P1] No coverage measurement.** `pytest-cov` is installed but unused in CI. Add a coverage
    report and a ratchet (fail-under that only rises). Target: ≥90% for `core/`, ≥85% overall.
18. **[P1] Warning hygiene.** 117 warnings in a green run, including matplotlib figure leaks
    ("More than 20 figures") from plotting tests. Add an autouse close-figures fixture, then turn
    on `filterwarnings = error` with explicit, commented exemptions.
19. **[P2] CI breadth.** Matrix ubuntu+macos × Python 3.11–3.13 for the base lane (full lane may
    stay pinned); publish the coverage artifact.
20. **[P2] Property-based tests.** Orientation algebra, symmetry reduction, index round-trips,
    and OR variant orbits are ideal Hypothesis targets (composition laws, involution of
    inverses, orbit invariance). Currently absent.
21. **[P2] Benchmarks are manifests without timings.** Stand up a runnable performance lane with
    pinned map-scale cases so vectorization claims are evidence-backed.
22. **[P2] Release engineering.** No CHANGELOG, no release process, version `0.1.0.dev0`
    untagged. Define the versioning and changelog policy now, before external users exist.

### 2.7 Documentation

23. **[P1] The docs index has drifted.** `docs/README.md` omits documents that exist and are
    normative: the VESTA parity matrix, the phase-transformation validation matrix, the API
    stability standard, the executable-examples standard, and every roadmap except one. Fixed
    alongside this review; keep the index complete as a review-gate obligation.
24. **[P1] Stale claims.** `docs/architecture/phase_transformation_foundation.md` still says
    variant doctrine is unimplemented and transformation code "should remain experimental" —
    contradicted by the shipped, tested stable surface. Foundation docs must be updated in the
    same change that lands the capability they describe.
25. **[P2] Theory notes lag the transformation subsystem.** `docs/site/theory/` has no OR/variant
    theory note while orientation, EBSD, and diffraction all have them. The OR foundation work
    must land with its LaTeX note, per standing policy.

## 3. Priority Program

The sequencing rule remains foundations before breadth, with OR analysis the flagship when
priorities compete within a horizon. The active program as of 2026-08-15 is:

1. **Governance and executable gates:** reconcile the roadmap and durable ledger, add Windows to
   the base CI lane, forbid growth in Sphinx warnings, and add a minimal real-browser Playwright
   lane for critical workbench behavior.
2. **Architecture and validation hardening:** remove the stable-core dependency on experimental
   parent scoring, finish `describe()`/JSON coverage for stable reports, and execute the prepared
   MTEX OR-fitting/reconstruction campaign before stronger parity claims.
3. **Five-feature delivery cycle, in dependency order:** measured powder-XRD I/O and comparison;
   random-standard defocus calibration; hex-grid EBSD; finite-thickness SAED; named-component ODF
   fitting.
4. **Feature completion rule:** every feature includes tests, an independently known numerical
   reference, `describe()`, portable JSON where appropriate, a worked example, theory/workflow
   documentation, a parity-ledger update, and a benchmark case. Synthetic measurements are
   permitted only when no redistributable real data is available and must be labelled as
   synthetic in the fixture metadata, worked example, and validation prose.
5. **After this cycle:** ghost/zero-range correction, measured OR/EBSD validation, and
   quantitative multiphase XRD fitting are the next scientific priorities; habit-plane/PTMC and
   full Rietveld remain longer-horizon programs.

## 4. Explainable-Results Doctrine (normative)

PyTex results must be explainable to a reader who did not run the code. From this review
forward:

1. Every stable report/result object gains a `describe()` method returning structured prose
   (a `str` or a small `Description` object with `text` and `citations`): what was computed,
   under which conventions (frames, symmetry, Euler convention), key numbers with units, and
   what the numbers mean scientifically — e.g. a `VariantSelectionReport.describe()` names the
   OR, the variant count, the dominant variants, and the residual statistics, and flags
   selection strength against a random-selection baseline.
2. `describe()` text is a tested surface: unit tests assert the presence of the governing
   convention statements and key computed numbers (not exact prose), so explanations cannot
   silently desynchronize from results.
3. Verbal output follows the terminology registry; symbols and terms used in prose are the
   registry's, and citations use the reference canon.
4. JSON contracts and `describe()` grow together: any object with a canonical JSON contract must
   also be describable, so machine and human interchange stay in lockstep.
5. Explanations state limitations honestly (kinematic-only intensities, unvalidated regimes,
   ambiguity flags) using the same language as the validation ledgers.

## 5. Standing Quality Bars (normative)

These consolidate and extend the existing standards; they are review-gate obligations for every
merge:

- **Scientific completeness:** a stable numerical feature ships with theory note (`docs/site/theory/`),
  concept/workflow page, executable worked example with independently-provenanced expected
  values, validation-ledger row, and `describe()` support. Partial landings must be staged in
  `experimental/`.
- **Algorithmic honesty:** every algorithm documents its assumptions, normalization, failure
  modes, and the regime in which it is validated; claims beyond the validated regime are
  forbidden in docs and docstrings.
- **Vectorization:** batch paths use array primitives end to end; a Python per-element loop on a
  path that scales with map or set size is a defect, not a style choice.
- **Semantics:** no naked arrays on stable surfaces where frame/phase/symmetry meaning exists —
  including inside existing objects (finding 6 shows internal drift is possible).
- **Hygiene:** new warnings are defects; coverage may not decrease; the docs index and
  foundation docs are updated in the same change as the capability.
- **Convention pinning:** wherever a literature convention exists (variant numbering, axis/angle
  representatives, section conventions), PyTex either conforms and pins it in tests or documents
  the deviation explicitly.

## 6. Review Cadence

- Re-run this repo-wide critical review every two development cycles or after any major
  subsystem lands; update this document in place with a dated changelog entry below.
- The subsystem scorecard in
  [Repo Review: 2026 Foundation Audit](../architecture/repo_review_2026_foundation_audit.md)
  remains the historical baseline; this guide supersedes its recommendations.

## Changelog

- **2026-08-15 reconciliation:** Re-verified repository scale and quality lanes; converted the
  July findings into an explicit disposition table; replaced the completed Cycle A/B program with
  the governance-plus-five-feature delivery cycle; recorded Windows, Sphinx-warning, browser-test,
  stable/experimental-boundary, report-contract, and external-validation debt without reopening
  capabilities that have already landed.
- **2026-07:** Initial review and guide; OR-analysis flagship program established;
  explainable-results doctrine adopted.
- **2026-07 (Cycle A executed):** Findings 1, 2, 6, 16, 17, 18, 23, 24 closed — index
  correspondence (F1–F3), OR misorientation + deviation (F5), parallelism finders (F4),
  `describe()` on all transformation reports, `SymmetrySpec` equality fix, unified phase
  matching, typed parallel directions, warnings-as-errors with a figure-close fixture (suite
  runs warning-free), and an 87% coverage ratchet in CI (measured 88%). Remaining Cycle A
  spillover: none. Next per §3: Cycle B (findings 3, 4, 7, 19, 20).
- **2026-07 (Cycle B executed):** Findings 3, 7, 19, 20 closed; finding 4 delivered at its
  planned experimental stage — variant hot paths vectorized, CI matrix (ubuntu+macos ×
  3.11–3.13), OR fitting (`fit_orientation_relationship`, recovers GT exactly from a KS
  nominal), first Hypothesis property suites, and experimental map-scale parent-grain
  reconstruction with quaternion-averaged parent refinement (synthetic planted-parent
  validation; literature fixtures and EBSD grain-graph wiring queued before stabilization).
  Next per §3: Cycle C+ (findings 5, 8–15, 21, 22), plus the queued ledger follow-ups.
- **2026-07 (finding 26, found and fixed):** The transformation stack composed predicted
  children as ``V @ P``, which contradicts the normative crystal-to-specimen orientation
  convention (correct: $g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}}$); internally self-consistent synthetic
  validations masked it, and it would have produced wrong variant assignments on real
  measured orientations (14–16 deg spurious residuals on canonically built data). Fixed
  across prediction, deviation, fitting, scoring, and reconstruction; a regression test now
  pins the specimen-space parallelism identity. Lesson recorded: synthetic round-trip tests
  must build inputs through the *canonical convention*, never through the code under test's
  own composition.

## References

### Normative

- `mission.md` and `specifications.md` (repository root)
- [Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md)
- [Development Principles](../standards/development_principles.md)
- [Engineering Governance](../standards/engineering_governance.md)

### Informative

- [World-Class Feature And Foundation Roadmap](world_class_feature_roadmap.md)
- [Implementation Roadmap](implementation_roadmap.md)
- [Repo Review: 2026 Foundation Audit](../architecture/repo_review_2026_foundation_audit.md)
