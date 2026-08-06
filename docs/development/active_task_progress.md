# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Current Task: TEM Tilt Navigation Program (TN) — IN PROGRESS (started 2026-08-06)

**Objective.** Implement, in full, the formulation approved in
[`docs/architecture/tem_tilt_navigation_foundation.md`](../architecture/tem_tilt_navigation_foundation.md):
given an indexed current zone axis and holder position, compute the holder $\alpha$/$\beta$ tilts
that bring a requested target zone axis onto the electron beam — with symmetry handling, ambiguity
classification, reachability, path planning, calibration, uncertainty, a publication-quality
annotated stereographic projection, tests, worked examples, docs and a runnable notebook.

**Governing documents.** `AGENTS.md`; the TN foundation document (normative for the science);
`docs/standards/notation_and_conventions.md` (frame domains — TN invents none);
`docs/standards/executable_examples.md`; `docs/standards/visualization_style_guide.md`.

**Branch policy for this task:** `main` only. No feature branches. Commit and push incrementally.

### Scope decisions pinned at the start (do not re-litigate on resume)

1. **No animation.** The user clarified: a *static* publication-quality stereographic projection
   showing the zone-axis trajectory as a series of dots from current to target, annotated with
   principal zone axes, the $\alpha$/$\beta$ tilt axes, and the reachable region. §12.5 of the
   foundation document is superseded on this point; the frame-sequence primitive is retained
   because the dotted trajectory *is* the sampled path.
2. **Holder frame is declared to be the pytex specimen-domain frame** (foundation §2.2). No new
   `FrameDomain` member is added. $\mathbf{U}$ (crystal→holder) is a plain `Orientation`.
3. **Mode B (two indexed zones) is the recommended reconstruction path** because it needs no
   diffraction-rotation calibration (foundation §5.2).

### Phase ledger

| Phase | Deliverable | Status |
| --- | --- | --- |
| TN0 | Foundation document, reviewed and approved | **done** (2026-08-06) |
| TN1 | `tem/stage.py`: kinematics, calibration record, envelopes | **done** (`bd0be70`) |
| TN2 | `tem/navigation.py`: closed form, branches, orbit, validation, ranking | **done** (`bd0be70`) |
| TN3 | `tem/reconstruction.py` + `tem/ambiguity.py` | **done** (`bd0be70`) |
| TN4 | `tem/path.py`: geodesic/Kikuchi-band paths, waypoints, backlash | **done** (`bd0be70`) |
| TN5 | `tem/calibration.py` + uncertainty propagation | **done** (`bd0be70`, YAML contract in `34755e1`) |
| TN6 | `plotting/tilt_stereogram.py`: the annotated figure | **done** (`4a16430`) |
| TN7 | `tests/unit/test_tem_tilt_navigation.py` (foundation §14 matrix) | **done** (`4a16430`, `34755e1`) |
| TN8 | Worked examples + schemas + top-level exports | **done** (`34755e1`) |
| TN9 | Notebook `24_tem_tilt_navigation.ipynb` | **done** (`bcdcfa1`) |
| TN10 | Theory note, site pages, indexes, registry, parity matrix | **done** |

### What TN10 delivered

| Artifact | Path |
| --- | --- |
| Canonical LaTeX derivation | `docs/tex/algorithms/tem_specimen_tilt_navigation.tex` |
| Rendered algorithm page | `docs/site/algorithms/tem_tilt_navigation.md` |
| Design record, marked delivered | `docs/architecture/tem_tilt_navigation_foundation.md` |
| Registered in | algorithms toctree, theory index, notebooks index, architecture index, `docs/README.md` |
| Terminology registry | 7 new terms (holder frame, crystal-to-holder orientation, holder tilt, tilt envelope, diffraction rotation, observation stabilizer, ambiguity family) |
| MTEX parity matrix | new row; **no parity claimed in either direction** — MTEX has no equivalent surface, this being instrument-operation geometry rather than texture analysis |

### Verification results (measured, not assumed)

| Gate | Result |
| --- | --- |
| `pytest tests/unit` | **3909 passed** before TN7; **3993 passed** after (84 new TN tests) |
| `ruff check` | clean across `src/pytex`, `tests`, `worked_examples` |
| `mypy` (strict) | clean, 105 source files |
| Worked-example gallery | 6 new TN examples, all computed == expected |
| Public exports | 428 → 463 (35 TEM surfaces) |

Numerical claims verified against closed form: moving-beta-axis cancellation
(exact); beam-direction formula (exact); closed-form tilt solution (max error
4.6e-16 over 20 000 random directions); all four solution branches; residual law
`2 asin(sin(dphi/2) sin theta)` (max deviation 1.3e-10 deg); the 180-degree
error negating both angles exactly; two-zone reconstruction (1.2e-6 deg over 200
trials); Laue-vs-proper enlargement over all 32 point groups (factor exactly 2
for the ten affected).

### Open issue to raise with the user

**Commit `34755e1` has a wider footprint than its message describes.** It was
staged with `git add -A src/pytex`, which swept in pre-existing uncommitted
docstring work from the completed Release-Readiness program that was sitting
unstaged in the working tree at session start. Nothing was lost or altered — the
full suite passes — but the commit message describes only the TN worked-example
and contract work. The history is already pushed to `main`, so correcting it
would need a force push, which was not done unilaterally. Roughly 26 further
pre-existing modified/untracked files (docs, tests, `docs/tex/algorithms/*.tex`)
remain uncommitted and untouched.

### Next actions

TN9: write and execute `docs/site/tutorials/notebooks/24_tem_tilt_navigation.ipynb`,
then TN10 documentation.

---

## Previous Task: Release-Readiness Program (RR) — COMPLETE (2026-08-05)

**Objective.** Identify and close the remaining gaps in foundational methods, classes, and
documentation across the whole repository, extend the library where a genuinely high-impact
capability is missing, and bring PyTex to a state suitable for a first public release. Every
change must add verified scientific or usability value; no filler.

**Governing documents.** `AGENTS.md`, `mission.md`, `specifications.md`,
[`docs/roadmap/critical_review_and_development_guide.md`](../roadmap/critical_review_and_development_guide.md)
(Cycle C+ obligations: findings 5, 8–15, 21, 22), and
[`docs/standards/documentation_architecture.md`](../standards/documentation_architecture.md).

### Verified baseline (measured 2026-08-05, not inherited from notes)

| Gate | At session start | At program end |
| --- | --- | --- |
| `pytest` (base lane) | 1473 passed | **3882 passed** |
| `ruff check .` | **16 errors** (11 × RUF059, 5 × RUF043) | clean |
| `mypy` (strict) | **2 errors** (`plotting/crystal3d.py` `to_hex`) | clean |
| Public docstring coverage | **33.1%** (443/1338) | **90.5%** (1250/1381) |
| Undocumented `pytex.__all__` surfaces | **574** | 0 |
| Public classes with only an auto dataclass signature | **118** | 0 |
| `pip install` then `import pytex` | **fails** (undeclared `pyyaml`) | verified in a clean venv |
| Version literals in the tree | **4**, free to drift | 1, ratcheted |
| Public exports | 411 | 428 |
| `hypothesis` in `.venv` | **missing**, broke collection | installed |

### The audit that set the program

Measured directly against the working tree with an AST/introspection sweep:

- **Docstring gap (the largest single documentation defect).** Of the 411 names in
  `pytex.__all__`, **112 carry no docstring at all**; a further **430 public class members**
  across 112 classes are undocumented — worst concentrations `CrystalMap` (25),
  `DiffractionGeometry` (19), `GrainSegmentation` (16), `OrientationSet` (16), `Rotation` (15),
  `SymmetrySpec` (15), `PointGroup` (14). Repo-wide public docstring coverage is **33.1%**
  (443/1338). This directly contradicts the standing rule in `AGENTS.md` that every documented
  method state purpose, when-to-use, inputs, and outputs, and it is the first thing an external
  user meets. **This is RR1.**
- **Scientific capability gaps confirmed absent from `src/` by grep**, all of them Cycle C+
  items the guide already names:
  - no gnomonic projection and no Kikuchi geometry surface anywhere (finding 14, the roadmap's
    highest-priority TEM item, and the shared geometric foundation of EBSD pattern formation);
  - no March–Dollase and no ODF-weighted preferred-orientation intensity correction
    (finding 13, the roadmap's identified *unique* PyTex strength — the texture core driving
    powder intensities);
  - no GND density (finding 12), although KAM/GROD/GOS/GAM all exist in `ebsd/models.py`;
  - no hex-grid `CrystalMap` and no h5ebsd-family reader (finding 11).
- Findings already closed by earlier programs and re-verified here: texture kernel breadth
  (dVP + Gaussian-SO(3) + Abel–Poisson all present in `texture/kernels.py`), grain scalar
  metrics, cleanup filters, CHANGELOG.

### Phase ledger

| Phase | Scope | State |
| --- | --- | --- |
| RR0 | Gate repair + baseline | **Complete** — ruff and mypy restored to clean; baseline recorded above |
| RR1 | Public-API docstring completion + regression test | **Complete** |
| RR2a | Kikuchi / gnomonic geometry foundation | **Complete** |
| RR2b | Powder-XRD preferred orientation (March–Dollase + ODF-weighted) | **Complete** |
| RR2c | GND density from orientation gradients | **Complete** |
| RR3 | Release engineering sweep | **Complete** |

Hex-grid `CrystalMap` and h5ebsd readers (finding 11) are **deliberately out of scope** for this
program: they are a data-ingestion project of their own and do not block a first release of the
computational surfaces. Recorded here so the decision is not re-litigated.

### RR0 — Complete

Eleven `RUF059` unused-unpacked-variable findings (underscore-prefixed) and five `RUF043`
unescaped-regex-metacharacter findings in `pytest.raises(match=...)` patterns (made raw and
escaped — these were silently weakening the assertions, since an unescaped `.` matched any
character). Two strict-mypy errors in `_separated_species_colors` fixed with an explicit
`_rgb_tuple` helper, because `matplotlib.colors.to_hex` does not accept an `ndarray`.

### RR1 — Complete

**Outcome.** Every name in `pytex.__all__` and every public member of every exported class now
carries a docstring: **574 members plus 118 class docstrings**, across all eight subpackages.
Repo-wide public docstring coverage rose from **33.1% to 90.2%** (1207/1338); the remaining 10%
is internal surfaces outside `__all__`.

| Gate | Before RR1 | After RR1 |
| --- | --- | --- |
| `pytest` | 1473 passed | 3614 passed (the ratchet adds ~2140 parametrized cases) |
| `ruff` / `mypy` | clean | clean |
| `pytex.__all__` surfaces undocumented | 574 | 0 |
| Public classes with only an auto-generated dataclass signature | 118 | 0 |

**Method.** Docstrings were inserted by AST location with a scratchpad helper that refuses to
touch an already-documented node, re-parses the file, re-checks that every target now carries a
docstring, and re-tokenizes — so a docstring cannot land inside the wrong function or inside a
string literal. Every batch was written after reading the implementation, then verified with
`ruff` plus the affected unit tests.

Depth was matched to surface importance rather than applied uniformly: full NumPy-style
Purpose/Parameters/Returns for scientific entry points, one accurate sentence for simple
accessors. A Parameters block for `Grain.size` would be filler, not documentation. Where an
algorithm has real limits — kinematic-only intensities, ill-posed PF-to-ODF inversion, centring-only
absences, `habit_plane_pairs` being descriptive rather than computed — the docstring states them,
per the algorithmic-honesty bar in the development guide.

**The ratchet.** [`tests/unit/test_public_api_docstrings.py`](../../tests/unit/test_public_api_docstrings.py)
parametrizes over every exported surface and asserts two structural properties: a docstring exists,
and its first line is a real summary rather than a section header or a two-word stub. It is
deliberately structural, so it cannot be satisfied by a placeholder and cannot go stale when prose
is improved.

**Two defects the ratchet found that the initial audit had missed:**

1. **118 public classes carried only Python's auto-generated dataclass signature**
   (`DeLaValleePoussinKernel(halfwidth_deg: 'float')`) as `__doc__`. An audit that only tests
   `__doc__` for emptiness scores these as documented; they are the *first* thing a user reads.
   All 118 now have real class docstrings.
2. Five `pytest.raises(match=...)` patterns contained unescaped regex metacharacters (RR0). An
   unescaped `.` matches any character, so those assertions were silently weaker than they read.

### RR2a — Complete

Closes critical-review finding 14, the roadmap's highest-priority TEM item, and supplies the
geometric layer that EBSD and TEM Kikuchi patterns share.

**New module** [`src/pytex/diffraction/kikuchi.py`](../../src/pytex/diffraction/kikuchi.py):
`GnomonicProjection`, `KikuchiBand`, `KikuchiZoneAxis`, `KikuchiPattern`,
`simulate_kikuchi_pattern`, plus `plot_kikuchi_pattern` in `plotting/diffraction.py`. Built on the
existing `DiffractionGeometry` rather than a private frame model, per the no-private-frames rule.

**The scientific point.** A Kikuchi band is the pair of curves where a lattice plane's two Kossel
cones meet the detector. In the gnomonic projection every great circle maps to a straight line, so
band *centre* lines are exactly straight whatever the detector tilt — which is why band detection
belongs in gnomonic coordinates. Band *edges* are cones, not great circles, so they are conics;
**PyTex does not make the usual small-angle approximation that draws them as straight parallel
lines**, and a test asserts an edge is measurably non-linear so a linearized implementation would
fail. Band angular width is exactly `2*theta_B`, making width a direct measurement of d-spacing.

**Verification — closed-form anchors, not recorded outputs:**

| Anchor | Result |
| --- | --- |
| Gnomonic radius = `tan(psi)` at 0/15/30/45/60 deg | exact to 1e-12; pins the projection scale with no crystallography |
| Band centre satisfies its own line equation, at 3 detector tilts incl. 25 deg | < 1e-14 relative |
| Both edges lie on their Kossel cones | < 1e-12 |
| Ni {111} band width at 20 kV, Bragg's law by hand | 2.4187 deg computed vs 2.4188 deg simulated |
| `[011]` at cube orientation → gnomonic radius | exactly 1.0 — pins the whole crystal→specimen→lab→detector chain |
| Analytic apparent width vs direct numerical measurement | < 1e-3 relative |

A **real defect was found and fixed** during this phase: `GnomonicProjection.contains` used strict
inequalities, so a point on the exact detector edge — which round-trips with a few units of
floating-point error — was reported as off-detector. `detector_corner_coordinates()` therefore
reported its own corners as outside the detector. Now takes an explicit `atol_px`.

**Deliverables:** 22 unit tests in
[`tests/unit/test_kikuchi_geometry.py`](../../tests/unit/test_kikuchi_geometry.py); two executable
worked examples (`diffraction-ni-111-kikuchi-band-width`,
`diffraction-gnomonic-zone-axis-radius`) with hand-derived expected values and citations; the LaTeX
theory note
[`docs/tex/algorithms/kikuchi_bands_and_gnomonic_projection.tex`](../tex/algorithms/kikuchi_bands_and_gnomonic_projection.tex);
the workflow page [`docs/site/workflows/kikuchi_geometry.md`](../site/workflows/kikuchi_geometry.md);
five new rows in the diffraction validation matrix; a registered plotting-validation case; and
`describe()` on all three report objects, tested.

Suite: **3686 passed**, ruff and mypy clean.

### RR2b — Complete

Closes critical-review finding 13 and delivers what the roadmap identified as the *unique* PyTex
strength: the texture core driving powder diffraction intensities.

**New module**
[`src/pytex/diffraction/preferred_orientation.py`](../../src/pytex/diffraction/preferred_orientation.py):
`march_dollase_factors`, `MarchDollaseModel`, `ODFPreferredOrientationModel`,
`PreferredOrientationModel` (a runtime-checkable protocol), `preferred_orientation_factor_table`;
plus `apply_preferred_orientation` and a `preferred_orientation=` argument on
`generate_xrd_pattern`. The dependency runs one way — `xrd` imports the models, not the reverse.

**Verification — closed-form anchors:**

| Anchor | Result |
| --- | --- |
| March distribution integrates to 1 over the sphere, for every `r` | equal-area quadrature at 5 coefficients |
| `P(0) = r^-3`, `P(pi/2) = r^(3/2)` | exact to 1e-12 |
| Cubic {111} family factor under a (111) axis at `r = 1/2`, hand-derived as `2 + 162/(65*sqrt(65))` | agrees to **1e-16** |
| Uniform ODF ⇒ factor 1 (the random-texture identity) | 0.9999 |
| Planted ND fibre ⇒ the planted pole is the enhanced one | confirmed |

**A real scientific defect found and fixed.** `ODF.evaluate_pole_density` returns a *kernel-weighted
response*, not a value in multiples of random — the smoothing kernel peaks at 1 rather than
integrating to 1, so a uniform ODF returns the kernel's spherical mean (~0.024 for a 15° dVP
kernel), not 1. The first draft of `ODFPreferredOrientationModel` claimed "multiples of random" and
would have scaled every intensity by ~1/41. Fixed by dividing by the kernel's spherical mean,
computed in closed form by Gauss–Legendre quadrature in `u = cos(omega)`. The uniform-ODF test is
what caught it, and it is now the pinned anchor. **Note for anyone reading pole densities
elsewhere: `evaluate_pole_density` is relative, not m.r.d.**

**Two documentation defects also found and fixed:** `S2Grid.equispaced`/`regular` docstrings (written
in RR1) listed the hemisphere values as `"upper"`/`"full"` when the accepted values are
`"upper"`/`"sphere"`; and the `pytex.__init__` export block was missing three imports that were
already listed in `__all__` — caught immediately by the RR1 docstring ratchet, which is the second
time that test has paid for itself.

`march_dollase_factors` returns exactly `1.0` at `r = 1` by an explicit branch rather than through
`cos^2 + sin^2`, so "no preferred orientation" leaves intensities bit-identical.

**Deliverables:** 32 unit tests in
[`tests/unit/test_preferred_orientation.py`](../../tests/unit/test_preferred_orientation.py); three
executable worked examples (`diffraction-march-dollase-family-factor`,
`diffraction-march-dollase-normalization`, `diffraction-odf-weighted-random-texture`) with
hand-derived values and March/Dollase/Von Dreele/Bunge citations; the LaTeX theory note
[`docs/tex/algorithms/preferred_orientation_in_powder_intensities.tex`](../tex/algorithms/preferred_orientation_in_powder_intensities.tex);
four new validation-matrix rows including an explicit `planned` row for the absorption, extinction
and surface-roughness effects that are *not* modelled; and `describe()` on both models, tested.

Suite: **3749 passed**, ruff and mypy clean.

### RR2c — Complete

Closes critical-review finding 12 and completes the local-misorientation family in `ebsd/`
(KAM, GROD, GOS, GAM → GND): those say *how much* the orientation varies, this says *what
dislocation content that variation implies*, in physical units.

**New module** [`src/pytex/ebsd/gnd.py`](../../src/pytex/ebsd/gnd.py):
`lattice_curvature_tensor`, `nye_dislocation_density_tensor`,
`geometrically_necessary_dislocation_density` (curvature/Nye route and the conventional
`2*theta/(b*u)` KAM route), plus `plot_gnd_density_map`.

**Verification — a planted gradient, not a recorded output.** The test map is built by rotating
the lattice about the specimen normal by a known 0.8°/µm through `Rotation.from_axis_angle`, so
the expectation is independent of the code under test:

| Anchor | Result |
| --- | --- |
| Planted `kappa_20` recovered | 13962.634015954637 rad/m, exact |
| No leakage into other components | < 1e-6 |
| `rho = (dtheta/dx)/b` for copper | 5.4627e13 m⁻², in the expected 1e13–1e14 range |
| Undeformed single crystal | exactly 0 |
| KAM route vs curvature route, pure tilt | agree to 1e-6 — two independent formulas cross-checked |
| Density linear in curvature, inverse in `b`, scales with step | all pinned |

**A real correctness bug found and fixed.** `alpha = kappa^T - tr(kappa)*I` was implemented with
the identity form. Because the out-of-plane curvature column is `NaN` (unmeasurable from a
surface), the trace is `NaN`, and **`NaN * 0` is `NaN` in IEEE** — so the identity term poisoned
every off-diagonal component, destroying exactly the information the Nye tensor exists to carry.
`alpha_02` came out `NaN` instead of the planted value. Fixed by subtracting the trace on the
diagonal directly; a regression test pins it.

**A second defect** caught by warnings-as-errors: `np.log10(x, where=...)` without `out=` reads
uninitialized memory. Fixed with the correct `out=`/`where=` idiom.

**Honesty built into the surface.** Unmeasurable components are `NaN`, never zero — a 2-D map
determines 6 of 9 curvature and 5 of 9 Nye components, and "measured to be zero" must stay
distinguishable from "not measured". Symmetry is deliberately *not* reduced between neighbours
(a disorientation reduction would replace the small physical rotation with an equivalent whose
axis is unrelated to the curvature), so grain boundaries produce artefacts and the docs say so.
Densities are documented as lower bounds and as resolution dependent, with the step-size scaling
pinned as a test so it cannot be silently normalized away.

**Deliverables:** 18 unit tests in
[`tests/unit/test_gnd_density.py`](../../tests/unit/test_gnd_density.py); a new `ebsd` worked-example
group with two examples (`ebsd-planted-lattice-curvature`, `ebsd-gnd-density-from-curvature`)
citing Nye, Pantleon and Kysar; the LaTeX theory note
[`docs/tex/algorithms/lattice_curvature_and_gnd_density.tex`](../tex/algorithms/lattice_curvature_and_gnd_density.tex);
and an MTEX parity row that explicitly **declines** to claim parity with MTEX `calcGND`, which
models dislocation systems rather than a scalar lower bound.

Suite: **3777 passed**, ruff and mypy clean.

### RR3 — Complete

The sweep found that **PyTex could not be installed and imported at all**, which no amount of
in-repo testing would have revealed: the development environment happens to have every optional
dependency present. Verified by building the wheel and importing it in a fresh virtual
environment.

| Defect | Impact | Fix |
| --- | --- | --- |
| `pyyaml` imported at module level but undeclared | **`pip install pytex` then `import pytex` raised `ModuleNotFoundError`** | declared as a runtime dependency; it backs the style themes and the measured-SAED format, both on the core import path |
| `matplotlib.colors` imported at module level in `crystal3d` / `scene3d` | the `plotting` extra was silently mandatory, contradicting the `_require_matplotlib()` guards used everywhere else | made lazy; the library now computes with numpy + scipy + pyyaml alone |
| version literal in **4** files | a release bump would leave manifests stamping the old version | `src/pytex/_version.py` is the single source; `pyproject.toml` reads it statically, `CITATION.cff` gained a matching `version:` |
| `get_phase_fixture` is public but its data is not shipped | a documented public surface failed for every installed user, with a confusing path in the error | added `phase_fixtures_available()`, an actionable error naming the alternatives, and docstring notes |

**Verified end to end in a clean environment** (numpy, scipy, pyyaml only, *no* matplotlib):
`import pytex` succeeds, all 428 exports resolve, `importlib.metadata` reports the right version,
the `pytex` console script runs, and a real Kikuchi simulation plus a March–Dollase-corrected XRD
pattern both compute.

**The ratchet:** [`tests/unit/test_release_metadata.py`](../../tests/unit/test_release_metadata.py)
statically parses every source file and fails on any module-level third-party import that is not a
declared dependency (allowing only `TYPE_CHECKING` blocks and `try`/`except ImportError` guards),
on any second version literal, and on citation metadata drifting out of step. It holds regardless
of what is installed in the environment running it.

Also: CHANGELOG `[Unreleased]` now carries a **Fixed** section leading with the seven correctness
and packaging defects this program found, per the repo rule that scientific behaviour changes be
stated explicitly.

## Release-readiness verdict

**Gates.** 3882 tests pass; ruff and strict mypy clean; wheel and sdist build; clean-environment
install and import verified.

**Known gaps, recorded deliberately rather than left to be rediscovered:**

1. **The phase-fixture corpus is not shipped in the wheel.** Relocating it into the package touches
   23 files including checksum-pinned manifests and `.gitattributes` line-ending rules — the exact
   configuration that previously broke every Windows clone. It should be done as its own change
   with the integrity script re-run, not folded into this program. Until then
   `phase_fixtures_available()` and the loader's error message make the constraint explicit, but
   several existing documentation examples (`saed_generation.md`, `phases_and_cif.md`) still call
   `get_phase_fixture(...)` and will fail for installed users. The new `kikuchi_geometry.md` page
   was written to build its phase explicitly instead.
2. **Version is `0.1.0.dev0` with classifier `Development Status :: 2 - Pre-Alpha`, and there are
   no git tags.** Choosing the first public version number and tagging it is a maintainer decision,
   not an agent one.
3. Hex-grid `CrystalMap` and h5ebsd readers (critical-review finding 11) remain out of scope, as
   recorded at the top of this note.
4. The open follow-ons from the earlier TX and reconstruction programs (measured-EBSD fixtures for
   OR determination, the MATLAB-side MTEX parity run, irregular grain geometry in the map sweep)
   are untouched by this program and still open.

## Previous Task: PyTex scientific capability presentation — COMPLETE (2026-08-04)

`docs/presentations/pytex_scientific_capabilities_2026.pptx` contains 24 slides built from the
17-slide source template, with real PyTex-generated evidence figures (NaCl crystal, HCP
reference frame, orientation/OR mapping, powder XRD, SAED and composite SAED), speaker-source
notes, and a completed render-and-inspect fidelity QA pass. The user asked for the white
scientific infographic style, so the image-generated dark transformation hero was excluded.

## Previous Goal — COMPLETE (2026-08-04)

**Transformation Crystallography And Composite Diffraction Program (TX)** — five user-facing
phase-transformation answers built on the existing OR and diffraction primitives. All phases
TX0–TX6 are delivered and committed to `main`.

| Ask | Delivered surface |
| --- | --- |
| **(a)** measured parent/child Euler angles ⇒ *what is the orientation relationship?* | `characterize_orientation_relationship`, `orientation_relationship_from_euler`, `describe_orientation_relationship` |
| **(b)** an OR + an arbitrary parent plane/direction ⇒ the parallel planes and directions in every product variant | `variant_correspondence_table` |
| **(c)** a parent zone axis ⇒ a composite kinematic SAED exportable as graphics **and** reflection tables | `composite_reflection_table`, `export_composite_saed`, `CompositeSAEDPattern.centering_audit` |
| **(d)** a product-variant zone axis ⇒ the same composite, matrix and siblings around it | `simulate_composite_saed_from_child_zone` |
| **(e)** a measured SAED pattern ⇒ solved, from spots picked interactively or listed in YAML | `solve_saed_pattern`, `MeasuredSAEDPattern`, `SAEDSpotPicker`, `assign_transformation_variant` |

Demonstrated end to end on Burgers β→α in the committed-executed notebook
`docs/site/tutorials/notebooks/23_transformation_crystallography_end_to_end.ipynb`.

### Where the durable records are

- **Specification (normative, now marked delivered):**
  [`docs/architecture/transformation_crystallography_and_diffraction_program.md`](../architecture/transformation_crystallography_and_diffraction_program.md)
- **Phase ledger with the full outcome, the defects found, and the open follow-ons:**
  [`docs/roadmap/working_notes_transformation_diffraction_program.md`](../roadmap/working_notes_transformation_diffraction_program.md)

### Open follow-ons from that program

Measured-EBSD fixtures for the OR determination (validation is synthetic and says so); JSON
round-trip contracts for the new report objects; canonical SVG figures and LaTeX theory notes
for the OR-statement extraction and the ratio/angle algorithm; structure-aware cubic-to-cubic
catalog dispatch. None blocks use of the delivered surfaces. Details in the ledger.

## Older Goal — archived

The reconstruction-stabilization program (Phases 1–5, commits `7dd77d7b` … `c8c6eb1b`) is
complete and its handoff record is archived at
[`docs/development/archive/reconstruction_stabilization_2026_07.md`](archive/reconstruction_stabilization_2026_07.md).
Its three open follow-ons — running the MTEX side of `or_transformation_v1` on a machine with
MATLAB, a measured-data reconstruction fixture, and irregular grain geometry in the map sweep
— are **still open** and remain blockers to moving parent-grain reconstruction out of
`experimental`.
