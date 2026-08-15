# Feature Capability Review (2026-08)

A workflow-oriented audit of the working tree against the nine niche capabilities
users actually come to PyTex for, plus the cross-cutting concerns that decide
whether those capabilities are usable. Scores are *usability by a researcher for
real work*, not lines of code.

Verified facts at the 2026-08-15 reconciliation: 138 Python source files, about
82.5k Python source lines, 558 symbols in `pytex.__all__`, 2,059 test functions
before parametrization, 31 executable tutorial notebooks, and 28 worked-example
sources. Strict mypy and Ruff lint are clean, `py.typed` ships, and the CI matrix
covers Ubuntu and macOS on Python 3.11–3.13 with an 87% coverage ratchet plus a
full-scientific lane. Windows and documentation-warning ratchets are the active
governance increment.

## Scorecard

| # | Capability | Score | One-line verdict |
| --- | --- | --- | --- |
| 1 | SAED indexing / solving | 8 | Complete geometric solver; blind to intensities and to raw images |
| 2 | Composite OR patterns | 9 | Best-in-class; not quantitative across phases |
| 3 | SAED for arbitrary orientation / zone axis | 8 | Correct kinematics; no thickness shape factor, no HOLZ, no dynamical |
| 4 | Tilt solving to a target zone axis | 9.5 | The strongest subsystem in the repo |
| 5 | ODF from XRDML / LaboTex / EBSD | 8 | Full stack minus ghost correction and a defocus model |
| 6 | Pole-figure arithmetic | 2 → **8** | Was the clearest hole; addressed in full, see below |
| 7 | XRD pattern computation | 7 | Excellent forward model, zero analysis of measured data |
| 8 | EBSD | 8.5 | Broad and deep; square plus `.ang` hex-grid topology, but no HDF5 readers |
| 9 | OR from two grains' Euler angles | 9 | Exactly the asked-for entry point, with honest ambiguity reporting |
| 10 | Visualization | 8 | Superb architecture, static-matplotlib ceiling |

---

## 1. SAED indexing / solving — 8/10

**Implemented.** `diffraction/solving.py::solve_saed_pattern` — classical
ratio/angle indexing from picked spots over multiple candidate phases, with
space-group absences, per-spot `hkl` assignment, residuals, deduplicated ranked
solutions, and an explicit `is_conclusive` verdict; the ±zone-axis ambiguity of a
centrosymmetric reflection set is stated rather than hidden.
`solve_saed_pattern_file` plus the `measured_saed_pattern` YAML contract and JSON
schema make a solve reproducible from a committed text file.
`diffraction/models.py::index_saed_pattern` is the calibrated-`DiffractionGeometry`
route; `estimate_zone_axis` and `FamilyIndexingReport` sit beside it.
`plotting/saed_picker.py` separates testable picking *logic* from the Matplotlib
event adapter. `tem/indexing.py` lifts an indexed pattern to a full 3-D
orientation (single and multi-pattern, with scatter statistics), and
`tem/ambiguity.py` enumerates Laue-coset ambiguity families and the experiments
that would discriminate them.

**Lacking.** Intensities are never used — no intensity-weighted tie-breaking and
no double-diffraction / forbidden-spot flagging. No spot *detection*: patterns
must be clicked or typed, so there is no image → `MeasuredSpot` path (no
centroiding, subpixel refinement, or background removal). No ring
(Debye-Scherrer) indexing for polycrystalline or nanocrystalline patterns. No
HOLZ use. No lattice-parameter refinement from a solved pattern, and no search
against a structure database when the phase is unknown.

**Next.**
1. Automated spot detection with subpixel centroids, feeding `MeasuredSAEDPattern` directly.
2. Intensity-aware ranking plus a double-diffraction flag on formally forbidden hits.
3. Ring-pattern indexing (radius histogram → phase/`d`-spacing match).

## 2. Composite patterns for a given OR — 9/10

**Implemented.** `diffraction/composite.py::simulate_composite_saed` — parent plus
any subset of variants on one shared, parent-anchored detector basis, with the
geometry conventions pinned in the module docstring; irrational child zone axes
handled exactly and rationalized with a reported deviation.
`simulate_composite_saed_from_child_zone` inverts the entry point,
`find_spot_coincidences` / `SpotCoincidenceReport` answer the question the
technique exists for, and `sweep_parent_zone_axes` finds the most discriminating
zone axis. Export is a first-class product: `composite_reflection_table`,
CSV/markdown, and a manifest with a JSON schema. `plotting/composite_saed.py`
renders it with label clustering and declutter.

**Lacking.** Parent and child intensities are each normalized to 1 per
sub-pattern — honest, but it means a composite cannot be compared to a
micrograph quantitatively. No double diffraction between phases (`g_parent +
g_child`), which is a real and often dominant feature of martensite and
precipitate patterns. No shape-factor anisotropy for plate/lath variants. No
objective-aperture / dark-field selection.

**Next.**
1. Cross-phase intensity scaling under a declared volume-fraction and thickness model.
2. Double-diffraction spot generation with provenance marking each such spot.
3. Relrod elongation from variant habit-plane morphology.

## 3. SAED for arbitrary crystal orientation and zone axis — 8/10

**Implemented.** `diffraction/kinematic.py::simulate_zone_axis_spots` is the
engine: fully vectorized, selects reflections by Ewald-sphere excitation error
`s_g = g_z - g^2*lambda/2`, applies centering absences, computes electron
structure factors with isotropic Debye-Waller damping, an exact optional plane-parallel
`sinc^2(t s_g)` finite-thickness factor (with the older Lorentzian retained as an explicitly
separate compatibility proxy), `g_max`/`max_index` cuts, and a deterministic right-handed
detector triad with `align_g` and in-plane rotation control. It accepts an
irrational `CrystalDirection`, so *any* beam direction — not just a rational zone
— is exact. `diffraction/saed.py::generate_saed_pattern` remains as the simple
zone-law route, and `DiffractionGeometry` / `KinematicSimulation` provide the
calibrated detector projection (pattern centre, non-square pixels, detector
tilts, acceptance masks).

**Lacking.** The plane-parallel finite-thickness shape factor is implemented and independently
pinned at its analytic landmarks, but no bending, thickness distribution, mosaic spread, or
absorption broadens it. Higher-order Laue zones are excluded by the default `s_g` window and
there is no Laue-zone index on spots; existing HOLZ-ring geometry is not integrated into each
spot table. The separate CBED/Bloch-wave foundations do not yet make this SAED surface a general
dynamical simulation, and precession is absent.

**Next.**
1. Add bending/thickness-distribution specimen envelopes without weakening the explicit slab model.
2. Tag every spot with its Laue-zone index and connect the existing HOLZ-ring geometry.
3. Expand two-beam / Bloch-wave intensity integration and external validation.

## 4. Tilts required to reach a target zone axis — 9.5/10

**Implemented.** This is the most complete subsystem. `tem/reconstruction.py::CurrentState`
builds the crystal-to-holder orientation from indexed patterns, and
`CurrentState.from_two_zone_axes` deliberately avoids the diffraction-rotation
calibration that most often sends operators tilting backwards. `tem/stage.py`
models double-tilt, tilt-rotate and single-tilt holders, non-ideal/oblique axes
(`GeneralStageAxes`, `StageCalibration`) and rectangular / elliptical / polygon /
masked envelopes. `solve_tilts_for_direction` gives the closed-form answer;
`plan_tilt_to_zone_axis` searches the full symmetry orbit and every ambiguity
family, re-derives each candidate through the calibrated forward model and
rejects on residual *before* ranking, then scores on travel, largest single-axis
move, envelope clearance, propagated uncertainty and conditioning. `tem/path.py`
plans the route (geodesic, two-leg, waypoints, `connecting_band`);
`tem/calibration.py` fits the stage from tilt excursions and persists it; and
`plotting/tilt_stereogram.py` renders the whole answer as a publication figure.
Reports carry `describe()` prose and name the orbit member actually landed on.

**Lacking.** No backlash / hysteresis model. No eucentric-height or beam-shift
coupling. No reachability *atlas* product (bulk sweeps are possible via
`include_paths=False` but there is no exportable "what can I reach from here"
table). Kikuchi bands are simulated but not wired into navigation, so band-guided
tilting is unavailable. No live microscope-control bridge (defensible, but should
be stated).

**Next.**
1. A reachability atlas: per phase and holder, an exportable table/figure of reachable zones.
2. Couple `simulate_kikuchi_pattern` into navigation for band-guided tilt paths.
3. Backlash-aware path strategy with an approach-direction convention.

## 5. ODF computation from XRDML, LaboTex and EBSD — 8/10

**Implemented.** Two ODF representations: discrete kernel-weighted `ODF` and
`HarmonicODF` (generalized spherical harmonics on a symmetry-projected,
weighted-orthonormalized basis) with texture index, entropy, pole-density
evaluation, pole-figure reconstruction and Tikhonov-regularized inversion. Kernel
library: de la Vallee Poussin, Gaussian SO(3), Abel-Poisson.
`texture/reconstruction.py::ODFReconstructionConfig` bundles correction,
algorithm, kernel, regularization and bandlimit into one declared object so the
correction cannot be skipped; `PoleFigureCorrectionSpec` applies background,
defocus and scale in a fixed, documented order;
`defocus_from_random_standard` derives same-reflection radial factors, retains
azimuthal-scatter diagnostics, and refuses extrapolation;
`residual_reports_for_pole_figures` closes the QC loop. Sections (`phi2_sections`,
`sigma_sections`), `volume_fraction`, `TextureComponent` / `Fibre` with standard
fcc/bcc rolling components, and `component_volume_fractions` cover the reporting
side. Ingest: `read_xrdml_pole_figure` / `invert_xrdml_pole_figures`,
`read_labotex_pole_figures` / `invert_labotex_pole_figures`, and
`CrystalMap.to_odf` for EBSD. `core/misorientation_distribution.py` provides the
MDF core.

**Lacking.** No ghost correction and no zero-range method — the odd part of the
ODF is unconstrained, and the docstrings say so. No inverse fitting of components
(Gauss/Bingham) to a measured ODF.
No uncertainty quantification (bootstrap on ODF or volume fractions). No `.epf`,
`.uxd` or popLA `.xpc` readers.

**Next.**
1. Ghost correction (zero-range / positivity) — the biggest scientific credibility gap in texture.
2. ~~A defocus model plus random-standard calibration~~ — implemented as
   `defocus_from_random_standard(...)`; broader experimental fixtures remain desirable.
3. Component fitting: recover Gauss components and volume fractions from a measured ODF.

## 6. Pole-figure arithmetic — 2/10 at review, **8/10 after the 2026-08-08 sprint**

**At review.** Essentially nothing, and blocked structurally rather than merely
unwritten. `PoleFigure` held scattered `sample_directions` and `intensities`
with only `from_orientations`, `project` and `histogram`; two figures generally
shared no sampling direction, so there was nothing to combine pointwise. No
resampling, no spherical interpolation. Measured intensities were normalized
only by `max` or `sum`, never to multiples of a random distribution, so even on
a shared support their magnitudes were not comparable. No symmetrization, no
rotation, no tilt-range mask, and no difference product — so ODF-inversion QC
had no visual counterpart either.

**Delivered** (`ee9591c`..`920e6c4`; see
`docs/development/active_task_progress.md` for the design record):

1. `PoleFigure.sampling` — records whether the intensities are per-pole weights
   or an evaluated density. The distinction was invisible until resampling and
   then decides the answer, so it is recorded rather than guessed.
2. `PoleFigure.on_grid` — kernel resampling onto any `S2Grid`, estimator chosen
   from the sampling tag. Supplies the shared support.
3. `PoleFigure.normalize_to_mrd` / `spherical_mean`, plus
   `raster_solid_angle_weights` and an `mrd` mode on the XRDML and LaboTex
   readers. Supplies the shared scale.
4. `__add__` / `__sub__` / `__mul__` / `__truediv__`, `difference`, `rotate`,
   `symmetrize`, `restrict_polar_range`, each raising on pole, frame, antipodal,
   family or support mismatch. Subtraction returns a signed
   `PoleFigureDifference`, because a pole density is non-negative and a
   difference is not.
5. `PoleFigureResidualReport.difference_figure()` and
   `plot_pole_figure_difference` — the residual figure ODF inversion never had.

**Still lacking.** The resampling kernel is one fixed von Mises-Fisher shape,
with no S2 kernel library matching the SO(3) one. m.r.d. over a partial figure
averages over the measured cap only. There is no contoured rendering of a
difference figure (it is drawn as a scatter, which is honest for a scattered
support but less readable than a contour on a dense grid). Ratio figures have no
masking helper for near-zero denominators beyond refusing outright.

## 7. XRD pattern computation — 7/10

**Implemented.** `diffraction/xrd.py` generates powder reflections and full
patterns with tabulated X-ray form factors (`diffraction/scattering.py` plus a
generated table, replacing the old `f = Z` proxy), structure factors,
multiplicities from the phase symmetry, Lorentz-polarization, Debye-Waller,
Gaussian / Lorentzian / pseudo-Voigt profiles, Caglioti `U,V,W` angular FWHM, and
Kalpha1/Kalpha2 doublets for Cu, Mo, Co, Cr and Fe with a neutron flag and an
explicit warning that X-ray form factors do not transfer to neutrons. Systematic
absences come from the declared space group. `apply_preferred_orientation`
supports March-Dollase *and* an ODF-driven model
(`ODFPreferredOrientationModel`) — texture driving powder intensities is a
genuine differentiator no comparable package offers so directly.

**Lacking.** The module is forward-only. There is no measured-pattern ingest at
all: the XRDML adapter reads pole figures, not scans, and there is no `.xy`,
`.raw` or `.xrdml`-scan reader and no `MeasuredPowderPattern` type. Consequently
there is no background model, no multi-phase summation with scale factors, no
peak search or single-peak fitting, no Le Bail / Pawley / Rietveld, no
quantitative phase analysis, no size-strain (Scherrer, Williamson-Hall), no
instrument-profile convolution, no sin^2(psi) residual stress, and no pattern
export.

**Next.**
1. Measured-pattern I/O and a `MeasuredPowderPattern` type — nothing downstream can start without it.
2. Chebyshev background plus multi-phase summation with refinable scale factors.
3. Le Bail / Pawley whole-pattern fitting as the stepping stone to Rietveld and QPA.

## 8. EBSD — 8/10

**Implemented.** `ebsd/models.py::CrystalMap` covers multiphase maps, arbitrary
property channels, phase masks and selection, point selection, and manifest
export. `GrainSegmentation` provides union-find segmentation, GOS/GAM/GROD maps,
fitted ellipses, aspect ratios, shape orientation, perimeter, area, shape factor
and bounding boxes; `GrainBoundaryNetwork` and `GrainGraph` give the boundary
topology. Analytics: KAM at arbitrary neighbour order, CSL classification with
Brandon criterion and twin laws (`ebsd/csl.py`), GND density via lattice
curvature → Nye tensor (`ebsd/gnd.py`), Schmid and Taylor factor maps. Cleanup:
`remove_wild_spikes`, `majority_smoothed`, `merge_small_grains`, property
thresholds. Texture: `to_odf`, `pole_figure`, `inverse_pole_figure`,
`texture_report`, plus the `EBSDTextureWorkflow`. Ingest: native `.ang` and `.ctf`
readers (`adapters/scan_files.py`) and kikuchipy / pyebsdindex / orix adapters
(`normalize_ebsd`, `index_hough`, `refine_orientations`). Plotting covers IPF,
KAM, GND, phase and property maps with boundary overlays. MTEX parity tests guard
rotations, fundamental regions and EBSD behaviour.

**Lacking.** `.ang` `HexGrid` import, six-neighbour KAM, segmentation, explainable metadata, and
portable contracts are now implemented with a labelled analytic fixture. Hexagonal curvature/GND
finite differences and pixel-cell perimeters remain deliberately unsupported. No HDF5 readers
(`.h5oina`, `.oh5`, Bruker H5), which is how modern systems
actually export. Native pattern indexing does not exist: Hough and dictionary
indexing are delegated to optional adapters, so the pure-Python-first claim does
not hold for indexing. No in-fill of non-indexed points and no pluggable filter
policy layer. Everything is in-RAM NumPy — no chunked or memory-mapped backing,
so large maps will fail rather than degrade. No 3D EBSD, no five-parameter
boundary character.

**Next.**
1. Add hexagonal finite-difference curvature/GND and a declared cell-boundary perimeter model.
2. Add one HDF5 core with `.h5oina` and `.oh5` readers.
3. Add non-indexed in-fill and a `CrystalMapFilter` policy layer; then native Hough indexing to close the pure-Python claim.

## 9. OR between two grains from Euler angles — 9/10

**Implemented.** `core/transformation.py::orientation_relationship_from_euler` is
literally the asked-for entry point: two columns of Euler angles from an EBSD
export, two phases, and it returns an `ORCharacterizationReport` via
`characterize_orientation_relationship` — catalog matching, plane and direction
parallelism statements with residuals, and a `describe()` prose summary.
`fit_orientation_relationship` performs a double-coset-seeded fit; `or_deviation`
quantifies departure from a nominal OR. Without any parent data,
`experimental/or_identification.py::identify_orientation_relationship` ranks
candidates from child-child boundary misorientations against intervariant
fingerprints and states plainly when the margin is too small to discriminate, and
`or_refinement.py` refines the rotation itself from boundaries. Catalogs ship
Bain, KS, NW, Greninger-Troiano, Pitsch, Burgers, Shoji-Nishiyama,
Pitsch-Schrader and ferrite-cementite. Downstream: variants, intervariant
misorientations and boundary fingerprints, `variant_correspondence_table` (with
`to_csv`), `map_direction_across_variants` / `map_plane_across_variants`,
`variant_pole_figure` and its plot, and parent-grain reconstruction at both
single-grain and map scale.

**Lacking.** No uncertainty on a fitted OR — no covariance or confidence interval
propagated from orientation noise, and the discriminability margin is reported
but not formalized as a test. Habit-plane determination is a declared slot
(`TransformationVariant.habit_plane_pairs`) with zero computation behind it; no
PTMC or invariant-line analysis. Parent/child pairs must be supplied already
matched — there is no automatic extraction of adjacent cross-phase grain pairs
from a `CrystalMap`.

**Next.**
1. Bootstrap confidence intervals on the fitted OR and a formal candidate-discriminability test.
2. Automatic parent/child pair extraction from the boundary network, closing the map → OR loop.
3. Habit-plane determination from trace analysis, then PTMC.

## 10. Visualization — 8/10

**Implemented.** The architecture is the strongest part: a declarative figure-spec
layer (`plotting/builders.py`, `_render.py`, typed layer dataclasses) rendered by
a Matplotlib runtime, a headless SVG path (`svg_primitives`, `_svg_text`, with
text-overflow and marker audits in `scripts/`), and a 3-D scene graph
(`CrystalScene`, `PrimitiveScene3D`, `WorldScene3D`). Products cover
stereographic and Wulff nets, crystal directions/planes/symmetry elements, pole
figures (scatter and contoured), IPF plus an IPF colour key, ODF phi2 and sigma
sections, orientation/rotation/Euler/quaternion clouds, EBSD IPF/KAM/GND/phase/
property maps, XRD, SAED and Kikuchi patterns, composite SAED with label
declutter, tilt stereograms, 3-D crystal structures with plane, direction and
polyhedron overlays, reference-frame and algorithm diagrams, elastic property
surfaces, panel grids, scale bars, YAML style themes, a purpose-built m.r.d.
colormap and publication export. Structural validation tests guard the figures.

**Lacking.** Static Matplotlib only — no interactive or web layer, which is what
"stunning" now implies; the figure-spec indirection exists precisely to allow a
second backend, and that promise is unredeemed. 3-D quality is capped by
Matplotlib's painter's algorithm, so VESTA-class crystal rendering is not
reachable on this path. No animation, though tilt paths and variant sweeps
already produce the required frame data. Pole-figure contouring goes through a
2-D `histogram2d`, which is blocky at low counts where a spherical kernel density
would be smooth. Pole-figure colorbar and RD/TD rim polish remains on the
roadmap.

**Next.**
1. Spherical kernel density for PF/IPF contouring — an immediate, visible quality jump.
2. An animation helper (tilt paths, variant sweeps, ODF sections) emitting GIF/MP4.
3. An optional interactive backend behind the existing figure-spec contract; the architecture already anticipates it.

---

## 11. Cross-cutting findings not in the original list

**Strengths worth protecting.** A large warning-strict test suite, strict mypy, Ruff, `py.typed`,
JSON schemas and manifests for every serialized contract, provenance records
threaded through the domain objects, a worked-example framework, MTEX/VESTA/
diffraction parity ledgers, a terminology registry and an enforced citation
policy. Very few academic packages have this governance layer; it is the
repository's real moat.

**Risk 1 — roadmap drift (governance increment in progress).**
`docs/roadmap/world_class_feature_roadmap.md` and
`docs/roadmap/critical_review_and_development_guide.md` (dated 2026-07) describe a
31.7k-line, 804-test repository with CI on Python 3.11/ubuntu only and no
coverage. The tree is now 66.1k lines with 4281 tests, an ubuntu+macos x
3.11-3.13 matrix, and an 87% coverage ratchet. Items those documents list as
missing have landed: named ORs, the populated catalog, intervariant tables,
variant selection, `or_deviation`, the kernel library, GND, sigma sections,
scattering tables, doublets, pseudo-Voigt, `py.typed`. A contributor reading them
today will be misdirected. Reconcile them, or mark them superseded by this
review.

**Risk 2 — API surface.** `pytex.__all__` exports 470 flat symbols. A stability
policy (`docs/standards/api_stability_and_deprecation.md`) and a `_deprecation`
helper exist, but a flat namespace this size is hard to learn and hard to evolve.
Surface stability tiers (stable / provisional / experimental) *in the API itself*
and add a task-oriented "how do I…" index keyed to the nine workflows above.

**Risk 3 — no Windows in CI (governance increment in progress).** The matrix is Ubuntu and macOS; primary
development happens on Windows. Path, encoding and line-ending regressions will
be found by the developer rather than by CI. Adding `windows-latest` to the base
lane is a one-line change.

**Risk 4 — performance has no evidence.** Everything is in-RAM NumPy and
`benchmarks/` holds evidence manifests rather than timings. There is no CI timing
lane, so a scale regression on map-scale EBSD or ODF inversion would go unnoticed.

**Risk 5 — adoption and discoverability trail the API.** Thirty-one notebooks and 28 worked
examples now provide substantial teaching coverage, but 558 flat stable exports remain difficult
to discover. A task-oriented index and stability tiers are now more valuable than simply adding
more examples by count.

## Recommended next sprint

Ordered by payoff per unit effort, not by ambition:

1. ~~**Pole-figure arithmetic sprint** (section 6)~~ — **done 2026-08-08.** Resampling, m.r.d. normalization, arithmetic and residual figures all landed; see section 6.
2. **Governance repair** — reconcile the roadmap/ledger, add `windows-latest`, ratchet Sphinx warnings, and automate critical browser behavior.
3. **Measured XRD pattern I/O** (section 7) — converts the forward model into an analysis capability.
4. **Random-standard defocus calibration** (section 5) — completes the correction input already accepted by `PoleFigureCorrectionSpec`; ghost correction follows as a separate, higher-risk algorithm.
5. **Hex-grid EBSD support** (section 8) — removes an explicit rejection that excludes real EDAX data.
6. **Finite-thickness SAED and the true shape factor** (section 3) — one parameter, materially better SAED intensities, and the prerequisite for later quantitative work.
7. **Named-component ODF fitting** — reuses the existing component library and non-negative SciPy solvers to turn an ODF into interpretable component fractions plus a residual.
