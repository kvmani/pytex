# World-Class Feature And Foundation Roadmap

This document is the cross-domain umbrella roadmap for making PyTex robust,
scientifically elegant, and world class across texture, EBSD, XRD, TEM, and
orientation-relationship (OR) analysis. It complements — and does not replace —
the MTEX Parity And EBSD Feature Roadmap
(`docs/roadmap/mtex_parity_and_ebsd_feature_roadmap.md`, the EBSD/texture
parity ledger) and the [Implementation Roadmap](implementation_roadmap.md)
(the validation-hardening program). Where those documents track parity and
evidence, this one tracks the feature program: what to build, why, and in what
order.

Strategic decisions governing this roadmap:

1. **Native implementations for heavy numerics.** Rietveld-class XRD
   refinement and dynamical electron diffraction will be built natively
   (pure-Python-first, per `mission.md`), validated against external tools
   (GSAS-II, py4DSTEM-class references) rather than delegated to them. The
   semantic-coherence differentiator applies to intensities, not just
   geometry.
2. **Foundations before breadth.** Every horizon pairs feature work with the
   engineering that keeps the library trustworthy: typed packaging, stability
   policy, re-enabled integrity gates, CI breadth, and performance evidence.

Size tags: **S** (one commit), **M** (one focused session), **L** (multi-session
program with its own working notes).

## 1. Immediate Horizon (next 1-2 cycles)

### 1.1 Engineering robustness

| # | Item | Size | Notes |
| --- | --- | --- | --- |
| I1 | Regenerate stale phase-fixture hashes; re-enable `test_phase_fixtures.py` + `test_repo_integrity.py` | S | `scripts/regenerate_phase_fixture_catalog_hashes.py` exists; removes the permanent 26-test deselect |
| I2 | Ship `py.typed` (PEP 561) | S | strict mypy already passes; downstream users currently get no types |
| I3 | API stability + deprecation policy (`docs/standards/`) with a warn-and-forward helper | S | the ~306-symbol flat `__all__` needs a stability contract before further growth |
| I4 | CI matrix: Python 3.11-3.13, ubuntu+macos base lane; coverage report | S | full-scientific no-skip lane stays pinned on 3.11 |
| I5 | Runnable performance benchmark lane with pinned map-scale cases | M | `benchmarks/` holds evidence manifests only; timing exists only in `scripts/benchmark_miller_vectorized.py` |

### 1.2 Orientation relationships (quick wins, `core/transformation.py`)

| # | Item | Size | Validation |
| --- | --- | --- | --- |
| I6 | Named ORs: Kurdjumov-Sachs, Greninger-Troiano, Pitsch (fcc<->bcc); Burgers (bcc<->hcp, lifts the cubic-only guard) | M | published variant counts: KS=24, NW=12, GT=24, Burgers=12 |
| I7 | Populated standard OR catalog (today `OrientationRelationshipCatalog` ships empty) | S | catalog resolves each named OR by key |
| I8 | Intervariant misorientation table (axis/angle between all variant pairs) | M | KS intervariant angle set vs Morito et al. |
| I9 | Variant-selection scoring: assign observed children to nearest variant; `VariantSelectionReport` gains its algorithm; variant frequency histograms | M | synthetic transformations recover planted variant indices |
| I10 | OR-deviation metric: mean angular deviation of parent/child pairs from a nominal OR | M | zero for exact synthetic data; bridges to MODF-based OR fitting (medium term) |

### 1.3 XRD realism (quick wins, `diffraction/xrd.py`, `physics.py`)

| # | Item | Size | Validation |
| --- | --- | --- | --- |
| I11 | Tabulated X-ray scattering factors (Waasmaier-Kirfel class) and electron factors (Doyle-Turner class), replacing the `f = Z` proxy | M | intensity ratios vs pinned pymatgen baselines (`fixtures/diffraction/`) |
| I12 | K-alpha1/K-alpha2 doublet + richer `RadiationSpec` (Co/Cr/Fe anodes, neutron flag) | S | doublet 2-theta splitting vs Bragg's law analytically |
| I13 | Pseudo-Voigt profile + Caglioti (U,V,W) angular FWHM beside the single-FWHM Gaussian | M | profile integrals and limiting cases (eta=0 Gaussian, eta=1 Lorentzian) |

### 1.4 Visualization (quick wins)

| # | Item | Size | Notes |
| --- | --- | --- | --- |
| I14 | Crystal viewer: periodic boundary atoms (VESTA behavior), crystallographic view presets (along a/b/c/[uvw]), per-species legend | M | `plotting/crystal3d.py`; reuses `_view_angles_from_direction` |
| I15 | ODF sigma sections to pair with phi2 sections | S | `texture/models.py::phi2_sections` is the template |

## 2. Medium Term (following 2-4 cycles)

### 2.1 Texture core

- Kernel library breadth: Gaussian (SO(3)), Abel-Poisson, von Mises-Fisher,
  bump/fibre kernels beside the lone `DeLaValleePoussinKernel`
  (`texture/kernels.py`); `KernelSpec.as_so3_kernel()` generalized.
- `SO3FunHarmonic` at full strength: quadrature from scattered orientations,
  convolution, Funk/PF projection; replace `HarmonicODF` internals.
- MTEX-grade PF->ODF inversion: zero-range method, ghost correction,
  defocusing models on `PoleFigureCorrectionSpec`, validated against the
  LaboTex/XRDML corpus.
- `BinghamODF` component + parameter estimation; MDF axis/angle marginals
  (MDF core exists in `core/misorientation_distribution.py`).

### 2.2 XRD toward quantitative analysis (native path)

- Background models (Chebyshev, spline) + multi-phase pattern summation with
  scale factors.
- Preferred-orientation intensity corrections: March-Dollase and ODF-weighted
  (a unique PyTex strength: the texture core drives powder intensities).
- Size/strain broadening (Scherrer, Williamson-Hall) and instrument-profile
  convolution.
- Pattern export (`.xy`, XRDML write) and CIF **export** (imports exist via
  pymatgen; no writers today).
- Whole-pattern least squares (Pawley/Le Bail) as the stepping stone to
  Rietveld.

### 2.3 TEM / SAED

- **Kikuchi/gnomonic geometry surface**: `KikuchiBand`/`KikuchiPattern`,
  gnomonic projection, pattern-center-aware overlays (scoped in
  `references/feature_opportunities.md` section 3; highest-priority TEM item).
- Ring/polycrystalline (Debye-Scherrer) SAED with texture-weighted ring
  intensities (couples the ODF into TEM).
- Double diffraction / forbidden-reflection flagging; relativistic wavelength
  correction surfaced in `DiffractionGeometry`.
- HOLZ-line geometry; dark-field/selected-spot simulation from the existing
  `KinematicSimulation.simulate_spots`.

### 2.4 OR / parent-grain reconstruction (the flagship)

- Map-scale **variant-graph parent reconstruction** (MTEX
  `parentGrainReconstructor` equivalent): bridge `core/transformation.py`
  variants with the `ebsd/models.py` grain-boundary network
  (`merge_by_csl` union-find is the reusable grouping machinery);
  boundary-misorientation voting, iterative parent-grain growth, per-grain
  fit/confidence maps; validated on martensite->austenite and alpha->beta Ti
  literature fixtures.
- MODF-based OR refinement: fit the operative OR from measured parent/child
  boundary misorientations.
- Twin-aware parent voting reusing `ebsd/csl.py` twin laws.
- Variant pole figures + variant-colored maps (plotting support).

### 2.5 EBSD and scale

- Hex-grid `CrystalMap` (offset rows, honeycomb neighbors); h5ebsd-family
  readers (EDAX/Oxford/Bruker H5) behind one HDF5 core.
- Denoising suite as `CrystalMapFilter` policy objects (spline and
  half-quadratic smoothing, Kuwahara, non-indexed in-fill).
- GND density from orientation gradients (curvature-tensor route).
- Out-of-core groundwork: chunked/memmapped `CrystalMap` backing for
  larger-than-RAM maps (everything is in-RAM NumPy today).

### 2.6 Visualization / GUI preparation

- Pole-figure/IPF publication upgrade on the `pytex.texture` colormap with
  m.r.d. colorbars and RD/TD rim annotations; IPF color-key legend panel.
- Interactive-layer groundwork: map exploration widgets (hover orientation,
  click-to-PF) on the Matplotlib runtime doctrine; the
  `CrystalScene`/figure-spec scene-graph split is the GUI contract.

## 3. Long Term (world-class differentiation)

- **Native Rietveld refinement** + quantitative phase analysis, built on the
  medium-term profile/background/Pawley stack; validated against GSAS-II.
- **Dynamical electron diffraction**: Bloch-wave (then multislice)
  intensities, CBED, thickness oscillations; validated against
  py4DSTEM/JEMS-class references.
- Habit-plane and interface crystallography: habit-plane determination from
  grain pairs (the `TransformationVariant.habit_plane_pairs` slot exists with
  zero computation today), invariant-line/PTMC analysis, five-parameter
  boundary character.
- 3D EBSD: `CrystalMap3d`, 3D grains and boundary meshes, serial-section
  registration.
- Spherical/dictionary indexing bridges; spot/band indexing beyond the
  current geometric scaffolding.
- Texture evolution simulation (Taylor-based rolling texture prediction; the
  `properties/taylor.py` LP solver is the seed).
- Advanced statistics: Bingham confidence regions, bootstrap uncertainty on
  ODFs and volume fractions.
- GPU/numba acceleration behind pure-Python fallbacks; streaming import for
  multi-GB vendor files.
- The generic GUI tool consuming the scene-graph/figure-spec layer.

## 4. Sequencing Rationale

The immediate horizon removes trust debt (hash gates, typing, stability
policy) while landing the highest value-per-cost science: named ORs and
variant analysis unlock steel/titanium transformation studies with machinery
that already exists (`generate_variants`); scattering tables and profile
shapes are prerequisites for every later quantitative XRD claim; the crystal
viewer and section plots finish surfaces users touch daily. The medium term
builds the two flagships (parent reconstruction; quantitative XRD path) plus
the TEM geometry layer, each on top of an immediate-horizon foundation. The
long term is differentiation that assumes those flagships are stable.

## References

### Normative

- `mission.md`
- `specifications.md`
- [Implementation Roadmap](implementation_roadmap.md)
- `docs/roadmap/mtex_parity_and_ebsd_feature_roadmap.md`

### Informative

- `references/feature_opportunities.md`
- `docs/roadmap/working_notes_plotting_foundation.md`
