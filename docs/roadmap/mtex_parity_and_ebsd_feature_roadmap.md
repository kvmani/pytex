# MTEX Parity And EBSD Feature Roadmap

This document is the detailed feature program for taking PyTex from its current validated
foundation to a comprehensive texture, microtexture, and crystallography analysis tool with MTEX
parity as the floor. It is organized in three horizons — immediate, medium term, and long term —
and it is deliberately foundational-class-first: each horizon names the concrete classes and
methods that must be created or upgraded, because downstream feature breadth is only sustainable
when the semantic spine underneath it is complete.

Companion documents:

- `implementation_roadmap.md` — governing phase discipline and release gates
- `../testing/mtex_parity_matrix.md` — the authoritative parity ledger
- `../../references/feature_opportunities.md` — reference-corpus-backed feature ideas

## 1. Current Foundation Assessment

### 1.1 What already exists and is strong

| Layer | Surface | Status |
| --- | --- | --- |
| Rotations and orientations | `Rotation`, `Orientation`, `Misorientation`, `OrientationSet` in `core/orientation.py`: Euler (Bunge/Matthies/αβγ), quaternion, matrix, axis-angle, Rodrigues/Rodrigues-Frank, Miller and plane-direction construction, SO(3) grids, symmetry-aware reduction and misorientation | strong foundation |
| Symmetry | `SymmetrySpec`, `FundamentalSector` in `core/symmetry.py`: proper point-group generation, orbit-based exact reduction, IPF sector logic | foundational, proper-rotation-only |
| Miller objects | `MillerPlane(Set)`, `MillerDirection(Set)`, Miller-Bravais wrappers, vectorized family expansion, angle/projection helpers in `core/miller.py` | strong foundation |
| Structure | `Phase`, `Lattice`, `SpaceGroupSpec`, `UnitCell`, metric tensors, CIF via pymatgen adapter in `core/lattice.py` | strong foundation |
| Batches | `VectorSet`, `EulerSet`, `QuaternionSet`, `RotationSet` in `core/batches.py` | foundational |
| Texture | `PoleFigure`, `InversePoleFigure`, discrete kernel `ODF`, `HarmonicODF` (band-limited Wigner basis), dictionary-based PF inversion in `texture/` | foundational |
| EBSD | `CrystalMap`, `CrystalMapPhase`, `GrainSegmentation`, `GrainBoundaryNetwork`, `GrainGraph`, KAM, GROD, cleanup, multiphase selection in `ebsd/models.py` | foundational, regular-grid-first |
| Diffraction | powder XRD, SAED, structure factors, zone-axis estimation, indexing scaffolding in `diffraction/` | foundational |
| Plotting | `IPFColorKey`, IPF/KAM maps, Wulff nets, stereographic surfaces, VESTA-like 3D crystal scenes, YAML styles | foundational |
| Interchange | JSON contracts, import/experiment/validation manifests, XRDML and LaboTex PF import, orix/KikuchiPy/PyEBSDIndex bridges | strong differentiator |

### 1.2 Structural gaps that block MTEX-class breadth

These are the gaps that matter most, because many MTEX features are thin wrappers over
capabilities PyTex does not yet have:

1. **Symmetry is proper-rotation-only.** `SymmetrySpec` maps the 32 point groups onto 11 proper
   groups. There is no first-class representation of improper operations (mirrors, inversion,
   rotoinversion), Laue classes as objects, or specimen symmetry. MTEX leans on full
   `crystalSymmetry`/`specimenSymmetry` objects for family expansion, antipodal logic, color keys,
   structure-factor extinctions, and ODF symmetrization.
2. **No S2 (sphere) function layer.** MTEX's `vector3d` + `S2Fun` pairing underlies pole figures,
   IPF keys, and defocusing corrections. PyTex has `VectorSet` but no spherical-coordinate
   semantics, no spherical harmonics, and no spherical kernel density surface.
3. **No SO(3) function algebra.** MTEX 5+ is architecturally an `SO3Fun` library: harmonic, RBF,
   Bingham, fibre, and uniform ODFs are all one algebra with `+`, `*`, rotation, convolution, and
   quadrature. PyTex has two disconnected ODF classes (discrete kernel and band-limited harmonic)
   and no component algebra.
4. **No fibre object.** Fibre textures (α, β, γ fibres) are first-class in MTEX (`fibre` class,
   `fibreODF`, `volume` along fibres) and central to rolling-texture analysis.
5. **No named texture-component vocabulary.** Cube, Goss, Brass, Copper, S, and named fibres are
   the daily language of texture analysis and are absent.
6. **Vendor-reader breadth remains limited.** Direct square/hex `.ang` and `.ctf` readers now
   exist; `.cpr/.crc`, `.osc`, and h5ebsd-family formats still require adapters or remain absent.
7. **Grain geometry remains segment-based.** Rectangular and hexagonal graph topology, shape
   descriptors, CSL/twin classification, and MDF now exist, but boundaries are still discrete
   measurement-center segments rather than reconstructed physical polylines; hex cell perimeters
   and curvature/GND stencils remain open.
8. **Kernel family is minimal.** MTEX's de la Vallée Poussin kernel with halfwidth ↔ bandwidth
   duality is the workhorse for both ODF estimation and PF inversion; `KernelSpec` does not yet
   model kernels with Chebyshev/Legendre coefficient expansions.
9. **No physical-property layer.** Tensors, Schmid/Taylor factors, and slip systems are entirely
   absent; these anchor the "texture to properties" half of MTEX.

## 2. Foundational Class And Method Program

This section is the heart of the roadmap: the classes and methods that must exist (or be
upgraded) for the horizons in Sections 3–5 to be buildable without semantic forks. Each item
lists the target module and the API shape at the level of intent, not final signature.

### 2.1 Symmetry completion (`core/symmetry.py`, new `core/point_groups.py`)

Create or update:

- `PointGroup` (new, or a major `SymmetrySpec` upgrade): full international (Hermann-Mauguin) and
  Schoenflies naming for all 32 crystallographic point groups; operator set including improper
  operations as 3×3 matrices with determinant −1; properties `proper_subgroup`, `laue_group`,
  `is_laue`, `order`, `is_proper`; iteration over `rotations`, `mirrors`, `inversion`.
- `LaueClass` semantics: `SymmetrySpec.laue()` returning the 11 Laue groups as first-class values
  so antipodal handling stops being an ad hoc boolean.
- `SpecimenSymmetry` (new): triclinic, monoclinic, orthotropic, axial, and isotropic specimen
  symmetry with the same operator interface, consumed by `Orientation.canonicalize`,
  ODF symmetrization, and PF plotting.
- `FundamentalSector` upgrade: exact closed-form sector polygons per Laue class (vertices, edge
  great circles) rather than wedge tests, so IPF color keys and sector plots share one geometry.
- `OrientationFundamentalRegion` upgrade (`core/orientation_geometry.py`): closed-form
  Rodrigues-space polytopes per proper group (the MacKenzie cells), `contains()`, `boundary()`
  mesh export for plotting, and maximum disorientation angles per group pair.

### 2.2 Spherical vector and function layer (new `core/sphere.py`, new `texture/s2fun.py`)

- `SphericalVectorSet` (upgrade of `VectorSet` or a sibling): polar constructors
  `from_polar(theta, rho)`, `theta`/`rho` properties, `antipodal` flag carried as semantics,
  `angle_to`, `cross`, `dot`, `mean` (spherical mean), `perp`, and region masks. This is the
  `vector3d` equivalent everything else consumes.
- `S2Grid` (new): equispaced and Gauss-Legendre sphere grids for quadrature and plotting,
  with hemisphere/sector restriction.
- `S2FunHarmonic` (new): real spherical-harmonic expansion with quadrature from scattered data,
  evaluation, `+`, `-`, `*`, rotation, and symmetrization by `PointGroup`. Pole figures,
  IPF densities, and defocusing corrections become instances of this class.
- `S2Kernel` family (new): Abel-Poisson, von Mises-Fisher, Gauss-Weierstrass kernels with
  Legendre-coefficient representations for spherical kernel density estimation.

### 2.3 SO(3) function algebra (new `texture/so3fun.py`, refactor of `texture/models.py` and `texture/harmonics.py`)

- `SO3Fun` (new abstract base): `evaluate(orientations)`, `mean()`, `max()` / `steepest_descent_modes()`,
  `entropy()`, `texture_index()`, `volume(center, radius)`, `fibre_volume(fibre, radius)`,
  arithmetic (`+`, `-`, scalar `*`), `rotate(rotation)`, `symmetrise()`, and quadrature to
  `SO3FunHarmonic`.
- `SO3FunHarmonic` (upgrade of `HarmonicODF`): generalized-spherical-harmonic (Wigner-D)
  coefficients up to configurable bandwidth, fast evaluation, convolution with `SO3Kernel` and
  `S2Kernel`, `to_pole_figure(h)` via the Funk/projection transform, `from_orientations` by
  kernel density quadrature.
- `SO3Kernel` family (upgrade of `KernelSpec`): de la Vallée Poussin (primary, MTEX default),
  Abel-Poisson, von Mises-Fisher, Gauss-Weierstrass, bump; each with Chebyshev coefficient
  expansion, `halfwidth` ↔ `bandwidth` conversion, and convolution semantics.
- `UnimodalODF`, `FibreODF`, `UniformODF`, `BinghamODF` (new components): weighted-mixture model
  `ODF = Σ wᵢ · componentᵢ` with the full `SO3Fun` interface. The existing discrete `ODF` becomes
  the `UnimodalODF`-mixture special case rather than a separate code path.
- `ODFEstimator` policy object: kernel choice, halfwidth selection (rule-of-thumb from grain
  counts, cross-validation later), resolution — one place where "EBSD to ODF" defaults live.

### 2.4 Fibres and named components (new `texture/fibres.py`, new `texture/components.py`)

- `Fibre` (new): crystal direction ∥ specimen direction with `orientations(n)` sampling,
  `angle_to(orientation)`, symmetrized membership tests; constructors for named fibres
  (`Fibre.alpha(phase)`, `.beta`, `.gamma`, `.eta`, `.tau`, `.theta` for bcc/fcc conventions).
- `TextureComponent` (new): named ideal orientation with phase, tolerance, and citation —
  registry with Cube, Goss, Brass, Copper, S, Rotated-Cube, P, Q, R and the standard hcp/bcc
  components; `component_volume_fractions(odf | orientations, components)` analysis function.

### 2.5 Orientation statistics (`core/orientation.py`, new `core/orientation_stats.py`)

- `OrientationSet.mean()` — symmetry-aware quaternion mean (eigenvector method) with convergence
  handling; `OrientationSet.std()` / spread; `OrientationSet.cluster()` (symmetry-aware
  hierarchical or ODF-mode clustering).
- `MisorientationDistribution` (new): uncorrelated and boundary-correlated MDF from orientation
  or boundary populations, axis and angle marginal distributions, MacKenzie random baseline.

### 2.6 EBSD model upgrades (`ebsd/models.py`, new `ebsd/io/`, new `ebsd/grains.py`)

- `CrystalMap` upgrades: hexagonal scan-grid support (offset-row indexing plus honeycomb
  neighbor semantics), per-point scalar properties as named channels (IQ, CI, BC, MAD, fit),
  rectangle/polygon subregion selection, line-profile extraction, `gridify` for scattered input.
- Vendor readers (new `ebsd/io/ang.py`, `ebsd/io/ctf.py`, then `ebsd/io/h5ebsd.py`): pure-Python
  readers producing `CrystalMap` plus an auto-generated `EBSDImportManifest`, with explicit
  vendor reference-frame normalization (the existing manifest doctrine becomes the provenance
  record for direct reads, not the entry barrier).
- `Grain` upgrade / `Grain2d` (new in `ebsd/grains.py`): boundary polyline geometry, area,
  perimeter, `equivalent_diameter`, `aspect_ratio`, fitted ellipse (centroid, semiaxes, angle),
  `shape_factor`, convexity/paris, grain mean orientation, GOS (grain orientation spread), GAM.
- `GrainBoundaryNetwork` upgrade: polyline segments with true geometric length and trace
  direction, misorientation attached per segment, `select(criterion)` filtering, special-boundary
  classification via `CSLBoundary` (Brandon-criterion Σ classification) and `TwinLaw` (named twin
  laws per phase, e.g. Σ3 60°⟨111⟩, hcp extension/contraction twins), `merge_by(law)` for
  twin-merged parent grains.
- `GNDEstimator` (new, medium-term): curvature-tensor-based geometrically-necessary-dislocation
  density from orientation gradients, per-phase dislocation-system aware.

### 2.7 Plotting and color-key upgrades (`plotting/`)

- `IPFColorKey` upgrade: exact per-Laue-class sector coloring parity with the MTEX/TSL keys for
  all 11 Laue classes, arbitrary specimen reference direction (X/Y/Z or any vector), key-legend
  figure builder, and vectorized map coloring throughput.
- `AxisAngleColorKey`, `PatalaColorKey` (new, later): misorientation/boundary coloring.
- Map plotting: `plot_phase_map`, `plot_property_map` (IQ/CI/GOS/any channel),
  `plot_grain_map(color_by=...)`, boundary overlays styled by misorientation class.
- ODF sections: `plot_phi2_sections`, `plot_sigma_sections`, generic `ODFSectionPlot` builder
  reusing one contouring engine; 3D Rodrigues/Euler-space views staged later.
- `plot_pole_figure(contour|contourf|scatter)` and IPF density plots driven by `S2FunHarmonic`.

### 2.8 Physical property layer (new `properties/`, long-term)

- `CrystalTensor` (new): rank-n tensors with crystal-frame semantics, symmetry enforcement,
  rotation, Voigt notation IO; `StiffnessTensor`, `ComplianceTensor`, `ThermalExpansionTensor`.
- `SlipSystem` (new): named systems per structure (fcc {111}⟨110⟩, bcc families, hcp basal/
  prismatic/pyramidal), Schmid matrix, `schmid_factor(orientation, stress_direction)` vectorized.
- ODF-weighted homogenization: Voigt/Reuss/Hill averages, directional property surfaces
  (Young's modulus sphere plots), later wave velocities and Taylor-factor maps.

## 3. Immediate Horizon (next 1–2 development cycles)

Priorities chosen to unblock the largest number of downstream features while staying inside the
current Phase 1 validation discipline. All items include tests, JSON contracts where stable, docs,
and parity-matrix updates per `AGENTS.md`.

### 3.1 Foundational classes (build first, in this order)

1. **Symmetry completion** (2.1): full 32-point-group `PointGroup` operators, Laue classes,
   `SpecimenSymmetry`, exact `FundamentalSector` polygons. Nearly every other item depends on it.
2. **Spherical vector semantics** (2.2, first half): `SphericalVectorSet` polar constructors and
   spherical ops, `S2Grid`. Defer `S2FunHarmonic` internals to medium term if needed, but fix the
   API now.
3. **Kernel family upgrade** (2.3, kernel part): de la Vallée Poussin kernel with Chebyshev
   coefficients and halfwidth ↔ bandwidth conversion, wired into the existing discrete `ODF` and
   `HarmonicODF` so both use one kernel vocabulary.
4. **`SO3Fun` base interface** (2.3, interface only): define the abstract algebra and retrofit
   the two existing ODF classes onto it — `texture_index`, `entropy`, `volume`, `mean`, `max`
   on both. This prevents the medium-term refactor from breaking users.
5. **Fibres and named components** (2.4): `Fibre`, `FibreODF` (kernel-smeared fibre component),
   `TextureComponent` registry, `component_volume_fractions`. High teaching and research value,
   low risk, immediately demonstrable.
6. **Orientation statistics** (2.5, first half): symmetry-aware `mean()` and spread on
   `OrientationSet`; required by grains (mean orientation, GOS) and by texture components.

### 3.2 EBSD features (immediate)

- **`.ang` and `.ctf` direct readers** producing `CrystalMap` + auto-manifest with vendor-frame
  normalization. This converts PyTex from "bring a manifest" to "open your data".
- **Per-point property channels** on `CrystalMap` (IQ/CI/BC/MAD as named channels) with
  `plot_property_map`.
- **IPF map hardening**: IPF-X/Y/Z maps from one call with per-Laue-class key parity for cubic
  and hexagonal first (`plot_ipf_map(direction="X"|"Y"|"Z"|vector)`), key-legend generation,
  non-indexed/masked-point rendering.
- **Phase map plotting** and multiphase legend support.
- **Grain scalar metrics** on the existing raster grains: grain mean orientation, GOS, GAM,
  equivalent diameter from pixel counts — the cheap subset of 2.6 that does not need polyline
  geometry.
- **Cleanup completion**: wild-spike removal and CI/IQ-threshold filtering alongside the existing
  small-grain merge and majority smoothing.

### 3.3 Validation and parity work (immediate, continuing Phase 1)

- Extend the MTEX campaign corpus: color-key RGB campaigns per Laue class, fundamental-region
  boundary campaigns, kernel-coefficient campaigns, component-volume-fraction campaigns.
- Land the pending cubic and hexagonal XRDML fixtures so the PF import → ODF path claims parity.
- Fixture-backed `.ang`/`.ctf` reader tests against small open datasets with pinned hashes.

## 4. Medium-Term Horizon (roughly the following 2–4 cycles)

### 4.1 Texture core

- **`SO3FunHarmonic` at full strength** (2.3): quadrature from scattered orientations,
  convolution, Funk projection to pole figures, symmetrization; replace the internal paths of
  `HarmonicODF` with it.
- **`S2FunHarmonic`** and spherical KDE; pole figures and IPF densities become S2 functions.
- **MTEX-grade PF-to-ODF inversion**: modified-least-squares solver with non-negativity,
  zero-range method, ghost-error handling, defocusing and background correction models on
  `PoleFigureCorrectionSpec`; validation against the existing LaboTex/XRDML corpus plus MTEX
  campaign results.
- **`BinghamODF`** component and Bingham parameter estimation.
- **ODF section plotting parity**: φ₂ and σ sections with shared contouring, annotation of named
  components on sections and pole figures.
- **`MisorientationDistribution` (MDF)** with MacKenzie baseline, axis/angle marginals, and
  boundary-correlated variants.

### 4.2 EBSD program

- **Grains 2.0** (2.6): polyline boundary geometry, hex-grid segmentation, shape descriptors,
  fitted ellipses, `smooth()` for boundaries; alpha-shape support for sparse data staged last.
- **Special boundaries**: CSL Σ classification with Brandon criterion, named twin laws, twin
  merging into parent grains, boundary maps colored by classification.
- **Denoising suite**: spline and half-quadratic orientation smoothing, Kuwahara filter,
  non-indexed in-fill, all as explicit `CrystalMapFilter` policy objects with provenance records.
- **GND density estimation** from orientation gradients (curvature tensor route).
- **h5ebsd-family readers** (EDAX, Bruker, Oxford H5, kikuchipy h5ebsd) behind one HDF5 reader
  core; `.osc`/`.cpr+.crc` best-effort readers documented as vendor-limited.
- **Per-grain texture workflows**: grain-resolved ODFs, intragranular misorientation profiling,
  line profiles across maps.
- **Parent-grain reconstruction stabilization**: promote `experimental/phase_transformation.py`
  plus `core/parent_reconstruction.py` toward a variant-graph reconstruction (MTEX
  `parentGrainReconstructor` equivalent) for martensite→austenite and α→β titanium, with
  literature-fixture validation.

### 4.3 Properties (start)

- `SlipSystem` and vectorized Schmid-factor maps over `CrystalMap` (huge applied demand,
  moderate cost once symmetry completion exists).
- `CrystalTensor` core with `StiffnessTensor` and rotation/averaging basics.

### 4.4 Interop and performance

- ORIX/kikuchipy round-trip depth: orientation + phase + property channel round trips with
  equivalence tests, dictionary-indexing results ingestion.
- Vectorization pass: quaternion-native symmetry reduction throughput for full maps (millions of
  points), benchmark manifests extended with map-scale cases.

## 5. Long-Term Horizon (world-class differentiation)

- **Full property/plasticity layer**: Voigt/Reuss/Hill ODF homogenization, directional stiffness
  surfaces, seismic/elastic wave velocities, Taylor and Sachs polycrystal models, texture
  evolution simulation, deformation-twin transfer analysis.
- **3D EBSD**: `CrystalMap3d`, 3D grains and boundary meshes, serial-section registration,
  triple-line and quadruple-junction topology statistics.
- **Kikuchi/pattern-space geometry**: gnomonic projection, `KikuchiBand`/`KikuchiPattern`
  geometry, pattern-center-aware overlays, spherical indexing bridges — connecting EBSD and TEM
  pedagogy (already scoped in `references/feature_opportunities.md`).
- **Advanced statistics**: Bingham confidence regions, bootstrap uncertainty on ODFs and volume
  fractions, orientation-clustering at map scale.
- **Diffraction depth**: full defocusing-correction library for lab XRD texture, neutron/
  synchrotron spectra-to-texture workflows, Rietveld-adjacent intensity corrections, and
  spot/band indexing beyond the current geometric scaffolding.
- **Habit-plane and interface crystallography**: habit-plane determination from grain pairs,
  interface-plane distributions (five-parameter boundary character where data supports it).
- **Performance and scale**: out-of-core crystal maps, optional GPU/numba acceleration behind
  pure-Python fallbacks, streaming import for multi-GB vendor files.
- **Interactive layer**: map exploration widgets (hover orientation, click-to-PF), notebook-first
  interactive IPF/PF linking, keeping the Matplotlib runtime doctrine intact.
- **Visual parity closure**: Patala and axis-angle misorientation color keys, MTEX-comparable
  default figure styling, publication-preset themes.

## 6. EBSD Analysis Feature Ledger

Consolidated EBSD-specific list across horizons (I = immediate, M = medium, L = long):

| Feature | Horizon | Depends on |
| --- | --- | --- |
| `.ang` / `.ctf` direct readers with frame normalization | I | manifest doctrine (exists) |
| IPF-X/Y/Z maps with per-Laue-class TSL-parity color keys | I | symmetry completion |
| Phase maps, property-channel maps (IQ, CI, BC, MAD) | I | `CrystalMap` channels |
| IPF key legend figures | I | `FundamentalSector` polygons |
| Grain mean orientation, GOS, GAM, equivalent diameter | I | orientation mean |
| Wild-spike removal, CI/IQ filtering | I | — |
| KAM with all orders/thresholds on hex grids | implemented | labelled analytic fixture; shared real MTEX fixture remains |
| Hex-grid segmentation and neighbor semantics | implemented | `CrystalMap.grid_kind` + `row_lengths` |
| Polyline grain boundaries, shape descriptors, fitted ellipses | M | Grains 2.0 |
| CSL/twin boundary classification and twin merging | M | `Misorientation` + laws |
| MDF and boundary character distributions | M | `MisorientationDistribution` |
| Orientation denoising (spline, half-quadratic, Kuwahara) | M | filter framework |
| GND density maps | M | gradients + slip systems |
| Grain-resolved ODF and intragranular analysis | M | `SO3Fun` algebra |
| Schmid-factor maps | M | `SlipSystem` |
| h5ebsd-family readers | M | HDF5 reader core |
| Parent-grain reconstruction (stable) | M | variant graph method |
| 3D EBSD maps and grains | L | `CrystalMap3d` |
| Kikuchi band geometry and pattern-space overlays | L | gnomonic layer |
| Interactive map exploration | L | plotting runtime |
| Five-parameter boundary character | L | 3D or stereology |

## 7. Sequencing Rationale

1. Symmetry completion is first because IPF keys, Miller families, fundamental regions,
   structure-factor extinctions, MDFs, and CSL classification all consume it; building any of
   those first would create the private symmetry semantics that `AGENTS.md` forbids.
2. Vendor readers are co-first because adoption and validation both improve when real datasets
   flow in directly; every later EBSD feature is easier to test with real `.ang`/`.ctf` corpora.
3. The `SO3Fun`/`S2Fun` algebra lands as an interface immediately and an implementation in the
   medium term, so the public API stabilizes one cycle before the heavy math replaces the
   internals.
4. Grains 2.0 precedes MDF/CSL/GND because those consume boundary geometry and grain statistics.
5. The property layer is deliberately after the texture core: Schmid/tensor analysis is only as
   trustworthy as the symmetry and orientation spine underneath it.
6. Each horizon keeps the Phase 1 rule: no parity claim without a campaign fixture or pinned
   external baseline, and docs/figures/tests land with the feature.

## References

### Normative

- [Implementation Roadmap](implementation_roadmap.md)
- [MTEX Parity Matrix](../testing/mtex_parity_matrix.md)
- [Reference Canon](../standards/reference_canon.md)
- [Canonical Data Model](../architecture/canonical_data_model.md)

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
- [Feature Opportunities From The Reference Corpus](../../references/feature_opportunities.md)
- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*
- Rowenhorst et al., *Consistent representations of and conversions between 3D rotations*, MSMSE 23 (2015)
- Randle, V. & Engler, O., *Introduction to Texture Analysis: Macrotexture, Microtexture and Orientation Mapping*
