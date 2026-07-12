# MTEX Parity Matrix

This document is the authoritative ledger for PyTex parity against public MTEX tests and documented example categories.

Reference baseline:

- MTEX `6.0.0`
- public repository `tests/` tree
- public documentation examples where they expose behavior not captured by tests

The campaign-based parity runner targets MTEX `6.0.0`, released in November 2024, so validation
does not depend on the newest MTEX release line. Regeneration scripts live under
`scripts/mtex_generators/`, and shared campaign inputs live under `fixtures/mtex_parity/campaigns/`.

## Status Keys

- `implemented`: PyTex has equivalent or stronger automated coverage
- `foundational`: PyTex has a correct foundational implementation, but not yet full parity with the MTEX behavior category
- `planned`: parity target is accepted but not implemented yet
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | MTEX reference examples/tests | PyTex status | Notes |
| --- | --- | --- | --- |
| Euler/quaternion conversions | `check_eulerquat.m`; `fixtures/mtex_parity/campaigns/orientation_core_cases.json` | implemented | Baseline conversion and normalization tests are present, and the public rotation surface now also covers vectorized axis-angle, matrix, and Rodrigues / Rodrigues--Frank conversions. The campaign runner adds shared JSON inputs and MTEX/PyTex result JSON for Euler, matrix, quaternion, axis-angle, Miller construction, composition, inverse, and misorientation cases. Misorientation now composes in the crystal frame (`inv(o1) * o2`, the MTEX convention) so crystal-symmetry reduction through fixed left/right operator products is exact; symmetry-equivalent orientations measure zero disorientation, with regression coverage. Symmetry-aware quaternion mean (`OrientationSet.mean_orientation`) and spread statistics are implemented. |
| Orientation construction and grids | `orientation.byEuler`, `orientation.byAxisAngle`, `orientation.byMiller`, `regularSO3Grid`, `equispacedSO3Grid` | foundational | Scalar and batch orientation constructors now cover Euler, axis-angle, matrix, quaternion, Rodrigues / Rodrigues-Frank, Miller plane-direction correspondences, regular Bunge SO3 grids, and deterministic quaternion SO3 grids with shared frame and symmetry checks. Campaign coverage exists for scalar construction and operations; external grid parity breadth remains ahead. |
| Miller planes, directions, and family expansion | MTEX `Miller` object workflows and indexing examples; `fixtures/mtex_parity/campaigns/miller_geometry_cases.json` | foundational | First-class `MillerPlane` / `MillerDirection` scalar and batch objects now cover vectorized reciprocal/direct vectors, d-spacings, hexagonal index conversion, angle relations, projection, symmetry-family expansion, JSON contracts, orix interop, and Miller-backed orientation construction from three- or four-index inputs. The campaign runner adds compact cubic and hexagonal geometry/family checks. Broader MTEX parity for plotting idioms and every crystal-class helper remains ahead. |
| Fundamental regions | `check_fundamentalRegion.m` | foundational | Exact orbit-based minimum-angle reduction is implemented and regression-tested across the supported proper point groups. `FundamentalSector` now carries exact Laue-reduced sector polygons per proper rotation group (ordered vertices plus inward edge-normal half-spaces), with vectorized containment, spherical boundary traces, and sector centers; Monte-Carlo area checks confirm each sector covers 1/(2·order) of the sphere and vector reduction is verified to land inside the exact sector for every class, fixing previously undersized cyclic-group (2, 3, 4, 6, 23) sectors and the trigonal 32 wedge placement. Closed-form Rodrigues-space orientation polytopes are still ahead. |
| Symmetry operators and SO(3) basics | `SO3FunTests`, `check_WignerD.m`, symmetry-related checks | foundational | All 32 crystallographic point groups are now first-class through `PointGroup` (full operator sets including mirrors, inversion, and rotoinversions; Hermann-Mauguin and Schoenflies naming; Laue-class and proper-subgroup bridges), with order, closure, Laue-assignment, mirror-count, and family-expansion regression coverage. `SymmetrySpec` accepts every group symbol, exposes Laue semantics, and gains a specimen-symmetry constructor (triclinic, monoclinic, orthorhombic/orthotropic). Exact orbit reduction, symmetry actions, and a first Wigner-basis harmonic reconstruction surface were already implemented; broader external parity is still ahead. |
| Spherical projections and stereonets | public spherical-projection examples and plotting workflows | foundational | Wulff-net plotting, stereographic direction/plane projection, great-circle traces, and rotational symmetry-axis symbols are implemented with deterministic regression coverage, but full MTEX visual-parity claims are still ahead. The core model now also provides `vector3d`-style semantic sphere primitives: `SphericalVectorSet` (polar constructors, antipodal-aware angles, spherical/orientation-tensor means, folding) and `S2Grid` (equispaced and regular sphere grids with normalized band-area quadrature weights), with the canonical polar-angle convention centralized in `core/sphere.py` and reused by the stereonet surface. |
| EBSD container basics | `check_ebsd.m` | implemented | `CrystalMap` plus fixture-backed regular-grid segmentation coverage are implemented. Per-point scalar property channels (IQ/CI/BC/MAD/fit) are now first-class on `CrystalMap` (`properties`, `get_property`, `property_map`, `with_properties`), carried through `select_phase`/`select_points` and populated directly by the `.ang`/`.ctf` readers. |
| KAM-related behavior | `testKAM2.m` | implemented | Fixture-backed regular-grid KAM support covers order, threshold, and max-style aggregation. |
| GROD and grain-local orientation metrics | public EBSD workflow examples | implemented | Fixture-backed GROD relative to a representative grain orientation is implemented, alongside per-grain scalar metrics on `GrainSegmentation`: symmetry-aware grain mean orientation, grain orientation spread (GOS), grain average misorientation (GAM), grain size, and equivalent circular diameter, each with per-point map broadcasts (`gos_map_deg`, `gam_map_deg`). Moment-based grain shape descriptors are also implemented: `FittedEllipse` second-moment equivalent ellipses (centroid, semi-axes, orientation angle), aspect ratios, and axis-aligned bounding boxes from member-pixel positions. True staircase grain perimeter (`grain_perimeters`, rectangular-step aware, including map-edge faces), grain area, and the compactness shape factor `P/(2 sqrt(pi A))` (`grain_shape_factors`) are now implemented on regular grids. Ordered polyline boundary tracing and convexity/paris descriptors remain ahead in Grains 2.0. |
| Grain-boundary and cleanup workflows | public EBSD segmentation examples | implemented | Fixture-backed boundary extraction and adjacency-based small-grain merging are implemented for regular grids, together with majority smoothing, wild-spike removal (`remove_wild_spikes`, neighborhood-mean infill of isolated spikes), and property-threshold cleanup (`property_threshold_mask`, `filter_by_property`, `select_points`) for CI/IQ-style quality gating. Special-boundary classification is now available for cubic maps: `pytex.ebsd.csl` provides the Sigma1-Sigma29 CSL registry, the Brandon criterion, symmetry-reduced misorientation-deviation matching (`classify_misorientations`), and named cubic twin laws (Sigma3 60 deg <111>), wired into `GrainBoundaryNetwork` via `classify_csl` / `csl_fraction` / `select_csl`. Twin/CSL merging into parent grains is now implemented (`merge_by_csl`, `twin_merge`), unioning grains joined by a chosen CSL boundary type into MTEX-style parent grains. Non-cubic CSL tables remain ahead. |
| IPF color coding | `checkIpfColorCoding.m`; `fixtures/mtex_parity/campaigns/ipf_color_cases.json` | foundational | `IPFColorKey` exists, is symmetry-aware, and now consumes the exact fundamental-sector geometry: key legend meshes and sector boundary outlines are available for all 11 Laue classes via `legend_mesh` / `boundary_points_2d` / `plot_ipf_key`, with low-symmetry classes falling back to documented octant color anchors. The campaign runner emits comparable RGB and reduced-direction JSON for cubic and hexagonal cases, but full MTEX color-key visual parity is not yet claimed. Map plotting now covers IPF maps for any sample direction plus an IPF-X/Y/Z triptych (`plot_ipf_xyz_maps`), per-point property maps (`plot_property_map`), and phase maps with named legends (`plot_phase_map`). |
| ODF and PF reconstruction | `check_FourierODF.m`, PF reconstruction examples; `fixtures/mtex_parity/campaigns/odf_discrete_cases.json` | foundational | Discrete/kernel ODF evaluation, PF/IPF synthesis, explicit dictionary-based PF inversion, and a band-limited harmonic ODF inversion surface are implemented. The de la Vallee Poussin kernel now exists as a first-class object (`DeLaValleePoussinKernel`) with the halfwidth/kappa relation, normalized evaluation, quadrature-based Chebyshev character coefficients, and bandwidth estimation. It is now wired into the discrete ODF path: `KernelSpec.as_so3_kernel()` returns the normalized dVP object, `KernelSpec.bandwidth()` exposes the halfwidth-to-bandwidth duality, and both `KernelSpec.evaluate` and `ODF.evaluate` accept `normalized=True` for m.r.d.-scaled densities. The full harmonic ODF-estimation replacement remains staged with the SO3Fun refactor. Named texture components (`TextureComponent`, standard fcc/bcc rolling registries) with symmetry-aware component volume fractions, and `Fibre` objects (named bcc alpha/gamma/eta/theta fibres) with fibre sampling, symmetry-aware fibre distances, and fibre volume fractions are implemented. The campaign runner currently covers sampled discrete ODF density from weighted Euler orientations. XRDML PF and reconstruction campaigns are defined as pending until cubic and hexagonal fixtures are supplied. |
| Interfaces/import-export | `checkInterfaces.m`; `fixtures/mtex_parity/campaigns/xrdml_pole_figure_cases.json` | foundational | Stable EBSD import manifests, manifest IO, object-backed vendor bridge adapters, CIF-backed phase creation, XRDML pole-figure import, LaboTex PPF/EPF pole-figure import, and lightweight orix / KikuchiPy bridge surfaces now exist. Direct pure-Python vendor EBSD scan readers are now available for EDAX/TSL `.ang` (square grids, single and multiphase, TSL symmetry-code mapping) and Oxford/HKL Channel 5 `.ctf` (multiphase, HKL Laue-code mapping, non-indexed point handling), producing `CrystalMap` objects plus auto-generated import manifests and per-point quality properties (IQ/CI/fit, bands/MAD/BC/BS). Hex-grid `.ang` scans and h5ebsd-family readers remain ahead, as do real-file external fixtures; current coverage is deterministic in-repo scan content. Shared XRDML parity campaign inputs are present but pending real cubic and hexagonal XRDML fixtures. |
| Misorientation distribution (MDF) | public MDF / MacKenzie-baseline workflows | foundational | `MisorientationDistribution` builds symmetry-reduced disorientation-angle distributions from orientation populations, both uncorrelated (all unique pairs) and correlated (explicit boundary/neighbour pairs), with angle statistics and histograms. `random_baseline` provides the random (MacKenzie-type) distribution via seeded Haar-uniform sampling with point-group reduction; validated against the known cubic bounds (max ~62.8 deg, mean and peak near 42-45 deg). Analytic closed-form MacKenzie densities and boundary-correlated MDF plotting remain ahead. |
| Physical properties (slip / Schmid) | public slip-system / Schmid-factor workflows | foundational | A first `pytex.properties` layer provides `SlipSystem`, symmetry-expanded `SlipSystemFamily` (12-system fcc `{111}<110>` and bcc `{110}<111>` families), vectorized absolute Schmid factors over orientation populations, and `CrystalMap.schmid_factor_map` for per-point maximum Schmid factor maps. Validated against the analytic `sqrt(6)/6` maximum Schmid factor for `[001]` uniaxial loading. A rank-4 elastic-tensor layer is now implemented too: `StiffnessTensor` / `ComplianceTensor` with 6x6 Voigt IO, cubic/hexagonal/isotropic constructors, 4th-rank rotation, compliance inversion, and directional Young's modulus, validated against the analytic cubic `1/E` formula and isotropy/rotation invariance. Orientation-weighted Voigt/Reuss/Hill polycrystal homogenization (`homogenize_elastic`) is implemented, validated by single-orientation reproduction of the crystal, near-isotropy of a random aggregate, and the Voigt > Hill > Reuss bound ordering. Taylor factors and full directional wave-velocity surfaces remain ahead. |

## PyTex-Only Extensions

The following categories must exceed MTEX coverage:

- vendor reference-frame normalization
- provenance retention and manifest integrity
- adapter interoperability with ORIX, KikuchiPy, PyEBSDIndex, pymatgen, and diffsims
- LaTeX/SVG documentation asset integrity
- workflow-level reproducibility from CLI entry points

## Current Hardening Priorities

The highest-value remaining MTEX-facing hardening work is concentrated in the foundational rows
that underpin multiple downstream subsystems:

- orientation construction and deterministic SO(3) grid coverage
- Miller plane and direction workflows beyond the current starter families
- broader class-by-class fundamental-region coverage
- richer symmetry and harmonic parity breadth
- broader ODF/PF reconstruction breadth
- stronger interface and import or export parity at the adapter boundary

## Campaign Regeneration Workflow

The detailed step-by-step workflow is documented in:

- `docs/site/validation/mtex_regeneration.md`
- `scripts/mtex_generators/README.md`

Minimal example on a MATLAB system with MTEX started:

```matlab
addpath("scripts/mtex_generators")
run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
                         "fixtures/mtex_parity/results/mtex")
```

Repeat for each campaign file, bring `fixtures/mtex_parity/results/mtex/` back to the PyTex
machine, then run:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
