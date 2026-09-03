<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->
# Worked Examples

This section is PyTex's *documentation-as-test* surface. Each worked example bundles a scientific scenario, a runnable snippet, the value that snippet computes from the live code, and an independently known reference value. The examples are the single source of truth for both this gallery and the regression test `tests/unit/test_worked_examples.py`.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Reference-value ledger

The complete set of computed-versus-expected values at a glance:

| Example | Computed (live) | Expected (reference) | Unit | Status |
| --- | --- | --- | --- | --- |
| `frame-axis-correspondence-angle` | 90.0000 | 90.0000 | deg | ✅ |
| `frame-axis-correspondence-components` | [0.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 0.0000] | &mdash; | ✅ |
| `frame-graph-multi-hop-angle` | 30.0000 | 30.0000 | deg | ✅ |
| `frame-round-trip-residual` | 0.0000 | 0.0000 | &mdash; | ✅ |
| `frame-right-handed-determinant` | 1.0000 | 1.0000 | &mdash; | ✅ |
| `reciprocal-frame-star-count` | [3, 0] | [3, 0] | &mdash; | ✅ |
| `miller-indices-carry-no-star` | 0.0000 | 0.0000 | &mdash; | ✅ |
| `cubic-angle-100-110` | 45.0000 | 45.0000 | deg | ✅ |
| `cubic-angle-100-111` | 54.7356 | 54.7356 | deg | ✅ |
| `cubic-angle-dir-110-111` | 35.2644 | 35.2644 | deg | ✅ |
| `cubic-dspacing-111` | 2.30940 | 2.30940 | angstrom | ✅ |
| `cubic-plane-multiplicity` | [3, 6, 4, 24] | [3, 6, 4, 24] | &mdash; | ✅ |
| `hex-angle-basal-prism` | 90.0000 | 90.0000 | deg | ✅ |
| `hex-angle-prism-prism` | 60.0000 | 60.0000 | deg | ✅ |
| `cubic-nearest-low-index-round-trip` | 0.0000 | 0.0000 | deg | ✅ |
| `hex-nearest-low-index-plane-vs-direction` | [1, 0, 1] | [1, 0, 1] | &mdash; | ✅ |
| `orientation-euler-matrix-roundtrip` | 0.0000 | 0.0000 | deg | ✅ |
| `orientation-sigma3-disorientation` | 60.0000 | 60.0000 | deg | ✅ |
| `core-orientation-equal-volume-charts-agree-on-the-so3-volume` | [1.000000000000, 1.000000000000, 1.000000000000] | [1.000000000000, 1.000000000000, 1.000000000000] | &mdash; | ✅ |
| `core-orientation-ideal-indices-invert-the-plane-direction-construction` | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | &mdash; | ✅ |
| `diffraction-ni-111-two-theta` | 44.496 | 44.496 | deg | ✅ |
| `diffraction-ni-111-kikuchi-band-width` | 2.4188 | 2.4187 | deg | ✅ |
| `diffraction-gnomonic-zone-axis-radius` | 1.000000000000 | 1.000000000000 | &mdash; | ✅ |
| `diffraction-ebsd-specimen-normal-radius` | 0.363970234266 | 0.363970234266 | &mdash; | ✅ |
| `diffraction-march-dollase-family-factor` | 2.309132723130 | 2.309132723130 | &mdash; | ✅ |
| `diffraction-march-dollase-normalization` | 1.000000000 | 1.000000000 | &mdash; | ✅ |
| `diffraction-odf-weighted-random-texture` | 0.9999 | 1.0000 | &mdash; | ✅ |
| `diffraction-powder-profile-affine-comparison` | [5.000000000000, 5.000000000000, 0.000000000000, 0.000000000000] | [5.000000000000, 5.000000000000, 0.000000000000, 0.000000000000] | &mdash; | ✅ |
| `diffraction-kikuchi-map-zone-axis-tilt-angles` | [45.000000, 54.735610, 35.264390] | [45.000000, 54.735610, 35.264390] | deg | ✅ |
| `ebsd-hex-grid-six-neighbor-kam` | [13.000000, 3.000000, 6.000000, 3.000000, 1.200000, 1.200000, 0.000000, 0.000000, 0.000000] | [13.000000, 3.000000, 6.000000, 3.000000, 1.200000, 1.200000, 0.000000, 0.000000, 0.000000] | degree (KAM entries) | ✅ |
| `ebsd-planted-lattice-curvature` | 13962.634016 | 13962.634016 | rad/m | ✅ |
| `ebsd-gnd-density-from-curvature` | 5.462689e+13 | 5.462689e+13 | 1/m^2 | ✅ |
| `composite-electron-wavelength-200kv` | 0.02508 | 0.02508 | angstrom | ✅ |
| `composite-ks-exact-child-zone` | 0.0000 | 0.0000 | deg | ✅ |
| `composite-burgers-exact-basal-zone` | 0.0000 | 0.0000 | deg | ✅ |
| `composite-burgers-110-0002-coincidence` | 0.15450 | 0.15450 | mm | ✅ |
| `composite-burgers-reflection-table-identities` | [0.0000, 0.0000, 0.0000, 2.3380] | [0.0000, 0.0000, 0.0000, 2.3380] | counts, angstrom, counts, angstrom | ✅ |
| `composite-child-anchored-geometry-consistency` | [0.0000, 0.0000] | [0.0000, 0.0000] | mm | ✅ |
| `solving-simulate-then-solve-closure` | [1.0000, 0.0000, 110.0000, 0.0000] | [1.0000, 0.0000, 110.0000, 0.0000] | fraction, 1/angstrom, family code, count | ✅ |
| `kinematic-silicon-double-diffraction-002` | 0.500000 | 0.500000 | dimensionless | ✅ |
| `saed-finite-thickness-rectangular-slab` | [1.000000000000, 0.405284734569, 0.000000000000] | [1.000000000000, 0.405284734569, 0.000000000000] | dimensionless | ✅ |
| `diffraction-cbed-aluminium-extinction-distances-at-100kv` | [555.2, 663.9, 1062.5] | [556.0, 673.0, 1057.0] | angstrom | ✅ |
| `diffraction-cbed-two-beam-thickness-inverts-the-fringe-relation` | [2000.000000, 500.000000] | [2000.000000, 500.000000] | angstrom | ✅ |
| `diffraction-dynamical-two-beam-limit-of-the-many-beam-solver` | < 1e-12 | 0.00e+00 | &mdash; | ✅ |
| `diffraction-dynamical-intensity-is-conserved-without-absorption` | < 1e-12 | 0.00e+00 | &mdash; | ✅ |
| `diffraction-holz-strain-and-wavelength-are-exactly-degenerate` | < 1e-15 | 0.00e+00 | radian | ✅ |
| `diffraction-groups-construction-yields-buxtons-thirty-one` | 31 | 31 | &mdash; | ✅ |
| `diffraction-groups-friedel-observation-splits-the-point-groups` | [21, 11] | [21, 11] | &mdash; | ✅ |
| `diffraction-groups-zincblende-down-001-loses-the-two-fold` | [1, 1, 1, 0, 1, 1, 1] | [1, 1, 1, 0, 1, 1, 1] | &mdash; | ✅ |
| `texture-gaussian-kernel-normalization-and-halfwidth` | [1.000000, 0.500000] | [1.000000, 0.500000] | &mdash; | ✅ |
| `texture-uniform-odf-pole-density-is-one-mrd` | [1.000013, 1.000248, 0.999713] | [1.000000, 1.000000, 1.000000] | m.r.d. | ✅ |
| `texture-odf-named-component-fit-recovers-exact-mixture` | [0.700000000, 0.300000000, 0.000000000] | [0.700000000, 0.300000000, 0.000000000] | volume fraction | ✅ |
| `ipf-cubic-sector-corners-are-primaries` | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | &mdash; | ✅ |
| `ipf-cubic-closed-form-colour-113` | [1.000000, 0.000000, 0.750000] | [1.000000, 0.000000, 0.750000] | &mdash; | ✅ |
| `ipf-symmetric-equivalents-share-one-colour` | < 1e-12 | 0.00e+00 | &mdash; | ✅ |
| `mdf-triclinic-mean-disorientation-angle` | 126.3581 | 126.4756 | deg | ✅ |
| `mdf-cubic-maximum-disorientation-angle` | 62.7994 | 62.7994 | deg | ✅ |
| `mdf-cubic-random-low-angle-fraction` | 0.0223 | 0.0223 | &mdash; | ✅ |
| `elastic-cubic-youngs-modulus-110-equals-112` | [66.6888, 130.3376, 130.3376, 191.1497] | [66.6888, 130.3376, 130.3376, 191.1497] | GPa | ✅ |
| `elastic-cubic-voigt-reuss-bulk-moduli-coincide` | < 1e-11 | 0.000e+00 | GPa | ✅ |
| `elastic-random-aggregate-matches-voigt-reuss-closed-form` | [54.6314, 40.0264] | [54.6400, 40.0339] | GPa | ✅ |
| `plasticity-fcc-cube-schmid-factor` | 0.408248 | 0.408248 | &mdash; | ✅ |
| `plasticity-fcc-cube-taylor-factor` | 2.449490 | 2.449490 | &mdash; | ✅ |
| `plasticity-random-fcc-taylor-factor` | 3.0546 | 3.0600 | &mdash; | ✅ |
| `pole-figure-raster-unweighted-mean-is-biased` | 0.500000 | 0.500000 | &mdash; | ✅ |
| `pole-figure-raster-bias-survives-refinement` | 1.000000000 | 1.000000000 | &mdash; | ✅ |
| `pole-figure-raster-weighted-mean-converges` | [0.31960, 0.32627] | [0.31960, 0.32627] | &mdash; | ✅ |
| `pole-figure-random-standard-defocus-calibration` | [1.000000, 0.800000, 0.500000, 20.000000] | [1.000000, 0.800000, 0.500000, 20.000000] | &mdash; | ✅ |
| `directional-orientation-tensor-limiting-eigenvalues` | [0.0000, 0.4996, 0.5004, 0.0000, 0.0000, 1.0000] | [0.0000, 0.5000, 0.5000, 0.0000, 0.0000, 1.0000] | &mdash; | ✅ |
| `directional-orientation-tensor-unit-trace` | 1.000000000000 | 1.000000000000 | &mdash; | ✅ |
| `directional-mean-axis-of-randomly-signed-axes` | 0.999982 | 1.000000 | &mdash; | ✅ |
| `ghost-pole-figure-is-centrosymmetric` | 1.000000 | 1.000000 | &mdash; | ✅ |
| `ghost-odd-degrees-are-half-the-harmonic-basis` | 0.46740 | 0.46739 | &mdash; | ✅ |
| `texture-pole-figure-mrd-unit-mean-density` | [1.000000000000, -0.000000000000] | [1.000000000000, 0.000000000000] | m.r.d. | ✅ |
| `texture-pole-figure-resampling-and-addition-identities` | [0.000000000000, 2.000000000000] | [0.000000000000, 2.000000000000] | m.r.d. | ✅ |
| `kearns-random-texture-is-one-third` | 0.001800 | 0.000000 | &mdash; | ✅ |
| `kearns-triad-sum-is-exactly-one` | 1.000000000000 | 1.000000000000 | &mdash; | ✅ |
| `kearns-ideal-basal-girdle` | [0.500000, 0.500000, 0.000000] | [0.500000, 0.500000, 0.000000] | &mdash; | ✅ |
| `kearns-1965-table-3-longitudinal-section` | 0.4879 | 0.4880 | &mdash; | ✅ |
| `tem-tilt-001-to-011-travel` | 45.000 | 45.000 | deg | ✅ |
| `tem-tilt-001-to-111-travel` | 54.7356 | 54.7356 | deg | ✅ |
| `tem-holder-accessible-solid-angle` | 1.04720 | 1.04720 | sr | ✅ |
| `tem-diffraction-rotation-residual` | 5.000000 | 5.000000 | deg | ✅ |
| `tem-observation-stabilizer-cubic-001` | 8 | 8 | operators | ✅ |
| `tem-symmetry-orbit-multiplicity` | 48 | 48 | directions | ✅ |
| `tem-indexed-orientation-identity` | 0.00e+00 | 0.00e+00 | dimensionless | ✅ |
| `tem-self-calibrated-diffraction-rotation` | 37.000000 | 37.000000 | deg | ✅ |
| `saed-practice-camera-constant-identity` | 4.95454 | 4.95454 | mm | ✅ |
| `saed-practice-roll-about-the-beam` | 30.000000 | 30.000000 | deg | ✅ |
| `saed-practice-hcp-prism-axial-ratio` | 1.08762 | 1.08762 | &mdash; | ✅ |
| `saed-practice-atlas-basal-to-prism` | 90.000000 | 90.000000 | deg | ✅ |
| `saed-lattice-fit-recovers-the-beam-centre` | 0.000000000 | 0.000000000 | px | ✅ |
| `saed-scoring-calibration-bias` | 0.0476190476 | 0.0476190476 | &mdash; | ✅ |
| `saed-practice-double-diffraction-forbidden-200` | 3.69422 | 3.69422 | mm | ✅ |
| `or-ks-plane-correspondence-identity` | [0.0000, 1.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 1.0000, 0.0000] | indices, deg | ✅ |
| `or-bain-direction-correspondence-identity` | [1.0000, 0.0000, 0.0000, 0.0000] | [1.0000, 0.0000, 0.0000, 0.0000] | indices, deg | ✅ |
| `or-ks-misorientation-representation` | [42.8478, 0.9679, 0.1776, 0.1776] | [42.8500, 0.9679, 0.1776, 0.1776] | deg, axis components | ✅ |
| `or-fit-recovers-gt-from-ks-nominal` | [2.4037, 0.0000] | [2.4037, 0.0000] | deg | ✅ |
| `or-ks-same-parent-boundary-fingerprint` | [0.0000, 0.0000] | [0.0000, 0.0000] | deg | ✅ |
| `or-ks-identified-from-measured-orientations` | [0.0000, 5.2644] | [0.0000, 5.2600] | deg | ✅ |
| `or-ks-parallelism-statement-from-rotation` | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | indices, deg | ✅ |
| `or-ks-variant-correspondence-packets` | [24, 4, 6, 6, 6] | [24, 4, 6, 6, 6] | counts | ✅ |
| `or-ks-variant-parallelisms-are-per-variant` | [0.0000, 4.0000] | [0.0000, 4.0000] | deg, count | ✅ |
| `or-dossier-agrees-with-its-sources` | [0.0000, 0.0000, 24.0000, 4.0000, 10.0000, 60.0000] | [0.0000, 0.0000, 24.0000, 4.0000, 10.0000, 60.0000] | angstrom^3, angstrom^3, counts, deg | ✅ |
| `or-rationalization-costs-the-ks-gt-separation` | [0.0000, 1.0000, 1.0000, 2.4037] | [0.0000, 1.0000, 1.0000, 2.4037] | deg, booleans, deg | ✅ |
| `or-fcc-twin-is-sigma-3` | [180.0000, 60.0000, 4.0000] | [180.0000, 60.0000, 4.0000] | deg, deg, count | ✅ |
| `or-cube-on-cube-is-the-identity` | [0.0000, 1.0000, 1.0000, 1.0000, 2.0000] | [0.0000, 1.0000, 1.0000, 1.0000, 2.0000] | deg, count, indices | ✅ |
| `or-custom-statement-reproduces-ks` | [0.0000, 24.0000, 42.8478] | [0.0000, 24.0000, 42.8500] | deg, count, deg | ✅ |
| `viz-transform-crystal-to-sample-consistency` | 0.0000 | 0.0000 | &mdash; | ✅ |
| `viz-or-parallel-direction-alignment` | 1.0000 | 1.0000 | &mdash; | ✅ |
| `viz-scene-bond-length-halite-identity` | 2.0000 | 2.0000 | angstrom | ✅ |
| `viz-or-stereogram-parallelism-coincides` | [0.0000, 0.0000] | [0.0000, 0.0000] | projection-plane units | ✅ |
| `workbench-ks-packet-size` | [6, 6, 6, 6] | [6, 6, 6, 6] | variants per packet | ✅ |
| `workbench-ks-intervariant-spectrum` | [10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00] | [10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00] | deg | ✅ |
| `workbench-mrd-mean-is-one` | [1.000, 1.000, 1.000] | [1.000, 1.000, 1.000] | m.r.d. | ✅ |
| `workbench-goss-pole-at-nd` | 0.0 | 0.0 | deg | ✅ |
| `workbench-crystal-viewer-goss-nd` | 45.0000 | 45.0000 | deg | ✅ |
| `workbench-crystal-viewer-euler-round-trip` | < 1e-08 | 0.00e+00 | deg | ✅ |

## Example groups

- {doc}`Reference Frames And Frame Transforms <generated/reference_frames>` &mdash; Creating standard frames, declaring frame relationships in words, and letting the frame graph compose multi-step chains — with the rotation angles, components, and invariants checked against exact analytic values. The last two examples pin the IUCr notation convention: the reciprocal star marks the basis, never the indices.
- {doc}`Crystal geometry: angles, spacings, and multiplicities <generated/crystal_geometry>` &mdash; Interplanar and interdirection angles, interplanar spacings, and symmetry multiplicities for cubic and hexagonal phases, and the naming of a direction or a plane that arrived as a Cartesian vector. Each result is checked against an analytic identity for the relevant crystal system.
- {doc}`Orientations and disorientation angles <generated/orientation>` &mdash; Round-trip consistency of orientation representations and symmetry-reduced disorientation angles, checked against exact identities and the Sigma 3 twin reference.
- {doc}`Orientation representations <generated/orientation-representations>` &mdash; The constants and identities behind the equal-volume charts of SO(3), and the inversion that names an orientation as a (hkl)[uvw] texture component.
- {doc}`Diffraction geometry <generated/diffraction>` &mdash; Powder scattering angles from PyTex interplanar spacings via Bragg's law, Kikuchi band and zone-axis geometry in the gnomonic projection, zone-axis routing on a stereographic Kikuchi map, the EBSD camera geometry, and preferred-orientation corrections to powder intensities — each checked against a standard reference value or a closed-form identity.
- {doc}`EBSD microstructure <generated/ebsd>` &mdash; Hex-grid KAM, lattice curvature, and geometrically necessary dislocation density recovered from analytically planted topology or orientation gradients.
- {doc}`Composite OR diffraction <generated/composite-diffraction>` &mdash; Numerical cornerstones of composite orientation-relationship SAED simulation: the relativistic electron wavelength against the standard 200 kV value, the exactness of the Kurdjumov-Sachs child-zone mapping, and the two defining Burgers beta->alpha signatures (exact basal zone and the {110}/(0002) near-coincidence), plus the identities the exported reflection table must satisfy, and the exact halfway position of the double-diffraction Si 002 spot, and the analytic landmarks of the finite-thickness relrod shape factor.
- {doc}`Convergent-beam diffraction <generated/convergent-beam-diffraction>` &mdash; The absolute scale of the two-beam extinction distance, checked against a published table, and the fringe analysis that measures a foil thickness without needing one.
- {doc}`Dynamical CBED and symmetry determination <generated/dynamical-cbed-and-symmetry>` &mdash; The exact limits that calibrate a many-beam calculation, the HOLZ degeneracy that makes voltage calibration mandatory, and the diffraction-group construction that determines a point group including its centre of symmetry.
- {doc}`Texture kernels <generated/texture>` &mdash; Analytic identities of the SO(3) kernel surface - normalization (A_0 = 1) and the halfwidth definition - together with the m.r.d. scale on which pole densities are reported, plus an exactly identifiable named-component mixture, all computed live.
- {doc}`Inverse-pole-figure colouring <generated/ipf-coloring>` &mdash; What an IPF colour actually is, checked against hand-derived values: the sector corners colour to exact primaries, a direction on the [001]-[111] edge colours to exactly (1, 0, 3/4) by the closed form, and every symmetric equivalent shares one colour.
- {doc}`Random disorientation baseline <generated/random-disorientation>` &mdash; The null hypothesis every MDF claim is measured against: the exact (1 - cos w)/pi density with no symmetry, the exact cubic maximum 62.7994 degrees from the Rodrigues zone vertex, and the 2.2 percent of random cubic boundaries that are low-angle by chance.
- {doc}`Elastic anisotropy and homogenization <generated/elastic-anisotropy>` &mdash; Cubic elasticity against closed forms: [110] and [112] are exactly equally stiff, the Voigt and Reuss bulk moduli of a cubic aggregate are identical so the entire bound gap is in the shear modulus, and a numerically homogenized random aggregate reproduces both shear bounds.
- {doc}`Schmid and Taylor plasticity factors <generated/schmid-and-taylor>` &mdash; Slip geometry against exact answers: eight fcc systems share a Schmid factor of 1/sqrt(6) under [001] tension, the cube orientation's full-constraint Taylor factor is exactly sqrt(6), and a random fcc texture averages Taylor's 1938 value of 3.06.
- {doc}`Pole-figure raster sampling <generated/pole-figure-sampling>` &mdash; Why a measured pole figure needs solid-angle weights: the unweighted mean of cos^2 over a tilt raster is exactly 1/2 against a true spherical mean of 1/3, a 50 percent bias that survives halving the raster step, while the weighted mean converges.
- {doc}`Directional statistics and mean axes <generated/directional-statistics>` &mdash; Averaging axes rather than vectors: the orientation tensor has unit trace, its eigenvalues take exact values at the girdle and cluster limits, and it recovers a fibre axis from randomly signed data where the vector resultant fails outright.
- {doc}`The ghost problem <generated/ghost-problem>` &mdash; What diffraction pole figures cannot determine: an asymmetric texture still gives a pole set closed under negation, and excluding the odd harmonic degrees that centrosymmetry annihilates discards nearly half the basis.
- {doc}`Pole-figure arithmetic <generated/pole-figure-arithmetic>` &mdash; Exact identities behind comparing two pole figures: the multiples-of-random scale, resampling onto a shared support, and the additivity of pole densities.
- {doc}`The Kearns parameter <generated/kearns-parameter>` &mdash; The scalar texture index the zirconium industry specifies components against, checked against the identities that calibrate it -- 1/3 for a random texture, an exact sum of 1 over any orthonormal triad, (1/2, 1/2, 0) for an ideal basal girdle -- and against the tabulated calculation in Kearns' own 1965 report.
- {doc}`TEM tilt navigation <generated/tem_tilt_navigation>` &mdash; Holder tilts that bring a target zone axis onto the electron beam: analytic interzonal travel for the standard cubic transitions, the closed-form solid angle a double-tilt holder commands, the cost of an uncalibrated diffraction rotation, and the group-order counts that decide whether a single indexed pattern leaves a real ambiguity.
- {doc}`Simulated SAED plates and the zone-axis atlas <generated/saed_practice_patterns>` &mdash; The geometry a practice diffraction pattern must reproduce if indexing it is to teach anything: the camera-constant identity that places every reflection, the hcp prism-zone aspect ratio that measures c/a without any calibration at all, and the basal-to-prism angle the zone-axis atlas has to report as exactly 90 degrees, the beam centre a lattice fit recovers from the spots, and the length bias a mis-set camera constant leaves in the scoring while the angles stay put, and the forbidden reflection that double diffraction puts on a real plate at exactly the radius a genuine one would occupy.
- {doc}`Orientation-relationship correspondence <generated/transformation>` &mdash; Index-correspondence identities for named orientation relationships: mapping parent planes and directions to their product-phase counterparts, with rationalized indices and angular residuals, the misorientation representation used for EBSD comparison, and the recovery of a relationship and its parallelism statement from measured parent/child orientation pairs.
- {doc}`Composable visualization primitives <generated/visualization>` &mdash; Geometric guarantees of the visualization layer: a placement transform that reproduces the crystal-to-sample map, the orientation-relationship placement that makes parallel directions coincide in one world frame, a scene bond-length measurement checked against the exact NaCl-type a/2 distance, and the OR stereogram plotting a parallelism as one point and one circle for every variant.
- {doc}`Workbench service layer <generated/workbench-service-layer>` &mdash; The three quantitative claims the workbench user guide makes, each checked against a value fixed independently of this code: the Kurdjumov-Sachs packet structure and intervariant spectrum from Morito et al., the closure of the m.r.d. scale as an exact identity, the assertion a Miller component label makes about where its poles land, and the crystal viewer's claim that its camera is an orientation.

```{toctree}
:maxdepth: 1
:hidden:

generated/reference_frames
generated/crystal_geometry
generated/orientation
generated/orientation-representations
generated/diffraction
generated/ebsd
generated/composite-diffraction
generated/convergent-beam-diffraction
generated/dynamical-cbed-and-symmetry
generated/texture
generated/ipf-coloring
generated/random-disorientation
generated/elastic-anisotropy
generated/schmid-and-taylor
generated/pole-figure-sampling
generated/directional-statistics
generated/ghost-problem
generated/pole-figure-arithmetic
generated/kearns-parameter
generated/tem_tilt_navigation
generated/saed_practice_patterns
generated/transformation
generated/visualization
generated/workbench-service-layer
```
