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
| `orientation-euler-matrix-roundtrip` | 0.0000 | 0.0000 | deg | ✅ |
| `orientation-sigma3-disorientation` | 60.0000 | 60.0000 | deg | ✅ |
| `core-orientation-equal-volume-charts-agree-on-the-so3-volume` | [1.000000000000, 1.000000000000, 1.000000000000] | [1.000000000000, 1.000000000000, 1.000000000000] | &mdash; | ✅ |
| `core-orientation-ideal-indices-invert-the-plane-direction-construction` | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | &mdash; | ✅ |
| `diffraction-ni-111-two-theta` | 44.496 | 44.496 | deg | ✅ |
| `diffraction-ni-111-kikuchi-band-width` | 2.4188 | 2.4187 | deg | ✅ |
| `diffraction-gnomonic-zone-axis-radius` | 1.000000000000 | 1.000000000000 | &mdash; | ✅ |
| `diffraction-march-dollase-family-factor` | 2.309132723130 | 2.309132723130 | &mdash; | ✅ |
| `diffraction-march-dollase-normalization` | 1.000000000 | 1.000000000 | &mdash; | ✅ |
| `diffraction-odf-weighted-random-texture` | 0.9999 | 1.0000 | &mdash; | ✅ |
| `diffraction-kikuchi-map-zone-axis-tilt-angles` | [45.000000, 54.735610, 35.264390] | [45.000000, 54.735610, 35.264390] | deg | ✅ |
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
| `diffraction-cbed-aluminium-extinction-distances-at-100kv` | [555.2, 663.9, 1062.5] | [556.0, 673.0, 1057.0] | angstrom | ✅ |
| `diffraction-cbed-two-beam-thickness-inverts-the-fringe-relation` | [2000.000000, 500.000000] | [2000.000000, 500.000000] | angstrom | ✅ |
| `diffraction-dynamical-two-beam-limit-of-the-many-beam-solver` | 1.14e-15 | 0.00e+00 | &mdash; | ✅ |
| `diffraction-dynamical-intensity-is-conserved-without-absorption` | 2.11e-15 | 0.00e+00 | &mdash; | ✅ |
| `diffraction-holz-strain-and-wavelength-are-exactly-degenerate` | 1.56e-17 | 0.00e+00 | radian | ✅ |
| `diffraction-groups-construction-yields-buxtons-thirty-one` | 31 | 31 | &mdash; | ✅ |
| `diffraction-groups-friedel-observation-splits-the-point-groups` | [21, 11] | [21, 11] | &mdash; | ✅ |
| `diffraction-groups-zincblende-down-001-loses-the-two-fold` | [1, 1, 1, 0, 1, 1, 1] | [1, 1, 1, 0, 1, 1, 1] | &mdash; | ✅ |
| `texture-gaussian-kernel-normalization-and-halfwidth` | [1.000000, 0.500000] | [1.000000, 0.500000] | &mdash; | ✅ |
| `texture-uniform-odf-pole-density-is-one-mrd` | [1.000013, 1.000248, 0.999713] | [1.000000, 1.000000, 1.000000] | m.r.d. | ✅ |
| `ipf-cubic-sector-corners-are-primaries` | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | &mdash; | ✅ |
| `ipf-cubic-closed-form-colour-113` | [1.000000, 0.000000, 0.750000] | [1.000000, 0.000000, 0.750000] | &mdash; | ✅ |
| `ipf-symmetric-equivalents-share-one-colour` | 7.08e-15 | 0.00e+00 | &mdash; | ✅ |
| `texture-pole-figure-mrd-unit-mean-density` | [1.000000000000, -0.000000000000] | [1.000000000000, 0.000000000000] | m.r.d. | ✅ |
| `texture-pole-figure-resampling-and-addition-identities` | [0.000000000000, 2.000000000000] | [0.000000000000, 2.000000000000] | m.r.d. | ✅ |
| `tem-tilt-001-to-011-travel` | 45.000 | 45.000 | deg | ✅ |
| `tem-tilt-001-to-111-travel` | 54.7356 | 54.7356 | deg | ✅ |
| `tem-holder-accessible-solid-angle` | 1.04720 | 1.04720 | sr | ✅ |
| `tem-diffraction-rotation-residual` | 5.000000 | 5.000000 | deg | ✅ |
| `tem-observation-stabilizer-cubic-001` | 8 | 8 | operators | ✅ |
| `tem-symmetry-orbit-multiplicity` | 48 | 48 | directions | ✅ |
| `tem-indexed-orientation-identity` | 0.00e+00 | 0.00e+00 | dimensionless | ✅ |
| `tem-self-calibrated-diffraction-rotation` | 37.000000 | 37.000000 | deg | ✅ |
| `or-ks-plane-correspondence-identity` | [0.0000, 1.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 1.0000, 0.0000] | indices, deg | ✅ |
| `or-bain-direction-correspondence-identity` | [1.0000, 0.0000, 0.0000, 0.0000] | [1.0000, 0.0000, 0.0000, 0.0000] | indices, deg | ✅ |
| `or-ks-misorientation-representation` | [42.8478, 0.9679, 0.1776, 0.1776] | [42.8500, 0.9679, 0.1776, 0.1776] | deg, axis components | ✅ |
| `or-fit-recovers-gt-from-ks-nominal` | [2.4037, 0.0000] | [2.4037, 0.0000] | deg | ✅ |
| `or-ks-same-parent-boundary-fingerprint` | [0.0000, 0.0000] | [0.0000, 0.0000] | deg | ✅ |
| `or-ks-identified-from-measured-orientations` | [0.0000, 5.2644] | [0.0000, 5.2600] | deg | ✅ |
| `or-ks-parallelism-statement-from-rotation` | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | indices, deg | ✅ |
| `or-ks-variant-correspondence-packets` | [24, 4, 6, 6, 6] | [24, 4, 6, 6, 6] | counts | ✅ |
| `viz-transform-crystal-to-sample-consistency` | 0.0000 | 0.0000 | &mdash; | ✅ |
| `viz-or-parallel-direction-alignment` | 1.0000 | 1.0000 | &mdash; | ✅ |
| `viz-scene-bond-length-halite-identity` | 2.0000 | 2.0000 | angstrom | ✅ |

## Example groups

- {doc}`Reference Frames And Frame Transforms <generated/reference_frames>` &mdash; Creating standard frames, declaring frame relationships in words, and letting the frame graph compose multi-step chains — with the rotation angles, components, and invariants checked against exact analytic values. The last two examples pin the IUCr notation convention: the reciprocal star marks the basis, never the indices.
- {doc}`Crystal geometry: angles, spacings, and multiplicities <generated/crystal_geometry>` &mdash; Interplanar and interdirection angles, interplanar spacings, and symmetry multiplicities for cubic and hexagonal phases. Each result is checked against an analytic identity for the relevant crystal system.
- {doc}`Orientations and disorientation angles <generated/orientation>` &mdash; Round-trip consistency of orientation representations and symmetry-reduced disorientation angles, checked against exact identities and the Sigma 3 twin reference.
- {doc}`Orientation representations <generated/orientation-representations>` &mdash; The constants and identities behind the equal-volume charts of SO(3), and the inversion that names an orientation as a (hkl)[uvw] texture component.
- {doc}`Diffraction geometry <generated/diffraction>` &mdash; Powder scattering angles from PyTex interplanar spacings via Bragg's law, Kikuchi band and zone-axis geometry in the gnomonic projection, zone-axis routing on a stereographic Kikuchi map, and preferred-orientation corrections to powder intensities — each checked against a standard reference value or a closed-form identity.
- {doc}`EBSD microstructure <generated/ebsd>` &mdash; Lattice curvature and geometrically necessary dislocation density recovered from a planted orientation gradient, checked against the closed-form Nye result.
- {doc}`Composite OR diffraction <generated/composite-diffraction>` &mdash; Numerical cornerstones of composite orientation-relationship SAED simulation: the relativistic electron wavelength against the standard 200 kV value, the exactness of the Kurdjumov-Sachs child-zone mapping, and the two defining Burgers beta->alpha signatures (exact basal zone and the {110}/(0002) near-coincidence), plus the identities the exported reflection table must satisfy, and the exact halfway position of the double-diffraction Si 002 spot.
- {doc}`Convergent-beam diffraction <generated/convergent-beam-diffraction>` &mdash; The absolute scale of the two-beam extinction distance, checked against a published table, and the fringe analysis that measures a foil thickness without needing one.
- {doc}`Dynamical CBED and symmetry determination <generated/dynamical-cbed-and-symmetry>` &mdash; The exact limits that calibrate a many-beam calculation, the HOLZ degeneracy that makes voltage calibration mandatory, and the diffraction-group construction that determines a point group including its centre of symmetry.
- {doc}`Texture kernels <generated/texture>` &mdash; Analytic identities of the SO(3) kernel surface - normalization (A_0 = 1) and the halfwidth definition - together with the m.r.d. scale on which pole densities are reported, all computed live.
- {doc}`Inverse-pole-figure colouring <generated/ipf-coloring>` &mdash; What an IPF colour actually is, checked against hand-derived values: the sector corners colour to exact primaries, a direction on the [001]-[111] edge colours to exactly (1, 0, 3/4) by the closed form, and every symmetric equivalent shares one colour.
- {doc}`Pole-figure arithmetic <generated/pole-figure-arithmetic>` &mdash; Exact identities behind comparing two pole figures: the multiples-of-random scale, resampling onto a shared support, and the additivity of pole densities.
- {doc}`TEM tilt navigation <generated/tem_tilt_navigation>` &mdash; Holder tilts that bring a target zone axis onto the electron beam: analytic interzonal travel for the standard cubic transitions, the closed-form solid angle a double-tilt holder commands, the cost of an uncalibrated diffraction rotation, and the group-order counts that decide whether a single indexed pattern leaves a real ambiguity.
- {doc}`Orientation-relationship correspondence <generated/transformation>` &mdash; Index-correspondence identities for named orientation relationships: mapping parent planes and directions to their product-phase counterparts, with rationalized indices and angular residuals, the misorientation representation used for EBSD comparison, and the recovery of a relationship and its parallelism statement from measured parent/child orientation pairs.
- {doc}`Composable visualization primitives <generated/visualization>` &mdash; Geometric guarantees of the visualization layer: a placement transform that reproduces the crystal-to-sample map, the orientation-relationship placement that makes parallel directions coincide in one world frame, and a scene bond-length measurement checked against the exact NaCl-type a/2 distance.

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
generated/pole-figure-arithmetic
generated/tem_tilt_navigation
generated/transformation
generated/visualization
```
