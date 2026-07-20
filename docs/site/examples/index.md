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
| `cubic-angle-100-110` | 45.0000 | 45.0000 | deg | ✅ |
| `cubic-angle-100-111` | 54.7356 | 54.7356 | deg | ✅ |
| `cubic-angle-dir-110-111` | 35.2644 | 35.2644 | deg | ✅ |
| `cubic-dspacing-111` | 2.30940 | 2.30940 | angstrom | ✅ |
| `cubic-plane-multiplicity` | [3, 6, 4, 24] | [3, 6, 4, 24] | &mdash; | ✅ |
| `hex-angle-basal-prism` | 90.0000 | 90.0000 | deg | ✅ |
| `hex-angle-prism-prism` | 60.0000 | 60.0000 | deg | ✅ |
| `orientation-euler-matrix-roundtrip` | 0.0000 | 0.0000 | deg | ✅ |
| `orientation-sigma3-disorientation` | 60.0000 | 60.0000 | deg | ✅ |
| `diffraction-ni-111-two-theta` | 44.496 | 44.496 | deg | ✅ |
| `composite-electron-wavelength-200kv` | 0.02508 | 0.02508 | angstrom | ✅ |
| `composite-ks-exact-child-zone` | 0.0000 | 0.0000 | deg | ✅ |
| `composite-burgers-exact-basal-zone` | 0.0000 | 0.0000 | deg | ✅ |
| `composite-burgers-110-0002-coincidence` | 0.15450 | 0.15450 | mm | ✅ |
| `texture-gaussian-kernel-normalization-and-halfwidth` | [1.000000, 0.500000] | [1.000000, 0.500000] | &mdash; | ✅ |
| `or-ks-plane-correspondence-identity` | [0.0000, 1.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 1.0000, 0.0000] | indices, deg | ✅ |
| `or-bain-direction-correspondence-identity` | [1.0000, 0.0000, 0.0000, 0.0000] | [1.0000, 0.0000, 0.0000, 0.0000] | indices, deg | ✅ |
| `or-ks-misorientation-representation` | [42.8478, 0.9679, 0.1776, 0.1776] | [42.8500, 0.9679, 0.1776, 0.1776] | deg, axis components | ✅ |
| `or-fit-recovers-gt-from-ks-nominal` | [2.4037, 0.0000] | [2.4037, 0.0000] | deg | ✅ |
| `viz-transform-crystal-to-sample-consistency` | 0.0000 | 0.0000 | &mdash; | ✅ |
| `viz-or-parallel-direction-alignment` | 1.0000 | 1.0000 | &mdash; | ✅ |
| `viz-scene-bond-length-halite-identity` | 2.0000 | 2.0000 | angstrom | ✅ |

## Example groups

- {doc}`Crystal geometry: angles, spacings, and multiplicities <generated/crystal_geometry>` &mdash; Interplanar and interdirection angles, interplanar spacings, and symmetry multiplicities for cubic and hexagonal phases. Each result is checked against an analytic identity for the relevant crystal system.
- {doc}`Orientations and disorientation angles <generated/orientation>` &mdash; Round-trip consistency of orientation representations and symmetry-reduced disorientation angles, checked against exact identities and the Sigma 3 twin reference.
- {doc}`Diffraction geometry <generated/diffraction>` &mdash; Powder scattering angles derived from PyTex interplanar spacings via Bragg's law, checked against a standard reference reflection position.
- {doc}`Composite OR diffraction <generated/composite-diffraction>` &mdash; Numerical cornerstones of composite orientation-relationship SAED simulation: the relativistic electron wavelength against the standard 200 kV value, the exactness of the Kurdjumov-Sachs child-zone mapping, and the two defining Burgers beta->alpha signatures (exact basal zone and the {110}/(0002) near-coincidence).
- {doc}`Texture kernels <generated/texture>` &mdash; Analytic identities of the SO(3) kernel surface: normalization (A_0 = 1) and the halfwidth definition, computed live.
- {doc}`Orientation-relationship correspondence <generated/transformation>` &mdash; Index-correspondence identities for named orientation relationships: mapping parent planes and directions to their product-phase counterparts, with rationalized indices and angular residuals, and the misorientation representation used for EBSD comparison.
- {doc}`Composable visualization primitives <generated/visualization>` &mdash; Geometric guarantees of the visualization layer: a placement transform that reproduces the crystal-to-sample map, the orientation-relationship placement that makes parallel directions coincide in one world frame, and a scene bond-length measurement checked against the exact NaCl-type a/2 distance.

```{toctree}
:maxdepth: 1
:hidden:

generated/crystal_geometry
generated/orientation
generated/diffraction
generated/composite-diffraction
generated/texture
generated/transformation
generated/visualization
```
