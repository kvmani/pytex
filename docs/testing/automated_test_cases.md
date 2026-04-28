# Automated Test Cases

PyTex uses automated tests as a scientific audit surface, not merely as a CI gate.

This page is written for domain experts who want to inspect whether a method, formula, or
convention is scientifically credible without first reading Python implementation details.

## What This Page Must Do

For every major crystallographic conversion, orientation convention, or diffraction
calculation, the validation surface should make five things visible in one place:

- the authoritative reference and exact page anchor
- the source notation and governing formula
- the algorithm currently implemented by PyTex
- the source-derived expected value for a worked example
- the actual value returned by the current code

The goal is not to decorate the test suite. The goal is to make the numerical pathway
auditable by a crystallographer, diffraction specialist, or EBSD expert who wants to
challenge assumptions, notation choices, or interpretation directly.

## Review Pattern

Each documented case on this page follows the same audit structure:

1. `Code surface`
   The method, class, and automated tests that define the behavior.
2. `Reference basis`
   Exact source and page numbers.
3. `Source formulation`
   The formula or convention written as close to the source notation as practical.
4. `Algorithm used by PyTex`
   A prose summary of the implementation path, without requiring code reading.
5. `Audit example`
   A worked example with source-derived expected values, current code outputs, and an
   interpretation column.

## Case 1: Euler Conversion Audit

### Scope

This case documents how PyTex handles Bunge Euler angles, matrix conversion, and
Matthies or ABG export.

### Code surface

- `src/pytex/core/orientation.py`
- `tests/unit/test_orientation.py`
- `tests/unit/test_batches.py`

### Reference basis

- Rowenhorst et al., `references/MathsOfrotations_RolletDegraef.pdf`, Appendix A.1 and A.5, pp. 21-22
- Randle and Engler, `references/Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 47-49

### Source formulation

Randle and Engler give the Bunge convention as the ordered triple

```text
(phi1, Phi, phi2)
```

with the three rotations transforming the specimen coordinate system into the crystal
coordinate system.

Rowenhorst Appendix A.1 gives the explicit Bunge `ZXZ` passive rotation matrix:

```text
g =
[ c1 c2 - s1 c s2      s1 c2 + c1 c s2      s s2 ]
[ -c1 s2 - s1 c c2     -s1 s2 + c1 c c2     s c2 ]
[ s1 s                  -c1 s                 c    ]
```

with

```text
c1 = cos(phi1), s1 = sin(phi1)
c2 = cos(phi2), s2 = sin(phi2)
c  = cos(Phi),  s  = sin(Phi)
```

Randle and Engler also state the conversion relations:

Bunge to Roe:

```text
Roe_1 = phi1 - 90 deg
Roe_2 = Phi
Roe_3 = phi2 + 90 deg
```

Bunge to Kocks:

```text
Kocks_1 = phi1 - 90 deg
Kocks_2 = Phi
Kocks_3 = 90 deg - phi2
```

### Algorithm used by PyTex

PyTex currently implements this path:

1. `Rotation.from_bunge_euler(phi1, Phi, phi2)` delegates to the generic Euler
   constructor with `convention="bunge"`.
2. The constructor converts degrees to radians and composes three axis-angle
   quaternions in `z-x-z` order.
3. The stable internal representation is a unit quaternion stored as `(w, x, y, z)`.
4. `Rotation.as_matrix()` converts that quaternion into a 3x3 matrix.
5. `Rotation.to_matthies_euler()` and `Rotation.to_abg_euler()` recover a `ZYZ`
   repeated-axis triple from the matrix.

### Audit example

Randle and Engler use the ideal orientation `(123)[634]` as a worked example and give
the Bunge angles:

```text
(phi1, Phi, phi2) = (301.0 deg, 36.7 deg, 26.63 deg)
```

From the reference matrix formula above, the source-derived passive matrix is:

```text
g_reference =
[[ 0.76844944, -0.58114553,  0.26787185],
 [ 0.38349746,  0.75334535,  0.53422887],
 [-0.51226473, -0.30779971,  0.80177564]]
```

Because PyTex documents `Orientation` as a crystal-frame to specimen-frame mapping, the
matrix surfaced by the current code should be compared against the transpose of the
textbook passive matrix:

```text
g_reference_transpose =
[[ 0.76844944,  0.38349746, -0.51226473],
 [-0.58114553,  0.75334535, -0.30779971],
 [ 0.26787185,  0.53422887,  0.80177564]]
```

The source-derived angle conversions are:

```text
Roe expected   = (211.0 deg, 36.7 deg, 116.63 deg)
Kocks expected = (211.0 deg, 36.7 deg,  63.37 deg)
```

### Audit table

| Quantity | Source-derived expected value | Current code output | Interpretation |
| --- | --- | --- | --- |
| Bunge input | `(301.0, 36.7, 26.63)` | accepted exactly | Matches the authoritative example. |
| Textbook passive matrix | `[[0.76844944, -0.58114553, 0.26787185], [0.38349746, 0.75334535, 0.53422887], [-0.51226473, -0.30779971, 0.80177564]]` | not returned directly by `as_matrix()` | The source formula is written for the passive matrix convention. |
| PyTex matrix under documented crystal-to-specimen mapping | `g_reference_transpose` | `[[0.76844944, 0.38349746, -0.51226473], [-0.58114553, 0.75334535, -0.30779971], [0.26787185, 0.53422887, 0.80177564]]` | Current PyTex output matches the transpose of the textbook passive matrix. This should be stated explicitly in the docs, not inferred by readers. |
| ABG or Matthies export | `Roe-style expected = (211.0, 36.7, 116.63)` | `(211.0, 36.7, 116.63)` | Current `abg` and `matthies` exports agree numerically with the Roe-style relation for this example. |
| Kocks export | `Kocks expected = (211.0, 36.7, 63.37)` | no dedicated public result | PyTex does not currently expose a distinct Kocks conversion. Domain experts should not interpret `abg` as Kocks unless that behavior is implemented and documented explicitly. |
| Quaternion storage | source texts often use `(q0, q1, q2, q3)` | `(-0.91153311, -0.23093746, 0.21396277, 0.26456609)` from the local run | PyTex uses `(w, x, y, z)` storage. Sign-equivalent quaternions represent the same rotation, so quaternion sign alone is not a correctness failure. |

### Automated assertion currently present

- `Rotation.from_bunge_euler(0, 0, 0)` yields identity.
- Matrix round-trips preserve the represented rotation.
- Euler round-trips preserve the represented rotation.

### Review prompts for domain experts

- Is the crystal-to-specimen matrix transpose relation the correct documented mapping for
  PyTex's public orientation meaning?
- Should PyTex surface a distinct Kocks conversion instead of only `matthies` and `abg`?
- Should the repository use one canonical sign convention for published quaternions in
  docs and examples?

## Case 2: Hexagonal 4-Index Audit

### Scope

This case documents the conversion between three-index and four-index notation for
hexagonal directions and planes.

### Code surface

- `src/pytex/core/hexagonal.py`
- `src/pytex/core/miller.py`
- `tests/unit/test_hexagonal.py`
- `tests/unit/test_miller_objects.py`

### Reference basis

- Donnay, `references/hexagnoal 4index mathematics.pdf`, pp. 1-4

### Source formulation

For planes, Donnay states the supernumerary index rule:

```text
i = -(h + k)
```

For zones and directions, the four-index constraint is:

```text
U + V + T = 0
```

The practical conversion used throughout PyTex is:

```text
U = (2u - v) / 3
V = (2v - u) / 3
T = -(u + v) / 3
W = w
```

Inverse:

```text
u = 2U + V
v = 2V + U
w = W
```

### Algorithm used by PyTex

PyTex currently:

1. applies the exact rational formulas above using `Fraction`
2. rescales to the smallest common integer tuple
3. reduces the resulting tuple by the common divisor
4. validates the defining constraints

This means the conversion is not implemented as floating-point rounding.

### Audit example

The minimal teaching examples used in the tests are:

```text
[1 0 0] -> [2 -1 -1 0]
[2 -1 -1 3] -> [1 0 1]
(1 0 2) -> (1 0 -1 2)
```

### Audit table

| Quantity | Source-derived expected value | Current code output | Interpretation |
| --- | --- | --- | --- |
| Direction `[1 0 0]` to four-index | `[2 -1 -1 0]` | `[2, -1, -1, 0]` | Matches the standard Miller-Bravais conversion. |
| Direction `[2 -1 -1 3]` back to three-index | `[1 0 1]` after reduction | `[1, 0, 1]` | The reduction step is essential and is handled correctly. |
| Plane `(1 0 2)` to four-index | `(1 0 -1 2)` | `[1, 0, -1, 2]` | Matches the supernumerary-index rule exactly. |
| Plane `(1 0 -1 2)` back to three-index | `(1 0 2)` | `[1, 0, 2]` | Matches the inverse plane rule. |
| Invalid four-index direction check | reject if `U + V + T != 0` | current code raises `ValueError` | Correct behavior for scientific data hygiene. |
| Invalid four-index plane check | reject if `i != -(h + k)` | current code raises `ValueError` | Correct behavior for scientific data hygiene. |

### Automated assertion currently present

- three-index to four-index direction conversion
- four-index to three-index round-trip
- three-index to four-index plane round-trip
- invalid four-index rejection

### Review prompts for domain experts

- Should PyTex keep three-index forms as the primary storage for stable public types, or
  should a first-class Miller-Bravais object surface be added?
- Are additional canonical teaching examples from Donnay pp. 3-4 worth promoting into
  the automated suite?

## Case 3: Reciprocal Vector And Interplanar Spacing Audit

### Scope

This case documents how PyTex computes `d_hkl` from the reciprocal basis and checks
cubic and hexagonal examples against analytic formulas.

### Code surface

- `src/pytex/core/lattice.py`
- `src/pytex/core/miller.py`
- `tests/unit/test_miller_objects.py`

### Reference basis

- De Graef, `references/crystallographY_calcualtions.pdf`, book pp. 10-14
- Kelly and Groves, `references/Kelly & Groves.pdf`, Appendix 3, pp. 469-472

### Source formulation

For the reciprocal vector:

```text
g_hkl = h a* + k b* + l c*
```

with

```text
|g_hkl| = 1 / d_hkl
```

For a hexagonal lattice:

```text
1 / d_hkl^2 = (4/3) * (h^2 + h k + k^2) / a^2 + l^2 / c^2
```

### Algorithm used by PyTex

PyTex currently:

1. constructs the reciprocal basis from the lattice
2. forms the reciprocal vector for each Miller triplet
3. computes `d_hkl = 1 / ||g_hkl||`
4. uses vectorized parity checks against the analytic cubic and hexagonal formulas

### Audit examples

The current automated checks use:

- cubic lattice with `a = 4.0 A`
- hexagonal lattice with `a = 2.95 A`, `c = 4.68 A`

The source-derived expected values are:

```text
cubic d_200 = 4.0 / 2 = 2.0 A
cubic d_111 = 4.0 / sqrt(3) = 2.3094010767585034 A
hexagonal d_100 = 2.5547749411640943 A
hexagonal d_101 = 2.2424130309239043 A
hexagonal d_002 = 2.34 A
```

### Audit table

| Quantity | Source-derived expected value | Current code output | Interpretation |
| --- | --- | --- | --- |
| Cubic `d_200` for `a = 4.0 A` | `2.0` | `2.0` | Exact match. |
| Cubic `d_111` for `a = 4.0 A` | `2.3094010767585034` | `2.3094010767585034` | Exact match to the printed precision. |
| Hexagonal `d_100` for `a = 2.95 A`, `c = 4.68 A` | `2.5547749411640943` | `2.5547749411640943` | Exact match to the printed precision. |
| Hexagonal `d_101` for `a = 2.95 A`, `c = 4.68 A` | `2.2424130309239043` | `2.2424130309239048` | Agreement within floating-point tolerance. |
| Hexagonal `d_002` for `a = 2.95 A`, `c = 4.68 A` | `2.34` | `2.34` | Exact match. |

### Automated assertion currently present

- cubic `d_hkl` parity against simple closed-form formulas
- hexagonal `d_hkl` parity against the Appendix 3 formula
- plane-normal and direction-vector checks for cubic examples

### Review prompts for domain experts

- Should a future test page also show the reciprocal vector itself before reducing to
  `d_hkl`?
- Should direct and reciprocal metric tensor surfaces be exposed explicitly in the public
  API rather than only through lattice helpers?

## Case 4: Powder XRD And Bragg-Law Audit

### Scope

This case documents the `2theta` and intensity pathway used by PyTex for the current
pedagogical powder-XRD generator.

### Code surface

- `src/pytex/diffraction/xrd.py`
- `tests/unit/test_xrd_saed_and_styles.py`
- `tests/unit/test_diffraction_external_baselines.py`

### Reference basis

- Williams and Carter, `references/williamsandcarter.pdf`, pp. 78-79

### Source formulation

The Bragg law used to determine the scattering angle is:

```text
n lambda = 2 d sin(theta_B)
```

PyTex stores and plots the full diffraction angle:

```text
two_theta = 2 theta_B
```

For the cubic test phase with `a = 3.0 A`, the `(111)` spacing is:

```text
d_111 = a / sqrt(3) = 1.7320508075688774 A
```

Using Cu Ka radiation:

```text
lambda = 1.5406 A
```

the source-derived expected angle is:

```text
two_theta = 2 * arcsin(lambda / (2 d_111))
          = 52.81250368742122 deg
```

### Algorithm used by PyTex

PyTex currently:

1. enumerates candidate `hkl` triplets
2. computes `d_hkl` from the reciprocal basis
3. applies the Bragg condition to obtain `2theta`
4. computes multiplicity from symmetry-equivalent reciprocal vectors
5. computes a pedagogical intensity using:
   structure-factor amplitude squared, multiplicity, and a Lorentz-polarization factor
6. broadens reflections onto a grid and normalizes the final intensity profile to a maximum
   of `1`

### Audit table

| Quantity | Source-derived expected value | Current code output | Interpretation |
| --- | --- | --- | --- |
| Cubic `(111)` spacing for `a = 3.0 A` | `1.7320508075688774 A` | `1.7320508075688774 A` | Exact match. |
| `(111)` Bragg angle for Cu Ka radiation | `52.81250368742122 deg` | `52.81250368742122 deg` | Exact match to the printed precision. |
| `(111)` multiplicity in cubic `m-3m` | `8` | `8` | Correct for the current symmetry handling. |
| Raw reflection intensity in current pedagogical model | model-dependent | `12084.825189587242` | Not a handbook intensity. It reflects the current atomic-number plus Lorentz-polarization model. |
| Peak-normalized intensity-grid maximum | `1.0` after normalization | `1.0` | Correct for the current grid normalization contract. |

### Automated assertion currently present

- the `(111)` reflection must be present in the generated pattern
- the `2theta` value must agree with Bragg-law expectation
- the normalized intensity grid must remain within the expected bounds

### Review prompts for domain experts

- Is the current pedagogical intensity model clearly labeled enough as an approximation
  rather than a full physically rigorous X-ray scattering model?
- Should future test documentation separate angle correctness from intensity-model
  correctness more sharply?

## Case 5: Harmonic ODF Reconstruction Audit

### Scope

This case documents the current harmonic PF-to-ODF reconstruction route, including the
band-limited basis, the crystal and specimen symmetry handling, and the even-degree policy
used for antipodal pole figures.

### Code surface

- `src/pytex/texture/harmonics.py`
- `tests/unit/test_harmonic_odf.py`
- `docs/site/workflows/harmonic_odf_reconstruction.md`

### Reference basis

- Randle and Engler, `references/Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 49-50
- Randle and Engler, `references/Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 105-107

### Source formulation

Randle and Engler describe the ODF series-expansion method through symmetrized generalized
spherical harmonics and note that PF and ODF coefficients are linked linearly after
truncation.

PyTex implements that doctrine in Bunge Euler space through real and imaginary parts of
Wigner `D` terms:

```text
D^l_mn(phi1, Phi, phi2) = exp(-i m phi1) d^l_mn(Phi) exp(-i n phi2)
```

and enforces the ODF invariance conditions

```text
f(g h) = f(g)   for crystal symmetry h in C
f(s g) = f(g)   for specimen symmetry s in S
```

through numerical group averaging of the basis.

For antipodal pole figures, odd harmonic orders are not observable uniquely. PyTex
therefore defaults to `even_degrees_only=True` unless the caller overrides it.

### Algorithm used by PyTex

PyTex currently:

1. builds a midpoint Bunge quadrature over full `SO(3)` with weights proportional to
   `sin(Phi)`
2. evaluates raw harmonic terms up to the requested `degree_bandlimit`
3. projects those terms through the crystal and specimen symmetry operators by left or
   right group averaging
4. orthonormalizes the symmetry-projected basis numerically with the weighted Gram matrix
5. builds the PF forward operator against the current pole-density kernel semantics
6. solves the retained harmonic coefficient system by Tikhonov-regularized least squares

### Audit example

The current harmonic regression case uses:

- cubic crystal symmetry `m-3m`
- no specimen symmetry
- antipodal pole figures, so even degrees only
- `degree_bandlimit = 4`
- Bunge quadrature steps `(60 deg, 60 deg, 60 deg)`
- a synthetic but source-governed harmonic ODF constructed from a random-texture baseline
  plus a small retained harmonic perturbation

The local run used the coefficient head:

```text
true coefficient head =
[0.0, 0.05, -0.0438932183, 0.0, 0.0]
```

and reconstructed pole figures for `(100)`, `(110)`, and `(111)` on a regular specimen
direction grid.

### Audit table

| Quantity | Source-derived expected value | Current code output | Interpretation |
| --- | --- | --- | --- |
| Harmonic truncation rule | use a finite `l_max` | `degree_bandlimit = 4` | Matches the classical series-expansion doctrine. |
| Antipodal PF policy | odd orders not uniquely observable | `even_degrees_only = True` | Correct conservative default for diffraction PFs. |
| Crystal symmetry handling | right action by the crystal group | crystal symmetry order `24` | Matches the `m-3m` proper rotation group used in the phase. |
| Specimen symmetry handling | left action by the specimen group | specimen symmetry order `1` | No specimen symmetry was imposed in this audit case. |
| Retained harmonic basis size | symmetry-reduced band-limited basis, not the raw term count | `basis_size = 10` from `raw_basis_size = 211` | The numerical basis is strongly compressed by symmetry and weighted orthonormalization. |
| Quadrature size | full `SO(3)` sampling with Haar weighting | `quadrature_size = 108` | Coarse but explicit quadrature used for the regression case. |
| Random-texture mean | mean ODF should remain close to `1` in this constructed example | `mean_density = 0.9999984758` | The reconstructed mean is effectively preserved. |
| Density bounds | small harmonic perturbation around random texture | `minimum_density = 0.9193180022`, `maximum_density = 1.0806786729` | The fitted ODF preserves the intended bounded perturbation. |
| Coefficient recovery | fitted coefficients should reproduce the synthetic coefficients within solver tolerance | fitted head `[-2.18e-17, 4.99996279e-02, -4.38925650e-02, 6.59e-08, -4.63e-07]` | The recovered coefficient head matches the synthetic source closely. |
| PF refit quality | reconstructed PFs should match the synthetic PF data to solver tolerance | `relative_residual_norm = 1.0201662e-06` | The harmonic inversion reproduces the generated PF data essentially exactly. |
| Linear solve conditioning | finite condition number and full retained rank | `matrix_rank = 10`, `condition_number = 28.9234296672` | The current regression case is well conditioned enough for a stable solve. |

### Automated assertion currently present

- harmonic ODF evaluation is invariant under crystal right actions and specimen left actions
- harmonic inversion reproduces synthetic PF data within tight residual tolerance
- antipodal PFs default to even harmonic degrees
- ODF section plotting accepts a harmonic ODF object directly

### Review prompts for domain experts

- Is the current explicit left/right symmetry averaging the clearest way to document the
  relationship between PyTex and the classical symmetrized harmonic formulation?
- Should PyTex eventually expose degree-selection diagnostics beyond basis size, rank, and
  condition number?
- Is the current kernel-regularized PF forward operator an acceptable stable bridge
  between measured PF sampling and harmonic ODF reconstruction, or should a closed-form
  coefficient transform be added later for parity and benchmarking?

## Foundation Feature Priority Audit

`Code surface`

- `src/pytex/texture/reconstruction.py`
- `src/pytex/core/orientation_geometry.py`
- `src/pytex/diffraction/physics.py`
- `src/pytex/ebsd/texture_workflow.py`
- `src/pytex/core/parent_reconstruction.py`
- `tests/unit/test_foundation_feature_priorities.py`

`Reference basis`

- Bunge and derivative Euler-convention relationships follow the notation discussion used in
  Randle and Engler, *Introduction to Texture Analysis*, pp. 47-49.
- The pole-density and harmonic reconstruction pathways follow the ODF/PF formulation summarized
  in `references/formulation_summary.md`.
- Structure-factor and reflection-condition behavior follows the crystallographic diffraction
  formulation summarized in `docs/tex/algorithms/foundation_feature_priorities.tex`.

`Formulas under automated review`

Pole-figure correction:

$$
I_i^{\mathrm{corr}} = s \left(\frac{I_i}{d_i} - b\right)
$$

Pole-figure residual:

$$
r_i = I_i^{\mathrm{fit}} - I_i^{\mathrm{obs}},
\qquad
\rho = \frac{\lVert r\rVert_2}{\max(\lVert I^{\mathrm{obs}}\rVert_2,\epsilon)}
$$

Bunge to Kocks:

$$
(\Psi,\Theta,\Phi)_{\mathrm{Kocks}}
=
(\phi_1 - 90^\circ,\Phi,90^\circ-\phi_2)
$$

Structure-factor intensity:

$$
I_{hkl} = m_{hkl}|F_{hkl}|^2 L_p(2\theta)
$$

EBSD texture weights:

$$
w_i^{\mathrm{norm}} = \frac{m_i w_i}{\sum_j m_j w_j}
$$

Parent candidate score:

$$
s_m = R(\Delta\theta_{m1},\Delta\theta_{m2},\ldots,\Delta\theta_{mn})
$$

`Current automated assertions`

| Area | Input | Expected output | Current code output |
| --- | --- | --- | --- |
| Texture V2 correction | PF intensities with `scale=2`, `background=0.25` | corrected intensities are finite and non-negative | passed |
| Texture V2 residual report | fitted ODF and corrected PF | residual count equals PF observation count | passed |
| Texture JSON contracts | `PoleFigure`, `ODF`, `ODFInversionReport` | round-trip payload equality | passed |
| Euler transform | Bunge `(301.0, 36.7, 26.63)` | Kocks `(211.0, 36.7, 63.37)` and Roe `(211.0, 36.7, 116.63)` | passed |
| IPF sector boundary | cubic `m-3m` symmetry | `[001]` lies inside the sector and equations are exposed | passed |
| Diffraction physics | fcc Ni-like phase | `(100)` forbidden, `(111)` allowed, positive intensity for `(111)` | passed |
| Diffraction guardrails | invalid scattering vector, angle, multiplicity, and structure-factor values | construction or evaluation raises `ValueError` before producing nonphysical intensity output | passed |
| EBSD texture workflow | three-orientation `CrystalMap` with weights `(0.2, 0.3, 0.5)` | normalized weights sum to `1.0`; one PF produced | passed |
| EBSD workflow guardrails | invalid segmentation threshold, empty sample-direction list, all-masked weights, and nonfinite result weights | construction or normalization raises `ValueError` | passed |
| EBSD JSON contracts | `CrystalMap`, `TextureReport` | round-trip payload equality | passed |
| Parent reconstruction | identity parent-child relationship | identity candidate has best score near `0 deg` | passed |
| Parent reconstruction guardrails | nonfinite ambiguity tolerance, invalid variant indices, and nonfinite report scores | construction raises `ValueError` | passed |

`Interpretation`

This audit confirms that the new foundation surfaces are wired, typed, and reconstructible. It does
not claim complete physical diffraction modeling, complete orientation-space boundary catalogs, or
full parent-reconstruction capability. The added guardrail cases confirm that the staged public
surfaces reject invalid metadata and nonphysical numerical inputs before downstream interpretation.
Algorithmic breadth remains an explicit later expansion target.

## Current Audit Findings

The documented cases above already show issues that are scientifically important even when
the current tests pass:

- The Bunge matrix in the reference texts and the matrix surfaced by PyTex differ by a
  transpose because they describe different mapping conventions. This is not a bug by
  itself, but it must be documented explicitly.
- The current `abg` or `matthies` export agrees numerically with the Roe-style relation for
  the worked orientation example. It should not be described as Kocks unless PyTex adds a
  distinct Kocks conversion and validates it.
- The powder-XRD angle pathway is source-auditable and strong. The intensity pathway is
  intentionally pedagogical and should stay labeled as such.

## Maintenance Rule

When a future task changes a formula, convention, or public numerical pathway:

- update the relevant automated tests
- update the corresponding audit case on this page
- refresh the source-derived expected values if the reference basis changes
- refresh the current code output block from a real local run
- describe any mismatch as a scientific interpretation issue, not only as a pass or fail
  test event

## References

### Normative

- <a href="testing_strategy.html">Testing Strategy</a>
- <a href="../standards/documentation_architecture.html">Documentation Architecture</a>
- <a href="../standards/development_principles.html">Development Principles</a>

### Informative

- `references/formulation_summary.md`
- `references/reference_index.md`
