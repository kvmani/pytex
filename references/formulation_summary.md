# Formulation Summary

Future tasks should consult this document before implementing or modifying formulas, examples, or figure geometry that rely on the PDFs in `references/`.

## Notation Crosswalk

| Topic | Source-preferred notation | PyTex-facing note |
| --- | --- | --- |
| Bunge Euler angles | `(phi1, Phi, phi2)` | Keep this ordering in docs and examples. |
| Rotation axis-angle | `(n_hat, omega)` | State whether the rotation is interpreted actively or passively. |
| Quaternion | `(q0, q1, q2, q3)` | PyTex stores quaternions in `(w, x, y, z)` order. |
| Harmonic ODF expansion | `f(g) = sum C_l^{mu nu} T_l^{mu nu}(g)` or Wigner-`D` equivalents | PyTex implements a band-limited real harmonic basis and symmetry-projects it numerically. |
| Hexagonal plane indices | `(h k i l)` | Enforce `i = -(h + k)`. |
| Hexagonal direction indices | `[U V T W]` or `[u v t w]` | PyTex currently exposes `UVTW`/`uvtw` forms but stores reduced 3-index forms internally. |
| Reciprocal lattice vector | `g_hkl` | In PyTex, this corresponds to reciprocal-basis components tied to a `Phase`. |

## Core Formulas

### 1. Reciprocal Basis

Source:
`crystallographY_calcualtions.pdf`, book pp. 10-11 (PDF pp. 20-21)

For direct basis vectors `a`, `b`, `c`:

```text
a* = (b x c) / (a . (b x c))
b* = (c x a) / (a . (b x c))
c* = (a x b) / (a . (b x c))
```

The defining orthogonality rule is:

```text
a_i . a_j* = delta_ij
```

Implementation consequence:

- Reciprocal-basis routines should preserve the explicit dual-basis semantics.
- Docs and tests should state that reciprocal vectors live in a different basis even when represented in the same Cartesian embedding.

### 2. Reciprocal Vector, Plane Normal, and Interplanar Spacing

Source:
`crystallographY_calcualtions.pdf`, book pp. 11-14 (PDF pp. 21-23)

For the reciprocal-lattice vector

```text
g_hkl = h a* + k b* + l c*
```

the key geometric relations are:

```text
g_hkl is normal to the plane (hkl)
|g_hkl| = 1 / d_hkl
```

Implementation consequence:

- `CrystalPlane.normal` and `CrystalPlane.d_spacing_angstrom` should be documented together, not separately.
- Validation examples should explicitly show the reciprocal-vector magnitude before converting to `d_hkl`.

### 3. Direct <-> Reciprocal Component Transforms

Source:
`crystallographY_calcualtions.pdf`, book pp. 16-18 (PDF pp. 25-27)

Using the direct metric tensor `g_ij` and reciprocal metric tensor `g*_ij`:

```text
p*_m = p_i g_im
p_i = p*_m g*_mi
```

Implementation consequence:

- Future vector APIs should expose direct-to-reciprocal and reciprocal-to-direct component transforms as named operations instead of leaving them as ad hoc matrix multiplications in user code.

### 4. Hexagonal Plane Conversion

Source:
`hexagnoal 4index mathematics.pdf`, pp. 1-2

For plane indices:

```text
i = -(h + k)
(h k l) -> (h k i l)
```

Implementation consequence:

- Plane conversions should preserve the supernumerary index in docs and examples, even if PyTex stores the canonical 3-index form internally.

### 5. Hexagonal Direction Conversion

Source:
`hexagnoal 4index mathematics.pdf`, pp. 2-3

From 3-index `[u v w]` to 4-index `[U V T W]`:

```text
U = (2u - v) / 3
V = (2v - u) / 3
T = -(u + v) / 3
W = w
```

Inverse transform:

```text
u = 2U + V
v = 2V + U
w = W
```

Constraint:

```text
U + V + T = 0
```

Implementation consequence:

- Converters should reduce to the smallest integer tuple after applying the formula.
- Test cases should include both divisible-by-3 and non-divisible-by-3 direction inputs.

### 6. Hexagonal Four-Index Zone Law

Source:
`hexagnoal 4index mathematics.pdf`, pp. 2-4

The four-index zone relation can be written as:

```text
h U + k V + i T + l W = 0
```

with

```text
i = -(h + k)
T = -(U + V)
```

Implementation consequence:

- PyTex currently has `ZoneAxis.contains_miller_index` for 3-index objects.
- A future Miller-Bravais API should expose this four-index zone law directly for hexagonal teaching and validation workflows.

### 7. Hexagonal Interplanar Spacing

Source:
`Kelly & Groves.pdf`, Appendix 3, pp. 469-472

For the hexagonal lattice:

```text
1 / d_hkl^2 = (4/3) * (h^2 + h k + k^2) / a^2 + l^2 / c^2
```

Implementation consequence:

- Hexagonal d-spacing tests should cite this formula explicitly instead of using only numeric parity checks.

### 8. Bragg Law

Source:
`williamsandcarter.pdf`, pp. 78-79

```text
n lambda = 2 d sin(theta_B)
```

Implementation consequence:

- Powder-XRD docs should distinguish `theta_B` from reported `2theta`.
- Any detector or plotting workflow should state whether its angular axis is `theta`, `2theta`, or detector-plane coordinates.

### 9. Bunge Euler Rotation Matrix

Sources:

- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 47-49
- `MathsOfrotations_RolletDegraef.pdf`, p. 21

PyTex docs should keep the Bunge tuple ordered as:

```text
(phi1, Phi, phi2)
```

and should present the matrix formula using the same symbol order whenever the docs explain Euler input/output semantics.

Implementation consequence:

- Public documentation and tests should avoid switching symbol names or axis order casually.
- When a page uses the matrix form, it should also state the mapping direction and frame meaning.

### 10. Axis-Angle and Quaternion

Source:
`MathsOfrotations_RolletDegraef.pdf`, pp. 6-7, 23-25

Axis-angle to quaternion:

```text
q = (cos(omega / 2), n_hat sin(omega / 2))
```

Implementation consequence:

- Quaternion docs should say explicitly that PyTex stores this as `(w, x, y, z)`.
- When code canonicalizes the sign of a quaternion, the docs should explain why equivalent `q` and `-q` need a canonical representative.

### 11. Kikuchi Pattern as Gnomonic Projection

Source:
`Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, p. 153

Implementation consequence:

- Kikuchi figures should be drawn as gnomonic-projection geometry, not as generic flat sketches.
- Pattern-center and band-normal diagrams should state the projection model explicitly.

### 12. Harmonic ODF Series Expansion

Source:
`Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`,
pp. 105-107

The classical texture-analysis statement is that the ODF can be expanded in symmetrized
generalized spherical harmonics and the PF coefficients are linked linearly to the ODF
coefficients after truncation.

Implementation consequence:

- Harmonic reconstruction docs should state the selected bandlimit explicitly.
- Crystal and specimen symmetry must be treated as invariance constraints on the basis, not
  as after-the-fact cosmetic reductions.
- The implementation should document clearly whether it uses a closed-form coefficient
  relation or a numerical forward operator built on the same harmonic basis.

### 13. Crystal And Specimen Symmetry Invariance For ODFs

Source basis:
texture-analysis convention summarized by the same series-expansion discussion together
with the repository frame model

The invariance relations are:

```text
f(g h) = f(g)   for crystal symmetry h
f(s g) = f(g)   for specimen symmetry s
```

Implementation consequence:

- Future harmonic or orientation-space algorithms should preserve the right action for the
  crystal group and the left action for the specimen group.
- Docs should never describe crystal and specimen symmetry as interchangeable in ODF
  reconstruction.

### 14. Antipodal Pole Figures And Even Harmonic Degrees

Source:
`Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`,
pp. 105-109

Implementation consequence:

- Antipodal diffraction PF workflows should default to even harmonic degrees only unless a
  task has a clear reason to recover odd terms.
- Validation docs should say explicitly when a degree-selection rule is driven by the
  missing odd-order information in antipodal PF data.

## Worked Examples Worth Reusing

| Topic | Source | Use in PyTex |
| --- | --- | --- |
| Four-index plane and zone examples | `hexagnoal 4index mathematics.pdf`, pp. 3-4 | Future test vectors and tutorial examples. |
| Reciprocal-basis and d-spacing derivations | `crystallographY_calcualtions.pdf`, book pp. 10-18 | Foundation for reciprocal-space docs and tests. |
| Orientation descriptor comparisons | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 34-50 | Basis for orientation tutorial tables and figure legends. |
| Harmonic ODF truncation and odd-order discussion | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 105-109 | Basis for harmonic inversion docs, test cases, and conservative antipodal defaults. |
| Rotation round-robin convention example | `MathsOfrotations_RolletDegraef.pdf`, pp. 1-3 | Good candidate for future PyTex parity tests across representations. |

## Immediate PyTex Relevance

- `src/pytex/core/hexagonal.py`
- `src/pytex/core/miller.py`
- `src/pytex/core/lattice.py`
- `src/pytex/core/orientation.py`
- `src/pytex/texture/harmonics.py`
- `src/pytex/diffraction/xrd.py`
- `src/pytex/diffraction/saed.py`
- `docs/figures/orientation_conventions.svg`
- `docs/figures/reference_frames.svg`
- `docs/figures/ipf_sector_geometry_detailed.svg`
- `docs/testing/automated_test_cases.md`
- `docs/site/workflows/harmonic_odf_reconstruction.md`

## References

### Normative

- [Notation And Conventions](../docs/standards/notation_and_conventions.md)
- [Terminology And Symbol Registry](../docs/standards/terminology_and_symbol_registry.md)

### Informative

- [Reference Index](reference_index.md)
- [Feature Opportunities](feature_opportunities.md)
