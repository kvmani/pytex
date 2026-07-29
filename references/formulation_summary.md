# Formulation Summary

Future tasks should consult this document before implementing or modifying formulas, examples, or figure geometry that rely on the PDFs in `references/`.

## Notation Crosswalk

| Topic | Source-preferred notation | PyTex-facing note |
| --- | --- | --- |
| Bunge Euler angles | $(\varphi_1, \Phi, \varphi_2)$ | Keep this ordering in docs and examples. |
| Rotation axis-angle | $(\hat{\mathbf{n}}, \omega)$ | State whether the rotation is interpreted actively or passively. |
| Quaternion | $(q_0, q_1, q_2, q_3)$ | PyTex stores quaternions in `(w, x, y, z)` order. |
| Harmonic ODF expansion | $f(g) = \sum C_{l}^{\mu\nu} T_{l}^{\mu\nu}(g)$ or Wigner-$D$ equivalents | PyTex implements a band-limited real harmonic basis and symmetry-projects it numerically. |
| Hexagonal plane indices | $(hkil)$ | Enforce $i = -(h + k)$. |
| Hexagonal direction indices | $[UVTW]$ or $[uvtw]$ | PyTex currently exposes `UVTW`/`uvtw` forms but stores reduced 3-index forms internally. |
| Reciprocal lattice vector | $\mathbf{g}_{hkl}$ | In PyTex, this corresponds to reciprocal-basis components tied to a `Phase`. |

## Core Formulas

### 1. Reciprocal Basis

Source:
`crystallographY_calcualtions.pdf`, book pp. 10-11 (PDF pp. 20-21)

For direct basis vectors $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$:

$$
\mathbf{a}^{*} = \frac{\mathbf{b} \times \mathbf{c}}{\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})},
\qquad
\mathbf{b}^{*} = \frac{\mathbf{c} \times \mathbf{a}}{\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})},
\qquad
\mathbf{c}^{*} = \frac{\mathbf{a} \times \mathbf{b}}{\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})}.
$$

The defining orthogonality rule is:

$$
\mathbf{a}_i \cdot \mathbf{a}_j^{*} = \delta_{ij}
$$

Implementation consequence:

- Reciprocal-basis routines should preserve the explicit dual-basis semantics.
- Docs and tests should state that reciprocal vectors live in a different basis even when represented in the same Cartesian embedding.

### 2. Reciprocal Vector, Plane Normal, and Interplanar Spacing

Source:
`crystallographY_calcualtions.pdf`, book pp. 11-14 (PDF pp. 21-23)

For the reciprocal-lattice vector

$$
\mathbf{g}_{hkl} = h\,\mathbf{a}^{*} + k\,\mathbf{b}^{*} + l\,\mathbf{c}^{*}
$$

the key geometric relations are:

$$
\mathbf{g}_{hkl} \perp (hkl),
\qquad
\lVert \mathbf{g}_{hkl} \rVert = \frac{1}{d_{hkl}}.
$$

Implementation consequence:

- `CrystalPlane.normal` and `CrystalPlane.d_spacing_angstrom` should be documented together, not separately.
- Validation examples should explicitly show the reciprocal-vector magnitude before converting to $d_{hkl}$.

### 3. Direct <-> Reciprocal Component Transforms

Source:
`crystallographY_calcualtions.pdf`, book pp. 16-18 (PDF pp. 25-27)

Using the direct metric tensor $g_{ij}$ and reciprocal metric tensor $g^{*}_{ij}$:

$$
p^{*}_{m} = p_i\, g_{im},
\qquad
p_i = p^{*}_{m}\, g^{*}_{mi}.
$$

Implementation consequence:

- Future vector APIs should expose direct-to-reciprocal and reciprocal-to-direct component transforms as named operations instead of leaving them as ad hoc matrix multiplications in user code.

### 4. Hexagonal Plane Conversion

Source:
`hexagnoal 4index mathematics.pdf`, pp. 1-2

For plane indices:

$$
i = -(h + k),
\qquad
(hkl) \rightarrow (hkil).
$$

Implementation consequence:

- Plane conversions should preserve the supernumerary index in docs and examples, even if PyTex stores the canonical 3-index form internally.

### 5. Hexagonal Direction Conversion

Source:
`hexagnoal 4index mathematics.pdf`, pp. 2-3

From 3-index $[uvw]$ to 4-index $[UVTW]$:

$$
U = \frac{2u - v}{3},
\qquad
V = \frac{2v - u}{3},
\qquad
T = -\frac{u + v}{3},
\qquad
W = w.
$$

Inverse transform:

$$
u = 2U + V,
\qquad
v = 2V + U,
\qquad
w = W.
$$

Constraint:

$$
U + V + T = 0
$$

Implementation consequence:

- Converters should reduce to the smallest integer tuple after applying the formula.
- Test cases should include both divisible-by-3 and non-divisible-by-3 direction inputs.

### 6. Hexagonal Four-Index Zone Law

Source:
`hexagnoal 4index mathematics.pdf`, pp. 2-4

The four-index zone relation can be written as:

$$
hU + kV + iT + lW = 0
$$

with

$$
i = -(h + k),
\qquad
T = -(U + V).
$$

Implementation consequence:

- PyTex currently has `ZoneAxis.contains_miller_index` for 3-index objects.
- A future Miller-Bravais API should expose this four-index zone law directly for hexagonal teaching and validation workflows.

### 7. Hexagonal Interplanar Spacing

Source:
`Kelly & Groves.pdf`, Appendix 3, pp. 469-472

For the hexagonal lattice:

$$
\frac{1}{d_{hkl}^{2}}
=
\frac{4}{3}\,\frac{h^{2} + hk + k^{2}}{a^{2}} + \frac{l^{2}}{c^{2}}
$$

Implementation consequence:

- Hexagonal d-spacing tests should cite this formula explicitly instead of using only numeric parity checks.

### 8. Bragg Law

Source:
`williamsandcarter.pdf`, pp. 78-79

$$
n\lambda = 2 d \sin\theta_{B}
$$

Implementation consequence:

- Powder-XRD docs should distinguish $\theta_{B}$ from the reported $2\theta$.
- Any detector or plotting workflow should state whether its angular axis is $\theta$, $2\theta$, or detector-plane coordinates.

### 9. Bunge Euler Rotation Matrix

Sources:

- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 47-49
- `MathsOfrotations_RolletDegraef.pdf`, p. 21

PyTex docs should keep the Bunge tuple ordered as:

$$
(\varphi_1, \Phi, \varphi_2)
$$

and should present the matrix formula using the same symbol order whenever the docs explain Euler input/output semantics.

Implementation consequence:

- Public documentation and tests should avoid switching symbol names or axis order casually.
- When a page uses the matrix form, it should also state the mapping direction and frame meaning.

### 10. Axis-Angle and Quaternion

Source:
`MathsOfrotations_RolletDegraef.pdf`, pp. 6-7, 23-25

Axis-angle to quaternion:

$$
q = \left(\cos\frac{\omega}{2},\; \hat{\mathbf{n}} \sin\frac{\omega}{2}\right)
$$

Implementation consequence:

- Quaternion docs should say explicitly that PyTex stores this as `(w, x, y, z)`.
- When code canonicalizes the sign of a quaternion, the docs should explain why the equivalent $q$ and $-q$ need a canonical representative.

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

$$
f(g h) = f(g) \quad \text{for crystal symmetry } h,
\qquad
f(s g) = f(g) \quad \text{for specimen symmetry } s.
$$

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
