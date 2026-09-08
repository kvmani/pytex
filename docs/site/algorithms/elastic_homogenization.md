# Elastic Homogenization and Directional Moduli

**Surface:** `pytex.properties.tensors.StiffnessTensor`, `ComplianceTensor`,
`homogenize_elastic`, `youngs_modulus_surface`,
`linear_compressibility_surface`, `shear_modulus_surface`,
`poisson_ratio_surface`, `DirectionalModulusSurface`.

Crystalline solids exhibit direction-dependent elastic stiffness governed by atomic
bonding anisotropy. In polycrystalline aggregates, the macroscopic elastic response
represents an orientation-weighted average of single-crystal stiffnesses modified
by crystallographic texture. Calculating the effective macroscopic elastic tensor
from single-crystal elasticity and an orientation distribution function (ODF) constitutes
the classical problem of **elastic homogenization**.

Because local stress and strain fields within polycrystalline aggregates vary across
grain boundaries, an orientation distribution alone does not uniquely specify the
internal field equilibrium. Consequently, elastic homogenization yields rigorous
upper and lower bounds rather than an isolated scalar prediction. This page presents
the fourth-rank tensor formulation, explains the conversion factors between tensor
and Voigt representations, details the Voigt, Reuss, and Hill homogenization bounds,
and formulates directional modulus surfaces on the unit sphere.

## 1. Tensor elasticity and matrix representations

### 1.1 Fourth-rank constitutive relations

Linear elasticity is governed by generalized Hooke's law relating the second-rank
Cauchy stress tensor $\sigma_{ij}$ to the second-rank infinitesimal strain tensor
$\varepsilon_{kl}$:

$$
\sigma_{ij} = C_{ijkl}\,\varepsilon_{kl}, \qquad \varepsilon_{ij} = S_{ijkl}\,\sigma_{kl},
$$

where $C_{ijkl}$ is the fourth-rank **elastic stiffness tensor** and $S_{ijkl}$ is
the fourth-rank **elastic compliance tensor**. Both tensors satisfy major and minor
thermodynamic symmetries:

$$
C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij}.
$$

The stiffness and compliance tensors are mutual inverses in the space of fourth-rank
symmetric tensors:

$$
C_{ijmn}\,S_{mnkl} = I^{\text{sym}}_{ijkl} = \frac{1}{2}\left(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}\right).
$$

### 1.2 The Voigt notation convention

In engineering literature, fourth-rank tensors are conventionally contracted into
$6 \times 6$ symmetric matrices using the Voigt index mapping:

$$
11 \to 1, \quad 22 \to 2, \quad 33 \to 3, \quad 23,32 \to 4, \quad 13,31 \to 5, \quad 12,21 \to 6.
$$

While stiffness components map directly between tensor and matrix forms
($C_{ijkl} \leftrightarrow C_{\alpha\beta}$), compliance components require metric factors
of 2 and 4 due to the distinction between tensor shear strain $\varepsilon_{ij}$ ($i \ne j$)
and engineering shear strain $\gamma_{\alpha} = 2\varepsilon_{ij}$:

$$
S_{\alpha\beta} = \begin{cases}
S_{ijkl}, & \alpha \le 3 \text{ and } \beta \le 3, \\
2\,S_{ijkl}, & \alpha \le 3, \beta > 3 \text{ or } \alpha > 3, \beta \le 3, \\
4\,S_{ijkl}, & \alpha > 3 \text{ and } \beta > 3.
\end{cases}
$$

While matrix inversion satisfies $[C_{\alpha\beta}]^{-1} = [S_{\alpha\beta}]$, coordinate
frame rotations cannot be performed on $6 \times 6$ Voigt matrices using standard orthogonal
transformation rules without introducing specialized Bond transformation matrices.

To prevent silent conversion errors, PyTex maintains the **fourth-rank Cartesian tensor**
$C_{ijkl}$ and $S_{ijkl}$ as the canonical representation. Coordinate rotations are
evaluated directly via fourth-rank tensor transformation:

$$
C'_{ijkl} = R_{ip}\,R_{jq}\,R_{kr}\,R_{ls}\,C_{pqrs}.
$$

## 2. Variational homogenization bounds

Because exact stress and strain distributions depend on intergranular boundary topology
and grain morphology, variational energy principles establish rigorous bounds on
effective polycrystalline elasticity.

### 2.1 The Voigt bound (uniform strain)

The Voigt model (Voigt, 1928) assumes an isostrain state: every grain experiences an
identical strain tensor equal to the macroscopic aggregate strain ($\boldsymbol{\varepsilon}(g) \equiv \bar{\boldsymbol{\varepsilon}}$).
Averaging the resulting local stresses over all crystal orientations yields the **Voigt effective stiffness**:

$$
\mathbf{C}^{\text{V}} = \langle \mathbf{C}(g) \rangle = \int_{\mathrm{SO}(3)} \mathbf{C}(g)\,f(g)\,\mathrm{d}g.
$$

By the principle of minimum potential energy, the assumption of uniform strain overconstrains
the internal degrees of freedom, rendering $\mathbf{C}^{\text{V}}$ a rigorous **upper bound**
on aggregate stiffness.

### 2.2 The Reuss bound (uniform stress)

The Reuss model (Reuss, 1929) assumes an isostress state: every crystallite carries an
identical stress tensor equal to the macroscopic aggregate stress ($\boldsymbol{\sigma}(g) \equiv \bar{\boldsymbol{\sigma}}$).
Averaging the local elastic strains over orientation space yields the **Reuss effective compliance**:

$$
\mathbf{S}^{\text{R}} = \langle \mathbf{S}(g) \rangle = \int_{\mathrm{SO}(3)} \mathbf{S}(g)\,f(g)\,\mathrm{d}g.
$$

The **Reuss effective stiffness** is obtained by tensor inversion of $\mathbf{S}^{\text{R}}$:

$$
\mathbf{C}^{\text{R}} = \left(\mathbf{S}^{\text{R}}\right)^{-1}.
$$

By the principle of minimum complementary energy, the assumption of uniform stress relaxes
intergranular displacement compatibility, rendering $\mathbf{C}^{\text{R}}$ a rigorous **lower bound**
on aggregate stiffness.

### 2.3 The Voigt–Reuss–Hill (VRH) approximation

Hill (1952) demonstrated that the true effective stiffness of a random polycrystalline aggregate
is bounded between the Voigt and Reuss limits:

$$
\mathbf{C}^{\text{R}} \le \mathbf{C}^{\text{effective}} \le \mathbf{C}^{\text{V}}.
$$

The **Voigt–Reuss–Hill (VRH) average** is defined as the arithmetic mean of the two bounds:

$$
\mathbf{C}^{\text{VRH}} = \frac{1}{2}\left(\mathbf{C}^{\text{V}} + \mathbf{C}^{\text{R}}\right).
$$

The VRH average is widely adopted as an engineering approximation that provides close
agreement with experimental measurements. The magnitude of the bound separation
$\lVert \mathbf{C}^{\text{V}} - \mathbf{C}^{\text{R}} \rVert$ quantifies the degree of
elastic anisotropy within the aggregate.

## 3. Homogenization algorithm

```text
Input : Single-crystal stiffness C (fourth-rank tensor),
        Discrete orientations R_n, Normalized weights w_n, Scheme

1  Transform single-crystal stiffness into specimen frame for each orientation:
       C'_n = einsum('nip,njq,nkr,nls,pqrs->nijkl', R, R, R, R, C)
2  Voigt stiffness average:
       C_V = sum_n w_n C'_n
3  If Scheme == "voigt": Return C_V
4  Compute single-crystal compliance S = inverse(C)
5  Transform compliance into specimen frame:
       S'_n = einsum('nip,njq,nkr,nls,pqrs->nijkl', R, R, R, R, S)
6  Reuss compliance average:
       S_mean = sum_n w_n S'_n
7  Reuss stiffness:
       C_R = inverse(S_mean)
8  If Scheme == "reuss": Return C_R
9  Hill average:
       C_VRH = 0.5 * (C_V + C_R)
   Return C_VRH
```

By expressing orientation rotations as batch Einstein summations over all grains
simultaneously, the algorithm eliminates per-orientation Python loops and executes
at vectorized linear algebra speeds.

## 4. Directional modulus surfaces

For any homogenized aggregate or single-crystal compliance tensor $\mathbf{S}$, directional
elastic moduli are evaluated along arbitrary unit vectors on the sphere $\mathbb{S}^2$:

| Property | Analytical Formulation | Angular Parameters |
| --- | --- | --- |
| **Young's modulus** $E(\mathbf{d})$ | $E(\mathbf{d}) = \frac{1}{S'_{1111}(\mathbf{d})} = \frac{1}{d_i d_j d_k d_l S_{ijkl}}$ | Single direction $\mathbf{d} \in \mathbb{S}^2$ |
| **Linear compressibility** $\beta(\mathbf{d})$ | $\beta(\mathbf{d}) = d_i d_j S_{ijkk}$ | Single direction $\mathbf{d} \in \mathbb{S}^2$ |
| **Shear modulus** $G(\mathbf{n}, \mathbf{m})$ | $G(\mathbf{n}, \mathbf{m}) = \frac{1}{4\,n_i m_j n_k m_l S_{ijkl}}$ | Normal $\mathbf{n}$, shear direction $\mathbf{m} \perp \mathbf{n}$ |
| **Poisson's ratio** $\nu(\mathbf{n}, \mathbf{m})$ | $\nu(\mathbf{n}, \mathbf{m}) = -\frac{n_i n_j m_k m_l S_{ijkl}}{n_p n_q n_r n_s S_{pqrs}}$ | Axial direction $\mathbf{n}$, transverse direction $\mathbf{m} \perp \mathbf{n}$ |

Because shear modulus and Poisson's ratio depend on a plane normal $\mathbf{n}$ and an
in-plane shear direction $\mathbf{m}$, `shear_modulus_surface` and `poisson_ratio_surface`
evaluate directional extrema ($\min_{\mathbf{m}} G$, $\max_{\mathbf{m}} G$, $\min_{\mathbf{m}} \nu$,
$\max_{\mathbf{m}} \nu$) across the orthogonal circle $\mathbf{m} \cdot \mathbf{n} = 0$.

> [!NOTE]
> In certain anisotropic cubic and low-symmetry crystals, Poisson's ratio can become
> negative along specific crystallographic axes. PyTex preserves negative Poisson's
> ratios (auxetic behavior) without artificial clamping.

## 5. Assumptions and limitations

| Domain Feature | Included in Formulation | Model Limitation |
| --- | --- | --- |
| Crystallographic texture | Yes (ODF-weighted tensor averaging) | Requires representative experimental orientation sampling |
| Single-crystal elastic anisotropy | Yes (Full 21-parameter anisotropic tensor) | Assumes linear elastic response without dislocation plasticity |
| Morphological texture | No | Does not account for non-equiaxed grain aspect ratios or alignment |
| Multi-phase composites | Single phase | Multi-phase composites require volume-fraction Mori–Tanaka or self-consistent extensions |
| Microstructural stress concentrations | No | Local boundary traction concentrations are not resolved |

## Verification

- `tests/unit/test_elastic.py`: Verifies Voigt-Reuss-Hill ordering ($\mathbf{C}^{\text{R}} \le \mathbf{C}^{\text{VRH}} \le \mathbf{C}^{\text{V}}$),
  exact isotropy recovery for uniform random ODFs, and fourth-rank tensor inversion accuracy.
- Executable worked examples:
  - {doc}`../examples/generated/elastic-anisotropy`

## See also

- {doc}`../theory/elastic_anisotropy_and_homogenization` — Derivation of variational energy principles and bounding theorems.
- {doc}`pole_figure_inversion` — Experimental ODF reconstruction providing homogenization weights.
- {doc}`schmid_and_taylor` — Plasticity analogues to elastic homogenization.

## References

### Normative

- Voigt, W. (1928). *Lehrbuch der Kristallphysik*. Teubner.
- Reuss, A. (1929). Berechnung der Fließgrenze von Mischkristallen. *Zeitschrift für Angewandte Mathematik und Mechanik* **9**, 49–58. <https://doi.org/10.1002/zamm.19290090104>
- Hill, R. (1952). The elastic behaviour of a crystalline aggregate. *Proceedings of the Physical Society A* **65**, 349–354. <https://doi.org/10.1088/0370-1298/65/5/307>

### Informative

- Hashin, Z. & Shtrikman, S. (1962). A variational approach to the theory of the elastic behaviour of polycrystals. *Journal of the Mechanics and Physics of Solids* **10**, 343–352. <https://doi.org/10.1016/0022-5096(62)90005-4>
- Nye, J. F. (1985). *Physical Properties of Crystals: Their Representation by Tensors and Matrices*. Oxford University Press.
- Kocks, U. F., Tomé, C. N. & Wenk, H.-R. (1998). *Texture and Anisotropy: Preferred Orientations in Polycrystals and their Effect on Materials Properties*. Cambridge University Press.
