# Schmid and Taylor Plasticity Metrics

**Surface:** `pytex.properties.slip.SlipSystemFamily.schmid_factors`,
`max_schmid_factor`, `pytex.properties.taylor.taylor_factors`,
`uniaxial_strain_tensor`, with `CrystalMap.schmid_factor_map` and
`taylor_factor_map` providing spatial microstructural field mapping.

In crystal plasticity, crystallographic orientation directly governs anisotropic
yielding and plastic work dissipation. PyTex implements two foundational, complementary
metrics:
1. **The Schmid factor ($m$):** A single-system resolved shear stress criterion under
   stress boundary conditions, identifying the most favorable slip system for initial yield.
2. **The Taylor factor ($M$):** A multi-system optimization metric under strain boundary
   conditions, calculating the minimum total slip activity required to accommodate an
   imposed macroscopic strain tensor.

While frequently compared in microstructural characterization, these metrics reflect
fundamentally distinct mechanical regimes: Schmid analysis assumes unconstrained single-crystal
stress states, whereas Taylor analysis enforces rigorous strain compatibility across a polycrystal.

## 1. The Schmid factor: stress-controlled yield initiation

### 1.1 Mathematical formulation

For a crystal possessing a slip system characterized by slip plane normal $\mathbf{n}$
and slip direction $\mathbf{b}$ (with unit vectors $\lVert\mathbf{n}\rVert = \lVert\mathbf{b}\rVert = 1$
and orthogonality $\mathbf{n} \cdot \mathbf{b} = 0$), subjected to a uniaxial stress
$\sigma$ along specimen unit direction $\mathbf{d}$, the resolved shear stress $\tau$
acting on the slip system is:

$$
\tau = \sigma \, m,
$$

where the **Schmid factor** $m$ is defined as:

$$
m = \cos\phi\,\cos\lambda = (\mathbf{n} \cdot \mathbf{d})\,(\mathbf{b} \cdot \mathbf{d}),
$$

with $\phi$ denoting the angle between the stress axis $\mathbf{d}$ and the slip plane
normal $\mathbf{n}$, and $\lambda$ denoting the angle between $\mathbf{d}$ and the slip
direction $\mathbf{b}$.

### 1.2 Theoretical bounds and kinematic properties

- **Theoretical upper bound ($m \le 0.5$):** Under the orthogonality constraint
  $\mathbf{n} \cdot \mathbf{b} = 0$, the maximum resolved shear stress occurs at
  $\phi = \lambda = 45^\circ$, where $m = \cos(45^\circ)\cos(45^\circ) = 0.5$.
  A calculated Schmid factor exceeding $0.5$ indicates non-orthogonal vectors or unnormalized inputs.
- **Absolute magnitude ($|m|$):** Dislocation glide operates along either forward or
  reverse slip directions ($\pm\mathbf{b}$) under reversed shear stress. Yield initiation
  depends on $|\tau| \ge \tau_{\text{CRSS}}$, where $\tau_{\text{CRSS}}$ is the critical
  resolved shear stress. PyTex evaluates $|m|$ to assess yield susceptibility irrespective
  of sign convention.

### 1.3 Vectorized evaluation algorithm

```text
Input : Slip system family (unit normals N, directions B in crystal frame),
        Crystal orientations R (crystal-to-specimen), Stress axis d

1  Normalize stress axis: d <- d / ||d||
2  Transform slip systems to specimen frame:
       n_spec = R n,   b_spec = R b
3  Compute direction cosines:
       cos_phi    = n_spec . d
       cos_lambda = b_spec . d
4  Evaluate Schmid factors:
       m = |cos_phi * cos_lambda|
5  Maximum Schmid factor:
       m_max = max_{s in family} m_s
```

All tensor contractions are evaluated via batch NumPy `einsum` operations across
arbitrary spatial pixel dimensions without Python loops.

### 1.4 Physical assumptions and limitations

1. **Uniaxial stress idealization:** The classical Schmid factor assumes a purely
   uniaxial stress tensor $\boldsymbol{\sigma} = \sigma(\mathbf{d}\otimes\mathbf{d})$.
   Under multiaxial stress states, the resolved shear stress requires the full tensor
   contraction $\tau_s = \boldsymbol{\sigma} : \mathbf{P}_s$, where
   $\mathbf{P}_s = \frac{1}{2}(\mathbf{n}_s \otimes \mathbf{b}_s + \mathbf{b}_s \otimes \mathbf{n}_s)$
   is the symmetric Schmid tensor.
2. **Asymmetric critical resolved shear stresses:** In low-symmetry materials (e.g.,
   HCP titanium, zirconium, magnesium), distinct slip families (basal, prismatic,
   first-order pyramidal, second-order pyramidal) exhibit CRSS values that differ by
   factors of 2 to 10. Ranking slip systems purely by geometric Schmid factor without
   accounting for CRSS ratios can misidentify the active yield mode. PyTex structures
   `SlipSystemFamily` instances independently, allowing CRSS-weighted stress evaluations.
3. **Intergranular constraint absence:** Individual grains embedded in a polycrystal
   experience internal constraint stresses imposed by surrounding grains, modifying local
   stress states away from the macroscopic axis.

## 2. The Taylor factor: strain-controlled polycrystal compatibility

### 2.1 Mechanical foundation and the von Mises criterion

In the full-constraint Taylor model (Taylor, 1938), every crystallite within a polycrystal
is assumed to undergo an identical macroscopic plastic strain tensor $\mathbf{E}$.
Plastic volume conservation requires isochoric deformation:

$$
\operatorname{tr}(\mathbf{E}) = E_{11} + E_{22} + E_{33} = 0.
$$

A symmetric, deviatoric strain tensor in $\mathbb{R}^{3 \times 3}$ possesses five
independent degrees of freedom. By the von Mises criterion (von Mises, 1928), a crystal
must activate at least five linearly independent slip systems to accommodate an
arbitrary isochoric strain state without intergranular void formation or overlap.

### 2.2 Linear programming formulation

The total macroscopic plastic strain is accommodated as a linear combination of shear
increments $\gamma_s$ across available slip systems $s = 1, \dots, S$:

$$
\mathbf{E} = \sum_{s=1}^S \gamma_s\,\mathbf{P}_s.
$$

Because glide can occur in positive or negative directions along $\mathbf{b}_s$, each
shear increment is decomposed into non-negative components:
$\gamma_s = \gamma_s^+ - \gamma_s^-$, with $\gamma_s^\pm \ge 0$.

The Taylor factor minimizes the total internal plastic dissipation per unit equivalent strain:

$$
\min_{\boldsymbol{\gamma}^+, \boldsymbol{\gamma}^-} \; \sum_{s=1}^S \left(\gamma_s^+ + \gamma_s^-\right),
$$

subject to the 5 independent strain accommodation equality constraints:

$$
\sum_{s=1}^S \left(\gamma_s^+ - \gamma_s^-\right)\mathbf{P}_s = \mathbf{E},
\qquad \gamma_s^+ \ge 0, \quad \gamma_s^- \ge 0.
$$

PyTex solves this constrained linear program using `scipy.optimize.linprog` with the
HiGHS interior-point and simplex solvers. The **Taylor factor** $M$ is then evaluated as:

$$
M = \frac{\sum_{s=1}^S |\gamma_s|}{\varepsilon_{\text{eq}}},
$$

where $\varepsilon_{\text{eq}} = \sqrt{\frac{2}{3}\mathbf{E}:\mathbf{E}}$ is the von Mises
equivalent plastic strain. For uniaxial tension along direction 1 with transverse isotropic
contraction ($E_{11} = \varepsilon, E_{22} = E_{33} = -\frac{1}{2}\varepsilon$),
$\varepsilon_{\text{eq}} = \varepsilon$.

### 2.3 Physical interpretation of solver infeasibility

When the set of available slip systems cannot span the 5-dimensional deviatoric strain
subspace for a given crystal orientation, no non-negative combination of shears can
satisfy the strain compatibility equations. In this scenario, the linear program is
primal infeasible.

PyTex returns $M = \infty$ for infeasible orientations. This is a direct physical result:
- In HCP metals, **basal slip $\{0001\}\langle 11\bar{2}0\rangle$ provides only 2 linearly
  independent slip systems**, and cannot accommodate strain along the crystal $\mathbf{c}$-axis.
  A basal-only slip family subjected to uniaxial extension along $[0001]$ is strictly
  incapable of accommodating the strain, correctly returning $M = \infty$.
- Polycrystalline deformation requires secondary non-basal mechanisms (prismatic slip,
  pyramidal slip, or deformation twinning) to fulfill the von Mises five-system requirement.

### 2.4 Model hierarchy and upper-bound characteristics

The full-constraint Taylor model represents a strict upper bound on polycrystalline
flow stress (Bishop & Hill, 1951), as it enforces rigid strain compatibility while
relaxing stress equilibrium across grain boundaries. Relaxed-constraint models
(e.g., pancake models for rolled sheet) and self-consistent viscoplastic formulations
(VPSC) predict intermediate flow stresses between the lower-bound Sachs (isostress)
limit and the Taylor (isostrain) limit.

## 3. Spatial microstructural mapping

PyTex applies both formulations across spatially resolved EBSD datasets via `CrystalMap`:

| Feature | `CrystalMap.schmid_factor_map` | `CrystalMap.taylor_factor_map` |
| --- | --- | --- |
| Controlling boundary condition | Uniaxial stress along specimen vector $\mathbf{d}$ | Deviatoric strain tensor $\mathbf{E}$ |
| Computational cost | $O(N_{\text{pix}} \cdot S)$ direct tensor inner products | $O(N_{\text{pix}})$ linear program solves |
| Performance scaling | Fully vectorized NumPy array operations | Parallelized batch linear programming via HiGHS |
| Mechanical interpretation | Initial yield susceptibility of surface grains | Hardness and plastic work dissipation in constrained bulk |

## 4. Method selection guide

| Engineering / Scientific Inquiry | Appropriate Metric | Key Considerations |
| --- | --- | --- |
| Slip trace identification on deformed surface | **Schmid factor** | Surface grains experience relaxed transverse constraints |
| Soft versus hard grain identification in tensile test | **Schmid factor** | Apply CRSS weighting when evaluating multi-family HCP alloys |
| Polycrystalline aggregate yield stress estimation | **Taylor factor** | Volume-average $M$ over the experimental ODF ($\bar{M} = \int M(g) f(g)\,\mathrm{d}g$) |
| Identification of deformation twinning requirements | **Taylor factor** | Orientations where slip-only Taylor factors diverge ($M \to \infty$) necessitate twinning |

## Verification

- `tests/unit/test_schmid.py` and `tests/unit/test_taylor.py`: Verify the theoretical bound
  $m \le 0.5$ at $\phi = \lambda = 45^\circ$, analytical Taylor factors for ideal FCC rolling
  orientations (Copper, Brass, Goss, Cube), and von Mises rank-deficiency detection.
- Executable worked examples:
  - {doc}`../examples/generated/schmid-and-taylor`

## See also

- {doc}`../theory/schmid_and_taylor_plasticity` — Detailed derivations of the von Mises criterion and Bishop–Hill stress polyhedra.
- {doc}`ipf_coloring` — Visualization of spatial orientation fields alongside plasticity maps.
- {doc}`ebsd_grains_and_local_misorientation` — Grain segmentation and misorientation gradient analysis.

## References

### Normative

- Schmid, E. & Boas, W. (1935). *Kristallplastizität*. Springer. <https://doi.org/10.1007/978-3-662-34532-0>
- Taylor, G. I. (1938). Plastic strain in metals. *Journal of the Institute of Metals* **62**, 307–324.
- von Mises, R. (1928). Mechanik der plastischen Formänderung von Kristallen. *Zeitschrift für Angewandte Mathematik und Mechanik* **8**, 161–185. <https://doi.org/10.1002/zamm.19280080302>

### Informative

- Bishop, J. F. W. & Hill, R. (1951). A theory of the plastic distortion of a polycrystalline aggregate under combined stresses. *Philosophical Magazine* **42**, 414–427. <https://doi.org/10.1080/14786445108561065>
- Kocks, U. F., Tomé, C. N. & Wenk, H.-R. (1998). *Texture and Anisotropy: Preferred Orientations in Polycrystals and their Effect on Materials Properties*. Cambridge University Press.
