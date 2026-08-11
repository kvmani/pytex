# Elastic Anisotropy And Polycrystal Homogenization

A single crystal is elastically anisotropic; a randomly textured aggregate of the same crystals is
not. Getting from one to the other involves a notation that hides two factors of two, a closed form
that most texts quote without deriving, and a pair of bounds whose gap for cubic symmetry lives
entirely in one modulus. This note fixes all three for `pytex.properties.tensors`.

## Voigt Notation, And The Factors That Break It

The elastic tensors are rank four with the symmetries
$C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij}$, leaving 21 independent constants. Voigt notation
compresses the symmetric index pairs onto $1\ldots6$ so the tensor becomes a $6\times6$ matrix.

The compression is **not** the same for stiffness and compliance, and this is the single most
common source of wrong elastic constants:

$$
C_{mn} = C_{ijkl},
\qquad
S_{mn} = f_{m} f_{n} S_{ijkl},
\qquad
f_{m} = \begin{cases} 1 & m \le 3 \\ 2 & m > 3 \end{cases}
$$ (eq-el-voigt)

Stiffness carries no factors; compliance carries a factor 2 on each shear index, hence 1, 2 or 4
depending on how many of the pair are shear. The factors exist so that the matrix inverse and the
tensor inverse agree — so that $\mathbf{S} = \mathbf{C}^{-1}$ in $6\times6$ form is the same object
as the tensor that contracts correctly in

$$
\frac{1}{E(\hat{\mathbf{n}})} = n_{i} n_{j} n_{k} n_{l} \, S_{ijkl} .
$$ (eq-el-young)

Dropping the factors leaves a matrix that still inverts, still looks symmetric, and gives shear
moduli wrong by 2 or 4. PyTex therefore stores the full rank-4 tensor as the primary
representation and treats the Voigt matrix as a view, which is why `StiffnessTensor` and
`ComplianceTensor` are separate types rather than one type with a flag: the factor rule is a
property of *which* tensor is being compressed, not of the compression.

Verified for copper ($C_{11}, C_{12}, C_{44} = 168.4,\ 121.4,\ 75.4$ GPa) against the closed-form
cubic inverse

$$
S_{11} = \frac{C_{11}+C_{12}}{(C_{11}-C_{12})(C_{11}+2C_{12})},
\quad
S_{12} = \frac{-C_{12}}{(C_{11}-C_{12})(C_{11}+2C_{12})},
\quad
S_{44} = \frac{1}{C_{44}},
$$ (eq-el-cubic-inverse)

which the implementation reproduces to $1.7\times10^{-18}$.

## Directional Young's Modulus Has A Closed Form

For cubic symmetry {eq}`eq-el-young` collapses. Only one anisotropy combination survives, and the
direction enters through a single scalar:

$$
\frac{1}{E(\hat{\mathbf{n}})}
= S_{11} - 2\left(S_{11} - S_{12} - \tfrac{1}{2}S_{44}\right) J(\hat{\mathbf{n}}),
\qquad
J = n_{1}^{2}n_{2}^{2} + n_{2}^{2}n_{3}^{2} + n_{3}^{2}n_{1}^{2} .
$$ (eq-el-cubic-young)

$J$ runs from 0 along $\langle 100 \rangle$ to $1/3$ along $\langle 111 \rangle$, so those two
directions are the extremes and everything else interpolates. For copper:

| Direction | $J$ | $E$ closed form (GPa) | $E$ computed (GPa) | Difference |
| --- | ---: | ---: | ---: | ---: |
| $[100]$ | $0$ | 66.6888 | 66.6888 | $0$ |
| $[110]$ | $1/4$ | 130.3376 | 130.3376 | $2.8\times10^{-14}$ |
| $[112]$ | $1/4$ | 130.3376 | 130.3376 | $1.7\times10^{-13}$ |
| $[111]$ | $1/3$ | 191.1497 | 191.1497 | $4.0\times10^{-13}$ |

**$[110]$ and $[112]$ have exactly the same stiffness**, which is not obvious from the indices and
is not a coincidence of copper: $J$ is $1/4$ for both, so the equality holds for every cubic
material. Directions of equal $J$ form contours on the sphere, and only $J$ — not the direction —
enters {eq}`eq-el-cubic-young`.

The degree of anisotropy is the Zener ratio

$$
A_{Z} = \frac{2C_{44}}{C_{11}-C_{12}},
$$ (eq-el-zener)

equal to 1 exactly when the bracket in {eq}`eq-el-cubic-young` vanishes, i.e. when
$S_{44} = 2(S_{11}-S_{12})$, in which case $E$ is direction-independent. Copper has
$A_{Z} = 3.21$ and $E_{[111]}/E_{[100]} = 2.87$ — nearly a factor of three between the soft and
stiff directions of the same crystal.

## Rotating The Tensor

Expressing the crystal tensor in the specimen frame is the rank-4 transformation

$$
C'_{ijkl} = R_{ip}R_{jq}R_{kr}R_{ls}\,C_{pqrs},
$$ (eq-el-rotate)

evaluated for a whole orientation population as one `einsum` contraction rather than a loop. This
is the step where a Voigt-matrix representation would have to be converted back to a tensor
anyway, which is the second reason for storing rank four.

## Voigt And Reuss Are Bounds, Hill Is Not

An aggregate's response depends on how stress and strain partition between grains, which a
texture alone does not determine. Two assumptions give the extremes:

- **Voigt**: uniform *strain* in every grain. Compatible but not equilibrated, so it is stiff — an
  upper bound. It averages stiffnesses.
- **Reuss**: uniform *stress* in every grain. Equilibrated but not compatible, so it is compliant —
  a lower bound. It averages compliances and inverts.

$$
\mathbf{C}_{V} = \bigl\langle \mathbf{C}' \bigr\rangle,
\qquad
\mathbf{C}_{R} = \bigl\langle \mathbf{S}' \bigr\rangle^{-1},
\qquad
\mathbf{C}_{H} = \tfrac{1}{2}\left(\mathbf{C}_{V} + \mathbf{C}_{R}\right).
$$ (eq-el-vrh)

The bounding property is variational — each follows from a minimum-energy principle — and the
ordering $\mathbf{C}_{R} \preceq \mathbf{C} \preceq \mathbf{C}_{V}$ is rigorous for any
microstructure with the given texture.

**The Hill average has no such status.** It is the arithmetic mean of an upper and a lower bound,
chosen because it is usually closer to measurement than either, and it is not itself a bound and
carries no variational justification. It is a reasonable default and should not be reported as
though it were derived. Where a genuine improvement is needed, the Hashin–Shtrikman bounds are
tighter and remain rigorous.

## The Cubic Bulk Modulus Is Not Bounded — It Is Exact

For a randomly textured cubic aggregate the closed forms are

$$
K_{V} = K_{R} = \frac{C_{11}+2C_{12}}{3},
\qquad
\mu_{V} = \frac{C_{11}-C_{12}+3C_{44}}{5},
\qquad
\mu_{R} = \frac{5}{4(S_{11}-S_{12}) + 3S_{44}} .
$$ (eq-el-cubic-vrh)

The first equality is the non-obvious one and it is exact, not approximate: **for cubic symmetry
the Voigt and Reuss bulk moduli coincide**, verified here to $8.5\times10^{-14}$. A cubic crystal's
response to hydrostatic pressure is isotropic — pressure produces the same dilatation whatever the
orientation — so uniform-stress and uniform-strain assumptions cannot disagree about it. The whole
Voigt–Reuss gap therefore lives in the shear modulus, and quoting a "Hill bulk modulus" for a cubic
aggregate suggests an uncertainty that does not exist.

For copper the gap is large: $\mu_{V} = 54.64$, $\mu_{R} = 40.03$ GPa, a spread of 14.6 GPa or
about 31% of the Hill value $47.34$ GPa. That is the honest uncertainty on a randomly textured
copper shear modulus from texture information alone.

### Agreement with the numerical route

Homogenizing over $4\times10^{4}$ Haar-random orientations, against {eq}`eq-el-cubic-vrh`:

| Scheme | $K$ numeric | $(C_{11}-C_{12})/2$ | $C_{44}$ | $\mu$ closed form |
| --- | ---: | ---: | ---: | ---: |
| Voigt | 137.02 | 54.62 | 54.63 | 54.64 |
| Reuss | 137.03 | 40.02 | 40.03 | 40.03 |
| Hill | 137.03 | 47.32 | 47.33 | 47.34 |

against the exact $K = 137.07$ GPa. Two checks are running at once. The columns
$(C_{11}-C_{12})/2$ and $C_{44}$ are independent measures of the *same* shear modulus for an
isotropic material, so their agreement to $0.01$ GPa shows the numerical aggregate really has
become isotropic; and both track the closed form to about $0.03\%$, the residual being finite-sample
texture in $4\times10^{4}$ orientations rather than an error in the averaging.

## Assumptions And Limits

- The bounds assume the texture is the only information. They know nothing about grain shape,
  spatial correlation, or connectivity, so a strongly banded microstructure can sit near a bound
  while an equiaxed one sits near Hill.
- Homogenization is orientation-weighted and treats grains as equal unless weights are supplied;
  volume weighting is the caller's responsibility.
- The single-crystal constants are inputs. Their temperature dependence is usually a larger error
  than the choice of averaging scheme, and none of the schemes can compensate for a wrong
  $C_{ijkl}$.
- Only linear elasticity is modelled: no pressure dependence, no anelasticity, no damage.

## References

### Normative

- IEEE/ANSI and IUCr conventions for elastic-constant reduction to Voigt form; see also
  *International Tables for Crystallography, Volume D: Physical Properties of Crystals*, IUCr.
  DOI: <https://doi.org/10.1107/97809553602060000105>. Fixes the tensor symmetries and the
  restrictions each crystal class imposes on $C_{ijkl}$.

### Informative

- J. F. Nye, *Physical Properties of Crystals*, Oxford University Press. The reference treatment of
  {eq}`eq-el-voigt` and the factor convention.
- R. Hill, *The elastic behaviour of a crystalline aggregate*, Proceedings of the Physical Society
  A **65** (1952) 349–354. DOI: <https://doi.org/10.1088/0370-1298/65/5/307>. The averaging scheme
  and, explicitly, its lack of variational standing.
- Z. Hashin and S. Shtrikman, *A variational approach to the theory of the elastic behaviour of
  polycrystals*, Journal of the Mechanics and Physics of Solids **10** (1962) 343–352.
  DOI: <https://doi.org/10.1016/0022-5096(62)90005-4>. The tighter rigorous bounds.
- G. Simmons and H. Wang, *Single Crystal Elastic Constants and Calculated Aggregate Properties*,
  MIT Press (1971). Source of the copper constants used above.

## See Also

- {doc}`orientation_space_and_disorientation` — the orientation populations being averaged over.
- {doc}`discrete_odf_and_pole_figures` — where a weighted orientation population comes from.
