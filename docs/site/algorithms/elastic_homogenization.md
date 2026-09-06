# Elastic Homogenisation And Directional Moduli

**Surface:** `pytex.properties.tensors.StiffnessTensor`, `ComplianceTensor`,
`homogenize_elastic`, `youngs_modulus_surface`,
`linear_compressibility_surface`, `shear_modulus_surface`,
`poisson_ratio_surface`, `DirectionalModulusSurface`.

A single crystal is elastically anisotropic; a textured polycrystal inherits a
weaker version of that anisotropy, and an untextured one is isotropic. Going
from the single-crystal stiffness plus an orientation distribution to the
aggregate response is **homogenisation**, and the answer is not a single number
but a bracket, because the exact result depends on grain-scale stress and strain
fields that an orientation distribution does not contain.

## 1. The tensors, and the Voigt trap

The stiffness $C_{ijkl}$ relates stress to strain, $\sigma_{ij} = C_{ijkl}
\varepsilon_{kl}$; the compliance $S_{ijkl}$ is its inverse *as a fourth-rank
tensor*.

Both are conventionally written as $6\times 6$ Voigt matrices, and this is the
single most common source of silent error in elasticity code, because
**stiffness and compliance do not use the same Voigt convention**:

$$
C_{ijkl} \leftrightarrow C_{mn} \text{ directly,}
\qquad
S_{ijkl} \leftrightarrow S_{mn} \text{ with factors of } 1, 2, 4 .
$$

The factors arise because engineering shear strain is twice the tensor shear
strain. Consequences that follow, and that PyTex handles by never letting the
$6\times6$ form be the source of truth:

- $[C_{mn}]^{-1} = [S_{mn}]$ **as matrices** — that identity does hold — but
  $S_{ijkl}$ recovered from $S_{mn}$ without the factors is wrong.
- Rotating a $6\times6$ matrix with a $6\times6$ Bond matrix is a different
  operation for stiffness and compliance.

`ElasticTensor` stores the **fourth-rank tensor**, and the Voigt matrix is a
view produced with the right factors on the way in and out. Rotation is then the
unambiguous tensor operation

$$
C'_{ijkl} = R_{ip}R_{jq}R_{kr}R_{ls}\,C_{pqrs},
$$

one `einsum` over all orientations at once.

## 2. The bounds: why there are two answers

An aggregate's true stiffness depends on how stress and strain distribute among
grains, which the ODF does not tell us. Two extreme assumptions give two exact
bounds.

### 2.1 Voigt — uniform strain

Assume every grain suffers the **same strain** as the aggregate. Then stresses
are averaged and

$$
\mathbf{C}^{\text{V}} = \bigl\langle \mathbf{C}(g) \bigr\rangle .
$$

Compatibility is satisfied everywhere (all grains deform alike) but equilibrium
is violated at grain boundaries, where tractions do not match. This
over-constrains, so **Voigt is an upper bound**.

### 2.2 Reuss — uniform stress

Assume every grain carries the **same stress**. Then compliances are averaged
and

$$
\mathbf{C}^{\text{R}} = \bigl\langle \mathbf{S}(g) \bigr\rangle^{-1}.
$$

Equilibrium is satisfied, compatibility is violated — grains would separate or
interpenetrate. This under-constrains, so **Reuss is a lower bound**.

Note the asymmetry the code respects: Reuss is the **inverse of the mean
compliance**, *not* the mean of the stiffnesses' inverses in any other order.
Averaging stiffness and inverting gives Voigt; inverting and averaging gives
Reuss, and they differ.

### 2.3 Hill — the average of the two

$$
\mathbf{C}^{\text{VRH}} = \tfrac{1}{2}\left(\mathbf{C}^{\text{V}} + \mathbf{C}^{\text{R}}\right)
$$

The Voigt-Reuss-Hill average is the default because it is usually closer to
measurement than either bound. It is worth being clear about what it is: an
**empirical midpoint, with no variational status**. Voigt and Reuss are rigorous
bounds; Hill is a useful convention. The gap between the bounds is the honest
statement of what the ODF alone can determine, and a narrow gap means the
aggregate is nearly isotropic, not that the model is precise.

Tighter bounds exist — Hashin-Shtrikman uses two-point statistics — and a
self-consistent scheme solves for a consistent effective medium. Neither is
implemented here, and the page says so rather than implying the bracket is the
last word.

## 3. The algorithm

```text
input : single-crystal stiffness C, orientations R_n, weights w_n, scheme

1  normalise the weights
2  rotate the stiffness into the sample frame for every orientation, at once:
       C'_n = einsum('nip,njq,nkr,nls,pqrs->nijkl', R, R, R, R, C)
3  Voigt  : C_V = sum_n w_n C'_n
4  if scheme is "voigt": return C_V
5  rotate the compliance the same way; S_mean = sum_n w_n S'_n
6  Reuss  : C_R = inverse(S_mean)          -- as a fourth-rank tensor
7  if scheme is "reuss": return C_R
8  Hill   : (C_V + C_R) / 2
```

The four-matrix `einsum` in step 2 is the whole cost, and it is done once over
all orientations rather than per grain. The weights come from the ODF, so a
homogenisation is only as good as the texture measurement behind it — see
{doc}`pole_figure_inversion`.

## 4. Directional surfaces

Once an aggregate (or single-crystal) stiffness is in hand, the directional
properties are sampled on a spherical grid:

| Function | Quantity | Depends on |
| --- | --- | --- |
| `youngs_modulus_surface` | $E(\mathbf{d}) = 1/S'_{1111}$ | one direction |
| `linear_compressibility_surface` | strain along $\mathbf{d}$ under hydrostatic pressure | one direction |
| `shear_modulus_surface` | $G$ on a plane, in a direction | **two** directions |
| `poisson_ratio_surface` | transverse contraction | **two** directions |

The last two need a second argument and are therefore not single-valued
functions of direction: for a given plane normal, $G$ and $\nu$ vary with the
in-plane direction, and the surface reports extrema over that in-plane freedom
rather than pretending one value exists. Poisson's ratio in particular can be
**negative** in some directions of some cubic crystals — auxetic behaviour that
a code assuming positivity would clip away.

## 5. What this does and does not model

| | |
| --- | --- |
| Modelled | orientation-weighted anisotropy of a single-phase aggregate |
| Bounds | Voigt (upper), Reuss (lower), rigorous |
| Hill | empirical midpoint, no variational status |
| Not modelled | grain shape and its own texture (morphological texture) |
| Not modelled | grain-boundary compliance, porosity, second phases |
| Not modelled | two-point statistics (Hashin-Shtrikman), self-consistent schemes |

## Verification

- The Voigt-Reuss ordering, isotropy of a random aggregate, and the Voigt
  convention factors, in {doc}`../examples/generated/elastic-anisotropy`.

## See also

- {doc}`../theory/elastic_anisotropy_and_homogenization` — the derivations and
  the bound proofs.
- {doc}`pole_figure_inversion` — where the orientation weights come from.
- {doc}`schmid_and_taylor` — the plastic counterpart, where Taylor plays the
  role Voigt plays here.

## References

### Normative

- Voigt, W. (1928). *Lehrbuch der Kristallphysik*. Teubner.
- Reuss, A. (1929). Berechnung der Fließgrenze von Mischkristallen.
  *Zeitschrift für Angewandte Mathematik und Mechanik* **9**, 49-58.
  <https://doi.org/10.1002/zamm.19290090104>
- Hill, R. (1952). The elastic behaviour of a crystalline aggregate.
  *Proceedings of the Physical Society A* **65**, 349-354.
  <https://doi.org/10.1088/0370-1298/65/5/307>

### Informative

- Nye, J. F. (1985). *Physical Properties of Crystals*. Oxford University Press.
- Hashin, Z. & Shtrikman, S. (1962). A variational approach to the theory of the
  elastic behaviour of polycrystals. *Journal of the Mechanics and Physics of
  Solids* **10**, 343-352.
  <https://doi.org/10.1016/0022-5096(62)90005-4>
- Kocks, U. F., Tomé, C. N. & Wenk, H.-R. (1998). *Texture and Anisotropy*.
  Cambridge University Press.
