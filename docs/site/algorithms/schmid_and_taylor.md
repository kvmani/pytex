# Schmid Factors And The Taylor Factor

**Surface:** `pytex.properties.slip.SlipSystemFamily.schmid_factors`,
`max_schmid_factor`, the family constructors (`fcc_octahedral_slip` and its
siblings), `pytex.properties.taylor.taylor_factors`,
`uniaxial_strain_tensor`, with `CrystalMap.schmid_factor_map` and
`taylor_factor_map` applying them across an EBSD map.

Two numbers connect a crystal's orientation to its plastic response, and they
answer opposite questions. The **Schmid factor** asks *which slip system yields
first* under a given stress; the **Taylor factor** asks *how much work it costs*
to impose a given strain. One is a single-system criterion under stress control;
the other is a multi-system optimisation under strain control. They are
routinely presented as interchangeable measures of "how hard this grain is" and
they are not.

## 1. The Schmid factor

### 1.1 Definition

For a slip system with plane normal $\mathbf{n}$ and slip direction
$\mathbf{b}$, under a uniaxial stress along $\mathbf{d}$, the resolved shear
stress is $\tau = \sigma\, m$ with

$$
m \;=\; \cos\phi\,\cos\lambda
    \;=\; (\mathbf{n}\cdot\mathbf{d})\,(\mathbf{b}\cdot\mathbf{d}),
$$

$\phi$ the angle between the stress axis and the plane normal, $\lambda$ that
between the stress axis and the slip direction. Slip begins on the system
reaching the critical resolved shear stress first, which under a single applied
stress is the system of largest $|m|$.

**$m$ is bounded above by exactly $1/2$**, attained at
$\phi = \lambda = 45^\circ$. A computed Schmid factor above $0.5$ is an error,
not a hard grain, and the bound is the cheapest available check on a
convention.

### 1.2 The algorithm

```text
input : slip family (plane normals N, slip directions B, in the crystal frame),
        orientations R (crystal-to-specimen), stress direction d

1  d <- d / |d|
2  map the family into the specimen frame:   n_spec = R n,   b_spec = R b
3  cos_phi    = n_spec . d
4  cos_lambda = b_spec . d
5  m = |cos_phi * cos_lambda|                 -- absolute: slip sense is immaterial
6  the grain's Schmid factor is max over the family
```

Steps 2-5 are one `einsum` over orientations and systems together, so a full
EBSD map is evaluated without a Python loop.

**Why the absolute value.** A slip system slips in either sense, so the sign of
$m$ carries no information about *whether* it slips. It does carry information
about *direction*, which matters for texture evolution — but not for the yield
criterion this factor expresses, and taking the signed maximum would report the
wrong system for half the grains.

### 1.3 What it assumes, and where that breaks

| Assumption | Consequence when false |
| --- | --- |
| Single applied stress, uniaxial | multiaxial stress needs the full tensor form, not $\cos\phi\cos\lambda$ |
| All systems share one CRSS | in hcp, basal, prismatic and pyramidal CRSS differ by factors of several, so the largest $m$ is often *not* the active system |
| Grain deforms freely | in a polycrystal it does not — which is exactly what section 2 addresses |

The hcp caveat is the one that bites in practice: ranking zirconium or magnesium
slip systems by Schmid factor alone, with no CRSS weighting, predicts the wrong
active system routinely. `SlipSystemFamily` is per-family precisely so families
are kept separable and can be weighted by their own CRSS rather than pooled.

## 2. The Taylor factor

### 2.1 The different question

A grain inside a polycrystal cannot deform freely: its neighbours impose a
strain. The **full-constraint (Taylor) model** requires each grain to accommodate
the *macroscopic* strain entirely by its own slip. Five independent shears are
needed for an arbitrary volume-conserving strain — the von Mises criterion —
so several systems must act together, and the question becomes which combination
does it most cheaply.

$$
M \;=\; \frac{\sum_s |\dot\gamma_s|}{\dot\varepsilon_{\text{eq}}},
$$

the total slip needed per unit equivalent strain. Higher $M$ means a harder
grain in this sense — more shear expended for the same shape change.

### 2.2 The algorithm: a linear program, not a formula

```text
input : slip family, orientations R, deviatoric strain tensor E (trace 0)

1  build each system's Schmid tensor in the specimen frame:
       P_s = sym( n_spec b_spec^T )
2  reduce to the 5 independent deviatoric components
3  solve, per orientation:
       minimise   sum_s gamma_s          (total slip)
       subject to sum_s gamma_s P_s = E  (the strain is exactly accommodated)
                  gamma_s >= 0
   with each system entered twice, as +P_s and -P_s, so gamma >= 0 covers both senses
4  M = (minimum total slip) / equivalent strain
5  infeasible  ->  M = inf
```

Step 3 is a genuine **linear program** — solved with `scipy.optimize.linprog`
(HiGHS) — and this is the substantive point. The Taylor factor is the *value of
an optimisation*, not a closed-form expression in the orientation. There is no
formula to inline; a code path that appears to compute a Taylor factor by
direct evaluation is computing something else.

**The doubling in step 3** is how a sign-free formulation covers both slip
senses: columns are $[+P_s \mid -P_s]$ with $\gamma \ge 0$, which is the
standard reformulation of a free-sign variable and keeps the problem linear.

### 2.3 Infeasibility is a result, not a failure

If the family cannot span the five independent strain components, no combination
of its systems can accommodate the imposed strain and the program is infeasible.
PyTex returns $\infty$, which is the honest answer: this family cannot do this
deformation.

This is not an edge case. **Basal slip alone in hcp provides only two
independent systems**, which is why magnesium and zirconium require prismatic,
pyramidal, or twinning contributions to deform polycrystalline at all — and why
a Taylor factor computed from a basal-only family returns $\infty$ rather than a
large number. That infinity is the von Mises criterion speaking.

### 2.4 Model limits

Full-constraint Taylor is an upper bound on strength: it over-constrains, since
real grains do not each accommodate the full macroscopic strain. Relaxed
-constraint and self-consistent models sit below it. PyTex implements the
full-constraint case and says so; the value is a bound with a known direction,
not an estimate with unknown error.

## 3. Maps

`CrystalMap.schmid_factor_map` and `taylor_factor_map` apply the two across an
orientation map, producing the per-pixel fields that are overlaid on
microstructure to correlate slip traces or hardness with orientation. Cost
differs sharply:

| | Cost per orientation |
| --- | --- |
| Schmid | two dot products per system; vectorised over the whole map at once |
| Taylor | **one linear program**, solved per orientation |

A Taylor map of a large scan is therefore orders of magnitude more expensive
than a Schmid map, and that is intrinsic to the definition rather than an
implementation shortcoming.

## 4. Choosing between them

| Question | Use |
| --- | --- |
| Which system yields first under an applied stress? | Schmid |
| Which grains are soft in a tensile test? | Schmid (with CRSS weighting for hcp) |
| How much slip does this imposed strain cost? | Taylor |
| Polycrystal flow stress from a texture | Taylor, averaged over the ODF |
| Texture evolution | neither alone; both feed a plasticity model |

## Verification

- The exact bound $m \le 1/2$ at $\phi = \lambda = 45^\circ$, and the Taylor
  factor of the standard fcc orientations, in
  {doc}`../examples/generated/schmid-and-taylor`.

## See also

- {doc}`../theory/schmid_and_taylor_plasticity` — the derivations, the von
  Mises five-system argument, and the model hierarchy.
- {doc}`ipf_coloring` — how a Schmid or Taylor map is displayed alongside
  orientation.
- {doc}`ebsd_grains_and_local_misorientation` — the map these fields are
  computed over.

## References

### Normative

- Schmid, E. & Boas, W. (1935). *Kristallplastizität*. Springer.
  <https://doi.org/10.1007/978-3-662-34532-0>
- Taylor, G. I. (1938). Plastic strain in metals. *Journal of the Institute of
  Metals* **62**, 307-324.
- von Mises, R. (1928). Mechanik der plastischen Formänderung von Kristallen.
  *Zeitschrift für Angewandte Mathematik und Mechanik* **8**, 161-185.
  <https://doi.org/10.1002/zamm.19280080302>

### Informative

- Kocks, U. F., Tomé, C. N. & Wenk, H.-R. (1998). *Texture and Anisotropy*.
  Cambridge University Press.
- Bishop, J. F. W. & Hill, R. (1951). A theory of the plastic distortion of a
  polycrystalline aggregate under combined stresses. *Philosophical Magazine*
  **42**, 414-427. <https://doi.org/10.1080/14786445108561065>
