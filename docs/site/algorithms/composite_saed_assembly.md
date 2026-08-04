# Composite Parent + Variant SAED Assembly

**Surface:** `pytex.diffraction.composite.simulate_composite_saed` and
`simulate_composite_saed_from_child_zone`, built on
`pytex.diffraction.kinematic.simulate_zone_axis_spots`.

Given an orientation relationship and a zone axis — of the **parent** or of a
**product variant** — this simulates the kinematic zone-axis pattern of every
requested phase on one shared detector, so the spots overlay exactly as they
would on the microscope screen.

```{figure} ../../figures/composite_saed_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 resolves the viewing direction from either a
  parent zone axis or a product-variant zone axis mapped back through the variant
  rotation, then builds the shared detector triad once. Lane 2 enumerates
  reflections, applies the centring condition, selects by excitation error,
  computes structure factors and projects. Lane 3 assembles the variant bases,
  finds coincidences and exports.
:width: 100%

The algorithm, with the constraint governing each stage.
```

## 1. Geometry: the Ewald construction and the excitation error

For a beam of wavelength $\lambda$ along $-\hat{\mathbf{z}}$, the Ewald sphere has
radius $1/\lambda$ and is centred at $-\hat{\mathbf{z}}/\lambda$. A reciprocal
lattice point $\mathbf{g}$ diffracts when it lies on that sphere. Real crystals
are thin, so the points are relrods and a reflection is excited when it lies
*near* the sphere. The signed distance along the beam is the **excitation
error**:

$$s_g = g_z - \frac{\lambda \lVert \mathbf{g} \rVert^2}{2},$$

with $g_z$ the component along the zone axis. Two consequences are worth
stating because both surprise readers:

- A reflection lying exactly in the zero-order Laue zone has $g_z = 0$, so
  $s_g = -\lambda\lVert\mathbf{g}\rVert^2/2$ — **not** zero. Exact Bragg
  condition is not the same as being on the zone.
- $\lambda$ is small ($0.025079$ Å at 200 kV, $0.019687$ Å at 300 kV,
  relativistically), so the sphere is nearly flat over the accessible
  $\lVert\mathbf{g}\rVert$ range and a zone-axis pattern looks like a planar
  section of the reciprocal lattice.

A reflection is kept when $\lvert s_g \rvert \le$
`max_excitation_error_inv_angstrom`. The default of 0.05 Å⁻¹ keeps every
zero-order-Laue-zone reflection within the default $\lVert\mathbf{g}\rVert$ range
while excluding higher-order zones.

The detector position follows from the camera constant $L\lambda$:

$$\mathbf{r}_{\text{mm}} = (L\lambda)\,
\begin{pmatrix} \mathbf{g}\cdot\hat{\mathbf{u}} \\
\mathbf{g}\cdot\hat{\mathbf{v}} \end{pmatrix},
\qquad d = \frac{1}{\lVert \mathbf{g} \rVert} .$$

Note the radius uses the **in-plane** part of $\mathbf{g}$ and $d$ uses the full
vector, so $r < (L\lambda)\lVert\mathbf{g}\rVert$ by exactly the out-of-plane
component the excitation error records.

## 2. The shared detector basis

Everything rests on one construction. Given a parent zone direction
$\mathbf{z}_p$, `zone_basis_from_axis` returns an orthonormal right-handed triad
$(\hat{\mathbf{u}}, \hat{\mathbf{v}}, \hat{\mathbf{z}})$ with
$\hat{\mathbf{z}} = \mathbf{z}_p / \lVert \mathbf{z}_p \rVert$ and
$\hat{\mathbf{u}} \times \hat{\mathbf{v}} = \hat{\mathbf{z}}$. An optional
in-plane reference $\mathbf{g}_{\text{align}}$ places a nominated reflection
along $+\hat{\mathbf{u}}$.

Each variant's basis is the parent's rotated into that child's frame:

$$\mathbf{B}_k = \mathbf{V}_k \mathbf{B}_{\text{parent}} .$$

This is algebraically identical to pulling child reciprocal vectors back into the
parent frame before projecting, so every sub-pattern is physically consistent on
one detector — that is what makes a composite pattern meaningful rather than a
collage.

The child zone axis is then $\mathbf{z}_c = \mathbf{V}_k \mathbf{z}_p$, which is
**generally irrational**. The exact direction drives the simulation and a nearest
rational label is reported with its angular deviation, so a label is never
mistaken for the geometry.

## 3. Anchoring on a product zone instead

The derivation's natural choice is a parent zone axis; the microscope's is a
low-index zone of the *product*. `simulate_composite_saed_from_child_zone` takes
the latter and maps it back through the anchor variant:

$$\mathbf{z}_p = \mathbf{R}_k^{\mathsf{T}} \mathbf{z}_c .$$

That parent direction is generally irrational, and is reported exactly alongside
its nearest rational label — the same honesty child zone axes receive, so
neither crystal is privileged in the output. The basis is then built by the same
`zone_basis_from_axis` call, which gives a **testable identity**:

:::{admonition} Identity: the two anchoring routes agree exactly
:class: tip

Anchoring on variant $k$'s image of a parent zone reproduces the parent-anchored
pattern for that zone **exactly** — measured to $10^{-13}$ mm on every
sub-pattern, for four anchor variants — because
$\mathbf{R}_k^{\mathsf{T}}(\mathbf{R}_k \mathbf{z}_p) = \mathbf{z}_p$ and both
routes then call one function. There is one detector-geometry definition, not
two that must be kept in step.
:::

## 4. What a composite pattern shows

### Hexagonal: Burgers down $\beta\,[110]$

The defining plane parallelism $\{110\}_{\beta} \parallel (0001)_{\alpha}$ means
that looking down a $\langle 110 \rangle_{\beta}$ axis looks straight down the
$c$-axis of the variants whose basal plane is *that* $\{110\}$. Simulating all
12 variants at 200 kV with a 180 mm·Å camera constant:

| variant | nearest child zone | deviation | reflections |
| --- | --- | --- | --- |
| 1 | $[0001]$ | **0.000°** | 54 |
| 10 | $[000\bar{1}]$ | **0.000°** | 54 |
| 5, 6, 8, 12 | $[1\,4\,\bar{5}\,\bar{3}]$ etc. | 1.020° | 21 |
| 3, 4, 7, 11 | $[\bar{5}\,5\,0\,3]$ etc. | 1.187° | 18 |
| 2, 9 | $[\bar{8}\,1\,7\,0]$ | 1.322° | 10 |

Two variants exactly on zone showing the full six-fold basal pattern, the rest
tilted off it and contributing fewer, weaker spots — the composite's
characteristic appearance, derived rather than asserted. The parent contributes
30 reflections, and the closest parent/child coincidence is
$(\bar{1}10)_{\beta}$ against $(000\bar{2})_{\alpha}$ at **0.15 mm**, because
$d_{110}^{\text{bcc}} = a/\sqrt{2}$ and $d_{0002}^{\text{hcp}} = c/2$ are nearly
equal. That near-superposition is the practical TEM signature of the Burgers
relationship.

### Cubic: Kurdjumov-Sachs down $\gamma\,[01\bar{1}]$

The defining direction parallelism
$\langle 10\bar{1}\rangle_{\gamma} \parallel \langle 11\bar{1}\rangle_{\alpha'}$
makes at least one child zone exactly rational: the minimum child-zone deviation
over the 24 variants is $0.0000°$, and the maximum is $5.264°$ — the
Kurdjumov-Sachs to Nishiyama-Wassermann separation reappearing as a zone-axis
deviation. 34 parent reflections, 338 spots in total.

## 5. Constraints

:::{admonition} Constraint: declare the space group, or absences are assumed away
:class: warning

Lattice centring is read from the **first letter of the space-group symbol**, and
`ReflectionCondition.from_phase` falls back to primitive when a phase carries
none. A body-centred phase supplied without its symbol is therefore simulated as
primitive and keeps reflections its real structure forbids — with nothing in the
spot list to say so.

`pattern.centering_audit()` reports, per phase, the centring applied and whether
it was **declared** or **assumed**; `describe()`, the reflection table and the
manifest carry the same statement, and an assumed centring produces an explicit
warning. If a simulated bcc pattern shows a $\{100\}$ reflection, this is why.

This was not hypothetical: the repository's own shared Burgers worked-example
setup declared no space groups and had been listing forbidden $\beta$
reflections until the audit was built.
:::

:::{admonition} Constraint: intensities are per sub-pattern, and cannot be shared
:class: important

Each sub-pattern's intensities are normalized to its own maximum. Kinematic
theory defines **no** intensity ratio between two different phases, so comparing
a $\beta$ spot's intensity with an $\alpha$ spot's is meaningless.

A shared normalization option was considered and **rejected**, not deferred: it
would manufacture a number the theory does not support. Compare within one
source only; `describe()` states this.
:::

:::{admonition} Constraint: sort keys are quantized
:class: warning

Spots sort by decreasing intensity, then detector radius, then lexicographic
$hkl$. Symmetry-equivalent reflections have mathematically equal intensity and
radius that differ by $\sim 10^{-14}$ depending on how the basis was built, so
raw keys let floating-point noise decide the order before the exact index
tie-break was ever reached — and the same pattern reached by the two anchoring
routes came out correctly positioned but **permuted**.

Both continuous keys are now quantized before sorting: 1 pm of detector radius
and $10^{-12}$ of full-scale intensity, far below anything physical and far above
the noise they suppress.
:::

## 6. Cost

The engine is fully vectorized: the reflection cube, the centring mask, the
excitation errors, the structure factors and the projection are each one array
operation over all reflections. Enumeration is
$\mathcal{O}\!\left((2 m + 1)^3\right)$ in `max_index`, filtered down long before
the per-reflection work, and the whole composite is
$\mathcal{O}(n_{\text{variants}})$ such passes. No Python loop runs per
reflection.

## 7. Limits

- **Kinematic only.** No dynamical (Bloch-wave, multi-beam) intensities. The
  intensities rank reflections; they do not predict what a plate of a given
  thickness will show.
- **Zero-order Laue zone only.** No HOLZ rings, and no double-diffraction spots
  — kinematically forbidden reflections excited via $\mathbf{g}_1+\mathbf{g}_2$
  paths are not modelled.
- Structure factors use an atomic-number electron scattering proxy, not a full
  parameterized scattering table. Adequate for ranking, not for quantitative
  intensity work.

## Verification

| Claim | Where it is checked |
| --- | --- |
| $d = 1/\lVert g\rVert$ and $r = (L\lambda)\lVert g_{\text{in-plane}}\rVert$ on every row | `tests/unit/test_composite_saed_export.py` |
| Friedel symmetry for a centrosymmetric phase | same |
| The two anchoring routes agree to $10^{-13}$ mm, four anchor variants | `tests/unit/test_composite_saed_child_anchor.py` |
| Sort order stable under a $10^{-15}$ perturbation | same |
| Declared vs assumed centring, and the warning | `test_composite_saed_export.py` |
| Burgers $\langle 110\rangle_{\beta}$ maps exactly onto $[0001]_{\alpha}$ | [worked examples](../examples/generated/composite-diffraction.md) |
| $\{110\}_{\beta}$ / $(0002)_{\alpha}$ separation 0.15450 mm, analytic | same |

## See also

- {doc}`variant_correspondence` — the variant machinery underneath
- {doc}`saed_pattern_indexing` — the inverse problem
- {doc}`../workflows/composite_or_diffraction` — how to drive it

## References

### Normative

- {doc}`../architecture/diffraction_foundation`
- {doc}`../standards/data_contracts_and_manifests`

### Informative

- De Graef, *Introduction to Conventional Transmission Electron Microscopy*,
  Cambridge University Press, 2003 — Ewald construction, excitation error,
  relativistic electron wavelength.
- Williams and Carter, *Transmission Electron Microscopy*, 2nd ed. — camera
  constant and SAED practice.
- Burgers, *Physica* 1 (1934) 561.
