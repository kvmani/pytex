# Classifying A Boundary: Coincidence-Site Lattices

**Surface:** `pytex.ebsd.csl.classify_misorientations`, `CSLType`, `CSLMatch`,
`TwinLaw`, `brandon_tolerance_deg`, `CUBIC_CSL_TYPES`, `CUBIC_TWIN_LAWS`, with
`GrainBoundaryNetwork` supplying the boundaries and the workbench operation
`ebsd.distribution` reporting the population.

Certain special grain boundary orientations exhibit low interfacial energy,
enhanced resistance to grain boundary migration, and improved resistance to
intergranular corrosion and cracking. The {ref}`coincidence-site lattice (CSL) <term-csl>`
provides the classical geometric framework for identifying these special
boundaries, serving as the quantitative foundation for grain-boundary engineering.
This page explains how a measured boundary is assigned a $\Sigma$ value, the
physical basis of angular tolerances, and the geometric limitations inherent
to boundary classification.

## 1. What $\Sigma$ counts

Superimpose the lattices of two grains, extended through each other, and rotate
one relative to the other. At special misorientations a fraction of lattice
points coincide, and those coincidences form a lattice of their own — the
coincidence-site lattice. $\Sigma$ is the reciprocal of that fraction:

$$
\Sigma \;=\; \frac{\text{volume of the CSL unit cell}}{\text{volume of the crystal unit cell}}.
$$

$\Sigma 3$ means one lattice site in three is shared; $\Sigma 29$, one in
twenty-nine. $\Sigma$ is always **odd** for cubic lattices, which is a
consequence of the lattice geometry and a useful check on any table.

The physical argument is that a boundary with many coincident sites needs less
distortion to build, so it costs less energy. That argument is geometric, and
its limits are the subject of section 5.

## 2. The registry

`CUBIC_CSL_TYPES` carries the standard cubic set as axis-angle pairs, $\Sigma 1$
through $\Sigma 29$. Some $\Sigma$ values admit **more than one** distinct
boundary, which the registry distinguishes with a variant suffix rather than
silently keeping one:

| $\Sigma$ | Angle | Axis | Note |
| --- | --- | --- | --- |
| 3 | $60.0^\circ$ | $\langle 111 \rangle$ | the coherent twin of fcc and bcc metals |
| 5 | $36.86^\circ$ | $\langle 100 \rangle$ | |
| 7 | $38.21^\circ$ | $\langle 111 \rangle$ | |
| 9 | $38.94^\circ$ | $\langle 110 \rangle$ | formed by adjacent $\Sigma 3$ variants |
| 11 | $50.47^\circ$ | $\langle 110 \rangle$ | |
| 13a / 13b | $22.62^\circ$ / $27.79^\circ$ | $\langle 100 \rangle$ / $\langle 111 \rangle$ | two distinct boundaries |
| … | | | through $\Sigma 29$a/b |

`CUBIC_TWIN_LAWS` names the $\Sigma 3$ boundary as the coherent twin, because
"twin" is a statement about a named law and not merely about a $\Sigma$ value.

## 3. The Brandon criterion

A measured boundary never sits exactly on an ideal misorientation. {ref}`Brandon's criterion <term-csl>` admits a deviation that tightens as $\Sigma$ rises:

$$
\Delta\theta_{\max}(\Sigma) \;=\; \frac{\theta_0}{\sqrt{\Sigma}},
\qquad \theta_0 = 15^\circ \text{ by default.}
$$

| $\Sigma$ | Tolerance |
| --- | --- |
| 3 | $8.66^\circ$ |
| 9 | $5.00^\circ$ |
| 29 | $2.79^\circ$ |

The $1/\sqrt{\Sigma}$ form comes from the spacing of the secondary dislocation
network that accommodates the deviation: a higher-$\Sigma$ boundary has a finer
CSL and can absorb less misfit before the dislocation cores overlap and the
special structure is destroyed.

$\theta_0$ is a **parameter, not a constant**. Brandon's $15^\circ$ is
conventional; Palumbo-Aust and other criteria give different values and
different exponents. A $\Sigma$ fraction quoted without its criterion is not
reproducible, and changing $\theta_0$ moves the number substantially.

## 4. The classification algorithm

```text
input : misorientation matrices M_i = inv(o1) @ o2, crystal operators G,
        registry, theta0, include_sigma1

1  drop Sigma1 from the candidate list unless include_sigma1   -- low-angle, not special
2  for each candidate CSL type c:
3      tol_c  <- theta0 / sqrt(sigma_c)
4      dev_i  <- symmetry-reduced angle between M_i and the ideal matrix of c
5      accept where dev_i <= tol_c
6      keep c for boundary i when it improves on the incumbent:
             smaller deviation, or
             equal deviation and smaller sigma
7  boundaries with no qualifying type return None -- "general", not "unclassified"
```

Two mathematical considerations govern classification:

**Step 4 requires symmetry reduction.** The angular deviation must be minimized
over the full bicrystal symmetry orbit, as formulated in
{doc}`misorientation_and_disorientation`. Comparing unreduced rotation matrices
evaluates only an arbitrary representative of the ideal boundary and neglects
symmetry equivalents (576 in cubic–cubic bicrystals).

**Step 6's tie-break prefers the lower $\Sigma$.** Because tolerance bands of
distinct CSL types can overlap, a boundary may fall within the acceptance threshold
of multiple types simultaneously. Preferring the smaller $\Sigma$ adheres to the
physical convention that the higher coincidence description governs interfacial
structure, while providing a deterministic classification.

**$\Sigma 1$ is excluded by default** because it corresponds to low-angle grain
boundaries. Classifying sub-boundaries as $\Sigma 1$ would conflate dislocation
cell walls with special coincidence structures.

## 5. Physical interpretation and geometric limitations

When reporting CSL distributions (such as $\Sigma 3$ fractions) in materials
characterization, several physical and stereological limitations must be observed:

- **Misorientation vs boundary plane.** A boundary's energy and mobility depend
  on five macroscopic degrees of freedom: three for misorientation and two for
  the boundary plane normal. The coherent $\Sigma 3$ boundary on $\{111\}$ exhibits
  exceptionally low interfacial energy, whereas an incoherent $\Sigma 3$ boundary
  with identical misorientation on an arbitrary plane behaves like a general
  high-angle boundary. Because standard 2D EBSD classifies boundaries by
  misorientation alone, the calculated $\Sigma 3$ fraction represents an upper bound
  on the coherent twin population.
- **CSL geometry vs interfacial energy.** While high coincidence (low $\Sigma$)
  correlates broadly with lower interfacial energy, coincidence is a geometric
  criterion rather than an energetic law; significant energy variations exist
  within any CSL category.
- **Angular tolerance dependence.** The classified fraction depends strongly on
  the chosen angular threshold $\Delta\theta_{\max}$ (Brandon, Palumbo–Aust, or
  custom criteria). Quoted fractions must always state the criterion employed.
- **Stereological sectioning effects.** Planar 2D sections sample boundary traces
  rather than true boundary surface areas, introducing stereological projection bias.

`CSLMatch` therefore carries `deviation_deg` alongside the $\Sigma$ value. A
boundary at $0.2^\circ$ from ideal and one at $8.5^\circ$ are both classified as "$\Sigma 3$",
and reporting the numerical deviation preserves this important structural distinction.

## 6. Cost

$O(n_{\text{boundaries}} \times n_{\text{types}} \times |G|)$, vectorised over
boundaries: each candidate type is tested against all boundaries at once, and
the incumbent is updated by masked comparison rather than per boundary.

## Verification

- The $\Sigma 3$ twin misorientation reproduced from the fcc twin
  correspondence, and the $\Sigma 9$ that two $\Sigma 3$ twins produce as a
  consequence the code was not told, in
  {doc}`../examples/generated/ebsd`.

## See also

- {doc}`misorientation_and_disorientation` — the reduction step 4 depends on,
  and the histogram whose $60^\circ$ spike leads here.
- {doc}`ebsd_grains_and_local_misorientation` — where the boundaries come from.
- {doc}`../concepts/ebsd_foundation` — the boundary network model.

## References

### Normative

- Brandon, D. G. (1966). The structure of high-angle grain boundaries. *Acta
  Metallurgica* **14**, 1479-1484.
  <https://doi.org/10.1016/0001-6160(66)90168-4>
- Grimmer, H., Bollmann, W. & Warrington, D. H. (1974). Coincidence-site
  lattices and complete pattern-shift in cubic crystals. *Acta Crystallographica
  A* **30**, 197-207. <https://doi.org/10.1107/S056773947400043X>

### Informative

- Randle, V. (2004). Twinning-related grain boundary engineering. *Acta
  Materialia* **52**, 4067-4081.
  <https://doi.org/10.1016/j.actamat.2004.05.031>
- Palumbo, G. & Aust, K. T. (1990). Structure-dependence of intergranular
  corrosion in high purity nickel. *Acta Metallurgica et Materialia* **38**,
  2343-2352. <https://doi.org/10.1016/0956-7151(90)90101-L>
