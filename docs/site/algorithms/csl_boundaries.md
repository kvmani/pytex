# Classifying A Boundary: Coincidence-Site Lattices

**Surface:** `pytex.ebsd.csl.classify_misorientations`, `CSLType`, `CSLMatch`,
`TwinLaw`, `brandon_tolerance_deg`, `CUBIC_CSL_TYPES`, `CUBIC_TWIN_LAWS`, with
`GrainBoundaryNetwork` supplying the boundaries and the workbench operation
`ebsd.distribution` reporting the population.

Some grain boundaries have low energy, resist migration, resist corrosion, and
resist cracking; most do not. The **coincidence-site lattice** is the classical
geometric criterion that separates them, and grain-boundary engineering is the
practice of deliberately increasing the fraction of boundaries that satisfy it.
This page states how a measured boundary is assigned a $\Sigma$ value, what the
tolerance means, and — because it is the part most often left unsaid — what the
classification does **not** establish.

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
| 9 | $38.94^\circ$ | $\langle 110 \rangle$ | what two $\Sigma 3$ twins make |
| 11 | $50.47^\circ$ | $\langle 110 \rangle$ | |
| 13a / 13b | $22.62^\circ$ / $27.79^\circ$ | $\langle 100 \rangle$ / $\langle 111 \rangle$ | two distinct boundaries |
| … | | | through $\Sigma 29$a/b |

`CUBIC_TWIN_LAWS` names the $\Sigma 3$ boundary as the coherent twin, because
"twin" is a statement about a named law and not merely about a $\Sigma$ value.

## 3. The Brandon criterion

A measured boundary never sits exactly on an ideal misorientation. Brandon's
criterion admits a deviation that tightens as $\Sigma$ rises:

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

Two details decide the answer and are easy to get wrong:

**Step 4 must be symmetry-reduced.** The deviation is the minimum over the
crystal symmetry orbit, exactly as in
{doc}`misorientation_and_disorientation`. Comparing raw matrices measures the
distance to one arbitrary representative of the ideal boundary and misses the
other 575.

**Step 6's tie-break prefers the lower $\Sigma$.** Tolerance bands overlap, so a
boundary can lie within tolerance of two types at once. Preferring the smaller
$\Sigma$ follows the convention that the more coincident description is the
operative one, and — more importantly — it is *deterministic*, so the same
boundary does not change class between runs.

**$\Sigma 1$ is excluded by default** because it is the low-angle case: every
boundary below the tolerance would classify as $\Sigma 1$ and swamp the
statistics with a category that says only "these grains are barely
misoriented".

## 5. What a $\Sigma$ value does not tell you

This is the section to read before quoting a $\Sigma 3$ fraction as a
materials-property result.

- **It is a misorientation criterion only.** A boundary's character depends on
  its **plane** as well as its misorientation. The coherent $\Sigma 3$ on
  $\{111\}$ has very low energy; the *incoherent* $\Sigma 3$, the same
  misorientation on a different plane, does not, and behaves like a general
  high-angle boundary. Classification here uses misorientation alone, so it
  cannot separate them. A $\Sigma 3$ fraction is an upper bound on the coherent
  twin fraction.
- **Low $\Sigma$ does not guarantee low energy.** The geometric argument is a
  correlation with exceptions, and the exceptions are not rare.
- **The tolerance is a convention.** See section 3.
- **A 2-D section undercounts.** Boundary planes intersecting a polished surface
  are sampled by their trace, not by their area, so a boundary-plane
  distribution from a single section is biased.

`CSLMatch` therefore carries `deviation_deg` alongside the $\Sigma$ value. A
boundary at $0.2^\circ$ from ideal and one at $8.5^\circ$ are both "$\Sigma 3$",
and reporting only the label discards the distinction.

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
