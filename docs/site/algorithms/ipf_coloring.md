# Colouring An Orientation: The IPF Colour Key

**Surface:** `pytex.plotting.ipf.IPFColorKey`, `ipf_color`, `ipf_colors`,
`plot_ipf_key`, resting on
`SymmetrySpec.fundamental_sector`,
`SymmetrySpec.reduce_vectors_to_fundamental_sector`, and
`pytex.core.sphere.FundamentalSector`, with the workbench operations
`ebsd.map` and `texture.inverse_pole_figure`.

An inverse-pole-figure map answers one question per pixel: **which crystal
direction is parallel to a chosen specimen direction**, and it answers it in
colour. This is the most-reproduced figure in the whole of EBSD, and it is also
the most-misread, because two IPF maps are comparable only when three separate
choices behind them agree. This page states the algorithm and then states those
choices, because a figure whose colour key is not declared is not a measurement.

## 1. What is being coloured

Given an orientation $g$ (crystal-to-specimen, Bunge) and a chosen specimen
direction $\mathbf{y}$ — ND, RD, TD, or anything else — the crystal direction
lying along $\mathbf{y}$ is

$$
\mathbf{h} \;=\; g^{-1}\,\mathbf{y},
$$

expressed in the crystal frame. **The colour encodes $\mathbf{h}$, not $g$.**
That is worth being blunt about: an IPF map does not show orientation. It shows
one *component* of orientation — the rotation of the crystal about $\mathbf{y}$
is discarded entirely. Two grains with identical IPF colour can differ by any
rotation about $\mathbf{y}$, which for a cubic phase can be up to $45^\circ$ of
misorientation. A grain boundary invisible in an IPF map is not thereby a small
-angle boundary, and the map is not a substitute for a misorientation figure.

## 2. Folding into the fundamental sector

Symmetrically equivalent directions are the same direction, so they must take the
same colour. The crystal point group partitions the sphere into equivalent
regions; one of them — the **fundamental sector** — holds exactly one
representative of each orbit. For $m\bar{3}m$ it is the familiar standard
stereographic triangle with corners $[001]$, $[101]$, $[111]$.

```text
h  <- normalize(g^-1 y)
h~ <- reduce_vectors_to_fundamental_sector(h, antipodal=...)
```

`antipodal` is a genuine choice, not a default to be ignored. A *direction* and
its reverse are physically distinct; a *plane normal*, and any quantity a
diffraction experiment produces, is not. Setting `antipodal=True` (the default)
folds $\mathbf{h}$ and $-\mathbf{h}$ together, which matches the Friedel
symmetry of diffraction and halves the sector.

The sector is derived from the declared symmetry, so **two maps of different
symmetries are not colour-comparable** — the same RGB triple means a different
direction in each. `IPFColorKey` carries its `crystal_symmetry` for exactly this
reason, and refuses a symmetry with no crystal-domain reference frame rather
than guessing one.

## 3. From a position in the sector to a colour

The sector is a spherical triangle with three corner directions
$\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3$. Colouring is barycentric in those
corners:

```text
1  B <- [v1 v2 v3] as columns
2  c <- solve(B, h~)              -- barycentric coordinates of the reduced direction
3  c <- max(c, 0)                 -- clip the small negatives of a boundary point
4  c <- c / sum(c)                -- normalise to the simplex
5  rgb <- c . [red, green, blue]  -- the corner colours
6  rgb <- rgb ** (1 / saturation_gamma)
7  rgb <- rgb / max(rgb)          -- rescale so the brightest channel saturates
```

Steps 1-5 are the mathematics: each sector corner is assigned a primary, and any
direction inside takes the mixture given by its barycentric position. The corners
therefore come out pure red, green and blue, which is what makes $[001]$, $[101]$
and $[111]$ readable at a glance on a cubic map.

Steps 6 and 7 are **presentation, with no crystallographic content**, and they
are where the colours a reader actually sees are decided:

- `saturation_gamma` (default $0.5$, so channels are squared) pushes mixtures
  away from grey. Without it the interior of the triangle is a wash of muted
  tones and grain contrast is poor. It is a contrast control, and it is
  registered as a separate symbol from the lattice angle $\gamma$ precisely
  because they are unrelated.
- The final division by the largest channel makes every colour fully saturated,
  so the map uses the whole gamut rather than a dim corner of it.

Neither step preserves distance: **colour distance is not misorientation**.
Nothing in this construction is a metric, and reading "these two grains look
similar" as "these two grains are close in orientation" is unsupported.

### Low-symmetry sectors

A sector with fewer than three corners — a hemisphere, or a wedge — has no
triangle to be barycentric in. The implementation then anchors the colour basis
on the reference octant rather than failing, which keeps low-symmetry maps
colourable at the cost of a key that is conventional rather than
symmetry-derived. Declare it when publishing such a map.

## 4. The legend is part of the figure

`legend_mesh` tiles the sector on a polar/azimuth grid at `resolution_deg`,
keeps the directions the sector actually contains, projects them
(stereographic by default), and colours them with the *same* function used for
the data. That sharing is the point: a legend drawn by any other route can drift
from the map it explains, and a drifted legend is worse than none.

`boundary_points_2d` returns the projected sector outline for the key's border.

**An IPF map published without its key is unreadable**, because the key is what
declares the symmetry, the specimen direction, and the sector — the three
choices of section 5.

## 5. The three declarations that make two maps comparable

| Choice | Where it lives | What changes if it differs |
| --- | --- | --- |
| **Crystal symmetry** | `crystal_symmetry` | the sector itself; the same RGB means a different direction |
| **Specimen direction** | `specimen_direction` (ND, RD, TD, …) | which component of orientation is shown at all |
| **Antipodal folding** | `antipodal` | the sector's size, hence the whole colour assignment |

Two maps agreeing on all three are comparable. Two maps differing in any one are
not, however similar they look — and looking similar is exactly the trap, since
the colour scheme is recognisable while the choices behind it are not visible in
the image.

`saturation_gamma` and the corner colours change appearance without changing
meaning, so they need declaring for reproduction but do not break comparability.

## 6. Cost and constraints

| | |
| --- | --- |
| Cost | one symmetry reduction and one $3\times3$ solve per pixel, fully vectorised over the map |
| Refusal | a reduced direction outside the sector cone (zero barycentric sum) raises rather than returning grey — it means the reduction and the sector disagree, which is a defect, not a data property |
| Refusal | a non-crystal reference frame, or a non-positive `saturation_gamma`, is rejected at construction |
| Legend | `resolution_deg` must lie in $(0, 15]$; coarser is not a legend |

## 7. How the rest of PyTex uses it

| Consumer | Uses the key for |
| --- | --- |
| `ebsd.map` (workbench) | the orientation map, greyed by any measured channel, boundaries over the top |
| `texture.inverse_pole_figure` | the scatter or density IPF of a whole orientation set |
| `plot_ipf_key` | the standalone legend that must accompany both |
| `pytex.plotting.ebsd` | grain-boundary overlays drawn on the coloured map |

Because the colour key is one object, a map and its legend cannot disagree, and
a change of symmetry propagates to both.

## Verification

- IPF colouring of the cubic sector corners, and the invariance of colour under
  symmetry operations, in {doc}`../examples/generated/ipf-coloring`.
- Fundamental-sector reduction itself, in
  {doc}`../examples/generated/crystal_geometry`.

## See also

- {doc}`../theory/ipf_color_keys` — the canonical derivation and the sector
  geometry for each Laue class.
- {doc}`../theory/fundamental_region_reduction` — reduction to the fundamental
  sector, and why it is well defined.
- {doc}`../concepts/symmetry_and_fundamental_regions` — what a fundamental
  region is.
- {doc}`pole_figure_inversion` — the other direction: from measured pole
  densities to a distribution.

## References

### Normative

- Nolze, G. & Hielscher, R. (2016). Orientations - perfectly colored. *Journal
  of Applied Crystallography* **49**, 1786-1802.
  <https://doi.org/10.1107/S1600576716012942>

### Informative

- Engler, O. & Randle, V. (2010). *Introduction to Texture Analysis:
  Macrotexture, Microtexture, and Orientation Mapping*, 2nd ed. CRC Press.
  <https://doi.org/10.1201/9781420063660>
- Schwartz, A. J., Kumar, M., Adams, B. L. & Field, D. P., eds. (2009).
  *Electron Backscatter Diffraction in Materials Science*, 2nd ed. Springer.
  <https://doi.org/10.1007/978-0-387-88136-2>
