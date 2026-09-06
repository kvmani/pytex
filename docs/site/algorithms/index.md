# Algorithms

These pages state **how** each scientific surface computes what it computes: the
mathematics, the algorithm as steps a reader could reimplement, the constraints
and tolerances with what each is calibrated against, the complexity, and the
failure modes.

They sit between the other documentation layers rather than duplicating them:

| Layer | Answers | Example |
| --- | --- | --- |
| concepts | what the objects mean | {doc}`../concepts/orientation_relationships` |
| workflows | how to get a result | {doc}`../workflows/composite_or_diffraction` |
| **algorithms** (here) | how the result is computed, and what constrains it | this section |
| worked examples | the numbers, computed live and checked | {doc}`../examples/index` |
| theory notes | the canonical LaTeX derivations | {doc}`../theory/index` |

## What is covered

The pages group by the question they answer, and each states the surface it
documents so the code and the description cannot drift apart.

**Texture and orientation**

| Page | The computation |
| --- | --- |
| {doc}`pole_figure_inversion` | measured pole figures to an ODF, by the discrete and the harmonic route |
| {doc}`ghost_correction` | recovering the odd part that a pole figure cannot measure |
| {doc}`ipf_coloring` | an orientation, a chosen specimen direction, and the colour that follows |
| {doc}`kearns_parameter` | the basal-pole second-moment tensor, by three routes that disagree predictably |
| {doc}`misorientation_and_disorientation` | the symmetry orbit, its canonical representative, and boundary statistics |

**Electron backscatter diffraction**

| Page | The computation |
| --- | --- |
| {doc}`ebsd_grains_and_local_misorientation` | grains from a point grid; KAM, GROD, GOS, GAM; GND density |
| {doc}`csl_boundaries` | assigning a Sigma value, and what it does not establish |

**Transmission electron microscopy and diffraction**

| Page | The computation |
| --- | --- |
| {doc}`saed_pattern_indexing` | phase, zone axis, orientation and indices from picked spots |
| {doc}`cbed_thickness_and_symmetry` | foil thickness from fringes; point group including the centre of symmetry |
| {doc}`tem_tilt_navigation` | the holder tilts that reach a target zone axis |
| {doc}`composite_saed_assembly` | a parent-plus-variant pattern from an orientation relationship |

**Mechanical response**

| Page | The computation |
| --- | --- |
| {doc}`schmid_and_taylor` | which system yields first, and what a strain costs |

**Orientation relationships**

| Page | The computation |
| --- | --- |
| {doc}`orientation_relationship_determination` | an OR from measured parent/child orientations |
| {doc}`variant_correspondence` | variant-resolved plane and direction correspondence |

**X-ray diffraction**

| Page | The computation |
| --- | --- |
| {doc}`precise_lattice_parameter_determination` | a cell from peak positions, with the systematic error extrapolated away |

## Conventions used throughout

Every page works two systems side by side so nothing hexagonal is left implicit:

- **cubic** — Kurdjumov-Sachs, fcc → bcc (austenite → martensite);
- **hexagonal** — Burgers, bcc → hcp (β → α titanium and zirconium).

Angles are in degrees on every public surface, lengths on a detector in
millimetres, reciprocal lengths in Å⁻¹. Orientation matrices are crystal-to-
specimen in the Bunge convention, so a parent/child pair shows the rotation
$\mathbf{V} = \mathbf{C}^{\mathsf{T}}\mathbf{P}$. Symbols follow
{doc}`../standards/terminology_and_symbol_registry`.

Each figure on these pages is **generated** by
`scripts/generate_algorithm_figures.py` rather than drawn by hand, so a diagram
cannot drift from the algorithm it illustrates, and each is held to the
repository's figure layout guards.

```{toctree}
:maxdepth: 1

orientation_relationship_determination
variant_correspondence
composite_saed_assembly
saed_pattern_indexing
tem_tilt_navigation
cbed_thickness_and_symmetry
precise_lattice_parameter_determination
pole_figure_inversion
ghost_correction
ipf_coloring
kearns_parameter
misorientation_and_disorientation
ebsd_grains_and_local_misorientation
csl_boundaries
schmid_and_taylor
```
