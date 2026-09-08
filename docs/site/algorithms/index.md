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

**Texture and orientation analysis**

| Page | Computational Core |
| --- | --- |
| {doc}`pole_figure_inversion` | Forward projection operator assembly, discrete simplex regularisation, and spherical harmonic expansion |
| {doc}`ghost_correction` | Group character projection, odd-harmonic basis construction, and non-negative profile regularisation |
| {doc}`ipf_coloring` | Directional projection $g^{-1}\mathbf{y}$, fundamental sector reduction, and barycentric RGB gamut mapping |
| {doc}`kearns_parameter` | Basal orientation tensor $\mathbf{A}$, principal axis spectral decomposition, and kernel deconvolution |
| {doc}`misorientation_and_disorientation` | Rotational symmetry cosets, Mackenzie statistical baseline, and canonical fundamental region reduction |

**Electron backscatter diffraction (EBSD)**

| Page | Computational Core |
| --- | --- |
| {doc}`ebsd_grains_and_local_misorientation` | Spatial flood-fill segmentation, misorientation gradient fields (KAM, GROD, GOS), and GND scaling |
| {doc}`csl_boundaries` | Brandon criterion tolerance boundaries, CSL rotation registry lookup, and metric tie-breaking |
| {doc}`kikuchi_band_geometry` | Gnomonic band projection, Bragg bandwidth calculation, and detector frame transformations |

**Transmission electron microscopy and diffraction**

| Page | Computational Core |
| --- | --- |
| {doc}`saed_pattern_indexing` | Invariant edge-ratio triangle matching, zone-axis determination, and dynamical intensity gating |
| {doc}`cbed_thickness_and_symmetry` | Two-beam dynamical fringe fitting, specimen thickness extraction, and diffraction group symmetry |
| {doc}`tem_tilt_navigation` | Double-tilt gimbal angle kinematics, rotation matrix decomposition, and shortest-path tilt planning |
| {doc}`composite_saed_assembly` | Orientation-relationship variant mapping, reciprocal lattice transformation, and composite pattern synthesis |

**Crystal mechanics and properties**

| Page | Computational Core |
| --- | --- |
| {doc}`schmid_and_taylor` | Single-system resolved shear stress optimization and multi-system full-constraint linear programming |
| {doc}`elastic_homogenization` | Fourth-rank Cartesian tensor rotation, Voigt–Reuss–Hill bounds, and directional modulus surfaces |

**Orientation relationships and microstructural transformations**

| Page | Computational Core |
| --- | --- |
| {doc}`orientation_relationship_determination` | Double coset symmetry filtering, variant cluster absorption, and habit plane/direction alignment |
| {doc}`variant_correspondence` | Covariant plane and contravariant direction metric transformations across non-cubic lattices |
| {doc}`parent_grain_reconstruction` | Graph voting, variant pair disorientation inversion, and median orientation voting |

**Powder X-ray diffraction**

| Page | Computational Core |
| --- | --- |
| {doc}`phase_identification` | Continuous wavelet peak detection, bipartite Hungarian matching, and multi-criteria figure of merit |
| {doc}`precise_lattice_parameter_determination` | Doublet-constrained profile fitting, iterative re-indexing, and generalized Cohen least-squares |
| {doc}`rietveld_refinement` | Whole-pattern profile least squares, Caglioti instrumental broadening, and March–Dollase texture fitting |

## Conventions and units

To ensure comprehensive crystallographic coverage across all crystal systems, algorithms
consistently address both high-symmetry cubic systems and non-orthogonal hexagonal systems:
- **Cubic systems:** Exemplified by Kurdjumov–Sachs and Nishiyama–Wassermann $\mathrm{fcc} \to \mathrm{bcc}$
  transformations (austenite $\to$ martensite/ferrite).
- **Hexagonal close-packed systems:** Exemplified by Burgers $\mathrm{bcc} \leftrightarrow \mathrm{hcp}$
  transformations ($\beta \leftrightarrow \alpha$ titanium and zirconium alloys).

Standard physical and geometrical units are adopted across all public interfaces:
- Angles: Degrees ($^\circ$) for orientation parameters and detector angles; radians in internal trigonometric routines.
- Reciprocal distances: Inverse ångströms ($\text{Å}^{-1}$) or inverse nanometers ($\text{nm}^{-1}$).
- Real-space lengths: Nanometers ($\text{nm}$) or ångströms ($\text{Å}$) for crystal unit cells; millimeters ($\text{mm}$) for detector dimensions.
- Orientations: Crystal-to-specimen rotation matrices in Bunge $(\phi_1, \Phi, \phi_2)$ Euler angle conventions.
- Symbols and notation: Governed strictly by {doc}`../standards/terminology_and_symbol_registry`.

All architectural diagrams in this section are generated directly from code specifications
via `scripts/generate_algorithm_figures.py`, guaranteeing permanent synchronization between
published diagrams and underlying software implementations.

```{toctree}
:maxdepth: 1

orientation_relationship_determination
variant_correspondence
composite_saed_assembly
saed_pattern_indexing
tem_tilt_navigation
cbed_thickness_and_symmetry
phase_identification
precise_lattice_parameter_determination
pole_figure_inversion
ghost_correction
ipf_coloring
kearns_parameter
misorientation_and_disorientation
ebsd_grains_and_local_misorientation
csl_boundaries
kikuchi_band_geometry
schmid_and_taylor
elastic_homogenization
parent_grain_reconstruction
rietveld_refinement
```
