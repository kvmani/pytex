# Plotting Semantic Primitives

PyTex now includes a reusable plotting subsystem for the implemented semantic primitives.

The output policy is intentionally split in two:

- runtime plotting APIs return ordinary Matplotlib figures for user code
- canonical repository-tracked documentation figures remain `.svg` assets

## Runtime Surface

The current runtime plotting surface includes:

- `plot_vector_set(...)`
- `plot_wulff_net(...)`
- `plot_crystal_directions(...)`
- `plot_crystal_planes(...)`
- `plot_symmetry_orbit(...)`
- `plot_symmetry_elements(...)`
- `plot_euler_set(...)`
- `plot_quaternion_set(...)`
- `plot_rotations(...)`
- `plot_orientations(...)`
- `plot_pole_figure(...)`
- `plot_inverse_pole_figure(...)`
- `plot_odf(...)`
- `plot_xrd_pattern(...)`
- `plot_saed_pattern(...)`
- `plot_crystal_structure_3d(...)`

These routines accept PyTex semantic objects and validate frame, symmetry, and convention meaning before rendering.

For spherical crystallographic plotting, PyTex now distinguishes between:

- semantic direction markers projected from `CrystalDirection`
- semantic plane poles and great-circle traces projected from `CrystalPlane`
- rotational symmetry-axis plots rendered with order-specific markers
- reusable Wulff-net overlays generated from explicit projected circles rather than raster backdrops

For texture plots, PyTex now distinguishes between:

- point or support scatter views for direct inspection of discrete data
- classical contour pole figures rendered from a smoothed projected density grid
- ODF views rendered either as Euler-space density contours or as kernel-smoothed classical Bunge sections

These texture plotting surfaces are the same ones used to inspect the outputs of
`ODF.reconstruct_pole_figure(...)` and `ODF.invert_pole_figures(...)`.

## Documentation Surface

For canonical docs assets, PyTex keeps `.svg` as the required tracked format. The helper
`save_documentation_figure_svg(...)` exists so documentation figures can be produced from the
same semantic plotting system without weakening the repository SVG rule.

The same plotting subsystem also now carries a shared YAML style system:

- built-in themes such as `journal`, `presentation`, and `dark`
- deep-merge application of base defaults, theme data, YAML files, and runtime overrides
- reuse across texture plots, diffraction plots, and 3D crystal visualization

The public helpers `list_style_themes(...)`, `load_style_theme(...)`, `read_style_yaml(...)`, and
`resolve_style(...)` are part of that stable plotting support.

The shared style system now also covers the spherical plotting surface, including Wulff-net line
weights, label backgrounds, marker palettes, and publication-oriented symmetry-element symbols.

## Example

```python
import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    plot_orientations,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

orientations = OrientationSet.from_euler_angles(
    np.array([[0.0, 0.0, 0.0], [45.0, 35.0, 10.0], [90.0, 0.0, 0.0]]),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
    phase=phase,
)

figure = plot_orientations(orientations, representation="euler")
figure.savefig("orientation_scatter.png", dpi=150)
```

## Texture Contours And ODF Sections

```python
from pytex import KernelSpec, ODF, plot_odf, plot_pole_figure
from pytex.core.lattice import CrystalPlane, MillerIndex

odf = ODF.from_orientations(
    orientations,
    weights=[4.0, 2.0, 1.0],
    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=12.0),
)
pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
pole_figure = odf.reconstruct_pole_figure(pole, include_symmetry_family=True, antipodal=True)

pf_figure = plot_pole_figure(pole_figure, kind="contour", levels=14)
pf_figure.savefig("pole_figure_contours.png", dpi=150)

odf_sections = plot_odf(
    odf,
    kind="sections",
    section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
    levels=12,
)
odf_sections.savefig("odf_bunge_sections.png", dpi=150)
```

## Stereographic Directions, Planes, And Symmetry

```python
from pytex import CrystalDirection, plot_crystal_directions, plot_crystal_planes, plot_symmetry_elements

direction_figure = plot_crystal_directions(
    (
        CrystalDirection(np.array([1.0, 0.0, 0.0]), phase=phase),
        CrystalDirection(np.array([1.0, 1.0, 1.0]), phase=phase),
    ),
    labels=((1, 0, 0), (1, 1, 1)),
)

plane_figure = plot_crystal_planes(
    (pole,),
    labels=((1, 1, 1),),
    render="both",
)

symmetry_figure = plot_symmetry_elements(
    phase.symmetry,
    annotate_axes=True,
)
```

## Design Notes

- Plotters are built on shared figure-spec builders rather than each routine owning private Matplotlib logic.
- Runtime plotting is not allowed to redefine frame or symmetry semantics.
- The plotting layer is part of the user-facing product surface, not a notebook-only convenience layer.
- The contour and section paths reuse the same shared density-grid and render infrastructure as the rest of the plotting subsystem.
- Diffraction and crystal-visualization plots reuse the same style-resolution layer even when they require different rendering primitives.
- Stereographic direction, plane, and symmetry-element plots use the same notation registry as the 3D crystal-scene overlays, so Miller-style labels render consistently across 2D and 3D surfaces.

## Related Material

- {doc}`../concepts/orientation_texture`
- {doc}`stereographic_projections`
- {doc}`../tutorials/notebooks`
- {doc}`ipf_colors`
