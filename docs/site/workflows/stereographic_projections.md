# Stereographic Projections And Wulff Nets

PyTex now exposes a publication-oriented spherical plotting surface for crystallographic directions, planes, and rotational symmetry axes.

Current runtime entry points:

- `plot_wulff_net(...)`
- `plot_crystal_directions(...)`
- `plot_crystal_planes(...)`
- `plot_symmetry_elements(...)`

These functions sit on the same explicit frame, phase, and notation model as the rest of PyTex. Directions are passed as `CrystalDirection` objects, planes as `CrystalPlane` objects, and symmetry content as `SymmetrySpec`.

## Current Scope

- stereographic and equal-area projection support for direction and plane data
- Wulff-net style grid overlays for publication-facing stereographic plots
- plane traces rendered as projected great circles
- plane poles rendered with Miller-notation annotations
- rotational symmetry-element plots with order-specific symbols
- shared YAML house-style support through the existing plotting theme system

`plot_symmetry_elements(...)` currently visualizes proper rotational symmetry axes only. Mirror planes, inversion centers, and nonsymmorphic element symbols are not yet part of this stable plotting surface.

## Direction And Plane Example

```python
import numpy as np

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    plot_crystal_directions,
    plot_crystal_planes,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

direction_figure = plot_crystal_directions(
    (
        CrystalDirection(np.array([1.0, 0.0, 0.0]), phase=phase),
        CrystalDirection(np.array([1.0, 1.0, 1.0]), phase=phase),
    ),
    labels=((1, 0, 0), (1, 1, 1)),
    theme="journal",
)

plane_figure = plot_crystal_planes(
    (
        CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
        CrystalPlane(MillerIndex([1, 0, 0], phase=phase), phase=phase),
    ),
    labels=((1, 1, 1), (1, 0, 0)),
    render="both",
    theme="journal",
)
```

## Symmetry-Element Example

```python
from pytex import plot_symmetry_elements

figure = plot_symmetry_elements(
    phase.symmetry,
    annotate_axes=True,
    theme="journal",
)
```

The rotational-order symbol mapping is intentionally semantic:

- 2-fold axes use a dyad marker
- 3-fold axes use a triangular marker
- 4-fold axes use a square marker
- 6-fold axes use a hexagonal marker

Annotations use the same Miller-style mathtext path as the crystal-scene overlays, so negative-index bars render consistently across 2D and 3D plotting surfaces.

## Interpretation Notes

- Wulff nets are currently generated from explicit projected great circles and small circles rather than from image backdrops.
- Plane labels are attached to pole locations even when only traces are rendered, because the pole remains the unambiguous annotation anchor.
- The spherical plotting surface is validated through deterministic structural checks over figure content and annotations, but PyTex does not yet claim pixel- or style-parity with MTEX.

## Related Material

- {doc}`plotting_primitives`
- {doc}`../concepts/symmetry_and_fundamental_regions`
- {doc}`../concepts/orientation_texture`
- `docs/testing/plotting_validation_matrix.md`

## References

### Normative

- `docs/standards/notation_and_conventions.md`
- `docs/standards/latex_and_figures.md`
- `docs/testing/plotting_validation_matrix.md`

### Informative

- MTEX documentation: [Spherical Projections](https://mtex-toolbox.github.io/SphericalProjections.html)
- Bunge, *Texture Analysis in Materials Science* (1982)
