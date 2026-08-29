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

Every one of these — and the TEM tilt stereograms, the Kikuchi maps, and the workbench's
interactive stereogram — takes its plane coordinates from the single helper
`pytex.core.sphere.project_directions(directions, method=..., antipodal=...)`, which states both
radial laws ($r = \tan(\rho/2)$ stereographic, $r = 2\sin(\rho/2)$ equal area) and applies the
antipodal folding first. `pytex.texture.projections` re-exports it under the name the texture
subsystem has always used. There is deliberately no second implementation: a projection that
differs by a factor between two figures is invisible until two figures are compared.

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

## The Orientation-Relationship Stereogram

`plot_or_stereogram(relationship, variant=k)` draws the figure by which orientation relationships
are read in the literature: one net, in the **parent crystal frame**, carrying

- the parent pole of each parallel pair as an **open** symbol and the child pole — carried back
  into the parent frame by the variant rotation — as a **filled** one, so a parallelism reads as
  two symbols on top of each other;
- a **tie-line** joining each pair along its great circle, labelled with the deviation in degrees;
- for plane pairs, the **great circles** of both planes (parent dashed, child solid), so a plane
  parallelism reads as two coincident circles rather than two coincident points.

`variant` takes a `TransformationVariant` or a one-based index; `None` means variant 1, which is
the relationship as stated. Each variant carries its *own* parallel pair, so the labels move with
the variant instead of repeating variant 1's indices — see
{doc}`../concepts/visualization_primitives`.

`or_stereogram_pairs(...)` returns the same pairs as data (`ORStereogramPair`: both vectors in the
parent frame, the label, the deviation) for callers that want the numbers rather than the figure.

**Read the deviation label precisely.** Nominating extra parent families through `parent_planes` /
`parent_directions` routes them through `find_parallel_planes` / `find_parallel_directions`, and
the number those report — the number this figure draws and labels — is the **rationalization
residual**: the *exact* child image of any parent plane is parallel to it by construction, so what
is being measured is the angle by which the nearest low-index child index misses that exact image.
A small tolerance therefore keeps the pairs for which a low-index child object really is parallel
and drops the parent members for which none is. The drawn tie-line is exactly that gap, so the
figure and its label agree.

One numerical detail is load-bearing rather than cosmetic. A stereogram folds antipodal directions
onto one point, and the fold rule must break ties for poles lying on the equator. The two ends of a
tie-line are therefore folded **together**, from the parent pole's decision: Kurdjumov-Sachs
variants 7 and 9 have a defining direction in the equatorial plane, and the variant rotation
returns the child copy at `z = -8e-16`. Folded independently, the two ends of a zero-deviation
tie-line land on opposite rims — a diameter drawn across a figure whose entire claim is that the
two poles coincide. The worked example
{doc}`../examples/generated/visualization` pins the coincidence over all 24 variants.

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
- `docs/standards/scientific_notes_and_figures.md`
- `docs/testing/plotting_validation_matrix.md`

### Informative

- MTEX documentation: [Spherical Projections](https://mtex-toolbox.github.io/SphericalProjections.html)
- Bunge, *Texture Analysis in Materials Science* (1982)
