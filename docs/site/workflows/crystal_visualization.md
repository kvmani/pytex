# 3D Crystal Visualization

PyTex now includes a VESTA-inspired 3D crystal viewer surface for unit cells, supercells, bond overlays, and crystallographic plane overlays.

![Crystal Visualization Example](../../figures/crystal_visualization_demo.svg)

## Scope

- unit-cell and repeated-supercell rendering
- atom glyphs with element-based color and radius scaling
- optional bond construction from covalent-radius heuristics
- lattice-edge rendering
- `(hkl)` plane overlays
- optional slab or section filtering
- camera alignment by view direction or manual angles
- publication-oriented rendering through shared YAML styles

## Scientific Model

The viewer is built from the canonical PyTex structure model:

- atoms begin in fractional coordinates from `UnitCell`
- fractional coordinates are mapped through the direct basis into Cartesian crystal-space coordinates
- supercell repetition is explicit in direct-lattice units
- planes are computed from the `(hkl)` plane equation over the supercell box

The viewer does not invent its own crystallography. It renders geometry already defined by the `Phase`, `Lattice`, and `UnitCell`.

## Example

```python
import numpy as np

from pytex import (
    AtomicSite,
    CrystalDirection,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    UnitCell,
    build_crystal_scene,
    plot_crystal_structure_3d,
)
from pytex.core.conventions import FrameDomain, Handedness

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
lattice = Lattice(5.6402, 5.6402, 5.6402, 90.0, 90.0, 90.0, crystal_frame=crystal)
unit_cell = UnitCell(
    lattice=lattice,
    sites=(
        AtomicSite("Na1", "Na", np.array([0.0, 0.0, 0.0])),
        AtomicSite("Cl1", "Cl", np.array([0.5, 0.5, 0.5])),
    ),
)
phase = Phase(
    "NaCl",
    lattice=lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=unit_cell,
)

scene = build_crystal_scene(
    phase,
    repeats=(2, 2, 2),
    plane_hkls=((1, 1, 1),),
    theme="presentation",
)
figure = plot_crystal_structure_3d(
    scene,
    view_direction=CrystalDirection(indices=np.array([1, 1, 0]), phase=phase),
)
figure.savefig("nacl_structure.png", dpi=220)
```

## Interpretation Notes

- `build_crystal_scene(...)` computes reusable scene geometry and can be inspected separately from rendering.
- `plot_crystal_structure_3d(...)` renders either a `Phase` directly or a precomputed `CrystalScene`.
- Plane overlays are rendering aides for geometric interpretation. They do not change the underlying structure semantics.

## Current Limits

- the current bond model is heuristic, not chemistry-complete
- slicing is currently a simple slab filter rather than a full interactive clipping engine
- the current 3D backend is Matplotlib-based and optimized for static publication-oriented views rather than interactive scene editing

## Related Material

- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`style_customization`
- {doc}`combined_structure_diffraction_visualization`
- [../../tex/theory/crystal_visualization_geometry.tex](../../tex/theory/crystal_visualization_geometry.tex)

## References

### Normative

- `../../standards/reference_canon.md`
- `../../standards/notation_and_conventions.md`

### Informative

- `../../testing/plotting_validation_matrix.md`
