<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Orientation representations

The constants and identities behind the equal-volume charts of SO(3), and the inversion that names an orientation as a (hkl)[uvw] texture component.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The homochoric ball and the cubochoric cube enclose the same volume, and the cube's corner lands on the ball

Uniform sampling of orientation space requires a chart whose volume element is the invariant measure of SO(3). Two such charts exist: the homochoric ball of radius (3*pi/4)^(1/3), and the cubochoric cube of edge pi^(2/3) that the equal-volume map carries onto it. Both must enclose the same volume, pi^2, and the map must send a corner of the cube exactly onto the surface of the ball. Either identity fails immediately for a wrong constant in the map, which is why they are checked here rather than trusted.

**Symbols**

- $\mathbf{h}$ &mdash; Homochoric vector of a rotation; the equal-volume chart of SO(3).
- $\mathbf{c}$ &mdash; Cubochoric coordinate; the equal-volume chart mapped onto a cube.
- $R_{1}$ &mdash; Radius of the homochoric ball, $(3\pi/4)^{1/3}$.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CUBOCHORIC_CUBE_EDGE,
    CUBOCHORIC_CUBE_HALF_EDGE,
    FrameDomain,
    Handedness,
    HOMOCHORIC_BALL_RADIUS,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    convert_orientations,
    ideal_orientation_indices,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
ball_volume = (4.0 / 3.0) * np.pi * HOMOCHORIC_BALL_RADIUS**3
cube_volume = CUBOCHORIC_CUBE_EDGE**3
corner = np.full((1, 3), CUBOCHORIC_CUBE_HALF_EDGE)
corner_radius = float(np.linalg.norm(
    convert_orientations(corner, source='cubochoric', target='homochoric')
))
result = np.array([
    ball_volume / np.pi**2,
    cube_volume / np.pi**2,
    corner_radius / HOMOCHORIC_BALL_RADIUS,
])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `core-orientation-equal-volume-charts-agree-on-the-so3-volume` | [1.000000000000, 1.000000000000, 1.000000000000] | [1.000000000000, 1.000000000000, 1.000000000000] | &mdash; | 2.22e-16 | 1e-12 | ✅ pass |

**Why this value**: The invariant measure on SO(3), (1 - cos w) dw dOmega / pi^2, gives the group total volume pi^2 before normalization. The homochoric radial function f(w) = [3(w - sin w)/4]^(1/3) reaches R1 = (3 pi/4)^(1/3) at w = pi, so (4/3) pi R1^3 = pi^2; the cubochoric edge is fixed as pi^(2/3) by requiring the same volume. A cube corner is at maximum distance from the centre in the cube, so the volume-preserving map must send it to the ball's surface.

**Citation**: Rosca, Morawiec and De Graef, Modelling Simul. Mater. Sci. Eng. 22 (2014) 075013, doi:10.1088/0965-0393/22/7/075013; Morawiec, Orientations and Rotations (Springer, 2004), for the invariant measure.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Core API <../../api/index>`

## Naming an orientation (hkl)[uvw] inverts the construction that built it

A rolling-texture component is named by the crystal plane lying in the sheet plane and the crystal direction along the rolling direction. Orientation.from_miller turns that name into an orientation; ideal_orientation_indices turns an orientation back into the name. Round-tripping the copper component {112}<111> must return the indices it was built from, and both deviation angles must vanish - the deviations being what distinguishes an exact component from a nearest label.

**Symbols**

- $(hkl)$ &mdash; Crystal plane lying in the specimen plane of a named texture component.
- $[uvw]$ &mdash; Crystal direction along the specimen reference direction of a component.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CUBOCHORIC_CUBE_EDGE,
    CUBOCHORIC_CUBE_HALF_EDGE,
    FrameDomain,
    Handedness,
    HOMOCHORIC_BALL_RADIUS,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    convert_orientations,
    ideal_orientation_indices,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
copper = Orientation.from_miller(
    (1, 1, 2), (1, 1, -1), phase=phase, specimen_frame=specimen
)
indices = ideal_orientation_indices(copper)
result = np.array([
    *indices.hkl,
    *indices.uvw,
    indices.plane_deviation_deg,
    indices.direction_deviation_deg,
])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `core-orientation-ideal-indices-invert-the-plane-direction-construction` | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | [1.000000000, 1.000000000, 2.000000000, 1.000000000, 1.000000000, -1.000000000, 0.000000000, 0.000000000] | &mdash; | 0.00e+00 | 1e-09 | ✅ pass |

**Why this value**: An exact inverse identity, not a fitted result: the plane normal aligned with ND and the direction aligned with RD are recovered by mapping those specimen axes back through g^T and expressing them in the reciprocal and direct bases respectively, so the integer indices and zero residual angles follow by construction.

**Citation**: Bunge, Texture Analysis in Materials Science (1982), for the (hkl)[uvw] naming of rolling-texture components; Hirsch et al., Electron Microscopy of Thin Crystals (1965), for the copper component.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Core API <../../api/index>`
