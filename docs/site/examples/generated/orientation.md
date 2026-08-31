<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Orientations and disorientation angles

Round-trip consistency of orientation representations and symmetry-reduced disorientation angles, checked against exact identities and the Sigma 3 twin reference.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Euler -> matrix -> Euler round trip returns the input orientation

Before trusting any orientation pipeline you must confirm that converting Bunge Euler angles to a rotation matrix and back is lossless. Here we take a non-degenerate orientation, convert to a matrix, rebuild an Orientation from that matrix, and measure the disorientation to the original. A correct implementation returns exactly zero.

**Symbols**

- $(\phi_1, \Phi, \phi_2)$ &mdash; Bunge Euler angles.
- $\mathbf{R}$ &mdash; Active rotation matrix.
- $\omega$ &mdash; Disorientation angle: minimal misorientation angle over symmetry.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cube = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
```

:::

**Compute**

```python
g = Orientation.from_euler(
    30.0, 40.0, 50.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
rebuilt = Orientation.from_matrix(
    g.as_matrix(), specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
result = float(np.degrees(g.distance_to(rebuilt)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `orientation-euler-matrix-roundtrip` | 0.0000 | 0.0000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: A representation round trip is an identity map; the disorientation must be 0 degrees.

**Citation**: Bunge, Texture Analysis in Materials Science, 1982, Chapter 2.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Orientation constructors and helpers <../../api/index>`

## Sigma 3 twin disorientation is 60 degrees about <111>

The coherent annealing twin in FCC metals is the Sigma 3 boundary, a 60-degree rotation about a <111> axis. Computing the symmetry-reduced disorientation between the cube orientation and its 60-degree/<111> partner is the standard validation that cubic crystal symmetry is applied correctly during misorientation reduction.

**Symbols**

- $\omega$ &mdash; Disorientation angle: minimal misorientation angle over symmetry.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cube = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
```

:::

**Compute**

```python
twin = Orientation.from_axis_angle(
    (1, 1, 1), np.radians(60.0),
    specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic,
)
result = float(np.degrees(cube.distance_to(twin)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `orientation-sigma3-disorientation` | 60.0000 | 60.0000 | deg | < 1e-08 | 1e-06 | ✅ pass |

**Why this value**: The Sigma 3 coincidence-site boundary is a 60-degree rotation about <111>.

**Citation**: Grimmer, Bollmann and Warrington, Acta Cryst. A30 (1974) 197; Randle, The Role of the CSL.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Symmetry and fundamental regions <../../concepts/symmetry_and_fundamental_regions>`
