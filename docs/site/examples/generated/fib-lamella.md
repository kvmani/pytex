<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# FIB lamella planning

The residual tilt of a vertically milled lamella checked on orientations built so it is known (0, 12.5 and 35.26 degrees), the crystal-to-specimen convention on a case that fails if it is read backwards, the lamella azimuth, the exact holder tilt for a mounted lamella against its small-angle approximation, a holder's guaranteed reach, the footprint clearance in a hand-measured grain, the uncertainty budget, and the recovery of a reversed FIB rotation from fiducials.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## A target built to rise 12.5 degrees out of the surface

A tetragonal crystal whose [001] axis is rotated 77.5 degrees about the sample X axis: [001] then points at (0, -sin 77.5, cos 77.5), which rises arcsin(cos 77.5) = 12.5 degrees out of the surface. Its orbit is only +/-[001], so no other member can come closer: eps* is exactly 12.5 degrees, the tilt the TEM holder must supply.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
```

:::

**Compute**

```python
orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))
result = lamella_geometry(about([1, 0, 0], 77.5), orbit).eps_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-constructed-residual` | 12.500000000 | 12.500000000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: arcsin(cos 77.5 deg) = 90 - 77.5 = 12.5 deg, by construction.

**Citation**: International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## A cube-oriented grain has <100> in the surface

At the identity orientation a cubic crystal's [100] and [010] lie in the surface, so a <100> lamella needs no residual tilt at all: eps* = 0, and the lamella normal points along X_s (the smallest canonical azimuth among the tied options).

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
```

:::

**Compute**

```python
orbit = target_orbit(phase('m-3m', 3.52, 3.52), (1, 0, 0))
result = lamella_geometry(np.eye(3), orbit).eps_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-in-plane-target` | 0.000000 | 0.000000 | deg | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: [100] . Z_s = 0 at the identity orientation.

**Citation**: International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## With [111] along the normal, every <001> is 35.26 degrees out

When a cubic grain has [111] along the surface normal, each <001> axis makes the same angle arccos(1/sqrt 3) with it, so each rises arcsin(1/sqrt 3) = 35.264 degrees out of the surface. That exceeds a +/-30 degree holder's guaranteed reach: no vertical lamella of this grain shows <001> for every mounting rotation.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
```

:::

**Compute**

```python
target = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
axis = np.cross(target, [0.0, 0.0, 1.0])
g = about(axis / np.linalg.norm(axis), math.degrees(math.acos(target[2])))
orbit = target_orbit(phase('m-3m', 3.52, 3.52), (0, 0, 1))
result = lamella_geometry(g, orbit).eps_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-111-along-normal` | 35.26438968 | 35.26438968 | deg | 2.75e-09 | 1e-07 | ✅ pass |

**Why this value**: arcsin(1/sqrt(3)) = 35.2643897 deg, by hand.

**Citation**: International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## The orientation matrix maps crystal to specimen

PyTex orientations map crystal to specimen, v_s = g v_c. For a triclinic crystal with g = Rx(30) Rz(40), g [100] = (cos 40, sin 40 cos 30, sin 40 sin 30), which rises arcsin(sin 40 sin 30) = 18.747 degrees out of the surface. Read backwards, g^T [100] would lie in the surface and report zero; the geometry must not do that.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
```

:::

**Compute**

```python
g = about([1, 0, 0], 30.0) @ about([0, 0, 1], 40.0)
orbit = target_orbit(phase('1', 4.0, 4.0), (1, 0, 0))
result = lamella_geometry(g, orbit).eps_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-orientation-convention` | 18.74723725 | 18.74723725 | deg | 1.04e-09 | 1e-07 | ✅ pass |

**Why this value**: arcsin(sin 40 deg x sin 30 deg) = arcsin(0.3213938) = 18.7472373 deg, by hand.

**Citation**: Morawiec, Orientations and Rotations, Springer, doi:10.1007/978-3-662-09156-2.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## The lamella normal follows the target's in-plane part

A tetragonal [001] placed at (cos 20 cos 35, cos 20 sin 35, sin 20) rises 20 degrees out of the surface with its in-plane part at 35 degrees from X_s. The lamella normal is that in-plane part, normalized, so the azimuth to cut at is theta_S = 35 degrees.

**Symbols**

- $\theta_{S}$ &mdash; Azimuth of the lamella normal in the sample frame.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
```

:::

**Compute**

```python
d = np.array([math.cos(math.radians(20)) * math.cos(math.radians(35)),
              math.cos(math.radians(20)) * math.sin(math.radians(35)),
              math.sin(math.radians(20))])
axis = np.cross([0.0, 0.0, 1.0], d)
g = about(axis / np.linalg.norm(axis), math.degrees(math.acos(d[2])))
orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))
result = lamella_geometry(g, orbit).theta_sample_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-lamella-azimuth` | 35.000000000 | 35.000000000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: atan2(sin 35, cos 35) = 35 deg, by construction.

**Citation**: International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`, {doc}`FIB lamella planning <../../workflows/fib_lamella_planning>`

## The alpha tilt for a lamella mounted at 45 degrees

A lamella whose target sits eps* = 20 degrees out of its plane, mounted on the grid at phi = 45 degrees, needs alpha = arcsin(cos phi sin eps*) = arcsin(0.70711 x 0.34202) = 13.995 degrees - not the small-angle eps* cos phi = 14.142 degrees. The value is read from the solve at phi = 45 degrees of the mounting sweep.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.
- $\phi$ &mdash; Unknown in-plane mounting rotation of the lamella on the grid.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit
def phase(point_group, a, c):
    frame = frame_catalog.crystal_frame()
    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),
                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
                 crystal_frame=frame)
def about(axis, degrees):
    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()
from pytex.fib import MountModel, sweep_mounting_rotation
```

:::

**Compute**

```python
orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))
g = about([1, 0, 0], 70.0)
geometry = lamella_geometry(g, orbit)
sweep = sweep_mounting_rotation(g, orbit, geometry,
                                mount=MountModel(phi_samples=8, solver='closed_form'))
result = abs(float(sweep.alpha_front_deg[1]))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-exact-mount-tilt` | 13.99545 | 13.99545 | deg | 4.64e-06 | 1e-05 | ✅ pass |

**Why this value**: arcsin(cos 45 deg x sin 20 deg) = arcsin(0.2418448) = 13.99545 deg, by hand; the double-tilt closed form of De Graef, doi:10.1017/CBO9780511615092.

**Citation**: De Graef, Introduction to Conventional TEM, doi:10.1017/CBO9780511615092.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## The residual a +/-30 x +/-25 degree holder reaches for every mount

On the exact tilt curve |alpha| and |beta| never exceed eps*, and each reaches it on an axis, so a symmetric rectangular holder guarantees exactly its smaller half-range: 25 degrees for alpha +/-30 and beta +/-25.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.
- $\phi$ &mdash; Unknown in-plane mounting rotation of the lamella on the grid.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.fib import guaranteed_radius_deg
from pytex.tem.stage import RectangularEnvelope
```

:::

**Compute**

```python
result = guaranteed_radius_deg(RectangularEnvelope(-30.0, 30.0, -25.0, 25.0))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-guaranteed-radius` | 25.0000 | 25.0000 | deg | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: min(30, 25) = 25 deg.

**Citation**: De Graef, Introduction to Conventional TEM, doi:10.1017/CBO9780511615092.

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## How far inside a 40 x 20 micrometre grain a lamella can sit

A rectangular grain 40 um long and 20 um tall (0.5 um step) holds a 15 x 2 um lamella whose normal is along the short side. The rectangle can grow by (20 - 2)/2 = 9 um across it before touching the boundary, and by (40 - 15)/2 = 12.5 um along it; the clearance is the smaller, 9 um.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.fib import fit_footprint
mask = np.zeros((60, 100), dtype=bool)
mask[10:50, 10:90] = True
```

:::

**Compute**

```python
result = fit_footprint(mask, step_um=(0.5, 0.5), theta_sample_deg=90.0,
                       length_um=15.0, width_um=2.0).margin_um
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-footprint-clearance` | 9.000 | 9.000 | um | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: min((20 - 2)/2, (40 - 15)/2) = 9 um, by hand.

**Citation**: Crow, Summed-area tables for texture mapping, SIGGRAPH 1984, doi:10.1145/964965.808600.

**See also**: {doc}`FIB lamella planning <../../workflows/fib_lamella_planning>`

## The expanded uncertainty of a residual tilt

EBSD accuracy 0.5 deg, grain spread 1.5 deg, an exact registration and a mount repeatability of 2 deg combine in quadrature to u = sqrt(0.25 + 2.25 + 0 + 4) = 2.5495 deg; with k = 2 the verdict is judged on eps* + 5.099 deg.

**Symbols**

- $\varepsilon^{*}$ &mdash; Residual tilt: the best member's angle out of the surface.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.fib import MountModel, UncertaintyInputs, plan_lamella
from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.symmetry import SymmetrySpec
frame = frame_catalog.crystal_frame()
nickel = Phase('nickel', lattice=Lattice(3.52, 3.52, 3.52, 90, 90, 90, crystal_frame=frame),
               symmetry=SymmetrySpec.from_point_group('m-3m', reference_frame=frame),
               crystal_frame=frame)
```

:::

**Compute**

```python
plan = plan_lamella(np.eye(3), (1, 0, 0), phase=nickel, grain_spread_deg=1.5,
                    uncertainty=UncertaintyInputs(0.5, 2.0, 2.0),
                    mount=MountModel(phi_samples=8, solver='closed_form'),
                    crosscheck=False)
result = plan.budget.expanded_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-uncertainty-budget` | 5.0990195 | 5.0990195 | deg | 1.36e-08 | 1e-06 | ✅ pass |

**Why this value**: 2 sqrt(0.5^2 + 1.5^2 + 0^2 + 2^2) = 2 sqrt(6.5) = 5.0990195 deg, by hand.

**Citation**: JCGM 100:2008, Guide to the expression of uncertainty in measurement (GUM).

**See also**: {doc}`FIB lamella geometry for a target zone axis <../../theory/fib_lamella_zone_axis_geometry>`

## Recovering a reversed FIB rotation sense and its offset

Fiducial trenches milled at pattern rotations 0, 30 and 60 degrees on an instrument with s_R = -1 and R0 = 90 degrees have their normals at 90, 60 and 30 degrees in the sample frame. The calibration recovers the offset R0 = 90 degrees (and the sense -1) from them.

:::{dropdown} Setup (imports and object construction)

```python
from pytex.fib import FiducialObservation, calibrate_chamber_from_fiducials
observations = [FiducialObservation(rho, (90.0 - rho) % 180.0) for rho in (0.0, 30.0, 60.0)]
```

:::

**Compute**

```python
chamber = calibrate_chamber_from_fiducials(observations)
result = chamber.rotation_sense * 1000.0 + (chamber.rotation_offset_deg % 180.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `fib-fiducial-calibration` | -910.000 | -910.000 | sense x 1000 + R0 (deg) | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: The sense -1 and offset 90 deg the observations were generated with: -1000 + 90.

**Citation**: Giannuzzi & Stevie, Micron 30 (1999) 197, doi:10.1016/S0968-4328(99)00005-0.

**See also**: {doc}`FIB lamella planning <../../workflows/fib_lamella_planning>`
