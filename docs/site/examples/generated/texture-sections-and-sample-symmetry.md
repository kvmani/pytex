<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Sample symmetry, ODF sections and component fractions

The texture-analysis surfaces behind the workbench's measured-texture panel, each checked against its definition: axial sample symmetry removes exactly the azimuthal part of a figure, a random cubic texture holds the Haar volume of 24 balls within a tolerance, and a hexagonal ODF under orthorhombic sample symmetry spans 90, 90 and 60 degrees.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Axial sample symmetry removes the azimuthal part of a pole figure, and nothing else

A drawn wire, an extruded rod or a tube read along its axis is taken to be axially symmetric: its texture is unchanged by any rotation about the axis. Imposing that symmetry on a measured pole figure replaces each value by the average over its ring of constant tilt. Take a figure measured on a goniometer raster - tilt 0 to 75 degrees in 5 degree steps, azimuth every 5 degrees - whose intensity is 2 + 0.6x + 0.3xy + 0.5z^2 in specimen coordinates. The x and xy terms depend on azimuth and must vanish; the z^2 term depends on tilt alone and must survive unchanged. The example returns the largest departure of the symmetrized figure from 2 + 0.5z^2.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np

from pytex import FrameDomain, Handedness, Lattice, Phase, ReferenceFrame, SymmetrySpec
from pytex.core.lattice import CrystalPlane
from pytex.texture import PoleFigure, euler_section_ranges, impose_sample_symmetry
from pytex.texture.components import random_component_fraction

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
SPECIMEN = ReferenceFrame(
    "sample_rd_td_nd", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
)
HEXAGONAL = SymmetrySpec.from_point_group("6/mmm", reference_frame=CRYSTAL)
CUBIC = SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL)
LATTICE = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=CRYSTAL)
ZIRCONIUM = Phase("alpha_zr", lattice=LATTICE, symmetry=HEXAGONAL, crystal_frame=CRYSTAL)
BASAL = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)
```

:::

**Compute**

```python
psi, phi = np.meshgrid(np.arange(0.0, 76.0, 5.0), np.arange(0.0, 360.0, 5.0), indexing="ij")
theta, azimuth = np.radians(psi.ravel()), np.radians(phi.ravel())
directions = np.column_stack(
    [np.sin(theta) * np.cos(azimuth), np.sin(theta) * np.sin(azimuth), np.cos(theta)]
)
x, y, z = directions.T
figure = PoleFigure(
    pole=BASAL,
    sample_directions=directions,
    intensities=2.0 + 0.6 * x + 0.3 * x * y + 0.5 * z**2,
    specimen_frame=SPECIMEN,
    sampling="sampled_density",
)
axial = impose_sample_symmetry(figure, "axial")
result = float(np.max(np.abs(np.asarray(axial.intensities) - (2.0 + 0.5 * z**2))))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `axial-sample-symmetry-keeps-the-polar-part` | < 1e-11 | 0.00e+00 | m.r.d. | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: Analytic: on a ring of constant tilt theta, x = sin(theta) cos(psi) and xy = sin^2(theta) sin(psi) cos(psi) have zero mean over the azimuth psi, while z = cos(theta) is constant. For azimuths equally spaced round the full circle the discrete means of cos(psi) and sin(2 psi) are exactly zero as well, so the residual is floating-point rather than quadrature error.

**Citation**: H.-J. Bunge, Texture Analysis in Materials Science (1982), section 4.2: a fibre (cylindrical) sample symmetry makes the pole figure independent of the azimuth about the fibre axis.

**See also**: {doc}`Texture analysis in the workbench <../../workflows/texture_analysis_workbench>`

## A random cubic texture holds 2.28% of its volume within 15 degrees of any orientation

A component volume fraction - '9% cube within 15 degrees' - is only meaningful beside what a texture-free specimen gives for the same tolerance. Under the invariant measure on rotations, a ball of radius w holds (w - sin w)/pi of orientation space, and a cubic ideal orientation has 24 symmetry-equivalent balls, which do not overlap below 45 degrees. The example evaluates the random reference the texture analysis prints beside every component.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np

from pytex import FrameDomain, Handedness, Lattice, Phase, ReferenceFrame, SymmetrySpec
from pytex.core.lattice import CrystalPlane
from pytex.texture import PoleFigure, euler_section_ranges, impose_sample_symmetry
from pytex.texture.components import random_component_fraction

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
SPECIMEN = ReferenceFrame(
    "sample_rd_td_nd", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
)
HEXAGONAL = SymmetrySpec.from_point_group("6/mmm", reference_frame=CRYSTAL)
CUBIC = SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL)
LATTICE = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=CRYSTAL)
ZIRCONIUM = Phase("alpha_zr", lattice=LATTICE, symmetry=HEXAGONAL, crystal_frame=CRYSTAL)
BASAL = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)
```

:::

**Compute**

```python
result = random_component_fraction(15.0, CUBIC.order)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `random-texture-fraction-within-15-degrees-cubic` | 0.022768 | 0.022768 | &mdash; | 1.41e-07 | 1e-06 | ✅ pass |

**Why this value**: Analytic: w = 15 deg = 0.2617993878 rad and sin 15 deg = (sqrt(6) - sqrt(2))/4 = 0.2588190451, so w - sin w = 0.0029803427; times |G| = 24 and divided by pi this is 0.0227680.

**Citation**: A. Morawiec, Orientations and Rotations: Computations in Crystallographic Textures, Springer (2004), chapter 2: the invariant (Haar) measure on SO(3), under which the rotation angle w of a uniform rotation has density (1 - cos w)/pi.

**See also**: {doc}`Texture analysis in the workbench <../../workflows/texture_analysis_workbench>`

## A hexagonal ODF with orthorhombic sample symmetry spans 90, 90 and 60 degrees

Before drawing ODF sections, decide how much of Euler space they must cover. Too little hides part of the texture; too much repeats it. For a hexagonal crystal the six-fold axis along c repeats every orientation after 60 degrees of phi2, and a rolled sheet's orthorhombic symmetry, with the two-fold axes perpendicular to c, folds phi1 and Phi into 90 degrees. The example returns the box (phi1, Phi, phi2) the texture analysis derives from the operators.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np

from pytex import FrameDomain, Handedness, Lattice, Phase, ReferenceFrame, SymmetrySpec
from pytex.core.lattice import CrystalPlane
from pytex.texture import PoleFigure, euler_section_ranges, impose_sample_symmetry
from pytex.texture.components import random_component_fraction

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
SPECIMEN = ReferenceFrame(
    "sample_rd_td_nd", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
)
HEXAGONAL = SymmetrySpec.from_point_group("6/mmm", reference_frame=CRYSTAL)
CUBIC = SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL)
LATTICE = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=CRYSTAL)
ZIRCONIUM = Phase("alpha_zr", lattice=LATTICE, symmetry=HEXAGONAL, crystal_frame=CRYSTAL)
BASAL = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)
```

:::

**Compute**

```python
ranges = euler_section_ranges(HEXAGONAL, "orthorhombic")
result = [ranges["phi1_max_deg"], ranges["big_phi_max_deg"], ranges["phi2_max_deg"]]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hexagonal-orthorhombic-odf-section-box` | [90.0, 90.0, 60.0] | [90.0, 90.0, 60.0] | deg | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: Bunge (1982) section 4.2: the asymmetric Euler region for crystal symmetry 6/mmm (proper group 622) with orthorhombic sample symmetry is 0 <= phi1 <= 90, 0 <= Phi <= 90, 0 <= phi2 <= 60 degrees.

**Citation**: H.-J. Bunge, Texture Analysis in Materials Science (1982), section 4.2 and its table of asymmetric units of Euler space.

**See also**: {doc}`Texture analysis in the workbench <../../workflows/texture_analysis_workbench>`
