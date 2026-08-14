<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Simulated SAED plates and the zone-axis atlas

The geometry a practice diffraction pattern must reproduce if indexing it is to teach anything: the camera-constant identity that places every reflection, the hcp prism-zone aspect ratio that measures c/a without any calibration at all, and the basal-to-prism angle the zone-axis atlas has to report as exactly 90 degrees.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Where the 200 reflection of aluminium lands on the detector

The single identity a diffraction pattern is calibrated by: a reflection sits at a distance from the transmitted beam equal to the camera constant divided by its d-spacing. Everything downstream — the indexed answer, the lattice parameter, the phase identification — inherits this one relation, which is why a camera constant taken from the wrong camera length produces a self-consistent pattern of the wrong material. Here the simulated plate is asked where its strongest reflection is, and the answer must be the one the definition gives.

**Symbols**

- $d$ &mdash; Interplanar spacing of a reflecting plane.
- $g$ &mdash; Reciprocal-lattice vector of a reflection; |g| = 1/d.
- $\lambda$ &mdash; Electron wavelength at the accelerating voltage.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Aluminium, a = 4.0495 A (Wyckoff, Crystal Structures Vol. 1).
aluminium = Phase(
    "aluminium-fcc",
    lattice=Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# Alpha zirconium, a = 3.232 A, c = 5.147 A, so c/a = 1.5925.
zirconium = Phase(
    "zirconium-hcp",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
# A 400 mm camera length at 200 kV, where lambda = 0.0250793 A: L*lambda rounded
# to 10.0317 mm.A. The camera constant is an input here, not the quantity under
# test, so it is written out rather than recomputed.
CAMERA_CONSTANT = 10.0317
RASTER = DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=0.024)
```

:::

**Compute**

```python
image = synthesize_saed_image(
    aluminium,
    ZoneAxis([0, 0, 1], phase=aluminium),
    camera_constant_mm_angstrom=CAMERA_CONSTANT,
    raster=RASTER,
)
spot = next(
    entry
    for entry in image.spots
    if sorted(abs(int(value)) for value in entry.miller_indices) == [0, 0, 2]
)
centre = np.asarray(image.centre_px)
radius_px = float(np.linalg.norm(np.asarray(spot.position_px) - centre))
result = radius_px * image.raster.pixel_size_mm
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-practice-camera-constant-identity` | 4.95454 | 4.95454 | mm | 0.00e+00 | 1e-09 | ✅ pass |

**Why this value**: r = (L*lambda) / d with d_200 = a/2 = 2.02475 A for a = 4.0495 A, giving r = 10.0317 / 2.02475 = 4.95454 mm. Analytic from the lattice parameter and the definition of the camera constant; no program output enters it.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., Springer, DOI: 10.1007/978-0-387-76501-3, chapter 18 (the camera equation R d = L lambda); Wyckoff, R. W. G., Crystal Structures Vol. 1 (1963) for a.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`

## Reading c/a off one hcp prism-zone pattern

The hcp [2-1-10] pattern is a rectangle whose two shortest vectors are 0002 along c* and 01-10 perpendicular to it. Their lengths are 2/c and 2/(sqrt(3) a), so their ratio is sqrt(3) a / c and depends on nothing else — not on the camera constant, not on the accelerating voltage, not on the exposure. That makes this one pattern a calibration-free measurement of the axial ratio, and the standard way to tell zirconium (1.0876) from titanium (1.0908) or magnesium (1.0668) on the microscope. Here the ratio is measured off the simulated plate exactly as it would be measured off a real one.

**Symbols**

- $d$ &mdash; Interplanar spacing of a reflecting plane.
- $g$ &mdash; Reciprocal-lattice vector of a reflection; |g| = 1/d.
- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Aluminium, a = 4.0495 A (Wyckoff, Crystal Structures Vol. 1).
aluminium = Phase(
    "aluminium-fcc",
    lattice=Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# Alpha zirconium, a = 3.232 A, c = 5.147 A, so c/a = 1.5925.
zirconium = Phase(
    "zirconium-hcp",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
# A 400 mm camera length at 200 kV, where lambda = 0.0250793 A: L*lambda rounded
# to 10.0317 mm.A. The camera constant is an input here, not the quantity under
# test, so it is written out rather than recomputed.
CAMERA_CONSTANT = 10.0317
RASTER = DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=0.024)
```

:::

**Compute**

```python
image = synthesize_saed_image(
    zirconium,
    ZoneAxis([1, 0, 0], phase=zirconium),
    camera_constant_mm_angstrom=CAMERA_CONSTANT,
    raster=RASTER,
)
by_indices = {
    tuple(int(value) for value in spot.miller_indices): spot for spot in image.spots
}
result = by_indices[(0, 0, 2)].g_inv_angstrom / by_indices[(0, 1, 0)].g_inv_angstrom
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-practice-hcp-prism-axial-ratio` | 1.08762 | 1.08762 | &mdash; | 0.00e+00 | 1e-09 | ✅ pass |

**Why this value**: |g_0002| / |g_01-10| = (2/c) / (2 / (sqrt(3) a)) = sqrt(3) a / c = sqrt(3) * 3.232 / 5.147 = 1.08762, from the hexagonal reciprocal metric alone. Independent of the camera constant, which cancels in the ratio.

**Citation**: Edington, J. W., Practical Electron Microscopy in Materials Science, Macmillan (1975), on interpreting hexagonal zone-axis patterns; lattice parameters from Wyckoff, R. W. G., Crystal Structures Vol. 1 (1963).

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`

## Basal to prism is 90 degrees for every hexagonal metal

The zone-axis atlas exists to answer 'where should I go next', and the first thing it has to get right is how far away each candidate is. The hexagonal basal-to-prism pair is the cleanest possible check: the c axis is perpendicular to every a axis by the definition of the hexagonal cell, so the angle is exactly 90 degrees whatever the axial ratio, whatever the metal. It is also the practical lesson the pair teaches, because 90 degrees is beyond any conventional double-tilt holder in one move.

**Symbols**

- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Aluminium, a = 4.0495 A (Wyckoff, Crystal Structures Vol. 1).
aluminium = Phase(
    "aluminium-fcc",
    lattice=Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# Alpha zirconium, a = 3.232 A, c = 5.147 A, so c/a = 1.5925.
zirconium = Phase(
    "zirconium-hcp",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
# A 400 mm camera length at 200 kV, where lambda = 0.0250793 A: L*lambda rounded
# to 10.0317 mm.A. The camera constant is an input here, not the quantity under
# test, so it is written out rather than recomputed.
CAMERA_CONSTANT = 10.0317
RASTER = DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=0.024)
```

:::

**Compute**

```python
from pytex.tem.atlas import zone_axis_atlas

atlas = zone_axis_atlas(
    zirconium,
    current_zone_axis=ZoneAxis([0, 0, 1], phase=zirconium),
    max_index=1,
)
prism = next(entry for entry in atlas.entries if entry.label == '[100]')
result = prism.angle_from_current_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-practice-atlas-basal-to-prism` | 90.000000 | 90.000000 | deg | 0.00e+00 | 1e-09 | ✅ pass |

**Why this value**: In a hexagonal lattice alpha = beta = 90 degrees, so c is orthogonal to both a axes by construction of the cell. The angle between [0001] and [2-1-10] is therefore exactly 90 degrees for any c/a.

**Citation**: International Tables for Crystallography, Volume A, on the hexagonal cell setting; Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., DOI: 10.1007/978-0-387-76501-3, chapter 18.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`TEM specimen tilt navigation <../../theory/tem_specimen_tilt_navigation>`
