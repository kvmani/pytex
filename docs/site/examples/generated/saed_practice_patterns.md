<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Simulated SAED plates and the zone-axis atlas

The geometry a practice diffraction pattern must reproduce if indexing it is to teach anything: the camera-constant identity that places every reflection, the hcp prism-zone aspect ratio that measures c/a without any calibration at all, and the basal-to-prism angle the zone-axis atlas has to report as exactly 90 degrees, the beam centre a lattice fit recovers from the spots, and the length bias a mis-set camera constant leaves in the scoring while the angles stay put, and the forbidden reflection that double diffraction puts on a real plate at exactly the radius a genuine one would occupy.

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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
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
| `saed-practice-camera-constant-identity` | 4.95454 | 4.95454 | mm | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: r = (L*lambda) / d with d_200 = a/2 = 2.02475 A for a = 4.0495 A, giving r = 10.0317 / 2.02475 = 4.95454 mm. Analytic from the lattice parameter and the definition of the camera constant; no program output enters it.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., Springer, DOI: 10.1007/978-0-387-76501-3, chapter 18 (the camera equation R d = L lambda); Wyckoff, R. W. G., Crystal Structures Vol. 1 (1963) for a.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`

## What rolling the crystal about the beam does to the pattern

A simulated plate states the orientation it was built from, as the rotation taking crystal vectors into the pattern frame. Anything else drawn on that pattern - Kikuchi bands, a calculated overlay, a stereogram - is placed with that matrix, so it has to mean exactly what the spots mean or every overlay is silently turned. The identity that tests it needs no reference table: rolling the crystal about the beam by an angle rotates the azimuth of every reflection on the plate by that same angle, because the roll is a rotation about the projection axis and the projection commutes with it. Here the 200 reflection of aluminium is projected through the matrix at two rolls thirty degrees apart, and the angle between the two positions is measured.

**Symbols**

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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
```

:::

**Compute**

```python
def projected(roll_deg):
    image = synthesize_saed_image(
        aluminium,
        ZoneAxis([0, 0, 1], phase=aluminium),
        camera_constant_mm_angstrom=CAMERA_CONSTANT,
        raster=RASTER,
        in_plane_rotation_deg=roll_deg,
    )
    reciprocal = aluminium.lattice.reciprocal_basis().matrix
    g_200 = reciprocal @ np.array([2.0, 0.0, 0.0])
    return image.crystal_to_pattern() @ g_200


start = projected(0.0)
rolled = projected(30.0)
result = float(
    np.degrees(
        np.arctan2(
            start[0] * rolled[1] - start[1] * rolled[0],
            start[0] * rolled[0] + start[1] * rolled[1],
        )
    )
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-practice-roll-about-the-beam` | 30.000000 | 30.000000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: Exact by construction rather than by measurement: the roll is a rotation about the beam, which is the projection axis, so it acts on the detector plane as a plane rotation through the same angle. Thirty degrees of roll must move every spot's azimuth by thirty degrees, whatever the lattice, the camera constant or the reflection.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., Springer, DOI: 10.1007/978-0-387-76501-3, chapter 18 - the pattern rotates rigidly with the specimen about the beam, which is why one pattern cannot fix that rotation.

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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
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
| `saed-practice-hcp-prism-axial-ratio` | 1.08762 | 1.08762 | &mdash; | < 1e-11 | 1e-09 | ✅ pass |

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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
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
| `saed-practice-atlas-basal-to-prism` | 90.000000 | 90.000000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: In a hexagonal lattice alpha = beta = 90 degrees, so c is orthogonal to both a axes by construction of the cell. The angle between [0001] and [2-1-10] is therefore exactly 90 degrees for any c/a.

**Citation**: International Tables for Crystallography, Volume A, on the hexagonal cell setting; Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., DOI: 10.1007/978-0-387-76501-3, chapter 18.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`TEM specimen tilt navigation <../../theory/tem_specimen_tilt_navigation>`

## A beam centre picked 30 pixels out, recovered from the spots

Picking the transmitted beam by eye is the largest avoidable error in the indexing workflow: it biases every d-spacing at once, and it does so while leaving the pattern self-consistent, so the result is a plausible answer for the wrong material rather than an obvious failure. But the spots of a zone-axis pattern lie on a plane lattice, and with four or more of them that constraint over-determines the centre. Here eight nodes of an exact square lattice are given with the centre deliberately misplaced by 30 pixels in each direction, and the fit is asked to put it back.

**Symbols**

- $g$ &mdash; Reciprocal-lattice vector of a reflection; |g| = 1/d.


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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
```

:::

**Compute**

```python
from pytex.diffraction.lattice_fit import fit_planar_lattice

basis = np.array([[100.0, 0.0], [0.0, 100.0]])
indices = np.array([[1, 0], [0, 1], [-1, 0], [0, -1],
                    [1, 1], [-1, -1], [2, 0], [0, 2]], dtype=float)
truth = np.array([512.0, 384.0])
nodes = truth + indices @ basis
fit = fit_planar_lattice(nodes, truth + np.array([30.0, 30.0]))
result = float(np.linalg.norm(fit.centre - truth))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-lattice-fit-recovers-the-beam-centre` | 0.000000000 | 0.000000000 | px | < 1e-08 | 1e-06 | ✅ pass |

**Why this value**: Exact. The eight points are exact nodes of the lattice about the true centre, so the least-squares problem for the centre with the integer assignment held fixed has that point as its exact solution: the residual is zero and the recovered centre is the generating one. Independent of the basis chosen and of the starting error, up to the half-spacing limit at which a fit would be relabelling which node the origin is.

**Citation**: Standard linear least squares on the lattice model p = c + m a + n b; Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., DOI: 10.1007/978-0-387-76501-3, chapter 18 on why the beam position governs every measured spacing.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`

## A camera constant five percent high, read back from the scoring

The one calibration error that does not announce itself. A camera constant taken from the wrong camera length rescales every measured spacing and leaves every measured angle untouched, so the pattern stays perfectly self-consistent while pointing at the wrong material. The scoring keeps lengths and angles apart for exactly this reason, and weights angles higher, because an angular disagreement is evidence about the crystallography while a length disagreement may only be evidence about the instrument.

**Symbols**

- $d$ &mdash; Interplanar spacing of a reflecting plane.
- $g$ &mdash; Reciprocal-lattice vector of a reflection; |g| = 1/d.


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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
```

:::

**Compute**

```python
from dataclasses import dataclass
from pytex.diffraction.solution_scoring import score_solution

@dataclass(frozen=True)
class Spot:
    measured_index: int
    hkl: tuple
    label: str
    predicted_g_inv_angstrom: tuple

@dataclass(frozen=True)
class Solution:
    solved_spots: tuple
    matched_fraction: float = 1.0

calculated = np.array([[0.5, 0.0], [0.0, 0.5], [0.5, 0.5], [1.0, 0.0]])
solution = Solution(tuple(
    Spot(index, (2, 0, 0), 'g', tuple(calculated[index]))
    for index in range(len(calculated))
))
score = score_solution(solution, 1.05 * calculated)
result = float(score.rms_relative_length_deviation)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-scoring-calibration-bias` | 0.0476190476 | 0.0476190476 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: d = 1/|g|, so measured g larger by a factor 1.05 makes every measured d smaller by 1/1.05. The relative deviation is 1/1.05 - 1 = -0.0476190476 on every spot, and the r.m.s. of a constant is that constant. Identical on every spot is the signature that distinguishes a calibration error from an indexing error.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., DOI: 10.1007/978-0-387-76501-3, chapter 18 (R d = L lambda).

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`

## A forbidden silicon 200 sits exactly where the lattice puts it

Silicon 200 is forbidden by the diamond glide, and a real [011] plate shows it anyway, produced by (1-1-1) + (111). The operationally important fact is not that the spot is there but *where* it is: at the same radius the camera equation gives for d = a/2, indistinguishable by position or by spacing from a genuine reflection. Nothing about the measurement reveals it, which is exactly why indexing a pattern on it silently yields the wrong cell — and why PyTex marks the reflection rather than leaving a reader to notice. Here the doubly diffracted spot is asked where it landed, and the answer must be the one the camera equation gives for a reflection that kinematic theory says is not there at all.

**Symbols**

- $d$ &mdash; Interplanar spacing of a reflecting plane.
- $g$ &mdash; Reciprocal-lattice vector of a reflection; |g| = 1/d.
- $\lambda$ &mdash; Electron wavelength at the accelerating voltage.
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
from pytex.core.lattice import AtomicSite, UnitCell
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

# Silicon, a = 5.43102 A (CODATA), with the full diamond-cubic basis. The motif
# is required here and not for the two phases above: a forbidden reflection is
# forbidden by its structure factor, and a phase carrying no atomic sites has no
# structure factor to vanish.
_DIAMOND_BASIS = (
    (0.0, 0.0, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.5, 0.5, 0.0),
    (0.25, 0.25, 0.25),
    (0.25, 0.75, 0.75),
    (0.75, 0.25, 0.75),
    (0.75, 0.75, 0.25),
)
_silicon_lattice = Lattice(
    5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0, crystal_frame=crystal
)
silicon = Phase(
    "silicon",
    lattice=_silicon_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=_silicon_lattice,
        sites=tuple(
            AtomicSite(label=f"Si{index + 1}", species="Si", fractional_coordinates=xyz)
            for index, xyz in enumerate(_DIAMOND_BASIS)
        ),
    ),
)
```

:::

**Compute**

```python
image = synthesize_saed_image(
    silicon,
    ZoneAxis([0, 1, -1], phase=silicon),
    camera_constant_mm_angstrom=CAMERA_CONSTANT,
    raster=RASTER,
    include_double_diffraction=True,
)
spot = next(
    entry
    for entry in image.spots
    if tuple(int(value) for value in entry.miller_indices) == (2, 0, 0)
)
# It is present only because the option was asked for, and it says so.
assert spot.is_double_diffraction
assert spot.double_diffraction_origin == '(1 -1 -1) + (111)'
centre = np.asarray(image.centre_px)
radius_px = float(np.linalg.norm(np.asarray(spot.position_px) - centre))
result = radius_px * image.raster.pixel_size_mm
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `saed-practice-double-diffraction-forbidden-200` | 3.69422 | 3.69422 | mm | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: r = (L*lambda) / d with d_200 = a/2 = 2.71551 A for a = 5.43102 A, giving r = 10.0317 / 2.71551 = 3.69420 mm. The same camera equation that places a genuine reflection, applied to a forbidden one: double diffraction changes which reflections are visible, never where they are. Analytic from the lattice parameter; no program output enters it.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., Springer, DOI: 10.1007/978-0-387-76501-3, chapter 16 (double diffraction; silicon 200 along [110] is the worked case) and chapter 18 (the camera equation R d = L lambda). CODATA lattice parameter of silicon.

**See also**: {doc}`TEM pattern indexing workflow <../../workflows/tem_pattern_indexing>`, {doc}`Ratio and angle indexing <../../theory/saed_ratio_angle_indexing>`
