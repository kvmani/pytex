<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Diffraction geometry

Powder scattering angles from PyTex interplanar spacings via Bragg's law, Kikuchi band and zone-axis geometry in the gnomonic projection, zone-axis routing on a stereographic Kikuchi map, the EBSD camera geometry, and preferred-orientation corrections to powder intensities — each checked against a standard reference value or a closed-form identity.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Ni(111) powder reflection angle for Cu K-alpha1

You are calibrating or interpreting a powder pattern and need to predict where the Ni(111) peak should appear with a copper source. PyTex supplies the interplanar spacing from the lattice metric; Bragg's law then gives the scattering angle. The result should land on the textbook Ni(111) position near 44.5 degrees for Cu K-alpha1.

**Symbols**

- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.
- $\lambda$ &mdash; Radiation wavelength.
- $\theta$ &mdash; Bragg half-angle.
- $2\theta$ &mdash; Powder-diffraction scattering angle.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom
```

:::

**Compute**

```python
d_111 = MillerPlane.from_hkl([1, 1, 1], phase=nickel).d_spacing_angstrom
theta = np.arcsin(cu_ka1 / (2.0 * d_111))
result = float(np.degrees(2.0 * theta))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-ni-111-two-theta` | 44.496 | 44.496 | deg | 1.43e-04 | 5e-03 | ✅ pass |

**Why this value**: d_111 = 3.52387 / sqrt(3) = 2.03451 angstrom; with lambda = 1.5406 angstrom, 2*theta = 2*arcsin(lambda / (2 d)) = 44.50 degrees, matching standard Ni powder data.

**Citation**: ICDD PDF 04-0850 (nickel); Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed.

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Ni{111} Kikuchi band width at 20 kV

You are reading an EBSD pattern and want to know how wide the strongest bands should be, either to check a detector calibration or to identify a phase from band widths alone. A Kikuchi band is bounded by the two Kossel cones of its lattice plane, so its angular width is exactly 2*theta_B and Bragg's law makes that a direct measurement of the interplanar spacing: wide bands mean large d-spacings. The widest band of nickel comes from {111}.

**Symbols**

- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.
- $\lambda$ &mdash; Radiation wavelength.
- $\theta$ &mdash; Bragg half-angle.
- $2\theta_B$ &mdash; Angular width of a Kikuchi band.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    DiffractionGeometry,
    GnomonicProjection,
    Orientation,
    simulate_kikuchi_pattern,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
detector = ReferenceFrame(
    name="detector",
    domain=FrameDomain.DETECTOR,
    axes=("u", "v", "n"),
    handedness=Handedness.RIGHT,
)
laboratory = ReferenceFrame(
    name="laboratory",
    domain=FrameDomain.LABORATORY,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
# A conventional 20 kV EBSD detector.
ebsd_geometry = DiffractionGeometry(
    detector_frame=detector,
    specimen_frame=specimen,
    laboratory_frame=laboratory,
    beam_energy_kev=20.0,
    camera_length_mm=15.0,
    pattern_center=np.array([0.5, 0.5, 0.6]),
    detector_pixel_size_um=(50.0, 50.0),
    detector_shape=(480, 640),
)
cube_orientation = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, phase=nickel
)
```

:::

**Compute**

```python
pattern = simulate_kikuchi_pattern(
    ebsd_geometry, nickel, cube_orientation, max_index=2
)
band = pattern.band_for_plane((1, 1, 1))
result = float(np.degrees(band.angular_width_rad))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-ni-111-kikuchi-band-width` | 2.4188 | 2.4187 | deg | 1.05e-04 | 2e-03 | ✅ pass |

**Why this value**: d_111 = 3.52387 / sqrt(3) = 2.034510 angstrom. The relativistic electron wavelength at 20 kV is 0.085883 angstrom, so sin(theta_B) = lambda / (2 d) = 0.0211066, giving theta_B = 1.20936 degrees and a band width of 2*theta_B = 2.4187 degrees. This is the familiar ~2.4 degree width of the strongest nickel bands in a 20 kV EBSD pattern.

**Citation**: Goldstein et al., Scanning Electron Microscopy and X-Ray Microanalysis, 4th ed. (electron wavelength table); Schwartz, Kumar, Adams and Field (eds.), Electron Backscatter Diffraction in Materials Science, 2nd ed.

**See also**: {doc}`EBSD foundation <../../concepts/ebsd_foundation>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Gnomonic radius of the [011] zone axis at the cube orientation

Zone axes are the landmarks of a Kikuchi pattern, and locating them is the first step in indexing one. The gnomonic projection places a direction at a radius equal to the tangent of its angle from the detector normal, so the geometry can be checked in closed form. With a cubic crystal at the cube orientation and an untilted detector, [011] lies exactly 45 degrees from the detector normal, and must therefore project to gnomonic radius tan(45 deg) = 1 exactly. This is an end-to-end check of the crystal to specimen to laboratory to detector chain: a transposed rotation anywhere along it moves the answer.

**Symbols**

- $r_g$ &mdash; Gnomonic radius, in units of the detector distance.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    DiffractionGeometry,
    GnomonicProjection,
    Orientation,
    simulate_kikuchi_pattern,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
detector = ReferenceFrame(
    name="detector",
    domain=FrameDomain.DETECTOR,
    axes=("u", "v", "n"),
    handedness=Handedness.RIGHT,
)
laboratory = ReferenceFrame(
    name="laboratory",
    domain=FrameDomain.LABORATORY,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
# A conventional 20 kV EBSD detector.
ebsd_geometry = DiffractionGeometry(
    detector_frame=detector,
    specimen_frame=specimen,
    laboratory_frame=laboratory,
    beam_energy_kev=20.0,
    camera_length_mm=15.0,
    pattern_center=np.array([0.5, 0.5, 0.6]),
    detector_pixel_size_um=(50.0, 50.0),
    detector_shape=(480, 640),
)
cube_orientation = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, phase=nickel
)
```

:::

**Compute**

```python
pattern = simulate_kikuchi_pattern(
    ebsd_geometry, nickel, cube_orientation, max_index=2
)
axis = next(
    zone for zone in pattern.zone_axes
    if tuple(int(value) for value in zone.indices) == (0, 1, 1)
)
result = float(np.hypot(*axis.coordinates))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-gnomonic-zone-axis-radius` | 1.000000000000 | 1.000000000000 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: The gnomonic projection maps a direction at angle psi from the detector normal to radius tan(psi). In a cubic crystal [011] makes an angle of 45 degrees with [001]; at the cube orientation with an untilted detector [001] is the detector normal, so the radius is tan(45 deg) = 1 exactly. The identity is exact, so the tolerance is numerical only.

**Citation**: Snyder, Map Projections: A Working Manual, USGS Professional Paper 1395 (gnomonic projection); Randle and Engler, Introduction to Texture Analysis, 2nd ed.

**See also**: {doc}`EBSD foundation <../../concepts/ebsd_foundation>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Where the specimen normal falls on a 70-degree EBSD screen

An EBSD geometry is stated in the terms the microscope is configured in - stage tilt, camera elevation, pattern centre - and a sign error in any of them produces a pattern that still looks like a plausible band network. One number checks the whole convention. With the beam as the laboratory z axis, the stage tilting the specimen normal towards the camera by sigma, and the camera axis raised by epsilon above the plane perpendicular to the beam, the specimen normal makes an angle of 90 - (sigma - epsilon) with the camera axis. The gnomonic projection therefore places it at radius tan(90 - sigma + epsilon) from the pattern centre: at the standard 70 degrees with the camera unelevated, tan(20 deg). That the value is well under one is the reason the specimen normal falls on a real screen at all.

**Symbols**

- $r_g$ &mdash; Gnomonic radius, in units of the detector distance.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    DiffractionGeometry,
    GnomonicProjection,
    Orientation,
    simulate_kikuchi_pattern,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
detector = ReferenceFrame(
    name="detector",
    domain=FrameDomain.DETECTOR,
    axes=("u", "v", "n"),
    handedness=Handedness.RIGHT,
)
laboratory = ReferenceFrame(
    name="laboratory",
    domain=FrameDomain.LABORATORY,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
# A conventional 20 kV EBSD detector.
ebsd_geometry = DiffractionGeometry(
    detector_frame=detector,
    specimen_frame=specimen,
    laboratory_frame=laboratory,
    beam_energy_kev=20.0,
    camera_length_mm=15.0,
    pattern_center=np.array([0.5, 0.5, 0.6]),
    detector_pixel_size_um=(50.0, 50.0),
    detector_shape=(480, 640),
)
cube_orientation = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, phase=nickel
)
```

:::

**Compute**

```python
geometry = DiffractionGeometry.for_ebsd(
    sample_tilt_deg=70.0, detector_elevation_deg=0.0
)
projection = GnomonicProjection(geometry=geometry)
normal_lab = geometry.specimen_vectors_to_lab(
    np.array([[0.0, 0.0, 1.0]])
)
coordinates, _ = projection.project_directions(normal_lab)
result = float(np.hypot(*coordinates[0]))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-ebsd-specimen-normal-radius` | 0.363970234266 | 0.363970234266 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: tan(90 deg - 70 deg + 0 deg) = tan(20 deg) = 0.36397023426620234, from the stated laboratory frame alone. The identity is exact, so the tolerance is numerical only.

**Citation**: Schwartz, Kumar, Adams and Field (eds.), Electron Backscatter Diffraction in Materials Science, 2nd ed.; Britton et al., Materials Characterization 117 (2016) 113, doi:10.1016/j.matchar.2016.04.008 (EBSD frame conventions).

**See also**: {doc}`EBSD foundation <../../concepts/ebsd_foundation>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## March-Dollase factor for cubic {111} under a (111) plate texture

You are refining a powder pattern from a specimen that will not pack randomly — a platy powder, or a rolled foil — and the measured {111} peak is far too strong. The March-Dollase model absorbs that into one parameter. Because every symmetry-equivalent plane of a family diffracts at the same angle, the factor is averaged over the whole family, and for cubic {111} with the preferred axis along (111) that average has a closed form: one member sits at 0 degrees to the axis and three at arccos(1/3).

**Symbols**

- $r$ &mdash; March coefficient of the preferred-orientation model.
- $P_{hkl}$ &mdash; Preferred-orientation intensity factor, in multiples of random.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    MarchDollaseModel,
    MillerPlane,
    ODF,
    ODFPreferredOrientationModel,
    OrientationSet,
    march_dollase_factors,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
```

:::

**Compute**

```python
model = MarchDollaseModel(
    preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=nickel),
    march_coefficient=0.5,
)
result = float(model.factors([MillerPlane.from_hkl([1, 1, 1], phase=nickel)])[0])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-march-dollase-family-factor` | 2.309132723130 | 2.309132723130 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: With r = 1/2 the March function is P(a) = (r^2 cos^2 a + sin^2 a / r)^(-3/2). The cubic {111} family has four members up to inversion: one at a = 0, giving P = r^-3 = 8, and three at cos^2 a = 1/9, where the bracket is 1/36 + 16/9 = 65/36 and P = (36/65)^(3/2) = 216 / (65 sqrt(65)). The family mean is therefore (8 + 3 * 216 / (65 sqrt(65))) / 4 = 2 + 162 / (65 sqrt(65)) = 2.3091327231300272, evaluated in exact decimal arithmetic independently of PyTex.

**Citation**: Dollase, W. A., J. Appl. Cryst. 19, 267-272 (1986), DOI: 10.1107/S0021889886089458; March, A., Z. Kristallogr. 81, 285-297 (1932).

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Texture foundation <../../concepts/texture_foundation>`

## The March distribution integrates to one over the sphere

Before trusting any preferred-orientation correction it is worth confirming that it *redistributes* diffracted intensity rather than inventing it. The March distribution has exactly that property: averaged over a uniform distribution of directions it is 1 for every March coefficient. This is the statement that makes a fitted r a description of texture rather than a free intensity scale.

**Symbols**

- $r$ &mdash; March coefficient of the preferred-orientation model.
- $P_{hkl}$ &mdash; Preferred-orientation intensity factor, in multiples of random.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    MarchDollaseModel,
    MillerPlane,
    ODF,
    ODFPreferredOrientationModel,
    OrientationSet,
    march_dollase_factors,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
```

:::

**Compute**

```python
u = np.linspace(-1.0, 1.0, 2_000_001)
factors = march_dollase_factors(np.arccos(u), 0.4)
result = float(np.trapezoid(factors, u) / 2.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-march-dollase-normalization` | 1.000000000 | 1.000000000 | &mdash; | 5.71e-11 | 1e-06 | ✅ pass |

**Why this value**: Substituting u = cos(a), the spherical average is the integral over u in [-1, 1] of ((r^2 - 1/r) u^2 + 1/r)^(-3/2), halved. Its antiderivative is u / (B sqrt(A u^2 + B)) with A = r^2 - 1/r and B = 1/r, so the average evaluates to 1 / (B sqrt(A + B)) = r / r = 1 for every positive r. The identity is exact; the tolerance is the quadrature error alone.

**Citation**: Dollase, W. A., J. Appl. Cryst. 19, 267-272 (1986), DOI: 10.1107/S0021889886089458.

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Texture foundation <../../concepts/texture_foundation>`

## ODF-weighted intensities reduce to the random powder

PyTex can drive powder intensities from a measured orientation distribution instead of a fitted parameter: the intensity of a reflection scales with the pole density along the scattering vector, and an ODF supplies exactly that. The check that makes the result interpretable is the limiting case — an untextured specimen must reproduce the random powder the uncorrected pattern already assumes, giving a factor of 1 for every reflection with no fitted parameter anywhere.

**Symbols**

- $P_{hkl}$ &mdash; Preferred-orientation intensity factor, in multiples of random.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    MarchDollaseModel,
    MillerPlane,
    ODF,
    ODFPreferredOrientationModel,
    OrientationSet,
    march_dollase_factors,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
```

:::

**Compute**

```python
grid = OrientationSet.from_equispaced_so3_grid(
    12.0,
    specimen_frame=specimen,
    phase=nickel,
    reduce_to_fundamental_region=False,
)
model = ODFPreferredOrientationModel(odf=ODF.from_orientations(grid))
result = float(model.factors([MillerPlane.from_hkl([1, 1, 1], phase=nickel)])[0])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-odf-weighted-random-texture` | 0.9999 | 1.0000 | &mdash; | 1.44e-04 | 1e-02 | ✅ pass |

**Why this value**: Pole density is defined in multiples of a random distribution, so a uniform orientation distribution has pole density 1 along every specimen direction and the correction is the identity. The tolerance reflects the finite SO(3) grid used to represent the uniform distribution, not any approximation in the correction itself.

**Citation**: Bunge, H.-J., Texture Analysis in Materials Science, DOI: 10.1016/C2013-0-11769-2; Von Dreele, R. B., J. Appl. Cryst. 30, 517-525 (1997), DOI: 10.1107/S0021889897005918 (texture in Rietveld refinement).

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Texture foundation <../../concepts/texture_foundation>`

## Measured powder-profile comparison recovers a known scale and background

You have imported a measured powder profile and want an auditable first comparison with a simulated profile before attempting any structural refinement. This deliberately synthetic validation case sets I_obs = 5 I_sim + 5 at five points with equal standard uncertainty. Weighted least squares must therefore recover scale 5 and background 5 exactly, while both IUCr profile residuals vanish.

**Symbols**

- $R_p$ &mdash; Unweighted whole-profile agreement factor.
- $R_{wp}$ &mdash; Weighted whole-profile agreement factor.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom

from pytex import (
    MeasuredPowderPattern,
    PowderPattern,
    compare_powder_patterns,
)

axis = np.arange(20.0, 25.0)
simulated_intensity = np.arange(1.0, 6.0)
simulated = PowderPattern(
    phase=nickel,
    radiation=RadiationSpec.cu_ka(),
    reflections=(),
    two_theta_grid_deg=axis,
    intensity_grid=simulated_intensity,
)
measured = MeasuredPowderPattern(
    name="synthetic affine validation profile",
    two_theta_deg=axis,
    intensity=5.0 * simulated_intensity + 5.0,
    standard_uncertainty=np.ones(5),
    intensity_unit="counts",
    radiation=RadiationSpec.cu_ka(),
    synthetic=True,
    metadata={"fixture_kind": "synthetic_validation"},
)
```

:::

**Compute**

```python
comparison = compare_powder_patterns(measured, simulated)
result = np.array([
    comparison.scale_factor,
    comparison.background_offset,
    comparison.profile_r_factor,
    comparison.weighted_profile_r_factor,
])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-powder-profile-affine-comparison` | [5.000000000000, 5.000000000000, 0.000000000000, 0.000000000000] | [5.000000000000, 5.000000000000, 0.000000000000, 0.000000000000] | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: The five observed values are constructed independently as 5*x + 5 from x = 1,...,5. The weighted design matrix therefore contains the exact affine solution (5, 5), every residual is zero, and the numerators of both R_p and R_wp are exactly zero.

**Citation**: IUCr pdCIF dictionary definitions _pd_proc_ls.prof_R_factor and _pd_proc_ls.prof_wR_factor; Young, The Rietveld Method (1993), Ch. 1.

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Powder XRD and SAED theory <../../theory/powder_xrd_and_saed>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Kikuchi-map routing reproduces the exact cubic zone-axis angles

Build the stereographic Kikuchi map of nickel and ask it for the tilt from [001] to [011], to [111], and to [112], with a leg budget large enough that each is a single hop along one band. The angles between low-index cubic directions are closed-form - 45 degrees, arccos(1/sqrt(3)), and arccos(2/sqrt(6)) - so the routed travel is checked against arithmetic rather than against a prior run. Getting these right exercises the whole chain: the direct basis, the map frame, the zone law that decides which bands join two axes, and the shortest-path search.

**Symbols**

- $\theta$ &mdash; Angle between two crystal zone axes; the stage travel of one routing leg.


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
    compute_kikuchi_map,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
kikuchi_map = compute_kikuchi_map(
    nickel,
    beam_energy_kev=200.0,
    max_index=4,
    zone_axis_max_index=3,
)
targets = ([0, 1, 1], [1, 1, 1], [1, 1, 2])
result = np.array(
    [
        kikuchi_map.route_to([0, 0, 1], target, max_leg_deg=90.0).total_tilt_deg
        for target in targets
    ]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-kikuchi-map-zone-axis-tilt-angles` | [45.000000, 54.735610, 35.264390] | [45.000000, 54.735610, 35.264390] | deg | < 1e-12 | 1e-06 | ✅ pass |

**Why this value**: Closed-form angles between cubic directions: arccos(1/sqrt(2)) = 45 deg for [001]-[011], arccos(1/sqrt(3)) = 54.735610 deg for [001]-[111], and arccos(2/sqrt(6)) = 35.264390 deg for [001]-[112]. The two [111] and [112] values are complementary, summing to 90 degrees, because [112] is the reflection of [001] in the plane perpendicular to [111].

**Citation**: Standard cubic interaxial angles; see International Tables for Crystallography Vol. C (1999) for the reciprocal-lattice conventions, and Williams and Carter, Transmission Electron Microscopy 2nd ed. (2009) Ch. 19 for Kikuchi-map tilting.

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Texture foundation <../../concepts/texture_foundation>`
