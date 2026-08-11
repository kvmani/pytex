# Kikuchi Band Geometry

PyTex models Kikuchi bands and the gnomonic projection as one geometric layer shared by EBSD and
TEM. This is the surface that answers *where* the bands of a known phase and orientation fall on a
given detector, and *how wide* they are.

## Scope

- gnomonic projection between laboratory directions, gnomonic coordinates, and detector pixels
- Kikuchi bands as pairs of Kossel-cone traces, with exact conic edges
- band angular width from Bragg's law, and the distinct *apparent* width in the projection
- zone axes as the points where the bands sharing them intersect
- plotting in either gnomonic or detector coordinates

## Why Gnomonic Coordinates

A Kikuchi band is not a spot. Electrons scattered inside the specimen travel in every direction;
those meeting a lattice plane at the Bragg angle diffract, and the diffracting directions form two
cones — the Kossel cones — of semi-angle $90^\circ - \theta_B$ about the plane normal. Their
intersections with the detector are the two edges of the band, and the plane's own trace runs
midway between them.

The gnomonic projection is central projection from the diffraction source onto the detector plane,
in units of the detector distance. Its defining property is that **a great circle maps to a
straight line**. A lattice-plane trace is a great circle, so:

> Band centre lines are exactly straight in gnomonic coordinates, whatever the detector tilt.

In raw detector pixels the same trace curves as soon as the detector is tilted. That is why band
detection, pattern-centre fitting, and orientation refinement are all done in gnomonic coordinates.

Band *edges*, by contrast, are **not** straight: the Kossel cones are not great circles, so the
edges are conics — hyperbolae at the small Bragg angles of electron diffraction. The common
textbook statement that the edges are straight and parallel to the centre is the small-angle
approximation; PyTex does not make it, and samples the exact cones instead.

## Band Width Measures The Lattice

The angular width of a band is exactly $2\theta_B$, and

$$\sin\theta_B = \frac{\lambda}{2d},$$

so a wide band means a large interplanar spacing. This is what makes band widths usable for phase
discrimination. For nickel at 20 kV the strongest bands, from $\{111\}$ with $d = 2.0345$ Å, are
about 2.42 degrees wide — a value the worked example
`diffraction-ni-111-kikuchi-band-width` computes live and checks against a hand derivation.

The *apparent* width in gnomonic units is a different quantity: the projection stretches with
distance from the pattern centre, so a band far from the centre looks wider than one near it at the
same $\theta_B$. `KikuchiBand.width_at_pattern_center` reports that apparent width separately, so
that converting a measured width to a $d$-spacing without accounting for position is a visible
choice rather than a silent error.

## Example

```python
import numpy as np

from pytex import (
    DiffractionGeometry,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    plot_kikuchi_pattern,
    simulate_kikuchi_pattern,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT)
detector = ReferenceFrame("detector", FrameDomain.DETECTOR, ("u", "v", "n"), Handedness.RIGHT)
laboratory = ReferenceFrame("laboratory", FrameDomain.LABORATORY, ("x", "y", "z"), Handedness.RIGHT)

# Nickel. Declaring the space group matters: without it the centring absences
# cannot be applied and forbidden reflections would be listed as bands.
phase = Phase(
    name="nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    space_group_symbol="Fm-3m",
    space_group_number=225,
)

geometry = DiffractionGeometry(
    detector_frame=detector,
    specimen_frame=specimen,
    laboratory_frame=laboratory,
    beam_energy_kev=20.0,
    camera_length_mm=15.0,
    pattern_center=np.array([0.5, 0.5, 0.6]),
    detector_pixel_size_um=(50.0, 50.0),
    detector_shape=(480, 640),
)
orientation = Orientation.from_euler(0.0, 0.0, 0.0, specimen_frame=specimen, phase=phase)

pattern = simulate_kikuchi_pattern(geometry, phase, orientation, max_index=2)
print(pattern.describe())

figure = plot_kikuchi_pattern(pattern, coordinates="gnomonic", max_bands=12)
figure.savefig("ni_fcc_kikuchi.png", dpi=200)
```

`describe()` states the phase, the orientation in Bunge Euler angles, the beam energy and
wavelength, the widest bands with their spacings and angular widths, and the kinematic-intensity
limitation — so a reader who did not run the code can still audit what was computed.

## Coordinate Semantics

Three coordinate meanings stay separate, and conversions between them are explicit:

| Coordinates | Origin | Unit | Use |
| --- | --- | --- | --- |
| laboratory directions | source | unit vectors | the physical rays |
| gnomonic | pattern centre | detector distances | geometry; centre lines are straight |
| detector pixels | image corner | pixels | what a camera records |

`GnomonicProjection` converts between all three. Because gnomonic coordinates are in units of the
detector distance, they are detector-size independent: the same crystal gives the same gnomonic
pattern on any detector, which is what makes them the right frame for comparing or fitting
patterns.

Directions travelling away from the detector have no intersection with it. They are reported as
`NaN` with a `False` validity flag rather than being given a coordinate, so an invalid projection
can never be mistaken for a position near the origin.

## Zone Axes

A direction $[uvw]$ lies in a plane $(hkl)$ when $hu + kv + lw = 0$. Every band whose plane
contains a zone axis therefore passes through that axis's projected point, so zone axes appear as
hubs where several bands meet — the features by which patterns are recognized. PyTex reports only
axes shared by at least two bands, since one band does not define a hub, and flags whether each
projects onto the physical detector.

The `[011]` axis of a cubic crystal at the cube orientation lies exactly 45 degrees from an
untilted detector normal, and so must project to gnomonic radius $\tan(45^\circ) = 1$ exactly. The
worked example `diffraction-gnomonic-zone-axis-radius` checks this, which pins the whole
crystal → specimen → laboratory → detector chain in closed form.

## Current Limits

- band positions and widths are exact for the stated geometry; **intensities are a kinematic
  $|F|^2$ proxy**
- the excess/deficiency asymmetry across a band comes from the dynamical theory and is not
  modelled, so simulated bands are symmetric
- higher-order Laue zone rings, background, and detector point-spread are not modelled
- orientation determination uses band *positions*, so the intensity limitation does not affect the
  geometric use of this surface

## Related Material

- {doc}`../concepts/ebsd_foundation`
- {doc}`../concepts/diffraction_foundation`
- {doc}`diffraction_geometry`
- {doc}`saed_generation`
- {doc}`orix_kikuchipy_interop`
- {doc}`/theory/kikuchi_bands_and_gnomonic_projection`

## References

### Normative

- `../../standards/reference_canon.md`
- `../../standards/notation_and_conventions.md`

### Informative

- `../../testing/diffraction_validation_matrix.md`
