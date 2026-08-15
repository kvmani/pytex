"""Worked examples: diffraction geometry from crystallographic spacings.

These examples connect the PyTex crystallographic core to a diffraction
observable. The interplanar spacing d_hkl is computed from the metric tensor by
PyTex, and Bragg's law lambda = 2 d sin(theta) then fixes the powder scattering
angle 2*theta. The reference is the well-known Ni(111) reflection position for
Cu K-alpha1 radiation, near 44.5 degrees.

See :doc:`../../workflows/xrd_generation` and the diffraction theory notes.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

NICKEL_SETUP = """
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
"""

_DSPACING = SymbolUse(r"d_{hkl}", "Interplanar spacing of the (hkl) family.")
_THETA = SymbolUse(r"\theta", "Bragg half-angle.")
_TWO_THETA = SymbolUse(r"2\theta", "Powder-diffraction scattering angle.")
_LAMBDA = SymbolUse(r"\lambda", "Radiation wavelength.")

_XRD_WORKFLOW = SeeAlso("Powder XRD generation", "../../workflows/xrd_generation")
_DIFF_CONCEPT = SeeAlso("Diffraction foundation", "../../concepts/diffraction_foundation")
_XRD_THEORY = SeeAlso("Powder XRD and SAED theory", "../../theory/powder_xrd_and_saed")


NI_111_TWO_THETA = WorkedExample(
    id="diffraction-ni-111-two-theta",
    title="Ni(111) powder reflection angle for Cu K-alpha1",
    domain="diffraction",
    scenario=(
        "You are calibrating or interpreting a powder pattern and need to predict where the Ni(111) "
        "peak should appear with a copper source. PyTex supplies the interplanar spacing from the "
        "lattice metric; Bragg's law then gives the scattering angle. The result should land on the "
        "textbook Ni(111) position near 44.5 degrees for Cu K-alpha1."
    ),
    setup=NICKEL_SETUP,
    code=(
        "d_111 = MillerPlane.from_hkl([1, 1, 1], phase=nickel).d_spacing_angstrom\n"
        "theta = np.arcsin(cu_ka1 / (2.0 * d_111))\n"
        "result = float(np.degrees(2.0 * theta))"
    ),
    expected=44.496,
    unit="deg",
    tolerance=5e-3,
    reference=(
        "d_111 = 3.52387 / sqrt(3) = 2.03451 angstrom; with lambda = 1.5406 angstrom, "
        "2*theta = 2*arcsin(lambda / (2 d)) = 44.50 degrees, matching standard Ni powder data."
    ),
    citation="ICDD PDF 04-0850 (nickel); Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed.",
    symbols=(_DSPACING, _LAMBDA, _THETA, _TWO_THETA),
    see_also=(_XRD_WORKFLOW, _DIFF_CONCEPT),
    result_format="{:.3f}",
)


KIKUCHI_SETUP = (
    NICKEL_SETUP
    + """
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
"""
)

_BAND_WIDTH = SymbolUse(r"2\theta_B", "Angular width of a Kikuchi band.")
_GNOMONIC = SymbolUse(r"r_g", "Gnomonic radius, in units of the detector distance.")

_EBSD_CONCEPT = SeeAlso("EBSD foundation", "../../concepts/ebsd_foundation")


NI_111_KIKUCHI_BAND_WIDTH = WorkedExample(
    id="diffraction-ni-111-kikuchi-band-width",
    title="Ni{111} Kikuchi band width at 20 kV",
    domain="diffraction",
    scenario=(
        "You are reading an EBSD pattern and want to know how wide the strongest bands should be, "
        "either to check a detector calibration or to identify a phase from band widths alone. A "
        "Kikuchi band is bounded by the two Kossel cones of its lattice plane, so its angular width "
        "is exactly 2*theta_B and Bragg's law makes that a direct measurement of the interplanar "
        "spacing: wide bands mean large d-spacings. The widest band of nickel comes from {111}."
    ),
    setup=KIKUCHI_SETUP,
    code=(
        "pattern = simulate_kikuchi_pattern(\n"
        "    ebsd_geometry, nickel, cube_orientation, max_index=2\n"
        ")\n"
        "band = pattern.band_for_plane((1, 1, 1))\n"
        "result = float(np.degrees(band.angular_width_rad))"
    ),
    expected=2.4187,
    unit="deg",
    tolerance=2e-3,
    reference=(
        "d_111 = 3.52387 / sqrt(3) = 2.034510 angstrom. The relativistic electron wavelength at "
        "20 kV is 0.085883 angstrom, so sin(theta_B) = lambda / (2 d) = 0.0211066, giving "
        "theta_B = 1.20936 degrees and a band width of 2*theta_B = 2.4187 degrees. This is the "
        "familiar ~2.4 degree width of the strongest nickel bands in a 20 kV EBSD pattern."
    ),
    citation=(
        "Goldstein et al., Scanning Electron Microscopy and X-Ray Microanalysis, 4th ed. "
        "(electron wavelength table); Schwartz, Kumar, Adams and Field (eds.), Electron "
        "Backscatter Diffraction in Materials Science, 2nd ed."
    ),
    symbols=(_DSPACING, _LAMBDA, _THETA, _BAND_WIDTH),
    see_also=(_EBSD_CONCEPT, _DIFF_CONCEPT),
    result_format="{:.4f}",
)


CUBE_ZONE_AXIS_GNOMONIC_RADIUS = WorkedExample(
    id="diffraction-gnomonic-zone-axis-radius",
    title="Gnomonic radius of the [011] zone axis at the cube orientation",
    domain="diffraction",
    scenario=(
        "Zone axes are the landmarks of a Kikuchi pattern, and locating them is the first step in "
        "indexing one. The gnomonic projection places a direction at a radius equal to the tangent "
        "of its angle from the detector normal, so the geometry can be checked in closed form. With "
        "a cubic crystal at the cube orientation and an untilted detector, [011] lies exactly 45 "
        "degrees from the detector normal, and must therefore project to gnomonic radius "
        "tan(45 deg) = 1 exactly. This is an end-to-end check of the crystal to specimen to "
        "laboratory to detector chain: a transposed rotation anywhere along it moves the answer."
    ),
    setup=KIKUCHI_SETUP,
    code=(
        "pattern = simulate_kikuchi_pattern(\n"
        "    ebsd_geometry, nickel, cube_orientation, max_index=2\n"
        ")\n"
        "axis = next(\n"
        "    zone for zone in pattern.zone_axes\n"
        "    if tuple(int(value) for value in zone.indices) == (0, 1, 1)\n"
        ")\n"
        "result = float(np.hypot(*axis.coordinates))"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "The gnomonic projection maps a direction at angle psi from the detector normal to radius "
        "tan(psi). In a cubic crystal [011] makes an angle of 45 degrees with [001]; at the cube "
        "orientation with an untilted detector [001] is the detector normal, so the radius is "
        "tan(45 deg) = 1 exactly. The identity is exact, so the tolerance is numerical only."
    ),
    citation=(
        "Snyder, Map Projections: A Working Manual, USGS Professional Paper 1395 (gnomonic "
        "projection); Randle and Engler, Introduction to Texture Analysis, 2nd ed."
    ),
    symbols=(_GNOMONIC,),
    see_also=(_EBSD_CONCEPT, _DIFF_CONCEPT),
    result_format="{:.12f}",
)


PREFERRED_ORIENTATION_SETUP = (
    NICKEL_SETUP
    + """
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
"""
)

_MARCH = SymbolUse("r", "March coefficient of the preferred-orientation model.")
_PO_FACTOR = SymbolUse(
    r"P_{hkl}", "Preferred-orientation intensity factor, in multiples of random."
)

_XRD_CONCEPT = SeeAlso("Texture foundation", "../../concepts/texture_foundation")


MARCH_DOLLASE_FAMILY_FACTOR = WorkedExample(
    id="diffraction-march-dollase-family-factor",
    title="March-Dollase factor for cubic {111} under a (111) plate texture",
    domain="diffraction",
    scenario=(
        "You are refining a powder pattern from a specimen that will not pack randomly — a "
        "platy powder, or a rolled foil — and the measured {111} peak is far too strong. The "
        "March-Dollase model absorbs that into one parameter. Because every symmetry-equivalent "
        "plane of a family diffracts at the same angle, the factor is averaged over the whole "
        "family, and for cubic {111} with the preferred axis along (111) that average has a "
        "closed form: one member sits at 0 degrees to the axis and three at arccos(1/3)."
    ),
    setup=PREFERRED_ORIENTATION_SETUP,
    code=(
        "model = MarchDollaseModel(\n"
        "    preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=nickel),\n"
        "    march_coefficient=0.5,\n"
        ")\n"
        "result = float(model.factors([MillerPlane.from_hkl([1, 1, 1], phase=nickel)])[0])"
    ),
    expected=2.3091327231300272,
    unit="",
    tolerance=1e-12,
    reference=(
        "With r = 1/2 the March function is P(a) = (r^2 cos^2 a + sin^2 a / r)^(-3/2). The cubic "
        "{111} family has four members up to inversion: one at a = 0, giving P = r^-3 = 8, and "
        "three at cos^2 a = 1/9, where the bracket is 1/36 + 16/9 = 65/36 and "
        "P = (36/65)^(3/2) = 216 / (65 sqrt(65)). The family mean is therefore "
        "(8 + 3 * 216 / (65 sqrt(65))) / 4 = 2 + 162 / (65 sqrt(65)) = 2.3091327231300272, "
        "evaluated in exact decimal arithmetic independently of PyTex."
    ),
    citation=(
        "Dollase, W. A., J. Appl. Cryst. 19, 267-272 (1986), "
        "DOI: 10.1107/S0021889886089458; March, A., Z. Kristallogr. 81, 285-297 (1932)."
    ),
    symbols=(_MARCH, _PO_FACTOR),
    see_also=(_XRD_WORKFLOW, _XRD_CONCEPT),
    result_format="{:.12f}",
)


MARCH_DOLLASE_NORMALIZATION = WorkedExample(
    id="diffraction-march-dollase-normalization",
    title="The March distribution integrates to one over the sphere",
    domain="diffraction",
    scenario=(
        "Before trusting any preferred-orientation correction it is worth confirming that it "
        "*redistributes* diffracted intensity rather than inventing it. The March distribution "
        "has exactly that property: averaged over a uniform distribution of directions it is 1 "
        "for every March coefficient. This is the statement that makes a fitted r a description "
        "of texture rather than a free intensity scale."
    ),
    setup=PREFERRED_ORIENTATION_SETUP,
    code=(
        "u = np.linspace(-1.0, 1.0, 2_000_001)\n"
        "factors = march_dollase_factors(np.arccos(u), 0.4)\n"
        "result = float(np.trapezoid(factors, u) / 2.0)"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-6,
    reference=(
        "Substituting u = cos(a), the spherical average is the integral over u in [-1, 1] of "
        "((r^2 - 1/r) u^2 + 1/r)^(-3/2), halved. Its antiderivative is "
        "u / (B sqrt(A u^2 + B)) with A = r^2 - 1/r and B = 1/r, so the average evaluates to "
        "1 / (B sqrt(A + B)) = r / r = 1 for every positive r. The identity is exact; the "
        "tolerance is the quadrature error alone."
    ),
    citation=(
        "Dollase, W. A., J. Appl. Cryst. 19, 267-272 (1986), "
        "DOI: 10.1107/S0021889886089458."
    ),
    symbols=(_MARCH, _PO_FACTOR),
    see_also=(_XRD_WORKFLOW, _XRD_CONCEPT),
    result_format="{:.9f}",
)


ODF_WEIGHTED_RANDOM_TEXTURE = WorkedExample(
    id="diffraction-odf-weighted-random-texture",
    title="ODF-weighted intensities reduce to the random powder",
    domain="diffraction",
    scenario=(
        "PyTex can drive powder intensities from a measured orientation distribution instead of "
        "a fitted parameter: the intensity of a reflection scales with the pole density along "
        "the scattering vector, and an ODF supplies exactly that. The check that makes the "
        "result interpretable is the limiting case — an untextured specimen must reproduce the "
        "random powder the uncorrected pattern already assumes, giving a factor of 1 for every "
        "reflection with no fitted parameter anywhere."
    ),
    setup=PREFERRED_ORIENTATION_SETUP,
    code=(
        "grid = OrientationSet.from_equispaced_so3_grid(\n"
        "    12.0,\n"
        "    specimen_frame=specimen,\n"
        "    phase=nickel,\n"
        "    reduce_to_fundamental_region=False,\n"
        ")\n"
        "model = ODFPreferredOrientationModel(odf=ODF.from_orientations(grid))\n"
        "result = float(model.factors([MillerPlane.from_hkl([1, 1, 1], phase=nickel)])[0])"
    ),
    expected=1.0,
    unit="",
    tolerance=0.01,
    reference=(
        "Pole density is defined in multiples of a random distribution, so a uniform orientation "
        "distribution has pole density 1 along every specimen direction and the correction is the "
        "identity. The tolerance reflects the finite SO(3) grid used to represent the uniform "
        "distribution, not any approximation in the correction itself."
    ),
    citation=(
        "Bunge, H.-J., Texture Analysis in Materials Science, "
        "DOI: 10.1016/C2013-0-11769-2; Von Dreele, R. B., J. Appl. Cryst. 30, 517-525 (1997), "
        "DOI: 10.1107/S0021889897005918 (texture in Rietveld refinement)."
    ),
    symbols=(_PO_FACTOR,),
    see_also=(_XRD_WORKFLOW, _XRD_CONCEPT),
    result_format="{:.4f}",
)


POWDER_PROFILE_AFFINE_COMPARISON = WorkedExample(
    id="diffraction-powder-profile-affine-comparison",
    title="Measured powder-profile comparison recovers a known scale and background",
    domain="diffraction",
    scenario=(
        "You have imported a measured powder profile and want an auditable first comparison with "
        "a simulated profile before attempting any structural refinement. This deliberately "
        "synthetic validation case sets I_obs = 5 I_sim + 5 at five points with equal standard "
        "uncertainty. Weighted least squares must therefore recover scale 5 and background 5 "
        "exactly, while both IUCr profile residuals vanish."
    ),
    setup=(
        NICKEL_SETUP
        + """
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
"""
    ),
    code=(
        "comparison = compare_powder_patterns(measured, simulated)\n"
        "result = np.array([\n"
        "    comparison.scale_factor,\n"
        "    comparison.background_offset,\n"
        "    comparison.profile_r_factor,\n"
        "    comparison.weighted_profile_r_factor,\n"
        "])"
    ),
    expected=[5.0, 5.0, 0.0, 0.0],
    unit="",
    tolerance=1e-12,
    reference=(
        "The five observed values are constructed independently as 5*x + 5 from x = 1,...,5. "
        "The weighted design matrix therefore contains the exact affine solution (5, 5), every "
        "residual is zero, and the numerators of both R_p and R_wp are exactly zero."
    ),
    citation=(
        "IUCr pdCIF dictionary definitions _pd_proc_ls.prof_R_factor and "
        "_pd_proc_ls.prof_wR_factor; Young, The Rietveld Method (1993), Ch. 1."
    ),
    symbols=(
        SymbolUse(r"R_p", "Unweighted whole-profile agreement factor."),
        SymbolUse(r"R_{wp}", "Weighted whole-profile agreement factor."),
    ),
    see_also=(_XRD_WORKFLOW, _XRD_THEORY, _DIFF_CONCEPT),
    result_format="{:.12f}",
)


KIKUCHI_MAP_SETUP = """
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
"""

_TILT = SymbolUse(
    r"\theta",
    "Angle between two crystal zone axes; the stage travel of one routing leg.",
)


KIKUCHI_MAP_ZONE_AXIS_TILT_ANGLES = WorkedExample(
    id="diffraction-kikuchi-map-zone-axis-tilt-angles",
    title="Kikuchi-map routing reproduces the exact cubic zone-axis angles",
    domain="diffraction",
    scenario=(
        "Build the stereographic Kikuchi map of nickel and ask it for the tilt "
        "from [001] to [011], to [111], and to [112], with a leg budget large "
        "enough that each is a single hop along one band. The angles between "
        "low-index cubic directions are closed-form - 45 degrees, "
        "arccos(1/sqrt(3)), and arccos(2/sqrt(6)) - so the routed travel is "
        "checked against arithmetic rather than against a prior run. Getting "
        "these right exercises the whole chain: the direct basis, the map frame, "
        "the zone law that decides which bands join two axes, and the "
        "shortest-path search."
    ),
    setup=KIKUCHI_MAP_SETUP,
    code=(
        "kikuchi_map = compute_kikuchi_map(\n"
        "    nickel,\n"
        "    beam_energy_kev=200.0,\n"
        "    max_index=4,\n"
        "    zone_axis_max_index=3,\n"
        ")\n"
        "targets = ([0, 1, 1], [1, 1, 1], [1, 1, 2])\n"
        "result = np.array(\n"
        "    [\n"
        "        kikuchi_map.route_to([0, 0, 1], target, max_leg_deg=90.0).total_tilt_deg\n"
        "        for target in targets\n"
        "    ]\n"
        ")"
    ),
    expected=[45.0, 54.735610317245346, 35.264389682754654],
    unit="deg",
    tolerance=1e-6,
    reference=(
        "Closed-form angles between cubic directions: arccos(1/sqrt(2)) = 45 deg "
        "for [001]-[011], arccos(1/sqrt(3)) = 54.735610 deg for [001]-[111], and "
        "arccos(2/sqrt(6)) = 35.264390 deg for [001]-[112]. The two [111] and "
        "[112] values are complementary, summing to 90 degrees, because [112] is "
        "the reflection of [001] in the plane perpendicular to [111]."
    ),
    citation=(
        "Standard cubic interaxial angles; see International Tables for "
        "Crystallography Vol. C (1999) for the reciprocal-lattice conventions, "
        "and Williams and Carter, Transmission Electron Microscopy 2nd ed. "
        "(2009) Ch. 19 for Kikuchi-map tilting."
    ),
    symbols=(_TILT,),
    see_also=(_XRD_WORKFLOW, _XRD_CONCEPT),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="diffraction",
    title="Diffraction geometry",
    summary=(
        "Powder scattering angles from PyTex interplanar spacings via Bragg's law, Kikuchi band "
        "and zone-axis geometry in the gnomonic projection, zone-axis routing on a "
        "stereographic Kikuchi map, and preferred-orientation corrections to powder "
        "intensities — each checked against a standard reference value or a closed-form "
        "identity."
    ),
    examples=(
        NI_111_TWO_THETA,
        NI_111_KIKUCHI_BAND_WIDTH,
        CUBE_ZONE_AXIS_GNOMONIC_RADIUS,
        MARCH_DOLLASE_FAMILY_FACTOR,
        MARCH_DOLLASE_NORMALIZATION,
        ODF_WEIGHTED_RANDOM_TEXTURE,
        POWDER_PROFILE_AFFINE_COMPARISON,
        KIKUCHI_MAP_ZONE_AXIS_TILT_ANGLES,
    ),
)

__all__ = ["GROUP"]
