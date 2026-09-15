"""Worked examples: sample symmetry, ODF section ranges and component volume fractions.

Three surfaces of the texture analysis, each checked against an answer that
follows from its definition rather than from a prior run of the code:

* imposing **axial** sample symmetry on a pole figure removes exactly the
  azimuth-dependent part and leaves the polar part untouched;
* the fraction of a **random** texture within a misorientation ball of an ideal
  orientation is the Haar volume of the ball times the number of
  symmetry-equivalent balls;
* the Euler box an ODF section must cover for a **hexagonal** crystal with
  orthorhombic sample symmetry is 90, 90 and 60 degrees.

See :doc:`../../workflows/texture_analysis_workbench`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, WorkedExample

SETUP = """
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
"""

_WORKFLOW = SeeAlso(
    "Texture analysis in the workbench",
    "../../workflows/texture_analysis_workbench",
)

AXIAL_AVERAGE_KEEPS_THE_POLAR_PART = WorkedExample(
    id="axial-sample-symmetry-keeps-the-polar-part",
    title="Axial sample symmetry removes the azimuthal part of a pole figure, and nothing else",
    domain="texture",
    scenario=(
        "A drawn wire, an extruded rod or a tube read along its axis is taken to be axially "
        "symmetric: its texture is unchanged by any rotation about the axis. Imposing that "
        "symmetry on a measured pole figure replaces each value by the average over its ring of "
        "constant tilt. Take a figure measured on a goniometer raster - tilt 0 to 75 degrees in "
        "5 degree steps, azimuth every 5 degrees - whose intensity is 2 + 0.6x + 0.3xy + 0.5z^2 in "
        "specimen coordinates. The x and xy terms depend on azimuth and must vanish; the z^2 term "
        "depends on tilt alone and must survive unchanged. The example returns the largest "
        "departure of the symmetrized figure from 2 + 0.5z^2."
    ),
    setup=SETUP,
    code=(
        "psi, phi = np.meshgrid(np.arange(0.0, 76.0, 5.0), np.arange(0.0, 360.0, 5.0), "
        'indexing="ij")\n'
        "theta, azimuth = np.radians(psi.ravel()), np.radians(phi.ravel())\n"
        "directions = np.column_stack(\n"
        "    [np.sin(theta) * np.cos(azimuth), np.sin(theta) * np.sin(azimuth), np.cos(theta)]\n"
        ")\n"
        "x, y, z = directions.T\n"
        "figure = PoleFigure(\n"
        "    pole=BASAL,\n"
        "    sample_directions=directions,\n"
        "    intensities=2.0 + 0.6 * x + 0.3 * x * y + 0.5 * z**2,\n"
        "    specimen_frame=SPECIMEN,\n"
        '    sampling="sampled_density",\n'
        ")\n"
        'axial = impose_sample_symmetry(figure, "axial")\n'
        "result = float(np.max(np.abs(np.asarray(axial.intensities) - (2.0 + 0.5 * z**2))))"
    ),
    expected=0.0,
    unit="m.r.d.",
    tolerance=1e-9,
    reference=(
        "Analytic: on a ring of constant tilt theta, x = sin(theta) cos(psi) and "
        "xy = sin^2(theta) sin(psi) cos(psi) have zero mean over the azimuth psi, while z = "
        "cos(theta) is constant. For azimuths equally spaced round the full circle the discrete "
        "means of cos(psi) and sin(2 psi) are exactly zero as well, so the residual is "
        "floating-point rather than quadrature error."
    ),
    citation=(
        "H.-J. Bunge, Texture Analysis in Materials Science (1982), section 4.2: a fibre "
        "(cylindrical) sample symmetry makes the pole figure independent of the azimuth about "
        "the fibre axis."
    ),
    see_also=(_WORKFLOW,),
    result_format="{:.2e}",
)

RANDOM_COMPONENT_FRACTION_CUBIC = WorkedExample(
    id="random-texture-fraction-within-15-degrees-cubic",
    title="A random cubic texture holds 2.28% of its volume within 15 degrees of any orientation",
    domain="texture",
    scenario=(
        "A component volume fraction - '9% cube within 15 degrees' - is only meaningful beside "
        "what a texture-free specimen gives for the same tolerance. Under the invariant measure "
        "on rotations, a ball of radius w holds (w - sin w)/pi of orientation space, and a cubic "
        "ideal orientation has 24 symmetry-equivalent balls, which do not overlap below 45 "
        "degrees. The example evaluates the random reference the texture analysis prints beside "
        "every component."
    ),
    setup=SETUP,
    code="result = random_component_fraction(15.0, CUBIC.order)",
    expected=0.022768,
    unit="",
    tolerance=1e-6,
    reference=(
        "Analytic: w = 15 deg = 0.2617993878 rad and sin 15 deg = (sqrt(6) - sqrt(2))/4 = "
        "0.2588190451, so w - sin w = 0.0029803427; times |G| = 24 and divided by pi this is "
        "0.0227680."
    ),
    citation=(
        "A. Morawiec, Orientations and Rotations: Computations in Crystallographic Textures, "
        "Springer (2004), chapter 2: the invariant (Haar) measure on SO(3), under which the "
        "rotation angle w of a uniform rotation has density (1 - cos w)/pi."
    ),
    see_also=(_WORKFLOW,),
    result_format="{:.6f}",
)

HEXAGONAL_SECTION_RANGE = WorkedExample(
    id="hexagonal-orthorhombic-odf-section-box",
    title="A hexagonal ODF with orthorhombic sample symmetry spans 90, 90 and 60 degrees",
    domain="texture",
    scenario=(
        "Before drawing ODF sections, decide how much of Euler space they must cover. Too little "
        "hides part of the texture; too much repeats it. For a hexagonal crystal the six-fold "
        "axis along c repeats every orientation after 60 degrees of phi2, and a rolled sheet's "
        "orthorhombic symmetry, with the two-fold axes perpendicular to c, folds phi1 and Phi "
        "into 90 degrees. The example returns the box (phi1, Phi, phi2) the texture analysis "
        "derives from the operators."
    ),
    setup=SETUP,
    code=(
        'ranges = euler_section_ranges(HEXAGONAL, "orthorhombic")\n'
        'result = [ranges["phi1_max_deg"], ranges["big_phi_max_deg"], ranges["phi2_max_deg"]]'
    ),
    expected=[90.0, 90.0, 60.0],
    unit="deg",
    tolerance=1e-12,
    reference=(
        "Bunge (1982) section 4.2: the asymmetric Euler region for crystal symmetry 6/mmm (proper "
        "group 622) with orthorhombic sample symmetry is 0 <= phi1 <= 90, 0 <= Phi <= 90, "
        "0 <= phi2 <= 60 degrees."
    ),
    citation=(
        "H.-J. Bunge, Texture Analysis in Materials Science (1982), section 4.2 and its table of "
        "asymmetric units of Euler space."
    ),
    see_also=(_WORKFLOW,),
    result_format="{:.1f}",
)

GROUP = ExampleGroup(
    slug="texture-sections-and-sample-symmetry",
    title="Sample symmetry, ODF sections and component fractions",
    summary=(
        "The texture-analysis surfaces behind the workbench's measured-texture panel, each "
        "checked against its definition: axial sample symmetry removes exactly the azimuthal "
        "part of a figure, a random cubic texture holds the Haar volume of 24 balls within a "
        "tolerance, and a hexagonal ODF under orthorhombic sample symmetry spans 90, 90 and 60 "
        "degrees."
    ),
    examples=(
        AXIAL_AVERAGE_KEEPS_THE_POLAR_PART,
        RANDOM_COMPONENT_FRACTION_CUBIC,
        HEXAGONAL_SECTION_RANGE,
    ),
)
