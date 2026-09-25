"""Worked examples: FIB lamella planning for a target zone axis.

Every expected value here has independent provenance: an orientation built so
the out-of-plane angle of the target is known by construction, the arcsine of a
product of sines worked by hand, the half-range of a rectangle, the side of a
rectangle less the lamella width halved, a sum of squares, or the sense and
offset a set of fiducials was generated with. None of them is a copied output of
this code.

See :doc:`../../theory/fib_lamella_zone_axis_geometry`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

_THEORY = SeeAlso(
    "FIB lamella geometry for a target zone axis", "../../theory/fib_lamella_zone_axis_geometry"
)
_WORKFLOW = SeeAlso("FIB lamella planning", "../../workflows/fib_lamella_planning")

_EPS = SymbolUse(r"\varepsilon^{*}", "Residual tilt: the best member's angle out of the surface.")
_THETA = SymbolUse(r"\theta_{S}", "Azimuth of the lamella normal in the sample frame.")
_PHI = SymbolUse(r"\phi", "Unknown in-plane mounting rotation of the lamella on the grid.")

_SETUP_PHASES = (
    "import math\n"
    "import numpy as np\n"
    "from pytex.core import frame_catalog\n"
    "from pytex.core.lattice import Lattice, Phase\n"
    "from pytex.core.orientation import Rotation\n"
    "from pytex.core.symmetry import SymmetrySpec\n"
    "from pytex.fib import lamella_geometry, target_orbit\n"
    "def phase(point_group, a, c):\n"
    "    frame = frame_catalog.crystal_frame()\n"
    "    return Phase(point_group, lattice=Lattice(a, a, c, 90, 90, 90, crystal_frame=frame),\n"
    "                 symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),\n"
    "                 crystal_frame=frame)\n"
    "def about(axis, degrees):\n"
    "    return Rotation.from_axis_angle(axis, math.radians(degrees)).as_matrix()\n"
)


CONSTRUCTED_RESIDUAL = WorkedExample(
    id="fib-constructed-residual",
    title="A target built to rise 12.5 degrees out of the surface",
    domain="fib-lamella",
    scenario=(
        "A tetragonal crystal whose [001] axis is rotated 77.5 degrees about the sample X axis: "
        "[001] then points at (0, -sin 77.5, cos 77.5), which rises arcsin(cos 77.5) = 12.5 "
        "degrees out of the surface. Its orbit is only +/-[001], so no other member can come "
        "closer: eps* is exactly 12.5 degrees, the tilt the TEM holder must supply."
    ),
    setup=_SETUP_PHASES,
    code=(
        "orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))\n"
        "result = lamella_geometry(about([1, 0, 0], 77.5), orbit).eps_deg"
    ),
    expected=12.5,
    unit="deg",
    tolerance=1e-9,
    reference="arcsin(cos 77.5 deg) = 90 - 77.5 = 12.5 deg, by construction.",
    citation=(
        "International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100."
    ),
    symbols=(_EPS,),
    see_also=(_THEORY,),
    result_format="{:.9f}",
)

IN_PLANE_TARGET = WorkedExample(
    id="fib-in-plane-target",
    title="A cube-oriented grain has <100> in the surface",
    domain="fib-lamella",
    scenario=(
        "At the identity orientation a cubic crystal's [100] and [010] lie in the surface, so a "
        "<100> lamella needs no residual tilt at all: eps* = 0, and the lamella normal points "
        "along X_s (the smallest canonical azimuth among the tied options)."
    ),
    setup=_SETUP_PHASES,
    code=(
        "orbit = target_orbit(phase('m-3m', 3.52, 3.52), (1, 0, 0))\n"
        "result = lamella_geometry(np.eye(3), orbit).eps_deg"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-12,
    reference="[100] . Z_s = 0 at the identity orientation.",
    citation=(
        "International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100."
    ),
    symbols=(_EPS,),
    see_also=(_THEORY,),
    result_format="{:.6f}",
)

CUBIC_111_ALONG_NORMAL = WorkedExample(
    id="fib-111-along-normal",
    title="With [111] along the normal, every <001> is 35.26 degrees out",
    domain="fib-lamella",
    scenario=(
        "When a cubic grain has [111] along the surface normal, each <001> axis makes the same "
        "angle arccos(1/sqrt 3) with it, so each rises arcsin(1/sqrt 3) = 35.264 degrees out of "
        "the surface. That exceeds a +/-30 degree holder's guaranteed reach: no vertical lamella "
        "of this grain shows <001> for every mounting rotation."
    ),
    setup=_SETUP_PHASES,
    code=(
        "target = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)\n"
        "axis = np.cross(target, [0.0, 0.0, 1.0])\n"
        "g = about(axis / np.linalg.norm(axis), math.degrees(math.acos(target[2])))\n"
        "orbit = target_orbit(phase('m-3m', 3.52, 3.52), (0, 0, 1))\n"
        "result = lamella_geometry(g, orbit).eps_deg"
    ),
    expected=35.26438968,
    unit="deg",
    tolerance=1e-7,
    reference="arcsin(1/sqrt(3)) = 35.2643897 deg, by hand.",
    citation=(
        "International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100."
    ),
    symbols=(_EPS,),
    see_also=(_THEORY,),
    result_format="{:.8f}",
)

CONVENTION = WorkedExample(
    id="fib-orientation-convention",
    title="The orientation matrix maps crystal to specimen",
    domain="fib-lamella",
    scenario=(
        "PyTex orientations map crystal to specimen, v_s = g v_c. For a triclinic crystal with "
        "g = Rx(30) Rz(40), g [100] = (cos 40, sin 40 cos 30, sin 40 sin 30), which rises "
        "arcsin(sin 40 sin 30) = 18.747 degrees out of the surface. Read backwards, g^T [100] "
        "would lie in the surface and report zero; the geometry must not do that."
    ),
    setup=_SETUP_PHASES,
    code=(
        "g = about([1, 0, 0], 30.0) @ about([0, 0, 1], 40.0)\n"
        "orbit = target_orbit(phase('1', 4.0, 4.0), (1, 0, 0))\n"
        "result = lamella_geometry(g, orbit).eps_deg"
    ),
    expected=18.74723725,
    unit="deg",
    tolerance=1e-7,
    reference="arcsin(sin 40 deg x sin 30 deg) = arcsin(0.3213938) = 18.7472373 deg, by hand.",
    citation="Morawiec, Orientations and Rotations, Springer, doi:10.1007/978-3-662-09156-2.",
    symbols=(_EPS,),
    see_also=(_THEORY,),
    result_format="{:.8f}",
)

LAMELLA_AZIMUTH = WorkedExample(
    id="fib-lamella-azimuth",
    title="The lamella normal follows the target's in-plane part",
    domain="fib-lamella",
    scenario=(
        "A tetragonal [001] placed at (cos 20 cos 35, cos 20 sin 35, sin 20) rises 20 degrees out "
        "of the surface with its in-plane part at 35 degrees from X_s. The lamella normal is that "
        "in-plane part, normalized, so the azimuth to cut at is theta_S = 35 degrees."
    ),
    setup=_SETUP_PHASES,
    code=(
        "d = np.array([math.cos(math.radians(20)) * math.cos(math.radians(35)),\n"
        "              math.cos(math.radians(20)) * math.sin(math.radians(35)),\n"
        "              math.sin(math.radians(20))])\n"
        "axis = np.cross([0.0, 0.0, 1.0], d)\n"
        "g = about(axis / np.linalg.norm(axis), math.degrees(math.acos(d[2])))\n"
        "orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))\n"
        "result = lamella_geometry(g, orbit).theta_sample_deg"
    ),
    expected=35.0,
    unit="deg",
    tolerance=1e-9,
    reference="atan2(sin 35, cos 35) = 35 deg, by construction.",
    citation=(
        "International Tables for Crystallography, Vol. A, doi:10.1107/97809553602060000100."
    ),
    symbols=(_THETA,),
    see_also=(_THEORY, _WORKFLOW),
    result_format="{:.9f}",
)

EXACT_MOUNT_TILT = WorkedExample(
    id="fib-exact-mount-tilt",
    title="The alpha tilt for a lamella mounted at 45 degrees",
    domain="fib-lamella",
    scenario=(
        "A lamella whose target sits eps* = 20 degrees out of its plane, mounted on the grid at "
        "phi = 45 degrees, needs alpha = arcsin(cos phi sin eps*) = arcsin(0.70711 x 0.34202) "
        "= 13.995 degrees - not the small-angle eps* cos phi = 14.142 degrees. The value is read "
        "from the solve at phi = 45 degrees of the mounting sweep."
    ),
    setup=_SETUP_PHASES + "from pytex.fib import MountModel, sweep_mounting_rotation\n",
    code=(
        "orbit = target_orbit(phase('4/mmm', 3.99, 4.04), (0, 0, 1))\n"
        "g = about([1, 0, 0], 70.0)\n"
        "geometry = lamella_geometry(g, orbit)\n"
        "sweep = sweep_mounting_rotation(g, orbit, geometry,\n"
        "                                mount=MountModel(phi_samples=8, solver='closed_form'))\n"
        "result = abs(float(sweep.alpha_front_deg[1]))"
    ),
    expected=13.99545,
    unit="deg",
    tolerance=1e-5,
    reference=(
        "arcsin(cos 45 deg x sin 20 deg) = arcsin(0.2418448) = 13.99545 deg, by hand; the "
        "double-tilt closed form of De Graef, doi:10.1017/CBO9780511615092."
    ),
    citation="De Graef, Introduction to Conventional TEM, doi:10.1017/CBO9780511615092.",
    symbols=(_EPS, _PHI),
    see_also=(_THEORY,),
    result_format="{:.5f}",
)

GUARANTEED_RADIUS = WorkedExample(
    id="fib-guaranteed-radius",
    title="The residual a +/-30 x +/-25 degree holder reaches for every mount",
    domain="fib-lamella",
    scenario=(
        "On the exact tilt curve |alpha| and |beta| never exceed eps*, and each reaches it on an "
        "axis, so a symmetric rectangular holder guarantees exactly its smaller half-range: 25 "
        "degrees for alpha +/-30 and beta +/-25."
    ),
    setup=(
        "from pytex.fib import guaranteed_radius_deg\n"
        "from pytex.tem.stage import RectangularEnvelope\n"
    ),
    code="result = guaranteed_radius_deg(RectangularEnvelope(-30.0, 30.0, -25.0, 25.0))",
    expected=25.0,
    unit="deg",
    tolerance=1e-12,
    reference="min(30, 25) = 25 deg.",
    citation="De Graef, Introduction to Conventional TEM, doi:10.1017/CBO9780511615092.",
    symbols=(_EPS, _PHI),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)

FOOTPRINT_CLEARANCE = WorkedExample(
    id="fib-footprint-clearance",
    title="How far inside a 40 x 20 micrometre grain a lamella can sit",
    domain="fib-lamella",
    scenario=(
        "A rectangular grain 40 um long and 20 um tall (0.5 um step) holds a 15 x 2 um lamella "
        "whose normal is along the short side. The rectangle can grow by (20 - 2)/2 = 9 um "
        "across it before touching the boundary, and by (40 - 15)/2 = 12.5 um along it; the "
        "clearance is the smaller, 9 um."
    ),
    setup=(
        "import numpy as np\n"
        "from pytex.fib import fit_footprint\n"
        "mask = np.zeros((60, 100), dtype=bool)\n"
        "mask[10:50, 10:90] = True\n"
    ),
    code=(
        "result = fit_footprint(mask, step_um=(0.5, 0.5), theta_sample_deg=90.0,\n"
        "                       length_um=15.0, width_um=2.0).margin_um"
    ),
    expected=9.0,
    unit="um",
    tolerance=1e-12,
    reference="min((20 - 2)/2, (40 - 15)/2) = 9 um, by hand.",
    citation="Crow, Summed-area tables for texture mapping, SIGGRAPH 1984, doi:10.1145/964965.808600.",
    see_also=(_WORKFLOW,),
    result_format="{:.3f}",
)

UNCERTAINTY_BUDGET = WorkedExample(
    id="fib-uncertainty-budget",
    title="The expanded uncertainty of a residual tilt",
    domain="fib-lamella",
    scenario=(
        "EBSD accuracy 0.5 deg, grain spread 1.5 deg, an exact registration and a mount "
        "repeatability of 2 deg combine in quadrature to u = sqrt(0.25 + 2.25 + 0 + 4) = "
        "2.5495 deg; with k = 2 the verdict is judged on eps* + 5.099 deg."
    ),
    setup=(
        "import numpy as np\n"
        "from pytex.fib import MountModel, UncertaintyInputs, plan_lamella\n"
        "from pytex.core import frame_catalog\n"
        "from pytex.core.lattice import Lattice, Phase\n"
        "from pytex.core.symmetry import SymmetrySpec\n"
        "frame = frame_catalog.crystal_frame()\n"
        "nickel = Phase('nickel', lattice=Lattice(3.52, 3.52, 3.52, 90, 90, 90, crystal_frame=frame),\n"
        "               symmetry=SymmetrySpec.from_point_group('m-3m', reference_frame=frame),\n"
        "               crystal_frame=frame)\n"
    ),
    code=(
        "plan = plan_lamella(np.eye(3), (1, 0, 0), phase=nickel, grain_spread_deg=1.5,\n"
        "                    uncertainty=UncertaintyInputs(0.5, 2.0, 2.0),\n"
        "                    mount=MountModel(phi_samples=8, solver='closed_form'),\n"
        "                    crosscheck=False)\n"
        "result = plan.budget.expanded_deg"
    ),
    expected=5.0990195,
    unit="deg",
    tolerance=1e-6,
    reference="2 sqrt(0.5^2 + 1.5^2 + 0^2 + 2^2) = 2 sqrt(6.5) = 5.0990195 deg, by hand.",
    citation="JCGM 100:2008, Guide to the expression of uncertainty in measurement (GUM).",
    symbols=(_EPS,),
    see_also=(_THEORY,),
    result_format="{:.7f}",
)

FIDUCIAL_CALIBRATION = WorkedExample(
    id="fib-fiducial-calibration",
    title="Recovering a reversed FIB rotation sense and its offset",
    domain="fib-lamella",
    scenario=(
        "Fiducial trenches milled at pattern rotations 0, 30 and 60 degrees on an instrument with "
        "s_R = -1 and R0 = 90 degrees have their normals at 90, 60 and 30 degrees in the sample "
        "frame. The calibration recovers the offset R0 = 90 degrees (and the sense -1) from them."
    ),
    setup=(
        "from pytex.fib import FiducialObservation, calibrate_chamber_from_fiducials\n"
        "observations = [FiducialObservation(rho, (90.0 - rho) % 180.0) for rho in (0.0, 30.0, 60.0)]\n"
    ),
    code=(
        "chamber = calibrate_chamber_from_fiducials(observations)\n"
        "result = chamber.rotation_sense * 1000.0 + (chamber.rotation_offset_deg % 180.0)"
    ),
    expected=-910.0,
    unit="sense x 1000 + R0 (deg)",
    tolerance=1e-9,
    reference="The sense -1 and offset 90 deg the observations were generated with: -1000 + 90.",
    citation="Giannuzzi & Stevie, Micron 30 (1999) 197, doi:10.1016/S0968-4328(99)00005-0.",
    see_also=(_WORKFLOW,),
    result_format="{:.3f}",
)


GROUP = ExampleGroup(
    slug="fib-lamella",
    title="FIB lamella planning",
    summary=(
        "The residual tilt of a vertically milled lamella checked on orientations built so it is "
        "known (0, 12.5 and 35.26 degrees), the crystal-to-specimen convention on a case that "
        "fails if it is read backwards, the lamella azimuth, the exact holder tilt for a mounted "
        "lamella against its small-angle approximation, a holder's guaranteed reach, the "
        "footprint clearance in a hand-measured grain, the uncertainty budget, and the recovery "
        "of a reversed FIB rotation from fiducials."
    ),
    examples=(
        CONSTRUCTED_RESIDUAL,
        IN_PLANE_TARGET,
        CUBIC_111_ALONG_NORMAL,
        CONVENTION,
        LAMELLA_AZIMUTH,
        EXACT_MOUNT_TILT,
        GUARANTEED_RADIUS,
        FOOTPRINT_CLEARANCE,
        UNCERTAINTY_BUDGET,
        FIDUCIAL_CALIBRATION,
    ),
)
