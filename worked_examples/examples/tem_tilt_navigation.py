"""Worked examples: TEM specimen-tilt navigation.

These examples check the tilt engine against quantities that are known
independently of any program run: analytic interzonal angles, the closed-form
solid angle of a tilt envelope, the residual law for a mis-calibrated
diffraction rotation, and the order of the observation stabilizer for a cubic
crystal down ``[001]``.

The last two are the ones worth reading. The residual law is the reason a
long excursion should be broken into short hops; the stabilizer order is the
reason the much-feared 180-degree ambiguity is *harmless* for a centrosymmetric
crystal, and the reason the engine stays silent for one.

See :doc:`../../architecture/tem_tilt_navigation_foundation`.
"""

from __future__ import annotations

import math

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

TEM_SETUP = """
import numpy as np
from pytex import (
    CurrentState,
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    plan_tilt_to_zone_axis,
)
from pytex.tem.reconstruction import HOLDER_FRAME

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
# A generous holder, so the geometry rather than the envelope decides.
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
# The crystal sits with [001] along the beam at zero tilt.
aligned = CurrentState.from_orientation(
    Orientation.from_matrix(
        np.eye(3), specimen_frame=HOLDER_FRAME, phase=nickel, crystal_frame=crystal
    ),
    StagePosition(0.0, 0.0),
    current_zone_axis=ZoneAxis([0, 0, 1], phase=nickel),
)
"""

_ALPHA = SymbolUse(r"\alpha", "Holder tilt about the rod axis.")
_BETA = SymbolUse(r"\beta", "Holder tilt about the cradle axis carried in the rod.")
_OMEGA = SymbolUse(r"\Omega", "Solid angle of beam directions a holder can reach.")
_ZONE = SymbolUse(r"[uvw]", "Lattice direction brought parallel to the beam.")
_THETA = SymbolUse(r"\theta", "Angle between the current and target zone axes.")

_FOUNDATION = SeeAlso(
    "TEM tilt navigation foundation", "../../architecture/tem_tilt_navigation_foundation"
)
_NOTEBOOK = SeeAlso(
    "TEM tilt navigation notebook", "../../tutorials/notebooks/24_tem_tilt_navigation"
)


TILT_TO_011_FROM_001 = WorkedExample(
    id="tem-tilt-001-to-011-travel",
    title="Crystal travel from [001] to [011] in a cubic crystal",
    domain="tem",
    scenario=(
        "You are down the [001] zone of an FCC metal and want [011]. The engine solves the holder "
        "angles and plans the path; the crystal travel along that path must equal the interplanar "
        "angle between the two zone axes, which for a cubic crystal is an analytic quantity "
        "independent of lattice parameter. This is the basic sanity check on the whole chain: "
        "orientation, closed-form solution, forward validation, and geodesic path planning."
    ),
    setup=TEM_SETUP,
    code=(
        "report = plan_tilt_to_zone_axis(\n"
        "    aligned, ZoneAxis([0, 1, 1], phase=nickel), wide_stage\n"
        ")\n"
        "result = float(report.best().path.total_travel_deg)"
    ),
    expected=45.0,
    unit="deg",
    tolerance=1e-3,
    reference=(
        "The angle between [001] and [011] in a cubic lattice is "
        "arccos(1/sqrt(2)) = 45 degrees exactly, from the dot product of the two "
        "directions divided by their lengths. Independent of the lattice parameter."
    ),
    citation=(
        "Edington, J. W., Practical Electron Microscopy in Materials Science, Macmillan; "
        "standard cubic interzonal-angle tables."
    ),
    symbols=(_ZONE, _ALPHA, _BETA),
    see_also=(_FOUNDATION, _NOTEBOOK),
    result_format="{:.3f}",
)


TILT_TO_111_FROM_001 = WorkedExample(
    id="tem-tilt-001-to-111-travel",
    title="Crystal travel from [001] to [111] in a cubic crystal",
    domain="tem",
    scenario=(
        "The [001] to [111] hop is the one every TEM course teaches, and it is the longest of the "
        "common cubic transitions — far enough that a typical +/-30 degree double-tilt holder "
        "cannot make it without help from a symmetry equivalent. Here a wide holder is used so "
        "that the geometry alone is tested."
    ),
    setup=TEM_SETUP,
    code=(
        "report = plan_tilt_to_zone_axis(\n"
        "    aligned, ZoneAxis([1, 1, 1], phase=nickel), wide_stage\n"
        ")\n"
        "result = float(report.best().path.total_travel_deg)"
    ),
    expected=54.735610317245346,
    unit="deg",
    tolerance=1e-3,
    reference=(
        "arccos(1/sqrt(3)) = 54.7356 degrees, the analytic angle between <001> and "
        "<111> in a cubic lattice, and the standard tetrahedral-angle complement."
    ),
    citation=(
        "Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, Springer, "
        "DOI: 10.1007/978-0-387-76501-3."
    ),
    symbols=(_ZONE, _THETA),
    see_also=(_FOUNDATION, _NOTEBOOK),
    result_format="{:.4f}",
)


HOLDER_SOLID_ANGLE = WorkedExample(
    id="tem-holder-accessible-solid-angle",
    title="Solid angle a +/-30 degree double-tilt holder reaches",
    domain="tem",
    scenario=(
        "Before asking whether a particular zone axis is reachable, it is worth knowing how much "
        "of orientation space the holder commands at all. Because the beam direction in holder "
        "coordinates is a spherical coordinate system whose pole is the beta axis, the Jacobian is "
        "cos(alpha) and the integral is elementary. The answer — a little over eight percent of "
        "all directions — is why symmetry equivalents matter so much in practice."
    ),
    setup=TEM_SETUP,
    code=(
        "envelope = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)\n"
        "result = float(envelope.accessible_solid_angle_sr())"
    ),
    expected=math.radians(60.0) * 2.0 * math.sin(math.radians(30.0)),
    unit="sr",
    tolerance=1e-9,
    reference=(
        "Omega = (beta_max - beta_min) * (sin alpha_max - sin alpha_min) = "
        "(pi/3) * (2 sin 30 deg) = pi/3 = 1.04720 sr, integrating the cos(alpha) "
        "Jacobian of the beam-direction map over the tilt rectangle."
    ),
    citation=(
        "Derived in section 10.2 of docs/architecture/tem_tilt_navigation_foundation.md; "
        "standard spherical-measure result."
    ),
    symbols=(_OMEGA, _ALPHA, _BETA),
    see_also=(_FOUNDATION,),
    result_format="{:.5f}",
)


ROTATION_ERROR_RESIDUAL = WorkedExample(
    id="tem-diffraction-rotation-residual",
    title="Angular miss from a 5 degree diffraction-rotation error over a 90 degree hop",
    domain="tem",
    scenario=(
        "The diffraction rotation is not recorded by instrument metadata and must be calibrated. "
        "This example quantifies what an uncalibrated value costs: the miss is "
        "2 asin(sin(dphi/2) sin(theta)), which grows with the length of the hop. The same 5 degree "
        "error costs 0.44 degrees over a 5 degree hop and the full 5 degrees over a 90 degree one "
        "— which is the argument for routing a long excursion through intermediate zones and "
        "re-indexing at each."
    ),
    setup="from pytex.tem.calibration import residual_from_rotation_error_deg",
    code="result = residual_from_rotation_error_deg(5.0, 90.0)",
    expected=5.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "For theta = 90 degrees, sin(theta) = 1, so the expression reduces to "
        "2 asin(sin(dphi/2)) = dphi exactly. The residual equals the calibration "
        "error itself when the target is perpendicular to the current zone axis."
    ),
    citation=(
        "Derived in section 8.2 of docs/architecture/tem_tilt_navigation_foundation.md "
        "and verified numerically over 3000 random orientations."
    ),
    symbols=(_THETA,),
    see_also=(_FOUNDATION,),
    result_format="{:.6f}",
)


OBSERVATION_STABILIZER_ORDER = WorkedExample(
    id="tem-observation-stabilizer-cubic-001",
    title="Order of the observation stabilizer for cubic m-3m down [001]",
    domain="tem",
    scenario=(
        "A single indexed SAED pattern determines the orientation only up to the rotations of the "
        "Laue class that map the zone plane to itself. Counting them answers the question that "
        "decides whether the classical 180-degree ambiguity matters: for cubic m-3m down [001] "
        "the stabilizer is the group 422, of order 8, and every one of its operators is already a "
        "crystal symmetry — so the ambiguity is entirely absorbed and nothing is left "
        "undetermined."
    ),
    setup="from pytex.tem.ambiguity import observation_stabilizer",
    code='result = len(observation_stabilizer("m-3m", [0.0, 0.0, 1.0]))',
    expected=8,
    unit="operators",
    tolerance=0,
    reference=(
        "The rotations of m-3m fixing the [001] axis line form the point group 422: "
        "the identity, three rotations about [001] (90, 180, 270 degrees), and four "
        "two-fold rotations about the in-plane <100> and <110> axes. Order 8, from "
        "International Tables Volume A."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, IUCr/Springer, "
        "DOI: 10.1107/97809553602060000100."
    ),
    symbols=(_ZONE,),
    see_also=(_FOUNDATION,),
    result_format="{:d}",
)


CUBIC_ORBIT_MULTIPLICITY = WorkedExample(
    id="tem-symmetry-orbit-multiplicity",
    title="Number of symmetry-equivalent targets for a general cubic direction",
    domain="tem",
    scenario=(
        "The user asks for one zone axis; the crystal offers a whole orbit, every member of which "
        "gives an identical diffraction pattern. The choice among them is therefore free, and the "
        "engine takes the cheapest reachable one — which is frequently the difference between a "
        "target being reachable and not. For a general direction in a cubic crystal the orbit has "
        "one member per proper operator, in both senses."
    ),
    setup=TEM_SETUP,
    code=(
        "report = plan_tilt_to_zone_axis(\n"
        "    aligned, ZoneAxis([1, 3, 5], phase=nickel), wide_stage, include_paths=False\n"
        ")\n"
        "result = int(report.orbit_size)"
    ),
    expected=48,
    unit="directions",
    tolerance=0,
    reference=(
        "The proper point group of m-3m is 432, of order 24. A general direction has a "
        "trivial stabilizer, so its orbit has 24 members; counting both senses of each "
        "gives 48 distinct directions."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, IUCr/Springer, "
        "DOI: 10.1107/97809553602060000100."
    ),
    symbols=(_ZONE,),
    see_also=(_FOUNDATION, _NOTEBOOK),
    result_format="{:d}",
)


GROUP = ExampleGroup(
    slug="tem_tilt_navigation",
    title="TEM tilt navigation",
    summary=(
        "Holder tilts that bring a target zone axis onto the electron beam: analytic interzonal "
        "travel for the standard cubic transitions, the closed-form solid angle a double-tilt "
        "holder commands, the cost of an uncalibrated diffraction rotation, and the group-order "
        "counts that decide whether a single indexed pattern leaves a real ambiguity."
    ),
    examples=(
        TILT_TO_011_FROM_001,
        TILT_TO_111_FROM_001,
        HOLDER_SOLID_ANGLE,
        ROTATION_ERROR_RESIDUAL,
        OBSERVATION_STABILIZER_ORDER,
        CUBIC_ORBIT_MULTIPLICITY,
    ),
)

__all__ = ["GROUP"]
