"""Worked examples: the equal-volume charts, and naming an orientation.

Two things are pinned here that nothing else in the registry pins. The first is
the pair of constants behind the homochoric ball and the cubochoric cube: both
enclose the volume of SO(3) under its invariant measure, and the cube's corner
must land exactly on the ball's surface. A mis-remembered constant in the
equal-volume map fails both at once. The second is that naming an orientation
``(hkl)[uvw]`` really does invert the construction that builds one from those
indices.

See :doc:`../../concepts/orientation_texture` for the surrounding orientation
doctrine.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

REPRESENTATION_SETUP = """
import numpy as np
from pytex import (
    CUBOCHORIC_CUBE_EDGE,
    CUBOCHORIC_CUBE_HALF_EDGE,
    FrameDomain,
    Handedness,
    HOMOCHORIC_BALL_RADIUS,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    convert_orientations,
    ideal_orientation_indices,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
"""

_HOMOCHORIC = SymbolUse(
    r"\mathbf{h}",
    "Homochoric vector of a rotation; the equal-volume chart of SO(3).",
)
_CUBOCHORIC = SymbolUse(
    r"\mathbf{c}",
    "Cubochoric coordinate; the equal-volume chart mapped onto a cube.",
)
_BALL_RADIUS = SymbolUse(
    r"R_{1}",
    r"Radius of the homochoric ball, $(3\pi/4)^{1/3}$.",
)
_MILLER_PLANE = SymbolUse(
    r"(hkl)",
    "Crystal plane lying in the specimen plane of a named texture component.",
)
_MILLER_DIRECTION = SymbolUse(
    r"[uvw]",
    "Crystal direction along the specimen reference direction of a component.",
)

_ORIENTATION_CONCEPT = SeeAlso("Orientations and texture", "../../concepts/orientation_texture")
_API = SeeAlso("Core API", "../../api/index")


EQUAL_VOLUME_CHARTS = WorkedExample(
    id="core-orientation-equal-volume-charts-agree-on-the-so3-volume",
    title="The homochoric ball and the cubochoric cube enclose the same volume, and the cube's corner lands on the ball",
    domain="core",
    scenario=(
        "Uniform sampling of orientation space requires a chart whose volume "
        "element is the invariant measure of SO(3). Two such charts exist: the "
        "homochoric ball of radius (3*pi/4)^(1/3), and the cubochoric cube of "
        "edge pi^(2/3) that the equal-volume map carries onto it. Both must "
        "enclose the same volume, pi^2, and the map must send a corner of the "
        "cube exactly onto the surface of the ball. Either identity fails "
        "immediately for a wrong constant in the map, which is why they are "
        "checked here rather than trusted."
    ),
    setup=REPRESENTATION_SETUP,
    code=(
        "ball_volume = (4.0 / 3.0) * np.pi * HOMOCHORIC_BALL_RADIUS**3\n"
        "cube_volume = CUBOCHORIC_CUBE_EDGE**3\n"
        "corner = np.full((1, 3), CUBOCHORIC_CUBE_HALF_EDGE)\n"
        "corner_radius = float(np.linalg.norm(\n"
        "    convert_orientations(corner, source='cubochoric', target='homochoric')\n"
        "))\n"
        "result = np.array([\n"
        "    ball_volume / np.pi**2,\n"
        "    cube_volume / np.pi**2,\n"
        "    corner_radius / HOMOCHORIC_BALL_RADIUS,\n"
        "])"
    ),
    expected=[1.0, 1.0, 1.0],
    unit="",
    tolerance=1e-12,
    reference=(
        "The invariant measure on SO(3), (1 - cos w) dw dOmega / pi^2, gives "
        "the group total volume pi^2 before normalization. The homochoric "
        "radial function f(w) = [3(w - sin w)/4]^(1/3) reaches "
        "R1 = (3 pi/4)^(1/3) at w = pi, so (4/3) pi R1^3 = pi^2; the "
        "cubochoric edge is fixed as pi^(2/3) by requiring the same volume. A "
        "cube corner is at maximum distance from the centre in the cube, so "
        "the volume-preserving map must send it to the ball's surface."
    ),
    citation=(
        "Rosca, Morawiec and De Graef, Modelling Simul. Mater. Sci. Eng. 22 "
        "(2014) 075013, doi:10.1088/0965-0393/22/7/075013; Morawiec, "
        "Orientations and Rotations (Springer, 2004), for the invariant measure."
    ),
    symbols=(_HOMOCHORIC, _CUBOCHORIC, _BALL_RADIUS),
    see_also=(_ORIENTATION_CONCEPT, _API),
    result_format="{:.12f}",
)


IDEAL_INDICES_INVERT = WorkedExample(
    id="core-orientation-ideal-indices-invert-the-plane-direction-construction",
    title="Naming an orientation (hkl)[uvw] inverts the construction that built it",
    domain="core",
    scenario=(
        "A rolling-texture component is named by the crystal plane lying in "
        "the sheet plane and the crystal direction along the rolling "
        "direction. Orientation.from_miller turns that name into an "
        "orientation; ideal_orientation_indices turns an orientation back into "
        "the name. Round-tripping the copper component {112}<111> must return "
        "the indices it was built from, and both deviation angles must vanish "
        "- the deviations being what distinguishes an exact component from a "
        "nearest label."
    ),
    setup=REPRESENTATION_SETUP,
    code=(
        "copper = Orientation.from_miller(\n"
        "    (1, 1, 2), (1, 1, -1), phase=phase, specimen_frame=specimen\n"
        ")\n"
        "indices = ideal_orientation_indices(copper)\n"
        "result = np.array([\n"
        "    *indices.hkl,\n"
        "    *indices.uvw,\n"
        "    indices.plane_deviation_deg,\n"
        "    indices.direction_deviation_deg,\n"
        "])"
    ),
    expected=[1.0, 1.0, 2.0, 1.0, 1.0, -1.0, 0.0, 0.0],
    unit="",
    tolerance=1e-9,
    reference=(
        "An exact inverse identity, not a fitted result: the plane normal "
        "aligned with ND and the direction aligned with RD are recovered by "
        "mapping those specimen axes back through g^T and expressing them in "
        "the reciprocal and direct bases respectively, so the integer indices "
        "and zero residual angles follow by construction."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982), for the "
        "(hkl)[uvw] naming of rolling-texture components; Hirsch et al., "
        "Electron Microscopy of Thin Crystals (1965), for the copper component."
    ),
    symbols=(_MILLER_PLANE, _MILLER_DIRECTION),
    see_also=(_ORIENTATION_CONCEPT, _API),
    result_format="{:.9f}",
)


GROUP = ExampleGroup(
    slug="orientation-representations",
    title="Orientation representations",
    summary=(
        "The constants and identities behind the equal-volume charts of "
        "SO(3), and the inversion that names an orientation as a "
        "(hkl)[uvw] texture component."
    ),
    examples=(EQUAL_VOLUME_CHARTS, IDEAL_INDICES_INVERT),
)
