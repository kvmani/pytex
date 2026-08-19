"""Worked examples: interplanar/interdirection angles, spacings, multiplicities.

Every value below is computed live from the public PyTex API and compared with
an analytic identity that holds for the relevant crystal system. These are the
canonical "sanity" quantities that any texture or diffraction practitioner
checks first, which makes them ideal executable documentation.

See the theory note :doc:`../../theory/index` and the concept page
:doc:`../../concepts/miller_planes_directions` for the underlying metric-tensor
formulation.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

# Shared preambles. The setup is rendered as a collapsible block in the docs and
# executed verbatim before each snippet, so the phase objects are real.

CUBIC_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
"""

HEX_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
hexagonal = Phase(
    "hcp-demo",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
"""

_PLANE_ANGLE = SymbolUse(r"\angle(\mathbf{n}_1, \mathbf{n}_2)", "Angle between two plane normals.")
_DIR_ANGLE = SymbolUse(r"\angle(\mathbf{d}_1, \mathbf{d}_2)", "Angle between two lattice directions.")
_DSPACING = SymbolUse(r"d_{hkl}", "Interplanar spacing of the (hkl) family.")
_MULT = SymbolUse(r"m_{\{hkl\}}", "Symmetry multiplicity of a plane family.")

_MILLER_CONCEPT = SeeAlso("Miller planes and directions", "../../concepts/miller_planes_directions")
_MILLER_API = SeeAlso("angle_plane_plane_rad / angle_dir_dir_rad", "../../api/index")


CUBIC_PLANE_100_110 = WorkedExample(
    id="cubic-angle-100-110",
    title="Angle between (100) and (110) in a cubic crystal",
    domain="core",
    scenario=(
        "You have indexed two poles as {100} and {110} on a cubic phase and want to confirm the "
        "geometry of a pole figure or a Kikuchi band intersection. In a cubic system the answer is "
        "exactly 45 degrees, independent of the lattice parameter, so this is the first check that "
        "your frame and symmetry wiring is correct."
    ),
    setup=CUBIC_SETUP,
    code=(
        "result = float(np.degrees(angle_plane_plane_rad(\n"
        "    MillerPlane.from_hkl([1, 0, 0], phase=cubic),\n"
        "    MillerPlane.from_hkl([1, 1, 0], phase=cubic),\n"
        ")))"
    ),
    expected=45.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "For cubic metrics the interplanar angle is arccos(h1 h2 + k1 k2 + l1 l2 over norms); "
        "for (100) and (110) this is arccos(1/sqrt(2)) = 45 degrees, independent of a."
    ),
    citation="Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.",
    symbols=(_PLANE_ANGLE,),
    see_also=(_MILLER_CONCEPT, _MILLER_API),
)

CUBIC_PLANE_100_111 = WorkedExample(
    id="cubic-angle-100-111",
    title="Angle between (100) and (111) in a cubic crystal",
    domain="core",
    scenario=(
        "The 100-to-111 angle sets the classic stereographic-triangle geometry and appears whenever "
        "you relate a rolling-plane normal to an octahedral slip plane. The exact value is "
        "arccos(1/sqrt(3))."
    ),
    setup=CUBIC_SETUP,
    code=(
        "result = float(np.degrees(angle_plane_plane_rad(\n"
        "    MillerPlane.from_hkl([1, 0, 0], phase=cubic),\n"
        "    MillerPlane.from_hkl([1, 1, 1], phase=cubic),\n"
        ")))"
    ),
    expected=54.735610317245346,
    unit="deg",
    tolerance=1e-9,
    reference="arccos(1/sqrt(3)) = 54.7356 degrees for cubic (100) vs (111).",
    citation="Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.",
    symbols=(_PLANE_ANGLE,),
    see_also=(_MILLER_CONCEPT, _MILLER_API),
)

CUBIC_DIR_110_111 = WorkedExample(
    id="cubic-angle-dir-110-111",
    title="Angle between [110] and [111] directions in a cubic crystal",
    domain="core",
    scenario=(
        "Slip-system and Schmid-factor calculations repeatedly need the angle between a slip "
        "direction such as [110] and a loading or plane-normal direction such as [111]. In cubic "
        "metrics the direction angle equals the same-index plane angle."
    ),
    setup=CUBIC_SETUP,
    code=(
        "result = float(np.degrees(angle_dir_dir_rad(\n"
        "    MillerDirection.from_uvw([1, 1, 0], phase=cubic),\n"
        "    MillerDirection.from_uvw([1, 1, 1], phase=cubic),\n"
        ")))"
    ),
    expected=35.264389682754654,
    unit="deg",
    tolerance=1e-9,
    reference="arccos(sqrt(2/3)) = 35.2644 degrees for cubic [110] vs [111].",
    citation="Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.",
    symbols=(_DIR_ANGLE,),
    see_also=(_MILLER_CONCEPT, _MILLER_API),
)

CUBIC_DSPACING_111 = WorkedExample(
    id="cubic-dspacing-111",
    title="Interplanar spacing of (111) in a cubic crystal (a = 4 angstrom)",
    domain="core",
    scenario=(
        "Interplanar spacing is the bridge from crystallography to diffraction: it fixes the Bragg "
        "angle for a reflection. For a cubic lattice d_hkl = a / sqrt(h^2 + k^2 + l^2), so d_111 = "
        "a / sqrt(3)."
    ),
    setup=CUBIC_SETUP,
    code="result = MillerPlane.from_hkl([1, 1, 1], phase=cubic).d_spacing_angstrom",
    expected=2.3094010767585034,
    unit="angstrom",
    tolerance=1e-9,
    reference="d_111 = a / sqrt(3) = 4 / sqrt(3) = 2.30940 angstrom.",
    citation="Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Eq. 3-10.",
    symbols=(_DSPACING,),
    see_also=(_MILLER_CONCEPT, SeeAlso("Diffraction geometry worked examples", "diffraction")),
    result_format="{:.5f}",
)

CUBIC_MULTIPLICITY = WorkedExample(
    id="cubic-plane-multiplicity",
    title="Symmetry multiplicity of {100}, {110}, {111}, {321} under m-3m",
    domain="core",
    scenario=(
        "Powder-diffraction intensities and pole-figure normalization both depend on how many "
        "symmetry-equivalent planes a family contains. PyTex treats plane families with antipodal "
        "equivalence (a plane and its opposite normal are the same plane), so the reduced "
        "multiplicities are half the full point-group orbit: 3, 6, 4, and 24."
    ),
    setup=CUBIC_SETUP,
    code=(
        "result = [\n"
        "    MillerPlane.from_hkl(hkl, phase=cubic).symmetry_equivalent_indices()[0].shape[1]\n"
        "    for hkl in ([1, 0, 0], [1, 1, 0], [1, 1, 1], [3, 2, 1])\n"
        "]"
    ),
    expected=[3, 6, 4, 24],
    unit="",
    tolerance=0.0,
    reference=(
        "Full m-3m orbit sizes are 6, 12, 8, 48; antipodal folding halves them to 3, 6, 4, 24."
    ),
    citation="Hahn (ed.), International Tables for Crystallography Vol. A, point group m-3m.",
    symbols=(_MULT,),
    see_also=(
        SeeAlso("Symmetry and fundamental regions", "../../concepts/symmetry_and_fundamental_regions"),
        _MILLER_API,
    ),
)

HEX_BASAL_PRISM = WorkedExample(
    id="hex-angle-basal-prism",
    title="Angle between (0001) and (10-10) in a hexagonal crystal",
    domain="core",
    scenario=(
        "In HCP metals the basal plane (0001) and the prismatic planes {10-10} are the dominant "
        "slip and texture planes. Their normals are exactly perpendicular for any c/a ratio, which "
        "makes this a robust convention check for the hexagonal metric and the four-index handling."
    ),
    setup=HEX_SETUP,
    code=(
        "result = float(np.degrees(angle_plane_plane_rad(\n"
        "    MillerPlane.from_hkl([0, 0, 1], phase=hexagonal),\n"
        "    MillerPlane.from_hkl([1, 0, 0], phase=hexagonal),\n"
        ")))"
    ),
    expected=90.0,
    unit="deg",
    tolerance=1e-9,
    reference="The basal-plane normal is c*; prismatic normals lie in the basal plane, so the angle is 90 degrees.",
    citation="Partridge, The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967).",
    symbols=(_PLANE_ANGLE,),
    see_also=(
        SeeAlso("Hexagonal and trigonal conventions", "../../standards/hexagonal_and_trigonal_conventions"),
        _MILLER_API,
    ),
)

HEX_PRISM_PRISM = WorkedExample(
    id="hex-angle-prism-prism",
    title="Angle between (10-10) and (01-10) in a hexagonal crystal",
    domain="core",
    scenario=(
        "Adjacent first-order prismatic planes in HCP are separated by 60 degrees. This confirms "
        "that the 120-degree gamma angle of the hexagonal cell is handled correctly when reasoning "
        "about prismatic slip variants."
    ),
    setup=HEX_SETUP,
    code=(
        "result = float(np.degrees(angle_plane_plane_rad(\n"
        "    MillerPlane.from_hkl([1, 0, 0], phase=hexagonal),\n"
        "    MillerPlane.from_hkl([0, 1, 0], phase=hexagonal),\n"
        ")))"
    ),
    expected=60.0,
    unit="deg",
    tolerance=1e-9,
    reference="First-order prismatic normals are separated by 60 degrees in the hexagonal basal plane.",
    citation="Partridge, The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967).",
    symbols=(_PLANE_ANGLE,),
    see_also=(
        SeeAlso("Hexagonal and trigonal conventions", "../../standards/hexagonal_and_trigonal_conventions"),
        _MILLER_API,
    ),
)



CUBIC_NAME_A_DIRECTION = WorkedExample(
    id="cubic-nearest-low-index-round-trip",
    title="Naming a direction that arrived as geometry recovers its own indices",
    domain="core",
    scenario=(
        "An inverse pole figure, a stereogram readout and an interactive viewer all face the same "
        "problem: the direction is known as a vector and the reader wants a label. The naming is "
        "a search over low-index triples, so the first thing to establish is that it is exact "
        "where it should be — a direction built from [3 2 1] must be named [3 2 1], with an "
        "angular residual of zero rather than of the search's step size."
    ),
    setup=CUBIC_SETUP,
    code=(
        "vector = MillerDirection.from_uvw([3, 2, 1], phase=cubic).unit_vector_cartesian\n"
        "indices, residual_deg = nearest_low_index_direction(vector, phase=cubic)\n"
        "result = float(residual_deg)"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "The vector is the Cartesian image of [3 2 1] under the direct basis, and gcd(3, 2, 1) = 1 "
        "so the triple is already primitive. The nearest primitive triple to a direction that is "
        "one is therefore itself, at zero angle."
    ),
    citation=(
        "International Tables for Crystallography, Volume A, section 1.3 "
        "(direct and reciprocal bases and their index conventions)."
    ),
    symbols=(_DIR_ANGLE,),
    see_also=(_MILLER_CONCEPT, _MILLER_API),
)

HEX_NAME_A_PLANE_NORMAL = WorkedExample(
    id="hex-nearest-low-index-plane-vs-direction",
    title="In a hexagonal phase a plane normal and the like-indexed direction are not parallel",
    domain="core",
    scenario=(
        "Miller indices are reciprocal-basis components and direction indices are direct-basis "
        "components, so naming a vector as a plane and naming it as a direction are different "
        "questions with different answers in every non-cubic phase. The pyramidal normal is the "
        "standard demonstration: name the (10-11) normal as a plane and it comes back exactly, "
        "name it as a direction and the angle to [10-11] is the metric difference itself."
    ),
    setup=HEX_SETUP,
    code=(
        "normal = MillerPlane.from_hkl([1, 0, 1], phase=hexagonal).normal_cartesian\n"
        "indices, residual_deg = nearest_low_index_plane(normal, phase=hexagonal)\n"
        "result = [int(value) for value in indices]"
    ),
    expected=[1, 0, 1],
    unit="",
    tolerance=0.0,
    reference=(
        "The vector is the unit normal of (10-11) by construction, and the plane search runs "
        "against the reciprocal basis, so it must return the indices the normal was built from. "
        "The direction search, which runs against the direct basis, returns a different triple "
        "for the same vector because c/a is not 1."
    ),
    citation=(
        "International Tables for Crystallography, Volume A, section 1.3; Partridge, "
        "The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967)."
    ),
    symbols=(_PLANE_ANGLE,),
    see_also=(
        SeeAlso(
            "Hexagonal and trigonal conventions",
            "../../standards/hexagonal_and_trigonal_conventions",
        ),
        _MILLER_API,
    ),
    result_format="{:.0f}",
)

GROUP = ExampleGroup(
    slug="crystal_geometry",
    title="Crystal geometry: angles, spacings, and multiplicities",
    summary=(
        "Interplanar and interdirection angles, interplanar spacings, and symmetry multiplicities "
        "for cubic and hexagonal phases, and the naming of a direction or a plane that arrived "
        "as a Cartesian vector. Each result is checked against an analytic identity for the "
        "relevant crystal system."
    ),
    examples=(
        CUBIC_PLANE_100_110,
        CUBIC_PLANE_100_111,
        CUBIC_DIR_110_111,
        CUBIC_DSPACING_111,
        CUBIC_MULTIPLICITY,
        HEX_BASAL_PRISM,
        HEX_PRISM_PRISM,
        CUBIC_NAME_A_DIRECTION,
        HEX_NAME_A_PLANE_NORMAL,
    ),
)

__all__ = ["GROUP"]
