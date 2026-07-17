"""Worked examples: orientation-relationship index correspondence.

These examples exercise the index-correspondence surface of
``OrientationRelationship``: mapping parent-phase Miller planes and directions
to their child-phase counterparts. The reference values are the defining
parallelisms of the named relationships themselves (analytic identities): the
Kurdjumov-Sachs construction fixes ``{111}_fcc || {011}_bcc``, and the Bain
correspondence fixes ``[110]_fcc || [100]_bcc``, so the mapped rational indices
and their angular residuals are known exactly.

See :doc:`../../concepts/orientation_relationships` and the theory note on
index correspondence.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

CORRESPONDENCE_SETUP = """
import numpy as np
from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationRelationship,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

parent_frame = ReferenceFrame(
    name="austenite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
child_frame = ReferenceFrame(
    name="ferrite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
ferrite = Phase(
    "ferrite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
"""

_M_DIRECT = SymbolUse(
    r"\mathbf{M}",
    "Direction-index correspondence matrix mapping parent [uvw] to child [uvw].",
)
_M_RECIP = SymbolUse(
    r"\mathbf{M}^{*}",
    "Plane-index correspondence matrix mapping parent (hkl) to child (hkl).",
)
_HKL = SymbolUse(r"(hkl)", "Miller plane indices.")
_UVW = SymbolUse(r"[uvw]", "Miller direction indices.")

_OR_CONCEPT = SeeAlso(
    "Orientation relationships", "../../concepts/orientation_relationships"
)
_API = SeeAlso("Transformation API", "../../api/index")


KS_PLANE_CORRESPONDENCE = WorkedExample(
    id="or-ks-plane-correspondence-identity",
    title="Kurdjumov-Sachs maps (111) austenite onto (011) ferrite exactly",
    domain="transformation",
    scenario=(
        "Given the Kurdjumov-Sachs relationship, find which ferrite plane "
        "corresponds to the austenite close-packed plane (111). Because "
        "{111}_fcc || {011}_bcc is the defining parallelism of the "
        "relationship, the mapped plane must rationalize to (011) with zero "
        "angular residual — the residual is the verifiable quantity."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=(
        "ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(\n"
        "    parent_phase=austenite, child_phase=ferrite\n"
        ")\n"
        "mapped = ks.map_plane_to_child(\n"
        "    CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=austenite), phase=austenite)\n"
        ")\n"
        "result = np.concatenate(\n"
        "    [mapped.rational_indices.astype(float), [mapped.angular_residual_deg]]\n"
        ")"
    ),
    expected=[0.0, 1.0, 1.0, 0.0],
    unit="indices, deg",
    tolerance=1e-9,
    reference=(
        "The Kurdjumov-Sachs relationship is constructed from the parallelism "
        "{111}_fcc || {011}_bcc, so mapping the defining parent plane must "
        "recover the defining child plane identically (analytic identity)."
    ),
    citation="Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.",
    symbols=(_M_RECIP, _HKL),
    see_also=(_OR_CONCEPT, _API),
)


BAIN_DIRECTION_CORRESPONDENCE = WorkedExample(
    id="or-bain-direction-correspondence-identity",
    title="Bain maps [110] austenite onto [100] ferrite exactly",
    domain="transformation",
    scenario=(
        "Given the Bain correspondence, find which ferrite direction "
        "corresponds to the austenite [110] direction. The Bain construction "
        "fixes [110]_fcc || [100]_bcc, so the mapped direction must "
        "rationalize to [100] with zero angular residual."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=(
        "bain = OrientationRelationship.from_bain_correspondence(\n"
        "    parent_phase=austenite, child_phase=ferrite\n"
        ")\n"
        "mapped = bain.map_direction_to_child(\n"
        "    CrystalDirection([1.0, 1.0, 0.0], phase=austenite)\n"
        ")\n"
        "result = np.concatenate(\n"
        "    [mapped.rational_indices.astype(float), [mapped.angular_residual_deg]]\n"
        ")"
    ),
    expected=[1.0, 0.0, 0.0, 0.0],
    unit="indices, deg",
    tolerance=1e-9,
    reference=(
        "The Bain correspondence is constructed from (001)_fcc || (001)_bcc "
        "with [110]_fcc || [100]_bcc, so mapping the defining parent direction "
        "must recover the defining child direction identically (analytic identity)."
    ),
    citation="Bain, Trans. AIME 70 (1924) 25.",
    symbols=(_M_DIRECT, _UVW),
    see_also=(_OR_CONCEPT, _API),
)


_OMEGA = SymbolUse(
    r"(\mathbf{n}, \omega)",
    "Axis-angle pair of the symmetry-reduced misorientation representative.",
)


KS_MISORIENTATION_REPRESENTATION = WorkedExample(
    id="or-ks-misorientation-representation",
    title="Kurdjumov-Sachs as a misorientation: 42.85 deg about <0.968 0.178 0.178>",
    domain="transformation",
    scenario=(
        "Express the Kurdjumov-Sachs relationship the way it is measured from "
        "EBSD boundary data: as the minimal-angle symmetry-reduced "
        "misorientation. The published representative is a rotation of "
        "42.85 deg about an axis with components <0.968 0.178 0.178>; the "
        "computed angle and sorted absolute axis components are compared "
        "against that tabulated value."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=(
        "ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(\n"
        "    parent_phase=austenite, child_phase=ferrite\n"
        ")\n"
        "misorientation = ks.misorientation()\n"
        "axis = np.sort(np.abs(misorientation.rotation.axis))[::-1]\n"
        "result = np.concatenate([[misorientation.angle_deg], axis])"
    ),
    expected=[42.85, 0.9679, 0.1776, 0.1776],
    unit="deg, axis components",
    tolerance=5e-3,
    reference=(
        "The Kurdjumov-Sachs disorientation representative is tabulated as a "
        "42.85 deg rotation about <0.968 0.178 0.178> in standard "
        "thermo-mechanical processing references."
    ),
    citation=(
        "Verlinden, Driver, Samajdar, Doherty, Thermo-Mechanical Processing of "
        "Metallic Materials (2007); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325."
    ),
    symbols=(_OMEGA,),
    see_also=(_OR_CONCEPT, _API),
    result_format="{:.4f}",
)


OR_FITTING_SETUP = CORRESPONDENCE_SETUP + """
from pytex import (
    Orientation,
    OrientationSet,
    fit_orientation_relationship,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
gt = OrientationRelationship.from_greninger_troiano_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
gt_variants = gt.generate_variants()
rng = np.random.default_rng(11)
eulers = rng.uniform(0.0, 60.0, size=(20, 3))
parents = OrientationSet.from_orientations(
    [
        Orientation.from_euler(
            *euler, specimen_frame=specimen, symmetry=austenite.symmetry, phase=austenite
        )
        for euler in eulers
    ]
)
picks = rng.integers(0, len(gt_variants), size=20)
children = OrientationSet(
    quaternions=np.stack(
        [
            parents[index]
            .rotation.compose(gt_variants[int(picks[index])].parent_to_child_rotation.inverse())
            .quaternion
            for index in range(20)
        ],
        axis=0,
    ),
    crystal_frame=ferrite.crystal_frame,
    specimen_frame=specimen,
    symmetry=ferrite.symmetry,
    phase=ferrite,
)
"""


OR_FITTING_RECOVERS_GT = WorkedExample(
    id="or-fit-recovers-gt-from-ks-nominal",
    title="OR fitting recovers Greninger-Troiano from a Kurdjumov-Sachs start",
    domain="transformation",
    scenario=(
        "Twenty parent/child pairs are generated with the Greninger-Troiano "
        "relationship, then fitted starting from a Kurdjumov-Sachs nominal. "
        "The fit must land on the operative relationship exactly (zero mean "
        "residual) while reporting the documented 2.40 deg separation between "
        "the fitted relationship and the assumed KS nominal."
    ),
    setup=OR_FITTING_SETUP,
    code=(
        "ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(\n"
        "    parent_phase=austenite, child_phase=ferrite\n"
        ")\n"
        "report = fit_orientation_relationship(parents, children, ks)\n"
        "result = np.array([report.deviation_from_nominal_deg, report.mean_residual_deg])"
    ),
    expected=[2.4037, 0.0],
    unit="deg",
    tolerance=5e-3,
    reference=(
        "Exact GT-generated pairs must refit GT identically (zero residual is "
        "an analytic identity), and the reported distance from the KS nominal "
        "is the documented KS-GT representative separation of 2.40 deg."
    ),
    citation=(
        "Greninger and Troiano, Trans. AIME 185 (1949) 590; "
        "Kurdjumov and Sachs, Z. Phys. 64 (1930) 325."
    ),
    symbols=(_OMEGA,),
    see_also=(_OR_CONCEPT, _API),
)


GROUP = ExampleGroup(
    slug="transformation",
    title="Orientation-relationship correspondence",
    summary=(
        "Index-correspondence identities for named orientation relationships: "
        "mapping parent planes and directions to their product-phase "
        "counterparts, with rationalized indices and angular residuals, and "
        "the misorientation representation used for EBSD comparison."
    ),
    examples=(
        KS_PLANE_CORRESPONDENCE,
        BAIN_DIRECTION_CORRESPONDENCE,
        KS_MISORIENTATION_REPRESENTATION,
        OR_FITTING_RECOVERS_GT,
    ),
)
