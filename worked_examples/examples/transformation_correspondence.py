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



_FINGERPRINT = SymbolUse(
    r"G_c \left(R G_p R^{\mathsf{T}}\right) G_c",
    "Same-parent boundary fingerprint: the admissible child-child "
    "misorientations of one parent grain.",
)


KS_SIGMA3_IS_AN_ADMISSIBLE_BOUNDARY = WorkedExample(
    id="or-ks-same-parent-boundary-fingerprint",
    title="The Sigma3 twin is an admissible Kurdjumov-Sachs same-parent boundary",
    domain="transformation",
    scenario=(
        "Deciding whether two neighbouring martensite grains descend from one "
        "austenite grain means asking whether their boundary misorientation is "
        "one the relationship can actually produce. That admissible set is "
        "``G_c (R G_p R^T) G_c``, because two children of one parent satisfy "
        "``C_i^T C_j = V_i V_j^T``. Two identities are checked: the published "
        "Kurdjumov-Sachs intervariant table contains a 60 deg rotation about "
        "<111> — the Sigma3 twin relation, Morito's V1-V20 pair — so the exact "
        "Sigma3 rotation must sit at zero distance from the fingerprint; and "
        "every one of the 276 distinct variant-pair boundaries of a common "
        "parent must sit at zero distance too, since they generate the set by "
        "construction."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=(
        "from pytex import (\n"
        "    Rotation,\n"
        "    boundary_fingerprint_distances_deg,\n"
        "    intervariant_boundary_fingerprint,\n"
        ")\n"
        "\n"
        "ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(\n"
        "    parent_phase=austenite, child_phase=ferrite\n"
        ")\n"
        "fingerprint = intervariant_boundary_fingerprint(ks)\n"
        "\n"
        "sigma3 = Rotation.from_axis_angle(\n"
        "    np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0), np.deg2rad(60.0)\n"
        ").as_matrix()\n"
        "sigma3_distance = float(\n"
        "    boundary_fingerprint_distances_deg(sigma3[None, :, :], fingerprint)[0]\n"
        ")\n"
        "\n"
        "variants = ks.generate_variants()\n"
        "children = np.stack(\n"
        "    [variant.parent_to_child_rotation.inverse().as_matrix() for variant in variants]\n"
        ")\n"
        "left, right = np.triu_indices(len(variants), k=1)\n"
        "boundaries = np.einsum(\n"
        "    'nji,njk->nik', children[left], children[right], optimize=True\n"
        ")\n"
        "worst_variant_pair = float(\n"
        "    boundary_fingerprint_distances_deg(boundaries, fingerprint).max()\n"
        ")\n"
        "result = [sigma3_distance, worst_variant_pair]"
    ),
    expected=[0.0, 0.0],
    unit="deg",
    tolerance=1e-5,
    reference=(
        "Both values are identities, not fitted numbers. The Kurdjumov-Sachs "
        "intervariant table published by Morito et al. lists a 60 deg / <111> "
        "variant pair (V1-V20), which is exactly the Sigma3 coincidence-site "
        "relation, so the Sigma3 rotation belongs to the admissible set. The "
        "variant-pair boundaries generate the set by construction, so their "
        "distance to it is identically zero. The 1e-5 deg tolerance is the "
        "arccos and quaternion/matrix round-trip noise floor, not a physical "
        "margin."
    ),
    citation=(
        "Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) "
        "1789 (KS intervariant table); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325."
    ),
    symbols=(_FINGERPRINT,),
    see_also=(_OR_CONCEPT, _API),
)

_R_FIT = SymbolUse(
    r"\mathbf{R}",
    "Parent-to-child rotation of an orientation relationship.",
)

_IDENTIFICATION_CODE = """
from pytex import (
    OrientationSet,
    Rotation,
    characterize_orientation_relationship,
    specimen_frame,
)

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
variants = ks.generate_variants()
parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
# Canonical crystal->specimen convention: C = P V^T.
child_matrices = np.stack(
    [
        parent_matrix @ variants[k].parent_to_child_rotation.as_matrix().T
        for k in (0, 4, 8, 13, 17, 22)
    ]
)
frame = specimen_frame()
parents = OrientationSet.from_matrices(
    np.stack([parent_matrix] * 6), specimen_frame=frame, phase=austenite
)
children = OrientationSet.from_matrices(
    child_matrices, specimen_frame=frame, phase=ferrite
)
report = characterize_orientation_relationship(parents, children)
deviations = dict(zip(report.catalog_names, report.catalog_deviations_deg, strict=True))
result = [deviations["kurdjumov_sachs"], deviations["nishiyama_wassermann"]]
""".strip()


KS_IDENTIFIED_FROM_MEASURED_ORIENTATIONS = WorkedExample(
    id="or-ks-identified-from-measured-orientations",
    title="Kurdjumov-Sachs recovered from measured parent/child orientation pairs",
    domain="transformation",
    scenario=(
        "The everyday EBSD question: a parent grain and several child grains "
        "were indexed, and the operative orientation relationship is wanted. "
        "Children are synthesized here through six known Kurdjumov-Sachs "
        "variants of one parent, and characterization runs with no nominal "
        "relationship supplied, so the answer comes from the data alone. Two "
        "quantities are checked: the deviation of the fitted rotation from "
        "catalog Kurdjumov-Sachs, and its deviation from Nishiyama-Wassermann. "
        "The first must be zero because the data were built from that "
        "relationship; the second must be the published separation between the "
        "two relationships, which is what makes them distinguishable at all."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=_IDENTIFICATION_CODE,
    expected=[0.0, 5.26],
    unit="deg",
    tolerance=0.01,
    reference=(
        "The first value is an analytic identity: the children were generated "
        "from exact Kurdjumov-Sachs variants, so the fitted rotation must "
        "coincide with the relationship it was built from. The second is the "
        "tabulated 5.26 deg angular separation between the Kurdjumov-Sachs and "
        "Nishiyama-Wassermann relationships. Neither number is copied from a "
        "previous program output."
    ),
    citation=(
        "Kurdjumov and Sachs, Z. Phys. 64 (1930) 325; Nishiyama, Sci. Rep. "
        "Tohoku Univ. 23 (1934) 637; Wassermann, Arch. Eisenhuettenwes. 16 "
        "(1933) 647."
    ),
    symbols=(_R_FIT, _OMEGA),
    see_also=(_OR_CONCEPT, _API),
)

_STATEMENT_CODE = """
from pytex import describe_orientation_relationship

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
planes, directions = describe_orientation_relationship(ks)
plane, direction = planes[0], directions[0]
result = np.concatenate(
    [
        np.sort(np.abs(plane.parent_indices)),
        np.sort(np.abs(plane.child_indices)),
        np.sort(np.abs(direction.parent_indices)),
        np.sort(np.abs(direction.child_indices)),
        [plane.deviation_deg, direction.deviation_deg],
    ]
)
""".strip()


KS_STATEMENT_IS_RECOVERED_FROM_THE_ROTATION = WorkedExample(
    id="or-ks-parallelism-statement-from-rotation",
    title="Reading the Kurdjumov-Sachs parallelisms back out of its rotation",
    domain="transformation",
    scenario=(
        "An orientation relationship is stored as a rotation, but the "
        "literature reports it as parallel planes and directions. This example "
        "recovers that statement from the rotation alone and checks it against "
        "the defining Kurdjumov-Sachs parallelisms: the parent plane must "
        "belong to {111} and its child image to {011}, the parent direction to "
        "<110> and its child image to <111>, all at zero angular deviation. "
        "Sorted absolute indices are compared because any member of a family is "
        "an equally correct statement of the same relationship."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=_STATEMENT_CODE,
    expected=[1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0.0, 0.0],
    unit="indices, deg",
    tolerance=1e-4,
    reference=(
        "The Kurdjumov-Sachs relationship is defined by {111}_fcc || {011}_bcc "
        "and <110>_fcc || <111>_bcc, so recovering the statement from the "
        "rotation must reproduce exactly those families at zero deviation "
        "(analytic identity). The 1e-4 deg tolerance is the matrix-quaternion "
        "round-trip noise floor, not a physical margin."
    ),
    citation="Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.",
    symbols=(_R_FIT, _HKL, _UVW),
    see_also=(_OR_CONCEPT, _API),
)

_VARIANT_TABLE_CODE = """
from pytex import variant_correspondence_table

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
table = variant_correspondence_table(
    ks, CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=austenite), phase=austenite)
)
exact = table.exact_rows()
group_sizes = {}
for row in table.rows:
    group_sizes[row.equivalence_group] = group_sizes.get(row.equivalence_group, 0) + 1
result = [
    len(table.rows),
    table.distinct_image_count((1, 1, 1)),
    len(exact),
    min(group_sizes.values()),
    max(group_sizes.values()),
]
""".strip()


KS_VARIANT_CORRESPONDENCE_TABLE = WorkedExample(
    id="or-ks-variant-correspondence-packets",
    title="The (111) variant correspondence table is the four Kurdjumov-Sachs packets",
    domain="transformation",
    scenario=(
        "Ask what one austenite plane becomes in every martensite variant. "
        "Mapping (111) through all 24 Kurdjumov-Sachs variants and grouping the "
        "images by index family must reproduce the packet structure of lath "
        "martensite: four crystallographically distinct answers, six variants "
        "each, of which exactly one group — six variants — carries (111) onto a "
        "{011} ferrite plane at zero residual. The computed values are the row "
        "count, the number of distinct images, the number of exactly parallel "
        "variants, and the smallest and largest group sizes."
    ),
    setup=CORRESPONDENCE_SETUP,
    code=_VARIANT_TABLE_CODE,
    expected=[24, 4, 6, 6, 6],
    unit="counts",
    tolerance=0,
    reference=(
        "Crystallographic identity, not a measured coincidence. Kurdjumov-Sachs "
        "has 24 variants; each carries exactly one member of the four-member "
        "{111} family onto its {011} close-packed child plane, so any nominated "
        "member is the close-packed plane of exactly 24/4 = 6 of them. Those six "
        "are one packet in the sense of Morito et al., and the parent symmetry "
        "acts transitively on the remaining images, so every group holds six."
    ),
    citation=(
        "Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) "
        "1789 (packet structure); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325."
    ),
    symbols=(_M_RECIP, _HKL),
    see_also=(_OR_CONCEPT, _API),
)

GROUP = ExampleGroup(
    slug="transformation",
    title="Orientation-relationship correspondence",
    summary=(
        "Index-correspondence identities for named orientation relationships: "
        "mapping parent planes and directions to their product-phase "
        "counterparts, with rationalized indices and angular residuals, the "
        "misorientation representation used for EBSD comparison, and the "
        "recovery of a relationship and its parallelism statement from measured "
        "parent/child orientation pairs."
    ),
    examples=(
        KS_PLANE_CORRESPONDENCE,
        BAIN_DIRECTION_CORRESPONDENCE,
        KS_MISORIENTATION_REPRESENTATION,
        OR_FITTING_RECOVERS_GT,
        KS_SIGMA3_IS_AN_ADMISSIBLE_BOUNDARY,
        KS_IDENTIFIED_FROM_MEASURED_ORIENTATIONS,
        KS_STATEMENT_IS_RECOVERED_FROM_THE_ROTATION,
        KS_VARIANT_CORRESPONDENCE_TABLE,
    ),
)
