"""Worked examples: the workbench service layer.

The application's service layer is JSON-in, JSON-out and knows nothing about
HTTP, so every capability the interface offers is callable directly and is
testable as a numerical surface like any other.

These examples validate the three quantitative claims the workbench user guide
makes, each against a value fixed independently of this code:

* the Kurdjumov-Sachs packet structure — 24 variants in 4 packets of 6, from
  Morito et al.;
* the ten distinct intervariant disorientations of the same relationship, from
  the same table;
* the closure of the m.r.d. scale — that the area-weighted mean of any correctly
  normalised pole figure is exactly 1, whatever the texture.

The third is the important one to have here rather than only in a unit test:
it is an *identity*, not a measurement, and it is the property that makes two
pole figures from different instruments comparable at all.

See ``docs/site/workflows/workbench_application.md`` and
``docs/architecture/application_platform.md``.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

REGISTRY_SETUP = """
from pytex.app import REGISTRY
"""

VARIANT_SETUP = (
    REGISTRY_SETUP
    + """
AUSTENITE_TO_FERRITE = {
    "phase": {"builtin": "austenite_fcc"},
    "child_phase": {"builtin": "fe_bcc"},
    "relationship": "kurdjumov_sachs",
    "packet_plane": [1, 1, 1],
}
"""
)

TEXTURE_SETUP = (
    REGISTRY_SETUP
    + """
def pole_figure(model):
    return REGISTRY.call(
        "texture.pole_figure",
        {
            "phase": {"builtin": "ni_fcc"},
            "model": model,
            "spread_deg": 10.0,
            "grain_count": 400,
            "halfwidth_deg": 10.0,
            "seed": 7,
            "pole": [1, 1, 1],
            "projection": "equal_area",
            "resolution_deg": 5.0,
        },
    )
"""
)

_MRD = SymbolUse("m.r.d.", "Multiples of a random distribution.")

_OR_CONCEPT = SeeAlso("Orientation relationships", "../../concepts/orientation_relationships")
_TEXTURE_CONCEPT = SeeAlso("Texture foundation", "../../concepts/texture_foundation")
_MRD_NOTE = SeeAlso(
    "Pole figure arithmetic and the m.r.d. scale",
    "../../theory/pole_figure_arithmetic_and_mrd",
)
_GUIDE = SeeAlso("The PyTex Workbench", "../../workflows/workbench_application")

_CITATION_MORITO = (
    "Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789, Table 2."
)


KS_PACKET_COUNT = WorkedExample(
    id="workbench-ks-packet-size",
    title="Kurdjumov-Sachs gives four packets of six variants",
    domain="application",
    scenario=(
        "One austenite grain transforming under Kurdjumov-Sachs produces 24 child orientations, "
        "and they are not unstructured: each variant carries exactly one member of the parent "
        "{111} family into exact parallelism with a child {110}, and the variants sharing that "
        "member form a packet. The {111} family has four members, so there are four packets, and "
        "24 variants divided among them evenly gives six each. That grouping is what a lath "
        "martensite micrograph shows as a block, and it is why one parent grain gives 24 "
        "orientations but only four apparent plate directions."
    ),
    setup=VARIANT_SETUP,
    code=(
        "response = REGISTRY.call(\n"
        "    'variants.pole_figure',\n"
        "    dict(AUSTENITE_TO_FERRITE, pole=[1, 0, 0], projection='stereographic',\n"
        "         include_parent=False),\n"
        ")\n"
        "sizes = response['data']['packet_sizes']\n"
        "result = sorted(sizes.values())"
    ),
    expected=[6, 6, 6, 6],
    unit="variants per packet",
    tolerance=0,
    reference=(
        "Morito et al. report the 24 Kurdjumov-Sachs variants of lath martensite as four packets "
        "of six, one packet per member of the parent {111} family."
    ),
    citation=_CITATION_MORITO,
    symbols=(),
    see_also=(_OR_CONCEPT, _GUIDE),
)


KS_DISORIENTATION_SPECTRUM = WorkedExample(
    id="workbench-ks-intervariant-spectrum",
    title="The 276 variant pairs fall on ten disorientations",
    domain="application",
    scenario=(
        "Two child grains that grew from the same parent cannot meet at an arbitrary "
        "misorientation: the admissible set is fixed by the relationship and the two point "
        "groups, and it is discrete. That discreteness is what makes a measured misorientation "
        "histogram a test — peaks away from these ten angles are boundaries between different "
        "parent grains, which is the reasoning parent-grain reconstruction rests on. The 24 "
        "variants give 276 unordered pairs, and they collapse onto ten values."
    ),
    setup=VARIANT_SETUP,
    code=(
        "response = REGISTRY.call(\n"
        "    'variants.intervariant_misorientations',\n"
        "    dict(AUSTENITE_TO_FERRITE, merge_equal_angles=True),\n"
        ")\n"
        "result = [row['angle_deg'] for row in response['table']['rows']]"
    ),
    expected=[10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00],
    unit="deg",
    tolerance=0.01,
    reference=(
        "The ten distinct intervariant disorientation angles of the Kurdjumov-Sachs relationship, "
        "as tabulated for lath martensite by Morito et al."
    ),
    citation=_CITATION_MORITO,
    symbols=(),
    see_also=(_OR_CONCEPT, _GUIDE),
    result_format="{:.2f}",
)


MRD_MEAN_IS_ONE = WorkedExample(
    id="workbench-mrd-mean-is-one",
    title="The area-weighted mean of any pole figure is 1 m.r.d.",
    domain="application",
    scenario=(
        "Multiples of a random distribution is not a display convenience but a normalisation with "
        "an exact consequence: since 1 m.r.d. is by definition what a texture-free material gives "
        "everywhere, the area-weighted mean over the hemisphere must be 1 for *every* texture, "
        "however sharp. The mean therefore carries no information, which is precisely why every "
        "feature of a pole figure is a departure from it — and why a figure whose mean is not 1 "
        "has not been normalised, and its numbers mean nothing outside itself.\n\n"
        "The identity is checked here on a sharp single component and on a random texture at "
        "once, because holding for both is the whole claim. Note that it is an *area*-weighted "
        "mean: the unweighted average of the same grid is biased, because an equispaced grid on a "
        "hemisphere is not an equal-area one."
    ),
    setup=TEXTURE_SETUP,
    code=(
        "result = [\n"
        "    pole_figure('random')['data']['mean_mrd'],\n"
        "    pole_figure('goss')['data']['mean_mrd'],\n"
        "    pole_figure('fcc_rolling')['data']['mean_mrd'],\n"
        "]"
    ),
    expected=[1.0, 1.0, 1.0],
    unit="m.r.d.",
    tolerance=0.01,
    reference=(
        "Definitional: a pole figure normalised to multiples of a random distribution integrates "
        "to the sphere's area, so its area-weighted mean is exactly 1 for any texture."
    ),
    citation=(
        "Randle and Engler, Introduction to Texture Analysis, 2nd ed., chapter 5 "
        "(pole figure normalisation)."
    ),
    symbols=(_MRD,),
    see_also=(_TEXTURE_CONCEPT, _MRD_NOTE, _GUIDE),
    result_format="{:.3f}",
)


GOSS_POLE_AT_CENTRE = WorkedExample(
    id="workbench-goss-pole-at-nd",
    title="A Miller label is a testable claim: Goss puts (011) on ND",
    domain="application",
    scenario=(
        "The Goss component is written {011}<100>, and the first half of that notation asserts "
        "that the {011} plane lies in the sheet plane. So the (011) pole of a Goss texture must "
        "point along the sheet normal — polar angle zero, the centre of the pole figure. Checking "
        "it needs no reference figure at all, only the notation, which makes it the sharpest "
        "available end-to-end test of the whole chain: the Euler convention, the crystal-to-"
        "specimen mapping, the symmetry family of the pole, and the projection."
    ),
    setup=TEXTURE_SETUP.replace("'model': model", "'model': model").replace(
        '"spread_deg": 10.0', '"spread_deg": 6.0'
    ),
    code=(
        "rows = REGISTRY.call(\n"
        "    'texture.pole_figure',\n"
        "    {\n"
        "        'phase': {'builtin': 'ni_fcc'}, 'model': 'goss', 'spread_deg': 6.0,\n"
        "        'grain_count': 400, 'halfwidth_deg': 10.0, 'seed': 7,\n"
        "        'pole': [0, 1, 1], 'projection': 'equal_area', 'resolution_deg': 5.0,\n"
        "    },\n"
        ")['table']['rows']\n"
        "result = max(rows, key=lambda row: row['mrd'])['polar_deg']"
    ),
    expected=0.0,
    unit="deg",
    tolerance=6.0,
    reference=(
        "The Goss component {011}<100> places {011} in the sheet plane, so its (011) pole lies "
        "along ND at a polar angle of zero. The tolerance is the grid spacing, not a fitted "
        "margin."
    ),
    citation=(
        "Randle and Engler, Introduction to Texture Analysis, 2nd ed., chapter 5 "
        "(ideal orientations and their Miller descriptions)."
    ),
    symbols=(),
    see_also=(_TEXTURE_CONCEPT, _GUIDE),
    result_format="{:.1f}",
)



CRYSTAL_VIEWER_GOSS_ND = WorkedExample(
    id="workbench-crystal-viewer-goss-nd",
    title="The crystal viewer's camera is an orientation: Goss puts ND 45 degrees from [001]",
    domain="application",
    scenario=(
        "The crystal viewer turns a structure by accumulating a rotation from the drag, and the "
        "orientation dock claims that rotation *is* an orientation in the crystal-to-specimen "
        "convention — so that entering Euler angles and turning the crystal by hand are the same "
        "act. The claim is testable from notation alone. Goss is {011}<100>, so its sheet normal "
        "is a <011>, and the angle between <011> and <001> is arccos(1/sqrt(2)). Were the camera "
        "a transpose away from an orientation, or the Euler convention the other handedness, "
        "this number would not be 45 degrees."
    ),
    setup=REGISTRY_SETUP,
    code=(
        "rows = REGISTRY.call(\n"
        "    'crystal.orientation',\n"
        "    {\n"
        "        'phase': {'builtin': 'cu_fcc'}, 'euler_convention': 'bunge',\n"
        "        'angle1': 0.0, 'angle2': 45.0, 'angle3': 0.0,\n"
        "    },\n"
        ")['table']['rows']\n"
        "result = next(row['polar_deg'] for row in rows if row['axis'] == 'ND')"
    ),
    expected=45.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "Goss is {011}<100>: the sheet normal is a <011> direction, and the angle between <011> "
        "and the crystal c axis <001> is arccos(1/sqrt(2)) = 45 degrees exactly, whatever the "
        "lattice parameter. The tolerance is machine precision, not a margin."
    ),
    citation=(
        "Randle and Engler, Introduction to Texture Analysis, 2nd ed., chapter 5 "
        "(ideal orientations and their Miller descriptions); Bunge, Texture Analysis in "
        "Materials Science (1982), for the Euler convention."
    ),
    symbols=(),
    see_also=(_TEXTURE_CONCEPT, _GUIDE),
    result_format="{:.4f}",
)


CRYSTAL_VIEWER_EULER_ROUND_TRIP = WorkedExample(
    id="workbench-crystal-viewer-euler-round-trip",
    title="The camera and the angle triple are the same orientation, both ways",
    domain="application",
    scenario=(
        "The dock reads Euler angles off the camera while you drag and writes them back when you "
        "type a triple, so the two conversions must be exact inverses: a triple that survives the "
        "trip through nine matrix entries and back is what makes 'set the view to brass' and "
        "'what view is this?' answers to the same question. Brass, {011}<211>, supplies the "
        "triple, because its first angle is the irrational arctan(1/sqrt(2)) rather than a round "
        "number — a decomposition that quietly snapped to a grid would show."
    ),
    setup=(
        "import numpy as np\n"
        "from pytex.app.services.crystal import (\n"
        "    camera_matrix_from_euler,\n"
        "    euler_from_camera_matrix,\n"
        ")\n"
    ),
    code=(
        "brass = (35.26438968275465, 45.0, 0.0)\n"
        "camera = camera_matrix_from_euler(*brass)\n"
        "recovered = np.asarray(euler_from_camera_matrix(camera))\n"
        "result = float(np.max(np.abs(recovered - np.asarray(brass))))"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-6,
    reference=(
        "A rotation matrix and its Euler decomposition name the same element of SO(3), so the "
        "round trip is the identity up to floating-point error. Brass is "
        "(arctan(1/sqrt(2)), 45, 0) degrees in Bunge angles. The tolerance is a micro-degree "
        "rather than machine epsilon because the reader deliberately rounds to a picodegree "
        "before wrapping into [0, 360), so that an angle landing a hair below zero is reported "
        "as zero rather than as a full turn."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science, Butterworths (1982), chapter 2 "
        "(the ZXZ Euler convention and its matrix)."
    ),
    symbols=(),
    see_also=(_TEXTURE_CONCEPT, _GUIDE),
    result_format="{:.2e}",
)

GROUP = ExampleGroup(
    slug="workbench-service-layer",
    title="Workbench service layer",
    summary=(
        "The three quantitative claims the workbench user guide makes, each checked against a "
        "value fixed independently of this code: the Kurdjumov-Sachs packet structure and "
        "intervariant spectrum from Morito et al., the closure of the m.r.d. scale as an exact "
        "identity, the assertion a Miller component label makes about where its poles land, and the "
        "crystal viewer's claim that its camera is an orientation."
    ),
    examples=(
        KS_PACKET_COUNT,
        KS_DISORIENTATION_SPECTRUM,
        MRD_MEAN_IS_ONE,
        GOSS_POLE_AT_CENTRE,
        CRYSTAL_VIEWER_GOSS_ND,
        CRYSTAL_VIEWER_EULER_ROUND_TRIP,
    ),
)

__all__ = ["GROUP"]
