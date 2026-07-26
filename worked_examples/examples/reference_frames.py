"""Worked examples: reference-frame creation, transformation, and resolution.

Reference frames are the foundation every other PyTex quantity rests on, so the
examples here check the properties a user actually relies on: that a declared
axis correspondence produces the rotation it claims, that the frame graph
composes a multi-hop chain correctly, and that a round trip through a chain
returns the original components.

Every expected value has independent provenance — an exact rotation-matrix
identity, a right-handed-basis convention from the International Tables, or an
angle fixed by the geometry of the declaration itself — never a copied prior
program output.

See the concept page :doc:`../../concepts/reference_frames_and_conventions` and
the architecture note on the reference-frame foundation for the underlying
model.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

# Shared preamble. Rendered as a collapsible block in the docs and executed
# verbatim before each snippet, so every object below is real.

FRAME_SETUP = """
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
"""

_ROTATION = SymbolUse(r"\mathbf{R}", "Rotation matrix mapping source-frame to target-frame components.")
_ANGLE = SymbolUse(r"\omega", "Rotation angle of a frame-to-frame transform.")

_FRAME_CONCEPT = SeeAlso(
    "Reference frames and conventions", "../../concepts/reference_frames_and_conventions"
)
_FRAME_API = SeeAlso("ReferenceFrame / FrameTransform / FrameGraph", "../../api/index")


AXIS_CORRESPONDENCE_ROTATION = WorkedExample(
    id="frame-axis-correspondence-angle",
    title="Rotation angle implied by a declared axis correspondence",
    domain="core",
    scenario=(
        "EBSD vendors and analysis tools disagree about which specimen axis is called what. Rather "
        "than hand-writing a permutation matrix, you declare the correspondence in words: specimen "
        "x is the sample TD axis, specimen y is the reversed RD axis, specimen z is ND. That "
        "declaration is a 90-degree rotation about the shared third axis, and this example checks "
        "that PyTex builds exactly that rotation."
    ),
    setup=FRAME_SETUP,
    code=(
        "transform = FrameTransform.from_axis_correspondence(\n"
        "    specimen, sample, {\"x\": \"TD\", \"y\": \"-RD\", \"z\": \"ND\"}\n"
        ")\n"
        "result = transform.rotation_angle_deg"
    ),
    expected=90.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "The declaration fixes R e_x = e_TD, R e_y = -e_RD, R e_z = e_ND, i.e. the signed "
        "permutation [[0,-1,0],[1,0,0],[0,0,1]]. Its trace is 1, and the rotation angle follows "
        "from cos(omega) = (trace - 1) / 2 = 0, so omega = 90 degrees exactly."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, "
        "DOI: 10.1107/97809553602060000100 (right-handed axis conventions); "
        "trace identity for a proper rotation, Bunge, Texture Analysis in Materials Science, "
        "DOI: 10.1016/C2013-0-11769-2."
    ),
    symbols=(_ROTATION, _ANGLE),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


AXIS_CORRESPONDENCE_COMPONENTS = WorkedExample(
    id="frame-axis-correspondence-components",
    title="Components of the specimen x axis in a relabelled sample frame",
    domain="core",
    scenario=(
        "The point of a typed frame transform is that it converts components, not just angles. "
        "Having declared that specimen x is the sample TD axis, a direction lying along specimen "
        "x must come back with sample components (0, 1, 0): purely TD, no RD, no ND. This is the "
        "check that catches a reversed or transposed convention immediately."
    ),
    setup=FRAME_SETUP,
    code=(
        "transform = FrameTransform.from_axis_correspondence(\n"
        "    specimen, sample, {\"x\": \"TD\", \"y\": \"-RD\", \"z\": \"ND\"}\n"
        ")\n"
        "result = np.asarray(transform.apply_to_directions(np.array([1.0, 0.0, 0.0])))"
    ),
    expected=(0.0, 1.0, 0.0),
    unit="",
    tolerance=1e-12,
    reference=(
        "By definition of a basis, the source frame's x axis has source components e_x. The "
        "declaration 'x is TD' therefore forces its target components to be the TD basis vector "
        "(0, 1, 0)."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, "
        "DOI: 10.1107/97809553602060000100."
    ),
    symbols=(_ROTATION,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


FRAME_GRAPH_MULTI_HOP = WorkedExample(
    id="frame-graph-multi-hop-angle",
    title="Composing a two-hop frame chain with the frame graph",
    domain="core",
    scenario=(
        "A rolled sheet mounted 30 degrees off the stage axis gives two declared relationships: "
        "the canonical Cartesian reference to the specimen frame, and the specimen frame to the "
        "RD/TD/ND sample frame. You never declared the direct Cartesian-to-sample relationship, "
        "but you need it. The frame graph composes the shortest declared chain for you, and the "
        "result must be the 30-degree mounting rotation."
    ),
    setup=FRAME_SETUP,
    code=(
        "graph = rolling_frame_graph(rd_offset_deg=30.0)\n"
        "transform = graph.transform_between(\"cartesian\", \"sample_rd_td_nd\")\n"
        "result = transform.rotation_angle_deg"
    ),
    expected=30.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "The Cartesian-to-specimen edge is the identity and the specimen-to-sample edge is the "
        "declared 30-degree mounting rotation about ND, so the composition R = R2 R1 = R2 has "
        "rotation angle 30 degrees exactly."
    ),
    citation=(
        "Bunge, H.-J., Texture Analysis in Materials Science: Mathematical Methods, "
        "DOI: 10.1016/C2013-0-11769-2 (specimen-frame conventions for rolled sheet)."
    ),
    symbols=(_ANGLE,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


FRAME_ROUND_TRIP = WorkedExample(
    id="frame-round-trip-residual",
    title="Round-tripping components through a frame chain and back",
    domain="core",
    scenario=(
        "Any chain of frame transforms must be exactly invertible: converting a direction from the "
        "specimen frame into the sample frame and back has to return the original components. This "
        "is the invariant that guarantees no convention is silently lost when data crosses several "
        "module boundaries, so the residual is checked against exact zero."
    ),
    setup=FRAME_SETUP,
    code=(
        "graph = rolling_frame_graph(rd_offset_deg=37.5)\n"
        "direction = np.array([0.3, -0.7, 0.5])\n"
        "forward = graph.convert(\n"
        "    direction, source=\"specimen\", target=\"sample_rd_td_nd\", directions=True\n"
        ")\n"
        "back = graph.convert(\n"
        "    forward, source=\"sample_rd_td_nd\", target=\"specimen\", directions=True\n"
        ")\n"
        "result = float(np.max(np.abs(np.asarray(back) - direction)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-14,
    reference=(
        "A frame transform is a proper rotation, so R^-1 = R^T and R^T R = I exactly. The "
        "round-trip residual is therefore zero up to floating-point rounding."
    ),
    citation=(
        "Orthogonality of proper rotations; Bunge, Texture Analysis in Materials Science, "
        "DOI: 10.1016/C2013-0-11769-2."
    ),
    symbols=(_ROTATION,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


FRAME_DETERMINANT = WorkedExample(
    id="frame-right-handed-determinant",
    title="A right-handed frame has axis-vector determinant +1",
    domain="core",
    scenario=(
        "PyTex refuses to build a frame whose declared handedness contradicts its axis geometry, "
        "because a silently mirrored frame turns every downstream chirality result — variant "
        "selection, twin sense, pole-figure handedness — inside out. This example shows the "
        "invariant being reported for the standard RD/TD/ND sample frame."
    ),
    setup=FRAME_SETUP,
    code="result = sample.determinant",
    expected=1.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "For a right-handed orthonormal triad the determinant of the matrix whose columns are the "
        "axis vectors is exactly +1; a left-handed triad gives -1."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, "
        "DOI: 10.1107/97809553602060000100 (right-handed axial-frame convention)."
    ),
    symbols=(_ROTATION,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


_STAR = SymbolUse(r"\mathbf{a}^{*}", "Reciprocal basis vector, dual to the direct basis vector a.")
_GVEC = SymbolUse(r"\mathbf{g}_{hkl}", "Reciprocal-lattice vector of the (hkl) reflection.")

_NOTATION_SETUP = """
import numpy as np
from pytex import (
    crystal_frame,
    format_plane_family_indices,
    format_plane_indices,
    reciprocal_frame_for,
)

crystal = crystal_frame()
reciprocal = reciprocal_frame_for(crystal)
"""


RECIPROCAL_STAR_COUNT = WorkedExample(
    id="reciprocal-frame-star-count",
    title="Every reciprocal-frame axis carries the IUCr star",
    domain="core",
    scenario=(
        "In a workflow holding both direct and reciprocal quantities, the single most valuable "
        "safeguard is that a reciprocal-space vector cannot be mistaken for a direct-space one. "
        "PyTex enforces that by starring every axis of a reciprocal-domain frame. This example "
        "counts the starred axes on both frames: exactly three on the reciprocal frame, because "
        "the star belongs to the basis, and none on the direct crystal frame."
    ),
    setup=_NOTATION_SETUP,
    code=(
        "starred = sum(1 for axis in reciprocal.axes if axis.endswith(\"*\"))\n"
        "direct_starred = sum(1 for axis in crystal.axes if axis.endswith(\"*\"))\n"
        "result = np.array([starred, direct_starred])"
    ),
    expected=(3, 0),
    unit="",
    tolerance=0.0,
    reference=(
        "The reciprocal basis of a three-dimensional lattice has exactly three vectors "
        "a*, b*, c*, each conventionally starred; the direct basis a, b, c is never starred."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, "
        "DOI: 10.1107/97809553602060000100; reciprocal-space definitions, International "
        "Tables Volume C."
    ),
    symbols=(_STAR,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


MILLER_INDICES_UNSTARRED = WorkedExample(
    id="miller-indices-carry-no-star",
    title="Miller indices are never starred, in any bracket form",
    domain="core",
    scenario=(
        "The star marks the basis, not the indices. Miller indices are already reciprocal-basis "
        "components by definition, so starring them would name a different quantity - a mistake "
        "that is easy to make when 'reciprocal quantities are starred' is applied too broadly. "
        "This example counts the stars produced by every bracket form, which must be zero."
    ),
    setup=_NOTATION_SETUP,
    code=(
        "forms = [\n"
        "    format_plane_indices((1, 1, 1), style=\"plain\"),\n"
        "    format_plane_family_indices((1, 1, 1), style=\"plain\"),\n"
        "    format_plane_indices((1, -1, 0), style=\"plain\"),\n"
        "]\n"
        "result = float(sum(text.count(\"*\") for text in forms))"
    ),
    expected=0.0,
    unit="",
    tolerance=0.0,
    reference=(
        "By definition g_hkl = h a* + k b* + l c*: the indices (h, k, l) are the scalar "
        "coefficients of the starred basis vectors, so the indices carry no star themselves."
    ),
    citation=(
        "Hahn, Th. (ed.), International Tables for Crystallography, Volume A, "
        "DOI: 10.1107/97809553602060000100."
    ),
    symbols=(_GVEC,),
    see_also=(_FRAME_CONCEPT, _FRAME_API),
)


GROUP = ExampleGroup(
    slug="reference_frames",
    title="Reference Frames And Frame Transforms",
    summary=(
        "Creating standard frames, declaring frame relationships in words, and letting the frame "
        "graph compose multi-step chains — with the rotation angles, components, and invariants "
        "checked against exact analytic values. The last two examples pin the IUCr notation "
        "convention: the reciprocal star marks the basis, never the indices."
    ),
    examples=(
        AXIS_CORRESPONDENCE_ROTATION,
        AXIS_CORRESPONDENCE_COMPONENTS,
        FRAME_GRAPH_MULTI_HOP,
        FRAME_ROUND_TRIP,
        FRAME_DETERMINANT,
        RECIPROCAL_STAR_COUNT,
        MILLER_INDICES_UNSTARRED,
    ),
)
