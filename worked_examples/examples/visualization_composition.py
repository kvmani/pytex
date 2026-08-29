"""Worked examples: composable 3D visualization primitives.

These examples validate the geometric core of the visualization-primitives
layer -- the part that must be numerically correct for every composite figure
built on top of it. The first checks that `Transform3D.from_orientation`
reproduces the canonical crystal-to-sample vector map, so a crystal drawn in the
sample frame lands exactly where the orientation says it should. The second
checks the placement used by `WorldScene3D.from_orientation_relationship`: after
the child crystal is placed by the inverse orientation-relationship rotation,
the relationship's parallel directions become exactly parallel in the shared
world frame -- the geometric statement of the relationship, shown directly.

Both reference values are exact identities: a definitional round trip (deviation
0) and the parallelism that *defines* a parallel-direction orientation
relationship (direction cosine 1).

See :doc:`../../concepts/visualization_primitives`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

TRANSFORM_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    Transform3D,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
grain = Orientation.from_euler(
    30.0, 40.0, 10.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
"""

RELATIONSHIP_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    Transform3D,
)
from pytex.core.transformation import OrientationRelationship

def cubic_phase(name, a):
    frame = ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    return Phase(
        name,
        lattice=Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
    )

fcc = cubic_phase("austenite", 3.60)
bcc = cubic_phase("ferrite", 2.87)
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=fcc, child_phase=bcc
)
"""

SCENE_MEASUREMENT_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    build_crystal_scene,
)
from pytex.core.lattice import AtomicSite, UnitCell

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
halite_like = Phase(
    "halite-like-pair",
    lattice=lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Na1", species="Na", fractional_coordinates=np.zeros(3)),
            AtomicSite(
                label="Cl1",
                species="Cl",
                fractional_coordinates=np.array([0.5, 0.0, 0.0]),
            ),
        ),
    ),
)
"""

_ROT = SymbolUse(r"\mathbf{R}", "Active rotation matrix.")
_TRANSFORM = SymbolUse(r"\mathbf{T}", "Rigid placement (rotation and translation) into the world frame.")
_LATTICE_A = SymbolUse(r"a", "Cubic lattice parameter (edge length).")

_VIZ_CONCEPT = SeeAlso("Visualization primitives", "../../concepts/visualization_primitives")
_API = SeeAlso("Plotting and visualization API", "../../api/index")


TRANSFORM_CRYSTAL_TO_SAMPLE = WorkedExample(
    id="viz-transform-crystal-to-sample-consistency",
    title="Transform3D.from_orientation reproduces the crystal-to-sample map",
    domain="visualization",
    scenario=(
        "To draw a grain's structure in the sample frame you place its crystal geometry with a "
        "Transform3D built from the grain orientation. That placement must agree exactly with the "
        "orientation's own crystal-to-sample vector map, or every arrow, plane, and atom would be "
        "drawn in the wrong direction. Here we map a crystal-frame vector both ways and measure the "
        "difference; a correct placement returns exactly zero."
    ),
    setup=TRANSFORM_SETUP,
    code=(
        "vector = np.array([1.0, 2.0, -1.0])\n"
        "placement = Transform3D.from_orientation(grain)\n"
        "placed = placement.apply_vector(vector)\n"
        "expected = grain.map_crystal_vector(vector)\n"
        "result = float(np.max(np.abs(placed - expected)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "Transform3D.from_orientation(o) is defined to apply the orientation matrix g, which is the "
        "crystal-to-sample map o.map_crystal_vector; the two are identical, so the deviation is 0."
    ),
    citation="Bunge, Texture Analysis in Materials Science, 1982, Chapter 2 (orientation matrix).",
    symbols=(_ROT, _TRANSFORM),
    see_also=(_VIZ_CONCEPT, _API),
)

OR_PARALLEL_DIRECTION_ALIGNMENT = WorkedExample(
    id="viz-or-parallel-direction-alignment",
    title="Placing the child crystal aligns the KS parallel directions",
    domain="visualization",
    scenario=(
        "A two-crystal orientation-relationship figure places the child crystal so the relationship "
        "holds. WorldScene3D.from_orientation_relationship uses the inverse parent-to-child rotation "
        "for that placement. We reproduce that placement and check the defining property of the "
        "Kurdjumov-Sachs relationship: the paired close-packed directions "
        "(<110> in fcc, <111> in bcc) become exactly parallel in the shared world frame. The "
        "direction cosine is therefore 1."
    ),
    setup=RELATIONSHIP_SETUP,
    code=(
        "child_placement = Transform3D.from_matrix(\n"
        "    ks.parent_to_child_rotation.inverse().as_matrix()\n"
        ")\n"
        "parent_direction, child_direction = ks.parallel_directions[0]\n"
        "placed = child_placement.apply_vector(child_direction.unit_vector)\n"
        "placed = placed / np.linalg.norm(placed)\n"
        "result = float(parent_direction.unit_vector @ placed)"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-9,
    reference=(
        "The Kurdjumov-Sachs relationship fixes <110>_fcc || <111>_bcc; placing the child by the "
        "inverse relationship rotation makes the paired directions collinear, so their cosine is 1."
    ),
    citation="Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.",
    symbols=(_ROT, _TRANSFORM),
    see_also=(
        _VIZ_CONCEPT,
        SeeAlso("Orientation relationships", "../../concepts/orientation_texture"),
    ),
)


SCENE_BOND_LENGTH_MEASUREMENT = WorkedExample(
    id="viz-scene-bond-length-halite-identity",
    title="Measured scene bond length equals the exact NaCl-type a/2 distance",
    domain="visualization",
    scenario=(
        "Crystal scenes are only trustworthy teaching and publication figures if the geometry they "
        "draw is exact. The programmatic distance readout (CrystalScene.bond_lengths_angstrom, the "
        "scriptable analog of VESTA's click measurement) makes that checkable: in an NaCl-type "
        "arrangement the nearest-neighbour cation-anion distance is exactly half the cubic lattice "
        "parameter, so a scene built at a = 4 angstrom must measure every bond at 2 angstrom."
    ),
    setup=SCENE_MEASUREMENT_SETUP,
    code=(
        "scene = build_crystal_scene(halite_like, include_boundary_atoms=False)\n"
        "result = float(scene.bond_lengths_angstrom()[0])"
    ),
    expected=2.0,
    unit="angstrom",
    tolerance=1e-12,
    reference=(
        "NaCl-type geometry: the cation at (0, 0, 0) and anion at (1/2, 0, 0) are separated by a/2 "
        "along the cube edge; with a = 4 angstrom the bond length is exactly 2 angstrom."
    ),
    citation=(
        "International Tables for Crystallography, Vol. A (rock-salt structure geometry); "
        "Momma and Izumi, J. Appl. Cryst. 44 (2011) 1272 (VESTA distance readout)."
    ),
    symbols=(_LATTICE_A,),
    see_also=(_VIZ_CONCEPT, _API),
)

_OR_STEREOGRAM_CODE = """
from pytex.plotting.spherical import build_or_stereogram_figure_spec

worst_pole_gap = 0.0
worst_circle_gap = 0.0
for index in range(1, 25):
    spec = build_or_stereogram_figure_spec(
        ks, variant=index, include_wulff_net=False, show_tie_lines=False
    )
    parent_markers, child_markers = spec.marker_layers
    worst_pole_gap = max(
        worst_pole_gap,
        float(np.max(np.linalg.norm(parent_markers.points - child_markers.points, axis=1))),
    )
    parent_circle, child_circle = spec.line_layers[0].points, spec.line_layers[1].points
    worst_circle_gap = max(
        worst_circle_gap,
        float(np.max(np.linalg.norm(parent_circle - child_circle, axis=1))),
    )
result = [worst_pole_gap, worst_circle_gap]
""".strip()


OR_STEREOGRAM_COINCIDENCE = WorkedExample(
    id="viz-or-stereogram-parallelism-coincides",
    title="The OR stereogram plots a parallelism as one point and one circle",
    domain="visualization",
    scenario=(
        "The orientation-relationship stereogram makes two visual claims: the "
        "parent pole and the child pole of a parallel pair land on the same point "
        "of the net, and the great circles of two parallel planes lie on top of "
        "each other. Both are checked here for all 24 Kurdjumov-Sachs variants, "
        "not only the one the relationship was written with, by measuring the "
        "worst separation in projection-plane units. A separation that is not "
        "zero would mean the figure draws two objects where the crystallography "
        "has one -- which is exactly what happened before the pair was folded onto "
        "the upper hemisphere as a pair: variants 7 and 9 have a defining "
        "direction in the equatorial plane, and their two ends landed on opposite "
        "rims, a full disc diameter apart."
    ),
    setup=RELATIONSHIP_SETUP,
    code=_OR_STEREOGRAM_CODE,
    expected=[0.0, 0.0],
    unit="projection-plane units",
    tolerance=1e-9,
    reference=(
        "An identity, not a measurement. Each transformation variant maps its own "
        "parent normal exactly onto its own child normal, so once the child pole "
        "is carried back into the parent frame the two are the same unit vector "
        "and every projection of them coincides. The 1e-9 tolerance is a "
        "floating-point floor, not a physical margin; the measured worst gap is "
        "of order 1e-15."
    ),
    citation="Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.",
    symbols=(_TRANSFORM,),
    see_also=(_VIZ_CONCEPT, _API),
)


GROUP = ExampleGroup(
    slug="visualization",
    title="Composable visualization primitives",
    summary=(
        "Geometric guarantees of the visualization layer: a placement transform that reproduces "
        "the crystal-to-sample map, the orientation-relationship placement that makes parallel "
        "directions coincide in one world frame, a scene bond-length measurement checked "
        "against the exact NaCl-type a/2 distance, and the OR stereogram plotting a "
        "parallelism as one point and one circle for every variant."
    ),
    examples=(
        TRANSFORM_CRYSTAL_TO_SAMPLE,
        OR_PARALLEL_DIRECTION_ALIGNMENT,
        SCENE_BOND_LENGTH_MEASUREMENT,
        OR_STEREOGRAM_COINCIDENCE,
    ),
)

__all__ = ["GROUP"]
