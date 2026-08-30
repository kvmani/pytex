from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import math

import numpy as np
import pytest

from pytex import (
    Arrow3D,
    AxisTriad3D,
    Handedness,
    Lattice,
    Phase,
    PlanePatch3D,
    PointCloud3D,
    PolyLine3D,
    PrimitiveScene3D,
    ReferenceFrame,
    SymmetrySpec,
    Transform3D,
    crystal_plane_patch,
    direction_arrow,
    lattice_point_cloud,
    plane_normal_arrow,
    reference_frame_triad,
    render_primitive_scene_3d,
    unit_cell_polylines,
    vector_arrow,
)
from pytex.core.conventions import FrameDomain
from pytex.core.lattice import AtomicSite, CrystalDirection, CrystalPlane, MillerIndex, UnitCell
from pytex.core.orientation import Orientation, Rotation


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def _cubic_phase(a: float = 4.0) -> Phase:
    crystal = _crystal_frame()
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(AtomicSite(label="A1", species="Fe", fractional_coordinates=np.zeros(3)),),
    )
    return Phase(
        "cubic",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )


# --------------------------------------------------------------------------- #
# Transform3D
# --------------------------------------------------------------------------- #


def test_transform_identity_is_a_no_op() -> None:
    transform = Transform3D.identity()
    points = np.array([[1.0, 2.0, 3.0], [-1.0, 0.0, 4.0]])
    np.testing.assert_allclose(transform.apply_points(points), points)
    assert transform.is_rigid


def test_transform_rotation_plus_translation_applies_in_order() -> None:
    rotation = Rotation.from_axis_angle([0.0, 0.0, 1.0], np.deg2rad(90.0))
    transform = Transform3D.from_rotation(rotation, translation=[1.0, 2.0, 3.0])
    # +90 deg about z sends x-hat to y-hat, then translate
    np.testing.assert_allclose(
        transform.apply_points([1.0, 0.0, 0.0]), [1.0, 3.0, 3.0], atol=1e-9
    )
    # vectors ignore translation
    np.testing.assert_allclose(transform.apply_vector([1.0, 0.0, 0.0]), [0.0, 1.0, 0.0], atol=1e-9)


def test_transform_inverse_and_compose_round_trip() -> None:
    rotation = Rotation.from_bunge_euler(35.0, 20.0, 15.0)
    transform = Transform3D.from_rotation(rotation, translation=[2.0, -1.0, 0.5])
    inverse = transform.inverse()
    point = np.array([0.3, -0.7, 1.1])
    round_trip = inverse.apply_points(transform.apply_points(point))
    np.testing.assert_allclose(round_trip, point, atol=1e-9)
    identity = transform.compose(inverse)
    np.testing.assert_allclose(identity.matrix, np.eye(3), atol=1e-9)
    np.testing.assert_allclose(identity.translation, 0.0, atol=1e-9)


def test_transform_from_orientation_matches_map_crystal_vector() -> None:
    crystal = _crystal_frame()
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    orientation = Orientation.from_euler(
        30.0, 40.0, 10.0, convention="bunge", crystal_frame=crystal, specimen_frame=specimen
    )
    transform = Transform3D.from_orientation(orientation)
    vector = np.array([1.0, 2.0, -1.0])
    np.testing.assert_allclose(
        transform.apply_vector(vector),
        orientation.map_crystal_vector(vector),
        atol=1e-12,
    )


def test_transform_apply_normal_is_covariant_under_scaling() -> None:
    # anisotropic scaling: a normal must stay perpendicular to its plane
    transform = Transform3D.from_matrix(np.diag([2.0, 1.0, 0.5]))
    normal = np.array([0.0, 0.0, 1.0])
    in_plane = np.array([1.0, 0.0, 0.0])
    mapped_normal = transform.apply_normal(normal)
    mapped_in_plane = transform.apply_vector(in_plane)
    assert abs(float(mapped_normal @ mapped_in_plane)) < 1e-9


def test_transform_rejects_singular_matrix() -> None:
    with pytest.raises(ValueError, match="invertible"):
        Transform3D.from_matrix(np.zeros((3, 3)))


# --------------------------------------------------------------------------- #
# Primitive dataclasses
# --------------------------------------------------------------------------- #


def test_arrow_requires_distinct_endpoints() -> None:
    with pytest.raises(ValueError, match="distinct"):
        Arrow3D(tail=[0.0, 0.0, 0.0], head=[0.0, 0.0, 0.0])


def test_plane_patch_normal_is_normalized_and_vertices_are_coplanar() -> None:
    patch = PlanePatch3D(
        vertices=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        normal=[0.0, 0.0, 3.0],
    )
    np.testing.assert_allclose(np.linalg.norm(patch.normal), 1.0)
    offsets = patch.vertices @ patch.normal
    np.testing.assert_allclose(offsets, offsets[0], atol=1e-12)


def test_primitive_scene_transform_and_bounds() -> None:
    scene = PrimitiveScene3D(
        arrows=(Arrow3D(tail=[0.0, 0.0, 0.0], head=[2.0, 0.0, 0.0]),),
        point_clouds=(PointCloud3D(points=[[0.0, 0.0, 0.0], [0.0, 3.0, 0.0]]),),
    )
    bounds = scene.bounds()
    np.testing.assert_allclose(bounds[0], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(bounds[1], [2.0, 3.0, 0.0])
    shifted = scene.transformed(Transform3D.from_matrix(np.eye(3), translation=[1.0, 1.0, 1.0]))
    np.testing.assert_allclose(shifted.bounds()[0], [1.0, 1.0, 1.0])
    assert not scene.is_empty()
    assert PrimitiveScene3D().is_empty()


def test_axis_triad_expands_to_three_labelled_arrows() -> None:
    triad = AxisTriad3D(labels=("a", "b", "c"))
    arrows = triad.arrows()
    assert len(arrows) == 3
    np.testing.assert_allclose(arrows[0].head, [1.0, 0.0, 0.0])
    labels = triad.tip_labels()
    assert tuple(label.text for label in labels) == ("a", "b", "c")


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #


def test_direction_arrow_points_along_unit_vector() -> None:
    phase = _cubic_phase()
    direction = CrystalDirection(np.array([1.0, 1.0, 0.0]), phase=phase)
    arrow = direction_arrow(direction, length=2.0)
    np.testing.assert_allclose(arrow.vector, 2.0 * direction.unit_vector, atol=1e-12)
    assert arrow.label == "$[110]$"


def test_plane_normal_arrow_and_patch_align_with_plane_normal() -> None:
    phase = _cubic_phase()
    plane = CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=phase), phase=phase)
    arrow = plane_normal_arrow(plane, length=1.0)
    np.testing.assert_allclose(arrow.vector, plane.normal, atol=1e-12)
    patch = crystal_plane_patch(plane, extent=1.5)
    np.testing.assert_allclose(patch.normal, plane.normal, atol=1e-12)
    assert patch.label == "$(111)$"


def test_unit_cell_polylines_have_twelve_edges() -> None:
    phase = _cubic_phase(a=3.0)
    edges = unit_cell_polylines(phase)
    assert len(edges) == 12
    assert all(isinstance(edge, PolyLine3D) for edge in edges)
    # every edge of a cubic cell has length a
    for edge in edges:
        np.testing.assert_allclose(np.linalg.norm(edge.points[1] - edge.points[0]), 3.0, atol=1e-9)


def test_lattice_point_cloud_counts_supercell_nodes() -> None:
    phase = _cubic_phase()
    cloud = lattice_point_cloud(phase, repeats=(2, 1, 1))
    # (2+1) * (1+1) * (1+1) = 12 Bravais nodes
    assert cloud.points.shape == (12, 3)


def test_reference_frame_triad_labels_from_frame_axes() -> None:
    triad = reference_frame_triad(_crystal_frame(), length=2.0)
    assert triad.labels == ("a", "b", "c")
    np.testing.assert_allclose(np.linalg.norm(triad.axes[:, 0]), 2.0)


def test_vector_arrow_scales_from_origin() -> None:
    arrow = vector_arrow([0.0, 0.0, 2.0], origin=[1.0, 0.0, 0.0], scale=0.5)
    np.testing.assert_allclose(arrow.tail, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(arrow.head, [1.0, 0.0, 1.0])


# --------------------------------------------------------------------------- #
# Renderer
# --------------------------------------------------------------------------- #


def test_render_primitive_scene_smoke() -> None:
    phase = _cubic_phase()
    scene = PrimitiveScene3D(
        arrows=(vector_arrow([1.0, 0.0, 0.0], color="#dc2626", label="v"),),
        polylines=unit_cell_polylines(phase),
        triads=(reference_frame_triad(_crystal_frame()),),
        point_clouds=(lattice_point_cloud(phase),),
    )
    figure = render_primitive_scene_3d(scene, title="primitives")
    assert len(figure.axes) == 1
    axis = figure.axes[0]
    assert axis.get_title() == "primitives"


class TestPlaneClipping:
    """A plane overlay is the lattice plane cut by the cell, not a free square.

    The claim under test is geometric and checkable without a reference image:
    every vertex satisfies the plane equation, every vertex lies in the box, and
    the polygon closes. The old square failed the second of those by
    construction — it was centred on the origin and sized by a scene-scale
    guess, so half of it lay outside the crystal it was drawn for.
    """

    @staticmethod
    def _phase(builtin: str):  # type: ignore[no-untyped-def]
        from pytex.app.phases import phase_from_request

        return phase_from_request({"builtin": builtin})[1]

    def test_the_cubic_110_is_the_diagonal_rectangle_of_the_cell(self) -> None:
        """(110) of a cubic cell enters at one edge and leaves at the opposite one.

        Its area is the exact ``a * a * sqrt(2)`` of the diagonal section, which
        is the whole of the claim: not a square of some convenient size with the
        right normal, but the plane the cell actually cuts.
        """

        from pytex.plotting.primitives import lattice_plane_polygon

        phase = self._phase("zr_bcc_beta")
        edge = float(phase.lattice.a)
        polygon = lattice_plane_polygon(phase, (1, 1, 0))
        assert polygon is not None
        assert polygon.shape == (4, 3)
        # Every vertex inside the cell box, to the floating-point floor.
        assert np.all(polygon >= -1e-9)
        assert np.all(polygon <= edge + 1e-9)
        # Every vertex on one plane: x + y is constant, and equal to the cell edge.
        assert np.allclose(polygon[:, 0] + polygon[:, 1], edge)
        area = float(np.linalg.norm(np.cross(polygon[1] - polygon[0], polygon[2] - polygon[1])))
        assert area == pytest.approx(edge * edge * math.sqrt(2.0), rel=1e-9)

    def test_the_offset_chosen_is_the_largest_cross_section(self) -> None:
        """The member a reader means by "the (110) plane of this cell".

        The family has members through the origin corner and through the far
        corner, and both are degenerate — they touch an edge and cut no area.
        Choosing by area rejects them without a special case.
        """

        from pytex.plotting.primitives import lattice_plane_polygon

        phase = self._phase("zr_bcc_beta")
        chosen = lattice_plane_polygon(phase, (1, 1, 0))
        assert chosen is not None
        best = 0.0
        for offset in range(0, 3):
            candidate = lattice_plane_polygon(phase, (1, 1, 0), offset=float(offset))
            if candidate is None:
                continue
            edges = np.diff(np.vstack([candidate, candidate[:1]]), axis=0)
            best = max(best, float(np.linalg.norm(np.cross(edges[0], edges[1]))))
        edges = np.diff(np.vstack([chosen, chosen[:1]]), axis=0)
        assert float(np.linalg.norm(np.cross(edges[0], edges[1]))) == pytest.approx(best)

    def test_a_hexagonal_basal_plane_is_the_cell_face(self) -> None:
        from pytex.plotting.primitives import lattice_plane_polygon

        phase = self._phase("zr_hcp")
        polygon = lattice_plane_polygon(phase, (0, 0, 1))
        assert polygon is not None
        # Flat in z, and at one of the two basal faces of the cell.
        assert np.allclose(polygon[:, 2], polygon[0, 2])
        assert float(polygon[0, 2]) == pytest.approx(float(phase.lattice.c))

    def test_the_patch_builder_clips_when_asked_and_squares_when_not(self) -> None:
        from pytex.core.lattice import CrystalPlane, MillerIndex
        from pytex.plotting.primitives import crystal_plane_patch

        phase = self._phase("zr_bcc_beta")
        plane = CrystalPlane(MillerIndex(np.array([1, 1, 0]), phase=phase), phase=phase)
        clipped = crystal_plane_patch(plane, cell_repeats=(1, 1, 1))
        square = crystal_plane_patch(plane, extent=2.0)
        assert np.all(clipped.vertices >= -1e-9)
        # The free square straddles the origin, which is exactly why it is not
        # what a plane overlay should be.
        assert np.any(square.vertices < -1e-9)

    def test_a_direction_in_the_plane_is_drawn_as_a_chord_of_it(self) -> None:
        """The arrow lies in the patch, which is the claim a pair makes.

        Both endpoints are on the polygon's boundary and every point between
        them is inside it, so the arrow cannot leave the plane or the cell.
        """

        from pytex.plotting.primitives import lattice_plane_polygon, segment_in_polygon

        phase = self._phase("zr_bcc_beta")
        polygon = lattice_plane_polygon(phase, (1, 1, 0))
        basis = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
        direction = basis @ np.array([-1.0, 1.0, 1.0])
        chord = segment_in_polygon(polygon, direction)
        assert chord is not None
        tail, head = chord
        plane_offset = float(polygon[0, 0] + polygon[0, 1])
        for point in (tail, head, 0.5 * (tail + head)):
            assert float(point[0] + point[1]) == pytest.approx(plane_offset)
            assert np.all(point >= -1e-9)
            assert np.all(point <= float(phase.lattice.a) + 1e-9)
        # And it really is along the direction asked for.
        along = (head - tail) / np.linalg.norm(head - tail)
        unit = direction / np.linalg.norm(direction)
        assert abs(float(np.dot(along, unit))) == pytest.approx(1.0)

    def test_a_direction_out_of_the_plane_is_refused_rather_than_projected(self) -> None:
        """Refused, because a projected arrow would assert a parallelism.

        The caller falls back to clipping against the cell, which keeps the
        arrow inside the crystal without claiming it lies in a plane it does
        not lie in.
        """

        from pytex.plotting.primitives import lattice_plane_polygon, segment_in_polygon

        phase = self._phase("zr_bcc_beta")
        polygon = lattice_plane_polygon(phase, (1, 1, 0))
        assert segment_in_polygon(polygon, np.array([1.0, 1.0, 0.0])) is None

    def test_the_cell_fallback_keeps_a_direction_inside_the_box(self) -> None:
        from pytex.plotting.primitives import segment_in_cell

        phase = self._phase("zr_bcc_beta")
        edge = float(phase.lattice.a)
        tail, head = segment_in_cell(phase, np.array([1.0, 0.0, 0.0]))
        for point in (tail, head):
            assert np.all(point >= -1e-9)
            assert np.all(point <= edge + 1e-9)
        assert float(np.linalg.norm(head - tail)) == pytest.approx(edge)
