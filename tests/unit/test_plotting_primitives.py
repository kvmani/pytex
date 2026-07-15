from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

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
