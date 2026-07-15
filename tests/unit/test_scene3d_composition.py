from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pytex import (
    Handedness,
    Lattice,
    Phase,
    PlacedCrystal,
    PrimitiveScene3D,
    ReferenceFrame,
    SymmetrySpec,
    Transform3D,
    WorldScene3D,
    build_crystal_scene,
    plot_stereographic_vectors,
    render_world_scene_3d,
    vector_arrow,
)
from pytex.core.conventions import FrameDomain
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.core.orientation import Rotation
from pytex.core.transformation import OrientationRelationship
from pytex.plotting.spherical import build_vector_stereogram_figure_spec


def _cubic_phase(name: str, a: float = 4.0, species: str = "Fe") -> Phase:
    crystal = ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    site = AtomicSite(label=f"{species}1", species=species, fractional_coordinates=np.zeros(3))
    unit_cell = UnitCell(lattice=lattice, sites=(site,))
    return Phase(
        name,
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )


def _ks_relationship() -> OrientationRelationship:
    parent = _cubic_phase("fcc", a=3.6, species="Ni")
    child = _cubic_phase("bcc", a=2.9, species="Fe")
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )


# --------------------------------------------------------------------------- #
# PlacedCrystal / WorldScene3D basics
# --------------------------------------------------------------------------- #


def test_placed_crystal_identity_returns_same_scene() -> None:
    scene = build_crystal_scene(_cubic_phase("fcc"), include_boundary_atoms=False)
    placed = PlacedCrystal(scene=scene)
    assert placed.placed_scene() is scene


def test_placed_crystal_translation_shifts_atoms() -> None:
    scene = build_crystal_scene(_cubic_phase("fcc"), include_boundary_atoms=False)
    transform = Transform3D.from_matrix(np.eye(3), translation=[10.0, 0.0, 0.0])
    placed = PlacedCrystal(scene=scene, transform=transform)
    moved = placed.placed_scene()
    np.testing.assert_allclose(
        moved.atoms[0].position_angstrom - scene.atoms[0].position_angstrom,
        [10.0, 0.0, 0.0],
        atol=1e-9,
    )


def test_world_scene_add_crystal_and_bounds() -> None:
    phase = _cubic_phase("fcc")
    world = WorldScene3D().add_crystal(phase, include_boundary_atoms=False)
    assert len(world.crystals) == 1
    shifted = WorldScene3D().add_crystal(
        phase,
        transform=Transform3D.from_matrix(np.eye(3), translation=[20.0, 0.0, 0.0]),
        include_boundary_atoms=False,
    )
    assert shifted.bounds()[1][0] > world.bounds()[1][0] + 15.0


def test_world_scene_is_immutable_under_add() -> None:
    phase = _cubic_phase("fcc")
    base = WorldScene3D()
    grown = base.add_crystal(phase, include_boundary_atoms=False)
    assert len(base.crystals) == 0
    assert len(grown.crystals) == 1


def test_crystal_scene_transformed_rejects_non_rigid() -> None:
    scene = build_crystal_scene(_cubic_phase("fcc"), include_boundary_atoms=False)
    with pytest.raises(ValueError, match="rigid"):
        scene.transformed(Transform3D.from_matrix(np.diag([2.0, 1.0, 1.0])))


# --------------------------------------------------------------------------- #
# Composite rendering: one globally depth-sorted mesh across crystals
# --------------------------------------------------------------------------- #


def test_two_plain_crystals_share_one_depth_sorted_mesh() -> None:
    phase = _cubic_phase("fcc")
    resolution = 10
    world = (
        WorldScene3D()
        .add_crystal(
            phase,
            include_boundary_atoms=False,
            show_bonds=False,
            style_overrides={"crystal": {"atom_surface_resolution": resolution}},
        )
        .add_crystal(
            phase,
            transform=Transform3D.from_matrix(np.eye(3), translation=[8.0, 0.0, 0.0]),
            include_boundary_atoms=False,
            show_bonds=False,
            style_overrides={"crystal": {"atom_surface_resolution": resolution}},
        )
    )
    figure = render_world_scene_3d(
        world, style_overrides={"crystal": {"atom_surface_resolution": resolution}}
    )
    figure.canvas.draw()
    meshes = [a for a in figure.axes[0].collections if isinstance(a, Poly3DCollection)]
    # both crystals' atom faces live in ONE collection -> correct global occlusion
    assert len(meshes) == 1
    assert len(meshes[0].get_paths()) == 2 * (resolution - 1) ** 2


def test_world_scene_renders_loose_primitives() -> None:
    phase = _cubic_phase("fcc")
    world = WorldScene3D().add_crystal(phase, include_boundary_atoms=False).add_primitives(
        PrimitiveScene3D(arrows=(vector_arrow([1.0, 0.0, 0.0], scale=5.0, label="load"),))
    )
    figure = render_world_scene_3d(world)
    assert len(figure.axes) == 1


# --------------------------------------------------------------------------- #
# Orientation-relationship composition (the north-star use case)
# --------------------------------------------------------------------------- #


def test_from_orientation_relationship_builds_two_crystals_and_primitives() -> None:
    relationship = _ks_relationship()
    world = WorldScene3D.from_orientation_relationship(relationship, repeats=(1, 1, 1))
    assert len(world.crystals) == 2
    assert world.crystals[0].label == "fcc"
    assert world.crystals[1].label == "bcc"
    assert len(world.primitives.arrows) == len(relationship.parallel_directions)
    assert len(world.primitives.patches) == len(relationship.parallel_planes)


def test_or_placement_makes_parallel_directions_and_planes_coincide() -> None:
    relationship = _ks_relationship()
    world = WorldScene3D.from_orientation_relationship(relationship)
    child_transform = world.crystals[1].transform
    for parent_direction, child_direction in relationship.parallel_directions:
        placed = child_transform.apply_vector(np.asarray(child_direction, dtype=np.float64))
        placed = placed / np.linalg.norm(placed)
        parent_unit = parent_direction / np.linalg.norm(parent_direction)
        assert float(parent_unit @ placed) == pytest.approx(1.0, abs=1e-6)
    for parent_plane, child_plane in relationship.parallel_planes:
        placed = child_transform.apply_normal(child_plane.normal)
        assert float(parent_plane.normal @ placed) == pytest.approx(1.0, abs=1e-6)


def test_or_child_translation_separates_the_crystals() -> None:
    relationship = _ks_relationship()
    overlapping = WorldScene3D.from_orientation_relationship(relationship)
    separated = WorldScene3D.from_orientation_relationship(
        relationship, child_translation=[30.0, 0.0, 0.0]
    )
    assert separated.bounds()[1][0] > overlapping.bounds()[1][0] + 20.0


def test_or_render_smoke_with_legend() -> None:
    relationship = _ks_relationship()
    world = WorldScene3D.from_orientation_relationship(relationship)
    figure = world.render(show_legend=True, title="KS OR")
    assert figure.axes[0].get_legend() is not None


# --------------------------------------------------------------------------- #
# Stereographic vector bridge (projection analog)
# --------------------------------------------------------------------------- #


def test_vector_stereogram_spec_poles_and_traces() -> None:
    vectors = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 1.0]])
    spec = build_vector_stereogram_figure_spec(vectors, labels=["a", "b", "c"], render="both")
    assert len(spec.marker_layers) == 1
    assert spec.marker_layers[0].points.shape[0] == 3
    # three great-circle traces on top of the Wulff-net lines
    net_only = build_vector_stereogram_figure_spec(vectors, render="pole")
    assert len(spec.line_layers) == len(net_only.line_layers) + 3


def test_vector_stereogram_rejects_zero_vector_and_bad_render() -> None:
    with pytest.raises(ValueError, match="zero vector"):
        build_vector_stereogram_figure_spec(np.array([[0.0, 0.0, 0.0]]))
    with pytest.raises(ValueError, match="render"):
        build_vector_stereogram_figure_spec(np.array([[1.0, 0.0, 0.0]]), render="bogus")


def test_plot_stereographic_vectors_overlays_two_crystals_in_common_frame() -> None:
    # parent [100] and a child direction rotated into the parent frame land on
    # one stereogram, coloured distinctly -- the projection analog of the 3D OR
    relationship = _ks_relationship()
    child_transform = Transform3D.from_matrix(
        relationship.parent_to_child_rotation.inverse().as_matrix()
    )
    parent_dir = np.array([1.0, 0.0, 0.0])
    child_dir_in_parent = child_transform.apply_vector(np.array([1.0, 0.0, 0.0]))
    figure = plot_stereographic_vectors(
        np.vstack([parent_dir, child_dir_in_parent]),
        labels=["fcc[100]", "bcc[100]"],
        colors=["#2563eb", "#dc2626"],
    )
    assert len(figure.axes) == 1


def test_transform_from_rotation_identity_matches_expected() -> None:
    transform = Transform3D.from_rotation(Rotation.identity())
    np.testing.assert_allclose(transform.matrix, np.eye(3))
