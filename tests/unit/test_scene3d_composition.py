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
    render_variant_contact_sheet,
    render_world_scene_3d,
    vector_arrow,
)
from pytex.core.conventions import FrameDomain
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.core.notation import format_miller_indices
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
    """Two crystals, and every parallel object drawn on **each** of them.

    Two per pair rather than one: a parallelism is a statement about both
    crystals, and each side is clipped to its own cell, so the parent's plane
    and the child's plane are the same physical plane with different outlines.
    One patch at the origin -- which is what this asserted before -- put the
    whole statement on one crystal and left the other unmarked.
    """

    relationship = _ks_relationship()
    world = WorldScene3D.from_orientation_relationship(relationship, repeats=(1, 1, 1))
    assert len(world.crystals) == 2
    assert world.crystals[0].label == "fcc"
    assert world.crystals[1].label == "bcc"
    assert len(world.primitives.arrows) == 2 * len(relationship.parallel_directions)
    assert len(world.primitives.patches) == 2 * len(relationship.parallel_planes)
    # The pair is named once, on the parent's copy; the child's copy is silent,
    # so one statement does not print itself twice over one figure.
    assert sum(1 for patch in world.primitives.patches if patch.label) == len(
        relationship.parallel_planes
    )
    assert sum(1 for arrow in world.primitives.arrows if arrow.label) == len(
        relationship.parallel_directions
    )


def test_or_placement_makes_parallel_directions_and_planes_coincide() -> None:
    relationship = _ks_relationship()
    world = WorldScene3D.from_orientation_relationship(relationship)
    child_transform = world.crystals[1].transform
    for parent_direction, child_direction in relationship.parallel_directions:
        placed = child_transform.apply_vector(child_direction.unit_vector)
        placed = placed / np.linalg.norm(placed)
        assert float(parent_direction.unit_vector @ placed) == pytest.approx(1.0, abs=1e-6)
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


# --------------------------------------------------------------------------- #
# Variant-aware composite scenes (F15)
# --------------------------------------------------------------------------- #


def _child_matrix(world: WorldScene3D) -> np.ndarray:
    return np.asarray(world.crystals[1].transform.matrix, dtype=np.float64)


def test_variant_scene_places_child_by_that_variants_rotation() -> None:
    relationship = _ks_relationship()
    variants = relationship.generate_variants()
    for index, variant in ((1, variants[0]), (17, variants[16])):
        world = WorldScene3D.from_orientation_relationship(relationship, variant=index)
        np.testing.assert_allclose(
            _child_matrix(world),
            variant.parent_to_child_rotation.inverse().as_matrix(),
            atol=1e-12,
        )


def test_every_variant_scene_makes_its_own_pair_coincide() -> None:
    # the M2 validation: for EVERY variant the world-frame images of the parent
    # and child defining normals (and directions) must coincide, not just V1
    relationship = _ks_relationship()
    for variant in relationship.generate_variants():
        world = WorldScene3D.from_orientation_relationship(relationship, variant=variant)
        child_transform = world.crystals[1].transform
        for parent_plane, child_plane in variant.parallel_planes:
            placed = child_transform.apply_normal(child_plane.normal)
            assert float(parent_plane.normal @ placed) == pytest.approx(1.0, abs=1e-12)
        for parent_direction, child_direction in variant.parallel_directions:
            placed = child_transform.apply_vector(child_direction.unit_vector)
            placed = placed / np.linalg.norm(placed)
            assert float(parent_direction.unit_vector @ placed) == pytest.approx(1.0, abs=1e-12)


def test_nominal_pair_would_be_wrong_on_a_later_variant() -> None:
    # the trap this feature exists to close: drawing the relationship's nominal
    # plane on variant 17 gives a figure that looks right and is wrong
    relationship = _ks_relationship()
    variant = relationship.generate_variants()[16]
    nominal_parent_plane = relationship.parallel_planes[0][0]
    child_transform = Transform3D.from_matrix(
        variant.parent_to_child_rotation.inverse().as_matrix()
    )
    nominal_child_normal = child_transform.apply_normal(relationship.parallel_planes[0][1].normal)
    assert abs(float(nominal_parent_plane.normal @ nominal_child_normal)) < 0.99
    # while the variant's own pair coincides exactly
    variant_parent_plane, variant_child_plane = variant.parallel_planes[0]
    assert float(
        variant_parent_plane.normal @ child_transform.apply_normal(variant_child_plane.normal)
    ) == pytest.approx(1.0, abs=1e-12)


def test_variant_scene_primitive_labels_name_the_variants_indices() -> None:
    relationship = _ks_relationship()
    variant = relationship.generate_variants()[16]
    world = WorldScene3D.from_orientation_relationship(relationship, variant=variant)
    parent_plane, child_plane = variant.parallel_planes[0]
    expected = "{} ∥ {}".format(
        format_miller_indices(parent_plane.miller.indices, family="plane", style="plain"),
        format_miller_indices(child_plane.miller.indices, family="plane", style="plain"),
    )
    assert world.primitives.patches[0].label == expected


def test_variant_scenes_reproduce_generate_variants_as_a_set() -> None:
    relationship = _ks_relationship()
    variants = relationship.generate_variants()
    scenes = WorldScene3D.variant_scenes(relationship)
    assert len(scenes) == len(variants)
    child_operators = np.asarray(relationship.child_phase.symmetry.operators, dtype=np.float64)
    matched: list[int] = []
    for scene in scenes:
        placement = _child_matrix(scene).T  # the variant rotation itself
        residuals = [
            float(
                np.min(
                    np.linalg.norm(
                        child_operators @ variant.parent_to_child_rotation.as_matrix()
                        - placement,
                        axis=(1, 2),
                    )
                )
            )
            for variant in variants
        ]
        best = int(np.argmin(residuals))
        assert residuals[best] < 1e-10
        matched.append(best)
    assert sorted(matched) == list(range(len(variants)))


def test_variant_scenes_accepts_a_subset_and_rejects_a_duplicate_variant_kwarg() -> None:
    relationship = _ks_relationship()
    variants = relationship.generate_variants()[:3]
    scenes = WorldScene3D.variant_scenes(relationship, variants=variants)
    assert len(scenes) == 3
    with pytest.raises(ValueError, match="supplies 'variant' itself"):
        WorldScene3D.variant_scenes(relationship, variant=1)


def test_from_orientation_relationship_rejects_out_of_range_variant_index() -> None:
    relationship = _ks_relationship()
    with pytest.raises(ValueError, match=r"one-based index in 1\.\.24"):
        WorldScene3D.from_orientation_relationship(relationship, variant=25)
    with pytest.raises(ValueError, match=r"one-based index in 1\.\.24"):
        WorldScene3D.from_orientation_relationship(relationship, variant=0)


def test_variant_contact_sheet_lays_out_one_panel_per_scene() -> None:
    relationship = _ks_relationship()
    variants = relationship.generate_variants()[:6]
    scenes = WorldScene3D.variant_scenes(relationship, variants=variants, repeats=(1, 1, 1))
    figure = render_variant_contact_sheet(
        scenes, variants=variants, columns=3, suptitle="Kurdjumov-Sachs"
    )
    assert len(figure.axes) == 6
    assert [axes.get_title() for axes in figure.axes] == [f"V{index}" for index in range(1, 7)]


def test_variant_contact_sheet_validates_its_inputs() -> None:
    relationship = _ks_relationship()
    variants = relationship.generate_variants()[:2]
    scenes = WorldScene3D.variant_scenes(relationship, variants=variants)
    with pytest.raises(ValueError, match="at least one scene"):
        render_variant_contact_sheet([])
    with pytest.raises(ValueError, match="columns must be strictly positive"):
        render_variant_contact_sheet(scenes, columns=0)
    with pytest.raises(ValueError, match="same length as scenes"):
        render_variant_contact_sheet(scenes, titles=["only one"])
    with pytest.raises(ValueError, match="controls 'title'"):
        render_variant_contact_sheet(scenes, title="nope")
