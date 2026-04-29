from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    ipf_color,
    ipf_colors,
    specimen_direction_vector,
)


def make_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
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
    symmetry = SymmetrySpec.from_point_group("432", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def test_specimen_direction_aliases_resolve_to_canonical_vectors() -> None:
    assert_allclose(specimen_direction_vector("RD"), [1.0, 0.0, 0.0])
    assert_allclose(specimen_direction_vector("td"), [0.0, 1.0, 0.0])
    assert_allclose(specimen_direction_vector("z"), [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="Specimen direction"):
        specimen_direction_vector("bad")


def test_so2_grid_spacing_and_axis_invariant() -> None:
    _, specimen, phase = make_context()
    grid = OrientationSet.from_so2_grid(
        "ND",
        90.0,
        specimen_frame=specimen,
        phase=phase,
    )
    assert len(grid) == 4
    assert grid.phase == phase
    assert_allclose(grid[0].as_matrix(), np.eye(3), atol=1e-8)
    mapped_axis = grid.map_crystal_directions([0.0, 0.0, 1.0])
    assert_allclose(mapped_axis, np.tile([0.0, 0.0, 1.0], (4, 1)), atol=1e-8)
    with pytest.raises(ValueError, match="spacing_deg"):
        OrientationSet.from_so2_grid("ND", 0.0, specimen_frame=specimen, phase=phase)


def test_regular_so3_grid_uses_bunge_ranges_and_preserves_context() -> None:
    _, specimen, phase = make_context()
    grid = OrientationSet.from_regular_so3_grid(
        90.0,
        specimen_frame=specimen,
        phase=phase,
        phi1_range_deg=(0.0, 180.0),
        big_phi_range_deg=(0.0, 90.0),
        phi2_range_deg=(0.0, 180.0),
    )
    assert len(grid) == 8
    assert grid.phase == phase
    assert grid.symmetry == phase.symmetry
    for matrix in grid.as_matrices():
        assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-8)
        assert np.isclose(np.linalg.det(matrix), 1.0, atol=1e-8)


def test_equispaced_so3_grid_is_deterministic_normalized_and_deduplicated() -> None:
    _, specimen, phase = make_context()
    left = OrientationSet.from_equispaced_so3_grid(
        120.0,
        specimen_frame=specimen,
        phase=phase,
    )
    right = OrientationSet.from_equispaced_so3_grid(
        120.0,
        specimen_frame=specimen,
        phase=phase,
    )
    assert_allclose(left.quaternions, right.quaternions, atol=1e-12)
    assert_allclose(np.linalg.norm(left.quaternions, axis=1), 1.0, atol=1e-12)
    assert np.all(left.quaternions[:, 0] >= -1e-12)
    keys = np.round(left.exact_fundamental_region_keys(), decimals=10)
    assert np.unique(keys, axis=0).shape[0] == len(left)


def test_ipf_color_alias_wrappers_are_finite_and_vectorized() -> None:
    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        [[0.0, 0.0, 0.0], [35.0, 25.0, 10.0]],
        specimen_frame=specimen,
        phase=phase,
    )
    for direction in ("RD", "TD", "ND"):
        colors = ipf_colors(orientations, direction=direction)
        assert colors.shape == (2, 3)
        assert np.all(np.isfinite(colors))
        assert np.all((colors >= 0.0) & (colors <= 1.0))
        assert_allclose(ipf_color(orientations[0], direction=direction), colors[0], atol=1e-12)


def test_orientation_convenience_constructors_match_existing_surfaces() -> None:
    _, specimen, phase = make_context()
    from_euler = Orientation.from_euler(
        35.0,
        25.0,
        10.0,
        specimen_frame=specimen,
        phase=phase,
    )
    assert_allclose(
        from_euler.as_matrix(),
        Rotation.from_bunge_euler(35.0, 25.0, 10.0).as_matrix(),
        atol=1e-8,
    )
    from_axis = Orientation.from_axis_angle(
        "ND",
        np.pi / 2.0,
        specimen_frame=specimen,
        phase=phase,
    )
    assert_allclose(from_axis.map_crystal_vector([0.0, 0.0, 1.0]), [0.0, 0.0, 1.0])
    assert_allclose(
        Orientation.from_matrix(
            from_euler.as_matrix(),
            specimen_frame=specimen,
            phase=phase,
        ).as_matrix(),
        from_euler.as_matrix(),
        atol=1e-8,
    )
    assert_allclose(
        Orientation.from_quaternion(
            from_euler.rotation.quaternion,
            specimen_frame=specimen,
            phase=phase,
        ).as_matrix(),
        from_euler.as_matrix(),
        atol=1e-8,
    )


def test_orientation_from_miller_accepts_objects_and_raw_indices() -> None:
    _, specimen, phase = make_context()
    plane = CrystalPlane(MillerIndex([0, 0, 1], phase=phase), phase=phase)
    direction = CrystalDirection([1.0, 0.0, 0.0], phase=phase)
    from_objects = Orientation.from_miller(
        plane,
        direction,
        specimen_frame=specimen,
        specimen_plane_normal="ND",
        specimen_direction="RD",
    )
    from_raw = Orientation.from_miller(
        [0, 0, 1],
        [1, 0, 0],
        specimen_frame=specimen,
        phase=phase,
        specimen_plane_normal="ND",
        specimen_direction="RD",
    )
    assert_allclose(from_objects.as_matrix(), np.eye(3), atol=1e-8)
    assert_allclose(from_raw.as_matrix(), np.eye(3), atol=1e-8)


@pytest.mark.parametrize(
    ("point_group", "specimen_point_group"),
    [
        ("1", None),
        ("2", None),
        ("222", None),
        ("3", None),
        ("32", None),
        ("4", None),
        ("422", None),
        ("6", None),
        ("622", None),
        ("23", None),
        ("432", "222"),
    ],
)
def test_orientation_projection_is_idempotent_across_supported_point_groups(
    point_group: str,
    specimen_point_group: str | None,
) -> None:
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
    crystal_symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    specimen_symmetry = (
        None
        if specimen_point_group is None
        else SymmetrySpec.from_point_group(specimen_point_group, reference_frame=specimen)
    )
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("phase", lattice=lattice, symmetry=crystal_symmetry, crystal_frame=crystal)
    orientation = Orientation.from_euler(
        35.0,
        25.0,
        10.0,
        specimen_frame=specimen,
        phase=phase,
    )
    projected = orientation.project_to_exact_fundamental_region(
        specimen_symmetry=specimen_symmetry
    )
    reprojection = projected.project_to_exact_fundamental_region(
        specimen_symmetry=specimen_symmetry
    )
    assert projected.is_in_fundamental_region(specimen_symmetry=specimen_symmetry)
    assert_allclose(reprojection.as_matrix(), projected.as_matrix(), atol=1e-10)
