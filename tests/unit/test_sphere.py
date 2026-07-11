from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    FrameDomain,
    Handedness,
    ReferenceFrame,
    Rotation,
    S2Grid,
    SphericalVectorSet,
    VectorSet,
)


def make_frame(name: str = "specimen") -> ReferenceFrame:
    return ReferenceFrame(
        name=name,
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )


def test_from_polar_matches_expected_directions() -> None:
    frame = make_frame()
    directions = SphericalVectorSet.from_polar(
        [0.0, 90.0, 90.0],
        [0.0, 0.0, 90.0],
        reference_frame=frame,
    )
    expected = np.array(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    assert_allclose(directions.as_array(), expected, atol=1e-12)


def test_polar_round_trip_in_degrees_and_radians() -> None:
    frame = make_frame()
    polar = np.array([15.0, 60.0, 120.0])
    azimuth = np.array([10.0, 200.0, 350.0])
    directions = SphericalVectorSet.from_polar(polar, azimuth, reference_frame=frame)
    polar_deg, azimuth_deg = directions.to_polar()
    assert_allclose(polar_deg, polar, atol=1e-10)
    assert_allclose(azimuth_deg, azimuth, atol=1e-10)
    polar_rad, azimuth_rad = directions.to_polar(degrees=False)
    assert_allclose(polar_rad, np.deg2rad(polar), atol=1e-12)
    assert_allclose(azimuth_rad, np.deg2rad(azimuth), atol=1e-12)


def test_construction_normalizes_and_rejects_zero_vectors() -> None:
    frame = make_frame()
    directions = SphericalVectorSet.from_vectors(
        [[2.0, 0.0, 0.0], [0.0, 0.0, -5.0]],
        reference_frame=frame,
    )
    assert_allclose(np.linalg.norm(directions.as_array(), axis=1), [1.0, 1.0])
    with pytest.raises(ValueError):
        SphericalVectorSet.from_vectors([[0.0, 0.0, 0.0]], reference_frame=frame)


def test_vector_set_round_trip_preserves_frame() -> None:
    frame = make_frame()
    vector_set = VectorSet(values=[[0.0, 3.0, 4.0]], reference_frame=frame)
    directions = SphericalVectorSet.from_vector_set(vector_set, antipodal=True)
    assert directions.reference_frame == frame
    assert directions.antipodal is True
    round_trip = directions.to_vector_set()
    assert round_trip.reference_frame == frame
    assert_allclose(round_trip.values, [[0.0, 0.6, 0.8]], atol=1e-12)


def test_getitem_and_subset_semantics() -> None:
    frame = make_frame()
    directions = SphericalVectorSet.from_polar(
        [0.0, 45.0, 90.0],
        [0.0, 90.0, 180.0],
        reference_frame=frame,
        antipodal=True,
    )
    single = directions[1]
    assert isinstance(single, np.ndarray)
    assert single.shape == (3,)
    sliced = directions[0:2]
    assert isinstance(sliced, SphericalVectorSet)
    assert len(sliced) == 2
    assert sliced.antipodal is True
    subset = directions.subset([2, 0])
    assert len(subset) == 2
    assert_allclose(subset.as_array()[1], [0.0, 0.0, 1.0], atol=1e-12)


def test_angles_respect_antipodal_equivalence() -> None:
    frame = make_frame()
    plus_z = SphericalVectorSet.from_vectors([[0.0, 0.0, 1.0]], reference_frame=frame)
    minus_z = SphericalVectorSet.from_vectors([[0.0, 0.0, -1.0]], reference_frame=frame)
    assert_allclose(plus_z.angles_to_rad(minus_z), [np.pi], atol=1e-12)
    minus_z_antipodal = SphericalVectorSet.from_vectors(
        [[0.0, 0.0, -1.0]],
        reference_frame=frame,
        antipodal=True,
    )
    assert_allclose(plus_z.angles_to_rad(minus_z_antipodal), [0.0], atol=1e-12)
    assert_allclose(plus_z.angles_to_deg(minus_z_antipodal), [0.0], atol=1e-12)


def test_binary_operations_broadcast_and_require_matching_frames() -> None:
    frame = make_frame()
    other_frame = make_frame("lab")
    x_axis = SphericalVectorSet.from_vectors([[1.0, 0.0, 0.0]], reference_frame=frame)
    axes = SphericalVectorSet.from_vectors(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        reference_frame=frame,
    )
    assert_allclose(x_axis.dot(axes), [0.0, 0.0], atol=1e-12)
    crossed = x_axis.cross(axes)
    assert_allclose(crossed.as_array(), [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], atol=1e-12)
    mismatched = SphericalVectorSet.from_vectors(
        [[1.0, 0.0, 0.0]],
        reference_frame=other_frame,
    )
    with pytest.raises(ValueError):
        x_axis.dot(mismatched)
    with pytest.raises(ValueError):
        x_axis.cross(x_axis)
    three = SphericalVectorSet.from_polar([10.0, 20.0, 30.0], 0.0, reference_frame=frame)
    two = SphericalVectorSet.from_polar([10.0, 20.0], 0.0, reference_frame=frame)
    with pytest.raises(ValueError):
        three.dot(two)


def test_fold_upper_hemisphere_requires_antipodal_semantics() -> None:
    frame = make_frame()
    directions = SphericalVectorSet.from_vectors(
        [[0.0, 0.0, -1.0], [0.0, 1.0, 0.5], [-1.0, 0.0, 0.0]],
        reference_frame=frame,
        antipodal=True,
    )
    folded = directions.fold_upper_hemisphere()
    assert np.all(folded.as_array()[:, 2] >= -1e-12)
    assert_allclose(folded.as_array()[0], [0.0, 0.0, 1.0], atol=1e-12)
    assert_allclose(folded.as_array()[2], [1.0, 0.0, 0.0], atol=1e-12)
    non_antipodal = SphericalVectorSet.from_vectors(
        [[0.0, 0.0, -1.0]],
        reference_frame=frame,
    )
    with pytest.raises(ValueError):
        non_antipodal.fold_upper_hemisphere()


def test_mean_direction_for_clustered_directions() -> None:
    frame = make_frame()
    clustered = SphericalVectorSet.from_polar(
        [5.0, 5.0, 5.0, 5.0],
        [0.0, 90.0, 180.0, 270.0],
        reference_frame=frame,
    )
    assert_allclose(clustered.mean_direction(), [0.0, 0.0, 1.0], atol=1e-12)
    balanced = SphericalVectorSet.from_vectors(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        reference_frame=frame,
    )
    with pytest.raises(ValueError):
        balanced.mean_direction()


def test_mean_direction_for_antipodal_axes_uses_orientation_tensor() -> None:
    frame = make_frame()
    axes = SphericalVectorSet.from_vectors(
        [[1.0, 0.0, 0.05], [-1.0, 0.0, 0.05], [1.0, 0.0, -0.05], [-1.0, 0.0, -0.05]],
        reference_frame=frame,
        antipodal=True,
    )
    mean_axis = axes.mean_direction()
    assert_allclose(np.abs(mean_axis), [1.0, 0.0, 0.0], atol=1e-12)
    tensor = axes.orientation_tensor()
    assert_allclose(tensor, tensor.T, atol=1e-12)
    assert_allclose(np.trace(tensor), 1.0, atol=1e-12)


def test_rotated_by_maps_directions_within_the_frame() -> None:
    frame = make_frame()
    directions = SphericalVectorSet.from_vectors(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        reference_frame=frame,
        antipodal=True,
    )
    rotation = Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 2.0)
    rotated = directions.rotated_by(rotation)
    assert rotated.reference_frame == frame
    assert rotated.antipodal is True
    assert_allclose(rotated.as_array(), [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], atol=1e-12)


def test_equispaced_grid_covers_hemisphere_with_normalized_weights() -> None:
    frame = make_frame()
    grid = S2Grid.equispaced(10.0, reference_frame=frame)
    assert grid.method == "equispaced"
    assert grid.hemisphere == "upper"
    assert len(grid) == len(grid.vectors)
    assert_allclose(np.linalg.norm(grid.vectors.as_array(), axis=1), 1.0, atol=1e-12)
    assert np.all(grid.vectors.as_array()[:, 2] >= -1e-12)
    assert_allclose(float(grid.weights.sum()), 1.0, atol=1e-12)
    assert np.all(grid.weights > 0.0)
    polar_deg, _ = grid.vectors.to_polar()
    assert np.isclose(polar_deg.min(), 0.0)
    assert np.isclose(polar_deg.max(), 90.0)


def test_equispaced_grid_full_sphere_reaches_south_pole() -> None:
    frame = make_frame()
    grid = S2Grid.equispaced(15.0, reference_frame=frame, hemisphere="sphere")
    polar_deg, _ = grid.vectors.to_polar()
    assert np.isclose(polar_deg.max(), 180.0)
    assert_allclose(float(grid.weights.sum()), 1.0, atol=1e-12)


def test_equispaced_grid_density_scales_with_resolution() -> None:
    frame = make_frame()
    coarse = S2Grid.equispaced(15.0, reference_frame=frame)
    fine = S2Grid.equispaced(5.0, reference_frame=frame)
    assert len(fine) > 4 * len(coarse)


def test_regular_grid_deduplicates_poles_and_weights_by_band() -> None:
    frame = make_frame()
    grid = S2Grid.regular(30.0, 45.0, reference_frame=frame)
    polar_deg, _ = grid.vectors.to_polar()
    pole_points = np.count_nonzero(np.isclose(polar_deg, 0.0))
    assert pole_points == 1
    ring_points = np.count_nonzero(np.isclose(polar_deg, 30.0))
    assert ring_points == 8
    assert_allclose(float(grid.weights.sum()), 1.0, atol=1e-12)
    near_equator = np.isclose(polar_deg, 60.0)
    near_pole = np.isclose(polar_deg, 30.0)
    assert grid.weights[near_equator][0] > grid.weights[near_pole][0]


def test_grid_validation_rejects_bad_inputs() -> None:
    frame = make_frame()
    with pytest.raises(ValueError):
        S2Grid.equispaced(0.0, reference_frame=frame)
    with pytest.raises(ValueError):
        S2Grid.equispaced(120.0, reference_frame=frame)
    with pytest.raises(ValueError):
        S2Grid.equispaced(10.0, reference_frame=frame, hemisphere="lower")
    with pytest.raises(ValueError):
        S2Grid.regular(30.0, 50.0, reference_frame=frame)
