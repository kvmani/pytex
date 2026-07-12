from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def _frames() -> tuple[ReferenceFrame, ReferenceFrame]:
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
    return crystal, specimen


def _single_phase_map(properties: dict[str, np.ndarray] | None = None) -> CrystalMap:
    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.5, 3.5, 3.5, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    quaternions = np.array(
        [
            Rotation.identity().quaternion,
            Rotation.from_bunge_euler(10.0, 0.0, 0.0).quaternion,
            Rotation.from_bunge_euler(0.0, 5.0, 0.0).quaternion,
            Rotation.from_bunge_euler(5.0, 5.0, 0.0).quaternion,
        ]
    )
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    return CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
        properties=properties,
    )


def test_property_channels_default_to_empty() -> None:
    crystal_map = _single_phase_map()
    assert crystal_map.property_names == ()


def test_property_channels_round_trip_and_map_reshape() -> None:
    iq = np.array([10.0, 20.0, 30.0, 40.0])
    crystal_map = _single_phase_map({"image_quality": iq})
    assert crystal_map.property_names == ("image_quality",)
    np.testing.assert_allclose(crystal_map.get_property("image_quality"), iq)
    np.testing.assert_allclose(
        crystal_map.property_map("image_quality"), [[10.0, 20.0], [30.0, 40.0]]
    )
    # channels are read-only
    with pytest.raises(ValueError):
        crystal_map.get_property("image_quality")[0] = 99.0


def test_property_channel_length_is_validated() -> None:
    with pytest.raises(ValueError, match="one value per map point"):
        _single_phase_map({"bad": np.array([1.0, 2.0])})


def test_property_channel_names_must_be_non_empty_strings() -> None:
    with pytest.raises(ValueError, match="non-empty strings"):
        _single_phase_map({"  ": np.zeros(4)})


def test_with_properties_merges_and_replaces() -> None:
    crystal_map = _single_phase_map({"a": np.zeros(4)})
    merged = crystal_map.with_properties({"b": np.ones(4)})
    assert set(merged.property_names) == {"a", "b"}
    replaced = crystal_map.with_properties({"b": np.ones(4)}, replace=True)
    assert replaced.property_names == ("b",)
    # original is untouched (immutability)
    assert crystal_map.property_names == ("a",)


def test_get_property_reports_available_channels() -> None:
    crystal_map = _single_phase_map({"a": np.zeros(4)})
    with pytest.raises(KeyError, match="Available channels: a"):
        crystal_map.get_property("missing")


def test_property_threshold_mask_and_filter() -> None:
    crystal_map = _single_phase_map({"ci": np.array([0.1, 0.6, 0.9, 0.2])})
    mask = crystal_map.property_threshold_mask("ci", minimum=0.5)
    assert mask.tolist() == [False, True, True, False]
    filtered = crystal_map.filter_by_property("ci", minimum=0.5)
    assert len(filtered.orientations) == 2
    # channels are carried forward, masked to the retained points
    np.testing.assert_allclose(filtered.get_property("ci"), [0.6, 0.9])
    # dropping points leaves graph mode (no grid)
    assert filtered.grid_shape is None


def test_property_threshold_mask_requires_a_bound() -> None:
    crystal_map = _single_phase_map({"ci": np.array([0.1, 0.6, 0.9, 0.2])})
    with pytest.raises(ValueError, match="minimum/maximum"):
        crystal_map.property_threshold_mask("ci")
