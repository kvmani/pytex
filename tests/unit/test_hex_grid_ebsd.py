from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pytex import CrystalMap, from_json_contract, to_json_contract
from pytex.adapters import read_ang
from pytex.plotting import plot_kam_map

_FIXTURE = Path(__file__).parents[2] / "fixtures" / "ebsd" / "synthetic_hex_grid.ang"

_FIRST_SHELL_PAIRS = np.array(
    [
        [0, 1],
        [0, 3],
        [1, 2],
        [1, 3],
        [1, 4],
        [2, 4],
        [3, 4],
        [3, 5],
        [3, 6],
        [4, 6],
        [4, 7],
        [5, 6],
        [6, 7],
    ],
    dtype=np.int64,
)


def test_hex_ang_import_preserves_ragged_topology_and_properties() -> None:
    result = read_ang(_FIXTURE)
    crystal_map = result.crystal_map

    assert crystal_map.grid_kind == "hexagonal"
    assert crystal_map.grid_shape is None
    assert crystal_map.row_lengths == (3, 2, 3)
    assert crystal_map.default_connectivity == 6
    assert_allclose(crystal_map.step_sizes, (1.0, np.sqrt(3.0) / 2.0), atol=5e-13)
    assert result.manifest.metadata["sampleid"] == "synthetic_hex_grid_not_experimental"
    image_quality = crystal_map.property_map("image_quality")
    assert image_quality.shape == (3, 3)
    assert np.isnan(image_quality[1, 2])
    assert_allclose(image_quality[~np.isnan(image_quality)], 60.0)


def test_hex_first_shell_and_kam_match_direct_analytic_count() -> None:
    crystal_map = read_ang(_FIXTURE).crystal_map
    graph = crystal_map.neighbor_graph()

    assert graph.mode == "hexagonal_grid"
    assert graph.connectivity == 6
    assert graph.order == 1
    assert_array_equal(graph.pairs, _FIRST_SHELL_PAIRS)
    assert_allclose(graph.distances, 1.0, atol=5e-10)
    second_shell = crystal_map.neighbor_graph(order=2)
    assert len(second_shell.pairs) == 26
    assert not np.any(np.all(second_shell.pairs == (0, 7), axis=1))

    # Only point 1 is rotated by 6 degrees. Its four incident edges are 6 deg;
    # every other edge is zero, so each KAM is 6 / local degree when incident.
    expected_kam = np.array([3.0, 6.0, 3.0, 1.2, 1.2, 0.0, 0.0, 0.0])
    assert_allclose(
        crystal_map.kernel_average_misorientation_deg(),
        expected_kam,
        atol=2e-7,
    )


def test_hex_segmentation_and_contract_keep_six_neighbor_semantics() -> None:
    crystal_map = read_ang(_FIXTURE).crystal_map
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)

    assert segmentation.connectivity == 6
    assert segmentation.grain_sizes() == {0: 7, 1: 1}
    assert segmentation.label_grid.shape == (3, 3)
    assert segmentation.label_grid[1, 2] == -1
    assert segmentation.grod_map_deg().shape == (8,)
    with pytest.raises(ValueError, match="rectangular pixel faces"):
        segmentation.grain_perimeters()

    restored = from_json_contract(to_json_contract(crystal_map))
    assert isinstance(restored, CrystalMap)
    assert restored.grid_kind == "hexagonal"
    assert restored.row_lengths == (3, 2, 3)
    assert_array_equal(restored.neighbor_pairs(), _FIRST_SHELL_PAIRS)
    assert_allclose(restored.get_property("confidence_index"), 0.95)
    assert "six-neighbour first-shell topology" in restored.describe()


def test_hex_topology_invariants_and_explicit_connectivity_errors() -> None:
    crystal_map = read_ang(_FIXTURE).crystal_map

    with pytest.raises(ValueError, match="exactly one entry"):
        replace(crystal_map, row_lengths=(3, 2))
    with pytest.raises(ValueError, match="requires connectivity=6"):
        crystal_map.neighbor_graph(connectivity=4)
    aligned_coordinates = np.array(crystal_map.coordinates, copy=True)
    aligned_coordinates[3:5, 0] = (0.0, 1.0)
    with pytest.raises(ValueError, match="do not form staggered"):
        replace(crystal_map, coordinates=aligned_coordinates)
    with pytest.raises(ValueError, match="require a 2D grid_shape"):
        crystal_map._require_regular_2d_grid()


def test_hex_kam_plot_uses_measured_coordinates() -> None:
    import matplotlib.pyplot as plt

    crystal_map = read_ang(_FIXTURE).crystal_map
    figure = plot_kam_map(crystal_map)
    try:
        axes = figure.axes[0]
        assert len(axes.collections) == 1
        assert_allclose(axes.collections[0].get_offsets(), crystal_map.coordinates)
    finally:
        plt.close(figure)
