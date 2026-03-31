from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    CrystalMapPhase,
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def _foundation() -> tuple[ReferenceFrame, ReferenceFrame]:
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


def _phase(name: str, point_group: str, a: float) -> Phase:
    crystal, _ = _foundation()
    symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=crystal)
    return Phase(name=name, lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def _multiphase_map() -> tuple[CrystalMap, Phase, Phase]:
    crystal, specimen = _foundation()
    fcc = _phase("fcc_demo", "m-3m", 3.5)
    bcc = _phase("bcc_demo", "m-3m", 2.9)
    orientations = OrientationSet.from_quaternions(
        np.array(
            [
                Rotation.identity().quaternion,
                Rotation.from_bunge_euler(10.0, 0.0, 0.0).quaternion,
                Rotation.identity().quaternion,
                Rotation.from_bunge_euler(5.0, 0.0, 0.0).quaternion,
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        phase_entries=(
            CrystalMapPhase(phase_id=0, name=fcc.name, symmetry=fcc.symmetry, phase=fcc),
            CrystalMapPhase(phase_id=1, name=bcc.name, symmetry=bcc.symmetry, phase=bcc),
        ),
        phase_ids=np.array([0, 0, 1, 1], dtype=np.int64),
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    return crystal_map, fcc, bcc


def test_multiphase_crystal_map_tracks_phase_assignments_and_selection() -> None:
    crystal_map, fcc, _ = _multiphase_map()
    assert crystal_map.is_multiphase
    assert crystal_map.phase_summary() == {"fcc_demo": 2, "bcc_demo": 2}
    assert crystal_map.phase_entry_for_index(0) is not None
    assert crystal_map.phase_entry_for_index(0).name == "fcc_demo"
    selected = crystal_map.select_phase("fcc_demo")
    assert len(selected.orientations) == 2
    assert selected.orientations.phase == fcc
    assert selected.grid_shape is None
    assert any("graph mode" in note for note in selected.validate())


def test_multiphase_texture_outputs_require_phase_selector() -> None:
    crystal_map, _, _ = _multiphase_map()
    with pytest.raises(ValueError):
        crystal_map.to_odf()
    odf = crystal_map.to_odf(phase="fcc_demo")
    assert odf.orientations.phase is not None
    assert odf.orientations.phase.name == "fcc_demo"


def test_multiphase_segmentation_does_not_merge_across_phase_boundaries() -> None:
    crystal_map, _, _ = _multiphase_map()
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=15.0,
        symmetry_aware=True,
        connectivity=4,
    )
    assert len(segmentation.grains) == 2
    assert np.array_equal(segmentation.labels, np.array([0, 0, 1, 1], dtype=np.int64))


def test_neighbor_graph_supports_irregular_coordinate_mode() -> None:
    crystal, specimen = _foundation()
    phase = _phase("single", "432", 3.0)
    orientations = OrientationSet.from_quaternions(
        np.array(
            [
                Rotation.identity().quaternion,
                Rotation.from_bunge_euler(2.0, 0.0, 0.0).quaternion,
                Rotation.from_bunge_euler(4.0, 0.0, 0.0).quaternion,
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [0.8, 0.0], [0.0, 0.9]], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
    )
    graph = crystal_map.neighbor_graph(connectivity=4, order=1)
    assert graph.mode == "coordinate_radius"
    assert graph.pairs.shape[1] == 2
    assert graph.distances.shape[0] == graph.pairs.shape[0]


def test_neighbor_graph_includes_diagonals_for_eight_connectivity() -> None:
    crystal, specimen = _foundation()
    phase = _phase("single", "432", 3.0)
    orientations = OrientationSet.from_quaternions(
        np.array(
            [
                Rotation.identity().quaternion,
                Rotation.identity().quaternion,
                Rotation.identity().quaternion,
                Rotation.from_bunge_euler(10.0, 0.0, 0.0).quaternion,
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    neighbor_pairs = crystal_map.neighbor_graph(connectivity=8, order=1).pairs
    pair_set = {tuple(pair.tolist()) for pair in neighbor_pairs}
    assert (0, 3) in pair_set
    assert (1, 2) in pair_set
    kam_4 = crystal_map.kernel_average_misorientation_deg(connectivity=4, symmetry_aware=False)
    kam_8 = crystal_map.kernel_average_misorientation_deg(connectivity=8, symmetry_aware=False)
    assert kam_4[0, 0] == pytest.approx(0.0)
    assert kam_8[0, 0] > 0.0
