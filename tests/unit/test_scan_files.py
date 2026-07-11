from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.adapters import read_ang, read_ctf

SINGLE_PHASE_ANG = """\
# TEM_PIXperUM          1.000000
# x-star                0.507048
# WorkinDistance        22.000000
#
# Phase 1
# MaterialName  	Nickel
# Formula     	Ni
# Symmetry              43
# LatticeConstants      3.520 3.520 3.520  90.000  90.000  90.000
#
# GRID: SqrGrid
# XSTEP: 1.000000
# YSTEP: 1.000000
# NCOLS_ODD: 2
# NCOLS_EVEN: 2
# NROWS: 2
# OPERATOR: 	tester
#
   0.00000   0.00000   0.00000      0.00000      0.00000  60.0  0.950  0  1  0.500
   0.78540   0.52360   0.26180      1.00000      0.00000  55.0  0.900  0  1  0.600
   1.57080   0.78540   0.52360      0.00000      1.00000  50.0  0.850  0  1  0.700
   0.26180   1.04720   1.57080      1.00000      1.00000  45.0  0.800  0  1  0.800
"""

MULTI_PHASE_ANG = """\
# Phase 1
# MaterialName  	Nickel
# Formula     	Ni
# Symmetry              43
# Phase 2
# MaterialName  	Zirconium
# Formula     	Zr
# Symmetry              62
# GRID: SqrGrid
# XSTEP: 1.000000
# YSTEP: 1.000000
   0.10000   0.20000   0.30000      0.00000      0.00000  60.0  0.950  1  1  0.500
   0.40000   0.50000   0.60000      1.00000      0.00000  55.0  0.900  2  1  0.600
   0.70000   0.80000   0.90000      0.00000      1.00000  50.0  0.850  1  1  0.700
   0.00000   0.00000   0.00000      1.00000      1.00000   5.0  0.010  0  1  0.900
"""

TWO_PHASE_CTF = (
    "Channel Text File\n"
    "Prj\ttest-project\n"
    "Author\ttester\n"
    "JobMode\tGrid\n"
    "XCells\t2\n"
    "YCells\t2\n"
    "XStep\t0.5000\n"
    "YStep\t0.5000\n"
    "AcqE1\t0\n"
    "Euler angles refer to Sample Coordinate system (CS0)!\tMag\t500\n"
    "Phases\t2\n"
    "3.524;3.524;3.524\t90;90;90\tNickel\t11\t225\t\tcomment\n"
    "3.232;3.232;5.147\t90;90;120\tZirconium\t9\t194\t\n"
    "Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n"
    "1\t0.0\t0.0\t10\t0\t10.5\t20.5\t30.5\t0.40\t180\t200\n"
    "2\t0.5\t0.0\t9\t0\t40.0\t50.0\t60.0\t0.50\t170\t190\n"
    "0\t0.0\t0.5\t0\t3\t0.0\t0.0\t0.0\t0.00\t10\t20\n"
    "1\t0.5\t0.5\t11\t0\t70.0\t80.0\t90.0\t0.30\t190\t210\n"
)


def write(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_read_ang_single_phase(tmp_path: Path) -> None:
    result = read_ang(write(tmp_path, "scan.ang", SINGLE_PHASE_ANG))
    crystal_map = result.crystal_map
    assert len(crystal_map.orientations) == 4
    assert crystal_map.grid_shape == (2, 2)
    assert crystal_map.step_sizes == (1.0, 1.0)
    assert crystal_map.orientations.phase is None
    assert crystal_map.orientations.symmetry is not None
    assert crystal_map.orientations.symmetry.point_group == "m-3m"
    # .ang angles are radians; the second point is (45, 30, 15) degrees.
    euler_deg = crystal_map.orientations.as_bunge_euler(degrees=True)
    assert_allclose(euler_deg[1], [45.0, 30.0, 15.0], atol=1e-3)
    assert_allclose(crystal_map.coordinates[2], [0.0, 1.0], atol=1e-12)

    manifest = result.manifest
    assert manifest.source_system == "tsl_ang"
    assert manifest.phase_name == "Nickel"
    assert manifest.point_group == "m-3m"
    assert manifest.angle_unit == "radian"
    assert manifest.metadata["dropped_unindexed_points"] == "0"
    assert manifest.metadata["operator"] == "tester"

    assert_allclose(result.properties["image_quality"], [60.0, 55.0, 50.0, 45.0])
    assert_allclose(result.properties["confidence_index"], [0.95, 0.9, 0.85, 0.8])
    assert_allclose(result.properties["fit"], [0.5, 0.6, 0.7, 0.8])
    assert result.header_metadata["GRID"] == "SqrGrid"


def test_read_ang_multiphase_drops_unindexed_rows(tmp_path: Path) -> None:
    result = read_ang(write(tmp_path, "scan.ang", MULTI_PHASE_ANG))
    crystal_map = result.crystal_map
    assert len(crystal_map.orientations) == 3
    assert crystal_map.is_multiphase
    names = {entry.phase_id: entry.name for entry in crystal_map.phase_entries}
    assert names == {1: "Nickel", 2: "Zirconium"}
    groups = {entry.phase_id: entry.point_group for entry in crystal_map.phase_entries}
    assert groups == {1: "m-3m", 2: "6/mmm"}
    assert crystal_map.phase_ids is not None
    assert crystal_map.phase_ids.tolist() == [1, 2, 1]
    assert result.manifest.metadata["dropped_unindexed_points"] == "1"
    assert len(result.properties["confidence_index"]) == 3
    # grid_shape must not be claimed once rows were dropped
    assert crystal_map.grid_shape is None


def test_read_ang_rejects_hex_grids_and_unknown_symmetry(tmp_path: Path) -> None:
    hex_content = SINGLE_PHASE_ANG.replace("SqrGrid", "HexGrid")
    with pytest.raises(ValueError, match="HexGrid"):
        read_ang(write(tmp_path, "hex.ang", hex_content))
    bad_symmetry = SINGLE_PHASE_ANG.replace("# Symmetry              43", "# Symmetry  99")
    with pytest.raises(ValueError, match="symmetry code"):
        read_ang(write(tmp_path, "bad.ang", bad_symmetry))
    with pytest.raises(ValueError, match="no data rows"):
        read_ang(write(tmp_path, "empty.ang", "# GRID: SqrGrid\n"))


def test_read_ctf_two_phases(tmp_path: Path) -> None:
    result = read_ctf(write(tmp_path, "scan.ctf", TWO_PHASE_CTF))
    crystal_map = result.crystal_map
    assert len(crystal_map.orientations) == 3
    assert crystal_map.is_multiphase
    groups = {entry.name: entry.point_group for entry in crystal_map.phase_entries}
    assert groups == {"Nickel": "m-3m", "Zirconium": "6/mmm"}
    assert crystal_map.phase_ids is not None
    assert crystal_map.phase_ids.tolist() == [1, 2, 1]
    assert crystal_map.step_sizes == (0.5, 0.5)
    assert crystal_map.grid_shape is None

    euler_deg = crystal_map.orientations.as_bunge_euler(degrees=True)
    assert_allclose(euler_deg[0], [10.5, 20.5, 30.5], atol=1e-6)

    manifest = result.manifest
    assert manifest.source_system == "hkl_ctf"
    assert manifest.angle_unit == "degree"
    assert manifest.metadata["dropped_unindexed_points"] == "1"
    assert manifest.metadata["jobmode"] == "Grid"

    assert_allclose(result.properties["bands"], [10.0, 9.0, 11.0])
    assert_allclose(result.properties["mean_angular_deviation"], [0.4, 0.5, 0.3])
    assert_allclose(result.properties["band_contrast"], [180.0, 170.0, 190.0])


def test_read_ctf_rejects_malformed_files(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Channel Text File"):
        read_ctf(write(tmp_path, "bad.ctf", "not a ctf\n"))
    truncated = TWO_PHASE_CTF.split("Phase\tX")[0]
    with pytest.raises(ValueError, match="data header"):
        read_ctf(write(tmp_path, "trunc.ctf", truncated))
    only_unindexed = TWO_PHASE_CTF.replace("\n1\t", "\n0\t").replace("\n2\t", "\n0\t")
    with pytest.raises(ValueError, match="non-indexed"):
        read_ctf(write(tmp_path, "empty.ctf", only_unindexed))


def test_readers_feed_existing_map_workflows(tmp_path: Path) -> None:
    result = read_ang(write(tmp_path, "scan.ang", SINGLE_PHASE_ANG))
    crystal_map = result.crystal_map
    kam = crystal_map.kernel_average_misorientation_deg()
    assert kam.shape == (2, 2)
    assert np.all(np.isfinite(kam))
    report = crystal_map.summary()
    assert report["point_count"] == 4
