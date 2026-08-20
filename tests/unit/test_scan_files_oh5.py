"""The EDAX OIM HDF5 scan reader (`.oh5` / `.h5`).

The fixtures here are written with `h5py` rather than tracked as binary blobs:
the layout under test *is* the assertion, so a file built in the test states it
where a reader can be checked against it, and nothing regenerable enters the
repository. The layout mirrors what OIM Analysis 8.6 writes, verified against a
real ``.oh5`` and the ``.ang`` OIM exported from the same scan.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.adapters import SCAN_FILE_SUFFIXES, read_oh5, read_scan, scan_reader_for

h5py = pytest.importorskip("h5py", reason="the .oh5/.h5 reader needs the optional 'hdf5' extra")

# A 2x3 square scan: three columns, two rows, 0.5 um steps.
SQUARE_EULER = np.array(
    [
        [0.10, 0.20, 0.30],
        [0.40, 0.50, 0.60],
        [0.70, 0.80, 0.90],
        [1.00, 1.10, 1.20],
        [1.30, 1.40, 1.50],
        [1.60, 1.70, 1.80],
    ],
    dtype=np.float64,
)
SQUARE_X = np.array([0.0, 0.5, 1.0, 0.0, 0.5, 1.0], dtype=np.float64)
SQUARE_Y = np.array([0.0, 0.0, 0.0, 0.5, 0.5, 0.5], dtype=np.float64)


def write_oh5(
    path: Path,
    *,
    scan_name: str = "Scan 1",
    grid_type: str = "SqrGrid",
    n_rows: int = 2,
    n_columns: int = 3,
    step_x: float = 0.5,
    step_y: float = 0.5,
    euler: np.ndarray = SQUARE_EULER,
    x_positions: np.ndarray = SQUARE_X,
    y_positions: np.ndarray = SQUARE_Y,
    phase_channel: np.ndarray | None = None,
    phases: tuple[tuple[str, str, int], ...] = (("Nickel", "Ni", 43),),
    channels: dict[str, np.ndarray] | None = None,
    extra_top_level: bool = True,
) -> Path:
    """Write a file shaped like an OIM Analysis HDF5 scan export.

    Parameters
    ----------
    phases : tuple of (str, str, int)
        Material name, formula, and ``LGsymID`` for each declared phase, in
        header order. OIM names the phase subgroups ``1``, ``2``, ...
    phase_channel : ndarray, optional
        The ``Phase`` per-point channel. A zero-based index into ``phases``
        with ``-1`` for unindexed points is written when omitted, which is what
        OIM 8 writes.
    """

    count = len(euler)
    if phase_channel is None:
        phase_channel = np.zeros(count, dtype=np.int8)
    with h5py.File(path, "w") as handle:
        if extra_top_level:
            handle.create_dataset("Manufacturer", data=np.bytes_(["EDAX"]))
            handle.create_dataset("Version", data=np.bytes_(["OIM Analysis 8.6.107 x64"]))
        header = handle.create_group(f"{scan_name}/EBSD/Header")
        data = handle.create_group(f"{scan_name}/EBSD/Data")
        header.create_dataset("Grid Type", data=np.bytes_([grid_type]))
        header.create_dataset("nRows", data=np.array([n_rows], dtype=np.int32))
        header.create_dataset("nColumns", data=np.array([n_columns], dtype=np.int32))
        header.create_dataset("Step X", data=np.array([step_x], dtype=np.float32))
        header.create_dataset("Step Y", data=np.array([step_y], dtype=np.float32))
        header.create_dataset("Sample Tilt", data=np.array([70.0], dtype=np.float64))
        header.create_dataset("Notes", data=np.bytes_(["written by a test"]))
        for index, (material, formula, symmetry_code) in enumerate(phases, start=1):
            group = header.create_group(f"Phase/{index}")
            group.create_dataset("MaterialName", data=np.bytes_([material]))
            group.create_dataset("Formula", data=np.bytes_([formula]))
            group.create_dataset("LGsymID", data=np.array([symmetry_code], dtype=np.int32))
        data.create_dataset("Phi1", data=euler[:, 0].astype(np.float32))
        data.create_dataset("Phi", data=euler[:, 1].astype(np.float32))
        data.create_dataset("Phi2", data=euler[:, 2].astype(np.float32))
        data.create_dataset("X Position", data=x_positions.astype(np.float32))
        data.create_dataset("Y Position", data=y_positions.astype(np.float32))
        data.create_dataset("Phase", data=np.asarray(phase_channel, dtype=np.int8))
        for name, values in (channels or {}).items():
            data.create_dataset(name, data=values)
    return path


def default_channels(count: int = 6) -> dict[str, np.ndarray]:
    return {
        "IQ": np.linspace(100.0, 600.0, count, dtype=np.float32),
        "CI": np.linspace(0.9, 0.4, count, dtype=np.float32),
        "Fit": np.linspace(0.5, 1.0, count, dtype=np.float32),
        "SEM Signal": np.arange(count, dtype=np.int32),
        "PRIAS Center Square": np.linspace(1.0, 2.0, count, dtype=np.float32),
        "Valid": np.zeros(count, dtype=np.int8),
    }


def test_read_oh5_single_phase_square_grid(tmp_path: Path) -> None:
    path = write_oh5(tmp_path / "scan.oh5", channels=default_channels())
    result = read_oh5(path)
    crystal_map = result.crystal_map

    assert len(crystal_map.orientations) == 6
    assert crystal_map.grid_shape == (2, 3)
    assert_allclose(crystal_map.step_sizes, (0.5, 0.5), atol=1e-6)
    assert crystal_map.orientations.symmetry is not None
    assert crystal_map.orientations.symmetry.point_group == "m-3m"
    assert_allclose(np.asarray(crystal_map.coordinates)[:, 0], SQUARE_X, atol=1e-6)
    # The angles are Bunge radians in the file, as in the .ang export, so they
    # arrive without a degree conversion anywhere in between.
    assert_allclose(
        crystal_map.orientations.as_bunge_euler(degrees=False),
        SQUARE_EULER,
        atol=1e-5,
    )


def test_read_oh5_names_channels_as_the_ang_reader_does(tmp_path: Path) -> None:
    path = write_oh5(tmp_path / "scan.oh5", channels=default_channels())
    result = read_oh5(path)

    # The four channels a .ang row also carries keep read_ang's names, so a
    # workflow cannot tell which of the two EDAX exports fed it...
    assert {"image_quality", "confidence_index", "detector_signal", "fit"} <= set(
        result.crystal_map.property_names
    )
    assert_allclose(result.properties["confidence_index"], np.linspace(0.9, 0.4, 6), atol=1e-6)
    assert_allclose(result.properties["image_quality"], np.linspace(100.0, 600.0, 6), atol=1e-4)
    # ...and the channels only the HDF5 export carries come through too, which
    # is the reason to prefer it.
    assert "prias_center_square" in result.crystal_map.property_names
    assert "valid" in result.crystal_map.property_names
    assert result.header_metadata["Notes"] == "written by a test"


def test_read_oh5_and_read_h5_are_one_format(tmp_path: Path) -> None:
    """`.oh5` and `.h5` are the same container under two extensions."""

    oh5 = write_oh5(tmp_path / "scan.oh5", channels=default_channels())
    h5 = write_oh5(tmp_path / "scan.h5", channels=default_channels())

    from_oh5 = read_scan(oh5).crystal_map
    from_h5 = read_scan(h5).crystal_map
    assert_allclose(
        from_oh5.orientations.as_matrices(),
        from_h5.orientations.as_matrices(),
        atol=0.0,
    )
    assert from_oh5.grid_shape == from_h5.grid_shape


def test_read_oh5_multiphase_drops_points_outside_the_declared_phases(tmp_path: Path) -> None:
    path = write_oh5(
        tmp_path / "two.oh5",
        phases=(("Nickel", "Ni", 43), ("Zirconium", "Zr", 62)),
        phase_channel=np.array([0, 1, 0, -1, 1, 0], dtype=np.int8),
    )
    result = read_oh5(path)
    crystal_map = result.crystal_map

    assert len(crystal_map.orientations) == 5
    assert result.manifest.metadata["dropped_unindexed_points"] == "1"
    # The zero-based Phase channel is resolved against the header's one-based
    # phase-group names, so a phase keeps the number the header gave it.
    assert [entry.phase_id for entry in crystal_map.phase_entries] == [1, 2]
    assert [entry.name for entry in crystal_map.phase_entries] == ["Nickel", "Zirconium"]
    assert [entry.point_group for entry in crystal_map.phase_entries] == ["m-3m", "6/mmm"]
    assert list(np.asarray(crystal_map.phase_ids)) == [1, 2, 1, 2, 1]
    # A row with a hole is no longer the vendor's complete topology.
    assert crystal_map.grid_shape is None


def test_read_oh5_reads_a_one_based_phase_channel(tmp_path: Path) -> None:
    """Older OIM writers number phases from one and mark unindexed points zero."""

    path = write_oh5(
        tmp_path / "one_based.oh5",
        phases=(("Nickel", "Ni", 43), ("Zirconium", "Zr", 62)),
        phase_channel=np.array([1, 2, 1, 0, 2, 1], dtype=np.int8),
    )
    result = read_oh5(path)

    assert result.manifest.metadata["phase_column_base"] == "1"
    assert len(result.crystal_map.orientations) == 5
    assert list(np.asarray(result.crystal_map.phase_ids)) == [1, 2, 1, 2, 1]


def test_read_oh5_keeps_every_point_of_a_single_phase_scan(tmp_path: Path) -> None:
    """Matching `read_ang`, so one scan imports the same way from either format."""

    path = write_oh5(
        tmp_path / "unindexed.oh5",
        phase_channel=np.array([0, 0, 0, 0, -1, -1], dtype=np.int8),
        channels=default_channels(),
    )
    result = read_oh5(path)

    assert len(result.crystal_map.orientations) == 6
    assert result.manifest.metadata["dropped_unindexed_points"] == "0"
    assert result.crystal_map.grid_shape == (2, 3)


def test_read_oh5_recovers_hexagonal_row_lengths(tmp_path: Path) -> None:
    """A HexGrid header states no odd/even split, so the rows are counted."""

    x_positions = np.array([0.0, 1.0, 2.0, 0.5, 1.5, 0.0, 1.0, 2.0], dtype=np.float64)
    y_positions = np.array([0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], dtype=np.float64)
    euler = np.tile(np.array([0.1, 0.2, 0.3]), (8, 1))
    path = write_oh5(
        tmp_path / "hex.oh5",
        grid_type="HexGrid",
        n_rows=3,
        n_columns=3,
        step_y=0.5,
        euler=euler,
        x_positions=x_positions,
        y_positions=y_positions,
    )
    crystal_map = read_oh5(path).crystal_map

    assert crystal_map.grid_kind == "hexagonal"
    assert crystal_map.row_lengths == (3, 2, 3)


def test_read_oh5_rejects_a_hexagonal_scan_with_uneven_rows(tmp_path: Path) -> None:
    y_positions = np.array([0.0, 0.0, 0.0, 0.5, 0.5, 1.0], dtype=np.float64)
    path = write_oh5(
        tmp_path / "ragged.oh5",
        grid_type="HexGrid",
        n_rows=3,
        step_y=0.5,
        y_positions=y_positions,
    )
    with pytest.raises(ValueError, match="distinct"):
        read_oh5(path)


def test_read_oh5_reads_a_named_scan_and_rejects_an_unknown_one(tmp_path: Path) -> None:
    path = write_oh5(tmp_path / "named.oh5", scan_name="Scan 7")
    assert len(read_oh5(path, scan="Scan 7").crystal_map.orientations) == 6
    assert read_oh5(path).manifest.metadata["scan"] == "Scan 7"
    with pytest.raises(ValueError, match="no EBSD scan named"):
        read_oh5(path, scan="Scan 9")


def test_read_oh5_rejects_files_it_cannot_interpret(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="not supported"):
        read_oh5(write_oh5(tmp_path / "grid.oh5", grid_type="RadialGrid"))
    with pytest.raises(ValueError, match="Unsupported TSL symmetry code"):
        read_oh5(write_oh5(tmp_path / "symmetry.oh5", phases=(("Mystery", "Xx", 99),)))

    empty = tmp_path / "empty.oh5"
    with h5py.File(empty, "w") as handle:
        handle.create_dataset("Manufacturer", data=np.bytes_(["EDAX"]))
    with pytest.raises(ValueError, match="no group holding an EBSD scan"):
        read_oh5(empty)

    unindexed = tmp_path / "unindexed.oh5"
    with h5py.File(unindexed, "w") as handle:
        header = handle.create_group("Scan 1/EBSD/Header")
        handle.create_group("Scan 1/EBSD/Data")
        header.create_dataset("nRows", data=np.array([1], dtype=np.int32))
        header.create_dataset("nColumns", data=np.array([1], dtype=np.int32))
        group = header.create_group("Phase/1")
        group.create_dataset("MaterialName", data=np.bytes_(["Nickel"]))
        group.create_dataset("LGsymID", data=np.array([43], dtype=np.int32))
    with pytest.raises(ValueError, match="not an indexed OIM EBSD scan"):
        read_oh5(unindexed)


def test_read_oh5_manifest_records_its_provenance(tmp_path: Path) -> None:
    path = write_oh5(tmp_path / "scan.oh5")
    manifest = read_oh5(path).manifest

    assert manifest.source_system == "edax_oh5"
    assert manifest.source_file == str(path)
    assert manifest.metadata["reader"] == "pytex.adapters.scan_files.read_oh5"
    assert manifest.metadata["manufacturer"] == "EDAX"
    assert manifest.metadata["version"].startswith("OIM Analysis")
    assert manifest.metadata["grid"] == "SqrGrid"
    assert manifest.metadata["frame_policy"] == "vendor specimen frame, no axis remapping applied"


def test_read_scan_dispatches_on_the_extension(tmp_path: Path) -> None:
    from pytex.adapters.scan_files import read_ang, read_ctf

    assert SCAN_FILE_SUFFIXES == (".ang", ".ctf", ".oh5", ".h5")
    assert scan_reader_for("a.ang") is read_ang
    assert scan_reader_for("a.CTF") is read_ctf
    assert scan_reader_for("a.oh5") is read_oh5
    assert scan_reader_for(tmp_path / "a.h5") is read_oh5
    with pytest.raises(ValueError, match="not a scan file PyTex reads"):
        scan_reader_for("a.txt")


def test_read_scan_passes_its_keywords_through(tmp_path: Path) -> None:
    """A keyword the chosen reader accepts reaches it through the dispatcher."""

    from pytex.adapters.scan_files import default_ebsd_frames
    from pytex.core.lattice import Lattice, Phase
    from pytex.core.symmetry import SymmetrySpec

    crystal_frame, _, _ = default_ebsd_frames()
    supplied: Any = Phase(
        "Nickel",
        lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal_frame),
        crystal_frame=crystal_frame,
    )
    result = read_scan(write_oh5(tmp_path / "scan.oh5"), phase=supplied)

    assert result.crystal_map.orientations.phase is supplied
