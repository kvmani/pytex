"""Tests for the standalone CIF-to-HRTEM-series driver."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from scripts.hrtem_cif_series import inclusive_range, load_cif_snapshot


def _write_cif(path: Path, *, gamma: float = 90.0) -> None:
    path.write_text(
        f"""data_md_snapshot
_cell_length_a 4.0
_cell_length_b 6.0
_cell_length_c 8.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma {gamma}
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0.25 0.50 0.75
Fe2 Fe 0.75 0.25 0.25
""",
        encoding="utf-8",
    )


def test_inclusive_range_includes_stop_without_overshoot() -> None:
    assert np.array_equal(inclusive_range(-50.0, 50.0, 25.0), [-50, -25, 0, 25, 50])


def test_load_cif_snapshot_puts_selected_cell_axis_along_beam(tmp_path: Path) -> None:
    path = tmp_path / "snapshot.cif"
    _write_cif(path)

    snapshot = load_cif_snapshot(path, "a")

    assert snapshot.species == ("Fe", "Fe")
    assert np.allclose(np.diag(snapshot.cell), (6.0, 8.0, 4.0))
    assert snapshot.periodicity == (True, True, False)
    assert np.all(snapshot.positions[:, :2] >= 0.0)
    assert np.all(snapshot.positions[:, :2] < np.diag(snapshot.cell)[:2])


def test_load_cif_snapshot_rejects_nonorthogonal_cell(tmp_path: Path) -> None:
    path = tmp_path / "triclinic.cif"
    _write_cif(path, gamma=85.0)

    with pytest.raises(ValueError, match="not orthogonal"):
        load_cif_snapshot(path, "c")
