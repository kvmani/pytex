from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from pytex import (
    FrameDomain,
    Handedness,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    from_orix_miller,
    from_orix_rotation,
    to_orix_miller_direction,
    to_orix_miller_plane,
    to_orix_rotation,
)

NACL_CIF = """
data_NaCl
_symmetry_space_group_name_H-M 'F m -3 m'
_cell_length_a 5.6402
_cell_length_b 5.6402
_cell_length_c 5.6402
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_Int_Tables_number 225
loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  Na Na 0.0 0.0 0.0
  Cl Cl 0.5 0.5 0.5
"""

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    paths = [str(SRC_ROOT), str(REPO_ROOT)]
    if existing:
        paths.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def test_python_m_pytex_info_runs() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "pytex", "info"],
        check=True,
        capture_output=True,
        env=_subprocess_env(),
        text=True,
    )
    assert "PyTex 0.1.0.dev0" in result.stdout
    assert "Canonical convention: pytex_canonical" in result.stdout


def test_phase_from_cif_string_optional_adapter_integration() -> None:
    pytest = __import__("pytest")
    pytest.importorskip("pymatgen.core")
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    phase = Phase.from_cif_string(NACL_CIF, crystal_frame=crystal, phase_name="nacl")
    assert phase.name == "nacl"
    assert phase.space_group_number == 225
    assert phase.unit_cell is not None


def test_orix_optional_adapter_boundary_preserves_core_semantics() -> None:
    pytest = __import__("pytest")
    pytest.importorskip("orix")
    pytest.importorskip("pymatgen.core")
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    phase = Phase.from_cif_string(NACL_CIF, crystal_frame=crystal, phase_name="nacl")
    rotation = Rotation.from_bunge_euler(35.0, 25.0, 10.0)
    recovered_rotation = from_orix_rotation(to_orix_rotation(rotation))
    assert recovered_rotation.quaternion.tolist() == pytest.approx(rotation.quaternion.tolist())

    plane = MillerPlane((1, 1, 1), phase=phase)
    recovered_plane = from_orix_miller(to_orix_miller_plane(plane), phase=phase)
    assert recovered_plane.phase.name == phase.name
    assert recovered_plane.indices.tolist() == [1, 1, 1]

    direction = MillerDirection((1, 0, 0), phase=phase)
    recovered_direction = from_orix_miller(to_orix_miller_direction(direction), phase=phase)
    assert recovered_direction.phase.name == phase.name
    assert recovered_direction.indices.tolist() == [1, 0, 0]
