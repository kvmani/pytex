from __future__ import annotations

import subprocess
import sys

from pytex import FrameDomain, Handedness, Phase, ReferenceFrame

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


def test_python_m_pytex_info_runs() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "pytex", "info"],
        check=True,
        capture_output=True,
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
