from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

NACL_CIF = """
data_NaCl
_symmetry_space_group_name_H-M 'F m -3 m'
_symmetry_Int_Tables_number 225
_cell_length_a 5.6402
_cell_length_b 5.6402
_cell_length_c 5.6402
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  Na1 Na 0.0 0.0 0.0
  Cl1 Cl 0.5 0.5 0.5
"""


def make_phase() -> Phase:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    lattice = Lattice(3.0, 3.0, 5.0, 90.0, 90.0, 120.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    return Phase(name="hcp-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def test_direct_and_reciprocal_bases_are_dual() -> None:
    phase = make_phase()
    direct = phase.lattice.direct_basis().matrix
    reciprocal = phase.lattice.reciprocal_basis().matrix
    assert_allclose(reciprocal.T @ direct, np.eye(3), atol=1e-8)


def test_crystal_direction_returns_unit_vector() -> None:
    phase = make_phase()
    direction = CrystalDirection(coordinates=np.array([1.0, 0.0, 0.0]), phase=phase)
    assert_allclose(np.linalg.norm(direction.unit_vector), 1.0)


def test_crystal_plane_normal_is_unit() -> None:
    phase = make_phase()
    plane = CrystalPlane(miller=MillerIndex(indices=np.array([0, 0, 1]), phase=phase), phase=phase)
    assert_allclose(np.linalg.norm(plane.normal), 1.0)


def test_lattice_rejects_non_positive_lengths() -> None:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    with pytest.raises(ValueError):
        Lattice(0.0, 3.0, 5.0, 90.0, 90.0, 120.0, crystal_frame=crystal)


def test_miller_index_rejects_zero_triplet() -> None:
    phase = make_phase()
    with pytest.raises(ValueError):
        MillerIndex(indices=np.array([0, 0, 0]), phase=phase)


def test_crystal_direction_rejects_zero_vector() -> None:
    phase = make_phase()
    with pytest.raises(ValueError):
        CrystalDirection(coordinates=np.array([0.0, 0.0, 0.0]), phase=phase)


def test_phase_rejects_mismatched_lattice_frame() -> None:
    phase = make_phase()
    other_crystal = ReferenceFrame(
        name="other_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    with pytest.raises(ValueError):
        Phase(
            name="bad-phase",
            lattice=phase.lattice,
            symmetry=phase.symmetry,
            crystal_frame=other_crystal,
        )


def test_phase_can_be_created_from_cif_string() -> None:
    pytest.importorskip("pymatgen.core")
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    phase = Phase.from_cif_string(NACL_CIF, crystal_frame=crystal, primitive=True)
    assert phase.name == "NaCl"
    assert phase.space_group_number == 225
    assert phase.space_group_symbol == "Fm-3m"
    assert phase.chemical_formula == "NaCl"
    assert phase.symmetry.point_group == "m-3m"
    assert phase.unit_cell is not None
    assert len(phase.unit_cell.sites) == 2
    assert phase.provenance is not None
    assert phase.provenance.source_system == "cif"


def test_phase_can_be_created_from_pymatgen_structure() -> None:
    structure_cls = pytest.importorskip("pymatgen.core").Structure
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    structure = structure_cls.from_str(NACL_CIF, fmt="cif")
    phase = Phase.from_pymatgen_structure(structure, crystal_frame=crystal, phase_name="rocksalt")
    assert phase.name == "rocksalt"
    assert phase.space_group_number == 225
    assert phase.space_group_symbol == "Fm-3m"
    assert phase.unit_cell is not None
    assert phase.unit_cell.lattice == phase.lattice


def test_phase_can_be_created_from_cif_path(tmp_path) -> None:
    pytest.importorskip("pymatgen.core")
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    cif_path = tmp_path / "nacl.cif"
    cif_path.write_text(NACL_CIF, encoding="utf-8")
    phase = Phase.from_cif(cif_path, crystal_frame=crystal, phase_name="nacl-phase")
    assert phase.name == "nacl-phase"
    assert phase.space_group_number == 225
    assert phase.provenance is not None
    assert phase.provenance.source_path == str(cif_path)
    assert phase.unit_cell is not None
    assert len(phase.unit_cell.sites) == 8


def test_unit_cell_can_be_created_from_pymatgen_structure() -> None:
    structure_cls = pytest.importorskip("pymatgen.core").Structure
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    structure = structure_cls.from_str(NACL_CIF, fmt="cif")
    lattice = Lattice.from_pymatgen_lattice(structure.lattice, crystal_frame=crystal)
    unit_cell = Phase.from_pymatgen_structure(structure, crystal_frame=crystal).unit_cell
    assert unit_cell is not None
    assert unit_cell.lattice.a == pytest.approx(lattice.a)
    assert unit_cell.lattice.alpha_deg == pytest.approx(lattice.alpha_deg)
