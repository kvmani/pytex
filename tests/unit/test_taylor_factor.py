from __future__ import annotations

import math

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    fcc_octahedral_slip,
    taylor_factors,
    uniaxial_strain_tensor,
)
from pytex.core.misorientation_distribution import _haar_uniform_quaternions


def _cubic_phase() -> tuple[Phase, ReferenceFrame, ReferenceFrame]:
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
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="fcc", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return phase, crystal, specimen


def test_uniaxial_strain_tensor_is_deviatoric_with_unit_equivalent() -> None:
    strain = uniaxial_strain_tensor((0, 0, 1))
    assert np.trace(strain) == pytest.approx(0.0, abs=1e-12)
    # von Mises equivalent of the normalized uniaxial strain is 1
    equivalent = math.sqrt(2.0 / 3.0 * np.tensordot(strain, strain))
    assert equivalent == pytest.approx(1.0)


def test_cube_orientation_taylor_factor_is_sqrt_six() -> None:
    phase, crystal, specimen = _cubic_phase()
    cube = Orientation(
        Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    factor = taylor_factors(fcc_octahedral_slip(phase), cube, tension_axis=(0, 0, 1))
    # the exact full-constraint Taylor factor of the cube orientation is sqrt(6)
    assert factor == pytest.approx(math.sqrt(6.0), abs=1e-6)


def test_random_fcc_mean_taylor_factor_matches_classic_value() -> None:
    phase, crystal, specimen = _cubic_phase()
    quaternions = _haar_uniform_quaternions(300, np.random.default_rng(0))
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    factors = taylor_factors(fcc_octahedral_slip(phase), orientations, tension_axis=(0, 0, 1))
    assert factors.shape == (300,)
    # the classic random-fcc average Taylor factor is ~3.06
    assert np.mean(factors) == pytest.approx(3.06, abs=0.06)
    assert np.all(np.isfinite(factors))


def test_taylor_factor_is_invariant_to_tension_axis_for_random_texture() -> None:
    phase, crystal, specimen = _cubic_phase()
    quaternions = _haar_uniform_quaternions(200, np.random.default_rng(1))
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    family = fcc_octahedral_slip(phase)
    along_z = np.mean(taylor_factors(family, orientations, tension_axis=(0, 0, 1)))
    along_x = np.mean(taylor_factors(family, orientations, tension_axis=(1, 0, 0)))
    # a random texture has no preferred axis, so the mean factor is isotropic
    assert along_z == pytest.approx(along_x, abs=0.1)


def test_taylor_factor_rejects_conflicting_and_non_deviatoric_inputs() -> None:
    phase, _, _ = _cubic_phase()
    family = fcc_octahedral_slip(phase)
    crystal = phase.crystal_frame
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    orientation = Orientation(
        Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    with pytest.raises(ValueError, match="only one"):
        taylor_factors(
            family, orientation, strain_tensor=np.eye(3) - np.eye(3), tension_axis=(0, 0, 1)
        )
    with pytest.raises(ValueError, match="deviatoric"):
        taylor_factors(family, orientation, strain_tensor=np.eye(3))


def test_crystal_map_taylor_factor_map_reshapes_to_grid() -> None:
    phase, crystal, specimen = _cubic_phase()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            )
            for _ in range(4)
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    taylor = crystal_map.taylor_factor_map(fcc_octahedral_slip(phase), tension_axis="z")
    assert taylor.shape == (2, 2)
    assert np.allclose(taylor, math.sqrt(6.0))
