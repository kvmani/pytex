from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    ComplianceTensor,
    FrameDomain,
    Handedness,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    StiffnessTensor,
    SymmetrySpec,
    homogenize_elastic,
)
from pytex.core.misorientation_distribution import _haar_uniform_quaternions


def _copper() -> StiffnessTensor:
    # Single-crystal copper elastic constants (GPa).
    return StiffnessTensor.cubic(168.4, 121.4, 75.4)


def _cubic_frames() -> tuple[ReferenceFrame, ReferenceFrame, SymmetrySpec]:
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
    return crystal, specimen, symmetry


def test_voigt_round_trip_preserves_stiffness() -> None:
    stiffness = _copper()
    reconstructed = StiffnessTensor.from_voigt(stiffness.voigt_matrix())
    np.testing.assert_allclose(reconstructed.voigt_matrix(), stiffness.voigt_matrix())


def test_stiffness_compliance_are_mutual_inverses() -> None:
    stiffness = _copper()
    round_trip = stiffness.compliance().stiffness()
    np.testing.assert_allclose(round_trip.voigt_matrix(), stiffness.voigt_matrix(), atol=1e-9)


def test_cubic_youngs_modulus_matches_analytic_formula() -> None:
    stiffness = _copper()
    voigt = stiffness.compliance().voigt_matrix()
    s11, s12, s44 = voigt[0, 0], voigt[0, 1], voigt[3, 3]

    def analytic(direction: tuple[float, float, float]) -> float:
        vector = np.asarray(direction, dtype=np.float64)
        vector = vector / np.linalg.norm(vector)
        anisotropy = (
            vector[0] ** 2 * vector[1] ** 2
            + vector[1] ** 2 * vector[2] ** 2
            + vector[2] ** 2 * vector[0] ** 2
        )
        return 1.0 / (s11 - 2.0 * (s11 - s12 - s44 / 2.0) * anisotropy)

    for direction in [(1, 0, 0), (1, 1, 0), (1, 1, 1)]:
        assert stiffness.youngs_modulus(direction) == pytest.approx(analytic(direction))


def test_isotropic_tensor_is_direction_independent_and_rotation_invariant() -> None:
    stiffness = StiffnessTensor.isotropic(youngs_modulus=200.0, poisson_ratio=0.3)
    for direction in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 3)]:
        assert stiffness.youngs_modulus(direction) == pytest.approx(200.0)
    rotation = np.linalg.qr(np.random.default_rng(0).standard_normal((3, 3)))[0]
    rotated = stiffness.rotate(rotation)
    assert rotated.youngs_modulus((1, 2, 3)) == pytest.approx(200.0)


def test_rotating_stiffness_permutes_cubic_axes_modulus() -> None:
    stiffness = _copper()
    # A 90-degree rotation about z maps [100] to [010]; the modulus along the
    # rotated [100] axis equals the original [010] modulus (both <100>).
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    rotated = stiffness.rotate(rotation)
    assert rotated.youngs_modulus((1, 0, 0)) == pytest.approx(
        stiffness.youngs_modulus((1, 0, 0))
    )


def test_compliance_from_voigt_matches_stiffness_inverse() -> None:
    stiffness = _copper()
    compliance = ComplianceTensor.from_voigt(np.linalg.inv(stiffness.voigt_matrix()))
    np.testing.assert_allclose(
        compliance.voigt_matrix(), stiffness.compliance().voigt_matrix(), atol=1e-12
    )


def test_homogenize_single_orientation_reproduces_crystal() -> None:
    crystal, specimen, symmetry = _cubic_frames()
    stiffness = _copper()
    orientations = OrientationSet.from_quaternions(
        Rotation.identity().quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    for scheme in ("voigt", "reuss", "hill"):
        aggregate = homogenize_elastic(stiffness, orientations, scheme=scheme)
        np.testing.assert_allclose(
            aggregate.voigt_matrix(), stiffness.voigt_matrix(), atol=1e-9
        )


def test_homogenize_random_population_is_nearly_isotropic_with_bounds() -> None:
    crystal, specimen, symmetry = _cubic_frames()
    stiffness = _copper()
    quaternions = _haar_uniform_quaternions(400, np.random.default_rng(0))
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    directions = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 3)]
    voigt = homogenize_elastic(stiffness, orientations, scheme="voigt")
    reuss = homogenize_elastic(stiffness, orientations, scheme="reuss")
    hill = homogenize_elastic(stiffness, orientations, scheme="hill")
    voigt_e = [voigt.youngs_modulus(d) for d in directions]
    # a random aggregate is nearly isotropic: modulus varies little with direction
    assert (max(voigt_e) - min(voigt_e)) / np.mean(voigt_e) < 0.05
    # Voigt (uniform strain) is stiffer than Reuss (uniform stress); Hill lies between
    for direction in directions:
        assert (
            voigt.youngs_modulus(direction)
            > hill.youngs_modulus(direction)
            > reuss.youngs_modulus(direction)
        )


def test_homogenize_rejects_unknown_scheme() -> None:
    crystal, specimen, symmetry = _cubic_frames()
    orientations = OrientationSet.from_quaternions(
        Rotation.identity().quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    with pytest.raises(ValueError, match="scheme"):
        homogenize_elastic(_copper(), orientations, scheme="bogus")


def test_hexagonal_stiffness_is_transversely_isotropic_in_basal_plane() -> None:
    # Titanium-like hexagonal constants (GPa).
    stiffness = StiffnessTensor.hexagonal(162.4, 92.0, 69.0, 180.7, 46.7)
    basal_x = stiffness.youngs_modulus((1, 0, 0))
    basal_y = stiffness.youngs_modulus((0, 1, 0))
    assert basal_x == pytest.approx(basal_y)
    # the c-axis modulus differs from the basal-plane modulus
    assert stiffness.youngs_modulus((0, 0, 1)) != pytest.approx(basal_x)
