from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    FrameDomain,
    Handedness,
    OrientationSet,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.texture import (
    STANDARD_BCC_ROLLING_COMPONENTS,
    STANDARD_FCC_ROLLING_COMPONENTS,
    DeLaValleePoussinKernel,
    Fibre,
    TextureComponent,
    component_volume_fractions,
)


def make_frames() -> tuple[ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    return crystal, specimen


def cubic(crystal: ReferenceFrame) -> SymmetrySpec:
    return SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)


def make_set(euler_deg: np.ndarray, *, symmetric: bool = True) -> OrientationSet:
    crystal, specimen = make_frames()
    return OrientationSet.from_euler_angles(
        euler_deg,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=cubic(crystal) if symmetric else None,
        degrees=True,
    )


def test_standard_component_registries_are_well_formed() -> None:
    names = [component.name for component in STANDARD_FCC_ROLLING_COMPONENTS]
    assert names == ["cube", "goss", "brass", "copper", "s"]
    assert all(component.miller_label for component in STANDARD_FCC_ROLLING_COMPONENTS)
    bcc_names = [component.name for component in STANDARD_BCC_ROLLING_COMPONENTS]
    assert "rotated_cube" in bcc_names


def test_component_orientation_construction_and_validation() -> None:
    crystal, specimen = make_frames()
    component = TextureComponent("cube", (0.0, 0.0, 0.0))
    orientation = component.orientation(
        specimen_frame=specimen,
        crystal_frame=crystal,
        symmetry=cubic(crystal),
    )
    assert_allclose(orientation.as_matrix(), np.eye(3), atol=1e-12)
    with pytest.raises(ValueError, match="non-empty"):
        TextureComponent("", (0.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="crystal_frame"):
        component.orientation(specimen_frame=specimen)


def test_volume_fractions_classify_a_cube_and_goss_mixture() -> None:
    euler = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 45.0, 0.0],
            [1.0, 46.0, 0.0],
            [0.0, 45.0, 90.0],
        ]
    )
    orientations = make_set(euler)
    fractions = component_volume_fractions(orientations, tolerance_deg=10.0)
    # 2 cube-like points and 3 goss-like points (the last is a symmetry
    # equivalent of goss, so symmetry-aware assignment must catch it).
    assert fractions["cube"] == pytest.approx(0.4)
    assert fractions["goss"] == pytest.approx(0.6)
    assert fractions["copper"] == pytest.approx(0.0)


def test_volume_fractions_respect_weights_and_validate_inputs() -> None:
    euler = np.array([[0.0, 0.0, 0.0], [0.0, 45.0, 0.0]])
    orientations = make_set(euler)
    fractions = component_volume_fractions(
        orientations,
        tolerance_deg=10.0,
        weights=[3.0, 1.0],
    )
    assert fractions["cube"] == pytest.approx(0.75)
    assert fractions["goss"] == pytest.approx(0.25)
    with pytest.raises(ValueError, match="tolerance_deg"):
        component_volume_fractions(orientations, tolerance_deg=0.0)
    with pytest.raises(ValueError, match="one value per orientation"):
        component_volume_fractions(orientations, weights=[1.0])


def test_fibre_sampling_lies_on_the_fibre() -> None:
    crystal, specimen = make_frames()
    fibre = Fibre.gamma_bcc()
    orientations = fibre.orientations(
        24,
        specimen_frame=specimen,
        crystal_frame=crystal,
        symmetry=cubic(crystal),
    )
    assert len(orientations) == 24
    angles = fibre.angles_to_deg(orientations)
    assert_allclose(angles, 0.0, atol=1e-8)
    assert fibre.volume_fraction(orientations, tolerance_deg=5.0) == pytest.approx(1.0)


def test_fibre_distance_measures_axis_tilt() -> None:
    crystal, specimen = make_frames()
    theta = Fibre.theta()  # <100> || ND
    # Cube orientation lies on the theta fibre; a 10-degree Phi tilt moves
    # the mapped <100> axis 10 degrees off ND.
    euler = np.array([[0.0, 0.0, 0.0], [0.0, 10.0, 0.0]])
    orientations = OrientationSet.from_euler_angles(
        euler,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=cubic(crystal),
        degrees=True,
    )
    angles = theta.angles_to_deg(orientations)
    assert angles[0] == pytest.approx(0.0, abs=1e-8)
    assert angles[1] == pytest.approx(10.0, abs=1e-6)


def test_fibre_symmetry_awareness_without_symmetry_uses_antipodal_only() -> None:
    crystal, specimen = make_frames()
    fibre = Fibre("custom", (0.0, 0.0, 1.0), "ND")
    euler = np.array([[0.0, 180.0, 0.0]])  # crystal +z mapped to -ND
    orientations = OrientationSet.from_euler_angles(
        euler,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=None,
        degrees=True,
    )
    angles = fibre.angles_to_deg(orientations)
    assert angles[0] == pytest.approx(0.0, abs=1e-8)


def test_fibre_validation() -> None:
    with pytest.raises(ValueError, match="non-zero"):
        Fibre("bad", (0.0, 0.0, 0.0), "ND")
    with pytest.raises(ValueError, match="Specimen direction"):
        Fibre("bad", (1.0, 0.0, 0.0), "QQ")
    crystal, specimen = make_frames()
    with pytest.raises(ValueError, match="at least one orientation"):
        Fibre.eta().orientations(0, specimen_frame=specimen, crystal_frame=crystal)


def test_de_la_vallee_poussin_kernel_normalization_and_halfwidth() -> None:
    kernel = DeLaValleePoussinKernel(halfwidth_deg=10.0)
    # Half maximum at the halfwidth by construction.
    ratio = float(
        kernel.evaluate_deg(np.array([10.0]))[0] / kernel.evaluate_deg(np.array([0.0]))[0]
    )
    assert ratio == pytest.approx(0.5, abs=1e-12)
    coefficients = kernel.chebyshev_coefficients(64)
    assert coefficients[0] == pytest.approx(1.0, abs=1e-6)
    assert np.all(coefficients >= -1e-9)
    assert coefficients[1] > coefficients[32] > coefficients[64]


def test_kernel_bandwidth_scales_inversely_with_halfwidth() -> None:
    sharp = DeLaValleePoussinKernel(halfwidth_deg=5.0)
    broad = DeLaValleePoussinKernel(halfwidth_deg=25.0)
    assert sharp.bandwidth(max_bandwidth=256) > broad.bandwidth(max_bandwidth=256)
    with pytest.raises(ValueError, match="halfwidth_deg"):
        DeLaValleePoussinKernel(halfwidth_deg=0.0)
    with pytest.raises(ValueError, match="threshold"):
        DeLaValleePoussinKernel(halfwidth_deg=10.0).bandwidth(threshold=2.0)
