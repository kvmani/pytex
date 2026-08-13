"""Tests for the Kearns orientation parameter and the pole orientation tensor.

The identities that make these tests meaningful are exact, not empirical: the
Kearns parameters along an orthonormal specimen triad sum to 1 for every
texture, a random texture gives 1/3 in every direction, and a single crystal
gives 1 along its own basal pole and 0 perpendicular to it. Where a literature
number is asserted it is the tabulated value of the defining report, not a
prior output of this code.
"""

from __future__ import annotations

import json
from itertools import pairwise

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation as SciRot

from pytex import (
    ODF,
    FrameDomain,
    Handedness,
    KernelSpec,
    Lattice,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.core.lattice import CrystalPlane
from pytex.core.sphere import raster_solid_angle_weights, spherical_angles_to_directions
from pytex.texture.kearns import (
    KEARNS_ISOTROPIC_VALUE,
    DiffractogramReflection,
    KearnsReport,
    basal_tilt_angle_deg,
    basal_tilt_profile,
    harris_texture_coefficients,
    kearns_from_diffractogram,
    kearns_from_odf,
    kearns_from_orientations,
    kearns_from_pole_figure,
    kearns_from_tilt_profile,
    kernel_axis_shrinkage,
    pole_orientation_tensor,
)

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
SPECIMEN = ReferenceFrame(
    "sample_rd_td_nd", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
)
SYMMETRY = SymmetrySpec.from_point_group("6/mmm", reference_frame=CRYSTAL)
# alpha-Zr; c/a = 1.5926, the ratio Kearns quotes as 1.59 for his tilt table.
LATTICE = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=CRYSTAL)
ZIRCONIUM = Phase("alpha_zr", lattice=LATTICE, symmetry=SYMMETRY, crystal_frame=CRYSTAL)
BASAL = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)

# Kearns (1965) WAPD-TM-472 Table 3, longitudinal section of the swaged
# Zircaloy-2 rod: the azimuthally averaged (0001) pole density at the midpoints
# of ten-degree tilt bins, and the f he tabulates from them.
KEARNS_TABLE3_LS_TILTS_DEG = np.arange(5.0, 90.0, 10.0)
KEARNS_TABLE3_LS_INTENSITY = np.array([3.27, 2.71, 1.69, 1.35, 1.17, 0.97, 0.73, 0.62, 0.55])
KEARNS_TABLE3_LS_F = 0.488

# Kearns (1965) Table 2, "angle of tilt to (0001)" for alpha-Zr at c/a = 1.59.
KEARNS_TABLE2_TILTS_DEG = {
    (0, 0, 0, 2): 0.0,
    (1, 0, -1, 5): 20.2,
    (1, 0, -1, 4): 24.7,
    (1, 0, -1, 3): 31.5,
    (2, 0, -2, 5): 36.3,
    (1, 0, -1, 2): 42.5,
    (2, 0, -2, 3): 50.7,
    (1, 1, -2, 2): 57.8,
    (1, 0, -1, 1): 61.4,
    (2, 1, -3, 2): 67.6,
    (3, 0, -3, 2): 70.1,
    (2, 0, -2, 1): 74.8,
    (2, 1, -3, 1): 78.4,
    (1, 0, -1, 0): 90.0,
    (1, 1, -2, 0): 90.0,
}


def orientations_from_matrices(matrices: np.ndarray) -> OrientationSet:
    return OrientationSet.from_matrices(
        np.asarray(matrices, dtype=np.float64).reshape(-1, 3, 3),
        crystal_frame=CRYSTAL,
        specimen_frame=SPECIMEN,
        phase=ZIRCONIUM,
        symmetry=SYMMETRY,
    )


def rotations_with_c_axes(c_axes: np.ndarray, *, seed: int) -> np.ndarray:
    """Rotation matrices whose crystal ``[0001]`` maps to the given specimen axes.

    The remaining degree of freedom — rotation about the c axis — is randomized,
    which is what makes the result a basal *fibre* rather than a component.
    """

    rng = np.random.default_rng(seed)
    axes = np.asarray(c_axes, dtype=np.float64)
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
    helper = np.tile(np.array([1.0, 0.0, 0.0]), (axes.shape[0], 1))
    degenerate = np.abs(axes @ np.array([1.0, 0.0, 0.0])) > 0.9
    helper[degenerate] = np.array([0.0, 1.0, 0.0])
    first = np.cross(helper, axes)
    first /= np.linalg.norm(first, axis=1, keepdims=True)
    second = np.cross(axes, first)
    spin = rng.uniform(0.0, 2.0 * np.pi, size=axes.shape[0])
    x_axis = first * np.cos(spin)[:, None] + second * np.sin(spin)[:, None]
    y_axis = np.cross(axes, x_axis)
    return np.stack([x_axis, y_axis, axes], axis=2)


def basal_fibre(axis: np.ndarray, *, spread_deg: float, count: int, seed: int) -> OrientationSet:
    """A basal fibre: c axes clustered about ``axis`` with a Gaussian pole density.

    ``spread_deg`` is the Gaussian halfwidth of the *pole density on the
    sphere*, so the tilt angles are drawn from a density proportional to
    ``exp(-phi^2 / 2 sigma^2) sin(phi)``. Drawing the tilt angle from the
    Gaussian directly would instead give a pole density diverging as
    ``1/sin(phi)`` at the fibre axis — a singular texture that no diffraction
    method can be fairly judged against.
    """

    rng = np.random.default_rng(seed)
    axis = np.asarray(axis, dtype=np.float64) / np.linalg.norm(axis)
    grid = np.linspace(0.0, np.pi / 2.0, 4001)
    density = np.exp(-0.5 * (grid / np.deg2rad(spread_deg)) ** 2) * np.sin(grid)
    cumulative = np.cumsum(density)
    cumulative /= cumulative[-1]
    tilt = np.interp(rng.uniform(size=count), cumulative, grid)
    azimuth = rng.uniform(0.0, 2.0 * np.pi, size=count)
    helper = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, axis)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    c_axes = np.cos(tilt)[:, None] * axis + np.sin(tilt)[:, None] * (
        np.cos(azimuth)[:, None] * u + np.sin(azimuth)[:, None] * v
    )
    return orientations_from_matrices(rotations_with_c_axes(c_axes, seed=seed + 1))


def random_orientations(count: int = 6000, *, seed: int = 5) -> OrientationSet:
    return orientations_from_matrices(SciRot.random(count, random_state=seed).as_matrix())


# --------------------------------------------------------------------------- #
# The tensor and the identities that follow from it
# --------------------------------------------------------------------------- #


def test_pole_orientation_tensor_of_an_orthonormal_triad_is_isotropic() -> None:
    tensor = pole_orientation_tensor(np.eye(3))
    assert_allclose(tensor, np.eye(3) / 3.0, atol=1e-15)
    assert_allclose(np.trace(tensor), 1.0, atol=1e-15)


def test_pole_orientation_tensor_always_has_unit_trace() -> None:
    directions = SciRot.random(200, random_state=1).apply(np.array([0.0, 0.0, 1.0]))
    weights = np.abs(np.random.default_rng(2).normal(size=200))
    assert_allclose(np.trace(pole_orientation_tensor(directions, weights)), 1.0, atol=1e-12)


def test_pole_orientation_tensor_is_read_only_and_symmetric() -> None:
    tensor = pole_orientation_tensor(SciRot.random(50, random_state=4).as_matrix()[:, :, 2])
    assert not tensor.flags.writeable
    assert_allclose(tensor, tensor.T, atol=1e-15)


@pytest.mark.parametrize(
    ("weights", "message"),
    [
        ([1.0, -1.0], "finite and non-negative"),
        ([0.0, 0.0], "positive value"),
    ],
)
def test_pole_orientation_tensor_rejects_bad_weights(weights, message) -> None:
    with pytest.raises(ValueError, match=message):
        pole_orientation_tensor(np.eye(3)[:2], weights)


def test_random_texture_gives_one_third_in_every_direction() -> None:
    report = kearns_from_orientations(random_orientations(), pole=BASAL)
    assert_allclose(report.values, KEARNS_ISOTROPIC_VALUE, atol=0.01)
    assert report.triad_sum == pytest.approx(1.0, abs=1e-12)


def test_triad_sum_is_exactly_one_for_a_strong_texture() -> None:
    report = kearns_from_orientations(
        basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=12.0, count=3000, seed=9),
        pole=BASAL,
    )
    assert report.triad_sum == pytest.approx(1.0, abs=1e-12)
    assert report.value("ND") > 0.9


def test_triad_sum_holds_for_a_rotated_triad_too() -> None:
    """The sum rule is a property of the trace, so any orthonormal triad obeys it."""

    triad = SciRot.from_euler("zyz", [37.0, 61.0, 13.0], degrees=True).as_matrix().T
    report = kearns_from_orientations(
        basal_fibre(np.array([1.0, 0.0, 2.0]), spread_deg=20.0, count=2000, seed=17),
        pole=BASAL,
        directions=triad,
        direction_labels=("u", "v", "w"),
    )
    assert report.is_orthonormal_triad
    assert report.triad_sum == pytest.approx(1.0, abs=1e-12)


def test_single_crystal_resolves_to_one_along_its_own_basal_pole() -> None:
    report = kearns_from_orientations(orientations_from_matrices(np.eye(3)), pole=BASAL)
    assert_allclose(report.values, [0.0, 0.0, 1.0], atol=1e-12)
    assert report.value("ND") == pytest.approx(1.0)


def test_an_ideal_basal_girdle_gives_one_half_in_the_girdle_plane() -> None:
    """c axes spread uniformly in the RD-TD plane give f = 1/2, 1/2, 0."""

    azimuth = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)
    c_axes = np.stack([np.cos(azimuth), np.sin(azimuth), np.zeros_like(azimuth)], axis=1)
    report = kearns_from_orientations(
        orientations_from_matrices(rotations_with_c_axes(c_axes, seed=3)), pole=BASAL
    )
    assert_allclose(report.values, [0.5, 0.5, 0.0], atol=1e-10)


def test_off_axis_direction_follows_the_quadratic_form() -> None:
    report = kearns_from_orientations(
        basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=25.0, count=2000, seed=21),
        pole=BASAL,
        directions=[[1.0, 0.0, 1.0]],
        direction_labels=("diagonal",),
    )
    triad = kearns_from_orientations(
        basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=25.0, count=2000, seed=21), pole=BASAL
    )
    tensor = np.asarray(triad.orientation_tensor)
    unit = np.array([1.0, 0.0, 1.0]) / np.sqrt(2.0)
    assert report.value("diagonal") == pytest.approx(float(unit @ tensor @ unit), abs=1e-12)


def test_grain_weights_change_the_answer() -> None:
    matrices = np.stack([np.eye(3), SciRot.from_euler("x", 90, degrees=True).as_matrix()])
    orientations = orientations_from_matrices(matrices)
    equal = kearns_from_orientations(orientations, pole=BASAL)
    weighted = kearns_from_orientations(orientations, pole=BASAL, weights=[3.0, 1.0])
    assert equal.value("ND") == pytest.approx(0.5)
    assert weighted.value("ND") == pytest.approx(0.75)


def test_orientations_route_rejects_a_foreign_crystal_frame() -> None:
    other = ReferenceFrame("other", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=other)
    phase = Phase(
        "other_zr",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=other),
        crystal_frame=other,
    )
    pole = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=phase)
    with pytest.raises(ValueError, match="crystal frame"):
        kearns_from_orientations(random_orientations(count=10), pole=pole)


# --------------------------------------------------------------------------- #
# The ODF route and the kernel shrinkage law
# --------------------------------------------------------------------------- #


def test_kernel_shrinkage_tends_to_one_for_a_narrow_kernel() -> None:
    assert kernel_axis_shrinkage(KernelSpec(halfwidth_deg=2.0)) == pytest.approx(1.0, abs=3e-3)


def test_kernel_shrinkage_falls_monotonically_with_halfwidth() -> None:
    values = [kernel_axis_shrinkage(KernelSpec(halfwidth_deg=h)) for h in (2, 5, 10, 20, 40)]
    assert all(later < earlier for earlier, later in pairwise(values))


def test_kernel_shrinkage_matches_direct_integration_over_so3() -> None:
    """The closed form against a Monte-Carlo integral of the kernel over SO(3).

    The derivation is analytic — Rodrigues plus the moments of ``(a . c)^2`` for
    a uniform axis — so this checks the algebra rather than tuning a constant.
    """

    sample = SciRot.random(120_000, random_state=13)
    angles = sample.magnitude()
    axes = sample.apply(np.array([0.0, 0.0, 1.0]))
    for halfwidth in (10.0, 25.0):
        kernel = KernelSpec(halfwidth_deg=halfwidth)
        weights = np.asarray(kernel.evaluate(angles), dtype=np.float64)
        monte_carlo = float(np.einsum("n,n,n->", weights, axes[:, 2], axes[:, 2]) / weights.sum())
        assert kernel_axis_shrinkage(kernel) == pytest.approx(monte_carlo, abs=5e-3)


def test_odf_support_reading_equals_the_orientation_route() -> None:
    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=15.0, count=1500, seed=31)
    odf = ODF.from_orientations(orientations, kernel=KernelSpec(halfwidth_deg=8.0))
    from_odf = kearns_from_odf(odf, pole=BASAL, deconvolve_kernel=True)
    from_orientations = kearns_from_orientations(orientations, pole=BASAL)
    assert_allclose(from_odf.values, from_orientations.values, atol=1e-12)


def test_odf_density_reading_is_the_support_reading_shrunk_toward_one_third() -> None:
    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=15.0, count=1500, seed=31)
    kernel = KernelSpec(halfwidth_deg=15.0)
    odf = ODF.from_orientations(orientations, kernel=kernel)
    support = kearns_from_odf(odf, pole=BASAL, deconvolve_kernel=True)
    density = kearns_from_odf(odf, pole=BASAL, deconvolve_kernel=False)
    shrinkage = (3.0 * kernel_axis_shrinkage(kernel) - 1.0) / 2.0
    expected = KEARNS_ISOTROPIC_VALUE + shrinkage * (
        np.asarray(support.values) - KEARNS_ISOTROPIC_VALUE
    )
    assert_allclose(density.values, expected, atol=1e-12)
    assert 0.0 < shrinkage < 1.0
    assert density.triad_sum == pytest.approx(1.0, abs=1e-12)


def test_odf_kernel_smoothing_costs_a_reported_amount() -> None:
    """A single crystal read as a density is short of 1 by exactly the kernel."""

    odf = ODF.from_orientations(
        orientations_from_matrices(np.eye(3)), kernel=KernelSpec(halfwidth_deg=10.0)
    )
    density = kearns_from_odf(odf, pole=BASAL)
    shrinkage = float(density.diagnostics["kernel_shrinkage_factor"])
    assert density.value("ND") == pytest.approx(
        KEARNS_ISOTROPIC_VALUE + shrinkage * (1.0 - KEARNS_ISOTROPIC_VALUE), abs=1e-12
    )
    assert density.value("ND") < 1.0
    assert kearns_from_odf(odf, pole=BASAL, deconvolve_kernel=True).value("ND") == pytest.approx(
        1.0
    )


# --------------------------------------------------------------------------- #
# The pole-figure route
# --------------------------------------------------------------------------- #


def test_pole_figure_from_orientations_agrees_with_the_orientation_route() -> None:
    orientations = basal_fibre(np.array([0.2, 0.0, 1.0]), spread_deg=18.0, count=2000, seed=41)
    figure = PoleFigure.from_orientations(orientations, BASAL)
    assert_allclose(
        kearns_from_pole_figure(figure).values,
        kearns_from_orientations(orientations, pole=BASAL).values,
        atol=1e-12,
    )


def raster_pole_figure(
    density, *, tilt_step_deg: float = 5.0, max_tilt_deg: float = 90.0
) -> tuple[PoleFigure, np.ndarray]:
    """A tilt/rotation raster carrying an analytic pole density, as measured data."""

    polar = np.arange(0.0, max_tilt_deg + 1e-9, tilt_step_deg)
    azimuth = np.arange(0.0, 360.0, tilt_step_deg)
    polar_grid, azimuth_grid = np.meshgrid(polar, azimuth, indexing="ij")
    polar_flat = polar_grid.ravel()
    directions = spherical_angles_to_directions(polar_flat, azimuth_grid.ravel())
    figure = PoleFigure(
        pole=BASAL,
        sample_directions=directions,
        intensities=density(directions),
        specimen_frame=SPECIMEN,
        sampling="sampled_density",
    )
    return figure, polar_flat


def test_uniform_raster_pole_figure_integrates_to_one_third() -> None:
    figure, _ = raster_pole_figure(lambda v: np.ones(v.shape[0]))
    report = kearns_from_pole_figure(figure)
    assert_allclose(report.values, KEARNS_ISOTROPIC_VALUE, atol=2e-3)
    assert report.triad_sum == pytest.approx(1.0, abs=1e-12)
    assert report.diagnostics["measured_solid_angle_fraction"] == pytest.approx(1.0, abs=1e-9)


def test_ignoring_solid_angle_weights_biases_a_raster_by_fifty_percent() -> None:
    """The bias the raster weights exist to remove, measured on a flat figure."""

    figure, polar = raster_pole_figure(lambda v: np.ones(v.shape[0]))
    unweighted = kearns_from_pole_figure(
        figure, integration_weights=np.ones(polar.size) / polar.size
    )
    assert unweighted.value("ND") == pytest.approx(0.5, abs=5e-3)
    assert kearns_from_pole_figure(figure).value("ND") == pytest.approx(
        KEARNS_ISOTROPIC_VALUE, abs=2e-3
    )


def test_explicit_raster_weights_reproduce_the_default() -> None:
    """The default is the hemisphere-capped raster quadrature, stated explicitly."""

    figure, polar = raster_pole_figure(lambda v: 1.0 + 3.0 * v[:, 2] ** 2)
    assert_allclose(
        kearns_from_pole_figure(
            figure,
            integration_weights=raster_solid_angle_weights(polar, polar_max_deg=90.0),
        ).values,
        kearns_from_pole_figure(figure).values,
        atol=1e-12,
    )


def test_capping_the_outermost_band_at_the_equator_matters() -> None:
    """Without the cap the equatorial ring claims twice the solid angle it owns.

    The analytic target is the spherical mean of ``cos^2``, exactly 1/3. The
    uncapped weighting is out by 4 percent on a 5 degree raster; the capped one
    by 0.06 percent. The difference is a quadrature defect, not sampling noise:
    the figure is exactly uniform.
    """

    figure, polar = raster_pole_figure(lambda v: np.ones(v.shape[0]))
    uncapped = kearns_from_pole_figure(
        figure, integration_weights=raster_solid_angle_weights(polar)
    ).value("ND")
    capped = kearns_from_pole_figure(figure).value("ND")
    assert uncapped == pytest.approx(0.3196, abs=1e-3)
    assert capped == pytest.approx(KEARNS_ISOTROPIC_VALUE, abs=1e-3)
    assert abs(capped - KEARNS_ISOTROPIC_VALUE) < abs(uncapped - KEARNS_ISOTROPIC_VALUE) / 10.0


def test_truncating_the_tilt_range_biases_a_basal_texture_upward() -> None:
    """Stopping at 75 degrees discards the poles that pull f down."""

    def density(v: np.ndarray) -> np.ndarray:
        return np.exp(6.0 * (v[:, 2] ** 2 - 1.0)) + 0.35

    complete, _ = raster_pole_figure(density, max_tilt_deg=90.0)
    truncated, _ = raster_pole_figure(density, max_tilt_deg=75.0)
    full_value = kearns_from_pole_figure(complete).value("ND")
    partial = kearns_from_pole_figure(truncated)
    assert partial.value("ND") > full_value
    assert partial.diagnostics["max_polar_deg"] == pytest.approx(75.0)
    assert partial.diagnostics["measured_solid_angle_fraction"] == pytest.approx(
        1.0 - np.cos(np.deg2rad(75.0)), abs=1e-9
    )
    assert any("pseudo-norm" in note for note in partial.notes)


def test_pole_figure_route_refuses_a_figure_with_no_intensity() -> None:
    figure, _ = raster_pole_figure(lambda v: np.zeros(v.shape[0]))
    with pytest.raises(ValueError, match="no positive weight"):
        kearns_from_pole_figure(figure)


def test_pole_figure_route_rejects_mismatched_explicit_weights() -> None:
    figure, _ = raster_pole_figure(lambda v: np.ones(v.shape[0]))
    with pytest.raises(ValueError):
        kearns_from_pole_figure(figure, integration_weights=np.ones(3))


# --------------------------------------------------------------------------- #
# The tilt profile and Kearns' own tabulated calculation
# --------------------------------------------------------------------------- #


def test_kearns_1965_table_3_longitudinal_section_is_reproduced() -> None:
    """Kearns (1965) Table 3: his tabulated intensities must give his tabulated f."""

    report = kearns_from_tilt_profile(
        KEARNS_TABLE3_LS_TILTS_DEG,
        KEARNS_TABLE3_LS_INTENSITY,
        pole=BASAL,
        specimen_frame=SPECIMEN,
    )
    assert report.value("ND") == pytest.approx(KEARNS_TABLE3_LS_F, abs=1e-3)
    assert report.orientation_tensor is None
    assert report.triad_sum is None


def test_a_flat_tilt_profile_is_the_random_value() -> None:
    report = kearns_from_tilt_profile(
        np.arange(2.5, 90.0, 5.0),
        np.ones(18),
        pole=BASAL,
        specimen_frame=SPECIMEN,
    )
    assert report.value("ND") == pytest.approx(KEARNS_ISOTROPIC_VALUE, abs=2e-3)


def test_tilt_profile_ignores_the_scale_of_the_intensities() -> None:
    base = kearns_from_tilt_profile(
        KEARNS_TABLE3_LS_TILTS_DEG, KEARNS_TABLE3_LS_INTENSITY, pole=BASAL, specimen_frame=SPECIMEN
    )
    scaled = kearns_from_tilt_profile(
        KEARNS_TABLE3_LS_TILTS_DEG,
        KEARNS_TABLE3_LS_INTENSITY * 137.0,
        pole=BASAL,
        specimen_frame=SPECIMEN,
    )
    assert scaled.value("ND") == pytest.approx(base.value("ND"), abs=1e-14)


@pytest.mark.parametrize(
    ("tilts", "intensity", "message"),
    [
        ([5.0, 15.0, 40.0], [1.0, 1.0, 1.0], "equally spaced"),
        ([5.0, 95.0], [1.0, 1.0], r"\[0, 90\]"),
        ([5.0, 15.0], [1.0, -1.0], "non-negative"),
        ([5.0, 15.0], [0.0, 0.0], "everywhere zero"),
        ([5.0], [1.0], "at least two"),
    ],
)
def test_tilt_profile_rejects_unusable_input(tilts, intensity, message) -> None:
    with pytest.raises(ValueError, match=message):
        kearns_from_tilt_profile(tilts, intensity, pole=BASAL, specimen_frame=SPECIMEN)


# --------------------------------------------------------------------------- #
# The diffractogram route
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("indices", "expected"), sorted(KEARNS_TABLE2_TILTS_DEG.items()))
def test_basal_tilt_angles_match_kearns_table_2(indices, expected) -> None:
    """Computed from the phase metric, against Kearns' hand-tabulated angles."""

    plane = CrystalPlane.from_miller_bravais(indices, phase=ZIRCONIUM)
    assert basal_tilt_angle_deg(plane) == pytest.approx(expected, abs=0.2)


def test_harris_texture_coefficients_average_to_one() -> None:
    coefficients = harris_texture_coefficients([0.2, 1.4, 3.1, 0.9])
    assert float(np.mean(coefficients)) == pytest.approx(1.0, abs=1e-14)
    assert not coefficients.flags.writeable


def test_harris_normalization_removes_an_unknown_scale() -> None:
    values = np.array([0.2, 1.4, 3.1, 0.9])
    assert_allclose(
        harris_texture_coefficients(values * 42.0),
        harris_texture_coefficients(values),
        atol=1e-14,
    )


def synthetic_reflections(
    orientations: OrientationSet, indices, *, axis: ArrayLike = (0.0, 0.0, 1.0)
) -> list[DiffractogramReflection]:
    """Reflection intensities a symmetric scan would record from this texture.

    The basal-pole density is binned in tilt from ``axis`` and divided by the
    random expectation of each bin, giving times-random units; each reflection
    then reads that density at its own basal tilt. This is the physical content
    of Kearns' assignment, so a route that recovers the texture's true ``f``
    from these numbers has been validated end to end.
    """

    poles = orientations.map_crystal_directions(BASAL.normal)
    vectors = np.asarray(getattr(poles, "values", poles), dtype=np.float64)
    cosine = np.abs(vectors @ np.asarray(axis, dtype=np.float64))
    edges = np.linspace(0.0, 90.0, 91)
    counts, _ = np.histogram(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))), bins=edges)
    band = np.cos(np.deg2rad(edges[:-1])) - np.cos(np.deg2rad(edges[1:]))
    density = (counts / counts.sum()) / band
    midpoints = 0.5 * (edges[:-1] + edges[1:])
    reflections = []
    for quadruple in indices:
        plane = CrystalPlane.from_miller_bravais(quadruple, phase=ZIRCONIUM)
        tilt = basal_tilt_angle_deg(plane)
        reflections.append(
            DiffractogramReflection(
                plane=plane,
                intensity=float(np.interp(tilt, midpoints, density)),
                random_intensity=1.0,
            )
        )
    return reflections


ALPHA_ZR_REFLECTIONS = (
    (0, 0, 0, 2),
    (1, 0, -1, 5),
    (1, 0, -1, 4),
    (1, 0, -1, 3),
    (2, 0, -2, 5),
    (1, 0, -1, 2),
    (2, 0, -2, 3),
    (1, 1, -2, 2),
    (1, 0, -1, 1),
    (2, 1, -3, 2),
    (2, 0, -2, 1),
    (2, 1, -3, 1),
    (1, 0, -1, 0),
)


def test_diffractogram_route_recovers_a_simulated_texture() -> None:
    """End to end: simulate the peak intensities, then recover f from them."""

    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=25.0, count=40_000, seed=61)
    truth = kearns_from_orientations(orientations, pole=BASAL).value("ND")
    recovered = kearns_from_diffractogram(
        synthetic_reflections(orientations, ALPHA_ZR_REFLECTIONS),
        specimen_frame=SPECIMEN,
    )
    assert recovered.value("ND") == pytest.approx(truth, abs=0.02)
    assert recovered.diagnostics["distinct_basal_tilt_count"] >= 10.0
    assert recovered.diagnostics["diffraction_vector_tilt_spread_deg"] == 0.0


def test_diffractogram_route_recovers_a_transverse_texture() -> None:
    """The hard case: almost no basal intensity near the section normal."""

    orientations = basal_fibre(np.array([0.0, 1.0, 0.0]), spread_deg=20.0, count=40_000, seed=71)
    truth = kearns_from_orientations(orientations, pole=BASAL).value("ND")
    recovered = kearns_from_diffractogram(
        synthetic_reflections(orientations, ALPHA_ZR_REFLECTIONS),
        specimen_frame=SPECIMEN,
    )
    assert truth < 0.15
    assert recovered.value("ND") == pytest.approx(truth, abs=0.05)


def test_diffractogram_route_flags_a_fixed_omega_scan() -> None:
    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=25.0, count=5000, seed=81)
    reflections = synthetic_reflections(orientations, ALPHA_ZR_REFLECTIONS)
    tilted = [
        DiffractogramReflection(
            plane=reflection.plane,
            intensity=reflection.intensity,
            random_intensity=reflection.random_intensity,
            specimen_tilt_deg=float(index) * 3.0,
        )
        for index, reflection in enumerate(reflections)
    ]
    report = kearns_from_diffractogram(tilted, specimen_frame=SPECIMEN)
    assert report.diagnostics["diffraction_vector_tilt_spread_deg"] > 2.0
    assert any("fixed-omega" in note for note in report.notes)


def test_diffractogram_route_always_requires_reference_intensities() -> None:
    """Raw peak intensities are not pole densities, whichever normalization is asked for."""

    plane_a = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)
    plane_b = CrystalPlane.from_miller_bravais((1, 0, -1, 0), phase=ZIRCONIUM)
    reflections = [
        DiffractogramReflection(plane=plane_a, intensity=3.0),
        DiffractogramReflection(plane=plane_b, intensity=1.0),
    ]
    for normalization in ("random_standard", "harris"):
        with pytest.raises(ValueError, match="random_intensity"):
            kearns_from_diffractogram(
                reflections, specimen_frame=SPECIMEN, normalization=normalization
            )


def test_the_two_normalizations_give_the_same_kearns_parameter() -> None:
    """f is a ratio, so rescaling every density by a common factor cannot move it.

    The Harris texture coefficient is exactly such a rescaling. It changes the
    absolute pole densities the profile reports -- Kearns measured that error at
    about 23 percent -- and changes f by nothing at all. Anyone expecting the
    normalization to repair a triad that misses 1 is expecting the wrong thing.
    """

    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=20.0, count=5000, seed=83)
    reflections = synthetic_reflections(orientations, ALPHA_ZR_REFLECTIONS)
    standard = kearns_from_diffractogram(reflections, specimen_frame=SPECIMEN)
    harris = kearns_from_diffractogram(reflections, specimen_frame=SPECIMEN, normalization="harris")
    assert harris.value("ND") == pytest.approx(standard.value("ND"), abs=1e-14)
    assert any("invariant under any common scale" in note for note in harris.notes)


def test_diffractogram_route_rejects_an_unknown_normalization() -> None:
    orientations = basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=20.0, count=500, seed=85)
    with pytest.raises(ValueError, match="normalization must be"):
        kearns_from_diffractogram(
            synthetic_reflections(orientations, ALPHA_ZR_REFLECTIONS),
            specimen_frame=SPECIMEN,
            normalization="none",
        )


def test_basal_tilt_profile_holds_the_edge_density_outside_the_measured_range() -> None:
    plane_a = CrystalPlane.from_miller_bravais((1, 0, -1, 1), phase=ZIRCONIUM)  # 61.5 deg
    plane_b = CrystalPlane.from_miller_bravais((1, 0, -1, 2), phase=ZIRCONIUM)  # 42.6 deg
    tilts, profile, notes = basal_tilt_profile(
        [
            DiffractogramReflection(plane=plane_a, intensity=2.0, random_intensity=1.0),
            DiffractogramReflection(plane=plane_b, intensity=4.0, random_intensity=1.0),
        ]
    )
    assert_allclose(tilts, np.arange(5.0, 90.0, 10.0))
    assert profile[0] == pytest.approx(4.0)  # below the measured range
    assert profile[-1] == pytest.approx(2.0)  # above it
    assert any("held constant" in note for note in notes)


def test_basal_tilt_profile_needs_two_distinct_tilts() -> None:
    plane = CrystalPlane.from_miller_bravais((1, 0, -1, 0), phase=ZIRCONIUM)
    other = CrystalPlane.from_miller_bravais((1, 1, -2, 0), phase=ZIRCONIUM)
    with pytest.raises(ValueError, match="same basal tilt"):
        basal_tilt_profile(
            [
                DiffractogramReflection(plane=plane, intensity=1.0, random_intensity=1.0),
                DiffractogramReflection(plane=other, intensity=2.0, random_intensity=1.0),
            ]
        )


def test_basal_tilt_profile_rejects_a_bin_width_that_does_not_divide_ninety() -> None:
    with pytest.raises(ValueError, match="divide 90"):
        basal_tilt_profile(
            synthetic_reflections(random_orientations(count=200), ALPHA_ZR_REFLECTIONS),
            bin_width_deg=7.0,
        )


@pytest.mark.parametrize(
    ("intensity", "random_intensity", "tilt"),
    [(-1.0, None, 0.0), (1.0, 0.0, 0.0), (1.0, None, np.nan)],
)
def test_diffractogram_reflection_validates_its_fields(intensity, random_intensity, tilt) -> None:
    plane = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)
    with pytest.raises(ValueError):
        DiffractogramReflection(
            plane=plane,
            intensity=intensity,
            random_intensity=random_intensity,
            specimen_tilt_deg=tilt,
        )


# --------------------------------------------------------------------------- #
# The report surface
# --------------------------------------------------------------------------- #


def test_describe_states_the_route_the_values_and_the_sum_rule() -> None:
    report = kearns_from_orientations(
        basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=15.0, count=1000, seed=91), pole=BASAL
    )
    prose = report.describe()
    assert "(0002)" in prose
    assert "f_ND" in prose and "f_RD" in prose
    assert "0.3333" in prose
    assert "Triad sum" in prose
    assert "Principal Kearns values" in prose


def test_json_contract_and_describe_agree_on_every_number() -> None:
    report = kearns_from_orientations(
        basal_fibre(np.array([0.0, 0.0, 1.0]), spread_deg=15.0, count=1000, seed=93), pole=BASAL
    )
    payload = report.to_json()
    assert json.loads(report.to_json_string()) == payload
    assert payload["method"] == "orientations"
    assert payload["pole"] == "(0002)"
    assert payload["triad_sum"] == pytest.approx(1.0)
    assert [entry["label"] for entry in payload["directions"]] == ["RD", "TD", "ND"]
    for entry in payload["directions"]:
        assert f"f_{entry['label']} = {entry['f']:.4f}" in report.describe()
    assert len(payload["orientation_tensor"]) == 3
    assert len(payload["principal_values"]) == 3


def test_report_value_lookup_by_label_and_its_error() -> None:
    report = kearns_from_orientations(random_orientations(count=500), pole=BASAL)
    assert report.value("ND") == pytest.approx(float(report.values[2]))
    with pytest.raises(KeyError, match="No direction labelled"):
        report.value("LD")


def test_report_rejects_values_outside_the_unit_interval() -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        KearnsReport(
            values=np.array([1.4]),
            directions=np.array([[0.0, 0.0, 1.0]]),
            direction_labels=("ND",),
            method="orientations",
            pole=BASAL,
            specimen_frame=SPECIMEN,
        )


def test_report_rejects_an_orientation_tensor_without_unit_trace() -> None:
    with pytest.raises(ValueError, match="unit trace"):
        KearnsReport(
            values=np.array([0.5]),
            directions=np.array([[0.0, 0.0, 1.0]]),
            direction_labels=("ND",),
            method="orientations",
            pole=BASAL,
            specimen_frame=SPECIMEN,
            orientation_tensor=np.eye(3),
        )


def test_report_rejects_a_crystal_frame_as_the_specimen_frame() -> None:
    with pytest.raises(ValueError, match="specimen domain"):
        KearnsReport(
            values=np.array([0.5]),
            directions=np.array([[0.0, 0.0, 1.0]]),
            direction_labels=("c",),
            method="orientations",
            pole=BASAL,
            specimen_frame=CRYSTAL,
        )


def test_direction_labels_must_match_the_directions() -> None:
    with pytest.raises(ValueError, match="one label per direction"):
        kearns_from_orientations(
            random_orientations(count=50),
            pole=BASAL,
            directions=[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            direction_labels=("only-one",),
        )


def test_labels_without_directions_are_rejected() -> None:
    with pytest.raises(ValueError, match="cannot be given without directions"):
        kearns_from_orientations(
            random_orientations(count=50), pole=BASAL, direction_labels=("a", "b", "c")
        )
