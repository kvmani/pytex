"""Pole-figure resampling, m.r.d. normalization and arithmetic.

These surfaces exist as one dependency chain: two pole figures cannot be
combined until they share a support (resampling) and a physical scale
(m.r.d.), so the tests are ordered the same way.
"""

from __future__ import annotations

from dataclasses import replace as dataclass_replace
from itertools import pairwise

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.contracts import from_json_contract, to_json_contract
from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    S2Grid,
    SymmetrySpec,
    raster_solid_angle_weights,
)
from pytex.plotting._render import ScatterLayer2D
from pytex.plotting.builders import build_pole_figure_difference_spec
from pytex.plotting.runtime import plot_pole_figure_difference
from pytex.texture import (
    ODF,
    HarmonicODF,
    KernelSpec,
    PoleFigure,
    PoleFigureResidualReport,
    residual_reports_for_pole_figures,
)
from pytex.texture.models import random_pole_density


def make_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
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
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def make_pole(phase: Phase) -> CrystalPlane:
    return CrystalPlane(MillerIndex((1, 1, 1), phase=phase), phase=phase)


def uniform_grid(specimen: ReferenceFrame, resolution_deg: float = 10.0) -> S2Grid:
    return S2Grid.equispaced(
        resolution_deg,
        reference_frame=specimen,
        hemisphere="upper",
        antipodal=True,
    )


def constant_density_figure(
    specimen: ReferenceFrame,
    phase: Phase,
    *,
    value: float = 1.0,
    resolution_deg: float = 10.0,
) -> PoleFigure:
    """A pole figure that is exactly ``value`` everywhere, sampled on a grid."""

    grid = uniform_grid(specimen, resolution_deg)
    directions = grid.vectors.values
    return PoleFigure(
        pole=make_pole(phase),
        sample_directions=directions,
        intensities=np.full(directions.shape[0], float(value)),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )


# --------------------------------------------------------------------------
# Sampling semantics
# --------------------------------------------------------------------------


def test_from_orientations_declares_the_pole_cloud_reading() -> None:
    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        np.zeros((4, 3)), specimen_frame=specimen, phase=phase
    )
    figure = PoleFigure.from_orientations(orientations, make_pole(phase))
    assert figure.sampling == "scattered_poles"


def test_sampled_density_reading_is_accepted_and_preserved() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    assert figure.sampling == "sampled_density"


def test_unknown_sampling_tag_raises() -> None:
    _, specimen, phase = make_context()
    with pytest.raises(ValueError, match="sampling must be"):
        PoleFigure(
            pole=make_pole(phase),
            sample_directions=np.array([[0.0, 0.0, 1.0]]),
            intensities=np.array([1.0]),
            specimen_frame=specimen,
            sampling="measured",  # type: ignore[arg-type]
        )


def test_sampling_survives_the_json_contract_round_trip() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    restored = from_json_contract(to_json_contract(figure))
    assert isinstance(restored, PoleFigure)
    assert restored.sampling == "sampled_density"


def test_payloads_predating_the_sampling_tag_read_as_pole_clouds() -> None:
    _, specimen, phase = make_context()
    payload = to_json_contract(constant_density_figure(specimen, phase))
    del payload["sampling"]
    restored = from_json_contract(payload)
    assert isinstance(restored, PoleFigure)
    assert restored.sampling == "scattered_poles"


def test_unknown_sampling_tag_in_a_payload_raises_at_the_contract_boundary() -> None:
    _, specimen, phase = make_context()
    payload = to_json_contract(constant_density_figure(specimen, phase))
    payload["sampling"] = "whatever"
    with pytest.raises(ValueError, match="must be 'scattered_poles' or 'sampled_density'"):
        from_json_contract(payload)


def test_correcting_a_figure_does_not_change_what_its_intensities_mean() -> None:
    from pytex.texture import PoleFigureCorrectionSpec

    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, value=2.0)
    corrected = PoleFigureCorrectionSpec(scale=0.5).apply(figure)
    assert corrected.sampling == "sampled_density"
    assert_allclose(corrected.intensities, 1.0)


# --------------------------------------------------------------------------
# Resampling
# --------------------------------------------------------------------------


def test_interpolation_reproduces_a_constant_field_exactly() -> None:
    """Nadaraya-Watson is a weighted mean, so a flat field must survive it.

    This is the partition-of-unity identity and it holds for any halfwidth and
    any target direction, which makes it the sharpest available check that the
    interpolating branch is a mean and not a sum.
    """

    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, value=3.0, resolution_deg=10.0)
    resampled = figure.on_grid(uniform_grid(specimen, 15.0), normalize=False)
    assert_allclose(resampled.intensities, 3.0, rtol=1e-12)


def uniform_pole_cloud(
    specimen: ReferenceFrame, phase: Phase, *, resolution_deg: float
) -> PoleFigure:
    """A pole cloud representing a genuinely random texture.

    Each pole carries its own solid angle as weight, so the cloud represents a
    uniform distribution rather than the ring structure of the grid that
    generated the directions.
    """

    grid = uniform_grid(specimen, resolution_deg)
    return PoleFigure(
        pole=make_pole(phase),
        sample_directions=grid.vectors.values,
        intensities=grid.weights,
        specimen_frame=specimen,
        antipodal=True,
        sampling="scattered_poles",
    )


def test_density_estimation_of_a_random_texture_converges_to_one_mrd() -> None:
    """A random texture must come back at 1 m.r.d. everywhere.

    The analytic normalization — divide by the total pole weight and by the
    kernel's spherical mean — is calibrated by exactly this case, so it is
    checked with ``normalize=False``; the grid-mean rescaling would otherwise
    make the identity true by construction and test nothing.

    Refining the cloud is the real assertion. The residual is the discretization
    error of the finite cloud, not a bias in the normalizing constant, so it
    must fall like the square of the spacing. A wrong constant would leave a
    floor that refinement could not remove.
    """

    _, specimen, phase = make_context()
    target = uniform_grid(specimen, 10.0)
    errors = [
        float(
            np.max(
                np.abs(
                    uniform_pole_cloud(specimen, phase, resolution_deg=resolution)
                    .on_grid(target, halfwidth_deg=10.0, normalize=False)
                    .intensities
                    - 1.0
                )
            )
        )
        for resolution in (8.0, 4.0, 2.0)
    ]
    assert errors[0] < 0.05
    # Second order: halving the spacing must quarter the error.
    for coarse, fine in pairwise(errors):
        assert fine < coarse / 3.0


def test_normalization_makes_the_grid_weighted_mean_exactly_one() -> None:
    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    figure = constant_density_figure(specimen, phase, value=7.5)
    resampled = figure.on_grid(grid, normalize=True)
    assert_allclose(float(np.sum(grid.weights * resampled.intensities)), 1.0, rtol=1e-12)


def test_resampling_declares_the_result_a_sampled_density() -> None:
    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        np.zeros((8, 3)), specimen_frame=specimen, phase=phase
    )
    figure = PoleFigure.from_orientations(orientations, make_pole(phase))
    resampled = figure.on_grid(uniform_grid(specimen))
    assert resampled.sampling == "sampled_density"
    assert resampled.sample_directions.shape[0] == len(uniform_grid(specimen))


def test_a_sharp_texture_resamples_to_a_peak_at_its_pole() -> None:
    """A single orientation must put its density maximum on its own pole."""

    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        np.zeros((1, 3)), specimen_frame=specimen, phase=phase
    )
    figure = PoleFigure.from_orientations(orientations, make_pole(phase))
    grid = uniform_grid(specimen, 5.0)
    resampled = figure.on_grid(grid, halfwidth_deg=8.0)
    peak = resampled.sample_directions[int(np.argmax(resampled.intensities))]
    expected = figure.sample_directions
    # The peak must coincide with one of the mapped {111} poles, up to the
    # antipodal identification the figure declares.
    separations_deg = np.degrees(np.arccos(np.clip(np.abs(expected @ peak), -1.0, 1.0)))
    assert float(np.min(separations_deg)) < 5.0
    assert float(np.max(resampled.intensities)) > 1.0


def test_resampling_rejects_a_grid_in_the_wrong_frame() -> None:
    _, specimen, phase = make_context()
    other = ReferenceFrame(
        name="other-specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    figure = constant_density_figure(specimen, phase)
    foreign = S2Grid.equispaced(15.0, reference_frame=other, hemisphere="upper", antipodal=True)
    with pytest.raises(ValueError, match="specimen frame"):
        figure.on_grid(foreign)


def test_resampling_rejects_an_unknown_estimator() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    with pytest.raises(ValueError, match="estimator must be"):
        figure.on_grid(uniform_grid(specimen), estimator="nearest")  # type: ignore[arg-type]


def test_resampling_rejects_an_impossible_halfwidth() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    with pytest.raises(ValueError, match="halfwidth_deg"):
        figure.on_grid(uniform_grid(specimen), halfwidth_deg=0.0)


def test_normalizing_an_everywhere_zero_field_raises_rather_than_dividing_by_zero() -> None:
    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    figure = PoleFigure(
        pole=make_pole(phase),
        sample_directions=grid.vectors.values,
        intensities=np.zeros(len(grid)),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    with pytest.raises(ValueError, match="zero everywhere"):
        figure.on_grid(grid, normalize=True)


# --------------------------------------------------------------------------
# Multiples of a random distribution
# --------------------------------------------------------------------------


def test_raster_weights_sum_to_one_and_never_vanish() -> None:
    polar = np.repeat(np.arange(0.0, 70.1, 5.0), 72)
    weights = raster_solid_angle_weights(polar)
    assert_allclose(float(np.sum(weights)), 1.0, rtol=1e-12)
    assert float(np.min(weights)) > 0.0


def test_raster_weights_integrate_a_known_function_to_the_analytic_cap_mean() -> None:
    """The weights must reproduce the solid-angle mean, not the naive mean.

    For ``f = cos(polar)`` over a cap of half-angle ``t``, the solid-angle mean
    is ``sin^2(t) / (2 * (1 - cos t))``, which is a strict identity. The raster
    covers rings from 0 to 70 degrees at 2.5 degree steps, so its bands extend
    to 71.25 degrees. Refinement must drive the error to zero; the unweighted
    mean does not converge to this value at all, which is the point.
    """

    step = 2.5
    polar = np.repeat(np.arange(0.0, 70.0 + 1e-9, step), 72)
    weights = raster_solid_angle_weights(polar)
    cap_rad = np.deg2rad(70.0 + step / 2.0)
    exact = float(np.sin(cap_rad) ** 2 / (2.0 * (1.0 - np.cos(cap_rad))))
    weighted = float(np.sum(weights * np.cos(np.deg2rad(polar))))
    assert abs(weighted - exact) < 5e-4
    naive = float(np.mean(np.cos(np.deg2rad(polar))))
    assert abs(naive - exact) > 1e-2


def test_raster_weights_favour_the_equator_over_the_pole() -> None:
    """The whole point of weighting: a raster over-samples near the pole."""

    rings = np.arange(0.0, 90.1, 10.0)
    polar = np.repeat(rings, 36)
    weights = raster_solid_angle_weights(polar).reshape(rings.size, 36)
    per_ring = weights.sum(axis=1)
    assert float(per_ring[-1]) > 5.0 * float(per_ring[0])


def test_normalize_to_mrd_with_supplied_weights_is_exact() -> None:
    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    rng = np.random.default_rng(20260808)
    figure = PoleFigure(
        pole=make_pole(phase),
        sample_directions=grid.vectors.values,
        intensities=rng.uniform(0.5, 4.0, size=len(grid)),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    normalized = figure.normalize_to_mrd(integration_weights=grid.weights)
    assert_allclose(normalized.spherical_mean(integration_weights=grid.weights), 1.0, rtol=1e-12)


def test_normalize_to_mrd_leaves_a_figure_already_in_mrd_alone() -> None:
    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    figure = constant_density_figure(specimen, phase, value=1.0)
    normalized = figure.normalize_to_mrd(integration_weights=grid.weights)
    assert_allclose(normalized.intensities, figure.intensities, rtol=1e-12)


def test_normalize_to_mrd_falls_back_to_an_estimated_spherical_mean() -> None:
    """Without weights the mean is estimated by resampling, not refused."""

    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, value=6.0, resolution_deg=8.0)
    assert_allclose(figure.spherical_mean(), 6.0, rtol=1e-9)
    assert_allclose(figure.normalize_to_mrd().intensities, 1.0, rtol=1e-9)


def test_spherical_mean_rejects_non_positive_weights() -> None:
    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    figure = constant_density_figure(specimen, phase)
    weights = np.array(grid.weights, copy=True)
    weights[0] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        figure.spherical_mean(integration_weights=weights)


def test_mrd_normalization_makes_two_differently_scaled_measurements_comparable() -> None:
    """The reason m.r.d. exists: detector scale must not survive normalization.

    Two measurements of the same texture that differ only by counting time —
    or by a `max` versus `sum` pre-normalization — must become numerically
    identical once both are on the m.r.d. scale.
    """

    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    rng = np.random.default_rng(4242)
    field = rng.uniform(0.5, 4.0, size=len(grid))
    common = {
        "pole": make_pole(phase),
        "sample_directions": grid.vectors.values,
        "specimen_frame": specimen,
        "antipodal": True,
        "sampling": "sampled_density",
    }
    long_count = PoleFigure(intensities=field * 9000.0, **common)  # type: ignore[arg-type]
    scaled_by_max = PoleFigure(intensities=field / field.max(), **common)  # type: ignore[arg-type]
    assert_allclose(
        long_count.normalize_to_mrd(integration_weights=grid.weights).intensities,
        scaled_by_max.normalize_to_mrd(integration_weights=grid.weights).intensities,
        rtol=1e-12,
    )


# --------------------------------------------------------------------------
# Arithmetic
# --------------------------------------------------------------------------


def two_figures_on_one_grid(
    specimen: ReferenceFrame, phase: Phase
) -> tuple[PoleFigure, PoleFigure, S2Grid]:
    """Two different textures resampled onto one shared support."""

    grid = uniform_grid(specimen, 10.0)
    sharp = PoleFigure.from_orientations(
        OrientationSet.from_euler_angles(
            np.zeros((1, 3)), specimen_frame=specimen, phase=phase
        ),
        make_pole(phase),
    ).on_grid(grid, halfwidth_deg=15.0)
    spread = PoleFigure.from_orientations(
        OrientationSet.from_euler_angles(
            np.array([[0.0, 0.0, 0.0], [30.0, 20.0, 10.0], [12.0, 45.0, 33.0]]),
            specimen_frame=specimen,
            phase=phase,
        ),
        make_pole(phase),
    ).on_grid(grid, halfwidth_deg=15.0)
    return sharp, spread, grid


def test_a_figure_minus_itself_is_identically_zero() -> None:
    _, specimen, phase = make_context()
    sharp, _, _ = two_figures_on_one_grid(specimen, phase)
    residual = sharp - sharp
    assert residual.max_absolute_deviation == 0.0
    assert residual.rms_deviation == 0.0


def test_subtraction_is_antisymmetric() -> None:
    _, specimen, phase = make_context()
    sharp, spread, _ = two_figures_on_one_grid(specimen, phase)
    assert_allclose((sharp - spread).values, -(spread - sharp).values, rtol=1e-12)


def test_addition_and_subtraction_are_consistent() -> None:
    """(a + b) - b must recover a, which ties the two operators together."""

    _, specimen, phase = make_context()
    sharp, spread, _ = two_figures_on_one_grid(specimen, phase)
    recovered = (sharp + spread) - spread
    assert_allclose(recovered.values, sharp.intensities, rtol=1e-12)


def test_sum_of_two_mrd_figures_has_mean_two() -> None:
    """Densities add, and so do their means; this fixes the scale convention."""

    _, specimen, phase = make_context()
    sharp, spread, grid = two_figures_on_one_grid(specimen, phase)
    total = sharp + spread
    assert_allclose(total.spherical_mean(integration_weights=grid.weights), 2.0, rtol=1e-9)


def test_subtracting_one_from_an_mrd_figure_gives_deviation_from_random() -> None:
    _, specimen, phase = make_context()
    sharp, _, grid = two_figures_on_one_grid(specimen, phase)
    deviation = sharp - 1.0
    # A normalized figure integrates to 1, so its deviation integrates to 0.
    assert_allclose(float(np.sum(grid.weights * deviation.values)), 0.0, atol=1e-9)
    # And it must be signed: a texture is above random somewhere and below it
    # elsewhere, which is precisely what a non-negative type could not express.
    assert float(np.max(deviation.values)) > 0.0
    assert float(np.min(deviation.values)) < 0.0


def test_ratio_of_a_figure_to_itself_is_one_everywhere() -> None:
    _, specimen, phase = make_context()
    sharp, _, _ = two_figures_on_one_grid(specimen, phase)
    assert_allclose((sharp / sharp).intensities, 1.0, rtol=1e-12)


def test_scaling_is_exactly_undone_by_dividing() -> None:
    _, specimen, phase = make_context()
    sharp, _, _ = two_figures_on_one_grid(specimen, phase)
    assert_allclose(((sharp * 3.5) / 3.5).intensities, sharp.intensities, rtol=1e-12)
    assert_allclose((2.0 * sharp).intensities, (sharp * 2.0).intensities, rtol=1e-12)


def test_scaling_moves_the_mean_proportionally() -> None:
    _, specimen, phase = make_context()
    sharp, _, grid = two_figures_on_one_grid(specimen, phase)
    assert_allclose(
        (sharp * 4.0).spherical_mean(integration_weights=grid.weights), 4.0, rtol=1e-9
    )


def test_difference_carries_labels_into_its_prose() -> None:
    _, specimen, phase = make_context()
    sharp, spread, _ = two_figures_on_one_grid(specimen, phase)
    text = sharp.difference(spread, left_label="measured", right_label="recalculated").describe()
    assert "measured minus recalculated" in text
    assert "m.r.d." in text
    assert "{111}" in text


def test_difference_of_scale_mismatched_figures_says_so() -> None:
    """The describe() must diagnose a normalization error, not merely report one."""

    _, specimen, phase = make_context()
    sharp, _, _ = two_figures_on_one_grid(specimen, phase)
    text = (sharp * 2.0).difference(sharp).describe()
    assert "not on the same normalization" in text


def test_weighted_rms_differs_from_the_unweighted_one() -> None:
    """On a non-equal-area support the unweighted RMS is the wrong summary."""

    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 10.0)
    sharp, spread, _ = two_figures_on_one_grid(specimen, phase)
    residual = sharp - spread
    assert residual.weighted_rms_deviation(grid.weights) != pytest.approx(
        residual.rms_deviation, rel=1e-6
    )


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        ("pole", "different poles"),
        ("frame", "different specimen frames"),
        ("antipodal", "antipodal conventions"),
        ("family", "single-plane"),
        ("support", "different supports"),
    ],
)
def test_combining_mismatched_figures_raises_with_the_reason(mutate: str, message: str) -> None:
    """Every mismatch that could otherwise produce a plausible wrong answer."""

    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    if mutate == "pole":
        other = dataclass_replace(
            figure, pole=CrystalPlane(MillerIndex((2, 0, 0), phase=phase), phase=phase)
        )
    elif mutate == "frame":
        other_frame = ReferenceFrame(
            name="rolled-specimen",
            domain=FrameDomain.SPECIMEN,
            axes=("rd", "td", "nd"),
            handedness=Handedness.RIGHT,
        )
        other = dataclass_replace(figure, specimen_frame=other_frame)
    elif mutate == "antipodal":
        other = dataclass_replace(figure, antipodal=False)
    elif mutate == "family":
        other = dataclass_replace(figure, includes_symmetry_family=False)
    else:
        other = dataclass_replace(
            figure,
            sample_directions=figure.sample_directions[:-1],
            intensities=figure.intensities[:-1],
        )
    with pytest.raises(ValueError, match=message):
        _ = figure - other


def test_dividing_by_a_figure_with_a_hole_refuses_rather_than_inventing_a_value() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    holed = PoleFigure(
        pole=figure.pole,
        sample_directions=figure.sample_directions,
        intensities=np.where(np.arange(len(figure.intensities)) == 3, 0.0, 1.0),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    with pytest.raises(ValueError, match="zero density"):
        _ = figure / holed


def test_adding_a_constant_that_would_go_below_zero_raises() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, value=0.5)
    with pytest.raises(ValueError, match="negative"):
        _ = figure + (-1.0)


def test_scaling_by_a_negative_factor_raises() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    with pytest.raises(ValueError, match="non-negative"):
        _ = figure * -1.0


# --------------------------------------------------------------------------
# Transforms
# --------------------------------------------------------------------------


def test_rotating_a_figure_moves_its_poles_and_keeps_its_densities() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, value=2.0)
    rotated = figure.rotate(Rotation.from_bunge_euler(30.0, 20.0, 10.0))
    assert_allclose(rotated.intensities, figure.intensities, rtol=1e-12)
    assert not np.allclose(rotated.sample_directions, figure.sample_directions)
    # A rotation is an isometry, so every angle to a fixed member is preserved.
    assert_allclose(
        rotated.sample_directions @ rotated.sample_directions[0],
        figure.sample_directions @ figure.sample_directions[0],
        atol=1e-12,
    )


def test_rotating_a_random_texture_leaves_it_random() -> None:
    """The invariance that makes the rotation trustworthy on real data."""

    _, specimen, phase = make_context()
    target = uniform_grid(specimen, 10.0)
    cloud = uniform_pole_cloud(specimen, phase, resolution_deg=4.0)
    before = cloud.on_grid(target, halfwidth_deg=12.0, normalize=False)
    after = cloud.rotate(Rotation.from_bunge_euler(45.0, 35.0, 25.0)).on_grid(
        target, halfwidth_deg=12.0, normalize=False
    )
    assert_allclose(after.intensities, before.intensities, atol=2e-2)


def test_symmetrize_records_the_assumption_and_replicates_the_support() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    orthorhombic = SymmetrySpec.from_point_group("222", reference_frame=specimen)
    symmetrized = figure.symmetrize(orthorhombic)
    assert symmetrized.sample_symmetry == orthorhombic
    expected = len(figure.intensities) * orthorhombic.operators.shape[0]
    assert symmetrized.sample_directions.shape[0] == expected


def test_symmetrizing_makes_the_field_invariant_under_the_group() -> None:
    """The point of symmetrizing: the resampled field becomes group-invariant.

    A single sharp component is deliberately asymmetric; after imposing
    orthorhombic specimen symmetry, rotating the figure by a group operator
    must leave the resampled field unchanged.
    """

    _, specimen, phase = make_context()
    orthorhombic = SymmetrySpec.from_point_group("222", reference_frame=specimen)
    target = uniform_grid(specimen, 10.0)
    figure = PoleFigure.from_orientations(
        OrientationSet.from_euler_angles(
            np.array([[17.0, 29.0, 41.0]]), specimen_frame=specimen, phase=phase
        ),
        make_pole(phase),
    )
    symmetrized = figure.symmetrize(orthorhombic)
    plain = symmetrized.on_grid(target, halfwidth_deg=12.0)
    two_fold_z = Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi)
    turned = symmetrized.rotate(two_fold_z).on_grid(target, halfwidth_deg=12.0)
    assert_allclose(turned.intensities, plain.intensities, atol=1e-9)


def test_symmetrize_rejects_a_symmetry_from_another_frame() -> None:
    crystal, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase)
    with pytest.raises(ValueError, match="specimen frame"):
        figure.symmetrize(SymmetrySpec.from_point_group("222", reference_frame=crystal))


def test_restricting_the_polar_range_drops_the_defocused_rim() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, resolution_deg=5.0)
    restricted = figure.restrict_polar_range(max_polar_deg=70.0)
    polar = np.degrees(
        np.arccos(np.clip(np.abs(restricted.sample_directions[:, 2]), -1.0, 1.0))
    )
    assert float(np.max(polar)) <= 70.0 + 1e-9
    assert len(restricted.intensities) < len(figure.intensities)


def test_restricting_to_an_empty_band_raises() -> None:
    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, resolution_deg=10.0)
    # The grid has rings every 10 degrees, so this band falls strictly between
    # two of them and contains nothing at all.
    with pytest.raises(ValueError, match="No sampled direction"):
        figure.restrict_polar_range(min_polar_deg=41.0, max_polar_deg=44.0)


def test_the_full_workflow_two_measurements_to_one_residual() -> None:
    """The end-to-end path this sprint exists to make possible.

    Two figures with different supports and incomparable scales — the normal
    situation, and previously a dead end — become one residual figure in m.r.d.
    """

    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        np.array([[0.0, 0.0, 0.0], [35.0, 20.0, 10.0]]), specimen_frame=specimen, phase=phase
    )
    measured = PoleFigure.from_orientations(orientations, make_pole(phase))
    coarse = uniform_grid(specimen, 12.0)
    modelled = (measured.on_grid(coarse, halfwidth_deg=18.0) * 850.0).normalize_to_mrd(
        integration_weights=coarse.weights
    )
    assert measured.sample_directions.shape != modelled.sample_directions.shape

    common = measured.integration_grid(resolution_deg=8.0)
    residual = measured.on_grid(common, halfwidth_deg=15.0).difference(
        modelled.on_grid(common, halfwidth_deg=15.0),
        left_label="measured",
        right_label="modelled",
    )
    assert len(residual) == len(common)
    # Both are normalized, so the residual must integrate to zero even where
    # the two figures disagree pointwise.
    assert_allclose(float(np.sum(common.weights * residual.values)), 0.0, atol=1e-9)
    assert residual.rms_deviation > 0.0
    assert "measured minus modelled" in residual.describe()


# --------------------------------------------------------------------------
# ODF residual QC product and its figure
# --------------------------------------------------------------------------


def odf_from_orientations(specimen: ReferenceFrame, phase: Phase) -> ODF:
    return ODF.from_orientations(
        OrientationSet.from_euler_angles(
            np.array([[0.0, 0.0, 0.0], [35.0, 20.0, 10.0], [12.0, 45.0, 33.0]]),
            specimen_frame=specimen,
            phase=phase,
        ),
        kernel=KernelSpec(halfwidth_deg=15.0),
    )


def test_residual_report_exposes_a_difference_figure_matching_its_residuals() -> None:
    _, specimen, phase = make_context()
    odf = odf_from_orientations(specimen, phase)
    grid = uniform_grid(specimen, 12.0)
    measured = odf.reconstruct_pole_figure(make_pole(phase)).on_grid(grid, halfwidth_deg=15.0)
    report = PoleFigureResidualReport.from_odf(odf, measured)

    figure = report.difference_figure()
    assert_allclose(figure.values, report.residuals, rtol=1e-12)
    assert figure.left_label == "recalculated"
    assert figure.right_label == "measured"
    assert len(figure) == report.observation_count
    assert figure.max_absolute_deviation == pytest.approx(report.max_absolute_error)


def test_a_perfect_reconstruction_has_an_identically_zero_residual_figure() -> None:
    """Comparing an ODF against the figure it itself predicts must be exact.

    This is the calibration of the whole QC path: any non-zero residual here
    would be a bug in the comparison, not a property of the data.

    Note the division by ``random_pole_density``. A discrete ODF's
    ``evaluate_pole_density`` returns a kernel-weighted response, so building a
    figure from it directly produces one in response units — and a pole figure
    is in multiples of random. The two must be put on the same scale before
    they can be subtracted at all, which is exactly what ``from_odf`` does
    internally.
    """

    _, specimen, phase = make_context()
    odf = odf_from_orientations(specimen, phase)
    pole = make_pole(phase)
    grid = uniform_grid(specimen, 12.0)
    predicted = PoleFigure(
        pole=pole,
        sample_directions=grid.vectors.values,
        intensities=odf.evaluate_pole_density(pole, grid.vectors.values)
        / random_pole_density(odf.kernel),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    report = PoleFigureResidualReport.from_odf(odf, predicted)
    assert report.difference_figure().max_absolute_deviation < 1e-12


def test_residual_report_describes_its_own_quality() -> None:
    _, specimen, phase = make_context()
    odf = odf_from_orientations(specimen, phase)
    grid = uniform_grid(specimen, 12.0)
    measured = odf.reconstruct_pole_figure(make_pole(phase)).on_grid(grid, halfwidth_deg=15.0)
    text = PoleFigureResidualReport.from_odf(odf, measured).describe()
    assert "relative residual norm" in text
    assert "recalculated minus measured" in text
    assert "acceptance test" in text


def test_residual_reports_for_pole_figures_yields_plottable_differences() -> None:
    _, specimen, phase = make_context()
    odf = odf_from_orientations(specimen, phase)
    grid = uniform_grid(specimen, 12.0)
    figures = [
        odf.reconstruct_pole_figure(
            CrystalPlane(MillerIndex(indices, phase=phase), phase=phase)
        ).on_grid(grid, halfwidth_deg=15.0)
        for indices in ((1, 1, 1), (2, 0, 0))
    ]
    reports = residual_reports_for_pole_figures(odf, figures)
    assert len(reports) == 2
    for report, figure in zip(reports, figures, strict=True):
        difference = report.difference_figure()
        assert difference.pole == figure.pole
        assert difference.project().shape == (len(figure.intensities), 2)


def test_difference_figure_spec_centres_the_diverging_scale_on_zero() -> None:
    """A diverging colormap is only honest with symmetric limits.

    Without them matplotlib centres the neutral colour on the data mean, so a
    residual that is positive everywhere would be drawn half in the colour that
    means negative.
    """

    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 15.0)
    figure = constant_density_figure(specimen, phase, value=1.0, resolution_deg=15.0)
    biased = PoleFigure(
        pole=figure.pole,
        sample_directions=grid.vectors.values,
        intensities=np.linspace(1.0, 3.0, len(grid)),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    spec = build_pole_figure_difference_spec(biased - figure)
    layer = spec.scatter_layers[0]
    assert layer.cmap == "pytex.diverging"
    assert layer.vmin == pytest.approx(-layer.vmax)
    assert layer.vmax == pytest.approx(2.0)
    assert "{111}" in spec.title


def test_difference_figure_spec_survives_an_identically_zero_residual() -> None:
    """A perfect fit must still render rather than divide by a zero range."""

    _, specimen, phase = make_context()
    figure = constant_density_figure(specimen, phase, resolution_deg=15.0)
    spec = build_pole_figure_difference_spec(figure - figure)
    layer = spec.scatter_layers[0]
    assert layer.vmin is not None and layer.vmax is not None
    assert layer.vmin < layer.vmax


def test_scatter_layer_rejects_half_specified_colour_limits() -> None:
    with pytest.raises(ValueError, match="together or not at all"):
        ScatterLayer2D(points=np.zeros((2, 2)), vmin=-1.0)
    with pytest.raises(ValueError, match="vmin < vmax"):
        ScatterLayer2D(points=np.zeros((2, 2)), vmin=1.0, vmax=1.0)


def test_plot_pole_figure_difference_renders() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, specimen, phase = make_context()
    sharp, spread, _ = two_figures_on_one_grid(specimen, phase)
    figure = plot_pole_figure_difference(
        sharp.difference(spread, left_label="sharp", right_label="spread")
    )
    assert isinstance(figure, matplotlib.figure.Figure)
    assert "sharp - spread" in figure.axes[0].get_title()
    plt.close(figure)


# --------------------------------------------------------------------------
# Discrete-ODF residuals are compared on the physical scale
# --------------------------------------------------------------------------


def test_random_pole_density_matches_the_response_of_a_random_texture() -> None:
    """The normalizing constant, checked against an actual random texture.

    A discrete ODF's pole density is a kernel-weighted response, and a random
    texture returns the kernel's spherical mean rather than 1. Sampling SO(3)
    uniformly and evaluating the response must reproduce that mean; the
    quadrature is then verified against the thing it stands for.
    """

    _, specimen, phase = make_context()
    crystal = phase.crystal_frame
    uniform = OrientationSet.from_equispaced_so3_grid(
        15.0,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    kernel = KernelSpec(halfwidth_deg=12.0)
    odf = ODF.from_orientations(uniform, kernel=kernel)
    grid = uniform_grid(specimen, 12.0)
    response = odf.evaluate_pole_density(make_pole(phase), grid.vectors.values)
    assert float(np.mean(response)) == pytest.approx(random_pole_density(kernel), rel=1e-3)
    # The scale error this guards against is one to two orders of magnitude,
    # not a rounding difference.
    assert random_pole_density(kernel) < 0.05


def test_a_discrete_odf_residual_is_small_against_the_texture_it_came_from() -> None:
    """The regression this exists to prevent.

    ``evaluate_pole_density`` returns a kernel response, not m.r.d., while a
    resampled pole figure is in m.r.d. Differencing them directly reported a
    relative residual of 0.99 for a *perfect* fit — a scale error masquerading
    as a total misfit, which would condemn every sound inversion.
    """

    _, specimen, phase = make_context()
    orientations = OrientationSet.from_euler_angles(
        np.array([[0.0, 0.0, 0.0], [35.0, 20.0, 10.0], [12.0, 45.0, 33.0]]),
        specimen_frame=specimen,
        phase=phase,
    )
    grid = uniform_grid(specimen, 9.0)
    odf = ODF.from_orientations(orientations, kernel=KernelSpec(halfwidth_deg=12.0))
    measured = PoleFigure.from_orientations(orientations, make_pole(phase)).on_grid(
        grid, halfwidth_deg=12.0
    )
    report = PoleFigureResidualReport.from_odf(odf, measured)
    assert report.relative_residual_norm < 0.1
    # Both sides must now sit on the same scale, which is the actual fix.
    assert float(np.mean(report.predicted_intensities)) == pytest.approx(
        float(np.mean(measured.intensities)), rel=0.05
    )


def test_a_harmonic_odf_residual_is_not_rescaled() -> None:
    """HarmonicODF pole densities are already m.r.d. and must be left alone.

    Applying the discrete correction here would introduce the very scale error
    it removes elsewhere, so the two paths are checked separately.
    """

    _, specimen, phase = make_context()
    grid = uniform_grid(specimen, 12.0)
    flat = PoleFigure(
        pole=make_pole(phase),
        sample_directions=grid.vectors.values,
        intensities=np.ones(len(grid)),
        specimen_frame=specimen,
        antipodal=True,
        sampling="sampled_density",
    )
    report = HarmonicODF.invert_pole_figures((flat,), degree_bandlimit=6, regularization=1e-6)
    residual = PoleFigureResidualReport.from_odf(report.odf, flat)
    assert float(np.mean(residual.predicted_intensities)) == pytest.approx(1.0, abs=0.1)
    assert residual.relative_residual_norm < 0.5
