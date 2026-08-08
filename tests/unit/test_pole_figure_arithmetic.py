"""Pole-figure resampling, m.r.d. normalization and arithmetic.

These surfaces exist as one dependency chain: two pole figures cannot be
combined until they share a support (resampling) and a physical scale
(m.r.d.), so the tests are ordered the same way.
"""

from __future__ import annotations

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
    S2Grid,
    SymmetrySpec,
)
from pytex.texture import PoleFigure


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
