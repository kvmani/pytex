"""Pole-figure resampling, m.r.d. normalization and arithmetic.

These surfaces exist as one dependency chain: two pole figures cannot be
combined until they share a support (resampling) and a physical scale
(m.r.d.), so the tests are ordered the same way.
"""

from __future__ import annotations

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
