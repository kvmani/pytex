"""Finite-thickness SAED shape-factor tests against the slab transform."""

from __future__ import annotations

import numpy as np
import pytest

from pytex import FiniteThicknessShapeFactor, from_json_contract, to_json_contract
from pytex.diffraction import KinematicSimulationConfig


def test_rectangular_slab_analytic_values() -> None:
    """A 100 A slab has its half-zero value at 1/(2t) and first zero at 1/t."""

    model = FiniteThicknessShapeFactor(100.0)
    excitation = np.array([0.0, 0.005, 0.01])

    assert model.intensity_factor(excitation) == pytest.approx(
        np.array([1.0, 4.0 / np.pi**2, 0.0]), abs=1e-15
    )
    assert model.first_zero_inv_angstrom == pytest.approx(0.01)


def test_shape_factor_is_even_and_scalar_safe() -> None:
    model = FiniteThicknessShapeFactor(75.0)

    assert model.intensity_factor(-0.004) == pytest.approx(model.intensity_factor(0.004))
    assert model.amplitude_factor(0.0) == 1.0


def test_shape_factor_contract_and_explanation_round_trip() -> None:
    model = FiniteThicknessShapeFactor(120.0)
    payload = to_json_contract(model)
    restored = from_json_contract(payload)

    assert restored == model
    assert payload["schema_id"] == "pytex.diffraction.finite_thickness_shape_factor"
    assert "sinc^2(t s_g)" in model.describe()
    assert "120 angstrom" in model.describe()


@pytest.mark.parametrize("thickness", [0.0, -1.0, np.inf, np.nan])
def test_invalid_thickness_raises(thickness: float) -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        FiniteThicknessShapeFactor(thickness)


def test_kinematic_config_rejects_ambiguous_shape_models() -> None:
    with pytest.raises(ValueError, match="mutually exclusive"):
        KinematicSimulationConfig(
            foil_thickness_angstrom=100.0,
            relrod_sigma_inv_angstrom=0.01,
        )
