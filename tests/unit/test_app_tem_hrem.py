"""The HRTEM workbench operations, checked against the physics they expose.

Nothing here compares against a recorded output of this code. The expectations
are algebraic identities of the wave aberration function — a two-fold
astigmatism along its own azimuth is a defocus offset, and its extremes are 90
degrees apart — or contract properties: a control the form declares must reach
the lens, and a term left at zero must not be reported as though it shaped the
result.
"""

from __future__ import annotations

import pytest

from pytex.app import REGISTRY

_SIMULATE = "tem.simulate_hrem"
_CTF = "tem.ctf_calculator"


def _ctf(**overrides: object) -> dict:
    request = {
        "beam_energy_kev": 300.0,
        "mode": "double_corrected",
        "defocus_angstrom": -50.0,
        "cs_um": 1.0,
        "focal_spread_angstrom": 10.0,
        "convergence_semiangle_mrad": 0.2,
        "aperture_cutoff_mrad": 0.0,
        "max_q_inv_angstrom": 2.5,
    }
    request.update(overrides)
    return REGISTRY.call(_CTF, request)


def _parameter_names(operation_id: str) -> set[str]:
    return {parameter.name for parameter in REGISTRY.get(operation_id).parameters}


class TestControlsReachTheLens:
    """Every aberration the core model carries is reachable from both forms."""

    @pytest.mark.parametrize("operation_id", [_SIMULATE, _CTF])
    def test_both_operations_declare_the_full_residual_aberration_set(
        self, operation_id: str
    ) -> None:
        names = _parameter_names(operation_id)
        for required in (
            "c5_mm",
            "astigmatism_angstrom",
            "astigmatism_angle_deg",
            "coma_angstrom",
            "coma_angle_deg",
            "trefoil_angstrom",
            "trefoil_angle_deg",
        ):
            assert required in names, f"{operation_id} cannot reach {required}"

    def test_the_simulation_no_longer_hard_codes_partial_coherence(self) -> None:
        """Focal spread was fixed by mode; the CTF form beside it exposed all three."""
        names = _parameter_names(_SIMULATE)
        assert {
            "focal_spread_angstrom",
            "convergence_semiangle_mrad",
            "aperture_cutoff_mrad",
        } <= names

    def test_the_ctf_declares_the_cut_it_reports(self) -> None:
        assert {"azimuth_deg", "max_q_inv_angstrom"} <= _parameter_names(_CTF)


class TestAzimuthalReporting:
    def test_a_round_lens_reports_isotropic_transfer_and_no_residual_terms(self) -> None:
        data = _ctf()["data"]

        assert data["has_azimuthal_aberrations"] is False
        assert data["residual_aberrations"] == []
        assert data["resolution_anisotropy_angstrom"] == pytest.approx(0.0, abs=1e-9)
        assert "isotropic" in data["azimuthal_description"]

    def test_astigmatism_along_its_azimuth_matches_the_equivalent_round_defocus(self) -> None:
        """C12 enters chi as a defocus offset of +C12 at its own azimuth."""
        astigmatic = _ctf(
            astigmatism_angstrom=20.0, astigmatism_angle_deg=30.0, azimuth_deg=30.0
        )["data"]
        equivalent = _ctf(defocus_angstrom=-30.0)["data"]

        assert astigmatic["point_resolution_angstrom"] == pytest.approx(
            equivalent["point_resolution_angstrom"], rel=1e-9
        )

    def test_the_reported_band_and_prose_state_the_anisotropy(self) -> None:
        result = _ctf(astigmatism_angstrom=20.0, astigmatism_angle_deg=30.0)
        data = result["data"]

        assert data["has_azimuthal_aberrations"] is True
        assert [term["symbol"] for term in data["residual_aberrations"]] == ["C12"]
        assert data["resolution_anisotropy_angstrom"] > 0.0
        # Two-fold astigmatism has period 180 degrees, so the coarsest direction
        # lies 90 degrees from the astigmatism azimuth.
        assert data["worst_azimuth_deg"] == pytest.approx(120.0, abs=5.0)
        assert len(data["azimuthal_transfer_min"]) == len(data["spatial_frequencies"])
        assert "C12" in result["summary"]
        assert "not the whole lens" in result["summary"]

    def test_the_band_brackets_the_cut_at_every_frequency(self) -> None:
        data = _ctf(coma_angstrom=25.0, coma_angle_deg=15.0, trefoil_angstrom=10.0)["data"]

        for low, value, high in zip(
            data["azimuthal_transfer_min"],
            data["transfer_function"],
            data["azimuthal_transfer_max"],
            strict=True,
        ):
            assert low - 1e-9 <= value <= high + 1e-9

    def test_the_frequency_range_control_sets_the_sampled_span(self) -> None:
        narrow = _ctf(max_q_inv_angstrom=1.5)["data"]
        wide = _ctf(max_q_inv_angstrom=4.0)["data"]

        assert narrow["spatial_frequencies"][-1] == pytest.approx(1.5)
        assert wide["spatial_frequencies"][-1] == pytest.approx(4.0)

    def test_a_zero_aperture_means_no_aperture_rather_than_a_closed_one(self) -> None:
        """A cutoff of zero must not mask the whole passband to nothing."""
        data = _ctf(aperture_cutoff_mrad=0.0)["data"]

        assert any(abs(value) > 1e-6 for value in data["transfer_function"])


class TestSimulationReportsWhatShapedTheImage:
    def test_residual_terms_appear_in_the_table_and_the_summary(self) -> None:
        result = REGISTRY.call(
            _SIMULATE,
            {
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 1,
                "beam_energy_kev": 200.0,
                "mode": "double_corrected",
                "defocus_angstrom": -30.0,
                "cs_um": 0.0,
                "astigmatism_angstrom": 15.0,
                "astigmatism_angle_deg": 45.0,
                "sampling_angstrom": 0.3,
            },
        )
        data = result["data"]
        metrics = {row["metric"] for row in result["table"]["rows"]}

        assert data["has_azimuthal_aberrations"] is True
        assert [term["symbol"] for term in data["residual_aberrations"]] == ["C12"]
        assert "Two-fold astigmatism C₁₂" in metrics
        assert "not resolved equally in every direction" in result["summary"]

    def test_a_round_run_does_not_recite_zero_valued_coefficients(self) -> None:
        result = REGISTRY.call(
            _SIMULATE,
            {
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 1,
                "beam_energy_kev": 200.0,
                "mode": "double_corrected",
                "defocus_angstrom": -30.0,
                "cs_um": 0.0,
                "sampling_angstrom": 0.3,
            },
        )
        metrics = {row["metric"] for row in result["table"]["rows"]}

        assert result["data"]["residual_aberrations"] == []
        assert not any(metric.startswith("Two-fold") for metric in metrics)
        assert not any(metric.startswith("Trefoil") for metric in metrics)
