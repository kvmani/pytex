"""Scientific and wire-contract tests for the CBED workbench module.

The claims worth pinning here are the ones a user would act on: which regime the
discs are in (because it decides whether a thickness can be measured at all),
that the rasterised image is registered to the disc geometry drawn over it, and
that the layer spacing read back from a HOLZ ring reproduces the lattice
parameter it came from.
"""

from __future__ import annotations

import base64

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError, UnsupportedRequestError

#: Silicon's cubic lattice parameter, CODATA-traceable via the silicon lattice
#: standard: a = 5.4310 A at 295 K. The [001] repeat is exactly this, which is
#: what the HOLZ ring radius has to reproduce.
SILICON_A_ANGSTROM = 5.4310


def pattern(**overrides: object) -> dict:
    request = {
        "phase": {"builtin": "si_diamond"},
        "zone_axis": [0, 0, 1],
        "beam_energy_kev": 200.0,
        "convergence_semi_angle_mrad": 3.0,
        "thickness_nm": 100.0,
        "method": "two-beam",
        "camera_constant_mm_angstrom": 180.0,
        "g_max_inv_angstrom": 1.2,
        "max_index": 4,
        "disc_samples": 41,
    }
    request.update(overrides)
    return REGISTRY.call("cbed.pattern", request)


def decode(image: dict) -> np.ndarray:
    """The transmitted greyscale bytes, as the panel reconstructs them."""

    raw = np.frombuffer(base64.b64decode(image["data"]), dtype=np.uint8)
    return raw.reshape((image["height"], image["width"]))


class TestRegime:
    """The convergence angle decides what the pattern can be used for."""

    def test_a_narrow_probe_keeps_the_discs_separated(self) -> None:
        result = pattern(convergence_semi_angle_mrad=3.0)
        assert result["data"]["regime"] == "kossel-moellenstedt"
        assert result["data"]["separated"] is True
        assert "independent rocking curve" in result["summary"]

    def test_a_wide_probe_overlaps_them(self) -> None:
        result = pattern(convergence_semi_angle_mrad=12.0)
        assert result["data"]["regime"] == "kossel"
        assert result["data"]["separated"] is False

    def test_the_regime_is_exactly_the_diameter_against_the_nearest_spacing(self) -> None:
        """Stated so a reader can check the claim rather than trust the label."""

        result = pattern(convergence_semi_angle_mrad=3.0)
        data = result["data"]
        assert 2 * data["disc_radius_mm"] < data["nearest_disc_separation_mm"]

    def test_an_overlapping_pattern_says_the_interference_is_not_modelled(self) -> None:
        result = pattern(convergence_semi_angle_mrad=12.0)
        assert any("interference" in note for note in result["notes"])


class TestRasterisedImage:
    """The picture and the vector overlay must describe the same pattern."""

    def test_the_image_is_square_greyscale_of_the_declared_size(self) -> None:
        image = pattern()["data"]["image"]
        assert image["encoding"] == "base64-gray8"
        assert image["width"] == image["height"]
        grid = decode(image)
        assert grid.shape == (image["height"], image["width"])

    def test_the_transmitted_disc_lands_at_the_centre_of_the_image(self) -> None:
        """The centre pixel must be lit, or the drawing is offset from the geometry."""

        image = pattern()["data"]["image"]
        grid = decode(image)
        middle = image["height"] // 2
        assert grid[middle, middle] > 0

    def test_every_disc_centre_falls_inside_the_stated_extent(self) -> None:
        data = pattern()["data"]
        extent = data["image"]["extent_mm"]
        for disc in data["discs"]:
            assert abs(disc["x_mm"]) <= extent
            assert abs(disc["y_mm"]) <= extent

    def test_a_disc_centre_maps_onto_lit_pixels(self) -> None:
        """The vector outline is drawn from `discs`; the intensity from `image`.

        If these two disagreed the panel would draw circles beside the discs
        rather than around them, which is the failure this test exists to catch.
        """

        data = pattern()
        image = data["data"]["image"]
        grid = decode(image)
        size = image["height"]
        extent = image["extent_mm"]
        for disc in data["data"]["discs"][:6]:
            column = round((disc["x_mm"] + extent) / (2 * extent) * (size - 1))
            row = round((extent - disc["y_mm"]) / (2 * extent) * (size - 1))
            assert grid[row, column] > 0, f"no intensity at the centre of {disc['hkl_label']}"

    def test_the_peak_intensity_is_reported_so_the_bytes_mean_something(self) -> None:
        image = pattern()["data"]["image"]
        assert image["peak_intensity"] > 0.0
        assert decode(image).max() == 255


class TestDiscTable:
    """The rows behind the figure."""

    def test_the_transmitted_disc_comes_first_and_is_named(self) -> None:
        rows = pattern()["table"]["rows"]
        assert "direct" in rows[0]["hkl_label"]
        assert rows[0]["g_inv_angstrom"] == pytest.approx(0.0)
        assert rows[0]["d_angstrom"] is None

    def test_hover_columns_are_the_export_table_columns(self) -> None:
        result = pattern()
        assert result["data"]["columns"] == result["table"]["columns"]

    def test_silicon_down_001_draws_only_hk0_reflections(self) -> None:
        """Every zeroth-Laue-zone reflection of [001] is perpendicular to it."""

        for disc in pattern()["data"]["discs"][1:]:
            assert disc["indices"][2] == 0

    def test_the_diamond_glide_extinguishes_200_through_its_structure_factor(self) -> None:
        """(200) is drawn, and is dark. Both halves of that matter.

        The centering filter removes only what F-centering forbids, and (200)
        is all-even so it survives it. What kills (200) in the diamond structure
        is the two-atom basis: the glide puts the second atom exactly out of
        phase, and `|F|` goes to zero. So the disc is enumerated and drawn with
        no intensity, which is the honest rendering — dropping it would claim a
        systematic absence where there is a cancellation.
        """

        discs = {tuple(disc["indices"]): disc for disc in pattern()["data"]["discs"]}
        allowed = [disc for key, disc in discs.items() if abs(key[0]) == 2 and abs(key[1]) == 2]
        extinguished = [
            disc for key, disc in discs.items() if abs(key[0]) == 2 and key[1] == 0 and key[2] == 0
        ]
        assert allowed, "(220) is the first allowed diamond-cubic reflection and must be drawn"
        assert all(disc["structure_factor_amplitude"] > 1.0 for disc in allowed)
        for disc in extinguished:
            assert disc["structure_factor_amplitude"] == pytest.approx(0.0, abs=1e-6)
            assert disc["mean_intensity"] == pytest.approx(0.0, abs=1e-9)

    def test_d_spacing_is_the_reciprocal_of_the_drawn_g(self) -> None:
        for disc in pattern()["data"]["discs"][1:]:
            assert disc["d_angstrom"] == pytest.approx(1.0 / disc["g_inv_angstrom"])


class TestRefusals:
    """A refusal must name the control that caused it."""

    def test_symmetry_determination_is_refused_on_a_two_beam_pattern(self) -> None:
        with pytest.raises(UnsupportedRequestError) as raised:
            pattern(method="two-beam", determine_point_group=True)
        assert raised.value.details["field"] == "determine_point_group"
        assert "excitation error" in raised.value.message

    def test_a_phase_without_atoms_is_refused_with_the_reason(self) -> None:
        with pytest.raises(InvalidInputError) as raised:
            pattern(
                phase={
                    "a": 4.0,
                    "b": 4.0,
                    "c": 4.0,
                    "alpha": 90.0,
                    "beta": 90.0,
                    "gamma": 90.0,
                    "point_group": "m-3m",
                }
            )
        assert "extinction distance" in raised.value.message

    def test_a_zero_convergence_angle_is_refused_by_the_declared_bound(self) -> None:
        """A parallel beam is SAED, not a degenerate CBED, so it never reaches the simulator."""

        with pytest.raises(InvalidInputError) as raised:
            pattern(convergence_semi_angle_mrad=0.0)
        assert raised.value.details["field"] == "convergence_semi_angle_mrad"
        assert "at least 0.1" in raised.value.message


class TestThicknessFromFringes:
    """The two-beam inversion, and the check it carries with it."""

    @staticmethod
    def fit(**request: object) -> dict:
        return REGISTRY.call("cbed.thickness_from_fringes", request)

    def test_three_minima_give_a_thickness_and_an_extinction_distance(self) -> None:
        result = self.fit(s1=0.0071, s2=0.0128, s3=0.0184)
        assert result["data"]["thickness_angstrom"] > 0.0
        assert result["data"]["extinction_distance_angstrom"] > 0.0
        assert result["data"]["thickness_nm"] == pytest.approx(
            result["data"]["thickness_angstrom"] / 10.0
        )

    def test_the_fit_recovers_a_thickness_the_minima_were_built_from(self) -> None:
        """The round trip: generate minima from known t and xi, then fit them.

        The expected value has independent provenance — it is the closed form
        `s_n = sqrt(n^2/t^2 - 1/xi^2)` inverted, not a previous program output.
        """

        thickness = 1500.0
        extinction = 400.0
        # The first *observable* minimum is the lowest n with n/t > 1/xi: below
        # it the closed form has no real root, because the effective deviation
        # parameter can never fall that low. Here t/xi = 3.75, so n starts at 4.
        first_order = 4
        minima = {
            f"s{slot}": float(
                np.sqrt((first_order + slot - 1) ** 2 / thickness**2 - 1.0 / extinction**2)
            )
            for slot in (1, 2, 3, 4)
        }
        result = self.fit(**minima)
        assert result["data"]["thickness_angstrom"] == pytest.approx(thickness, rel=1e-6)
        assert result["data"]["extinction_distance_angstrom"] == pytest.approx(
            extinction, rel=1e-6
        )

    def test_the_linearised_points_are_the_columns_the_plot_draws(self) -> None:
        result = self.fit(s1=0.0071, s2=0.0128, s3=0.0184)
        rows = result["data"]["fit_points"]
        assert len(rows) == 3
        for row in rows:
            assert row["inverse_order_squared"] == pytest.approx(1.0 / row["order"] ** 2)
            assert row["s_over_n_squared"] == pytest.approx(
                (row["excitation_error_inv_angstrom"] / row["order"]) ** 2
            )

    def test_two_minima_are_enough_and_use_exactly_two_points(self) -> None:
        """The optional minima must stay optional.

        `s3` deliberately carries no default: an optional parameter that has one
        is substituted whenever it is omitted, so a default there would make a
        two-minimum fit impossible to ask for — the form would silently add a
        third measurement the user never made.
        """

        result = self.fit(s1=0.0071, s2=0.0128)
        assert len(result["data"]["fit_points"]) == 2
        assert result["data"]["thickness_angstrom"] > 0.0

    def test_minima_beyond_the_third_are_optional_too(self) -> None:
        result = self.fit(s1=0.0071, s2=0.0128, s3=0.0184, s4=0.0240)
        assert len(result["data"]["fit_points"]) == 4


class TestHOLZRings:
    """The one dimension the zone-axis pattern cannot see."""

    @staticmethod
    def rings(**overrides: object) -> dict:
        request = {
            "phase": {"builtin": "si_diamond"},
            "zone_axis": [0, 0, 1],
            "beam_energy_kev": 200.0,
            "orders": 2,
            "camera_constant_mm_angstrom": 180.0,
        }
        request.update(overrides)
        return REGISTRY.call("cbed.holz_rings", request)

    def test_the_layer_spacing_reproduces_the_silicon_lattice_parameter(self) -> None:
        """`H = 1/|r_uvw|`, and down [001] that repeat is `a` exactly.

        This is the whole point of the measurement, so it is checked against the
        cited lattice parameter rather than against a stored output.
        """

        data = self.rings()["data"]
        assert data["real_space_repeat_angstrom"] == pytest.approx(SILICON_A_ANGSTROM, abs=5e-3)
        assert data["layer_spacing_inv_angstrom"] == pytest.approx(
            1.0 / SILICON_A_ANGSTROM, rel=1e-3
        )

    def test_ring_radii_follow_the_square_root_of_the_order(self) -> None:
        """`G_n = sqrt(2nH/lambda)`, so the ratio of the first two is sqrt(2)."""

        rings = self.rings(orders=2)["data"]["rings"]
        assert len(rings) >= 2
        first, second = rings[0], rings[1]
        expected = np.sqrt(second["order"] / first["order"])
        assert second["radius_inv_angstrom"] / first["radius_inv_angstrom"] == pytest.approx(
            expected, rel=1e-9
        )

    def test_the_millimetre_radius_is_the_camera_constant_times_the_reciprocal_one(self) -> None:
        rings = self.rings(camera_constant_mm_angstrom=250.0)["data"]["rings"]
        for ring in rings:
            assert ring["radius_mm"] == pytest.approx(ring["radius_inv_angstrom"] * 250.0)

    def test_a_shorter_wavelength_pushes_the_rings_outward(self) -> None:
        """`G_n` scales as `lambda^{-1/2}`, so a higher voltage widens the ring."""

        low = self.rings(beam_energy_kev=100.0)["data"]["rings"][0]
        high = self.rings(beam_energy_kev=300.0)["data"]["rings"][0]
        assert high["radius_inv_angstrom"] > low["radius_inv_angstrom"]


class TestNarration:
    """Every CBED operation reports into the centralized log."""

    def test_the_pattern_reports_its_regime_and_disc_count(self) -> None:
        from pytex.app.contracts import execute

        envelope, status = execute(
            "cbed.pattern",
            {
                "phase": {"builtin": "si_diamond"},
                "zone_axis": [0, 0, 1],
                "beam_energy_kev": 200.0,
                "convergence_semi_angle_mrad": 3.0,
                "thickness_nm": 100.0,
                "method": "two-beam",
                "camera_constant_mm_angstrom": 180.0,
                "g_max_inv_angstrom": 1.2,
                "max_index": 4,
                "disc_samples": 21,
            },
        )
        assert status == 200
        messages = [record["message"] for record in envelope["log"]]
        assert any("discs, separated" in message for message in messages)
        assert any("Drawing the convergent-beam discs" in message for message in messages)
