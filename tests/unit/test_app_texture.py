"""The texture panel, checked against answers crystallography fixes in advance.

This is why the panel builds a model texture rather than loading a file: a model
has a known answer and a data set does not.

Three independent checks anchor everything here.

* **A random texture is 1 m.r.d. everywhere.** Not approximately, and not in
  arbitrary units — exactly 1, in every view, for every symmetry. It is the
  definition of the scale, and it caught the defect this module was written
  around: the ODF sections were 24 times too large and labelled "m.r.d.".
* **The area-weighted mean of any pole figure is 1 m.r.d.**, whatever the
  texture, because that is what normalising to a random distribution means.
* **A named component puts its poles where its Miller label says.** Goss is
  {011}<100>, so the (011) poles of a Goss texture sit at the centre of the
  figure, which is ND. No reference figure is needed to check that; the notation
  asserts it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY

NICKEL = {"builtin": "ni_fcc"}
ZIRCONIUM = {"builtin": "zr_hcp"}

#: A texture model shared by the tests, with the scatter fixed by a seed.
BASE: dict[str, object] = {
    "phase": NICKEL,
    "model": "fcc_rolling",
    "spread_deg": 10.0,
    "grain_count": 200,
    "halfwidth_deg": 10.0,
    "seed": 7,
}


def call(operation: str, **request: object) -> dict:
    return REGISTRY.call(operation, request)


def pole_figure(**overrides: object) -> dict:
    request = dict(BASE, pole=[1, 1, 1], projection="equal_area", resolution_deg=5.0)
    request.update(overrides)
    return call("texture.pole_figure", **request)


def odf_sections(**overrides: object) -> dict:
    request = dict(BASE, section_resolution_deg=5.0)
    request.update(overrides)
    return call("texture.odf_sections", **request)


def inverse_pole_figure(**overrides: object) -> dict:
    request = dict(BASE, sample_direction="nd", projection="equal_area")
    request.update(overrides)
    return call("texture.inverse_pole_figure", **request)


def section_values(result: dict) -> np.ndarray:
    return np.array(
        [
            value
            for section in result["data"]["sections"]
            for row in section["densities"]
            for value in row
        ]
    )


class TestTheRandomBaseline:
    """The definition of the m.r.d. scale, checked as such."""

    @pytest.mark.parametrize("phase", [NICKEL, ZIRCONIUM], ids=["cubic", "hexagonal"])
    def test_a_random_texture_is_one_mrd_across_the_pole_figure(self, phase: dict) -> None:
        """Flat at 1, to within the noise of a finite number of grains.

        The tolerance is on the *spread*, not on the mean: the mean is 1 by
        construction and is checked exactly elsewhere. What this asserts is that
        a random texture has no features — that nothing in the pipeline invents
        a peak where the material has none.
        """

        data = pole_figure(phase=phase, model="random", grain_count=600)["data"]
        assert data["mean_mrd"] == pytest.approx(1.0, abs=0.01)
        assert data["max_mrd"] < 1.35, "a random texture must not show a peak"
        assert data["min_mrd"] > 0.65, "a random texture must not show a hole"

    @pytest.mark.parametrize("phase", [NICKEL, ZIRCONIUM], ids=["cubic", "hexagonal"])
    def test_a_random_texture_reads_one_mrd_in_every_view(self, phase: dict) -> None:
        """The ODF sections must be on the same scale as the pole figure.

        They were not. `ODF.evaluate(normalized=True)` normalises the kernel
        over the whole of SO(3), while a symmetry-aware evaluation folds every
        query into the fundamental zone — which is 1/|G| of it. So a random
        texture read |G|: 23.9 for m-3m and 11.9 for 6/mmm, against operator
        counts of 24 and 12, in a column labelled "m.r.d.".

        Both symmetries are checked because a single one cannot distinguish a
        correct normalisation from a constant fudge factor.
        """

        values = section_values(odf_sections(phase=phase, model="random", grain_count=600))
        assert float(values.mean()) == pytest.approx(1.0, abs=0.05)
        assert float(values.max()) < 1.6, "a random texture must not peak in Euler space either"


class TestTheMrdScaleIsClosed:
    """The arithmetic identity that makes two pole figures comparable."""

    @pytest.mark.parametrize(
        "model", ["random", "cube", "goss", "brass", "copper", "s", "fcc_rolling"]
    )
    def test_the_area_weighted_mean_is_one_whatever_the_texture(self, model: str) -> None:
        """1 m.r.d. on average is what normalising to a random distribution means.

        It holds for a sharp single component exactly as it does for a random
        one, which is the whole point: the mean carries no information, so every
        feature of a pole figure is a departure from it. A figure whose mean is
        not 1 has not been normalised and its numbers mean nothing outside
        itself.
        """

        assert pole_figure(model=model)["data"]["mean_mrd"] == pytest.approx(1.0, abs=0.01)

    def test_the_mean_is_area_weighted_and_not_a_raw_row_average(self) -> None:
        """The unweighted mean of the same grid is biased, and must differ.

        An equispaced grid on a hemisphere is not equal-area, so averaging the
        rows of the exported table without weights gives the wrong answer. If
        the two ever agreed, the weighting would have been dropped.
        """

        result = pole_figure(model="goss")
        rows = np.array([row["mrd"] for row in result["table"]["rows"]])
        assert result["data"]["mean_mrd"] == pytest.approx(1.0, abs=0.01)
        assert float(rows.mean()) != pytest.approx(result["data"]["mean_mrd"], abs=1e-3)


class TestComponentsLandWhereTheirLabelSays:
    """A Miller label is a testable claim about where the poles go."""

    def test_goss_puts_its_011_at_the_centre(self) -> None:
        """Goss is {011}<100>: the {011} plane lies in the sheet plane.

        So the (011) pole is along ND, which is the centre of the figure — polar
        angle 0. This checks the entire chain at once: Euler convention, the
        crystal-to-specimen mapping, the pole family, and the projection.
        """

        result = pole_figure(model="goss", pole=[0, 1, 1], spread_deg=6.0)
        peak = max(result["table"]["rows"], key=lambda row: row["mrd"])
        assert peak["polar_deg"] == pytest.approx(0.0, abs=6.0)
        assert result["data"]["max_mrd"] > 3.0, "a sharp component must show a real peak"

    def test_cube_puts_its_001_at_the_centre(self) -> None:
        """Cube is {001}<100>, so its (001) pole is along ND too."""

        result = pole_figure(model="cube", pole=[0, 0, 1], spread_deg=6.0)
        peak = max(result["table"]["rows"], key=lambda row: row["mrd"])
        assert peak["polar_deg"] == pytest.approx(0.0, abs=6.0)

    def test_cube_and_goss_are_not_the_same_figure(self) -> None:
        """Both peak at the centre on their own plane, and differ on {111}.

        Without this, the two tests above would both pass on a pipeline that
        ignored the component entirely and always peaked at the centre.
        """

        cube = np.array(
            [row["mrd"] for row in pole_figure(model="cube", pole=[1, 1, 1])["table"]["rows"]]
        )
        goss = np.array(
            [row["mrd"] for row in pole_figure(model="goss", pole=[1, 1, 1])["table"]["rows"]]
        )
        assert not np.allclose(cube, goss, atol=0.1)

    def test_a_sharper_component_gives_a_higher_peak(self) -> None:
        """Halving the scatter must concentrate the poles, not just move them."""

        tight = pole_figure(model="goss", pole=[0, 1, 1], spread_deg=4.0)["data"]["max_mrd"]
        loose = pole_figure(model="goss", pole=[0, 1, 1], spread_deg=20.0)["data"]["max_mrd"]
        assert tight > loose


class TestProjections:
    """Both projections must share a rim, or two figures cannot be compared."""

    @pytest.mark.parametrize("method", ["equal_area", "stereographic"])
    def test_every_point_lies_inside_the_unit_disc(self, method: str) -> None:
        for row in pole_figure(projection=method)["table"]["rows"]:
            assert math.hypot(row["x"], row["y"]) <= 1.0 + 1e-9

    @pytest.mark.parametrize("method", ["equal_area", "stereographic"])
    def test_the_rim_is_ninety_degrees_from_the_centre(self, method: str) -> None:
        rim = [
            row
            for row in pole_figure(projection=method)["table"]["rows"]
            if abs(row["polar_deg"] - 90.0) < 1e-6
        ]
        assert rim, "the hemisphere grid must reach the equator"
        for row in rim:
            assert math.hypot(row["x"], row["y"]) == pytest.approx(1.0, abs=1e-9)

    def test_the_two_projections_are_not_the_same(self) -> None:
        radius = {}
        for method in ("equal_area", "stereographic"):
            rows = pole_figure(projection=method)["table"]["rows"]
            middle = next(row for row in rows if 40.0 < row["polar_deg"] < 50.0)
            radius[method] = math.hypot(middle["x"], middle["y"])
        assert radius["equal_area"] != pytest.approx(radius["stereographic"], abs=1e-6)


class TestInversePoleFigure:
    """One point per grain, folded into the fundamental sector."""

    def test_every_grain_appears_once(self) -> None:
        result = inverse_pole_figure()
        # Five components at 200 grains each.
        assert result["data"]["grain_count"] == 5 * 200
        assert len(result["table"]["rows"]) == 5 * 200

    def test_the_sector_has_three_corners(self) -> None:
        """The cubic standard triangle: [001], [101], [111]."""

        vertices = inverse_pole_figure()["data"]["sector_vertices"]
        assert len(vertices) == 3
        for x, y in vertices:
            assert math.hypot(x, y) <= 1.0 + 1e-9

    def test_every_point_is_inside_the_disc(self) -> None:
        for row in inverse_pole_figure()["table"]["rows"]:
            assert math.hypot(row["x"], row["y"]) <= 1.0 + 1e-9

    def test_a_cube_texture_puts_001_along_nd(self) -> None:
        """Cube is {001}<100>: [001] is along ND, so every grain is at a corner.

        The polar angle reported is measured from the crystal [001], so for a
        cube texture looked at along ND it must be near zero for every grain.
        """

        rows = inverse_pole_figure(model="cube", spread_deg=5.0)["table"]["rows"]
        polar = np.array([row["polar_deg"] for row in rows])
        assert float(np.median(polar)) < 12.0

    def test_the_three_specimen_axes_give_different_answers(self) -> None:
        """A texture is not isotropic, so ND, RD and TD cannot agree."""

        seen = {
            axis: inverse_pole_figure(sample_direction=axis)["table"]["rows"][0]["polar_deg"]
            for axis in ("nd", "rd", "td")
        }
        assert len(set(round(value, 6) for value in seen.values())) == 3


class TestOdfSections:
    """Three sections, and the coordinate degeneracy they cannot hide."""

    def test_three_sections_are_returned(self) -> None:
        sections = odf_sections()["data"]["sections"]
        assert [section["phi2_deg"] for section in sections] == [0.0, 45.0, 65.0]

    def test_a_sharp_component_gives_a_strong_peak_in_the_published_range(self) -> None:
        """A tight cube reaches tens of m.r.d., not hundreds or fractions.

        The number is what the defect this module found made wrong: before the
        fundamental-zone normalisation it read 737 m.r.d. for this case, which
        is not a texture strength any material has. The bound is loose on
        purpose — it is checking an order of magnitude, which is exactly what
        was wrong.
        """

        peak = odf_sections(model="cube", spread_deg=8.0)["data"]["max_mrd"]
        assert 5.0 < peak < 200.0

    def test_the_table_carries_every_section_grid_point(self) -> None:
        result = odf_sections()
        sections = result["data"]["sections"]
        expected = sum(
            len(section["phi1_deg"]) * len(section["big_phi_deg"]) for section in sections
        )
        assert len(result["table"]["rows"]) == expected

    def test_a_rolling_texture_is_stronger_than_random_but_not_a_single_crystal(self) -> None:
        random_peak = odf_sections(model="random", grain_count=600)["data"]["max_mrd"]
        rolling_peak = odf_sections(model="fcc_rolling")["data"]["max_mrd"]
        assert rolling_peak > 2.0 * random_peak


class TestVolumeFractions:
    """The quantitative reading a figure only supports qualitatively."""

    def test_a_single_component_texture_is_mostly_that_component(self) -> None:
        fractions = pole_figure(model="goss", spread_deg=8.0)["data"]["component_fractions"]
        assert len(fractions) == 1
        assert fractions[0]["component"] == "goss"
        assert fractions[0]["fraction"] > 0.5

    def test_the_rolling_texture_lists_all_five_components(self) -> None:
        fractions = pole_figure(model="fcc_rolling")["data"]["component_fractions"]
        assert [entry["component"] for entry in fractions] == [
            "cube",
            "goss",
            "brass",
            "copper",
            "s",
        ]
        assert all(0.0 <= entry["fraction"] <= 1.0 for entry in fractions)

    def test_a_random_texture_has_no_named_components(self) -> None:
        assert pole_figure(model="random")["data"]["component_fractions"] == []


class TestReproducibility:
    """The seed must actually fix the scatter, or nothing above is repeatable."""

    def test_the_same_seed_gives_the_same_figure(self) -> None:
        first = pole_figure(seed=11)["data"]["max_mrd"]
        second = pole_figure(seed=11)["data"]["max_mrd"]
        assert first == second

    def test_a_different_seed_gives_a_different_sample(self) -> None:
        assert pole_figure(seed=11)["data"]["max_mrd"] != pole_figure(seed=12)["data"]["max_mrd"]


class TestResultContract:
    """What the panel and the export both read."""

    def test_the_declared_columns_match_the_table(self) -> None:
        for result in (pole_figure(), inverse_pole_figure(), odf_sections()):
            declared = [column["key"] for column in result["data"]["columns"]]
            assert declared == [column["key"] for column in result["table"]["columns"]]

    def test_every_result_explains_itself(self) -> None:
        for result in (pole_figure(), inverse_pole_figure(), odf_sections()):
            assert len(result["summary"]) > 80
            assert result["notes"]
            assert result["citations"]

    def test_the_pole_figure_table_and_the_plotted_points_agree(self) -> None:
        result = pole_figure()
        assert len(result["table"]["rows"]) == len(result["data"]["points"])
        assert result["table"]["rows"][0]["x"] == result["data"]["points"][0]["x"]
