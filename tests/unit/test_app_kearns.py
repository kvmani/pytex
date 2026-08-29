"""The Kearns panel: the three routes, and the agreement that justifies them.

The point of this suite is not that each operation returns a number. It is that
the three operations describe **one synthetic specimen whose true ``f`` is
known**, and all three recover it. A route that agrees with itself proves
nothing; three independent routes agreeing with an independently known answer is
the claim the panel makes to a user, so it is the claim under test.

The truth comes from the exact route on a large model texture, and the defaults
of the two diffraction routes were generated from that same texture. If someone
regenerates the defaults, these tests fail until the truth is regenerated with
them — which is the intended coupling.
"""

from __future__ import annotations

from typing import Any

import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.registry import REGISTRY

# The Kearns parameter of the synthetic specimen the panel's defaults describe:
# a basal fibre of 30 degree spread about ND, evaluated by the exact route on
# 200,000 orientations. Every default in `pytex.app.services.kearns` was
# generated from this texture.
TRUE_F_ND = 0.6423

#: How far a route may sit from the truth before the panel is misleading. The
#: diffraction routes carry binning and interpolation error the exact route does
#: not, and 0.01 in f is far below the difference between any two real specimens.
ROUTE_TOLERANCE = 0.01


def defaults(operation_id: str) -> dict[str, Any]:
    """The request the form sends when nothing has been touched."""

    spec = next(item for item in REGISTRY.operations() if item.id == operation_id)
    return {
        parameter.name: parameter.default
        for parameter in spec.parameters
        if parameter.default is not None
    }


def run(operation_id: str, **overrides: Any) -> dict[str, Any]:
    request = defaults(operation_id)
    request.update(overrides)
    return REGISTRY.call(operation_id, request)


def value_of(result: dict[str, Any], label: str) -> float:
    entry = next(item for item in result["data"]["directions"] if item["label"] == label)
    return float(entry["f"])


class TestTheThreeRoutesAgree:
    """The claim the panel exists to support."""

    def test_exact_route_recovers_the_model_truth(self) -> None:
        result = run("kearns.from_orientations")
        assert value_of(result, "ND") == pytest.approx(TRUE_F_ND, abs=ROUTE_TOLERANCE)

    def test_diffractogram_route_recovers_the_same_number(self) -> None:
        """Thirteen peak heights, no orientation data, and the same answer."""

        result = run("kearns.from_diffractogram")
        assert value_of(result, "ND") == pytest.approx(TRUE_F_ND, abs=ROUTE_TOLERANCE)

    def test_tilt_profile_route_recovers_the_same_number(self) -> None:
        result = run("kearns.from_tilt_profile")
        assert value_of(result, "ND") == pytest.approx(TRUE_F_ND, abs=ROUTE_TOLERANCE)

    def test_the_three_agree_with_each_other(self) -> None:
        """Stated separately because it is what a user is told to check."""

        exact = value_of(run("kearns.from_orientations"), "ND")
        diffractogram = value_of(run("kearns.from_diffractogram"), "ND")
        profile = value_of(run("kearns.from_tilt_profile"), "ND")
        assert diffractogram == pytest.approx(exact, abs=ROUTE_TOLERANCE)
        assert profile == pytest.approx(exact, abs=ROUTE_TOLERANCE)


class TestTheTriadAndItsClosure:
    def test_the_exact_route_closes_the_triad_but_says_it_proves_nothing(self) -> None:
        """The sum is 1, and the summary must not present that as a passed check.

        These three values are the diagonal of one pole orientation tensor whose
        trace is identically 1, so they close whatever the data were. Claiming
        the closure as evidence of a sound measurement is false reassurance, and
        `TestTheUploadRoutes.test_closure_passes_while_the_answer_is_wrong`
        shows a case where it is closed and wrong by more than half.
        """

        result = run("kearns.from_orientations")
        assert result["data"]["triad_sum"] == pytest.approx(1.0, abs=1e-9)
        assert result["data"]["is_orthonormal_triad"] is True
        assert "closure by construction, not evidence" in result["summary"]

    def test_a_fibre_texture_is_transversely_isotropic(self) -> None:
        """A fibre about ND must give equal f along RD and TD, by construction."""

        result = run("kearns.from_orientations", spread_deg=30.0, count=50000.0)
        assert value_of(result, "RD") == pytest.approx(value_of(result, "TD"), abs=0.01)

    def test_widening_the_fibre_weakens_it_towards_random(self) -> None:
        """f falls monotonically towards 1/3 as the spread grows — but does not reach it.

        The model cannot produce a random texture, and the panel says so. A
        Gaussian truncated to the quadrant is still concentrated towards the
        axis at any width, so the widest fibre the control allows lands near
        0.37 rather than at 1/3. Pinned here because it would otherwise look
        like a defect to the next reader, and because the help text makes the
        claim in exactly these terms.
        """

        values = [
            value_of(run("kearns.from_orientations", spread_deg=spread, count=80000.0), "ND")
            for spread in (30.0, 60.0, 90.0)
        ]
        assert values[0] > values[1] > values[2] > 1.0 / 3.0
        assert values[2] == pytest.approx(0.374, abs=0.01)

    def test_a_sharp_fibre_concentrates_along_its_axis(self) -> None:
        result = run("kearns.from_orientations", spread_deg=5.0, count=20000.0)
        assert value_of(result, "ND") > 0.97

    def test_the_fibre_axis_is_honoured(self) -> None:
        """A fibre about RD must put the large value on RD, not on ND."""

        result = run("kearns.from_orientations", fibre_axis="RD", spread_deg=20.0)
        assert value_of(result, "RD") > value_of(result, "ND")
        assert value_of(result, "RD") > 0.7

    def test_a_single_section_route_reports_no_closure(self) -> None:
        """One section gives one number, and the panel must not imply otherwise."""

        result = run("kearns.from_diffractogram")
        assert result["data"]["triad_sum"] is None
        assert "no closure check" in result["summary"]


class TestTheEvidenceTravelsWithTheNumber:
    def test_the_orientation_tensor_has_unit_trace(self) -> None:
        tensor = run("kearns.from_orientations")["data"]["orientation_tensor"]
        assert sum(tensor[index][index] for index in range(3)) == pytest.approx(1.0, abs=1e-9)

    def test_the_diffractogram_reports_every_reflection_at_its_tilt(self) -> None:
        result = run("kearns.from_diffractogram")
        rows = result["table"]["rows"]
        assert len(rows) == 13
        # Sorted by tilt, and spanning the full range: the basal pole at 0 and
        # the prism pole at 90 are what make the integral determinate.
        assert rows[0]["basal_tilt_deg"] == pytest.approx(0.0, abs=1e-6)
        assert rows[-1]["basal_tilt_deg"] == pytest.approx(90.0, abs=1e-6)
        assert rows[0]["plane"] == "(0002)"

    def test_the_diffractogram_returns_the_profile_it_integrated(self) -> None:
        profile = run("kearns.from_diffractogram")["data"]["profile"]
        assert len(profile["polar_deg"]) == len(profile["intensity"]) >= 5

    def test_the_quadrature_contributions_sum_to_f(self) -> None:
        """The tilt-profile table must add up to the number above it."""

        result = run("kearns.from_tilt_profile")
        total = sum(row["contribution"] for row in result["table"]["rows"])
        assert total == pytest.approx(value_of(result, "ND"), abs=1e-9)

    def test_the_harris_coefficients_average_to_one(self) -> None:
        rows = run("kearns.from_diffractogram")["table"]["rows"]
        mean = sum(row["harris"] for row in rows) / len(rows)
        assert mean == pytest.approx(1.0, abs=1e-9)

    def test_the_report_carries_its_own_prose(self) -> None:
        described = run("kearns.from_orientations")["data"]["describe"]
        assert "Kearns parameter f" in described
        assert "Triad sum" in described
        # The random baseline must be stated, or the number has no scale.
        assert "0.3333" in described

    def test_the_phase_name_is_not_shadowed_by_the_specification(self) -> None:
        """The report's own JSON contract key must survive the panel's additions."""

        data = run("kearns.from_orientations")["data"]
        assert isinstance(data["phase"], str)
        assert isinstance(data["phase_spec"], dict)


class TestTheInputsAreCheckedWhereTheyAreWrong:
    def test_a_cubic_phase_is_refused_by_name(self) -> None:
        """The parameter is defined on a basal pole; a cubic phase has none."""

        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_orientations", phase={"builtin": "ni_fcc"})
        assert "hexagonal" in str(excinfo.value)

    def test_a_reflection_line_with_the_wrong_column_count_names_its_line(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_diffractogram", reflections="0 0 0 2 1.0\n1 0 -1 3 2.0 1.0 9.0\n")
        assert "Line 2" in str(excinfo.value)

    def test_miller_bravais_indices_must_close(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_diffractogram", reflections="0 0 0 2  1.0  1.0\n1 1 -1 3  2.0  1.0\n")
        assert "h + k + i = 0" in str(excinfo.value)

    def test_one_reflection_is_refused_with_the_reason(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_diffractogram", reflections="0 0 0 2  1.0  1.0\n")
        assert "at least two reflections" in str(excinfo.value)

    def test_a_zero_random_intensity_is_refused_before_it_divides(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_diffractogram", reflections="0 0 0 2  1.0  0.0\n1 0 -1 0  1.0  1.0\n")
        assert "strictly positive" in str(excinfo.value)

    def test_unequally_spaced_profile_nodes_are_refused(self) -> None:
        """The common bin width is what cancels; unequal nodes silently break that."""

        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_tilt_profile", profile="5 1.0\n15 0.8\n40 0.2\n")
        assert "equally spaced" in str(excinfo.value)

    def test_a_tilt_outside_the_quadrant_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_tilt_profile", profile="5 1.0\n50 0.8\n95 0.2\n")
        assert "[0, 90]" in str(excinfo.value)

    def test_comments_and_blank_lines_are_ignored_not_parsed(self) -> None:
        result = run(
            "kearns.from_tilt_profile",
            profile="# a comment\n\n5 1.0\n\n45 0.5  # trailing note\n85 0.1\n",
        )
        assert len(result["table"]["rows"]) == 3

    def test_an_all_comment_block_says_so_rather_than_failing_obscurely(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_tilt_profile", profile="# nothing here\n\n")
        assert "No data rows" in str(excinfo.value)
        # The hint is where the reason lives: the user is looking at a box that
        # is visibly not empty, so "no rows" alone would read as a bug.
        assert "blank or a comment" in (excinfo.value.hint or "")

    def test_pasted_euler_angles_replace_the_model(self) -> None:
        """The EBSD route: three orientations in, a triad out, closing as always."""

        result = run(
            "kearns.from_orientations",
            euler_angles="0 0 0\n0 90 0\n0 45 0\n",
        )
        assert result["data"]["triad_sum"] == pytest.approx(1.0, abs=1e-9)
        assert result["inputs"]["source"] == "pasted"

    def test_a_malformed_euler_line_names_its_line(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_orientations", euler_angles="0 0 0\n10 20\n")
        assert "Line 2" in str(excinfo.value)


class TestTheExamplesRunAsAdvertised:
    @pytest.mark.parametrize(
        "example",
        [item for item in REGISTRY.examples() if item.panel == "kearns"],
        ids=lambda example: example.id,
    )
    def test_every_example_runs(self, example: Any) -> None:
        request = defaults(example.operation)
        request.update(example.request)
        result = REGISTRY.call(example.operation, request)
        assert result["summary"]

    def test_the_examples_describe_one_specimen(self) -> None:
        """All three examples must land on the same f, or the lesson is false."""

        values = []
        for example in (item for item in REGISTRY.examples() if item.panel == "kearns"):
            request = defaults(example.operation)
            request.update(example.request)
            values.append(value_of(REGISTRY.call(example.operation, request), "ND"))
        assert max(values) - min(values) < ROUTE_TOLERANCE


class TestTheUploadRoutes:
    """The two routes that read a measurement, on the repository's XRDML fixture.

    The fixture is a *random standard* truncated at 60 degrees of tilt, and that
    truncation is what makes it valuable here. A uniform pole figure integrated
    over the whole hemisphere gives exactly 1/3 in every direction; integrated
    over the measured cap alone it does not, because the cap is biased towards
    the pole. So the fixture is a ready-made demonstration of this route's real
    systematic error, and of the fact that the triad sum cannot detect it.
    """

    FIXTURE = "fixtures/xrdml/synthetic_random_standard.xrdml"

    def files(self, count: int = 1) -> dict[str, Any]:
        from pathlib import Path

        text = Path(self.FIXTURE).read_text(encoding="utf-8")
        return {"items": [{"name": f"zr-{index}.xrdml", "text": text} for index in range(count)]}

    def test_the_pole_figure_route_integrates_the_opened_figure(self) -> None:
        result = run(
            "kearns.from_pole_figure", files=self.files(), poles=((0, 0, 0, 2),)
        )
        assert result["data"]["diagnostics"]["sampled_point_count"] == 12.0
        assert value_of(result, "RD") == pytest.approx(value_of(result, "TD"), abs=1e-9)

    def test_closure_passes_while_the_answer_is_wrong(self) -> None:
        """The reason the panel no longer presents closure as a passed check.

        This figure is uniform, so the true f is 1/3 in every direction. Read
        over the measured cap alone it comes out at 0.518 — wrong by more than
        half the value — and the triad still sums to 1.0000, because the three
        numbers are the diagonal of one unit-trace tensor and always will be.
        """

        result = run("kearns.from_pole_figure", files=self.files(), poles=((0, 0, 0, 2),))
        assert result["data"]["triad_sum"] == pytest.approx(1.0, abs=1e-9)
        assert value_of(result, "ND") == pytest.approx(0.518, abs=0.005)
        assert value_of(result, "ND") != pytest.approx(1.0 / 3.0, abs=0.1)

    def test_the_incomplete_figure_is_reported_as_such(self) -> None:
        result = run("kearns.from_pole_figure", files=self.files(), poles=((0, 0, 0, 2),))
        assert result["data"]["diagnostics"]["max_polar_deg"] == pytest.approx(60.0, abs=0.1)
        assert result["data"]["diagnostics"]["measured_solid_angle_fraction"] == pytest.approx(
            0.5, abs=0.01
        )
        assert "covers 50% of the hemisphere" in result["summary"]
        assert "not detectable from the triad sum" in result["summary"]

    def test_the_summary_says_closure_is_by_construction(self) -> None:
        """Every single-tensor route must say so, or the check misleads."""

        for operation in ("kearns.from_orientations",):
            assert "closure by construction" in run(operation)["summary"]

    def test_a_non_basal_pole_is_not_called_a_kearns_parameter(self) -> None:
        result = run("kearns.from_pole_figure", files=self.files(), poles=((1, 0, -1, 1),))
        assert any("not the Kearns parameter" in note for note in result["notes"])

    def test_a_basal_pole_carries_no_such_warning(self) -> None:
        result = run("kearns.from_pole_figure", files=self.files(), poles=((0, 0, 0, 2),))
        assert not any("not the Kearns parameter" in note for note in result["notes"])

    def test_the_pole_figure_route_refuses_more_than_one_figure(self) -> None:
        """One figure is what it integrates; several is the ODF route's question."""

        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_pole_figure", files=self.files(3), poles=((0, 0, 0, 2),))
        assert "integrates one figure" in str(excinfo.value)

    def test_no_opened_file_is_reported_as_such(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_pole_figure", files={"items": []}, poles=((0, 0, 0, 2),))
        assert "No pole-figure file has been opened" in str(excinfo.value)

    def test_an_unreadable_file_names_the_file(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run(
                "kearns.from_pole_figure",
                files={"items": [{"name": "broken.xrdml", "text": "<not>xrdml</not>"}]},
                poles=((0, 0, 0, 2),),
            )
        assert "broken.xrdml" in str(excinfo.value)

    def test_miller_bravais_indices_must_close_on_the_upload_routes_too(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_pole_figure", files=self.files(), poles=((1, 1, -1, 1),))
        assert "h + k + i = 0" in str(excinfo.value)

    def test_the_odf_route_inverts_and_resolves(self) -> None:
        result = run(
            "kearns.from_odf",
            files=self.files(3),
            poles=((0, 0, 0, 2), (1, 0, -1, 1), (1, 0, -1, 2)),
            dictionary_count=300,
        )
        assert result["data"]["pole_figure_count"] == 3
        assert result["data"]["residual"] >= 0.0
        assert result["data"]["triad_sum"] == pytest.approx(1.0, abs=1e-9)
        assert "closes by construction" in result["summary"]

    def test_the_support_tensor_is_more_anisotropic_than_the_density_tensor(self) -> None:
        """The kernel shrinks anisotropy, so deconvolving it must increase the spread.

        This is the whole content of the choice the operation exposes, and it is
        a closed-form relation rather than a numerical accident: the two tensors
        differ by A_density = I/3 + beta (A_support - I/3) with beta < 1.
        """

        common = {
            "files": self.files(3),
            "poles": ((0, 0, 0, 2), (1, 0, -1, 1), (1, 0, -1, 2)),
            "dictionary_count": 300,
        }
        density = run("kearns.from_odf", tensor="density", **common)
        support = run("kearns.from_odf", tensor="support", **common)
        assert 0.0 < density["data"]["kernel_shrinkage"] < 1.0

        def spread(result: dict[str, Any]) -> float:
            values = [value_of(result, label) for label in ("RD", "TD", "ND")]
            return max(values) - min(values)

        assert spread(support) > spread(density)

    def test_the_odf_route_refuses_a_cubic_phase(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run("kearns.from_odf", files=self.files(3), phase={"builtin": "ni_fcc"})
        assert "hexagonal" in str(excinfo.value)
