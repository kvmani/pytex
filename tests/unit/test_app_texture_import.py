"""Measured pole figures: reading XRDML files, and the scale they are read on.

The panel's model textures are checkable because the answer is known before the
calculation runs. A measurement is not, so the claims worth testing here are
different ones: that the file is read, that several figures land on *one* scale
so they can be compared, that the contour levels are the ones asked for, and
that the normalisation means what it says.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.services.texture import _contour_levels

FIXTURE = Path("fixtures/xrdml/synthetic_random_standard.xrdml")


def measured(names_and_poles: list[tuple[str, list[int]]], **overrides: object) -> dict:
    text = FIXTURE.read_text(encoding="utf-8")
    request: dict[str, object] = {
        "files": {"items": [{"name": name, "text": text} for name, _ in names_and_poles]},
        "phase": {"builtin": "ni_fcc"},
        "poles": [pole for _, pole in names_and_poles],
    }
    request.update(overrides)
    return REGISTRY.call("texture.measured_pole_figures", request)


def test_an_xrdml_file_becomes_a_readable_pole_figure() -> None:
    result = measured([("random.xrdml", [1, 1, 1])])
    figure = result["data"]["figures"][0]

    assert figure["file"] == "random.xrdml"
    # A measured figure collects the whole family, so it is written {111}.
    assert figure["label"] == "{111}"
    assert figure["count"] == len(figure["points"])
    assert figure["count"] > 0
    # Every point lands inside the projection disc.
    for point in figure["points"]:
        assert (point["x"] ** 2 + point["y"] ** 2) <= 1.0 + 1e-9
        assert 0.0 <= point["polar_deg"] <= 90.0 + 1e-9


def test_the_random_standard_reads_near_one_mrd() -> None:
    """The fixture is a random standard, so m.r.d. must put it near 1.

    This is the one thing a measurement *can* be checked against: a specimen
    with no texture reads 1 m.r.d. everywhere, and a normalisation that does not
    reproduce that is not a normalisation.
    """

    figure = measured([("random.xrdml", [1, 1, 1])])["data"]["figures"][0]
    assert figure["mean"] == pytest.approx(1.0, abs=0.25)


def test_several_figures_share_one_scale_so_they_can_be_compared() -> None:
    result = measured([("a.xrdml", [1, 1, 1]), ("b.xrdml", [2, 0, 0])])
    figures = result["data"]["figures"]
    assert [figure["label"] for figure in figures] == ["{111}", "{200}"]
    assert figures[0]["scale"] == figures[1]["scale"]
    assert figures[0]["scale"] == result["data"]["scale"]

    apart = measured(
        [("a.xrdml", [1, 1, 1]), ("b.xrdml", [2, 0, 0])],
        shared_scale=False,
    )["data"]["figures"]
    assert apart[0]["scale"]["maximum"] == apart[0]["maximum"]


def test_the_contour_levels_are_the_ones_asked_for() -> None:
    result = measured([("a.xrdml", [1, 1, 1])], contour_levels="1, 2, 4, 7, 10")
    assert result["data"]["levels"] == [1.0, 2.0, 4.0, 7.0, 10.0]
    # Whitespace separation, duplicates and disorder are all accepted.
    assert _contour_levels("4 1 2 2", count=3, minimum=0.0, maximum=5.0) == [1.0, 2.0, 4.0]


def test_automatic_levels_span_the_data_without_touching_its_ends() -> None:
    """A contour at the maximum draws a dot; one at the minimum draws the rim."""

    levels = _contour_levels("", count=4, minimum=1.0, maximum=6.0)
    assert len(levels) == 4
    assert levels[0] > 1.0
    assert levels[-1] < 6.0
    assert levels == sorted(levels)


def test_a_level_that_is_not_a_number_is_refused_rather_than_dropped() -> None:
    """Dropping it silently would draw fewer contours with no way to notice."""

    with pytest.raises(InvalidInputError) as raised:
        measured([("a.xrdml", [1, 1, 1])], contour_levels="1, 2, strong")
    assert raised.value.details["field"] == "contour_levels"
    assert "strong" in str(raised.value)


def test_the_normalisation_choice_changes_what_the_numbers_mean() -> None:
    peak = measured([("a.xrdml", [1, 1, 1])], intensity_normalization="max")
    assert peak["data"]["unit"] == "peak = 1"
    assert peak["data"]["figures"][0]["maximum"] == pytest.approx(1.0)

    raw = measured([("a.xrdml", [1, 1, 1])], intensity_normalization="none")
    assert raw["data"]["unit"] == "counts"
    assert raw["data"]["figures"][0]["maximum"] > 1.0


def test_the_result_says_the_plane_came_from_the_file_order() -> None:
    """XRDML records the diffraction angle, not the reflection.

    Assigning the wrong plane to a file is the mistake this operation makes
    easiest, so the result must say where the assignment came from.
    """

    notes = " ".join(measured([("a.xrdml", [1, 1, 1])])["notes"])
    assert "order the files were opened" in notes
    assert "Defocusing" in notes


def test_opening_nothing_is_refused_beside_its_own_control() -> None:
    with pytest.raises(InvalidInputError) as raised:
        REGISTRY.call(
            "texture.measured_pole_figures",
            {"files": {"items": []}, "phase": {"builtin": "ni_fcc"}},
        )
    assert raised.value.details["field"] == "files"


def test_a_file_that_is_not_a_pole_figure_is_refused_with_a_reason() -> None:
    with pytest.raises(InvalidInputError) as raised:
        REGISTRY.call(
            "texture.measured_pole_figures",
            {
                "files": {"items": [{"name": "scan.xrdml", "text": "<xrdMeasurements/>"}]},
                "phase": {"builtin": "ni_fcc"},
            },
        )
    assert raised.value.details["field"] == "files"
    assert "could not be read" in str(raised.value)


def test_fewer_planes_than_files_reuses_the_last_one() -> None:
    """Three files and one plane is a common way to open a set by mistake.

    Reusing the last plane is the recoverable behaviour: the figures draw, they
    are visibly labelled with the same plane, and the user fixes the list.
    """

    text = FIXTURE.read_text(encoding="utf-8")
    result = REGISTRY.call(
        "texture.measured_pole_figures",
        {
            "files": {"items": [{"name": f"{index}.xrdml", "text": text} for index in range(3)]},
            "phase": {"builtin": "ni_fcc"},
            "poles": [[1, 1, 1]],
        },
    )
    assert [figure["label"] for figure in result["data"]["figures"]] == ["{111}"] * 3


def test_the_odf_is_reconstructed_from_the_opened_figures() -> None:
    """The inversion runs, is sliced, and reports what it cost.

    The claim is deliberately modest, because the inversion is ill-posed: the
    sections exist, they are on the m.r.d. scale, and the residual travels with
    them. A test asserting a *particular* ODF from three copies of one random
    standard would be asserting the regularization.
    """

    result = measured(
        [("a.xrdml", [1, 1, 1]), ("b.xrdml", [2, 0, 0]), ("c.xrdml", [2, 2, 0])],
        reconstruct_odf=True,
        dictionary_count=300,
    )
    odf = result["data"]["odf"]

    assert odf["pole_figure_count"] == 3
    assert odf["dictionary_count"] == 300
    # The three sections every fcc texture paper prints.
    assert [section["phi2_deg"] for section in odf["sections"]] == [0.0, 45.0, 65.0]
    for section in odf["sections"]:
        assert len(section["densities"]) == len(section["phi1_deg"])
        assert len(section["densities"][0]) == len(section["big_phi_deg"])
        assert all(value >= 0.0 for row in section["densities"] for value in row)
    assert odf["max_mrd"] > 0.0
    assert odf["residual"] >= 0.0
    # The summary and the notes both say the inversion is ill-posed, because a
    # number offered without that caveat would be read as a measurement.
    assert "residual" in result["summary"]
    assert any("ill-posed" in note for note in result["notes"])


def test_no_odf_is_computed_unless_it_is_asked_for() -> None:
    """It is the expensive branch, and most sessions only want the figures."""

    assert measured([("a.xrdml", [1, 1, 1])])["data"]["odf"] is None


def test_a_plate_of_figures_is_labelled_by_something_that_distinguishes_them() -> None:
    """A label repeated on every panel identifies nothing.

    The sample name inside the file is the better identifier when it varies, and
    the file name is the only one that always does. Every fixture here is the
    same file, so the fallback is exactly the case under test.
    """

    figures = measured(
        [("cold-rolled.xrdml", [1, 1, 1]), ("annealed.xrdml", [2, 0, 0])]
    )["data"]["figures"]
    labels = [figure["sample_label"] for figure in figures]
    assert labels == ["cold-rolled", "annealed"]
    assert len(set(labels)) == len(labels)


@pytest.fixture(scope="module")
def harmonic_reconstruction() -> dict:
    """One harmonic inversion of the fixture, with ghost correction requested."""

    return measured(
        [("a.xrdml", [1, 1, 1])],
        reconstruct_odf=True,
        odf_method="harmonic",
        odf_bandlimit=6,
        ghost_correction="positivity",
    )


def test_the_harmonic_route_says_what_it_solved_for(harmonic_reconstruction: dict) -> None:
    odf = harmonic_reconstruction["data"]["odf"]
    assert odf["method"] == "harmonic"
    assert "harmonic series to degree 6" in odf["method_label"]
    assert odf["coefficient_count"] > 0
    assert odf["observation_count"] > 0


def test_an_under_determined_harmonic_fit_says_so(harmonic_reconstruction: dict) -> None:
    """Fewer measured intensities than coefficients means the regularization decides.

    Stated rather than blocked: an under-determined fit is legitimate when the
    smoothing is doing the work knowingly, and misleading when it is not. The
    reader is the one who can tell the difference, so the reader is told.
    """

    odf = harmonic_reconstruction["data"]["odf"]
    assert odf["coefficient_count"] > odf["observation_count"]
    notes = " ".join(harmonic_reconstruction["notes"])
    assert "fewer data than unknowns" in notes
    assert "regularization, not the specimen" in notes


def test_ghost_correction_on_a_cubic_material_reports_that_it_could_not_act(
    harmonic_reconstruction: dict,
) -> None:
    """The rotation group 432 has no odd-degree invariant below degree 9.

    So a cubic ODF expanded to degree 6 has no ghost part to correct, and the
    application must say that rather than report a correction of size zero as
    though a correction had been made.
    """

    ghost = harmonic_reconstruction["data"]["odf"]["ghost"]
    assert ghost is not None
    assert ghost["odd_basis_size"] == 0
    assert ghost["amplitude_ratio"] == 0.0
    assert "no odd-degree harmonic term" in harmonic_reconstruction["summary"]
    assert "first odd invariant is at degree 9" in " ".join(harmonic_reconstruction["notes"])


def test_the_dictionary_route_has_no_odd_part_to_correct() -> None:
    """Ghost correction is defined on the expansion, not on a cloud of weights."""

    odf = measured(
        [("a.xrdml", [1, 1, 1])],
        reconstruct_odf=True,
        dictionary_count=120,
    )["data"]["odf"]
    assert odf["method"] == "dictionary"
    assert odf["ghost"] is None
    assert "dictionary of 120 orientations" in odf["method_label"]
