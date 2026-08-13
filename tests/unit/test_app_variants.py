"""The variant panel, checked against counts and angles the literature fixes.

Nothing here is compared against a previously recorded output of this code. The
variant counts (24 for Kurdjumov-Sachs, 12 for Nishiyama-Wassermann and Burgers,
3 for Bain), the packet structure (4 packets of 6 on the parent {111}; 6 of 2 on
the parent {110}) and the ten distinct intervariant disorientations of
Kurdjumov-Sachs are all published values, listed in the citations the operations
carry. The projection geometry is checked against exact identities: a pole on
the projection axis lands at the origin, one on the equator lands on the unit
circle, and the two projections agree there and nowhere else.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

AUSTENITE = {"builtin": "austenite_fcc"}
FERRITE = {"builtin": "fe_bcc"}
ZIRCONIUM = {"builtin": "zr_hcp"}

#: The ten distinct Kurdjumov-Sachs intervariant disorientations, in degrees.
#:
#: Morito et al., Acta Mater. 51 (2003) 1789, Table 2 — the reference table for
#: lath martensite variant analysis. Quoted there to two decimals, so the
#: tolerance below is set by the quotation, not by the arithmetic.
MORITO_ANGLES_DEG = (
    10.53,
    14.88,
    20.61,
    21.06,
    47.11,
    49.47,
    50.51,
    51.73,
    57.21,
    60.00,
)

#: The disorientations available *within* a Kurdjumov-Sachs packet.
#:
#: The six variants sharing a parent {111} meet at only three of the ten angles
#: (Morito Table 2, the V1-V2 through V1-V6 rows). This is the sharper claim:
#: the whole spectrum being discrete is one thing, the within-packet subset
#: being three values out of ten is what makes a packet identifiable.
MORITO_WITHIN_PACKET_DEG = (10.53, 49.47, 60.00)


def call(operation: str, **request: object) -> dict:
    return REGISTRY.call(operation, request)


def pole_figure(**overrides: object) -> dict:
    request: dict[str, object] = {
        "phase": AUSTENITE,
        "child_phase": FERRITE,
        "relationship": "kurdjumov_sachs",
        "pole": [1, 0, 0],
        "packet_plane": [1, 1, 1],
        "projection": "stereographic",
        "include_parent": True,
    }
    request.update(overrides)
    return call("variants.pole_figure", **request)


class TestVariantCounts:
    """The multiplicities every textbook quotes."""

    @pytest.mark.parametrize(
        ("relationship", "parent", "child", "expected"),
        [
            ("kurdjumov_sachs", AUSTENITE, FERRITE, 24),
            ("nishiyama_wassermann", AUSTENITE, FERRITE, 12),
            ("bain", AUSTENITE, FERRITE, 3),
            ("burgers", FERRITE, ZIRCONIUM, 12),
        ],
    )
    def test_the_variant_count_is_the_published_one(
        self, relationship: str, parent: dict, child: dict, expected: int
    ) -> None:
        result = pole_figure(relationship=relationship, phase=parent, child_phase=child)
        assert result["data"]["variant_count"] == expected


class TestPackets:
    """Packet structure: the grouping a micrograph shows as a block."""

    def test_kurdjumov_sachs_gives_four_packets_of_six(self) -> None:
        """Morito et al.: 24 variants, 4 packets, 6 variants each.

        Four because the parent {111} family has four members, and each variant
        carries exactly one of them into parallelism with a child {110}.
        """

        data = pole_figure()["data"]
        assert data["packet_count"] == 4
        assert sorted(data["packet_sizes"].values()) == [6, 6, 6, 6]

    def test_nishiyama_wassermann_keeps_four_packets_with_three_each(self) -> None:
        """Halving the variants halves the packet size, not the packet count.

        Packets are counted by the parent family, which has four members
        whatever relationship sits on it.
        """

        data = pole_figure(relationship="nishiyama_wassermann")["data"]
        assert data["packet_count"] == 4
        assert sorted(data["packet_sizes"].values()) == [3, 3, 3, 3]

    def test_burgers_gives_six_packets_of_two_on_the_parent_110(self) -> None:
        data = pole_figure(
            relationship="burgers",
            phase=FERRITE,
            child_phase=ZIRCONIUM,
            pole=[0, 0, 1],
            packet_plane=[1, 1, 0],
        )["data"]
        assert data["packet_count"] == 6
        assert sorted(data["packet_sizes"].values()) == [2] * 6

    def test_an_uneven_grouping_says_so_rather_than_colouring_silently(self) -> None:
        """Pitsch is defined on {100}, so grouping it by {111} is not a packet.

        The grouping is still computed — it is the nearest-parallel parent plane
        of each variant, which is a real answer to a real question — but a
        figure whose colours carry no meaning must not pretend otherwise.
        """

        result = pole_figure(relationship="pitsch")
        assert len({*result["data"]["packet_sizes"].values()}) > 1
        assert "not the family" in result["notes"][0]

    def test_an_even_grouping_says_that_too(self) -> None:
        result = pole_figure()
        assert "Every packet holds 6 variants" in result["notes"][0]


class TestProjectionGeometry:
    """The projection is checked against identities, not against itself."""

    def test_a_pole_on_the_axis_lands_at_the_centre(self) -> None:
        """The parent (001) is the projection axis, so it projects to the origin.

        True of both projections, and the one point at which they agree
        trivially. It is included because it is the check that catches an
        inverted or transposed rotation immediately.
        """

        for method in ("stereographic", "equal_area"):
            parent = pole_figure(projection=method)["data"]["parent_poles"]
            centre = [row for row in parent if row["polar_deg"] < 1e-9]
            assert len(centre) == 1, f"{method}: expected exactly one pole on the axis"
            assert math.hypot(centre[0]["x"], centre[0]["y"]) == pytest.approx(0.0, abs=1e-12)

    def test_an_equatorial_pole_lands_on_the_unit_circle(self) -> None:
        """A pole at 90 degrees is on the rim, at radius 1 in either projection.

        The two projections differ everywhere between the centre and the rim —
        stereographic goes as tan(theta/2), equal-area as sqrt(2)sin(theta/2) —
        and coincide at both ends only once each is divided by its own
        equatorial radius. That division is the service's doing: the library
        returns each projection at its natural scale, where equal-area reaches
        √2, and drawing the two in the same circle without rescaling puts 41% of
        the equal-area figure outside the rim.
        """

        for method in ("stereographic", "equal_area"):
            parent = pole_figure(projection=method)["data"]["parent_poles"]
            rim = [row for row in parent if abs(row["polar_deg"] - 90.0) < 1e-6]
            assert rim, f"{method}: the (100) family must reach the equator"
            for row in rim:
                assert math.hypot(row["x"], row["y"]) == pytest.approx(1.0, abs=1e-9)

    def test_the_two_projections_differ_where_they_must(self) -> None:
        """Away from the centre and the rim, the radii must not agree.

        A projection control that silently did nothing would pass every other
        test in this class.
        """

        radii = {}
        for method in ("stereographic", "equal_area"):
            rows = pole_figure(projection=method)["data"]["poles"]
            middle = [row for row in rows if 20.0 < row["polar_deg"] < 70.0]
            assert middle
            radii[method] = math.hypot(middle[0]["x"], middle[0]["y"])
        assert radii["stereographic"] != pytest.approx(radii["equal_area"], abs=1e-6)

    def test_every_pole_lies_inside_the_disc(self) -> None:
        for method in ("stereographic", "equal_area"):
            result = pole_figure(projection=method)
            for row in result["table"]["rows"]:
                assert math.hypot(row["x"], row["y"]) <= 1.0 + 1e-9

    def test_the_reported_angles_describe_the_plotted_point(self) -> None:
        """Polar angle and radius must agree, or the hover card lies.

        The table is what the CSV writes and what a hover shows, so a row whose
        stated polar angle does not match where the point was drawn is a defect
        no plot can reveal. Checked through the stereographic closed form,
        r = tan(theta/2).
        """

        for row in pole_figure(projection="stereographic")["table"]["rows"]:
            expected = math.tan(math.radians(row["polar_deg"]) / 2.0)
            assert math.hypot(row["x"], row["y"]) == pytest.approx(expected, abs=1e-9)


class TestBainIsTheReferenceCase:
    """Bain has an answer that can be written down without computing it."""

    def test_every_bain_pole_lies_on_a_parent_100_or_110(self) -> None:
        """The Bain correspondence puts [100] of the child along [110] of the parent.

        That is the whole content of the Bain path: the child cell is the parent
        cell described on axes turned 45 degrees about one cube axis, with no
        rotation beyond what the correspondence itself requires. So each child
        {100} pole must lie exactly along a parent ⟨100⟩ or a parent ⟨110⟩ — at
        0, 45 or 90 degrees from the projection axis and at a multiple of 45
        degrees in azimuth, with nothing in between.

        The claim is exact, not approximate, which makes it the sharpest
        available check on the whole mapping chain: a transposed rotation, a
        swapped variant index or a projection applied in the wrong frame all
        move a pole off these angles.
        """

        rows = pole_figure(relationship="bain")["data"]["poles"]
        assert len(rows) == 3 * 3
        for row in rows:
            assert row["polar_deg"] == pytest.approx(
                min((0.0, 45.0, 90.0), key=lambda value: abs(value - row["polar_deg"])),
                abs=1e-9,
            ), f"{row['pole']} of variant {row['variant']} is not on a ⟨100⟩ or ⟨110⟩"
            # Azimuth is undefined at the pole itself, where the point is the
            # origin and any angle describes it.
            if row["polar_deg"] > 1e-9:
                # Distance to the nearest multiple of 45, not the remainder:
                # an azimuth a hair under 360 has remainder 45, not 0.
                offset = row["azimuth_deg"] % 45.0
                assert min(offset, 45.0 - offset) == pytest.approx(0.0, abs=1e-9)

    def test_the_parent_axis_is_shared_by_exactly_one_bain_variant(self) -> None:
        """Each Bain variant keeps one cube axis; only one of the three keeps [001]."""

        rows = pole_figure(relationship="bain")["data"]["poles"]
        on_axis = [row for row in rows if row["polar_deg"] < 1e-9]
        assert len(on_axis) == 1
        assert on_axis[0]["variant"] == 1


class TestIntervariantSpectrum:
    """Checked against Morito et al., Table 2."""

    @staticmethod
    def spectrum(**overrides: object) -> dict:
        request: dict[str, object] = {
            "phase": AUSTENITE,
            "child_phase": FERRITE,
            "relationship": "kurdjumov_sachs",
            "packet_plane": [1, 1, 1],
            "merge_equal_angles": True,
        }
        request.update(overrides)
        return call("variants.intervariant_misorientations", **request)

    def test_the_pair_count_is_every_unordered_pair(self) -> None:
        data = self.spectrum()["data"]
        assert data["pair_count"] == 24 * 23 // 2

    def test_the_distinct_angles_are_the_published_ten(self) -> None:
        rows = self.spectrum()["table"]["rows"]
        angles = [float(row["angle_deg"]) for row in rows]
        assert len(angles) == len(MORITO_ANGLES_DEG)
        assert np.allclose(angles, MORITO_ANGLES_DEG, atol=0.01)

    def test_the_multiplicities_account_for_every_pair(self) -> None:
        rows = self.spectrum()["table"]["rows"]
        assert sum(int(row["pairs"]) for row in rows) == 24 * 23 // 2

    def test_within_a_packet_only_three_of_the_ten_angles_occur(self) -> None:
        """The sharper claim, and the one that makes a packet identifiable."""

        rows = self.spectrum()["table"]["rows"]
        within = [float(row["angle_deg"]) for row in rows if int(row["same_packet"])]
        assert np.allclose(within, MORITO_WITHIN_PACKET_DEG, atol=0.01)

    def test_the_same_packet_pairs_are_exactly_the_within_packet_combinations(self) -> None:
        """Four packets of six give 4 * C(6,2) = 60 pairs, and no more."""

        rows = self.spectrum()["table"]["rows"]
        assert sum(int(row["same_packet"]) for row in rows) == 4 * (6 * 5 // 2)

    def test_the_unmerged_table_lists_every_pair_once(self) -> None:
        result = self.spectrum(merge_equal_angles=False)
        rows = result["table"]["rows"]
        assert len(rows) == 24 * 23 // 2
        pairs = {(row["variant_a"], row["variant_b"]) for row in rows}
        assert len(pairs) == len(rows)
        assert all(row["variant_a"] < row["variant_b"] for row in rows)

    def test_the_pairs_are_sorted_by_angle(self) -> None:
        rows = self.spectrum(merge_equal_angles=False)["table"]["rows"]
        angles = [float(row["angle_deg"]) for row in rows]
        assert angles == sorted(angles)

    def test_every_axis_is_a_unit_vector(self) -> None:
        rows = self.spectrum(merge_equal_angles=False)["table"]["rows"]
        for row in rows:
            norm = math.sqrt(row["axis_x"] ** 2 + row["axis_y"] ** 2 + row["axis_z"] ** 2)
            assert norm == pytest.approx(1.0, abs=1e-9)

    def test_the_sixty_degree_pairs_are_exactly_about_111_or_110(self) -> None:
        """The 60° rows must land on a rational axis with no deviation at all.

        Morito Table 2 gives the V1-V2 pair as 60° about [11-1]. That is the
        Σ3 twin relation of the bcc child, and it is exact — not a low-index
        approximation to an irrational axis, which most of the other pairs are.
        So this is the one row set where a zero deviation is a *requirement*,
        and it checks the axis labelling and the symmetry reduction together:
        a wrong basis, a missed symmetry operator or a sign error all show up
        here as a nonzero residual or a different axis.
        """

        rows = self.spectrum(merge_equal_angles=False)["table"]["rows"]
        sixty = [row for row in rows if abs(float(row["angle_deg"]) - 60.0) < 1e-6]
        assert sixty, "the 60° pairs are the ones with an exactly rational axis"
        for row in sixty:
            assert float(row["axis_deviation_deg"]) == pytest.approx(0.0, abs=1e-6), (
                f"the 60° axis {row['axis']} is not exact"
            )
            assert row["axis"] in {"[111]", "[011]"}, row["axis"]

    def test_an_irrational_axis_is_labelled_with_its_deviation(self) -> None:
        """Most intervariant axes are irrational, and the table must admit it.

        A rational parent axis maps to an irrational child axis under a real
        orientation relationship. Labelling those with the nearest low-index
        direction is useful; labelling them *silently* would claim an exactness
        that is not there, so every label carries how far off it is.
        """

        rows = self.spectrum(merge_equal_angles=False)["table"]["rows"]
        deviations = [float(row["axis_deviation_deg"]) for row in rows]
        assert max(deviations) > 0.1, "no pair is irrational; the labelling is suspiciously exact"
        assert max(deviations) < 5.0, (
            "some axis is more than 5° from every low-index direction, which is a label "
            "worth nothing; raise the index limit or report it as irrational"
        )

    def test_an_axis_and_its_reverse_are_not_listed_separately(self) -> None:
        """An axis is a line: [1 -1 1] and [-1 1 -1] must not both appear.

        The sign that comes out of a symmetry reduction is arbitrary, and
        without canonicalising it the same physical axis appears under two
        labels in one table.
        """

        rows = self.spectrum(merge_equal_angles=False)["table"]["rows"]
        for label in {row["axis"] for row in rows}:
            # Plain style writes a negative index with a minus sign, and
            # separates components with spaces whenever one is negative or
            # multi-digit — so the first character after the bracket is a minus
            # exactly when the first index is negative.
            assert not label.startswith("[-"), (
                f"{label} is the reverse of the canonical representative"
            )

    def test_burgers_has_its_own_spectrum(self) -> None:
        data = self.spectrum(
            phase=FERRITE,
            child_phase=ZIRCONIUM,
            relationship="burgers",
            packet_plane=[1, 1, 0],
        )["data"]
        assert data["variant_count"] == 12
        assert data["pair_count"] == 12 * 11 // 2


class TestErrorPaths:
    """What a user who asks for something impossible is told."""

    def test_an_inapplicable_relationship_names_its_field(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            pole_figure(relationship="kurdjumov_sachs", phase=ZIRCONIUM, child_phase=ZIRCONIUM)
        assert excinfo.value.details["field"] == "relationship"

    def test_bain_has_no_variant_pairs_worth_less_than_a_message(self) -> None:
        """Three variants do have pairs; the guard is for a one-variant case.

        Kept as a test of the message rather than of the guard, because the
        guard should not fire for any relationship the picker offers — if it
        ever does, that is worth knowing.
        """

        result = call(
            "variants.intervariant_misorientations",
            phase=AUSTENITE,
            child_phase=FERRITE,
            relationship="bain",
            packet_plane=[1, 1, 1],
            merge_equal_angles=True,
        )
        assert result["data"]["pair_count"] == 3

    def test_a_packet_plane_that_cannot_group_names_its_own_field(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            pole_figure(packet_plane=[0, 0, 0])
        assert excinfo.value.details["field"] in {"packet_plane", "pole"}


class TestPublicationFigure:
    """The published figure must be of the same poles the table exports."""

    @staticmethod
    def render(**overrides: object) -> dict:
        request: dict[str, object] = {
            "phase": AUSTENITE,
            "child_phase": FERRITE,
            "relationship": "kurdjumov_sachs",
            "pole": [1, 0, 0],
            "packet_plane": [1, 1, 1],
            "projection": "stereographic",
            "include_parent": True,
            "format": "svg",
            "dpi": 600,
        }
        request.update(overrides)
        return call("variants.render", **request)

    def test_an_svg_figure_is_produced_as_text(self) -> None:
        pytest.importorskip("matplotlib")
        result = self.render()
        assert result["data"]["encoding"] == "text"
        assert result["data"]["image"].lstrip().startswith(("<?xml", "<svg"))
        assert result["data"]["bytes"] > 0

    def test_a_png_figure_is_produced_as_base64(self) -> None:
        pytest.importorskip("matplotlib")
        import base64

        result = self.render(format="png")
        assert result["data"]["encoding"] == "base64"
        # The PNG magic number, so this is an image rather than an error page.
        assert base64.b64decode(result["data"]["image"])[:4] == b"\x89PNG"

    def test_the_figure_carries_every_pole_the_table_does(self) -> None:
        """One marker per row, or the figure is not of what was exported."""

        pytest.importorskip("matplotlib")
        table = pole_figure()
        rendered = self.render()
        assert rendered["data"]["pole_count"] == len(table["table"]["rows"])

    def test_no_figure_is_left_open(self) -> None:
        """A leaked figure is a defect here and a memory leak in the server."""

        matplotlib = pytest.importorskip("matplotlib")
        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        plt.close("all")
        self.render()
        assert plt.get_fignums() == []

    def test_the_packet_colours_are_the_same_in_both_renderers(self) -> None:
        """The screen and the figure must colour packet 3 the same.

        The colours exist twice — as hex in Python for matplotlib, as HSL in
        CSS for the browser — because neither renderer can read the other's
        constant. That duplication is only safe if it is checked, so the two
        lists are compared by converting the CSS ones and matching them within
        the rounding of an 8-bit channel.
        """

        import colorsys
        import re

        from pytex.app.server import STATIC_ROOT
        from pytex.app.services.variants import _PACKET_COLORS

        source = (STATIC_ROOT / "js" / "panels" / "variants.js").read_text(encoding="utf-8")
        block = re.search(r"const PACKET_COLORS = \[(.*?)\];", source, flags=re.S)
        assert block is not None, "PACKET_COLORS has been renamed; revisit this invariant"
        css = re.findall(r"hsl\(([\d.]+) ([\d.]+)% ([\d.]+)%\)", block.group(1))
        assert len(css) == len(_PACKET_COLORS), (
            "the two packet palettes have different lengths, so they cannot agree"
        )
        for (hue, saturation, lightness), expected in zip(css, _PACKET_COLORS, strict=True):
            red, green, blue = colorsys.hls_to_rgb(
                float(hue) / 360.0, float(lightness) / 100.0, float(saturation) / 100.0
            )
            actual = f"#{round(red * 255):02x}{round(green * 255):02x}{round(blue * 255):02x}"
            assert actual == expected, (
                f"the browser draws {actual} where the figure draws {expected}"
            )


class TestResultContract:
    """The panel and the export both read these, so both are checked."""

    def test_the_pole_table_and_the_plotted_data_are_the_same_rows(self) -> None:
        """A hover card and a CSV must not be able to disagree.

        The table is the concatenation of the variant poles and the parent
        poles, and the panel draws from `data`, so the two must have the same
        length and the same coordinates.
        """

        result = pole_figure()
        data = result["data"]
        assert len(result["table"]["rows"]) == len(data["poles"]) + len(data["parent_poles"])
        assert result["table"]["rows"][0]["x"] == data["poles"][0]["x"]

    def test_the_parent_poles_can_be_switched_off(self) -> None:
        assert pole_figure(include_parent=False)["data"]["parent_poles"] == []

    def test_the_declared_columns_match_the_table(self) -> None:
        result = pole_figure()
        declared = [column["key"] for column in result["data"]["columns"]]
        assert declared == [column["key"] for column in result["table"]["columns"]]

    def test_every_result_says_what_it_means(self) -> None:
        for result in (pole_figure(), TestIntervariantSpectrum.spectrum()):
            assert len(result["summary"]) > 60
            assert result["notes"]
            assert result["citations"]
