"""TX2: variant-resolved parallel planes and directions.

Expected values are the crystallography of the named relationships themselves —
the packet structure of Kurdjumov-Sachs, the zone law, and the closure of the
forward and reverse index maps — not copied program output.
"""

from __future__ import annotations

import csv
import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalDirection,
    CrystalPlane,
    MillerIndex,
    OrientationRelationship,
    Phase,
    VariantCorrespondenceTable,
    variant_correspondence_table,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases
from tests.unit.test_transformation import make_phases


def _fcc_bcc() -> tuple[Phase, Phase]:
    _, _, parent, child = make_phases()
    return parent, child


def _plane(indices: tuple[int, int, int], phase: Phase) -> CrystalPlane:
    return CrystalPlane(MillerIndex(np.asarray(indices, dtype=np.int64), phase=phase), phase=phase)


def _ks() -> OrientationRelationship:
    parent, child = _fcc_bcc()
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )


class TestKurdjumovSachsPacketStructure:
    """The (111) table is the packet structure of lath martensite, in table form."""

    def test_six_of_twentyfour_variants_map_111_onto_011_exactly(self) -> None:
        """Each KS variant carries exactly one {111} member onto a {011} child plane.

        There are four {111} members and 24 variants, so any one member is the
        close-packed plane of exactly six of them — Morito's four packets of
        six. This is the defining crystallography, not a measured coincidence.
        """

        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        exact = table.exact_rows()
        assert len(table.rows) == 24
        assert len(exact) == 6
        for row in exact:
            assert sorted(np.abs(row.indices).tolist()) == [0, 1, 1]
            assert row.residual_deg < 1e-9

    def test_the_exact_variants_are_exactly_one_packet(self) -> None:
        """The six exact variants must be the packet `variant_close_packed_groups` finds."""

        from pytex.core import variant_close_packed_groups

        relationship = _ks()
        plane = _plane((1, 1, 1), relationship.parent_phase)
        table = variant_correspondence_table(relationship, plane)
        labels = variant_close_packed_groups(relationship, plane)
        variants = relationship.generate_variants()
        exact_indices = {row.variant_index for row in table.exact_rows()}
        # The packet label of every exactly-parallel variant must be the same,
        # and no variant outside the table's exact set may share it.
        packets = {
            int(labels[position])
            for position, variant in enumerate(variants)
            if variant.variant_index in exact_indices
        }
        assert len(packets) == 1
        packet = packets.pop()
        packet_members = {
            variant.variant_index
            for position, variant in enumerate(variants)
            if int(labels[position]) == packet
        }
        assert packet_members == exact_indices

    def test_grouping_reports_four_distinct_images(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        assert table.distinct_image_count((1, 1, 1)) == 4
        # Every group must be populated by the same number of variants: the
        # 24 variants split evenly because the parent symmetry acts transitively.
        counts: dict[int, int] = {}
        for row in table.rows:
            counts[row.equivalence_group] = counts.get(row.equivalence_group, 0) + 1
        assert set(counts.values()) == {6}

    def test_a_specific_close_packed_direction_is_exact_in_four_variants(self) -> None:
        """[1-10] lies in two of the four {111} planes, and each contributes two variants.

        KS pairs one <110> parent direction with one <111> child direction per
        variant, so a nominated direction is exactly parallel in
        2 planes x 2 variants = 4 of the 24.
        """

        relationship = _ks()
        direction = CrystalDirection([1.0, -1.0, 0.0], phase=relationship.parent_phase)
        table = variant_correspondence_table(relationship, direction)
        exact = table.exact_rows()
        assert len(exact) == 4
        for row in exact:
            assert sorted(np.abs(row.indices).tolist()) == [1, 1, 1]


class TestReverseSense:
    def test_every_variant_maps_the_child_close_packed_plane_back_onto_111(self) -> None:
        """(011) child is the close-packed plane of *every* KS variant.

        The forward direction is selective — one {111} member per variant — but
        the reverse is not: each variant's own {011} image comes from some {111}
        member, so mapping (011) back lands on a {111} plane in all 24.
        """

        relationship = _ks()
        table = variant_correspondence_table(
            relationship,
            _plane((0, 1, 1), relationship.child_phase),
            sense="child_to_parent",
        )
        assert len(table.rows) == 24
        assert len(table.exact_rows()) == 24
        assert table.distinct_image_count((0, 1, 1)) == 1
        for row in table.rows:
            assert sorted(np.abs(row.indices).tolist()) == [1, 1, 1]

    def test_forward_then_reverse_returns_the_source(self) -> None:
        """Round trip through one variant must recover the input indices."""

        relationship = _ks()
        variants = relationship.generate_variants()
        source = _plane((1, 1, 1), relationship.parent_phase)
        forward = variant_correspondence_table(relationship, source, variants=variants[:1])
        image = forward.rows[0]
        back = variant_correspondence_table(
            relationship,
            _plane(
                (int(image.indices[0]), int(image.indices[1]), int(image.indices[2])),
                relationship.child_phase,
            ),
            sense="child_to_parent",
            variants=variants[:1],
        )
        assert_allclose(back.rows[0].indices, source.miller.indices, atol=0)
        assert back.rows[0].residual_deg < 1e-9

    def test_reverse_sense_rejects_parent_phase_objects(self) -> None:
        relationship = _ks()
        with pytest.raises(ValueError, match="must belong to the child phase"):
            variant_correspondence_table(
                relationship,
                _plane((1, 1, 1), relationship.parent_phase),
                sense="child_to_parent",
            )


class TestHexagonalLabels:
    def test_burgers_maps_011_beta_onto_the_basal_plane_in_two_variants(self) -> None:
        """Burgers has six packets of two, so one {110} member is basal in exactly two."""

        parent, child = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=child
        )
        table = variant_correspondence_table(relationship, _plane((0, 1, 1), parent))
        assert len(table.rows) == 12
        exact = table.exact_rows()
        assert len(exact) == 2
        for row in exact:
            # Rationalization is sign-sensitive (it matches the exact image's
            # direction, which matters for a diffraction vector g), so the basal
            # plane may be reported as (0001) or its antiparallel (000-1).
            assert row.image_label.replace(" ", "") in {"(0001)", "(000-1)"}
            assert sorted(np.abs(row.indices).tolist()) == [0, 0, 1]

    def test_hexagonal_images_are_labeled_in_four_index_form(self) -> None:
        parent, child = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=child
        )
        table = variant_correspondence_table(
            relationship, CrystalDirection([1.0, -1.0, 1.0], phase=parent)
        )
        # Four-index labels carry three separators; three-index labels carry two
        # (or none, when every component is a single digit).
        assert all(row.image_label.count(" ") == 3 for row in table.exact_rows())
        # The parent side is cubic and stays three-index.
        assert table.rows[0].source_label.count(" ") == 2


class TestMultipleSources:
    def test_several_planes_are_tabulated_together_with_independent_grouping(self) -> None:
        relationship = _ks()
        parent = relationship.parent_phase
        table = variant_correspondence_table(
            relationship, [_plane((1, 1, 1), parent), _plane((1, 0, 0), parent)]
        )
        assert len(table.rows) == 48
        assert table.source_indices == ((1, 1, 1), (1, 0, 0))
        assert len(table.rows_for((1, 1, 1))) == 24
        assert len(table.rows_for_variant(1)) == 2
        # Grouping is per source object, so both start their numbering at 0.
        assert min(row.equivalence_group for row in table.rows_for((1, 1, 1))) == 0
        assert min(row.equivalence_group for row in table.rows_for((1, 0, 0))) == 0

    def test_mixing_planes_and_directions_is_rejected(self) -> None:
        relationship = _ks()
        parent = relationship.parent_phase
        with pytest.raises(ValueError, match="all objects to be of one kind"):
            variant_correspondence_table(
                relationship,
                [_plane((1, 1, 1), parent), CrystalDirection([1.0, 0.0, 0.0], phase=parent)],
            )

    def test_a_single_object_need_not_be_wrapped_in_a_sequence(self) -> None:
        relationship = _ks()
        single = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        wrapped = variant_correspondence_table(
            relationship, [_plane((1, 1, 1), relationship.parent_phase)]
        )
        assert [row.image_label for row in single.rows] == [
            row.image_label for row in wrapped.rows
        ]


class TestRationalization:
    def test_raising_the_index_bound_never_worsens_the_residual(self) -> None:
        """A larger candidate set contains the smaller one, so the best fit cannot get worse."""

        relationship = _ks()
        plane = _plane((1, 1, 1), relationship.parent_phase)
        coarse = variant_correspondence_table(relationship, plane, max_index=4)
        fine = variant_correspondence_table(relationship, plane, max_index=12)
        for low, high in zip(coarse.rows, fine.rows, strict=True):
            assert high.residual_deg <= low.residual_deg + 1e-12

    def test_exact_correspondences_are_independent_of_the_index_bound(self) -> None:
        relationship = _ks()
        plane = _plane((1, 1, 1), relationship.parent_phase)
        coarse = variant_correspondence_table(relationship, plane, max_index=3)
        fine = variant_correspondence_table(relationship, plane, max_index=17)
        assert {row.variant_index for row in coarse.exact_rows()} == {
            row.variant_index for row in fine.exact_rows()
        }

    def test_the_exact_components_and_indices_agree_where_the_residual_is_zero(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        for row in table.exact_rows():
            exact = row.exact_components / np.linalg.norm(row.exact_components)
            rational = row.indices.astype(float) / np.linalg.norm(row.indices.astype(float))
            assert_allclose(np.abs(exact @ rational), 1.0, atol=1e-9)


class TestExports:
    def test_csv_round_trips_the_indices_exactly(self, tmp_path) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        path = table.to_csv(tmp_path / "table.csv")
        with path.open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle))
        assert len(records) == len(table.rows)
        for record, row in zip(records, table.rows, strict=True):
            assert record["image_indices"] == " ".join(
                str(int(value)) for value in row.indices
            )
            assert float(record["residual_deg"]) == pytest.approx(row.residual_deg, abs=1e-12)
            assert int(record["variant"]) == row.variant_index

    def test_markdown_has_one_row_per_entry_plus_a_header(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        lines = table.to_markdown().strip().splitlines()
        assert len(lines) == len(table.rows) + 2
        assert lines[0].startswith("| variant |")

    def test_json_payload_is_serializable_and_matches_the_rows(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        payload = json.loads(json.dumps(table.to_json_dict()))
        assert payload["schema"] == "pytex.variant_correspondence_table/1"
        assert payload["kind"] == "plane"
        assert payload["sense"] == "parent_to_child"
        assert len(payload["rows"]) == len(table.rows)

    def test_describe_states_the_distinct_and_exact_counts(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        text = table.describe()
        assert "4 crystallographically distinct image(s)" in text
        assert "exactly parallel in 6 of them" in text
        assert "(011)" in text


class TestValidation:
    def test_rejects_an_empty_object_list(self) -> None:
        with pytest.raises(ValueError, match="at least one object"):
            variant_correspondence_table(_ks(), [])

    def test_rejects_an_unknown_sense(self) -> None:
        relationship = _ks()
        with pytest.raises(ValueError, match="sense must be"):
            variant_correspondence_table(
                relationship,
                _plane((1, 1, 1), relationship.parent_phase),
                sense="sideways",
            )

    def test_rejects_an_empty_table(self) -> None:
        with pytest.raises(ValueError, match="at least one row"):
            VariantCorrespondenceTable(
                relationship_name="x",
                kind="plane",
                sense="parent_to_child",
                source_phase_name="a",
                image_phase_name="b",
                rows=(),
                max_index=3,
            )

    def test_row_arrays_are_read_only(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        with pytest.raises(ValueError):
            table.rows[0].indices[0] = 9

    def test_unknown_source_raises_rather_than_returning_zero(self) -> None:
        relationship = _ks()
        table = variant_correspondence_table(
            relationship, _plane((1, 1, 1), relationship.parent_phase)
        )
        with pytest.raises(KeyError):
            table.distinct_image_count((2, 3, 5))
