"""TX3: reflection tables, figure and manifest exports for composite SAED.

Expected values are analytic identities of the diffraction geometry (d = 1/|g|,
r = camera constant x |g|, Friedel symmetry of a centrosymmetric structure) or
structural properties of the export contract — never a stored spot list.
"""

from __future__ import annotations

import csv
import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.adapters import composite_saed_manifest_schema_path
from pytex.core import (
    FrameDomain,
    Handedness,
    Lattice,
    OrientationRelationship,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.diffraction import (
    COMPOSITE_REFLECTION_TABLE_SCHEMA,
    COMPOSITE_SAED_MANIFEST_SCHEMA,
    REFLECTION_TABLE_COLUMNS,
    CompositeSAEDPattern,
    composite_reflection_table,
    composite_saed_manifest,
    export_composite_saed,
    phase_centering_is_declared,
    simulate_composite_saed,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases, make_fcc_bcc_phases


@pytest.fixture(scope="module")
def burgers_composite() -> CompositeSAEDPattern:
    parent, child = make_bcc_hcp_phases()
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=child
    )
    return simulate_composite_saed(
        relationship,
        ZoneAxis(np.array([1, 1, 0]), phase=parent),
        variant_indices=(1, 2, 3),
    )


class TestReflectionTableGeometry:
    def test_every_row_satisfies_the_bragg_and_camera_identities(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        """d = 1/|g|, and detector radius = camera constant x the in-plane part of g.

        Both are definitions, not fits, so they hold to machine precision. The
        radius uses the *projected* g: a reflection sits slightly off the exact
        Bragg condition, and that out-of-plane component is what the excitation
        error records, so `camera_constant * |g|` would be a small overestimate.
        """

        table = composite_reflection_table(burgers_composite)
        assert len(table) == burgers_composite.spot_count()
        camera_constant = table.camera_constant_mm_angstrom
        in_plane = np.vstack(
            [
                spots.g_detector_inv_angstrom
                for _, spots in burgers_composite.iter_spot_tables()
                if len(spots)
            ]
        )
        assert in_plane.shape[0] == len(table)
        for row, projected in zip(table.rows, in_plane, strict=True):
            assert_allclose(row.d_angstrom, 1.0 / row.g_inv_angstrom, rtol=1e-12)
            assert_allclose(
                row.detector_radius_mm,
                camera_constant * float(np.linalg.norm(projected)),
                rtol=1e-12,
            )
            # The projection can only shorten the vector.
            assert row.detector_radius_mm <= camera_constant * row.g_inv_angstrom + 1e-9

    def test_the_table_reproduces_the_engine_arrays_exactly(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        """The table must be a view of the simulation, not a parallel computation."""

        table = composite_reflection_table(burgers_composite)
        parent_rows = table.rows_for_source("parent")
        spots = burgers_composite.parent_spots
        assert spots is not None
        assert len(parent_rows) == len(spots)
        for position, row in enumerate(parent_rows):
            assert row.hkl == tuple(int(value) for value in spots.hkl[position])
            assert_allclose(row.detector_mm, spots.detector_mm[position], atol=0.0)
            assert row.relative_intensity == float(spots.intensity[position])

    def test_row_order_follows_parent_then_variants_in_composite_order(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        table = composite_reflection_table(burgers_composite)
        sources = [
            (row.source, row.variant_index)
            for row in table.rows
        ]
        # Every parent row precedes every variant row.
        first_variant = next(i for i, entry in enumerate(sources) if entry[0] == "variant")
        assert all(entry[0] == "parent" for entry in sources[:first_variant])
        assert table.variant_indices == burgers_composite.variant_indices

    def test_friedel_pairs_appear_with_equal_radius_and_intensity(self) -> None:
        """A centrosymmetric structure gives I(g) = I(-g); the table must show it.

        Friedel's law for a kinematic, centrosymmetric crystal without anomalous
        scattering — an identity of the theory, so any asymmetry in the exported
        table would be a defect in the export, not physics.
        """

        parent, child = make_fcc_bcc_phases()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        pattern = simulate_composite_saed(
            relationship, ZoneAxis(np.array([0, 0, 1]), phase=parent), variant_indices=(1,)
        )
        table = composite_reflection_table(pattern)
        by_hkl = {row.hkl: row for row in table.rows_for_source("parent")}
        assert by_hkl
        for hkl, row in by_hkl.items():
            opposite = (-hkl[0], -hkl[1], -hkl[2])
            assert opposite in by_hkl, f"Friedel partner of {hkl} missing"
            partner = by_hkl[opposite]
            assert_allclose(partner.d_angstrom, row.d_angstrom, rtol=1e-12)
            assert_allclose(
                partner.detector_radius_mm, row.detector_radius_mm, rtol=1e-9
            )
            assert_allclose(
                partner.relative_intensity, row.relative_intensity, rtol=1e-12
            )
            assert_allclose(partner.detector_mm, -np.asarray(row.detector_mm), atol=1e-9)


class TestIntensityThreshold:
    def test_thresholding_removes_only_weak_rows(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        full = composite_reflection_table(burgers_composite)
        filtered = composite_reflection_table(burgers_composite, intensity_threshold=0.5)
        assert len(filtered) < len(full)
        assert all(row.relative_intensity >= 0.5 for row in filtered.rows)
        assert "were excluded from this table" in filtered.describe()

    def test_the_pattern_itself_is_untouched(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        before = burgers_composite.spot_count()
        composite_reflection_table(burgers_composite, intensity_threshold=0.9)
        assert burgers_composite.spot_count() == before

    def test_rejects_an_out_of_range_threshold(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        with pytest.raises(ValueError, match=r"interval \[0, 1\)"):
            composite_reflection_table(burgers_composite, intensity_threshold=1.0)


class TestCenteringAudit:
    def test_declared_space_groups_are_reported_as_declared(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        audit = burgers_composite.centering_audit()
        assert audit == (
            ("beta-titanium", "I", True),
            ("alpha-titanium", "P", True),
        )
        assert "ASSUMED" not in burgers_composite.describe()

    def test_a_phase_without_a_space_group_is_flagged_as_assumed(self) -> None:
        """The silent failure this audit exists to catch.

        A body-centred phase supplied without a space-group symbol is simulated
        as primitive, keeping reflections its real symmetry forbids. Nothing in
        the pattern shows that, so the report has to.
        """

        frame = ReferenceFrame(
            "bare_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
        )
        lattice = Lattice(3.3, 3.3, 3.3, 90.0, 90.0, 90.0, crystal_frame=frame)
        bare = Phase(
            "undeclared",
            lattice=lattice,
            symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
            crystal_frame=frame,
        )
        assert not phase_centering_is_declared(bare)

        _, child = make_fcc_bcc_phases()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=bare, child_phase=child
        )
        pattern = simulate_composite_saed(
            relationship, ZoneAxis(np.array([0, 0, 1]), phase=bare), variant_indices=(1,)
        )
        audit = dict((name, declared) for name, _, declared in pattern.centering_audit())
        assert audit["undeclared"] is False
        text = pattern.describe()
        assert "ASSUMED" in text
        assert "forbidden reflections may be present" in text
        assert "WARNING" in composite_reflection_table(pattern).describe()

    def test_disabling_absences_is_stated_explicitly(self) -> None:
        from pytex.diffraction import KinematicSimulationConfig

        parent, child = make_fcc_bcc_phases()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        pattern = simulate_composite_saed(
            relationship,
            ZoneAxis(np.array([0, 0, 1]), phase=parent),
            variant_indices=(1,),
            config=KinematicSimulationConfig(apply_centering_absences=False),
        )
        assert "were NOT applied" in pattern.describe()


class TestExportFiles:
    def test_writes_table_figure_coincidences_and_manifest(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        export = export_composite_saed(burgers_composite, tmp_path, stem="burgers")
        assert export.reflection_table_path.exists()
        assert export.manifest_path.exists()
        assert export.coincidence_table_path is not None
        assert export.coincidence_table_path.exists()
        assert len(export.figure_paths) == 1
        assert export.figure_paths[0].suffix == ".svg"
        assert export.figure_paths[0].stat().st_size > 0
        assert all(path.exists() for path in export.paths())
        assert "burgers" in export.describe()

    def test_skipping_figures_avoids_importing_matplotlib_paths(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        export = export_composite_saed(
            burgers_composite, tmp_path, figure_formats=(), include_coincidences=False
        )
        assert export.figure_paths == ()
        assert export.coincidence_table_path is None
        assert export.reflection_table_path.exists()

    def test_rendering_does_not_leak_open_figures(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        import matplotlib.pyplot as plt

        plt.close("all")
        export_composite_saed(burgers_composite, tmp_path, figure_formats=("svg",))
        assert plt.get_fignums() == []

    def test_the_directory_is_created_when_missing(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        target = tmp_path / "nested" / "output"
        export = export_composite_saed(burgers_composite, target, figure_formats=())
        assert target.is_dir()
        assert export.directory == target


class TestCsvContract:
    def test_csv_columns_are_the_declared_contract(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        table = composite_reflection_table(burgers_composite)
        path = table.to_csv(tmp_path / "reflections.csv")
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            assert tuple(reader.fieldnames or ()) == REFLECTION_TABLE_COLUMNS
            records = list(reader)
        assert len(records) == len(table)

    def test_csv_round_trips_indices_and_positions(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        table = composite_reflection_table(burgers_composite)
        path = table.to_csv(tmp_path / "reflections.csv")
        with path.open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle))
        for record, row in zip(records, table.rows, strict=True):
            assert (int(record["h"]), int(record["k"]), int(record["l"])) == row.hkl
            assert float(record["detector_x_mm"]) == pytest.approx(
                row.detector_mm[0], abs=1e-9
            )
            assert float(record["detector_y_mm"]) == pytest.approx(
                row.detector_mm[1], abs=1e-9
            )

    def test_markdown_truncation_states_what_it_omitted(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        table = composite_reflection_table(burgers_composite)
        text = table.to_markdown(max_rows=5)
        assert text.count("\n") >= 7
        assert f"{len(table) - 5} further row(s) omitted" in text

    def test_json_payload_declares_its_schema_and_columns(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        table = composite_reflection_table(burgers_composite)
        payload = json.loads(json.dumps(table.to_json_dict()))
        assert payload["schema"] == COMPOSITE_REFLECTION_TABLE_SCHEMA
        assert tuple(payload["columns"]) == REFLECTION_TABLE_COLUMNS
        assert len(payload["rows"]) == len(table)


class TestManifest:
    def test_manifest_validates_against_its_schema(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        jsonschema = pytest.importorskip("jsonschema")
        export = export_composite_saed(burgers_composite, tmp_path, figure_formats=("svg",))
        payload = json.loads(export.manifest_path.read_text(encoding="utf-8"))
        schema = json.loads(
            composite_saed_manifest_schema_path().read_text(encoding="utf-8")
        )
        jsonschema.validate(payload, schema)
        assert payload["schema"] == COMPOSITE_SAED_MANIFEST_SCHEMA

    def test_manifest_records_every_variant_zone_axis_exactly_and_rationally(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        manifest = composite_saed_manifest(burgers_composite)
        assert len(manifest["variants"]) == len(burgers_composite.variant_patterns)
        for entry, pattern in zip(
            manifest["variants"], burgers_composite.variant_patterns, strict=True
        ):
            assert entry["variant_index"] == pattern.variant_index
            assert_allclose(
                entry["child_zone_axis_exact"],
                pattern.zone_axis_child.coordinates,
                atol=0.0,
            )
            assert entry["child_zone_axis_deviation_deg"] == pytest.approx(
                pattern.nearest_zone_axis.deviation_deg
            )

    def test_manifest_states_the_theory_level_and_normalization(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        manifest = composite_saed_manifest(burgers_composite)
        assert manifest["theory_level"] == "kinematic"
        assert manifest["intensity_normalization"] == "per_sub_pattern_max"

    def test_manifest_file_inventory_matches_what_was_written(
        self, burgers_composite: CompositeSAEDPattern, tmp_path
    ) -> None:
        export = export_composite_saed(burgers_composite, tmp_path, stem="case")
        payload = json.loads(export.manifest_path.read_text(encoding="utf-8"))
        listed = set(payload["files"].values())
        written = {path.name for path in export.paths()} - {export.manifest_path.name}
        assert listed == written
        for name in listed:
            assert (export.directory / name).exists()

    def test_manifest_reflection_counts_agree_with_the_pattern(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        manifest = composite_saed_manifest(burgers_composite)
        total = manifest["parent_reflection_count"] + sum(
            entry["reflection_count"] for entry in manifest["variants"]
        )
        assert total == manifest["total_reflection_count"]
        assert total == burgers_composite.spot_count()
