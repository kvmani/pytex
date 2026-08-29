"""The OR dossier, checked against the functions it is supposed to be calling.

The rule the dossier is built on is that it aggregates and never recomputes, so
these tests do not compare its numbers against recorded outputs. They compare
each one against the function a reader would check it against — the cell against
the lattice, the correspondence against `OrientationRelationship`, the spectrum
against `intervariant_misorientation_angles_deg` — plus the published values
that fix those functions in turn: 24 Kurdjumov-Sachs variants in 4 packets, the
42.85 deg disorientation, and Morito's ten intervariant angles.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.core.or_dossier import (
    OR_DOSSIER_SCHEMA_ID,
    OR_DOSSIER_SCHEMA_VERSION,
    ORDossierParallelism,
    or_dossier,
    or_dossier_schema_path,
)
from pytex.core.transformation import (
    OrientationRelationship,
    intervariant_misorientation_angles_deg,
)

#: Morito et al., Acta Mater. 51 (2003) 1789, Table 2 — the ten distinct
#: Kurdjumov-Sachs intervariant disorientations, in degrees.
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


def _cubic(name: str, a: float, species: str) -> Phase:
    frame = ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=frame)
    site = AtomicSite(label=f"{species}1", species=species, fractional_coordinates=np.zeros(3))
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=(site,)),
    )


def _hexagonal(name: str, a: float, c: float, species: str) -> Phase:
    frame = ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(a, a, c, 90.0, 90.0, 120.0, crystal_frame=frame)
    site = AtomicSite(label=f"{species}1", species=species, fractional_coordinates=np.zeros(3))
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=(site,)),
    )


def _ks() -> OrientationRelationship:
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=_cubic("austenite", 3.6, "Fe"),
        child_phase=_cubic("ferrite", 2.87, "Fe"),
    )


# --------------------------------------------------------------------------- #
# It aggregates; it does not recompute
# --------------------------------------------------------------------------- #


class TestAgreesWithItsSources:
    """Every block against the function it came from."""

    def test_the_lattice_block_is_the_lattice(self) -> None:
        relationship = _ks()
        dossier = or_dossier(relationship)
        for block, phase in (
            (dossier.parent, relationship.parent_phase),
            (dossier.child, relationship.child_phase),
        ):
            lattice = phase.lattice
            assert block.phase_name == phase.name
            assert block.cell_lengths_angstrom == pytest.approx(
                (lattice.a, lattice.b, lattice.c)
            )
            assert block.volume_angstrom3 == pytest.approx(lattice.volume_angstrom3())
            np.testing.assert_allclose(block.direct_basis, lattice.direct_basis().matrix)
            np.testing.assert_allclose(block.metric_tensor, lattice.metric_tensor())

    def test_the_cubic_volume_is_the_cube_of_the_edge(self) -> None:
        """An identity, so the shared volume definition cannot drift unnoticed."""

        dossier = or_dossier(_ks())
        assert dossier.parent.volume_angstrom3 == pytest.approx(3.6**3, rel=1e-12)
        assert dossier.child.volume_angstrom3 == pytest.approx(2.87**3, rel=1e-12)

    def test_the_correspondence_matrices_are_the_relationships_own(self) -> None:
        relationship = _ks()
        variants = relationship.generate_variants()
        dossier = or_dossier(relationship, variant=17)
        np.testing.assert_allclose(
            dossier.transformation.correspondence_direct,
            relationship.correspondence_direct(variant=variants[16]),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            dossier.transformation.correspondence_reciprocal,
            relationship.correspondence_reciprocal(variant=variants[16]),
            atol=1e-12,
        )

    def test_the_two_correspondences_are_inverse_transposes(self) -> None:
        """The property that preserves the zone law across the mapping."""

        dossier = or_dossier(_ks())
        direct = dossier.transformation.correspondence_direct
        reciprocal = dossier.transformation.correspondence_reciprocal
        np.testing.assert_allclose(direct.T @ reciprocal, np.eye(3), atol=1e-10)

    def test_the_spectrum_is_the_published_one(self) -> None:
        dossier = or_dossier(_ks())
        angles = dossier.misorientation.intervariant_angles_deg
        assert len(angles) == len(MORITO_ANGLES_DEG)
        np.testing.assert_allclose(sorted(angles), MORITO_ANGLES_DEG, atol=5e-3)

    def test_the_spectrum_is_the_functions_own(self) -> None:
        relationship = _ks()
        dossier = or_dossier(relationship)
        matrix = intervariant_misorientation_angles_deg(relationship)
        upper = matrix[np.triu_indices(matrix.shape[0], k=1)]
        assert set(dossier.misorientation.intervariant_angles_deg.tolist()) == set(
            np.unique(np.round(upper, 2)).tolist()
        )

    def test_the_disorientation_is_the_published_representative(self) -> None:
        dossier = or_dossier(_ks())
        assert dossier.misorientation.angle_deg == pytest.approx(42.85, abs=5e-3)
        axis = np.sort(np.abs(dossier.misorientation.axis))[::-1]
        np.testing.assert_allclose(axis, [0.9679, 0.1776, 0.1776], atol=5e-3)

    def test_the_packets_are_four_of_six(self) -> None:
        dossier = or_dossier(_ks())
        assert dossier.misorientation.variant_count == 24
        assert dossier.misorientation.packet_count == 4
        labels = dossier.misorientation.packet_labels
        assert sorted(np.bincount(labels)[1:].tolist()) == [6, 6, 6, 6]

    def test_the_packet_plane_defaults_to_the_defining_family(self) -> None:
        """A user should not have to know that Kurdjumov-Sachs packets on {111}."""

        dossier = or_dossier(_ks())
        assert dossier.misorientation.packet_plane == (1, 1, 1)


# --------------------------------------------------------------------------- #
# The parallelism block
# --------------------------------------------------------------------------- #


class TestParallelisms:
    def test_the_defining_pair_is_the_variants_own(self) -> None:
        """Not the relationship's nominal pair, which is only variant 1's."""

        relationship = _ks()
        planes = set()
        for index in range(1, 25):
            defining = or_dossier(relationship, variant=index).parallelism.defining
            plane = next(pair for pair in defining if pair.kind == "plane")
            assert plane.deviation_deg == 0.0
            planes.add(plane.parent_indices)
        assert len(planes) == 4

    def test_a_nominated_family_adds_rows_without_repeating_the_declared_one(self) -> None:
        relationship = _ks()
        parent = relationship.parent_phase
        dossier = or_dossier(
            relationship,
            planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        )
        assert dossier.parallelism.discovered
        declared = {
            (pair.kind, pair.parent_indices, pair.child_indices)
            for pair in dossier.parallelism.defining
        }
        for pair in dossier.parallelism.discovered:
            assert (pair.kind, pair.parent_indices, pair.child_indices) not in declared

    def test_nominated_directions_are_accepted(self) -> None:
        relationship = _ks()
        parent = relationship.parent_phase
        dossier = or_dossier(
            relationship,
            directions=[CrystalDirection([1.0, 1.0, 0.0], phase=parent)],
            tolerance_deg=3.0,
        )
        assert any(pair.kind == "direction" for pair in dossier.parallelism.discovered)
        assert all(pair.deviation_deg <= 3.0 for pair in dossier.parallelism.discovered)

    def test_labels_use_publication_notation(self) -> None:
        pair = ORDossierParallelism(
            kind="plane",
            origin="defining",
            parent_indices=(1, -1, 1),
            child_indices=(0, 1, 1),
            deviation_deg=0.0,
        )
        assert pair.parent_label == "(1 -1 1)"
        assert pair.child_label == "(011)"

    def test_a_pair_rejects_a_bad_kind_origin_or_deviation(self) -> None:
        for field, value in (
            ("kind", "bogus"),
            ("origin", "bogus"),
            ("deviation_deg", -1.0),
        ):
            payload = {
                "kind": "plane",
                "origin": "defining",
                "parent_indices": (1, 1, 1),
                "child_indices": (0, 1, 1),
                "deviation_deg": 0.0,
            }
            payload[field] = value
            with pytest.raises(ValueError):
                ORDossierParallelism(**payload)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# describe(), to_json() and the schema
# --------------------------------------------------------------------------- #


class TestExplainableAndSerializable:
    def test_describe_states_the_conventions_and_the_key_numbers(self) -> None:
        text = or_dossier(_ks(), variant=17).describe()
        assert "variant 17" in text
        assert "Angles are in degrees" in text
        assert "42.8" in text  # the disorientation
        assert "4 packets" in text
        assert "Morito" in text

    def test_describe_says_the_interface_was_not_analysed(self) -> None:
        """A missing section reads as an oversight; a stated absence does not."""

        text = or_dossier(_ks()).describe()
        assert "Not analysed" in text
        assert or_dossier(_ks()).to_json()["interface"] is None

    def test_describe_names_what_the_discovered_deviation_measures(self) -> None:
        relationship = _ks()
        parent = relationship.parent_phase
        text = or_dossier(
            relationship,
            planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        ).describe()
        assert "rationalization residual" in text

    def test_to_json_round_trips_through_json(self) -> None:
        payload = or_dossier(_ks(), variant=3).to_json()
        assert json.loads(json.dumps(payload)) == payload
        assert payload["schema_id"] == OR_DOSSIER_SCHEMA_ID
        assert payload["schema_version"] == OR_DOSSIER_SCHEMA_VERSION
        assert payload["variant_index"] == 3

    def test_the_schema_file_matches_what_to_json_writes(self) -> None:
        """The contract and the writer cannot drift apart silently."""

        schema = json.loads(or_dossier_schema_path().read_text(encoding="utf-8"))
        assert schema["$id"] == OR_DOSSIER_SCHEMA_ID
        assert schema["properties"]["schema_version"]["const"] == OR_DOSSIER_SCHEMA_VERSION
        payload = or_dossier(_ks()).to_json()
        assert set(schema["required"]) == set(payload)
        assert set(schema["properties"]) == set(payload)
        for key, block in (
            ("lattice_block", payload["parent"]),
            ("transformation_block", payload["transformation"]),
            ("misorientation_block", payload["misorientation"]),
            ("parallelism_block", payload["parallelism"]),
        ):
            assert set(schema["$defs"][key]["required"]) == set(block)
        pair_schema = schema["$defs"]["parallelism"]
        assert set(pair_schema["required"]) == set(payload["parallelism"]["pairs"][0])

    def test_the_reported_name_follows_the_relationship(self) -> None:
        relationship = _ks()
        assert or_dossier(relationship).relationship_name == relationship.name


# --------------------------------------------------------------------------- #
# export()
# --------------------------------------------------------------------------- #


class TestExport:
    def test_the_bundle_holds_numbers_tables_prose_and_figures(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        written = or_dossier(_ks(), variant=5).export(tmp_path)
        names = {path.name for path in written}
        assert names == {
            "or_dossier.json",
            "describe.md",
            "parallelisms.csv",
            "parallelisms.md",
            "intervariant_angles.csv",
            "or_stereogram.svg",
            "variant_contact_sheet.svg",
        }
        for path in written:
            assert path.exists()
            assert path.stat().st_size > 0

    def test_the_written_json_is_the_objects_json(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        dossier = or_dossier(_ks(), variant=5)
        dossier.export(tmp_path, figures=False)
        written = json.loads((tmp_path / "or_dossier.json").read_text(encoding="utf-8"))
        assert written == dossier.to_json()

    def test_the_csv_carries_every_parallelism_row(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        relationship = _ks()
        parent = relationship.parent_phase
        dossier = or_dossier(
            relationship,
            planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        )
        dossier.export(tmp_path, figures=False)
        lines = (tmp_path / "parallelisms.csv").read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) == len(dossier.parallelism.pairs) + 1  # header
        assert lines[0].split(",")[0] == "kind"

    def test_figures_can_be_skipped_without_losing_the_numbers(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        written = or_dossier(_ks()).export(tmp_path, figures=False)
        assert not any(path.suffix == ".svg" for path in written)
        assert (tmp_path / "or_dossier.json").exists()

    def test_export_leaves_no_open_figures(self, tmp_path) -> None:  # type: ignore[no-untyped-def]
        import matplotlib.pyplot as plt

        plt.close("all")
        or_dossier(_ks()).export(tmp_path)
        assert plt.get_fignums() == []


# --------------------------------------------------------------------------- #
# Other relationships, and the failure modes
# --------------------------------------------------------------------------- #


def test_burgers_dossier_reports_twelve_variants_in_six_packets() -> None:
    """The bcc-to-hcp path: the parent {110} family has six members, so the
    twelve variants fall into six packets of two rather than four of six."""

    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=_cubic("beta_zr", 3.57, "Zr"),
        child_phase=_hexagonal("alpha_zr", 3.23, 5.15, "Zr"),
    )
    dossier = or_dossier(relationship)
    assert dossier.misorientation.variant_count == 12
    assert dossier.misorientation.packet_count == 6
    assert sorted(np.bincount(dossier.misorientation.packet_labels)[1:].tolist()) == [2] * 6


def test_bain_dossier_reports_no_residual_polar_rotation() -> None:
    """Bain is the pure correspondence distortion, so its polar rotation is zero
    while a Kurdjumov-Sachs variant's is not — the classic distinction."""

    parent, child = _cubic("austenite", 3.6, "Fe"), _cubic("ferrite", 2.87, "Fe")
    bain = or_dossier(
        OrientationRelationship.from_bain_correspondence(parent_phase=parent, child_phase=child)
    )
    assert bain.transformation.polar_rotation_deg == pytest.approx(0.0, abs=1e-6)
    assert or_dossier(_ks()).transformation.polar_rotation_deg > 1.0


def test_an_out_of_range_variant_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match=r"one-based index in 1\.\.24"):
        or_dossier(_ks(), variant=25)


def test_a_relationship_without_a_plane_parallelism_groups_nothing() -> None:
    """Packets are counted on a plane family; without one, the dossier reports no
    grouping rather than inventing a family to group on."""

    relationship = _ks()
    bare = OrientationRelationship(
        name="bare",
        parent_phase=relationship.parent_phase,
        child_phase=relationship.child_phase,
        parent_to_child_rotation=relationship.parent_to_child_rotation,
    )
    dossier = or_dossier(bare)
    assert dossier.misorientation.packet_plane is None
    assert dossier.misorientation.packet_count == 0
    assert "no packet grouping was requested" in dossier.describe()
