from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import pytest

from pytex.core import (
    FrameDomain,
    Handedness,
    Phase,
    ReferenceFrame,
    get_phase_fixture,
    list_phase_fixtures,
    phase_fixture_catalog_path,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "fixtures/phases/catalog.json"
AUDIT_SUMMARY_PATH = REPO_ROOT / "benchmarks/structure_import/phase_fixture_audit_summary.json"


def _catalog() -> dict:
    return json.loads(CATALOG_PATH.read_text(encoding="utf-8"))


def _fixture_entries() -> list[dict]:
    return list(_catalog()["fixtures"])


def _metadata(entry: dict) -> dict:
    metadata_path = REPO_ROOT / entry["metadata_path"]
    return json.loads(metadata_path.read_text(encoding="utf-8"))


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _audit_summary() -> dict:
    return json.loads(AUDIT_SUMMARY_PATH.read_text(encoding="utf-8"))


def _require_pymatgen() -> None:
    pytest.importorskip(
        "pymatgen.core",
        reason="CIF-backed phase fixture loading requires the optional pymatgen dependency.",
    )


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["fixture_id"])
def test_phase_fixture_metadata_is_complete(entry: dict) -> None:
    metadata = _metadata(entry)
    required_fields = {
        "fixture_id",
        "display_name",
        "phase_name",
        "chemical_formula",
        "source_family",
        "source_record_id",
        "source_url",
        "authoring_mode",
        "redistribution",
        "citation",
        "crystal_system",
        "space_group_symbol",
        "space_group_number",
        "point_group",
        "lattice_parameters_angstrom",
        "lattice_angles_deg",
        "expected_conventional_cell_site_count",
        "expected_primitive_site_count",
        "intended_uses",
    }
    assert required_fields <= set(metadata)
    assert metadata["fixture_id"] == entry["fixture_id"]
    assert metadata["source_family"] == "Crystallography Open Database"
    assert metadata["source_url"].startswith("https://www.crystallography.net/cod/")
    assert "CC0" in metadata["redistribution"]
    assert metadata["intended_uses"]
    assert "structure_import_validation" in metadata["intended_uses"]


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["fixture_id"])
def test_phase_fixture_catalog_is_hash_pinned(entry: dict) -> None:
    artifact_path = REPO_ROOT / entry["artifact_path"]
    metadata_path = REPO_ROOT / entry["metadata_path"]
    assert _sha256(artifact_path) == entry["artifact_sha256"]
    assert _sha256(metadata_path) == entry["metadata_sha256"]


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["fixture_id"])
def test_phase_fixture_cif_round_trips_into_phase(entry: dict) -> None:
    _require_pymatgen()
    metadata = _metadata(entry)
    fixture_path = REPO_ROOT / entry["artifact_path"]
    crystal = _crystal_frame()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        warnings.filterwarnings("ignore", message="dict interface is deprecated.*")
        conventional = Phase.from_cif(fixture_path, crystal_frame=crystal)
        primitive = Phase.from_cif(fixture_path, crystal_frame=crystal, primitive=True)

    assert conventional.chemical_formula == metadata["chemical_formula"]
    assert conventional.space_group is not None
    assert conventional.space_group_symbol == metadata["space_group_symbol"]
    assert conventional.space_group_number == metadata["space_group_number"]
    assert conventional.symmetry.point_group == metadata["point_group"]
    assert conventional.unit_cell is not None
    assert len(conventional.unit_cell.sites) == metadata["expected_conventional_cell_site_count"]
    assert primitive.unit_cell is not None
    assert len(primitive.unit_cell.sites) == metadata["expected_primitive_site_count"]

    lengths = metadata["lattice_parameters_angstrom"]
    angles = metadata["lattice_angles_deg"]
    assert conventional.lattice.a == pytest.approx(lengths["a"], abs=1e-5)
    assert conventional.lattice.b == pytest.approx(lengths["b"], abs=1e-5)
    assert conventional.lattice.c == pytest.approx(lengths["c"], abs=1e-5)
    assert conventional.lattice.alpha_deg == pytest.approx(angles["alpha"], abs=1e-8)
    assert conventional.lattice.beta_deg == pytest.approx(angles["beta"], abs=1e-8)
    assert conventional.lattice.gamma_deg == pytest.approx(angles["gamma"], abs=1e-8)


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["fixture_id"])
def test_phase_fixture_cif_string_round_trips_into_phase(entry: dict) -> None:
    _require_pymatgen()
    metadata = _metadata(entry)
    cif_text = (REPO_ROOT / entry["artifact_path"]).read_text(encoding="utf-8")
    crystal = _crystal_frame()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        warnings.filterwarnings("ignore", message="dict interface is deprecated.*")
        conventional = Phase.from_cif_string(cif_text, crystal_frame=crystal)
        primitive = Phase.from_cif_string(cif_text, crystal_frame=crystal, primitive=True)

    assert conventional.chemical_formula == metadata["chemical_formula"]
    assert conventional.space_group_symbol == metadata["space_group_symbol"]
    assert conventional.space_group_number == metadata["space_group_number"]
    assert conventional.symmetry.point_group == metadata["point_group"]
    assert conventional.unit_cell is not None
    assert primitive.unit_cell is not None
    assert len(conventional.unit_cell.sites) == metadata["expected_conventional_cell_site_count"]
    assert len(primitive.unit_cell.sites) == metadata["expected_primitive_site_count"]


def test_phase_fixture_loader_api_exposes_catalog_and_loaders() -> None:
    _require_pymatgen()
    records = list_phase_fixtures()
    assert phase_fixture_catalog_path() == CATALOG_PATH
    assert {record.fixture_id for record in records} == {
        entry["fixture_id"] for entry in _fixture_entries()
    }

    record = get_phase_fixture("zr_hcp")
    assert record.display_name == "Zirconium (HCP)"
    assert record.chemical_formula == "Zr"
    crystal = _crystal_frame()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        warnings.filterwarnings("ignore", message="dict interface is deprecated.*")
        phase = record.load_phase(crystal_frame=crystal)

    assert phase.name == "zirconium-hcp"
    assert phase.space_group_symbol == "P6_3/mmc"
    assert phase.space_group_number == 194


def test_phase_fixtures_are_referenced_by_benchmark_or_validation_manifests() -> None:
    referenced_fixture_ids: set[str] = set()
    for manifest_path in sorted((REPO_ROOT / "benchmarks").glob("**/*.json")):
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        for fixture_id in payload.get("fixture_ids", []):
            referenced_fixture_ids.add(fixture_id)
        for artifact_path in payload.get("artifact_paths", []):
            assert (REPO_ROOT / artifact_path).exists(), artifact_path

    for entry in _fixture_entries():
        assert entry["fixture_id"] in referenced_fixture_ids


def test_structure_import_manifests_cover_full_phase_fixture_catalog() -> None:
    benchmark = json.loads(
        (REPO_ROOT / "benchmarks/structure_import/foundation_benchmark_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    validation = json.loads(
        (REPO_ROOT / "benchmarks/validation/structure_import_validation_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    expected_ids = sorted(entry["fixture_id"] for entry in _fixture_entries())
    assert sorted(benchmark["fixture_ids"]) == expected_ids
    assert sorted(validation["fixture_ids"]) == expected_ids


def test_structure_import_fixture_audit_summary_covers_full_phase_fixture_catalog() -> None:
    payload = _audit_summary()
    assert payload["schema_id"] == "pytex.structure_import_fixture_audit"
    assert payload["schema_version"] == "1.0.0"
    fixture_rows = payload["fixtures"]
    expected_ids = sorted(entry["fixture_id"] for entry in _fixture_entries())
    assert sorted(row["fixture_id"] for row in fixture_rows) == expected_ids
    for row in fixture_rows:
        assert (REPO_ROOT / row["artifact_path"]).exists()
        assert (REPO_ROOT / row["metadata_path"]).exists()
        assert set(row["conventional"]) >= {
            "space_group_symbol",
            "space_group_number",
            "point_group",
            "chemical_formula",
            "site_count",
            "lattice_parameters_angstrom",
            "lattice_angles_deg",
        }
        assert set(row["primitive"]) >= {
            "space_group_symbol",
            "space_group_number",
            "point_group",
            "chemical_formula",
            "site_count",
            "lattice_parameters_angstrom",
            "lattice_angles_deg",
        }


def test_structure_import_fixture_audit_summary_matches_loaded_phase_semantics() -> None:
    _require_pymatgen()
    fixture_rows = {row["fixture_id"]: row for row in _audit_summary()["fixtures"]}
    crystal = _crystal_frame()
    for entry in _fixture_entries():
        metadata = _metadata(entry)
        row = fixture_rows[entry["fixture_id"]]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
            warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
            warnings.filterwarnings("ignore", message="dict interface is deprecated.*")
            conventional = Phase.from_cif(REPO_ROOT / entry["artifact_path"], crystal_frame=crystal)
            primitive = Phase.from_cif(
                REPO_ROOT / entry["artifact_path"],
                crystal_frame=crystal,
                primitive=True,
            )

        assert row["conventional"]["space_group_symbol"] == conventional.space_group_symbol
        assert row["conventional"]["space_group_number"] == conventional.space_group_number
        assert row["conventional"]["point_group"] == conventional.symmetry.point_group
        assert row["conventional"]["chemical_formula"] == metadata["chemical_formula"]
        assert (
            row["conventional"]["site_count"] == metadata["expected_conventional_cell_site_count"]
        )
        assert row["primitive"]["space_group_symbol"] == primitive.space_group_symbol
        assert row["primitive"]["space_group_number"] == primitive.space_group_number
        assert row["primitive"]["point_group"] == primitive.symmetry.point_group
        assert row["primitive"]["chemical_formula"] == metadata["chemical_formula"]
        assert row["primitive"]["site_count"] == metadata["expected_primitive_site_count"]
