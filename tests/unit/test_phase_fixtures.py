from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

from pytex.core import FrameDomain, Handedness, Phase, ReferenceFrame

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "fixtures/phases/catalog.json"


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


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["fixture_id"])
def test_phase_fixture_cif_round_trips_into_phase(entry: dict) -> None:
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
