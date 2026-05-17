from __future__ import annotations

import json
from pathlib import Path

from scripts.compare_parity_results import compare_result_roots
from scripts.generate_pytex_parity_campaign import generate_campaign

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ROOT = REPO_ROOT / "fixtures" / "mtex_parity" / "campaigns"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_parity_campaign_files_use_shared_contract() -> None:
    campaign_files = sorted(CAMPAIGN_ROOT.glob("*.json"))
    assert campaign_files
    for path in campaign_files:
        payload = _load_json(path)
        assert payload["schema_id"] == "pytex.parity_case_campaign"
        assert payload["schema_version"] == "0.1.0"
        assert payload["target_baseline"]["system"] == "mtex"
        assert payload["target_baseline"]["version"] == "6.0.0"
        assert payload["conventions"]["quaternion_order"] == "wxyz"
        assert payload["conventions"]["matrix_storage"] == "row_major"
        case_ids = [case["case_id"] for case in payload["cases"]]
        assert len(case_ids) == len(set(case_ids))
        for case in payload["cases"]:
            assert case["status"] in {"active", "pending"}
            assert case["operation_family"]
            assert case["operation"]
            assert case["requested_outputs"]
            assert case["tolerances"]
            if case["status"] == "pending":
                assert case["reason_pending"]


def test_pytex_generator_emits_results_for_active_and_pending_cases(tmp_path: Path) -> None:
    written = generate_campaign(
        CAMPAIGN_ROOT / "orientation_core_cases.json",
        tmp_path,
        created_utc="2000-01-01T00:00:00Z",
        include_pending=True,
    )
    assert len(written) == 6
    result = _load_json(tmp_path / "orientation_core_v1" / "cubic_bunge_euler_001.json")
    assert result["schema_id"] == "pytex.parity_result"
    assert result["case_status"] == "active"
    assert result["producer"]["system"] == "pytex"
    assert result["producer"]["created_utc"] == "2000-01-01T00:00:00Z"
    assert result["results"]["quaternion_wxyz"]
    assert result["results"]["rotation_matrix"]
    assert result["results"]["mapped_crystal_vectors"]


def test_pytex_generator_writes_pending_xrdml_cases_as_skipped(tmp_path: Path) -> None:
    written = generate_campaign(
        CAMPAIGN_ROOT / "xrdml_pole_figure_cases.json",
        tmp_path,
        created_utc="2000-01-01T00:00:00Z",
        include_pending=True,
    )
    assert len(written) == 2
    result = _load_json(tmp_path / "xrdml_pole_figure_v1" / "cubic_xrdml_pole_figure_pending.json")
    assert result["case_status"] == "skipped"
    assert result["results"]["reason_pending"]


def test_parity_comparator_passes_identical_pytex_results(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    generate_campaign(
        CAMPAIGN_ROOT / "ipf_color_cases.json",
        left,
        created_utc="2000-01-01T00:00:00Z",
        include_pending=True,
    )
    generate_campaign(
        CAMPAIGN_ROOT / "ipf_color_cases.json",
        right,
        created_utc="2000-01-01T00:00:00Z",
        include_pending=True,
    )
    assert compare_result_roots(left, right) == []

