from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

REQUIRED_PATHS = [
    "README.md",
    "mission.md",
    "specifications.md",
    "AGENTS.md",
    "LICENSE",
    "pyproject.toml",
    "src/pytex/__init__.py",
    "src/pytex/core/__init__.py",
    "src/pytex/core/fixtures.py",
    "docs/README.md",
    "docs/architecture/overview.md",
    "docs/architecture/canonical_data_model.md",
    "docs/architecture/orientation_and_texture_foundation.md",
    "docs/architecture/ebsd_foundation.md",
    "docs/architecture/diffraction_foundation.md",
    "docs/architecture/multimodal_characterization_foundation.md",
    "docs/architecture/phase_transformation_foundation.md",
    "docs/architecture/repo_review_2026_foundation_audit.md",
    "docs/testing/strategy.md",
    "docs/testing/mtex_parity_matrix.md",
    "docs/testing/diffraction_validation_matrix.md",
    "docs/testing/structure_validation_matrix.md",
    "docs/standards/notation_and_conventions.md",
    "docs/standards/documentation_architecture.md",
    "docs/standards/latex_and_figures.md",
    "docs/standards/scientific_citation_policy.md",
    "docs/standards/benchmark_and_tolerance_governance.md",
    "docs/standards/hexagonal_and_trigonal_conventions.md",
    "docs/standards/development_principles.md",
    "docs/standards/data_contracts_and_manifests.md",
    "docs/standards/reference_canon.md",
    "docs/tex/README.md",
    "docs/tex/theory/orientation_space_and_disorientation.tex",
    "docs/tex/theory/euler_convention_handling.tex",
    "docs/tex/theory/fundamental_region_reduction.tex",
    "docs/tex/algorithms/discrete_odf_and_pole_figures.tex",
    "docs/tex/algorithms/ebsd_kam_parameterization.tex",
    "docs/tex/algorithms/ebsd_local_misorientation.tex",
    "docs/tex/algorithms/ebsd_grain_segmentation_and_grod.tex",
    "docs/tex/algorithms/ebsd_boundaries_and_cleanup.tex",
    "docs/tex/theory/hexagonal_conventions.tex",
    "docs/figures/reference_frames.svg",
    "docs/figures/hcp_reference_frame.svg",
    "docs/figures/pole_figure_construction.svg",
    "docs/site/README.md",
    "docs/site/conf.py",
    "docs/site/index.md",
    "docs/site/workflows/ebsd_grains.md",
    "fixtures/mtex_parity/README.md",
    "fixtures/phases/README.md",
    "fixtures/phases/catalog.json",
    "fixtures/phases/fe_bcc/phase.cif",
    "fixtures/phases/fe_bcc/metadata.json",
    "fixtures/phases/zr_hcp/phase.cif",
    "fixtures/phases/zr_hcp/metadata.json",
    "fixtures/phases/ni_fcc/phase.cif",
    "fixtures/phases/ni_fcc/metadata.json",
    "fixtures/phases/nicl/phase.cif",
    "fixtures/phases/nicl/metadata.json",
    "fixtures/phases/diamond/phase.cif",
    "fixtures/phases/diamond/metadata.json",
    "fixtures/mtex_parity/rotation/euler_quaternion_cases.json",
    "fixtures/mtex_parity/fundamental_region/vector_cases.json",
    "fixtures/mtex_parity/ebsd/kam_cases.json",
    "tests/parity/test_mtex_rotation_parity.py",
    "tests/parity/test_mtex_fundamental_region_parity.py",
    "tests/parity/test_mtex_ebsd_parity.py",
    "tests/integration/test_cli_and_optional_adapters.py",
    "schemas/README.md",
    "tests/unit/test_frames.py",
    "tests/unit/test_hexagonal.py",
    "tests/unit/test_reference_policy.py",
    "tests/unit/test_manifests.py",
    "tests/unit/test_symmetry.py",
    "tests/unit/test_transformation.py",
    "tests/unit/test_texture.py",
    "benchmarks/ebsd/foundation_benchmark_manifest.json",
    "benchmarks/diffraction/foundation_benchmark_manifest.json",
    "benchmarks/structure_import/foundation_benchmark_manifest.json",
    "benchmarks/structure_import/foundation_workflow_result_manifest.json",
    "benchmarks/validation/README.md",
    "benchmarks/validation/diffraction_validation_manifest.json",
    "benchmarks/validation/structure_import_validation_manifest.json",
    "schemas/benchmark_manifest.schema.json",
    "schemas/experiment_manifest.schema.json",
    "schemas/validation_manifest.schema.json",
    "schemas/workflow_result_manifest.schema.json",
]


def _read_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must decode to a JSON object.")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _check_phase_fixture_catalog(repo_root: Path) -> list[str]:
    issues: list[str] = []
    catalog = _read_json(repo_root / "fixtures/phases/catalog.json")
    required_catalog_fields = {"catalog_id", "schema_version", "source_policy", "fixtures"}
    missing_catalog_fields = sorted(required_catalog_fields - set(catalog))
    if missing_catalog_fields:
        return [
            "INVALID: fixtures/phases/catalog.json missing fields: "
            + ", ".join(missing_catalog_fields)
        ]

    fixtures = catalog.get("fixtures")
    if not isinstance(fixtures, list) or not fixtures:
        return ["INVALID: fixtures/phases/catalog.json must contain a non-empty fixtures list."]

    required_metadata_fields = {
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
    seen_ids: set[str] = set()
    for entry in fixtures:
        if not isinstance(entry, dict):
            issues.append("INVALID: fixtures/phases/catalog.json contains a non-object fixture entry.")
            continue
        fixture_id = entry.get("fixture_id")
        artifact_path = entry.get("artifact_path")
        metadata_path = entry.get("metadata_path")
        artifact_sha256 = entry.get("artifact_sha256")
        metadata_sha256 = entry.get("metadata_sha256")
        if not isinstance(fixture_id, str) or not fixture_id:
            issues.append("INVALID: fixtures/phases/catalog.json fixture entries must define fixture_id.")
            continue
        if fixture_id in seen_ids:
            issues.append(f"INVALID: duplicate phase fixture id '{fixture_id}'.")
            continue
        seen_ids.add(fixture_id)
        artifact_full_path: Path | None = None
        if not isinstance(artifact_path, str) or not artifact_path:
            issues.append(f"INVALID: phase fixture '{fixture_id}' must define artifact_path.")
        else:
            artifact_full_path = repo_root / artifact_path
        if artifact_full_path is not None and not artifact_full_path.exists():
            issues.append(f"MISSING: {artifact_path}")
        elif artifact_full_path is not None:
            if not isinstance(artifact_sha256, str) or len(artifact_sha256) != 64:
                issues.append(
                    f"INVALID: phase fixture '{fixture_id}' must define a 64-character artifact_sha256."
                )
            elif _sha256(artifact_full_path) != artifact_sha256:
                issues.append(
                    f"INVALID: {artifact_path} sha256 does not match catalog entry for '{fixture_id}'."
                )
        if not isinstance(metadata_path, str) or not metadata_path:
            issues.append(f"INVALID: phase fixture '{fixture_id}' must define metadata_path.")
            continue
        metadata_full_path = repo_root / metadata_path
        if not metadata_full_path.exists():
            issues.append(f"MISSING: {metadata_path}")
            continue
        if not isinstance(metadata_sha256, str) or len(metadata_sha256) != 64:
            issues.append(
                f"INVALID: phase fixture '{fixture_id}' must define a 64-character metadata_sha256."
            )
        elif _sha256(metadata_full_path) != metadata_sha256:
            issues.append(
                f"INVALID: {metadata_path} sha256 does not match catalog entry for '{fixture_id}'."
            )
        metadata = _read_json(metadata_full_path)
        missing_metadata_fields = sorted(required_metadata_fields - set(metadata))
        if missing_metadata_fields:
            issues.append(
                f"INVALID: {metadata_path} missing fields: {', '.join(missing_metadata_fields)}"
            )
        if metadata.get("fixture_id") != fixture_id:
            issues.append(
                f"INVALID: {metadata_path} fixture_id does not match catalog entry '{fixture_id}'."
            )
        intended_uses = metadata.get("intended_uses")
        if not isinstance(intended_uses, list) or not intended_uses:
            issues.append(f"INVALID: {metadata_path} intended_uses must be a non-empty list.")
        elif "structure_import_validation" not in intended_uses:
            issues.append(
                f"INVALID: {metadata_path} must include 'structure_import_validation' in intended_uses."
            )
    return issues


def _check_manifest_artifacts_and_fixtures(repo_root: Path) -> list[str]:
    issues: list[str] = []
    catalog = _read_json(repo_root / "fixtures/phases/catalog.json")
    fixture_ids = {entry["fixture_id"] for entry in catalog["fixtures"]}
    for manifest_path in sorted(repo_root.glob("benchmarks/**/*.json")):
        payload = _read_json(manifest_path)
        relative_manifest_path = manifest_path.relative_to(repo_root)
        manifest_fixture_ids = payload.get("fixture_ids", [])
        if not isinstance(manifest_fixture_ids, list):
            issues.append(f"INVALID: {relative_manifest_path} fixture_ids must be a list.")
            manifest_fixture_ids = []
        for fixture_id in manifest_fixture_ids:
            if fixture_id not in fixture_ids and not (repo_root / fixture_id).exists():
                issues.append(
                    f"INVALID: {relative_manifest_path} references unknown fixture id '{fixture_id}'."
                )
        artifact_paths = payload.get("artifact_paths", [])
        if not isinstance(artifact_paths, list):
            issues.append(f"INVALID: {relative_manifest_path} artifact_paths must be a list.")
            continue
        for artifact_path in artifact_paths:
            if not isinstance(artifact_path, str) or not artifact_path:
                issues.append(f"INVALID: {relative_manifest_path} contains an empty artifact path.")
                continue
            if not (repo_root / artifact_path).exists():
                issues.append(f"MISSING: {artifact_path}")
    structure_fixture_ids = tuple(
        _read_json(repo_root / "benchmarks/structure_import/foundation_benchmark_manifest.json").get(
            "fixture_ids", []
        )
    )
    validation_fixture_ids = tuple(
        _read_json(repo_root / "benchmarks/validation/structure_import_validation_manifest.json").get(
            "fixture_ids", []
        )
    )
    expected_fixture_ids = tuple(sorted(fixture_ids))
    if tuple(sorted(structure_fixture_ids)) != expected_fixture_ids:
        issues.append(
            "INVALID: structure_import foundation benchmark manifest fixture_ids must match the pinned phase fixture catalog."
        )
    if tuple(sorted(validation_fixture_ids)) != expected_fixture_ids:
        issues.append(
            "INVALID: structure_import validation manifest fixture_ids must match the pinned phase fixture catalog."
        )
    return issues


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    issues = [path for path in REQUIRED_PATHS if not (repo_root / path).exists()]
    issues.extend(_check_phase_fixture_catalog(repo_root))
    issues.extend(_check_manifest_artifacts_and_fixtures(repo_root))
    if issues:
        for issue in issues:
            print(issue if issue.startswith(("MISSING:", "INVALID:")) else f"MISSING: {issue}")
        return 1

    print("Repository integrity check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
