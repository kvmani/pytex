from __future__ import annotations

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
    "tests/unit/test_symmetry.py",
    "tests/unit/test_texture.py",
    "benchmarks/ebsd/foundation_benchmark_manifest.json",
    "benchmarks/diffraction/foundation_benchmark_manifest.json",
]


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    missing = [path for path in REQUIRED_PATHS if not (repo_root / path).exists()]
    if missing:
        for path in missing:
            print(f"MISSING: {path}")
        return 1

    print("Repository integrity check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
