from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_foundational_docs_declare_normative_and_informative_references() -> None:
    foundational_docs = [
        "mission.md",
        "specifications.md",
        "docs/architecture/overview.md",
        "docs/architecture/canonical_data_model.md",
        "docs/architecture/orientation_and_texture_foundation.md",
        "docs/architecture/ebsd_foundation.md",
        "docs/architecture/diffraction_foundation.md",
        "docs/architecture/multimodal_characterization_foundation.md",
        "docs/architecture/phase_transformation_foundation.md",
        "docs/architecture/repo_review_2026_foundation_audit.md",
        "docs/standards/engineering_governance.md",
        "docs/standards/notation_and_conventions.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/latex_and_figures.md",
        "docs/standards/terminology_and_symbol_registry.md",
        "docs/standards/scientific_citation_policy.md",
        "docs/standards/development_principles.md",
        "docs/standards/data_contracts_and_manifests.md",
        "docs/standards/reference_canon.md",
        "docs/testing/strategy.md",
        "docs/testing/automated_test_cases.md",
        "docs/testing/mtex_parity_matrix.md",
        "docs/testing/diffraction_validation_matrix.md",
        "docs/testing/structure_validation_matrix.md",
        "docs/testing/plotting_validation_matrix.md",
        "docs/roadmap/implementation_roadmap.md",
    ]
    for path in foundational_docs:
        content = _read(path)
        assert "### Normative" in content, path
        assert "### Informative" in content, path


def test_space_group_spec_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    canonical_data_model = _read("docs/architecture/canonical_data_model.md")
    assert "SpaceGroupSpec" in specifications
    assert "SpaceGroupSpec" in api_guide
    assert "SpaceGroupSpec" in canonical_data_model


def test_multimodal_primitives_are_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    core_model = _read("docs/site/concepts/core_model.md")
    for primitive in (
        "AcquisitionGeometry",
        "CalibrationRecord",
        "MeasurementQuality",
        "ScatteringSetup",
    ):
        assert primitive in specifications
        assert primitive in api_guide
        assert primitive in core_model


def test_texture_inversion_surface_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    canonical_data_model = _read("docs/architecture/canonical_data_model.md")
    assert "ODFInversionReport" in specifications
    assert "ODFInversionReport" in api_guide
    assert "ODFInversionReport" in canonical_data_model
    assert "HarmonicODF" in specifications
    assert "HarmonicODF" in api_guide
    assert "HarmonicODF" in canonical_data_model
    assert "HarmonicODFReconstructionReport" in specifications
    assert "HarmonicODFReconstructionReport" in api_guide
    assert "HarmonicODFReconstructionReport" in canonical_data_model


def test_plotting_surface_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    plotting_workflow = _read("docs/site/workflows/plotting_primitives.md")
    for symbol in (
        "plot_vector_set",
        "plot_wulff_net",
        "plot_crystal_directions",
        "plot_crystal_planes",
        "plot_symmetry_orbit",
        "plot_symmetry_elements",
        "plot_euler_set",
        "plot_quaternion_set",
        "plot_rotations",
        "plot_orientations",
        "plot_pole_figure",
        "plot_inverse_pole_figure",
        "plot_odf",
        "plot_xrd_pattern",
        "plot_saed_pattern",
        "plot_crystal_structure_3d",
        "resolve_style",
        "save_documentation_figure_svg",
    ):
        assert symbol in specifications
        assert symbol in api_guide
        assert symbol in plotting_workflow


def test_xrdml_texture_surface_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    texture_workflow = _read("docs/site/workflows/texture_odf_inversion.md")
    xrdml_workflow = _read("docs/site/workflows/xrdml_texture_import.md")
    for symbol in (
        "XRDMLPoleFigureMeasurement",
        "read_xrdml_pole_figure",
        "load_xrdml_pole_figure",
        "invert_xrdml_pole_figures",
        "LaboTexPoleFigureMeasurement",
        "read_labotex_pole_figures",
        "load_labotex_pole_figures",
        "invert_labotex_pole_figures",
    ):
        assert symbol in specifications
        assert symbol in api_guide
        assert symbol in xrdml_workflow or symbol in _read(
            "docs/site/workflows/labotex_open_pole_figures.md"
        )
    assert "ODF.evaluate_pole_density" in texture_workflow


def test_xrd_and_saed_surface_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    xrd_workflow = _read("docs/site/workflows/xrd_generation.md")
    saed_workflow = _read("docs/site/workflows/saed_generation.md")
    for symbol in (
        "RadiationSpec",
        "PowderReflection",
        "PowderPattern",
        "SAEDSpot",
        "SAEDPattern",
    ):
        assert symbol in specifications
        assert symbol in api_guide
    assert "AtomicSite" in specifications
    assert "AtomicSite" in api_guide
    assert "RadiationSpec" in xrd_workflow
    assert "PowderPattern" in xrd_workflow
    assert "SAEDPattern" in saed_workflow


def test_terminology_registry_and_glossary_are_present() -> None:
    standards_doc = _read("docs/standards/terminology_and_symbol_registry.md")
    glossary_page = _read("docs/site/concepts/technical_glossary_and_symbols.md")
    assert "symbol" in standards_doc.lower()
    assert "glossary" in glossary_page.lower()
    assert "zone axis" in glossary_page.lower()
    assert "2\\theta" in glossary_page or "2θ" in glossary_page


def test_transformation_primitives_are_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    api_guide = _read("docs/site/api/index.md")
    canonical_data_model = _read("docs/architecture/canonical_data_model.md")
    for primitive in (
        "OrientationRelationship",
        "TransformationVariant",
        "PhaseTransformationRecord",
    ):
        assert primitive in specifications
        assert primitive in api_guide
        assert primitive in canonical_data_model


def test_manifest_family_is_part_of_stable_documented_surface() -> None:
    specifications = _read("specifications.md")
    canonical_data_model = _read("docs/architecture/canonical_data_model.md")
    manifest_policy = _read("docs/standards/data_contracts_and_manifests.md")
    for primitive in (
        "ExperimentManifest",
        "BenchmarkManifest",
        "ValidationManifest",
        "WorkflowResultManifest",
    ):
        assert primitive in specifications
        assert primitive in canonical_data_model
        assert primitive in manifest_policy


def test_reference_working_corpus_assets_exist_and_are_linked_from_core_docs() -> None:
    reference_index = _read("references/reference_index.md")
    formulation_summary = _read("references/formulation_summary.md")
    feature_opportunities = _read("references/feature_opportunities.md")
    testing_strategy = _read("docs/testing/strategy.md")
    development_principles = _read("docs/standards/development_principles.md")
    reference_canon = _read("docs/standards/reference_canon.md")

    assert "Aspect" in reference_index
    assert "Future tasks should consult this document" in formulation_summary
    assert "High-Value Opportunities" in feature_opportunities
    assert "Automated Test Cases" in testing_strategy
    assert "formulation_summary" in development_principles
    assert "Repository Working Corpus" in reference_canon


def test_validation_site_exposes_automated_test_case_docs() -> None:
    validation_index = _read("docs/site/validation/index.md")
    site_case_page = _read("docs/site/validation/automated_test_cases.md")
    canonical_docs = _read("docs/site/reference/canonical_docs.md")

    assert "automated_test_cases" in validation_index
    assert "automated_test_cases.md" in site_case_page
    assert "validation/automated_test_cases" in canonical_docs
