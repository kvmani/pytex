from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_foundational_docs_agree_on_layered_documentation_policy() -> None:
    """Every foundational doc must name the same three documentation pillars.

    Sphinx is the browsable surface, the scientific notes are canonical prose
    and mathematics, and SVG is canonical for figures. The notes were LaTeX
    until 2026-08-11; they are now MyST under docs/site/theory/ so that the
    derivations render on the site instead of downloading as .tex sources.
    """
    foundational_docs = [
        "README.md",
        "mission.md",
        "specifications.md",
        "AGENTS.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/terminology_and_symbol_registry.md",
        "docs/standards/visualization_style_guide.md",
    ]
    for path in foundational_docs:
        content = _read(path).lower()
        assert "sphinx" in content, path
        assert "notes" in content, path
        assert "svg" in content, path


def test_no_foundational_doc_still_calls_latex_canonical() -> None:
    """The retired policy must not survive anywhere in the governing docs."""
    for path in (
        "README.md",
        "mission.md",
        "specifications.md",
        "AGENTS.md",
        "docs/README.md",
        "docs/standards/scientific_notes_and_figures.md",
        "docs/standards/documentation_architecture.md",
    ):
        content = _read(path).lower()
        assert "docs/tex/" not in content or "until 2026-08-11" in content, path
        assert "latex is the canonical" not in content, path
        assert "canonically in latex" not in content, path


def test_docs_site_placeholder_encodes_public_entry_point() -> None:
    content = _read("docs/site/README.md").lower()
    assert "public documentation entry point" in content
    assert "concept" in content
    assert "tutorial" in content
    assert "api" in content
    assert "notebook" in content


def test_notes_standard_points_back_to_documentation_architecture() -> None:
    content = _read("docs/standards/scientific_notes_and_figures.md").lower()
    assert "documentation_architecture.md" in content
    assert "visualization_style_guide.md" in content
    assert "sphinx" in content


def test_visualization_style_guide_defines_canonical_svg_policy() -> None:
    content = _read("docs/standards/visualization_style_guide.md").lower()
    for token in ("#07122f", "#2563eb", "#7c3aed", "mermaid", "title", "desc"):
        assert token in content


def test_foundational_docs_encode_notebook_policy() -> None:
    foundational_docs = [
        "mission.md",
        "specifications.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/development_principles.md",
    ]
    for path in foundational_docs:
        content = _read(path).lower()
        assert "notebook" in content, path


def test_foundational_docs_encode_dual_plotting_output_policy() -> None:
    policy_docs = [
        "specifications.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/scientific_notes_and_figures.md",
    ]
    for path in policy_docs:
        content = _read(path).lower()
        assert "matplotlib" in content, path
        assert "svg" in content, path


def test_installation_guide_covers_platforms_and_docs_builds() -> None:
    content = _read("docs/site/tutorials/installation_and_build.md").lower()
    for token in ("windows", "macos", "linux", "sphinx", "latex", "jupyter"):
        assert token in content


def test_site_index_references_glossary_and_installation_pages() -> None:
    content = _read("docs/site/index.md")
    assert "concepts/technical_glossary_and_symbols" in content
    assert "tutorials/installation_and_build" in content


def test_site_index_and_readme_expose_active_roadmap() -> None:
    site_index = _read("docs/site/index.md").lower()
    site_readme = _read("docs/site/README.md").lower()
    assert "roadmap/implementation_roadmap" in site_index
    assert "implementation roadmap" in site_readme


def test_governing_roadmaps_point_to_the_current_delivery_cycle() -> None:
    critical = _read("docs/roadmap/critical_review_and_development_guide.md")
    feature_review = _read("docs/roadmap/feature_capability_review_2026_08.md")
    world_class = _read("docs/roadmap/world_class_feature_roadmap.md")
    for content in (critical, feature_review, world_class):
        lowered = " ".join(content.lower().split())
        for tokens in (
            ("measured", "xrd"),
            ("defocus",),
            ("hex-grid", "ebsd"),
            ("finite-thickness", "saed"),
            ("component", "odf", "fitting"),
        ):
            assert all(token in lowered for token in tokens), tokens
    assert "804 tests pass" not in critical
    assert "i1" in world_class.lower()
    assert "i15 all landed" in world_class.lower()


def test_completed_application_ledger_has_no_open_status() -> None:
    ledger = _read("docs/development/active_task_progress.md")
    section = ledger.split(
        "## Application Platform: Desktop + Intranet Workbench — COMPLETE", maxsplit=1
    )[1].split("\n## ", maxsplit=1)[0]
    assert "| in progress |" not in section.lower()
    assert "all twelve steps are done" in section.lower()


def test_ci_covers_windows_and_ratchets_sphinx_warnings() -> None:
    workflow = _read(".github/workflows/ci.yml")
    strategy = _read("docs/testing/strategy.md")
    assert "windows-latest" in workflow
    assert "scripts/check_sphinx_warnings.py --max-warnings 0" in workflow
    assert "scripts/check_sphinx_warnings.py --max-warnings 0" in strategy


def test_ci_and_strategy_define_the_critical_browser_lane() -> None:
    workflow = _read(".github/workflows/ci.yml")
    strategy = _read("docs/testing/strategy.md")

    assert 'node-version: "22"' in workflow
    assert "npx playwright install --with-deps chromium" in workflow
    assert "npm run test:browser" in workflow
    assert "all seven workspaces" in strategy
    assert "test-only" in strategy


def test_executable_examples_standard_is_encoded_and_cross_linked() -> None:
    standard = _read("docs/standards/executable_examples.md").lower()
    for token in ("worked example", "computed", "reference value", "tolerance", "citation"):
        assert token in standard, token

    agents = _read("AGENTS.md")
    assert "docs/standards/executable_examples.md" in agents
    assert "worked_examples/" in agents

    doc_arch = _read("docs/standards/documentation_architecture.md").lower()
    assert "executable worked example" in doc_arch
    assert "docstring contract" in doc_arch

    site_index = _read("docs/site/index.md")
    assert "examples/index" in site_index

    standards_index = _read("docs/site/standards/index.md")
    assert "executable_examples" in standards_index


def test_foundational_docs_encode_executable_example_policy() -> None:
    for path in ("mission.md", "specifications.md", "AGENTS.md"):
        content = _read(path).lower()
        assert "worked example" in content, path


def test_interop_docs_state_validation_boundaries_explicitly() -> None:
    interop = _read("docs/site/workflows/orix_kikuchipy_interop.md").lower()
    ebsd = _read("docs/site/workflows/ebsd_import_normalization.md").lower()
    assert "what is executable today" in interop
    assert "what is not being claimed" in interop
    assert "current limits" in interop
    assert "what is verified today" in ebsd
    assert "interpretation rule" in ebsd


def test_workbench_guide_quotes_numbers_the_code_actually_produces() -> None:
    """The user guide's numbers are computed here, not trusted.

    `AGENTS.md`: "Documentation numbers must not be hand-transcribed." The
    workbench guide makes three quantitative claims a reader would act on — that
    Kurdjumov-Sachs gives 24 variants in 4 packets of 6, that the 276 variant
    pairs fall on ten specific disorientations, and that only three of those
    occur within a packet. Each is recomputed from the service layer and matched
    against the text, so a change in either that leaves them disagreeing is a
    failure rather than a silently wrong page.
    """

    from pytex.app import REGISTRY

    guide = _read("docs/site/workflows/workbench_application.md")

    pole_figure = REGISTRY.call(
        "variants.pole_figure",
        {
            "phase": {"builtin": "austenite_fcc"},
            "child_phase": {"builtin": "fe_bcc"},
            "relationship": "kurdjumov_sachs",
            "pole": [1, 0, 0],
            "packet_plane": [1, 1, 1],
            "projection": "stereographic",
            "include_parent": True,
        },
    )["data"]
    assert f"**{pole_figure['variant_count']} variants**" in guide
    assert f"**{pole_figure['packet_count']} groups of " in guide
    assert set(pole_figure["packet_sizes"].values()) == {6}, (
        "the guide says 4 groups of 6; the computation no longer agrees"
    )

    spectrum = REGISTRY.call(
        "variants.intervariant_misorientations",
        {
            "phase": {"builtin": "austenite_fcc"},
            "child_phase": {"builtin": "fe_bcc"},
            "relationship": "kurdjumov_sachs",
            "packet_plane": [1, 1, 1],
            "merge_equal_angles": True,
        },
    )
    rows = spectrum["table"]["rows"]
    assert f"The {spectrum['data']['pair_count']} variant pairs" in guide
    for row in rows:
        # The guide quotes each to two decimals, which is how the source table
        # quotes them.
        assert f"{float(row['angle_deg']):.2f}°" in guide, (
            f"the guide does not quote the {row['angle_deg']}° disorientation"
        )
    within = [row for row in rows if int(row["same_packet"])]
    assert len(within) == 3, "the guide says three of the ten occur within a packet"

    # And the m.r.d. claim the texture section rests on.
    texture = REGISTRY.call(
        "texture.pole_figure",
        {
            "phase": {"builtin": "ni_fcc"},
            "model": "random",
            "spread_deg": 10.0,
            "grain_count": 600,
            "halfwidth_deg": 10.0,
            "seed": 7,
            "pole": [1, 1, 1],
            "projection": "equal_area",
            "resolution_deg": 5.0,
        },
    )
    assert abs(texture["data"]["mean_mrd"] - 1.0) < 0.01
    assert "1.000" in guide, "the guide quotes the random baseline mean as 1.000"


def test_site_mathematics_renders_without_internet_access() -> None:
    """PyTex is read on a closed office intranet.

    The default MathJax CDN URL is unreachable there, which silently degraded
    every equation on the site to raw TeX. The bundle is vendored under
    _static/mathjax instead, and the config must keep pointing at it.
    """
    conf = _read("docs/site/conf.py")

    assert 'mathjax_path = "mathjax/tex-chtml-full.js"' in conf
    assert "cdn.jsdelivr.net" not in conf

    vendored = REPO_ROOT / "docs" / "site" / "_static" / "mathjax"
    bundle = vendored / "tex-chtml-full.js"
    assert bundle.is_file(), "the offline MathJax bundle is missing"
    # "tex-chtml-full" embeds every TeX extension, so nothing is fetched lazily.
    assert bundle.stat().st_size > 500_000
    assert (vendored / "LICENSE").is_file()

    fonts = sorted((vendored / "output" / "chtml" / "fonts" / "woff-v2").glob("*.woff"))
    assert len(fonts) >= 15, "MathJax web fonts are missing"


def test_vendored_site_assets_are_tracked_by_git() -> None:
    """Assets present only on a developer's disk would not reach a deployment."""
    import subprocess

    tracked = subprocess.run(
        ["git", "ls-files", "docs/site/_static/mathjax"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()

    assert any(name.endswith("tex-chtml-full.js") for name in tracked)
    assert sum(1 for name in tracked if name.endswith(".woff")) >= 15


def test_site_does_not_register_the_unused_mermaid_extension() -> None:
    """Mermaid is already forbidden on canonical visual pages in favour of SVG.

    Registering the extension anyway injected Mermaid and D3 script tags from
    jsdelivr onto every generated page, which simply fail on the intranet while
    the site contains no Mermaid diagram at all.
    """
    conf = _read("docs/site/conf.py")
    assert '"sphinxcontrib.mermaid"' not in conf

    sources = list((REPO_ROOT / "docs").rglob("*.md"))
    assert sources, "no documentation sources found"
    for path in sources:
        content = path.read_text(encoding="utf-8", errors="ignore")
        assert "```{mermaid}" not in content, path
        assert "```mermaid" not in content, path
