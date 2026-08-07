"""Release-metadata gates: the package must be installable and importable.

These tests exist because both of the conditions they check were violated at
once, and neither was visible from inside the development environment — where
every optional dependency happens to be present.

1. **Undeclared runtime dependency.** ``pytex.diffraction.solving`` and
   ``pytex.plotting.styles`` import ``yaml`` at module level, but ``pyyaml`` was
   not a declared dependency. A fresh ``pip install pytex`` therefore failed on
   the very first ``import pytex``.
2. **Import-time coupling to an optional stack.** ``pytex.plotting.crystal3d``
   and ``pytex.plotting.scene3d`` imported ``matplotlib.colors`` at module
   level, silently making the ``plotting`` extra mandatory and contradicting the
   ``_require_matplotlib()`` guards used everywhere else.
3. **A version literal in four files.** Bumping the release version would have
   left the manifest writers stamping the old one.

The checks are static, so they hold regardless of what happens to be installed
in the environment running them.
"""

from __future__ import annotations

import ast
import sys
import tomllib
from pathlib import Path

import pytest

import pytex

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPO_ROOT / "src" / "pytex"
PYPROJECT = REPO_ROOT / "pyproject.toml"

#: Distribution names that supply an importable module under a different name.
_IMPORT_NAME_BY_DISTRIBUTION = {
    "pyyaml": "yaml",
}


def _declared_runtime_imports() -> set[str]:
    """Top-level module names the declared runtime dependencies provide."""

    payload = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    names: set[str] = set()
    for requirement in payload["project"]["dependencies"]:
        distribution = (
            requirement.split(">")[0].split("<")[0].split("=")[0].split("[")[0].strip().lower()
        )
        names.add(_IMPORT_NAME_BY_DISTRIBUTION.get(distribution, distribution))
    return names


def _unconditional_module_level_imports(tree: ast.Module) -> set[str]:
    """Third-party top-level modules imported when the module is imported.

    Excludes two legitimate patterns that do not execute an import at runtime or
    that guard it explicitly: ``if TYPE_CHECKING:`` blocks, and ``try: import x
    except ImportError:`` optional-dependency guards. Everything else at module
    level runs on ``import pytex`` and therefore must be a declared dependency.
    """

    found: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Try):
            continue
        if isinstance(node, ast.If):
            test = node.test
            is_type_checking = (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
                isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
            )
            if is_type_checking:
                continue
        if isinstance(node, ast.Import):
            found.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            found.add(node.module.split(".")[0])
    return {
        name
        for name in found
        if name not in sys.stdlib_module_names and name != "pytex" and not name.startswith("_")
    }


_SOURCE_FILES = sorted(
    path for path in SOURCE_ROOT.rglob("*.py") if "__pycache__" not in str(path)
)


@pytest.mark.parametrize(
    "source_path",
    _SOURCE_FILES,
    ids=[str(path.relative_to(SOURCE_ROOT)) for path in _SOURCE_FILES],
)
def test_module_level_imports_are_declared_dependencies(source_path: Path) -> None:
    """Importing PyTex must not require anything beyond its declared dependencies.

    An optional stack may only be imported inside a function or behind a
    ``try``/``except ImportError`` guard, so that ``import pytex`` succeeds on a
    minimal install and the missing extra is reported with a useful message at
    the point of use.
    """

    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imported = _unconditional_module_level_imports(tree)
    undeclared = imported - _declared_runtime_imports()
    assert not undeclared, (
        f"{source_path.relative_to(REPO_ROOT)} imports {sorted(undeclared)} at module level, "
        "but they are not declared runtime dependencies in pyproject.toml. Either declare "
        "them, or import them lazily inside the functions that need them (see "
        "pytex.plotting.crystal3d._to_hex for the pattern)."
    )


def test_importing_pytex_needs_only_declared_dependencies() -> None:
    """The whole public package, not just its top module, must be minimal-safe."""

    declared = _declared_runtime_imports()
    offenders: dict[str, set[str]] = {}
    for source_path in _SOURCE_FILES:
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        undeclared = _unconditional_module_level_imports(tree) - declared
        if undeclared:
            offenders[str(source_path.relative_to(REPO_ROOT))] = undeclared
    assert not offenders, offenders


# --------------------------------------------------------------------------- #
# Version single-sourcing
# --------------------------------------------------------------------------- #


def test_the_version_has_exactly_one_literal_in_the_source_tree() -> None:
    """``src/pytex/_version.py`` must be the only place the version is written.

    A second literal drifts silently: the manifest writers previously carried
    their own copy, so a release bump would have stamped every exported manifest
    with the previous version.
    """

    version = pytex.__version__
    carriers = [
        path
        for path in _SOURCE_FILES
        if f'"{version}"' in path.read_text(encoding="utf-8")
        or f"'{version}'" in path.read_text(encoding="utf-8")
    ]
    assert [path.name for path in carriers] == ["_version.py"], (
        "The version literal must appear only in src/pytex/_version.py; found it in "
        f"{[str(p.relative_to(REPO_ROOT)) for p in carriers]}."
    )


def test_packaging_metadata_reads_the_single_version_source() -> None:
    """``pyproject.toml`` must derive the version rather than restate it."""

    payload = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    assert "version" not in payload["project"], (
        "pyproject.toml must not hardcode a version; it is derived from "
        "pytex._version.__version__."
    )
    assert payload["project"]["dynamic"] == ["version"]
    dynamic = payload["tool"]["setuptools"]["dynamic"]
    assert dynamic["version"] == {"attr": "pytex._version.__version__"}


def test_citation_metadata_matches_the_package_version() -> None:
    """A citation file naming the wrong version misattributes results."""

    citation = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert f"version: {pytex.__version__}" in citation, (
        "CITATION.cff must state the same version as the package."
    )


def test_manifest_writers_stamp_the_package_version() -> None:
    """Manifests record the version that produced them, not a stale copy."""

    from pytex.adapters import ebsd as ebsd_adapter
    from pytex.adapters import manifests as manifest_adapter

    assert ebsd_adapter._PYTEX_VERSION == pytex.__version__
    assert manifest_adapter._PYTEX_VERSION == pytex.__version__


def test_the_phase_fixture_corpus_reports_its_own_availability() -> None:
    """A repository-only asset must be detectable, not discovered by exception.

    ``get_phase_fixture`` is a public export whose data ships with the
    repository rather than with the wheel. Code that may run from either must be
    able to ask, and the failure message must say what to do instead.
    """

    from pytex.core.fixtures import phase_fixture_catalog_path, phase_fixtures_available

    assert phase_fixtures_available() is phase_fixture_catalog_path().is_file()
    # The tests run from a source checkout, so the corpus must be present here.
    assert phase_fixtures_available()
