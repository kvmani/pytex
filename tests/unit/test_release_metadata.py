"""Release-metadata gates: the package must be installable and importable.

These tests exist because both of the conditions they check were violated at
once, and neither was visible from inside the development environment — where
every optional dependency happens to be present.

1. **Undeclared runtime dependency.** ``pytex.diffraction.solving`` and
   ``pytex.plotting.styles`` import ``yaml`` at module level, but ``pyyaml`` was
   not a declared dependency. A fresh ``pip install pytex`` therefore failed on
   the very first ``import pytex``.
2. **Import-time coupling to a heavy stack.** ``pytex.plotting.crystal3d`` and
   ``pytex.plotting.scene3d`` imported ``matplotlib.colors`` at module level,
   making every ``import pytex`` pay for matplotlib. It is a required
   dependency now, so this is a cost rather than a failure -- but the lazy
   imports the rest of the plotting layer uses are still the contract.
3. **A version literal in four files.** Bumping the release version would have
   left the manifest writers stamping the old one.

The checks are static, so they hold regardless of what happens to be installed
in the environment running them.
"""

from __future__ import annotations

import ast
import importlib
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

#: Declared dependencies whose top-level module PyTex never imports itself.
#:
#: They are still required, and this file still checks they import: the adapter
#: surfaces are written against them and read the caller's object structurally,
#: so `pytex.adapters.index_hough` calls `signal.hough_indexing(...)` on an
#: object only KikuchiPy makes. An environment without them cannot run those
#: surfaces, which is the definition of a runtime dependency.
_DEPENDENCIES_PYTEX_DOES_NOT_IMPORT = frozenset({"diffsims", "kikuchipy"})


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

    Anything imported at module level runs on ``import pytex``, so it must be a
    declared dependency. A package that is *not* declared -- a test-only tool, a
    truly optional extra like ``pywebview`` -- may only be imported inside a
    function or behind a ``try``/``except ImportError`` guard.
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


@pytest.mark.parametrize("module_name", sorted(_declared_runtime_imports()))
def test_every_declared_dependency_actually_imports(module_name: str) -> None:
    """A declared dependency that is not installed is a broken environment.

    This is the executable form of the 0.5.0 packaging decision. PyTex used to
    keep its scientific stack behind an optional extra, so the same call could
    read a CIF on one machine and raise `ModuleNotFoundError` on another -- and
    the machine it raised on was the deployed one. Making the stack required is
    only worth something if "required" is checked, so every distribution named
    in `[project] dependencies` is imported here.

    It also catches the quieter failure: an import name that does not match its
    distribution name and was never added to `_IMPORT_NAME_BY_DISTRIBUTION`,
    which would otherwise make this check silently vacuous.
    """

    importlib.import_module(module_name)


def test_the_dependencies_pytex_does_not_import_are_still_declared() -> None:
    """The adapter-facing packages must not be quietly dropped.

    Nothing in `src/pytex` contains the string `import kikuchipy`, so the
    module-level import checks above can never notice if KikuchiPy leaves the
    dependency list. The bridges would keep importing and keep type-checking,
    and `index_hough` would fail at the call site in a user's session instead.
    """

    declared = _declared_runtime_imports()
    missing = _DEPENDENCIES_PYTEX_DOES_NOT_IMPORT - declared
    assert not missing, (
        f"{sorted(missing)} are required by the adapter surfaces but are no longer declared "
        "in pyproject.toml. See the comment beside them there."
    )


# --------------------------------------------------------------------------- #
# Version single-sourcing
# --------------------------------------------------------------------------- #


def test_the_version_has_exactly_one_literal_in_the_source_tree() -> None:
    """``src/pytex/_version.py`` must be the only place the version is written.

    Source *and* tests. A second literal drifts silently: the manifest writers
    previously carried their own copy, so a release bump would have stamped every
    exported manifest with the previous version. A literal in a test fails
    differently and just as uselessly — cutting 0.1.1 broke the `python -m pytex
    info` case, which had transcribed `0.1.0.dev0` and therefore tested the
    release number rather than the command.
    """

    version = pytex.__version__
    # The tests are scanned too. A version literal in a *test* is the same defect
    # wearing a different hat: it does not break the package, it breaks the suite
    # on the release commit, and it says nothing about the behaviour it covers.
    # `python -m pytex info` carried one and failed exactly that way.
    scanned = [
        *_SOURCE_FILES,
        *(
            path
            for path in (REPO_ROOT / "tests").rglob("*.py")
            if "__pycache__" not in str(path)
        ),
    ]
    carriers = [
        path
        for path in scanned
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


def test_the_documentation_build_reads_the_version_rather_than_restating_it() -> None:
    """A hard-coded ``release`` in ``conf.py`` prints the wrong number on every page.

    It is the same defect the source-tree check above exists for, one directory
    outside its reach: ``docs/site/conf.py`` carried its own literal, so the
    built site kept claiming the previous version indefinitely — nothing failed,
    and every page said so.
    """

    conf = (REPO_ROOT / "docs" / "site" / "conf.py").read_text(encoding="utf-8")
    assert f'"{pytex.__version__}"' not in conf, (
        "docs/site/conf.py must not restate the version literal."
    )
    assert "from pytex._version import __version__" in conf
    assert "release = __version__" in conf


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
