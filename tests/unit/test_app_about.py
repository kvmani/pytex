"""The About document is the program's statement of what it is.

It names a version, an author and a licence, and each of those is checked here
against its own source of truth rather than against a copy: the version against
``pytex.__version__``, the licence against ``pyproject.toml`` and ``LICENSE``,
and the whole document against the manifest the browser actually receives.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

import pytest

from pytex import __version__
from pytex.app import REGISTRY
from pytex.app.about import LICENSE_SPDX, about_document
from pytex.app.contracts import dumps

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_about_reports_the_running_version() -> None:
    document = about_document()
    assert document["name"] == "PyTex"
    assert document["version"] == __version__


def test_about_names_the_author_with_affiliation_and_contacts() -> None:
    author = about_document()["author"]
    assert author["name"] == "Dr K V Mani Krishna"
    assert "Materials Group" in author["affiliation"]
    assert "Bhabha Atomic Research Centre" in author["affiliation"]
    assert author["emails"] == ["kvmani@barc.gov.in", "kvmani@gmail.com"]


def test_about_licence_matches_the_packaging_metadata() -> None:
    """The displayed licence must be the licence the package is distributed under.

    A program that shows one licence in its About panel and ships another in its
    metadata is making a legal claim it cannot support, so this reads both.
    """

    pyproject = tomllib.loads((REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["license"] == LICENSE_SPDX

    licence_text = (REPOSITORY_ROOT / "LICENSE").read_text(encoding="utf-8")
    assert f"SPDX-License-Identifier: {LICENSE_SPDX}" in licence_text

    licence = about_document()["license"]
    assert licence["spdx"] == LICENSE_SPDX
    assert "GNU General Public License" in licence["name"]
    assert licence["url"].startswith("https://")
    # GPL-3.0 section 15: an interactive program displays the warranty
    # disclaimer, so the notice must actually contain it.
    assert "WITHOUT ANY WARRANTY" in licence["notice"]


def test_about_describes_the_program_in_prose() -> None:
    document = about_document()
    assert len(document["description"]) > 200
    assert len(document["tagline"].strip()) > 10
    for topic in ("texture", "diffraction", "EBSD"):
        assert topic in document["description"], f"the description never mentions {topic}"


@pytest.mark.parametrize("link", about_document()["links"], ids=lambda link: link["label"])
def test_about_links_are_labelled_absolute_urls(link: dict[str, str]) -> None:
    assert link["label"].strip()
    assert link["url"].startswith("https://")


def test_manifest_carries_the_about_document() -> None:
    """The frontend reads About off the manifest, so it must travel with it.

    Serialised through the application's own encoder, because a document the
    browser cannot receive is not published however well formed it is in Python.
    """

    manifest = REGISTRY.manifest()
    assert manifest["about"] == about_document()
    delivered = json.loads(dumps(manifest))
    assert delivered["about"]["author"]["name"] == "Dr K V Mani Krishna"
    assert delivered["about"]["version"] == __version__


def test_about_page_is_reachable_from_the_shell_markup() -> None:
    """The masthead must offer the panel; a document nobody can open is not one."""

    markup = (REPOSITORY_ROOT / "src" / "pytex" / "app" / "static" / "index.html").read_text(
        encoding="utf-8"
    )
    assert 'id="open-about"' in markup
    assert 'id="about-drawer"' in markup
    assert 'id="about-body"' in markup
