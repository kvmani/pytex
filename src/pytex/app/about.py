"""What the application says about itself.

Purpose
-------
One place answers "what is this program, which version am I running, who wrote
it, and under what licence may I use it". The workbench masthead has an *About*
entry, the desktop shell shows the same panel, and both read this document, so
the answer cannot differ between shells or drift from the package it describes.

The version is never written here: it is imported from :mod:`pytex._version`, the
single source of truth the package metadata and every written manifest also read.
The licence identifier is the SPDX expression declared in ``pyproject.toml`` and
carried by ``LICENSE``.

When and where to use it
------------------------
Call :func:`about_document` when a shell needs the identity block. It is embedded
in the application manifest under ``about``, so the browser already has it after
startup and needs no extra request.

Returns
-------
A JSON-serialisable mapping; see :func:`about_document`.
"""

from __future__ import annotations

from typing import Any

#: SPDX expression for the licence this program is distributed under. It matches
#: the ``license`` field of ``pyproject.toml`` and the header of ``LICENSE``.
LICENSE_SPDX = "GPL-3.0-or-later"

LICENSE_NAME = "GNU General Public License, version 3 or later"

LICENSE_URL = "https://www.gnu.org/licenses/gpl-3.0.txt"

#: The warranty disclaimer GPL-3.0 section 15 asks an interactive program to
#: display. Kept as prose rather than a link, because a user offline on an
#: instrument PC must still be able to read it.
LICENSE_NOTICE = (
    "PyTex is free software: you may redistribute it and modify it under the "
    "terms of the GNU General Public License as published by the Free Software "
    "Foundation, either version 3 of the licence, or (at your option) any later "
    "version. It is distributed in the hope that it will be useful, but WITHOUT "
    "ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or "
    "FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for "
    "more details."
)

DESCRIPTION = (
    "PyTex is a pure-Python library and workbench for crystallographic texture "
    "and diffraction. It carries one canonical data model — reference frames, "
    "symmetry, orientations and phases named explicitly rather than passed as "
    "bare arrays — and builds every workspace on it: crystal structure viewing, "
    "electron and X-ray diffraction simulation, TEM spot-pattern indexing and "
    "tilting, CBED, EBSD map analysis, variant and orientation-relationship "
    "analysis, and pole-figure and ODF texture analysis."
)

AUDIENCE = (
    "It is written for research and for teaching at the same standard: every "
    "result states the convention it was computed under, and the workbench is "
    "the same library a script would call."
)

AUTHOR_NAME = "Dr K V Mani Krishna"

AUTHOR_AFFILIATION = "Materials Group, Bhabha Atomic Research Centre (BARC)"

AUTHOR_EMAILS = ("kvmani@barc.gov.in", "kvmani@gmail.com")


def about_document() -> dict[str, Any]:
    """The identity block every shell displays in its About panel.

    When and where to use it
    ------------------------
    Called by :meth:`pytex.app.registry.ServiceRegistry.manifest` so the block
    travels with the manifest the frontend already fetches at startup. Call it
    directly from a shell that needs the same facts without a manifest.

    Returns
    -------
    dict
        ``name``, ``version``, ``tagline``, ``description``, ``audience``,
        ``author`` (``name``, ``affiliation``, ``emails``), ``license``
        (``spdx``, ``name``, ``url``, ``notice``), and ``links``: a list of
        ``{label, url}`` pairs. Every value is a string, a list of strings, or a
        mapping of those, so the document serialises without a custom encoder.

    Examples
    --------
    >>> from pytex import __version__
    >>> document = about_document()
    >>> document["version"] == __version__
    True
    >>> document["author"]["name"]
    'Dr K V Mani Krishna'
    >>> document["license"]["spdx"]
    'GPL-3.0-or-later'
    """

    from pytex import __version__

    return {
        "name": "PyTex",
        "version": __version__,
        "tagline": "Crystallographic texture and diffraction, for research and teaching.",
        "description": DESCRIPTION,
        "audience": AUDIENCE,
        "author": {
            "name": AUTHOR_NAME,
            "affiliation": AUTHOR_AFFILIATION,
            "emails": list(AUTHOR_EMAILS),
        },
        "license": {
            "spdx": LICENSE_SPDX,
            "name": LICENSE_NAME,
            "url": LICENSE_URL,
            "notice": LICENSE_NOTICE,
        },
        "links": [
            {"label": "Documentation", "url": "/docs/index.html"},
            {"label": "Source repository", "url": "https://github.com/kvmani/pytex"},
            {"label": "Licence text", "url": LICENSE_URL},
        ],
    }


__all__ = [
    "AUTHOR_AFFILIATION",
    "AUTHOR_EMAILS",
    "AUTHOR_NAME",
    "DESCRIPTION",
    "LICENSE_NAME",
    "LICENSE_NOTICE",
    "LICENSE_SPDX",
    "LICENSE_URL",
    "about_document",
]
