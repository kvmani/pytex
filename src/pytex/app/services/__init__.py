"""Every user-visible capability of the application, one module per domain.

Importing this package registers all operations in
:data:`pytex.app.registry.REGISTRY`. Nothing here imports a web or GUI
framework: a service takes a validated parameter mapping and returns a
:class:`~pytex.app.results.AppResult`, which is as true in a notebook as it is
behind an HTTP request.
"""

from __future__ import annotations

from pytex.app.services import calculator as calculator
from pytex.app.services import cbed as cbed
from pytex.app.services import crystal as crystal
from pytex.app.services import diffraction as diffraction
from pytex.app.services import tem as tem
from pytex.app.services import texture as texture
from pytex.app.services import variants as variants
from pytex.app.services import xrd as xrd

__all__ = [
    "calculator",
    "cbed",
    "crystal",
    "diffraction",
    "tem",
    "texture",
    "variants",
    "xrd",
]
