"""The PyTex application: one workbench, a desktop shell and an intranet shell.

This package is the application layer over the PyTex library. Its design is
recorded in ``docs/architecture/application_platform.md``; the short version is
that every user-visible capability is implemented once as a JSON-in/JSON-out
*service*, the user interface is generated from a self-describing operation
manifest, and the desktop and web shells differ only in how a window opens and
where files are written.

Importing this package registers every service in
:data:`pytex.app.registry.REGISTRY`, and imports nothing from a web or GUI
framework — the shells do that, lazily.

Entry points
------------
``python -m pytex.app serve``
    Start the intranet server.
``python -m pytex.app desktop``
    Open the desktop window.
"""

from __future__ import annotations

# Registers the operations. Kept last, and imported for its side effect, because
# the registry must be populated before anyone reads the manifest.
from pytex.app import services as services
from pytex.app.contracts import execute, success_envelope, to_jsonable
from pytex.app.errors import (
    DependencyMissingError,
    InvalidInputError,
    ServiceError,
    UnknownOperationError,
    UnsupportedRequestError,
)
from pytex.app.phases import BUILTIN_PHASES, PhaseSpec, SiteSpec, builtin_phase
from pytex.app.registry import REGISTRY, OperationSpec, ServiceRegistry

__all__ = [
    "BUILTIN_PHASES",
    "REGISTRY",
    "DependencyMissingError",
    "InvalidInputError",
    "OperationSpec",
    "PhaseSpec",
    "ServiceError",
    "ServiceRegistry",
    "SiteSpec",
    "UnknownOperationError",
    "UnsupportedRequestError",
    "builtin_phase",
    "execute",
    "services",
    "success_envelope",
    "to_jsonable",
]
