"""The single source of truth for the PyTex version.

Purpose
-------
The version appears in the installed package metadata, in ``pytex.__version__``,
and stamped into every manifest PyTex writes. Those must agree: a manifest that
records a version the code was not built from is worse than one recording no
version at all, because it looks authoritative.

This module holds the only literal. ``pyproject.toml`` reads it statically for
the package metadata, :mod:`pytex` re-exports it, and the manifest writers import
it. It deliberately imports nothing, so build-time static reading and runtime
import are both trivially safe.

Versioning follows the pre-1.0 policy in
``docs/standards/api_stability_and_deprecation.md``: minor versions may break
with deprecation warnings, patch versions may not.
"""

from __future__ import annotations

__version__ = "0.4.0"

__all__ = ["__version__"]
