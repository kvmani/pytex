"""Azimuthal projections for texture plotting.

The implementation lives in :mod:`pytex.core.sphere`, which owns the canonical
unit-sphere semantics; this module re-exports it under the name the texture and
plotting subsystems have always used, so there is exactly one projection
formula in the repository.
"""

from __future__ import annotations

from pytex.core.sphere import fold_upper_hemisphere, project_directions

__all__ = ["fold_upper_hemisphere", "project_directions"]
