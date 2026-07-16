"""Concrete worked-example groups, one module per scientific domain."""

from __future__ import annotations

from ..framework import ExampleGroup
from . import (
    core_crystal_geometry,
    diffraction_geometry,
    orientation_and_misorientation,
    transformation_correspondence,
    visualization_composition,
)

GROUPS: tuple[ExampleGroup, ...] = (
    core_crystal_geometry.GROUP,
    orientation_and_misorientation.GROUP,
    diffraction_geometry.GROUP,
    transformation_correspondence.GROUP,
    visualization_composition.GROUP,
)

__all__ = ["GROUPS"]
