"""Concrete worked-example groups, one module per scientific domain."""

from __future__ import annotations

from ..framework import ExampleGroup
from . import (
    composite_diffraction,
    core_crystal_geometry,
    diffraction_geometry,
    ebsd_microstructure,
    orientation_and_misorientation,
    pole_figure_arithmetic,
    reference_frames,
    tem_tilt_navigation,
    texture_kernels,
    transformation_correspondence,
    visualization_composition,
)

GROUPS: tuple[ExampleGroup, ...] = (
    reference_frames.GROUP,
    core_crystal_geometry.GROUP,
    orientation_and_misorientation.GROUP,
    diffraction_geometry.GROUP,
    ebsd_microstructure.GROUP,
    composite_diffraction.GROUP,
    texture_kernels.GROUP,
    pole_figure_arithmetic.GROUP,
    tem_tilt_navigation.GROUP,
    transformation_correspondence.GROUP,
    visualization_composition.GROUP,
)

__all__ = ["GROUPS"]
