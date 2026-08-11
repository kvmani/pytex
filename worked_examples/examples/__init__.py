"""Concrete worked-example groups, one module per scientific domain."""

from __future__ import annotations

from ..framework import ExampleGroup
from . import (
    composite_diffraction,
    convergent_beam_diffraction,
    core_crystal_geometry,
    diffraction_geometry,
    dynamical_cbed,
    ebsd_microstructure,
    elastic_anisotropy,
    ipf_coloring,
    orientation_and_misorientation,
    orientation_representations,
    pole_figure_arithmetic,
    random_disorientation,
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
    orientation_representations.GROUP,
    diffraction_geometry.GROUP,
    ebsd_microstructure.GROUP,
    composite_diffraction.GROUP,
    convergent_beam_diffraction.GROUP,
    dynamical_cbed.GROUP,
    texture_kernels.GROUP,
    ipf_coloring.GROUP,
    random_disorientation.GROUP,
    elastic_anisotropy.GROUP,
    pole_figure_arithmetic.GROUP,
    tem_tilt_navigation.GROUP,
    transformation_correspondence.GROUP,
    visualization_composition.GROUP,
)

__all__ = ["GROUPS"]
