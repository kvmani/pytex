"""Concrete worked-example groups, one module per scientific domain."""

from __future__ import annotations

from ..framework import ExampleGroup
from . import (
    composite_diffraction,
    convergent_beam_diffraction,
    core_crystal_geometry,
    diffraction_geometry,
    directional_statistics,
    dynamical_cbed,
    ebsd_microstructure,
    elastic_anisotropy,
    ghost_problem,
    ipf_coloring,
    orientation_and_misorientation,
    orientation_representations,
    pole_figure_arithmetic,
    pole_figure_sampling,
    random_disorientation,
    reference_frames,
    schmid_and_taylor,
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
    schmid_and_taylor.GROUP,
    pole_figure_sampling.GROUP,
    directional_statistics.GROUP,
    ghost_problem.GROUP,
    pole_figure_arithmetic.GROUP,
    tem_tilt_navigation.GROUP,
    transformation_correspondence.GROUP,
    visualization_composition.GROUP,
)

__all__ = ["GROUPS"]
