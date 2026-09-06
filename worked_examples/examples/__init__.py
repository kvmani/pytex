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
    kearns_parameter,
    orientation_and_misorientation,
    orientation_representations,
    phase_identification,
    pole_figure_arithmetic,
    pole_figure_sampling,
    precise_lattice_parameters,
    random_disorientation,
    reference_frames,
    saed_practice_patterns,
    schmid_and_taylor,
    tem_tilt_navigation,
    texture_kernels,
    transformation_correspondence,
    visualization_composition,
    workbench_service_layer,
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
    kearns_parameter.GROUP,
    tem_tilt_navigation.GROUP,
    saed_practice_patterns.GROUP,
    precise_lattice_parameters.GROUP,
    phase_identification.GROUP,
    transformation_correspondence.GROUP,
    visualization_composition.GROUP,
    workbench_service_layer.GROUP,
)

__all__ = ["GROUPS"]
