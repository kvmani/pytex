"""Lattice curvature and geometrically necessary dislocation density.

Expected values are analytic. The map is built by planting a *known* orientation
gradient through the canonical constructors — never through the code under
test — so the expectation is independent of the implementation:

* a lattice rotated about the specimen z axis by an angle linear in x has
  curvature ``kappa_20 = d(theta)/dx`` and no other component;
* Nye's relation then gives ``alpha_02 = kappa_20`` and nothing else, so
  ``rho_GND = (d(theta)/dx) / b`` in closed form;
* the KAM estimate ``rho = 2*theta/(b*u)`` must reproduce the same value for a
  pure single-axis tilt, which cross-checks two independent formulas.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import CRYSTAL_FRAME, MAP_FRAME, SPECIMEN_FRAME
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.ebsd.gnd import (
    geometrically_necessary_dislocation_density,
    lattice_curvature_tensor,
    nye_dislocation_density_tensor,
)
from pytex.ebsd.models import CrystalMap

_ROWS, _COLS = 9, 11
_STEP_UM = 0.5
#: Copper, from the FCC lattice parameter: b = a / sqrt(2) for a <110>/2 vector.
_BURGERS_NM = 0.2556


def _copper() -> Phase:
    return Phase(
        name="copper",
        lattice=Lattice(3.615, 3.615, 3.615, 90.0, 90.0, 90.0, crystal_frame=CRYSTAL_FRAME),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
    )


def _tilt_map(gradient_deg_per_um: float, *, axis: str = "x") -> CrystalMap:
    """A map whose lattice rotates about specimen z linearly along one axis."""

    phase = _copper()
    coordinates: list[tuple[float, float]] = []
    quaternions: list[np.ndarray] = []
    for row in range(_ROWS):
        for col in range(_COLS):
            x_um = col * _STEP_UM
            y_um = row * _STEP_UM
            coordinates.append((x_um, y_um))
            position = x_um if axis == "x" else y_um
            angle = np.deg2rad(gradient_deg_per_um * position)
            quaternions.append(Rotation.from_axis_angle([0.0, 0.0, 1.0], angle).quaternion)
    orientations = OrientationSet.from_quaternions(
        np.asarray(quaternions),
        specimen_frame=SPECIMEN_FRAME,
        phase=phase,
    )
    return CrystalMap(
        coordinates=np.asarray(coordinates, dtype=np.float64),
        orientations=orientations,
        map_frame=MAP_FRAME,
        grid_shape=(_ROWS, _COLS),
        step_sizes=(_STEP_UM, _STEP_UM),
    )


def _expected_curvature_rad_per_m(gradient_deg_per_um: float) -> float:
    """The planted gradient in radians per metre."""

    return float(np.deg2rad(gradient_deg_per_um) / 1e-6)


# --------------------------------------------------------------------------- #
# Curvature
# --------------------------------------------------------------------------- #


def test_planted_tilt_gradient_is_recovered_exactly() -> None:
    """The curvature component the gradient was planted in must match it."""

    gradient = 0.8
    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(gradient)))
    assert curvature.shape == (_ROWS, _COLS, 3, 3)
    expected = _expected_curvature_rad_per_m(gradient)
    # kappa_20 = d(omega_z)/dx for a rotation about specimen z varying along x.
    assert curvature[4, 5, 2, 0] == pytest.approx(expected, rel=1e-9)


def test_every_other_measurable_curvature_component_vanishes() -> None:
    """A single planted component must not leak into the others."""

    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(0.8)))
    interior = curvature[1:-1, 1:-1]
    for i in range(3):
        for j in range(2):
            if (i, j) == (2, 0):
                continue
            assert np.allclose(interior[:, :, i, j], 0.0, atol=1e-6), (i, j)


def test_the_out_of_plane_gradient_is_nan_not_zero() -> None:
    """A surface map cannot measure the depth gradient, and must say so.

    Reporting it as zero would silently assert that the lattice is unbent in
    depth, which no 2-D measurement can establish.
    """

    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(0.8)))
    assert np.all(np.isnan(curvature[..., 2]))


def test_gradient_along_y_lands_in_the_other_curvature_column() -> None:
    """Planting the gradient along y must move it to ``kappa_21``, not ``kappa_20``."""

    gradient = 0.6
    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(gradient, axis="y")))
    expected = _expected_curvature_rad_per_m(gradient)
    assert curvature[4, 5, 2, 1] == pytest.approx(expected, rel=1e-9)
    assert curvature[4, 5, 2, 0] == pytest.approx(0.0, abs=1e-6)


def test_curvature_scales_inversely_with_the_step_size_unit() -> None:
    """``step_scale_m`` is a unit declaration, so it must scale the result."""

    crystal_map = _tilt_map(0.8)
    micrometres = np.asarray(lattice_curvature_tensor(crystal_map, step_scale_m=1e-6))
    millimetres = np.asarray(lattice_curvature_tensor(crystal_map, step_scale_m=1e-3))
    assert millimetres[4, 5, 2, 0] == pytest.approx(micrometres[4, 5, 2, 0] / 1000.0, rel=1e-12)


def test_curvature_requires_a_regular_grid_with_step_sizes() -> None:
    phase = _copper()
    orientations = OrientationSet.from_quaternions(
        np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (6, 1)),
        specimen_frame=SPECIMEN_FRAME,
        phase=phase,
    )
    without_steps = CrystalMap(
        coordinates=np.array([[float(i), 0.0] for i in range(6)]),
        orientations=orientations,
        map_frame=MAP_FRAME,
        grid_shape=(2, 3),
    )
    with pytest.raises(ValueError, match="requires a CrystalMap with 2-D step_sizes"):
        lattice_curvature_tensor(without_steps)


# --------------------------------------------------------------------------- #
# Nye tensor
# --------------------------------------------------------------------------- #


def test_nye_relation_maps_the_planted_curvature_to_its_component() -> None:
    """``alpha_ij = kappa_ji`` off the diagonal, so ``alpha_02 = kappa_20``."""

    gradient = 0.8
    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(gradient)))
    nye = np.asarray(nye_dislocation_density_tensor(curvature))
    assert nye[4, 5, 0, 2] == pytest.approx(_expected_curvature_rad_per_m(gradient), rel=1e-9)


def test_unmeasurable_trace_does_not_poison_the_measurable_components() -> None:
    """A regression guard: ``NaN * 0`` is ``NaN``.

    Subtracting the trace as ``trace * identity`` propagates the NaN of the
    unmeasurable ``kappa_22`` into every off-diagonal component, destroying
    exactly the information the Nye tensor exists to carry. The trace must be
    applied to the diagonal alone. This was a real defect, caught by the planted
    gradient giving ``alpha_02 = NaN``.
    """

    curvature = np.asarray(lattice_curvature_tensor(_tilt_map(0.8)))
    nye = np.asarray(nye_dislocation_density_tensor(curvature))
    interior = nye[1:-1, 1:-1]
    assert np.all(np.isfinite(interior[:, :, 0, 1]))
    assert np.all(np.isfinite(interior[:, :, 0, 2]))
    assert np.all(np.isfinite(interior[:, :, 1, 0]))
    assert np.all(np.isfinite(interior[:, :, 1, 2]))
    # The diagonal genuinely is unmeasurable and must stay NaN.
    for index in range(3):
        assert np.all(np.isnan(interior[:, :, index, index]))


def test_nye_rejects_a_wrongly_shaped_curvature() -> None:
    with pytest.raises(ValueError, match=r"trailing shape \(3, 3\)"):
        nye_dislocation_density_tensor(np.zeros((4, 2)))


# --------------------------------------------------------------------------- #
# GND density
# --------------------------------------------------------------------------- #


def test_gnd_density_matches_the_closed_form_for_a_planted_tilt() -> None:
    """``rho = (d(theta)/dx) / b`` exactly, for a pure single-axis tilt."""

    gradient = 0.8
    density = np.asarray(
        geometrically_necessary_dislocation_density(
            _tilt_map(gradient), burgers_vector_nm=_BURGERS_NM
        )
    )
    expected = _expected_curvature_rad_per_m(gradient) / (_BURGERS_NM * 1e-9)
    assert density[4, 5] == pytest.approx(expected, rel=1e-9)
    # A 0.8 deg/um gradient in copper is a lightly deformed microstructure, so
    # the density must land in the physically sensible 1e13-1e14 range.
    assert 1e13 < density[4, 5] < 1e14


def test_an_undeformed_single_crystal_has_exactly_zero_density() -> None:
    density = np.asarray(
        geometrically_necessary_dislocation_density(
            _tilt_map(0.0), burgers_vector_nm=_BURGERS_NM
        )
    )
    assert float(np.nanmax(density)) == 0.0


def test_density_is_linear_in_the_planted_curvature() -> None:
    single = np.asarray(
        geometrically_necessary_dislocation_density(
            _tilt_map(0.8), burgers_vector_nm=_BURGERS_NM
        )
    )
    double = np.asarray(
        geometrically_necessary_dislocation_density(
            _tilt_map(1.6), burgers_vector_nm=_BURGERS_NM
        )
    )
    assert double[4, 5] == pytest.approx(2.0 * single[4, 5], rel=1e-9)


def test_density_is_inversely_proportional_to_the_burgers_vector() -> None:
    crystal_map = _tilt_map(0.8)
    small = np.asarray(
        geometrically_necessary_dislocation_density(crystal_map, burgers_vector_nm=0.25)
    )
    large = np.asarray(
        geometrically_necessary_dislocation_density(crystal_map, burgers_vector_nm=0.50)
    )
    assert large[4, 5] == pytest.approx(small[4, 5] / 2.0, rel=1e-12)


def test_the_kam_and_curvature_routes_agree_for_a_pure_tilt() -> None:
    """Two independent formulas must coincide where both are valid.

    For a single-axis tilt the 4-connectivity KAM at an interior point is half
    the per-step misorientation, so ``rho = 2*theta/(b*u)`` reduces to
    ``(d(theta)/dx)/b`` — exactly the curvature result. Agreement here checks
    the two implementations against each other, not against a recorded value.
    """

    crystal_map = _tilt_map(0.8)
    curvature_route = np.asarray(
        geometrically_necessary_dislocation_density(
            crystal_map, burgers_vector_nm=_BURGERS_NM, method="curvature"
        )
    )
    kam_route = np.asarray(
        geometrically_necessary_dislocation_density(
            crystal_map, burgers_vector_nm=_BURGERS_NM, method="kam"
        )
    )
    assert kam_route[4, 5] == pytest.approx(curvature_route[4, 5], rel=1e-6)


def test_density_rises_as_the_step_size_falls_at_fixed_curvature() -> None:
    """GND density from a map is resolution dependent, and must be seen to be.

    Declaring the same map to have a ten-times smaller physical step is
    declaring a ten-times steeper gradient, so the density must rise by ten.
    This is the documented resolution dependence, pinned so it cannot be
    silently normalized away.
    """

    crystal_map = _tilt_map(0.8)
    coarse = np.asarray(
        geometrically_necessary_dislocation_density(
            crystal_map, burgers_vector_nm=_BURGERS_NM, step_scale_m=1e-6
        )
    )
    fine = np.asarray(
        geometrically_necessary_dislocation_density(
            crystal_map, burgers_vector_nm=_BURGERS_NM, step_scale_m=1e-7
        )
    )
    assert fine[4, 5] == pytest.approx(10.0 * coarse[4, 5], rel=1e-9)


def test_density_rejects_invalid_arguments() -> None:
    crystal_map = _tilt_map(0.8)
    with pytest.raises(ValueError, match="burgers_vector_nm must be finite"):
        geometrically_necessary_dislocation_density(crystal_map, burgers_vector_nm=0.0)
    with pytest.raises(ValueError, match="method must be either"):
        geometrically_necessary_dislocation_density(
            crystal_map,
            burgers_vector_nm=_BURGERS_NM,
            method="nye",  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="step_scale_m must be finite"):
        lattice_curvature_tensor(crystal_map, step_scale_m=0.0)


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #


def test_gnd_density_map_renders_with_the_units_it_plots() -> None:
    """The colorbar must name the scale actually used, log or linear.

    A GND map on a linear scale is dominated by boundary artefacts, so the
    logarithmic default is the useful one — and the label has to say which was
    drawn, or the numbers on it are uninterpretable.
    """

    pytest.importorskip("matplotlib")
    from pytex.plotting import plot_gnd_density_map

    crystal_map = _tilt_map(0.8)
    figure = plot_gnd_density_map(crystal_map, burgers_vector_nm=_BURGERS_NM)
    axes = figure.axes[0]
    assert axes.get_title() == "Geometrically Necessary Dislocation Density"
    assert "log" in figure.axes[-1].get_ylabel()
    assert "m$^{-2}$" in figure.axes[-1].get_ylabel()

    linear = plot_gnd_density_map(crystal_map, burgers_vector_nm=_BURGERS_NM, log_scale=False)
    assert "log" not in linear.axes[-1].get_ylabel()


def test_gnd_density_map_plots_the_computed_field() -> None:
    """The rendered image must be the density, not a re-derived quantity."""

    pytest.importorskip("matplotlib")
    from pytex.plotting import plot_gnd_density_map

    crystal_map = _tilt_map(0.8)
    expected = np.log10(
        np.asarray(
            geometrically_necessary_dislocation_density(
                crystal_map, burgers_vector_nm=_BURGERS_NM
            )
        )
    )
    figure = plot_gnd_density_map(crystal_map, burgers_vector_nm=_BURGERS_NM)
    drawn = figure.axes[0].get_images()[0].get_array()
    assert np.allclose(np.asarray(drawn), expected, equal_nan=True)
