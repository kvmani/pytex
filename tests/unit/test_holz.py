"""Tests for `pytex.diffraction.holz`.

HOLZ line geometry is closed-form, so it can be checked against the thing it is
supposed to describe rather than against stored numbers. Three checks carry the
module:

1. **A point on the line is at the Bragg condition.** The excitation errors are
   computed by `pytex.diffraction.dynamical`, which arrived at them
   independently, so agreeing to machine precision means the two modules share
   one geometry rather than two consistent-looking ones.
2. **The strain/wavelength degeneracy is exact.** A fractional lattice strain and
   a fractional wavelength change of the same size cancel at every reflection
   simultaneously. This is the reason HOLZ metrology needs a calibrated
   accelerating voltage, and it is asserted rather than described.
3. **The analytic strain sensitivity is the derivative of the exact offset.**
   Checked by central differences, which catches a wrong power of ``(1 + eps)``
   that a plausibility check would not.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import crystal_frame
from pytex.core.lattice import AtomicSite, Lattice, Phase, SpaceGroupSpec, UnitCell, ZoneAxis
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.dynamical import beam_set_from_indices
from pytex.diffraction.holz import HOLZ_LINE_PATTERN_SCHEMA, holz_line_pattern
from pytex.diffraction.kinematic import electron_wavelength_angstrom


@pytest.fixture(scope="module")
def nickel() -> Phase:
    """FCC nickel. Its first-order Laue zone down ``[001]`` is rich in lines."""

    frame = crystal_frame()
    parameter = 3.5239
    lattice = Lattice(parameter, parameter, parameter, 90.0, 90.0, 90.0, crystal_frame=frame)
    sites = tuple(
        AtomicSite(
            label=f"Ni{index}",
            species="Ni",
            fractional_coordinates=np.asarray(position, dtype=np.float64),
        )
        for index, position in enumerate(
            [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
        )
    )
    return Phase(
        "nickel-fcc",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=frame),
    )


@pytest.fixture(scope="module")
def pattern(nickel: Phase):  # type: ignore[no-untyped-def]
    """HOLZ lines of nickel down ``[001]`` at 200 kV with an 8 mrad probe."""

    return holz_line_pattern(
        nickel,
        ZoneAxis(indices=(0, 0, 1), phase=nickel),
        convergence_semi_angle_mrad=8.0,
        max_index=24,
        g_max_inv_angstrom=6.0,
    )


# --------------------------------------------------------------------------- #
# The line is the Bragg locus
# --------------------------------------------------------------------------- #


def test_a_point_on_the_line_is_exactly_at_the_bragg_condition(nickel: Phase, pattern) -> None:  # type: ignore[no-untyped-def]
    """``s_g = 0`` along the line, as computed by the *dynamical* module.

    The two modules derive the excitation error independently — one to place a
    line, the other to build a structure matrix — so this is a cross-check on the
    shared geometry rather than a restatement of one formula.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    assert pattern.lines
    for line in pattern.lines[:6]:
        chord = line.chord_tilt_rad()
        assert chord is not None
        midpoint = chord.mean(axis=0)
        samples = np.stack([chord[0], midpoint, chord[1]])
        beams = beam_set_from_indices(nickel, zone, line.miller_indices[None, :])
        errors = beams.excitation_errors(samples)[:, 1]
        assert errors == pytest.approx(0.0, abs=1e-15)


def test_the_line_normal_is_the_in_plane_direction_of_g(pattern) -> None:  # type: ignore[no-untyped-def]
    """A Bragg locus runs perpendicular to ``g``; its normal is ``g_perp``."""

    for line in pattern.lines[:10]:
        assert float(np.linalg.norm(line.normal_tilt)) == pytest.approx(1.0, abs=1e-12)
        chord = line.chord_tilt_rad()
        assert chord is not None
        direction = chord[1] - chord[0]
        assert float(np.dot(direction, line.normal_tilt)) == pytest.approx(0.0, abs=1e-15)


def test_chord_endpoints_sit_on_the_cone_boundary(pattern) -> None:  # type: ignore[no-untyped-def]
    """The visible part of a line ends where the illumination cone does."""

    alpha = pattern.convergence_semi_angle_mrad * 1e-3
    for line in pattern.bright_field_lines[:10]:
        chord = line.chord_tilt_rad()
        assert chord is not None
        assert np.linalg.norm(chord, axis=1) == pytest.approx(alpha, rel=1e-12)


def test_the_excess_line_is_the_deficiency_line_translated_to_its_own_disc(pattern) -> None:  # type: ignore[no-untyped-def]
    """Same Bragg condition, same incident tilts, different disc.

    Recognising the dark/bright pair is how a HOLZ line is distinguished from a
    bend contour in a real pattern, so the geometry has to state the relation
    rather than leave it implied.
    """

    for line in pattern.bright_field_lines[:5]:
        deficiency = line.deficiency_chord_mm()
        excess = line.excess_chord_mm()
        assert deficiency is not None and excess is not None
        assert excess - deficiency == pytest.approx(
            np.broadcast_to(line.disc_centre_mm, excess.shape), rel=1e-12
        )


def test_a_line_outside_the_cone_has_no_chord(nickel: Phase) -> None:
    """A margin admits lines that a wider probe would show, and marks them as such."""

    wide = holz_line_pattern(
        nickel,
        ZoneAxis(indices=(0, 0, 1), phase=nickel),
        convergence_semi_angle_mrad=4.0,
        max_index=24,
        g_max_inv_angstrom=6.0,
        offset_margin_rad=6e-3,
    )
    outside = [line for line in wide.lines if not line.crosses_bright_field]
    assert outside, "the margin should admit lines beyond the cone"
    for line in outside[:5]:
        assert line.chord_tilt_rad() is None
        assert line.deficiency_chord_mm() is None
        assert line.excess_chord_mm() is None
        assert "outside" in line.describe()


# --------------------------------------------------------------------------- #
# Metrology
# --------------------------------------------------------------------------- #


def test_lattice_strain_and_wavelength_change_are_exactly_degenerate(pattern) -> None:  # type: ignore[no-untyped-def]
    """The reason a HOLZ measurement is a voltage measurement until it is calibrated.

    Scaling the lattice by ``1 + eps`` and the wavelength by the same factor
    returns every line to its original position simultaneously. Nothing in the
    pattern can separate the two, so a lattice parameter quoted from an
    uncalibrated microscope is really a statement about its high-tension supply.
    """

    wavelength = pattern.wavelength_angstrom
    for strain in (1e-4, 1e-3, 2e-2):
        deviations = [
            abs(
                line.offset_at(
                    lattice_strain=strain, wavelength_angstrom=wavelength * (1.0 + strain)
                )
                - line.offset_rad
            )
            for line in pattern.lines
        ]
        assert max(deviations) < 1e-16


def test_strain_sensitivity_is_the_derivative_of_the_exact_offset(pattern) -> None:  # type: ignore[no-untyped-def]
    """Central differences against the closed form catch a wrong power of ``1 + eps``."""

    step = 1e-6
    for line in pattern.lines[:8]:
        numerical = (
            line.offset_at(lattice_strain=step) - line.offset_at(lattice_strain=-step)
        ) / (2.0 * step)
        assert numerical == pytest.approx(line.strain_sensitivity_rad, rel=1e-6)


def test_offset_at_zero_strain_reproduces_the_stored_offset(pattern) -> None:  # type: ignore[no-untyped-def]
    """The closed form and the constructor must agree, or one of them is wrong.

    They differ only in the order of the two divisions, so the tolerance is
    floating-point rounding rather than physics.
    """

    for line in pattern.lines:
        assert line.offset_at() == pytest.approx(line.offset_rad, rel=1e-14, abs=1e-18)


def test_line_width_falls_as_one_over_thickness(pattern) -> None:  # type: ignore[no-untyped-def]
    """Sharper lines in a thicker foil: the reason HOLZ metrology wants thick specimens.

    The resolvable strain improves in exact proportion, which is the practical
    consequence and the opposite of the usual thin-foil instinct.
    """

    line = pattern.lines[0]
    assert line.angular_width_rad(2000.0) == pytest.approx(
        0.5 * line.angular_width_rad(1000.0), rel=1e-12
    )
    assert line.resolvable_strain(2000.0) == pytest.approx(
        0.5 * line.resolvable_strain(1000.0), rel=1e-12
    )
    assert line.angular_width_rad(1000.0) == pytest.approx(
        1.0 / (1000.0 * line.g_perp_inv_angstrom), rel=1e-12
    )


def test_a_near_parallel_intersection_moves_faster_than_either_line(pattern) -> None:  # type: ignore[no-untyped-def]
    """Why HOLZ measurements read intersections rather than single lines.

    The crossing of two lines meeting at angle ``phi`` moves as ``1/sin(phi)``
    times faster than the lines do. For this pattern the best crossing improves
    the resolvable strain by more than an order of magnitude over the best single
    line, which is exactly the gap between what a line alone can do and the
    ``1e-4`` that the technique is known for.
    """

    intersections = pattern.intersections()
    assert intersections
    best = intersections[0]
    assert best.inside_bright_field
    assert best.angle_deg <= 90.0
    single = max(line.strain_sensitivity_rad for line in pattern.bright_field_lines)
    assert best.strain_sensitivity_rad > 10.0 * single


def test_an_intersection_lies_on_both_of_its_lines(pattern) -> None:  # type: ignore[no-untyped-def]
    """The crossing must satisfy both line equations, or it is not a crossing."""

    for intersection in pattern.intersections()[:8]:
        first = pattern.line_for(intersection.first_indices)
        second = pattern.line_for(intersection.second_indices)
        for line in (first, second):
            assert float(
                np.dot(intersection.position_tilt_rad, line.normal_tilt)
            ) == pytest.approx(line.offset_rad, abs=1e-12)


def test_accelerating_voltage_moves_every_line(nickel: Phase, pattern) -> None:  # type: ignore[no-untyped-def]
    """A one-percent voltage error is a one-percent-scale lattice error.

    The shift predicted by the closed form must match a full recomputation at the
    other voltage, so a caller can trust ``offset_at`` for the sensitivity study
    without re-running the enumeration.
    """

    other_energy = 202.0
    other_wavelength = electron_wavelength_angstrom(other_energy)
    recomputed = holz_line_pattern(
        nickel,
        ZoneAxis(indices=(0, 0, 1), phase=nickel),
        beam_energy_kev=other_energy,
        convergence_semi_angle_mrad=8.0,
        max_index=24,
        g_max_inv_angstrom=6.0,
    )
    moved = 0
    for line in pattern.lines[:12]:
        predicted = line.offset_at(wavelength_angstrom=other_wavelength)
        actual = recomputed.line_for(line.miller_indices).offset_rad
        assert predicted == pytest.approx(actual, abs=1e-15)
        if abs(predicted - line.offset_rad) > 1e-9:
            moved += 1
    assert moved == 12


# --------------------------------------------------------------------------- #
# Guardrails and reporting
# --------------------------------------------------------------------------- #


def test_zeroth_laue_zone_is_refused_with_the_reason(nickel: Phase) -> None:
    """A ZOLZ reflection has ``g_z = 0`` and produces no line inside the pattern."""

    with pytest.raises(ValueError, match="Laue zone 0 has no HOLZ lines"):
        holz_line_pattern(nickel, ZoneAxis(indices=(0, 0, 1), phase=nickel), laue_zones=(0, 1))


def test_a_parallel_beam_is_refused(nickel: Phase) -> None:
    """One incident direction samples one point, so it shows no line at all."""

    with pytest.raises(ValueError, match="parallel beam"):
        holz_line_pattern(
            nickel,
            ZoneAxis(indices=(0, 0, 1), phase=nickel),
            convergence_semi_angle_mrad=0.0,
        )


def test_holz_lines_do_not_need_a_unit_cell(nickel: Phase) -> None:
    """Line positions are lattice geometry; the atoms decide only visibility.

    A dynamical calculation refuses a bare lattice because the potential
    coefficients need atom positions. This one must not, because the question it
    answers does not.
    """

    bare = Phase(
        "bare-fcc",
        lattice=nickel.lattice,
        symmetry=nickel.symmetry,
        crystal_frame=nickel.crystal_frame,
        space_group=nickel.space_group,
    )
    result = holz_line_pattern(
        bare,
        ZoneAxis(indices=(0, 0, 1), phase=bare),
        convergence_semi_angle_mrad=8.0,
        max_index=24,
        g_max_inv_angstrom=6.0,
    )
    assert result.lines


def test_describe_and_json_stay_in_lockstep(pattern) -> None:  # type: ignore[no-untyped-def]
    """The explainable-results contract, including the degeneracy warning."""

    prose = pattern.describe()
    assert "[001]" in prose
    assert "accelerating voltage" in prose
    assert "Contrast is not" in prose

    payload = pattern.to_json_dict()
    assert payload["schema"] == HOLZ_LINE_PATTERN_SCHEMA
    assert payload["line_count"] == len(pattern.lines)
    assert payload["bright_field_line_count"] == len(pattern.bright_field_lines)
    assert len(payload["lines"]) == len(pattern.lines)

    indices = pattern.miller_indices
    assert indices.shape == (len(pattern.lines), 3)
    assert np.all(indices @ np.array([0, 0, 1]) != 0), "every line comes from a HOLZ reflection"

    with pytest.raises(KeyError, match="no HOLZ line here"):
        pattern.line_for([2, 2, 0])
