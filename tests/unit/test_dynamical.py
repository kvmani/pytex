"""Tests for `pytex.diffraction.dynamical`.

A many-beam calculation is unusually easy to get subtly wrong and unusually easy
to check, because it has three exact properties that a wrong implementation
cannot fake:

1. **With two beams it must reproduce the closed form** of
   `pytex.diffraction.cbed.two_beam_rocking_curve` to machine precision. A wrong
   scale factor, a missing ``cos(theta_B)``, or a factor of ``pi`` in the wrong
   place all survive a plausibility check and all fail this one.
2. **Without absorption the propagator is unitary**, so the beam intensities sum
   to one at every thickness and every incident direction. Assuming the
   eigenvectors are orthogonal — the classic Bloch-wave implementation error —
   breaks this.
3. **Normal absorption is exactly a scalar.** It sits on every diagonal element,
   so it factors out of the exponential and cannot change a relative intensity.

The physical content is then pinned by two theorems rather than by stored
numbers: the Hashimoto-Howie-Whelan result that absorption makes the
bright-field rocking curve asymmetric while leaving the dark-field one
symmetric, and the propagator-symmetry statement that Friedel's law holds if and
only if the sampled structure is centrosymmetric.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import crystal_frame
from pytex.core.lattice import AtomicSite, Lattice, Phase, SpaceGroupSpec, UnitCell, ZoneAxis
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.cbed import extinction_distance_angstrom, two_beam_rocking_curve
from pytex.diffraction.dynamical import (
    BLOCH_WAVE_SOLUTION_SCHEMA,
    AbsorptionModel,
    beam_set_for_zone,
    beam_set_from_indices,
    potential_coefficients_inv_angstrom,
    solve_bloch_waves,
    structure_matrix,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom

_FCC_SITES = ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0))


def _cubic_phase(
    name: str,
    parameter_angstrom: float,
    species_and_offsets: tuple[tuple[str, tuple[float, float, float]], ...],
    point_group: str,
    space_group: tuple[str, int],
) -> Phase:
    """An FCC-derived cubic phase: one FCC sublattice per (species, offset) pair."""

    frame = crystal_frame()
    lattice = Lattice(
        parameter_angstrom,
        parameter_angstrom,
        parameter_angstrom,
        90.0,
        90.0,
        90.0,
        crystal_frame=frame,
    )
    sites = tuple(
        AtomicSite(
            label=f"{species}{index}",
            species=species,
            fractional_coordinates=np.asarray(base, dtype=np.float64)
            + np.asarray(offset, dtype=np.float64),
        )
        for species, offset in species_and_offsets
        for index, base in enumerate(_FCC_SITES)
    )
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group=SpaceGroupSpec(
            symbol=space_group[0], number=space_group[1], reference_frame=frame
        ),
    )


@pytest.fixture(scope="module")
def nickel() -> Phase:
    """FCC nickel: the centrosymmetric elemental reference case."""

    return _cubic_phase(
        "nickel-fcc", 3.5239, (("Ni", (0.0, 0.0, 0.0)),), "m-3m", ("Fm-3m", 225)
    )


@pytest.fixture(scope="module")
def gallium_arsenide() -> Phase:
    """Zincblende GaAs: the classic **non**-centrosymmetric CBED test case.

    Two FCC sublattices offset by ``(1/4, 1/4, 1/4)``. Its point group ``-43m``
    has no operation that reverses ``[111]``, which is exactly why CBED down
    ``[111]`` can determine the polarity of a GaAs crystal.
    """

    return _cubic_phase(
        "gallium-arsenide-zincblende",
        5.6535,
        (("Ga", (0.0, 0.0, 0.0)), ("As", (0.25, 0.25, 0.25))),
        "-43m",
        ("F-43m", 216),
    )


@pytest.fixture(scope="module")
def rocksalt_control() -> Phase:
    """The same two elements on the same lattice, but **centrosymmetric**.

    The offset is ``(1/2, 1/2, 1/2)`` instead of ``(1/4, 1/4, 1/4)``, which adds
    a centre of symmetry while changing nothing else — same cell, same species,
    same scattering strengths. It is the control that shows the Friedel test
    responds to the centre of symmetry and not to the chemistry.
    """

    return _cubic_phase(
        "rocksalt-control",
        5.6535,
        (("Ga", (0.0, 0.0, 0.0)), ("As", (0.5, 0.5, 0.5))),
        "m-3m",
        ("Fm-3m", 225),
    )


def _tilts_for_excitation_errors(
    beams: object, beam_index: int, targets: np.ndarray
) -> np.ndarray:
    """Incident tilts placing one beam at prescribed excitation errors.

    The excitation error is affine in the tilt, so walking along the in-plane
    direction of ``g`` sweeps ``s`` linearly; inverting that gives the tilts a
    two-beam comparison needs.
    """

    g_zone = np.asarray(beams.g_zone[beam_index], dtype=np.float64)  # type: ignore[attr-defined]
    in_plane = float(np.linalg.norm(g_zone[:2]))
    wavelength = electron_wavelength_angstrom(beams.beam_energy_kev)  # type: ignore[attr-defined]
    magnitude = float(beams.g_magnitude_inv_angstrom[beam_index])  # type: ignore[attr-defined]
    zero_tilt = g_zone[2] - 0.5 * wavelength * magnitude * magnitude
    distances = (zero_tilt - np.asarray(targets, dtype=np.float64)) / in_plane
    return distances[:, None] * (g_zone[:2] / in_plane)[None, :]


# --------------------------------------------------------------------------- #
# The absolute scale: two beams must reproduce the closed form
# --------------------------------------------------------------------------- #


def test_two_beam_limit_reproduces_the_closed_form(nickel: Phase) -> None:
    """With one reflection the solver must equal ``two_beam_rocking_curve`` exactly.

    This is the calibration of the whole module. The closed form is
    ``sin^2(pi t s_eff) / (xi_g s_eff)^2``, and reproducing it to machine
    precision pins the diagonal convention (``2 s_g``), the off-diagonal scale
    (``|nu_g| = 1/xi_g``, including the ``cos(theta_B)`` factor), and the ``i
    pi`` in the propagator simultaneously. Any one of them wrong gives a curve
    of the right general shape and the wrong fringe spacing.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_from_indices(nickel, zone, [[2, 2, 0]])
    assert beams.size == 2
    assert np.array_equal(beams.miller_indices[0], np.zeros(3, dtype=np.int64))

    targets = np.linspace(-0.02, 0.02, 41)
    tilts = _tilts_for_excitation_errors(beams, 1, targets)
    assert beams.excitation_errors(tilts)[:, 1] == pytest.approx(targets, abs=1e-15)

    thickness = 800.0
    solution = solve_bloch_waves(beams, tilts, thickness_angstrom=thickness)
    expected = two_beam_rocking_curve(
        targets,
        thickness_angstrom=thickness,
        extinction_distance_angstrom=float(
            extinction_distance_angstrom(nickel, [[2, 2, 0]])[0]
        ),
    )
    assert solution.intensities[:, 1] == pytest.approx(expected, abs=1e-12)
    assert solution.transmitted_intensity == pytest.approx(1.0 - expected, abs=1e-12)


def test_potential_coefficient_modulus_is_the_reciprocal_extinction_distance(
    nickel: Phase,
) -> None:
    """``|nu_g| = 1/xi_g`` ties the many-beam scale to the published two-beam one.

    The extinction distances themselves are validated against Williams and
    Carter in ``test_cbed.py``; this test is what makes that validation carry
    over to the dynamical module instead of being a second, independent scale.
    """

    reflections = [[1, 1, 1], [2, 0, 0], [2, 2, 0], [3, 1, 1]]
    coefficients = potential_coefficients_inv_angstrom(nickel, reflections)
    distances = extinction_distance_angstrom(nickel, reflections)
    assert np.abs(coefficients) == pytest.approx(1.0 / distances, rel=1e-12)


def test_potential_coefficients_are_real_only_for_a_centrosymmetric_structure(
    rocksalt_control: Phase, gallium_arsenide: Phase
) -> None:
    """The phase of ``nu_g`` is the centre of symmetry, made numerical.

    With the origin on a centre, every coefficient is real. Displace one
    sublattice to ``(1/4, 1/4, 1/4)`` and the odd-index coefficients acquire a
    phase — which is the only thing that distinguishes the two structures
    dynamically, and therefore the only thing a symmetry determination can be
    reading.
    """

    odd = [[1, 1, 1], [3, 1, 1]]
    assert np.max(np.abs(np.imag(
        potential_coefficients_inv_angstrom(rocksalt_control, odd)
    ))) < 1e-15
    assert np.max(np.abs(np.imag(
        potential_coefficients_inv_angstrom(gallium_arsenide, odd)
    ))) > 1e-4


# --------------------------------------------------------------------------- #
# Conservation, absorption, and the structure matrix
# --------------------------------------------------------------------------- #


def test_intensity_is_conserved_without_absorption(nickel: Phase) -> None:
    """Unitarity: the elastic structure matrix is Hermitian, so the beams sum to one.

    This is the check that the coupled solution is being recombined correctly.
    The eigenvectors of a complex matrix are not orthogonal, so obtaining the
    excitation amplitudes by projection rather than by solving would break this
    identity while still producing plausible rocking curves.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    assert beams.size > 8

    tilts = np.random.default_rng(0).uniform(-6e-3, 6e-3, size=(64, 2))
    for thickness in (100.0, 700.0, 2500.0):
        solution = solve_bloch_waves(beams, tilts, thickness_angstrom=thickness)
        assert solution.total_intensity == pytest.approx(1.0, abs=1e-12)


def test_structure_matrix_is_hermitian_without_absorption(nickel: Phase) -> None:
    """The elastic matrix must be Hermitian; the absorptive one must not be."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    tilts = np.array([[0.0, 0.0], [2e-3, -1e-3]])

    elastic = structure_matrix(beams, tilts)
    assert elastic == pytest.approx(np.conj(np.swapaxes(elastic, 1, 2)), abs=1e-14)

    absorbing = structure_matrix(beams, tilts, absorption=AbsorptionModel())
    assert np.max(np.abs(absorbing - np.conj(np.swapaxes(absorbing, 1, 2)))) > 1e-6


def test_structure_matrix_is_symmetric_exactly_when_centrosymmetric(
    rocksalt_control: Phase, gallium_arsenide: Phase
) -> None:
    """Symmetry of ``A`` — not Hermiticity — is the centre of symmetry.

    ``A`` is always Hermitian for a real potential. It is *symmetric* only when
    every included coefficient is real, and that is the propagator property the
    Friedel test reads.
    """

    tilts = np.array([[1e-3, 2e-3]])
    for phase, expect_symmetric in ((rocksalt_control, True), (gallium_arsenide, False)):
        zone = ZoneAxis(indices=(1, 1, 1), phase=phase)
        beams = beam_set_for_zone(
            phase,
            zone,
            laue_zones=(0, 1, -1),
            max_index=18,
            g_max_inv_angstrom=3.2,
            max_excitation_error_inv_angstrom=0.008,
            convergence_semi_angle_mrad=8.0,
        )
        matrix = structure_matrix(beams, tilts)
        asymmetry = float(np.max(np.abs(matrix - np.swapaxes(matrix, 1, 2))))
        if expect_symmetric:
            assert asymmetry < 1e-14
        else:
            assert asymmetry > 1e-5


def test_normal_absorption_is_exactly_a_scalar(nickel: Phase) -> None:
    """``exp(-2 pi t / xi'_0)`` multiplies everything and changes no relative intensity.

    The mean absorptive coefficient sits on every diagonal element of ``A``, so
    it commutes out of the matrix exponential. Dividing it out must recover the
    absorption-free result exactly when the anomalous term is switched off —
    which is the statement that the phenomenological ``mean_ratio`` cannot
    contaminate any conclusion about shape, position, or symmetry.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    tilts = np.random.default_rng(1).uniform(-6e-3, 6e-3, size=(24, 2))

    free = solve_bloch_waves(beams, tilts, thickness_angstrom=900.0)
    normal_only = solve_bloch_waves(
        beams,
        tilts,
        thickness_angstrom=900.0,
        absorption=AbsorptionModel(mean_ratio=0.12, reflection_ratio=0.0),
    )
    factor = normal_only.normal_absorption_factor
    assert 0.0 < factor < 1.0
    assert normal_only.intensities / factor == pytest.approx(free.intensities, rel=1e-10)


def test_absorption_makes_bright_field_asymmetric_and_leaves_dark_field_symmetric(
    nickel: Phase,
) -> None:
    """The Hashimoto-Howie-Whelan theorem, as a numerical assertion.

    In a two-beam calculation with absorption the dark-field rocking curve stays
    an even function of the excitation error while the bright-field one does
    not. The asymmetry is anomalous absorption: the two Bloch waves are absorbed
    at different rates, and which of them is preferentially excited depends on
    the sign of ``s``. Without absorption both curves are even, so the pair of
    assertions isolates the effect.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_from_indices(nickel, zone, [[2, 2, 0]])
    targets = np.linspace(-0.03, 0.03, 121)
    tilts = _tilts_for_excitation_errors(beams, 1, targets)

    free = solve_bloch_waves(beams, tilts, thickness_angstrom=1500.0)
    assert free.transmitted_intensity == pytest.approx(
        free.transmitted_intensity[::-1], abs=1e-12
    )

    absorbed = solve_bloch_waves(
        beams, tilts, thickness_angstrom=1500.0, absorption=AbsorptionModel()
    )
    bright = absorbed.transmitted_intensity
    dark = absorbed.intensity_of([2, 2, 0])
    assert np.max(np.abs(dark - dark[::-1])) / dark.max() < 1e-10
    assert np.max(np.abs(bright - bright[::-1])) / bright.max() > 0.1


def test_absorption_removes_intensity_and_never_adds_it(nickel: Phase) -> None:
    """Absorption must be a loss at every tilt, thickness and beam."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    tilts = np.random.default_rng(2).uniform(-6e-3, 6e-3, size=(32, 2))
    for thickness in (200.0, 1200.0):
        solution = solve_bloch_waves(
            beams, tilts, thickness_angstrom=thickness, absorption=AbsorptionModel()
        )
        assert np.all(solution.total_intensity < 1.0)
        assert np.all(solution.total_intensity > 0.0)


def test_absorption_model_rejects_a_gaining_bloch_wave() -> None:
    """A reflection ratio above the mean ratio would make a Bloch wave grow."""

    with pytest.raises(ValueError, match="gains intensity"):
        AbsorptionModel(mean_ratio=0.05, reflection_ratio=0.2)
    with pytest.raises(ValueError, match="non-negative"):
        AbsorptionModel(mean_ratio=-0.1, reflection_ratio=0.0)
    assert AbsorptionModel.none().is_absorbing is False
    assert AbsorptionModel(mean_ratio=0.2, reflection_ratio=0.2).is_absorbing is True


# --------------------------------------------------------------------------- #
# Friedel's law: the propagator-symmetry theorem
# --------------------------------------------------------------------------- #


def _maximum_friedel_violation(phase: Phase, laue_zones: tuple[int, ...]) -> float:
    """Largest ``|I_g(theta) - I_-g(-theta)|`` over a random tilt sample."""

    zone = ZoneAxis(indices=(1, 1, 1), phase=phase)
    beams = beam_set_for_zone(
        phase,
        zone,
        laue_zones=laue_zones,
        max_index=18 if len(laue_zones) > 1 else 4,
        g_max_inv_angstrom=3.2 if len(laue_zones) > 1 else 1.6,
        max_excitation_error_inv_angstrom=0.008 if len(laue_zones) > 1 else 0.03,
        convergence_semi_angle_mrad=8.0,
    )
    tilts = np.random.default_rng(3).uniform(-8e-3, 8e-3, size=(24, 2))
    solution = solve_bloch_waves(
        beams, np.concatenate([tilts, -tilts]), thickness_angstrom=600.0
    )
    half = tilts.shape[0]
    worst = 0.0
    for index in range(1, beams.size):
        try:
            partner = beams.index_of(-beams.miller_indices[index])
        except KeyError:
            continue
        difference = np.abs(
            solution.intensities[:half, index] - solution.intensities[half:, partner]
        )
        worst = max(worst, float(np.max(difference)))
    return worst


def test_zolz_only_calculation_cannot_see_the_lack_of_a_centre(
    gallium_arsenide: Phase,
) -> None:
    """A projection calculation reports Friedel's law even for zincblende.

    Every zeroth-Laue-zone coefficient of zincblende down ``[111]`` is real, so
    the structure matrix is symmetric and the two discs of a ``+-g`` pair are
    identical to machine precision. This is not a bug: it is why a CBED symmetry
    determination that ignores higher-order Laue zones is worthless, and it is
    stated on `pytex.diffraction.dynamical.BeamSet.holz_mask` for that reason.
    """

    assert _maximum_friedel_violation(gallium_arsenide, (0,)) < 1e-12


def test_holz_interaction_reveals_the_missing_centre_of_symmetry(
    gallium_arsenide: Phase, rocksalt_control: Phase
) -> None:
    """With HOLZ beams admitted, only the non-centrosymmetric structure breaks Friedel.

    The two phases share a lattice parameter, a pair of species and a cell
    content; they differ only by where the second sublattice sits, and therefore
    only by the presence of a centre of symmetry. The Friedel violation
    separates them by more than three orders of magnitude, which is the
    measurement CBED point-group determination rests on.
    """

    polar = _maximum_friedel_violation(gallium_arsenide, (0, 1, -1))
    centric = _maximum_friedel_violation(rocksalt_control, (0, 1, -1))
    assert polar > 0.1
    assert centric < 1e-3
    assert polar > 100.0 * centric


# --------------------------------------------------------------------------- #
# Beam sets
# --------------------------------------------------------------------------- #


def test_beam_set_puts_the_transmitted_beam_first_and_sorts_by_g(nickel: Phase) -> None:
    """Determinism of the beam ordering, which every returned array is indexed by."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    assert np.array_equal(beams.miller_indices[0], np.zeros(3, dtype=np.int64))
    assert beams.index_of([0, 0, 0]) == 0
    magnitudes = beams.g_magnitude_inv_angstrom
    assert np.all(np.diff(magnitudes) >= -1e-12)
    assert np.all(beams.laue_zone == 0)
    assert not beams.holz_mask.any()


def test_convergence_cone_is_what_admits_holz_reflections(nickel: Phase) -> None:
    """A HOLZ beam is far from Bragg on axis and exactly at Bragg inside the disc.

    Selecting on the zero-tilt excitation error alone discards every HOLZ
    reflection, and with it every HOLZ line. The selection therefore uses the
    minimum of ``|s_g|`` over the illumination cone, and this test is what
    documents that the difference is not cosmetic.
    """

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    common = {
        "laue_zones": (0, 1),
        "max_index": 20,
        "g_max_inv_angstrom": 5.2,
        "max_excitation_error_inv_angstrom": 0.004,
    }
    parallel = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=0.0, **common)
    convergent = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=8.0, **common)
    assert int(np.count_nonzero(convergent.holz_mask)) > int(
        np.count_nonzero(parallel.holz_mask)
    )
    assert convergent.holz_mask.any()


def test_beam_set_from_indices_takes_exactly_what_it_is_given(nickel: Phase) -> None:
    """No enumeration, no cut-off: the beams named plus the transmitted one."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_from_indices(nickel, zone, [[2, 2, 0], [-2, -2, 0], [2, 2, 0]])
    assert beams.size == 3
    assert beams.index_of([-2, -2, 0]) == 2
    with pytest.raises(KeyError, match="not in this beam set"):
        beams.index_of([1, 1, 1])


def test_beam_set_requires_a_unit_cell(nickel: Phase) -> None:
    """A bare lattice has no potential coefficients, so it cannot be propagated."""

    bare = Phase(
        "bare",
        lattice=nickel.lattice,
        symmetry=nickel.symmetry,
        crystal_frame=nickel.crystal_frame,
    )
    with pytest.raises(ValueError, match="carries no unit cell"):
        beam_set_for_zone(bare, ZoneAxis(indices=(0, 0, 1), phase=bare))


# --------------------------------------------------------------------------- #
# Reporting surface
# --------------------------------------------------------------------------- #


def test_thickness_series_reuses_the_eigenbasis(nickel: Phase) -> None:
    """``intensities_at`` must agree with a fresh solve at the same thickness."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    tilts = np.random.default_rng(4).uniform(-6e-3, 6e-3, size=(16, 2))
    kept = solve_bloch_waves(
        beams,
        tilts,
        thickness_angstrom=500.0,
        absorption=AbsorptionModel(),
        keep_eigenbasis=True,
    )
    for thickness in (250.0, 1750.0):
        fresh = solve_bloch_waves(
            beams, tilts, thickness_angstrom=thickness, absorption=AbsorptionModel()
        )
        assert kept.intensities_at(thickness) == pytest.approx(fresh.intensities, rel=1e-9)
    series = kept.intensities_at(np.array([300.0, 600.0]))
    assert series.shape == (2, tilts.shape[0], beams.size)


def test_thickness_series_refuses_without_the_eigenbasis(nickel: Phase) -> None:
    """Silently re-solving would hide a cost the caller chose not to pay."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    solution = solve_bloch_waves(beams, [[0.0, 0.0]], thickness_angstrom=500.0)
    with pytest.raises(ValueError, match="keep_eigenbasis=True"):
        solution.intensities_at(700.0)


def test_describe_and_json_report_the_limits_and_stay_in_lockstep(nickel: Phase) -> None:
    """The explainable-results contract: prose and payload must agree."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    solution = solve_bloch_waves(
        beams,
        np.random.default_rng(5).uniform(-6e-3, 6e-3, size=(8, 2)),
        thickness_angstrom=800.0,
        absorption=AbsorptionModel(),
    )
    prose = solution.describe()
    assert "Bloch-wave" in prose
    assert "[001]" in prose
    assert "absorbed" in prose

    payload = solution.to_json_dict()
    assert payload["schema"] == BLOCH_WAVE_SOLUTION_SCHEMA
    assert payload["beam_count"] == beams.size
    assert payload["holz_beam_count"] == 0
    assert payload["absorption"]["mean_ratio"] == pytest.approx(0.1)
    assert len(payload["beams"]) == beams.size

    assert "projection" in beams.describe()
    assert "phenomenological" in AbsorptionModel().describe()
    assert "unitary" in AbsorptionModel.none().describe()


def test_solver_rejects_a_non_positive_thickness(nickel: Phase) -> None:
    """A zero-thickness foil is not a limiting case worth silently accepting."""

    zone = ZoneAxis(indices=(0, 0, 1), phase=nickel)
    beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
    with pytest.raises(ValueError, match="strictly positive"):
        solve_bloch_waves(beams, [[0.0, 0.0]], thickness_angstrom=0.0)
