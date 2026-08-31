"""Double diffraction in the zone-law SAED engine: which spots appear, and why.

Every assertion here has a source outside the code under test.

* **Silicon down [011].** ``{200}`` is forbidden by the diamond glide and is the
  textbook illustration of double diffraction, produced by ``111 + (-1)11``.
  Williams & Carter, *Transmission Electron Microscopy*, 2nd ed., ch. 16.
* **hcp down [11-20].** ``(0001)`` is forbidden by the ``6_3`` screw axis and
  appears through ``(-1101) + (1-100)``. Edington, *Electron Diffraction in the
  Electron Microscope*, ch. 2.
* **Centring absences are never revived.** Lattice-centring conditions define a
  *sublattice* of reciprocal space, and a sublattice is closed under addition,
  so no sum of two centring-allowed reflections can land on a centring absence.
  This is a group-theoretic fact, not a property of the implementation, and bcc
  and fcc are checked against it rather than against stored numbers.

The intensity of a doubly diffracted spot is a model choice and is declared as
indicative, so nothing here pins its value. What is pinned is the *selection* —
which reflections appear, which do not, and that the marking says so.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.core.lattice import ZoneAxis
from pytex.diffraction.saed import generate_saed_pattern


def pattern(phase_id: str, zone: tuple[int, int, int], **overrides):
    phase = builtin_phase(phase_id).to_phase()
    options = {"max_index": 3, "include_double_diffraction": True}
    options.update(overrides)
    return generate_saed_pattern(phase, ZoneAxis(indices=zone, phase=phase), **options)


def marked(result) -> set[tuple[int, ...]]:
    return {
        tuple(int(value) for value in spot.miller_indices)
        for spot in result.spots
        if spot.is_double_diffraction
    }


def test_silicon_down_011_revives_the_forbidden_200_reflection() -> None:
    """The textbook case: {200} is absent kinematically and present on a plate."""

    assert (2, 0, 0) not in marked(
        pattern("si_diamond", (0, 1, -1), include_double_diffraction=False)
    )
    revived = marked(pattern("si_diamond", (0, 1, -1)))
    assert (2, 0, 0) in revived
    assert (-2, 0, 0) in revived


def test_the_revived_silicon_reflection_names_a_pair_of_real_reflections() -> None:
    """The parents must be genuine 111-type reflections that sum to the spot."""

    spot = next(
        spot
        for spot in pattern("si_diamond", (0, 1, -1)).spots
        if tuple(int(value) for value in spot.miller_indices) == (2, 0, 0)
    )
    parents = np.asarray(spot.double_diffraction_parents)
    assert parents.shape == (2, 3)
    assert tuple(int(value) for value in parents.sum(axis=0)) == (2, 0, 0)
    # Each parent is a {111} reflection, which is what the diamond lattice
    # actually has strongly excited in this zone.
    assert all(sorted(abs(int(value)) for value in parent) == [1, 1, 1] for parent in parents)
    assert spot.double_diffraction_origin_label() == "(1 -1 -1) + (111)"


def test_hcp_down_11_2bar_0_revives_the_forbidden_0001_reflection() -> None:
    """The 6_3 screw axis forbids (0001); a real plate shows it anyway.

    The zone is written [110] in the three-index basis the engine takes, which
    is the four-index [11-20] direction.
    """

    for phase_id in ("ti_hcp", "zr_hcp", "mg_hcp"):
        revived = marked(pattern(phase_id, (1, 1, 0)))
        assert (0, 0, 1) in revived, phase_id
        assert (0, 0, -1) in revived, phase_id


def test_a_centring_absence_is_never_revived() -> None:
    """A centred reciprocal lattice is closed under addition, so nothing appears.

    bcc and fcc absences are *centring* absences, and this is the property that
    distinguishes them from the basis absences of diamond and hcp. If this ever
    fails, the selection rule has stopped being a sum over the excited set.
    """

    assert marked(pattern("fe_bcc", (0, 0, 1))) == set()
    assert marked(pattern("fe_bcc", (1, 1, 0))) == set()
    assert marked(pattern("ni_fcc", (0, 1, -1))) == set()
    assert marked(pattern("ni_fcc", (0, 0, 1))) == set()


def test_the_option_moves_no_spot_and_relabels_no_allowed_reflection() -> None:
    """Switching it on is additive: the kinematic pattern inside it is unchanged.

    This is what makes the toggle a fair comparison. If enabling the option also
    nudged positions or re-weighted allowed reflections, a user flipping it could
    not attribute a difference to double diffraction.
    """

    off = pattern("si_diamond", (0, 1, -1), include_double_diffraction=False)
    on = pattern("si_diamond", (0, 1, -1))
    by_index = {tuple(int(value) for value in spot.miller_indices): spot for spot in on.spots}
    assert len(off.spots) == len(on.spots)
    for spot in off.spots:
        key = tuple(int(value) for value in spot.miller_indices)
        partner = by_index[key]
        assert np.allclose(partner.detector_coordinates, spot.detector_coordinates)
        if not partner.is_double_diffraction:
            assert partner.intensity == pytest.approx(spot.intensity)


def test_a_marked_reflection_is_weaker_than_every_genuine_one() -> None:
    """An observability estimate must never outrank a real reflection.

    A doubly diffracted spot on a real plate is faint. If the coupling model ever
    put one above a genuine reflection, the pattern would teach the opposite of
    what a plate shows.
    """

    result = pattern("si_diamond", (0, 1, -1))
    genuine = [spot.intensity for spot in result.spots if not spot.is_double_diffraction]
    forbidden = [spot.intensity for spot in result.spots if spot.is_double_diffraction]
    assert forbidden
    assert max(forbidden) < max(genuine)


def test_the_coupling_scales_the_estimate_without_changing_the_selection() -> None:
    """Coupling sets how prominent these spots look, not whether they are there."""

    weak = pattern("si_diamond", (0, 1, -1), double_diffraction_coupling=0.01)
    strong = pattern("si_diamond", (0, 1, -1), double_diffraction_coupling=0.5)
    assert marked(weak) == marked(strong)

    def intensity_of(result, hkl):
        return next(
            spot.intensity
            for spot in result.spots
            if tuple(int(value) for value in spot.miller_indices) == hkl
        )

    assert intensity_of(strong, (2, 0, 0)) > intensity_of(weak, (2, 0, 0))


def test_an_out_of_range_coupling_is_refused() -> None:
    phase = builtin_phase("si_diamond").to_phase()
    axis = ZoneAxis(indices=(0, 1, -1), phase=phase)
    for bad in (0.0, -0.1, 1.5, float("nan")):
        with pytest.raises(ValueError, match="double_diffraction_coupling"):
            generate_saed_pattern(
                phase, axis, include_double_diffraction=True, double_diffraction_coupling=bad
            )


def test_parents_that_do_not_sum_to_the_reflection_are_refused() -> None:
    """A construction-time invariant, per the repository's preference for them."""

    from pytex.diffraction.saed import SAEDSpot

    with pytest.raises(ValueError, match="must sum to the reflection"):
        SAEDSpot(
            miller_indices=(2, 0, 0),
            reciprocal_vector_crystal=(1.0, 0.0, 0.0),
            reciprocal_vector_detector=(1.0, 0.0, 0.0),
            detector_coordinates=(1.0, 0.0),
            intensity=0.1,
            excitation_error_inv_angstrom=0.0,
            double_diffraction_parents=((1, 1, 1), (1, 1, 1)),
        )
