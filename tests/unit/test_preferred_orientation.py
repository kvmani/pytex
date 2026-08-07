"""Preferred-orientation corrections for powder intensities.

The expected values are analytic, not recorded outputs. Four independent anchors
are used:

1. **The exact normalization of the March distribution.** Its average over a
   uniform distribution of directions is 1 for every March coefficient, because

   .. math:: \\int_0^1 \\bigl((r^2 - r^{-1})u^2 + r^{-1}\\bigr)^{-3/2} du = 1 .

   This is checked by equal-area quadrature on the sphere, and it is the
   statement that texture *redistributes* intensity rather than creating it.
2. **Closed-form limits**: ``P(0) = r^{-3}`` and ``P(pi/2) = r^{3/2}``.
3. **The random-texture identity**: a uniform ODF has pole density 1 in every
   direction, so the ODF-weighted correction must be 1 for every reflection.
4. **A planted single-crystal-like texture**: the pole that was planted along
   the scattering direction must be the enhanced one.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import CRYSTAL_FRAME, SPECIMEN_FRAME
from pytex.core.lattice import Lattice, Phase
from pytex.core.miller import MillerPlane
from pytex.core.orientation import OrientationSet
from pytex.core.sphere import S2Grid
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.preferred_orientation import (
    MarchDollaseModel,
    ODFPreferredOrientationModel,
    PreferredOrientationModel,
    march_dollase_factors,
    preferred_orientation_factor_table,
)
from pytex.diffraction.xrd import (
    apply_preferred_orientation,
    generate_powder_reflections,
    generate_xrd_pattern,
)
from pytex.texture.models import ODF, KernelSpec

_NICKEL_A = 3.52387


def _cubic_phase() -> Phase:
    lattice = Lattice(
        a=_NICKEL_A,
        b=_NICKEL_A,
        c=_NICKEL_A,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=CRYSTAL_FRAME,
    )
    return Phase(
        name="test-cubic",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
        space_group_symbol="Fm-3m",
        space_group_number=225,
    )


# --------------------------------------------------------------------------- #
# The March distribution
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("march_coefficient", [0.4, 0.7, 1.0, 1.5, 2.5])
def test_march_distribution_averages_to_one_over_the_sphere(march_coefficient: float) -> None:
    """The exact normalization identity, by equal-area quadrature.

    Preferred orientation redistributes diffracted intensity; it does not create
    or destroy it. That statement *is* this integral, so it is the primary
    anchor for the whole model.
    """

    grid = S2Grid.equispaced(2.0, reference_frame=SPECIMEN_FRAME, hemisphere="sphere")
    axis = np.array([0.0, 0.0, 1.0])
    angles = np.arccos(np.clip(np.asarray(grid.vectors.values) @ axis, -1.0, 1.0))
    factors = march_dollase_factors(angles, march_coefficient)
    weights = np.asarray(grid.weights, dtype=np.float64)
    average = float(np.sum(factors * weights) / np.sum(weights))
    assert average == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize("march_coefficient", [0.25, 0.5, 0.9, 1.0, 1.3, 3.0])
def test_march_distribution_matches_its_closed_form_limits(march_coefficient: float) -> None:
    """``P(0) = r^-3`` and ``P(pi/2) = r^(3/2)``, exactly."""

    values = march_dollase_factors([0.0, 0.5 * np.pi], march_coefficient)
    assert values[0] == pytest.approx(march_coefficient**-3.0, rel=1e-12)
    assert values[1] == pytest.approx(march_coefficient**1.5, rel=1e-12)


def test_random_march_coefficient_is_exactly_the_identity() -> None:
    """``r = 1`` must give exactly 1 at every angle, with no numerical drift."""

    angles = np.linspace(0.0, np.pi, 257)
    assert np.all(march_dollase_factors(angles, 1.0) == 1.0)


def test_plate_and_needle_habits_act_in_opposite_directions() -> None:
    """``r < 1`` enhances the preferred pole; ``r > 1`` suppresses it."""

    plate = march_dollase_factors([0.0], 0.6)
    needle = march_dollase_factors([0.0], 1.6)
    assert float(plate[0]) > 1.0
    assert float(needle[0]) < 1.0


def test_march_distribution_rejects_a_non_positive_coefficient() -> None:
    with pytest.raises(ValueError, match="march_coefficient must be finite"):
        march_dollase_factors(0.0, 0.0)
    with pytest.raises(ValueError, match="march_coefficient must be finite"):
        march_dollase_factors(0.0, -1.0)


# --------------------------------------------------------------------------- #
# MarchDollaseModel
# --------------------------------------------------------------------------- #


def test_random_model_leaves_every_reflection_untouched() -> None:
    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=1.0,
    )
    assert model.is_random
    planes = [
        MillerPlane.from_hkl(indices, phase=phase)
        for indices in ([1, 1, 1], [0, 0, 2], [0, 2, 2], [1, 1, 3])
    ]
    assert np.all(model.factors(planes) == 1.0)


def test_preferred_reflection_is_the_one_enhanced_by_a_plate_texture() -> None:
    """A plate texture on (111) must enhance {111} and suppress something else.

    The {111} family contains the preferred pole itself, so it always has a
    member at zero angle; a family that does not is pushed the other way. This
    checks the sign of the physics, not a recorded number.
    """

    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.6,
    )
    preferred = MillerPlane.from_hkl([1, 1, 1], phase=phase)
    other = MillerPlane.from_hkl([0, 0, 2], phase=phase)
    factors = model.factors([preferred, other])
    assert factors[0] > 1.0
    assert factors[1] < 1.0


def test_family_averaging_makes_the_factor_representative_independent() -> None:
    """Any member of a symmetry family must receive the same factor.

    Which member of ``{200}`` the enumeration happened to emit is an artefact,
    so a correction that depended on it would be wrong. Family averaging is what
    removes that dependence.
    """

    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([0, 0, 1], phase=phase),
        march_coefficient=0.7,
    )
    representatives = [
        MillerPlane.from_hkl(indices, phase=phase)
        for indices in ([2, 0, 0], [0, 2, 0], [0, 0, 2], [-2, 0, 0])
    ]
    factors = model.factors(representatives)
    assert np.allclose(factors, factors[0], rtol=0.0, atol=1e-12)


def test_family_averaged_factors_respect_the_normalization_on_average() -> None:
    """Averaged over many families the correction stays near unity.

    A model that systematically scaled every reflection up or down would be
    creating intensity rather than redistributing it. Over a spread of families
    the multiplicity-weighted mean must stay close to 1.
    """

    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 0, 0], phase=phase),
        march_coefficient=0.75,
    )
    reflections = generate_powder_reflections(phase, two_theta_range_deg=(20.0, 150.0))
    planes = [MillerPlane(indices=r.miller_indices, phase=phase) for r in reflections]
    factors = np.asarray(model.factors(planes))
    multiplicities = np.array([r.multiplicity for r in reflections], dtype=np.float64)
    weighted = float(np.sum(factors * multiplicities) / np.sum(multiplicities))
    assert 0.8 < weighted < 1.25


def test_model_rejects_planes_from_another_phase() -> None:
    phase = _cubic_phase()
    other = Phase(
        name="other",
        lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=CRYSTAL_FRAME),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
    )
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.8,
    )
    with pytest.raises(ValueError, match="share the preferred-orientation phase"):
        model.factors([MillerPlane.from_hkl([1, 1, 1], phase=other)])


def test_model_rejects_a_non_positive_coefficient() -> None:
    phase = _cubic_phase()
    with pytest.raises(ValueError, match="march_coefficient must be finite"):
        MarchDollaseModel(
            preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
            march_coefficient=0.0,
        )


def test_march_dollase_model_satisfies_the_protocol() -> None:
    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.8,
    )
    assert isinstance(model, PreferredOrientationModel)


# --------------------------------------------------------------------------- #
# ODF-weighted correction
# --------------------------------------------------------------------------- #


def _uniform_odf(phase: Phase, *, count: int = 4000) -> ODF:
    orientations = OrientationSet.from_equispaced_so3_grid(
        12.0,
        specimen_frame=SPECIMEN_FRAME,
        phase=phase,
        reduce_to_fundamental_region=False,
    )
    return ODF.from_orientations(orientations, kernel=KernelSpec(halfwidth_deg=15.0))


def test_a_uniform_texture_gives_no_correction() -> None:
    """The random-texture identity: pole density is 1 everywhere, so factors are 1.

    This is the anchor that makes the ODF correction interpretable — it must
    reduce to the random-powder case that the uncorrected pattern assumes.
    """

    phase = _cubic_phase()
    model = ODFPreferredOrientationModel(odf=_uniform_odf(phase))
    planes = [
        MillerPlane.from_hkl(indices, phase=phase)
        for indices in ([1, 1, 1], [0, 0, 2], [0, 2, 2], [1, 1, 3])
    ]
    factors = np.asarray(model.factors(planes))
    assert np.allclose(factors, 1.0, atol=0.06), factors


def test_a_planted_fibre_texture_enhances_the_pole_it_was_planted_on() -> None:
    """Plant {001} along ND and {001} must be the enhanced reflection.

    Built through the canonical constructors rather than through the code under
    test, so the expectation does not depend on the correction's own algebra.
    """

    phase = _cubic_phase()
    # A fibre about ND with (001) along ND: rotations about the specimen normal
    # leave the (001) pole on ND, so the {001} pole density there is high.
    orientations = OrientationSet.from_so2_grid(
        "ND", 5.0, specimen_frame=SPECIMEN_FRAME, phase=phase
    )
    odf = ODF.from_orientations(orientations, kernel=KernelSpec(halfwidth_deg=10.0))
    model = ODFPreferredOrientationModel(odf=odf)

    planted = MillerPlane.from_hkl([0, 0, 2], phase=phase)
    away = MillerPlane.from_hkl([1, 1, 1], phase=phase)
    factors = model.factors([planted, away])
    assert factors[0] > factors[1]
    assert factors[0] > 1.0


def test_odf_model_defaults_to_the_specimen_normal_and_normalizes_it() -> None:
    phase = _cubic_phase()
    model = ODFPreferredOrientationModel(odf=_uniform_odf(phase))
    assert np.allclose(model.scattering_direction, (0.0, 0.0, 1.0))

    tilted = ODFPreferredOrientationModel(
        odf=_uniform_odf(phase), scattering_direction=np.array([0.0, 0.0, 5.0])
    )
    assert np.allclose(tilted.scattering_direction, (0.0, 0.0, 1.0))


def test_odf_model_satisfies_the_protocol_and_names_its_frame() -> None:
    phase = _cubic_phase()
    model = ODFPreferredOrientationModel(odf=_uniform_odf(phase))
    assert isinstance(model, PreferredOrientationModel)
    assert model.specimen_frame_name == SPECIMEN_FRAME.name


# --------------------------------------------------------------------------- #
# Applying a correction to a pattern
# --------------------------------------------------------------------------- #


def test_applying_a_random_model_changes_no_intensity() -> None:
    phase = _cubic_phase()
    reflections = generate_powder_reflections(phase, two_theta_range_deg=(20.0, 120.0))
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=1.0,
    )
    corrected = apply_preferred_orientation(reflections, model, phase=phase)
    assert [r.intensity for r in corrected] == [r.intensity for r in reflections]


def test_correction_preserves_every_field_except_intensity() -> None:
    """A texture correction is about intensity; nothing else may move."""

    phase = _cubic_phase()
    reflections = generate_powder_reflections(phase, two_theta_range_deg=(20.0, 120.0))
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.55,
    )
    corrected = apply_preferred_orientation(reflections, model, phase=phase)
    assert len(corrected) == len(reflections)
    for before, after in zip(reflections, corrected, strict=True):
        assert np.array_equal(after.miller_indices, before.miller_indices)
        assert after.two_theta_deg == before.two_theta_deg
        assert after.d_spacing_angstrom == before.d_spacing_angstrom
        assert after.multiplicity == before.multiplicity
        assert after.structure_factor_amplitude == before.structure_factor_amplitude
    assert any(a.intensity != b.intensity for a, b in zip(corrected, reflections, strict=True))


def test_empty_reflection_lists_are_handled() -> None:
    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.8,
    )
    assert apply_preferred_orientation((), model, phase=phase) == ()


def test_pattern_generation_accepts_and_applies_a_model() -> None:
    """The pattern and its reflection list must both carry the correction."""

    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.5,
    )
    random_pattern = generate_xrd_pattern(phase, two_theta_range_deg=(20.0, 120.0))
    textured = generate_xrd_pattern(
        phase, two_theta_range_deg=(20.0, 120.0), preferred_orientation=model
    )

    random_by_hkl = {
        tuple(int(v) for v in r.miller_indices): r.intensity for r in random_pattern.reflections
    }
    textured_by_hkl = {
        tuple(int(v) for v in r.miller_indices): r.intensity for r in textured.reflections
    }
    assert random_by_hkl.keys() == textured_by_hkl.keys()

    # {111} is enhanced by a (111) plate texture, so its share of the pattern
    # rises relative to {200}. The family representative the enumeration emits
    # is an artefact, so both families are located by their index multiset.
    def _family(table: dict[tuple[int, int, int], float], want: tuple[int, ...]) -> float:
        key = next(k for k in table if tuple(sorted(abs(v) for v in k)) == want)
        return table[key]

    ratio_random = _family(random_by_hkl, (1, 1, 1)) / _family(random_by_hkl, (0, 0, 2))
    ratio_textured = _family(textured_by_hkl, (1, 1, 1)) / _family(textured_by_hkl, (0, 0, 2))
    assert ratio_textured > ratio_random

    # Peak positions must not move: texture changes intensities, not angles.
    assert np.allclose(random_pattern.two_theta_grid_deg, textured.two_theta_grid_deg)


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #


def test_factor_table_labels_reflections_as_families() -> None:
    """A powder reflection is a family, so it takes ``{hkl}`` brackets."""

    phase = _cubic_phase()
    model = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([1, 1, 1], phase=phase),
        march_coefficient=0.7,
    )
    planes = [
        MillerPlane.from_hkl(indices, phase=phase) for indices in ([1, 1, 1], [0, 0, 2])
    ]
    table = preferred_orientation_factor_table(planes, model)
    assert [label for label, _ in table] == ["{111}", "{002}"]
    assert table[0][1] > 1.0


def test_march_describe_states_the_habit_and_the_fibre_assumption() -> None:
    phase = _cubic_phase()
    plate = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([0, 0, 1], phase=phase),
        march_coefficient=0.6,
    ).describe()
    assert "plate-like" in plate
    assert "{001}" in plate
    assert "fibre texture" in plate
    assert "multiples of the random-powder intensity" in plate

    needle = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([0, 0, 1], phase=phase),
        march_coefficient=1.7,
    ).describe()
    assert "needle-like" in needle

    random_text = MarchDollaseModel(
        preferred_orientation=MillerPlane.from_hkl([0, 0, 1], phase=phase),
        march_coefficient=1.0,
    ).describe()
    assert "random powder" in random_text


def test_odf_describe_states_the_direction_frame_and_geometry_assumption() -> None:
    phase = _cubic_phase()
    text = ODFPreferredOrientationModel(odf=_uniform_odf(phase)).describe()
    assert "pole density" in text
    assert "multiples of a random distribution" in text
    assert "Bragg-Brentano" in text
    assert SPECIMEN_FRAME.name in text
