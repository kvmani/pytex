"""Validation of the indexed-pattern to crystal-orientation bridge.

The tests are built by *inverting* the composition: given a known orientation,
stage position and diffraction rotation, compute the crystal-to-pattern rotation
that indexing would have reported, then check the bridge recovers the
orientation. That direction matters — it exercises the composition against an
independently constructed input rather than against its own output.

Two cases carry most of the weight. The self-calibration test recovers a planted
diffraction rotation from the patterns alone, which is the claim that makes the
multi-pattern path worth having. And the mirrored-storage test uses the
*physical* manifestation of a mirrored pattern — an in-plane mirror together
with a reversed pattern normal, which is a proper rotation — rather than an
improper matrix no real solver would emit.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.core._angles import rotation_angle_from_matrix_rad
from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frame_catalog import sample_frame
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice, Phase, ZoneAxis
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.tem.indexing import (
    INDEXED_ORIENTATION_SCHEMA,
    IndexedPatternObservation,
    orientation_from_indexed_pattern,
    orientation_from_indexed_patterns,
    orientations_from_pattern_report,
)
from pytex.tem.navigation import solve_tilts_for_direction
from pytex.tem.reconstruction import HOLDER_FRAME
from pytex.tem.stage import (
    DoubleTiltStage,
    RectangularEnvelope,
    StageCalibration,
    StagePosition,
    rotation_x,
    rotation_z,
)

WIDE = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def nickel(point_group: str = "m-3m") -> Phase:
    frame = _crystal_frame()
    return Phase(
        "nickel-fcc",
        lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=frame),
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


def random_orientations(count: int, seed: int) -> list[np.ndarray]:
    generator = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        matrix = np.linalg.qr(generator.normal(size=(3, 3)))[0]
        if np.linalg.det(matrix) < 0:
            matrix[:, 0] *= -1
        out.append(matrix)
    return out


def pattern_rotation_for(
    crystal_to_holder: np.ndarray,
    position: StagePosition,
    diffraction_rotation_deg: float,
    stage: DoubleTiltStage = WIDE,
) -> np.ndarray:
    """The crystal-to-pattern rotation indexing would report, by inversion.

    From ``U = R_stage^T Rz(phi) R``, so ``R = Rz(-phi) R_stage U``.
    """

    stage_matrix = stage.rotation_matrix(position.alpha_deg, position.beta_deg)
    return rotation_z(math.radians(-diffraction_rotation_deg)) @ stage_matrix @ crystal_to_holder


def position_for(crystal_to_holder: np.ndarray, zone: ZoneAxis) -> StagePosition:
    """The stage position that puts ``zone`` on the beam under this orientation."""

    return StagePosition(*solve_tilts_for_direction(crystal_to_holder @ zone.unit_vector)[0])


def symmetry_reduced_angle_deg(
    first: np.ndarray, second: np.ndarray, phase: Phase
) -> float:
    """Smallest angle between two orientations under crystal symmetry.

    The angle comes from the skew part of each relative rotation rather than
    from its trace: this helper exists to certify that a recovery is *exact*,
    and the trace form cannot report better than about ``1e-06`` degrees for a
    near-identity rotation whatever the recovery actually achieved. See
    `pytex.core._angles`.
    """

    operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
    relative = np.einsum(
        "ij,njk->nik", first, np.einsum("nij,jk->nik", operators, second.T)
    )
    angles = np.degrees(rotation_angle_from_matrix_rad(relative))
    return float(np.min(angles))


class TestSinglePattern:
    def test_recovers_a_known_orientation_exactly(self) -> None:
        """The whole point: indexing plus stage angles gives back the orientation."""

        phase = nickel()
        zone = ZoneAxis([0, 0, 1], phase=phase)
        generator = np.random.default_rng(5)
        worst = 0.0
        for truth in random_orientations(200, seed=5):
            rotation_deg = float(generator.uniform(-180.0, 180.0))
            stage = DoubleTiltStage(
                envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
                calibration=StageCalibration(diffraction_rotation_deg=rotation_deg),
            )
            position = position_for(truth, zone)
            result = orientation_from_indexed_pattern(
                pattern_rotation_for(truth, position, rotation_deg),
                zone,
                position,
                stage,
            )
            worst = max(worst, float(np.max(np.abs(result.matrix - truth))))
        assert worst < 1e-12

    def test_identity_case_pins_the_sign_conventions(self) -> None:
        """[001] on the beam at zero tilt, zero rotation, identity pattern.

        The trivial case, and the one that fixes every sign in the composition:
        anything transposed or mis-ordered moves the answer off the identity.
        """

        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        result = orientation_from_indexed_pattern(
            np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), stage
        )
        assert np.allclose(result.matrix, np.eye(3), atol=1e-12)
        assert result.euler_bunge_deg == pytest.approx((0.0, 0.0, 0.0), abs=1e-9)

    def test_forward_validation_places_the_zone_axis_on_the_beam(self) -> None:
        phase = nickel()
        zone = ZoneAxis([1, 1, 2], phase=phase)
        for truth in random_orientations(30, seed=9):
            stage = DoubleTiltStage(
                envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
                calibration=StageCalibration(diffraction_rotation_deg=15.0),
            )
            position = position_for(truth, zone)
            result = orientation_from_indexed_pattern(
                pattern_rotation_for(truth, position, 15.0), zone, position, stage
            )
            assert result.zone_axis_residual_deg < 1e-9

    def test_refuses_without_a_diffraction_rotation_calibration(self) -> None:
        """Guessing that constant rotates the orientation bodily about the beam."""

        phase = nickel()
        with pytest.raises(ValueError, match="calibrated diffraction rotation"):
            orientation_from_indexed_pattern(
                np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), WIDE
            )

    def test_refuses_an_improper_pattern_matrix(self) -> None:
        """Indexing never returns a reflection; one here is a caller error."""

        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        with pytest.raises(ValueError, match="improper"):
            orientation_from_indexed_pattern(
                np.diag([1.0, -1.0, 1.0]),
                ZoneAxis([0, 0, 1], phase=phase),
                StagePosition(0.0, 0.0),
                stage,
            )

    def test_wrong_parity_makes_the_orientation_improper(self) -> None:
        phase = nickel()
        stage = DoubleTiltStage(
            calibration=StageCalibration(
                diffraction_rotation_deg=0.0, pattern_is_mirrored=True
            )
        )
        with pytest.raises(ValueError, match="improper"):
            orientation_from_indexed_pattern(
                np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), stage
            )

    def test_a_rotation_error_rotates_the_orientation_about_the_beam(self) -> None:
        """The error is not absorbed by symmetry; it is a bodily rotation.

        This is why the single-pattern path insists on a calibration: the wrong
        constant yields a clean, self-consistent, wrong orientation.
        """

        phase = nickel()
        zone = ZoneAxis([0, 0, 1], phase=phase)
        position = StagePosition(0.0, 0.0)
        # A rotation about the crystal c axis keeps [001] on the beam at zero
        # tilt, so the beam axis and the holder z axis coincide and the error's
        # axis can be checked directly.
        truth = np.asarray(
            Rotation.from_axis_angle([0.0, 0.0, 1.0], 0.3).as_matrix(), dtype=np.float64
        )
        pattern = pattern_rotation_for(truth, position, 0.0)
        for error_deg in (10.0, 180.0):
            stage = DoubleTiltStage(
                calibration=StageCalibration(diffraction_rotation_deg=error_deg)
            )
            result = orientation_from_indexed_pattern(pattern, zone, position, stage)
            relative = result.matrix @ truth.T
            angle = math.degrees(
                math.acos(max(-1.0, min(1.0, (float(np.trace(relative)) - 1.0) / 2.0)))
            )
            assert angle == pytest.approx(error_deg, abs=1e-9)
            # and the axis of that error is the beam axis
            assert np.allclose(relative @ np.array([0.0, 0.0, 1.0]), [0.0, 0.0, 1.0], atol=1e-9)

    def test_euler_angles_are_bunge_and_round_trip(self) -> None:
        phase = nickel()
        stage = DoubleTiltStage(
            envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
            calibration=StageCalibration(diffraction_rotation_deg=20.0),
        )
        zone = ZoneAxis([0, 1, 1], phase=phase)
        truth = random_orientations(1, seed=17)[0]
        position = position_for(truth, zone)
        result = orientation_from_indexed_pattern(
            pattern_rotation_for(truth, position, 20.0), zone, position, stage
        )
        phi1, capital_phi, phi2 = result.euler_bunge_deg
        rebuilt = Rotation.from_bunge_euler(phi1, capital_phi, phi2, degrees=True)
        assert np.allclose(rebuilt.as_matrix(), result.matrix, atol=1e-9)

    def test_reports_the_ambiguity_and_stays_quiet_when_there_is_none(self) -> None:
        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        result = orientation_from_indexed_pattern(
            np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), stage
        )
        assert result.is_unique
        assert result.equivalent_orientations == ()

    def test_non_centrosymmetric_phase_reports_the_alternative(self) -> None:
        """GaAs down [001] genuinely admits two orientations; both are returned."""

        phase = nickel("-43m")
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        result = orientation_from_indexed_pattern(
            np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), stage
        )
        assert not result.is_unique
        assert len(result.equivalent_orientations) == 1
        assert "hypotheses" in result.describe()

    def test_can_be_declared_against_a_named_sample_frame(self) -> None:
        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        rd_td_nd = sample_frame()
        result = orientation_from_indexed_pattern(
            np.eye(3),
            ZoneAxis([0, 0, 1], phase=phase),
            StagePosition(0.0, 0.0),
            stage,
            specimen_frame=rd_td_nd,
        )
        assert result.orientation.specimen_frame == rd_td_nd
        # Relabelling must not rotate anything.
        assert np.allclose(result.matrix, np.eye(3), atol=1e-12)
        assert np.allclose(result.in_frame(HOLDER_FRAME).as_matrix(), result.matrix)

    def test_describe_and_json_stay_in_lockstep(self) -> None:
        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=12.0))
        result = orientation_from_indexed_pattern(
            np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0), stage
        )
        text = result.describe()
        payload = result.to_json_dict()
        assert "Bunge" in text
        assert payload["schema"] == INDEXED_ORIENTATION_SCHEMA
        assert payload["diffraction_rotation_deg"] == 12.0
        assert payload["euler_bunge_deg"]["phi1"] == pytest.approx(
            result.euler_bunge_deg[0]
        )
        assert payload["is_unique"] is result.is_unique


class TestMultiPatternSelfCalibration:
    @pytest.mark.parametrize("rotation_deg", [0.0, 37.0, -110.0, 175.0])
    def test_recovers_orientation_and_rotation_together(self, rotation_deg: float) -> None:
        """The claim that makes this path worth having: no calibration needed.

        The diffraction rotation is *determined* from the patterns, so a stage
        carrying no calibration at all still yields both the orientation and the
        instrument constant.
        """

        phase = nickel()
        truth = random_orientations(1, seed=23)[0]
        observations = []
        for indices in ([0, 0, 1], [0, 1, 1], [1, 1, 1]):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    pattern_rotation_for(truth, position, rotation_deg),
                    zone,
                    position,
                    label=str(indices),
                )
            )
        result = orientation_from_indexed_patterns(observations, WIDE)

        assert result.is_consistent
        assert symmetry_reduced_angle_deg(result.matrix, truth, phase) < 1e-9
        wrapped = ((result.diffraction_rotation_deg - rotation_deg + 180.0) % 360.0) - 180.0
        assert wrapped == pytest.approx(0.0, abs=1e-9)
        assert result.diffraction_rotation_scatter_deg < 1e-9
        assert result.beam_deviation_deg < 1e-9
        assert result.interzonal_residual_deg < 1e-9

    def test_two_patterns_are_enough(self) -> None:
        phase = nickel()
        truth = random_orientations(1, seed=29)[0]
        observations = []
        for indices in ([0, 0, 1], [0, 1, 1]):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    pattern_rotation_for(truth, position, 42.0), zone, position
                )
            )
        result = orientation_from_indexed_patterns(observations, WIDE)
        assert result.is_consistent
        assert symmetry_reduced_angle_deg(result.matrix, truth, phase) < 1e-9
        assert result.diffraction_rotation_deg == pytest.approx(42.0, abs=1e-9)

    def test_mirrored_storage_is_detected_by_the_beam_deviation(self) -> None:
        """The physical manifestation of a mirrored pattern, not an improper matrix.

        Flipping the stored y axis mirrors the in-plane pattern *and* reverses the
        derived pattern normal. That pair is a 180-degree rotation about the
        pattern x axis — proper, and therefore exactly what a solver building
        right-handed triads would return. It is caught by the residual component
        that no diffraction rotation can absorb.
        """

        phase = nickel()
        truth = random_orientations(1, seed=2)[0]
        mirror = rotation_x(math.pi)
        observations = []
        for indices in ([0, 0, 1], [0, 1, 1], [1, 1, 1]):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    mirror @ pattern_rotation_for(truth, position, 25.0),
                    ZoneAxis([-v for v in indices], phase=phase),
                    position,
                )
            )
        result = orientation_from_indexed_patterns(observations, WIDE)
        assert not result.is_consistent
        assert result.beam_deviation_deg > 5.0
        assert "NOT SELF-CONSISTENT" in result.describe()

    def test_an_inconsistent_fit_refuses_to_become_a_calibration(self) -> None:
        """Propagating a contradicted constant is worse than having none."""

        phase = nickel()
        truth = random_orientations(1, seed=2)[0]
        mirror = rotation_x(math.pi)
        observations = [
            IndexedPatternObservation(
                mirror @ pattern_rotation_for(truth, position_for(truth, zone), 25.0),
                ZoneAxis([-v for v in zone.indices], phase=phase),
                position_for(truth, zone),
            )
            for zone in (
                ZoneAxis([0, 0, 1], phase=phase),
                ZoneAxis([0, 1, 1], phase=phase),
                ZoneAxis([1, 1, 1], phase=phase),
            )
        ]
        result = orientation_from_indexed_patterns(observations, WIDE)
        with pytest.raises(ValueError, match="not self-consistent"):
            result.as_calibration()

    def test_a_consistent_fit_yields_a_usable_calibration(self) -> None:
        """Closing the loop: two patterns calibrate, then singles work directly."""

        phase = nickel()
        truth = random_orientations(1, seed=31)[0]
        observations = []
        for indices in ([0, 0, 1], [0, 1, 1]):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    pattern_rotation_for(truth, position, -63.0), zone, position
                )
            )
        fitted = orientation_from_indexed_patterns(observations, WIDE)
        calibration = fitted.as_calibration()
        assert calibration.is_rotation_calibrated
        assert calibration.diffraction_rotation_deg == pytest.approx(-63.0, abs=1e-9)

        # A third, previously unused pattern now converts directly.
        stage = DoubleTiltStage(
            envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
            calibration=calibration,
        )
        zone = ZoneAxis([1, 1, 1], phase=phase)
        position = position_for(truth, zone)
        single = orientation_from_indexed_pattern(
            pattern_rotation_for(truth, position, -63.0), zone, position, stage
        )
        assert symmetry_reduced_angle_deg(single.matrix, truth, phase) < 1e-9

    def test_zone_axis_senses_are_resolved_not_assumed(self) -> None:
        """A single pattern does not fix the sense; the fit chooses it."""

        phase = nickel()
        truth = random_orientations(1, seed=37)[0]
        observations = []
        for sign, indices in ((1, [0, 0, 1]), (-1, [0, 1, 1])):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    pattern_rotation_for(truth, position, 0.0),
                    ZoneAxis([sign * v for v in indices], phase=phase),
                    position,
                )
            )
        result = orientation_from_indexed_patterns(observations, WIDE)
        assert result.zone_axis_signs[1] == -1
        assert result.interzonal_residual_deg < 1e-9

    def test_requires_two_patterns(self) -> None:
        phase = nickel()
        zone = ZoneAxis([0, 0, 1], phase=phase)
        with pytest.raises(ValueError, match="at least two indexed patterns"):
            orientation_from_indexed_patterns(
                [IndexedPatternObservation(np.eye(3), zone, StagePosition(0.0, 0.0))],
                WIDE,
            )

    def test_rejects_parallel_zone_axes(self) -> None:
        phase = nickel()
        observations = [
            IndexedPatternObservation(
                np.eye(3), ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0)
            ),
            IndexedPatternObservation(
                np.eye(3), ZoneAxis([0, 0, 2], phase=phase), StagePosition(0.0, 0.0)
            ),
        ]
        with pytest.raises(ValueError, match="parallel"):
            orientation_from_indexed_patterns(observations, WIDE)

    def test_rejects_mixed_phases(self) -> None:
        first, second = nickel(), nickel()
        other = Phase(
            "other",
            lattice=Lattice(
                4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=second.crystal_frame
            ),
            symmetry=second.symmetry,
            crystal_frame=second.crystal_frame,
        )
        observations = [
            IndexedPatternObservation(
                np.eye(3), ZoneAxis([0, 0, 1], phase=first), StagePosition(0.0, 0.0)
            ),
            IndexedPatternObservation(
                np.eye(3),
                ZoneAxis([0, 1, 1], phase=other),
                StagePosition(10.0, 0.0),
            ),
        ]
        with pytest.raises(ValueError, match="one phase"):
            orientation_from_indexed_patterns(observations, WIDE)

    def test_describe_and_json_stay_in_lockstep(self) -> None:
        phase = nickel()
        truth = random_orientations(1, seed=41)[0]
        observations = []
        for indices in ([0, 0, 1], [0, 1, 1]):
            zone = ZoneAxis(indices, phase=phase)
            position = position_for(truth, zone)
            observations.append(
                IndexedPatternObservation(
                    pattern_rotation_for(truth, position, 5.0), zone, position
                )
            )
        result = orientation_from_indexed_patterns(observations, WIDE)
        text = result.describe()
        payload = result.to_json_dict()
        assert "determined" in text
        assert payload["pattern_count"] == 2
        assert payload["beam_deviation_deg"] == pytest.approx(result.beam_deviation_deg)
        assert payload["is_consistent"] is True


class TestReportBridge:
    def test_maps_over_ranked_solutions(self) -> None:
        """Competing indexings imply competing orientations; all are returned."""

        phase = nickel()
        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))

        class _Solution:
            def __init__(self, indices: list[int]) -> None:
                self.orientation = Rotation.identity()
                self.zone_axis = ZoneAxis(indices, phase=phase)

        class _Report:
            solutions = (_Solution([0, 0, 1]), _Solution([0, 1, 1]))

        results = orientations_from_pattern_report(
            _Report(), StagePosition(0.0, 0.0), stage
        )
        assert len(results) == 2
        assert [tuple(r.zone_axis.indices) for r in results] == [(0, 0, 1), (0, 1, 1)]

    def test_empty_report_yields_no_orientations(self) -> None:
        class _Empty:
            solutions = ()

        stage = DoubleTiltStage(calibration=StageCalibration(diffraction_rotation_deg=0.0))
        assert orientations_from_pattern_report(_Empty(), StagePosition(0.0, 0.0), stage) == ()


class TestPublicSurface:
    def test_exported_from_the_package_root(self) -> None:
        import pytex

        for name in (
            "orientation_from_indexed_pattern",
            "orientation_from_indexed_patterns",
            "IndexedPatternObservation",
            "IndexedOrientation",
            "MultiPatternOrientation",
        ):
            assert name in pytex.__all__, name
            assert hasattr(pytex, name), name
