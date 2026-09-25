"""The normative geometry and the frames of FIB lamella planning.

The lanes are section 12 of ``docs/architecture/fib_lamella_planning_foundation.md``:

- *analytic* --- orientations built so the out-of-plane angle is known in closed
  form: exactly in-plane (0), exactly along the normal (90), and a constructed
  12.5 degrees matched to 1e-9;
- *convention* --- a case whose answer changes if ``g`` is read as
  sample-to-crystal instead of crystal-to-sample, so the orientation convention
  of section 6 is pinned by a test that fails under the opposite reading;
- *invariance* --- applying any symmetry operator to ``g`` leaves ``eps*``,
  ``n_L`` and ``t_L`` unchanged;
- *round trip* --- the frame chain composed with its inverse is the identity;
- *registration* --- identity is a no-op, a known rotation and flip round-trip,
  and an affine fit reports its residual.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import (
    FIB_SAMPLE_FRAME,
    LAMELLA_FRAME,
    ChamberGeometry,
    FiducialObservation,
    ImageRegistration,
    SurfaceGeometry,
    calibrate_chamber_from_fiducials,
    lamella_frame_graph,
    lamella_geometry,
    lamella_geometry_batch,
    lamella_options,
    target_orbit,
    wrap_azimuth_deg,
)


def _frame() -> ReferenceFrame:
    return ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def _phase(name: str, point_group: str, a: float, c: float | None = None, gamma: float = 90.0):
    frame = _frame()
    lattice = Lattice(a, a, c if c is not None else a, 90.0, 90.0, gamma, crystal_frame=frame)
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


@pytest.fixture(scope="module")
def cubic() -> Phase:
    return _phase("nickel", "m-3m", 3.52)


@pytest.fixture(scope="module")
def tetragonal() -> Phase:
    return _phase("tetragonal", "4/mmm", 3.99, 4.04)


@pytest.fixture(scope="module")
def triclinic() -> Phase:
    return _phase("triclinic", "1", 4.0)


def _rx(degrees: float) -> np.ndarray:
    return Rotation.from_axis_angle([1.0, 0.0, 0.0], math.radians(degrees)).as_matrix()


def _rz(degrees: float) -> np.ndarray:
    return Rotation.from_axis_angle([0.0, 0.0, 1.0], math.radians(degrees)).as_matrix()


class TestOrbit:
    @pytest.mark.parametrize(
        ("target", "size"), [((1, 0, 0), 6), ((1, 1, 0), 12), ((1, 1, 1), 8), ((1, 1, 2), 24)]
    )
    def test_cubic_orbit_sizes_count_both_senses(self, cubic, target, size) -> None:
        assert len(target_orbit(cubic, target)) == size

    def test_one_sense_only_halves_a_non_centrosymmetric_orbit(self) -> None:
        polar = _phase("polar", "4mm", 3.9, 4.1)
        assert len(target_orbit(polar, (0, 0, 1))) == 2
        assert len(target_orbit(polar, (0, 0, 1), both_senses=False)) == 1

    def test_members_are_unit_images_of_their_indices(self, tetragonal) -> None:
        orbit = target_orbit(tetragonal, (1, 0, 1))
        direct = np.asarray(tetragonal.lattice.direct_basis().matrix)
        for indices, vector in zip(orbit.indices, orbit.cartesian, strict=True):
            expected = direct @ indices
            np.testing.assert_allclose(vector, expected / np.linalg.norm(expected), atol=1e-12)

    def test_the_target_is_written_as_a_family(self, cubic) -> None:
        assert target_orbit(cubic, (0, 1, 1)).family_text == "<011>"

    @pytest.mark.parametrize("bad", [(0, 0, 0), (1, 0), (0.5, 0, 1)])
    def test_bad_targets_are_refused(self, cubic, bad) -> None:
        with pytest.raises(ValueError):
            target_orbit(cubic, bad)


class TestAnalytic:
    def test_a_target_lying_in_the_surface_has_zero_residual(self, cubic) -> None:
        geometry = lamella_geometry(np.eye(3), target_orbit(cubic, (1, 0, 0)))
        assert geometry.eps_deg == pytest.approx(0.0, abs=1e-12)
        np.testing.assert_allclose(geometry.normal_sample, [1.0, 0.0, 0.0], atol=1e-12)
        np.testing.assert_allclose(geometry.long_axis_sample, [0.0, 1.0, 0.0], atol=1e-12)
        assert geometry.theta_sample_deg == pytest.approx(0.0, abs=1e-12)

    def test_a_target_along_the_normal_is_ninety_degrees_and_degenerate(self, tetragonal) -> None:
        geometry = lamella_geometry(np.eye(3), target_orbit(tetragonal, (0, 0, 1)))
        assert geometry.eps_deg == pytest.approx(90.0, abs=1e-12)
        assert geometry.degenerate

    def test_a_constructed_twelve_and_a_half_degrees_is_recovered_to_1e_9(self, tetragonal) -> None:
        # [001] rotated 77.5 degrees about x rises 12.5 degrees out of the surface.
        geometry = lamella_geometry(_rx(77.5), target_orbit(tetragonal, (0, 0, 1)))
        assert abs(geometry.eps_deg - 12.5) < 1e-9
        # Its in-plane part is along -y, so the canonical normal is +y (theta 90).
        assert geometry.theta_sample_deg == pytest.approx(90.0, abs=1e-9)
        assert not geometry.degenerate

    def test_the_lamella_frame_is_right_handed_with_the_normal_in_the_surface(self, cubic) -> None:
        g = Rotation.from_axis_angle([0.3, -0.5, 0.8], 1.1).as_matrix()
        geometry = lamella_geometry(g, target_orbit(cubic, (1, 1, 2)))
        n, t = geometry.normal_sample, geometry.long_axis_sample
        assert n[2] == pytest.approx(0.0, abs=1e-15)
        np.testing.assert_allclose(np.cross([0.0, 0.0, 1.0], n), t, atol=1e-12)
        assert np.linalg.det(geometry.sample_to_lamella_matrix()) == pytest.approx(1.0)
        # d* has no component along t_L and makes eps* with n_L.
        d = geometry.direction_sample
        assert float(d @ t) == pytest.approx(0.0, abs=1e-12)
        assert math.degrees(math.acos(float(d @ n))) == pytest.approx(geometry.eps_deg, abs=1e-9)


class TestOrientationConvention:
    """PyTex orientations map crystal to specimen; the geometry must use g, not g^T."""

    def test_the_answer_is_g_times_u_and_not_its_transpose(self, triclinic) -> None:
        # g = Rx(30) Rz(40): g [100] rises asin(sin40 sin30) = 18.75 deg out of the
        # surface, while the transposed reading g^T [100] = (cos40, -sin40, 0) lies
        # in it. A sign or transpose error in the geometry flips which one is found.
        g = _rx(30.0) @ _rz(40.0)
        orbit = target_orbit(triclinic, (1, 0, 0))
        expected = math.degrees(math.asin(math.sin(math.radians(40.0)) * 0.5))
        assert lamella_geometry(g, orbit).eps_deg == pytest.approx(expected, abs=1e-9)
        assert lamella_geometry(g.T, orbit).eps_deg == pytest.approx(0.0, abs=1e-9)

    def test_it_agrees_with_the_library_mapping_of_crystal_directions(self, cubic) -> None:
        euler = np.array([[35.0, 50.0, 12.0], [120.0, 20.0, 75.0]])
        orientations = OrientationSet.from_euler_angles(
            euler,
            crystal_frame=cubic.crystal_frame,
            specimen_frame=SurfaceGeometry().frame_transform().source,
            phase=cubic,
            symmetry=cubic.symmetry,
        )
        orbit = target_orbit(cubic, (1, 1, 1))
        batch = lamella_geometry_batch(orientations.as_matrices(), orbit)
        mapped = orientations.map_crystal_directions(orbit.cartesian[batch.member_index[0]])
        first = np.asarray(mapped)[0] if np.asarray(mapped).ndim == 2 else np.asarray(mapped)
        np.testing.assert_allclose(first, batch.direction_sample[0], atol=1e-12)


class TestInvariance:
    def test_a_symmetry_operator_on_g_changes_nothing(self, cubic) -> None:
        rng = np.random.default_rng(7)
        orbit = target_orbit(cubic, (1, 1, 2))
        operators = np.asarray(cubic.symmetry.operators)
        for _ in range(5):
            axis = rng.normal(size=3)
            g = Rotation.from_axis_angle(
                axis / np.linalg.norm(axis), rng.uniform(0, math.pi)
            ).as_matrix()
            reference = lamella_geometry(g, orbit)
            for operator in operators:
                other = lamella_geometry(g @ operator, orbit)
                assert other.eps_deg == pytest.approx(reference.eps_deg, abs=1e-9)
                assert other.theta_sample_deg == pytest.approx(reference.theta_sample_deg, abs=1e-7)
                np.testing.assert_allclose(other.normal_sample, reference.normal_sample, atol=1e-8)
                np.testing.assert_allclose(
                    other.long_axis_sample, reference.long_axis_sample, atol=1e-8
                )

    def test_the_batch_equals_the_single_solve(self, cubic) -> None:
        rng = np.random.default_rng(11)
        rotations = [
            Rotation.from_axis_angle(v / np.linalg.norm(v), a).as_matrix()
            for v, a in zip(rng.normal(size=(20, 3)), rng.uniform(0, math.pi, 20), strict=True)
        ]
        orbit = target_orbit(cubic, (0, 1, 1))
        batch = lamella_geometry_batch(np.stack(rotations), orbit)
        for index, g in enumerate(rotations):
            single = lamella_geometry(g, orbit)
            assert batch.eps_deg[index] == pytest.approx(single.eps_deg, abs=1e-12)
            assert batch.theta_sample_deg[index] == pytest.approx(single.theta_sample_deg)

    def test_the_options_start_with_the_chosen_lamella(self, cubic) -> None:
        g = Rotation.from_axis_angle([0.2, 0.4, 0.9], 0.7).as_matrix()
        orbit = target_orbit(cubic, (0, 1, 1))
        options = lamella_options(g, orbit)
        geometry = lamella_geometry(g, orbit)
        assert options[0].eps_deg == pytest.approx(geometry.eps_deg)
        assert options[0].theta_sample_deg == pytest.approx(geometry.theta_sample_deg)
        assert [item.eps_deg for item in options] == sorted(item.eps_deg for item in options)
        assert len({round(item.theta_sample_deg, 6) for item in options}) == len(options)


class TestFramesAndRoundTrip:
    @pytest.mark.parametrize("normal", [1, -1])
    @pytest.mark.parametrize("rows", [1, -1])
    def test_every_surface_convention_gives_a_proper_rotation(self, normal, rows) -> None:
        surface = SurfaceGeometry(normal_sign=normal, scan_y_sign=rows)
        matrix = surface.specimen_to_sample_matrix()
        assert np.linalg.det(matrix) == pytest.approx(1.0)
        np.testing.assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-15)
        point = np.array([[2.0, 3.0]])
        np.testing.assert_allclose(
            surface.sample_to_scan_xy(surface.scan_to_sample_xy(point)), point
        )

    def test_an_inward_normal_reverses_the_out_of_plane_sense(self, tetragonal) -> None:
        orbit = target_orbit(tetragonal, (0, 0, 1))
        up = lamella_geometry(_rx(77.5), orbit, SurfaceGeometry(normal_sign=1))
        down = lamella_geometry(_rx(77.5), orbit, SurfaceGeometry(normal_sign=-1))
        assert up.eps_deg == pytest.approx(down.eps_deg)
        # Y_s reverses with Z_s, so the in-plane azimuth is mirrored.
        assert down.theta_sample_deg == pytest.approx(
            float(wrap_azimuth_deg(-up.theta_sample_deg, period=180.0)), abs=1e-9
        )

    def test_the_rigid_chain_round_trips_to_1e_9(self, cubic) -> None:
        g = Rotation.from_axis_angle([0.6, -0.2, 0.3], 0.9).as_matrix()
        geometry = lamella_geometry(g, target_orbit(cubic, (1, 1, 1)))
        from pytex.fib import MountModel

        graph = lamella_frame_graph(
            SurfaceGeometry(normal_sign=-1),
            sample_to_lamella=geometry.sample_to_lamella_matrix(),
            lamella_to_holder=MountModel.lamella_to_holder(37.0, -1),
        )
        vectors = np.random.default_rng(3).normal(size=(10, 3))
        there = graph.convert(vectors, source="specimen", target="holder", directions=True)
        back = graph.convert(there, source="holder", target="specimen", directions=True)
        np.testing.assert_allclose(back, vectors, atol=1e-9)
        assert graph.has_frame(FIB_SAMPLE_FRAME.name) and graph.has_frame(LAMELLA_FRAME.name)

    def test_the_chosen_member_sits_at_eps_in_the_lamella_frame(self, cubic) -> None:
        g = Rotation.from_axis_angle([0.6, -0.2, 0.3], 0.9).as_matrix()
        geometry = lamella_geometry(g, target_orbit(cubic, (1, 1, 1)))
        local = geometry.sample_to_lamella_matrix() @ geometry.direction_sample
        eps = math.radians(geometry.eps_deg)
        np.testing.assert_allclose(
            local, [math.cos(eps), 0.0, geometry.out_of_plane_sign * math.sin(eps)], atol=1e-12
        )

    @pytest.mark.parametrize("flip", [1, -1])
    def test_mounting_rotations_are_proper_and_send_the_normal_to_the_beam(self, flip) -> None:
        from pytex.fib import MountModel

        matrices = MountModel.lamella_to_holder(np.arange(0.0, 360.0, 15.0), flip)
        np.testing.assert_allclose(np.linalg.det(matrices), 1.0)
        np.testing.assert_allclose(
            matrices[:, :, 0], [[0.0, 0.0, flip]] * len(matrices), atol=1e-15
        )


class TestRegistration:
    def test_identity_is_a_proven_no_op(self) -> None:
        reg = ImageRegistration.identity()
        points = np.random.default_rng(1).uniform(-50, 50, size=(25, 2))
        np.testing.assert_array_equal(reg.to_image(points), points)
        np.testing.assert_array_equal(reg.to_scan(points), points)
        for theta in (0.0, 17.5, 90.0, 133.0):
            assert reg.image_azimuth_deg(theta, SurfaceGeometry()) == pytest.approx(theta)
        assert reg.is_identity and reg.residual_deg == 0.0

    def test_a_known_rotation_and_flip_round_trip(self) -> None:
        reg = ImageRegistration.similarity(
            rotation_deg=30.0, scale=2.5, flip_y=True, translation=(4, -3)
        )
        points = np.random.default_rng(2).uniform(-10, 10, size=(12, 2))
        np.testing.assert_allclose(reg.to_scan(reg.to_image(points)), points, atol=1e-12)
        assert reg.is_mirror
        # A flip mirrors the azimuth, then the rotation adds 30 degrees.
        assert reg.image_azimuth_deg(10.0, SurfaceGeometry()) == pytest.approx(20.0)

    def test_an_affine_fit_recovers_the_map_and_reports_its_residual(self) -> None:
        truth = ImageRegistration(
            matrix=np.array([[1.9, 0.2], [-0.1, 2.1]]), translation=np.array([5.0, 7.0])
        )
        scan = np.array([[0.0, 0.0], [40.0, 0.0], [0.0, 30.0], [40.0, 30.0], [20.0, 15.0]])
        exact = ImageRegistration.from_control_points(scan, truth.to_image(scan))
        np.testing.assert_allclose(exact.matrix, truth.matrix, atol=1e-10)
        assert exact.residual_um == pytest.approx(0.0, abs=1e-9)
        noisy = truth.to_image(scan) + np.array(
            [[0.3, -0.2], [-0.1, 0.4], [0.2, 0.1], [-0.4, 0.0], [0.1, -0.3]]
        )
        fitted = ImageRegistration.from_control_points(scan, noisy)
        assert fitted.residual_um > 0.0 and fitted.residual_deg > 0.0
        assert "RMS residual" in fitted.describe()
        assert fitted.to_json_dict()["residual_deg"] == pytest.approx(fitted.residual_deg)

    def test_three_points_are_exact_and_say_so(self) -> None:
        scan = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]])
        reg = ImageRegistration.from_control_points(scan, scan * 2.0)
        assert reg.residual_um == 0.0
        assert "no redundancy" in reg.describe()

    def test_collinear_or_too_few_points_are_refused(self) -> None:
        with pytest.raises(ValueError, match="collinear"):
            ImageRegistration.from_control_points(
                [[0, 0], [1, 1], [2, 2]], [[0, 0], [1, 1], [2, 2]]
            )
        with pytest.raises(ValueError, match="three"):
            ImageRegistration.from_control_points([[0, 0], [1, 0]], [[0, 0], [1, 0]])


class TestChamber:
    def test_the_ion_azimuth_follows_the_documented_formula(self) -> None:
        chamber = ChamberGeometry(
            rotation_sense=-1, rotation_offset_deg=12.0, stage_rotation_deg=5.0
        )
        assert chamber.ion_azimuth_deg(40.0) == pytest.approx((-(40.0 + 5.0) + 12.0) % 360.0)
        assert not chamber.calibrated
        assert "UNCALIBRATED" in chamber.calibration_caveat()

    @pytest.mark.parametrize("angle", [0.0, 90.0, -5.0])
    def test_the_column_angle_must_be_between_zero_and_ninety(self, angle) -> None:
        with pytest.raises(ValueError):
            ChamberGeometry(column_angle_deg=angle)

    def test_the_default_is_fifty_four_and_fifty_two_is_standard(self) -> None:
        assert ChamberGeometry().column_angle_deg == 54.0
        assert ChamberGeometry(column_angle_deg=52.0).is_standard_column_angle
        assert not ChamberGeometry(column_angle_deg=45.0).is_standard_column_angle

    def test_foreshortening_is_cos_t_across_the_tilt_axis(self) -> None:
        chamber = ChamberGeometry(column_angle_deg=54.0)
        assert chamber.apparent_sem_length(10.0, 0.0) == pytest.approx(10.0)
        assert chamber.apparent_sem_length(10.0, 90.0) == pytest.approx(
            10.0 * math.cos(math.radians(54.0))
        )

    @pytest.mark.parametrize(
        ("sense", "offset", "stage"), [(1, 17.0, 0.0), (-1, 90.0, 0.0), (-1, 33.0, 20.0)]
    )
    def test_the_fiducial_calibration_recovers_sense_and_offset(self, sense, offset, stage) -> None:
        truth = ChamberGeometry(rotation_sense=sense, rotation_offset_deg=offset)
        observations = []
        for rho in (0.0, 30.0, 60.0, 100.0):
            # The trench normal azimuth that pattern rotation rho produces.
            theta = (sense * (rho - offset)) - stage
            observations.append(FiducialObservation(rho, theta % 180.0, stage_rotation_deg=stage))
        fitted = calibrate_chamber_from_fiducials(observations)
        assert fitted.calibrated and fitted.rotation_sense == sense
        assert fitted.calibration_residual_deg == pytest.approx(0.0, abs=1e-9)
        # The fitted offset reproduces every milled rotation modulo 180 degrees.
        for item in observations:
            predicted = (
                fitted.rotation_sense * (item.sample_azimuth_deg + stage)
                + fitted.rotation_offset_deg
            )
            assert (predicted - item.pattern_rotation_deg) % 180.0 == pytest.approx(
                0.0, abs=1e-6
            ) or (predicted - item.pattern_rotation_deg) % 180.0 == pytest.approx(180.0, abs=1e-6)
        assert truth.rotation_sense == fitted.rotation_sense

    def test_uninformative_fiducials_are_refused(self) -> None:
        with pytest.raises(ValueError, match="at least two"):
            calibrate_chamber_from_fiducials([FiducialObservation(0.0, 0.0)])
        with pytest.raises(ValueError, match="senses"):
            calibrate_chamber_from_fiducials(
                [FiducialObservation(0.0, 0.0), FiducialObservation(90.0, 90.0)]
            )
