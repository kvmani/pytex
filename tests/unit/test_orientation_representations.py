"""Tests for `pytex.core.representations`.

The equal-volume maps are tested against their *defining properties* rather than
against transcribed reference numbers: a volume-preserving bijection of the cube
onto the ball is pinned by its Jacobian determinant, by where the cube corner and
face centre land, and by round-tripping. Those checks would fail for any wrong
constant, which a table of copied outputs would not.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import crystal_frame, specimen_frame
from pytex.core.orientation import Orientation, Rotation
from pytex.core.representations import (
    CUBOCHORIC_CUBE_EDGE,
    CUBOCHORIC_CUBE_HALF_EDGE,
    HOMOCHORIC_BALL_RADIUS,
    ORIENTATION_REPRESENTATIONS_SCHEMA,
    OrientationRepresentationSet,
    RepresentationKind,
    convert_orientations,
    cubochoric_from_homochoric,
    cubochoric_from_quaternions,
    homochoric_from_cubochoric,
    homochoric_from_quaternions,
    ideal_orientation_indices,
    orientation_representations,
    quaternions_from_cubochoric,
    quaternions_from_euler_angles,
    quaternions_from_homochoric,
    quaternions_to_euler_angles,
    rotation_representations,
)


@pytest.fixture(scope="module")
def ni_fcc_phase() -> object:
    """Nickel, the canonical FCC case, from the pinned CIF fixture."""

    pytest.importorskip(
        "pymatgen.core",
        reason="CIF-backed phase fixtures require the optional pymatgen dependency.",
    )
    from pytex.core.fixtures import get_phase_fixture

    return get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal_frame())


@pytest.fixture(scope="module")
def zr_hcp_phase() -> object:
    """Zirconium, the canonical HCP case, from the pinned CIF fixture."""

    pytest.importorskip(
        "pymatgen.core",
        reason="CIF-backed phase fixtures require the optional pymatgen dependency.",
    )
    from pytex.core.fixtures import get_phase_fixture

    return get_phase_fixture("zr_hcp").load_phase(crystal_frame=crystal_frame())


def _random_quaternions(count: int, seed: int = 20260809) -> np.ndarray:
    rng = np.random.default_rng(seed)
    quaternions = rng.normal(size=(count, 4))
    quaternions /= np.linalg.norm(quaternions, axis=1)[:, None]
    return np.where(quaternions[:, :1] < 0.0, -quaternions, quaternions)


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #


def test_ball_and_cube_enclose_the_same_volume() -> None:
    """The equal-volume map exists only because these two volumes agree."""

    ball = (4.0 / 3.0) * np.pi * HOMOCHORIC_BALL_RADIUS**3
    cube = CUBOCHORIC_CUBE_EDGE**3
    assert ball == pytest.approx(np.pi**2, rel=1e-14)
    assert cube == pytest.approx(np.pi**2, rel=1e-14)
    assert CUBOCHORIC_CUBE_HALF_EDGE == pytest.approx(CUBOCHORIC_CUBE_EDGE / 2.0)


# --------------------------------------------------------------------------- #
# Homochoric
# --------------------------------------------------------------------------- #


def test_homochoric_radius_follows_the_analytic_law() -> None:
    """``|h| = [3/4 (omega - sin omega)]^(1/3)``, checked at three known angles."""

    angles_rad = np.array([np.pi / 3.0, np.pi / 2.0, np.pi])
    axis = np.tile(np.array([0.0, 0.0, 1.0]), (3, 1))
    quaternions = np.column_stack(
        [np.cos(angles_rad / 2.0), axis * np.sin(angles_rad / 2.0)[:, None]]
    )
    radii = np.linalg.norm(homochoric_from_quaternions(quaternions), axis=1)
    expected = np.cbrt(0.75 * (angles_rad - np.sin(angles_rad)))
    assert radii == pytest.approx(expected, abs=1e-12)
    # A rotation by pi sits exactly on the bounding sphere.
    assert radii[-1] == pytest.approx(HOMOCHORIC_BALL_RADIUS, abs=1e-12)


def test_identity_maps_to_the_origin_in_both_equal_volume_charts() -> None:
    identity = np.array([[1.0, 0.0, 0.0, 0.0]])
    assert homochoric_from_quaternions(identity) == pytest.approx(np.zeros((1, 3)))
    assert cubochoric_from_quaternions(identity) == pytest.approx(np.zeros((1, 3)))


def test_homochoric_round_trips_through_quaternions() -> None:
    quaternions = _random_quaternions(512)
    recovered = quaternions_from_homochoric(homochoric_from_quaternions(quaternions))
    assert recovered == pytest.approx(quaternions, abs=1e-12)


def test_homochoric_rejects_a_vector_outside_the_ball() -> None:
    """Rodrigues vectors share the shape and are unbounded; the mix-up must raise."""

    with pytest.raises(ValueError, match="outside the homochoric ball"):
        quaternions_from_homochoric([[2.0, 0.0, 0.0]])


# --------------------------------------------------------------------------- #
# Cubochoric: the defining properties of the equal-volume map
# --------------------------------------------------------------------------- #


def test_cube_corner_lands_on_the_ball_surface() -> None:
    """A corner of the cube is a rotation by pi, so it must land on the sphere."""

    corner = np.full((1, 3), CUBOCHORIC_CUBE_HALF_EDGE)
    assert np.linalg.norm(homochoric_from_cubochoric(corner)) == pytest.approx(
        HOMOCHORIC_BALL_RADIUS, abs=1e-12
    )


def test_face_centre_maps_along_its_own_axis_to_the_ball_surface() -> None:
    """The nested-surface law: the face centre of the full cube is the pole."""

    face = np.array([[0.0, 0.0, CUBOCHORIC_CUBE_HALF_EDGE]])
    assert homochoric_from_cubochoric(face)[0] == pytest.approx(
        np.array([0.0, 0.0, HOMOCHORIC_BALL_RADIUS]), abs=1e-12
    )


def test_nested_cubes_map_onto_nested_spheres() -> None:
    """Every point of the sub-cube of half-edge ``t`` lands on one sphere.

    This is the first of the two conditions that determine the map, and it is
    what makes the radial law ``r = t (6 / pi)^(1/3)`` rather than anything else.
    """

    rng = np.random.default_rng(11)
    half = 0.63 * CUBOCHORIC_CUBE_HALF_EDGE
    faces = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            points = rng.uniform(-half, half, size=(200, 3))
            points[:, axis] = sign * half
            faces.append(points)
    radii = np.linalg.norm(homochoric_from_cubochoric(np.vstack(faces)), axis=1)
    expected = half * (6.0 / np.pi) ** (1.0 / 3.0)
    assert radii == pytest.approx(np.full(radii.shape, expected), abs=1e-12)


def test_cube_to_ball_map_preserves_volume_pointwise() -> None:
    """The Jacobian determinant is 1 everywhere, which is the whole point.

    Computed by central differences at random interior points, plus points just
    off the face diagonal and just off the pyramid boundary. The probes stay
    *off* those surfaces by more than the difference step, because the map has a
    kink there: a centred difference straddling a kink averages two different
    derivatives and would report a determinant belonging to neither side.
    Continuity across those same surfaces is tested separately.
    """

    rng = np.random.default_rng(5)
    step = 1e-6
    probes = np.vstack(
        [
            rng.uniform(-0.95, 0.95, size=(12, 3)) * CUBOCHORIC_CUBE_HALF_EDGE,
            np.array([[0.41, 0.39, 0.90], [0.52, 0.47, 0.55], [0.70, 0.10, 0.73]])
            * CUBOCHORIC_CUBE_HALF_EDGE,
        ]
    )
    for probe in probes:
        jacobian = np.empty((3, 3))
        for column in range(3):
            offset = np.zeros(3)
            offset[column] = step
            forward = homochoric_from_cubochoric((probe + offset)[None, :])[0]
            backward = homochoric_from_cubochoric((probe - offset)[None, :])[0]
            jacobian[:, column] = (forward - backward) / (2.0 * step)
        assert abs(float(np.linalg.det(jacobian))) == pytest.approx(1.0, abs=1e-6)


def test_cube_to_ball_map_is_continuous_across_the_branch_boundaries() -> None:
    """Straddle the pyramid boundary and the face diagonal, where branches swap."""

    half = CUBOCHORIC_CUBE_HALF_EDGE
    gap = 1e-8
    pairs = [
        # |x| = |z|: the pyramid boundary.
        (np.array([0.6 * half - gap, 0.2 * half, 0.6 * half]),
         np.array([0.6 * half + gap, 0.2 * half, 0.6 * half])),
        # |x| = |y|: the face diagonal, where the wedge map swaps its roles.
        (np.array([0.5 * half, 0.5 * half - gap, 0.9 * half]),
         np.array([0.5 * half, 0.5 * half + gap, 0.9 * half])),
        # z = 0 with x dominant: the map through the cube centre.
        (np.array([0.4 * half, 0.1 * half, -gap]),
         np.array([0.4 * half, 0.1 * half, gap])),
    ]
    for left, right in pairs:
        images = homochoric_from_cubochoric(np.vstack([left, right]))
        assert np.max(np.abs(images[0] - images[1])) < 1e-7


def test_cubochoric_round_trips_both_ways() -> None:
    rng = np.random.default_rng(2024)
    cube = rng.uniform(-CUBOCHORIC_CUBE_HALF_EDGE, CUBOCHORIC_CUBE_HALF_EDGE, (2048, 3))
    assert cubochoric_from_homochoric(homochoric_from_cubochoric(cube)) == pytest.approx(
        cube, abs=1e-11
    )
    quaternions = _random_quaternions(512, seed=31)
    assert quaternions_from_cubochoric(
        cubochoric_from_quaternions(quaternions)
    ) == pytest.approx(quaternions, abs=1e-11)


def test_uniform_cubochoric_sampling_reproduces_the_so3_angle_law() -> None:
    """Uniform in the cube must mean uniform in SO(3).

    The invariant measure gives rotation angles the density
    ``(1 - cos w) / pi`` on ``[0, pi]``, whose mean is
    ``pi / 2 + 2 / pi`` radians = 126.4756 deg. That number is analytic, not a
    prior program output, and it is the sharpest available test that the map is
    the equal-volume one.
    """

    rng = np.random.default_rng(99)
    cube = rng.uniform(
        -CUBOCHORIC_CUBE_HALF_EDGE, CUBOCHORIC_CUBE_HALF_EDGE, (200_000, 3)
    )
    quaternions = quaternions_from_cubochoric(cube)
    angles_deg = np.degrees(2.0 * np.arccos(np.clip(quaternions[:, 0], -1.0, 1.0)))
    analytic_mean_deg = np.degrees(np.pi / 2.0 + 2.0 / np.pi)
    assert float(angles_deg.mean()) == pytest.approx(analytic_mean_deg, abs=0.3)


def test_cubochoric_rejects_a_coordinate_outside_the_cube() -> None:
    with pytest.raises(ValueError, match="outside the cubochoric cube"):
        homochoric_from_cubochoric([[CUBOCHORIC_CUBE_HALF_EDGE * 1.5, 0.0, 0.0]])


# --------------------------------------------------------------------------- #
# Vectorized Euler conversions
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("convention", ["bunge", "matthies"])
def test_vectorized_euler_matches_the_object_path_exactly(convention: str) -> None:
    """The array path must not be a second, subtly different implementation."""

    rng = np.random.default_rng(4)
    angles = np.column_stack(
        [
            rng.uniform(0.0, 360.0, 64),
            rng.uniform(0.0, 180.0, 64),
            rng.uniform(0.0, 360.0, 64),
        ]
    )
    vectorized = quaternions_from_euler_angles(angles, convention=convention)
    per_object = np.stack(
        [
            Rotation.from_euler(*triple, convention=convention).quaternion
            for triple in angles
        ]
    )
    assert vectorized == pytest.approx(per_object, abs=1e-14)

    back = quaternions_to_euler_angles(vectorized, convention=convention)
    per_object_back = np.stack(
        [
            Rotation(quaternion=quaternion).to_euler(convention=convention)
            for quaternion in vectorized
        ]
    )
    assert back == pytest.approx(per_object_back, abs=1e-9)


def test_vectorized_euler_handles_the_gimbal_degenerate_row() -> None:
    """At ``Phi = 0`` only ``phi1 + phi2`` is determined; the rotation must survive."""

    angles = np.array([[30.0, 0.0, 20.0]])
    quaternions = quaternions_from_euler_angles(angles)
    recovered = quaternions_to_euler_angles(quaternions)
    assert recovered[0, 1] == pytest.approx(0.0, abs=1e-9)
    assert recovered[0, 2] == pytest.approx(0.0, abs=1e-9)
    assert quaternions_from_euler_angles(recovered) == pytest.approx(
        quaternions, abs=1e-12
    )


# --------------------------------------------------------------------------- #
# The generic converter
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", list(RepresentationKind))
def test_every_representation_round_trips_through_the_converter(
    kind: RepresentationKind,
) -> None:
    quaternions = _random_quaternions(128, seed=77)
    values = convert_orientations(
        quaternions, source=RepresentationKind.QUATERNION, target=kind
    )
    recovered = convert_orientations(
        values, source=kind, target=RepresentationKind.QUATERNION
    )
    # q and -q denote the same rotation, and the matrix and Euler routes do not
    # promise to return the same representative, so compare the rotations.
    alignment = np.abs(np.sum(recovered * quaternions, axis=1))
    assert alignment == pytest.approx(np.ones(quaternions.shape[0]), abs=1e-9)


def test_converter_accepts_plain_strings_for_the_kinds() -> None:
    cube = convert_orientations(
        [[0.0, 0.0, 0.0]], source="euler_bunge", target="cubochoric"
    )
    assert cube == pytest.approx(np.zeros((1, 3)), abs=1e-12)


def test_converter_rejects_an_unknown_kind() -> None:
    with pytest.raises(ValueError):
        convert_orientations([[1.0, 0.0, 0.0, 0.0]], source="quaternion", target="wxyz")


# --------------------------------------------------------------------------- #
# Ideal (hkl)[uvw] recovery
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("plane", "direction"),
    [((1, 1, 2), (1, 1, -1)), ((1, 1, 0), (1, -1, 2)), ((0, 0, 1), (1, 0, 0))],
)
def test_ideal_indices_invert_from_miller_exactly(
    ni_fcc_phase: object, plane: tuple[int, int, int], direction: tuple[int, int, int]
) -> None:
    """Copper, brass, and cube: build from the indices and read them back."""

    orientation = Orientation.from_miller(
        plane, direction, phase=ni_fcc_phase, specimen_frame=specimen_frame()
    )
    indices = ideal_orientation_indices(orientation)
    assert tuple(int(value) for value in indices.hkl) == plane
    assert tuple(int(value) for value in indices.uvw) == direction
    assert indices.is_exact
    assert indices.plane_deviation_deg == pytest.approx(0.0, abs=1e-9)
    assert indices.direction_deviation_deg == pytest.approx(0.0, abs=1e-9)


def test_ideal_indices_report_the_deviation_of_a_perturbed_orientation(
    ni_fcc_phase: object,
) -> None:
    """A label without its deviation would be a claim the data does not support."""

    phase = ni_fcc_phase
    exact = Orientation.from_miller(
        (1, 1, 0), (1, -1, 2), phase=phase, specimen_frame=specimen_frame()
    )
    tilt = Rotation.from_axis_angle((0.0, 1.0, 0.0), np.radians(4.0))
    perturbed = Orientation.from_matrix(
        tilt.as_matrix() @ exact.as_matrix(),
        specimen_frame=specimen_frame(),
        phase=phase,
        crystal_frame=phase.crystal_frame,
    )
    indices = ideal_orientation_indices(perturbed, max_index=4)
    assert not indices.is_exact
    assert indices.plane_deviation_deg > 0.5
    assert "nearest label" in indices.describe()


def test_hexagonal_phase_gets_four_index_labels(zr_hcp_phase: object) -> None:
    """The literature writes hexagonal indices with four components; so must we."""

    orientation = Orientation.from_miller(
        (0, 0, 1), (1, 0, 0), phase=zr_hcp_phase, specimen_frame=specimen_frame()
    )
    indices = ideal_orientation_indices(orientation)
    assert indices.hkil is not None
    assert indices.uvtw is not None
    assert tuple(int(value) for value in indices.hkil) == (0, 0, 0, 1)
    # The zone-law identity h + k + i = 0 must hold for the plane indices.
    assert int(indices.hkil[0] + indices.hkil[1] + indices.hkil[2]) == 0
    # And u + v + t = 0 for the direction indices.
    assert int(indices.uvtw[0] + indices.uvtw[1] + indices.uvtw[2]) == 0


def test_ideal_indices_require_a_phase() -> None:
    orientation = Orientation.from_matrix(
        np.eye(3), specimen_frame=specimen_frame(), crystal_frame=crystal_frame()
    )
    with pytest.raises(ValueError, match="needs the phase"):
        ideal_orientation_indices(orientation)


# --------------------------------------------------------------------------- #
# The reports
# --------------------------------------------------------------------------- #


def test_all_representations_denote_the_same_rotation(ni_fcc_phase: object) -> None:
    """Every field must reconstruct the one rotation the report is about."""

    orientation = Orientation.from_miller(
        (1, 1, 2), (1, 1, -1), phase=ni_fcc_phase, specimen_frame=specimen_frame()
    )
    report = orientation_representations(orientation)
    reference = orientation.as_matrix()

    assert report.matrix == pytest.approx(reference, abs=1e-12)
    for kind, values in (
        (RepresentationKind.QUATERNION, report.quaternion[None, :]),
        (
            RepresentationKind.AXIS_ANGLE,
            np.concatenate([report.axis, [np.radians(report.angle_deg)]])[None, :],
        ),
        (RepresentationKind.RODRIGUES, report.rodrigues[None, :]),
        (RepresentationKind.RODRIGUES_FRANK, report.rodrigues_frank[None, :]),
        (RepresentationKind.EULER_BUNGE, report.euler_bunge_deg[None, :]),
        (RepresentationKind.EULER_MATTHIES, report.euler_matthies_deg[None, :]),
        (RepresentationKind.HOMOCHORIC, report.homochoric[None, :]),
        (RepresentationKind.CUBOCHORIC, report.cubochoric[None, :]),
    ):
        rebuilt = convert_orientations(
            values, source=kind, target=RepresentationKind.MATRIX
        )[0]
        assert rebuilt == pytest.approx(reference, abs=1e-9), kind


def test_report_describe_and_json_stay_in_lockstep(ni_fcc_phase: object) -> None:
    orientation = Orientation.from_miller(
        (1, 1, 0), (1, -1, 2), phase=ni_fcc_phase, specimen_frame=specimen_frame()
    )
    report = orientation_representations(orientation)
    payload = report.to_json_dict()
    prose = report.describe()

    assert payload["schema"] == ORIENTATION_REPRESENTATIONS_SCHEMA
    assert payload["phase"] == ni_fcc_phase.name  # type: ignore[attr-defined]
    assert payload["euler_bunge_deg"] == pytest.approx(
        report.euler_bunge_deg.tolist(), abs=1e-12
    )
    assert payload["ideal_orientation"] is not None
    assert payload["ideal_orientation"]["label"] == report.ideal_indices.label  # type: ignore[union-attr]
    assert "crystal-to-specimen" in prose
    assert "Bunge" in prose and "ZYZ" in prose
    assert report.ideal_indices is not None
    assert report.ideal_indices.label in prose
    table = report.to_table()
    assert "cubochoric" in table and "Rodrigues-Frank" in table


def test_bare_rotation_report_carries_no_crystal_claims() -> None:
    """A stage tilt has no phase, so it must not be given an (hkl)[uvw] name."""

    report = rotation_representations(
        Rotation.from_axis_angle((0.0, 0.0, 1.0), np.radians(30.0))
    )
    assert report.ideal_indices is None
    assert report.phase_name is None
    assert report.angle_deg == pytest.approx(30.0, abs=1e-12)
    assert "not an orientation" in report.describe()
    assert report.to_json_dict()["ideal_orientation"] is None


def test_representation_set_matches_the_single_orientation_reports() -> None:
    quaternions = _random_quaternions(32, seed=12)
    batch = OrientationRepresentationSet.from_quaternions(quaternions)
    assert len(batch) == 32
    for index in (0, 7, 31):
        single = rotation_representations(Rotation(quaternion=quaternions[index]))
        row = batch.row(index)
        assert row.cubochoric == pytest.approx(single.cubochoric, abs=1e-12)
        assert row.euler_bunge_deg == pytest.approx(single.euler_bunge_deg, abs=1e-12)
        assert row.angle_deg == pytest.approx(single.angle_deg, abs=1e-12)
    assert "ten representations" in batch.describe()


def test_reports_use_the_canonical_quaternion_sign() -> None:
    """``q`` and ``-q`` are one rotation; a *report* must pick one of them.

    Without this, a batch row and the single-orientation report for the same
    rotation could differ by a global sign, and a componentwise comparison of
    two identical rotations could report a difference of 2.
    """

    quaternions = _random_quaternions(64, seed=5)
    flipped = -quaternions
    batch = OrientationRepresentationSet.from_quaternions(flipped)
    assert np.all(batch.quaternions[:, 0] >= 0.0)
    assert batch.quaternions == pytest.approx(quaternions, abs=1e-12)

    single = rotation_representations(Rotation(quaternion=flipped[3]))
    assert single.quaternion == pytest.approx(batch.quaternions[3], abs=1e-12)


def test_canonical_sign_resolves_the_180_degree_tie() -> None:
    """At 180 degrees both signs have ``w = 0``, so the tie needs its own rule."""

    from pytex.core.representations import canonical_quaternions

    half_turns = np.array([[0.0, -1.0, 0.0, 0.0], [0.0, 0.0, 0.0, -1.0]])
    canonical = canonical_quaternions(half_turns)
    assert canonical == pytest.approx(np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]))
    # Idempotent: canonicalizing an already-canonical batch changes nothing.
    assert canonical_quaternions(canonical) == pytest.approx(canonical)


def test_representation_set_can_be_built_from_any_single_form() -> None:
    angles = np.array([[35.0, 45.0, 0.0], [59.0, 37.0, 63.0]])
    batch = OrientationRepresentationSet.from_values(
        angles, source=RepresentationKind.EULER_BUNGE
    )
    assert batch.euler_bunge_deg == pytest.approx(angles, abs=1e-9)
    assert batch.matrices.shape == (2, 3, 3)


def test_empty_representation_set_describes_itself_without_crashing() -> None:
    batch = OrientationRepresentationSet.from_quaternions(np.zeros((0, 4)))
    assert len(batch) == 0
    assert "empty" in batch.describe()
