from __future__ import annotations

from dataclasses import dataclass, field
from typing import TypedDict

import numpy as np

from pytex.core._arrays import (
    as_float_array,
    as_int_array,
    is_rotation_matrix,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, MillerIndex, Phase, ReciprocalLatticeVector, ZoneAxis
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord

_DETECTOR_PROJECTION_EPSILON = 1e-12
_BRAGG_ARGUMENT_TOLERANCE = 1e-12
_ZONE_AXIS_ORTHOGONALITY_ATOL = 1e-8
_ZONE_AXIS_ORTHOGONALITY_RTOL = 1e-8
_INTENSITY_EPSILON = 1e-12


class _CandidateSpot(TypedDict):
    miller_indices: np.ndarray
    reciprocal_vector_lab: np.ndarray
    outgoing_direction_lab: np.ndarray
    detector_coordinates_px: np.ndarray
    excitation_error_inv_angstrom: float
    intensity: float
    two_theta_rad: float
    azimuth_rad: float
    on_detector: bool
    accepted_by_mask: bool
    family_key: tuple[float, ...]


def _rotation_matrix_x(angle_rad: float) -> np.ndarray:
    cos_angle = float(np.cos(angle_rad))
    sin_angle = float(np.sin(angle_rad))
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, cos_angle, -sin_angle],
            [0.0, sin_angle, cos_angle],
        ],
        dtype=np.float64,
    )


def _rotation_matrix_y(angle_rad: float) -> np.ndarray:
    cos_angle = float(np.cos(angle_rad))
    sin_angle = float(np.sin(angle_rad))
    return np.array(
        [
            [cos_angle, 0.0, sin_angle],
            [0.0, 1.0, 0.0],
            [-sin_angle, 0.0, cos_angle],
        ],
        dtype=np.float64,
    )


def _rotation_matrix_z(angle_rad: float) -> np.ndarray:
    cos_angle = float(np.cos(angle_rad))
    sin_angle = float(np.sin(angle_rad))
    return np.array(
        [
            [cos_angle, -sin_angle, 0.0],
            [sin_angle, cos_angle, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )


def _kinematic_intensity(
    reciprocal_vector_lab: np.ndarray,
    excitation_error_inv_angstrom: float,
    *,
    model: str,
    excitation_sigma_inv_angstrom: float,
) -> float:
    if model == "unit":
        return 1.0
    if model != "kinematic_proxy":
        raise ValueError("intensity_model must be either 'unit' or 'kinematic_proxy'.")
    reciprocal_magnitude = float(np.linalg.norm(reciprocal_vector_lab))
    excitation_ratio = excitation_error_inv_angstrom / excitation_sigma_inv_angstrom
    excitation_weight = 1.0 / (1.0 + excitation_ratio * excitation_ratio)
    resolution_weight = 1.0 / (1.0 + reciprocal_magnitude * reciprocal_magnitude)
    return float(excitation_weight * resolution_weight)


def _reflection_family_key(miller_indices: np.ndarray, phase: Phase) -> tuple[float, ...]:
    reciprocal_vector = ReciprocalLatticeVector.from_miller_index(
        MillerIndex(miller_indices, phase=phase)
    ).cartesian_vector
    canonical_direction = phase.symmetry.canonicalize_vector(reciprocal_vector, antipodal=False)
    magnitude = float(np.linalg.norm(reciprocal_vector))
    rounded_direction = tuple(float(value) for value in np.round(canonical_direction, decimals=8))
    return (*rounded_direction, float(np.round(magnitude, decimals=8)))


@dataclass(frozen=True, slots=True)
class DiffractionGeometry:
    detector_frame: ReferenceFrame
    specimen_frame: ReferenceFrame
    laboratory_frame: ReferenceFrame
    beam_energy_kev: float
    camera_length_mm: float
    pattern_center: np.ndarray
    detector_pixel_size_um: tuple[float, float]
    detector_shape: tuple[int, int]
    beam_direction_lab: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 0.0, 1.0], dtype=np.float64)
    )
    specimen_to_lab_matrix: np.ndarray = field(default_factory=lambda: np.eye(3, dtype=np.float64))
    tilt_degrees: tuple[float, float, float] = (0.0, 0.0, 0.0)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.detector_frame.domain is not FrameDomain.DETECTOR:
            raise ValueError(
                "DiffractionGeometry.detector_frame must belong to the detector domain."
            )
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError(
                "DiffractionGeometry.specimen_frame must belong to the specimen domain."
            )
        if self.laboratory_frame.domain is not FrameDomain.LABORATORY:
            raise ValueError(
                "DiffractionGeometry.laboratory_frame must belong to the laboratory domain."
            )
        if self.beam_energy_kev <= 0.0:
            raise ValueError("DiffractionGeometry.beam_energy_kev must be strictly positive.")
        if self.camera_length_mm <= 0.0:
            raise ValueError("DiffractionGeometry.camera_length_mm must be strictly positive.")
        if len(self.detector_pixel_size_um) != 2:
            raise ValueError("DiffractionGeometry.detector_pixel_size_um must have length 2.")
        if len(self.detector_shape) != 2:
            raise ValueError("DiffractionGeometry.detector_shape must have length 2.")
        if any(value <= 0.0 for value in self.detector_pixel_size_um):
            raise ValueError(
                "DiffractionGeometry.detector_pixel_size_um values must be strictly positive."
            )
        if any(value <= 0 for value in self.detector_shape):
            raise ValueError("DiffractionGeometry.detector_shape values must be strictly positive.")
        if len(self.tilt_degrees) != 3:
            raise ValueError("DiffractionGeometry.tilt_degrees must contain exactly three values.")
        object.__setattr__(self, "pattern_center", as_float_array(self.pattern_center, shape=(3,)))
        if np.any(~np.isfinite(self.pattern_center)):
            raise ValueError("DiffractionGeometry.pattern_center must be finite.")
        if np.any((self.pattern_center[:2] < 0.0) | (self.pattern_center[:2] > 1.0)):
            raise ValueError(
                "DiffractionGeometry.pattern_center x and y components must lie in [0, 1]."
            )
        if self.pattern_center[2] <= 0.0:
            raise ValueError("DiffractionGeometry.pattern_center z component must be positive.")
        object.__setattr__(self, "beam_direction_lab", normalize_vector(self.beam_direction_lab))
        object.__setattr__(
            self,
            "specimen_to_lab_matrix",
            as_float_array(self.specimen_to_lab_matrix, shape=(3, 3)),
        )
        if not is_rotation_matrix(self.specimen_to_lab_matrix):
            raise ValueError(
                "DiffractionGeometry.specimen_to_lab_matrix must be a rotation matrix."
            )
        object.__setattr__(
            self,
            "detector_pixel_size_um",
            tuple(float(value) for value in self.detector_pixel_size_um),
        )
        object.__setattr__(
            self,
            "detector_shape",
            tuple(int(value) for value in self.detector_shape),
        )
        object.__setattr__(
            self,
            "tilt_degrees",
            tuple(float(value) for value in self.tilt_degrees),
        )

    @property
    def electron_wavelength_angstrom(self) -> float:
        voltage = self.beam_energy_kev * 1000.0
        numerator = 12.2639
        denominator = np.sqrt(voltage * (1.0 + 0.97845e-6 * voltage))
        return float(numerator / denominator)

    @property
    def ewald_sphere_radius_inv_angstrom(self) -> float:
        return float(1.0 / self.electron_wavelength_angstrom)

    @property
    def incident_wavevector_lab(self) -> np.ndarray:
        wavevector = self.beam_direction_lab / self.electron_wavelength_angstrom
        wavevector = np.ascontiguousarray(wavevector)
        wavevector.setflags(write=False)
        return wavevector

    @property
    def pattern_center_px(self) -> np.ndarray:
        width_px = self.detector_shape[1]
        height_px = self.detector_shape[0]
        center = np.array(
            [
                self.pattern_center[0] * (width_px - 1),
                self.pattern_center[1] * (height_px - 1),
            ],
            dtype=np.float64,
        )
        center.setflags(write=False)
        return center

    @property
    def detector_basis_lab(self) -> np.ndarray:
        beam = self.beam_direction_lab
        trial_up = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        if np.isclose(abs(float(np.dot(beam, trial_up))), 1.0, atol=1e-8):
            trial_up = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        u_axis = normalize_vector(np.cross(trial_up, beam))
        v_axis = normalize_vector(np.cross(beam, u_axis))
        local_basis = np.column_stack([u_axis, v_axis, beam])
        tilt_rad = np.deg2rad(np.asarray(self.tilt_degrees, dtype=np.float64))
        local_rotation = (
            _rotation_matrix_x(float(tilt_rad[0]))
            @ _rotation_matrix_y(float(tilt_rad[1]))
            @ _rotation_matrix_z(float(tilt_rad[2]))
        )
        basis = local_basis @ local_rotation
        basis = np.ascontiguousarray(basis)
        basis.setflags(write=False)
        return basis

    def specimen_vectors_to_lab(self, vectors: np.ndarray) -> np.ndarray:
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Specimen vectors must end with dimension 3.")
        mapped = array @ self.specimen_to_lab_matrix.T
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def lab_vectors_to_specimen(self, vectors: np.ndarray) -> np.ndarray:
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Laboratory vectors must end with dimension 3.")
        mapped = array @ self.specimen_to_lab_matrix
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def detector_coordinates_mm(self, coordinates_px: np.ndarray) -> np.ndarray:
        detector_pixels = as_float_array(coordinates_px, shape=(None, 2))
        pixel_size_mm = np.array(self.detector_pixel_size_um, dtype=np.float64) / 1000.0
        offsets_mm = (detector_pixels - self.pattern_center_px[None, :]) * pixel_size_mm[None, :]
        offsets_mm = np.ascontiguousarray(offsets_mm)
        offsets_mm.setflags(write=False)
        return offsets_mm

    def detector_points_lab_mm(self, coordinates_px: np.ndarray) -> np.ndarray:
        offsets_mm = self.detector_coordinates_mm(coordinates_px)
        basis = self.detector_basis_lab
        center_lab = basis[:, 2] * self.camera_length_mm
        points_lab = (
            center_lab[None, :]
            + offsets_mm[:, [0]] * basis[:, 0][None, :]
            + offsets_mm[:, [1]] * basis[:, 1][None, :]
        )
        points_lab = np.ascontiguousarray(points_lab)
        points_lab.setflags(write=False)
        return points_lab

    def outgoing_directions_lab(self, coordinates_px: np.ndarray) -> np.ndarray:
        directions = normalize_vectors(self.detector_points_lab_mm(coordinates_px))
        directions = np.ascontiguousarray(directions)
        directions.setflags(write=False)
        return directions

    def project_directions_to_detector_px(
        self, directions_lab: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        directions = normalize_vectors(directions_lab)
        basis = self.detector_basis_lab
        detector_normal = basis[:, 2]
        denominators = directions @ detector_normal
        valid = denominators > _DETECTOR_PROJECTION_EPSILON
        scales = np.full(directions.shape[0], np.nan, dtype=np.float64)
        scales[valid] = self.camera_length_mm / denominators[valid]
        intersection = directions * scales[:, None]
        u_offsets = intersection @ basis[:, 0]
        v_offsets = intersection @ basis[:, 1]
        pixel_size_mm = np.array(self.detector_pixel_size_um, dtype=np.float64) / 1000.0
        coordinates_px = np.column_stack(
            [
                u_offsets / pixel_size_mm[0] + self.pattern_center_px[0],
                v_offsets / pixel_size_mm[1] + self.pattern_center_px[1],
            ]
        )
        coordinates_px = np.ascontiguousarray(coordinates_px)
        coordinates_px.setflags(write=False)
        valid = np.ascontiguousarray(valid)
        valid.setflags(write=False)
        return coordinates_px, valid

    def scattering_vectors_lab(self, coordinates_px: np.ndarray) -> np.ndarray:
        outgoing = self.outgoing_directions_lab(coordinates_px)
        wavelength = self.electron_wavelength_angstrom
        incident_wavevector = self.beam_direction_lab[None, :] / wavelength
        scattered_wavevectors = outgoing / wavelength
        scattering = scattered_wavevectors - incident_wavevector
        scattering = np.ascontiguousarray(scattering)
        scattering.setflags(write=False)
        return scattering

    def two_theta_rad(self, coordinates_px: np.ndarray) -> np.ndarray:
        outgoing = self.outgoing_directions_lab(coordinates_px)
        cos_angles = np.clip(outgoing @ self.beam_direction_lab, -1.0, 1.0)
        angles = np.arccos(cos_angles)
        angles = np.ascontiguousarray(angles)
        angles.setflags(write=False)
        return angles

    def azimuth_rad(self, coordinates_px: np.ndarray) -> np.ndarray:
        outgoing = self.outgoing_directions_lab(coordinates_px)
        basis = self.detector_basis_lab
        u_coords = outgoing @ basis[:, 0]
        v_coords = outgoing @ basis[:, 1]
        azimuth = np.arctan2(v_coords, u_coords)
        azimuth = np.ascontiguousarray(azimuth)
        azimuth.setflags(write=False)
        return azimuth

    def bragg_two_theta_rad(self, d_spacing_angstrom: float) -> float:
        if d_spacing_angstrom <= 0.0:
            raise ValueError("d_spacing_angstrom must be strictly positive.")
        argument = self.electron_wavelength_angstrom / (2.0 * d_spacing_angstrom)
        if argument > 1.0 + _BRAGG_ARGUMENT_TOLERANCE:
            raise ValueError("Bragg condition cannot be satisfied for the given spacing.")
        return float(2.0 * np.arcsin(np.clip(argument, -1.0, 1.0)))

    def ring_radius_mm(self, two_theta_rad: float) -> float:
        if two_theta_rad < 0.0:
            raise ValueError("two_theta_rad must be non-negative.")
        return float(self.camera_length_mm * np.tan(two_theta_rad))

    def ring_radius_mm_for_d_spacing(self, d_spacing_angstrom: float) -> float:
        return self.ring_radius_mm(self.bragg_two_theta_rad(d_spacing_angstrom))

    def ring_radius_mm_for_plane(self, plane: CrystalPlane) -> float:
        return self.ring_radius_mm_for_d_spacing(plane.d_spacing_angstrom)


@dataclass(frozen=True, slots=True)
class DetectorAcceptanceMask:
    inset_px: tuple[float, float] = (0.0, 0.0)
    max_radius_px: float | None = None

    def __post_init__(self) -> None:
        if len(self.inset_px) != 2:
            raise ValueError("DetectorAcceptanceMask.inset_px must have length 2.")
        inset = tuple(float(value) for value in self.inset_px)
        if any(not np.isfinite(value) or value < 0.0 for value in inset):
            raise ValueError(
                "DetectorAcceptanceMask.inset_px values must be finite and non-negative."
            )
        if self.max_radius_px is not None:
            radius = float(self.max_radius_px)
            if not np.isfinite(radius) or radius <= 0.0:
                raise ValueError(
                    "DetectorAcceptanceMask.max_radius_px must be finite "
                    "and positive when provided."
                )
            object.__setattr__(self, "max_radius_px", radius)
        object.__setattr__(self, "inset_px", inset)

    def contains(self, geometry: DiffractionGeometry, coordinates_px: np.ndarray) -> np.ndarray:
        coordinates = as_float_array(coordinates_px, shape=(None, 2))
        finite = np.all(np.isfinite(coordinates), axis=1)
        inset_u, inset_v = self.inset_px
        min_u = inset_u
        max_u = float(geometry.detector_shape[1] - 1) - inset_u
        min_v = inset_v
        max_v = float(geometry.detector_shape[0] - 1) - inset_v
        within_rectangle = (
            finite
            & (coordinates[:, 0] >= min_u)
            & (coordinates[:, 0] <= max_u)
            & (coordinates[:, 1] >= min_v)
            & (coordinates[:, 1] <= max_v)
        )
        if self.max_radius_px is None:
            accepted = within_rectangle
        else:
            radial_distance = np.linalg.norm(
                coordinates - geometry.pattern_center_px[None, :], axis=1
            )
            accepted = within_rectangle & (radial_distance <= self.max_radius_px)
        accepted = np.ascontiguousarray(accepted)
        accepted.setflags(write=False)
        return accepted


@dataclass(frozen=True, slots=True)
class ReflectionFamily:
    representative_miller_indices: np.ndarray
    member_miller_indices: np.ndarray
    representative_spot_index: int
    spot_indices: np.ndarray
    multiplicity: int
    total_intensity: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "representative_miller_indices",
            as_int_array(self.representative_miller_indices, shape=(3,)),
        )
        object.__setattr__(
            self,
            "member_miller_indices",
            as_int_array(self.member_miller_indices, shape=(None, 3)),
        )
        object.__setattr__(self, "spot_indices", as_int_array(self.spot_indices, shape=(None,)))
        if self.multiplicity <= 0:
            raise ValueError("ReflectionFamily.multiplicity must be positive.")
        if (
            self.spot_indices.shape[0] == 0
            or self.spot_indices.shape[0] > self.member_miller_indices.shape[0]
        ):
            raise ValueError(
                "ReflectionFamily.spot_indices must contain between one and "
                "the full set of member indices."
            )
        if not np.isfinite(self.total_intensity) or self.total_intensity < 0.0:
            raise ValueError("ReflectionFamily.total_intensity must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class KinematicSpot:
    miller_indices: np.ndarray
    reciprocal_vector_lab: np.ndarray
    outgoing_direction_lab: np.ndarray
    detector_coordinates_px: np.ndarray
    excitation_error_inv_angstrom: float
    intensity: float
    two_theta_rad: float
    azimuth_rad: float
    on_detector: bool
    accepted_by_mask: bool = True
    family_id: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        object.__setattr__(
            self,
            "reciprocal_vector_lab",
            as_float_array(self.reciprocal_vector_lab, shape=(3,)),
        )
        object.__setattr__(
            self,
            "outgoing_direction_lab",
            as_float_array(self.outgoing_direction_lab, shape=(3,)),
        )
        object.__setattr__(
            self,
            "detector_coordinates_px",
            as_float_array(self.detector_coordinates_px, shape=(2,)),
        )
        if np.any(~np.isfinite(self.reciprocal_vector_lab)):
            raise ValueError("KinematicSpot.reciprocal_vector_lab must be finite.")
        if np.any(~np.isfinite(self.outgoing_direction_lab)):
            raise ValueError("KinematicSpot.outgoing_direction_lab must be finite.")
        if not np.isfinite(self.excitation_error_inv_angstrom):
            raise ValueError("KinematicSpot.excitation_error_inv_angstrom must be finite.")
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("KinematicSpot.intensity must be finite and non-negative.")
        if not np.isfinite(self.two_theta_rad) or self.two_theta_rad < 0.0:
            raise ValueError("KinematicSpot.two_theta_rad must be finite and non-negative.")
        if not np.isfinite(self.azimuth_rad):
            raise ValueError("KinematicSpot.azimuth_rad must be finite.")
        if self.accepted_by_mask and not self.on_detector:
            raise ValueError(
                "KinematicSpot.accepted_by_mask cannot be true for off-detector spots."
            )
        if self.on_detector:
            if np.any(~np.isfinite(self.detector_coordinates_px)):
                raise ValueError(
                    "KinematicSpot.detector_coordinates_px must be finite for on-detector spots."
                )


@dataclass(frozen=True, slots=True)
class KinematicSimulation:
    geometry: DiffractionGeometry
    phase: Phase
    spots: tuple[KinematicSpot, ...]
    reflection_families: tuple[ReflectionFamily, ...] = ()
    orientation: Orientation | None = None
    zone_axis: ZoneAxis | None = None
    provenance: ProvenanceRecord | None = None

    def accepted_spots(self) -> tuple[KinematicSpot, ...]:
        return tuple(spot for spot in self.spots if spot.accepted_by_mask)

    @classmethod
    def rank_orientation_candidates(
        cls,
        geometry: DiffractionGeometry,
        phase: Phase,
        pattern: DiffractionPattern,
        miller_indices: np.ndarray,
        candidate_orientations: list[Orientation] | OrientationSet,
        *,
        zone_axis: ZoneAxis | None = None,
        max_excitation_error_inv_angstrom: float = 5e-2,
        intensity_model: str = "kinematic_proxy",
        excitation_sigma_inv_angstrom: float = 5e-2,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        max_distance_px: float = 10.0,
        cluster_radius_px: float = 5.0,
        use_only_accepted: bool = True,
    ) -> tuple[OrientationIndexingCandidate, ...]:
        if isinstance(candidate_orientations, OrientationSet):
            orientations = [
                candidate_orientations[index] for index in range(len(candidate_orientations))
            ]
        else:
            orientations = list(candidate_orientations)
        candidates: list[OrientationIndexingCandidate] = []
        for orientation_index, orientation in enumerate(orientations):
            simulation = cls.simulate_spots(
                geometry,
                phase,
                miller_indices,
                orientation=orientation,
                zone_axis=zone_axis,
                max_excitation_error_inv_angstrom=max_excitation_error_inv_angstrom,
                intensity_model=intensity_model,
                excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
                acceptance_mask=acceptance_mask,
            )
            indexing = simulation.associate_to_pattern(
                pattern,
                max_distance_px=max_distance_px,
                cluster_radius_px=cluster_radius_px,
                use_only_accepted=use_only_accepted,
            )
            candidates.append(
                OrientationIndexingCandidate(
                    orientation_index=orientation_index,
                    orientation=orientation,
                    indexing=indexing,
                )
            )
        candidates = sorted(
            candidates,
            key=lambda candidate: (
                -candidate.score,
                candidate.indexing.mean_residual_px,
                -candidate.indexing.match_fraction,
            ),
        )
        return tuple(candidates)

    @classmethod
    def refine_orientation_candidate(
        cls,
        geometry: DiffractionGeometry,
        phase: Phase,
        pattern: DiffractionPattern,
        miller_indices: np.ndarray,
        initial_orientation: Orientation,
        *,
        zone_axis: ZoneAxis | None = None,
        max_excitation_error_inv_angstrom: float = 5e-2,
        intensity_model: str = "kinematic_proxy",
        excitation_sigma_inv_angstrom: float = 5e-2,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        max_distance_px: float = 10.0,
        cluster_radius_px: float = 5.0,
        use_only_accepted: bool = True,
        search_half_width_deg: float = 2.0,
        step_deg: float = 1.0,
        iterations: int = 2,
    ) -> OrientationRefinementResult:
        if search_half_width_deg < 0.0:
            raise ValueError("search_half_width_deg must be non-negative.")
        if step_deg <= 0.0:
            raise ValueError("step_deg must be strictly positive.")
        if iterations <= 0:
            raise ValueError("iterations must be strictly positive.")
        best_orientation = initial_orientation
        best_simulation = cls.simulate_spots(
            geometry,
            phase,
            miller_indices,
            orientation=best_orientation,
            zone_axis=zone_axis,
            max_excitation_error_inv_angstrom=max_excitation_error_inv_angstrom,
            intensity_model=intensity_model,
            excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
            acceptance_mask=acceptance_mask,
        )
        best_indexing = best_simulation.associate_to_pattern(
            pattern,
            max_distance_px=max_distance_px,
            cluster_radius_px=cluster_radius_px,
            use_only_accepted=use_only_accepted,
        )
        evaluated_candidates = 1
        current_half_width = float(search_half_width_deg)
        current_step = float(step_deg)
        for _ in range(iterations):
            center_phi1, center_phi, center_phi2 = best_orientation.rotation.to_bunge_euler()
            offsets = np.arange(
                -current_half_width,
                current_half_width + 0.5 * current_step,
                current_step,
                dtype=np.float64,
            )
            candidates: list[OrientationIndexingCandidate] = []
            for delta_phi1 in offsets:
                for delta_phi in offsets:
                    for delta_phi2 in offsets:
                        candidate_orientation = Orientation(
                            rotation=Rotation.from_bunge_euler(
                                center_phi1 + float(delta_phi1),
                                center_phi + float(delta_phi),
                                center_phi2 + float(delta_phi2),
                            ),
                            crystal_frame=initial_orientation.crystal_frame,
                            specimen_frame=initial_orientation.specimen_frame,
                            symmetry=initial_orientation.symmetry,
                            phase=initial_orientation.phase,
                            provenance=initial_orientation.provenance,
                        )
                        candidate_simulation = cls.simulate_spots(
                            geometry,
                            phase,
                            miller_indices,
                            orientation=candidate_orientation,
                            zone_axis=zone_axis,
                            max_excitation_error_inv_angstrom=max_excitation_error_inv_angstrom,
                            intensity_model=intensity_model,
                            excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
                            acceptance_mask=acceptance_mask,
                        )
                        candidate_indexing = candidate_simulation.associate_to_pattern(
                            pattern,
                            max_distance_px=max_distance_px,
                            cluster_radius_px=cluster_radius_px,
                            use_only_accepted=use_only_accepted,
                        )
                        candidates.append(
                            OrientationIndexingCandidate(
                                orientation_index=-1,
                                orientation=candidate_orientation,
                                indexing=candidate_indexing,
                            )
                        )
            evaluated_candidates += len(candidates)
            best_candidate = sorted(
                candidates,
                key=lambda candidate: (
                    -candidate.score,
                    candidate.indexing.mean_residual_px,
                    -candidate.indexing.match_fraction,
                ),
            )[0]
            best_orientation = best_candidate.orientation
            best_indexing = best_candidate.indexing
            current_half_width *= 0.5
            current_step *= 0.5
        return OrientationRefinementResult(
            seed_orientation=initial_orientation,
            refined_candidate=OrientationIndexingCandidate(
                orientation_index=-1,
                orientation=best_orientation,
                indexing=best_indexing,
            ),
            evaluated_candidates=evaluated_candidates,
            iterations=iterations,
            initial_search_half_width_deg=float(search_half_width_deg),
            final_step_deg=float(current_step),
        )

    def associate_to_pattern(
        self,
        pattern: DiffractionPattern,
        *,
        max_distance_px: float = 10.0,
        cluster_radius_px: float = 5.0,
        use_only_accepted: bool = True,
    ) -> IndexingCandidate:
        if pattern.geometry != self.geometry:
            raise ValueError("pattern.geometry must match simulation.geometry.")
        if pattern.phase != self.phase:
            raise ValueError("pattern.phase must match simulation.phase.")
        if max_distance_px <= 0.0:
            raise ValueError("max_distance_px must be strictly positive.")
        clusters = pattern.cluster_observations(max_distance_px=cluster_radius_px)
        candidate_spots = self.accepted_spots() if use_only_accepted else self.spots
        remaining_spots = set(range(len(candidate_spots)))
        matches: list[SpotAssignment] = []
        for cluster in clusters:
            if not remaining_spots:
                break
            best_spot_index = min(
                remaining_spots,
                key=lambda spot_index: float(
                    np.linalg.norm(
                        candidate_spots[spot_index].detector_coordinates_px - cluster.center_px
                    )
                ),
            )
            residual = float(
                np.linalg.norm(
                    candidate_spots[best_spot_index].detector_coordinates_px - cluster.center_px
                )
            )
            if residual <= max_distance_px:
                remaining_spots.remove(best_spot_index)
                matches.append(
                    SpotAssignment(
                        observed_cluster_id=cluster.cluster_id,
                        simulated_spot_index=best_spot_index,
                        residual_px=residual,
                        family_id=candidate_spots[best_spot_index].family_id,
                    )
                )
        matched_cluster_ids = {match.observed_cluster_id for match in matches}
        unmatched_clusters = np.array(
            [
                cluster.cluster_id
                for cluster in clusters
                if cluster.cluster_id not in matched_cluster_ids
            ],
            dtype=np.int64,
        )
        unmatched_spots = np.array(sorted(remaining_spots), dtype=np.int64)
        return IndexingCandidate(
            pattern=pattern,
            simulation=self,
            observation_clusters=tuple(clusters),
            matches=tuple(matches),
            unmatched_observed_cluster_ids=unmatched_clusters,
            unmatched_simulated_spot_indices=unmatched_spots,
        )

    @classmethod
    def simulate_spots(
        cls,
        geometry: DiffractionGeometry,
        phase: Phase,
        miller_indices: np.ndarray,
        *,
        orientation: Orientation | None = None,
        zone_axis: ZoneAxis | None = None,
        max_excitation_error_inv_angstrom: float = 5e-2,
        intensity_model: str = "kinematic_proxy",
        excitation_sigma_inv_angstrom: float = 5e-2,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        deduplicate_families: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> KinematicSimulation:
        miller_array = as_float_array(miller_indices, shape=(None, 3))
        if max_excitation_error_inv_angstrom < 0.0:
            raise ValueError("max_excitation_error_inv_angstrom must be non-negative.")
        if excitation_sigma_inv_angstrom <= 0.0:
            raise ValueError("excitation_sigma_inv_angstrom must be strictly positive.")
        rounded_miller = np.rint(miller_array)
        if not np.allclose(miller_array, rounded_miller, atol=1e-12):
            raise ValueError("miller_indices must contain integer-valued triplets.")
        if orientation is not None:
            if orientation.phase is not None and orientation.phase != phase:
                raise ValueError("orientation.phase must match phase when specified.")
            if orientation.crystal_frame != phase.crystal_frame:
                raise ValueError("orientation.crystal_frame must match phase.crystal_frame.")
            if orientation.specimen_frame != geometry.specimen_frame:
                raise ValueError("orientation.specimen_frame must match geometry.specimen_frame.")
        if zone_axis is not None and zone_axis.phase != phase:
            raise ValueError("zone_axis.phase must match phase.")

        incident = geometry.incident_wavevector_lab
        incident_magnitude = float(np.linalg.norm(incident))
        zone_axis_vector = zone_axis.unit_vector if zone_axis is not None else None
        candidate_spots: list[_CandidateSpot] = []
        for miller_triplet in miller_array:
            miller_triplet_int = np.rint(miller_triplet).astype(np.int64)
            reciprocal_vector = ReciprocalLatticeVector.from_miller_index(
                MillerIndex(miller_triplet_int, phase=phase)
            ).cartesian_vector
            if orientation is not None:
                reciprocal_vector_specimen = orientation.map_crystal_vector(reciprocal_vector)
            else:
                reciprocal_vector_specimen = reciprocal_vector
            if zone_axis_vector is not None:
                zone_axis_dot = float(np.dot(zone_axis_vector, reciprocal_vector_specimen))
                zone_axis_scale = max(1.0, float(np.linalg.norm(reciprocal_vector_specimen)))
                if not np.isclose(
                    zone_axis_dot,
                    0.0,
                    atol=_ZONE_AXIS_ORTHOGONALITY_ATOL * zone_axis_scale,
                    rtol=_ZONE_AXIS_ORTHOGONALITY_RTOL,
                ):
                    continue
            reciprocal_vector_lab = geometry.specimen_vectors_to_lab(
                reciprocal_vector_specimen[None, :]
            )[0]
            outgoing_wavevector = incident + reciprocal_vector_lab
            excitation_error = float(np.linalg.norm(outgoing_wavevector) - incident_magnitude)
            if abs(excitation_error) > max_excitation_error_inv_angstrom:
                continue
            outgoing_direction = normalize_vector(outgoing_wavevector)
            coordinates_px, valid = geometry.project_directions_to_detector_px(
                outgoing_direction[None, :]
            )
            two_theta = float(
                np.arccos(
                    np.clip(
                        float(np.dot(outgoing_direction, geometry.beam_direction_lab)), -1.0, 1.0
                    )
                )
            )
            basis = geometry.detector_basis_lab
            azimuth = float(
                np.arctan2(
                    float(np.dot(outgoing_direction, basis[:, 1])),
                    float(np.dot(outgoing_direction, basis[:, 0])),
                )
            )
            finite_coordinates = bool(np.all(np.isfinite(coordinates_px[0])))
            on_detector = (
                bool(valid[0])
                and finite_coordinates
                and 0.0 <= coordinates_px[0, 0] <= geometry.detector_shape[1] - 1
                and 0.0 <= coordinates_px[0, 1] <= geometry.detector_shape[0] - 1
            )
            if acceptance_mask is None:
                accepted_by_mask = on_detector
            else:
                accepted_by_mask = (
                    bool(acceptance_mask.contains(geometry, coordinates_px)[0]) and on_detector
                )
            intensity = _kinematic_intensity(
                reciprocal_vector_lab,
                excitation_error,
                model=intensity_model,
                excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
            )
            candidate_spots.append(
                {
                    "miller_indices": miller_triplet_int,
                    "reciprocal_vector_lab": reciprocal_vector_lab,
                    "outgoing_direction_lab": outgoing_direction,
                    "detector_coordinates_px": coordinates_px[0],
                    "excitation_error_inv_angstrom": excitation_error,
                    "intensity": intensity,
                    "two_theta_rad": two_theta,
                    "azimuth_rad": azimuth,
                    "on_detector": on_detector,
                    "accepted_by_mask": accepted_by_mask,
                    "family_key": _reflection_family_key(miller_triplet_int, phase),
                }
            )
        family_index_by_key: dict[tuple[float, ...], int] = {}
        family_members: list[list[int]] = []
        for index, spot in enumerate(candidate_spots):
            family_key = spot["family_key"]
            if family_key not in family_index_by_key:
                family_index_by_key[family_key] = len(family_members)
                family_members.append([])
            family_members[family_index_by_key[family_key]].append(index)

        family_id_by_index = np.full(len(candidate_spots), -1, dtype=np.int64)
        reflection_families: list[ReflectionFamily] = []
        representative_indices: list[int] = []
        for family_id, spot_indices in enumerate(family_members):
            member_rows = [candidate_spots[index] for index in spot_indices]
            for index in spot_indices:
                family_id_by_index[index] = family_id
            representative_index = max(
                spot_indices,
                key=lambda idx: (
                    float(candidate_spots[idx]["intensity"]),
                    -abs(float(candidate_spots[idx]["excitation_error_inv_angstrom"])),
                    -float(np.linalg.norm(candidate_spots[idx]["miller_indices"])),
                ),
            )
            representative_indices.append(representative_index)
            reflection_families.append(
                ReflectionFamily(
                    representative_miller_indices=candidate_spots[representative_index][
                        "miller_indices"
                    ],
                    member_miller_indices=np.stack(
                        [row["miller_indices"] for row in member_rows],
                        axis=0,
                    ),
                    representative_spot_index=representative_index,
                    spot_indices=np.array(spot_indices, dtype=np.int64),
                    multiplicity=len(spot_indices),
                    total_intensity=sum(float(row["intensity"]) for row in member_rows),
                )
            )
        instantiated_spots = [
            KinematicSpot(
                miller_indices=spot["miller_indices"],
                reciprocal_vector_lab=spot["reciprocal_vector_lab"],
                outgoing_direction_lab=spot["outgoing_direction_lab"],
                detector_coordinates_px=spot["detector_coordinates_px"],
                excitation_error_inv_angstrom=spot["excitation_error_inv_angstrom"],
                intensity=spot["intensity"],
                two_theta_rad=spot["two_theta_rad"],
                azimuth_rad=spot["azimuth_rad"],
                on_detector=spot["on_detector"],
                accepted_by_mask=spot["accepted_by_mask"],
                family_id=int(family_id_by_index[index]),
            )
            for index, spot in enumerate(candidate_spots)
        ]
        if deduplicate_families:
            selected = set(representative_indices)
            spots = tuple(
                spot for index, spot in enumerate(instantiated_spots) if index in selected
            )
            reflection_families = [
                ReflectionFamily(
                    representative_miller_indices=family.representative_miller_indices,
                    member_miller_indices=family.member_miller_indices,
                    representative_spot_index=family_id,
                    spot_indices=np.array([family_id], dtype=np.int64),
                    multiplicity=family.multiplicity,
                    total_intensity=family.total_intensity,
                )
                for family_id, family in enumerate(reflection_families)
            ]
        else:
            spots = tuple(instantiated_spots)
        return cls(
            geometry=geometry,
            phase=phase,
            spots=spots,
            reflection_families=tuple(reflection_families),
            orientation=orientation,
            zone_axis=zone_axis,
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class DiffractionPattern:
    coordinates_px: np.ndarray
    intensities: np.ndarray
    geometry: DiffractionGeometry
    phase: Phase
    orientation: Orientation | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        coordinates = as_float_array(self.coordinates_px, shape=(None, 2))
        intensities = as_float_array(self.intensities, shape=(coordinates.shape[0],))
        if np.any(~np.isfinite(intensities)) or np.any(intensities < 0.0):
            raise ValueError("DiffractionPattern intensities must be finite and non-negative.")
        if self.orientation is not None:
            if self.orientation.specimen_frame != self.geometry.specimen_frame:
                raise ValueError(
                    "DiffractionPattern.orientation.specimen_frame must match "
                    "geometry.specimen_frame."
                )
            if self.orientation.crystal_frame != self.phase.crystal_frame:
                raise ValueError(
                    "DiffractionPattern.orientation.crystal_frame must match phase.crystal_frame."
                )
            if self.orientation.phase is not None and self.orientation.phase != self.phase:
                raise ValueError("DiffractionPattern.orientation.phase must match phase.")
        object.__setattr__(self, "coordinates_px", coordinates)
        object.__setattr__(self, "intensities", intensities)

    def detector_coordinates_mm(self) -> np.ndarray:
        return self.geometry.detector_coordinates_mm(self.coordinates_px)

    def outgoing_directions_lab(self) -> np.ndarray:
        return self.geometry.outgoing_directions_lab(self.coordinates_px)

    def scattering_vectors_lab(self) -> np.ndarray:
        return self.geometry.scattering_vectors_lab(self.coordinates_px)

    def two_theta_rad(self) -> np.ndarray:
        return self.geometry.two_theta_rad(self.coordinates_px)

    def azimuth_rad(self) -> np.ndarray:
        return self.geometry.azimuth_rad(self.coordinates_px)

    def cluster_observations(
        self, *, max_distance_px: float = 5.0
    ) -> tuple[DetectedSpotCluster, ...]:
        if max_distance_px <= 0.0:
            raise ValueError("max_distance_px must be strictly positive.")
        unassigned = set(range(self.coordinates_px.shape[0]))
        clusters: list[DetectedSpotCluster] = []
        cluster_id = 0
        while unassigned:
            seed = min(unassigned)
            members = {seed}
            frontier = [seed]
            unassigned.remove(seed)
            while frontier:
                current = frontier.pop()
                current_coordinate = self.coordinates_px[current]
                close_indices = [
                    candidate
                    for candidate in list(unassigned)
                    if np.linalg.norm(self.coordinates_px[candidate] - current_coordinate)
                    <= max_distance_px
                ]
                for candidate in close_indices:
                    unassigned.remove(candidate)
                    members.add(candidate)
                    frontier.append(candidate)
            member_indices = np.array(sorted(members), dtype=np.int64)
            weights = self.intensities[member_indices]
            if np.isclose(float(np.sum(weights)), 0.0):
                center = np.mean(self.coordinates_px[member_indices], axis=0)
            else:
                center = np.average(self.coordinates_px[member_indices], axis=0, weights=weights)
            clusters.append(
                DetectedSpotCluster(
                    cluster_id=cluster_id,
                    member_indices=member_indices,
                    center_px=center,
                    total_intensity=float(np.sum(weights)),
                )
            )
            cluster_id += 1
        return tuple(clusters)


@dataclass(frozen=True, slots=True)
class DetectedSpotCluster:
    cluster_id: int
    member_indices: np.ndarray
    center_px: np.ndarray
    total_intensity: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "member_indices", as_int_array(self.member_indices, shape=(None,)))
        object.__setattr__(self, "center_px", as_float_array(self.center_px, shape=(2,)))
        if self.member_indices.size == 0:
            raise ValueError("DetectedSpotCluster.member_indices must be non-empty.")
        if not np.isfinite(self.total_intensity) or self.total_intensity < 0.0:
            raise ValueError("DetectedSpotCluster.total_intensity must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class SpotAssignment:
    observed_cluster_id: int
    simulated_spot_index: int
    residual_px: float
    family_id: int | None = None

    def __post_init__(self) -> None:
        if self.observed_cluster_id < 0:
            raise ValueError("SpotAssignment.observed_cluster_id must be non-negative.")
        if self.simulated_spot_index < 0:
            raise ValueError("SpotAssignment.simulated_spot_index must be non-negative.")
        if not np.isfinite(self.residual_px) or self.residual_px < 0.0:
            raise ValueError("SpotAssignment.residual_px must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class IndexingCandidate:
    pattern: DiffractionPattern
    simulation: KinematicSimulation
    observation_clusters: tuple[DetectedSpotCluster, ...]
    matches: tuple[SpotAssignment, ...]
    unmatched_observed_cluster_ids: np.ndarray
    unmatched_simulated_spot_indices: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "observation_clusters", tuple(self.observation_clusters))
        object.__setattr__(self, "matches", tuple(self.matches))
        object.__setattr__(
            self,
            "unmatched_observed_cluster_ids",
            as_int_array(self.unmatched_observed_cluster_ids, shape=(None,)),
        )
        object.__setattr__(
            self,
            "unmatched_simulated_spot_indices",
            as_int_array(self.unmatched_simulated_spot_indices, shape=(None,)),
        )

    @property
    def match_fraction(self) -> float:
        if not self.observation_clusters:
            return 0.0
        return float(len(self.matches) / len(self.observation_clusters))

    @property
    def mean_residual_px(self) -> float:
        if not self.matches:
            return np.inf
        return float(np.mean([match.residual_px for match in self.matches]))

    @property
    def score(self) -> float:
        if not self.matches:
            return 0.0
        residual_penalty = 1.0 / (1.0 + self.mean_residual_px)
        return float(self.match_fraction * residual_penalty)

    def family_reports(self) -> tuple[FamilyIndexingReport, ...]:
        if not self.simulation.reflection_families:
            return ()
        matches_by_family: dict[int, list[SpotAssignment]] = {}
        for match in self.matches:
            if match.family_id is None:
                continue
            matches_by_family.setdefault(match.family_id, []).append(match)
        reports: list[FamilyIndexingReport] = []
        for family_id, family in enumerate(self.simulation.reflection_families):
            family_matches = matches_by_family.get(family_id, [])
            observed_cluster_ids = np.array(
                [match.observed_cluster_id for match in family_matches],
                dtype=np.int64,
            )
            mean_residual_px = (
                float(np.mean([match.residual_px for match in family_matches]))
                if family_matches
                else np.inf
            )
            simulated_spot_count = int(family.spot_indices.shape[0])
            reports.append(
                FamilyIndexingReport(
                    family_id=family_id,
                    representative_miller_indices=family.representative_miller_indices,
                    multiplicity=family.multiplicity,
                    simulated_spot_count=simulated_spot_count,
                    matched_spot_count=len(family_matches),
                    matched_fraction=float(
                        len(family_matches) / simulated_spot_count if simulated_spot_count else 0.0
                    ),
                    total_family_intensity=family.total_intensity,
                    mean_residual_px=mean_residual_px,
                    observed_cluster_ids=observed_cluster_ids,
                )
            )
        return tuple(reports)


@dataclass(frozen=True, slots=True)
class OrientationIndexingCandidate:
    orientation_index: int
    orientation: Orientation
    indexing: IndexingCandidate

    @property
    def score(self) -> float:
        return self.indexing.score


@dataclass(frozen=True, slots=True)
class FamilyIndexingReport:
    family_id: int
    representative_miller_indices: np.ndarray
    multiplicity: int
    simulated_spot_count: int
    matched_spot_count: int
    matched_fraction: float
    total_family_intensity: float
    mean_residual_px: float
    observed_cluster_ids: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "representative_miller_indices",
            as_int_array(self.representative_miller_indices, shape=(3,)),
        )
        object.__setattr__(
            self,
            "observed_cluster_ids",
            as_int_array(self.observed_cluster_ids, shape=(None,)),
        )
        if self.family_id < 0:
            raise ValueError("FamilyIndexingReport.family_id must be non-negative.")
        if self.multiplicity <= 0:
            raise ValueError("FamilyIndexingReport.multiplicity must be positive.")
        if self.simulated_spot_count <= 0:
            raise ValueError("FamilyIndexingReport.simulated_spot_count must be positive.")
        if self.matched_spot_count < 0 or self.matched_spot_count > self.simulated_spot_count:
            raise ValueError(
                "FamilyIndexingReport.matched_spot_count must lie between "
                "zero and simulated_spot_count."
            )
        if not np.isfinite(self.matched_fraction) or not (0.0 <= self.matched_fraction <= 1.0):
            raise ValueError("FamilyIndexingReport.matched_fraction must lie in [0, 1].")
        if not np.isfinite(self.total_family_intensity) or self.total_family_intensity < 0.0:
            raise ValueError(
                "FamilyIndexingReport.total_family_intensity must be finite and non-negative."
            )
        if not np.isfinite(self.mean_residual_px) and self.matched_spot_count > 0:
            raise ValueError(
                "FamilyIndexingReport.mean_residual_px must be finite when matches exist."
            )


@dataclass(frozen=True, slots=True)
class OrientationRefinementResult:
    seed_orientation: Orientation
    refined_candidate: OrientationIndexingCandidate
    evaluated_candidates: int
    iterations: int
    initial_search_half_width_deg: float
    final_step_deg: float

    def __post_init__(self) -> None:
        if self.evaluated_candidates <= 0:
            raise ValueError("OrientationRefinementResult.evaluated_candidates must be positive.")
        if self.iterations <= 0:
            raise ValueError("OrientationRefinementResult.iterations must be positive.")
        if self.initial_search_half_width_deg < 0.0:
            raise ValueError(
                "OrientationRefinementResult.initial_search_half_width_deg must be non-negative."
            )
        if self.final_step_deg <= 0.0:
            raise ValueError("OrientationRefinementResult.final_step_deg must be positive.")
