from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, TypedDict

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    as_int_array,
    is_rotation_matrix,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.acquisition import (
    AcquisitionGeometry,
    CalibrationRecord,
    MeasurementQuality,
    ScatteringSetup,
)
from pytex.core.conventions import FrameDomain
from pytex.core.frames import FrameTransform, ReferenceFrame
from pytex.core.lattice import (
    CrystalPlane,
    MillerIndex,
    Phase,
    ReciprocalLatticeVector,
    ZoneAxis,
)
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.shape_factors import FiniteThicknessShapeFactor

if TYPE_CHECKING:
    from pytex.adapters import ExperimentManifest

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
    foil_thickness_angstrom: float | None,
) -> float:
    if model not in {"unit", "kinematic_proxy"}:
        raise ValueError("intensity_model must be either 'unit' or 'kinematic_proxy'.")
    if foil_thickness_angstrom is not None:
        excitation_weight = float(
            FiniteThicknessShapeFactor(foil_thickness_angstrom).intensity_factor(
                excitation_error_inv_angstrom
            )
        )
    elif model == "unit":
        return 1.0
    else:
        excitation_ratio = excitation_error_inv_angstrom / excitation_sigma_inv_angstrom
        excitation_weight = 1.0 / (1.0 + excitation_ratio * excitation_ratio)
    if model == "unit":
        return excitation_weight
    reciprocal_magnitude = float(np.linalg.norm(reciprocal_vector_lab))
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
    """The complete geometry of a diffraction experiment.

    Purpose
    -------
    The single authority on where a diffracted beam lands and what angle it
    corresponds to. Detector position, orientation, pixel size, beam energy,
    and the specimen-to-laboratory relationship are held together so that no
    projection has to re-derive them, and so that frame domains stay
    separate: the detector, specimen, and laboratory frames are distinct and
    are checked as such at construction.

    Attributes
    ----------
    detector_frame, specimen_frame, laboratory_frame : ReferenceFrame
        Must belong to their respective domains; enforced.
    beam_energy_kev : float
        Accelerating voltage, from which the relativistically corrected
        electron wavelength follows.
    camera_length_mm : float
        Specimen-to-detector distance along the detector normal.
    pattern_center : np.ndarray
        ``(x, y, z)``; the in-plane components are fractions of the detector
        extent in ``[0, 1]``, and ``z`` is a positive distance fraction.
    detector_pixel_size_um : tuple of float
        Per-axis pixel size, so non-square pixels are handled correctly.
    detector_shape : tuple of int
        ``(height, width)`` in pixels.
    beam_direction_lab : np.ndarray
        Incident beam direction; normalized on construction.
    specimen_to_lab_matrix : np.ndarray
        Must be a proper rotation; enforced.
    tilt_degrees : tuple of float
        Detector tilt as an X-Y-Z rotation sequence.
    acquisition_geometry, calibration_record, measurement_quality,
    scattering_setup, provenance : optional
        Context records. When supplied they are cross-checked against the
        explicit fields above, so the two cannot disagree silently.
    """

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
    acquisition_geometry: AcquisitionGeometry | None = None
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    scattering_setup: ScatteringSetup | None = None
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
        if self.acquisition_geometry is not None:
            if self.acquisition_geometry.specimen_frame != self.specimen_frame:
                raise ValueError(
                    "DiffractionGeometry.acquisition_geometry.specimen_frame must match "
                    "DiffractionGeometry.specimen_frame."
                )
            if self.acquisition_geometry.detector_frame != self.detector_frame:
                raise ValueError(
                    "DiffractionGeometry.acquisition_geometry.detector_frame must match "
                    "DiffractionGeometry.detector_frame."
                )
            if self.acquisition_geometry.laboratory_frame != self.laboratory_frame:
                raise ValueError(
                    "DiffractionGeometry.acquisition_geometry.laboratory_frame must match "
                    "DiffractionGeometry.laboratory_frame."
                )
            if self.acquisition_geometry.specimen_to_laboratory is not None:
                transform = self.acquisition_geometry.specimen_to_laboratory
                if not np.allclose(transform.translation_vector, np.zeros(3), atol=1e-8):
                    raise ValueError(
                        "DiffractionGeometry currently requires a zero-translation "
                        "specimen_to_laboratory transform."
                    )
                if not np.allclose(
                    transform.rotation_matrix,
                    self.specimen_to_lab_matrix,
                    atol=1e-8,
                ):
                    raise ValueError(
                        "DiffractionGeometry.specimen_to_lab_matrix must match the "
                        "acquisition geometry specimen_to_laboratory transform."
                    )
            if (
                self.calibration_record is not None
                and self.acquisition_geometry.calibration_record is not None
                and self.calibration_record != self.acquisition_geometry.calibration_record
            ):
                raise ValueError(
                    "DiffractionGeometry.calibration_record must match the acquisition "
                    "geometry calibration record when both are provided."
                )
            if (
                self.measurement_quality is not None
                and self.acquisition_geometry.measurement_quality is not None
                and self.measurement_quality != self.acquisition_geometry.measurement_quality
            ):
                raise ValueError(
                    "DiffractionGeometry.measurement_quality must match the acquisition "
                    "geometry measurement quality when both are provided."
                )
        if self.scattering_setup is not None:
            if self.scattering_setup.laboratory_frame != self.laboratory_frame:
                raise ValueError(
                    "DiffractionGeometry.scattering_setup.laboratory_frame must match "
                    "DiffractionGeometry.laboratory_frame."
                )
            if self.scattering_setup.radiation_type not in {"electron", "generic"}:
                raise ValueError(
                    "DiffractionGeometry currently supports only electron or generic "
                    "scattering setups."
                )
            if not np.allclose(
                self.scattering_setup.incident_beam_direction,
                self.beam_direction_lab,
                atol=1e-8,
            ):
                raise ValueError(
                    "DiffractionGeometry.beam_direction_lab must match the scattering "
                    "setup incident beam direction."
                )
            if self.scattering_setup.beam_energy_kev is not None and not np.isclose(
                self.scattering_setup.beam_energy_kev,
                self.beam_energy_kev,
                atol=1e-8,
            ):
                raise ValueError(
                    "DiffractionGeometry.beam_energy_kev must match the scattering setup "
                    "beam energy when both are provided."
                )

    @classmethod
    def for_ebsd(
        cls,
        *,
        beam_energy_kev: float = 20.0,
        sample_tilt_deg: float = 70.0,
        detector_elevation_deg: float = 0.0,
        detector_azimuth_deg: float = 0.0,
        pattern_center: ArrayLike = (0.5, 0.5, 0.65),
        detector_shape: tuple[int, int] = (480, 640),
        detector_pixel_size_um: tuple[float, float] = (50.0, 50.0),
        detector_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame | None = None,
        laboratory_frame: ReferenceFrame | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> DiffractionGeometry:
        r"""Build the geometry of a backscatter-diffraction (EBSD) camera.

        Purpose
        -------
        State an EBSD setup in the terms it is actually configured in — how far
        the stage is tilted, where the camera sits, and where the pattern centre
        falls — and get back the same
        :class:`DiffractionGeometry` every other diffraction surface consumes.
        Without this the caller has to compose the stage rotation and the camera
        placement by hand, which is precisely the kind of private frame model
        this class exists to prevent.

        Convention
        ----------
        The laboratory frame is fixed by the column and the stage, not by the
        camera:

        - :math:`\hat{\mathbf{z}}_{\text{lab}}` is the **beam**, travelling
          from the gun to the specimen.
        - :math:`\hat{\mathbf{x}}_{\text{lab}}` is the **stage tilt axis**.
        - :math:`\hat{\mathbf{y}}_{\text{lab}}` completes the right-handed
          triad and points **away from the camera**, so the camera is on
          :math:`-\hat{\mathbf{y}}_{\text{lab}}`.

        An untilted specimen has its normal along
        :math:`-\hat{\mathbf{z}}_{\text{lab}}`, facing into the beam.
        ``sample_tilt_deg`` rotates the specimen about the tilt axis so that its
        normal tips *towards* the camera — the 70 degrees of a standard EBSD
        stage. ``detector_elevation_deg`` raises the camera's own axis above the
        plane perpendicular to the beam, and ``detector_azimuth_deg`` swings it
        about the beam, so a camera port that is not on the nominal side is
        described rather than approximated.

        Two consequences are worth stating because they are what the geometry is
        checked against:

        - the specimen normal makes an angle of
          :math:`90^\circ - (\sigma - \epsilon)` with the camera axis, so it
          projects at gnomonic radius
          :math:`	an(90^\circ - \sigma + \epsilon)` — about 0.364 at the
          standard 70 degrees, which is why the specimen normal falls on a real
          screen at all;
        - the beam makes :math:`90^\circ - \epsilon` with the camera axis, so
          at zero elevation it is parallel to the screen and has no gnomonic
          image, which is correct: an EBSD screen never sees the beam.

        Parameters
        ----------
        beam_energy_kev : float
            Accelerating voltage. 20 kV is the usual EBSD working point, an
            order of magnitude below TEM, and it is what sets the band widths.
        sample_tilt_deg : float
            Stage tilt, in degrees, about the tilt axis. Must lie in
            ``[0, 90)``: at 90 degrees the specimen surface contains the beam.
        detector_elevation_deg : float
            Elevation of the camera axis above the plane perpendicular to the
            beam, in degrees. Must lie in ``(-90, 90)``.
        detector_azimuth_deg : float
            Rotation of the camera about the beam, in degrees, from the nominal
            port.
        pattern_center : ArrayLike
            ``(x*, y*, z*)``. The in-plane pair are fractions of the detector
            extent, as elsewhere in this class; ``z*`` is the camera distance as
            a fraction of the detector **width**, which is the quantity every
            EBSD calibration reports. Vendors differ over the origin and the
            axis the fractions are taken along, so a value copied from a vendor
            file may need converting; PyTex takes them in its own stated frame
            rather than guessing which vendor wrote them.
        detector_shape : tuple of int
            ``(height, width)`` in pixels.
        detector_pixel_size_um : tuple of float
            Per-axis pixel pitch.
        detector_frame, specimen_frame, laboratory_frame : ReferenceFrame, optional
            Default to the catalogue frames.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        DiffractionGeometry
            Ready for :func:`pytex.diffraction.kikuchi.simulate_kikuchi_pattern`.

        Raises
        ------
        ValueError
            For a tilt or elevation outside its range, or a non-positive
            camera distance.

        See Also
        --------
        pytex.diffraction.kikuchi.simulate_kikuchi_pattern : What consumes it.
        """

        from pytex.core.frame_catalog import (
            DETECTOR_FRAME,
            LABORATORY_FRAME,
            SPECIMEN_FRAME,
        )

        tilt = float(sample_tilt_deg)
        elevation = float(detector_elevation_deg)
        azimuth = float(detector_azimuth_deg)
        if not 0.0 <= tilt < 90.0:
            raise ValueError("sample_tilt_deg must lie in [0, 90).")
        if not -90.0 < elevation < 90.0:
            raise ValueError("detector_elevation_deg must lie in (-90, 90).")
        centre = as_float_array(pattern_center, shape=(3,))
        shape = (int(detector_shape[0]), int(detector_shape[1]))
        width_mm = shape[1] * float(detector_pixel_size_um[1]) / 1000.0
        distance_mm = float(centre[2]) * width_mm
        if distance_mm <= 0.0:
            raise ValueError("pattern_center z component must be strictly positive.")
        return cls(
            detector_frame=detector_frame or DETECTOR_FRAME,
            specimen_frame=specimen_frame or SPECIMEN_FRAME,
            laboratory_frame=laboratory_frame or LABORATORY_FRAME,
            beam_energy_kev=float(beam_energy_kev),
            camera_length_mm=distance_mm,
            pattern_center=centre,
            detector_pixel_size_um=(
                float(detector_pixel_size_um[0]),
                float(detector_pixel_size_um[1]),
            ),
            detector_shape=shape,
            # The specimen normal faces the beam and tips towards the camera:
            # an untilted specimen is the half-turn about the tilt axis, and the
            # stage tilt takes it back from there.
            specimen_to_lab_matrix=_rotation_matrix_x(float(np.deg2rad(180.0 - tilt))),
            # The camera axis starts along the beam and is laid down onto the
            # port: 90 degrees less the elevation, then swung about the beam.
            tilt_degrees=(90.0 - elevation, 0.0, azimuth),
            provenance=provenance,
        )

    def to_experiment_manifest(
        self,
        *,
        source_system: str = "pytex",
        phase: Phase | None = None,
        referenced_files: tuple[str, ...] = (),
        metadata: dict[str, str] | None = None,
    ) -> ExperimentManifest:
        """Export this geometry as a schema-validated experiment manifest.

        Synthesizes a minimal acquisition geometry when none is attached, and
        records camera length, detector shape, and pixel size as metadata. The
        manifest states explicitly that the specimen-to-detector transform was
        taken as the nominal identity when no measured transform was available,
        so a downstream reader is never misled about calibration provenance.
        """

        from pytex.adapters import ExperimentManifest

        acquisition_geometry = self.acquisition_geometry
        if acquisition_geometry is None:
            acquisition_geometry = AcquisitionGeometry(
                specimen_frame=self.specimen_frame,
                modality="tem" if self.scattering_setup is not None else "generic",
                detector_frame=self.detector_frame,
                laboratory_frame=self.laboratory_frame,
                specimen_to_detector=FrameTransform(
                    source=self.specimen_frame,
                    target=self.detector_frame,
                    rotation_matrix=np.eye(3),
                    provenance=self.provenance,
                ),
                specimen_to_laboratory=FrameTransform(
                    source=self.specimen_frame,
                    target=self.laboratory_frame,
                    rotation_matrix=self.specimen_to_lab_matrix,
                    provenance=self.provenance,
                ),
                calibration_record=self.calibration_record,
                measurement_quality=self.measurement_quality,
                provenance=self.provenance,
            )
        merged_metadata = {
            "camera_length_mm": f"{self.camera_length_mm:g}",
            "detector_shape": "x".join(str(value) for value in self.detector_shape),
            "detector_pixel_size_um": ",".join(
                f"{value:g}" for value in self.detector_pixel_size_um
            ),
            "detector_alignment_contract": (
                "nominal_identity_when_no_explicit_specimen_to_detector_transform_is_available"
            ),
        }
        if metadata is not None:
            merged_metadata.update(metadata)
        return ExperimentManifest.from_acquisition_geometry(
            acquisition_geometry,
            source_system=source_system,
            phase=phase,
            scattering_setup=self.scattering_setup,
            referenced_files=referenced_files,
            metadata=merged_metadata,
        )

    @property
    def electron_wavelength_angstrom(self) -> float:
        """Relativistically corrected electron wavelength, in angstroms.

        Evaluates ``lambda = 12.2639 / sqrt(V (1 + 0.97845e-6 V))`` with ``V``
        the accelerating voltage in volts. The correction matters: at 200 kV the
        non-relativistic value is about 6 percent too large, which would
        misplace every simulated spot.
        """

        voltage = self.beam_energy_kev * 1000.0
        numerator = 12.2639
        denominator = np.sqrt(voltage * (1.0 + 0.97845e-6 * voltage))
        return float(numerator / denominator)

    @property
    def ewald_sphere_radius_inv_angstrom(self) -> float:
        """Radius of the Ewald sphere, ``1 / lambda``, in inverse angstroms.

        The scale that decides how nearly flat the sphere is over the pattern —
        the reason electron patterns show a whole zone at once while X-ray
        patterns do not.
        """

        return float(1.0 / self.electron_wavelength_angstrom)

    @property
    def incident_wavevector_lab(self) -> np.ndarray:
        """Incident wavevector ``k_0`` in the laboratory frame, in inverse angstroms.

        Magnitude ``1 / lambda`` along the beam direction. Read-only.
        """

        wavevector = self.beam_direction_lab / self.electron_wavelength_angstrom
        wavevector = np.ascontiguousarray(wavevector)
        wavevector.setflags(write=False)
        return wavevector

    @property
    def pattern_center_px(self) -> np.ndarray:
        """Pattern centre in pixel coordinates, as ``(u, v)``.

        Converted from the stored fractional pattern centre using the detector
        shape. Read-only.
        """

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
        """Detector axes expressed in the laboratory frame, as columns ``(u, v, n)``.

        Purpose
        -------
        The single definition of detector orientation that every projection and
        back-projection in this class goes through, so no other code needs to
        re-derive it.

        The third column is the detector normal. The in-plane axes are
        constructed orthogonal to the beam and then rotated by the configured
        X-Y-Z tilt sequence, so a tilted detector is handled without the caller
        composing rotations by hand. Read-only.
        """

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
        """Map specimen-frame vectors into the laboratory frame.

        Applies the stored specimen-to-laboratory rotation. Accepts any array
        ending in dimension 3 and returns a read-only array.
        """

        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Specimen vectors must end with dimension 3.")
        mapped = array @ self.specimen_to_lab_matrix.T
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def lab_vectors_to_specimen(self, vectors: np.ndarray) -> np.ndarray:
        """Map laboratory-frame vectors into the specimen frame.

        The inverse of :meth:`specimen_vectors_to_lab`.
        """

        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Laboratory vectors must end with dimension 3.")
        mapped = array @ self.specimen_to_lab_matrix
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def detector_coordinates_mm(self, coordinates_px: np.ndarray) -> np.ndarray:
        """Detector-plane offsets from the pattern centre, in millimetres.

        Converts pixel coordinates using the per-axis pixel size, so
        non-square pixels are handled correctly. Returns ``(n, 2)``, read-only.
        """

        detector_pixels = as_float_array(coordinates_px, shape=(None, 2))
        pixel_size_mm = np.array(self.detector_pixel_size_um, dtype=np.float64) / 1000.0
        offsets_mm = (detector_pixels - self.pattern_center_px[None, :]) * pixel_size_mm[None, :]
        offsets_mm = np.ascontiguousarray(offsets_mm)
        offsets_mm.setflags(write=False)
        return offsets_mm

    def detector_points_lab_mm(self, coordinates_px: np.ndarray) -> np.ndarray:
        """Laboratory-frame positions of detector pixels, in millimetres.

        Places each pixel on the detector plane at the configured camera length
        along the detector normal, offset by its in-plane coordinates. Returns
        ``(n, 3)``, read-only.
        """

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
        """Unit scattered-beam directions for the given detector pixels.

        The direction from the specimen to each pixel — the geometric input to
        :meth:`scattering_vectors_lab` and :meth:`two_theta_rad`.
        """

        directions = normalize_vectors(self.detector_points_lab_mm(coordinates_px))
        directions = np.ascontiguousarray(directions)
        directions.setflags(write=False)
        return directions

    def project_directions_to_detector_px(
        self, directions_lab: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Project laboratory directions onto the detector, in pixels.

        Purpose
        -------
        The forward projection every simulation uses: given a diffracted beam
        direction, where does it land on the detector?

        Parameters
        ----------
        directions_lab : np.ndarray
            ``(n, 3)`` directions in the laboratory frame; normalized
            internally.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(n, 2)`` pixel coordinates and an ``(n,)`` validity mask. Beams
            travelling away from the detector, or parallel to its plane, cannot
            intersect it: those rows are ``NaN`` and are flagged ``False``
            rather than silently projected to a spurious position.
        """

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
        """Scattering vectors ``q = k - k_0`` for the given pixels.

        In the laboratory frame and in inverse angstroms. The quantity compared
        against reciprocal-lattice vectors when a pattern is indexed.
        """

        outgoing = self.outgoing_directions_lab(coordinates_px)
        wavelength = self.electron_wavelength_angstrom
        incident_wavevector = self.beam_direction_lab[None, :] / wavelength
        scattered_wavevectors = outgoing / wavelength
        scattering = scattered_wavevectors - incident_wavevector
        scattering = np.ascontiguousarray(scattering)
        scattering.setflags(write=False)
        return scattering

    def two_theta_rad(self, coordinates_px: np.ndarray) -> np.ndarray:
        """Scattering angle ``2*theta`` at the given pixels, in radians.

        Measured from the incident beam direction.
        """

        outgoing = self.outgoing_directions_lab(coordinates_px)
        cos_angles = np.clip(outgoing @ self.beam_direction_lab, -1.0, 1.0)
        angles = np.arccos(cos_angles)
        angles = np.ascontiguousarray(angles)
        angles.setflags(write=False)
        return angles

    def azimuth_rad(self, coordinates_px: np.ndarray) -> np.ndarray:
        """Azimuthal angle of the given pixels about the beam, in radians.

        Measured in the detector plane from the detector ``u`` axis, in
        ``(-pi, pi]``. Together with :meth:`two_theta_rad` this gives the polar
        coordinates used for ring integration and texture-resolved analysis.
        """

        outgoing = self.outgoing_directions_lab(coordinates_px)
        basis = self.detector_basis_lab
        u_coords = outgoing @ basis[:, 0]
        v_coords = outgoing @ basis[:, 1]
        azimuth = np.arctan2(v_coords, u_coords)
        azimuth = np.ascontiguousarray(azimuth)
        azimuth.setflags(write=False)
        return azimuth

    def bragg_two_theta_rad(self, d_spacing_angstrom: float) -> float:
        """Bragg angle ``2*theta`` for a given interplanar spacing, in radians.

        Solves ``lambda = 2 d sin(theta)`` at first order using this geometry's
        electron wavelength. Raises when the spacing is too small to satisfy the
        Bragg condition at this wavelength, rather than returning a clipped
        angle that would look plausible and be wrong.
        """

        if d_spacing_angstrom <= 0.0:
            raise ValueError("d_spacing_angstrom must be strictly positive.")
        argument = self.electron_wavelength_angstrom / (2.0 * d_spacing_angstrom)
        if argument > 1.0 + _BRAGG_ARGUMENT_TOLERANCE:
            raise ValueError("Bragg condition cannot be satisfied for the given spacing.")
        return float(2.0 * np.arcsin(np.clip(argument, -1.0, 1.0)))

    def ring_radius_mm(self, two_theta_rad: float) -> float:
        """Detector radius of a diffraction ring at a given ``2*theta``.

        ``L tan(2 theta)`` for camera length ``L``. Note the tangent: the common
        small-angle form ``L * 2 theta`` is an approximation, and this is not
        it.
        """

        if two_theta_rad < 0.0:
            raise ValueError("two_theta_rad must be non-negative.")
        return float(self.camera_length_mm * np.tan(two_theta_rad))

    def ring_radius_mm_for_d_spacing(self, d_spacing_angstrom: float) -> float:
        """Detector ring radius for a given interplanar spacing, in millimetres.

        Composes :meth:`bragg_two_theta_rad` with :meth:`ring_radius_mm` — the
        calculation behind camera-constant calibration, since the product of
        ring radius and d-spacing is the camera constant.
        """

        return self.ring_radius_mm(self.bragg_two_theta_rad(d_spacing_angstrom))

    def ring_radius_mm_for_plane(self, plane: CrystalPlane) -> float:
        """Detector ring radius for a given crystal plane, in millimetres.

        The typed form of :meth:`ring_radius_mm_for_d_spacing`.
        """

        return self.ring_radius_mm_for_d_spacing(plane.d_spacing_angstrom)


@dataclass(frozen=True, slots=True)
class DetectorAcceptanceMask:
    """Which region of a detector counts as usable.

    Purpose
    -------
    A simulated reflection may project outside the physical detector, onto
    its unusable border, or beyond a useful radius. Declaring acceptance once
    means simulation, indexing, and scoring all agree on which spots were
    observable in principle — the distinction between "not predicted" and
    "predicted but unobservable".

    Attributes
    ----------
    inset_px : tuple of float
        Border to exclude on each axis, in pixels.
    max_radius_px : float, optional
        Radial cut-off from the pattern centre. ``None`` means no radial
        limit.
    """

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
        """Which detector coordinates the mask accepts.

        Purpose
        -------
        A simulated reflection can project outside the physical detector, onto
        its unusable border, or beyond the useful radius. This decides
        acceptance once so simulation, indexing, and scoring all agree on which
        spots were observable in principle.

        Parameters
        ----------
        geometry : DiffractionGeometry
            Supplies detector shape and pattern centre.
        coordinates_px : np.ndarray
            ``(n, 2)`` pixel coordinates.

        Returns
        -------
        np.ndarray
            ``(n,)`` boolean acceptance mask, read-only. Non-finite coordinates
            — the output of a failed projection — are rejected.
        """

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
    """A group of symmetry-equivalent reflections in a simulation.

    Purpose
    -------
    Symmetry-equivalent reflections carry the same structure factor and
    d-spacing, so they succeed or fail together. Grouping them lets an
    indexing report say *which family* went unmatched, which distinguishes a
    wrong centring or wrong phase from a wrong orientation.

    Attributes
    ----------
    representative_miller_indices : np.ndarray
        The family's representative ``(hkl)``.
    member_miller_indices : np.ndarray
        ``(m, 3)`` all members of the family.
    representative_spot_index : int
    spot_indices : np.ndarray
        Indices of the simulated spots belonging to this family.
    multiplicity : int
        Family size; strictly positive.
    total_intensity : float
        Summed intensity over the family's spots.
    """

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
    """One simulated reflection, with its geometry and intensity.

    Attributes
    ----------
    miller_indices : np.ndarray
        The reflection ``(hkl)``.
    reciprocal_vector_lab : np.ndarray
        ``g`` in the laboratory frame, in inverse angstroms.
    outgoing_direction_lab : np.ndarray
        Unit diffracted-beam direction.
    detector_coordinates_px : np.ndarray
        Where the beam meets the detector.
    excitation_error_inv_angstrom : float
        Distance from the Ewald sphere. Zero is exact Bragg; the further from
        zero, the weaker the reflection.
    intensity : float
        Kinematic intensity, relative and not absolute.
    two_theta_rad, azimuth_rad : float
        Detector polar coordinates of the spot.
    on_detector : bool
        Whether the beam intersects the detector plane at all.
    accepted_by_mask : bool
        Whether it also falls inside the acceptance mask.
    family_id : int, optional
        Index into the simulation's reflection families.
    """

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
    """A simulated kinematic diffraction pattern with its full context.

    Purpose
    -------
    The forward model that pattern indexing inverts: given geometry, phase,
    and orientation, the reflections that appear and where. It also carries
    the matching and orientation-search entry points, so simulation and
    comparison stay bound to the same conventions.

    Limits
    ------
    Kinematic (single-scattering) throughout. Multiple scattering, dynamical
    intensity transfer, and absorption are not modelled, so intensities are
    indicative; relative intensities within one zone are the meaningful
    output.

    Attributes
    ----------
    geometry : DiffractionGeometry
    phase : Phase
    spots : tuple of KinematicSpot
        Every simulated reflection, including those off the detector. Use
        :meth:`accepted_spots` when comparing against a measurement.
    reflection_families : tuple of ReflectionFamily
    orientation : Orientation, optional
    zone_axis : ZoneAxis, optional
    provenance : ProvenanceRecord, optional
    """

    geometry: DiffractionGeometry
    phase: Phase
    spots: tuple[KinematicSpot, ...]
    reflection_families: tuple[ReflectionFamily, ...] = ()
    orientation: Orientation | None = None
    zone_axis: ZoneAxis | None = None
    provenance: ProvenanceRecord | None = None

    def accepted_spots(self) -> tuple[KinematicSpot, ...]:
        """The simulated spots that fell inside the detector acceptance mask.

        Use these, not :attr:`spots`, when comparing against a measurement: the
        rejected ones were never observable.
        """

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
        foil_thickness_angstrom: float | None = None,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        max_distance_px: float = 10.0,
        cluster_radius_px: float = 5.0,
        use_only_accepted: bool = True,
    ) -> tuple[OrientationIndexingCandidate, ...]:
        """Score a list of candidate orientations against an observed pattern.

        Purpose
        -------
        Orientation determination by exhaustive comparison: simulate the pattern
        for each candidate, match it to the observation, and rank by agreement.
        This is the geometric ancestor of dictionary indexing.

        Parameters
        ----------
        geometry : DiffractionGeometry
            Beam, specimen, and detector geometry used for every simulation.
        phase : Phase
            Crystal structure and symmetry used to construct reciprocal vectors.
        miller_indices : np.ndarray
            Candidate reflection indices passed to :meth:`simulate_spots`.
        pattern : DiffractionPattern
            The observation to match against.
        candidate_orientations : list of Orientation or OrientationSet
            The search space. Its resolution bounds the achievable angular
            accuracy; refine the winner with
            :meth:`refine_orientation_candidate`.
        zone_axis : ZoneAxis or None
            Optional zone-axis constraint passed to :meth:`simulate_spots`.
        max_excitation_error_inv_angstrom : float
            Maximum excitation error retained by :meth:`simulate_spots`.
        intensity_model : str
            Reflection intensity model used by :meth:`simulate_spots`.
        excitation_sigma_inv_angstrom : float
            Excitation-error width used by the kinematic intensity proxy.
        foil_thickness_angstrom : float or None
            Optional finite-thickness shape-factor input.
        acceptance_mask : DetectorAcceptanceMask or None
            Optional detector region in which simulated spots are observable.
        max_distance_px : float
            Largest observed-to-simulated spot separation accepted as a match.
        cluster_radius_px : float
            Radius used to cluster nearby observed spots before association.
        use_only_accepted : bool
            If true, compare only observations accepted by the detector mask.

        Returns
        -------
        tuple of OrientationIndexingCandidate
            Sorted best first by :attr:`IndexingCandidate.score`.
        """

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
                foil_thickness_angstrom=foil_thickness_angstrom,
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
        foil_thickness_angstrom: float | None = None,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        max_distance_px: float = 10.0,
        cluster_radius_px: float = 5.0,
        use_only_accepted: bool = True,
        search_half_width_deg: float = 2.0,
        step_deg: float = 1.0,
        iterations: int = 2,
    ) -> OrientationRefinementResult:
        """Locally refine an orientation estimate against an observed pattern.

        Purpose
        -------
        Improve an orientation beyond the resolution of the coarse candidate
        list, by searching a shrinking neighbourhood around the current best.

        Method
        ------
        Iterative local grid search: at each iteration, orientations within
        ``search_half_width_deg`` are sampled at ``step_deg``, scored as in
        :meth:`rank_orientation_candidates`, and the best becomes the next
        centre with a reduced step. This is a local search — it improves a good
        estimate and cannot recover from a wrong one.

        Parameters
        ----------
        initial_orientation : Orientation
            Starting estimate, typically the winner from
            :meth:`rank_orientation_candidates`.
        search_half_width_deg : float
            Half-width of the first search neighbourhood.
        step_deg : float
            Initial angular step.
        iterations : int
            Number of refinement rounds.
        Other parameters : see :meth:`simulate_spots` and
            :meth:`associate_to_pattern`.

        Returns
        -------
        OrientationRefinementResult
            The refined orientation and the score history, so convergence can be
            inspected rather than assumed.
        """

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
            foil_thickness_angstrom=foil_thickness_angstrom,
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
                            foil_thickness_angstrom=foil_thickness_angstrom,
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
        """Match this simulation's spots against an observed pattern.

        Purpose
        -------
        The assignment step of indexing: pair each observed spot cluster with
        the nearest unused simulated reflection, and report what was left over
        on both sides.

        Parameters
        ----------
        pattern : DiffractionPattern
            The observation. Its geometry and phase must equal this
            simulation's; a mismatch raises rather than producing a meaningless
            match.
        max_distance_px : float
            Largest accepted separation between an observation and its assigned
            reflection.
        cluster_radius_px : float
            Radius used to merge raw detections into spot clusters first.
        use_only_accepted : bool
            Match only against detector-accepted spots (default).

        Returns
        -------
        IndexingCandidate
            Matches, unmatched observed clusters, and unmatched simulated spots
            — the unmatched sets are as diagnostic as the matches.
        """

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
        foil_thickness_angstrom: float | None = None,
        acceptance_mask: DetectorAcceptanceMask | None = None,
        deduplicate_families: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> KinematicSimulation:
        """Simulate a kinematic diffraction pattern on the detector.

        Purpose
        -------
        Given a geometry, a phase, and a set of reflections, compute where each
        reflection lands and how strong it is — the forward model that pattern
        indexing inverts.

        Method and limits
        -----------------
        Kinematic (single-scattering) theory throughout: intensities are
        ``|F|^2`` modulated by an excitation-error term. Multiple scattering,
        dynamical intensity transfer, and absorption are not modelled, so
        absolute intensities are indicative and relative intensities within a
        zone are the meaningful output.

        Parameters
        ----------
        geometry : DiffractionGeometry
            Detector, beam, and frame definitions.
        phase : Phase
            Supplies the lattice and scattering factors.
        miller_indices : np.ndarray
            ``(n, 3)`` integer-valued reflections to consider.
        orientation : Orientation, optional
            Crystal orientation. Must be consistent with ``phase`` when its
            phase is set.
        zone_axis : ZoneAxis, optional
            Zone axis to align to, for zone-axis pattern simulation.
        max_excitation_error_inv_angstrom : float
            Reflections further than this from the Ewald sphere are dropped.
            This is the relrod cut-off; widening it admits weaker, more
            off-axis reflections.
        intensity_model : str
            Intensity weighting model, ``"kinematic_proxy"`` by default.
        excitation_sigma_inv_angstrom : float
            Width of the legacy Lorentzian excitation-error falloff. Used only
            when ``foil_thickness_angstrom`` is ``None``.
        foil_thickness_angstrom : float, optional
            Positive plane-parallel foil thickness. When supplied, the exact
            normalized ``sinc^2(t s_g)`` shape factor replaces the legacy
            Lorentzian excitation-error proxy.
        acceptance_mask : DetectorAcceptanceMask, optional
            Marks which spots would be observable.
        deduplicate_families : bool
            Collapse symmetry-equivalent reflections to one representative each.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        KinematicSimulation
            Spots with Miller indices, detector coordinates, intensities, and
            acceptance flags, plus reflection-family grouping.
        """

        miller_array = as_float_array(miller_indices, shape=(None, 3))
        if max_excitation_error_inv_angstrom < 0.0:
            raise ValueError("max_excitation_error_inv_angstrom must be non-negative.")
        if excitation_sigma_inv_angstrom <= 0.0:
            raise ValueError("excitation_sigma_inv_angstrom must be strictly positive.")
        if foil_thickness_angstrom is not None and (
            not np.isfinite(foil_thickness_angstrom) or foil_thickness_angstrom <= 0.0
        ):
            raise ValueError(
                "foil_thickness_angstrom must be finite and strictly positive when set."
            )
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
        # Vectorised over all candidate reflections. Only zone-axis
        # orthogonality and excitation error remove reflections; on-detector /
        # acceptance are stored fields, not filters. Geometry is batched; the
        # trivial per-reflection intensity and family key are kept in the
        # survivor assembly loop for an exact match with the scalar path.
        candidate_spots: list[_CandidateSpot] = []
        miller_int_all = np.rint(miller_array).astype(np.int64)
        reciprocal_basis = phase.lattice.reciprocal_basis().matrix
        reciprocal_crystal = miller_int_all.astype(np.float64) @ reciprocal_basis.T
        if orientation is not None:
            reciprocal_specimen = reciprocal_crystal @ orientation.rotation.as_matrix().T
        else:
            reciprocal_specimen = reciprocal_crystal
        if zone_axis_vector is not None:
            zone_dot = reciprocal_specimen @ zone_axis_vector
            zone_scale = np.maximum(1.0, np.linalg.norm(reciprocal_specimen, axis=1))
            zone_ok = np.abs(zone_dot) <= _ZONE_AXIS_ORTHOGONALITY_ATOL * zone_scale
        else:
            zone_ok = np.ones(miller_int_all.shape[0], dtype=bool)
        reciprocal_lab_all = geometry.specimen_vectors_to_lab(reciprocal_specimen)
        outgoing_wavevectors = incident[None, :] + reciprocal_lab_all
        excitation_all = np.linalg.norm(outgoing_wavevectors, axis=1) - incident_magnitude
        excitation_ok = np.abs(excitation_all) <= max_excitation_error_inv_angstrom
        survivors = np.flatnonzero(zone_ok & excitation_ok)
        if survivors.size:
            outgoing_directions = normalize_vectors(outgoing_wavevectors[survivors])
            coordinates_px_all, valid_all = geometry.project_directions_to_detector_px(
                outgoing_directions
            )
            basis = geometry.detector_basis_lab
            two_theta_all = np.arccos(
                np.clip(outgoing_directions @ geometry.beam_direction_lab, -1.0, 1.0)
            )
            azimuth_all = np.arctan2(
                outgoing_directions @ basis[:, 1], outgoing_directions @ basis[:, 0]
            )
            finite_all = np.all(np.isfinite(coordinates_px_all), axis=1)
            on_detector_all = (
                valid_all
                & finite_all
                & (coordinates_px_all[:, 0] >= 0.0)
                & (coordinates_px_all[:, 0] <= geometry.detector_shape[1] - 1)
                & (coordinates_px_all[:, 1] >= 0.0)
                & (coordinates_px_all[:, 1] <= geometry.detector_shape[0] - 1)
            )
            if acceptance_mask is None:
                accepted_all = on_detector_all
            else:
                accepted_all = (
                    np.asarray(acceptance_mask.contains(geometry, coordinates_px_all), dtype=bool)
                    & on_detector_all
                )
            reciprocal_lab_survivors = reciprocal_lab_all[survivors]
            for local, global_index in enumerate(survivors):
                miller_triplet_int = miller_int_all[global_index]
                excitation_error = float(excitation_all[global_index])
                candidate_spots.append(
                    {
                        "miller_indices": miller_triplet_int,
                        "reciprocal_vector_lab": reciprocal_lab_survivors[local],
                        "outgoing_direction_lab": outgoing_directions[local],
                        "detector_coordinates_px": coordinates_px_all[local],
                        "excitation_error_inv_angstrom": excitation_error,
                        "intensity": _kinematic_intensity(
                            reciprocal_lab_survivors[local],
                            excitation_error,
                            model=intensity_model,
                            excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
                            foil_thickness_angstrom=foil_thickness_angstrom,
                        ),
                        "two_theta_rad": float(two_theta_all[local]),
                        "azimuth_rad": float(azimuth_all[local]),
                        "on_detector": bool(on_detector_all[local]),
                        "accepted_by_mask": bool(accepted_all[local]),
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
    """A measured diffraction pattern: spot positions and intensities.

    Purpose
    -------
    The observation side of indexing. Because it carries its geometry and
    phase, the detector coordinates can be converted to scattering vectors
    and angles without the caller restating the experiment.

    Attributes
    ----------
    coordinates_px : np.ndarray
        ``(n, 2)`` detected positions in detector pixels.
    intensities : np.ndarray
        ``(n,)`` intensities; must be finite and non-negative.
    geometry : DiffractionGeometry
    phase : Phase
    orientation : Orientation, optional
        When known; its frames and phase are checked against the geometry
        and phase rather than assumed compatible.
    provenance : ProvenanceRecord, optional
    """

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
        """Detector-plane offsets of the observed spots, in millimetres."""

        return self.geometry.detector_coordinates_mm(self.coordinates_px)

    def outgoing_directions_lab(self) -> np.ndarray:
        """Unit scattered-beam directions of the observed spots, in the lab frame."""

        return self.geometry.outgoing_directions_lab(self.coordinates_px)

    def scattering_vectors_lab(self) -> np.ndarray:
        """Scattering vectors ``q`` of the observed spots, in inverse angstroms.

        These are the measured reciprocal-space positions that indexing compares
        against computed reciprocal-lattice vectors.
        """

        return self.geometry.scattering_vectors_lab(self.coordinates_px)

    def two_theta_rad(self) -> np.ndarray:
        """Scattering angles ``2*theta`` of the observed spots, in radians."""

        return self.geometry.two_theta_rad(self.coordinates_px)

    def azimuth_rad(self) -> np.ndarray:
        """Azimuthal angles of the observed spots about the beam, in radians."""

        return self.geometry.azimuth_rad(self.coordinates_px)

    def cluster_observations(
        self, *, max_distance_px: float = 5.0
    ) -> tuple[DetectedSpotCluster, ...]:
        """Group nearby detections into spot clusters.

        Purpose
        -------
        A measured spot usually appears as several detections. Indexing needs
        one position per physical reflection, so detections within
        ``max_distance_px`` of one another are merged by single-linkage
        clustering and reduced to an intensity-weighted centroid — sub-pixel
        accurate, and robust to asymmetric spot shapes.

        Parameters
        ----------
        max_distance_px : float
            Linkage radius. Too large merges genuinely distinct reflections;
            too small splits one spot into several.

        Returns
        -------
        tuple of DetectedSpotCluster
            Each carrying its member detection indices, weighted centre, and
            total intensity. Clusters of zero total intensity fall back to the
            unweighted centroid.
        """

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
    """A group of detections merged into one observed reflection.

    Purpose
    -------
    A measured spot usually appears as several detections. Indexing needs one
    position per physical reflection, so detections are clustered and reduced
    to an intensity-weighted centroid — sub-pixel accurate and robust to
    asymmetric spot shapes.

    Attributes
    ----------
    cluster_id : int
    member_indices : np.ndarray
        Indices of the detections merged; must be non-empty.
    center_px : np.ndarray
        Intensity-weighted centre.
    total_intensity : float
        Summed intensity; finite and non-negative.
    """

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
    """A match between one observed spot cluster and one simulated reflection.

    Attributes
    ----------
    observed_cluster_id : int
    simulated_spot_index : int
    residual_px : float
        Detector-plane distance between the two; finite and non-negative.
        This is the per-spot accuracy of the indexing solution.
    family_id : int, optional
        The reflection family the matched spot belongs to.
    """

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
    """One indexing solution, with everything needed to judge it.

    Purpose
    -------
    Reports not only the matches but both unmatched sets — observed clusters
    with no reflection, and predicted reflections with no observation. Those
    are as diagnostic as the matches: unexplained observations suggest a
    wrong phase or a second grain, while unobserved predictions suggest wrong
    centring or an intensity cut-off set too low.

    Attributes
    ----------
    pattern : DiffractionPattern
        The observation.
    simulation : KinematicSimulation
        The prediction being tested.
    observation_clusters : tuple of DetectedSpotCluster
    matches : tuple of SpotAssignment
    unmatched_observed_cluster_ids : np.ndarray
    unmatched_simulated_spot_indices : np.ndarray
    """

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
        """Fraction of observed spot clusters that received an assignment.

        In ``[0, 1]``; ``0.0`` when there were no observations. The recall side
        of indexing quality.
        """

        if not self.observation_clusters:
            return 0.0
        return float(len(self.matches) / len(self.observation_clusters))

    @property
    def mean_residual_px(self) -> float:
        """Mean detector-plane residual of the assigned spots, in pixels.

        ``inf`` when nothing matched, so an empty solution never wins a
        smallest-residual comparison. The precision side of indexing quality.
        """

        if not self.matches:
            return np.inf
        return float(np.mean([match.residual_px for match in self.matches]))

    @property
    def score(self) -> float:
        """Combined indexing quality in ``[0, 1]``, larger is better.

        ``match_fraction / (1 + mean_residual_px)``: a solution must both explain
        many observed spots and place them accurately. A high match fraction
        with poor residuals, or tight residuals on a handful of spots, both
        score low. This is a heuristic ranking score, not a likelihood.
        """

        if not self.matches:
            return 0.0
        residual_penalty = 1.0 / (1.0 + self.mean_residual_px)
        return float(self.match_fraction * residual_penalty)

    def family_reports(self) -> tuple[FamilyIndexingReport, ...]:
        """Per-reflection-family breakdown of this indexing solution.

        Purpose
        -------
        An overall score hides *which* reflections were explained. Reporting per
        symmetry family exposes the diagnostic pattern — for example a family
        systematically unmatched, which points to a wrong centring assumption or
        a wrong phase rather than to a wrong orientation.

        Returns
        -------
        tuple of FamilyIndexingReport
            One entry per simulated family, with matched counts, matched
            fraction, total family intensity, and mean residual. Empty when the
            simulation carried no family grouping.
        """

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
    """One candidate orientation together with its indexing solution.

    The unit that :meth:`KinematicSimulation.rank_orientation_candidates`
    sorts, pairing the trial orientation with the quality of the pattern
    match it produced.

    Attributes
    ----------
    orientation_index : int
        Position in the candidate list.
    orientation : Orientation
    indexing : IndexingCandidate
    """

    orientation_index: int
    orientation: Orientation
    indexing: IndexingCandidate

    @property
    def score(self) -> float:
        """Indexing quality of this orientation candidate; see
        :attr:`IndexingCandidate.score`.
        """

        return self.indexing.score


@dataclass(frozen=True, slots=True)
class FamilyIndexingReport:
    """How one reflection family fared in an indexing solution.

    Purpose
    -------
    The per-family breakdown that an overall score hides. A family
    systematically unmatched — high simulated count, zero matches — points to
    a wrong centring assumption or a wrong phase, not to a wrong orientation.

    Attributes
    ----------
    family_id : int
    representative_miller_indices : np.ndarray
    multiplicity : int
    simulated_spot_count : int
    matched_spot_count : int
    matched_fraction : float
    total_family_intensity : float
        Strong families that go unmatched are the informative ones.
    mean_residual_px : float
        ``inf`` when nothing in the family matched.
    observed_cluster_ids : np.ndarray
    """

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
    """The outcome of a local orientation refinement, with its search history.

    Purpose
    -------
    Records not just the refined orientation but how it was reached — how
    many candidates were evaluated, over how many iterations, and how the
    search window shrank — so convergence can be inspected rather than
    assumed. Refinement is a local search: it improves a good estimate and
    cannot recover from a wrong one.

    Attributes
    ----------
    seed_orientation : Orientation
        The starting estimate.
    refined_candidate : OrientationIndexingCandidate
        The best solution found.
    evaluated_candidates : int
    iterations : int
    initial_search_half_width_deg : float
    final_step_deg : float
        The final angular step, which bounds the achieved precision.
    """

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


def _candidate_zone_axes(max_index: int) -> np.ndarray:
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    rows = [
        (u, v, w)
        for u in range(-max_index, max_index + 1)
        for v in range(-max_index, max_index + 1)
        for w in range(-max_index, max_index + 1)
        if not (u == 0 and v == 0 and w == 0)
    ]
    values = np.asarray(rows, dtype=np.int64)
    gcds = np.gcd.reduce(np.abs(values), axis=1)
    values = values // gcds[:, None]
    signs = np.sign(values[np.arange(values.shape[0]), np.argmax(values != 0, axis=1)])
    values = values * signs[:, None]
    values = np.unique(values, axis=0)
    return as_int_array(values, shape=(values.shape[0], 3))


def estimate_zone_axis(
    pattern: DiffractionPattern,
    *,
    max_index: int = 4,
    min_intensity_fraction: float = 0.0,
) -> ZoneAxis:
    """Estimate the zone axis of a measured diffraction pattern.

    Purpose
    -------
    Recover ``[uvw]`` from spot positions alone, without an orientation
    estimate — the first step in indexing an unknown pattern.

    Method
    ------
    All reflections of a zone lie in the plane perpendicular to the zone
    axis, so their scattering vectors are perpendicular to it. Every
    candidate direction up to ``max_index`` is scored by the mean absolute
    normalized projection of the measured scattering vectors onto it, and
    the smallest-residual candidate wins. Enumeration is exhaustive, so the
    result is deterministic and independent of any starting guess.

    Parameters
    ----------
    pattern : DiffractionPattern
        The observation. Must contain at least one spot.
    max_index : int
        Largest absolute index of the enumerated candidates. Raising it
        admits higher-order zone axes at cubic cost.
    min_intensity_fraction : float
        Discard spots weaker than this fraction of the strongest, in
        ``[0, 1]``. Use it to keep noise out of the fit.

    Returns
    -------
    ZoneAxis
        The best candidate, on the pattern's phase.

    Notes
    -----
    Reports the best candidate unconditionally; it does not decide whether
    the pattern really is a zone-axis pattern. Check the resulting indexing
    residuals before trusting the axis.
    """

    if not 0.0 <= min_intensity_fraction <= 1.0:
        raise ValueError("min_intensity_fraction must lie in [0, 1].")
    scattering_vectors = pattern.geometry.lab_vectors_to_specimen(pattern.scattering_vectors_lab())
    intensities = pattern.intensities
    if intensities.size == 0:
        raise ValueError("pattern must contain at least one spot.")
    if min_intensity_fraction > 0.0:
        threshold = float(np.max(intensities)) * min_intensity_fraction
        mask = intensities >= threshold
        if not np.any(mask):
            raise ValueError("No pattern spots remain after intensity filtering.")
        scattering_vectors = scattering_vectors[mask]
    norms = np.linalg.norm(scattering_vectors, axis=1)
    valid = norms > _INTENSITY_EPSILON
    if not np.any(valid):
        raise ValueError("pattern scattering vectors must contain a non-zero vector.")
    scattering_vectors = scattering_vectors[valid]
    norms = norms[valid]
    direct_basis = pattern.phase.lattice.direct_basis().matrix
    candidates = _candidate_zone_axes(max_index)
    best_axis: np.ndarray | None = None
    best_score = np.inf
    for candidate in candidates:
        direction = normalize_vector(direct_basis @ candidate.astype(np.float64))
        residuals = np.abs(scattering_vectors @ direction) / np.maximum(norms, 1.0)
        score = float(np.mean(residuals))
        if score < best_score:
            best_score = score
            best_axis = candidate
    if best_axis is None:
        raise ValueError("No zone-axis candidate could be estimated.")
    return ZoneAxis(best_axis, phase=pattern.phase)


def index_saed_pattern(
    pattern: DiffractionPattern,
    miller_indices: np.ndarray,
    *,
    orientation: Orientation | None = None,
    zone_axis: ZoneAxis | None = None,
    max_excitation_error_inv_angstrom: float = 5e-2,
    intensity_model: str = "kinematic_proxy",
    excitation_sigma_inv_angstrom: float = 5e-2,
    foil_thickness_angstrom: float | None = None,
    acceptance_mask: DetectorAcceptanceMask | None = None,
    max_distance_px: float = 10.0,
    cluster_radius_px: float = 5.0,
    use_only_accepted: bool = True,
) -> IndexingCandidate:
    """Index a measured SAED pattern against a phase.

    Purpose
    -------
    The single-call route from an observed spot list to an indexed solution:
    determine the zone axis if it was not supplied, simulate the expected
    pattern, and match it to the observation.

    Parameters
    ----------
    pattern : DiffractionPattern
        The observation.
    miller_indices : np.ndarray
        ``(n, 3)`` candidate reflections to simulate.
    orientation : Orientation, optional
        Fixes the in-plane rotation when known.
    zone_axis : ZoneAxis, optional
        Skips the estimation step; supply it when the axis is known.
    max_excitation_error_inv_angstrom : float
        Excitation-error limit; see :meth:`KinematicSimulation.simulate_spots`.
    intensity_model : str
        Spot-intensity policy; see :meth:`KinematicSimulation.simulate_spots`.
    excitation_sigma_inv_angstrom : float
        Excitation envelope width; see :meth:`KinematicSimulation.simulate_spots`.
    foil_thickness_angstrom : float, optional
        Exact finite-thickness shape-factor input; see
        :meth:`KinematicSimulation.simulate_spots`.
    acceptance_mask : DetectorAcceptanceMask, optional
        Detector acceptance; see :meth:`KinematicSimulation.simulate_spots`.
    max_distance_px : float
        Matching radius; see :meth:`KinematicSimulation.associate_to_pattern`.
    cluster_radius_px : float
        Observed-spot clustering radius; see
        :meth:`KinematicSimulation.associate_to_pattern`.
    use_only_accepted : bool
        Whether to match only accepted simulated spots; see
        :meth:`KinematicSimulation.associate_to_pattern`.

    Returns
    -------
    IndexingCandidate
        The solution, with matched and unmatched spots on both sides and a
        per-family breakdown available through
        :meth:`IndexingCandidate.family_reports`.
    """

    selected_zone_axis = zone_axis if zone_axis is not None else estimate_zone_axis(pattern)
    simulation = KinematicSimulation.simulate_spots(
        pattern.geometry,
        pattern.phase,
        miller_indices,
        orientation=orientation,
        zone_axis=selected_zone_axis,
        max_excitation_error_inv_angstrom=max_excitation_error_inv_angstrom,
        intensity_model=intensity_model,
        excitation_sigma_inv_angstrom=excitation_sigma_inv_angstrom,
        foil_thickness_angstrom=foil_thickness_angstrom,
        acceptance_mask=acceptance_mask,
    )
    return simulation.associate_to_pattern(
        pattern,
        max_distance_px=max_distance_px,
        cluster_radius_px=cluster_radius_px,
        use_only_accepted=use_only_accepted,
    )
