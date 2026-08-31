from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array, normalize_vector
from pytex.core._chemistry import atomic_number
from pytex.core.conventions import FrameDomain
from pytex.core.frame_catalog import detector_frame as catalog_detector_frame
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import MillerIndex, Phase, ReciprocalLatticeVector, ZoneAxis
from pytex.core.notation import format_plane_indices
from pytex.core.provenance import ProvenanceRecord

#: Relative intensity at or below which a zone reflection counts as
#: kinematically forbidden.
#:
#: A forbidden reflection's structure factor is exactly zero in theory and a few
#: parts in :math:`10^{16}` of the strongest reflection in floating point. One
#: part in :math:`10^{4}` is the same threshold
#: :class:`~pytex.diffraction.kinematic.KinematicSimulationConfig` uses for
#: ``min_relative_intensity``, so the two engines agree on which reflections are
#: absent. The gap between a genuinely weak reflection and a forbidden one spans
#: many orders of magnitude, so the exact value is not delicate.
FORBIDDEN_RELATIVE_INTENSITY = 1e-4


def _choose_zone_basis(zone_axis: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(trial, zone_axis))), 1.0, atol=1e-8):
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u_axis = normalize_vector(np.cross(zone_axis, trial))
    v_axis = normalize_vector(np.cross(zone_axis, u_axis))
    basis = np.column_stack([u_axis, v_axis, zone_axis])
    basis = np.ascontiguousarray(basis)
    basis.setflags(write=False)
    return basis


def _structure_factor_electron(phase: Phase, hkl: np.ndarray) -> complex:
    if phase.unit_cell is None or not phase.unit_cell.sites:
        return complex(1.0, 0.0)
    reciprocal = phase.lattice.reciprocal_basis().matrix
    g_cart = reciprocal @ hkl.astype(np.float64)
    g_sq = float(np.dot(g_cart, g_cart))
    total = 0.0j
    for site in phase.unit_cell.sites:
        z = float(atomic_number(site.species))
        occupancy = float(site.occupancy)
        b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
        damping = np.exp(-(b_iso * g_sq) / max(16.0 * np.pi * np.pi, 1e-12))
        phase_factor = np.exp(
            2.0j * np.pi * float(np.dot(hkl.astype(np.float64), site.fractional_coordinates))
        )
        total += occupancy * z * damping * phase_factor
    return complex(total)


def _enumerate_zone_hkls(max_index: int) -> np.ndarray:
    values = range(-max_index, max_index + 1)
    hkls = [
        (h, k, ell)
        for h in values
        for k in values
        for ell in values
        if not (h == 0 and k == 0 and ell == 0)
    ]
    array = np.asarray(hkls, dtype=np.int64)
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class SAEDSpot:
    """One reflection in a simulated selected-area electron diffraction pattern.

    Attributes
    ----------
    miller_indices : np.ndarray
        The reflection ``(hkl)``.
    detector_coordinates : np.ndarray
        Position on the detector in millimetres, at radius proportional to
        ``|g|`` through the camera constant.
    g_magnitude_inv_angstrom : float
    d_spacing_angstrom : float
    intensity : float
        Kinematic and relative, unless the spot is marked as double diffraction,
        in which case it is an indicative observability estimate instead. See
        :attr:`double_diffraction_parents`.
    double_diffraction_parents : np.ndarray, optional
        ``(2, 3)`` integer pair ``(g1, g2)`` summing to :attr:`miller_indices`,
        present exactly for a reflection that is kinematically forbidden and
        appears only because ``include_double_diffraction`` was requested.
        ``None`` for an ordinary kinematic reflection.
    label : str
        Rendered index label, empty when the spot is beyond the label limit.
    """

    miller_indices: np.ndarray
    reciprocal_vector_crystal: np.ndarray
    reciprocal_vector_detector: np.ndarray
    detector_coordinates: np.ndarray
    intensity: float
    excitation_error_inv_angstrom: float
    label: str | None = None
    double_diffraction_parents: np.ndarray | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        object.__setattr__(
            self,
            "reciprocal_vector_crystal",
            as_float_array(self.reciprocal_vector_crystal, shape=(3,)),
        )
        object.__setattr__(
            self,
            "reciprocal_vector_detector",
            as_float_array(self.reciprocal_vector_detector, shape=(3,)),
        )
        object.__setattr__(
            self,
            "detector_coordinates",
            as_float_array(self.detector_coordinates, shape=(2,)),
        )
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("SAEDSpot.intensity must be finite and non-negative.")
        if not np.isfinite(self.excitation_error_inv_angstrom):
            raise ValueError("SAEDSpot.excitation_error_inv_angstrom must be finite.")
        if self.double_diffraction_parents is not None:
            parents = as_int_array(self.double_diffraction_parents, shape=(2, 3))
            if not np.array_equal(parents.sum(axis=0), self.miller_indices):
                raise ValueError(
                    "SAEDSpot.double_diffraction_parents must sum to the reflection they produce."
                )
            object.__setattr__(self, "double_diffraction_parents", parents)

    @property
    def is_double_diffraction(self) -> bool:
        """Whether this reflection is present only through double diffraction.

        ``True`` marks a reflection whose structure factor is (near) zero and
        which is drawn because two excited reflections sum to it. Read it as an
        observability statement, never as a kinematic intensity.
        """

        return self.double_diffraction_parents is not None

    def double_diffraction_origin_label(self) -> str:
        """``"(g1) + (g2)"`` for a marked reflection, or the empty string.

        The path that puts the spot on the plate, written in the repository's
        plane notation so it reads the way the literature writes it.
        """

        if self.double_diffraction_parents is None:
            return ""
        first, second = np.asarray(self.double_diffraction_parents, dtype=int)
        return (
            f"{format_plane_indices(tuple(int(v) for v in first), style='plain')} + "
            f"{format_plane_indices(tuple(int(v) for v in second), style='plain')}"
        )


@dataclass(frozen=True, slots=True)
class SAEDPattern:
    """A simulated zone-axis electron diffraction pattern.

    Purpose
    -------
    The spot pattern seen in a TEM with the beam along a given zone axis,
    scaled by the camera constant so radial position is proportional to
    ``|g|``.

    Limits
    ------
    Kinematic and geometric: relative intensities are indicative rather than
    quantitative, and double diffraction — which can make a formally
    forbidden reflection appear — is not modelled here. For that, use the
    vectorized engine
    :func:`pytex.diffraction.kinematic.simulate_zone_axis_spots` with
    ``KinematicSimulationConfig(include_double_diffraction=True)``, which
    adds those reflections and flags them.

    Attributes
    ----------
    phase : Phase
    zone_axis : ZoneAxis
        The beam direction in crystal indices.
    spots : tuple of SAEDSpot
    camera_constant_mm_angstrom : float
        The product ``L*lambda`` setting the pattern scale.
    provenance : ProvenanceRecord, optional
    """

    phase: Phase
    zone_axis: ZoneAxis
    detector_frame: ReferenceFrame
    reciprocal_frame: ReferenceFrame
    camera_constant_mm_angstrom: float
    spots: tuple[SAEDSpot, ...]
    zone_basis_crystal: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "spots", tuple(self.spots))
        object.__setattr__(
            self, "zone_basis_crystal", as_float_array(self.zone_basis_crystal, shape=(3, 3))
        )
        if self.detector_frame.domain is not FrameDomain.DETECTOR:
            raise ValueError("SAEDPattern.detector_frame must belong to the detector domain.")
        if self.reciprocal_frame.domain is not FrameDomain.RECIPROCAL:
            raise ValueError("SAEDPattern.reciprocal_frame must belong to the reciprocal domain.")
        if self.camera_constant_mm_angstrom <= 0.0:
            raise ValueError("SAEDPattern.camera_constant_mm_angstrom must be positive.")

    def detector_extent_mm(self) -> float:
        """A suggested plot half-extent covering all spots, in millimetres.

        The largest spot radius plus 15 percent margin, so no spot sits on the
        frame edge. Returns ``1.0`` for an empty pattern.
        """

        if not self.spots:
            return 1.0
        radii = [float(np.linalg.norm(spot.detector_coordinates)) for spot in self.spots]
        return max(radii) * 1.15


def _apply_double_diffraction(spots: list[SAEDSpot], *, coupling: float) -> list[SAEDSpot]:
    """Give the reachable forbidden reflections of a zone an indicative intensity.

    A beam diffracted by ``g1`` is itself an incident beam inside the crystal, so
    diffracting it again by ``g2`` sends it out along ``g1 + g2``. The reachable
    set is therefore the pairwise algebraic sums of the excited reflections —
    the shared rule in :func:`pytex.diffraction.kinematic.double_diffraction_sums`,
    reused here rather than restated.

    Every reflection of the zone is already enumerated by the caller, forbidden
    ones included, at (near) zero intensity. So this replaces the forbidden spots
    that the rule can reach with copies carrying an observability intensity and
    the parent pair that produced them, and leaves every other spot untouched. No
    position changes, and an allowed reflection is never re-labelled as double
    diffraction even though the rule also reaches it — its own structure factor
    already explains it.
    """

    from pytex.diffraction.kinematic import double_diffraction_sums

    if not spots:
        return spots
    intensities = np.asarray([spot.intensity for spot in spots], dtype=np.float64)
    peak = float(intensities.max())
    if peak <= 0.0:
        return spots
    relative = intensities / peak

    excited = relative > FORBIDDEN_RELATIVE_INTENSITY
    forbidden = ~excited
    if not bool(excited.any()) or not bool(forbidden.any()):
        return spots

    indices = np.asarray([spot.miller_indices for spot in spots], dtype=np.int64)
    sums, weight, parents = double_diffraction_sums(indices[excited], relative[excited])
    if sums.shape[0] == 0:
        return spots

    # The two-step amplitude scales as F(g1) F(g2), so the intensity scales as
    # the product of the parent intensities; paths are summed because their
    # relative phases are outside what a kinematic treatment can supply. Clipping
    # at the strongest kinematic spot keeps the pattern's intensity scale intact.
    reachable = {
        tuple(int(value) for value in row): (float(total), pair)
        for row, total, pair in zip(sums, weight, parents, strict=True)
    }

    updated = list(spots)
    for row in np.flatnonzero(forbidden):
        spot = spots[int(row)]
        match = reachable.get(tuple(int(value) for value in spot.miller_indices))
        if match is None:
            continue
        total, pair = match
        boosted = min(1.0, coupling * total) * peak
        if boosted <= spot.intensity:
            continue
        updated[int(row)] = SAEDSpot(
            miller_indices=spot.miller_indices,
            reciprocal_vector_crystal=spot.reciprocal_vector_crystal,
            reciprocal_vector_detector=spot.reciprocal_vector_detector,
            detector_coordinates=spot.detector_coordinates,
            intensity=boosted,
            excitation_error_inv_angstrom=spot.excitation_error_inv_angstrom,
            label=spot.label,
            double_diffraction_parents=pair,
        )
    return updated


def generate_saed_pattern(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    camera_constant_mm_angstrom: float = 180.0,
    max_index: int = 6,
    max_g_inv_angstrom: float | None = None,
    zone_tolerance_inv_angstrom: float = 1e-6,
    intensity_model: Literal["electron_atomic_number", "unit"] = "electron_atomic_number",
    label_limit: int = 20,
    include_double_diffraction: bool = False,
    double_diffraction_coupling: float = 0.05,
    provenance: ProvenanceRecord | None = None,
) -> SAEDPattern:
    """Simulate a selected-area electron diffraction pattern down a zone axis.

    Purpose
    -------
    The spot pattern seen in a TEM with the beam along ``[uvw]``: which
    reflections of the zone appear, where they land on the detector, and how
    strong they are.

    Method and limits
    -----------------
    Kinematic and geometric. Reflections belonging to the zone are placed by
    the camera constant ``L*lambda``, so radial position is proportional to
    ``|g|``. Intensities come from electron structure factors with no
    dynamical scattering, so relative intensities within a zone are
    indicative rather than quantitative.

    Double diffraction — a diffracted beam re-diffracting inside the crystal,
    which puts a spot at ``g1 + g2`` and can therefore make a formally forbidden
    reflection appear — is off by default and enabled with
    ``include_double_diffraction``. The selection rule is the shared one,
    :func:`pytex.diffraction.kinematic.double_diffraction_sums`. Because this
    function already enumerates the whole ``hkl`` cube of the zone and keeps
    forbidden reflections at (near) zero intensity, enabling the option
    **re-weights and marks reflections that are already present** rather than
    appending rows, which is how the vectorized engine
    :func:`pytex.diffraction.kinematic.simulate_zone_axis_spots` expresses the
    same physics. No spot moves: only intensity and the marking change.

    A centring absence is never revived. Lattice centring conditions define a
    *sublattice* of reciprocal space, a sublattice is closed under addition, and
    so no sum of two excited reflections can land on one. Only a basis absence —
    from a glide plane, a screw axis, or the motif — can be revived, which is
    the physically correct outcome and the reason the hcp ``(0001)`` reflection
    appears on a real plate.

    Parameters
    ----------
    phase : Phase
        Must be the same phase the zone axis is defined on.
    zone_axis : ZoneAxis
        The beam direction in crystal indices.
    camera_constant_mm_angstrom : float
        The product ``L*lambda`` setting the pattern scale; the quantity a
        real microscope is calibrated for.
    max_index : int
        Largest absolute Miller index enumerated.
    max_g_inv_angstrom : float, optional
        Radial cut-off in reciprocal space, equivalent to the recorded
        detector extent.
    zone_tolerance_inv_angstrom : float
        Tolerance on the zone-law condition, standing in for the finite
        thickness of the Ewald-sphere intersection.
    intensity_model : str
        ``"electron_atomic_number"`` (default) or ``"unit"``.
    label_limit : int
        Maximum number of spots to label, so a dense pattern stays legible.
    include_double_diffraction : bool
        Whether to give kinematically forbidden reflections reachable as
        ``g1 + g2`` an indicative intensity and mark them. Default ``False``,
        which leaves the pattern purely kinematic.
    double_diffraction_coupling : float
        Scale applied to the two-step weight, in ``(0, 1]``. It carries
        everything kinematic theory cannot supply — beam coupling strength and
        specimen thickness — so a marked reflection's intensity is an
        observability estimate, not a measurement. Ignored when
        ``include_double_diffraction`` is ``False``.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    SAEDPattern
        Spots with indices, detector coordinates, and intensities.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if camera_constant_mm_angstrom <= 0.0:
        raise ValueError("camera_constant_mm_angstrom must be strictly positive.")
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if max_g_inv_angstrom is not None and max_g_inv_angstrom <= 0.0:
        raise ValueError("max_g_inv_angstrom must be strictly positive when provided.")
    if zone_tolerance_inv_angstrom < 0.0:
        raise ValueError("zone_tolerance_inv_angstrom must be non-negative.")
    if not np.isfinite(double_diffraction_coupling) or not 0.0 < double_diffraction_coupling <= 1.0:
        raise ValueError("double_diffraction_coupling must lie in the interval (0, 1].")

    zone_vector = zone_axis.unit_vector
    zone_basis = _choose_zone_basis(zone_vector)
    reciprocal_basis = phase.lattice.reciprocal_basis()
    reciprocal_frame = reciprocal_basis.frame
    detector_frame = catalog_detector_frame(
        f"{phase.name}_saed_detector",
        description="Detector plane for kinematic SAED plotting.",
        provenance=provenance,
    )
    spots: list[SAEDSpot] = []
    for hkl in _enumerate_zone_hkls(max_index):
        reciprocal_vector = ReciprocalLatticeVector.from_miller_index(
            MillerIndex(hkl, phase=phase)
        ).cartesian_vector
        zone_projection = float(np.dot(reciprocal_vector, zone_vector))
        if abs(zone_projection) > zone_tolerance_inv_angstrom:
            continue
        g_magnitude = float(np.linalg.norm(reciprocal_vector))
        if np.isclose(g_magnitude, 0.0):
            continue
        if max_g_inv_angstrom is not None and g_magnitude > max_g_inv_angstrom:
            continue
        detector_coords_mm = camera_constant_mm_angstrom * np.array(
            [
                float(np.dot(reciprocal_vector, zone_basis[:, 0])),
                float(np.dot(reciprocal_vector, zone_basis[:, 1])),
            ],
            dtype=np.float64,
        )
        intensity = (
            1.0
            if intensity_model == "unit"
            else float(abs(_structure_factor_electron(phase, hkl)) ** 2)
            / (1.0 + g_magnitude * g_magnitude)
        )
        spots.append(
            SAEDSpot(
                miller_indices=hkl,
                reciprocal_vector_crystal=reciprocal_vector,
                reciprocal_vector_detector=zone_basis.T @ reciprocal_vector,
                detector_coordinates=detector_coords_mm,
                intensity=intensity,
                excitation_error_inv_angstrom=zone_projection,
                # The repository notation standard renders a specific plane as
                # (hkl) with overbarred negatives. The old space-joined digits
                # were both off-standard and wide enough that adjacent spot
                # labels ran together in a dense zone.
                label=format_plane_indices(tuple(int(value) for value in hkl), style="mathtext"),
            )
        )
    if include_double_diffraction:
        spots = _apply_double_diffraction(spots, coupling=float(double_diffraction_coupling))
    # A *total* order. Intensity and radius alone leave the symmetry-equivalent
    # reflections of one ring tied, and a tie is broken by whatever order they
    # were generated in -- which is not guaranteed to be the same on another
    # machine. The tracked test pattern's sidecar listed its spots in a
    # different order on Linux than on Windows for exactly that reason, and the
    # label limit below is applied by position, so the tie decided which spots
    # got named. The indices settle it.
    spots.sort(
        key=lambda spot: (
            -spot.intensity,
            float(np.linalg.norm(spot.detector_coordinates)),
            tuple(int(value) for value in spot.miller_indices),
        )
    )
    limited_spots = []
    for index, spot in enumerate(spots):
        if index >= label_limit:
            limited_spots.append(
                SAEDSpot(
                    miller_indices=spot.miller_indices,
                    reciprocal_vector_crystal=spot.reciprocal_vector_crystal,
                    reciprocal_vector_detector=spot.reciprocal_vector_detector,
                    detector_coordinates=spot.detector_coordinates,
                    intensity=spot.intensity,
                    excitation_error_inv_angstrom=spot.excitation_error_inv_angstrom,
                    label=None,
                    double_diffraction_parents=spot.double_diffraction_parents,
                )
            )
        else:
            limited_spots.append(spot)
    return SAEDPattern(
        phase=phase,
        zone_axis=zone_axis,
        detector_frame=detector_frame,
        reciprocal_frame=reciprocal_frame,
        camera_constant_mm_angstrom=camera_constant_mm_angstrom,
        spots=tuple(limited_spots),
        zone_basis_crystal=zone_basis,
        provenance=provenance,
    )
