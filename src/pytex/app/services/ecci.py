"""The ECCI workflow: from a measured EBSD orientation to a two-beam tilt.

Electron channelling contrast imaging (ECCI) images dislocations by tilting the
specimen so that one strong reflection is close to the Bragg condition while its
neighbours are suppressed — a two-beam condition. The tilt an operator needs is
most conveniently found starting from an orientation an EBSD system has already
measured: index the point, choose the crystallographic direction you want on the
beam, and solve for the stage move.

This module answers three questions from one EBSD-measured orientation:

1. What does the EBSD camera see right now? Reused wholesale from
   :func:`pytex.diffraction.kikuchi.simulate_kikuchi_pattern`, exactly as
   :mod:`pytex.app.services.ebsd_pattern` computes it.
2. What would an **on-axis** detector see — a BSE detector on the beam axis,
   looking back along it the way a TEM screen does in SAED mode? This is a
   zone-axis pattern of whichever low-index direction the beam is currently
   nearest to, so it is built with the same on-axis machinery TEM patterns use:
   :func:`pytex.tem.synthetic.synthesize_saed_image` down the nearest zone axis,
   with the excitation error of every reflection recomputed against the
   *actual*, continuous beam direction rather than the nominal zone axis — see
   :func:`_on_axis_pattern`.
3. What stage tilt and rotation bring a chosen crystal direction onto the beam?

Why this does not build a second on-axis ``DiffractionGeometry``
------------------------------------------------------------------
:class:`~pytex.diffraction.models.DiffractionGeometry` earns its keep when a
detector sits off the beam axis at a stated elevation and azimuth — that is what
an EBSD screen is. An on-axis detector has no such freedom: by definition its
normal *is* the beam direction, whatever the specimen tilt. That is exactly the
geometry :func:`~pytex.diffraction.saed.generate_saed_pattern` (via
:func:`~pytex.tem.synthetic.synthesize_saed_image`) already assumes for a
zone-axis pattern, so reusing it is simpler than constructing a
``DiffractionGeometry`` at ``detector_elevation_deg`` near 90 and would still
have to special-case the excitation-error computation, which
``DiffractionGeometry`` does not provide.

Why the stage solver is new code, not a reuse of `pytex.tem.navigation`
-------------------------------------------------------------------------
:func:`pytex.tem.navigation.plan_tilt_to_zone_axis` solves a **double-tilt TEM
holder**: two independent tilts, alpha about a fixed axis and beta about a
second axis fixed in the tilted frame. An SEM/ECCI stage is not that stage: it
has one mechanical tilt about a fixed laboratory axis (the same axis
:class:`DiffractionGeometry.for_ebsd` calls the tilt axis) and a rotation about
the specimen's own normal, applied *before* the tilt, exactly as a eucentric SEM
stage is driven. The closed-form geometry is therefore different — derived in
:func:`_stage_branches` — though it is built the same way
``pytex.tem.navigation`` builds its own: an exact closed form, forward-validated
by re-deriving the achieved alignment through the same matrices used to search
for it, never by trusting the algebra that produced it.
"""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import direction_label, phase_parameter, plane_label
from pytex.app.services.traces import clipped_runs

__all__: tuple[str, ...] = ()

_CITATION_ECCI = (
    "Picard, Twigg, Caldwell, Eddy, Neudeck & Trexler, 'Electron channeling contrast imaging "
    "of atomic scale defects', in Springer Handbook of Microscopy (2019), and Crimp, "
    "Microsc. Microanal. 12 (2006) 102, doi:10.1017/S1431927606060438."
)
_CITATION_WILLIAMS = (
    "Williams & Carter, Transmission Electron Microscopy, 2nd ed., chapters 12, 16."
)
_CITATION_SCHWARTZ = (
    "Schwartz, Kumar, Adams & Field (eds.), Electron Backscatter Diffraction in Materials "
    "Science, 2nd ed., Springer 2009, Chs. 1-5."
)

#: Points sampled along each Kikuchi trace before clipping, as in ebsd_pattern.py.
_TRACE_SAMPLES = 721

#: The tilt range DiffractionGeometry.for_ebsd accepts. A solved branch landing
#: at or past this is not a move the geometry (or, in practice, a stage rated to
#: [0, 90)) can make.
_MAX_TILT_DEG = 89.9


def _rotation_matrix_x(angle_rad: float) -> np.ndarray:
    cosine, sine = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, cosine, -sine], [0.0, sine, cosine]], dtype=np.float64)


def _rotation_matrix_z(angle_rad: float) -> np.ndarray:
    cosine, sine = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


def _wrap180(angle_deg: float) -> float:
    """Wrap an angle to ``(-180, 180]``, so rotation deltas read as the short way round."""

    wrapped = (angle_deg + 180.0) % 360.0 - 180.0
    return 180.0 if wrapped == -180.0 else wrapped


def _stage_to_lab_matrix(tilt_deg: float, rotation_deg: float) -> np.ndarray:
    """Specimen-to-laboratory rotation of a tilt-plus-rotation SEM/ECCI stage.

    Extends the convention of :meth:`DiffractionGeometry.for_ebsd` — beam along
    laboratory z, tilt about laboratory x, untilted specimen normal facing the
    beam — with the second stage degree of freedom a real eucentric stage has:
    rotation about the specimen's own normal, applied before the tilt. At
    ``rotation_deg = 0`` this is exactly ``for_ebsd``'s own matrix, so an
    EBSD-measured orientation and geometry pass through unchanged.
    """

    return np.asarray(
        _rotation_matrix_x(math.radians(180.0 - tilt_deg))
        @ _rotation_matrix_z(math.radians(rotation_deg)),
        dtype=np.float64,
    )


def _beam_direction_specimen(tilt_deg: float, rotation_deg: float) -> np.ndarray:
    """The specimen-frame direction that lies along the beam at this stage state."""

    matrix = _stage_to_lab_matrix(tilt_deg, rotation_deg)
    return matrix.T @ np.array([0.0, 0.0, 1.0], dtype=np.float64)


def _beam_direction_crystal(orientation: Any, tilt_deg: float, rotation_deg: float) -> np.ndarray:
    """The crystal-frame direction that lies along the beam at this stage state."""

    beam_specimen = _beam_direction_specimen(tilt_deg, rotation_deg)
    crystal_to_specimen = np.asarray(orientation.rotation.as_matrix(), dtype=np.float64)
    return np.asarray(crystal_to_specimen.T @ beam_specimen, dtype=np.float64)


#: Step used to report which way a stage control moves the target, in degrees.
#:
#: Small enough that the reported displacement is the local derivative for any
#: practical purpose, large enough to stay well clear of floating-point noise.
_STAGE_PROBE_DEG = 1.0


def _target_in_lab(
    crystal_to_specimen: np.ndarray,
    target_crystal: np.ndarray,
    tilt_deg: float,
    rotation_deg: float,
) -> np.ndarray:
    """The target direction in the laboratory frame, where the beam is ``+z``.

    This is the frame the operator actually reasons in: the beam comes down
    ``+z``, so the target's ``(x, y)`` say how far off the beam it is and in
    which direction, and its ``z`` says whether it is pointing towards the gun
    or away from it.
    """

    stage = _stage_to_lab_matrix(tilt_deg, rotation_deg)
    direction = np.asarray(stage @ crystal_to_specimen @ target_crystal, dtype=np.float64)
    norm = float(np.linalg.norm(direction))
    return direction / norm if norm > 0.0 else direction


def _angle_from_beam_deg(target_crystal: np.ndarray, beam_crystal: np.ndarray) -> float:
    """How far the target direction lies from the beam, in degrees.

    Folded over the sense of both directions, because a zone axis and its
    opposite are the same axis: a target reached by tilting the crystal through
    the beam is on axis, not 180 degrees away from it.
    """

    beam_norm = float(np.linalg.norm(beam_crystal))
    if beam_norm <= 0.0:
        return 0.0
    cosine = abs(float(np.dot(target_crystal, beam_crystal)) / beam_norm)
    return float(math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0)))))


def _stage_view_payload(
    orientation: Any,
    target_crystal: np.ndarray,
    tilt_deg: float,
    rotation_deg: float,
) -> dict[str, Any]:
    """Where the target sits relative to the beam, and which way each control moves it.

    Purpose
    -------
    Turn "the target is 12 degrees off the beam" into something an operator can
    act on. The angle alone says how far there is to go but not which way, so a
    user moves a control, watches the number, and moves it back — which is the
    slow way to find a two-beam condition and the reason a solver exists at all.

    What it adds is the direction. The target's position is reported in the
    laboratory frame, where the beam is ``+z`` and the origin of the ``(x, y)``
    plane is the beam itself, so the browser can draw the target as a point that
    has to be walked to the centre. Alongside it are the displacements that one
    degree of tilt and one degree of rotation produce, which is what makes the
    two controls legible: they say, at this stage state, which way each knob
    pushes the target and how far.

    The displacements are finite differences over
    :data:`_STAGE_PROBE_DEG` rather than an analytic derivative. The stage
    geometry is a product of two rotations and differentiating it in closed form
    would be a second derivation to keep in step with the first; a difference of
    the same function the rest of the panel uses cannot drift away from it.
    """

    crystal_to_specimen = np.asarray(orientation.rotation.as_matrix(), dtype=np.float64)
    here = _target_in_lab(crystal_to_specimen, target_crystal, tilt_deg, rotation_deg)
    # A zone axis and its opposite are the same axis, and the deviation reported
    # everywhere else folds them. So does this: the target is always drawn in the
    # hemisphere facing the gun, or the marker would jump to the far side of the
    # plot the moment the direction passed through the beam, for no physical
    # reason at all. The probes are folded against the same base sense, so an
    # arrow can never be reversed by a fold that happened to it and not to it.
    if float(here[2]) < 0.0:
        here = -here
    probes = []
    for probe_tilt, probe_rotation in (
        (tilt_deg + _STAGE_PROBE_DEG, rotation_deg),
        (tilt_deg, rotation_deg + _STAGE_PROBE_DEG),
    ):
        probe = _target_in_lab(crystal_to_specimen, target_crystal, probe_tilt, probe_rotation)
        if float(np.dot(probe, here)) < 0.0:
            probe = -probe
        probes.append(probe)
    tilted, rotated = probes
    # The specimen normal is specimen +z carried into the laboratory frame. It
    # is what the schematic draws the holder around.
    normal = np.asarray(
        _stage_to_lab_matrix(tilt_deg, rotation_deg) @ np.array([0.0, 0.0, 1.0]),
        dtype=np.float64,
    )
    azimuth_deg = float(math.degrees(math.atan2(float(here[1]), float(here[0]))))
    return {
        "probe_step_deg": _STAGE_PROBE_DEG,
        "target_lab": [float(value) for value in here],
        "target_azimuth_deg": azimuth_deg,
        "specimen_normal_lab": [float(value) for value in normal],
        # The tilt axis is laboratory x by construction of the stage matrix.
        "tilt_axis_lab": [1.0, 0.0, 0.0],
        "per_tilt_degree": [float(tilted[0] - here[0]), float(tilted[1] - here[1])],
        "per_rotation_degree": [float(rotated[0] - here[0]), float(rotated[1] - here[1])],
    }


def _stage_branches(direction_specimen: np.ndarray) -> list[tuple[float, float]]:
    """Every ``(tilt_deg, rotation_deg)`` bringing a specimen direction onto the beam.

    Purpose
    -------
    The closed form behind the stage solver: given a direction already expressed
    in the specimen frame, the exact stage angles that place it along
    laboratory +z, for the ``Rx(180 - tilt) @ Rz(rotation)`` stage of
    :func:`_stage_to_lab_matrix`. Pure geometry — no envelope, no forward
    validation, no choice among branches; :func:`_solve_stage_for_direction`
    does that.

    Derivation
    ----------
    Write ``v = (v0, v1, v2)`` for the unit direction and ``rho = hypot(v0, v1)``.
    ``Rz(phi) @ v`` has first component zero exactly when
    ``phi = atan2(v0, v1)`` or that plus pi, giving second component
    ``+rho`` or ``-rho`` respectively (third component unchanged at ``v2``).
    ``Rx(theta)`` then zeroes the second component when
    ``theta = atan2(w1, v2)``, and because ``w1**2 + v2**2 = rho**2 + v2**2 = 1``
    for a unit ``v``, the third component that survives is exactly ``+1`` for
    both branches: both bring ``v`` onto ``+z``, not ``-z``. ``theta`` is stored
    as a stage tilt through ``theta = radians(180 - tilt_deg)``.

    When ``rho`` is degenerate (``v`` already along specimen z) the rotation is
    indeterminate and is reported as zero.

    Returns
    -------
    list of (float, float)
        ``(tilt_deg, rotation_deg)`` pairs, unfiltered — a caller checks the
        physical tilt range.
    """

    v0, v1, v2 = (float(value) for value in direction_specimen)
    rho = math.hypot(v0, v1)
    if rho < 1e-9:
        theta = math.atan2(0.0, v2)
        tilt_deg = 180.0 - math.degrees(theta)
        return [(tilt_deg, 0.0)]
    base_phi = math.atan2(v0, v1)
    branches: list[tuple[float, float]] = []
    for phi in (base_phi, base_phi + math.pi):
        sine_phi, cosine_phi = math.sin(phi), math.cos(phi)
        w1 = v0 * sine_phi + v1 * cosine_phi
        theta = math.atan2(w1, v2)
        tilt_deg = 180.0 - math.degrees(theta)
        rotation_deg = math.degrees(phi)
        branches.append((tilt_deg, _wrap180(rotation_deg)))
    return branches


def _forward_residual_deg(
    orientation: Any, direction_crystal_unit: np.ndarray, tilt_deg: float, rotation_deg: float
) -> float:
    """Re-derive the achieved beam direction and measure its angle from the target.

    Deliberately recomputed from the crystal direction through the full stage
    chain rather than reusing anything :func:`_stage_branches` produced, so an
    algebra error in the closed form cannot validate itself.
    """

    achieved_specimen = np.asarray(orientation.rotation.as_matrix(), dtype=np.float64) @ (
        direction_crystal_unit / np.linalg.norm(direction_crystal_unit)
    )
    achieved_lab = _stage_to_lab_matrix(tilt_deg, rotation_deg) @ achieved_specimen
    cosine = float(np.clip(achieved_lab[2] / np.linalg.norm(achieved_lab), -1.0, 1.0))
    return float(math.degrees(math.acos(cosine)))


def _solve_stage_for_direction(
    orientation: Any,
    direction_crystal: np.ndarray,
    *,
    current_tilt_deg: float,
    current_rotation_deg: float,
    allow_reverse: bool,
) -> list[dict[str, Any]]:
    """Rank the reachable stage moves that bring a crystal direction onto the beam.

    Every branch of :func:`_stage_branches`, for the direction and (if
    ``allow_reverse``) its antiparallel sense, is forward-validated and kept
    only if it lands within the tilt range a real stage — and
    ``DiffractionGeometry.for_ebsd`` — accepts. Survivors are sorted by total
    travel from the current stage state, shortest first.
    """

    direction_unit = direction_crystal / np.linalg.norm(direction_crystal)
    crystal_to_specimen = np.asarray(orientation.rotation.as_matrix(), dtype=np.float64)
    senses = (1.0, -1.0) if allow_reverse else (1.0,)
    solutions: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for sense in senses:
        direction_specimen = crystal_to_specimen @ (sense * direction_unit)
        for tilt_deg, rotation_deg in _stage_branches(direction_specimen):
            if not (0.0 <= tilt_deg < _MAX_TILT_DEG):
                continue
            key = (round(tilt_deg * 100.0), round(rotation_deg * 100.0))
            if key in seen:
                continue
            seen.add(key)
            residual_deg = _forward_residual_deg(
                orientation, sense * direction_unit, tilt_deg, rotation_deg
            )
            delta_tilt = tilt_deg - current_tilt_deg
            delta_rotation = _wrap180(rotation_deg - current_rotation_deg)
            travel_deg = math.hypot(delta_tilt, delta_rotation)
            solutions.append(
                {
                    "tilt_deg": tilt_deg,
                    "rotation_deg": rotation_deg,
                    "delta_tilt_deg": delta_tilt,
                    "delta_rotation_deg": delta_rotation,
                    "travel_deg": travel_deg,
                    "residual_deg": residual_deg,
                    "is_reversed": sense < 0.0,
                }
            )
    solutions.sort(key=lambda entry: entry["travel_deg"])
    return solutions


def _nearest_zone_axis(
    phase: Any, direction_crystal_unit: np.ndarray, *, max_index: int
) -> tuple[tuple[int, int, int], float]:
    """The low-index real-space direction closest to a crystal-frame direction.

    A small, local duplicate of the technique in
    :func:`pytex.app.services.tem._nearest_zone_axis`: it is a dozen lines
    behind a private name in another service's module, not something worth an
    inter-service import for. Both senses of a direction give the same zone-axis
    pattern, so they are folded together.
    """

    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=np.float64)
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid != 0, axis=1)]
    divisor = np.gcd(np.gcd(np.abs(grid[:, 0]), np.abs(grid[:, 1])), np.abs(grid[:, 2]))
    grid = grid // divisor[:, None]
    grid = np.unique(grid, axis=0)
    cartesian = grid.astype(np.float64) @ direct.T
    norms = np.linalg.norm(cartesian, axis=1)
    keep = norms > 1e-12
    cartesian = cartesian[keep] / norms[keep, None]
    grid = grid[keep]
    beam = direction_crystal_unit / np.linalg.norm(direction_crystal_unit)
    cosines = np.clip(np.abs(cartesian @ beam), -1.0, 1.0)
    best = int(np.argmax(cosines))
    deviation_deg = float(math.degrees(math.acos(cosines[best])))
    indices = (int(grid[best, 0]), int(grid[best, 1]), int(grid[best, 2]))
    if float(cartesian[best] @ beam) < 0.0:
        indices = (-indices[0], -indices[1], -indices[2])
    return indices, deviation_deg


def _on_axis_pattern(
    phase: Any,
    beam_crystal_unit: np.ndarray,
    *,
    zone_search_max_index: int,
    camera_length_mm: float,
    beam_energy_kev: float,
    detector_px: int,
    pixel_size_mm: float,
    max_index: int,
) -> dict[str, Any]:
    """Simulate what an on-axis detector shows, at the actual continuous beam direction.

    The pattern drawn is the exact zone-axis pattern of the nearest low-index
    direction — a spot pattern exists only on a rational zone, exactly as
    :mod:`pytex.app.services.tem`'s SAED simulator explains — but every spot's
    excitation error is recomputed against the *actual* beam direction rather
    than the nominal zone axis, using the same
    ``dot(g, beam_direction)`` definition
    :class:`~pytex.diffraction.saed.SAEDSpot` already uses for
    ``excitation_error_inv_angstrom`` (dot with the zone axis it was built
    from). This is the number that goes to zero, for every reflection of the
    zone, exactly when the stage solver's target lands on the beam — which is
    what demonstrates the two-beam condition has actually been reached.
    """

    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.kinematic import electron_wavelength_angstrom
    from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

    indices, deviation_deg = _nearest_zone_axis(
        phase, beam_crystal_unit, max_index=zone_search_max_index
    )
    axis = ZoneAxis(indices=np.asarray(indices, dtype=int), phase=phase)
    wavelength = float(electron_wavelength_angstrom(beam_energy_kev))
    camera_constant = float(camera_length_mm * wavelength)
    image = synthesize_saed_image(
        phase,
        axis,
        camera_constant_mm_angstrom=camera_constant,
        raster=DetectorRaster(
            width_px=detector_px, height_px=detector_px, pixel_size_mm=pixel_size_mm
        ),
        in_plane_rotation_deg=0.0,
        max_index=max_index,
    )
    if not image.spots:
        raise InvalidInputError(
            "No reflection of the nearest zone lands on the on-axis detector at this "
            "camera length.",
            field="on_axis_camera_length_mm",
            hint="Shorten the on-axis camera length, or enlarge the on-axis detector.",
        )
    payload = image.to_json()
    reciprocal = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=np.float64)
    beam_unit = beam_crystal_unit / np.linalg.norm(beam_crystal_unit)
    max_abs_excitation = 0.0
    for spot_payload, spot in zip(payload["spots"], image.spots, strict=True):
        g_crystal = reciprocal @ spot.miller_indices.astype(np.float64)
        excitation = float(np.dot(g_crystal, beam_unit))
        spot_payload["excitation_error_inv_angstrom"] = excitation
        max_abs_excitation = max(max_abs_excitation, abs(excitation))
    payload["nearest_zone_axis"] = [int(value) for value in indices]
    payload["nearest_zone_axis_deviation_deg"] = deviation_deg
    payload["max_abs_excitation_error_inv_angstrom"] = max_abs_excitation
    return payload


def _kikuchi_payload(
    phase: Any,
    spec: Any,
    orientation: Any,
    geometry: Any,
    *,
    max_bands: int,
    max_index: int,
    zone_axis_max_index: int,
) -> dict[str, Any]:
    """The EBSD Kikuchi pattern in the wire form ``ebsdkikuchi.js`` already draws.

    A lean re-derivation of the payload
    :func:`pytex.app.services.ebsd_pattern._simulate_kikuchi_pattern` builds,
    trimmed to what an embedded panel needs (no result table): both call the
    same :func:`~pytex.diffraction.kikuchi.simulate_kikuchi_pattern`, so the two
    payloads describe the identical geometry whenever the inputs agree.
    """

    from pytex.diffraction.kikuchi import simulate_kikuchi_pattern

    width_px, height_px = geometry.detector_shape[1], geometry.detector_shape[0]
    pixel_size_mm = geometry.detector_pixel_size_um[0] / 1000.0
    try:
        pattern = simulate_kikuchi_pattern(
            geometry,
            phase,
            orientation,
            max_index=max_index,
            max_bands=max_bands,
            zone_axis_max_index=zone_axis_max_index,
        )
    except ValueError as error:
        raise InvalidInputError(
            f"No EBSD pattern could be simulated at this stage state: {error}",
            field="max_index",
            hint="Raise the index limit, or check the phase carries an atomic basis.",
        ) from error

    projection = pattern.projection
    to_px = projection.to_detector_px
    bands: list[dict[str, Any]] = []
    for band in pattern.bands:
        indices = tuple(int(value) for value in band.plane.indices)
        centre = clipped_runs(
            np.asarray(to_px(band.center_trace(projection, samples=_TRACE_SAMPLES)), dtype=float),
            width_px,
            height_px,
        )
        edges = [
            clipped_runs(np.asarray(to_px(edge), dtype=float), width_px, height_px)
            for edge in band.edge_traces(projection, samples=_TRACE_SAMPLES)
        ]
        if not centre and not any(edges):
            continue
        gnomonic_width = band.width_at_pattern_center(projection)
        width_in_px = (
            float(gnomonic_width * geometry.camera_length_mm / pixel_size_mm)
            if math.isfinite(gnomonic_width)
            else float("inf")
        )
        bands.append(
            {
                "hkl": list(indices),
                "label": plane_label(indices, spec=spec),
                "d_angstrom": float(band.d_spacing_angstrom),
                "bragg_angle_deg": float(np.degrees(band.bragg_angle_rad)),
                "width_deg": float(np.degrees(band.angular_width_rad)),
                "width_px": width_in_px,
                "intensity": float(band.intensity),
                "centre": centre,
                "edges": edges,
            }
        )
    axes: list[dict[str, Any]] = []
    for axis in pattern.zone_axes:
        point = np.asarray(to_px(np.asarray(axis.coordinates, dtype=float).reshape(1, 2)))[0]
        indices = tuple(int(value) for value in axis.indices)
        axes.append(
            {
                "uvw": list(indices),
                "label": direction_label(indices, spec=spec),
                "x": float(point[0]),
                "y": float(point[1]),
                "order": int(axis.band_count),
                "on_detector": bool(axis.on_detector),
            }
        )
    if not bands:
        raise InvalidInputError(
            "Every simulated EBSD band misses the screen at this geometry.",
            field="detector_distance",
            hint="Reduce z*, or move the pattern centre back towards the middle of the screen.",
        )
    return {
        "width_px": width_px,
        "height_px": height_px,
        "pattern_centre_px": [
            float(geometry.pattern_center_px[0]),
            float(geometry.pattern_center_px[1]),
        ],
        "camera_length_mm": float(geometry.camera_length_mm),
        "pixel_size_mm": pixel_size_mm,
        "wavelength_angstrom": float(pattern.wavelength_angstrom),
        "bands": bands,
        "zone_axes": axes,
    }


def _build_geometry(request: dict[str, Any], *, tilt_deg: float, rotation_deg: float) -> Any:
    """The EBSD ``DiffractionGeometry`` at a given stage tilt and rotation."""

    from pytex.diffraction.models import DiffractionGeometry

    base = DiffractionGeometry.for_ebsd(
        beam_energy_kev=float(request["beam_energy_kev"]),
        sample_tilt_deg=tilt_deg,
        detector_elevation_deg=float(request["detector_elevation_deg"]),
        detector_azimuth_deg=float(request["detector_azimuth_deg"]),
        pattern_center=(
            float(request["pattern_centre_x"]),
            float(request["pattern_centre_y"]),
            float(request["detector_distance"]),
        ),
        detector_shape=(int(request["detector_height_px"]), int(request["detector_width_px"])),
        detector_pixel_size_um=(float(request["pixel_size_um"]), float(request["pixel_size_um"])),
    )
    return dataclasses.replace(
        base, specimen_to_lab_matrix=_stage_to_lab_matrix(tilt_deg, rotation_deg)
    )


def _build_orientation(request: dict[str, Any], phase: Any) -> Any:
    from pytex.core.frame_catalog import SPECIMEN_FRAME
    from pytex.core.orientation import Orientation

    return Orientation.from_euler(
        float(request["phi1_deg"]),
        float(request["Phi_deg"]),
        float(request["phi2_deg"]),
        degrees=True,
        specimen_frame=SPECIMEN_FRAME,
        phase=phase,
    )


def _target_direction_crystal(request: dict[str, Any], phase: Any) -> np.ndarray:
    from pytex.core.lattice import CrystalDirection

    indices = tuple(int(value) for value in request["target_zone_axis"])
    direction = CrystalDirection(coordinates=np.asarray(indices, dtype=float), phase=phase)
    return direction.unit_vector


def _shared_parameters(*, tilt_default: float = 70.0, rotation_default: float = 0.0) -> tuple:
    return (
        phase_parameter(help_text="The phase whose bands and reflections are simulated."),
        NumberParameter(
            name="phi1_deg",
            label="First Bunge angle",
            symbol="phi_1",
            help_text=(
                "First Bunge angle of the EBSD-measured orientation, about the specimen z axis."
            ),
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="Phi_deg",
            label="Second Bunge angle",
            symbol="Phi",
            help_text="Second Bunge angle of the EBSD-measured orientation, about the new x axis.",
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="phi2_deg",
            label="Third Bunge angle",
            symbol="phi_2",
            help_text="Third Bunge angle of the EBSD-measured orientation, about the new z axis.",
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="stage_tilt_deg",
            label="Stage tilt",
            help_text=(
                "Current stage tilt about the tilt axis — 70 deg for an as-measured EBSD point."
            ),
            units="deg",
            default=tilt_default,
            minimum=0.0,
            maximum=89.0,
            group="Stage (current)",
        ),
        NumberParameter(
            name="stage_rotation_deg",
            label="Stage rotation",
            help_text=(
                "Current stage rotation about the specimen normal, applied before the tilt. Most "
                "EBSD acquisitions are taken at the stage's zero-rotation reading."
            ),
            units="deg",
            default=rotation_default,
            minimum=-180.0,
            maximum=180.0,
            group="Stage (current)",
        ),
        NumberParameter(
            name="detector_elevation_deg",
            label="EBSD camera elevation",
            help_text=(
                "How far the EBSD camera axis is raised above the plane perpendicular to the beam."
            ),
            units="deg",
            default=0.0,
            minimum=-40.0,
            maximum=40.0,
            group="EBSD camera",
        ),
        NumberParameter(
            name="detector_azimuth_deg",
            label="EBSD camera azimuth",
            help_text="Rotation of the EBSD camera about the beam, from the nominal port.",
            units="deg",
            default=0.0,
            minimum=-180.0,
            maximum=180.0,
            group="EBSD camera",
        ),
        NumberParameter(
            name="pattern_centre_x",
            label="Pattern centre x*",
            help_text=(
                "Where the EBSD camera axis meets the screen, as a fraction of the screen width."
            ),
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            group="EBSD camera",
        ),
        NumberParameter(
            name="pattern_centre_y",
            label="Pattern centre y*",
            help_text=(
                "Where the EBSD camera axis meets the screen, as a fraction of the screen height."
            ),
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            group="EBSD camera",
        ),
        NumberParameter(
            name="detector_distance",
            label="EBSD camera distance z*",
            help_text="EBSD specimen-to-screen distance, as a fraction of the screen width.",
            default=0.65,
            minimum=0.1,
            maximum=3.0,
            group="EBSD camera",
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text=(
                "Sets the electron wavelength for both the EBSD pattern and the on-axis pattern."
            ),
            units="kV",
            default=20.0,
            minimum=5.0,
            maximum=40.0,
            group="EBSD camera",
        ),
        IntegerParameter(
            name="detector_width_px",
            label="EBSD screen width",
            units="px",
            default=640,
            minimum=64,
            maximum=2048,
            group="EBSD camera",
            advanced=True,
            help_text="EBSD camera width in pixels.",
        ),
        IntegerParameter(
            name="detector_height_px",
            label="EBSD screen height",
            units="px",
            default=480,
            minimum=64,
            maximum=2048,
            group="EBSD camera",
            advanced=True,
            help_text="EBSD camera height in pixels.",
        ),
        NumberParameter(
            name="pixel_size_um",
            label="EBSD pixel size",
            units="um",
            default=50.0,
            minimum=1.0,
            maximum=500.0,
            group="EBSD camera",
            advanced=True,
            help_text="EBSD screen pixel pitch.",
        ),
        IntegerParameter(
            name="max_bands",
            label="EBSD bands drawn",
            default=24,
            minimum=1,
            maximum=120,
            help_text="Keep the strongest this many EBSD bands.",
            group="EBSD camera",
        ),
        IntegerParameter(
            name="max_index",
            label="EBSD index limit",
            default=4,
            minimum=1,
            maximum=8,
            advanced=True,
            group="EBSD camera",
            help_text="Largest |h|, |k| or |l| enumerated for the EBSD pattern.",
        ),
        IntegerParameter(
            name="zone_axis_max_index",
            label="EBSD zone-axis limit",
            default=3,
            minimum=1,
            maximum=6,
            advanced=True,
            group="EBSD camera",
            help_text="Largest |u|, |v| or |w| enumerated for the EBSD pattern's zone axes.",
        ),
        NumberParameter(
            name="on_axis_camera_length_mm",
            label="On-axis camera length",
            help_text=(
                "Effective specimen-to-detector distance for the on-axis BSE/SAED-mode view. Sets "
                "the display scale only — it does not change which reflections are near the beam."
            ),
            units="mm",
            default=100.0,
            minimum=1.0,
            maximum=2000.0,
            group="On-axis detector",
        ),
        IntegerParameter(
            name="on_axis_detector_px",
            label="On-axis detector size",
            units="px",
            default=512,
            minimum=64,
            maximum=2048,
            group="On-axis detector",
            advanced=True,
            help_text="On-axis detector raster size, in pixels, square.",
        ),
        NumberParameter(
            name="on_axis_pixel_size_mm",
            label="On-axis pixel size",
            units="mm",
            default=0.05,
            minimum=0.001,
            maximum=1.0,
            group="On-axis detector",
            advanced=True,
            help_text="On-axis detector pixel pitch.",
        ),
        IntegerParameter(
            name="on_axis_max_index",
            label="On-axis reflection limit",
            default=6,
            minimum=1,
            maximum=10,
            advanced=True,
            group="On-axis detector",
            help_text="Largest |h|, |k| or |l| enumerated for the on-axis pattern.",
        ),
        IntegerParameter(
            name="zone_search_max_index",
            label="Zone search limit",
            default=4,
            minimum=1,
            maximum=8,
            advanced=True,
            group="On-axis detector",
            help_text=(
                "How far to search for the low-index zone axis nearest the actual beam direction."
            ),
        ),
        IndicesParameter(
            name="target_zone_axis",
            label="Target direction [uvw]",
            width=3,
            default=(0, 0, 1),
            help_text=(
                "The crystallographic direction to bring onto the beam for the two-beam condition. "
                "Its proximity to the current beam direction is always reported."
            ),
        ),
    )


@REGISTRY.operation(
    "ecci.solve_workflow",
    title="Solve the ECCI two-beam tilt",
    summary=(
        "From an EBSD orientation: the EBSD pattern, the on-axis view, and the tilt to a target."
    ),
    help_text=(
        "Starts from an orientation an EBSD system has already indexed and answers the question "
        "an ECCI operator actually asks: *what stage move brings this direction onto the beam?*\n\n"
        "**Three views of one crystal.** The EBSD Kikuchi pattern, exactly as the EBSD camera "
        "records it at the stated stage tilt and rotation. The on-axis view — what a BSE detector "
        "on the beam axis sees, in TEM-style on-axis diffraction mode, down whichever low-index "
        "direction the beam currently sits nearest. And the stage tilt and rotation that bring a "
        "chosen direction onto the beam exactly, ranked by how far the stage has to move.\n\n"
        "**The stage, stated explicitly.** One mechanical tilt about a fixed laboratory axis, plus "
        "a rotation about the specimen's own normal applied before the tilt — the two degrees of "
        "freedom a eucentric SEM/ECCI stage actually has, which is not the double-tilt holder a "
        "TEM uses.\n\n"
        "**Forward-validated.** Every candidate move is checked by re-deriving the direction it "
        "actually places on the beam, so the residual reported is never the solver trusting its "
        "own algebra.\n\n"
        "**What it is not.** Kinematic and geometric, like every pattern in this workspace: "
        "two-beam "
        "*contrast* — which reflection dominates, how strong the channelling signal is — is not "
        "modelled, only the geometric condition that makes a two-beam setup possible."
    ),
    parameters=(
        *_shared_parameters(),
        BooleanParameter(
            name="allow_reverse",
            label="Allow the opposite sense",
            help_text=(
                "Treat [uvw] and [-u-v-w] as the same target; usually halves the stage travel."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns=(
        "One row per candidate stage move; the EBSD and on-axis patterns at the current state "
        "under `data`."
    ),
    panel="ecci",
    citations=(_CITATION_ECCI, _CITATION_WILLIAMS, _CITATION_SCHWARTZ),
    tags=("EBSD", "ECCI", "channelling", "two-beam", "tilt", "on-axis", "zone axis", "dislocation"),
)
def _solve_workflow(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    orientation = _build_orientation(request, phase)
    tilt_deg = float(request["stage_tilt_deg"])
    rotation_deg = float(request["stage_rotation_deg"])

    geometry = _build_geometry(request, tilt_deg=tilt_deg, rotation_deg=rotation_deg)
    kikuchi = _kikuchi_payload(
        phase,
        spec,
        orientation,
        geometry,
        max_bands=int(request["max_bands"]),
        max_index=int(request["max_index"]),
        zone_axis_max_index=int(request["zone_axis_max_index"]),
    )

    beam_crystal = _beam_direction_crystal(orientation, tilt_deg, rotation_deg)
    on_axis = _on_axis_pattern(
        phase,
        beam_crystal,
        zone_search_max_index=int(request["zone_search_max_index"]),
        camera_length_mm=float(request["on_axis_camera_length_mm"]),
        beam_energy_kev=float(request["beam_energy_kev"]),
        detector_px=int(request["on_axis_detector_px"]),
        pixel_size_mm=float(request["on_axis_pixel_size_mm"]),
        max_index=int(request["on_axis_max_index"]),
    )

    target_direction = _target_direction_crystal(request, phase)
    target_indices = tuple(int(value) for value in request["target_zone_axis"])
    target_label = direction_label(target_indices, spec=spec)

    solutions = _solve_stage_for_direction(
        orientation,
        target_direction,
        current_tilt_deg=tilt_deg,
        current_rotation_deg=rotation_deg,
        allow_reverse=bool(request["allow_reverse"]),
    )
    if not solutions:
        raise InvalidInputError(
            f"No stage tilt within [0, {_MAX_TILT_DEG:g}) deg brings {target_label} onto the beam "
            "from this orientation.",
            field="target_zone_axis",
            hint="Choose a different target direction, or allow the opposite sense.",
        )

    rows = [
        {
            "rank": index + 1,
            "tilt_deg": float(solution["tilt_deg"]),
            "rotation_deg": float(solution["rotation_deg"]),
            "delta_tilt_deg": float(solution["delta_tilt_deg"]),
            "delta_rotation_deg": float(solution["delta_rotation_deg"]),
            "travel_deg": float(solution["travel_deg"]),
            "residual_deg": float(solution["residual_deg"]),
            "sense": "reversed" if solution["is_reversed"] else "as given",
        }
        for index, solution in enumerate(solutions)
    ]
    best = solutions[0]

    proximity_indices, proximity_deviation = _nearest_zone_axis(
        phase, beam_crystal, max_index=int(request["zone_search_max_index"])
    )
    proximity_label = direction_label(proximity_indices, spec=spec)

    result = AppResult(
        title=f"ECCI tilt to {target_label} from {spec.name} at Bunge "
        f"({request['phi1_deg']:.1f}, {request['Phi_deg']:.1f}, {request['phi2_deg']:.1f}) deg",
        summary=(
            f"At the current stage state (tilt {tilt_deg:.1f} deg, rotation "
            f"{rotation_deg:.1f} deg) "
            f"the beam sits {proximity_deviation:.2f} deg from {proximity_label}, the nearest "
            f"low-index direction. Bringing {target_label} onto the beam for a two-beam condition "
            f"needs a stage tilt of {best['tilt_deg']:.2f} deg and a rotation of "
            f"{best['rotation_deg']:.2f} deg — a move of {best['delta_tilt_deg']:+.2f} deg in tilt "
            f"and {best['delta_rotation_deg']:+.2f} deg in rotation, forward-validated to "
            f"{best['residual_deg']:.4f} deg off the beam. {len(rows)} reachable stage "
            "move(s) were "
            "found within the [0, 90) deg tilt range."
        ),
        table=ResultTable(
            columns=(
                Column("rank", "Rank", numeric=True),
                Column("sense", "Sense"),
                Column("tilt_deg", "Tilt", units="deg", numeric=True, digits=2),
                Column("rotation_deg", "Rotation", units="deg", numeric=True, digits=2),
                Column("delta_tilt_deg", "Delta tilt", units="deg", numeric=True, digits=2),
                Column("delta_rotation_deg", "Delta rotation", units="deg", numeric=True, digits=2),
                Column("travel_deg", "Total travel", units="deg", numeric=True, digits=2),
                Column(
                    "residual_deg",
                    "Off axis after move",
                    units="deg",
                    numeric=True,
                    digits=4,
                    help_text=(
                        "Forward-validated angle between the target direction and the beam after "
                        "the move."
                    ),
                ),
            ),
            rows=tuple(rows),
            caption=f"Stage moves bringing {target_label} onto the beam, shortest travel first.",
        ),
        data={
            "kikuchi": kikuchi,
            "on_axis": on_axis,
            "proximity": {
                "indices": list(proximity_indices),
                "label": proximity_label,
                "deviation_deg": proximity_deviation,
            },
            "target": {
                "indices": list(target_indices),
                "label": target_label,
                # The same deviation the live operation reports, so the stage
                # console is populated by the first solve rather than staying
                # blank until something is moved.
                "angle_from_beam_deg": _angle_from_beam_deg(target_direction, beam_crystal),
            },
            "solution": rows[0],
            "current": {"tilt_deg": tilt_deg, "rotation_deg": rotation_deg},
            "stage_view": _stage_view_payload(
                orientation, target_direction, tilt_deg, rotation_deg
            ),
        },
        inputs={
            "phase": spec.name,
            "euler_deg": [
                float(request["phi1_deg"]),
                float(request["Phi_deg"]),
                float(request["phi2_deg"]),
            ],
            "stage_tilt_deg": tilt_deg,
            "stage_rotation_deg": rotation_deg,
            "target_zone_axis": list(target_indices),
        },
        notes=(
            "Kinematic and geometric, as every simulated pattern in this workspace is: band and "
            "spot positions are exact for the stated geometry, intensities are a kinematic proxy, "
            "and dynamical two-beam contrast is not modelled.",
            "The on-axis pattern shows the zone axis nearest the actual beam direction; its "
            "reflections' excitation errors are computed against the actual continuous beam "
            "direction, not the nominal zone, so they read near zero exactly when the beam is "
            "genuinely close to that axis.",
            "The stage solved for has one tilt about a fixed axis and one rotation about the "
            "specimen normal applied before the tilt - a eucentric SEM/ECCI stage, not a TEM "
            "double-tilt holder.",
        ),
        citations=(_CITATION_ECCI, _CITATION_WILLIAMS, _CITATION_SCHWARTZ),
    )
    return result.to_json()


@REGISTRY.operation(
    "ecci.resimulate",
    title="Re-simulate the ECCI patterns at a given tilt and rotation",
    summary="The EBSD and on-axis patterns at an explicit stage state, without re-solving.",
    help_text=(
        "The live counterpart to solving the workflow: given the same crystal and camera "
        "description, plus an explicit stage tilt and rotation, recompute both patterns and the "
        "zone-axis proximity at that state. This is what a tilt/rotation slider calls on every "
        "move, and what a notebook calls to check a solved tilt actually reaches the two-beam "
        "condition it was solved for."
    ),
    parameters=_shared_parameters(),
    returns="The EBSD and on-axis patterns and the zone-axis proximity at the stated stage state.",
    panel="ecci",
    citations=(_CITATION_ECCI, _CITATION_WILLIAMS),
    tags=("EBSD", "ECCI", "channelling", "two-beam", "on-axis", "zone axis", "live"),
)
def _resimulate(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    orientation = _build_orientation(request, phase)
    tilt_deg = float(request["stage_tilt_deg"])
    rotation_deg = float(request["stage_rotation_deg"])

    geometry = _build_geometry(request, tilt_deg=tilt_deg, rotation_deg=rotation_deg)
    kikuchi = _kikuchi_payload(
        phase,
        spec,
        orientation,
        geometry,
        max_bands=int(request["max_bands"]),
        max_index=int(request["max_index"]),
        zone_axis_max_index=int(request["zone_axis_max_index"]),
    )

    beam_crystal = _beam_direction_crystal(orientation, tilt_deg, rotation_deg)
    on_axis = _on_axis_pattern(
        phase,
        beam_crystal,
        zone_search_max_index=int(request["zone_search_max_index"]),
        camera_length_mm=float(request["on_axis_camera_length_mm"]),
        beam_energy_kev=float(request["beam_energy_kev"]),
        detector_px=int(request["on_axis_detector_px"]),
        pixel_size_mm=float(request["on_axis_pixel_size_mm"]),
        max_index=int(request["on_axis_max_index"]),
    )

    proximity_indices, proximity_deviation = _nearest_zone_axis(
        phase, beam_crystal, max_index=int(request["zone_search_max_index"])
    )
    proximity_label = direction_label(proximity_indices, spec=spec)

    target_indices = tuple(int(value) for value in request["target_zone_axis"])
    target_direction = _target_direction_crystal(request, phase)
    target_angle_deg = _angle_from_beam_deg(target_direction, beam_crystal)

    result = AppResult(
        title=(
            f"ECCI patterns of {spec.name} at tilt {tilt_deg:.2f} deg, rotation "
            f"{rotation_deg:.2f} deg"
        ),
        summary=(
            f"At tilt {tilt_deg:.2f} deg and rotation {rotation_deg:.2f} deg the beam sits "
            f"{proximity_deviation:.2f} deg from {proximity_label} and {target_angle_deg:.2f} deg "
            f"from the target direction {direction_label(target_indices, spec=spec)}. The on-axis "
            f"pattern's largest reflection excitation error is "
            f"{on_axis['max_abs_excitation_error_inv_angstrom']:.4f} inverse angstrom."
        ),
        data={
            "kikuchi": kikuchi,
            "on_axis": on_axis,
            "proximity": {
                "indices": list(proximity_indices),
                "label": proximity_label,
                "deviation_deg": proximity_deviation,
            },
            "target": {
                "indices": list(target_indices),
                "label": direction_label(target_indices, spec=spec),
                "angle_from_beam_deg": target_angle_deg,
            },
            "state": {"tilt_deg": tilt_deg, "rotation_deg": rotation_deg},
            "stage_view": _stage_view_payload(
                orientation, target_direction, tilt_deg, rotation_deg
            ),
        },
        inputs={
            "phase": spec.name,
            "euler_deg": [
                float(request["phi1_deg"]),
                float(request["Phi_deg"]),
                float(request["phi2_deg"]),
            ],
            "stage_tilt_deg": tilt_deg,
            "stage_rotation_deg": rotation_deg,
        },
        notes=(
            "Kinematic and geometric: band and spot positions are exact for the stated geometry, "
            "intensities are a kinematic proxy.",
        ),
        citations=(_CITATION_ECCI, _CITATION_WILLIAMS),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="ecci.example.nickel_cube_to_111",
            title="Nickel: cube orientation, tilt to [111]",
            panel="ecci",
            summary=(
                "A cube-oriented fcc grain, EBSD-measured at 70 deg, tilted to a <111> "
                "two-beam condition."
            ),
            teaches=(
                "The stage move needed is not the crystal misorientation from [001] to [111] "
                "(54.7 deg): with the specimen already tilted 70 deg for EBSD, the stage tilt that "
                "lands [111] on the beam is a very different, non-obvious number, which is exactly "
                "why this is solved rather than estimated by eye."
            ),
            operation="ecci.solve_workflow",
            request={
                "phase": {"builtin": "ni_fcc"},
                "phi1_deg": 0.0,
                "Phi_deg": 0.0,
                "phi2_deg": 0.0,
                "stage_tilt_deg": 70.0,
                "stage_rotation_deg": 0.0,
                "target_zone_axis": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="ecci.example.ferrite_near_001",
            title="Ferrite near [001], tilt back on-axis",
            panel="ecci",
            summary="A bcc grain close to cube, brought exactly onto its own [001] for ECCI.",
            teaches=(
                "The nearest-zone-axis readout catches a grain that is already close to a "
                "low-index "
                "pole: the on-axis pattern shows the [001] zone from the first response, and the "
                "solved tilt is a small correction rather than a large excursion."
            ),
            operation="ecci.solve_workflow",
            request={
                "phase": {"builtin": "fe_bcc"},
                "phi1_deg": 5.0,
                "Phi_deg": 3.0,
                "phi2_deg": 2.0,
                "stage_tilt_deg": 70.0,
                "stage_rotation_deg": 0.0,
                "target_zone_axis": [0, 0, 1],
            },
        ),
        ExampleScenario(
            id="ecci.example.zirconium_basal",
            title="Zirconium hcp: tilt to the basal pole",
            panel="ecci",
            summary="An hcp grain, EBSD-measured near the basal orientation, tilted onto [0001].",
            teaches=(
                "The same solver and the same on-axis machinery work unchanged on a non-cubic "
                "phase: the target direction is resolved through the phase's own direct lattice "
                "basis, so a hexagonal c-axis is handled exactly as a cubic pole is."
            ),
            operation="ecci.solve_workflow",
            request={
                "phase": {"builtin": "zr_hcp"},
                "phi1_deg": 0.0,
                "Phi_deg": 10.0,
                "phi2_deg": 0.0,
                "stage_tilt_deg": 70.0,
                "stage_rotation_deg": 0.0,
                "target_zone_axis": [0, 0, 1],
            },
        ),
    )
)
