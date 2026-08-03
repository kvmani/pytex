"""Solving a measured SAED pattern from a calibrated list of spot positions.

The question this module answers is the one asked in front of a microscope:
*here are the spots I picked — which phase is this, down which zone axis, and
what is every spot?*

The input is deliberately minimal: spot positions relative to the transmitted
beam, plus enough calibration to turn them into reciprocal-space lengths. That
is all a printed pattern gives you, and it is all this needs. Spots may be
clicked interactively (see `pytex.plotting.saed_picker`) or listed in a YAML
file; the file is the contract, so a solved pattern is reproducible from a
committed text file.

Relation to `pytex.diffraction.models.index_saed_pattern`
--------------------------------------------------------
That surface solves a different problem. It starts from a `DiffractionGeometry`
— a calibrated detector with a known distance, pattern centre, tilt and pixel
size — and works in detector pixels, which is what an automated indexing
pipeline attached to an instrument has. This module starts from a bare list of
spot coordinates and a camera constant, which is what a person reading a
micrograph has. Use `index_saed_pattern` when the detector model is known; use
`solve_saed_pattern` when only the pattern is.

Algorithm
---------
Classical ratio/angle indexing (Edington; Williams and Carter): two
non-collinear reflections fix the zone, so the two shortest measured vectors
seed the solution. A candidate assignment is admissible when both lengths match
a calculated reflection within a relative tolerance and their interplanar angle
matches within an angular tolerance. The zone axis follows from the cross
product, the crystal-to-detector rotation from aligning the calculated pair onto
the observed pair, and every remaining spot is then indexed by projection and
scored. Solutions are ranked by matched fraction, then by residual.

Intensities are never used: a kinematic intensity model is not reliable enough
to index against, and a printed pattern rarely carries calibrated intensities at
all. Geometry alone decides.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pytex.core._arrays import as_float_array
from pytex.core.hexagonal import is_hexagonal_phase, plane_hkl_to_hkil
from pytex.core.lattice import CrystalDirection, Phase, ZoneAxis
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.orientation import Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import rationalize_zone_axis
from pytex.diffraction.kinematic import (
    centering_allowed_mask,
    electron_wavelength_angstrom,
)
from pytex.diffraction.physics import ReflectionCondition

#: Schema identifier of the measured-pattern YAML contract.
MEASURED_PATTERN_SCHEMA = "pytex.measured_saed_pattern/1"

#: Schema identifier of the solution-report payload.
PATTERN_SOLUTION_SCHEMA = "pytex.pattern_solution_report/1"

#: Coordinate units a measured pattern may be expressed in.
COORDINATE_UNITS: tuple[str, ...] = ("px", "mm", "reciprocal_angstrom")

#: Default relative tolerance on a reflection length.
#:
#: Three percent is the working figure for spot positions read off a printed or
#: digitized pattern: it comfortably covers the centring error of a hand-picked
#: spot at typical camera constants while still separating, for example, the
#: {111} and {200} rings of an fcc metal, whose lengths differ by 15 percent.
DEFAULT_LENGTH_TOLERANCE_RELATIVE = 0.03

#: Default tolerance on an interplanar angle, in degrees.
DEFAULT_ANGLE_TOLERANCE_DEG = 2.0

#: How many seed pairs, and how many admissible index assignments per seed pair,
#: are carried into full verification.
#:
#: Symmetry makes many assignments equivalent, and verifying all of them would
#: re-derive the same solution dozens of times. These bounds keep the cost
#: predictable; solutions are deduplicated afterwards, so raising them does not
#: change the answer on a well-posed pattern, only the time taken.
_MAX_SEED_PAIRS = 6
_MAX_ASSIGNMENTS_PER_SEED = 64


@dataclass(frozen=True, slots=True)
class PatternCalibration:
    """How to turn picked spot coordinates into reciprocal-space vectors.

    ``units`` names what the coordinates are: detector pixels, millimetres on
    the detector, or reciprocal angstroms (already calibrated). Pixels need
    ``pixel_size_mm``; pixels and millimetres both need a camera constant
    ``L·lambda``, supplied directly as ``camera_constant_mm_angstrom`` or
    derived from ``camera_length_mm`` and ``beam_energy_kev``.

    ``centre`` is the transmitted-beam position **in the same units as the
    coordinates**, subtracted before anything else. A pattern already measured
    relative to the beam leaves it at the origin.
    """

    units: str = "reciprocal_angstrom"
    camera_constant_mm_angstrom: float | None = None
    camera_length_mm: float | None = None
    beam_energy_kev: float | None = None
    pixel_size_mm: float | None = None
    centre: tuple[float, float] = (0.0, 0.0)

    def __post_init__(self) -> None:
        if self.units not in COORDINATE_UNITS:
            raise ValueError(
                "PatternCalibration.units must be one of " + ", ".join(COORDINATE_UNITS) + "."
            )
        if self.units == "px" and (
            self.pixel_size_mm is None or self.pixel_size_mm <= 0.0
        ):
            raise ValueError("Pixel coordinates require a positive pixel_size_mm.")
        if self.camera_length_mm is not None and self.camera_length_mm <= 0.0:
            raise ValueError("camera_length_mm must be strictly positive when given.")
        if self.beam_energy_kev is not None and self.beam_energy_kev <= 0.0:
            raise ValueError("beam_energy_kev must be strictly positive when given.")
        if (
            self.camera_constant_mm_angstrom is not None
            and self.camera_constant_mm_angstrom <= 0.0
        ):
            raise ValueError(
                "camera_constant_mm_angstrom must be strictly positive when given."
            )
        if self.units != "reciprocal_angstrom":
            # Fail here rather than at the first spot conversion: an
            # uncalibrated length is not a recoverable state.
            self._require_camera_constant()
        object.__setattr__(
            self, "centre", (float(self.centre[0]), float(self.centre[1]))
        )

    def _require_camera_constant(self) -> float:
        if self.camera_constant_mm_angstrom is not None:
            return float(self.camera_constant_mm_angstrom)
        if self.camera_length_mm is not None and self.beam_energy_kev is not None:
            return float(self.camera_length_mm) * electron_wavelength_angstrom(
                float(self.beam_energy_kev)
            )
        raise ValueError(
            f"Coordinates in '{self.units}' need a camera constant: supply "
            "camera_constant_mm_angstrom, or camera_length_mm together with "
            "beam_energy_kev."
        )

    @property
    def effective_camera_constant_mm_angstrom(self) -> float | None:
        """The camera constant ``L·lambda`` in mm·angstrom, if one applies."""

        if self.units == "reciprocal_angstrom":
            return None
        return self._require_camera_constant()

    def to_reciprocal_angstrom(self, coordinates: Any) -> np.ndarray:
        """Convert raw picked coordinates to in-plane ``g`` vectors in 1/angstrom.

        The centre is subtracted first, then pixels are scaled to millimetres,
        then millimetres are divided by the camera constant. A vanishing camera
        constant or pixel size cannot occur — both are validated at construction.
        """

        raw = as_float_array(
            np.asarray(coordinates, dtype=np.float64).reshape(-1, 2), shape=(None, 2)
        )
        centred = raw - np.asarray(self.centre, dtype=np.float64)
        if self.units == "reciprocal_angstrom":
            return np.ascontiguousarray(centred)
        if self.units == "px":
            # Construction guarantees a positive pixel size for pixel units.
            assert self.pixel_size_mm is not None
            millimetres = centred * float(self.pixel_size_mm)
        else:
            millimetres = centred
        return np.ascontiguousarray(millimetres / self._require_camera_constant())

    def describe(self) -> str:
        """Prose summary of how coordinates become reciprocal-space lengths."""

        if self.units == "reciprocal_angstrom":
            return (
                "Spot coordinates are already reciprocal-space vectors in 1/angstrom, "
                f"measured from a transmitted beam at {self.centre}."
            )
        constant = self._require_camera_constant()
        source = (
            "given directly"
            if self.camera_constant_mm_angstrom is not None
            else (
                f"derived from a {self.camera_length_mm:g} mm camera length at "
                f"{self.beam_energy_kev:g} kV"
            )
        )
        scale = (
            f" Pixels are scaled by {self.pixel_size_mm:g} mm/px first."
            if self.units == "px"
            else ""
        )
        return (
            f"Spot coordinates are in {self.units}, measured from a transmitted beam at "
            f"{self.centre}.{scale} Reciprocal lengths follow from |g| = r / (L*lambda) "
            f"with a camera constant of {constant:.4f} mm*angstrom ({source})."
        )

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {"units": self.units, "centre": list(self.centre)}
        for name in (
            "camera_constant_mm_angstrom",
            "camera_length_mm",
            "beam_energy_kev",
            "pixel_size_mm",
        ):
            value = getattr(self, name)
            if value is not None:
                payload[name] = float(value)
        return payload


@dataclass(frozen=True, slots=True)
class MeasuredSpot:
    """One picked diffraction spot, in the pattern's own coordinate units."""

    position: tuple[float, float]
    intensity: float | None = None
    label: str | None = None

    def __post_init__(self) -> None:
        position = (float(self.position[0]), float(self.position[1]))
        if not all(np.isfinite(position)):
            raise ValueError("MeasuredSpot.position must be finite.")
        if self.intensity is not None and (
            not np.isfinite(self.intensity) or self.intensity < 0.0
        ):
            raise ValueError("MeasuredSpot.intensity must be finite and non-negative.")
        object.__setattr__(self, "position", position)


@dataclass(frozen=True, slots=True)
class MeasuredSAEDPattern:
    """A measured pattern: picked spots plus the calibration that scales them.

    The transmitted beam is **not** a spot: it is the calibration's ``centre``,
    and every position is taken relative to it. A spot at the centre would have
    no direction and is rejected.
    """

    name: str
    spots: tuple[MeasuredSpot, ...]
    calibration: PatternCalibration
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        spots = tuple(self.spots)
        if len(spots) < 2:
            raise ValueError(
                "A pattern needs at least two spots: two non-collinear reflections are "
                "what fix a zone axis."
            )
        object.__setattr__(self, "spots", spots)
        if not self.name.strip():
            raise ValueError("MeasuredSAEDPattern.name must be non-empty.")
        magnitudes = np.linalg.norm(self.g_vectors_inv_angstrom(), axis=1)
        if np.any(magnitudes <= 0.0):
            raise ValueError(
                "A spot coincides with the transmitted beam; the beam is the calibration "
                "centre, not a reflection."
            )

    def __len__(self) -> int:
        return len(self.spots)

    def g_vectors_inv_angstrom(self) -> np.ndarray:
        """In-plane ``(n, 2)`` reciprocal vectors, beam-centred, in 1/angstrom."""

        return self.calibration.to_reciprocal_angstrom(
            [spot.position for spot in self.spots]
        )

    def g_magnitudes_inv_angstrom(self) -> np.ndarray:
        return np.asarray(
            np.linalg.norm(self.g_vectors_inv_angstrom(), axis=1), dtype=np.float64
        )

    def d_spacings_angstrom(self) -> np.ndarray:
        """Observed ``d`` per spot, ``d = 1/|g|`` — the quantity to compare with tables."""

        return 1.0 / self.g_magnitudes_inv_angstrom()

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> MeasuredSAEDPattern:
        """Build from a decoded YAML/JSON mapping following the documented schema."""

        schema = payload.get("schema")
        if schema is not None and schema != MEASURED_PATTERN_SCHEMA:
            raise ValueError(
                f"Unsupported measured-pattern schema '{schema}'; expected "
                f"'{MEASURED_PATTERN_SCHEMA}'."
            )
        calibration_payload = dict(payload.get("calibration") or {})
        centre = calibration_payload.pop("centre", (0.0, 0.0))
        calibration = PatternCalibration(
            centre=(float(centre[0]), float(centre[1])), **calibration_payload
        )
        spots = tuple(
            MeasuredSpot(
                position=(float(entry["x"]), float(entry["y"])),
                intensity=(
                    None if entry.get("intensity") is None else float(entry["intensity"])
                ),
                label=entry.get("label"),
            )
            for entry in payload["spots"]
        )
        return cls(
            name=str(payload.get("name", "measured_pattern")),
            spots=spots,
            calibration=calibration,
        )

    @classmethod
    def from_yaml(cls, path: str | Path) -> MeasuredSAEDPattern:
        """Read a measured pattern from its YAML file.

        The file is the reproducibility boundary: a pattern solved from a
        committed YAML gives the same answer on any machine, whether the spots
        were originally clicked or typed.
        """

        decoded = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if not isinstance(decoded, dict):
            raise ValueError("A measured-pattern YAML file must decode to a mapping.")
        return cls.from_dict(decoded)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": MEASURED_PATTERN_SCHEMA,
            "name": self.name,
            "calibration": self.calibration.to_dict(),
            "spots": [
                {
                    "x": spot.position[0],
                    "y": spot.position[1],
                    **({} if spot.intensity is None else {"intensity": spot.intensity}),
                    **({} if spot.label is None else {"label": spot.label}),
                }
                for spot in self.spots
            ],
        }

    def to_yaml(self, path: str | Path) -> Path:
        """Write the pattern to a YAML file and return the path."""

        output = Path(path)
        output.write_text(
            yaml.safe_dump(self.to_dict(), sort_keys=False), encoding="utf-8"
        )
        return output

    def describe(self) -> str:
        """Prose summary: how many spots, their calibration, and their d-range."""

        spacings = self.d_spacings_angstrom()
        return (
            f"Measured SAED pattern '{self.name}': {len(self)} spot(s). "
            f"{self.calibration.describe()} Observed d-spacings range from "
            f"{float(np.min(spacings)):.4f} to {float(np.max(spacings)):.4f} angstrom "
            f"(|g| from {1.0 / float(np.max(spacings)):.4f} to "
            f"{1.0 / float(np.min(spacings)):.4f} 1/angstrom)."
        )


@dataclass(frozen=True, slots=True)
class SolvedSpot:
    """One measured spot with the reflection assigned to it."""

    measured_index: int
    hkl: tuple[int, int, int]
    label: str
    predicted_g_inv_angstrom: tuple[float, float]
    residual_inv_angstrom: float

    def __post_init__(self) -> None:
        if self.measured_index < 0:
            raise ValueError("measured_index must be non-negative.")
        if not np.isfinite(self.residual_inv_angstrom) or self.residual_inv_angstrom < 0.0:
            raise ValueError("residual_inv_angstrom must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class PatternSolution:
    """One candidate answer: a phase, a zone axis, and an indexing of the spots.

    ``orientation`` maps crystal Cartesian vectors into the pattern frame, whose
    ``x`` and ``y`` axes are the picked coordinates' axes and whose ``z`` axis
    points along the zone axis toward the viewer.

    A single SAED pattern cannot distinguish a zone axis from its reverse when
    the reflection set is centrosymmetric — inverting the crystal through the
    origin leaves the pattern unchanged. `PatternSolutionReport` reports that
    explicitly rather than presenting one of the two as the answer.
    """

    phase_name: str
    zone_axis: ZoneAxis
    zone_axis_label: str
    orientation: Rotation
    solved_spots: tuple[SolvedSpot, ...]
    unindexed_spot_indices: tuple[int, ...]
    measured_spot_count: int
    seed_spot_indices: tuple[int, int]
    variant_index: int | None = None
    variant_deviation_deg: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "solved_spots", tuple(self.solved_spots))
        object.__setattr__(
            self, "unindexed_spot_indices", tuple(self.unindexed_spot_indices)
        )
        if self.measured_spot_count <= 0:
            raise ValueError("measured_spot_count must be positive.")
        if len(self.solved_spots) + len(self.unindexed_spot_indices) != (
            self.measured_spot_count
        ):
            raise ValueError(
                "Every measured spot must be either solved or listed as unindexed."
            )

    @property
    def matched_fraction(self) -> float:
        return len(self.solved_spots) / self.measured_spot_count

    @property
    def mean_residual_inv_angstrom(self) -> float:
        if not self.solved_spots:
            return float("inf")
        return float(
            np.mean([spot.residual_inv_angstrom for spot in self.solved_spots])
        )

    @property
    def max_residual_inv_angstrom(self) -> float:
        if not self.solved_spots:
            return float("inf")
        return float(np.max([spot.residual_inv_angstrom for spot in self.solved_spots]))

    @property
    def score(self) -> tuple[float, float]:
        """Ranking key: more spots indexed first, then smaller residuals.

        Matched fraction dominates deliberately. A solution that explains every
        spot with moderate residuals is a better answer than one that explains
        half of them perfectly, which is usually a coincidence on a sub-lattice.
        """

        return (-self.matched_fraction, self.mean_residual_inv_angstrom)

    def describe(self) -> str:
        """Prose summary: phase, zone, indexing quality, and the assignments."""

        lines = [
            f"Phase '{self.phase_name}' viewed along zone axis {self.zone_axis_label}: "
            f"{len(self.solved_spots)} of {self.measured_spot_count} measured spot(s) "
            f"indexed ({100.0 * self.matched_fraction:.1f}%), mean residual "
            f"{self.mean_residual_inv_angstrom:.5f} 1/angstrom, max "
            f"{self.max_residual_inv_angstrom:.5f} 1/angstrom. Seeded from measured "
            f"spots {self.seed_spot_indices[0]} and {self.seed_spot_indices[1]}."
        ]
        if self.variant_index is not None:
            lines.append(
                f"Best-matching transformation variant: {self.variant_index} "
                f"(deviation {self.variant_deviation_deg:.3f} deg)."
            )
        assignments = ", ".join(
            f"{spot.label}"
            for spot in sorted(self.solved_spots, key=lambda item: item.measured_index)[:12]
        )
        lines.append(f"Assignments: {assignments}.")
        if self.unindexed_spot_indices:
            lines.append(
                f"Unindexed spot(s): {list(self.unindexed_spot_indices)} — not explained "
                "by this phase and zone at the given tolerances. A spot whose indices "
                "exceed the solver's max_index is simply never offered a match, so raise "
                "it before widening the tolerances."
            )
        return " ".join(lines)

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "phase": self.phase_name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "zone_axis_label": self.zone_axis_label,
            "orientation_quaternion": [
                float(value) for value in self.orientation.quaternion
            ],
            "matched_fraction": self.matched_fraction,
            "mean_residual_inv_angstrom": self.mean_residual_inv_angstrom,
            "max_residual_inv_angstrom": self.max_residual_inv_angstrom,
            "seed_spot_indices": list(self.seed_spot_indices),
            "variant_index": self.variant_index,
            "variant_deviation_deg": self.variant_deviation_deg,
            "spots": [
                {
                    "measured_index": spot.measured_index,
                    "hkl": list(spot.hkl),
                    "label": spot.label,
                    "predicted_g_inv_angstrom": list(spot.predicted_g_inv_angstrom),
                    "residual_inv_angstrom": spot.residual_inv_angstrom,
                }
                for spot in self.solved_spots
            ],
            "unindexed_spot_indices": list(self.unindexed_spot_indices),
        }


@dataclass(frozen=True, slots=True)
class PatternSolutionReport:
    """Ranked solutions for one measured pattern, with an honest verdict.

    ``is_conclusive`` requires the best solution to index every spot and to lead
    the runner-up in matched fraction, *or* to be the only solution found. Two
    solutions that differ only by the zone-axis sense are not counted against
    each other — that ambiguity is intrinsic to a single pattern, not a failure
    to discriminate, and `describe()` names it.
    """

    pattern_name: str
    solutions: tuple[PatternSolution, ...]
    considered_phase_names: tuple[str, ...]
    measured_spot_count: int
    length_tolerance_relative: float
    angle_tolerance_deg: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "solutions", tuple(self.solutions))
        object.__setattr__(
            self, "considered_phase_names", tuple(self.considered_phase_names)
        )
        if not self.considered_phase_names:
            raise ValueError("At least one candidate phase must have been considered.")

    def __len__(self) -> int:
        return len(self.solutions)

    def best(self) -> PatternSolution:
        """The highest-ranked solution, or raise if the pattern was not solved."""

        if not self.solutions:
            raise ValueError(
                f"No solution was found for pattern '{self.pattern_name}'. Widen the "
                "tolerances, raise max_index, or add the correct candidate phase."
            )
        return self.solutions[0]

    @property
    def is_conclusive(self) -> bool:
        if not self.solutions:
            return False
        best = self.solutions[0]
        if best.matched_fraction < 1.0:
            return False
        distinct = [
            solution
            for solution in self.solutions[1:]
            if solution.phase_name != best.phase_name
            or not _same_zone_family(solution.zone_axis, best.zone_axis)
        ]
        return not distinct or distinct[0].matched_fraction < 1.0

    def describe(self) -> str:
        """Prose summary: the ranking, the verdict, and the intrinsic ambiguity."""

        phases = ", ".join(self.considered_phase_names)
        if not self.solutions:
            return (
                f"Measured SAED pattern '{self.pattern_name}' ({self.measured_spot_count} "
                f"spot(s)) could not be solved against {phases} at a "
                f"{100.0 * self.length_tolerance_relative:.1f}% length tolerance and "
                f"{self.angle_tolerance_deg:g} deg angle tolerance. No candidate phase "
                "supplies two reflections matching the two seed spots in both length and "
                "interplanar angle."
            )
        best = self.solutions[0]
        verdict = (
            "The solution is unambiguous among the candidates."
            if self.is_conclusive
            else (
                "The solution is NOT unambiguous: another candidate explains the pattern "
                "equally well. Tilt to a second zone axis, or add spots at larger |g|, "
                "before trusting it."
            )
        )
        ranking = "; ".join(
            f"{solution.phase_name} {solution.zone_axis_label} "
            f"({100.0 * solution.matched_fraction:.0f}%, "
            f"{solution.mean_residual_inv_angstrom:.4f} 1/A)"
            for solution in self.solutions[:5]
        )
        return (
            f"Measured SAED pattern '{self.pattern_name}' ({self.measured_spot_count} "
            f"spot(s)) solved against {phases}: best is {best.describe()} {verdict} "
            f"Ranking: {ranking}. Note that a single SAED pattern cannot distinguish a "
            "zone axis from its reverse when the reflection set is centrosymmetric, "
            "because inverting the crystal leaves the pattern unchanged; the reported "
            "sense is one of the two equally valid descriptions."
        )

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "schema": PATTERN_SOLUTION_SCHEMA,
            "pattern": self.pattern_name,
            "considered_phases": list(self.considered_phase_names),
            "measured_spot_count": self.measured_spot_count,
            "length_tolerance_relative": self.length_tolerance_relative,
            "angle_tolerance_deg": self.angle_tolerance_deg,
            "is_conclusive": self.is_conclusive,
            "solutions": [solution.to_json_dict() for solution in self.solutions],
        }


def _same_zone_family(left: ZoneAxis, right: ZoneAxis) -> bool:
    """Whether two zone axes are the same axis up to sign.

    Sign is deliberately ignored: the two senses of one zone axis are the
    intrinsic single-pattern ambiguity, not two competing answers.
    """

    a = np.asarray(left.indices, dtype=np.int64)
    b = np.asarray(right.indices, dtype=np.int64)
    return bool(np.array_equal(a, b) or np.array_equal(a, -b))


def _allowed_reflections(phase: Phase, max_index: int) -> tuple[np.ndarray, np.ndarray]:
    """Allowed ``hkl`` and their Cartesian ``g`` vectors, up to ``max_index``.

    Reuses the same centering test the simulation engine applies, so a
    reflection the simulator would not draw is never offered as an index here.
    """

    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1)
    hkl = grid.reshape(-1, 3)
    hkl = hkl[np.any(hkl != 0, axis=1)]
    hkl = hkl[centering_allowed_mask(hkl, ReflectionCondition.from_phase(phase))]
    reciprocal = phase.lattice.reciprocal_basis().matrix
    return hkl, hkl.astype(np.float64) @ reciprocal.T


def _seed_pairs(g_observed: np.ndarray, magnitudes: np.ndarray) -> list[tuple[int, int]]:
    """Index pairs of measured spots that can seed a solution.

    The shortest vectors are preferred because they are the best-determined
    relative to picking error, and pairs that are (anti)parallel are skipped
    because two collinear reflections do not define a plane.
    """

    order = np.argsort(magnitudes)
    pairs: list[tuple[int, int]] = []
    for first in range(len(order)):
        for second in range(first + 1, len(order)):
            a, b = int(order[first]), int(order[second])
            cross = abs(
                g_observed[a, 0] * g_observed[b, 1] - g_observed[a, 1] * g_observed[b, 0]
            )
            if cross <= 1e-9 * magnitudes[a] * magnitudes[b]:
                continue
            pairs.append((a, b))
            if len(pairs) >= _MAX_SEED_PAIRS:
                return pairs
    return pairs


def _rotation_from_pair(
    g1_crystal: np.ndarray,
    g2_crystal: np.ndarray,
    g1_observed: np.ndarray,
    g2_observed: np.ndarray,
) -> np.ndarray | None:
    """Crystal-to-pattern rotation carrying a calculated pair onto an observed pair.

    Both pairs are turned into right-handed orthonormal triads by Gram-Schmidt —
    first vector, then the in-plane part of the second, then their cross product
    — and the rotation is the product of the two triads. Returns ``None`` when
    either pair is degenerate, which the caller treats as an inadmissible
    assignment rather than an error.
    """

    def _triad(first: np.ndarray, second: np.ndarray) -> np.ndarray | None:
        e1 = first / np.linalg.norm(first)
        residual = second - float(np.dot(second, e1)) * e1
        norm = float(np.linalg.norm(residual))
        if norm < 1e-12:
            return None
        e2 = residual / norm
        return np.column_stack([e1, e2, np.cross(e1, e2)])

    crystal = _triad(g1_crystal, g2_crystal)
    observed = _triad(
        np.array([g1_observed[0], g1_observed[1], 0.0]),
        np.array([g2_observed[0], g2_observed[1], 0.0]),
    )
    if crystal is None or observed is None:
        return None
    return np.asarray(observed @ crystal.T, dtype=np.float64)


def _zone_axis_from_rotation(
    rotation: np.ndarray, phase: Phase, max_index: int
) -> tuple[ZoneAxis, float] | None:
    """The direct-lattice zone axis whose Cartesian image is the pattern normal."""

    normal_crystal = rotation.T @ np.array([0.0, 0.0, 1.0])
    direction = CrystalDirection.from_cartesian(normal_crystal, phase=phase)
    rationalized = rationalize_zone_axis(direction, max_index=max_index)
    return (
        ZoneAxis(rationalized.indices, phase=phase),
        float(rationalized.deviation_deg),
    )


def solve_saed_pattern(
    pattern: MeasuredSAEDPattern,
    phases: Sequence[Phase],
    *,
    max_index: int = 4,
    length_tolerance_relative: float = DEFAULT_LENGTH_TOLERANCE_RELATIVE,
    angle_tolerance_deg: float = DEFAULT_ANGLE_TOLERANCE_DEG,
    max_solutions: int = 5,
    provenance: ProvenanceRecord | None = None,
) -> PatternSolutionReport:
    """Solve a measured SAED pattern: which phase, which zone axis, which spots.

    Purpose: answers the question asked in front of a microscope. Given picked
    spot positions and enough calibration to scale them, it determines the
    phase, the zone axis, the crystal orientation in the pattern frame, and the
    Miller indices of every spot — with residuals, alternatives, and an explicit
    verdict on whether the answer is unambiguous.

    When to use: interpreting a measured or printed diffraction pattern, or
    checking a simulated one. Use
    `pytex.diffraction.models.index_saed_pattern` instead when a calibrated
    `DiffractionGeometry` is available: that path knows the detector, this one
    only needs the pattern.

    Algorithm: classical ratio/angle indexing. The two shortest non-collinear
    measured vectors seed the solution; a calculated pair is admissible when
    both lengths match within ``length_tolerance_relative`` and their
    interplanar angle within ``angle_tolerance_deg``. The zone axis follows from
    the pair, the crystal-to-pattern rotation from aligning calculated onto
    observed, and every remaining spot is indexed by projection. Intensities are
    never used — geometry alone decides.

    Inputs: the pattern; the candidate ``phases``; ``max_index``, the reflection
    index bound enumerated per phase (systematic absences from each phase's
    space group are applied, so an undeclared space group means a primitive
    lattice is assumed); the two tolerances; and how many ranked solutions to
    return.

    Output: a `PatternSolutionReport` — read its ``describe()``. It may hold no
    solutions, which is a legitimate answer meaning no candidate phase explains
    the pattern at these tolerances.

    A single SAED pattern cannot distinguish a zone axis from its reverse for a
    centrosymmetric reflection set. The report says so and does not count the two
    senses as competing answers.

    See also
    --------
    `MeasuredSAEDPattern.from_yaml` : the file contract for picked spots.
    `solve_saed_pattern_file` : read and solve in one call.
    """

    if not phases:
        raise ValueError("solve_saed_pattern requires at least one candidate phase.")
    if max_index < 1:
        raise ValueError("max_index must be at least 1.")
    if not 0.0 < length_tolerance_relative < 1.0:
        raise ValueError("length_tolerance_relative must lie in (0, 1).")
    if not 0.0 < angle_tolerance_deg < 180.0:
        raise ValueError("angle_tolerance_deg must lie in (0, 180).")

    g_observed = pattern.g_vectors_inv_angstrom()
    magnitudes = pattern.g_magnitudes_inv_angstrom()
    match_radius = length_tolerance_relative * magnitudes
    seeds = _seed_pairs(g_observed, magnitudes)

    solutions: list[PatternSolution] = []
    phase_by_name = {phase.name: phase for phase in phases}
    reflections: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for phase in phases:
        hkl, g_crystal = _allowed_reflections(phase, max_index)
        reflections[phase.name] = (hkl, g_crystal)
        if hkl.shape[0] == 0:
            continue
        calculated = np.linalg.norm(g_crystal, axis=1)
        units = g_crystal / calculated[:, None]
        for first, second in seeds:
            target_a, target_b = magnitudes[first], magnitudes[second]
            pool_a = np.flatnonzero(
                np.abs(calculated - target_a) <= length_tolerance_relative * target_a
            )
            pool_b = np.flatnonzero(
                np.abs(calculated - target_b) <= length_tolerance_relative * target_b
            )
            if pool_a.size == 0 or pool_b.size == 0:
                continue
            observed_cos = float(
                np.dot(g_observed[first], g_observed[second]) / (target_a * target_b)
            )
            cosines = units[pool_a] @ units[pool_b].T
            observed_angle = float(np.arccos(np.clip(observed_cos, -1.0, 1.0)))
            calculated_angles = np.arccos(np.clip(cosines, -1.0, 1.0))
            admissible = np.argwhere(
                np.abs(calculated_angles - observed_angle)
                <= np.radians(angle_tolerance_deg)
            )
            for row in admissible[:_MAX_ASSIGNMENTS_PER_SEED]:
                index_a = int(pool_a[int(row[0])])
                index_b = int(pool_b[int(row[1])])
                rotation = _rotation_from_pair(
                    g_crystal[index_a],
                    g_crystal[index_b],
                    g_observed[first],
                    g_observed[second],
                )
                if rotation is None:
                    continue
                zone = _zone_axis_from_rotation(rotation, phase, max_index)
                if zone is None:
                    continue
                zone_axis, _ = zone
                solved, unindexed = _assign_spots(
                    rotation=rotation,
                    hkl=hkl,
                    g_crystal=g_crystal,
                    g_observed=g_observed,
                    match_radius=match_radius,
                    phase=phase,
                )
                if not solved:
                    continue
                solutions.append(
                    PatternSolution(
                        phase_name=phase.name,
                        zone_axis=zone_axis,
                        zone_axis_label=format_direction_indices(
                            tuple(int(value) for value in zone_axis.indices), style="plain"
                        ),
                        orientation=Rotation.from_matrix(rotation).canonicalized(),
                        solved_spots=tuple(solved),
                        unindexed_spot_indices=tuple(unindexed),
                        measured_spot_count=len(pattern),
                        seed_spot_indices=(first, second),
                    )
                )
    ranked = _deduplicate_solutions(
        sorted(
            solutions,
            key=lambda item: (item.score, _zone_axis_preference(item.zone_axis)),
        ),
        phases,
        # One extra so `is_conclusive` can see whether a genuine competitor
        # exists even when the caller only wants the top few reported.
        limit=max_solutions + 1,
    )
    ranked = [
        _canonicalize_description(
            solution,
            phase=phase_by_name[solution.phase_name],
            hkl=reflections[solution.phase_name][0],
            g_crystal=reflections[solution.phase_name][1],
            g_observed=g_observed,
            match_radius=match_radius,
            max_index=max_index,
        )
        for solution in ranked
    ]
    return PatternSolutionReport(
        pattern_name=pattern.name,
        solutions=tuple(ranked[:max_solutions]),
        considered_phase_names=tuple(phase.name for phase in phases),
        measured_spot_count=len(pattern),
        length_tolerance_relative=float(length_tolerance_relative),
        angle_tolerance_deg=float(angle_tolerance_deg),
        provenance=provenance,
    )


def _assign_spots(
    *,
    rotation: np.ndarray,
    hkl: np.ndarray,
    g_crystal: np.ndarray,
    g_observed: np.ndarray,
    match_radius: np.ndarray,
    phase: Phase,
) -> tuple[list[SolvedSpot], list[int]]:
    """Index every measured spot against the reflections this rotation predicts."""

    projected = g_crystal @ rotation.T
    tolerance = float(np.max(match_radius))
    in_zone = np.abs(projected[:, 2]) <= tolerance
    plane = projected[in_zone][:, :2]
    candidates = hkl[in_zone]
    bravais = is_hexagonal_phase(phase)
    solved: list[SolvedSpot] = []
    unindexed: list[int] = []
    used: set[int] = set()
    for index in range(g_observed.shape[0]):
        if plane.shape[0] == 0:
            unindexed.append(index)
            continue
        distances = np.linalg.norm(plane - g_observed[index], axis=1)
        order = np.argsort(distances)
        chosen: int | None = None
        for candidate in order:
            position = int(candidate)
            if position in used:
                continue
            if distances[position] <= match_radius[index]:
                chosen = position
            break
        if chosen is None:
            unindexed.append(index)
            continue
        used.add(chosen)
        indices = tuple(int(value) for value in candidates[chosen])
        label_indices = (
            tuple(int(value) for value in plane_hkl_to_hkil(candidates[chosen]))
            if bravais
            else indices
        )
        solved.append(
            SolvedSpot(
                measured_index=index,
                hkl=indices,  # type: ignore[arg-type]
                label=format_plane_indices(label_indices, style="plain"),
                predicted_g_inv_angstrom=(
                    float(plane[chosen, 0]),
                    float(plane[chosen, 1]),
                ),
                residual_inv_angstrom=float(distances[chosen]),
            )
        )
    return solved, unindexed


#: How close two crystal orientations must be to count as the same solution.
#:
#: Generous on purpose. Distinct indexings of one pattern are separated by whole
#: symmetry operations — tens of degrees — while the same solution reached
#: through two different seed pairs agrees to well under a degree. Anything in
#: between would mean the seed assignments themselves disagreed, which the
#: length and angle tolerances already exclude.
_SOLUTION_EQUIVALENCE_TOLERANCE_DEG = 1.0


def _zone_axis_preference(zone_axis: ZoneAxis) -> tuple[int, tuple[int, ...]]:
    """Sort key preferring the conventional member of an equivalent zone family.

    Symmetry-equivalent solutions are the same answer written differently, so
    the one reported should be the one a crystallographer would write: fewest
    negative components, then the lowest indices. For a cubic pattern down a
    cube axis this reports ``[001]`` rather than the equally valid ``[0-10]``.
    """

    indices = tuple(int(value) for value in zone_axis.indices)
    return (sum(1 for value in indices if value < 0), indices)


def _orientations_equivalent(
    left: np.ndarray, right: np.ndarray, operators: np.ndarray
) -> bool:
    """Whether two crystal-to-pattern rotations describe the same orientation.

    They do when they differ by a crystal symmetry operation, which is exactly
    what makes two seed assignments that look different produce one physical
    answer.
    """

    candidates = np.einsum("ij,sjk->sik", right, operators, optimize=True)
    relative = np.einsum("sij,kj->sik", candidates, left, optimize=True)
    traces = np.trace(relative, axis1=-2, axis2=-1)
    best = float(np.max(traces))
    angle = float(np.degrees(np.arccos(np.clip((best - 1.0) * 0.5, -1.0, 1.0))))
    return angle <= _SOLUTION_EQUIVALENCE_TOLERANCE_DEG


def _canonicalize_description(
    solution: PatternSolution,
    *,
    phase: Phase,
    hkl: np.ndarray,
    g_crystal: np.ndarray,
    g_observed: np.ndarray,
    match_radius: np.ndarray,
    max_index: int,
) -> PatternSolution:
    """Rewrite a solution using the most conventional equivalent description.

    A crystal symmetry operation relabels the crystal without changing anything
    physical, so a cubic pattern down a cube axis is equally correctly reported
    as ``[001]``, ``[010]`` or ``[0-10]``. Which one comes out of the solver
    depends on which seed assignment happened to be tried first, which is not a
    fact about the specimen. This picks the description a crystallographer would
    write — fewest negative indices, then lowest — and reindexes every spot in
    it, so the reported zone axis and the reported spot indices stay consistent.
    """

    operators = (
        phase.symmetry.operators
        if phase.symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    base = solution.orientation.as_matrix()
    best = solution
    best_key = _zone_axis_preference(solution.zone_axis)
    for operator in operators:
        rotation = np.ascontiguousarray(base @ operator)
        zone = _zone_axis_from_rotation(rotation, phase, max_index)
        if zone is None:
            continue
        zone_axis, deviation = zone
        if deviation > _SOLUTION_EQUIVALENCE_TOLERANCE_DEG:
            continue
        key = _zone_axis_preference(zone_axis)
        if key >= best_key:
            continue
        solved, unindexed = _assign_spots(
            rotation=rotation,
            hkl=hkl,
            g_crystal=g_crystal,
            g_observed=g_observed,
            match_radius=match_radius,
            phase=phase,
        )
        # A relabelling must not change how much of the pattern is explained.
        # If it does, the operator is not a symmetry of this solution and the
        # original description stands.
        if len(solved) != len(solution.solved_spots):
            continue
        best_key = key
        best = PatternSolution(
            phase_name=solution.phase_name,
            zone_axis=zone_axis,
            zone_axis_label=format_direction_indices(
                tuple(int(value) for value in zone_axis.indices), style="plain"
            ),
            orientation=Rotation.from_matrix(rotation).canonicalized(),
            solved_spots=tuple(solved),
            unindexed_spot_indices=tuple(unindexed),
            measured_spot_count=solution.measured_spot_count,
            seed_spot_indices=solution.seed_spot_indices,
            variant_index=solution.variant_index,
            variant_deviation_deg=solution.variant_deviation_deg,
        )
    return best


def _deduplicate_solutions(
    solutions: Sequence[PatternSolution],
    phases: Sequence[Phase],
    *,
    limit: int,
) -> list[PatternSolution]:
    """Collapse solutions that describe the same crystal orientation.

    Many seed assignments are related by a crystal symmetry operation and give
    the same physical answer through different bookkeeping; reporting them all
    would make an unambiguous solve look contested. Solutions arrive already
    sorted by score and then by zone-axis conventionality, so the survivor of
    each group is the one worth showing.

    Comparison stops once ``limit`` distinct solutions are held, which bounds
    the cost: each candidate is compared only against the few kept, not against
    every other candidate.
    """

    operators = {
        phase.name: (
            phase.symmetry.operators
            if phase.symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        for phase in phases
    }
    unique: list[PatternSolution] = []
    matrices: list[np.ndarray] = []
    for solution in solutions:
        matrix = solution.orientation.as_matrix()
        group = operators[solution.phase_name]
        if any(
            kept.phase_name == solution.phase_name
            and _orientations_equivalent(existing, matrix, group)
            for kept, existing in zip(unique, matrices, strict=True)
        ):
            continue
        unique.append(solution)
        matrices.append(matrix)
        if len(unique) >= limit:
            break
    return unique


def assign_transformation_variant(
    solution: PatternSolution,
    relationship: OrientationRelationship,
    parent_orientation: Rotation,
) -> PatternSolution:
    """Name which transformation variant a solved child pattern belongs to.

    Purpose: completes the composite story. Once a product-phase pattern has
    been solved, and the parent's orientation in the same pattern frame is
    known — from solving the parent's own spots in that pattern, or from EBSD —
    the variant follows: the child orientation a variant predicts is
    ``P V_k^T`` in the canonical crystal-to-specimen convention, and the variant
    whose prediction sits closest to the solved orientation is the answer.

    Inputs: a solved child-phase `PatternSolution`; the relationship, whose
    child phase must be the solution's phase; and the parent's crystal-to-pattern
    rotation.

    Output: a copy of the solution carrying ``variant_index`` and
    ``variant_deviation_deg``. The deviation is symmetry-reduced under the child
    point group, so it is the true disorientation from the prediction; a large
    value means the pattern does not belong to this relationship at all, which
    the caller should check rather than assume.

    See also
    --------
    `pytex.core.transformation.OrientationRelationship.generate_variants`.
    """

    if relationship.child_phase.name != solution.phase_name:
        raise ValueError(
            f"The solution is for phase '{solution.phase_name}' but the relationship's "
            f"child phase is '{relationship.child_phase.name}'."
        )
    variants = relationship.generate_variants()
    parent_matrix = parent_orientation.as_matrix()
    solved = solution.orientation.as_matrix()
    child_symmetry = relationship.child_phase.symmetry
    operators = (
        child_symmetry.operators
        if child_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    best_index = 0
    best_angle = float("inf")
    for variant in variants:
        # Canonical convention: C = P V^T.
        predicted = parent_matrix @ variant.parent_to_child_rotation.as_matrix().T
        candidates = np.einsum("sij,jk->sik", operators, predicted, optimize=True)
        relative = np.einsum("sij,kj->sik", candidates, solved, optimize=True)
        traces = np.trace(relative, axis1=-2, axis2=-1)
        angle = float(
            np.degrees(np.arccos(np.clip((float(np.max(traces)) - 1.0) * 0.5, -1.0, 1.0)))
        )
        if angle < best_angle:
            best_angle = angle
            best_index = variant.variant_index
    return PatternSolution(
        phase_name=solution.phase_name,
        zone_axis=solution.zone_axis,
        zone_axis_label=solution.zone_axis_label,
        orientation=solution.orientation,
        solved_spots=solution.solved_spots,
        unindexed_spot_indices=solution.unindexed_spot_indices,
        measured_spot_count=solution.measured_spot_count,
        seed_spot_indices=solution.seed_spot_indices,
        variant_index=best_index,
        variant_deviation_deg=best_angle,
    )


def solve_saed_pattern_file(
    path: str | Path,
    phases: Sequence[Phase],
    **kwargs: Any,
) -> PatternSolutionReport:
    """Read a measured-pattern YAML file and solve it.

    The convenience form of `solve_saed_pattern` for the common case where the
    spots live in a file. All keyword arguments are forwarded.
    """

    return solve_saed_pattern(MeasuredSAEDPattern.from_yaml(path), phases, **kwargs)

__all__ = [
    "COORDINATE_UNITS",
    "DEFAULT_ANGLE_TOLERANCE_DEG",
    "DEFAULT_LENGTH_TOLERANCE_RELATIVE",
    "MEASURED_PATTERN_SCHEMA",
    "PATTERN_SOLUTION_SCHEMA",
    "MeasuredSAEDPattern",
    "MeasuredSpot",
    "PatternCalibration",
    "PatternSolution",
    "PatternSolutionReport",
    "SolvedSpot",
    "assign_transformation_variant",
    "solve_saed_pattern",
    "solve_saed_pattern_file",
]
