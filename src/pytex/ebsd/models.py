from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import normalize_vector
from pytex.core.acquisition import AcquisitionGeometry, CalibrationRecord, MeasurementQuality
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, MillerIndex, Phase
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec

if TYPE_CHECKING:
    from pytex.adapters import ExperimentManifest
    from pytex.texture import ODF, InversePoleFigure, KernelSpec, PoleFigure


def _specimen_direction_vector(
    direction: str | ArrayLike,
    specimen_frame: ReferenceFrame,
) -> np.ndarray:
    if isinstance(direction, str):
        normalized = direction.strip().lower()
        axis_lookup = {label.lower(): index for index, label in enumerate(specimen_frame.axes)}
        axis_lookup.update({"x": 0, "y": 1, "z": 2})
        if normalized not in axis_lookup:
            raise ValueError(
                "Sample direction labels must be one of the specimen-frame axis labels or "
                "'x', 'y', 'z'."
            )
        vector = np.zeros(3, dtype=np.float64)
        vector[axis_lookup[normalized]] = 1.0
        return vector
    return normalize_vector(direction)


def _coerce_pole(
    pole: CrystalPlane | ArrayLike,
    *,
    phase: Phase | None,
) -> CrystalPlane:
    if isinstance(pole, CrystalPlane):
        return pole
    if phase is None:
        raise ValueError(
            "CrystalMap.pole_figure() requires OrientationSet.phase when poles are passed as "
            "raw Miller indices."
        )
    indices = np.asarray(pole, dtype=np.int64)
    if indices.shape != (3,):
        raise ValueError("Raw pole indices must have shape (3,).")
    return CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase)


def _coerce_pole_sequence(
    poles: (
        CrystalPlane
        | ArrayLike
        | tuple[CrystalPlane | ArrayLike, ...]
        | list[CrystalPlane | ArrayLike]
    ),
    *,
    phase: Phase | None,
) -> tuple[CrystalPlane, ...]:
    if isinstance(poles, CrystalPlane):
        return (poles,)
    if isinstance(poles, np.ndarray) and poles.shape == (3,):
        return (_coerce_pole(poles, phase=phase),)
    if isinstance(poles, (list, tuple)):
        return tuple(_coerce_pole(pole, phase=phase) for pole in poles)
    return (_coerce_pole(poles, phase=phase),)


def _coerce_sample_direction_sequence(
    sample_directions: str
    | ArrayLike
    | tuple[str | ArrayLike, ...]
    | list[str | ArrayLike],
    specimen_frame: ReferenceFrame,
) -> tuple[np.ndarray, ...]:
    if isinstance(sample_directions, str):
        return (_specimen_direction_vector(sample_directions, specimen_frame),)
    if isinstance(sample_directions, np.ndarray) and sample_directions.shape == (3,):
        return (_specimen_direction_vector(sample_directions, specimen_frame),)
    if isinstance(sample_directions, (list, tuple)):
        if len(sample_directions) == 3 and not isinstance(sample_directions[0], str):
            candidate = np.asarray(sample_directions, dtype=np.float64)
            if candidate.shape == (3,):
                return (_specimen_direction_vector(candidate, specimen_frame),)
        return tuple(
            _specimen_direction_vector(direction, specimen_frame)
            for direction in sample_directions
        )
    return (_specimen_direction_vector(sample_directions, specimen_frame),)


def _canonical_phase_name(name: str) -> str:
    normalized = name.strip()
    if not normalized:
        raise ValueError("Phase names must be non-empty strings.")
    return normalized


def _readonly_int_array(
    array: ArrayLike,
    *,
    shape: tuple[int | None, ...],
    name: str,
) -> np.ndarray:
    values = np.asarray(array, dtype=np.int64)
    if values.ndim != len(shape):
        raise ValueError(f"{name} must have shape {shape}.")
    for axis, expected in enumerate(shape):
        if expected is not None and values.shape[axis] != expected:
            raise ValueError(f"{name} must have shape {shape}.")
    values = np.ascontiguousarray(values)
    values.setflags(write=False)
    return values


def _readonly_float_array(
    array: ArrayLike, *, shape: tuple[int | None, ...], name: str
) -> np.ndarray:
    values = np.asarray(array, dtype=np.float64)
    if values.ndim != len(shape):
        raise ValueError(f"{name} must have shape {shape}.")
    for axis, expected in enumerate(shape):
        if expected is not None and values.shape[axis] != expected:
            raise ValueError(f"{name} must have shape {shape}.")
    values = np.ascontiguousarray(values)
    values.setflags(write=False)
    return values


def _rotation_angles_from_matrices(matrices: np.ndarray) -> np.ndarray:
    traces = np.trace(matrices, axis1=1, axis2=2)
    cos_theta = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    return np.asarray(np.arccos(cos_theta), dtype=np.float64)


def _relative_rotation_matrices(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.asarray(
        np.einsum("nij,nkj->nik", right, left, optimize=True),
        dtype=np.float64,
    )


def _disorientation_angles_from_relative_matrices(
    relative_matrices: np.ndarray,
    *,
    left_symmetry: SymmetrySpec | None,
    right_symmetry: SymmetrySpec | None,
) -> np.ndarray:
    left_ops = (
        left_symmetry.operators
        if left_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    right_ops = (
        right_symmetry.operators
        if right_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    left_applied = np.einsum("aij,njk->naik", left_ops, relative_matrices, optimize=True)
    candidates = np.einsum(
        "naij,bkj->nabik",
        left_applied,
        right_ops,
        optimize=True,
    )
    candidate_matrices = candidates.reshape(-1, 3, 3)
    angles = _rotation_angles_from_matrices(candidate_matrices).reshape(
        relative_matrices.shape[0],
        left_ops.shape[0] * right_ops.shape[0],
    )
    return np.asarray(np.min(angles, axis=1), dtype=np.float64)


def _pairwise_distances(coordinates: np.ndarray) -> np.ndarray:
    deltas = coordinates[:, None, :] - coordinates[None, :, :]
    return np.asarray(np.linalg.norm(deltas, axis=2), dtype=np.float64)


def _inferred_base_spacing(coordinates: np.ndarray, step_sizes: tuple[float, ...] | None) -> float:
    if step_sizes is not None:
        return float(min(step_sizes))
    distances = _pairwise_distances(coordinates)
    positive = distances[distances > 1e-12]
    if positive.size == 0:
        raise ValueError(
            "CrystalMap requires at least two distinct coordinates for graph workflows."
        )
    return float(np.min(positive))


def _vectorized_regular_grid_pairs(
    rows: int,
    cols: int,
    *,
    connectivity: int,
    order: int,
) -> np.ndarray:
    if connectivity not in {4, 8}:
        raise ValueError("connectivity must be either 4 or 8.")
    if order <= 0:
        raise ValueError("order must be strictly positive.")
    grid = np.arange(rows * cols, dtype=np.int64).reshape(rows, cols)
    offsets: list[tuple[int, int]] = []
    for drow in range(-order, order + 1):
        for dcol in range(-order, order + 1):
            if drow == 0 and dcol == 0:
                continue
            if connectivity == 4 and abs(drow) + abs(dcol) > order:
                continue
            if connectivity == 8 and max(abs(drow), abs(dcol)) > order:
                continue
            if drow < 0 or (drow == 0 and dcol <= 0):
                continue
            offsets.append((drow, dcol))
    pair_blocks: list[np.ndarray] = []
    for drow, dcol in offsets:
        row_from = slice(0, rows - drow)
        row_to = slice(drow, rows)
        if dcol >= 0:
            col_from = slice(0, cols - dcol)
            col_to = slice(dcol, cols)
        else:
            col_from = slice(-dcol, cols)
            col_to = slice(0, cols + dcol)
        sources = grid[row_from, col_from].reshape(-1)
        targets = grid[row_to, col_to].reshape(-1)
        if sources.size == 0:
            continue
        pair_blocks.append(np.column_stack([sources, targets]))
    if not pair_blocks:
        return np.empty((0, 2), dtype=np.int64)
    pairs = np.concatenate(pair_blocks, axis=0)
    pairs = np.ascontiguousarray(pairs, dtype=np.int64)
    pairs.setflags(write=False)
    return pairs


@dataclass(frozen=True, slots=True)
class CrystalMapPhase:
    phase_id: int
    name: str
    symmetry: SymmetrySpec
    phase: Phase | None = None
    aliases: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        normalized_name = _canonical_phase_name(self.name)
        if self.phase_id < 0:
            raise ValueError("CrystalMapPhase.phase_id must be non-negative.")
        if self.symmetry.reference_frame is None:
            raise ValueError(
                "CrystalMapPhase.symmetry.reference_frame must be set for phase-resolved maps."
            )
        if self.phase is not None:
            if self.phase.name != normalized_name and normalized_name not in self.phase.aliases:
                raise ValueError(
                    "CrystalMapPhase.name must match phase.name or one of phase.aliases."
                )
            if self.phase.symmetry != self.symmetry:
                raise ValueError(
                    "CrystalMapPhase.phase.symmetry must match CrystalMapPhase.symmetry."
                )
            if self.phase.crystal_frame != self.symmetry.reference_frame:
                raise ValueError(
                    "CrystalMapPhase.phase.crystal_frame must match symmetry.reference_frame."
                )
        object.__setattr__(self, "name", normalized_name)
        object.__setattr__(
            self,
            "aliases",
            tuple(_canonical_phase_name(alias) for alias in self.aliases),
        )

    @property
    def crystal_frame(self) -> ReferenceFrame:
        if self.phase is not None:
            return self.phase.crystal_frame
        return cast(ReferenceFrame, self.symmetry.reference_frame)

    @property
    def point_group(self) -> str:
        return self.symmetry.point_group

    def matches(self, selector: int | str | Phase | CrystalMapPhase) -> bool:
        if isinstance(selector, CrystalMapPhase):
            return selector.phase_id == self.phase_id
        if isinstance(selector, Phase):
            if self.phase is not None:
                return self.phase == selector
            return selector.name == self.name or selector.name in self.aliases
        if isinstance(selector, int):
            return self.phase_id == selector
        normalized = selector.strip()
        return normalized == self.name or normalized in self.aliases


@dataclass(frozen=True, slots=True)
class CoordinateNeighborGraph:
    pairs: np.ndarray
    distances: np.ndarray
    connectivity: int
    order: int
    mode: str
    max_distance: float | None = None

    def __post_init__(self) -> None:
        pairs = _readonly_int_array(
            self.pairs,
            shape=(None, 2),
            name="CoordinateNeighborGraph.pairs",
        )
        distances = _readonly_float_array(
            self.distances,
            shape=(pairs.shape[0],),
            name="CoordinateNeighborGraph.distances",
        )
        if self.connectivity not in {4, 8}:
            raise ValueError("CoordinateNeighborGraph.connectivity must be either 4 or 8.")
        if self.order <= 0:
            raise ValueError("CoordinateNeighborGraph.order must be strictly positive.")
        if self.max_distance is not None and self.max_distance <= 0.0:
            raise ValueError("CoordinateNeighborGraph.max_distance must be positive when provided.")
        if self.mode not in {"regular_grid", "coordinate_radius"}:
            raise ValueError(
                "CoordinateNeighborGraph.mode must be 'regular_grid' or 'coordinate_radius'."
            )
        object.__setattr__(self, "pairs", pairs)
        object.__setattr__(self, "distances", distances)


@dataclass(frozen=True, slots=True)
class TextureReport:
    odf: ODF
    pole_figures: tuple[PoleFigure, ...] = ()
    inverse_pole_figures: tuple[InversePoleFigure, ...] = ()
    odf_figure: Any | None = None
    pole_figure_figures: tuple[Any, ...] = ()
    inverse_pole_figure_figures: tuple[Any, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "pole_figures", tuple(self.pole_figures))
        object.__setattr__(self, "inverse_pole_figures", tuple(self.inverse_pole_figures))
        object.__setattr__(self, "pole_figure_figures", tuple(self.pole_figure_figures))
        object.__setattr__(
            self,
            "inverse_pole_figure_figures",
            tuple(self.inverse_pole_figure_figures),
        )


@dataclass(frozen=True, slots=True)
class Grain:
    grain_id: int
    member_indices: np.ndarray
    mean_coordinate: np.ndarray
    reference_orientation_index: int
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        member_indices = np.asarray(self.member_indices, dtype=np.int64)
        if member_indices.ndim != 1 or member_indices.size == 0:
            raise ValueError("Grain.member_indices must be a non-empty 1D array.")
        if self.reference_orientation_index not in set(int(value) for value in member_indices):
            raise ValueError("Grain.reference_orientation_index must belong to member_indices.")
        mean_coordinate = np.asarray(self.mean_coordinate, dtype=np.float64)
        if mean_coordinate.ndim != 1:
            raise ValueError("Grain.mean_coordinate must be a 1D array.")
        member_indices = np.ascontiguousarray(member_indices)
        member_indices.setflags(write=False)
        mean_coordinate = np.ascontiguousarray(mean_coordinate)
        mean_coordinate.setflags(write=False)
        object.__setattr__(self, "member_indices", member_indices)
        object.__setattr__(self, "mean_coordinate", mean_coordinate)

    @property
    def size(self) -> int:
        return int(self.member_indices.size)


@dataclass(frozen=True, slots=True)
class GrainBoundarySegment:
    left_index: int
    right_index: int
    left_grain_id: int
    right_grain_id: int
    misorientation_deg: float
    length: float
    midpoint: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.left_index == self.right_index:
            raise ValueError("GrainBoundarySegment endpoints must be distinct.")
        if self.left_grain_id == self.right_grain_id:
            raise ValueError("GrainBoundarySegment must connect two distinct grains.")
        if self.misorientation_deg < 0.0:
            raise ValueError("GrainBoundarySegment.misorientation_deg must be non-negative.")
        if self.length < 0.0:
            raise ValueError("GrainBoundarySegment.length must be non-negative.")
        midpoint = np.asarray(self.midpoint, dtype=np.float64)
        if midpoint.ndim != 1:
            raise ValueError("GrainBoundarySegment.midpoint must be a 1D array.")
        midpoint = np.ascontiguousarray(midpoint)
        midpoint.setflags(write=False)
        object.__setattr__(self, "midpoint", midpoint)

    def classify(self, *, high_angle_threshold_deg: float = 15.0) -> str:
        if high_angle_threshold_deg < 0.0:
            raise ValueError("high_angle_threshold_deg must be non-negative.")
        return "high_angle" if self.misorientation_deg >= high_angle_threshold_deg else "low_angle"


@dataclass(frozen=True, slots=True)
class GrainBoundaryNetwork:
    segmentation: GrainSegmentation
    segments: tuple[GrainBoundarySegment, ...]
    min_misorientation_deg: float
    high_angle_threshold_deg: float = 15.0
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.min_misorientation_deg < 0.0:
            raise ValueError("GrainBoundaryNetwork.min_misorientation_deg must be non-negative.")
        if self.high_angle_threshold_deg < 0.0:
            raise ValueError("GrainBoundaryNetwork.high_angle_threshold_deg must be non-negative.")
        object.__setattr__(self, "segments", tuple(self.segments))

    @property
    def count(self) -> int:
        return len(self.segments)

    @property
    def mean_misorientation_deg(self) -> float:
        if not self.segments:
            return 0.0
        return float(np.mean([segment.misorientation_deg for segment in self.segments]))

    @property
    def total_length(self) -> float:
        return float(np.sum([segment.length for segment in self.segments]))

    @property
    def high_angle_count(self) -> int:
        return int(
            sum(
                segment.classify(high_angle_threshold_deg=self.high_angle_threshold_deg)
                == "high_angle"
                for segment in self.segments
            )
        )

    def grain_graph(self) -> GrainGraph:
        edge_groups: dict[tuple[int, int], list[tuple[int, GrainBoundarySegment]]] = {}
        for index, segment in enumerate(self.segments):
            left_grain_id, right_grain_id = sorted((segment.left_grain_id, segment.right_grain_id))
            edge_key = (left_grain_id, right_grain_id)
            edge_groups.setdefault(edge_key, []).append((index, segment))
        edges: list[GrainGraphEdge] = []
        for edge_key, members in sorted(edge_groups.items()):
            segment_indices = np.array(
                [member_index for member_index, _ in members], dtype=np.int64
            )
            segments = [segment for _, segment in members]
            total_length = float(np.sum([segment.length for segment in segments]))
            mean_misorientation = float(
                np.mean([segment.misorientation_deg for segment in segments])
            )
            high_angle_fraction = float(
                np.mean(
                    [
                        segment.classify(high_angle_threshold_deg=self.high_angle_threshold_deg)
                        == "high_angle"
                        for segment in segments
                    ]
                )
            )
            edges.append(
                GrainGraphEdge(
                    left_grain_id=edge_key[0],
                    right_grain_id=edge_key[1],
                    segment_indices=segment_indices,
                    total_length=total_length,
                    mean_misorientation_deg=mean_misorientation,
                    high_angle_fraction=high_angle_fraction,
                )
            )
        node_grain_ids = np.array(
            [grain.grain_id for grain in self.segmentation.grains], dtype=np.int64
        )
        return GrainGraph(
            segmentation=self.segmentation,
            edges=tuple(edges),
            node_grain_ids=node_grain_ids,
            high_angle_threshold_deg=self.high_angle_threshold_deg,
            provenance=self.provenance,
        )


@dataclass(frozen=True, slots=True)
class GrainGraphEdge:
    left_grain_id: int
    right_grain_id: int
    segment_indices: np.ndarray
    total_length: float
    mean_misorientation_deg: float
    high_angle_fraction: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.left_grain_id == self.right_grain_id:
            raise ValueError("GrainGraphEdge must connect distinct grains.")
        object.__setattr__(
            self,
            "segment_indices",
            np.ascontiguousarray(np.asarray(self.segment_indices, dtype=np.int64)),
        )
        if self.segment_indices.ndim != 1 or self.segment_indices.size == 0:
            raise ValueError("GrainGraphEdge.segment_indices must be a non-empty 1D array.")
        if self.total_length < 0.0:
            raise ValueError("GrainGraphEdge.total_length must be non-negative.")
        if self.mean_misorientation_deg < 0.0:
            raise ValueError("GrainGraphEdge.mean_misorientation_deg must be non-negative.")
        if not 0.0 <= self.high_angle_fraction <= 1.0:
            raise ValueError("GrainGraphEdge.high_angle_fraction must lie in [0, 1].")

    @property
    def grain_pair(self) -> tuple[int, int]:
        return (self.left_grain_id, self.right_grain_id)


@dataclass(frozen=True, slots=True)
class GrainGraph:
    segmentation: GrainSegmentation
    edges: tuple[GrainGraphEdge, ...]
    node_grain_ids: np.ndarray
    high_angle_threshold_deg: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "edges", tuple(self.edges))
        node_ids = np.ascontiguousarray(np.asarray(self.node_grain_ids, dtype=np.int64))
        if node_ids.ndim != 1:
            raise ValueError("GrainGraph.node_grain_ids must be a 1D array.")
        object.__setattr__(self, "node_grain_ids", node_ids)
        if self.high_angle_threshold_deg < 0.0:
            raise ValueError("GrainGraph.high_angle_threshold_deg must be non-negative.")

    @property
    def edge_count(self) -> int:
        return len(self.edges)

    @property
    def adjacency_matrix(self) -> np.ndarray:
        node_index = {int(grain_id): idx for idx, grain_id in enumerate(self.node_grain_ids)}
        matrix = np.zeros((len(self.node_grain_ids), len(self.node_grain_ids)), dtype=np.int64)
        for edge in self.edges:
            left = node_index[int(edge.left_grain_id)]
            right = node_index[int(edge.right_grain_id)]
            matrix[left, right] = 1
            matrix[right, left] = 1
        matrix = np.ascontiguousarray(matrix)
        matrix.setflags(write=False)
        return matrix

    def neighbors(self, grain_id: int) -> np.ndarray:
        neighbors = [
            edge.right_grain_id if edge.left_grain_id == grain_id else edge.left_grain_id
            for edge in self.edges
            if grain_id in edge.grain_pair
        ]
        array = np.asarray(sorted(neighbors), dtype=np.int64)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        return array


@dataclass(frozen=True, slots=True)
class GrainSegmentation:
    crystal_map: CrystalMap
    labels: np.ndarray
    grains: tuple[Grain, ...]
    max_misorientation_deg: float
    connectivity: int
    symmetry_aware: bool
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        labels = np.asarray(self.labels, dtype=np.int64)
        if labels.shape != (len(self.crystal_map.orientations),):
            raise ValueError(
                "GrainSegmentation.labels must have one entry per CrystalMap orientation."
            )
        if self.max_misorientation_deg < 0.0:
            raise ValueError("GrainSegmentation.max_misorientation_deg must be non-negative.")
        if self.connectivity not in {4, 8}:
            raise ValueError("GrainSegmentation.connectivity must be either 4 or 8.")
        labels = np.ascontiguousarray(labels)
        labels.setflags(write=False)
        object.__setattr__(self, "labels", labels)
        object.__setattr__(self, "grains", tuple(self.grains))

    @property
    def label_grid(self) -> np.ndarray:
        rows, cols = self.crystal_map._require_regular_2d_grid()
        labels = np.ascontiguousarray(self.labels.reshape((rows, cols)))
        labels.setflags(write=False)
        return labels

    def reference_orientation(self, grain: Grain) -> Orientation:
        return self.crystal_map.orientations[grain.reference_orientation_index]

    def grod_map_deg(self) -> np.ndarray:
        rows, cols = self.crystal_map._require_regular_2d_grid()
        deviations = np.zeros(len(self.crystal_map.orientations), dtype=np.float64)
        for grain in self.grains:
            reference = self.reference_orientation(grain)
            for index in grain.member_indices:
                deviations[int(index)] = float(
                    np.rad2deg(
                        self.crystal_map.orientations[int(index)].distance_to(
                            reference,
                            symmetry_aware=self.symmetry_aware,
                        )
                    )
                )
        deviations = np.ascontiguousarray(deviations.reshape((rows, cols)))
        deviations.setflags(write=False)
        return deviations

    def boundary_network(
        self,
        *,
        min_misorientation_deg: float = 0.0,
        high_angle_threshold_deg: float = 15.0,
    ) -> GrainBoundaryNetwork:
        if min_misorientation_deg < 0.0:
            raise ValueError("min_misorientation_deg must be non-negative.")
        segments: list[GrainBoundarySegment] = []
        neighbor_pairs = self.crystal_map.neighbor_pairs(connectivity=self.connectivity)
        neighbor_pairs = neighbor_pairs[self.crystal_map._same_phase_pair_mask(neighbor_pairs)]
        for left_index, right_index in neighbor_pairs:
            left_label = int(self.labels[int(left_index)])
            right_label = int(self.labels[int(right_index)])
            if left_label == right_label:
                continue
            misorientation_deg = float(
                np.rad2deg(
                    self.crystal_map.orientations[int(left_index)].distance_to(
                        self.crystal_map.orientations[int(right_index)],
                        symmetry_aware=self.symmetry_aware,
                    )
                )
            )
            if misorientation_deg < min_misorientation_deg:
                continue
            segments.append(
                GrainBoundarySegment(
                    left_index=int(left_index),
                    right_index=int(right_index),
                    left_grain_id=left_label,
                    right_grain_id=right_label,
                    misorientation_deg=misorientation_deg,
                    length=float(
                        np.linalg.norm(
                            self.crystal_map.coordinates[int(left_index)]
                            - self.crystal_map.coordinates[int(right_index)]
                        )
                    ),
                    midpoint=np.mean(
                        self.crystal_map.coordinates[[int(left_index), int(right_index)]],
                        axis=0,
                    ),
                    provenance=self.provenance,
                )
            )
        return GrainBoundaryNetwork(
            segmentation=self,
            segments=tuple(segments),
            min_misorientation_deg=min_misorientation_deg,
            high_angle_threshold_deg=high_angle_threshold_deg,
            provenance=self.provenance,
        )

    def grain_graph(
        self,
        *,
        min_misorientation_deg: float = 0.0,
        high_angle_threshold_deg: float = 15.0,
    ) -> GrainGraph:
        return self.boundary_network(
            min_misorientation_deg=min_misorientation_deg,
            high_angle_threshold_deg=high_angle_threshold_deg,
        ).grain_graph()

    def majority_smoothed(
        self,
        *,
        iterations: int = 1,
        min_neighbor_votes: int = 3,
    ) -> GrainSegmentation:
        if iterations <= 0:
            raise ValueError("iterations must be strictly positive.")
        if min_neighbor_votes <= 0:
            raise ValueError("min_neighbor_votes must be strictly positive.")
        self.crystal_map._require_regular_2d_grid()
        labels = np.array(self.labels, copy=True)
        neighbor_pairs = self.crystal_map.neighbor_pairs(connectivity=self.connectivity)
        neighbor_pairs = neighbor_pairs[self.crystal_map._same_phase_pair_mask(neighbor_pairs)]
        adjacency: dict[int, list[int]] = {index: [] for index in range(len(labels))}
        for left_index, right_index in neighbor_pairs:
            adjacency[int(left_index)].append(int(right_index))
            adjacency[int(right_index)].append(int(left_index))
        for _ in range(iterations):
            updated = labels.copy()
            for index in range(len(labels)):
                neighbor_labels = [int(labels[neighbor]) for neighbor in adjacency[index]]
                if not neighbor_labels:
                    continue
                unique, counts = np.unique(
                    np.asarray(neighbor_labels, dtype=np.int64), return_counts=True
                )
                best_position = int(np.argmax(counts))
                best_label = int(unique[best_position])
                best_count = int(counts[best_position])
                if best_count >= min_neighbor_votes and best_label != int(labels[index]):
                    updated[index] = best_label
            labels = updated
        smoothed = self.crystal_map._segmentation_from_labels(
            labels,
            max_misorientation_deg=self.max_misorientation_deg,
            symmetry_aware=self.symmetry_aware,
            connectivity=self.connectivity,
        )
        return smoothed

    def merge_small_grains(
        self,
        *,
        min_size: int,
        until_stable: bool = True,
        max_iterations: int | None = None,
    ) -> GrainSegmentation:
        if min_size <= 0:
            raise ValueError("min_size must be strictly positive.")
        if min_size <= 1:
            return self
        labels = np.array(self.labels, copy=True)
        iterations = 0
        while True:
            if max_iterations is not None and iterations >= max_iterations:
                return self.crystal_map._segmentation_from_labels(
                    labels,
                    max_misorientation_deg=self.max_misorientation_deg,
                    symmetry_aware=self.symmetry_aware,
                    connectivity=self.connectivity,
                )
            current = self.crystal_map._segmentation_from_labels(
                labels,
                max_misorientation_deg=self.max_misorientation_deg,
                symmetry_aware=self.symmetry_aware,
                connectivity=self.connectivity,
            )
            small_grains = sorted(
                (grain for grain in current.grains if grain.size < min_size),
                key=lambda grain: (grain.size, grain.grain_id),
            )
            if not small_grains:
                return current
            grain = small_grains[0]
            adjacency: dict[int, tuple[int, float]] = {}
            neighbor_pairs = current.crystal_map.neighbor_pairs(connectivity=current.connectivity)
            neighbor_pairs = neighbor_pairs[
                current.crystal_map._same_phase_pair_mask(neighbor_pairs)
            ]
            for left_index, right_index in neighbor_pairs:
                left_label = int(current.labels[int(left_index)])
                right_label = int(current.labels[int(right_index)])
                if left_label == right_label:
                    continue
                if left_label == grain.grain_id:
                    target_label = right_label
                    source_index = int(left_index)
                    target_index = int(right_index)
                elif right_label == grain.grain_id:
                    target_label = left_label
                    source_index = int(right_index)
                    target_index = int(left_index)
                else:
                    continue
                misorientation_deg = float(
                    np.rad2deg(
                        current.crystal_map.orientations[source_index].distance_to(
                            current.crystal_map.orientations[target_index],
                            symmetry_aware=current.symmetry_aware,
                        )
                    )
                )
                count, total = adjacency.get(target_label, (0, 0.0))
                adjacency[target_label] = (count + 1, total + misorientation_deg)
            if not adjacency:
                return current
            target_label = min(
                adjacency,
                key=lambda label: (
                    -adjacency[label][0],
                    adjacency[label][1] / adjacency[label][0],
                    label,
                ),
            )
            labels[current.labels == grain.grain_id] = target_label
            iterations += 1
            if not until_stable:
                return self.crystal_map._segmentation_from_labels(
                    labels,
                    max_misorientation_deg=self.max_misorientation_deg,
                    symmetry_aware=self.symmetry_aware,
                    connectivity=self.connectivity,
                )


@dataclass(frozen=True, slots=True)
class CrystalMap:
    coordinates: np.ndarray
    orientations: OrientationSet
    map_frame: ReferenceFrame
    phase_entries: tuple[CrystalMapPhase, ...] = ()
    phase_ids: np.ndarray | None = None
    grid_shape: tuple[int, ...] | None = None
    step_sizes: tuple[float, ...] | None = None
    acquisition_geometry: AcquisitionGeometry | None = None
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] not in (2, 3):
            raise ValueError("CrystalMap.coordinates must have shape (n, 2) or (n, 3).")
        if coordinates.shape[0] != len(self.orientations):
            raise ValueError("CrystalMap coordinates and orientations must have matching lengths.")
        if self.map_frame.domain not in {FrameDomain.MAP, FrameDomain.SPECIMEN}:
            raise ValueError("CrystalMap.map_frame must belong to the map or specimen domain.")
        phase_entries = tuple(self.phase_entries)
        if phase_entries:
            phase_ids_seen: set[int] = set()
            phase_names_seen: set[str] = set()
            for entry in phase_entries:
                if entry.phase_id in phase_ids_seen:
                    raise ValueError("CrystalMap.phase_entries phase_id values must be unique.")
                if entry.name in phase_names_seen:
                    raise ValueError("CrystalMap.phase_entries names must be unique.")
                if entry.crystal_frame != self.orientations.crystal_frame:
                    raise ValueError(
                        "CrystalMap.phase_entries must share OrientationSet.crystal_frame."
                    )
                phase_ids_seen.add(entry.phase_id)
                phase_names_seen.add(entry.name)
            if len(phase_entries) > 1 and (
                self.orientations.phase is not None or self.orientations.symmetry is not None
            ):
                raise ValueError(
                    "Multiphase CrystalMap instances require OrientationSet.phase and "
                    "OrientationSet.symmetry to be None so phase semantics remain attached "
                    "explicitly through CrystalMap.phase_entries."
                )
        phase_ids = None
        if self.phase_ids is not None:
            phase_ids = np.asarray(self.phase_ids, dtype=np.int64)
            if phase_ids.shape != (len(self.orientations),):
                raise ValueError("CrystalMap.phase_ids must have one entry per orientation.")
            if np.any(phase_ids < 0):
                raise ValueError("CrystalMap.phase_ids must be non-negative.")
            if not phase_entries:
                raise ValueError(
                    "CrystalMap.phase_ids requires CrystalMap.phase_entries to be provided."
                )
            available_ids = {entry.phase_id for entry in phase_entries}
            if any(int(value) not in available_ids for value in np.unique(phase_ids)):
                raise ValueError("CrystalMap.phase_ids must refer only to declared phase_ids.")
            phase_ids = np.ascontiguousarray(phase_ids)
            phase_ids.setflags(write=False)
        elif phase_entries:
            if len(phase_entries) == 1:
                phase_ids = np.zeros(len(self.orientations), dtype=np.int64)
                phase_ids.setflags(write=False)
            else:
                raise ValueError(
                    "Multiphase CrystalMap instances require CrystalMap.phase_ids to be provided."
                )
        if (
            self.orientations.phase is not None
            and phase_entries
            and phase_entries[0].phase is not None
            and self.orientations.phase != phase_entries[0].phase
        ):
            raise ValueError(
                "Single-phase CrystalMap.phase_entries[0].phase must match OrientationSet.phase."
            )
        if self.acquisition_geometry is not None:
            if self.acquisition_geometry.specimen_frame != self.orientations.specimen_frame:
                raise ValueError(
                    "CrystalMap.acquisition_geometry.specimen_frame must match "
                    "CrystalMap.orientations.specimen_frame."
                )
            if self.acquisition_geometry.map_frame is not None:
                if self.map_frame != self.acquisition_geometry.map_frame:
                    raise ValueError(
                        "CrystalMap.map_frame must match "
                        "CrystalMap.acquisition_geometry.map_frame when provided."
                    )
            elif self.map_frame != self.orientations.specimen_frame:
                raise ValueError(
                    "CrystalMap.map_frame must equal the specimen frame when no acquisition "
                    "map_frame is provided."
                )
            if (
                self.calibration_record is not None
                and self.acquisition_geometry.calibration_record is not None
                and self.calibration_record != self.acquisition_geometry.calibration_record
            ):
                raise ValueError(
                    "CrystalMap.calibration_record must match the acquisition geometry "
                    "calibration record when both are provided."
                )
            if (
                self.measurement_quality is not None
                and self.acquisition_geometry.measurement_quality is not None
                and self.measurement_quality != self.acquisition_geometry.measurement_quality
            ):
                raise ValueError(
                    "CrystalMap.measurement_quality must match the acquisition geometry "
                    "measurement quality when both are provided."
                )
        coordinates = np.ascontiguousarray(coordinates)
        coordinates.setflags(write=False)
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "phase_entries", phase_entries)
        object.__setattr__(self, "phase_ids", phase_ids)
        coordinate_dims = int(coordinates.shape[1])
        if self.grid_shape is not None:
            if len(self.grid_shape) != coordinate_dims:
                raise ValueError("CrystalMap.grid_shape must match the coordinate dimensionality.")
            if any(size <= 0 for size in self.grid_shape):
                raise ValueError("CrystalMap.grid_shape entries must be strictly positive.")
        if self.step_sizes is not None:
            if len(self.step_sizes) != coordinate_dims:
                raise ValueError("CrystalMap.step_sizes must match the coordinate dimensionality.")
            step_sizes = tuple(float(step) for step in self.step_sizes)
            if any(step <= 0.0 for step in step_sizes):
                raise ValueError("CrystalMap.step_sizes entries must be strictly positive.")
            object.__setattr__(self, "step_sizes", step_sizes)

    @property
    def is_multiphase(self) -> bool:
        return len(self.phase_entries) > 1

    @property
    def has_phase_assignments(self) -> bool:
        return bool(self.phase_entries) or self.orientations.phase is not None

    @property
    def phase_id_array(self) -> np.ndarray | None:
        if self.phase_ids is not None:
            return self.phase_ids
        if self.orientations.phase is None and self.orientations.symmetry is None:
            return None
        values = np.zeros(len(self.orientations), dtype=np.int64)
        values.setflags(write=False)
        return values

    @property
    def resolved_phase_entries(self) -> tuple[CrystalMapPhase, ...]:
        if self.phase_entries:
            return self.phase_entries
        if self.orientations.phase is None and self.orientations.symmetry is None:
            return ()
        resolved_symmetry = (
            self.orientations.phase.symmetry
            if self.orientations.phase is not None
            else self.orientations.symmetry
        )
        if resolved_symmetry is None:
            raise ValueError(
                "CrystalMap requires phase-resolved symmetry when phase entries are synthesized."
            )
        return (
            CrystalMapPhase(
                phase_id=0,
                name=(
                    self.orientations.phase.name
                    if self.orientations.phase is not None
                    else "unresolved_phase"
                ),
                symmetry=resolved_symmetry,
                phase=self.orientations.phase,
                provenance=self.provenance,
            ),
        )

    @property
    def primary_phase(self) -> Phase | None:
        if self.orientations.phase is not None and not self.phase_entries:
            return self.orientations.phase
        if len(self.resolved_phase_entries) == 1:
            return self.resolved_phase_entries[0].phase
        return None

    def phase_summary(self) -> dict[str, int]:
        phase_entries = self.resolved_phase_entries
        phase_ids = self.phase_id_array
        if not phase_entries or phase_ids is None:
            return {}
        return {
            entry.name: int(np.count_nonzero(phase_ids == entry.phase_id))
            for entry in phase_entries
        }

    def summary(self) -> dict[str, Any]:
        return {
            "point_count": len(self.orientations),
            "coordinate_dimensions": int(self.coordinates.shape[1]),
            "grid_shape": (
                None
                if self.grid_shape is None
                else tuple(int(value) for value in self.grid_shape)
            ),
            "step_sizes": (
                None
                if self.step_sizes is None
                else tuple(float(value) for value in self.step_sizes)
            ),
            "is_multiphase": self.is_multiphase,
            "phases": self.phase_summary(),
            "map_frame": self.map_frame.name,
            "specimen_frame": self.orientations.specimen_frame.name,
        }

    def validate(self) -> tuple[str, ...]:
        notes: list[str] = []
        if self.is_multiphase:
            unresolved = [
                entry.name for entry in self.resolved_phase_entries if entry.phase is None
            ]
            if unresolved:
                notes.append(
                    "Full Phase objects are not attached for phases: " + ", ".join(unresolved)
                )
        if self.grid_shape is None:
            notes.append("Map is operating in graph mode rather than regular-grid mode.")
        return tuple(notes)

    def _phase_entry_by_id(self) -> dict[int, CrystalMapPhase]:
        return {entry.phase_id: entry for entry in self.resolved_phase_entries}

    def _resolve_phase_entry(
        self, selector: int | str | Phase | CrystalMapPhase
    ) -> CrystalMapPhase:
        for entry in self.resolved_phase_entries:
            if entry.matches(selector):
                return entry
        raise ValueError(f"Unknown phase selector: {selector!r}.")

    def phase_mask(self, selector: int | str | Phase | CrystalMapPhase) -> np.ndarray:
        phase_ids = self.phase_id_array
        if phase_ids is None:
            raise ValueError("CrystalMap does not carry explicit phase assignments.")
        entry = self._resolve_phase_entry(selector)
        mask = np.ascontiguousarray(phase_ids == entry.phase_id)
        mask.setflags(write=False)
        return mask

    def phase_entry_for_index(self, index: int) -> CrystalMapPhase | None:
        phase_ids = self.phase_id_array
        if phase_ids is None:
            return None
        return self._phase_entry_by_id()[int(phase_ids[int(index)])]

    def select_phase(self, selector: int | str | Phase | CrystalMapPhase) -> CrystalMap:
        mask = np.asarray(self.phase_mask(selector), dtype=bool)
        if not np.any(mask):
            raise ValueError("Selected phase has no points in this CrystalMap.")
        entry = self._resolve_phase_entry(selector)
        provenance = self.provenance
        phase = entry.phase
        if phase is not None:
            orientations = OrientationSet.from_quaternions(
                self.orientations.quaternions[mask],
                crystal_frame=self.orientations.crystal_frame,
                specimen_frame=self.orientations.specimen_frame,
                phase=phase,
                provenance=provenance,
            )
        else:
            orientations = OrientationSet.from_quaternions(
                self.orientations.quaternions[mask],
                crystal_frame=self.orientations.crystal_frame,
                specimen_frame=self.orientations.specimen_frame,
                symmetry=entry.symmetry,
                provenance=provenance,
            )
        full_selection = bool(np.all(mask))
        return CrystalMap(
            coordinates=self.coordinates[mask],
            orientations=orientations,
            map_frame=self.map_frame,
            grid_shape=self.grid_shape if full_selection else None,
            step_sizes=self.step_sizes,
            acquisition_geometry=self.acquisition_geometry,
            calibration_record=self.calibration_record,
            measurement_quality=self.measurement_quality,
            provenance=provenance,
        )

    def _same_phase_pair_mask(self, pairs: np.ndarray) -> np.ndarray:
        phase_ids = self.phase_id_array
        if phase_ids is None:
            return np.ones(pairs.shape[0], dtype=bool)
        return np.asarray(phase_ids[pairs[:, 0]] == phase_ids[pairs[:, 1]], dtype=bool)

    def _pair_misorientation_rad(
        self,
        pairs: np.ndarray,
        *,
        symmetry_aware: bool,
    ) -> np.ndarray:
        if pairs.size == 0:
            return np.empty(0, dtype=np.float64)
        matrices = self.orientations.as_matrices()
        relative = _relative_rotation_matrices(matrices[pairs[:, 0]], matrices[pairs[:, 1]])
        if not symmetry_aware:
            return _rotation_angles_from_matrices(relative)
        phase_ids = self.phase_id_array
        if phase_ids is None:
            return _disorientation_angles_from_relative_matrices(
                relative,
                left_symmetry=self.orientations.symmetry,
                right_symmetry=self.orientations.symmetry,
            )
        angles = np.full(pairs.shape[0], np.nan, dtype=np.float64)
        phase_lookup = self._phase_entry_by_id()
        left_ids = phase_ids[pairs[:, 0]]
        right_ids = phase_ids[pairs[:, 1]]
        same_phase = left_ids == right_ids
        if not np.any(same_phase):
            return angles
        for phase_id in np.unique(left_ids[same_phase]):
            mask = same_phase & (left_ids == phase_id)
            entry = phase_lookup[int(phase_id)]
            angles[mask] = _disorientation_angles_from_relative_matrices(
                relative[mask],
                left_symmetry=entry.symmetry,
                right_symmetry=entry.symmetry,
            )
        return angles

    def _require_regular_2d_grid(self) -> tuple[int, int]:
        if self.grid_shape is None or len(self.grid_shape) != 2:
            raise ValueError("CrystalMap regular-grid workflows require a 2D grid_shape.")
        rows, cols = self.grid_shape
        if rows * cols != len(self.orientations):
            raise ValueError(
                "CrystalMap.grid_shape must match the number of orientations for regular-grid "
                "workflows."
            )
        return int(rows), int(cols)

    def to_experiment_manifest(
        self,
        *,
        source_system: str = "pytex",
        referenced_files: tuple[str, ...] = (),
        metadata: dict[str, str] | None = None,
    ) -> ExperimentManifest:
        from pytex.adapters import ExperimentManifest

        acquisition_geometry = self.acquisition_geometry
        if acquisition_geometry is None:
            acquisition_geometry = AcquisitionGeometry(
                specimen_frame=self.orientations.specimen_frame,
                modality="ebsd",
                map_frame=self.map_frame if self.map_frame.domain is FrameDomain.MAP else None,
                calibration_record=self.calibration_record,
                measurement_quality=self.measurement_quality,
                provenance=self.provenance,
            )
        merged_metadata: dict[str, str] = {}
        if self.grid_shape is not None:
            merged_metadata["grid_shape"] = "x".join(str(value) for value in self.grid_shape)
        if self.step_sizes is not None:
            merged_metadata["step_sizes"] = ",".join(f"{value:g}" for value in self.step_sizes)
        if self.is_multiphase:
            merged_metadata["phase_names"] = ",".join(
                entry.name for entry in self.resolved_phase_entries
            )
        if metadata is not None:
            merged_metadata.update(metadata)
        return ExperimentManifest.from_acquisition_geometry(
            acquisition_geometry,
            source_system=source_system,
            phase=self.primary_phase,
            phases=tuple(
                entry.phase for entry in self.resolved_phase_entries if entry.phase is not None
            ),
            referenced_files=referenced_files,
            metadata=merged_metadata,
        )

    def _phase_resolved_view(
        self,
        *,
        phase: int | str | Phase | CrystalMapPhase | None,
        operation: str,
    ) -> CrystalMap:
        if phase is None:
            if self.is_multiphase:
                raise ValueError(
                    f"CrystalMap.{operation}() requires a phase selector for multiphase maps."
                )
            return self
        return self.select_phase(phase)

    def to_odf(
        self,
        *,
        phase: int | str | Phase | CrystalMapPhase | None = None,
        weights: ArrayLike | None = None,
        kernel: KernelSpec | None = None,
        specimen_symmetry: SymmetrySpec | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> ODF:
        from pytex.texture import ODF

        phase_view = self._phase_resolved_view(phase=phase, operation="to_odf")
        return ODF.from_orientations(
            phase_view.orientations,
            weights=weights,
            kernel=kernel,
            specimen_symmetry=specimen_symmetry,
            provenance=phase_view.provenance if provenance is None else provenance,
        )

    def pole_figure(
        self,
        pole: CrystalPlane | ArrayLike,
        *,
        phase: int | str | Phase | CrystalMapPhase | None = None,
        weights: ArrayLike | None = None,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
        sample_symmetry: SymmetrySpec | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigure:
        from pytex.texture import PoleFigure

        phase_view = self._phase_resolved_view(phase=phase, operation="pole_figure")
        return PoleFigure.from_orientations(
            phase_view.orientations,
            _coerce_pole(pole, phase=phase_view.orientations.phase),
            weights=weights,
            include_symmetry_family=include_symmetry_family,
            antipodal=antipodal,
            sample_symmetry=sample_symmetry,
            provenance=phase_view.provenance if provenance is None else provenance,
        )

    def inverse_pole_figure(
        self,
        sample_direction: str | ArrayLike = "z",
        *,
        phase: int | str | Phase | CrystalMapPhase | None = None,
        weights: ArrayLike | None = None,
        reduce_by_symmetry: bool = True,
        antipodal: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> InversePoleFigure:
        from pytex.texture import InversePoleFigure

        phase_view = self._phase_resolved_view(phase=phase, operation="inverse_pole_figure")
        return InversePoleFigure.from_orientations(
            phase_view.orientations,
            _specimen_direction_vector(sample_direction, phase_view.orientations.specimen_frame),
            weights=weights,
            reduce_by_symmetry=reduce_by_symmetry,
            antipodal=antipodal,
            provenance=phase_view.provenance if provenance is None else provenance,
        )

    def texture_report(
        self,
        *,
        poles: CrystalPlane
        | ArrayLike
        | tuple[CrystalPlane | ArrayLike, ...]
        | list[CrystalPlane | ArrayLike] = (),
        sample_directions: str
        | ArrayLike
        | tuple[str | ArrayLike, ...]
        | list[str | ArrayLike] = ("x", "y", "z"),
        phase: int | str | Phase | CrystalMapPhase | None = None,
        weights: ArrayLike | None = None,
        kernel: KernelSpec | None = None,
        specimen_symmetry: SymmetrySpec | None = None,
        include_symmetry_family: bool = True,
        reduce_by_symmetry: bool = True,
        antipodal: bool = True,
        plot: bool = False,
        provenance: ProvenanceRecord | None = None,
        pole_figure_plot_kwargs: dict[str, Any] | None = None,
        inverse_pole_figure_plot_kwargs: dict[str, Any] | None = None,
        odf_plot_kwargs: dict[str, Any] | None = None,
    ) -> TextureReport:
        phase_view = self._phase_resolved_view(phase=phase, operation="texture_report")
        report_provenance = phase_view.provenance if provenance is None else provenance
        odf = phase_view.to_odf(
            weights=weights,
            kernel=kernel,
            specimen_symmetry=specimen_symmetry,
            provenance=report_provenance,
        )
        if isinstance(poles, (list, tuple)) and len(poles) == 0:
            pole_sequence: tuple[CrystalPlane, ...] = ()
        else:
            pole_sequence = _coerce_pole_sequence(poles, phase=phase_view.orientations.phase)
        direction_sequence = _coerce_sample_direction_sequence(
            sample_directions,
            phase_view.orientations.specimen_frame,
        )
        pole_figures = tuple(
            phase_view.pole_figure(
                pole,
                weights=weights,
                include_symmetry_family=include_symmetry_family,
                antipodal=antipodal,
                sample_symmetry=specimen_symmetry,
                provenance=report_provenance,
            )
            for pole in pole_sequence
        )
        inverse_pole_figures = tuple(
            phase_view.inverse_pole_figure(
                sample_direction=direction,
                weights=weights,
                reduce_by_symmetry=reduce_by_symmetry,
                antipodal=antipodal,
                provenance=report_provenance,
            )
            for direction in direction_sequence
        )
        odf_figure: Any | None = None
        pole_figure_figures: tuple[Any, ...] = ()
        inverse_pole_figure_figures: tuple[Any, ...] = ()
        if plot:
            from pytex.plotting import plot_inverse_pole_figure, plot_odf, plot_pole_figure

            odf_figure = plot_odf(odf, **(odf_plot_kwargs or {}))
            pole_figure_figures = tuple(
                plot_pole_figure(pole_figure, **(pole_figure_plot_kwargs or {}))
                for pole_figure in pole_figures
            )
            inverse_pole_figure_figures = tuple(
                plot_inverse_pole_figure(
                    inverse_pole_figure,
                    **(inverse_pole_figure_plot_kwargs or {}),
                )
                for inverse_pole_figure in inverse_pole_figures
            )
        return TextureReport(
            odf=odf,
            pole_figures=pole_figures,
            inverse_pole_figures=inverse_pole_figures,
            odf_figure=odf_figure,
            pole_figure_figures=pole_figure_figures,
            inverse_pole_figure_figures=inverse_pole_figure_figures,
            provenance=report_provenance,
        )

    def neighbor_graph(
        self,
        *,
        connectivity: int = 4,
        order: int = 1,
        max_distance: float | None = None,
    ) -> CoordinateNeighborGraph:
        if connectivity not in {4, 8}:
            raise ValueError("connectivity must be either 4 or 8.")
        if order <= 0:
            raise ValueError("order must be strictly positive.")
        if self.grid_shape is not None and len(self.grid_shape) == 2:
            rows, cols = self._require_regular_2d_grid()
            if max_distance is None:
                pairs = _vectorized_regular_grid_pairs(
                    rows,
                    cols,
                    connectivity=connectivity,
                    order=order,
                )
                distances = np.linalg.norm(
                    self.coordinates[pairs[:, 0]] - self.coordinates[pairs[:, 1]],
                    axis=1,
                )
                return CoordinateNeighborGraph(
                    pairs=pairs,
                    distances=distances,
                    connectivity=connectivity,
                    order=order,
                    mode="regular_grid",
                )
        radius = max_distance
        if radius is None:
            base_spacing = _inferred_base_spacing(self.coordinates, self.step_sizes)
            radius_scale = float(order) * (np.sqrt(2.0) if connectivity == 8 else 1.0)
            radius = base_spacing * radius_scale + 1e-9
        distances = _pairwise_distances(self.coordinates)
        upper_mask = np.triu(np.ones_like(distances, dtype=bool), k=1)
        pair_mask = upper_mask & (distances <= radius)
        indices = np.column_stack(np.nonzero(pair_mask)).astype(np.int64)
        pair_distances = np.asarray(distances[pair_mask], dtype=np.float64)
        return CoordinateNeighborGraph(
            pairs=indices,
            distances=pair_distances,
            connectivity=connectivity,
            order=order,
            mode="coordinate_radius",
            max_distance=radius,
        )

    def neighbor_pairs(self, *, connectivity: int = 4) -> np.ndarray:
        return self.neighbor_graph(connectivity=connectivity, order=1).pairs

    def kam_neighbor_pairs(self, *, order: int = 1) -> np.ndarray:
        return self.neighbor_graph(connectivity=4, order=order).pairs

    def kernel_average_misorientation_deg(
        self,
        *,
        symmetry_aware: bool = True,
        connectivity: int = 4,
        order: int = 1,
        threshold_deg: float | None = None,
        statistic: str = "mean",
        segmentation: GrainSegmentation | None = None,
    ) -> np.ndarray:
        if connectivity not in {4, 8}:
            raise ValueError("connectivity must be either 4 or 8.")
        if threshold_deg is not None and threshold_deg < 0.0:
            raise ValueError("threshold_deg must be non-negative when provided.")
        if statistic not in {"mean", "max"}:
            raise ValueError("statistic must be either 'mean' or 'max'.")
        if segmentation is not None and segmentation.crystal_map is not self:
            raise ValueError("segmentation.crystal_map must be this CrystalMap instance.")
        graph = self.neighbor_graph(connectivity=connectivity, order=order)
        neighbor_pairs = graph.pairs
        valid_mask = self._same_phase_pair_mask(neighbor_pairs)
        if segmentation is not None:
            valid_mask &= (
                segmentation.labels[neighbor_pairs[:, 0]]
                == segmentation.labels[neighbor_pairs[:, 1]]
            )
        filtered_pairs = neighbor_pairs[valid_mask]
        sums = np.zeros(len(self.orientations), dtype=np.float64)
        counts = np.zeros(len(self.orientations), dtype=np.int64)
        maxima = np.zeros(len(self.orientations), dtype=np.float64)
        if filtered_pairs.size:
            angle_deg = np.rad2deg(
                self._pair_misorientation_rad(filtered_pairs, symmetry_aware=symmetry_aware)
            )
            finite_mask = np.isfinite(angle_deg)
            if threshold_deg is not None:
                finite_mask &= angle_deg <= threshold_deg
            filtered_pairs = filtered_pairs[finite_mask]
            angle_deg = angle_deg[finite_mask]
            if filtered_pairs.size:
                pair_indices = np.concatenate([filtered_pairs[:, 0], filtered_pairs[:, 1]])
                pair_values = np.concatenate([angle_deg, angle_deg])
                sums = np.asarray(
                    np.bincount(
                        pair_indices,
                        weights=pair_values,
                        minlength=len(self.orientations),
                    ),
                    dtype=np.float64,
                )
                counts = np.asarray(
                    np.bincount(pair_indices, minlength=len(self.orientations)),
                    dtype=np.int64,
                )
                np.maximum.at(maxima, filtered_pairs[:, 0], angle_deg)
                np.maximum.at(maxima, filtered_pairs[:, 1], angle_deg)
        if statistic == "mean":
            with np.errstate(divide="ignore", invalid="ignore"):
                values = np.divide(
                    sums,
                    counts,
                    out=np.zeros_like(sums),
                    where=counts > 0,
                )
        else:
            values = maxima
        if self.grid_shape is not None and len(self.grid_shape) == 2:
            rows, cols = self._require_regular_2d_grid()
            values = values.reshape((rows, cols))
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def _representative_orientation_index(
        self,
        member_indices: np.ndarray,
        *,
        symmetry_aware: bool,
    ) -> int:
        if member_indices.size == 1:
            return int(member_indices[0])
        best_index = int(member_indices[0])
        best_score = float("inf")
        for candidate in member_indices:
            candidate_orientation = self.orientations[int(candidate)]
            score = 0.0
            for other in member_indices:
                score += candidate_orientation.distance_to(
                    self.orientations[int(other)],
                    symmetry_aware=symmetry_aware,
                )
            if score < best_score:
                best_index = int(candidate)
                best_score = score
        return best_index

    def _segmentation_from_labels(
        self,
        labels: np.ndarray,
        *,
        max_misorientation_deg: float,
        symmetry_aware: bool,
        connectivity: int,
    ) -> GrainSegmentation:
        labels_array = np.asarray(labels, dtype=np.int64)
        if labels_array.shape != (len(self.orientations),):
            raise ValueError("labels must contain one entry per orientation.")
        unique_labels = sorted(int(label) for label in np.unique(labels_array))
        relabeled = np.empty_like(labels_array)
        grains: list[Grain] = []
        for grain_id, old_label in enumerate(unique_labels):
            member_indices = np.flatnonzero(labels_array == old_label).astype(np.int64)
            relabeled[member_indices] = grain_id
            reference_index = self._representative_orientation_index(
                member_indices,
                symmetry_aware=symmetry_aware,
            )
            grains.append(
                Grain(
                    grain_id=grain_id,
                    member_indices=member_indices,
                    mean_coordinate=np.mean(self.coordinates[member_indices], axis=0),
                    reference_orientation_index=reference_index,
                    provenance=self.provenance,
                )
            )
        return GrainSegmentation(
            crystal_map=self,
            labels=relabeled,
            grains=tuple(grains),
            max_misorientation_deg=max_misorientation_deg,
            connectivity=connectivity,
            symmetry_aware=symmetry_aware,
            provenance=self.provenance,
        )

    def segment_grains(
        self,
        *,
        max_misorientation_deg: float,
        symmetry_aware: bool = True,
        connectivity: int = 4,
    ) -> GrainSegmentation:
        if max_misorientation_deg < 0.0:
            raise ValueError("max_misorientation_deg must be non-negative.")
        neighbor_pairs = self.neighbor_graph(connectivity=connectivity, order=1).pairs
        same_phase = self._same_phase_pair_mask(neighbor_pairs)
        neighbor_pairs = neighbor_pairs[same_phase]
        parent = np.arange(len(self.orientations), dtype=np.int64)

        def find(index: int) -> int:
            root = index
            while parent[root] != root:
                root = int(parent[root])
            while parent[index] != index:
                next_index = int(parent[index])
                parent[index] = root
                index = next_index
            return root

        def union(left: int, right: int) -> None:
            left_root = find(left)
            right_root = find(right)
            if left_root != right_root:
                parent[right_root] = left_root

        if neighbor_pairs.size:
            angles_deg = np.rad2deg(
                self._pair_misorientation_rad(neighbor_pairs, symmetry_aware=symmetry_aware)
            )
            valid_pairs = neighbor_pairs[
                np.isfinite(angles_deg) & (angles_deg <= max_misorientation_deg)
            ]
            for left_index, right_index in valid_pairs:
                union(int(left_index), int(right_index))

        component_map: dict[int, list[int]] = {}
        for index in range(len(self.orientations)):
            component_map.setdefault(find(index), []).append(index)

        labels = np.empty(len(self.orientations), dtype=np.int64)
        for grain_id, member_list in enumerate(sorted(component_map.values(), key=lambda x: x[0])):
            member_indices = np.asarray(member_list, dtype=np.int64)
            labels[member_indices] = grain_id

        return self._segmentation_from_labels(
            labels,
            max_misorientation_deg=max_misorientation_deg,
            symmetry_aware=symmetry_aware,
            connectivity=connectivity,
        )
