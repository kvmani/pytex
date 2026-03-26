from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.orientation import Orientation
from pytex.core.orientation import OrientationSet
from pytex.core.provenance import ProvenanceRecord


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
            edge_key = tuple(sorted((segment.left_grain_id, segment.right_grain_id)))
            edge_groups.setdefault(edge_key, []).append((index, segment))
        edges: list[GrainGraphEdge] = []
        for edge_key, members in sorted(edge_groups.items()):
            segment_indices = np.array([member_index for member_index, _ in members], dtype=np.int64)
            segments = [segment for _, segment in members]
            total_length = float(np.sum([segment.length for segment in segments]))
            mean_misorientation = float(np.mean([segment.misorientation_deg for segment in segments]))
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
        node_grain_ids = np.array([grain.grain_id for grain in self.segmentation.grains], dtype=np.int64)
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
        object.__setattr__(self, "segment_indices", np.ascontiguousarray(np.asarray(self.segment_indices, dtype=np.int64)))
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
        deviations = np.ascontiguousarray(deviations.reshape(self.crystal_map.grid_shape))
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
        for left_index, right_index in self.crystal_map.neighbor_pairs(connectivity=self.connectivity):
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
                unique, counts = np.unique(np.asarray(neighbor_labels, dtype=np.int64), return_counts=True)
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
            for left_index, right_index in current.crystal_map.neighbor_pairs(
                connectivity=current.connectivity
            ):
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
    grid_shape: tuple[int, ...] | None = None
    step_sizes: tuple[float, ...] | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] not in (2, 3):
            raise ValueError("CrystalMap.coordinates must have shape (n, 2) or (n, 3).")
        if coordinates.shape[0] != len(self.orientations):
            raise ValueError("CrystalMap coordinates and orientations must have matching lengths.")
        if self.map_frame.domain not in {FrameDomain.MAP, FrameDomain.SPECIMEN}:
            raise ValueError("CrystalMap.map_frame must belong to the map or specimen domain.")
        coordinates = np.ascontiguousarray(coordinates)
        coordinates.setflags(write=False)
        object.__setattr__(self, "coordinates", coordinates)
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

    def neighbor_pairs(self, *, connectivity: int = 4) -> np.ndarray:
        rows, cols = self._require_regular_2d_grid()
        if connectivity not in {4, 8}:
            raise ValueError("connectivity must be either 4 or 8.")
        offsets = [(0, 1), (1, 0)]
        if connectivity == 8:
            offsets.extend([(1, 1), (1, -1)])
        pairs: list[tuple[int, int]] = []
        for row in range(rows):
            for col in range(cols):
                source = row * cols + col
                for drow, dcol in offsets:
                    next_row = row + drow
                    next_col = col + dcol
                    if 0 <= next_row < rows and 0 <= next_col < cols:
                        pairs.append((source, next_row * cols + next_col))
        array = np.asarray(pairs, dtype=np.int64)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        return array

    def kam_neighbor_pairs(self, *, order: int = 1) -> np.ndarray:
        rows, cols = self._require_regular_2d_grid()
        if order <= 0:
            raise ValueError("order must be strictly positive.")
        pairs: list[tuple[int, int]] = []
        for row in range(rows):
            for col in range(cols):
                source = row * cols + col
                for drow in range(-order, order + 1):
                    for dcol in range(-order, order + 1):
                        if drow == 0 and dcol == 0:
                            continue
                        if abs(drow) + abs(dcol) > order:
                            continue
                        next_row = row + drow
                        next_col = col + dcol
                        if not (0 <= next_row < rows and 0 <= next_col < cols):
                            continue
                        target = next_row * cols + next_col
                        if source < target:
                            pairs.append((source, target))
        array = np.asarray(pairs, dtype=np.int64)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        return array

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
        neighbor_pairs = self.kam_neighbor_pairs(order=order)
        sums = np.zeros(len(self.orientations), dtype=np.float64)
        counts = np.zeros(len(self.orientations), dtype=np.int64)
        maxima = np.zeros(len(self.orientations), dtype=np.float64)
        _, cols = self._require_regular_2d_grid()
        for left_index, right_index in neighbor_pairs:
            if segmentation is not None and segmentation.labels[int(left_index)] != segmentation.labels[int(right_index)]:
                continue
            if connectivity == 8:
                row_left, col_left = divmod(int(left_index), cols)
                row_right, col_right = divmod(int(right_index), cols)
                if max(abs(row_left - row_right), abs(col_left - col_right)) > order:
                    continue
            angle_deg = self.orientations[left_index].distance_to(
                self.orientations[right_index],
                symmetry_aware=symmetry_aware,
            )
            angle_deg = float(np.rad2deg(angle_deg))
            if threshold_deg is not None and angle_deg > threshold_deg:
                continue
            sums[left_index] += angle_deg
            sums[right_index] += angle_deg
            counts[left_index] += 1
            counts[right_index] += 1
            maxima[left_index] = max(maxima[left_index], angle_deg)
            maxima[right_index] = max(maxima[right_index], angle_deg)
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
        values = np.ascontiguousarray(values.reshape(self.grid_shape))
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
        neighbor_pairs = self.neighbor_pairs(connectivity=connectivity)
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

        for left_index, right_index in neighbor_pairs:
            angle_deg = float(
                np.rad2deg(
                    self.orientations[int(left_index)].distance_to(
                        self.orientations[int(right_index)],
                        symmetry_aware=symmetry_aware,
                    )
                )
            )
            if angle_deg <= max_misorientation_deg:
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
