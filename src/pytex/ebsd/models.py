from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import csr_matrix, triu

from pytex.core._arrays import normalize_vector
from pytex.core.acquisition import AcquisitionGeometry, CalibrationRecord, MeasurementQuality
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, MillerIndex, Phase
from pytex.core.orientation import (
    Orientation,
    OrientationSet,
    _disorientation_medoid_index,
    _reduced_pair_disorientation_angles,
)
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
    if isinstance(poles, list | tuple):
        return tuple(_coerce_pole(pole, phase=phase) for pole in poles)
    return (_coerce_pole(poles, phase=phase),)


def _coerce_sample_direction_sequence(
    sample_directions: str | ArrayLike | tuple[str | ArrayLike, ...] | list[str | ArrayLike],
    specimen_frame: ReferenceFrame,
) -> tuple[np.ndarray, ...]:
    if isinstance(sample_directions, str):
        return (_specimen_direction_vector(sample_directions, specimen_frame),)
    if isinstance(sample_directions, np.ndarray) and sample_directions.shape == (3,):
        return (_specimen_direction_vector(sample_directions, specimen_frame),)
    if isinstance(sample_directions, list | tuple):
        if len(sample_directions) == 3 and not isinstance(sample_directions[0], str):
            candidate = np.asarray(sample_directions, dtype=np.float64)
            if candidate.shape == (3,):
                return (_specimen_direction_vector(candidate, specimen_frame),)
        return tuple(
            _specimen_direction_vector(direction, specimen_frame) for direction in sample_directions
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


def _freeze_property_channels(
    properties: Mapping[str, ArrayLike] | None,
    *,
    point_count: int,
) -> MappingProxyType[str, np.ndarray]:
    if properties is None:
        return MappingProxyType({})
    frozen: dict[str, np.ndarray] = {}
    for name, values in dict(properties).items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("CrystalMap property channel names must be non-empty strings.")
        array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        if array.shape != (point_count,):
            raise ValueError(
                f"CrystalMap property channel '{name}' must have one value per map point "
                f"(expected shape ({point_count},), got {array.shape})."
            )
        array.setflags(write=False)
        frozen[name] = array
    return MappingProxyType(frozen)


def _rotation_angles_from_matrices(matrices: np.ndarray) -> np.ndarray:
    traces = np.trace(matrices, axis1=1, axis2=2)
    cos_theta = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    return np.asarray(np.arccos(cos_theta), dtype=np.float64)


def _relative_rotation_matrices(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    # Crystal-frame relative rotation (left^T @ right): crystal symmetry then
    # reduces through fixed left/right operator products in
    # _disorientation_angles_from_relative_matrices.
    return np.asarray(
        np.einsum("nji,njk->nik", left, right, optimize=True),
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
    # Shares the core quaternion kernel: one dense product per memory-bounded
    # block, instead of expanding an (n, |S_l|, |S_r|, 3, 3) candidate array
    # that grows without limit in the number of pairs.
    return _reduced_pair_disorientation_angles(relative_matrices, left_ops, right_ops)


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


def _vectorized_hexagonal_grid_pairs(
    coordinates: np.ndarray,
    row_lengths: tuple[int, ...],
    *,
    order: int,
) -> np.ndarray:
    """Unique pairs within ``order`` graph steps on a staggered hexagonal grid."""

    if order <= 0:
        raise ValueError("order must be strictly positive.")
    row_starts = np.concatenate(
        [np.array([0], dtype=np.int64), np.cumsum(row_lengths, dtype=np.int64)]
    )
    pair_blocks: list[np.ndarray] = []
    for row, length in enumerate(row_lengths):
        if length > 1:
            start = int(row_starts[row])
            left = np.arange(start, start + length - 1, dtype=np.int64)
            pair_blocks.append(np.column_stack([left, left + 1]))
    for row in range(len(row_lengths) - 1):
        upper = np.arange(row_starts[row], row_starts[row + 1], dtype=np.int64)
        lower = np.arange(row_starts[row + 1], row_starts[row + 2], dtype=np.int64)
        x_delta = np.abs(coordinates[upper, 0, None] - coordinates[lower, 0][None, :])
        tolerance = max(float(np.max(np.abs(coordinates[:, 0]))) * 1e-12, 1e-12)
        nearest_from_upper = x_delta <= np.min(x_delta, axis=1, keepdims=True) + tolerance
        nearest_from_lower = x_delta <= np.min(x_delta, axis=0, keepdims=True) + tolerance
        upper_positions, lower_positions = np.nonzero(nearest_from_upper | nearest_from_lower)
        cross_row_pairs = np.column_stack(
            [upper[upper_positions], lower[lower_positions]]
        )
        expected_cross_pairs = (
            2 * min(upper.size, lower.size)
            if upper.size != lower.size
            else max(2 * upper.size - 1, 1)
        )
        if cross_row_pairs.shape[0] != expected_cross_pairs:
            raise ValueError(
                "CrystalMap hexagonal coordinates do not form staggered adjacent rows "
                f"{row} and {row + 1}."
            )
        pair_blocks.append(cross_row_pairs)
    if not pair_blocks:
        return np.empty((0, 2), dtype=np.int64)
    first_shell = np.unique(np.concatenate(pair_blocks, axis=0), axis=0)
    if order == 1:
        pairs = np.ascontiguousarray(first_shell, dtype=np.int64)
        pairs.setflags(write=False)
        return pairs

    point_count = int(coordinates.shape[0])
    sources = np.concatenate([first_shell[:, 0], first_shell[:, 1]])
    targets = np.concatenate([first_shell[:, 1], first_shell[:, 0]])
    adjacency = csr_matrix(
        (np.ones(sources.size, dtype=np.int8), (sources, targets)),
        shape=(point_count, point_count),
    )
    frontier = adjacency.copy()
    reachable = adjacency.copy()
    for _ in range(1, order):
        frontier = frontier @ adjacency
        frontier.data[:] = 1
        frontier.setdiag(0)
        frontier.eliminate_zeros()
        reachable = reachable.maximum(frontier)
    upper_triangle = triu(reachable, k=1, format="coo")
    pairs = np.column_stack([upper_triangle.row, upper_triangle.col])
    pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    pairs = np.ascontiguousarray(pairs, dtype=np.int64)
    pairs.setflags(write=False)
    return pairs


@dataclass(frozen=True, slots=True)
class CrystalMapPhase:
    """One phase declared in a crystal map.

    Purpose
    -------
    Vendor files routinely name a phase and give its point group without a
    full crystal structure. This type accommodates that: symmetry is
    required, a full :class:`~pytex.core.lattice.Phase` is optional, and
    either is enough to resolve the crystal frame.

    Attributes
    ----------
    phase_id : int
        Non-negative identifier, as used in the file.
    name : str
        Canonical phase name.
    symmetry : SymmetrySpec
        Required, and its reference frame must be set — phase-resolved maps
        cannot reduce orientations without it.
    phase : Phase, optional
        The full structure, when available. Its name and symmetry are checked
        against this entry rather than assumed to agree.
    aliases : tuple of str
        Alternative names, for matching against other systems' naming.
    provenance : ProvenanceRecord, optional
    """

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
        """The crystal-domain reference frame of this map phase.

        Taken from the full :class:`~pytex.core.lattice.Phase` when one is
        attached, otherwise from the symmetry specification. A map phase can
        exist without a full phase — vendor files often name a phase and give
        its point group but no lattice — so this property is the single place
        that resolves which frame applies.
        """

        if self.phase is not None:
            return self.phase.crystal_frame
        return cast(ReferenceFrame, self.symmetry.reference_frame)

    @property
    def point_group(self) -> str:
        """The Hermann-Mauguin point-group symbol of this map phase."""

        return self.symmetry.point_group

    def matches(self, selector: int | str | Phase | CrystalMapPhase) -> bool:
        """Whether this map phase is the one a selector refers to.

        Purpose
        -------
        EBSD workflows identify a phase in whichever way is convenient at the
        call site — by the integer phase id in the file, by name, by a
        :class:`~pytex.core.lattice.Phase`, or by another
        :class:`CrystalMapPhase`. This resolves all four to one boolean so
        callers need not branch on the selector type.

        Parameters
        ----------
        selector : int, str, Phase, or CrystalMapPhase
            An integer matches ``phase_id``; a string matches the canonical
            name or any registered alias; a ``Phase`` matches the attached phase
            when there is one, and otherwise falls back to name/alias matching;
            a ``CrystalMapPhase`` matches on ``phase_id``.
        """

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
    """Which measurement points count as neighbours, and how far apart they are.

    Purpose
    -------
    The neighbourhood definition that KAM, grain segmentation, boundary
    extraction, and smoothing all share, so the definition lives in one place
    rather than being re-derived per metric.

    Attributes
    ----------
    pairs : np.ndarray
        ``(m, 2)`` unique unordered index pairs.
    distances : np.ndarray
        ``(m,)`` distances between the paired points.
    connectivity : int
        ``4`` (square edges), ``6`` (hexagonal first shell), or ``8``
        (square edges and corners).
    order : int
        Neighbour shell; ``1`` is nearest neighbours.
    mode : str
        How the pairs were built: ``"regular_grid"`` for the vectorized grid
        construction, ``"hexagonal_grid"`` for staggered logical rows, or
        ``"coordinate_radius"`` for the distance-based fallback used on
        irregular point sets.
    max_distance : float, optional
        Radius cut-off, when the coordinate-radius path was used.
    """

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
        if self.connectivity not in {4, 6, 8}:
            raise ValueError("CoordinateNeighborGraph.connectivity must be 4, 6, or 8.")
        if self.order <= 0:
            raise ValueError("CoordinateNeighborGraph.order must be strictly positive.")
        if self.max_distance is not None and self.max_distance <= 0.0:
            raise ValueError("CoordinateNeighborGraph.max_distance must be positive when provided.")
        if self.mode not in {"regular_grid", "hexagonal_grid", "coordinate_radius"}:
            raise ValueError(
                "CoordinateNeighborGraph.mode must be 'regular_grid', 'hexagonal_grid', "
                "or 'coordinate_radius'."
            )
        object.__setattr__(self, "pairs", pairs)
        object.__setattr__(self, "distances", distances)


@dataclass(frozen=True, slots=True)
class TextureReport:
    """A complete texture characterization computed under one set of choices.

    Purpose
    -------
    Bundles the ODF with its pole figures and inverse pole figures, all
    produced from the same weighting, kernel, and symmetry assumptions — so
    the parts of a texture description cannot drift apart in conventions.

    Attributes
    ----------
    odf : ODF
    pole_figures : tuple of PoleFigure
    inverse_pole_figures : tuple of InversePoleFigure
    odf_figure, pole_figure_figures, inverse_pole_figure_figures : optional
        Rendered Matplotlib figures, present only when plotting was
        requested; the numerical products never depend on them.
    provenance : ProvenanceRecord, optional
    """

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
class FittedEllipse:
    """Second-moment equivalent ellipse of a grain's member pixel positions.

    ``semi_axes`` are ordered major-first (``a >= b``), computed from the
    covariance eigenvalues as ``2 * sqrt(lambda)`` so a uniformly filled
    ellipse recovers its own semi-axes. ``angle_deg`` is the major-axis
    orientation measured counter-clockwise from the map x-axis in [0, 180).
    """

    grain_id: int
    centroid: np.ndarray
    semi_axes: tuple[float, float]
    angle_deg: float
    aspect_ratio: float

    def __post_init__(self) -> None:
        centroid = np.ascontiguousarray(np.asarray(self.centroid, dtype=np.float64))
        centroid.setflags(write=False)
        object.__setattr__(self, "centroid", centroid)
        if self.semi_axes[0] < self.semi_axes[1]:
            raise ValueError("FittedEllipse.semi_axes must be ordered major-axis first.")


@dataclass(frozen=True, slots=True)
class Grain:
    """One segmented grain: its member points and its representative orientation.

    Attributes
    ----------
    grain_id : int
    member_indices : np.ndarray
        Indices of the measurement points in this grain; non-empty.
    mean_coordinate : np.ndarray
        Centroid of the member points.
    reference_orientation_index : int
        The member point chosen as representative. Must belong to
        ``member_indices``; this is the reference GROD is measured against,
        and it is a measured orientation rather than an average.
    provenance : ProvenanceRecord, optional
    """

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
        """Number of measurement points belonging to this grain.

        A pixel count, not a physical area. For area or diameter in specimen
        units use :meth:`GrainSegmentation.grain_areas` or
        :meth:`GrainSegmentation.grain_equivalent_diameters`, which apply the
        map step sizes.
        """

        return int(self.member_indices.size)


@dataclass(frozen=True, slots=True)
class GrainBoundarySegment:
    """One pixel face separating two grains, with its misorientation.

    Attributes
    ----------
    left_index, right_index : int
        The two measurement points; must be distinct.
    left_grain_id, right_grain_id : int
        Their grains; must be distinct, since a segment by definition
        straddles a boundary.
    misorientation_deg : float
        Non-negative; symmetry-reduced when the segmentation was.
    length : float
        Face length in map coordinate units.
    midpoint : np.ndarray
        Face centre, for plotting the boundary network.
    provenance : ProvenanceRecord, optional
    """

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
        """Label this segment ``"high_angle"`` or ``"low_angle"``.

        The 15-degree default threshold is the conventional dividing line
        between low-angle (dislocation-wall) and high-angle boundaries in the
        Read-Shockley picture; it is a convention, not a physical constant, and
        is exposed so it can be set to whatever a study uses.
        """

        if high_angle_threshold_deg < 0.0:
            raise ValueError("high_angle_threshold_deg must be non-negative.")
        return "high_angle" if self.misorientation_deg >= high_angle_threshold_deg else "low_angle"


@dataclass(frozen=True, slots=True)
class GrainBoundaryNetwork:
    """Every boundary segment of a segmentation, with classification thresholds.

    Purpose
    -------
    The raw material for boundary-character statistics: segment
    misorientations, lengths, and the high-angle fraction. Collapse it with
    :meth:`grain_graph` when grain-scale adjacency is what is wanted.

    Attributes
    ----------
    segmentation : GrainSegmentation
        The segmentation the segments came from.
    segments : tuple of GrainBoundarySegment
    min_misorientation_deg : float
        The threshold below which segments were discarded.
    high_angle_threshold_deg : float
        Dividing line between low- and high-angle boundaries; 15 degrees by
        convention, and a convention rather than a physical constant.
    provenance : ProvenanceRecord, optional
    """

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
        """Number of boundary segments in the network."""

        return len(self.segments)

    @property
    def mean_misorientation_deg(self) -> float:
        """Mean segment misorientation in degrees; ``0.0`` for an empty network.

        Unweighted over segments. Because segments are pixel-face sized, this
        is close to a length-weighted mean on a regular grid, but it is not the
        same statistic on an irregular one.
        """

        if not self.segments:
            return 0.0
        return float(np.mean([segment.misorientation_deg for segment in self.segments]))

    @property
    def total_length(self) -> float:
        """Summed length of all boundary segments, in map coordinate units."""

        return float(np.sum([segment.length for segment in self.segments]))

    @property
    def high_angle_count(self) -> int:
        """Number of segments at or above the network's high-angle threshold."""

        return int(
            sum(
                segment.classify(high_angle_threshold_deg=self.high_angle_threshold_deg)
                == "high_angle"
                for segment in self.segments
            )
        )

    def grain_graph(self) -> GrainGraph:
        """Collapse the segment list into a grain-adjacency graph.

        Purpose
        -------
        The segment network is per pixel face; most grain-scale reasoning —
        neighbour queries, CSL classification, parent-grain reconstruction —
        wants one edge per grain pair instead. This groups segments by grain
        pair and summarizes each pair's total length, mean misorientation, and
        high-angle fraction.

        Returns
        -------
        GrainGraph
            Nodes are grain ids, edges are grain pairs, in sorted pair order so
            the graph is reproducible.
        """

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

    def _segment_misorientation_matrices(self) -> np.ndarray:
        matrices = self.segmentation.crystal_map.orientations.as_matrices()
        left = np.array(
            [matrices[segment.left_index] for segment in self.segments], dtype=np.float64
        )
        right = np.array(
            [matrices[segment.right_index] for segment in self.segments], dtype=np.float64
        )
        if left.size == 0:
            return np.empty((0, 3, 3), dtype=np.float64)
        return _relative_rotation_matrices(left, right)

    def _cubic_operators(self) -> np.ndarray:
        crystal_map = self.segmentation.crystal_map
        entries = crystal_map.resolved_phase_entries
        if len(entries) != 1:
            raise ValueError("CSL classification currently supports single-phase cubic maps only.")
        symmetry = entries[0].symmetry
        if symmetry is None or symmetry.point_group not in {"m-3m", "m-3", "432", "23"}:
            raise ValueError(
                "CSL classification requires a cubic crystal symmetry "
                "(point group 23, m-3, 432, or m-3m)."
            )
        return symmetry.operators

    def classify_csl(
        self,
        *,
        theta0_deg: float = 15.0,
        include_sigma1: bool = False,
    ) -> tuple[Any, ...]:
        """Classify each boundary segment against the cubic CSL registry.

        Returns one entry per segment: a ``CSLMatch`` (with the assigned Sigma
        and deviation) or ``None`` when no CSL type fits within the Brandon
        criterion. Requires a single-phase cubic map.
        """

        from pytex.ebsd.csl import classify_misorientations

        matrices = self._segment_misorientation_matrices()
        if matrices.shape[0] == 0:
            return ()
        return tuple(
            classify_misorientations(
                matrices,
                operators=self._cubic_operators(),
                theta0_deg=theta0_deg,
                include_sigma1=include_sigma1,
            )
        )

    def csl_fraction(self, sigma: int, *, theta0_deg: float = 15.0) -> float:
        """Fraction of boundary length classified as the given CSL Sigma."""

        matches = self.classify_csl(theta0_deg=theta0_deg)
        if not matches:
            return 0.0
        total_length = self.total_length
        if total_length <= 0.0:
            return 0.0
        matched_length = sum(
            segment.length
            for segment, match in zip(self.segments, matches, strict=True)
            if match is not None and match.sigma == sigma
        )
        return float(matched_length / total_length)

    def select_csl(
        self, sigma: int, *, theta0_deg: float = 15.0
    ) -> tuple[GrainBoundarySegment, ...]:
        """Return the boundary segments classified as the given CSL Sigma."""

        matches = self.classify_csl(theta0_deg=theta0_deg)
        return tuple(
            segment
            for segment, match in zip(self.segments, matches, strict=True)
            if match is not None and match.sigma == sigma
        )

    def merge_by_csl(self, sigma: int, *, theta0_deg: float = 15.0) -> GrainSegmentation:
        """Merge grains joined by the given CSL boundary type into parent grains.

        Grains connected across a boundary classified as the requested CSL Sigma
        are unioned; a new `GrainSegmentation` is returned in which each such
        connected group carries a single label (an MTEX-style twin-merged parent
        grain). The underlying crystal map is unchanged.
        """

        segmentation = self.segmentation
        grain_ids = [grain.grain_id for grain in segmentation.grains]
        parent = {grain_id: grain_id for grain_id in grain_ids}

        def find(grain_id: int) -> int:
            root = grain_id
            while parent[root] != root:
                root = parent[root]
            while parent[grain_id] != root:
                parent[grain_id], grain_id = root, parent[grain_id]
            return root

        matches = self.classify_csl(theta0_deg=theta0_deg)
        for segment, match in zip(self.segments, matches, strict=True):
            if match is not None and match.sigma == sigma:
                left_root = find(segment.left_grain_id)
                right_root = find(segment.right_grain_id)
                if left_root != right_root:
                    parent[right_root] = left_root
        merged_labels = np.array(
            [find(int(label)) for label in segmentation.labels], dtype=np.int64
        )
        return segmentation.crystal_map._segmentation_from_labels(
            merged_labels,
            max_misorientation_deg=segmentation.max_misorientation_deg,
            symmetry_aware=segmentation.symmetry_aware,
            connectivity=segmentation.connectivity,
        )

    def twin_merge(self, *, theta0_deg: float = 15.0) -> GrainSegmentation:
        """Merge Sigma3-twin-related grains into parent grains."""

        return self.merge_by_csl(3, theta0_deg=theta0_deg)


@dataclass(frozen=True, slots=True)
class GrainGraphEdge:
    """The adjacency between two grains, summarizing their shared boundary.

    Attributes
    ----------
    left_grain_id, right_grain_id : int
        The connected grains; must be distinct.
    segment_indices : np.ndarray
        The boundary segments forming this edge; non-empty.
    total_length : float
        Shared boundary length. This is the natural weight for voting
        algorithms — a long shared boundary is stronger evidence than a
        single shared pixel face.
    mean_misorientation_deg : float
    high_angle_fraction : float
        Fraction of the shared boundary above the high-angle threshold, in
        ``[0, 1]``.
    provenance : ProvenanceRecord, optional
    """

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
        """The ``(left_grain_id, right_grain_id)`` pair this edge connects."""

        return (self.left_grain_id, self.right_grain_id)


@dataclass(frozen=True, slots=True)
class GrainGraph:
    """The grain-adjacency graph of a segmentation.

    Purpose
    -------
    Grain-scale connectivity: one edge per grain pair rather than one per
    pixel face. This is the structure region-growing algorithms iterate on —
    parent-grain reconstruction, twin merging, small-grain absorption.

    Attributes
    ----------
    segmentation : GrainSegmentation
    edges : tuple of GrainGraphEdge
        In sorted grain-pair order, so the graph is reproducible.
    node_grain_ids : np.ndarray
        Grain ids in adjacency-matrix row order. Note that row ``i``
        corresponds to ``node_grain_ids[i]``, not to grain id ``i``.
    high_angle_threshold_deg : float
    provenance : ProvenanceRecord, optional
    """

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
        """Number of grain-pair edges in the graph."""

        return len(self.edges)

    @property
    def adjacency_matrix(self) -> np.ndarray:
        """Symmetric 0/1 adjacency matrix over :attr:`node_grain_ids`.

        Row and column ``i`` correspond to ``node_grain_ids[i]``, *not* to grain
        id ``i``. Returned read-only.
        """

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
        """Sorted grain ids adjacent to ``grain_id``.

        Returns an empty array for an isolated grain. This is the primitive that
        region-growing algorithms — parent-grain reconstruction, twin merging,
        small-grain absorption — iterate on.
        """

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
    """A crystal map partitioned into grains, with the settings that produced it.

    Purpose
    -------
    The bridge from point-wise orientations to microstructural statistics:
    grain sizes, shapes, boundaries, averages, and the local-misorientation
    family (GROD, GOS, GAM) all derive from here.

    The segmentation parameters are stored, not just applied, so every
    derived metric inherits the same conventions and the result stays
    reproducible.

    Attributes
    ----------
    crystal_map : CrystalMap
        The map segmented; derived metrics query it directly.
    labels : np.ndarray
        One grain id per measurement point.
    grains : tuple of Grain
    max_misorientation_deg : float
        The grain-boundary criterion used. It determines whether subgrains
        are resolved as separate grains, so it must be reported with any
        grain-size result.
    connectivity : int
        ``4`` or ``8`` neighbourhood.
    symmetry_aware : bool
        Whether disorientation rather than raw rotation angle was used.
    provenance : ProvenanceRecord, optional
    """

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
        if self.connectivity not in {4, 6, 8}:
            raise ValueError("GrainSegmentation.connectivity must be 4, 6, or 8.")
        labels = np.ascontiguousarray(labels)
        labels.setflags(write=False)
        object.__setattr__(self, "labels", labels)
        object.__setattr__(self, "grains", tuple(self.grains))

    @property
    def label_grid(self) -> np.ndarray:
        """Per-pixel grain labels reshaped to the map grid.

        Rectangular maps return their native shape. Hexagonal maps return a
        padded ``(n_rows, max_row_length)`` array with ``-1`` in positions that
        are not measurement points. The non-negative values are grain ids.
        """

        if self.crystal_map.grid_kind == "hexagonal":
            assert self.crystal_map.row_lengths is not None
            labels = np.full(
                (len(self.crystal_map.row_lengths), max(self.crystal_map.row_lengths)),
                -1,
                dtype=np.int64,
            )
            offset = 0
            for row, length in enumerate(self.crystal_map.row_lengths):
                labels[row, :length] = self.labels[offset : offset + length]
                offset += length
        else:
            rows, cols = self.crystal_map._require_regular_2d_grid()
            labels = self.labels.reshape((rows, cols))
        labels = np.ascontiguousarray(labels)
        labels.setflags(write=False)
        return labels

    def reference_orientation(self, grain: Grain) -> Orientation:
        """The orientation of a grain's representative measurement point.

        This is the *measured* orientation at the point the segmentation chose
        as representative, not an average. For the averaged quantity use
        :meth:`grain_mean_orientation`; GROD is measured against this reference.
        """

        return self.crystal_map.orientations[grain.reference_orientation_index]

    def grod_map_deg(self) -> np.ndarray:
        """Grain Reference Orientation Deviation map, in degrees.

        Purpose
        -------
        Per pixel, the misorientation between that point and its own grain's
        reference orientation. GROD reveals intragranular orientation gradients
        — the signature of stored plastic strain, subgrain structure, and
        recovery — which a grain-average number hides.

        Returns
        -------
        np.ndarray
            Degrees, read-only. Rectangular maps return ``(rows, cols)``;
            hexagonal and unstructured maps return one value per point.
            Symmetry-aware iff the segmentation was built that way.

        See Also
        --------
        gos_map_deg : One number per grain instead of per pixel.
        kernel_average_misorientation_deg : Local (neighbour) rather than
            grain-referenced deviation.
        """

        point_count = len(self.crystal_map.orientations)
        # Per-point reference-orientation index: each point takes its grain's
        # representative orientation. Vectorised over all points at once.
        reference_index = np.empty(point_count, dtype=np.int64)
        for grain in self.grains:
            reference_index[grain.member_indices] = grain.reference_orientation_index
        matrices = self.crystal_map.orientations.as_matrices()
        relative = _relative_rotation_matrices(matrices, matrices[reference_index])
        if self.symmetry_aware:
            angles_rad = _disorientation_angles_from_relative_matrices(
                relative,
                left_symmetry=self.crystal_map.orientations.symmetry,
                right_symmetry=self.crystal_map.orientations.symmetry,
            )
        else:
            angles_rad = _rotation_angles_from_matrices(relative)
        deviations = np.rad2deg(angles_rad)
        if self.crystal_map.grid_shape is not None:
            rows, cols = self.crystal_map._require_regular_2d_grid()
            deviations = deviations.reshape((rows, cols))
        deviations = np.ascontiguousarray(deviations)
        deviations.setflags(write=False)
        return deviations

    def grain_mean_orientation(self, grain: Grain) -> Orientation:
        """Symmetry-aware mean orientation of one grain's member points.

        Uses quaternion eigenvector averaging with per-member symmetry-branch
        selection; see :meth:`~pytex.core.orientation.OrientationSet.mean_orientation`.
        """

        return self.crystal_map.orientations.subset(grain.member_indices).mean_orientation()

    def grain_mean_orientations(self) -> dict[int, Orientation]:
        """Mean orientation of every grain, keyed by grain id."""

        return {grain.grain_id: self.grain_mean_orientation(grain) for grain in self.grains}

    def grain_orientation_spread_deg(self) -> dict[int, float]:
        """Grain Orientation Spread (GOS): mean deviation of members from the grain mean."""

        spreads: dict[int, float] = {}
        for grain in self.grains:
            member_set = self.crystal_map.orientations.subset(grain.member_indices)
            mean = member_set.mean_orientation()
            angles = member_set.spread_angles_deg(
                reference=mean,
                symmetry_aware=self.symmetry_aware,
            )
            spreads[grain.grain_id] = float(np.mean(angles))
        return spreads

    def _broadcast_grain_values_to_points(self, values: dict[int, float]) -> np.ndarray:
        per_point = np.zeros(len(self.crystal_map.orientations), dtype=np.float64)
        for grain in self.grains:
            per_point[grain.member_indices] = values[grain.grain_id]
        return per_point

    def _reshape_to_grid_if_regular(self, per_point: np.ndarray) -> np.ndarray:
        if self.crystal_map.grid_shape is not None and len(self.crystal_map.grid_shape) == 2:
            rows, cols = self.crystal_map._require_regular_2d_grid()
            per_point = per_point.reshape((rows, cols))
        per_point = np.ascontiguousarray(per_point)
        per_point.setflags(write=False)
        return per_point

    def gos_map_deg(self) -> np.ndarray:
        """Grain Orientation Spread broadcast to every pixel of its grain.

        Each pixel carries its grain's GOS value, so the map shows which grains
        are deformed rather than where within a grain the deformation sits.
        Returned as a grid for a regular 2-D map and per point otherwise.
        """

        return self._reshape_to_grid_if_regular(
            self._broadcast_grain_values_to_points(self.grain_orientation_spread_deg())
        )

    def grain_average_misorientation_deg(self) -> dict[int, float]:
        """Grain Average Misorientation (GAM): mean intragranular KAM per grain."""

        kam = np.asarray(
            self.crystal_map.kernel_average_misorientation_deg(
                symmetry_aware=self.symmetry_aware,
                connectivity=self.connectivity,
                segmentation=self,
            ),
            dtype=np.float64,
        ).ravel()
        return {grain.grain_id: float(np.mean(kam[grain.member_indices])) for grain in self.grains}

    def gam_map_deg(self) -> np.ndarray:
        """Grain Average Misorientation broadcast to every pixel of its grain.

        GAM averages the *local* (kernel) misorientation inside a grain, so it
        responds to short-wavelength gradients — geometrically necessary
        dislocation content — where GOS responds to the total spread. The two
        differ sharply for a grain with a single sharp subgrain wall.
        """

        return self._reshape_to_grid_if_regular(
            self._broadcast_grain_values_to_points(self.grain_average_misorientation_deg())
        )

    def grain_equivalent_diameters(self) -> dict[int, float]:
        """Equivalent circular diameter per grain from member pixel area."""

        if self.crystal_map.step_sizes is None or len(self.crystal_map.step_sizes) != 2:
            raise ValueError(
                "grain_equivalent_diameters requires a CrystalMap with 2-D step_sizes."
            )
        dx, dy = self.crystal_map.step_sizes
        pixel_area = float(dx) * float(dy)
        return {
            grain.grain_id: float(2.0 * np.sqrt(grain.size * pixel_area / np.pi))
            for grain in self.grains
        }

    def grain_sizes(self) -> dict[int, int]:
        """Member pixel count of every grain, keyed by grain id."""

        return {grain.grain_id: grain.size for grain in self.grains}

    def _fitted_ellipse(self, grain: Grain) -> FittedEllipse:
        coordinates = self.crystal_map.coordinates[grain.member_indices, :2]
        centroid = np.mean(coordinates, axis=0)
        if grain.size < 2:
            return FittedEllipse(
                grain_id=grain.grain_id,
                centroid=centroid,
                semi_axes=(0.0, 0.0),
                angle_deg=0.0,
                aspect_ratio=1.0,
            )
        covariance = np.cov(coordinates, rowvar=False, bias=True)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        # eigh returns ascending eigenvalues; take the major axis last.
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = np.clip(eigenvalues[order], 0.0, None)
        major_vector = eigenvectors[:, order[0]]
        semi_major = float(2.0 * np.sqrt(eigenvalues[0]))
        semi_minor = float(2.0 * np.sqrt(eigenvalues[1]))
        angle = float(np.degrees(np.arctan2(major_vector[1], major_vector[0])) % 180.0)
        aspect_ratio = (
            float(np.sqrt(eigenvalues[0] / eigenvalues[1]))
            if eigenvalues[1] > 1e-12
            else float("inf")
        )
        return FittedEllipse(
            grain_id=grain.grain_id,
            centroid=centroid,
            semi_axes=(semi_major, semi_minor),
            angle_deg=angle,
            aspect_ratio=aspect_ratio,
        )

    def grain_fitted_ellipse(self, grain: Grain) -> FittedEllipse:
        """Best-fit ellipse of one grain's member-pixel cloud.

        Fitted from the eigen-decomposition of the pixel-coordinate covariance,
        with semi-axes set to twice the eigenvalue square roots. This is the
        standard second-moment shape descriptor: it captures elongation and
        alignment, and it is not a boundary fit, so a concave grain is still
        described by a convex ellipse.
        """

        return self._fitted_ellipse(grain)

    def grain_fitted_ellipses(self) -> dict[int, FittedEllipse]:
        """Best-fit ellipse of every grain, keyed by grain id."""

        return {grain.grain_id: self._fitted_ellipse(grain) for grain in self.grains}

    def grain_aspect_ratios(self) -> dict[int, float]:
        """Major-to-minor axis ratio of every grain's fitted ellipse.

        The elongation measure behind grain-morphology statistics for rolled and
        drawn material. Degenerate (single-line) grains report infinity.
        """

        return {grain.grain_id: self._fitted_ellipse(grain).aspect_ratio for grain in self.grains}

    def grain_shape_orientations_deg(self) -> dict[int, float]:
        """Angle of each grain's major axis, in degrees within ``[0, 180)``.

        Measured from the map x-axis. Together with
        :meth:`grain_aspect_ratios` this gives the morphological-texture
        counterpart of crystallographic texture.
        """

        return {grain.grain_id: self._fitted_ellipse(grain).angle_deg for grain in self.grains}

    def _grid_step_sizes(self) -> tuple[float, float]:
        if self.crystal_map.step_sizes is not None and len(self.crystal_map.step_sizes) == 2:
            return float(self.crystal_map.step_sizes[0]), float(self.crystal_map.step_sizes[1])
        coordinates = self.crystal_map.coordinates
        x_values = np.unique(coordinates[:, 0])
        y_values = np.unique(coordinates[:, 1])
        dx = float(np.min(np.diff(x_values))) if x_values.size > 1 else 1.0
        dy = float(np.min(np.diff(y_values))) if y_values.size > 1 else 1.0
        return dx, dy

    def grain_perimeters(self) -> dict[int, float]:
        """True staircase perimeter of each grain on a regular 2-D grid.

        Sums the lengths of the pixel faces on each grain's boundary, including
        faces on the map edge. Horizontal (column) neighbours share a vertical
        face of length ``dy``; vertical (row) neighbours share a face of length
        ``dx``; step sizes are honoured for rectangular grids.
        """

        if self.crystal_map.grid_kind == "hexagonal":
            raise ValueError(
                "grain_perimeters is currently defined only for rectangular pixel faces; "
                "hexagonal center sampling requires a declared cell-boundary model."
            )
        label_grid = self.label_grid
        rows, cols = label_grid.shape
        dx, dy = self._grid_step_sizes()
        # A face is a boundary when the neighbour is off-grid or a different
        # grain. Map-edge faces default to True; interior faces compare shifted
        # label slices. Column neighbours give dy faces, row neighbours dx faces.
        left_boundary = np.ones((rows, cols), dtype=bool)
        left_boundary[:, 1:] = label_grid[:, 1:] != label_grid[:, :-1]
        right_boundary = np.ones((rows, cols), dtype=bool)
        right_boundary[:, :-1] = label_grid[:, :-1] != label_grid[:, 1:]
        up_boundary = np.ones((rows, cols), dtype=bool)
        up_boundary[1:, :] = label_grid[1:, :] != label_grid[:-1, :]
        down_boundary = np.ones((rows, cols), dtype=bool)
        down_boundary[:-1, :] = label_grid[:-1, :] != label_grid[1:, :]
        per_cell = dy * (left_boundary.astype(np.float64) + right_boundary) + dx * (
            up_boundary.astype(np.float64) + down_boundary
        )
        grain_count = len(self.grains)
        totals = np.bincount(label_grid.ravel(), weights=per_cell.ravel(), minlength=grain_count)
        return {grain.grain_id: float(totals[grain.grain_id]) for grain in self.grains}

    def grain_areas(self) -> dict[int, float]:
        """Area of each grain (member pixel count times the pixel area)."""

        dx, dy = self._grid_step_sizes()
        pixel_area = dx * dy
        return {grain.grain_id: grain.size * pixel_area for grain in self.grains}

    def grain_shape_factors(self) -> dict[int, float]:
        """Shape factor ``P / (2 sqrt(pi A))`` (1 for a circle, larger otherwise)."""

        perimeters = self.grain_perimeters()
        areas = self.grain_areas()
        return {
            grain_id: float(perimeters[grain_id] / (2.0 * np.sqrt(np.pi * areas[grain_id])))
            for grain_id in perimeters
        }

    def grain_bounding_boxes(self) -> dict[int, tuple[float, float]]:
        """Axis-aligned (width, height) of each grain's member-pixel extent."""

        boxes: dict[int, tuple[float, float]] = {}
        for grain in self.grains:
            coordinates = self.crystal_map.coordinates[grain.member_indices, :2]
            extent = np.ptp(coordinates, axis=0)
            boxes[grain.grain_id] = (float(extent[0]), float(extent[1]))
        return boxes

    def boundary_network(
        self,
        *,
        min_misorientation_deg: float = 0.0,
        high_angle_threshold_deg: float = 15.0,
    ) -> GrainBoundaryNetwork:
        """Build the grain-boundary segment network of this segmentation.

        Purpose
        -------
        Extract every pixel face that straddles two grains, with its
        misorientation, length, and midpoint — the raw material for
        boundary-character statistics, CSL classification, and boundary
        plotting.

        Parameters
        ----------
        min_misorientation_deg : float
            Discard segments below this misorientation. Use it to suppress the
            noise floor of the indexing.
        high_angle_threshold_deg : float
            The threshold the resulting network uses to classify segments;
            15 degrees by convention.

        Returns
        -------
        GrainBoundaryNetwork
            Segments only between points of the same phase; interphase
            interfaces are not boundary segments in this model.
        """

        if min_misorientation_deg < 0.0:
            raise ValueError("min_misorientation_deg must be non-negative.")
        segments: list[GrainBoundarySegment] = []
        neighbor_pairs = self.crystal_map.neighbor_pairs(connectivity=self.connectivity)
        neighbor_pairs = neighbor_pairs[self.crystal_map._same_phase_pair_mask(neighbor_pairs)]
        # Keep only pairs that straddle a grain boundary, then compute their
        # misorientations and geometry vectorised in one shot.
        labels = self.labels
        boundary_mask = labels[neighbor_pairs[:, 0]] != labels[neighbor_pairs[:, 1]]
        boundary_pairs = neighbor_pairs[boundary_mask]
        if boundary_pairs.size:
            # Mirror Orientation.distance_to exactly (reduction by the set-level
            # symmetry) so results are identical to the previous scalar path.
            matrices = self.crystal_map.orientations.as_matrices()
            relative = _relative_rotation_matrices(
                matrices[boundary_pairs[:, 0]], matrices[boundary_pairs[:, 1]]
            )
            if self.symmetry_aware:
                angles_rad = _disorientation_angles_from_relative_matrices(
                    relative,
                    left_symmetry=self.crystal_map.orientations.symmetry,
                    right_symmetry=self.crystal_map.orientations.symmetry,
                )
            else:
                angles_rad = _rotation_angles_from_matrices(relative)
            misorientation_deg = np.rad2deg(angles_rad)
            coordinates = self.crystal_map.coordinates
            left_coordinates = coordinates[boundary_pairs[:, 0]]
            right_coordinates = coordinates[boundary_pairs[:, 1]]
            lengths = np.linalg.norm(left_coordinates - right_coordinates, axis=1)
            midpoints = 0.5 * (left_coordinates + right_coordinates)
            keep = misorientation_deg >= min_misorientation_deg
            for position in np.flatnonzero(keep):
                left_index = int(boundary_pairs[position, 0])
                right_index = int(boundary_pairs[position, 1])
                segments.append(
                    GrainBoundarySegment(
                        left_index=left_index,
                        right_index=right_index,
                        left_grain_id=int(labels[left_index]),
                        right_grain_id=int(labels[right_index]),
                        misorientation_deg=float(misorientation_deg[position]),
                        length=float(lengths[position]),
                        midpoint=midpoints[position],
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
        """The grain-adjacency graph of this segmentation.

        Convenience for ``boundary_network(...).grain_graph()``; see
        :meth:`GrainBoundaryNetwork.grain_graph`.
        """

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
        """Reassign isolated points to the label most of their neighbours carry.

        Purpose
        -------
        Clean up salt-and-pepper labelling left by mis-indexed points, without
        moving real boundaries. A point changes label only when at least
        ``min_neighbor_votes`` of its neighbours agree on a different one, so a
        genuine boundary — where votes are split — is left alone.

        Parameters
        ----------
        iterations : int
            Number of smoothing sweeps. Each sweep uses the labels from the end
            of the previous one, so the operation is deterministic.
        min_neighbor_votes : int
            Votes required to overrule a point's own label. Raising it makes the
            filter more conservative.

        Returns
        -------
        GrainSegmentation
            A new segmentation on the same map, re-derived from the smoothed
            labels. Requires a regular 2-D grid.
        """

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
        """Absorb grains below a size threshold into their best neighbour.

        Purpose
        -------
        Remove sub-resolution "grains" that are really indexing noise, so that
        grain-size statistics are not dominated by single-pixel artefacts.

        Method
        ------
        Repeatedly takes the smallest grain below ``min_size`` and merges it
        into the neighbour with which it shares the most boundary faces, using
        the lower mean misorientation and then the lower grain id to break ties.
        Deterministic by construction. Merging can push a grain over the
        threshold, which is why the default runs to a fixed point.

        Parameters
        ----------
        min_size : int
            Minimum member-pixel count a grain must have to survive.
        until_stable : bool
            Repeat until no grain is below the threshold (default). ``False``
            performs a single merge.
        max_iterations : int, optional
            Hard cap on merge steps, as a guard on pathological maps.

        Returns
        -------
        GrainSegmentation
            A new segmentation with the surviving grains relabelled.
        """

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
    """An EBSD orientation map: one orientation per measurement point.

    Purpose
    -------
    The central EBSD object. It holds the measured orientations together with
    their spatial coordinates, phase assignments, grid geometry, acquisition
    context, and any auxiliary quality channels — everything needed to
    interpret the orientations as a microstructure rather than as a bare
    list.

    Immutable by construction: derived maps (a phase subset, a filtered
    selection, an added property channel) are produced as new maps, so a
    downstream calculation cannot change the data another calculation is
    reading.

    Grid versus graph mode
    ----------------------
    With ``grid_shape`` set, the map is a regular rectangular raster and
    grid-shaped outputs are available. A staggered hexagonal scan instead uses
    ``grid_kind="hexagonal"`` and ``row_lengths`` because alternating rows are
    ragged. Without either topology the map operates in coordinate-graph mode.

    Attributes
    ----------
    coordinates : np.ndarray
        ``(n, 2)`` or ``(n, 3)`` point positions in the map frame.
    orientations : OrientationSet
        One orientation per point, in matching order.
    map_frame : ReferenceFrame
        Must belong to the map or specimen domain.
    phase_entries : tuple of CrystalMapPhase
        Declared phases. A multiphase map additionally requires per-point
        ``phase_ids`` and forbids phase and symmetry on the shared
        ``OrientationSet`` — so one phase's symmetry can never be applied to
        another phase's points.
    phase_ids : np.ndarray, optional
        Per-point phase assignment.
    grid_shape : tuple of int, optional
        Rectangular raster shape.
    grid_kind : {"square", "hexagonal"}, optional
        Logical scan topology. Existing rectangular maps infer ``"square"``
        from ``grid_shape``.
    row_lengths : tuple of int, optional
        Number of points in each staggered row of a hexagonal scan.
    step_sizes : tuple of float, optional
        Physical step per axis; required for areas and diameters in specimen
        units.
    acquisition_geometry, calibration_record, measurement_quality : optional
        Acquisition context, cross-checked against the orientations' frames.
    properties : Mapping[str, ArrayLike], optional
        Auxiliary per-point channels such as confidence index or image
        quality. They travel with every derived map.
    provenance : ProvenanceRecord, optional
    """

    coordinates: np.ndarray
    orientations: OrientationSet
    map_frame: ReferenceFrame
    phase_entries: tuple[CrystalMapPhase, ...] = ()
    phase_ids: np.ndarray | None = None
    grid_shape: tuple[int, ...] | None = None
    grid_kind: Literal["square", "hexagonal"] | None = None
    row_lengths: tuple[int, ...] | None = None
    step_sizes: tuple[float, ...] | None = None
    acquisition_geometry: AcquisitionGeometry | None = None
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    properties: Mapping[str, ArrayLike] | None = None
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
        object.__setattr__(
            self,
            "properties",
            _freeze_property_channels(self.properties, point_count=len(self.orientations)),
        )
        coordinate_dims = int(coordinates.shape[1])
        grid_shape = None
        if self.grid_shape is not None:
            grid_shape = tuple(int(size) for size in self.grid_shape)
            if len(grid_shape) != coordinate_dims:
                raise ValueError("CrystalMap.grid_shape must match the coordinate dimensionality.")
            if any(size <= 0 for size in grid_shape):
                raise ValueError("CrystalMap.grid_shape entries must be strictly positive.")
            if int(np.prod(grid_shape)) != len(self.orientations):
                raise ValueError(
                    "CrystalMap.grid_shape must contain exactly one cell per orientation."
                )
            object.__setattr__(self, "grid_shape", grid_shape)
        grid_kind = self.grid_kind
        if grid_kind is None and grid_shape is not None:
            grid_kind = "square"
        if grid_kind not in {None, "square", "hexagonal"}:
            raise ValueError("CrystalMap.grid_kind must be 'square', 'hexagonal', or None.")
        row_lengths = None
        if self.row_lengths is not None:
            row_lengths = tuple(int(length) for length in self.row_lengths)
            if not row_lengths or any(length <= 0 for length in row_lengths):
                raise ValueError("CrystalMap.row_lengths entries must be strictly positive.")
            if sum(row_lengths) != len(self.orientations):
                raise ValueError(
                    "CrystalMap.row_lengths must contain exactly one entry per orientation."
                )
            if any(abs(left - right) > 1 for left, right in pairwise(row_lengths)):
                raise ValueError(
                    "Adjacent CrystalMap hexagonal row lengths may differ by at most one point."
                )
        if grid_kind == "square":
            if grid_shape is None:
                raise ValueError("CrystalMap square topology requires grid_shape.")
            if row_lengths is not None:
                raise ValueError("CrystalMap square topology does not use row_lengths.")
        elif grid_kind == "hexagonal":
            if coordinate_dims != 2:
                raise ValueError("CrystalMap hexagonal topology requires 2-D coordinates.")
            if grid_shape is not None:
                raise ValueError(
                    "CrystalMap hexagonal topology uses row_lengths, not rectangular grid_shape."
                )
            if row_lengths is None:
                raise ValueError("CrystalMap hexagonal topology requires row_lengths.")
            _vectorized_hexagonal_grid_pairs(coordinates, row_lengths, order=1)
        elif row_lengths is not None:
            raise ValueError("CrystalMap.row_lengths requires grid_kind='hexagonal'.")
        object.__setattr__(self, "grid_kind", grid_kind)
        object.__setattr__(self, "row_lengths", row_lengths)
        if self.step_sizes is not None:
            if len(self.step_sizes) != coordinate_dims:
                raise ValueError("CrystalMap.step_sizes must match the coordinate dimensionality.")
            step_sizes = tuple(float(step) for step in self.step_sizes)
            if any(step <= 0.0 for step in step_sizes):
                raise ValueError("CrystalMap.step_sizes entries must be strictly positive.")
            object.__setattr__(self, "step_sizes", step_sizes)

    @property
    def is_multiphase(self) -> bool:
        """Whether the map declares more than one phase.

        Multiphase maps must carry per-point ``phase_ids`` and must leave phase
        and symmetry off the shared ``OrientationSet``, so that no calculation
        can silently apply one phase's symmetry to another's points.
        """

        return len(self.phase_entries) > 1

    @property
    def has_phase_assignments(self) -> bool:
        """Whether any phase information is attached, in either supported form.

        True when explicit ``phase_entries`` exist, or when the underlying
        ``OrientationSet`` carries a single phase.
        """

        return bool(self.phase_entries) or self.orientations.phase is not None

    @property
    def phase_id_array(self) -> np.ndarray | None:
        """Per-point phase ids, or ``None`` when the map has no phase information.

        For a single-phase map without explicit ids this synthesizes an
        all-zero array, so downstream code can treat both cases uniformly.
        """

        if self.phase_ids is not None:
            return self.phase_ids
        if self.orientations.phase is None and self.orientations.symmetry is None:
            return None
        values = np.zeros(len(self.orientations), dtype=np.int64)
        values.setflags(write=False)
        return values

    @property
    def resolved_phase_entries(self) -> tuple[CrystalMapPhase, ...]:
        """The map's phase entries, synthesizing one for a single-phase map.

        A map built from an ``OrientationSet`` that carries a phase has no
        explicit entries; this presents that case in the same form as an
        explicitly multiphase map. Returns an empty tuple when no phase
        information exists at all.
        """

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
        """The single :class:`~pytex.core.lattice.Phase` of the map, when unambiguous.

        Returns ``None`` for a genuinely multiphase map, and also when the one
        phase present is known only by symmetry and has no full ``Phase``
        attached — a common state for vendor imports.
        """

        if self.orientations.phase is not None and not self.phase_entries:
            return self.orientations.phase
        if len(self.resolved_phase_entries) == 1:
            return self.resolved_phase_entries[0].phase
        return None

    @property
    def property_names(self) -> tuple[str, ...]:
        """Names of the auxiliary per-point channels carried with the map.

        Typical channels are the vendor's confidence index, image quality, band
        contrast, or fit. They travel with every derived map, so a
        quality-filtered sub-map keeps its quality channel.
        """

        properties = cast("Mapping[str, np.ndarray]", self.properties)
        return tuple(properties.keys())

    def get_property(self, name: str) -> np.ndarray:
        """The raw per-point values of one auxiliary channel.

        Raises ``KeyError`` naming the available channels when the requested one
        does not exist. Returns the flat per-point array; use
        :meth:`property_map` for the grid-shaped view.
        """

        properties = cast("Mapping[str, np.ndarray]", self.properties)
        if name not in properties:
            available = ", ".join(sorted(properties)) or "<none>"
            raise KeyError(
                f"CrystalMap has no property channel '{name}'. Available channels: {available}."
            )
        return properties[name]

    def property_map(self, name: str) -> np.ndarray:
        """One auxiliary channel reshaped to the map grid.

        Rectangular maps return their native shape. Hexagonal maps use a
        padded array with ``NaN`` where the shorter rows have no measurement.
        Use the flat :meth:`get_property` values for numerical reductions.
        """

        raw = self.get_property(name)
        if self.grid_kind == "hexagonal":
            assert self.row_lengths is not None
            values = np.full(
                (len(self.row_lengths), max(self.row_lengths)),
                np.nan,
                dtype=np.float64,
            )
            offset = 0
            for row, length in enumerate(self.row_lengths):
                values[row, :length] = raw[offset : offset + length]
                offset += length
        else:
            rows, cols = self._require_regular_2d_grid()
            values = raw.reshape((rows, cols))
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def with_properties(
        self,
        properties: Mapping[str, ArrayLike],
        *,
        replace: bool = False,
    ) -> CrystalMap:
        """A copy of the map with additional or replaced auxiliary channels.

        Purpose
        -------
        ``CrystalMap`` is immutable, so derived per-point quantities — a KAM
        map, a Schmid-factor map, a computed mask — are attached by producing a
        new map rather than mutating this one.

        Parameters
        ----------
        properties : Mapping[str, ArrayLike]
            Channels to attach; each must have one value per point.
        replace : bool
            When ``False`` (default) the new channels are merged over the
            existing ones. When ``True`` the existing channels are discarded.

        Returns
        -------
        CrystalMap
            A new map sharing coordinates, orientations, grid, and provenance.
        """

        existing = dict(cast("Mapping[str, np.ndarray]", self.properties))
        merged: dict[str, ArrayLike] = {} if replace else dict(existing)
        merged.update(dict(properties))
        return CrystalMap(
            coordinates=self.coordinates,
            orientations=self.orientations,
            map_frame=self.map_frame,
            phase_entries=self.phase_entries,
            phase_ids=self.phase_ids,
            grid_shape=self.grid_shape,
            grid_kind=self.grid_kind,
            row_lengths=self.row_lengths,
            step_sizes=self.step_sizes,
            acquisition_geometry=self.acquisition_geometry,
            calibration_record=self.calibration_record,
            measurement_quality=self.measurement_quality,
            properties=merged,
            provenance=self.provenance,
        )

    def phase_summary(self) -> dict[str, int]:
        """Point count per phase name.

        Returns an empty dict when the map carries no phase assignments. This
        is the phase-fraction count by *points*, which equals the area fraction
        only on a regular grid with uniform step size.
        """

        phase_entries = self.resolved_phase_entries
        phase_ids = self.phase_id_array
        if not phase_entries or phase_ids is None:
            return {}
        return {
            entry.name: int(np.count_nonzero(phase_ids == entry.phase_id))
            for entry in phase_entries
        }

    def summary(self) -> dict[str, Any]:
        """A compact machine-readable description of the map.

        Reports point count, coordinate dimensionality, grid shape and step
        sizes (``None`` in graph mode), multiphase flag, per-phase point counts,
        and the map and specimen frame names. Intended for logging, manifests,
        and quick inspection rather than for numerical work.
        """

        return {
            "point_count": len(self.orientations),
            "coordinate_dimensions": int(self.coordinates.shape[1]),
            "grid_kind": self.grid_kind,
            "grid_shape": (
                None if self.grid_shape is None else tuple(int(value) for value in self.grid_shape)
            ),
            "row_lengths": self.row_lengths,
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

    @property
    def default_connectivity(self) -> int:
        """Natural first-shell connectivity of the declared topology.

        Hexagonal scans use six neighbours. Rectangular and unstructured maps
        retain the historical four-neighbour default.
        """

        return 6 if self.grid_kind == "hexagonal" else 4

    def describe(self) -> str:
        """Convention-explicit scientific prose describing this EBSD map.

        The description names the scan topology, neighbourhood convention,
        phase state, output-shape behavior, and the rectangular-only limit on
        curvature and pixel-face perimeter calculations.
        """

        if self.grid_kind == "hexagonal":
            assert self.row_lengths is not None
            topology = (
                f"a staggered hexagonal scan with {len(self.row_lengths)} rows "
                f"({', '.join(str(value) for value in self.row_lengths)} points) and "
                "six-neighbour first-shell topology"
            )
            limits = (
                "Local metrics return one value per measured point; padded display arrays mark "
                "missing row positions. Curvature/GND and pixel-face perimeter calculations "
                "remain rectangular-grid-only."
            )
        elif self.grid_shape is not None:
            topology = (
                "a rectangular " + " x ".join(str(value) for value in self.grid_shape) + " grid"
            )
            limits = "Grid-shaped local metrics use the declared row-major raster."
        else:
            topology = "an unstructured coordinate graph"
            limits = "Local metrics return one value per point and require coordinate neighbours."
        phases = self.phase_summary()
        phase_text = (
            ", ".join(f"{name}: {count}" for name, count in phases.items())
            if phases
            else "no explicit phase counts"
        )
        return (
            f"CrystalMap contains {len(self.orientations)} orientations in the "
            f"{self.map_frame.name} frame on {topology}. Phase assignments: {phase_text}. "
            f"{limits}"
        )

    def validate(self) -> tuple[str, ...]:
        """Advisory notes about limitations of this map, as human-readable strings.

        Purpose
        -------
        Report conditions that are legal but will restrict what can be computed
        — phases known only by symmetry with no full ``Phase`` attached, or a
        map in graph mode where grid-shaped outputs are unavailable. An empty
        tuple means no such limitation was found.

        This is an advisory surface, not a validator: structural invariants are
        enforced at construction time and raise there.
        """

        notes: list[str] = []
        if self.is_multiphase:
            unresolved = [
                entry.name for entry in self.resolved_phase_entries if entry.phase is None
            ]
            if unresolved:
                notes.append(
                    "Full Phase objects are not attached for phases: " + ", ".join(unresolved)
                )
        if self.grid_shape is None and self.grid_kind is None:
            notes.append("Map is operating in graph mode rather than regular-grid mode.")
        if self.grid_kind == "hexagonal":
            notes.append(
                "Hexagonal topology supports neighbourhood metrics, KAM, and segmentation; "
                "curvature/GND and pixel-face perimeters require a rectangular grid."
            )
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
        """Boolean per-point mask selecting one phase.

        Accepts a phase id, name, ``Phase``, or ``CrystalMapPhase``. Raises when
        the map carries no phase assignments, rather than returning an
        all-true mask that would quietly conflate phases.
        """

        phase_ids = self.phase_id_array
        if phase_ids is None:
            raise ValueError("CrystalMap does not carry explicit phase assignments.")
        entry = self._resolve_phase_entry(selector)
        mask = np.ascontiguousarray(phase_ids == entry.phase_id)
        mask.setflags(write=False)
        return mask

    def phase_entry_for_index(self, index: int) -> CrystalMapPhase | None:
        """The phase entry of a single measurement point, or ``None``.

        ``None`` means the map carries no phase assignments at all.
        """

        phase_ids = self.phase_id_array
        if phase_ids is None:
            return None
        return self._phase_entry_by_id()[int(phase_ids[int(index)])]

    def select_phase(self, selector: int | str | Phase | CrystalMapPhase) -> CrystalMap:
        """The sub-map containing only the points of one phase.

        Purpose
        -------
        Most texture and grain calculations are single-phase by nature. This
        extracts one phase's points and rebuilds their ``OrientationSet`` with
        that phase's symmetry attached, so the sub-map is fully self-describing.

        Parameters
        ----------
        selector : int, str, Phase, or CrystalMapPhase
            Which phase to keep; see :meth:`CrystalMapPhase.matches`.

        Returns
        -------
        CrystalMap
            The sub-map, carrying masked property channels. Its ``grid_shape``
            is dropped to ``None`` — graph mode — unless every point was
            selected, because a phase subset of a grid is generally not a grid.

        Raises
        ------
        ValueError
            When the phase has no points in this map.
        """

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
        properties = {
            name: values[mask]
            for name, values in cast("Mapping[str, np.ndarray]", self.properties).items()
        }
        return CrystalMap(
            coordinates=self.coordinates[mask],
            orientations=orientations,
            map_frame=self.map_frame,
            grid_shape=self.grid_shape if full_selection else None,
            grid_kind=self.grid_kind if full_selection else None,
            row_lengths=self.row_lengths if full_selection else None,
            step_sizes=self.step_sizes,
            acquisition_geometry=self.acquisition_geometry,
            calibration_record=self.calibration_record,
            measurement_quality=self.measurement_quality,
            properties=properties or None,
            provenance=provenance,
        )

    def select_points(self, mask: ArrayLike) -> CrystalMap:
        """Return the sub-map of points where ``mask`` is True (graph mode).

        Phase assignments and property channels are carried forward, masked to
        the retained points. The result drops ``grid_shape`` unless every point
        is selected, since an arbitrary mask is generally not a full raster.
        """

        mask_array = np.asarray(mask, dtype=bool)
        if mask_array.shape != (len(self.orientations),):
            raise ValueError("select_points mask must have one entry per orientation.")
        if not np.any(mask_array):
            raise ValueError("select_points mask must retain at least one point.")
        full_selection = bool(np.all(mask_array))
        indices = np.flatnonzero(mask_array)
        orientations = self.orientations.subset(indices)
        phase_ids = None if self.phase_ids is None else self.phase_ids[mask_array]
        properties = {
            name: values[mask_array]
            for name, values in cast("Mapping[str, np.ndarray]", self.properties).items()
        }
        return CrystalMap(
            coordinates=self.coordinates[mask_array],
            orientations=orientations,
            map_frame=self.map_frame,
            phase_entries=self.phase_entries,
            phase_ids=phase_ids,
            grid_shape=self.grid_shape if full_selection else None,
            grid_kind=self.grid_kind if full_selection else None,
            row_lengths=self.row_lengths if full_selection else None,
            step_sizes=self.step_sizes,
            acquisition_geometry=self.acquisition_geometry,
            calibration_record=self.calibration_record,
            measurement_quality=self.measurement_quality,
            properties=properties or None,
            provenance=self.provenance,
        )

    def property_threshold_mask(
        self,
        name: str,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
    ) -> np.ndarray:
        """Boolean keep-mask for a property channel within ``[minimum, maximum]``."""

        if minimum is None and maximum is None:
            raise ValueError("property_threshold_mask requires at least one of minimum/maximum.")
        values = self.get_property(name)
        keep = np.ones(values.shape, dtype=bool)
        if minimum is not None:
            keep &= values >= minimum
        if maximum is not None:
            keep &= values <= maximum
        keep = np.ascontiguousarray(keep)
        keep.setflags(write=False)
        return keep

    def filter_by_property(
        self,
        name: str,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
    ) -> CrystalMap:
        """Drop points whose ``name`` channel falls outside ``[minimum, maximum]``."""

        return self.select_points(
            self.property_threshold_mask(name, minimum=minimum, maximum=maximum)
        )

    def remove_wild_spikes(
        self,
        *,
        threshold_deg: float,
        symmetry_aware: bool = True,
        connectivity: int | None = None,
    ) -> CrystalMap:
        """Replace isolated wild-spike points with their neighborhood mean orientation.

        A point is a wild spike when its minimum disorientation to every one of
        its same-phase neighbors exceeds ``threshold_deg``; such points are
        overwritten with the symmetry-aware mean orientation of those neighbors.
        The map geometry, phase assignments, and property channels are preserved.
        """

        if threshold_deg <= 0.0:
            raise ValueError("threshold_deg must be strictly positive.")
        resolved_connectivity = (
            (6 if self.grid_kind == "hexagonal" else 8) if connectivity is None else connectivity
        )
        neighbor_pairs = self.neighbor_graph(connectivity=resolved_connectivity, order=1).pairs
        neighbor_pairs = neighbor_pairs[self._same_phase_pair_mask(neighbor_pairs)]
        point_count = len(self.orientations)
        adjacency: dict[int, list[int]] = {index: [] for index in range(point_count)}
        for left_index, right_index in neighbor_pairs:
            adjacency[int(left_index)].append(int(right_index))
            adjacency[int(right_index)].append(int(left_index))
        # Detect spikes in one batched reduction over all boundary edges rather
        # than a separate misorientation call per point. The misorientation
        # angle is symmetric in the pair, so each undirected edge contributes to
        # the incident minimum of both endpoints.
        quaternions = np.array(self.orientations.quaternions, copy=True)
        min_incident = np.full(point_count, np.inf, dtype=np.float64)
        has_neighbor = np.zeros(point_count, dtype=bool)
        if neighbor_pairs.size:
            edge_angles_deg = np.rad2deg(
                self._pair_misorientation_rad(neighbor_pairs, symmetry_aware=symmetry_aware)
            )
            finite = np.isfinite(edge_angles_deg)
            finite_pairs = neighbor_pairs[finite]
            finite_angles = edge_angles_deg[finite]
            np.minimum.at(min_incident, finite_pairs[:, 0], finite_angles)
            np.minimum.at(min_incident, finite_pairs[:, 1], finite_angles)
            has_neighbor[finite_pairs[:, 0]] = True
            has_neighbor[finite_pairs[:, 1]] = True
        spike_mask = has_neighbor & (min_incident > threshold_deg)
        if not np.any(spike_mask):
            return self
        for index in np.flatnonzero(spike_mask):
            neighbor_indices = np.asarray(adjacency[int(index)], dtype=np.int64)
            neighbor_mean = self.orientations.subset(neighbor_indices).mean_orientation()
            quaternions[int(index)] = neighbor_mean.rotation.quaternion
        orientations = OrientationSet.from_quaternions(
            quaternions,
            crystal_frame=self.orientations.crystal_frame,
            specimen_frame=self.orientations.specimen_frame,
            symmetry=self.orientations.symmetry,
            phase=self.orientations.phase,
            provenance=self.provenance,
        )
        return CrystalMap(
            coordinates=self.coordinates,
            orientations=orientations,
            map_frame=self.map_frame,
            phase_entries=self.phase_entries,
            phase_ids=self.phase_ids,
            grid_shape=self.grid_shape,
            grid_kind=self.grid_kind,
            row_lengths=self.row_lengths,
            step_sizes=self.step_sizes,
            acquisition_geometry=self.acquisition_geometry,
            calibration_record=self.calibration_record,
            measurement_quality=self.measurement_quality,
            properties=cast("Mapping[str, np.ndarray]", self.properties) or None,
            provenance=self.provenance,
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
        """Export the map's acquisition context as a machine-readable manifest.

        Purpose
        -------
        Cross-tool workflows must carry frames, calibration, and measurement
        quality with the data. This produces the schema-validated
        :class:`~pytex.adapters.ExperimentManifest`, synthesizing a minimal EBSD
        acquisition geometry when the map does not already carry one, and adding
        grid shape, step sizes, and phase names as metadata.

        Parameters
        ----------
        source_system : str
            Name recorded as the producing system.
        referenced_files : tuple of str
            Paths recorded alongside the manifest.
        metadata : dict, optional
            Extra key/value pairs, merged over the derived ones.
        """

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
        """Estimate the orientation distribution function from the map.

        Purpose
        -------
        The bridge from a measured EBSD map to quantitative texture: every
        indexed point contributes one orientation, kernel-smoothed into a
        continuous ODF from which volume fractions and texture strength follow.

        Parameters
        ----------
        phase : int, str, Phase, or CrystalMapPhase, optional
            Required for a multiphase map; an ODF mixes phases meaninglessly.
        weights : ArrayLike, optional
            One weight per point — use a confidence-index or image-quality
            channel to down-weight poorly indexed points. Uniform when omitted.
        kernel : KernelSpec, optional
            Smoothing kernel and half-width. The half-width controls the
            bias-variance trade-off and should reflect the measurement's angular
            resolution, not be chosen for appearance.
        specimen_symmetry : SymmetrySpec, optional
            Statistical sample symmetry to impose (orthorhombic for rolled
            sheet, for example). Imposing it is an assumption about the process,
            so it is off by default.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        ODF
        """

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
        """Pole figure of one crystal plane, computed from the map.

        Purpose
        -------
        Where a given crystal plane normal points in the specimen frame across
        all measured points — the standard texture representation, here computed
        directly from EBSD orientations rather than measured by diffraction.

        Parameters
        ----------
        pole : CrystalPlane or ArrayLike
            The plane whose normal is plotted, as a typed plane or ``(hkl)``.
        phase : optional
            Required for a multiphase map.
        weights : ArrayLike, optional
            One weight per point.
        include_symmetry_family : bool
            Plot the whole ``{hkl}`` family (default), as measured pole figures
            inevitably do, rather than the single variant.
        antipodal : bool
            Fold each pole onto one hemisphere (default), the usual convention.
        sample_symmetry : SymmetrySpec, optional
            Statistical specimen symmetry to impose.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        PoleFigure
        """

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
        """Inverse pole figure of a specimen direction, computed from the map.

        Purpose
        -------
        Which crystal directions align with a chosen specimen axis — the
        representation behind IPF colouring and behind fibre-texture statements
        such as "a strong <111> fibre along ND".

        Parameters
        ----------
        sample_direction : str or ArrayLike
            Specimen axis, by name (``"x"``, ``"y"``, ``"z"``, ``"RD"``,
            ``"TD"``, ``"ND"``) or as a vector. Defaults to ``"z"``.
        phase : optional
            Required for a multiphase map.
        weights : ArrayLike, optional
            One weight per point.
        reduce_by_symmetry : bool
            Fold directions into the crystal-symmetry fundamental sector
            (default), which is what makes the standard triangle standard.
        antipodal : bool
            Treat a direction and its reverse as equivalent (default).
        provenance : ProvenanceRecord, optional

        Returns
        -------
        InversePoleFigure
        """

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
        sample_directions: str | ArrayLike | tuple[str | ArrayLike, ...] | list[str | ArrayLike] = (
            "x",
            "y",
            "z",
        ),
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
        """Compute a complete texture characterization of the map in one call.

        Purpose
        -------
        The teaching- and reporting-oriented entry point: it builds the ODF, the
        requested pole figures, and the inverse pole figures from one consistent
        set of choices, so that the pieces of a texture description cannot drift
        apart in conventions, weighting, or symmetry assumptions.

        Parameters
        ----------
        poles : CrystalPlane, ArrayLike, or sequence thereof
            Planes to produce pole figures for. Empty by default.
        sample_directions : str, ArrayLike, or sequence thereof
            Specimen axes to produce inverse pole figures for; ``("x", "y", "z")``
            by default.
        phase : optional
            Required for a multiphase map.
        weights : ArrayLike, optional
            One weight per point, applied consistently to every product.
        kernel : KernelSpec, optional
            ODF smoothing kernel.
        specimen_symmetry : SymmetrySpec, optional
            Statistical specimen symmetry, applied to both ODF and pole figures.
        include_symmetry_family, reduce_by_symmetry, antipodal : bool
            Conventions passed through to the pole and inverse pole figures.
        plot : bool
            Also render Matplotlib figures for each product and attach them to
            the report. Off by default so the computation has no plotting
            dependency.
        provenance : ProvenanceRecord, optional
        pole_figure_plot_kwargs, inverse_pole_figure_plot_kwargs, odf_plot_kwargs : dict, optional
            Forwarded to the corresponding plotting functions when ``plot`` is
            set.

        Returns
        -------
        TextureReport
            Carrying the ODF, the figures, and (when requested) the rendered
            Matplotlib figures.
        """

        phase_view = self._phase_resolved_view(phase=phase, operation="texture_report")
        report_provenance = phase_view.provenance if provenance is None else provenance
        odf = phase_view.to_odf(
            weights=weights,
            kernel=kernel,
            specimen_symmetry=specimen_symmetry,
            provenance=report_provenance,
        )
        if isinstance(poles, list | tuple) and len(poles) == 0:
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
        connectivity: int | None = None,
        order: int = 1,
        max_distance: float | None = None,
    ) -> CoordinateNeighborGraph:
        """The point-adjacency graph underlying every local-misorientation metric.

        Purpose
        -------
        Decide, once and explicitly, which measurement points count as
        neighbours. KAM, grain segmentation, boundary extraction, and smoothing
        all consume this graph, so the neighbourhood definition lives in one
        place rather than being re-derived per metric.

        Parameters
        ----------
        connectivity : int
            ``4`` (edges) or ``8`` (edges and corners) on a rectangular grid;
            ``6`` on a hexagonal grid. The default follows the declared
            topology.
        order : int
            Neighbour shell. ``1`` is nearest neighbours; higher orders reach
            further and are what "KAM of order n" means.
        max_distance : float, optional
            Radius cut-off. Forces the coordinate-radius path even on a regular
            grid, and is required for irregular point clouds.

        Returns
        -------
        CoordinateNeighborGraph
            Unique unordered pairs and their distances. ``mode`` records which
            path produced them: ``"regular_grid"`` for the fast vectorized grid
            construction, ``"hexagonal_grid"`` for staggered logical rows, or
            ``"coordinate_radius"`` for the distance-based fallback.
        """

        resolved_connectivity = self.default_connectivity if connectivity is None else connectivity
        if resolved_connectivity not in {4, 6, 8}:
            raise ValueError("connectivity must be 4, 6, or 8.")
        if self.grid_kind == "hexagonal" and resolved_connectivity != 6:
            raise ValueError("Hexagonal CrystalMap topology requires connectivity=6.")
        if self.grid_kind == "square" and resolved_connectivity == 6:
            raise ValueError("Rectangular CrystalMap topology requires connectivity=4 or 8.")
        if order <= 0:
            raise ValueError("order must be strictly positive.")
        if self.grid_kind == "hexagonal" and max_distance is None:
            assert self.row_lengths is not None
            pairs = _vectorized_hexagonal_grid_pairs(
                self.coordinates,
                self.row_lengths,
                order=order,
            )
            distances = np.linalg.norm(
                self.coordinates[pairs[:, 0]] - self.coordinates[pairs[:, 1]],
                axis=1,
            )
            return CoordinateNeighborGraph(
                pairs=pairs,
                distances=distances,
                connectivity=6,
                order=order,
                mode="hexagonal_grid",
            )
        if self.grid_shape is not None and len(self.grid_shape) == 2:
            rows, cols = self._require_regular_2d_grid()
            if max_distance is None:
                pairs = _vectorized_regular_grid_pairs(
                    rows,
                    cols,
                    connectivity=resolved_connectivity,
                    order=order,
                )
                distances = np.linalg.norm(
                    self.coordinates[pairs[:, 0]] - self.coordinates[pairs[:, 1]],
                    axis=1,
                )
                return CoordinateNeighborGraph(
                    pairs=pairs,
                    distances=distances,
                    connectivity=resolved_connectivity,
                    order=order,
                    mode="regular_grid",
                )
        radius = max_distance
        if radius is None:
            base_spacing = _inferred_base_spacing(self.coordinates, self.step_sizes)
            radius_scale = float(order) * (np.sqrt(2.0) if resolved_connectivity == 8 else 1.0)
            radius = base_spacing * radius_scale + 1e-9
        distances = _pairwise_distances(self.coordinates)
        upper_mask = np.triu(np.ones_like(distances, dtype=bool), k=1)
        pair_mask = upper_mask & (distances <= radius)
        indices = np.column_stack(np.nonzero(pair_mask)).astype(np.int64)
        pair_distances = np.asarray(distances[pair_mask], dtype=np.float64)
        return CoordinateNeighborGraph(
            pairs=indices,
            distances=pair_distances,
            connectivity=resolved_connectivity,
            order=order,
            mode="coordinate_radius",
            max_distance=radius,
        )

    def neighbor_pairs(self, *, connectivity: int | None = None) -> np.ndarray:
        """First-shell neighbour index pairs, as an ``(m, 2)`` array.

        Shorthand for ``neighbor_graph(connectivity=..., order=1).pairs``.
        """

        return self.neighbor_graph(connectivity=connectivity, order=1).pairs

    def kam_neighbor_pairs(self, *, order: int = 1) -> np.ndarray:
        """Neighbour pairs for a KAM calculation at the given shell order.

        Uses four-connectivity on a square grid and six-connectivity on a
        hexagonal grid.
        """

        return self.neighbor_graph(order=order).pairs

    def kernel_average_misorientation_deg(
        self,
        *,
        symmetry_aware: bool = True,
        connectivity: int | None = None,
        order: int = 1,
        threshold_deg: float | None = None,
        statistic: str = "mean",
        segmentation: GrainSegmentation | None = None,
    ) -> np.ndarray:
        """Kernel Average Misorientation, in degrees.

        Purpose
        -------
        Per point, the average misorientation to its neighbours — the standard
        local-deformation measure. KAM responds to the geometrically necessary
        dislocation content stored in short-wavelength orientation gradients,
        which is why it maps deformation structure where a grain-average number
        cannot.

        Parameters
        ----------
        symmetry_aware : bool
            Reduce each pair misorientation by crystal symmetry (default).
        connectivity : int
            ``4`` or ``8`` on a rectangular grid, ``6`` on a hexagonal grid;
            omitted to use the natural topology.
        order : int
            Neighbour shell; higher orders average over a wider kernel.
        threshold_deg : float, optional
            Exclude neighbour pairs above this misorientation. This is the
            standard way to keep grain boundaries out of an intragranular KAM:
            without it, boundary pixels report the boundary misorientation
            rather than the local gradient.
        statistic : str
            ``"mean"`` (default) or ``"max"`` over the neighbourhood.
        segmentation : GrainSegmentation, optional
            Restrict averaging to same-grain pairs. Stricter and more physical
            than a threshold, and required for the GAM definition. Must be a
            segmentation of this same map instance.

        Returns
        -------
        np.ndarray
            Degrees, shaped to the grid for a regular 2-D map and per point
            otherwise; read-only. Points with no admissible neighbour report
            zero. Pairs that cross a phase boundary are never included, since a
            misorientation between different phases is not defined here.
        """

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

    def schmid_factor_map(
        self,
        family: Any,
        stress_direction: str | ArrayLike,
    ) -> np.ndarray:
        """Per-point maximum Schmid factor for a slip-system family.

        ``family`` is a `pytex.properties.SlipSystemFamily`; ``stress_direction``
        is a specimen-frame axis (label such as ``"x"`` or an explicit vector).
        The result is reshaped to the map grid when a regular 2-D grid exists.
        """

        stress = _specimen_direction_vector(stress_direction, self.orientations.specimen_frame)
        values = np.asarray(
            family.max_schmid_factor(self.orientations, stress),
            dtype=np.float64,
        )
        if self.grid_shape is not None and len(self.grid_shape) == 2:
            rows, cols = self._require_regular_2d_grid()
            values = values.reshape((rows, cols))
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def taylor_factor_map(
        self,
        family: Any,
        *,
        tension_axis: str | ArrayLike = "z",
        strain_tensor: ArrayLike | None = None,
    ) -> np.ndarray:
        """Per-point full-constraint Taylor factor for a slip-system family.

        ``family`` is a `pytex.properties.SlipSystemFamily`. Provide either a
        ``tension_axis`` (specimen-frame label or vector, default ``"z"``) for
        uniaxial tension or an explicit deviatoric ``strain_tensor``. The result
        is reshaped to the map grid when a regular 2-D grid exists.
        """

        from pytex.properties.taylor import taylor_factors

        if strain_tensor is not None:
            values = np.asarray(
                taylor_factors(family, self.orientations, strain_tensor=strain_tensor),
                dtype=np.float64,
            )
        else:
            axis = _specimen_direction_vector(tension_axis, self.orientations.specimen_frame)
            values = np.asarray(
                taylor_factors(family, self.orientations, tension_axis=axis),
                dtype=np.float64,
            )
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
        # Representative = member with the least total disorientation to the
        # others. The row sums are accumulated in blocks, so a large grain never
        # allocates an (n, n) matrix, and members tied to within floating-point
        # noise resolve to the lowest index rather than to whatever the
        # summation order happened to produce.
        best = _disorientation_medoid_index(
            self.orientations.subset(member_indices),
            symmetry_aware=symmetry_aware,
        )
        return int(member_indices[best])

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
        connectivity: int | None = None,
    ) -> GrainSegmentation:
        """Group measurement points into grains by neighbour misorientation.

        Purpose
        -------
        Turn a point-wise orientation map into the grain objects that
        microstructural statistics need — sizes, shapes, boundaries, and grain
        averages all follow from this segmentation.

        Method
        ------
        Union-find over the first-shell neighbour graph: two neighbouring points
        join the same grain when their misorientation is at or below
        ``max_misorientation_deg``. This is the standard flood-fill
        segmentation; it merges points connected by a chain of small steps, so a
        grain with a continuous gradient can exceed the threshold end to end.
        Pairs across a phase boundary are never joined.

        Parameters
        ----------
        max_misorientation_deg : float
            The grain-boundary criterion. Conventionally 5-15 degrees; the value
            determines whether subgrains are resolved as separate grains.
        symmetry_aware : bool
            Use disorientation rather than raw rotation angle (default).
        connectivity : int
            ``4`` or ``8`` on a rectangular grid, ``6`` on a hexagonal grid;
            omitted to use the natural topology.

        Returns
        -------
        GrainSegmentation
            Grain labels, grain objects, and the settings used, so downstream
            metrics inherit the same conventions.
        """

        if max_misorientation_deg < 0.0:
            raise ValueError("max_misorientation_deg must be non-negative.")
        graph = self.neighbor_graph(connectivity=connectivity, order=1)
        neighbor_pairs = graph.pairs
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
            connectivity=graph.connectivity,
        )
