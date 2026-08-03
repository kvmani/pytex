"""Publication-grade rendering for composite OR SAED patterns.

Purpose: a typed, layered, highly configurable matplotlib renderer for
:class:`pytex.diffraction.composite.CompositeSAEDPattern`. Every visual
aspect — marker shape, fill, size mapping, color, opacity, edge, z-order,
variant subset, axes units, legend, transmitted beam — is controlled by an
explicit frozen configuration object with defaults tuned for publication
figures (open indigo parent circles over filled colored variant markers on a
white background, equal-aspect detector axes).

Rendering is structural: each sub-pattern becomes one scatter collection
tagged with a machine-readable ``gid`` (``pytex-composite:<label>``), which
is also how the regression tests verify the figure semantically (per repo
policy, no image-byte baselines).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, Literal

import numpy as np

from pytex.core.hexagonal import plane_hkl_to_hkil
from pytex.diffraction.composite import (
    CompositeSAEDPattern,
    VariantZonePattern,
    is_hexagonal_phase,
)
from pytex.diffraction.kinematic import SpotTable
from pytex.plotting.frames import add_frame_indicator

SizeModeName = Literal["intensity_area", "intensity_radius", "constant"]
AxesUnitsName = Literal["mm", "inv_angstrom"]
IndexFormatName = Literal["plain", "overline"]

GID_PREFIX = "pytex-composite"

#: Colorblind-aware categorical palette for variant sub-patterns (12 hues:
#: Paul Tol's muted scheme extended with two high-contrast extras).
VARIANT_COLOR_PALETTE: tuple[str, ...] = (
    "#cc6677",
    "#332288",
    "#ddcc77",
    "#117733",
    "#88ccee",
    "#882255",
    "#44aa99",
    "#999933",
    "#aa4499",
    "#dd7788",
    "#6699cc",
    "#661100",
)

#: Marker cycle for variant sub-patterns; combined with the color palette the
#: default styling distinguishes up to 96 variants before repeating.
VARIANT_MARKER_CYCLE: tuple[str, ...] = ("s", "^", "D", "v", "P", "X", "h", "*")


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt


@dataclass(frozen=True, slots=True)
class SpotStyle:
    """Appearance of one sub-pattern's diffraction spots.

    ``size_mode`` maps normalized intensity to the matplotlib scatter size
    (points^2): ``"intensity_area"`` makes marker **area** proportional to
    intensity (the perceptually honest default), ``"intensity_radius"`` makes
    marker radius proportional to intensity (dramatizes strong spots), and
    ``"constant"`` ignores intensity. ``min_size_pt2`` keeps weak reflections
    visible. ``filled=False`` renders hollow markers (edge in ``color``),
    the classic convention for parent-phase overlays.
    """

    marker: str = "o"
    color: str = "#3f51b5"
    filled: bool = True
    size_scale: float = 90.0
    size_mode: SizeModeName = "intensity_area"
    min_size_pt2: float = 6.0
    alpha: float = 0.9
    edge_color: str = "#111111"
    edge_width: float = 0.6
    zorder: float = 3.0

    def __post_init__(self) -> None:
        if not self.marker:
            raise ValueError("SpotStyle.marker must be non-empty.")
        if not np.isfinite(self.size_scale) or self.size_scale <= 0.0:
            raise ValueError("SpotStyle.size_scale must be finite and strictly positive.")
        if self.size_mode not in {"intensity_area", "intensity_radius", "constant"}:
            raise ValueError(
                "SpotStyle.size_mode must be 'intensity_area', 'intensity_radius' or "
                "'constant'."
            )
        if not np.isfinite(self.min_size_pt2) or self.min_size_pt2 <= 0.0:
            raise ValueError("SpotStyle.min_size_pt2 must be finite and strictly positive.")
        if not 0.0 < self.alpha <= 1.0:
            raise ValueError("SpotStyle.alpha must lie in the interval (0, 1].")
        if not np.isfinite(self.edge_width) or self.edge_width < 0.0:
            raise ValueError("SpotStyle.edge_width must be finite and non-negative.")

    def marker_sizes_pt2(self, intensity: np.ndarray) -> np.ndarray:
        """Scatter sizes (points^2) for max-normalized intensities."""

        values = np.asarray(intensity, dtype=np.float64)
        if self.size_mode == "constant":
            sizes = np.full(values.shape, self.size_scale)
        elif self.size_mode == "intensity_area":
            sizes = self.size_scale * values
        else:
            sizes = self.size_scale * values**2
        return np.asarray(np.maximum(sizes, self.min_size_pt2), dtype=np.float64)


def format_hkl(
    hkl: Any, *, index_format: IndexFormatName = "overline", bravais: bool = False
) -> str:
    """Format Miller indices as a TEM-convention reflection label.

    Purpose: consistent spot labels across composite figures. ``"plain"``
    yields ``(1 1 -1)``; ``"overline"`` yields crystallographic mathtext with
    overlined negative indices — compact for single-digit indices
    (``$(11\\bar{1})$``) and thin-space separated when any index has more
    than one digit (``$(12\\;\\bar{1}\\;1)$``).

    With ``bravais`` the three-index ``(hkl)`` is expanded to the four-index
    Miller-Bravais ``(h k i l)`` form with ``i = -(h + k)``, the convention for
    hexagonal phases such as the alpha-hcp product of the Burgers relationship.
    """

    values = [int(value) for value in np.asarray(hkl, dtype=np.int64).reshape(3)]
    if bravais:
        values = [int(value) for value in plane_hkl_to_hkil(values)]
    if index_format == "plain":
        return "(" + " ".join(str(value) for value in values) + ")"
    if index_format != "overline":
        raise ValueError("index_format must be 'plain' or 'overline'.")
    compact = all(abs(value) <= 9 for value in values)
    tokens = [
        rf"\bar{{{abs(value)}}}" if value < 0 else str(abs(value)) for value in values
    ]
    separator = "" if compact else r"\;"
    return "$(" + separator.join(tokens) + ")$"


@dataclass(frozen=True, slots=True)
class SpotAnnotationConfig:
    """Configuration for spot labeling with crowding avoidance.

    ``merge_coincident`` collapses reflections from different sub-patterns
    that land within ``coincidence_tolerance_mm`` on the detector into one
    multi-line label (each line tagged ``p`` for the parent or ``Vn`` for a
    variant) — the standard way composite OR patterns are annotated.
    ``max_labels`` caps the number of label *clusters*, keeping the densest
    composites readable; clusters are prioritized by intensity, then radius.
    With ``avoid_overlap`` labels are placed greedily on a two-ring compass
    of candidate anchors and dropped (never overlapped) when no free
    position exists; ``leader_lines`` draws a thin connector when a label
    lands on the outer ring.
    """

    enabled: bool = True
    max_labels: int = 24
    min_intensity: float = 0.05
    index_format: IndexFormatName = "overline"
    merge_coincident: bool = True
    coincidence_tolerance_mm: float = 2.0
    offset_pt: float = 7.0
    font_size: float = 7.0
    text_color: str = "#111111"
    label_color_follows_spot: bool = True
    bbox_alpha: float = 0.65
    avoid_overlap: bool = True
    leader_lines: bool = False

    def __post_init__(self) -> None:
        if self.max_labels < 1:
            raise ValueError("max_labels must be at least 1.")
        if not 0.0 <= self.min_intensity <= 1.0:
            raise ValueError("min_intensity must lie in the interval [0, 1].")
        if self.index_format not in {"plain", "overline"}:
            raise ValueError("index_format must be 'plain' or 'overline'.")
        if (
            not np.isfinite(self.coincidence_tolerance_mm)
            or self.coincidence_tolerance_mm < 0.0
        ):
            raise ValueError("coincidence_tolerance_mm must be finite and non-negative.")
        if not np.isfinite(self.offset_pt) or self.offset_pt <= 0.0:
            raise ValueError("offset_pt must be finite and strictly positive.")
        if not np.isfinite(self.font_size) or self.font_size <= 0.0:
            raise ValueError("font_size must be finite and strictly positive.")
        if not 0.0 <= self.bbox_alpha <= 1.0:
            raise ValueError("bbox_alpha must lie in the interval [0, 1].")


@dataclass(frozen=True, slots=True)
class AnnotationResult:
    """Outcome of the annotation pass: what was labeled and what was dropped.

    ``texts`` and ``positions_data`` describe the placed labels (positions
    are label anchor points in data units); ``skipped_count`` counts label
    clusters dropped because no collision-free anchor existed within the
    candidate rings. ``cluster_count`` is the number of label clusters that
    passed the intensity floor and budget (placed + skipped).
    """

    texts: tuple[str, ...]
    positions_data: np.ndarray
    cluster_count: int
    placed_count: int
    skipped_count: int
    merged_cluster_count: int

    def __post_init__(self) -> None:
        positions = as_float_array_2d(self.positions_data, rows=len(self.texts))
        positions.setflags(write=False)
        object.__setattr__(self, "positions_data", positions)
        object.__setattr__(self, "texts", tuple(self.texts))
        if self.placed_count != len(self.texts):
            raise ValueError("placed_count must equal the number of placed texts.")
        if self.placed_count + self.skipped_count != self.cluster_count:
            raise ValueError("cluster_count must equal placed_count + skipped_count.")
        if self.merged_cluster_count < 0 or self.merged_cluster_count > self.cluster_count:
            raise ValueError("merged_cluster_count must lie in [0, cluster_count].")

    def describe(self) -> str:
        """Prose summary of annotation coverage and crowding decisions."""

        return (
            f"Spot annotation: {self.placed_count} label(s) placed out of "
            f"{self.cluster_count} candidate cluster(s) ({self.merged_cluster_count} "
            f"merged coincident-reflection label(s); {self.skipped_count} dropped for "
            "lack of collision-free space). Labels use TEM reflection notation; lines "
            "in a merged label are tagged 'p' (parent) or 'Vn' (variant n) and refer "
            "to reflections landing within the coincidence tolerance on the shared "
            "detector."
        )


def as_float_array_2d(values: Any, *, rows: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0:
        return np.zeros((0, 2), dtype=np.float64)
    if array.shape != (rows, 2):
        raise ValueError(f"positions_data must have shape ({rows}, 2).")
    return np.ascontiguousarray(array)


_DEFAULT_PARENT_STYLE = SpotStyle(
    marker="o",
    color="#3f51b5",
    filled=False,
    size_scale=140.0,
    edge_width=1.4,
    alpha=1.0,
    zorder=5.0,
)
_DEFAULT_CHILD_STYLE = SpotStyle(
    marker="s",
    color="#f57c00",
    filled=True,
    size_scale=80.0,
    edge_width=0.5,
    alpha=0.85,
    zorder=3.0,
)


@dataclass(frozen=True, slots=True)
class CompositeSAEDPlotConfig:
    """Full rendering configuration for composite SAED figures.

    ``variant_styles`` overrides the palette/marker cycling per 1-based
    variant index; unlisted variants fall back to ``child_base_style`` with
    the color palette and marker cycle applied in rendering order.
    ``variant_indices`` restricts rendering to a subset without re-simulating.
    ``axes_units`` switches between calibrated detector millimetres and
    reciprocal angstrom. The transmitted beam is drawn as a distinct central
    marker (never part of any spot table).

    ``show_frame_indicator`` adds a small gizmo showing where the **parent
    crystal axes** point on this detector, projected through the pattern's
    parent-anchored detector basis. It answers "which way is the crystal
    oriented in this pattern?" without a prose caption; off by default so
    existing figures are unchanged.
    """

    parent_style: SpotStyle = field(default_factory=lambda: _DEFAULT_PARENT_STYLE)
    child_base_style: SpotStyle = field(default_factory=lambda: _DEFAULT_CHILD_STYLE)
    variant_styles: Mapping[int, SpotStyle] | None = None
    variant_color_palette: tuple[str, ...] = VARIANT_COLOR_PALETTE
    variant_marker_cycle: tuple[str, ...] = VARIANT_MARKER_CYCLE
    variant_indices: tuple[int, ...] | None = None
    show_parent: bool = True
    show_transmitted_beam: bool = True
    transmitted_beam_size_pt2: float = 240.0
    transmitted_beam_color: str = "#111111"
    axes_units: AxesUnitsName = "mm"
    show_legend: bool = True
    legend_max_entries: int = 13
    legend_outside: bool = True
    show_title: bool = True
    title: str | None = None
    background: str = "#ffffff"
    figsize: tuple[float, float] = (7.0, 7.0)
    dpi: int = 200
    limit_padding_fraction: float = 0.08
    annotation: SpotAnnotationConfig = field(default_factory=SpotAnnotationConfig)
    show_frame_indicator: bool = False
    frame_indicator_loc: str = "lower left"

    def __post_init__(self) -> None:
        if not self.variant_color_palette:
            raise ValueError("variant_color_palette must contain at least one color.")
        if not self.variant_marker_cycle:
            raise ValueError("variant_marker_cycle must contain at least one marker.")
        if (
            not np.isfinite(self.transmitted_beam_size_pt2)
            or self.transmitted_beam_size_pt2 <= 0.0
        ):
            raise ValueError("transmitted_beam_size_pt2 must be finite and strictly positive.")
        if self.axes_units not in {"mm", "inv_angstrom"}:
            raise ValueError("axes_units must be 'mm' or 'inv_angstrom'.")
        if self.legend_max_entries < 1:
            raise ValueError("legend_max_entries must be at least 1.")
        if not 0.0 <= self.limit_padding_fraction < 1.0:
            raise ValueError("limit_padding_fraction must lie in the interval [0, 1).")
        if self.dpi < 1:
            raise ValueError("dpi must be at least 1.")

    def style_for_variant(self, variant_index: int, render_position: int) -> SpotStyle:
        """Resolved style for a variant: explicit override or cycled default."""

        if self.variant_styles is not None and variant_index in self.variant_styles:
            return self.variant_styles[variant_index]
        color = self.variant_color_palette[
            render_position % len(self.variant_color_palette)
        ]
        marker = self.variant_marker_cycle[
            render_position % len(self.variant_marker_cycle)
        ]
        return replace(self.child_base_style, color=color, marker=marker)


def _spot_coordinates(table: SpotTable, axes_units: AxesUnitsName) -> np.ndarray:
    if axes_units == "mm":
        return table.detector_mm
    return table.g_detector_inv_angstrom


def _scatter_kwargs(style: SpotStyle) -> dict[str, Any]:
    if style.filled:
        return {
            "facecolors": style.color,
            "edgecolors": style.edge_color,
            "linewidths": style.edge_width,
        }
    return {
        "facecolors": "none",
        "edgecolors": style.color,
        "linewidths": max(style.edge_width, 0.8),
    }


def _scatter_sub_pattern(
    ax: Any,
    table: SpotTable,
    style: SpotStyle,
    *,
    label: str,
    gid: str,
    axes_units: AxesUnitsName,
) -> None:
    coordinates = _spot_coordinates(table, axes_units)
    collection = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        s=style.marker_sizes_pt2(table.intensity),
        marker=style.marker,
        alpha=style.alpha,
        zorder=style.zorder,
        label=label,
        **_scatter_kwargs(style),
    )
    collection.set_gid(gid)


@dataclass(frozen=True, slots=True)
class _LabelSpot:
    """One reflection queued for labeling (internal)."""

    coordinates: np.ndarray
    intensity: float
    line: str
    sub_pattern_order: int
    color: str


def _collect_label_spots(
    rendered: CompositeSAEDPattern,
    plot_config: CompositeSAEDPlotConfig,
) -> list[_LabelSpot]:
    annotation = plot_config.annotation
    spots: list[_LabelSpot] = []
    order = 0
    parent_bravais = is_hexagonal_phase(rendered.relationship.parent_phase)
    child_bravais = is_hexagonal_phase(rendered.relationship.child_phase)
    parent_table = rendered.parent_spots
    if plot_config.show_parent and parent_table is not None:
        coordinates = _spot_coordinates(parent_table, plot_config.axes_units)
        for row in range(len(parent_table)):
            label = format_hkl(
                parent_table.hkl[row],
                index_format=annotation.index_format,
                bravais=parent_bravais,
            )
            spots.append(
                _LabelSpot(
                    coordinates=coordinates[row],
                    intensity=float(parent_table.intensity[row]),
                    line=f"{label} p",
                    sub_pattern_order=order,
                    color=plot_config.parent_style.color,
                )
            )
    order += 1
    for position, variant_pattern in enumerate(rendered.variant_patterns):
        style = plot_config.style_for_variant(variant_pattern.variant_index, position)
        table = variant_pattern.spots
        coordinates = _spot_coordinates(table, plot_config.axes_units)
        for row in range(len(table)):
            label = format_hkl(
                table.hkl[row],
                index_format=annotation.index_format,
                bravais=child_bravais,
            )
            spots.append(
                _LabelSpot(
                    coordinates=coordinates[row],
                    intensity=float(table.intensity[row]),
                    line=f"{label} V{variant_pattern.variant_index}",
                    sub_pattern_order=order + position,
                    color=style.color,
                )
            )
    return spots


def _cluster_label_spots(
    spots: list[_LabelSpot], tolerance_units: float
) -> list[list[int]]:
    """Union-find clustering of spot indices within the coincidence tolerance."""

    count = len(spots)
    parent_index = list(range(count))

    def find(index: int) -> int:
        while parent_index[index] != index:
            parent_index[index] = parent_index[parent_index[index]]
            index = parent_index[index]
        return index

    if count and tolerance_units > 0.0:
        from scipy.spatial import cKDTree

        coordinates = np.vstack([spot.coordinates for spot in spots])
        tree = cKDTree(coordinates)
        for left, right in tree.query_pairs(tolerance_units):
            root_left, root_right = find(int(left)), find(int(right))
            if root_left != root_right:
                parent_index[root_right] = root_left
    clusters: dict[int, list[int]] = {}
    for index in range(count):
        clusters.setdefault(find(index), []).append(index)
    return list(clusters.values())


def _annotate_composite(
    fig: Any,
    axes: Any,
    rendered: CompositeSAEDPattern,
    plot_config: CompositeSAEDPlotConfig,
) -> AnnotationResult:
    annotation = plot_config.annotation
    spots = _collect_label_spots(rendered, plot_config)
    if not spots:
        return AnnotationResult(
            texts=(),
            positions_data=np.zeros((0, 2)),
            cluster_count=0,
            placed_count=0,
            skipped_count=0,
            merged_cluster_count=0,
        )

    tolerance_units = annotation.coincidence_tolerance_mm
    if plot_config.axes_units == "inv_angstrom":
        tolerance_units /= rendered.config.camera_constant_mm_angstrom
    if annotation.merge_coincident:
        clusters = _cluster_label_spots(spots, tolerance_units)
    else:
        clusters = [[index] for index in range(len(spots))]

    ranked: list[tuple[float, float, float, float, list[int]]] = []
    for members in clusters:
        max_intensity = max(spots[index].intensity for index in members)
        if max_intensity < annotation.min_intensity:
            continue
        anchor = np.mean(
            np.vstack([spots[index].coordinates for index in members]), axis=0
        )
        radius = float(np.linalg.norm(anchor))
        ranked.append((-max_intensity, radius, float(anchor[0]), float(anchor[1]), members))
    ranked.sort(key=lambda item: item[:4])
    ranked = ranked[: annotation.max_labels]

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    display_transform = axes.transData

    spot_display = display_transform.transform(
        np.vstack([spot.coordinates for spot in spots])
    )
    placed_boxes: list[tuple[float, float, float, float]] = []
    texts: list[str] = []
    positions: list[np.ndarray] = []
    skipped = 0
    merged_count = 0
    directions = np.array(
        [(1, 1), (-1, 1), (1, -1), (-1, -1), (0, 1), (0, -1), (1, 0), (-1, 0)],
        dtype=np.float64,
    )
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    rings = (1.0, 2.4)

    def boxes_overlap(
        left: tuple[float, float, float, float],
        right: tuple[float, float, float, float],
    ) -> bool:
        return not (
            left[2] <= right[0]
            or right[2] <= left[0]
            or left[3] <= right[1]
            or right[3] <= left[1]
        )

    for _, _, _, _, members in ranked:
        ordered = sorted(
            members,
            key=lambda index: (spots[index].sub_pattern_order, -spots[index].intensity),
        )
        lines = [spots[index].line for index in ordered]
        text_value = "\n".join(lines)
        anchor = np.mean(np.vstack([spots[index].coordinates for index in ordered]), axis=0)
        single = len({spots[index].sub_pattern_order for index in ordered}) == 1
        color = (
            spots[ordered[0]].color
            if annotation.label_color_follows_spot and single
            else annotation.text_color
        )

        artist = axes.annotate(
            text_value,
            xy=(float(anchor[0]), float(anchor[1])),
            xytext=(annotation.offset_pt, annotation.offset_pt),
            textcoords="offset points",
            fontsize=annotation.font_size,
            color=color,
            zorder=7.0,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "#ffffff",
                "alpha": annotation.bbox_alpha,
                "edgecolor": "none",
            },
        )
        placed = False
        used_ring = 0.0
        for ring in rings:
            for direction in directions:
                offset = direction * annotation.offset_pt * ring
                artist.set_position((float(offset[0]), float(offset[1])))
                artist.set_horizontalalignment("left" if direction[0] > 0 else
                                               "right" if direction[0] < 0 else "center")
                artist.set_verticalalignment("bottom" if direction[1] > 0 else
                                             "top" if direction[1] < 0 else "center")
                extent = artist.get_window_extent(renderer=renderer)
                box = (
                    float(extent.x0) - 1.0,
                    float(extent.y0) - 1.0,
                    float(extent.x1) + 1.0,
                    float(extent.y1) + 1.0,
                )
                if any(boxes_overlap(box, other) for other in placed_boxes):
                    continue
                inside = (
                    (spot_display[:, 0] > box[0])
                    & (spot_display[:, 0] < box[2])
                    & (spot_display[:, 1] > box[1])
                    & (spot_display[:, 1] < box[3])
                )
                for index in ordered:
                    inside[index] = False
                if bool(np.any(inside)):
                    continue
                placed = True
                used_ring = ring
                placed_boxes.append(box)
                break
            if placed:
                break
        if not placed:
            if annotation.avoid_overlap:
                artist.remove()
                skipped += 1
                continue
            used_ring = rings[0]
            extent = artist.get_window_extent(renderer=renderer)
            placed_boxes.append(
                (float(extent.x0), float(extent.y0), float(extent.x1), float(extent.y1))
            )
        artist.set_gid(f"{GID_PREFIX}:annotation:{len(texts)}")
        if len(members) > 1:
            merged_count += 1
        if annotation.leader_lines and used_ring > 1.0:
            offset_now = artist.get_position()
            leader = axes.annotate(
                "",
                xy=(float(anchor[0]), float(anchor[1])),
                xytext=offset_now,
                textcoords="offset points",
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#888888",
                    "linewidth": 0.6,
                    "shrinkA": 0.0,
                    "shrinkB": 2.0,
                },
                zorder=6.5,
            )
            leader.set_gid(f"{GID_PREFIX}:leader:{len(texts)}")
        texts.append(text_value)
        positions.append(anchor)

    return AnnotationResult(
        texts=tuple(texts),
        positions_data=np.vstack(positions) if positions else np.zeros((0, 2)),
        cluster_count=len(ranked),
        placed_count=len(texts),
        skipped_count=skipped,
        merged_cluster_count=merged_count,
    )


def render_composite_saed(
    pattern: CompositeSAEDPattern,
    *,
    config: CompositeSAEDPlotConfig | None = None,
    ax: Any | None = None,
    return_annotations: bool = False,
) -> Any:
    """Render a composite OR SAED pattern to a matplotlib figure.

    Purpose: the presentation layer over
    :func:`pytex.diffraction.composite.simulate_composite_saed`. Children are
    drawn first (palette/marker cycling or explicit per-variant styles), the
    parent on top as hollow markers, then the transmitted beam at the origin;
    axes are equal-aspect with units chosen by the configuration.

    Inputs: ``config`` — a :class:`CompositeSAEDPlotConfig` (defaults are
    publication-ready); ``ax`` — optional existing matplotlib axes (its
    figure is reused), otherwise a new figure is created with the configured
    size, dpi and background.

    Output: the matplotlib ``Figure``, or ``(figure, AnnotationResult)``
    when ``return_annotations`` is true. Every sub-pattern scatter carries
    the gid ``pytex-composite:<label>``, every label
    ``pytex-composite:annotation:<i>`` (leaders
    ``pytex-composite:leader:<i>``), for structural inspection and testing.
    Spot labels follow ``config.annotation`` (see
    :class:`SpotAnnotationConfig`): coincident reflections merge into
    phase-tagged multi-line labels placed collision-free.
    """

    plot_config = CompositeSAEDPlotConfig() if config is None else config
    plt = _require_matplotlib()

    rendered = pattern
    if plot_config.variant_indices is not None:
        rendered = pattern.select_variants(list(plot_config.variant_indices))

    if ax is None:
        fig, axes = plt.subplots(
            figsize=plot_config.figsize,
            dpi=plot_config.dpi,
            facecolor=plot_config.background,
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(plot_config.background)

    child_phase_name = rendered.relationship.child_phase.name
    variant_patterns: tuple[VariantZonePattern, ...] = rendered.variant_patterns
    for position, variant_pattern in enumerate(variant_patterns):
        style = plot_config.style_for_variant(variant_pattern.variant_index, position)
        label = f"{child_phase_name} {variant_pattern.label()}"
        _scatter_sub_pattern(
            axes,
            variant_pattern.spots,
            style,
            label=label,
            gid=f"{GID_PREFIX}:variant:{variant_pattern.variant_index}",
            axes_units=plot_config.axes_units,
        )

    parent_table = rendered.parent_spots
    if plot_config.show_parent and parent_table is not None:
        parent_zone = rendered.parent_zone_axis_label()
        parent_label = f"{rendered.relationship.parent_phase.name} {parent_zone}"
        _scatter_sub_pattern(
            axes,
            parent_table,
            plot_config.parent_style,
            label=parent_label,
            gid=f"{GID_PREFIX}:parent",
            axes_units=plot_config.axes_units,
        )

    if plot_config.show_transmitted_beam:
        beam = axes.scatter(
            [0.0],
            [0.0],
            s=plot_config.transmitted_beam_size_pt2,
            marker="o",
            facecolors=plot_config.transmitted_beam_color,
            edgecolors="none",
            zorder=6.0,
            label="transmitted beam",
        )
        beam.set_gid(f"{GID_PREFIX}:transmitted-beam")

    coordinates = rendered.all_detector_coordinates()
    if plot_config.axes_units == "inv_angstrom":
        coordinates = (
            coordinates / rendered.config.camera_constant_mm_angstrom
            if coordinates.size
            else coordinates
        )
    if coordinates.size:
        extent = float(np.max(np.abs(coordinates)))
    else:
        extent = 1.0
    extent = extent * (1.0 + plot_config.limit_padding_fraction) or 1.0
    axes.set_xlim(-extent, extent)
    axes.set_ylim(-extent, extent)
    axes.set_aspect("equal", adjustable="box")

    if plot_config.axes_units == "mm":
        axes.set_xlabel("detector u (mm)")
        axes.set_ylabel("detector v (mm)")
    else:
        # Reciprocal-space axes: g is the reciprocal-lattice vector, bolded per
        # IUCr vector convention, and the unit is inverse angstrom.
        axes.set_xlabel(r"$\mathbf{g}\cdot\hat{u}$ ($\mathrm{\AA}^{-1}$)")
        axes.set_ylabel(r"$\mathbf{g}\cdot\hat{v}$ ($\mathrm{\AA}^{-1}$)")

    if plot_config.show_title:
        if plot_config.title is not None:
            axes.set_title(plot_config.title)
        else:
            parent_zone = rendered.parent_zone_axis_label()
            axes.set_title(
                f"{rendered.relationship.name}: composite SAED, parent {parent_zone} zone"
            )

    if plot_config.show_legend:
        handles, labels = axes.get_legend_handles_labels()
        if len(handles) > plot_config.legend_max_entries:
            handles = handles[: plot_config.legend_max_entries]
            labels = labels[: plot_config.legend_max_entries]
        if handles:
            legend_kwargs: dict[str, Any] = {"fontsize": 8, "framealpha": 0.9}
            if plot_config.legend_outside:
                legend_kwargs.update({"loc": "upper left", "bbox_to_anchor": (1.02, 1.0)})
            axes.legend(handles, labels, **legend_kwargs)

    if plot_config.show_frame_indicator:
        # zone_basis_parent has the orthonormal detector basis (u, v, zone) as
        # columns in parent-crystal Cartesian components, so its transpose maps
        # a parent-crystal vector into detector components. Its columns are
        # therefore the parent crystal axes as seen on this detector.
        add_frame_indicator(
            axes,
            pattern.relationship.parent_phase.crystal_frame,
            loc=plot_config.frame_indicator_loc,
            basis=np.asarray(pattern.zone_basis_parent, dtype=np.float64).T,
            elev_deg=90.0,
            azim_deg=-90.0,
            label_frame=True,
        )

    # Layout must be final before label placement: the collision checks
    # measure display-space extents, which tight_layout would invalidate.
    fig.tight_layout()
    annotation_result: AnnotationResult | None = None
    if plot_config.annotation.enabled:
        annotation_result = _annotate_composite(fig, axes, rendered, plot_config)

    if return_annotations:
        if annotation_result is None:
            annotation_result = AnnotationResult(
                texts=(),
                positions_data=np.zeros((0, 2)),
                cluster_count=0,
                placed_count=0,
                skipped_count=0,
                merged_cluster_count=0,
            )
        return fig, annotation_result
    return fig


__all__ = [
    "GID_PREFIX",
    "VARIANT_COLOR_PALETTE",
    "VARIANT_MARKER_CYCLE",
    "AnnotationResult",
    "CompositeSAEDPlotConfig",
    "SpotAnnotationConfig",
    "SpotStyle",
    "format_hkl",
    "render_composite_saed",
]
