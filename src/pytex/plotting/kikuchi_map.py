r"""Drawing the stereographic Kikuchi map.

The map is the operator's road atlas, so the drawing has one job: make the
network readable at a glance. Three encodings carry the information, and each
answers a question an operator actually asks.

**Line weight is intensity.** A strong band is a band you will see on the screen;
a weak one you will not. Weight is the kinematic :math:`|F_g|^2` of the
reflection, so the heavy lines on the plot are the heavy lines in the microscope.

**The gap between an edge pair is the band.** Centre lines are drawn as thin
curves and the two Kossel-cone edges as lighter ones on either side, exactly as
:func:`pytex.plotting.plot_kikuchi_pattern` does for a detector pattern. The gap
is the true angular width, :math:`2\theta_B`, which for electrons is of order one
degree — so at map scale it is a hairline. ``width_scale`` exaggerates it for
teaching, and any figure that uses it must say so in its caption.

**Marker area is zone-axis order.** The number of bands crossing an axis is how
conspicuous it is: a four-band crossing is unmistakable on the screen, a two-band
one is a guess. Area scales with that count, so the large markers are the axes
worth aiming for.

A route drawn over the map is a heavy arc per leg, with its waypoint axes
highlighted, because that is what the operator will be tracking.

See Also
--------
pytex.diffraction.kikuchi_map : The map itself, and the routing.
pytex.plotting.tilt_stereogram : The complementary picture — the same
    stereographic plane, but showing the *stage's* reachable region rather than
    the crystal's band network.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pytex.core.hexagonal import direction_uvw_to_uvtw, is_hexagonal_phase
from pytex.core.lattice import Phase
from pytex.core.notation import format_direction_indices
from pytex.diffraction.kikuchi_map import (
    KikuchiMapBand,
    KikuchiRoute,
    StereographicKikuchiMap,
    projected_trace_runs,
)
from pytex.diffraction.stereonets import (
    generate_stereonet_grid,
    projection_boundary_radius,
)
from pytex.plotting.styles import resolve_style
from pytex.texture.projections import project_directions

__all__ = ["plot_kikuchi_map"]

#: Points sampled along each band trace.
#:
#: A band edge is a small circle of angular radius near 90 degrees, so its
#: projection is a long curve; at fewer samples than this the arcs visibly
#: polygonize at publication size.
DEFAULT_TRACE_SAMPLES = 721


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt


def _project_runs(directions: np.ndarray, *, method: str) -> list[np.ndarray]:
    """Projected polylines of a direction sequence, split at equator crossings.

    The splitting itself lives in :func:`pytex.diffraction.projected_trace_runs`,
    because the workbench draws the same map in the browser and must break the
    same curves in the same places.
    """

    return list(projected_trace_runs(directions, method=method))


def _great_circle_arc(start: np.ndarray, end: np.ndarray, samples: int) -> np.ndarray:
    """The short great-circle arc between two unit directions.

    The route legs are arcs, not straight lines in the projection plane, because
    the stage turns the crystal about an axis and the beam direction traces a
    great circle. Drawing them straight would misrepresent the path an operator
    follows.
    """

    first = np.asarray(start, dtype=np.float64)
    second = np.asarray(end, dtype=np.float64)
    if float(np.dot(first, second)) < 0.0:
        second = -second
    angle = float(np.arccos(np.clip(float(np.dot(first, second)), -1.0, 1.0)))
    if angle < 1e-12:
        return first[None, :]
    fractions = np.linspace(0.0, 1.0, max(samples, 2))
    sine = np.sin(angle)
    weights_first = np.sin((1.0 - fractions) * angle) / sine
    weights_second = np.sin(fractions * angle) / sine
    arc = weights_first[:, None] * first[None, :] + weights_second[:, None] * second[None, :]
    return np.asarray(arc / np.linalg.norm(arc, axis=1)[:, None], dtype=np.float64)


def _axis_label(phase: Phase, indices: tuple[int, ...]) -> str:
    values = [int(value) for value in indices]
    if is_hexagonal_phase(phase):
        values = [int(value) for value in direction_uvw_to_uvtw(values)]
    return format_direction_indices(values, style="mathtext")


def plot_kikuchi_map(
    kikuchi_map: StereographicKikuchiMap,
    *,
    route: KikuchiRoute | None = None,
    method: str = "stereographic",
    max_bands: int | None = None,
    show_edges: bool = True,
    width_scale: float = 1.0,
    min_label_order: int = 4,
    max_labels: int | None = 14,
    include_net: bool = True,
    include_minor_net: bool = False,
    samples: int = DEFAULT_TRACE_SAMPLES,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    r"""Draw a stereographic Kikuchi map, optionally with a planned route on it.

    Purpose
    -------
    Turn the computed band network into the picture a TEM operator navigates by:
    which bands are strong, which zone axes are conspicuous, and — with ``route``
    — which band to follow from here to there.

    When to use
    -----------
    After :func:`pytex.diffraction.compute_kikuchi_map`, for planning a tilt
    series or as a reference atlas beside the microscope. For the pattern on a
    physical detector at one known orientation, use
    :func:`pytex.plotting.plot_kikuchi_pattern` instead.

    Reading the figure
    ------------------
    Line weight is kinematic intensity; the gap between an edge pair is the band's
    true angular width; marker area is the number of bands crossing a zone axis.
    The bounding circle is the equator of the crystal sphere — directions at
    90 degrees from the map centre.

    Parameters
    ----------
    kikuchi_map : StereographicKikuchiMap
    route : KikuchiRoute, optional
        Drawn as one heavy arc per leg, with the endpoints and waypoints marked.
        An unreachable route is drawn dashed, so it cannot be mistaken for a plan.
    method : str
        ``"stereographic"`` (default, angle-preserving — the right choice when the
        numbers read off the map are angles the stage must turn through) or
        ``"equal_area"``.
    max_bands : int, optional
        Draw only the strongest this many bands. A full map of a cubic phase to
        index 4 has of order a hundred, which is more than a printed figure can
        carry.
    show_edges : bool
        Draw the Kossel-cone edges as well as the centre lines.
    width_scale : float
        Multiplies each band's Bragg angle for *drawing only*. The true width is
        of order one degree, so at map scale a band is a hairline; a value above 1
        makes the width legible at the cost of no longer being to scale, and a
        caption must say which was used. Does not affect any returned geometry.
    min_label_order : int
        Label zone axes crossed by at least this many bands.
    max_labels : int, optional
        Keep at most this many labels, the most conspicuous first — highest band
        count, then closest to the map centre. A cubic map to index 4 has of order
        a hundred axes above any useful ``min_label_order``, and labelling them all
        produces an unreadable figure; ``None`` labels every qualifying axis.
    include_net : bool
        Draw the Wulff-net graticule behind the bands.
    include_minor_net : bool
        Add the fine graticule as well. Off by default: the style's minor step is
        two degrees, which suits a pole figure carrying a handful of poles and
        reads as a grey smudge behind a band network of a hundred curves.
    samples : int
        Points along each trace.
    title : str, optional
        Overrides the generated title.
    theme, style_path, style_overrides :
        Styling, resolved through :func:`pytex.plotting.styles.resolve_style`.
    ax : matplotlib Axes, optional
        Draw into an existing axes; a new figure is created when omitted.

    Returns
    -------
    matplotlib.figure.Figure
        The figure. The caller owns it and is responsible for closing it.

    Raises
    ------
    ValueError
        For an unknown projection method, fewer than two samples, a non-positive
        ``max_bands``, or a non-positive ``width_scale``.
    """

    if method not in {"stereographic", "equal_area"}:
        raise ValueError("method must be either 'stereographic' or 'equal_area'.")
    if samples < 2:
        raise ValueError("samples must be at least two.")
    if max_bands is not None and max_bands <= 0:
        raise ValueError("max_bands must be strictly positive when provided.")
    if not np.isfinite(width_scale) or width_scale <= 0.0:
        raise ValueError("width_scale must be finite and strictly positive.")
    if min_label_order < 2:
        raise ValueError("min_label_order must be at least 2: a crossing needs two bands.")
    if max_labels is not None and max_labels <= 0:
        raise ValueError("max_labels must be strictly positive when provided.")

    plt = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    spherical = style["spherical"]

    if ax is None:
        figure, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
        )
    else:
        axes = ax
        figure = axes.figure

    radius = projection_boundary_radius(method)
    if include_net:
        net = generate_stereonet_grid(
            method=method,
            major_step_deg=float(spherical["major_step_deg"]),
            minor_step_deg=(
                float(spherical["minor_step_deg"])
                if include_minor_net and bool(spherical["show_minor_grid"])
                else None
            ),
        )
        for line in net.minor_lines:
            axes.plot(
                line[:, 0],
                line[:, 1],
                color=spherical["net_minor_color"],
                linewidth=float(spherical["net_minor_linewidth"]),
                alpha=float(spherical["net_alpha"]),
                zorder=1,
            )
        for line in net.major_lines:
            axes.plot(
                line[:, 0],
                line[:, 1],
                color=spherical["net_major_color"],
                linewidth=float(spherical["net_major_linewidth"]),
                alpha=float(spherical["net_alpha"]),
                zorder=1,
            )
    axes.add_patch(
        plt.Circle(
            (0.0, 0.0),
            radius,
            fill=False,
            edgecolor=spherical["boundary_color"],
            linewidth=float(spherical["boundary_linewidth"]),
            zorder=2,
        )
    )

    bands = kikuchi_map.bands if max_bands is None else kikuchi_map.bands[:max_bands]
    trace_colour = spherical["plane_colors"][0]
    for band in bands:
        # Intensity spans orders of magnitude, so the cube root compresses it into
        # a line-weight range the eye can rank without the weakest bands vanishing.
        weight = float(spherical["plane_trace_linewidth"]) * (
            0.45 + 0.95 * float(band.relative_intensity) ** (1.0 / 3.0)
        )
        for polyline in _project_runs(_band_centre_directions(band, samples), method=method):
            axes.plot(
                polyline[:, 0],
                polyline[:, 1],
                color=trace_colour,
                linewidth=weight,
                alpha=float(spherical["plane_trace_alpha"]),
                zorder=3,
                solid_capstyle="round",
            )
        if not show_edges:
            continue
        for offset in (-1.0, 1.0):
            directions = _band_edge_directions(band, offset, width_scale, samples)
            for polyline in _project_runs(directions, method=method):
                axes.plot(
                    polyline[:, 0],
                    polyline[:, 1],
                    color=trace_colour,
                    linewidth=weight * 0.55,
                    alpha=float(spherical["plane_trace_alpha"]) * 0.55,
                    zorder=3,
                )

    axes_points = np.asarray(
        [axis.projected(method=method) for axis in kikuchi_map.zone_axes], dtype=np.float64
    )
    orders = np.asarray([axis.order for axis in kikuchi_map.zone_axes], dtype=np.float64)
    if axes_points.size:
        axes.scatter(
            axes_points[:, 0],
            axes_points[:, 1],
            s=float(spherical["direction_size"]) * (0.10 + 0.06 * orders),
            marker=spherical["direction_marker"],
            color=spherical["direction_colors"][0],
            edgecolor=spherical["direction_edgecolor"],
            linewidth=float(spherical["direction_linewidth"]),
            zorder=5,
        )
        labelled = [
            (axis, point)
            for axis, point in zip(kikuchi_map.zone_axes, axes_points, strict=True)
            if axis.order >= min_label_order
        ]
        # zone_axes is already ordered by decreasing band count then increasing
        # polar angle, so truncation keeps the most conspicuous axes.
        if max_labels is not None:
            labelled = labelled[:max_labels]
        for axis, point in labelled:
            axes.annotate(
                _axis_label(kikuchi_map.phase, axis.indices),
                (float(point[0]), float(point[1])),
                textcoords="offset points",
                # label_offset in the style is a data-space fraction, which suits a
                # pole figure of fixed radius; on a map whose limits the caller may
                # change, a fixed offset in points keeps labels clear of markers of
                # any size.
                xytext=(5, 4),
                fontsize=float(spherical["label_fontsize"]),
                color=spherical["label_color"],
                zorder=6,
            )

    if route is not None:
        _draw_route(axes, kikuchi_map, route, method=method, samples=samples, style=style)

    limit = radius * 1.08
    axes.set_xlim(-limit, limit)
    axes.set_ylim(-limit, limit)
    axes.set_aspect("equal", adjustable="box")
    axes.set_axis_off()
    if title is None:
        centre = _axis_label(kikuchi_map.phase, kikuchi_map.centre_indices)
        scale_note = "" if width_scale == 1.0 else f", widths x{width_scale:g}"
        title = (
            f"Kikuchi map of {kikuchi_map.phase.name} at "
            f"{kikuchi_map.beam_energy_kev:.0f} kV, centred on {centre}{scale_note}"
        )
    axes.set_title(title)
    figure.tight_layout()
    return figure


def _band_centre_directions(band: KikuchiMapBand, samples: int) -> np.ndarray:
    """Unit directions along a band's centre line, before projection."""

    return np.asarray(band.centre_directions(samples=samples), dtype=np.float64)


def _band_edge_directions(
    band: KikuchiMapBand, offset: float, width_scale: float, samples: int
) -> np.ndarray:
    """Unit directions along one band edge, before projection.

    ``offset`` selects the edge: ``-1`` the narrow side, ``+1`` the far one.
    ``width_scale`` exaggerates the Bragg angle for drawing only.
    """

    narrow, far = band.edge_directions(samples=samples, width_scale=width_scale)
    return np.asarray(narrow if offset < 0.0 else far, dtype=np.float64)


def _draw_route(
    axes: Any,
    kikuchi_map: StereographicKikuchiMap,
    route: KikuchiRoute,
    *,
    method: str,
    samples: int,
    style: dict[str, Any],
) -> None:
    """Overlay a planned route: one arc per leg, with the landmarks marked."""

    spherical = style["spherical"]
    # The route is the one thing on the figure that is not a property of the
    # crystal, so it takes the accent colour rather than any of the band or
    # direction palettes.
    colour = style["common"]["colors"]["accent"]
    linestyle = "-" if route.reachable else "--"
    for index, leg in enumerate(route.legs):
        # The leg's own oriented directions, not the map's canonical representative
        # of each axis: a zone axis is a line, and the two senses project to
        # opposite sides of the disc. Using the canonical sense would draw a path
        # that jumps across the projection between legs.
        arc = _great_circle_arc(
            np.asarray(leg.start_direction), np.asarray(leg.end_direction), samples
        )
        for polyline in _project_runs(arc, method=method):
            # A pale halo under the route: the band network beneath it is dense by
            # nature, and a single coloured line of any weight disappears into it.
            axes.plot(
                polyline[:, 0],
                polyline[:, 1],
                color=style["common"]["colors"]["background"],
                linewidth=float(spherical["plane_trace_linewidth"]) * 5.0,
                alpha=0.85,
                zorder=7,
                solid_capstyle="round",
            )
            axes.plot(
                polyline[:, 0],
                polyline[:, 1],
                color=colour,
                linewidth=float(spherical["plane_trace_linewidth"]) * 2.4,
                linestyle=linestyle,
                alpha=0.95,
                zorder=8,
                solid_capstyle="round",
            )
        midpoint = np.asarray(
            project_directions(arc[len(arc) // 2][None, :], method=method, antipodal=True)
        )[0]
        axes.annotate(
            f"{index + 1}",
            (float(midpoint[0]), float(midpoint[1])),
            textcoords="offset points",
            xytext=(5, -11),
            fontsize=float(spherical["label_fontsize"]),
            color=colour,
            fontweight="bold",
            zorder=9,
        )
        for indices in leg.waypoint_indices:
            waypoint = kikuchi_map.zone_axis_for_direction(indices)
            if waypoint is None:
                continue
            # Same reason: take the sense the route passes through, which is the
            # sense the reported indices carry.
            sense = np.asarray(waypoint.direction_map, dtype=np.float64)
            if not np.array_equal(
                np.sign(np.asarray(waypoint.indices)), np.sign(np.asarray(indices))
            ):
                sense = -sense
            point = np.asarray(
                project_directions(sense[None, :], method=method, antipodal=True)
            )[0]
            axes.scatter(
                [float(point[0])],
                [float(point[1])],
                s=float(spherical["direction_size"]) * 0.55,
                marker="x",
                color=colour,
                linewidth=float(spherical["direction_linewidth"]) * 2.0,
                zorder=8,
            )
    for direction, indices, marker in (
        (route.start_direction, route.start_indices, "o"),
        (route.target_direction, route.target_indices, "*"),
    ):
        point = np.asarray(
            project_directions(
                np.asarray(direction, dtype=np.float64)[None, :],
                method=method,
                antipodal=True,
            )
        )[0]
        # Label the endpoints in the sense the route actually uses. The map's own
        # label for the same axis may sit on the opposite side of the disc, since
        # the two senses of a zone axis are the same axis and the map labels only
        # one of them.
        axes.annotate(
            _axis_label(kikuchi_map.phase, indices),
            (float(point[0]), float(point[1])),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=float(spherical["label_fontsize"]),
            color=colour,
            fontweight="bold",
            zorder=10,
        )
        axes.scatter(
            [float(point[0])],
            [float(point[1])],
            s=float(spherical["direction_size"]) * (2.2 if marker == "*" else 1.4),
            marker=marker,
            facecolor="none" if marker == "o" else colour,
            edgecolor=colour,
            linewidth=float(spherical["direction_linewidth"]) * 2.0,
            zorder=9,
        )
