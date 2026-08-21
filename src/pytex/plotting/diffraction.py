from __future__ import annotations

import re
from typing import Any

import numpy as np

from pytex.core.notation import (
    format_direction_indices,
    format_plane_family_indices,
    format_plane_indices,
)
from pytex.diffraction.kikuchi import GnomonicProjection, KikuchiPattern
from pytex.diffraction.saed import SAEDPattern
from pytex.diffraction.xrd import PowderPattern
from pytex.plotting.frames import add_frame_indicator
from pytex.plotting.styles import resolve_style


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, object


def plot_xrd_pattern(
    pattern: PowderPattern,
    *,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot a simulated powder X-ray diffraction pattern.

    Parameters
    ----------
    pattern : PowderPattern
        Intensity against ``2*theta``, from
        :func:`~pytex.diffraction.generate_xrd_pattern`.
    Remaining parameters control peak labelling, styling, and the target
    axes.

    Returns
    -------
    Any
        The Matplotlib axes.
    """

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    xrd_style = style["xrd"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    axes.plot(
        pattern.two_theta_grid_deg,
        pattern.intensity_grid,
        color=xrd_style["line_color"],
        linewidth=float(common["line"]["width"]),
    )
    axes.fill_between(
        pattern.two_theta_grid_deg,
        pattern.intensity_grid,
        color=xrd_style["fill_color"],
        alpha=0.55,
    )
    axes.set_xlabel(r"2$\theta$ (deg)")
    axes.set_ylabel("normalized intensity")
    if bool(xrd_style.get("show_title", True)):
        # The rotated {hkl} peak labels stand in the strip above the axes, which
        # is exactly where an unpadded title sits: the two printed on top of
        # each other. Reserve that strip when the labels are drawn.
        title_pad = (
            float(xrd_style.get("label_title_pad", 34.0))
            if xrd_style.get("annotate_peaks", True)
            else None
        )
        axes.set_title(
            f"{pattern.phase.name} Powder XRD ({pattern.radiation.name})", pad=title_pad
        )
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    if xrd_style.get("annotate_peaks", True):
        ranked = sorted(
            pattern.reflections,
            key=lambda reflection: reflection.intensity,
            reverse=True,
        )
        low, high = pattern.two_theta_grid_deg[0], pattern.two_theta_grid_deg[-1]
        span = float(high - low)
        # Distinct families can sit within a fraction of a degree of each other,
        # and one label was drawn per reflection with no regard for position, so
        # in a dense pattern several rotated {hkl} strings printed on the same
        # spot as an unreadable stack. Keep the strongest of any crowded group.
        separation = float(xrd_style.get("label_min_separation_fraction", 0.018)) * span
        labelled: list[float] = []
        for reflection in ranked[: int(xrd_style.get("max_labels", 12))]:
            two_theta = float(reflection.two_theta_deg)
            if any(abs(two_theta - placed) < separation for placed in labelled):
                continue
            labelled.append(two_theta)
            # A powder reflection is the whole symmetry-related family (that is
            # what its multiplicity counts), so the family brackets {hkl} are the
            # correct notation, not the single-plane form (hkl).
            label = format_plane_family_indices(
                tuple(int(value) for value in reflection.miller_indices), style="plain"
            )
            axes.axvline(
                reflection.two_theta_deg,
                color=xrd_style["peak_color"],
                alpha=0.25,
                linewidth=0.8,
            )
            axes.text(
                reflection.two_theta_deg,
                1.015,
                label,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=float(common["font"]["size"]) - 1.0,
                color=xrd_style["peak_color"],
                transform=axes.get_xaxis_transform(),
            )
    fig.tight_layout()
    return fig


def plot_saed_pattern(
    pattern: SAEDPattern,
    *,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    show_frame_indicator: bool = False,
    frame_indicator_loc: str = "lower right",
    ax: Any | None = None,
) -> Any:
    """Render a kinematic SAED pattern in detector coordinates.

    Parameters
    ----------
    pattern:
        The `SAEDPattern` to draw. Its spot positions are detector-plane
        coordinates in millimetres.
    theme, style_path, style_overrides:
        Plot styling, resolved through `pytex.plotting.styles.resolve_style`.
    show_frame_indicator:
        Draw the pattern's own detector frame as a small gizmo in a corner, so
        the figure states which way ``u`` and ``v`` point rather than relying on
        the axis labels alone. The detector normal is omitted because it points
        at the viewer. Off by default so existing figures are unchanged.
    frame_indicator_loc:
        Which corner the gizmo occupies; see
        `pytex.plotting.frames.add_frame_indicator`.
    ax:
        An existing axes to draw into. A new figure is created when omitted.

    Returns
    -------
    matplotlib.figure.Figure
        The figure holding the pattern. The caller owns it.
    """

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    saed_style = style["saed"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=saed_style["background"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(saed_style["background"])
    if pattern.spots:
        coordinates = np.vstack([spot.detector_coordinates for spot in pattern.spots])
        intensities = np.array([spot.intensity for spot in pattern.spots], dtype=np.float64)
        if np.max(intensities) > 0.0:
            sizes = float(saed_style["spot_scale"]) * intensities / np.max(intensities)
        else:
            sizes = np.full_like(intensities, float(saed_style["spot_scale"]))
        axes.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=sizes,
            color=saed_style["spot_color"],
            edgecolors="white",
            linewidths=0.35,
            alpha=0.95,
        )
        if saed_style.get("annotate_spots", True):
            span = float(pattern.detector_extent_mm())
            # Labels sat exactly on their spots, so in a dense zone every index
            # string overprinted its neighbours into an unreadable smear. Offset
            # each label clear of its own spot and drop any that would land on a
            # label already placed.
            offset = float(saed_style.get("label_offset_fraction", 0.022)) * span
            separation = float(saed_style.get("label_min_separation_fraction", 0.075)) * span
            # An index string is several times wider than it is tall, so a
            # circular exclusion zone still lets neighbouring labels run into
            # each other horizontally. The zone is widened per character.
            placed: list[tuple[np.ndarray, float]] = []
            for spot in pattern.spots[: int(saed_style.get("max_labels", 16))]:
                if spot.label is None:
                    continue
                position = np.asarray(spot.detector_coordinates, dtype=np.float64)
                # mathtext markup is not drawn, so measure the glyphs the
                # label actually renders as.
                glyphs = max(1, len(re.sub(r"\\[a-zA-Z]+|[${}]", "", spot.label)))
                width = separation * (0.6 + 0.55 * glyphs)
                if any(
                    abs(float(position[0] - other[0])) < 0.5 * (width + other_width)
                    and abs(float(position[1] - other[1])) < separation
                    for other, other_width in placed
                ):
                    continue
                placed.append((position, width))
                axes.text(
                    float(position[0]) + offset,
                    float(position[1]) + offset,
                    spot.label,
                    fontsize=float(common["font"]["size"]) - 1.5,
                    color=saed_style["label_color"],
                    ha="left",
                    va="bottom",
                )
    extent = pattern.detector_extent_mm()
    axes.axhline(0.0, color=saed_style["ring_color"], linewidth=0.8, alpha=0.45)
    axes.axvline(0.0, color=saed_style["ring_color"], linewidth=0.8, alpha=0.45)
    axes.set_xlim(-extent, extent)
    axes.set_ylim(-extent, extent)
    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("detector u (mm)")
    axes.set_ylabel("detector v (mm)")
    if bool(saed_style.get("show_title", True)):
        axes.set_title(
            "SAED Pattern "
            + "("
            + " ".join(str(int(value)) for value in pattern.zone_axis.indices)
            + " zone axis)"
        )
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    if show_frame_indicator:
        # Viewed down the detector normal, so u and v lie in the page exactly as
        # the plotted coordinates do. The normal n is omitted: it points at the
        # viewer and would project to a point on top of the origin.
        add_frame_indicator(
            axes,
            pattern.detector_frame,
            loc=frame_indicator_loc,
            axis_subset=("u", "v"),
            elev_deg=90.0,
            azim_deg=-90.0,
            label_frame=True,
        )
    fig.tight_layout()
    return fig


def _label_along_band(
    axes: Any,
    trace: np.ndarray,
    label: str,
    *,
    bounds: tuple[float, float, float, float],
    color: str,
    fontsize: float,
) -> None:
    """Write a band's name along the band, a fifth of the way in from one end.

    Only the samples inside ``bounds`` are candidates. A centre trace is sampled
    across the whole hemisphere of directions, and most of it is nowhere near
    the detector, so a fraction taken along the whole trace puts the name off
    the picture — where the reader never sees it and nothing reports that it is
    missing.

    Within the visible stretch, the middle is where the bands of a zone cross,
    so a name placed there lands on the hub the pattern is read by, and the very
    end is under the frame. A fifth of the way in is clear of both.

    The angle is folded into the readable half-turn, because a line has no
    direction and text upside down is not a label. ``rotation_mode="anchor"``
    with ``transform_rotates_text`` keeps the text on the band when the axes are
    rescaled, rather than at an angle that was right for one figure size.
    """

    left, right, bottom, top = bounds
    inside = (
        (trace[:, 0] >= left)
        & (trace[:, 0] <= right)
        & (trace[:, 1] >= bottom)
        & (trace[:, 1] <= top)
    )
    # The *longest contiguous* visible stretch, not every visible sample. A
    # trace can leave the picture and come back, and a baseline measured across
    # that gap is not a direction: it is a chord between two different parts of
    # the band, and the name would be written at an angle the band never takes.
    longest: np.ndarray | None = None
    start: int | None = None
    for position, visible in enumerate([*inside, False]):
        if visible and start is None:
            start = position
        elif not visible and start is not None:
            run = trace[start:position]
            if longest is None or run.shape[0] > longest.shape[0]:
                longest = run
            start = None
    if longest is None or longest.shape[0] < 2:
        return
    trace = longest
    index = max(1, round(trace.shape[0] * 0.2))
    span = max(1, trace.shape[0] // 20)
    start = trace[max(0, index - span)]
    point = trace[index]
    delta = point - start
    if not np.any(np.abs(delta) > 0.0):
        return
    angle = float(np.degrees(np.arctan2(delta[1], delta[0])))
    if angle > 90.0:
        angle -= 180.0
    elif angle <= -90.0:
        angle += 180.0
    axes.text(
        float(point[0]),
        float(point[1]),
        label,
        rotation=angle,
        rotation_mode="anchor",
        transform_rotates_text=True,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        color=color,
        zorder=6,
        # Clipped to the axes, and not merely for tidiness: an unclipped text
        # artist counts towards the figure's tight bounding box, so a band named
        # near the edge would push the layout out until matplotlib reported that
        # it could not fit the decorations.
        clip_on=True,
    )


def plot_kikuchi_pattern(
    pattern: KikuchiPattern,
    *,
    coordinates: str = "gnomonic",
    show_edges: bool = True,
    show_zone_axes: bool = True,
    max_bands: int | None = None,
    label_zone_axes: bool = True,
    label_bands: bool = False,
    samples: int = 361,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Draw a simulated Kikuchi pattern as band traces and zone axes.

    Purpose
    -------
    Make the geometry legible: each band is drawn as its centre line with the
    two Kossel-cone edges that bound it, and the zone axes are marked where
    bands intersect. Useful for teaching how a pattern is built, and for
    overlaying a simulation on a measured pattern to check an indexing solution.

    Parameters
    ----------
    pattern : KikuchiPattern
        The simulated pattern.
    coordinates : str
        ``"gnomonic"`` (default) plots in the projection where band centres are
        exactly straight — the right frame for judging geometry. ``"detector"``
        plots in detector pixels, which is what a measured pattern looks like
        and where centre lines curve on a tilted detector.
    show_edges : bool
        Draw the Kossel-cone edges as well as the centre lines. The visible gap
        between an edge pair *is* the band.
    show_zone_axes : bool
        Mark the zone axes.
    max_bands : int, optional
        Draw only the strongest ``max_bands`` bands, for legibility.
    label_zone_axes : bool
        Annotate zone axes with their ``[uvw]`` indices.
    label_bands : bool
        Write each band's ``(hkl)`` **along** the band rather than beside it. A
        band is identified by which line it is, so a horizontal caption belongs
        — visually — to whichever line happens to be nearest the text, and on a
        zone-axis pattern several bands cross within a few pixels. Off by
        default, because a pattern carrying tens of bands becomes a page of
        text; lower ``max_bands`` before turning it on.
    samples : int
        Points sampled along each trace.
    theme, style_path, style_overrides :
        Plot styling, resolved through `pytex.plotting.styles.resolve_style`.
    ax : matplotlib Axes, optional
        An existing axes to draw into. A new figure is created when omitted.

    Returns
    -------
    matplotlib.figure.Figure
        The figure holding the pattern. The caller owns it.
    """

    if coordinates not in {"gnomonic", "detector"}:
        raise ValueError("coordinates must be either 'gnomonic' or 'detector'.")
    if samples < 2:
        raise ValueError("samples must be at least two.")
    if max_bands is not None and max_bands <= 0:
        raise ValueError("max_bands must be strictly positive when provided.")

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    saed_style = style["saed"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=saed_style["background"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(saed_style["background"])

    projection = GnomonicProjection(geometry=pattern.geometry)
    in_detector = coordinates == "detector"

    def _place(points: np.ndarray) -> np.ndarray:
        if points.shape[0] == 0:
            return points
        return np.asarray(projection.to_detector_px(points)) if in_detector else points

    bands = pattern.bands if max_bands is None else pattern.bands[:max_bands]
    # Names are written after the limits are set, because where a name belongs
    # depends on which stretch of the band is on the picture.
    named: list[tuple[np.ndarray, str]] = []
    for band in bands:
        centre = _place(band.center_trace(projection, samples=samples))
        if centre.shape[0] > 1:
            axes.plot(
                centre[:, 0],
                centre[:, 1],
                color=saed_style["spot_color"],
                linewidth=1.1,
                alpha=0.9,
            )
            if label_bands:
                named.append(
                    (
                        centre,
                        format_plane_indices(
                            tuple(int(value) for value in band.plane.indices), style="mathtext"
                        ),
                    )
                )
        if not show_edges:
            continue
        for edge in band.edge_traces(projection, samples=samples):
            placed = _place(edge)
            if placed.shape[0] > 1:
                axes.plot(
                    placed[:, 0],
                    placed[:, 1],
                    color=saed_style["ring_color"],
                    linewidth=0.7,
                    alpha=0.55,
                )

    if show_zone_axes:
        visible = [axis for axis in pattern.zone_axes if axis.on_detector]
        if visible:
            points = _place(np.vstack([axis.coordinates for axis in visible]))
            axes.scatter(
                points[:, 0],
                points[:, 1],
                s=float(saed_style["spot_scale"]) * 0.25,
                color=saed_style["label_color"],
                zorder=5,
            )
            if label_zone_axes:
                for axis, point in zip(visible, points, strict=True):
                    indices = tuple(int(value) for value in axis.indices)
                    axes.annotate(
                        format_direction_indices(indices, style="mathtext"),
                        (float(point[0]), float(point[1])),
                        textcoords="offset points",
                        xytext=(4, 4),
                        fontsize=float(common["font"]["size"]) - 1.5,
                        color=saed_style["label_color"],
                    )

    if in_detector:
        height, width = pattern.geometry.detector_shape
        axes.set_xlim(0.0, float(width - 1))
        # Detector row 0 is conventionally the top of the image.
        axes.set_ylim(float(height - 1), 0.0)
        axes.set_xlabel("detector u (px)")
        axes.set_ylabel("detector v (px)")
    else:
        corners = np.asarray(projection.detector_corner_coordinates())
        axes.set_xlim(float(corners[:, 0].min()), float(corners[:, 0].max()))
        axes.set_ylim(float(corners[:, 1].min()), float(corners[:, 1].max()))
        axes.set_xlabel("gnomonic x (detector distances)")
        axes.set_ylabel("gnomonic y (detector distances)")
    axes.set_aspect("equal", adjustable="box")
    if named:
        left, right = sorted(axes.get_xlim())
        bottom, top = sorted(axes.get_ylim())
        for trace, label in named:
            _label_along_band(
                axes,
                trace,
                label,
                bounds=(left, right, bottom, top),
                color=saed_style["label_color"],
                fontsize=float(common["font"]["size"]) - 1.5,
            )
    if bool(saed_style.get("show_title", True)):
        axes.set_title(f"Kikuchi pattern ({pattern.phase.name}, {coordinates} coordinates)")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig
