from __future__ import annotations

import re
from typing import Any

import numpy as np

from pytex.core.notation import format_plane_family_indices
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
