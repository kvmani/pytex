"""The figures of FIB lamella planning: seven views of one plan.

Purpose
-------
Decision D9 of the FIB lamella foundation document: a numeric answer without its
figure is an incomplete deliverable. Each function here draws one of the figures
of section 9 into a caller-supplied :class:`matplotlib.figure.Figure`, so the
workbench (which renders through :func:`pytex.app.figures.render_figure`) and a
script (:func:`fib_figure`) draw exactly the same picture.

The figures
-----------
1. The frame-convention figure is a canonical SVG asset, written by
   ``scripts/generate_fib_lamella_figures.py`` into ``docs/figures/``.
2. :func:`draw_plan_view` --- the SEM or IQ image with the grain outline, the
   placement rectangle, the azimuth, a scale bar and the lift-out end.
3. :func:`draw_section_schematic` --- the trench, slab, depth and thickness seen
   end-on, with the lamella normal in the surface and the target axis ``eps*``
   out of it.
4. :func:`draw_stereogram` --- the signature figure: every orbit member on an
   upper-hemisphere stereogram of the sample frame, whose primitive circle *is*
   the great circle of achievable beam directions, with ``eps*`` marked.
5. :func:`draw_phi_feasibility` --- a polar plot of the mounting rotation ``phi``:
   the largest residual the holder removes at each ``phi``, the guaranteed disc,
   and the reachable arcs at ``eps*``.
6. :func:`draw_preparability_map` --- ``eps*`` at every point, sequential colour
   map with an explicit scale, the ranked sites on top.
7. :func:`draw_predicted_saed` --- the kinematic pattern at the achieved axis,
   from :func:`pytex.diffraction.saed.generate_saed_pattern`.

Conventions
-----------
Plan views and maps are drawn with scan ``y`` increasing **downwards**, the way
an SEM or EBSD image is displayed; the sample axes are drawn as a gizmo so the
orientation of ``X_s`` and ``Y_s`` on screen is never implied. Every figure takes
``theme="light"`` or ``"dark"``; colours come from the canonical tokens of
``docs/standards/visualization_style_guide.md`` and symbols from
:mod:`pytex.core.symbols`.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from pytex.core.notation import format_direction_family_indices, format_direction_indices
from pytex.core.symbols import symbol_latex

__all__ = [
    "FIB_FIGURE_KINDS",
    "FIB_THEMES",
    "draw_phi_feasibility",
    "draw_plan_view",
    "draw_predicted_saed",
    "draw_preparability_map",
    "draw_section_schematic",
    "draw_stereogram",
    "fib_figure",
    "reach_curve_deg",
]

#: Colours by role, per theme. Light follows the style-guide tokens; dark keeps
#: the hues and lifts them for contrast on the dark paper.
FIB_THEMES: dict[str, dict[str, str]] = {
    "light": {
        "paper": "#fbfdff",
        "ink": "#07122f",
        "muted": "#40506f",
        "grid": "#cbd5e1",
        "primary": "#2563eb",
        "accent": "#0f9f9f",
        "caution": "#f59e0b",
        "risk": "#e11d48",
        "ok": "#16a34a",
        "violet": "#7c3aed",
        "surface": "#e2e8f0",
        "slab": "#93c5fd",
    },
    "dark": {
        "paper": "#0b1220",
        "ink": "#e5e7eb",
        "muted": "#94a3b8",
        "grid": "#334155",
        "primary": "#60a5fa",
        "accent": "#2dd4bf",
        "caution": "#fbbf24",
        "risk": "#fb7185",
        "ok": "#4ade80",
        "violet": "#a78bfa",
        "surface": "#1e293b",
        "slab": "#1d4ed8",
    },
}


def _palette(theme: str) -> dict[str, str]:
    try:
        return FIB_THEMES[theme]
    except KeyError:
        raise ValueError(f"theme must be one of {sorted(FIB_THEMES)}; got {theme!r}.") from None


def _math(name: str) -> str:
    return f"${symbol_latex(name)}$"


def _member_math(plan: Any) -> str:
    """The chosen member as publication notation, negative indices overbarred."""

    return format_direction_indices([int(value) for value in plan.geometry.member_indices])


def _family_math(plan: Any) -> str:
    return format_direction_family_indices([int(value) for value in plan.orbit.target_indices])


def _style_axes(axes: Any, colors: dict[str, str]) -> None:
    axes.set_facecolor(colors["paper"])
    for spine in axes.spines.values():
        spine.set_color(colors["muted"])
    axes.tick_params(colors=colors["muted"], labelcolor=colors["ink"])
    axes.xaxis.label.set_color(colors["ink"])
    axes.yaxis.label.set_color(colors["ink"])
    axes.title.set_color(colors["ink"])


def _prepare(figure: Any, colors: dict[str, str]) -> None:
    figure.set_facecolor(colors["paper"])


def _nice_length(span: float) -> float:
    """A round scale-bar length near a fifth of ``span``."""

    target = max(span / 5.0, 1e-9)
    exponent = math.floor(math.log10(target))
    for mantissa in (1.0, 2.0, 5.0, 10.0):
        if mantissa * 10**exponent >= target:
            return float(mantissa * 10**exponent)
    return float(10 ** (exponent + 1))


def _sample_gizmo(axes: Any, plan: Any, colors: dict[str, str]) -> None:
    """Draw ``X_s`` and ``Y_s`` as they point on the (y-down) scan display.

    Drawn in axes coordinates on a paper-coloured card, so it stays legible
    over any image in either theme.
    """

    from matplotlib.patches import FancyBboxPatch

    axes.add_patch(
        FancyBboxPatch(
            (0.02, 0.02),
            0.2,
            0.2,
            boxstyle="round,pad=0.01",
            transform=axes.transAxes,
            fc=colors["paper"],
            ec=colors["muted"],
            alpha=0.85,
            zorder=5,
        )
    )
    origin = np.array([0.12, 0.12])
    surface = plan.surface
    for label, vector in (("X_s", (1.0, 0.0)), ("Y_s", (0.0, 1.0))):
        scan = surface.sample_to_scan_xy(np.array(vector))
        # Scan y runs down the screen, axes-fraction y runs up.
        screen = np.array([scan[0], -scan[1]])
        tip = origin + 0.07 * screen
        axes.annotate(
            "",
            xy=tuple(tip),
            xytext=tuple(origin),
            xycoords="axes fraction",
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "-|>", "color": colors["ink"], "lw": 1.0},
            zorder=6,
        )
        label_at = origin + 0.093 * screen
        axes.text(
            label_at[0],
            label_at[1],
            f"${label}$",
            transform=axes.transAxes,
            color=colors["ink"],
            ha="center",
            va="center",
            fontsize=8,
            zorder=6,
        )


# --------------------------------------------------------------------------- #
# 2. Plan view
# --------------------------------------------------------------------------- #


def draw_plan_view(
    figure: Any,
    plan: Any,
    *,
    background: np.ndarray | None = None,
    extent: tuple[float, float, float, float] | None = None,
    grain_mask: np.ndarray | None = None,
    background_label: str = "image quality",
    theme: str = "light",
) -> None:
    """The site on the SEM or IQ image: outline, rectangle, azimuth, scale, lift-out end.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
    plan : LamellaPlan
        Must carry a placement.
    background : (rows, cols) array, optional
        Greyscale image on the scan grid (an IQ channel, or a registered SEM
        image resampled to it).
    extent : (x_min, x_max, y_min, y_max), optional
        Pixel-edge extent of ``background`` in scan micrometres.
    grain_mask : (rows, cols) bool array, optional
        The planned grain, outlined.
    background_label : str
    theme : {"light", "dark"}
    """

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1)
    _style_axes(axes, colors)
    placement = plan.placement
    if placement is None:
        raise ValueError("draw_plan_view needs a plan with a footprint placement.")
    if background is not None and extent is not None:
        x0, x1, y0, y1 = extent
        axes.imshow(
            np.asarray(background, dtype=np.float64),
            cmap="gray",
            extent=(x0, x1, y1, y0),
            origin="upper",
            interpolation="nearest",
        )
        if grain_mask is not None:
            rows, cols = grain_mask.shape
            xs = np.linspace(x0, x1, cols + 1)[:-1] + (x1 - x0) / cols / 2
            ys = np.linspace(y0, y1, rows + 1)[:-1] + (y1 - y0) / rows / 2
            axes.contour(
                xs,
                ys,
                grain_mask.astype(float),
                levels=[0.5],
                colors=[colors["caution"]],
                linewidths=1.2,
            )
    corners = np.vstack([placement.corners_scan_um, placement.corners_scan_um[:1]])
    axes.fill(corners[:, 0], corners[:, 1], color=colors["primary"], alpha=0.35, lw=0)
    axes.plot(corners[:, 0], corners[:, 1], color=colors["primary"], lw=1.6)
    cx, cy = placement.center_scan_um
    theta = math.radians(plan.geometry.theta_sample_deg)
    surface = plan.surface
    normal = surface.sample_to_scan_xy(np.array([math.cos(theta), math.sin(theta)]))
    long_axis = surface.sample_to_scan_xy(np.array([-math.sin(theta), math.cos(theta)]))
    # A window on the site: the grain's extent, but never less than 2.5 L across,
    # and never beyond the map.
    half = 1.25 * plan.spec.length_um
    if grain_mask is not None and extent is not None and grain_mask.any():
        rows, cols = grain_mask.shape
        x0, x1, y0, y1 = extent
        rr, cc = np.nonzero(grain_mask)
        gx = x0 + np.array([cc.min(), cc.max() + 1]) * (x1 - x0) / cols
        gy = y0 + np.array([rr.min(), rr.max() + 1]) * (y1 - y0) / rows
        half = max(half, 0.6 * float(gx[1] - gx[0]), 0.6 * float(gy[1] - gy[0]))
    window = [cx - half, cx + half, cy - half, cy + half]
    if extent is not None:
        window = [
            max(window[0], extent[0]),
            min(window[1], extent[1]),
            max(window[2], extent[2]),
            min(window[3], extent[3]),
        ]
    axes.set_xlim(window[0], window[1])
    axes.set_ylim(window[3], window[2])
    reach = 0.45 * plan.spec.length_um
    axes.annotate(
        "",
        xy=(cx + reach * normal[0], cy + reach * normal[1]),
        xytext=(cx, cy),
        arrowprops={"arrowstyle": "-|>", "color": colors["risk"], "lw": 1.4},
    )
    axes.text(
        cx + 1.2 * reach * normal[0],
        cy + 1.2 * reach * normal[1],
        "$n_L$",
        color=colors["risk"],
        fontsize=8,
        ha="center",
        va="center",
    )
    end = np.array([cx, cy]) + 0.5 * plan.spec.length_um * long_axis
    axes.plot([end[0]], [end[1]], "D", ms=4, color=colors["caution"])
    fit_text = f", clearance {placement.margin_um:.1f} µm" if placement.fits else ", does NOT fit"
    axes.text(
        0.98,
        0.98,
        f"{_math('lamella_azimuth')} = {plan.geometry.theta_sample_deg:.1f}°   "
        f"{_math('lamella_image_azimuth')} = {plan.theta_image_deg:.1f}°\n"
        f"{_math('lamella_length')} x {_math('lamella_width')} = "
        f"{plan.spec.length_um:g} x {plan.spec.width_um:g} µm{fit_text}\n"
        "diamond: lift-out end (+$t_L$)",
        transform=axes.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        color=colors["ink"],
        bbox={
            "boxstyle": "round,pad=0.3",
            "fc": colors["paper"],
            "ec": colors["muted"],
            "alpha": 0.85,
        },
    )
    x_lo, x_hi = sorted(axes.get_xlim())
    y_lo, y_hi = sorted(axes.get_ylim())
    span = x_hi - x_lo
    from pytex.plotting.figure import add_scale_bar

    bar = add_scale_bar(
        axes,
        _nice_length(span),
        units="µm",
        location="lower right",
        color=colors["ink"],
        fontsize=7,
    )
    bar.patch.set_facecolor(colors["paper"])
    del y_lo, y_hi
    _sample_gizmo(axes, plan, colors)
    axes.set_aspect("equal")
    axes.set_xlabel("scan x (µm)")
    axes.set_ylabel("scan y (µm)")
    title = f"Grain {plan.grain_id}: " if plan.grain_id is not None else ""
    axes.set_title(f"{title}lamella site on the {background_label} map", fontsize=10)


# --------------------------------------------------------------------------- #
# 3. Section schematic
# --------------------------------------------------------------------------- #


def draw_section_schematic(figure: Any, plan: Any, *, theme: str = "light") -> None:
    """The lamella seen along its long axis: trenches, slab, depth, thickness, eps*.

    The horizontal axis is ``n_L`` (in the surface), the vertical axis ``Z_s``.
    The target member ``d*`` lies exactly in this plane, at ``eps*`` from
    ``n_L``, which is the whole of section 2 in one picture.
    """

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1)
    _style_axes(axes, colors)
    spec = plan.spec
    width, depth = spec.width_um, spec.depth_um
    trench = max(1.2 * depth * 0.6, 3.0 * width)
    thin = max(spec.final_thickness_nm / 1000.0, 0.02)
    # Surface and material.
    left, right = -trench - width, trench + width
    axes.fill_between([left, right], [-1.35 * depth] * 2, [0.0, 0.0], color=colors["surface"], lw=0)
    # The two trenches, sloped as ion milling leaves them.
    for side in (-1, 1):
        inner = side * width / 2.0
        outer = side * (width / 2.0 + trench)
        axes.fill(
            [inner, outer, outer, inner],
            [0.0, 0.0, -0.05 * depth, -depth],
            color=colors["paper"],
            lw=0,
        )
    # The slab, then its thinned window.
    axes.fill(
        [-width / 2, width / 2, width / 2, -width / 2],
        [0, 0, -depth, -depth],
        color=colors["slab"],
        lw=1.0,
        ec=colors["primary"],
    )
    axes.fill(
        [-thin / 2, thin / 2, thin / 2, -thin / 2],
        [-0.1 * depth, -0.1 * depth, -0.7 * depth, -0.7 * depth],
        color=colors["accent"],
        lw=0,
    )
    axes.plot([left, right], [0, 0], color=colors["ink"], lw=1.0)
    # Ion beam down, beam direction along n_L.
    axes.annotate(
        "",
        xy=(0.0, 0.05 * depth),
        xytext=(0.0, 0.55 * depth),
        arrowprops={"arrowstyle": "-|>", "color": colors["caution"], "lw": 1.6},
    )
    axes.text(
        0.12 * trench, 0.4 * depth, "ion beam, mills along $-Z_s$", color=colors["ink"], fontsize=8
    )
    y_axis = -0.4 * depth
    axes.annotate(
        "",
        xy=(0.95 * right, y_axis),
        xytext=(0.0, y_axis),
        arrowprops={"arrowstyle": "-|>", "color": colors["risk"], "lw": 1.4},
    )
    axes.text(
        0.95 * right,
        y_axis - 0.08 * depth,
        "$n_L$ (TEM beam)",
        color=colors["risk"],
        fontsize=8,
        ha="right",
        va="top",
    )
    eps = math.radians(plan.eps_deg)
    sign = plan.geometry.out_of_plane_sign
    reach = 0.8 * right
    axes.annotate(
        "",
        xy=(reach * math.cos(eps), y_axis + sign * reach * math.sin(eps)),
        xytext=(0.0, y_axis),
        arrowprops={"arrowstyle": "-|>", "color": colors["violet"], "lw": 1.4},
    )
    axes.text(
        reach * math.cos(eps),
        y_axis + sign * reach * math.sin(eps) + sign * 0.05 * depth,
        f"{_member_math(plan)}, {_math('lamella_residual')} = {plan.eps_deg:.1f}°",
        color=colors["violet"],
        fontsize=8,
        ha="right",
        va="bottom" if sign > 0 else "top",
    )
    # Dimension labels.
    axes.annotate(
        "",
        xy=(-width / 2, -1.12 * depth),
        xytext=(width / 2, -1.12 * depth),
        arrowprops={"arrowstyle": "<->", "color": colors["ink"], "lw": 0.8},
    )
    axes.text(
        0.0,
        -1.22 * depth,
        f"{_math('lamella_width')} = {width:g} µm",
        ha="center",
        va="top",
        color=colors["ink"],
        fontsize=8,
    )
    axes.annotate(
        "",
        xy=(left * 0.85, 0.0),
        xytext=(left * 0.85, -depth),
        arrowprops={"arrowstyle": "<->", "color": colors["ink"], "lw": 0.8},
    )
    axes.text(
        left * 0.83,
        -0.5 * depth,
        f"{_math('mill_depth')} = {depth:g} µm",
        color=colors["ink"],
        fontsize=8,
        va="center",
    )
    axes.text(
        thin,
        -0.78 * depth,
        f"thinned to {spec.final_thickness_nm:g} nm",
        color=colors["ink"],
        fontsize=7,
        ha="left",
    )
    axes.set_xlim(left, right)
    axes.set_ylim(-1.4 * depth, 0.7 * depth)
    axes.set_aspect("equal")
    axes.set_xlabel("along $n_L$ (µm)")
    axes.set_ylabel("along $Z_s$ (µm)")
    axes.set_title("Section across the lamella (viewed along $t_L$)", fontsize=10)


# --------------------------------------------------------------------------- #
# 4. The stereogram
# --------------------------------------------------------------------------- #


def _stereo(vectors: np.ndarray) -> np.ndarray:
    """Upper-hemisphere stereographic projection about ``Z_s``."""

    vectors = np.asarray(vectors, dtype=np.float64)
    return vectors[..., :2] / (1.0 + vectors[..., 2:3])


def draw_stereogram(
    figure: Any, plan: Any, *, guaranteed_radius_deg: float | None = None, theme: str = "light"
) -> None:
    """Every orbit member against the circle of achievable beam directions.

    Upper-hemisphere stereographic projection of the sample frame with ``Z_s``
    at the centre, so the **primitive circle is the surface plane** --- the
    great circle every vertically milled lamella's normal lies on. A member's
    distance from the primitive is its ``eps``; the band the holder can close
    for every mounting rotation is shaded, and the chosen member's ``eps*`` is
    drawn as the short radial arc to the circle.
    """

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1)
    _style_axes(axes, colors)
    axes.set_aspect("equal")
    axes.axis("off")
    t = np.linspace(0.0, 2.0 * math.pi, 361)
    if guaranteed_radius_deg is not None and guaranteed_radius_deg > 0.0:
        inner = math.tan(math.radians(90.0 - guaranteed_radius_deg) / 2.0)
        ring_x = np.concatenate([np.cos(t), inner * np.cos(t[::-1])])
        ring_y = np.concatenate([np.sin(t), inner * np.sin(t[::-1])])
        axes.fill(ring_x, ring_y, color=colors["ok"], alpha=0.18, lw=0)
        axes.plot(inner * np.cos(t), inner * np.sin(t), color=colors["ok"], lw=0.8, ls="--")
        axes.text(
            0.0,
            -inner + 0.03,
            f"guaranteed band, {_math('lamella_residual')} ≤ {guaranteed_radius_deg:.0f}°",
            color=colors["ok"],
            fontsize=7,
            ha="center",
            va="bottom",
        )
    for eps in (30.0, 60.0):
        radius = math.tan(math.radians(90.0 - eps) / 2.0)
        axes.plot(radius * np.cos(t), radius * np.sin(t), color=colors["grid"], lw=0.5)
    axes.plot(np.cos(t), np.sin(t), color=colors["primary"], lw=1.6)
    axes.text(
        0.72,
        0.74,
        "surface plane\n(achievable beam\ndirections)",
        color=colors["primary"],
        fontsize=7,
        ha="left",
    )
    # Orbit members on the upper hemisphere.
    directions = _member_directions(plan)
    if plan.orbit.both_senses:
        # Every member has its reverse in the orbit, so the upper hemisphere
        # already shows each lamella option once.
        shown = directions[directions[:, 2] >= -1e-12]
        hollow = np.zeros(len(shown), dtype=bool)
    else:
        hollow = directions[:, 2] < -1e-12
        shown = np.where(hollow[:, None], -directions, directions)
    for (x, y), is_hollow in zip(_stereo(shown), hollow, strict=True):
        axes.plot(
            x,
            y,
            "o",
            ms=4.5,
            mec=colors["ink"],
            mfc=colors["paper"] if is_hollow else colors["ink"],
            mew=0.8,
        )
    chosen = np.asarray(plan.geometry.direction_sample, dtype=np.float64)
    if chosen[2] < 0:
        chosen = -chosen
    cx, cy = _stereo(chosen)
    axes.plot(cx, cy, "o", ms=9, mfc="none", mec=colors["risk"], mew=1.6)
    rim = np.array([chosen[0], chosen[1]]) / max(float(np.hypot(chosen[0], chosen[1])), 1e-12)
    axes.plot([cx, rim[0]], [cy, rim[1]], color=colors["risk"], lw=1.6)
    theta = math.radians(plan.geometry.theta_sample_deg)
    for angle, label, color in (
        (theta, "$n_L$", colors["risk"]),
        (theta + math.pi / 2, "$t_L$", colors["muted"]),
    ):
        axes.plot([math.cos(angle)], [math.sin(angle)], "s", ms=5, color=color)
        axes.text(
            1.12 * math.cos(angle),
            1.12 * math.sin(angle),
            label,
            color=color,
            fontsize=8,
            ha="center",
            va="center",
        )
    axes.text(
        -1.28,
        -1.28,
        f"selected {_member_math(plan)}\n{_math('lamella_residual')} = {plan.eps_deg:.2f}° "
        "from the surface plane",
        color=colors["risk"],
        fontsize=8,
        ha="left",
        va="bottom",
    )
    if not plan.orbit.both_senses:
        axes.text(
            1.28,
            -1.28,
            "open: lower hemisphere,\nshown reversed",
            color=colors["muted"],
            fontsize=7,
            ha="right",
            va="bottom",
        )
    for label, (x, y) in (("$X_s$", (1.22, 0.0)), ("$Y_s$", (0.0, 1.22))):
        axes.text(x, y, label, color=colors["ink"], fontsize=8, ha="center", va="center")
    axes.plot(0, 0, "+", color=colors["ink"], ms=7)
    axes.text(0.03, 0.03, "$Z_s$", color=colors["ink"], fontsize=8)
    axes.set_xlim(-1.32, 1.32)
    axes.set_ylim(-1.32, 1.32)
    axes.set_title(
        f"Orbit of {_family_math(plan)}, upper hemisphere of the sample frame", fontsize=10
    )


def _member_directions(plan: Any) -> np.ndarray:
    """Sample-frame directions of every orbit member, as the plan computed them."""

    directions = plan.member_directions_sample
    if directions is None:
        raise ValueError("The plan does not carry its orbit's sample-frame directions.")
    return np.asarray(directions, dtype=np.float64)


# --------------------------------------------------------------------------- #
# 5. The phi feasibility plot
# --------------------------------------------------------------------------- #


def reach_curve_deg(
    envelope: Any,
    phi_deg: np.ndarray,
    *,
    rise_sign: int = 1,
    flip: int = 1,
    resolution_deg: float = 0.1,
) -> np.ndarray:
    """Largest residual the holder removes at each mounting rotation, exactly.

    Scans the residual upwards in ``resolution_deg`` steps and records, for each
    ``phi``, the last value whose exact tilt stays inside the envelope.
    """

    from pytex.fib.planning import _curve_inside

    phi = np.asarray(phi_deg, dtype=np.float64)
    reach = np.zeros(phi.shape)
    alive = np.ones(phi.shape, dtype=bool)
    for eps in np.arange(resolution_deg, 89.9, resolution_deg):
        inside = _curve_inside(envelope, float(eps), rise_sign, phi, flip)
        alive &= inside
        reach = np.where(alive, eps, reach)
        if not alive.any():
            break
    return reach


def draw_phi_feasibility(figure: Any, plan: Any, *, envelope: Any, theme: str = "light") -> None:
    """Reachable and unreachable arcs of the unknown mounting rotation.

    Polar axes: angle = ``phi``, radius = residual tilt. The shaded region is
    what the holder can remove at each ``phi`` (the envelope's exact reach),
    the dashed circle the guaranteed disc, the solid circle ``eps*`` and the
    dotted circle its upper bound; the thick arcs at ``eps*`` show where the
    front and back branches actually reach the target.
    """

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1, projection="polar")
    axes.set_facecolor(colors["paper"])
    phi = np.linspace(0.0, 360.0, 721)
    reach = reach_curve_deg(envelope, phi, rise_sign=plan.geometry.out_of_plane_sign)
    theta = np.radians(phi)
    axes.fill_between(theta, 0.0, reach, color=colors["accent"], alpha=0.2, lw=0)
    axes.plot(theta, reach, color=colors["accent"], lw=1.0, label="holder reach")
    guaranteed = float(reach.min())
    axes.plot(
        theta,
        np.full_like(theta, guaranteed),
        color=colors["ok"],
        ls="--",
        lw=1.0,
        label=f"guaranteed disc ({guaranteed:.1f}°)",
    )
    eps = plan.eps_deg
    upper = plan.budget.eps_upper_deg
    axes.plot(
        theta,
        np.full_like(theta, eps),
        color=colors["ink"],
        lw=0.8,
        label=f"{_math('lamella_residual')} = {eps:.1f}°",
    )
    axes.plot(
        theta,
        np.full_like(theta, upper),
        color=colors["muted"],
        lw=0.8,
        ls=":",
        label=f"upper bound {upper:.1f}°",
    )
    sweep = plan.sweep
    sweep_theta = np.radians(np.asarray(sweep.phi_deg))
    for mask, offset, color, label in (
        (np.asarray(sweep.reachable_front), 0.0, colors["primary"], "front: reachable"),
        (np.asarray(sweep.reachable_back), 1.2, colors["violet"], "back: reachable"),
    ):
        radius = np.where(mask, eps + offset, np.nan)
        axes.plot(sweep_theta, radius, color=color, lw=3.0, label=label, solid_capstyle="butt")
    top = max(float(reach.max()), upper) * 1.15 + 2.0
    axes.set_ylim(0.0, top)
    axes.set_theta_zero_location("E")
    axes.set_theta_direction(1)
    axes.tick_params(colors=colors["muted"], labelcolor=colors["ink"], labelsize=7)
    axes.grid(color=colors["grid"], lw=0.5)
    axes.spines["polar"].set_color(colors["muted"])
    axes.set_title(
        f"Mounting rotation {_math('mount_rotation')}: reachable for "
        f"{100 * sweep.fraction:.0f}% ({plan.feasibility.value})",
        fontsize=10,
        color=colors["ink"],
    )
    legend = axes.legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), fontsize=7, frameon=False)
    for text in legend.get_texts():
        text.set_color(colors["ink"])


# --------------------------------------------------------------------------- #
# 6. The preparability map
# --------------------------------------------------------------------------- #


def draw_preparability_map(
    figure: Any,
    raster: Any,
    *,
    plans: Sequence[Any] = (),
    guaranteed_radius_deg: float | None = None,
    theme: str = "light",
) -> None:
    """``eps*`` at every point of the map, with the ranked sites on top.

    Sequential colour map (``viridis``, perceptually uniform and legible in
    greyscale) with an explicit colour scale in degrees; low is good. The
    contour marks the guaranteed radius of the holder.
    """

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1)
    _style_axes(axes, colors)
    x0, x1, y0, y1 = raster.extent_um
    finite = raster.eps_deg[np.isfinite(raster.eps_deg)]
    vmax = float(max(finite.max(), 1.0)) if finite.size else 45.0
    image = axes.imshow(
        raster.eps_deg,
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        extent=(x0, x1, y1, y0),
        origin="upper",
        interpolation="nearest",
    )
    if (
        guaranteed_radius_deg is not None
        and finite.size
        and finite.min() < guaranteed_radius_deg < finite.max()
    ):
        rows, cols = raster.eps_deg.shape
        xs = np.linspace(x0, x1, cols + 1)[:-1] + (x1 - x0) / cols / 2
        ys = np.linspace(y0, y1, rows + 1)[:-1] + (y1 - y0) / rows / 2
        axes.contour(
            xs,
            ys,
            np.nan_to_num(raster.eps_deg, nan=90.0),
            levels=[guaranteed_radius_deg],
            colors=[colors["paper"]],
            linewidths=0.8,
            linestyles="--",
        )
    for rank, plan in enumerate(plans, start=1):
        placement = plan.placement
        if placement is None:
            continue
        corners = np.vstack([placement.corners_scan_um, placement.corners_scan_um[:1]])
        color = colors["risk"] if rank == 1 else colors["caution"]
        axes.plot(corners[:, 0], corners[:, 1], color=color, lw=1.4 if rank == 1 else 0.9)
        cx, cy = placement.center_scan_um
        axes.text(
            cx, cy, str(rank), color=color, fontsize=8, ha="center", va="center", fontweight="bold"
        )
    bar = figure.colorbar(image, ax=axes, fraction=0.046, pad=0.03)
    bar.set_label(f"{_math('lamella_residual')} (°)", color=colors["ink"])
    bar.ax.tick_params(colors=colors["muted"], labelcolor=colors["ink"])
    bar.outline.set_edgecolor(colors["muted"])
    axes.set_aspect("equal")
    axes.set_xlim(x0, x1)
    axes.set_ylim(y1, y0)
    axes.set_xlabel("scan x (µm)")
    axes.set_ylabel("scan y (µm)")
    axes.set_title(
        f"Preparability for {raster.target_text}: residual tilt at every point", fontsize=10
    )


# --------------------------------------------------------------------------- #
# 7. Predicted SAED
# --------------------------------------------------------------------------- #


def draw_predicted_saed(
    figure: Any, plan: Any, *, theme: str = "light", max_index: int = 4
) -> None:
    """The kinematic SAED pattern at the achieved zone axis (the chosen member)."""

    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.saed import generate_saed_pattern

    colors = _palette(theme)
    _prepare(figure, colors)
    axes = figure.add_subplot(1, 1, 1)
    _style_axes(axes, colors)
    phase = plan.orbit.phase
    zone = ZoneAxis(np.asarray(plan.geometry.member_indices), phase)
    pattern = generate_saed_pattern(phase, zone, max_index=max_index, label_limit=12)
    spots = list(pattern.spots)
    if spots:
        coordinates = np.vstack([spot.detector_coordinates for spot in spots])
        intensities = np.array([spot.intensity for spot in spots], dtype=np.float64)
        peak = float(intensities.max()) if intensities.max() > 0 else 1.0
        sizes = 8.0 + 90.0 * intensities / peak
        axes.scatter(coordinates[:, 0], coordinates[:, 1], s=sizes, color=colors["ink"], lw=0)
        axes.scatter([0.0], [0.0], s=120.0, color=colors["primary"], lw=0)
        extent = float(np.abs(coordinates).max()) * 1.15 + 1.0
        offset = 0.035 * extent
        order = np.argsort(-intensities)[:12]
        for index in order:
            spot = spots[int(index)]
            if not np.any(spot.miller_indices):
                continue
            x, y = spot.detector_coordinates
            axes.text(
                x,
                y + offset,
                spot.label,
                color=colors["muted"],
                fontsize=6,
                ha="center",
                va="bottom",
            )
        axes.set_xlim(-extent, extent)
        axes.set_ylim(-extent, extent)
    axes.set_aspect("equal")
    axes.set_xlabel("detector u (mm)")
    axes.set_ylabel("detector v (mm)")
    axes.set_title(f"Predicted SAED down {_member_math(plan)} (kinematic)", fontsize=10)


# --------------------------------------------------------------------------- #
# Script entry point
# --------------------------------------------------------------------------- #

#: The figure kinds :func:`fib_figure` draws.
FIB_FIGURE_KINDS: tuple[str, ...] = (
    "plan_view",
    "section",
    "stereogram",
    "phi_feasibility",
    "preparability",
    "saed",
)


def fib_figure(
    kind: str,
    draw_arguments: dict[str, Any],
    *,
    theme: str = "light",
    size_in: tuple[float, float] = (6.4, 4.8),
) -> Any:
    """Draw one figure into a new bare :class:`matplotlib.figure.Figure`.

    No ``pyplot`` figure manager is involved, so nothing is left open.

    Parameters
    ----------
    kind : str
        One of :data:`FIB_FIGURE_KINDS`.
    draw_arguments : dict
        Arguments of the matching ``draw_*`` function after ``figure``.
    theme : {"light", "dark"}
    size_in : (width, height)

    Returns
    -------
    matplotlib.figure.Figure
    """

    from matplotlib.figure import Figure

    functions: dict[str, Callable[..., None]] = {
        "plan_view": draw_plan_view,
        "section": draw_section_schematic,
        "stereogram": draw_stereogram,
        "phi_feasibility": draw_phi_feasibility,
        "preparability": draw_preparability_map,
        "saed": draw_predicted_saed,
    }
    if kind not in functions:
        raise ValueError(f"kind must be one of {FIB_FIGURE_KINDS}; got {kind!r}.")
    figure = Figure(figsize=size_in, layout="constrained")
    functions[kind](figure, theme=theme, **draw_arguments)
    return figure
