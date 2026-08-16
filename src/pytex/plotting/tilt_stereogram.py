"""The tilt-navigation stereogram: where the beam is, where it is going, and how.

An annotated stereographic projection that answers, at a glance, the questions an
operator cannot check by eye:

- where the beam sits in the crystal now, and where the requested target is;
- **the trajectory, drawn as a series of dots** — one per sampled stage position,
  so the spacing itself shows how fast the beam moves through the crystal, and
  the dots grow toward the target so the direction of travel survives a
  greyscale reprint;
- which region of the crystal the holder can reach at all, drawn as the exact
  curve boundary rather than a sampled blob;
- every symmetry-equivalent version of the target, distinguished as reachable or
  not, which is usually the surprise: a target that looks hopeless is routine
  because an equivalent sits inside the envelope;
- the principal zone axes, so the figure reads as a crystallographic stereogram
  and not merely as a plot;
- and the sense of each tilt knob, drawn from the engine's own forward model.

Why two panels
--------------
A double-tilt holder reaches roughly eight percent of the sphere, so everything
interesting happens inside a small patch of a full stereogram. A single panel
must then choose between crystallographic context and legible detail. It does not
have to: the default figure shows both, an overview that places the move among
the low-index poles and a detail view zoomed to the reachable region where the
stage angles and the tilt-axis senses are readable.

Why the reachable region draws exactly
--------------------------------------
The beam direction in holder coordinates is
``(-cos a sin b, sin a, cos a cos b)`` — spherical coordinates whose pole is the
beta axis. Curves of constant beta are therefore great circles through that pole
and curves of constant alpha are small circles about it, so the boundary of the
accessible region is four exact circular arcs, not a sampled outline.

Independence from the engine
----------------------------
The renderer consumes `TiltPath` samples and does **no kinematics of its own**;
it plots what the engine produced. A drawing that agreed with the engine because
it re-implemented the same formula would be evidence about nothing. Independence
is obtained in the test suite, which recomputes the plotted trajectory by
accumulating small rotations and asserts agreement.

The module follows the repository's builder/renderer split: `build_*_figure_spec`
returns a declarative `FigureSpec2D` that can be asserted against structurally,
and `plot_*` renders it. Tests check the spec, never the pixels.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal

import numpy as np

from pytex.core._arrays import normalize_vector
from pytex.core.lattice import Phase
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.sphere import project_directions
from pytex.diffraction.stereonets import projection_boundary_radius
from pytex.plotting._render import (
    FigureSpec2D,
    LineLayer2D,
    MarkerLayer2D,
    TextLayer2D,
    render_figure_spec,
)
from pytex.plotting.spherical import build_wulff_net_figure_spec
from pytex.tem.navigation import TiltPlanReport, TiltSolution
from pytex.tem.path import TiltPath
from pytex.tem.stage import StageModel

__all__ = [
    "TILT_STEREOGRAM_COLORS",
    "build_tilt_stereogram_figure_spec",
    "plot_tilt_stereogram",
]

#: Colours for the stereogram, drawn from the repository visualization tokens.
#:
#: Semantics rather than decoration: teal marks where the specimen *is* (a
#: measurement), green a reachable destination, rose an unreachable one, core
#: blue the trajectory, violet the crystallographic scaffolding, and amber the
#: instrument axes. Keeping the assignment here rather than inline is what lets a
#: caller restyle the figure without re-deriving which colour meant what.
TILT_STEREOGRAM_COLORS: dict[str, str] = {
    "ink": "#07122f",
    "muted": "#40506f",
    "current": "#0f9f9f",
    "target": "#16a34a",
    "unreachable": "#e11d48",
    "trajectory": "#2563eb",
    "poles": "#7c3aed",
    "axes": "#f59e0b",
    "region": "#2563eb",
}

#: Low-index directions labelled on the overview panel by default.
#:
#: The families an operator navigates by. Every phase gets the same list expanded
#: through its own symmetry, which is the right behaviour: the poles worth
#: labelling are the low-index ones of *that* lattice.
_DEFAULT_POLE_FAMILIES: tuple[tuple[int, int, int], ...] = (
    (1, 0, 0),
    (1, 1, 0),
    (1, 1, 1),
    (2, 1, 0),
    (2, 1, 1),
    (3, 1, 1),
)

#: Wulff-net styling for this figure: coarse and pale.
#:
#: A full two-degree minor net is right for reading angles off a stereogram by
#: hand, and wrong here — it competes with the trajectory, which is the subject.
#: The net's job in this figure is to establish that the plot *is* a stereogram
#: and to give the eye a coarse angular reference.
_NET_STYLE: dict[str, Any] = {
    "show_minor_grid": False,
    "major_step_deg": 15.0,
    "net_major_color": "#d5dbe6",
    "net_major_linewidth": 0.55,
    "net_alpha": 0.9,
}

#: Minimum separation, in projected units, between two labelled poles.
#:
#: Pole labels crowd badly near the rim, where the projection compresses. Rather
#: than shrink the font until nothing is legible, labels are dropped where they
#: would collide; the marker stays, so no pole disappears from the figure.
_LABEL_SEPARATION = 0.085

#: Which panel a spec describes.
StereogramView = Literal["overview", "detail"]


def _radial_offset(point: np.ndarray, distance: float) -> np.ndarray:
    """Offset pointing away from the projection centre.

    Pushing a label outward keeps it clear of the marker it names and of the
    dense middle of the stereogram, and it degrades gracefully at the origin,
    where any direction is as good as another.
    """

    norm = float(np.linalg.norm(point))
    if norm < 1e-9:
        return np.array([0.0, distance])
    return np.asarray(point / norm * distance, dtype=np.float64)


def _separated_offsets(
    first: np.ndarray, second: np.ndarray, distance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Label offsets that push two nearby points apart along their own axis.

    The start and target markers are frequently close — a short hop is the common
    case — so offsetting both radially would stack their labels. Pushing each one
    *away from the other* guarantees separation wherever they sit, and reduces to
    a sensible vertical split when they coincide.
    """

    separation = second - first
    norm = float(np.linalg.norm(separation))
    if norm < 1e-9:
        return np.array([0.0, -distance]), np.array([0.0, distance])
    unit = separation / norm
    return -unit * distance, unit * distance


def _pole_directions(
    phase: Phase, families: Sequence[tuple[int, int, int]], max_poles: int
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Symmetry-expanded low-index poles with their formatted labels.

    Labels come from `pytex.core.notation` alone, per the repository rule that
    crystallographic notation is never formatted inline.
    """

    direct = phase.lattice.direct_basis().matrix
    operators = np.asarray(phase.symmetry.operators, dtype=np.float64)

    directions: list[np.ndarray] = []
    labels: list[str] = []
    for family in families:
        base = np.asarray(family, dtype=np.float64) @ direct.T
        base = base / float(np.linalg.norm(base))
        for image in np.einsum("nij,j->ni", operators, base):
            if image[2] < 0.0:
                image = -image
            if any(float(np.linalg.norm(image - kept)) < 1e-6 for kept in directions):
                continue
            directions.append(np.asarray(image, dtype=np.float64))
            labels.append(_family_label(family, image, direct))
            if len(directions) >= max_poles:
                break
        if len(directions) >= max_poles:
            break
    if not directions:
        return np.zeros((0, 3), dtype=np.float64), ()
    return np.stack(directions), tuple(labels)


def _family_label(
    family: tuple[int, int, int], image: np.ndarray, direct: np.ndarray
) -> str:
    """The specific ``[uvw]`` label of one symmetry image of a family."""

    indices = np.linalg.solve(direct, image)
    scale = float(np.max(np.abs(indices)))
    if scale <= 0.0:
        return format_direction_indices(list(family))
    scaled = indices / scale * max(abs(value) for value in family)
    rounded = [round(float(value)) for value in scaled]
    if all(value == 0 for value in rounded):
        return format_direction_indices(list(family))
    return format_direction_indices(rounded)


def _beam_grid(stage: StageModel, alphas: np.ndarray, betas: np.ndarray) -> np.ndarray:
    """Beam directions in holder coordinates along a sequence of stage positions."""

    return np.stack(
        [
            stage.beam_direction(float(alpha), float(beta))
            for alpha, beta in zip(alphas, betas, strict=True)
        ]
    )


def _to_crystal(points_holder: np.ndarray, crystal_to_holder: np.ndarray) -> np.ndarray:
    """Holder-frame rows into crystal-frame rows: ``(U^T v)^T = v^T U``."""

    return np.asarray(points_holder @ crystal_to_holder, dtype=np.float64)


def _split_on_wrap(projected: np.ndarray, threshold: float = 0.5) -> list[np.ndarray]:
    """Break a projected trace where it crosses the hemisphere boundary.

    Without this, a trace that leaves the upper hemisphere and re-enters draws a
    straight chord across the disc, which is not a curve on the sphere at all.
    """

    if projected.shape[0] < 2:
        return [projected]
    breaks = np.nonzero(np.linalg.norm(np.diff(projected, axis=0), axis=1) > threshold)[0]
    return [segment for segment in np.split(projected, breaks + 1) if segment.shape[0] >= 2]


def _reachable_region_layers(
    stage: StageModel, crystal_to_holder: np.ndarray, method: str, color: str
) -> tuple[LineLayer2D, ...]:
    """The envelope boundary as exact arcs on the crystal stereogram."""

    alpha_min, alpha_max, beta_min, beta_max = stage.envelope.bounds()
    layers: list[LineLayer2D] = []
    samples = 181
    edges: list[tuple[np.ndarray, np.ndarray]] = []
    for alpha_edge in (alpha_min, alpha_max):
        betas = np.linspace(beta_min, beta_max, samples)
        edges.append((np.full_like(betas, alpha_edge), betas))
    for beta_edge in (beta_min, beta_max):
        alphas = np.linspace(alpha_min, alpha_max, samples)
        edges.append((alphas, np.full_like(alphas, beta_edge)))

    for index, (alphas, betas) in enumerate(edges):
        crystal = _to_crystal(_beam_grid(stage, alphas, betas), crystal_to_holder)
        for segment in _split_on_wrap(project_directions(crystal, method=method)):
            layers.append(
                LineLayer2D(
                    points=segment,
                    label="reachable with this holder" if index == 0 else None,
                    color=color,
                    linewidth=1.5,
                    linestyle="-",
                    alpha=0.9,
                )
            )
    return tuple(layers)


def _tilt_axis_layers(
    stage: StageModel,
    crystal_to_holder: np.ndarray,
    current_alpha: float,
    current_beta: float,
    method: str,
    color: str,
) -> tuple[tuple[LineLayer2D, ...], tuple[TextLayer2D, ...]]:
    """Short arcs showing which way each knob moves the beam, from the current point.

    This answers "what does positive alpha actually do to my pattern?" directly on
    the figure. Each arc is the locus traced by moving one axis alone, generated
    through the engine's forward model rather than from an assumed sense, so a
    reversed sign convention is visible rather than hidden.
    """

    span = 15.0
    lines: list[LineLayer2D] = []
    texts: list[TextLayer2D] = []
    excursions = (
        ("+α", np.linspace(current_alpha, current_alpha + span, 40), np.full(40, current_beta)),
        ("+β", np.full(40, current_alpha), np.linspace(current_beta, current_beta + span, 40)),
    )
    for index, (label, alphas, betas) in enumerate(excursions):
        crystal = _to_crystal(_beam_grid(stage, alphas, betas), crystal_to_holder)
        projected = project_directions(crystal, method=method)
        lines.append(
            LineLayer2D(
                points=projected,
                label="tilt-axis sense" if index == 0 else None,
                color=color,
                linewidth=2.6,
                linestyle="-",
                alpha=0.9,
            )
        )
        step = projected[-1] - projected[-2]
        norm = float(np.linalg.norm(step))
        direction = step / norm if norm > 1e-12 else np.array([1.0, 0.0])
        texts.append(
            TextLayer2D(
                position=projected[-1] + direction * 0.05,
                text=label,
                color=color,
                fontsize=12.0,
                bbox_facecolor="#ffffff",
                bbox_edgecolor=color,
                bbox_alpha=0.95,
                zorder=8.0,
            )
        )
    return tuple(lines), tuple(texts)


def _equivalent_points(
    report: TiltPlanReport,
    stage: StageModel,
    crystal_to_holder: np.ndarray,
    method: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Projected symmetry equivalents, split into reachable and out of range."""

    from pytex.tem.navigation import solve_tilts_for_direction

    operators = np.asarray(report.current.phase.symmetry.operators, dtype=np.float64)
    base = normalize_vector(report.target.unit_vector)
    images = np.einsum("nij,j->ni", operators, base)
    images = np.vstack([images, -images])

    reachable: list[np.ndarray] = []
    unreachable: list[np.ndarray] = []
    seen: list[np.ndarray] = []
    for image in images:
        if any(float(np.linalg.norm(image - kept)) < 1e-8 for kept in seen):
            continue
        seen.append(image)
        branches = solve_tilts_for_direction(crystal_to_holder @ image, allow_reverse=True)
        inside = any(stage.envelope.contains(alpha, beta) for alpha, beta in branches)
        (reachable if inside else unreachable).append(image)

    def stack(items: list[np.ndarray]) -> np.ndarray:
        if not items:
            return np.zeros((0, 2), dtype=np.float64)
        return project_directions(np.stack(items), method=method)

    return stack(reachable), stack(unreachable)


def _trajectory_annotations(
    path: TiltPath, trajectory: np.ndarray, count: int
) -> list[TextLayer2D]:
    """Stage-angle labels at a few points along the trajectory.

    Enough to read the schedule off the figure, few enough not to bury it. Each
    label is offset perpendicular to the local direction of travel so it sits
    beside the dotted arc rather than on the next dot.
    """

    total = trajectory.shape[0]
    if total < 5 or count < 1:
        return []
    indices = sorted({round(f * (total - 1)) for f in np.linspace(0.2, 0.8, count)})
    annotations: list[TextLayer2D] = []
    for index in indices:
        sample = path.samples[index]
        previous = trajectory[max(index - 1, 0)]
        following = trajectory[min(index + 1, total - 1)]
        step = following - previous
        norm = float(np.linalg.norm(step))
        normal = (
            np.array([-step[1], step[0]]) / norm if norm > 1e-9 else np.array([0.0, 1.0])
        )
        annotations.append(
            TextLayer2D(
                position=trajectory[index] + normal * 0.055,
                text=(
                    f"α {sample.position.alpha_deg:+.0f}°, "
                    f"β {sample.position.beta_deg:+.0f}°"
                ),
                color=TILT_STEREOGRAM_COLORS["trajectory"],
                fontsize=8.0,
                ha="center",
                bbox_facecolor="#ffffff",
                bbox_edgecolor="none",
                bbox_alpha=0.8,
                zorder=7.0,
            )
        )
    return annotations


def build_tilt_stereogram_figure_spec(
    report: TiltPlanReport,
    stage: StageModel,
    *,
    solution: TiltSolution | None = None,
    view: StereogramView = "overview",
    method: str = "stereographic",
    title: str | None = None,
    show_wulff_net: bool = True,
    show_reachable_region: bool = True,
    show_equivalents: bool = True,
    show_poles: bool = True,
    show_tilt_axes: bool = True,
    max_poles: int = 26,
    pole_families: Sequence[tuple[int, int, int]] = _DEFAULT_POLE_FAMILIES,
    theme: str = "journal",
) -> FigureSpec2D:
    """Assemble the declarative spec for one panel of the tilt stereogram.

    Purpose
    -------
    Builds the figure as data, so its scientific content can be asserted in tests
    — that the trajectory has one dot per path sample, that the target marker
    sits where the engine put it, that unreachable equivalents carry the
    rejection colour — without rendering anything or comparing images.

    When to use
    -----------
    Call this when a panel must be inspected, restyled, or embedded in a larger
    layout. Call :func:`plot_tilt_stereogram` for the finished two-panel figure.

    Parameters
    ----------
    report : TiltPlanReport
        The result of `pytex.tem.plan_tilt_to_zone_axis`. Supplies the current
        state, the target, the symmetry orbit and the ranked solutions.
    stage : StageModel
        Used for the reachable region and the tilt-axis arcs. The trajectory
        itself comes from the report's path, never from re-running kinematics
        here.
    solution : TiltSolution, optional
        Which solution to draw. Defaults to the best; pass an alternative to show
        a competing ambiguity hypothesis.
    view : {"overview", "detail"}, default "overview"
        ``"overview"`` fills the whole projection disc and labels the principal
        poles, placing the move in crystallographic context. ``"detail"`` zooms
        to the reachable region, where the stage angles and tilt-axis senses are
        legible. A double-tilt holder reaches only a few percent of the sphere,
        so one panel cannot do both jobs well.
    method : {"stereographic", "equal_area"}, default "stereographic"
        Equal-angle is the crystallographic default and preserves angles, which
        is what makes a stereogram readable as geometry.
    title : str, optional
    show_wulff_net, show_reachable_region, show_equivalents, show_poles : bool
    show_tilt_axes : bool, default True
        Draw the arcs showing which way each knob moves the beam. Honoured on
        the detail view, where they are legible.
    max_poles : int, default 26
        Cap on labelled poles, so a low-symmetry phase does not produce a figure
        that is more text than geometry.
    pole_families : sequence of (int, int, int)
    theme : str, default "journal"

    Returns
    -------
    FigureSpec2D

    Notes
    -----
    The trajectory is drawn as **dots, one per sampled stage position**, not as a
    continuous line: their spacing carries information. Where they crowd, the
    beam is moving slowly through the crystal for a given change of stage angle,
    which is what happens as the beta-axis pole is approached.
    """

    if solution is None and report.solutions:
        solution = report.best()
    if solution is None:
        solution = report.nearest_approach
    if solution is None:
        raise ValueError(
            "The report contains neither a reachable solution nor a nearest approach, "
            "so there is no trajectory to draw. That happens only when the holder "
            "envelope is empty; check the stage model."
        )

    phase = report.current.phase
    crystal_to_holder = report.current.matrix @ solution.family.operator
    radius = projection_boundary_radius(method)
    is_detail = view == "detail"

    current_point = project_directions(
        report.current.beam_direction_crystal(stage), method=method
    )
    target_point = project_directions(solution.orbit_member, method=method)
    path: TiltPath | None = solution.path
    trajectory: np.ndarray | None = (
        project_directions(path.beam_directions_crystal(), method=method)
        if path is not None
        else None
    )

    # The detail limits must be known *before* the annotation layers are built:
    # the renderer expands the axes to fit every text layer, so a label for a
    # pole far outside the zoom would silently undo the zoom.
    if is_detail:
        xlim, ylim = _detail_limits(
            stage, crystal_to_holder, method, current_point, target_point, trajectory
        )
        boundary_radius = None
    else:
        xlim = (-radius - 0.14, radius + 0.14)
        ylim = (-radius - 0.14, radius + 0.14)
        boundary_radius = radius

    def in_view(points: np.ndarray) -> np.ndarray:
        """Mask of projected points inside the panel limits."""

        if points.shape[0] == 0:
            return np.zeros(0, dtype=bool)
        return np.asarray(
            (points[:, 0] >= xlim[0])
            & (points[:, 0] <= xlim[1])
            & (points[:, 1] >= ylim[0])
            & (points[:, 1] <= ylim[1])
        )

    base = (
        build_wulff_net_figure_spec(
            method=method, theme=theme, style_overrides={"spherical": dict(_NET_STYLE)}
        )
        if show_wulff_net
        else FigureSpec2D()
    )
    line_layers: list[LineLayer2D] = list(base.line_layers)
    marker_layers: list[MarkerLayer2D] = []
    text_layers: list[TextLayer2D] = []

    if show_reachable_region:
        line_layers.extend(
            _reachable_region_layers(
                stage, crystal_to_holder, method, TILT_STEREOGRAM_COLORS["region"]
            )
        )

    if show_poles:
        poles, labels = _pole_directions(phase, pole_families, max_poles)
        if poles.shape[0]:
            projected = project_directions(poles, method=method)
            visible = in_view(projected)
            projected = projected[visible]
            labels = tuple(
                label for label, keep in zip(labels, visible, strict=True) if keep
            )
        if poles.shape[0] and projected.shape[0]:
            marker_layers.append(
                MarkerLayer2D(
                    points=projected,
                    marker="+",
                    label="low-index zone axes",
                    facecolors=TILT_STEREOGRAM_COLORS["poles"],
                    edgecolors=None,
                    sizes=30.0,
                    alpha=0.8,
                    linewidths=1.0,
                )
            )
            placed: list[np.ndarray] = []
            for point, label in zip(projected, labels, strict=True):
                anchor = point + _radial_offset(point, 0.045)
                if any(
                    float(np.linalg.norm(anchor - other)) < _LABEL_SEPARATION
                    for other in placed
                ):
                    continue
                placed.append(anchor)
                text_layers.append(
                    TextLayer2D(
                        position=anchor,
                        text=label,
                        color=TILT_STEREOGRAM_COLORS["poles"],
                        fontsize=8.5 if is_detail else 8.0,
                        zorder=5.0,
                    )
                )

    if show_equivalents:
        reachable, unreachable = _equivalent_points(
            report, stage, crystal_to_holder, method
        )
        unreachable = unreachable[in_view(unreachable)]
        reachable = reachable[in_view(reachable)]
        if unreachable.shape[0]:
            marker_layers.append(
                MarkerLayer2D(
                    points=unreachable,
                    marker="o",
                    label="equivalent target, out of range",
                    facecolors="none",
                    edgecolors=TILT_STEREOGRAM_COLORS["unreachable"],
                    sizes=54.0,
                    alpha=0.85,
                    linewidths=1.1,
                )
            )
        if reachable.shape[0]:
            marker_layers.append(
                MarkerLayer2D(
                    points=reachable,
                    marker="o",
                    label="equivalent target, reachable",
                    facecolors="none",
                    edgecolors=TILT_STEREOGRAM_COLORS["target"],
                    # Larger than the target star so it reads as a halo around
                    # the selected one rather than disappearing beneath it.
                    sizes=430.0,
                    alpha=0.95,
                    linewidths=1.8,
                )
            )

    if path is not None and trajectory is not None:
        sizes = (
            np.linspace(10.0, 60.0, trajectory.shape[0])
            if is_detail
            else np.full(trajectory.shape[0], 12.0)
        )
        marker_layers.append(
            MarkerLayer2D(
                points=trajectory,
                marker="o",
                label="beam direction during the tilt (dots grow toward the target)",
                facecolors=TILT_STEREOGRAM_COLORS["trajectory"],
                edgecolors="#ffffff",
                sizes=sizes,
                alpha=0.95,
                linewidths=0.6,
            )
        )
        if is_detail:
            text_layers.extend(_trajectory_annotations(path, trajectory, count=1))

    separation = float(np.linalg.norm(target_point[0] - current_point[0]))
    start_offset, target_offset = _separated_offsets(
        current_point[0], target_point[0], max(0.35 * separation, 0.06)
    )
    reached = solution.verdict.value in ("exact", "within_tolerance")
    target_color = (
        TILT_STEREOGRAM_COLORS["target"] if reached else TILT_STEREOGRAM_COLORS["unreachable"]
    )

    marker_layers.append(
        MarkerLayer2D(
            points=current_point,
            marker="o",
            label="current beam direction",
            facecolors=TILT_STEREOGRAM_COLORS["current"],
            edgecolors=TILT_STEREOGRAM_COLORS["ink"],
            sizes=130.0,
            alpha=1.0,
            linewidths=1.3,
        )
    )
    marker_layers.append(
        MarkerLayer2D(
            points=target_point,
            marker="*",
            label="selected target",
            facecolors=target_color,
            edgecolors=TILT_STEREOGRAM_COLORS["ink"],
            sizes=300.0,
            alpha=1.0,
            linewidths=1.2,
        )
    )

    if is_detail:
        current_label = (
            format_direction_indices(
                [int(value) for value in report.current.current_zone_axis.indices]
            )
            if report.current.current_zone_axis is not None
            else "current"
        )
        target_label = (
            format_direction_indices([int(v) for v in solution.orbit_member_indices])
            if solution.orbit_member_indices is not None
            else format_direction_indices([int(v) for v in report.target.indices])
        )
        text_layers.append(
            TextLayer2D(
                position=current_point[0] + start_offset,
                text=f"start {current_label}",
                color=TILT_STEREOGRAM_COLORS["current"],
                fontsize=10.5,
                bbox_facecolor="#ffffff",
                bbox_edgecolor=TILT_STEREOGRAM_COLORS["current"],
                bbox_alpha=0.95,
                zorder=9.0,
            )
        )
        text_layers.append(
            TextLayer2D(
                position=target_point[0] + target_offset,
                text=f"target {target_label}",
                color=target_color,
                fontsize=10.5,
                bbox_facecolor="#ffffff",
                bbox_edgecolor=target_color,
                bbox_alpha=0.95,
                zorder=9.0,
            )
        )
        if show_tilt_axes:
            axis_lines, axis_texts = _tilt_axis_layers(
                stage,
                crystal_to_holder,
                report.current.position.alpha_deg,
                report.current.position.beta_deg,
                method,
                TILT_STEREOGRAM_COLORS["axes"],
            )
            line_layers.extend(axis_lines)
            text_layers.extend(axis_texts)

    return FigureSpec2D(
        title=title or _panel_title(report, solution, view),
        xlabel="",
        ylabel="",
        xlim=xlim,
        ylim=ylim,
        equal_aspect=True,
        grid=False,
        show_axes=False,
        boundary_circle_radius=boundary_radius,
        boundary_circle_color=TILT_STEREOGRAM_COLORS["ink"],
        boundary_circle_linewidth=1.3,
        boundary_circle_linestyle="-",
        marker_layers=tuple(marker_layers),
        line_layers=tuple(line_layers),
        text_layers=tuple(text_layers),
    )


def _detail_limits(
    stage: StageModel,
    crystal_to_holder: np.ndarray,
    method: str,
    current_point: np.ndarray,
    target_point: np.ndarray,
    trajectory: np.ndarray | None,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Square limits enclosing the reachable region, the path and both endpoints.

    Framed on the reachable region rather than on the trajectory alone, so the
    detail view always shows *how much room is left* around the destination —
    which is the question the operator asks next.
    """

    alpha_min, alpha_max, beta_min, beta_max = stage.envelope.bounds()
    alphas = np.linspace(alpha_min, alpha_max, 40)
    betas = np.linspace(beta_min, beta_max, 40)
    grid_alpha, grid_beta = np.meshgrid(alphas, betas, indexing="ij")
    region = project_directions(
        _to_crystal(
            _beam_grid(stage, grid_alpha.ravel(), grid_beta.ravel()), crystal_to_holder
        ),
        method=method,
    )
    points = [region, current_point, target_point]
    if trajectory is not None:
        points.append(trajectory)
    stacked = np.vstack(points)

    centre = 0.5 * (stacked.max(axis=0) + stacked.min(axis=0))
    half = 0.5 * float(np.max(stacked.max(axis=0) - stacked.min(axis=0)))
    half = max(half, 0.05) * 1.22
    return (
        (float(centre[0] - half), float(centre[0] + half)),
        (float(centre[1] - half), float(centre[1] + half)),
    )


def _panel_title(
    report: TiltPlanReport, solution: TiltSolution, view: StereogramView
) -> str:
    """Panel title: what the panel shows, and the move it describes."""

    target = format_direction_indices([int(value) for value in report.target.indices])
    if view == "detail":
        return (
            f"Detail: α {solution.position.alpha_deg:+.1f}°, "
            f"β {solution.position.beta_deg:+.1f}°"
        )
    return f"{report.current.phase.name}: tilt to {target}"


def _caption(report: TiltPlanReport, solution: TiltSolution) -> str:
    """Two lines stating what was found and what remains undetermined."""

    landed = (
        format_direction_indices([int(v) for v in solution.orbit_member_indices])
        if solution.orbit_member_indices is not None
        else "the requested direction"
    )
    travel = (
        f"{solution.path.total_travel_deg:.1f}° of crystal travel"
        if solution.path is not None
        else f"{solution.travel_deg:.1f}° of specimen rotation"
    )
    band = ""
    if solution.path is not None and solution.path.connecting_band is not None:
        indices = [int(v) for v in solution.path.connecting_band.indices]
        band = f"  Follow the {format_plane_indices(indices)} Kikuchi band."
    ambiguity = (
        "Orientation unambiguous."
        if report.ambiguity.is_unique
        else (
            f"{len(report.ambiguity.families)} orientation hypotheses; "
            f"this is family {solution.family.index}."
        )
    )
    return (
        f"Drive to α {solution.position.alpha_deg:+.2f}°, "
        f"β {solution.position.beta_deg:+.2f}° to place {landed} on the beam, "
        f"after {travel}; forward-validated residual "
        f"{solution.residual_deg:.3f}°.{band}\n"
        f"{report.reachable_orbit_size} of {report.orbit_size} symmetry-equivalent "
        f"targets are reachable with this holder.  {ambiguity}"
    )


def plot_tilt_stereogram(
    report: TiltPlanReport,
    stage: StageModel,
    *,
    solution: TiltSolution | None = None,
    view: StereogramView | Literal["both"] = "both",
    method: str = "stereographic",
    title: str | None = None,
    show_wulff_net: bool = True,
    show_reachable_region: bool = True,
    show_equivalents: bool = True,
    show_poles: bool = True,
    show_tilt_axes: bool = True,
    max_poles: int = 26,
    pole_families: Sequence[tuple[int, int, int]] = _DEFAULT_POLE_FAMILIES,
    theme: str = "journal",
    figsize: tuple[float, float] | None = None,
    ax: Any | None = None,
) -> Any:
    """Render the tilt-navigation stereogram: overview and detail, with a caption.

    Purpose
    -------
    The one figure that makes a tilt plan checkable by eye. It shows the current
    beam direction, the requested target and all its symmetry equivalents, the
    region the holder can reach, the trajectory as a series of dots, the
    principal zone axes and the sense of each tilt knob — all in the crystal's
    own stereographic frame, and all from the same numbers the engine produced.

    When to use
    -----------
    Immediately after `pytex.tem.plan_tilt_to_zone_axis`, before driving the
    stage. It is also the fastest way to spot a 180-degree diffraction-rotation
    error: the trajectory visibly heads away from the target rather than toward
    it.

    Parameters
    ----------
    report : TiltPlanReport
    stage : StageModel
    solution : TiltSolution, optional
        Defaults to the best-ranked solution.
    view : {"both", "overview", "detail"}, default "both"
        ``"both"`` draws the two panels side by side with a shared legend and
        caption, which is the publication form. The single-panel values are for
        embedding.
    method, title, show_* , max_poles, pole_families, theme
        As :func:`build_tilt_stereogram_figure_spec`.
    figsize : tuple of float, optional
        Overrides the default, which is sized for a two-column journal figure.
    ax : matplotlib axis, optional
        Draw a single panel onto an existing axis. Requires ``view`` to name one
        panel, since two panels cannot share one axis; the caller then owns the
        legend and caption.

    Returns
    -------
    matplotlib figure

    Examples
    --------
    Plan a move and draw it::

        report = plan_tilt_to_zone_axis(current, target, stage)
        figure = plot_tilt_stereogram(report, stage)

    See Also
    --------
    pytex.tem.navigation.plan_tilt_to_zone_axis : produces the report drawn here.
    build_tilt_stereogram_figure_spec : the declarative form, for testing.
    """

    import matplotlib.pyplot as plt

    def spec_for(panel: StereogramView) -> FigureSpec2D:
        return build_tilt_stereogram_figure_spec(
            report,
            stage,
            solution=solution,
            view=panel,
            method=method,
            title=title if view != "both" else None,
            show_wulff_net=show_wulff_net,
            show_reachable_region=show_reachable_region,
            show_equivalents=show_equivalents,
            show_poles=show_poles,
            show_tilt_axes=show_tilt_axes,
            max_poles=max_poles,
            pole_families=pole_families,
            theme=theme,
        )

    if ax is not None:
        if view == "both":
            raise ValueError(
                'Two panels cannot share one axis. Pass view="overview" or '
                'view="detail" when supplying ax, or omit ax to get the full figure.'
            )
        render_figure_spec(spec_for(view), ax=ax)
        return ax.figure

    chosen = report.best() if (solution is None and report.solutions) else solution
    if chosen is None:
        chosen = report.nearest_approach

    if view == "both":
        figure, axes = plt.subplots(1, 2, figsize=figsize or (13.0, 7.6))
        for panel, axis in zip(("overview", "detail"), axes, strict=True):
            render_figure_spec(spec_for(panel), ax=axis)  # type: ignore[arg-type]
        legend_axis = axes[0]
        _finish_layout(figure, axes, legend_axis, report, chosen, title)
        return figure

    figure, axis = plt.subplots(figsize=figsize or (7.4, 8.0))
    render_figure_spec(spec_for(view), ax=axis)
    _finish_layout(figure, (axis,), axis, report, chosen, title)
    return figure


def _finish_layout(
    figure: Any,
    axes: Sequence[Any],
    legend_axis: Any,
    report: TiltPlanReport,
    solution: TiltSolution | None,
    title: str | None,
) -> None:
    """Move the legend out of the projection and place the caption beneath it.

    The renderer's default ``loc="best"`` legend lands inside the disc and covers
    the trajectory, which is the one thing the figure exists to show.
    """

    handles: list[Any] = []
    labels: list[str] = []
    for axis in axes:
        existing = axis.get_legend()
        if existing is not None:
            existing.remove()
        for handle, label in zip(*axis.get_legend_handles_labels(), strict=True):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        figure.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.085),
            ncol=3 if len(axes) > 1 else 2,
            frameon=False,
            fontsize=9.5,
            handletextpad=0.5,
            columnspacing=1.6,
        )
    if solution is not None:
        figure.text(
            0.5,
            0.018,
            _caption(report, solution),
            ha="center",
            va="bottom",
            fontsize=9.5,
            color=TILT_STEREOGRAM_COLORS["muted"],
        )
    if title is not None and len(axes) > 1:
        figure.suptitle(title, fontsize=14.0, color=TILT_STEREOGRAM_COLORS["ink"])
    figure.tight_layout(rect=(0.0, 0.155, 1.0, 0.94 if title else 0.97))
