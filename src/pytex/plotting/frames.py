"""Drawing reference frames: 3D triads, embeddable gizmos, and documentation SVG.

A reference frame is only as useful as a reader's ability to see it. This module
renders the *same* `pytex.core.frames.ReferenceFrame` three ways, from one
geometry computation, so a frame looks identical wherever it appears:

1. **Scene primitives** — `frame_triad` and `frame_triad_primitives` produce
   `pytex.plotting.primitives.AxisTriad3D` / `PrimitiveScene3D` objects that drop
   straight into the existing 3D crystal and world-scene renderers.
2. **Embeddable gizmo** — `add_frame_indicator` draws a small orientation
   indicator into a corner of *any* existing matplotlib axes: a SAED
   diffractogram, a pole figure, an IPF map, a crystal viewer panel. This is how
   a figure states which way ``RD``, ``ND``, or the detector ``u`` axis points
   without a separate legend.
3. **Standalone SVG** — `reference_frame_svg` and `frame_catalog_svg` emit
   style-guide-compliant SVG text with no matplotlib involved at all, for the
   documentation system and for `docs/figures/`.

Projection
----------

The 2D renderers use an orthographic projection specified the way matplotlib's
3D axes specify a view: ``elev_deg`` above the ``XY`` plane and ``azim_deg``
around it from ``+X`` toward ``+Y``. The default view (elevation 22 degrees,
azimuth 34 degrees) matches the default used by the 3D crystal renderer, so a
gizmo and the scene it annotates agree.

Colors
------

Axis colors come from the repository triad palette already fixed in
`pytex.plotting.primitives` (blue, green, vermillion — Okabe-Ito derived, so
axis identity survives grayscale printing and common color-vision deficiencies).

See also
--------
`pytex.core.frames` : the frame model being drawn.
`pytex.core.frame_catalog` : the standard frames.
`docs/standards/visualization_style_guide.md` : the canonical SVG tokens.
"""

from __future__ import annotations

import zlib
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any
from xml.sax.saxutils import escape

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.frames import FrameTransform, ReferenceFrame
from pytex.plotting.primitives import (
    TRIAD_AXIS_COLORS,
    AxisTriad3D,
    Label3D,
    PrimitiveScene3D,
    reference_frame_triad,
)

__all__ = [
    "DEFAULT_VIEW_AZIM_DEG",
    "DEFAULT_VIEW_ELEV_DEG",
    "FrameTriad",
    "add_frame_indicator",
    "frame_catalog_svg",
    "frame_triad",
    "frame_triad_primitives",
    "plot_frame_relationship",
    "plot_reference_frame",
    "project_orthographic",
    "reference_frame_svg",
]

#: Default viewing elevation in degrees, matching the 3D crystal renderer.
DEFAULT_VIEW_ELEV_DEG = 22.0

#: Default viewing azimuth in degrees, matching the 3D crystal renderer.
DEFAULT_VIEW_AZIM_DEG = 34.0

# Canonical documentation tokens (docs/standards/visualization_style_guide.md).
_INK = "#07122f"
_MUTED = "#40506f"
_PAPER = "#fbfdff"
_PANEL = "#ffffff"
_PANEL_STROKE = "#d7e0ef"
_SVG_FONT = "Arial, Helvetica, sans-serif"


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment-dependent branch
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt


# --------------------------------------------------------------------------- #
# Projection
# --------------------------------------------------------------------------- #


def project_orthographic(
    vectors: ArrayLike,
    *,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
) -> tuple[np.ndarray, np.ndarray]:
    """Project 3D vectors onto a 2D viewing plane, orthographically.

    What it does
        Builds the camera basis for the requested view and returns both the
        on-screen coordinates and the depth along the view direction, so callers
        can draw far-to-near and keep occlusion order right.

    When to use it
        Whenever frame geometry has to appear in a 2D figure — the corner gizmo,
        the documentation SVG, or any custom 2D annotation of 3D axes.

    Parameters
    ----------
    vectors:
        Array-like of shape ``(n, 3)`` or ``(3,)`` in canonical Cartesian
        components.
    elev_deg:
        Elevation of the camera above the ``XY`` plane, in degrees.
    azim_deg:
        Azimuth of the camera measured from ``+X`` toward ``+Y``, in degrees.

    Returns
    -------
    screen : numpy.ndarray
        ``(n, 2)`` screen coordinates, right-handed with ``+y`` up.
    depth : numpy.ndarray
        ``(n,)`` signed distance along the view direction; larger is nearer the
        camera.
    """

    array = np.asarray(vectors, dtype=np.float64)
    single = array.ndim == 1
    stacked = array.reshape(1, 3) if single else array
    if stacked.ndim != 2 or stacked.shape[1] != 3:
        raise ValueError("project_orthographic expects vectors of shape (n, 3) or (3,).")

    elevation = float(np.deg2rad(elev_deg))
    azimuth = float(np.deg2rad(azim_deg))
    view = np.array(
        [
            np.cos(elevation) * np.cos(azimuth),
            np.cos(elevation) * np.sin(azimuth),
            np.sin(elevation),
        ]
    )
    right = np.array([-np.sin(azimuth), np.cos(azimuth), 0.0])
    up = np.cross(view, right)

    screen = np.column_stack([stacked @ right, stacked @ up])
    depth = stacked @ view
    return np.ascontiguousarray(screen), np.ascontiguousarray(depth)


# --------------------------------------------------------------------------- #
# Triad model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class FrameTriad:
    """A drawable triad computed once from a `ReferenceFrame`.

    What it does
        Resolves a frame's axis geometry, labels, colors, origin, and display
        length into a single immutable object that every renderer in this module
        consumes. Computing it once is what keeps the 3D scene triad, the corner
        gizmo, and the documentation SVG in exact agreement.

    When to use it
        Directly, when a caller needs the resolved axis endpoints (for a custom
        renderer or a layout calculation). Most callers should use `frame_triad`,
        `add_frame_indicator`, or `reference_frame_svg`, which build one
        internally.

    Parameters
    ----------
    frame:
        The frame being drawn.
    origin:
        Where the triad is anchored, in canonical Cartesian components.
    length:
        Display length of each axis arrow.
    colors:
        Three color strings, one per axis.
    normalize:
        When ``True`` (the default) each axis is drawn at unit length before
        scaling, so an oblique or scaled frame still yields a legible gizmo. Set
        ``False`` to draw the frame's true axis vectors.
    basis:
        Optional ``(3, 3)`` override whose **columns** are the axis vectors to
        draw, replacing the frame's own geometry. Use it when the axes must be
        expressed in some other frame's coordinates — for instance drawing a
        source frame's axes inside its target frame via
        `pytex.core.frames.FrameTransform.source_axes_in_target`. The frame is
        still used for labels, name, and domain.
    """

    frame: ReferenceFrame
    origin: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    length: float = 1.0
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS
    normalize: bool = True
    basis: np.ndarray | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin", as_float_array(self.origin, shape=(3,)))
        if len(self.colors) != 3:
            raise ValueError("FrameTriad.colors must contain exactly three colors.")
        if not np.isfinite(self.length) or self.length <= 0.0:
            raise ValueError("FrameTriad.length must be finite and strictly positive.")
        if self.basis is not None:
            object.__setattr__(self, "basis", as_float_array(self.basis, shape=(3, 3)))

    @property
    def labels(self) -> tuple[str, str, str]:
        """The frame's three axis labels."""

        return self.frame.axes

    def axis_matrix(self) -> np.ndarray:
        """``(3, 3)`` matrix whose columns are the drawn axis vectors."""

        if self.basis is not None:
            matrix = np.asarray(self.basis, dtype=np.float64)
            if self.normalize:
                norms = np.linalg.norm(matrix, axis=0, keepdims=True)
                matrix = matrix / np.where(norms == 0.0, 1.0, norms)
        else:
            matrix = np.asarray(
                self.frame.unit_axis_matrix() if self.normalize else self.frame.basis_matrix,
                dtype=np.float64,
            )
        return as_float_array(float(self.length) * matrix, shape=(3, 3))

    def endpoints(self) -> np.ndarray:
        """``(3, 3)`` array whose row ``i`` is the tip of axis ``i``."""

        return as_float_array(self.origin[None, :] + self.axis_matrix().T, shape=(3, 3))

    def describe(self) -> str:
        """Prose summary of what a viewer of this triad is looking at."""

        tips = self.endpoints()
        parts = [
            f"{label} to [{tips[index][0]:.3f} {tips[index][1]:.3f} {tips[index][2]:.3f}] "
            f"in {self.colors[index]}"
            for index, label in enumerate(self.labels)
        ]
        scaling = (
            "axes drawn at equal display length"
            if self.normalize
            else "axes drawn at their true relative lengths"
        )
        return (
            f"Triad for reference frame '{self.frame.name}' ({self.frame.domain.value} domain) "
            f"anchored at [{self.origin[0]:.3f} {self.origin[1]:.3f} {self.origin[2]:.3f}], "
            f"{scaling}: {'; '.join(parts)}."
        )


def frame_triad(
    frame: ReferenceFrame,
    *,
    origin: ArrayLike = (0.0, 0.0, 0.0),
    length: float = 1.0,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    normalize: bool = True,
    **triad_kwargs: Any,
) -> AxisTriad3D:
    """Build an `AxisTriad3D` scene primitive for a reference frame.

    What it does
        Turns a frame into the same primitive the 3D renderers already draw, with
        the frame's axis labels as tip labels and its axis vectors as directions.

    When to use it
        To put a frame gizmo into a `pytex.plotting.primitives.PrimitiveScene3D`
        or a crystal/world scene — for instance to show the specimen ``RD/TD/ND``
        triad beside a rendered grain.

    Parameters
    ----------
    frame:
        The frame to draw.
    origin:
        Anchor point in the scene's world frame (angstrom for crystal scenes).
    length:
        Axis arrow length in the same units as ``origin``.
    colors:
        Per-axis colors; defaults to the repository triad palette.
    normalize:
        Draw unit-length axes (default) or the frame's true axis vectors.
    **triad_kwargs:
        Passed through to `AxisTriad3D` (``linewidth``, ``fontsize``, ...).

    Returns
    -------
    AxisTriad3D
        Ready to add to a scene's ``triads``.
    """

    triad = FrameTriad(
        frame=frame,
        origin=as_float_array(origin, shape=(3,)),
        length=length,
        colors=colors,
        normalize=normalize,
    )
    return reference_frame_triad(
        frame,
        basis=triad.axis_matrix(),
        origin=triad.origin,
        length=1.0,
        colors=colors,
        orthonormalize=False,
        **triad_kwargs,
    )


def frame_triad_primitives(
    frame: ReferenceFrame,
    *,
    origin: ArrayLike = (0.0, 0.0, 0.0),
    length: float = 1.0,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    normalize: bool = True,
    caption: str | None = None,
    **triad_kwargs: Any,
) -> PrimitiveScene3D:
    """Build a `PrimitiveScene3D` holding a frame's triad and an optional caption.

    Use this rather than `frame_triad` when the frame should be captioned in the
    scene (for example labelling two crystals in an orientation relationship);
    the returned scene merges into any other scene with
    `PrimitiveScene3D.merge`.

    Parameters
    ----------
    caption:
        Text anchored just below the triad origin. Pass ``None`` for no caption,
        or omit it and pass the frame name explicitly.
    """

    triad = frame_triad(
        frame,
        origin=origin,
        length=length,
        colors=colors,
        normalize=normalize,
        **triad_kwargs,
    )
    labels: tuple[Label3D, ...] = ()
    if caption is not None:
        anchor = as_float_array(origin, shape=(3,)) - np.array([0.0, 0.0, 0.22 * float(length)])
        labels = (Label3D(position=anchor, text=str(caption), color=_INK, fontsize=11.0),)
    return PrimitiveScene3D(triads=(triad,), labels=labels)


# --------------------------------------------------------------------------- #
# Matplotlib figures
# --------------------------------------------------------------------------- #


def plot_reference_frame(
    frame: ReferenceFrame,
    *,
    length: float = 1.0,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    normalize: bool = True,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    show_reference_box: bool = True,
    title: str | None = None,
    ax: Any | None = None,
) -> tuple[Any, Any]:
    """Draw one reference frame as a labelled 3D triad.

    What it does
        Renders the frame's three axes as colored arrows with their labels, on a
        3D axes, optionally inside a faint unit box that gives the eye a sense of
        the canonical ``X, Y, Z`` reference the axis vectors are quoted in.

    When to use it
        For teaching material and quick inspection — "where does this frame
        actually point?" — and as the panel builder behind multi-frame figures.

    Parameters
    ----------
    frame:
        The frame to draw.
    length:
        Axis arrow length.
    colors:
        Per-axis colors.
    normalize:
        Draw unit-length axes (default) or the frame's true axis vectors.
    elev_deg, azim_deg:
        Viewing angles, matching matplotlib's 3D convention.
    show_reference_box:
        Draw the faint canonical-Cartesian reference cube.
    title:
        Axes title; defaults to the frame name and domain.
    ax:
        An existing 3D axes to draw into. A new figure is created when omitted.

    Returns
    -------
    tuple
        ``(figure, axes)``. The caller owns the figure and is responsible for
        closing it.
    """

    plt = _require_matplotlib()
    if ax is None:
        figure = plt.figure(figsize=(4.2, 4.2))
        axes = figure.add_subplot(111, projection="3d")
    else:
        axes = ax
        figure = axes.figure

    triad = FrameTriad(
        frame=frame,
        length=length,
        colors=colors,
        normalize=normalize,
    )
    _draw_triad_3d(axes, triad)

    if show_reference_box:
        _draw_reference_box(axes, extent=float(length))

    span = 1.35 * float(length)
    axes.set_xlim(-span, span)
    axes.set_ylim(-span, span)
    axes.set_zlim(-span, span)
    axes.set_box_aspect((1.0, 1.0, 1.0))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    axes.set_axis_off()
    axes.set_title(
        f"{frame.name} ({frame.domain.value})" if title is None else title,
        fontsize=11.0,
        color=_INK,
    )
    return figure, axes


def plot_frame_relationship(
    transform: FrameTransform,
    *,
    length: float = 1.0,
    source_colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    target_colors: tuple[str, str, str] = ("#93b4f5", "#8fd6bd", "#f2a0a0"),
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    title: str | None = None,
    annotate: bool = True,
    ax: Any | None = None,
) -> tuple[Any, Any]:
    """Draw both frames of a `FrameTransform` in one 3D view.

    What it does
        Draws both frames **in the target frame's coordinates**, where the target
        frame is the identity triad (pale, dashed) and the source frame's axes
        are the columns of the transform's rotation matrix (saturated, solid).
        The rotation between them is therefore visible rather than merely stated.

    When to use it
        To document or debug a frame relationship: a vendor axis convention, a
        specimen mounted at an angle, a detector rotated about the beam. It is
        the visual companion to `FrameTransform.describe`.

    Parameters
    ----------
    transform:
        The relationship to draw.
    length:
        Axis arrow length.
    source_colors, target_colors:
        Per-axis colors for each triad. The defaults pair a saturated source
        triad with a pale target triad so the two stay distinguishable.
    elev_deg, azim_deg:
        Viewing angles.
    title:
        Axes title; defaults to naming both frames and the rotation angle.
    annotate:
        Add a caption stating the rotation angle and axis.
    ax:
        An existing 3D axes to draw into.

    Returns
    -------
    tuple
        ``(figure, axes)``. The caller owns the figure.
    """

    plt = _require_matplotlib()
    if ax is None:
        figure = plt.figure(figsize=(4.6, 4.6))
        axes = figure.add_subplot(111, projection="3d")
    else:
        axes = ax
        figure = axes.figure

    # Both triads are drawn in the *target* frame's coordinates: the target frame
    # is the identity triad there by definition, and the source frame's axes are
    # the columns of the transform's rotation matrix.
    target_triad = FrameTriad(
        frame=transform.target,
        length=length,
        colors=target_colors,
        basis=np.eye(3),
    )
    source_triad = FrameTriad(
        frame=transform.source,
        length=length,
        colors=source_colors,
        basis=transform.source_axes_in_target(),
    )

    _draw_triad_3d(axes, target_triad, linestyle="--", alpha=0.75)
    _draw_triad_3d(axes, source_triad)

    span = 1.35 * float(length)
    axes.set_xlim(-span, span)
    axes.set_ylim(-span, span)
    axes.set_zlim(-span, span)
    axes.set_box_aspect((1.0, 1.0, 1.0))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    axes.set_axis_off()

    if title is None:
        title = (
            f"{transform.source.name} -> {transform.target.name} "
            f"({transform.rotation_angle_deg:.1f} deg)"
        )
    axes.set_title(title, fontsize=11.0, color=_INK)
    if annotate:
        axis = transform.rotation_axis
        axes.text2D(
            0.02,
            0.02,
            f"solid: {transform.source.name} axes seen in {transform.target.name}\n"
            f"dashed: {transform.target.name} axes\n"
            f"rotation {transform.rotation_angle_deg:.2f} deg about "
            f"[{axis[0]:.3f} {axis[1]:.3f} {axis[2]:.3f}]",
            transform=axes.transAxes,
            fontsize=8.0,
            color=_MUTED,
            va="bottom",
        )
    return figure, axes


def add_frame_indicator(
    ax: Any,
    frame: ReferenceFrame,
    *,
    loc: str = "lower right",
    size: float = 0.18,
    pad: float = 0.02,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    axis_subset: Sequence[str | int] | None = None,
    basis: ArrayLike | None = None,
    fontsize: float = 8.0,
    linewidth: float = 1.6,
    label_frame: bool = False,
) -> Any:
    """Embed a small frame-orientation gizmo into the corner of an existing plot.

    What it does
        Adds an inset axes in the requested corner and draws the frame's axes,
        projected orthographically, as short labelled arrows. Axes pointing away
        from the viewer are drawn thinner and paler so the sense of depth reads
        correctly.

    When to use it
        On any 2D figure whose orientation would otherwise be ambiguous: a
        simulated SAED diffractogram (show the detector ``u/v`` axes or the zone
        axis), a pole figure (show ``RD``/``TD``), an IPF or KAM map (show the
        map axes), or a projected crystal-viewer panel. It is the standard way
        PyTex figures state their frame without a prose caption.

    Parameters
    ----------
    ax:
        The matplotlib axes to annotate. Any axes type works, including polar.
    frame:
        The frame to indicate.
    loc:
        One of ``"upper left"``, ``"upper right"``, ``"lower left"``,
        ``"lower right"``.
    size:
        Gizmo size as a fraction of the host axes.
    pad:
        Padding from the host axes edge, as a fraction.
    elev_deg, azim_deg:
        Viewing angles for the projection. Match these to the host figure's own
        view when annotating a projected 3D scene.
    colors:
        Per-axis colors.
    axis_subset:
        Draw only these axes (labels or indices), e.g. ``("RD", "TD")`` for an
        in-plane pole figure where ``ND`` points at the viewer and would
        degenerate to a point. Defaults to all three.
    basis:
        Optional ``(3, 3)`` override whose **columns** are the axis vectors to
        draw, replacing the frame's own geometry while keeping its labels. Use
        it when the axes must be expressed in the host figure's coordinates
        rather than the canonical Cartesian reference — for example drawing a
        crystal frame's axes on a detector, where the detector basis supplies
        the mapping.
    fontsize:
        Axis-label font size.
    linewidth:
        Arrow line width for axes pointing toward the viewer.
    label_frame:
        Add the frame name beneath the gizmo.

    Returns
    -------
    matplotlib.axes.Axes
        The inset axes, so callers can restyle it further.

    Raises
    ------
    ValueError
        If ``loc`` is not one of the four supported corners.
    """

    _require_matplotlib()
    corners = {
        "lower left": (pad, pad),
        "lower right": (1.0 - size - pad, pad),
        "upper left": (pad, 1.0 - size - pad),
        "upper right": (1.0 - size - pad, 1.0 - size - pad),
    }
    key = str(loc).strip().lower()
    if key not in corners:
        raise ValueError(
            f"loc must be one of {', '.join(sorted(corners))}; received '{loc}'."
        )
    left, bottom = corners[key]
    inset = ax.inset_axes((left, bottom, size, size))
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_facecolor("none")
    for spine in inset.spines.values():
        spine.set_visible(False)

    indices = (
        [frame.axis_index(item) for item in axis_subset]
        if axis_subset is not None
        else [0, 1, 2]
    )
    resolved_basis = None if basis is None else as_float_array(basis, shape=(3, 3))
    endpoints = FrameTriad(
        frame=frame, length=1.0, colors=colors, basis=resolved_basis
    ).endpoints()
    screen, depth = project_orthographic(endpoints, elev_deg=elev_deg, azim_deg=azim_deg)

    # Draw far axes first so nearer ones overlay them.
    for index in sorted(indices, key=lambda item: float(depth[item])):
        toward_viewer = float(depth[index]) >= 0.0
        tip = screen[index]
        inset.annotate(
            "",
            xy=(float(tip[0]), float(tip[1])),
            xytext=(0.0, 0.0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": colors[index],
                "linewidth": linewidth if toward_viewer else 0.6 * linewidth,
                "alpha": 1.0 if toward_viewer else 0.55,
                "shrinkA": 0.0,
                "shrinkB": 0.0,
            },
        )
        # Place the label a fixed distance beyond the tip along the axis's own
        # screen direction, so a strongly foreshortened axis (one pointing at the
        # viewer) still gets a label clear of the origin instead of on top of it.
        label_at = _label_anchor(tip)
        inset.text(
            float(label_at[0]),
            float(label_at[1]),
            frame.axes[index],
            color=colors[index],
            fontsize=fontsize,
            ha="center",
            va="center",
            alpha=1.0 if toward_viewer else 0.55,
        )

    if label_frame:
        inset.text(
            0.0,
            -1.85,
            frame.name,
            color=_MUTED,
            fontsize=0.85 * fontsize,
            ha="center",
            va="center",
        )
        inset.set_ylim(-2.05, 1.45)
    else:
        inset.set_ylim(-1.45, 1.45)
    inset.set_xlim(-1.45, 1.45)
    inset.set_aspect("equal")
    return inset


def _label_anchor(tip: ArrayLike, *, offset: float = 0.30, minimum: float = 0.45) -> np.ndarray:
    """Where to put an axis label given the axis tip's projected position.

    The label sits ``offset`` beyond the tip along the axis's screen direction,
    but never closer than ``minimum`` to the origin, so an axis pointing almost
    straight at the viewer (which projects to a near-zero-length arrow) still
    gets a legible, non-overlapping label. A degenerate zero-length projection
    falls back to placing the label straight up.
    """

    point = np.asarray(tip, dtype=np.float64).reshape(2)
    radius = float(np.hypot(point[0], point[1]))
    if radius < 1e-9:
        return np.array([0.0, minimum])
    return point * (max(radius + offset, minimum) / radius)


def _draw_triad_3d(
    axes: Any,
    triad: FrameTriad,
    *,
    linestyle: str = "-",
    alpha: float = 1.0,
) -> None:
    """Draw a `FrameTriad` onto a matplotlib 3D axes."""

    origin = triad.origin
    tips = triad.endpoints()
    for index, label in enumerate(triad.labels):
        vector = tips[index] - origin
        axes.quiver(
            origin[0],
            origin[1],
            origin[2],
            vector[0],
            vector[1],
            vector[2],
            color=triad.colors[index],
            linewidth=2.0,
            linestyle=linestyle,
            alpha=alpha,
            arrow_length_ratio=0.16,
        )
        tip = origin + 1.12 * vector
        axes.text(
            tip[0],
            tip[1],
            tip[2],
            label,
            color=triad.colors[index],
            fontsize=10.0,
            alpha=alpha,
            ha="center",
            va="center",
        )


def _draw_reference_box(axes: Any, *, extent: float) -> None:
    """Draw the faint canonical-Cartesian reference cube behind a triad."""

    corners = np.array(
        [
            [0.0, 0.0, 0.0],
            [extent, 0.0, 0.0],
            [extent, extent, 0.0],
            [0.0, extent, 0.0],
            [0.0, 0.0, extent],
            [extent, 0.0, extent],
            [extent, extent, extent],
            [0.0, extent, extent],
        ]
    )
    edges = (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    )
    for start, end in edges:
        segment = corners[[start, end]]
        axes.plot(
            segment[:, 0],
            segment[:, 1],
            segment[:, 2],
            color=_PANEL_STROKE,
            linewidth=0.7,
            zorder=0,
        )


# --------------------------------------------------------------------------- #
# Standalone SVG (no matplotlib)
# --------------------------------------------------------------------------- #


def _svg_triad_group(
    triad: FrameTriad,
    *,
    centre: tuple[float, float],
    scale: float,
    elev_deg: float,
    azim_deg: float,
    fontsize: float,
    marker_id: str,
    label_scale: float = 1.0,
) -> str:
    """Emit the SVG fragment for one projected triad.

    ``label_scale`` pushes the axis labels further from the origin, so a second
    triad drawn in the same panel can label its axes clear of the first one's.
    """

    screen, depth = project_orthographic(
        triad.endpoints(), elev_deg=elev_deg, azim_deg=azim_deg
    )
    parts: list[str] = []
    for index in sorted(range(3), key=lambda item: float(depth[item])):
        # SVG y grows downward, so the screen y axis is negated.
        x = centre[0] + scale * float(screen[index, 0])
        y = centre[1] - scale * float(screen[index, 1])
        toward_viewer = float(depth[index]) >= 0.0
        opacity = "1" if toward_viewer else "0.55"
        width = 3.4 if toward_viewer else 2.2
        colour = triad.colors[index]
        parts.append(
            f'    <line x1="{centre[0]:.2f}" y1="{centre[1]:.2f}" x2="{x:.2f}" y2="{y:.2f}" '
            f'stroke="{colour}" stroke-width="{width:.1f}" stroke-linecap="round" '
            f'opacity="{opacity}" marker-end="url(#{marker_id}-{index})"/>'
        )
        anchor = _label_anchor(screen[index]) * float(label_scale)
        label_x = centre[0] + scale * float(anchor[0])
        label_y = centre[1] - scale * float(anchor[1])
        parts.append(
            f'    <text x="{label_x:.2f}" y="{label_y:.2f}" font-size="{fontsize:.0f}" '
            f'font-family="{_SVG_FONT}" fill="{colour}" opacity="{opacity}" '
            f'text-anchor="middle" dominant-baseline="middle">'
            f"{escape(triad.labels[index])}</text>"
        )
    parts.append(
        f'    <circle cx="{centre[0]:.2f}" cy="{centre[1]:.2f}" r="3.4" fill="{_INK}"/>'
    )
    return "\n".join(parts)


def _svg_arrow_markers(marker_id: str, colors: Sequence[str]) -> str:
    """Emit one arrowhead marker per axis color.

    ``markerUnits="userSpaceOnUse"`` is essential: the default
    ``strokeWidth`` units would scale the head by the line width, so a triad
    drawn with a 3.4-unit stroke would sprout 30-unit arrowheads.
    """

    return "\n".join(
        f'    <marker id="{marker_id}-{index}" markerUnits="userSpaceOnUse" '
        f'markerWidth="13" markerHeight="10" refX="12" refY="5" orient="auto">\n'
        f'      <path d="M0,0 L13,5 L0,10 z" fill="{colour}"/>\n'
        f"    </marker>"
        for index, colour in enumerate(colors)
    )


def reference_frame_svg(
    frame: ReferenceFrame,
    *,
    width: float = 340.0,
    height: float = 320.0,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    title: str | None = None,
    subtitle: str | None = None,
    normalize: bool = True,
) -> str:
    """Render one reference frame as a standalone, self-describing SVG document.

    What it does
        Projects the frame's axes orthographically and writes a complete SVG
        document — arrowheads, labels, title, and an accessible ``<desc>`` — in
        pure Python. **No matplotlib is involved**, so documentation figures can
        be generated in a minimal environment.

    When to use it
        For documentation assets: Sphinx pages, `docs/figures/` SVGs, LaTeX
        includes, and anywhere a frame diagram should be a crisp vector figure
        rather than a rasterized plot. Use `plot_reference_frame` instead when
        the frame belongs inside a larger matplotlib figure.

    Parameters
    ----------
    frame:
        The frame to draw.
    width, height:
        Document size in SVG user units.
    elev_deg, azim_deg:
        Viewing angles for the orthographic projection.
    colors:
        Per-axis colors.
    title:
        Heading text and SVG ``<title>``; defaults to the frame name and domain.
    subtitle:
        Caption under the heading; defaults to the axis labels and their long
        names when the frame carries `ReferenceFrame.axis_descriptions`.
    normalize:
        Draw unit-length axes (default) or the frame's true axis vectors.

    Returns
    -------
    str
        A complete SVG document, ready to write to a ``.svg`` file.

    Notes
    -----
    Output follows `docs/standards/visualization_style_guide.md`: Arial-family
    text, the canonical ink/paper tokens, and mandatory ``<title>``/``<desc>``
    elements.
    """

    heading = title if title is not None else f"{frame.name} ({frame.domain.value} frame)"
    if subtitle is None:
        if frame.axis_descriptions:
            subtitle = "; ".join(
                f"{label} = {frame.axis_descriptions[index]}"
                for index, label in enumerate(frame.axes)
            )
        else:
            subtitle = f"axes {', '.join(frame.axes)}, {frame.handedness.value}-handed"

    triad = FrameTriad(frame=frame, length=1.0, colors=colors, normalize=normalize)
    # Pushed below the mid-line so the vertical axis label clears the subtitle.
    centre = (width / 2.0, height / 2.0 + 34.0)
    scale = 0.30 * min(width, height)
    # A *stable* digest, not builtins.hash: Python randomizes string hashing per
    # process, so a hash-derived id changed the marker names on every run and
    # regenerating a committed figure always produced a diff. A generated asset
    # whose bytes move for no reason cannot be checked for drift, which is the
    # whole reason these figures are generated rather than drawn.
    digest = zlib.crc32(frame.name.encode("utf-8")) % 100000
    marker_id = f"pytex-frame-arrow-{digest}"
    triad_svg = _svg_triad_group(
        triad,
        centre=centre,
        scale=scale,
        elev_deg=elev_deg,
        azim_deg=azim_deg,
        fontsize=16.0,
        marker_id=marker_id,
    )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img">\n'
        f"  <title>{escape(heading)}</title>\n"
        f"  <desc>{escape(frame.describe())}</desc>\n"
        f"  <defs>\n{_svg_arrow_markers(marker_id, colors)}\n  </defs>\n"
        f'  <rect width="{width:.0f}" height="{height:.0f}" fill="{_PAPER}"/>\n'
        f'  <text x="{width / 2.0:.0f}" y="34" font-size="18" font-family="{_SVG_FONT}" '
        f'fill="{_INK}" text-anchor="middle">{escape(heading)}</text>\n'
        f'  <text x="{width / 2.0:.0f}" y="56" font-size="12" font-family="{_SVG_FONT}" '
        f'fill="{_MUTED}" text-anchor="middle">{escape(subtitle)}</text>\n'
        f"  <g>\n"
        f"{triad_svg}\n"
        f"  </g>\n"
        f"</svg>\n"
    )


def frame_catalog_svg(
    frames: Sequence[ReferenceFrame],
    *,
    columns: int = 3,
    panel_width: float = 280.0,
    panel_height: float = 260.0,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    title: str = "PyTex Reference Frame Catalog",
    subtitle: str = (
        "Each panel shows one standard frame's axes as components in the canonical "
        "Cartesian reference."
    ),
) -> str:
    """Render several reference frames as one multi-panel documentation SVG.

    What it does
        Lays the supplied frames out on a grid of captioned panels, each showing
        the frame's projected triad, its domain, and its axis labels. Like
        `reference_frame_svg` this is pure Python with no matplotlib dependency.

    When to use it
        To generate the canonical catalog figure for the documentation
        (``docs/figures/reference_frame_catalog.svg``), or any teaching figure
        that has to contrast several frames side by side — crystal versus
        specimen versus map versus detector.

    Parameters
    ----------
    frames:
        The frames to lay out, in reading order.
    columns:
        Panels per row.
    panel_width, panel_height:
        Panel size in SVG user units.
    elev_deg, azim_deg:
        Viewing angles, shared by every panel so the panels are comparable.
    colors:
        Per-axis colors.
    title, subtitle:
        Figure heading and caption.

    Returns
    -------
    str
        A complete SVG document.

    Raises
    ------
    ValueError
        If no frames are supplied or ``columns`` is not positive.
    """

    if not frames:
        raise ValueError("frame_catalog_svg requires at least one frame.")
    if columns <= 0:
        raise ValueError("columns must be a positive integer.")

    count = len(frames)
    column_count = min(columns, count)
    row_count = int(np.ceil(count / column_count))
    margin = 26.0
    header = 92.0
    width = margin * 2.0 + column_count * panel_width
    height = header + margin + row_count * panel_height

    marker_id = "pytex-catalog-arrow"
    panels: list[str] = []
    for index, frame in enumerate(frames):
        row, column = divmod(index, column_count)
        left = margin + column * panel_width
        top = header + row * panel_height
        # Offset below the panel mid-line so the vertical axis label clears the caption.
        centre = (left + panel_width / 2.0, top + panel_height / 2.0 + 24.0)
        scale = 0.26 * min(panel_width, panel_height)
        triad = FrameTriad(frame=frame, length=1.0, colors=colors)
        triad_svg = _svg_triad_group(
            triad,
            centre=centre,
            scale=scale,
            elev_deg=elev_deg,
            azim_deg=azim_deg,
            fontsize=13.0,
            marker_id=marker_id,
        )
        panels.append(
            f'  <g>\n'
            f'    <rect x="{left + 6:.1f}" y="{top + 6:.1f}" '
            f'width="{panel_width - 12:.1f}" height="{panel_height - 12:.1f}" rx="8" '
            f'fill="{_PANEL}" stroke="{_PANEL_STROKE}" stroke-width="1"/>\n'
            f'    <text x="{centre[0]:.1f}" y="{top + 32:.1f}" font-size="15" '
            f'font-family="{_SVG_FONT}" fill="{_INK}" text-anchor="middle">'
            f"{escape(frame.name)}</text>\n"
            f'    <text x="{centre[0]:.1f}" y="{top + 50:.1f}" font-size="11" '
            f'font-family="{_SVG_FONT}" fill="{_MUTED}" text-anchor="middle">'
            f"{escape(frame.domain.value)} domain &#183; {escape('/'.join(frame.axes))}"
            f"</text>\n"
            f"{triad_svg}\n"
            f"  </g>"
        )

    description = "; ".join(
        f"{frame.name} ({frame.domain.value}): axes {', '.join(frame.axes)}"
        for frame in frames
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img">\n'
        f"  <title>{escape(title)}</title>\n"
        f"  <desc>{escape(subtitle)} Frames shown: {escape(description)}.</desc>\n"
        f"  <defs>\n{_svg_arrow_markers(marker_id, colors)}\n  </defs>\n"
        f'  <rect width="{width:.0f}" height="{height:.0f}" fill="{_PAPER}"/>\n'
        f'  <text x="{margin:.0f}" y="46" font-size="24" font-family="{_SVG_FONT}" '
        f'fill="{_INK}">{escape(title)}</text>\n'
        f'  <text x="{margin:.0f}" y="72" font-size="13" font-family="{_SVG_FONT}" '
        f'fill="{_MUTED}">{escape(subtitle)}</text>\n'
        + "\n".join(panels)
        + "\n</svg>\n"
    )
