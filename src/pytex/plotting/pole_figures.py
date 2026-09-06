"""Publication-quality pole-figure drawings, and the settings behind them.

Purpose
-------
A pole figure in a paper is a contoured density on an equal-area disc, with
levels the reader can name, a scale shared with every figure it is compared
against, and enough identification on the drawing itself to survive being cut
out and pasted into a montage. This module produces that drawing, and exposes
every choice it makes as a declared setting rather than a hidden default.

Three objects carry the settings:

* :class:`ContourSpec` — which levels, on which scale, over which range. Level
  choice is where a pole figure most easily misleads: two figures contoured on
  their own maxima look alike however different they are, which is why a shared
  range is a first-class option rather than something a caller reconstructs.
* :class:`PoleFigureStyle` — how the drawing is composed: raster resolution,
  smoothing halfwidth, what is annotated, whether the axes and grid of a
  diagnostic plot are shown at all.
* :class:`PoleFigureSet` — several figures drawn together on one scale, which
  is the comparison case.

The density itself is evaluated on the sphere at the inverse projections of the
raster points, never binned in the drawing plane; see
:meth:`pytex.texture.PoleFigure.density_on_directions`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from itertools import pairwise

import numpy as np

from pytex.core.notation import format_plane_family_indices, format_plane_indices
from pytex.core.sphere import unproject_plane_points
from pytex.plotting._render import (
    ContourLayer2D,
    FigureSpec2D,
    MarkerLayer2D,
    MultiFigureSpec2D,
    TextLayer2D,
)
from pytex.plotting.colormaps import register_pytex_colormaps
from pytex.texture.models import PoleFigure, ResamplingEstimator

#: Level ladders a contoured density can use.
CONTOUR_SCALES = ("linear", "geometric", "explicit")

#: The density of an untextured specimen, in multiples of a random
#: distribution. A contour at this value separates orientations that are
#: over-represented from those that are under-represented, which is the only
#: level on a pole figure with an absolute meaning.
RANDOM_LEVEL_MRD = 1.0


@dataclass(frozen=True, slots=True)
class ContourSpec:
    """Which contour levels a density map is drawn at.

    Purpose
    -------
    The level set is the quantitative content of a contoured pole figure: it is
    what turns a coloured blob into a statement that some orientation is four
    times more common than random. Leaving it to the plotting library means the
    levels change whenever the data do, so no two figures can be compared and
    none can be quoted.

    When and where
    --------------
    Pass one to :func:`build_pole_figure_contour_spec` or to
    :class:`PoleFigureStyle`. Use :meth:`shared_across` when several figures
    must be read against each other — the whole point of a comparison plate.

    Attributes
    ----------
    scale : str
        ``"linear"`` (equal steps), ``"geometric"`` (equal ratios, which is how
        a sharp texture is usually shown because its dynamic range is large),
        or ``"explicit"`` (use :attr:`values` verbatim).
    count : int
        Number of intervals for the generated ladders; ignored when explicit.
    values : tuple of float, optional
        The levels themselves, for ``scale="explicit"``. Required then, ignored
        otherwise.
    vmin, vmax : float, optional
        Range the ladder spans. ``None`` takes the value from the data, which
        makes the figure self-scaled; setting them is what makes two figures
        comparable.
    filled : bool
        Draw filled bands (default). Turning this off leaves line contours
        alone, which reproduces well in a single-colour journal.
    lines : bool
        Draw the level lines over the bands (default).
    label_lines : bool
        Write the level value onto each line. Off by default because it
        crowds a small panel; on, it makes the drawing readable without a
        colour bar.
    cmap : str
        Colormap name for the filled bands.
    include_random_level : bool
        Insert a level at exactly 1 m.r.d. when it lies inside the range
        (default). It is the one level with an absolute meaning, and a ladder
        that steps over it hides the boundary between over- and
        under-representation.

    See Also
    --------
    PoleFigureStyle : the rest of the drawing's settings.
    """

    scale: str = "linear"
    count: int = 8
    values: tuple[float, ...] | None = None
    vmin: float | None = None
    vmax: float | None = None
    filled: bool = True
    lines: bool = True
    label_lines: bool = False
    cmap: str = "viridis"
    include_random_level: bool = True

    def __post_init__(self) -> None:
        if self.scale not in CONTOUR_SCALES:
            raise ValueError(f"ContourSpec.scale must be one of {CONTOUR_SCALES!r}.")
        if self.count < 2:
            raise ValueError("ContourSpec.count must be at least 2.")
        if self.scale == "explicit":
            if self.values is None or len(self.values) < 2:
                raise ValueError("ContourSpec.scale='explicit' requires at least two values.")
            ordered = tuple(float(value) for value in self.values)
            if any(later <= earlier for earlier, later in pairwise(ordered)):
                raise ValueError("ContourSpec.values must be strictly increasing.")
            object.__setattr__(self, "values", ordered)
        if self.vmin is not None and self.vmax is not None and not self.vmin < self.vmax:
            raise ValueError("ContourSpec requires vmin < vmax.")
        if not self.filled and not self.lines:
            raise ValueError("A contour spec that draws neither bands nor lines draws nothing.")

    def range_for(self, values: np.ndarray) -> tuple[float, float]:
        """The level range this spec implies for one data set.

        A declared bound always wins; an undeclared one is taken from the data.
        A field with no spread at all is given a nominal unit range so that a
        uniform pole figure still renders instead of failing.
        """

        finite = np.asarray(values, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            raise ValueError("Contour levels need at least one finite density value.")
        low = float(np.min(finite)) if self.vmin is None else float(self.vmin)
        high = float(np.max(finite)) if self.vmax is None else float(self.vmax)
        if not high > low:
            high = low + 1.0
        return low, high

    def levels_for(self, values: np.ndarray) -> np.ndarray:
        """The contour levels for one data set, as a strictly increasing array.

        Parameters
        ----------
        values : np.ndarray
            The density field, which may contain NaN outside the projection
            disc; those are ignored.

        Returns
        -------
        np.ndarray
            The level values, in m.r.d.

        Raises
        ------
        ValueError
            If a geometric ladder is asked for over a range that reaches zero
            or below. A ratio ladder has no meaning there, and silently
            shifting the range would misstate every level; declare ``vmin``
            instead.
        """

        if self.scale == "explicit":
            assert self.values is not None  # guaranteed by __post_init__
            return np.asarray(self.values, dtype=np.float64)
        low, high = self.range_for(values)
        if self.scale == "geometric":
            if low <= 0.0:
                raise ValueError(
                    "A geometric contour ladder needs a strictly positive lower bound: a density "
                    "of zero has no ratio to anything. Set ContourSpec.vmin to the smallest "
                    "density worth drawing, or use scale='linear'."
                )
            levels = np.geomspace(low, high, self.count + 1)
        else:
            levels = np.linspace(low, high, self.count + 1)
        if self.include_random_level and low < RANDOM_LEVEL_MRD < high:
            # Insert the random level in place of any neighbour it would crowd:
            # a ladder carrying both 1.00 and 1.07 draws a band too thin to see
            # and a colour bar with an unreadable pair of ticks.
            spacing = float(np.min(np.diff(levels))) if levels.size > 1 else 1.0
            keep = np.abs(levels - RANDOM_LEVEL_MRD) > 0.25 * spacing
            levels = np.union1d(
                levels[keep], np.array([RANDOM_LEVEL_MRD], dtype=np.float64)
            )
        return np.asarray(np.unique(levels), dtype=np.float64)

    def shared_across(self, value_sets: Sequence[np.ndarray]) -> ContourSpec:
        """This spec with its range widened to cover every data set given.

        Purpose
        -------
        Two pole figures contoured on their own maxima look equally strong
        however different they are; a plate of them is a picture of nothing. A
        shared range is what makes the comparison quantitative, and it is the
        default behaviour of :class:`PoleFigureSet`.

        Declared bounds are preserved: if a caller has already fixed ``vmax``,
        pooling the data must not silently override it.
        """

        lows: list[float] = []
        highs: list[float] = []
        for values in value_sets:
            low, high = self.range_for(values)
            lows.append(low)
            highs.append(high)
        if not lows:
            raise ValueError("A shared contour range needs at least one data set.")
        return replace(
            self,
            vmin=self.vmin if self.vmin is not None else min(lows),
            vmax=self.vmax if self.vmax is not None else max(highs),
        )


@dataclass(frozen=True, slots=True)
class PoleFigureStyle:
    """How a pole figure is drawn, apart from its contour levels.

    Purpose
    -------
    Separates the two kinds of setting a figure has: what the numbers are
    (:class:`ContourSpec`) and how the drawing presents them. The defaults are
    the publication ones — no Cartesian axes, no grid, a solid boundary, the
    specimen axes named on the rim and the extrema stated — because a figure
    that needs five keyword arguments before it can be published is a figure
    most people will publish unadjusted.

    Attributes
    ----------
    method : str
        ``"equal_area"`` (default) or ``"stereographic"``. Equal area is the
        correct choice for a density: it is the projection under which a random
        distribution is uniform, so equal areas on the drawing carry equal
        expected numbers of poles.
    resolution : int
        Raster points across the diameter. Cost is quadratic in it; 241 is
        smooth at full page width.
    halfwidth_deg : float
        Kernel halfwidth for the density estimate. The consequential physical
        choice — set it from the angular resolution of the measurement, not
        until the figure looks right.
    estimator : str, optional
        Override the estimator implied by the figure's ``sampling``.
    contours : ContourSpec
        The level settings.
    show_frame_axes : bool
        Name the specimen axes on the rim (default).
    frame_axis_labels : tuple of str, optional
        Override the names taken from the specimen frame.
    rotation_deg : float
        Rotate the drawing in its own plane. Zero (the default) draws the
        specimen frame's first axis to the right, which is what the projection
        actually does; 90 puts it at the top, which is the usual rolling-plane
        convention with RD up. The rim labels rotate with the data, so the
        drawing stays self-describing either way.
    sample_label : str, optional
        An identifier drawn on the figure itself. On a comparison plate this is
        what lets a reader attribute a panel without counting rows.
    show_pole_label : bool
        Draw the ``{hkl}`` family the figure belongs to (default).
    show_extrema : bool
        State the maximum and minimum density under the drawing (default). A
        contoured figure without them cannot be read quantitatively at all.
    mark_maximum : bool
        Put a marker at the strongest point of the density (default).
    show_colorbar : bool
        Attach a colour bar (default). Turned off for the panels of a shared
        plate, which carries one bar for the set.
    show_axes, show_grid : bool
        The diagnostic Cartesian frame. Off by default.
    title : str, optional
        Overrides the generated title.
    show_title : bool
        Draw a title at all (default).

    See Also
    --------
    build_pole_figure_contour_spec : the drawing this describes.
    PoleFigureSet : several of them on one scale.
    """

    method: str = "equal_area"
    resolution: int = 241
    halfwidth_deg: float = 10.0
    estimator: ResamplingEstimator | None = None
    contours: ContourSpec = field(default_factory=ContourSpec)
    show_frame_axes: bool = True
    frame_axis_labels: tuple[str, str] | None = None
    rotation_deg: float = 0.0
    sample_label: str | None = None
    show_pole_label: bool = True
    show_extrema: bool = True
    mark_maximum: bool = True
    show_colorbar: bool = True
    show_axes: bool = False
    show_grid: bool = False
    title: str | None = None
    show_title: bool = True

    def __post_init__(self) -> None:
        if self.method not in {"equal_area", "stereographic"}:
            raise ValueError("PoleFigureStyle.method must be 'equal_area' or 'stereographic'.")
        if self.resolution < 16:
            raise ValueError("PoleFigureStyle.resolution must be at least 16 points.")
        if not 0.0 < self.halfwidth_deg < 180.0:
            raise ValueError("PoleFigureStyle.halfwidth_deg must lie in (0, 180).")


def _projection_radius(method: str) -> float:
    return float(np.sqrt(2.0)) if method == "equal_area" else 1.0


def pole_figure_label(pole_figure: PoleFigure) -> str:
    """Bracket notation matching what the figure actually plots.

    A pole figure normally shows the whole symmetry-related orbit of its pole,
    which is written ``{hkl}``; only when family expansion was switched off does
    it show the single plane ``(hkl)``.
    """

    indices = tuple(int(value) for value in pole_figure.pole.miller.indices)
    if getattr(pole_figure, "includes_symmetry_family", True):
        return format_plane_family_indices(indices, style="plain")
    return format_plane_indices(indices, style="plain")


def pole_figure_density_grid(
    pole_figure: PoleFigure,
    *,
    style: PoleFigureStyle | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a pole figure's density on the raster of its drawing.

    Purpose
    -------
    The field a contour plot needs, computed the right way round: each raster
    point is inverse-projected to the direction it stands for, and the density
    is estimated there on the sphere. Nothing is binned into drawing pixels and
    nothing is smoothed in the projection plane, where the distortion is
    largest exactly where pole figures are most crowded.

    Parameters
    ----------
    pole_figure : PoleFigure
    style : PoleFigureStyle, optional

    Returns
    -------
    (x, y, density)
        Raster coordinates and the ``(len(y), len(x))`` density in m.r.d., with
        NaN outside the projection boundary so the contour stops at the rim.
    """

    settings = PoleFigureStyle() if style is None else style
    radius = _projection_radius(settings.method)
    axis = np.linspace(-radius, radius, settings.resolution, dtype=np.float64)
    grid_x, grid_y = np.meshgrid(axis, axis, indexing="xy")
    inside = (grid_x * grid_x + grid_y * grid_y) <= radius * radius + 1e-12
    points = np.column_stack([grid_x[inside], grid_y[inside]])
    if settings.rotation_deg:
        # Rotating the *drawing* means asking each raster point for the
        # direction that would land on it after rotation, which is the inverse
        # rotation applied to the point before unprojecting.
        angle = np.deg2rad(settings.rotation_deg)
        rotation = np.array(
            [[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]], dtype=np.float64
        )
        points = points @ rotation
    directions = unproject_plane_points(points, method=settings.method)
    values = pole_figure.density_on_directions(
        directions,
        halfwidth_deg=settings.halfwidth_deg,
        estimator=settings.estimator,
    )
    density = np.full(grid_x.shape, np.nan, dtype=np.float64)
    density[inside] = values
    density.setflags(write=False)
    return axis, axis, density


def _rim_label_layers(
    pole_figure: PoleFigure,
    settings: PoleFigureStyle,
    *,
    radius: float,
) -> tuple[TextLayer2D, ...]:
    """Name the specimen axes where they meet the rim."""

    if not settings.show_frame_axes:
        return ()
    axes = settings.frame_axis_labels
    if axes is None:
        frame_axes = pole_figure.specimen_frame.axes
        axes = (str(frame_axes[0]), str(frame_axes[1]))
    angle = np.deg2rad(settings.rotation_deg)
    offset = 1.09 * radius
    directions = {
        axes[0]: (np.cos(angle), np.sin(angle)),
        axes[1]: (-np.sin(angle), np.cos(angle)),
    }
    return tuple(
        TextLayer2D(
            position=np.array([offset * x, offset * y], dtype=np.float64),
            text=name,
            fontsize=10.0,
            ha="center",
            va="center",
        )
        for name, (x, y) in directions.items()
    )


def _annotation_layers(
    pole_figure: PoleFigure,
    settings: PoleFigureStyle,
    density: np.ndarray,
    *,
    radius: float,
) -> tuple[tuple[TextLayer2D, ...], tuple[MarkerLayer2D, ...]]:
    """Sample identity, pole family, and the extrema that make it quantitative."""

    texts: list[TextLayer2D] = []
    markers: list[MarkerLayer2D] = []
    finite = density[np.isfinite(density)]
    maximum = float(np.max(finite)) if finite.size else float("nan")
    minimum = float(np.min(finite)) if finite.size else float("nan")
    if settings.sample_label:
        texts.append(
            TextLayer2D(
                position=np.array([-1.05 * radius, 1.12 * radius], dtype=np.float64),
                text=settings.sample_label,
                fontsize=11.0,
                ha="left",
                va="center",
                bbox_facecolor="white",
                bbox_edgecolor="0.35",
                bbox_alpha=0.85,
            )
        )
    if settings.show_pole_label:
        texts.append(
            TextLayer2D(
                position=np.array([1.05 * radius, 1.12 * radius], dtype=np.float64),
                text=pole_figure_label(pole_figure),
                fontsize=11.0,
                ha="right",
                va="center",
            )
        )
    if settings.show_extrema and finite.size:
        texts.append(
            TextLayer2D(
                position=np.array([0.0, -1.16 * radius], dtype=np.float64),
                text=f"max {maximum:.2f} / min {minimum:.2f} m.r.d.",
                fontsize=9.0,
                ha="center",
                va="center",
            )
        )
    return tuple(texts), tuple(markers)


def _maximum_marker(
    x: np.ndarray, y: np.ndarray, density: np.ndarray
) -> tuple[MarkerLayer2D, ...]:
    if not np.any(np.isfinite(density)):
        return ()
    flat = int(np.nanargmax(density))
    row, column = np.unravel_index(flat, density.shape)
    return (
        MarkerLayer2D(
            points=np.array([[x[column], y[row]]], dtype=np.float64),
            marker="+",
            # An unfilled marker takes its colour from the face, and passing an
            # edge colour as well makes Matplotlib warn and ignore one of them.
            facecolors="black",
            edgecolors=None,
            sizes=70.0,
            linewidths=1.1,
        ),
    )


def build_pole_figure_contour_spec(
    pole_figure: PoleFigure,
    *,
    style: PoleFigureStyle | None = None,
    levels: np.ndarray | None = None,
) -> FigureSpec2D:
    """A publication-quality contoured pole figure.

    Purpose
    -------
    The drawing a paper needs, from a measured or reconstructed pole figure and
    a declared set of settings: contoured density on the projection disc, the
    specimen axes named, the sample identified, and the extrema stated so the
    contours can be read as numbers.

    When and where
    --------------
    Use it for any figure that will be shown to somebody. The scatter and
    histogram builders in :mod:`pytex.plotting.builders` remain the diagnostic
    views — they show where the data actually are, which a contour hides.

    Parameters
    ----------
    pole_figure : PoleFigure
    style : PoleFigureStyle, optional
        Defaults to the publication settings.
    levels : np.ndarray, optional
        Pre-resolved levels, which is how a comparison plate gives every panel
        the same ladder. Normally left to the style's contour spec.

    Returns
    -------
    FigureSpec2D

    See Also
    --------
    PoleFigureSet.build : the multi-sample form.
    """

    settings = PoleFigureStyle() if style is None else style
    register_pytex_colormaps()
    radius = _projection_radius(settings.method)
    x, y, density = pole_figure_density_grid(pole_figure, style=settings)
    resolved = settings.contours.levels_for(density) if levels is None else np.asarray(levels)
    texts, markers = _annotation_layers(pole_figure, settings, density, radius=radius)
    texts = texts + _rim_label_layers(pole_figure, settings, radius=radius)
    if settings.mark_maximum:
        markers = markers + _maximum_marker(x, y, density)
    title: str | None = settings.title
    if title is None and settings.show_title:
        title = f"Pole Figure {pole_figure_label(pole_figure)}"
    if not settings.show_title:
        title = None
    return FigureSpec2D(
        title=title,
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        show_axes=settings.show_axes,
        grid=settings.show_grid,
        boundary_circle_radius=radius,
        boundary_circle_color="black",
        boundary_circle_linestyle="-",
        boundary_circle_linewidth=1.1,
        contour_layers=(
            ContourLayer2D(
                x=x,
                y=y,
                values=density,
                levels=resolved,
                cmap=settings.contours.cmap,
                filled=settings.contours.filled,
                label_lines=settings.contours.label_lines,
                line_color="black" if settings.contours.lines else None,
                colorbar_label="m.r.d." if settings.show_colorbar else None,
            ),
        ),
        marker_layers=markers,
        text_layers=texts,
    )


@dataclass(frozen=True, slots=True)
class PoleFigureSet:
    """Several pole figures drawn together, on one scale, each identified.

    Purpose
    -------
    The comparison plate: rows of samples against columns of poles, every panel
    contoured at the same levels so that a difference on the page is a
    difference in the specimens. Contouring each panel on its own maximum —
    which is what happens when nothing says otherwise — makes every sample look
    equally textured, and is the single most common way a pole-figure plate
    misleads.

    When and where
    --------------
    Whenever more than one measurement is being shown at once: several samples,
    several poles of one sample, or a measured figure beside its
    reconstruction. For one figure alone, use
    :func:`build_pole_figure_contour_spec`.

    Attributes
    ----------
    figures : tuple of tuple of PoleFigure
        One inner tuple per sample; each holds that sample's pole figures in
        the order they should appear across the row. Rows must have equal
        length, so that a column means one pole.
    sample_labels : tuple of str
        One identifier per row, drawn on every panel of that row.
    style : PoleFigureStyle
        Applied to every panel; its ``sample_label`` is overridden per row and
        its colour bar is suppressed in favour of the plate's own.
    shared_scale : bool
        Contour every panel at the same levels (default). Switching it off
        gives each panel its own ladder, which is occasionally what a reader
        wants and is never a comparison.

    See Also
    --------
    ContourSpec.shared_across : the range computation behind ``shared_scale``.
    """

    figures: tuple[tuple[PoleFigure, ...], ...]
    sample_labels: tuple[str, ...]
    style: PoleFigureStyle = field(default_factory=PoleFigureStyle)
    shared_scale: bool = True

    def __post_init__(self) -> None:
        if not self.figures:
            raise ValueError("A pole figure set needs at least one sample.")
        widths = {len(row) for row in self.figures}
        if len(widths) != 1 or 0 in widths:
            raise ValueError(
                "Every sample in a pole figure set must contribute the same, non-zero number of "
                "pole figures, so that a column of the plate means one pole."
            )
        if len(self.sample_labels) != len(self.figures):
            raise ValueError("A pole figure set needs exactly one label per sample.")

    @property
    def sample_count(self) -> int:
        """Number of samples, i.e. rows of the plate."""

        return len(self.figures)

    @property
    def pole_count(self) -> int:
        """Number of pole figures per sample, i.e. columns of the plate."""

        return len(self.figures[0])

    def shared_levels(self) -> np.ndarray | None:
        """The level ladder every panel is drawn at, or ``None`` if not shared.

        Computing it requires the density of every panel, so this evaluates the
        whole plate; :meth:`build` reuses that work rather than repeating it.
        """

        if not self.shared_scale:
            return None
        densities = [
            pole_figure_density_grid(figure, style=self.style)[2]
            for row in self.figures
            for figure in row
        ]
        return self.style.contours.shared_across(densities).levels_for(
            np.concatenate([density[np.isfinite(density)] for density in densities])
        )

    def build(self, *, suptitle: str | None = None) -> MultiFigureSpec2D:
        """The plate: one panel per pole figure, rows of samples.

        Parameters
        ----------
        suptitle : str, optional

        Returns
        -------
        MultiFigureSpec2D
            With ``shared_colorbar`` set when the scale is shared, so the plate
            carries one bar for every panel rather than one bar each — which
            would suggest, wrongly, that each panel had its own scale.
        """

        levels = self.shared_levels()
        panels: list[FigureSpec2D] = []
        for label, row in zip(self.sample_labels, self.figures, strict=True):
            for column, figure in enumerate(row):
                style = replace(
                    self.style,
                    sample_label=label,
                    show_colorbar=not self.shared_scale,
                    show_title=False,
                    # Only the first column needs the specimen axes named; the
                    # rest of the row shares them, and repeating them on every
                    # panel crowds a plate without adding anything.
                    show_frame_axes=self.style.show_frame_axes and column == 0,
                )
                panels.append(
                    build_pole_figure_contour_spec(figure, style=style, levels=levels)
                )
        return MultiFigureSpec2D(
            panels=tuple(panels),
            ncols=self.pole_count,
            # Discs are square and carry a caption line, so the generic panel
            # size would leave a plate mostly white space.
            figsize=(3.5 * self.pole_count + 1.4, 3.9 * self.sample_count),
            suptitle=suptitle,
            shared_colorbar_label="m.r.d." if self.shared_scale else None,
        )


__all__ = [
    "CONTOUR_SCALES",
    "RANDOM_LEVEL_MRD",
    "ContourSpec",
    "PoleFigureSet",
    "PoleFigureStyle",
    "build_pole_figure_contour_spec",
    "pole_figure_density_grid",
    "pole_figure_label",
]
