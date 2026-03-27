from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from pytex.core._arrays import as_float_array

if TYPE_CHECKING:
    from matplotlib.figure import Figure


@dataclass(frozen=True, slots=True)
class ScatterLayer2D:
    points: np.ndarray
    label: str | None = None
    colors: str | np.ndarray | None = None
    values: np.ndarray | None = None
    sizes: float | np.ndarray = 30.0
    cmap: str = "viridis"
    alpha: float = 0.9
    edgecolors: str | None = "black"
    linewidths: float = 0.4
    colorbar_label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "points", as_float_array(self.points, shape=(None, 2)))
        if self.values is not None:
            object.__setattr__(
                self,
                "values",
                as_float_array(self.values, shape=(self.points.shape[0],)),
            )
        if isinstance(self.sizes, np.ndarray):
            object.__setattr__(
                self,
                "sizes",
                as_float_array(self.sizes, shape=(self.points.shape[0],)),
            )


@dataclass(frozen=True, slots=True)
class LineLayer2D:
    points: np.ndarray
    label: str | None = None
    color: str = "black"
    linewidth: float = 1.5
    linestyle: str = "-"
    closed: bool = False
    alpha: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "points", as_float_array(self.points, shape=(None, 2)))


@dataclass(frozen=True, slots=True)
class ImageLayer2D:
    image: np.ndarray
    extent: tuple[float, float, float, float]
    cmap: str = "viridis"
    alpha: float = 0.95
    colorbar_label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "image", as_float_array(self.image, shape=(None, None)))


@dataclass(frozen=True, slots=True)
class ContourLayer2D:
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    levels: int | np.ndarray = 12
    cmap: str = "viridis"
    alpha: float = 0.95
    colorbar_label: str | None = None
    line_color: str | None = "black"
    line_width: float = 0.45

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", as_float_array(self.x, shape=(None,)))
        object.__setattr__(self, "y", as_float_array(self.y, shape=(None,)))
        object.__setattr__(self, "values", as_float_array(self.values, shape=(None, None)))
        if self.values.shape != (self.y.shape[0], self.x.shape[0]):
            raise ValueError("ContourLayer2D.values must have shape (len(y), len(x)).")
        if isinstance(self.levels, np.ndarray):
            object.__setattr__(self, "levels", as_float_array(self.levels, shape=(None,)))


@dataclass(frozen=True, slots=True)
class TextLayer2D:
    position: np.ndarray
    text: str
    color: str = "black"
    fontsize: float = 10.0
    ha: str = "center"
    va: str = "center"

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", as_float_array(self.position, shape=(2,)))


@dataclass(frozen=True, slots=True)
class FigureSpec2D:
    title: str | None = None
    xlabel: str = "x"
    ylabel: str = "y"
    xlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    equal_aspect: bool = True
    grid: bool = True
    boundary_circle_radius: float | None = None
    scatter_layers: tuple[ScatterLayer2D, ...] = field(default_factory=tuple)
    line_layers: tuple[LineLayer2D, ...] = field(default_factory=tuple)
    image_layers: tuple[ImageLayer2D, ...] = field(default_factory=tuple)
    contour_layers: tuple[ContourLayer2D, ...] = field(default_factory=tuple)
    text_layers: tuple[TextLayer2D, ...] = field(default_factory=tuple)


@dataclass(frozen=True, slots=True)
class MultiFigureSpec2D:
    panels: tuple[FigureSpec2D, ...]
    ncols: int = 3
    figsize: tuple[float, float] | None = None
    suptitle: str | None = None

    def __post_init__(self) -> None:
        if not self.panels:
            raise ValueError("MultiFigureSpec2D requires at least one panel.")
        if self.ncols <= 0:
            raise ValueError("MultiFigureSpec2D.ncols must be strictly positive.")


@dataclass(frozen=True, slots=True)
class ScatterLayer3D:
    points: np.ndarray
    label: str | None = None
    colors: str | np.ndarray | None = None
    values: np.ndarray | None = None
    sizes: float | np.ndarray = 36.0
    cmap: str = "viridis"
    alpha: float = 0.9
    depthshade: bool = True
    colorbar_label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "points", as_float_array(self.points, shape=(None, 3)))
        if self.values is not None:
            object.__setattr__(
                self,
                "values",
                as_float_array(self.values, shape=(self.points.shape[0],)),
            )
        if isinstance(self.sizes, np.ndarray):
            object.__setattr__(
                self,
                "sizes",
                as_float_array(self.sizes, shape=(self.points.shape[0],)),
            )


@dataclass(frozen=True, slots=True)
class ArrowLayer3D:
    vectors: np.ndarray
    labels: tuple[str, ...] | None = None
    color: str = "tab:blue"
    linewidth: float = 1.4
    normalize: bool = False
    alpha: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "vectors", as_float_array(self.vectors, shape=(None, 3)))
        if self.labels is not None and len(self.labels) != self.vectors.shape[0]:
            raise ValueError("ArrowLayer3D.labels must match the number of vectors.")


@dataclass(frozen=True, slots=True)
class FigureSpec3D:
    title: str | None = None
    xlabel: str = "x"
    ylabel: str = "y"
    zlabel: str = "z"
    scatter_layers: tuple[ScatterLayer3D, ...] = field(default_factory=tuple)
    arrow_layers: tuple[ArrowLayer3D, ...] = field(default_factory=tuple)
    unit_sphere_radius: float | None = None
    grid: bool = True


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
    except ImportError as exc:  # pragma: no cover - environment-dependent branch
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, Figure


def _set_equal_aspect_3d(ax: Any, all_points: np.ndarray) -> None:
    mins = np.min(all_points, axis=0)
    maxs = np.max(all_points, axis=0)
    center = 0.5 * (mins + maxs)
    radius = 0.5 * float(np.max(maxs - mins))
    radius = max(radius, 1.0)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1.0, 1.0, 1.0))


def render_figure_spec(
    spec: FigureSpec2D | FigureSpec3D | MultiFigureSpec2D,
    *,
    ax: Any | None = None,
) -> Any:
    if isinstance(spec, MultiFigureSpec2D):
        return render_multi_figure_spec_2d(spec)
    if isinstance(spec, FigureSpec3D):
        return render_figure_spec_3d(spec, ax=ax)
    return render_figure_spec_2d(spec, ax=ax)


def render_figure_spec_2d(spec: FigureSpec2D, *, ax: Any | None = None) -> Any:
    plt, _ = _require_matplotlib()
    if ax is None:
        fig, axes = plt.subplots(figsize=(6.2, 6.2))
    else:
        axes = ax
        fig = axes.figure
    image_artist: Any | None = None
    for image_layer in spec.image_layers:
        image_artist = axes.imshow(
            image_layer.image,
            origin="lower",
            extent=image_layer.extent,
            cmap=image_layer.cmap,
            alpha=image_layer.alpha,
            aspect="auto",
        )
        if image_layer.colorbar_label is not None:
            colorbar = fig.colorbar(image_artist, ax=axes)
            colorbar.set_label(image_layer.colorbar_label)
    contour_artist: Any | None = None
    for contour_layer in spec.contour_layers:
        xx, yy = np.meshgrid(contour_layer.x, contour_layer.y)
        contour_artist = axes.contourf(
            xx,
            yy,
            contour_layer.values,
            levels=contour_layer.levels,
            cmap=contour_layer.cmap,
            alpha=contour_layer.alpha,
        )
        if contour_layer.line_color is not None:
            axes.contour(
                xx,
                yy,
                contour_layer.values,
                levels=contour_layer.levels,
                colors=contour_layer.line_color,
                linewidths=contour_layer.line_width,
                alpha=min(1.0, contour_layer.alpha + 0.05),
            )
        if contour_layer.colorbar_label is not None:
            colorbar = fig.colorbar(contour_artist, ax=axes)
            colorbar.set_label(contour_layer.colorbar_label)
    scatter_artist: Any | None = None
    for scatter_layer in spec.scatter_layers:
        color_input = (
            scatter_layer.colors if scatter_layer.values is None else scatter_layer.values
        )
        scatter_artist = axes.scatter(
            scatter_layer.points[:, 0],
            scatter_layer.points[:, 1],
            c=color_input,
            s=scatter_layer.sizes,
            cmap=scatter_layer.cmap if scatter_layer.values is not None else None,
            alpha=scatter_layer.alpha,
            edgecolors=scatter_layer.edgecolors,
            linewidths=scatter_layer.linewidths,
            label=scatter_layer.label,
        )
        if scatter_layer.values is not None and scatter_layer.colorbar_label is not None:
            colorbar = fig.colorbar(scatter_artist, ax=axes)
            colorbar.set_label(scatter_layer.colorbar_label)
    for line_layer in spec.line_layers:
        points = line_layer.points
        if line_layer.closed:
            points = np.vstack([points, points[0]])
        axes.plot(
            points[:, 0],
            points[:, 1],
            color=line_layer.color,
            linewidth=line_layer.linewidth,
            linestyle=line_layer.linestyle,
            alpha=line_layer.alpha,
            label=line_layer.label,
        )
    for text_layer in spec.text_layers:
        axes.text(
            float(text_layer.position[0]),
            float(text_layer.position[1]),
            text_layer.text,
            color=text_layer.color,
            fontsize=text_layer.fontsize,
            ha=text_layer.ha,
            va=text_layer.va,
        )
    if spec.boundary_circle_radius is not None:
        from matplotlib.patches import Circle

        axes.add_patch(
            Circle(
                (0.0, 0.0),
                radius=spec.boundary_circle_radius,
                fill=False,
                linewidth=1.2,
                linestyle="--",
                color="dimgray",
            )
        )
    axes.set_xlabel(spec.xlabel)
    axes.set_ylabel(spec.ylabel)
    if spec.title is not None:
        axes.set_title(spec.title)
    if spec.equal_aspect:
        axes.set_aspect("equal", adjustable="box")
    if spec.xlim is not None:
        axes.set_xlim(*spec.xlim)
    if spec.ylim is not None:
        axes.set_ylim(*spec.ylim)
    axes.grid(spec.grid, alpha=0.3)
    if any(layer.label is not None for layer in spec.scatter_layers) or any(
        layer.label is not None for layer in spec.line_layers
    ):
        axes.legend(loc="best")
    fig.tight_layout()
    return fig


def render_multi_figure_spec_2d(spec: MultiFigureSpec2D) -> Any:
    plt, _ = _require_matplotlib()
    n_panels = len(spec.panels)
    ncols = min(spec.ncols, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    figsize = spec.figsize or (6.0 * ncols, 5.4 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)
    flat_axes = axes.ravel()
    for axis, panel in zip(flat_axes, spec.panels, strict=False):
        render_figure_spec_2d(panel, ax=axis)
    for axis in flat_axes[n_panels:]:
        axis.set_visible(False)
    if spec.suptitle is not None:
        fig.suptitle(spec.suptitle, y=0.995)
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.975))
    else:
        fig.tight_layout()
    return fig


def render_figure_spec_3d(spec: FigureSpec3D, *, ax: Any | None = None) -> Any:
    plt, _ = _require_matplotlib()
    if ax is None:
        fig = plt.figure(figsize=(7.0, 6.2))
        axes = fig.add_subplot(111, projection="3d")
    else:
        axes = ax
        fig = axes.figure
    scatter_artist: Any | None = None
    all_points: list[np.ndarray] = []
    for scatter_layer in spec.scatter_layers:
        color_input = (
            scatter_layer.colors if scatter_layer.values is None else scatter_layer.values
        )
        scatter_artist = axes.scatter(
            scatter_layer.points[:, 0],
            scatter_layer.points[:, 1],
            scatter_layer.points[:, 2],
            c=color_input,
            s=scatter_layer.sizes,
            cmap=scatter_layer.cmap if scatter_layer.values is not None else None,
            alpha=scatter_layer.alpha,
            depthshade=scatter_layer.depthshade,
            label=scatter_layer.label,
        )
        all_points.append(scatter_layer.points)
        if scatter_layer.values is not None and scatter_layer.colorbar_label is not None:
            colorbar = fig.colorbar(scatter_artist, ax=axes, shrink=0.75)
            colorbar.set_label(scatter_layer.colorbar_label)
    for arrow_layer in spec.arrow_layers:
        vectors = arrow_layer.vectors
        all_points.append(vectors)
        axes.quiver(
            np.zeros(vectors.shape[0], dtype=np.float64),
            np.zeros(vectors.shape[0], dtype=np.float64),
            np.zeros(vectors.shape[0], dtype=np.float64),
            vectors[:, 0],
            vectors[:, 1],
            vectors[:, 2],
            normalize=arrow_layer.normalize,
            color=arrow_layer.color,
            linewidth=arrow_layer.linewidth,
            alpha=arrow_layer.alpha,
        )
        if arrow_layer.labels is not None:
            for vector, label in zip(vectors, arrow_layer.labels, strict=True):
                axes.text(
                    float(vector[0]),
                    float(vector[1]),
                    float(vector[2]),
                    label,
                    color=arrow_layer.color,
                )
    if spec.unit_sphere_radius is not None:
        u = np.linspace(0.0, 2.0 * np.pi, 48)
        v = np.linspace(0.0, np.pi, 24)
        x = spec.unit_sphere_radius * np.outer(np.cos(u), np.sin(v))
        y = spec.unit_sphere_radius * np.outer(np.sin(u), np.sin(v))
        z = spec.unit_sphere_radius * np.outer(np.ones_like(u), np.cos(v))
        axes.plot_wireframe(x, y, z, color="lightgray", linewidth=0.5, alpha=0.35)
        all_points.append(
            np.array(
                [
                    [-spec.unit_sphere_radius, 0.0, 0.0],
                    [spec.unit_sphere_radius, 0.0, 0.0],
                    [0.0, -spec.unit_sphere_radius, 0.0],
                    [0.0, spec.unit_sphere_radius, 0.0],
                    [0.0, 0.0, -spec.unit_sphere_radius],
                    [0.0, 0.0, spec.unit_sphere_radius],
                ],
                dtype=np.float64,
            )
        )
    axes.set_xlabel(spec.xlabel)
    axes.set_ylabel(spec.ylabel)
    axes.set_zlabel(spec.zlabel)
    if spec.title is not None:
        axes.set_title(spec.title)
    axes.grid(spec.grid, alpha=0.3)
    if all_points:
        _set_equal_aspect_3d(axes, np.vstack(all_points))
    if any(layer.label is not None for layer in spec.scatter_layers):
        axes.legend(loc="best")
    fig.tight_layout()
    return fig


def save_documentation_figure_svg(figure: Figure, path: str | Path) -> Path:
    destination = Path(path)
    if destination.suffix.lower() != ".svg":
        raise ValueError("Documentation figures must be saved with a .svg suffix.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, format="svg", bbox_inches="tight")
    return destination
