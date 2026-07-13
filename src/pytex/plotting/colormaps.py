"""Scientific colormap and palette foundation for PyTex plotting.

This module is the single source of color semantics for the plotting layer:

- ``pytex.texture``: sequential white-to-dark ramp for pole-figure / ODF
  intensity in multiples of random density (m.r.d.). White anchors zero so
  an untextured background reads as paper, following texture-literature
  convention; lightness decreases monotonically so magnitude is legible in
  grayscale print and to color-vision-deficient readers.
- ``pytex.misorientation``: single-hue sequential ramp (white to deep blue)
  for KAM / GROD / misorientation-magnitude maps.
- ``pytex.diverging``: blue-to-red diverging ramp with a neutral midpoint for
  signed quantities (e.g. strain, deviation-from-mean).
- Okabe-Ito categorical palette for identity coloring (phase maps, legends):
  the standard color-vision-deficiency-safe set, assigned in fixed order.

All colormaps are plain matplotlib ``LinearSegmentedColormap`` objects
registered under their PyTex names by `register_pytex_colormaps` (idempotent),
so every plotting entry point and end users can request them by name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from matplotlib.colors import Colormap


def _require_matplotlib_colors() -> Any:
    try:
        import matplotlib
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return matplotlib


@dataclass(frozen=True, slots=True)
class ColormapSpec:
    """Declarative colormap definition: name, ordered color stops, and kind."""

    name: str
    colors: tuple[str, ...]
    kind: str

    def __post_init__(self) -> None:
        if self.kind not in {"sequential", "diverging"}:
            raise ValueError("ColormapSpec.kind must be 'sequential' or 'diverging'.")
        if len(self.colors) < 2:
            raise ValueError("ColormapSpec needs at least two color stops.")


# Texture intensity (m.r.d.): white zero anchor, then a monotone-lightness
# yellow-orange-brown ramp (ColorBrewer YlOrBr anchors).
TEXTURE_INTENSITY_SPEC = ColormapSpec(
    name="pytex.texture",
    colors=(
        "#ffffff",
        "#fff7bc",
        "#fee391",
        "#fec44f",
        "#fe9929",
        "#ec7014",
        "#cc4c02",
        "#993404",
        "#662506",
    ),
    kind="sequential",
)

# Misorientation magnitude: single-hue white-to-deep-blue (ColorBrewer Blues).
MISORIENTATION_SPEC = ColormapSpec(
    name="pytex.misorientation",
    colors=(
        "#ffffff",
        "#deebf7",
        "#c6dbef",
        "#9ecae1",
        "#6baed6",
        "#4292c6",
        "#2171b5",
        "#08519c",
        "#08306b",
    ),
    kind="sequential",
)

# Signed quantities: cool-to-warm through a neutral near-white midpoint
# (ColorBrewer RdBu reversed to blue -> red reading order).
DIVERGING_SPEC = ColormapSpec(
    name="pytex.diverging",
    colors=(
        "#2166ac",
        "#4393c3",
        "#92c5de",
        "#d1e5f0",
        "#f7f7f7",
        "#fddbc7",
        "#f4a582",
        "#d6604d",
        "#b2182b",
    ),
    kind="diverging",
)

PYTEX_COLORMAP_SPECS: tuple[ColormapSpec, ...] = (
    TEXTURE_INTENSITY_SPEC,
    MISORIENTATION_SPEC,
    DIVERGING_SPEC,
)

# Okabe-Ito categorical palette (Wong 2011, Nature Methods 8:441): the standard
# CVD-safe identity set. Fixed assignment order; black last (outline-like, used
# only when seven hues are exhausted).
OKABE_ITO_COLORS: tuple[str, ...] = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
    "#F0E442",
    "#000000",
)


def categorical_colors(count: int) -> tuple[str, ...]:
    """Return ``count`` identity colors in the fixed Okabe-Ito order.

    Colors are assigned in a fixed order (never re-shuffled by count) so the
    same phase keeps the same color as maps gain or lose phases. Counts beyond
    the palette cycle with a repeat; prefer regrouping into fewer classes when
    that happens.
    """

    if count < 0:
        raise ValueError("count must be non-negative.")
    palette = OKABE_ITO_COLORS
    return tuple(palette[index % len(palette)] for index in range(count))


def _build_colormap(spec: ColormapSpec) -> Colormap:
    matplotlib = _require_matplotlib_colors()
    colormap: Colormap = matplotlib.colors.LinearSegmentedColormap.from_list(
        spec.name, list(spec.colors), N=256
    )
    return colormap


def register_pytex_colormaps() -> tuple[str, ...]:
    """Register all PyTex colormaps with matplotlib (idempotent).

    Returns the tuple of registered colormap names. Safe to call repeatedly;
    already-registered names are left untouched.
    """

    matplotlib = _require_matplotlib_colors()
    registered = []
    for spec in PYTEX_COLORMAP_SPECS:
        if spec.name not in matplotlib.colormaps:
            matplotlib.colormaps.register(_build_colormap(spec), name=spec.name)
        registered.append(spec.name)
    return tuple(registered)


def get_pytex_colormap(name: str) -> Colormap:
    """Return a registered PyTex colormap by name, registering all on demand."""

    names = register_pytex_colormaps()
    if name not in names:
        raise ValueError(f"Unknown PyTex colormap {name!r}; available: {names}.")
    matplotlib = _require_matplotlib_colors()
    colormap: Colormap = matplotlib.colormaps[name]
    return colormap


def srgb_to_lightness(colors: np.ndarray) -> np.ndarray:
    """CIELAB L* for an ``(n, 3)`` array of sRGB colors in [0, 1].

    Used to verify that sequential ramps are perceptually ordered (monotone
    lightness); exposed publicly so downstream palettes can run the same check.
    """

    rgb = np.asarray(colors, dtype=np.float64)
    if rgb.ndim != 2 or rgb.shape[1] != 3:
        raise ValueError("colors must have shape (n, 3).")
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    # Y (relative luminance) row of the sRGB -> XYZ (D65) matrix.
    luminance = linear @ np.array([0.2126729, 0.7151522, 0.0721750])
    epsilon = (6.0 / 29.0) ** 3
    kappa = (29.0 / 3.0) ** 3
    scaled = np.where(
        luminance > epsilon, np.cbrt(luminance) * 116.0 - 16.0, luminance * kappa
    )
    return np.ascontiguousarray(scaled)


__all__ = [
    "DIVERGING_SPEC",
    "MISORIENTATION_SPEC",
    "OKABE_ITO_COLORS",
    "PYTEX_COLORMAP_SPECS",
    "TEXTURE_INTENSITY_SPEC",
    "ColormapSpec",
    "categorical_colors",
    "get_pytex_colormap",
    "register_pytex_colormaps",
    "srgb_to_lightness",
]
