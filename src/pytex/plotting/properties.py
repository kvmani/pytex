from __future__ import annotations

from typing import Any

from pytex.plotting.styles import resolve_style
from pytex.properties.tensors import (
    ComplianceTensor,
    DirectionalModulusSurface,
    StiffnessTensor,
    youngs_modulus_surface,
)


def _matplotlib() -> Any:
    import matplotlib.pyplot as plt

    return plt


def plot_youngs_modulus_surface(
    tensor: StiffnessTensor | ComplianceTensor | DirectionalModulusSurface,
    *,
    n_theta: int = 90,
    n_phi: int = 180,
    cmap: str = "viridis",
    unit: str | None = "auto",
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> Any:
    """Render the directional Young's-modulus surface as a 3D shape.

    ``tensor`` may be a `StiffnessTensor`/`ComplianceTensor` (sampled here) or a
    pre-computed `DirectionalModulusSurface`. The radius in each direction is the
    modulus along that direction; the colour also encodes the modulus, so elastic
    anisotropy reads directly off the shape. ``unit="auto"`` labels the colorbar
    with the conventional unit for the property (GPa for moduli, 1/GPa for
    compressibility, none for Poisson's ratio); pass a string or ``None`` to
    override.
    """

    plt = _matplotlib()
    if isinstance(tensor, DirectionalModulusSurface):
        surface = tensor
    else:
        surface = youngs_modulus_surface(tensor, n_theta=n_theta, n_phi=n_phi)
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    x, y, z = surface.cartesian_surface()
    values = surface.values
    normalized = (values - surface.minimum) / max(surface.maximum - surface.minimum, 1e-12)
    colormap = plt.get_cmap(cmap)
    fig = plt.figure(
        figsize=tuple(common["figure"]["figsize"]),
        dpi=int(common["figure"]["dpi"]),
        facecolor=common["figure"]["facecolor"],
    )
    axes = fig.add_subplot(111, projection="3d")
    axes.plot_surface(
        x,
        y,
        z,
        facecolors=colormap(normalized),
        rstride=1,
        cstride=1,
        linewidth=0,
        antialiased=True,
        shade=False,
    )
    limit = float(surface.maximum) * 1.05
    axes.set_xlim(-limit, limit)
    axes.set_ylim(-limit, limit)
    axes.set_zlim(-limit, limit)
    axes.set_box_aspect((1, 1, 1))
    display_names = {
        "youngs_modulus": "Young's modulus",
        "linear_compressibility": "linear compressibility",
        "shear_modulus_min": "min shear modulus",
        "shear_modulus_max": "max shear modulus",
        "poisson_ratio_min": "min Poisson's ratio",
        "poisson_ratio_max": "max Poisson's ratio",
    }
    default_units = {
        "youngs_modulus": "GPa",
        "linear_compressibility": "1/GPa",
        "shear_modulus_min": "GPa",
        "shear_modulus_max": "GPa",
    }
    label = display_names.get(surface.property_name, surface.property_name.replace("_", " "))
    resolved_unit = default_units.get(surface.property_name) if unit == "auto" else unit
    colorbar_label = f"{label} ({resolved_unit})" if resolved_unit else label
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_zlabel("z")
    axes.set_title(f"{label} surface  (anisotropy {surface.anisotropy_ratio:.2f}x)")
    mappable = plt.cm.ScalarMappable(cmap=colormap)
    mappable.set_array(values)
    colorbar = fig.colorbar(mappable, ax=axes, shrink=0.6, pad=0.1)
    colorbar.set_label(colorbar_label)
    fig.tight_layout()
    return fig


__all__ = ["plot_youngs_modulus_surface"]
