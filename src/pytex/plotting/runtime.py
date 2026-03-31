from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from numpy.typing import ArrayLike

from pytex.core.batches import EulerSet, QuaternionSet, RotationSet, VectorSet
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalDirection, CrystalPlane
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.plotting._render import render_figure_spec
from pytex.plotting.builders import (
    build_euler_figure_spec,
    build_inverse_pole_figure_spec,
    build_odf_figure_spec,
    build_orientation_figure_spec,
    build_pole_figure_spec,
    build_quaternion_figure_spec,
    build_rotation_figure_spec,
    build_symmetry_orbit_figure_spec,
    build_vector_figure_spec,
    coerce_orientation_set,
    coerce_rotation_set,
    coerce_vector_set,
)
from pytex.plotting.spherical import (
    plot_crystal_directions as _plot_crystal_directions,
)
from pytex.plotting.spherical import (
    plot_crystal_planes as _plot_crystal_planes,
)
from pytex.plotting.spherical import (
    plot_symmetry_elements as _plot_symmetry_elements,
)
from pytex.plotting.spherical import (
    plot_wulff_net as _plot_wulff_net,
)
from pytex.texture.models import ODF, InversePoleFigure, PoleFigure


def plot_vector_set(
    vectors: VectorSet | ArrayLike,
    *,
    reference_frame: ReferenceFrame | None = None,
    normalize: bool = False,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    vector_set = coerce_vector_set(vectors, reference_frame=reference_frame)
    return render_figure_spec(
        build_vector_figure_spec(vector_set, normalize=normalize, title=title),
        ax=ax,
    )


def plot_symmetry_orbit(
    symmetry: SymmetrySpec,
    seed_vector: ArrayLike | VectorSet,
    *,
    antipodal: bool = True,
    method: str = "equal_area",
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_symmetry_orbit_figure_spec(
            symmetry,
            seed_vector,
            antipodal=antipodal,
            method=method,
            title=title,
        ),
        ax=ax,
    )


def plot_symmetry_elements(
    symmetry: SymmetrySpec,
    *,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    annotate_axes: bool = False,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return _plot_symmetry_elements(
        symmetry,
        method=method,
        include_wulff_net=include_wulff_net,
        annotate_axes=annotate_axes,
        title=title,
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
        ax=ax,
    )


def plot_wulff_net(
    *,
    method: str = "stereographic",
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return _plot_wulff_net(
        method=method,
        title=title,
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
        ax=ax,
    )


def plot_crystal_directions(
    directions: CrystalDirection | Sequence[CrystalDirection],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return _plot_crystal_directions(
        directions,
        labels=labels,
        method=method,
        include_wulff_net=include_wulff_net,
        title=title,
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
        ax=ax,
    )


def plot_crystal_planes(
    planes: CrystalPlane | Sequence[CrystalPlane],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    render: str = "trace",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return _plot_crystal_planes(
        planes,
        labels=labels,
        method=method,
        render=render,
        include_wulff_net=include_wulff_net,
        title=title,
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
        ax=ax,
    )


def plot_euler_set(
    euler_set: EulerSet,
    *,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(build_euler_figure_spec(euler_set, title=title), ax=ax)


def plot_quaternion_set(
    quaternions: QuaternionSet,
    *,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(build_quaternion_figure_spec(quaternions, title=title), ax=ax)


def plot_rotations(
    rotations: Rotation | RotationSet | QuaternionSet | EulerSet,
    *,
    representation: str = "axis_angle",
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    rotation_set = coerce_rotation_set(rotations)
    if representation == "axis_angle":
        return render_figure_spec(build_rotation_figure_spec(rotation_set, title=title), ax=ax)
    if representation == "euler":
        return render_figure_spec(
            build_euler_figure_spec(
                rotation_set.as_euler_set(convention="bunge", degrees=True),
                title=title or "Rotations In Euler Space",
            ),
            ax=ax,
        )
    raise ValueError("representation must be 'axis_angle' or 'euler'.")


def plot_orientations(
    orientations: Orientation | OrientationSet,
    *,
    representation: str = "euler",
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    orientation_set = coerce_orientation_set(orientations)
    if representation == "axis_angle":
        return render_figure_spec(
            build_orientation_figure_spec(orientation_set, title=title),
            ax=ax,
        )
    if representation == "euler":
        return render_figure_spec(
            build_euler_figure_spec(
                orientation_set.as_euler_set(convention="bunge", degrees=True),
                title=title or "Orientations In Euler Space",
            ),
            ax=ax,
        )
    raise ValueError("representation must be 'axis_angle' or 'euler'.")


def plot_pole_figure(
    pole_figure: PoleFigure,
    *,
    method: str = "equal_area",
    kind: str = "scatter",
    bins: int = 72,
    sigma_bins: float = 1.25,
    levels: int = 12,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_pole_figure_spec(
            pole_figure,
            method=method,
            kind=kind,
            bins=bins,
            sigma_bins=sigma_bins,
            levels=levels,
            title=title,
        ),
        ax=ax,
    )


def plot_inverse_pole_figure(
    inverse_pole_figure: InversePoleFigure,
    *,
    method: str = "equal_area",
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_inverse_pole_figure_spec(inverse_pole_figure, method=method, title=title),
        ax=ax,
    )


def plot_odf(
    odf: ODF,
    *,
    kind: str = "scatter",
    bins: int = 72,
    sigma_bins: float = 1.25,
    levels: int = 12,
    section_phi2_deg: tuple[float, ...] = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
    section_phi1_steps: int = 181,
    section_big_phi_steps: int = 91,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_odf_figure_spec(
            odf,
            kind=kind,
            bins=bins,
            sigma_bins=sigma_bins,
            levels=levels,
            section_phi2_deg=section_phi2_deg,
            section_phi1_steps=section_phi1_steps,
            section_big_phi_steps=section_big_phi_steps,
            title=title,
        ),
        ax=ax,
    )
