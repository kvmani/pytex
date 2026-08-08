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
    build_pole_figure_difference_spec,
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
    plot_stereographic_vectors as _plot_stereographic_vectors,
)
from pytex.plotting.spherical import (
    plot_symmetry_elements as _plot_symmetry_elements,
)
from pytex.plotting.spherical import (
    plot_wulff_net as _plot_wulff_net,
)
from pytex.texture.harmonics import HarmonicODF
from pytex.texture.models import (
    ODF,
    InversePoleFigure,
    PoleFigure,
    PoleFigureDifference,
)


def plot_vector_set(
    vectors: VectorSet | ArrayLike,
    *,
    reference_frame: ReferenceFrame | None = None,
    normalize: bool = False,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot a set of vectors on a spherical projection.

    Parameters
    ----------
    vectors : VectorSet or ArrayLike
        A typed vector set, or a raw ``(n, 3)`` array together with
        ``reference_frame``.
    reference_frame : ReferenceFrame, optional
        Required when raw arrays are passed, so the plot states which frame
        the directions live in.
    title : str, optional
    ax : matplotlib Axes, optional
        Draw into an existing axes instead of creating a figure.

    Returns
    -------
    Any
        The Matplotlib axes.
    """

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
    """Plot the symmetry orbit of a direction on a spherical projection.

    Purpose
    -------
    Show every direction the point group makes equivalent to the given one —
    the picture behind a ``<uvw>`` family, and the clearest way to see how
    many equivalents a symmetry produces.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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
    """Plot the symmetry elements of a point group stereographically.

    Purpose
    -------
    The classical stereogram: rotation axes and mirror traces of the group,
    drawn in the conventional crystallographic layout for teaching and for
    checking that a symmetry has been declared as intended.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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
    """Draw a Wulff (or Schmidt) net graticule.

    The reference net against which angles are read off a stereographic
    projection by hand. Useful as a background layer under a pole figure.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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
    """Plot crystal directions on a spherical projection.

    Accepts one direction or a sequence. Directions are resolved through the
    direct basis, so the plotted positions are correct in non-cubic lattices
    where index and direction do not coincide.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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
    """Plot crystal planes as poles on a spherical projection.

    Accepts one plane or a sequence. Planes are plotted as their normals,
    resolved through the reciprocal basis — the only correct route outside
    the cubic system.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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


def plot_stereographic_vectors(
    vectors: Any,
    *,
    labels: Sequence[str | None] | None = None,
    colors: Sequence[str] | None = None,
    method: str = "stereographic",
    render: str = "pole",
    antipodal: bool = True,
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot arbitrary vectors stereographically, with labels and colours.

    The general-purpose spherical plotting entry point that the more
    specific functions delegate to. Use it for overlays that mix sources —
    measured poles against computed ones, for example.

    Parameters
    ----------
    vectors : Any
        ``(n, 3)`` directions or a typed vector set.
    method : str
        Projection method; see the note below.
    antipodal : bool
        Fold onto one hemisphere (default), the usual pole-figure
        convention.
    title : str, optional
    ax : matplotlib Axes, optional

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

    return _plot_stereographic_vectors(
        vectors,
        labels=labels,
        colors=colors,
        method=method,
        render=render,
        antipodal=antipodal,
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
    """Plot Euler angles as a scatter in Euler space.

    Shows where a set of orientations sits in the ``(phi1, Phi, phi2)``
    coordinates the texture literature sections. Note that Euler space is not
    metrically faithful — equal angular differences do not correspond to
    equal distances, and the space is strongly distorted near ``Phi = 0``.

    Parameters
    ----------
    euler_set : EulerSet
        The typed set, so the convention is known rather than assumed.
    title : str, optional
    ax : matplotlib Axes, optional
    """

    return render_figure_spec(build_euler_figure_spec(euler_set, title=title), ax=ax)


def plot_quaternion_set(
    quaternions: QuaternionSet,
    *,
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot a quaternion set in orientation space.

    Parameters
    ----------
    quaternions : QuaternionSet
    title : str, optional
    ax : matplotlib Axes, optional
    """

    return render_figure_spec(build_quaternion_figure_spec(quaternions, title=title), ax=ax)


def plot_rotations(
    rotations: Rotation | RotationSet | QuaternionSet | EulerSet,
    *,
    representation: str = "axis_angle",
    title: str | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot rotations in axis-angle or Euler representation.

    Parameters
    ----------
    rotations : Rotation, RotationSet, QuaternionSet, or EulerSet
        Any rotation representation; coerced internally.
    representation : str
        ``"axis_angle"`` or ``"euler"``. Axis-angle shows the rotation axis
        distribution, which reveals fibres directly; Euler matches the
        literature's sectioning convention.
    title : str, optional
    ax : matplotlib Axes, optional

    Raises
    ------
    ValueError
        For an unrecognized representation.
    """

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
    """Plot crystal orientations in Euler or axis-angle representation.

    Parameters
    ----------
    orientations : Orientation or OrientationSet
    representation : str
        ``"euler"`` (default) or ``"axis_angle"``; see :func:`plot_rotations`.
    title : str, optional
    ax : matplotlib Axes, optional
    """

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


def plot_variant_pole_figure(
    variant_poles: Any,
    *,
    method: str = "stereographic",
    antipodal: bool = True,
    include_wulff_net: bool = True,
    label_variants: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot a ``VariantPoleFigure`` as a color-per-variant stereographic overlay.

    Purpose: renders the predicted specimen-frame child-plane poles of every
    transformation variant (from ``pytex.variant_pole_figure``) on a Wulff
    net, cycling one color per variant so the overlay can be read against a
    measured child pole figure. With ``label_variants`` the first pole of
    each variant carries a ``V<k>`` label.

    Output: the Matplotlib axes from ``plot_stereographic_vectors``.
    """

    indices = variant_poles.variant_indices
    palette = [
        "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b",
        "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#aec7e8", "#98df8a",
    ]
    unique = {int(index): position for position, index in enumerate(dict.fromkeys(indices))}
    colors = [palette[unique[int(index)] % len(palette)] for index in indices]
    labels: list[str | None] | None = None
    if label_variants:
        seen: set[int] = set()
        labels = []
        for index in indices:
            value = int(index)
            labels.append(f"V{value}" if value not in seen else None)
            seen.add(value)
    resolved_title = title or (
        f"{variant_poles.relationship_name}: variant poles"
    )
    return _plot_stereographic_vectors(
        variant_poles.poles.values,
        labels=labels,
        colors=colors,
        method=method,
        render="pole",
        antipodal=antipodal,
        include_wulff_net=include_wulff_net,
        title=resolved_title,
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
        ax=ax,
    )


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
    """Plot a pole figure.

    Parameters
    ----------
    pole_figure : PoleFigure
        The figure to render, from a measurement or reconstructed from an
        ODF.
    title : str, optional
    ax : matplotlib Axes, optional

    Returns
    -------
    Any
        The Matplotlib axes.

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

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


def plot_pole_figure_difference(
    difference: PoleFigureDifference,
    *,
    method: str = "equal_area",
    title: str | None = None,
    symmetric_limits: bool = True,
    ax: Any | None = None,
) -> Any:
    """Plot a signed pole-figure residual on a diverging colour scale.

    Purpose
    -------
    The QC figure for a PF-to-ODF inversion, and the standard way to compare
    two measurements. A residual norm says *how badly* two figures disagree;
    only this says *where*, and where is what identifies the cause — a
    systematic miss over one region of the specimen sphere points at an
    unmodelled component or an uncorrected defocusing loss, while noise spread
    evenly over the whole figure points at counting statistics.

    Parameters
    ----------
    difference : PoleFigureDifference
        From ``a.difference(b)`` or ``PoleFigureResidualReport.difference_figure()``.
    method : str
        Projection; see :func:`plot_pole_figure`.
    title : str, optional
        Defaults to a title naming the pole and both operands.
    symmetric_limits : bool
        Centre the colour scale on zero (default). Leave this on: the diverging
        colormap's neutral colour marks the zero crossing only when the limits
        are symmetric, so turning it off will colour part of a one-signed
        residual as though it had the other sign.
    ax : matplotlib Axes, optional

    Returns
    -------
    Any
        The Matplotlib axes.
    """

    return render_figure_spec(
        build_pole_figure_difference_spec(
            difference,
            method=method,
            title=title,
            symmetric_limits=symmetric_limits,
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
    """Plot an inverse pole figure in the standard triangle.

    The fundamental-sector outline is drawn when the figure carries crystal
    symmetry, so the plot reads as the conventional standard triangle rather
    than as an unbounded scatter.

    Parameters
    ----------
    inverse_pole_figure : InversePoleFigure
    title : str, optional
    ax : matplotlib Axes, optional

    The ``method`` argument selects the projection: ``"equal_area"``
    (Schmidt) preserves area, so densities are comparable across the figure,
    while ``"stereographic"`` (Wulff) preserves angles and is the right
    choice for angle-measuring constructions.
    """

    return render_figure_spec(
        build_inverse_pole_figure_spec(inverse_pole_figure, method=method, title=title),
        ax=ax,
    )


def plot_odf_phi2_sections(
    sections: Any,
    *,
    cmap: str = "pytex.texture",
    levels: int = 12,
    max_cols: int = 3,
    panel_labels: bool = True,
    colorbar_label: str = "ODF density (m.r.d.)",
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> Any:
    """Render a pre-computed `ODFSectionData` as a publication section panel.

    Complements ``plot_odf(kind="sections")`` by plotting an already-sampled
    `ODF.phi2_sections(...)` result, so the density grids can be computed once
    and then analysed or plotted without recomputation. Sections lay out on a
    `PanelGrid` (at most ``max_cols`` per row) under the requested theme, with
    Bunge-Euler math labels, optional (a), (b), (c) panel labels, a shared
    density scale across sections, and an m.r.d. colorbar. The default
    colormap is the white-anchored ``pytex.texture`` intensity ramp.
    """

    from pytex.plotting.colormaps import register_pytex_colormaps
    from pytex.plotting.figure import PanelGrid, label_panels

    register_pytex_colormaps()
    count = int(sections.section_count)
    cols = max(1, min(count, int(max_cols)))
    rows = (count + cols - 1) // cols
    grid = PanelGrid(
        rows,
        cols,
        panel_size=(3.2, 3.0),
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    vmax = sections.max_density
    contour = None
    section_symbol = (
        r"\sigma" if getattr(sections, "section_kind", "phi2") == "sigma" else r"\varphi_2"
    )
    for index in range(count):
        axis = grid.axes_flat[index]
        contour = axis.contourf(
            sections.phi1_deg,
            sections.big_phi_deg,
            sections.densities[index],
            levels=levels,
            cmap=cmap,
            vmin=0.0,
            vmax=vmax,
        )
        axis.set_title(rf"${section_symbol} = {sections.phi2_deg[index]:.0f}^\circ$")
        if index // cols == rows - 1 or index + cols >= count:
            axis.set_xlabel(r"$\varphi_1$ (deg)")
        if index % cols == 0:
            axis.set_ylabel(r"$\Phi$ (deg)")
        axis.set_aspect("equal", adjustable="box")
        axis.invert_yaxis()
    grid.hide_unused(count)
    if panel_labels and count > 1:
        label_panels(grid.axes_flat[:count])
    if contour is not None:
        grid.shared_colorbar(
            contour,
            axes=grid.axes_flat[:count],
            label=colorbar_label,
        )
    if title is not None:
        grid.figure.suptitle(title)
    return grid.figure


def plot_odf(
    odf: ODF | HarmonicODF,
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
    """Plot an orientation distribution function as Euler-space sections.

    Purpose
    -------
    The conventional presentation of an ODF: density contoured on
    constant-``phi2`` sections through Bunge Euler space, in multiples of a
    random distribution.

    Parameters
    ----------
    odf : ODF or HarmonicODF
        Either ODF representation.
    title : str, optional
    ax : matplotlib Axes, optional
    """

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
