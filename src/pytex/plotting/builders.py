from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import normalize_vector
from pytex.core.batches import EulerSet, QuaternionSet, RotationSet, VectorSet
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.plotting._render import (
    ArrowLayer3D,
    ContourLayer2D,
    FigureSpec2D,
    FigureSpec3D,
    ImageLayer2D,
    LineLayer2D,
    MultiFigureSpec2D,
    ScatterLayer2D,
    ScatterLayer3D,
    TextLayer2D,
)
from pytex.texture.harmonics import HarmonicODF
from pytex.texture.models import ODF, InversePoleFigure, PoleFigure


def _frame_axis_labels(frame_axes: tuple[str, str, str]) -> tuple[str, str, str]:
    return (frame_axes[0], frame_axes[1], frame_axes[2])


def _projection_radius(method: str) -> float:
    if method == "equal_area":
        return float(np.sqrt(2.0))
    if method == "stereographic":
        return 1.0
    raise ValueError("Projection method must be 'equal_area' or 'stereographic'.")


def _projection_title(noun: str, method: str) -> str:
    return f"{noun} ({method.replace('_', ' ').title()})"


def _miller_label(indices: np.ndarray) -> str:
    return "(" + " ".join(str(int(value)) for value in indices) + ")"


def _operator_axis_and_angle(operator: np.ndarray) -> tuple[np.ndarray, float]:
    rotation = Rotation.from_matrix(operator)
    axis = rotation.axis
    nonzero = np.flatnonzero(np.abs(axis) > 1e-12)
    if nonzero.size > 0 and axis[int(nonzero[0])] < 0.0:
        axis = -axis
    return axis, rotation.angle_deg


def _angles_to_sizes(
    values: np.ndarray,
    *,
    min_size: float = 36.0,
    max_size: float = 180.0,
) -> np.ndarray:
    if values.size == 0:
        return np.empty(0, dtype=np.float64)
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if np.isclose(vmin, vmax):
        sizes = np.full(values.shape, 0.5 * (min_size + max_size), dtype=np.float64)
    else:
        sizes = min_size + (values - vmin) * (max_size - min_size) / (vmax - vmin)
    sizes.setflags(write=False)
    return sizes


def _weights_to_sizes(weights: np.ndarray) -> np.ndarray:
    normalized = weights / np.max(weights)
    sizes = 40.0 + 180.0 * normalized
    sizes = np.ascontiguousarray(sizes)
    sizes.setflags(write=False)
    return sizes


def _grid_centers(edges: np.ndarray) -> np.ndarray:
    centers = 0.5 * (edges[:-1] + edges[1:])
    centers = np.ascontiguousarray(centers, dtype=np.float64)
    centers.setflags(write=False)
    return centers


def _gaussian_kernel1d(sigma_bins: float) -> np.ndarray:
    radius = max(1, int(np.ceil(3.0 * sigma_bins)))
    offsets = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (offsets / max(sigma_bins, 1e-12)) ** 2)
    kernel /= np.sum(kernel)
    kernel = np.ascontiguousarray(kernel)
    kernel.setflags(write=False)
    return kernel


def _smooth_histogram(histogram: np.ndarray, *, sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0.0:
        smoothed = np.ascontiguousarray(histogram, dtype=np.float64)
        smoothed.setflags(write=False)
        return smoothed
    kernel = _gaussian_kernel1d(sigma_bins)
    padded_x = np.pad(histogram, ((kernel.size // 2, kernel.size // 2), (0, 0)), mode="edge")
    smooth_x = np.apply_along_axis(
        lambda column: np.convolve(column, kernel, mode="valid"),
        0,
        padded_x,
    )
    padded_y = np.pad(smooth_x, ((0, 0), (kernel.size // 2, kernel.size // 2)), mode="edge")
    smooth_xy = np.apply_along_axis(
        lambda row: np.convolve(row, kernel, mode="valid"),
        1,
        padded_y,
    )
    smooth_xy = np.ascontiguousarray(smooth_xy, dtype=np.float64)
    smooth_xy.setflags(write=False)
    return smooth_xy


def _masked_projection_density(
    histogram: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    *,
    radius: float,
    sigma_bins: float,
) -> np.ndarray:
    smoothed = _smooth_histogram(histogram, sigma_bins=sigma_bins)
    xcenters = _grid_centers(xedges)
    ycenters = _grid_centers(yedges)
    xx, yy = np.meshgrid(xcenters, ycenters, indexing="xy")
    mask = (xx * xx + yy * yy) <= radius * radius + 1e-12
    density = np.where(mask, smoothed.T, np.nan)
    density = np.ascontiguousarray(density, dtype=np.float64)
    density.setflags(write=False)
    return density


def _bunge_section_grid(
    odf: ODF | HarmonicODF,
    *,
    phi2_deg: float,
    phi1_steps: int,
    big_phi_steps: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    phi1 = np.linspace(0.0, 360.0, phi1_steps, dtype=np.float64)
    big_phi = np.linspace(0.0, 90.0, big_phi_steps, dtype=np.float64)
    phi1_mesh, big_phi_mesh = np.meshgrid(phi1, big_phi, indexing="xy")
    if isinstance(odf, HarmonicODF):
        angles_deg = np.column_stack(
            [
                phi1_mesh.reshape(-1),
                big_phi_mesh.reshape(-1),
                np.full(phi1_mesh.size, float(phi2_deg), dtype=np.float64),
            ]
        )
        section_orientations = OrientationSet.from_euler_angles(
            angles_deg,
            crystal_frame=odf.crystal_frame,
            specimen_frame=odf.specimen_frame,
            symmetry=odf.crystal_symmetry,
            phase=odf.phase,
            convention="bunge",
            degrees=True,
            provenance=odf.provenance,
        )
        values = np.asarray(odf.evaluate(section_orientations), dtype=np.float64).reshape(
            big_phi_steps,
            phi1_steps,
        )
        values = np.ascontiguousarray(values, dtype=np.float64)
        values.setflags(write=False)
        return phi1, big_phi, values
    support = odf.orientations.as_euler_set(convention="bunge", degrees=True).angles
    support_phi1 = support[:, 0][None, None, :]
    support_big_phi = support[:, 1][None, None, :]
    support_phi2 = support[:, 2][None, None, :]
    delta_phi1 = np.minimum(
        np.abs(phi1_mesh[:, :, None] - support_phi1),
        360.0 - np.abs(phi1_mesh[:, :, None] - support_phi1),
    )
    delta_big_phi = np.abs(big_phi_mesh[:, :, None] - support_big_phi)
    delta_phi2 = np.minimum(
        np.abs(float(phi2_deg) - support_phi2),
        360.0 - np.abs(float(phi2_deg) - support_phi2),
    )
    angular_distance_deg = np.sqrt(delta_phi1**2 + delta_big_phi**2 + delta_phi2**2)
    angular_distance_rad = np.deg2rad(angular_distance_deg)
    kernel_values = odf.kernel.evaluate(angular_distance_rad)
    values = np.tensordot(kernel_values, odf.normalized_weights, axes=([2], [0]))
    values = np.ascontiguousarray(values)
    values.setflags(write=False)
    return phi1, big_phi, values


def build_vector_figure_spec(
    vectors: VectorSet,
    *,
    normalize: bool = False,
    title: str | None = None,
) -> FigureSpec3D:
    plotted = vectors.normalized() if normalize else vectors
    xlabel, ylabel, zlabel = _frame_axis_labels(plotted.reference_frame.axes)
    labels = tuple(f"v{i}" for i in range(plotted.values.shape[0]))
    return FigureSpec3D(
        title=title or f"Vectors In {plotted.reference_frame.name}",
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel,
        arrow_layers=(ArrowLayer3D(vectors=plotted.values, labels=labels),),
        unit_sphere_radius=1.0 if normalize else None,
    )


def build_symmetry_orbit_figure_spec(
    symmetry: SymmetrySpec,
    seed_vector: ArrayLike | VectorSet,
    *,
    antipodal: bool = True,
    method: str = "equal_area",
    title: str | None = None,
) -> FigureSpec2D:
    if isinstance(seed_vector, VectorSet):
        if seed_vector.reference_frame != symmetry.reference_frame:
            raise ValueError("VectorSet.reference_frame must match SymmetrySpec.reference_frame.")
        if len(seed_vector) != 1:
            raise ValueError("Symmetry orbit plotting requires exactly one seed vector.")
        vector = seed_vector.values[0]
    else:
        vector = normalize_vector(seed_vector)
    orbit = symmetry.equivalent_vectors(vector)
    orbit_array = orbit.values if isinstance(orbit, VectorSet) else orbit
    reduced = symmetry.reduce_vector_to_fundamental_sector(vector, antipodal=antipodal)
    reduced_array = reduced.values[0] if isinstance(reduced, VectorSet) else reduced
    from pytex.texture.projections import project_directions

    orbit_points = project_directions(orbit_array, method=method, antipodal=antipodal)
    reduced_point = project_directions(reduced_array[None, :], method=method, antipodal=antipodal)
    sector_vertices = project_directions(
        symmetry.fundamental_sector(antipodal=antipodal).vertices,
        method=method,
        antipodal=antipodal,
    )
    radius = _projection_radius(method)
    return FigureSpec2D(
        title=title or _projection_title(f"{symmetry.point_group} Orbit", method),
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        boundary_circle_radius=radius,
        scatter_layers=(
            ScatterLayer2D(points=orbit_points, label="orbit", colors="tab:blue", sizes=60.0),
            ScatterLayer2D(points=reduced_point, label="reduced", colors="tab:red", sizes=90.0),
        ),
        line_layers=(LineLayer2D(points=sector_vertices, label="fundamental sector", closed=True),),
    )


def build_symmetry_elements_figure_spec(
    symmetry: SymmetrySpec,
    *,
    title: str | None = None,
) -> FigureSpec2D:
    from pytex.plotting.spherical import build_symmetry_elements_figure_spec as _build

    return _build(symmetry, title=title)


def build_euler_figure_spec(
    euler_set: EulerSet,
    *,
    values: np.ndarray | None = None,
    sizes: np.ndarray | None = None,
    title: str | None = None,
) -> FigureSpec2D:
    angles = euler_set.angles if euler_set.degrees else np.rad2deg(euler_set.angles)
    x = angles[:, 0]
    y = angles[:, 1]
    labels = {
        "bunge": (r"$\phi_1$ (deg)", r"$\Phi$ (deg)", r"$\phi_2$ (deg)"),
        "matthies": (r"$\alpha$ (deg)", r"$\beta$ (deg)", r"$\gamma$ (deg)"),
        "abg": (r"$\alpha$ (deg)", r"$\beta$ (deg)", r"$\gamma$ (deg)"),
    }
    xlabel, ylabel, color_label = labels[euler_set.convention]
    if values is None:
        values = angles[:, 2]
    return FigureSpec2D(
        title=title or f"Euler Scatter ({euler_set.convention})",
        xlabel=xlabel,
        ylabel=ylabel,
        xlim=(0.0, 360.0),
        ylim=(0.0, 180.0 if euler_set.convention == "bunge" else 360.0),
        scatter_layers=(
            ScatterLayer2D(
                points=np.column_stack([x, y]),
                values=values,
                sizes=36.0 if sizes is None else sizes,
                colorbar_label=color_label,
            ),
        ),
    )


def build_quaternion_figure_spec(
    quaternions: QuaternionSet,
    *,
    title: str | None = None,
) -> FigureSpec3D:
    vector_part = quaternions.quaternions[:, 1:]
    scalar_part = quaternions.quaternions[:, 0]
    return FigureSpec3D(
        title=title or "Quaternion Scatter",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        scatter_layers=(
            ScatterLayer3D(
                points=vector_part,
                values=scalar_part,
                colorbar_label="w",
                sizes=60.0,
                label="quaternions",
            ),
        ),
        unit_sphere_radius=1.0,
    )


def build_rotation_figure_spec(
    rotations: RotationSet,
    *,
    title: str | None = None,
) -> FigureSpec3D:
    axes = np.stack([Rotation(quaternion).axis for quaternion in rotations.quaternions], axis=0)
    angles = np.asarray(
        [Rotation(quaternion).angle_deg for quaternion in rotations.quaternions],
        dtype=np.float64,
    )
    return FigureSpec3D(
        title=title or "Rotation Axes",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        scatter_layers=(
            ScatterLayer3D(
                points=axes,
                values=angles,
                sizes=_angles_to_sizes(angles),
                colorbar_label="rotation angle (deg)",
                label="rotations",
            ),
        ),
        unit_sphere_radius=1.0,
    )


def build_orientation_figure_spec(
    orientations: OrientationSet,
    *,
    title: str | None = None,
) -> FigureSpec3D:
    rotation_set = orientations.as_rotation_set()
    xlabel, ylabel, zlabel = _frame_axis_labels(orientations.crystal_frame.axes)
    spec = build_rotation_figure_spec(rotation_set, title=title or "Orientation Axes")
    return FigureSpec3D(
        title=spec.title,
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel,
        scatter_layers=spec.scatter_layers,
        arrow_layers=spec.arrow_layers,
        unit_sphere_radius=spec.unit_sphere_radius,
    )


def build_pole_figure_spec(
    pole_figure: PoleFigure,
    *,
    method: str = "equal_area",
    kind: str = "scatter",
    bins: int = 72,
    sigma_bins: float = 1.25,
    levels: int = 12,
    title: str | None = None,
) -> FigureSpec2D:
    projected = pole_figure.project(method=method)
    radius = _projection_radius(method)
    if kind == "scatter":
        sizes = _weights_to_sizes(pole_figure.intensities)
        return FigureSpec2D(
            title=title or f"Pole Figure {_miller_label(pole_figure.pole.miller.indices)}",
            xlabel="projection x",
            ylabel="projection y",
            xlim=(-radius, radius),
            ylim=(-radius, radius),
            boundary_circle_radius=radius,
            scatter_layers=(
                ScatterLayer2D(
                    points=projected,
                    values=pole_figure.intensities,
                    sizes=sizes,
                    colorbar_label="intensity",
                    label="poles",
                ),
            ),
        )
    if kind == "histogram":
        histogram, xedges, yedges = pole_figure.histogram(bins=bins, method=method)
        return FigureSpec2D(
            title=title
            or f"Pole Figure Histogram {_miller_label(pole_figure.pole.miller.indices)}",
            xlabel="projection x",
            ylabel="projection y",
            xlim=(-radius, radius),
            ylim=(-radius, radius),
            boundary_circle_radius=radius,
            image_layers=(
                ImageLayer2D(
                    image=histogram.T,
                    extent=(xedges[0], xedges[-1], yedges[0], yedges[-1]),
                    colorbar_label="intensity",
                ),
            ),
        )
    if kind == "contour":
        histogram, xedges, yedges = pole_figure.histogram(bins=bins, method=method)
        return FigureSpec2D(
            title=title or f"Pole Figure Contours {_miller_label(pole_figure.pole.miller.indices)}",
            xlabel="projection x",
            ylabel="projection y",
            xlim=(-radius, radius),
            ylim=(-radius, radius),
            boundary_circle_radius=radius,
            contour_layers=(
                ContourLayer2D(
                    x=_grid_centers(xedges),
                    y=_grid_centers(yedges),
                    values=_masked_projection_density(
                        histogram,
                        xedges,
                        yedges,
                        radius=radius,
                        sigma_bins=sigma_bins,
                    ),
                    levels=levels,
                    colorbar_label="intensity",
                ),
            ),
        )
    raise ValueError("Pole figure kind must be 'scatter', 'histogram', or 'contour'.")


def build_inverse_pole_figure_spec(
    inverse_pole_figure: InversePoleFigure,
    *,
    method: str = "equal_area",
    title: str | None = None,
) -> FigureSpec2D:
    projected = inverse_pole_figure.project(method=method)
    radius = _projection_radius(method)
    line_layers: tuple[LineLayer2D, ...] = ()
    text_layers: tuple[TextLayer2D, ...] = ()
    sector_vertices = inverse_pole_figure.project_sector_vertices(method=method)
    if sector_vertices is not None:
        line_layers = (
            LineLayer2D(points=sector_vertices, label="fundamental sector", closed=True),
        )
        text_layers = tuple(
            TextLayer2D(position=vertex, text=axis)
            for vertex, axis in zip(
                sector_vertices,
                inverse_pole_figure.crystal_frame.axes,
                strict=True,
            )
        )
    return FigureSpec2D(
        title=title or "Inverse Pole Figure",
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        boundary_circle_radius=radius,
        scatter_layers=(
            ScatterLayer2D(
                points=projected,
                values=inverse_pole_figure.intensities,
                sizes=_weights_to_sizes(inverse_pole_figure.intensities),
                colorbar_label="intensity",
                label="directions",
            ),
        ),
        line_layers=line_layers,
        text_layers=text_layers,
    )


def build_odf_figure_spec(
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
) -> FigureSpec2D | MultiFigureSpec2D:
    if isinstance(odf, HarmonicODF):
        euler_set = odf.quadrature_orientations.as_euler_set(convention="bunge", degrees=True)
        weights = odf.quadrature_densities
    else:
        euler_set = odf.orientations.as_euler_set(convention="bunge", degrees=True)
        weights = odf.normalized_weights
    if kind == "scatter":
        return build_euler_figure_spec(
            euler_set,
            values=weights,
            sizes=(
                _angles_to_sizes(np.abs(weights), min_size=20.0, max_size=120.0)
                if isinstance(odf, HarmonicODF)
                else _weights_to_sizes(weights)
            ),
            title=title
            or (
                "Harmonic ODF Samples In Euler Space"
                if isinstance(odf, HarmonicODF)
                else "ODF Support In Euler Space"
            ),
        )
    if kind == "contour":
        angles = euler_set.angles if euler_set.degrees else np.rad2deg(euler_set.angles)
        histogram, xedges, yedges = np.histogram2d(
            angles[:, 0],
            angles[:, 1],
            bins=bins,
            range=[[0.0, 360.0], [0.0, 180.0]],
            weights=weights,
        )
        return FigureSpec2D(
            title=title or "ODF Density In Euler Space",
            xlabel=r"$\phi_1$ (deg)",
            ylabel=r"$\Phi$ (deg)",
            xlim=(0.0, 360.0),
            ylim=(0.0, 180.0),
            equal_aspect=False,
            contour_layers=(
                ContourLayer2D(
                    x=_grid_centers(xedges),
                    y=_grid_centers(yedges),
                    values=_smooth_histogram(histogram.T, sigma_bins=sigma_bins),
                    levels=levels,
                    colorbar_label="normalized density",
                ),
            ),
        )
    if kind == "sections":
        panels = []
        for phi2_deg in section_phi2_deg:
            phi1, big_phi, values = _bunge_section_grid(
                odf,
                phi2_deg=phi2_deg,
                phi1_steps=section_phi1_steps,
                big_phi_steps=section_big_phi_steps,
            )
            panels.append(
                FigureSpec2D(
                    title=rf"$\phi_2 = {phi2_deg:.0f}^\circ$",
                    xlabel=r"$\phi_1$ (deg)",
                    ylabel=r"$\Phi$ (deg)",
                    xlim=(0.0, 360.0),
                    ylim=(0.0, 90.0),
                    equal_aspect=False,
                    contour_layers=(
                        ContourLayer2D(
                            x=phi1,
                            y=big_phi,
                            values=values,
                            levels=levels,
                            colorbar_label="density",
                        ),
                    ),
                )
            )
        return MultiFigureSpec2D(
            panels=tuple(panels),
            ncols=3,
            suptitle=title or "ODF Bunge Sections",
        )
    raise ValueError("ODF kind must be 'scatter', 'contour', or 'sections'.")


def coerce_vector_set(
    vectors: VectorSet | ArrayLike,
    *,
    reference_frame: object | None = None,
) -> VectorSet:
    if isinstance(vectors, VectorSet):
        return vectors
    if reference_frame is None:
        raise ValueError("reference_frame is required when plotting raw vectors.")
    from pytex.core.frames import ReferenceFrame

    frame = reference_frame
    if not isinstance(frame, ReferenceFrame):
        raise TypeError("reference_frame must be a ReferenceFrame.")
    return VectorSet(values=np.asarray(vectors, dtype=np.float64), reference_frame=frame)


def coerce_rotation_set(
    rotations: Rotation | RotationSet | QuaternionSet | EulerSet,
) -> RotationSet:
    if isinstance(rotations, RotationSet):
        return rotations
    if isinstance(rotations, Rotation):
        return RotationSet.from_rotations([rotations])
    if isinstance(rotations, QuaternionSet):
        return rotations.to_rotation_set()
    return rotations.to_rotation_set()


def coerce_orientation_set(orientations: Orientation | OrientationSet) -> OrientationSet:
    if isinstance(orientations, OrientationSet):
        return orientations
    return OrientationSet.from_orientations([orientations])
