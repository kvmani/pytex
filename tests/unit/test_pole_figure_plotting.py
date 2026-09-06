"""Contour control and multi-sample pole-figure comparison.

The tests are structural and semantic rather than pixel comparisons, per the
repository's plotting-validation rule: what matters is that the levels are the
ones asked for, that every panel of a comparison plate is drawn at the *same*
levels, and that the density being contoured is the density on the sphere.
"""

from __future__ import annotations

from dataclasses import replace

import matplotlib
import numpy as np
import pytest
from numpy.testing import assert_allclose

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.core.orientation import OrientationSet
from pytex.plotting.pole_figures import (
    RANDOM_LEVEL_MRD,
    ContourSpec,
    PoleFigureSet,
    PoleFigureStyle,
    build_pole_figure_contour_spec,
    pole_figure_density_grid,
)
from pytex.plotting.runtime import (
    plot_pole_figure_comparison,
    plot_pole_figure_contours,
)
from pytex.texture import ODF, KernelSpec, PoleFigure


def _context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("RD", "TD", "ND"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def _pole_figure(centre: tuple[float, float, float], spread: float, seed: int) -> PoleFigure:
    crystal, specimen, phase = _context()
    generator = np.random.default_rng(seed)
    angles = np.asarray(centre, dtype=np.float64) + spread * generator.standard_normal((250, 3))
    orientations = OrientationSet.from_euler_angles(
        angles,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF(
        orientations=orientations,
        weights=np.ones(len(orientations)),
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=12.0),
    )
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    return odf.reconstruct_pole_figure(pole)


def _uniform_pole_figure(count: int = 4000, seed: int = 3) -> PoleFigure:
    """A pole figure from an isotropic specimen: 1 m.r.d. everywhere."""

    crystal, specimen, phase = _context()
    generator = np.random.default_rng(seed)
    quaternions = generator.standard_normal((count, 4))
    quaternions /= np.linalg.norm(quaternions, axis=1, keepdims=True)
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF(
        orientations=orientations,
        weights=np.ones(count),
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=15.0),
    )
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    return odf.reconstruct_pole_figure(pole)


def test_linear_levels_span_the_data_and_carry_the_random_level() -> None:
    values = np.array([0.0, 2.0, 8.0])
    levels = ContourSpec(scale="linear", count=4).levels_for(values)
    assert levels[0] == pytest.approx(0.0)
    assert levels[-1] == pytest.approx(8.0)
    assert RANDOM_LEVEL_MRD in levels


def test_the_random_level_replaces_a_neighbour_rather_than_crowding_it() -> None:
    """A ladder carrying both 1.00 and 1.07 draws a band nobody can see."""

    values = np.array([0.0, 8.538])
    levels = ContourSpec(scale="linear", count=8).levels_for(values)
    assert RANDOM_LEVEL_MRD in levels
    spacing = np.diff(levels)
    assert float(np.min(spacing)) > 0.25 * float(np.max(spacing))


def test_geometric_levels_have_a_constant_ratio() -> None:
    levels = ContourSpec(
        scale="geometric", count=4, vmin=0.5, vmax=8.0, include_random_level=False
    ).levels_for(np.array([1.0]))
    ratios = levels[1:] / levels[:-1]
    assert_allclose(ratios, np.full(ratios.shape, 2.0), rtol=1e-12)


def test_a_geometric_ladder_over_zero_is_refused_with_a_usable_message() -> None:
    """Silently shifting the range would misstate every level on the figure."""

    spec = ContourSpec(scale="geometric")
    with pytest.raises(ValueError, match="strictly positive lower bound"):
        spec.levels_for(np.array([0.0, 4.0]))


def test_explicit_levels_are_used_verbatim_and_must_increase() -> None:
    spec = ContourSpec(scale="explicit", values=(0.5, 1.0, 2.0, 4.0))
    assert_allclose(spec.levels_for(np.array([0.0, 99.0])), np.array([0.5, 1.0, 2.0, 4.0]))
    with pytest.raises(ValueError, match="strictly increasing"):
        ContourSpec(scale="explicit", values=(2.0, 1.0))


def test_a_shared_range_covers_every_data_set_but_respects_declared_bounds() -> None:
    weak = np.array([0.4, 1.2])
    strong = np.array([0.1, 9.0])
    shared = ContourSpec().shared_across([weak, strong])
    assert shared.vmin == pytest.approx(0.1)
    assert shared.vmax == pytest.approx(9.0)
    pinned = ContourSpec(vmax=4.0).shared_across([weak, strong])
    assert pinned.vmax == pytest.approx(4.0)


def test_a_uniform_specimen_contours_at_one_mrd_everywhere() -> None:
    """The calibration of the density raster, as an identity rather than a tolerance.

    An isotropic specimen puts poles uniformly on the sphere, so its pole figure
    is flat at 1 multiple of a random distribution. Anything else means the
    raster is being evaluated on the wrong directions or normalized on the wrong
    domain — the two errors this pipeline could plausibly make.
    """

    figure = _uniform_pole_figure()
    x, _, density = pole_figure_density_grid(
        figure, style=PoleFigureStyle(resolution=41, halfwidth_deg=20.0)
    )
    inside = np.isfinite(density)
    assert x.shape == (41,)
    assert density.shape == (41, 41)
    # The corners of the raster square lie outside the projection disc.
    assert not inside[0, 0]
    assert_allclose(density[inside], np.ones(int(np.sum(inside))), atol=0.15)


def test_rotating_the_drawing_rotates_the_density_and_nothing_else() -> None:
    """A quarter turn of the drawing is a quarter turn of the raster, exactly.

    The rotation is applied to the raster points before they are unprojected,
    so it is a property of the drawing rather than of the data. On a square
    raster that makes the identity checkable to machine precision.
    """

    figure = _pole_figure((35.0, 45.0, 0.0), 10.0, seed=11)
    style = PoleFigureStyle(resolution=61, halfwidth_deg=15.0)
    _, _, plain = pole_figure_density_grid(figure, style=style)
    _, _, turned = pole_figure_density_grid(figure, style=replace(style, rotation_deg=90.0))
    assert_allclose(turned, np.rot90(plain), atol=1e-12, equal_nan=True)


def test_a_publication_figure_states_what_it_shows() -> None:
    figure = _pole_figure((0.0, 0.0, 0.0), 8.0, seed=5)
    style = PoleFigureStyle(resolution=61, halfwidth_deg=15.0, sample_label="Sample A")
    spec = build_pole_figure_contour_spec(figure, style=style)
    texts = {layer.text for layer in spec.text_layers}
    assert "Sample A" in texts
    assert "{111}" in texts
    assert any(text.startswith("max ") and "m.r.d." in text for text in texts)
    # The specimen axes are named where they meet the rim.
    assert {"RD", "TD"} <= texts
    # A publication figure carries no Cartesian axes or grid.
    assert spec.show_axes is False
    assert spec.grid is False
    assert spec.boundary_circle_radius == pytest.approx(np.sqrt(2.0))
    assert len(spec.contour_layers) == 1
    assert spec.contour_layers[0].values.shape == (61, 61)
    # And a marker on the strongest point.
    assert len(spec.marker_layers) == 1


def test_every_panel_of_a_shared_plate_is_drawn_at_the_same_levels() -> None:
    """The guarantee that makes a comparison plate a comparison.

    Contouring each panel on its own maximum makes a weak texture and a strong
    one look alike, which is the most common way a pole-figure plate misleads.
    """

    strong = _pole_figure((0.0, 0.0, 0.0), 7.0, seed=1)
    weak = _pole_figure((0.0, 0.0, 0.0), 40.0, seed=2)
    style = PoleFigureStyle(resolution=41, halfwidth_deg=18.0)
    plate = PoleFigureSet(
        figures=((strong,), (weak,)),
        sample_labels=("strong", "weak"),
        style=style,
    )
    spec = plate.build(suptitle="two samples")
    assert len(spec.panels) == 2
    assert spec.ncols == 1
    assert spec.shared_colorbar_label == "m.r.d."
    first, second = (panel.contour_layers[0] for panel in spec.panels)
    assert_allclose(first.levels, second.levels)
    # The shared ladder reaches the stronger sample's maximum, which a
    # self-scaled weak panel never would.
    weak_density = second.values[np.isfinite(second.values)]
    assert float(np.max(second.levels)) > 1.5 * float(np.max(weak_density))
    # No panel carries its own colour bar; the plate carries one.
    assert all(panel.contour_layers[0].colorbar_label is None for panel in spec.panels)
    labels = [
        {layer.text for layer in panel.text_layers} for panel in spec.panels
    ]
    assert "strong" in labels[0]
    assert "weak" in labels[1]


def test_an_unshared_plate_scales_every_panel_to_itself_and_says_so() -> None:
    strong = _pole_figure((0.0, 0.0, 0.0), 7.0, seed=1)
    weak = _pole_figure((0.0, 0.0, 0.0), 40.0, seed=2)
    plate = PoleFigureSet(
        figures=((strong,), (weak,)),
        sample_labels=("strong", "weak"),
        style=PoleFigureStyle(resolution=41, halfwidth_deg=18.0),
        shared_scale=False,
    )
    spec = plate.build()
    assert spec.shared_colorbar_label is None
    assert plate.shared_levels() is None
    first, second = (panel.contour_layers[0] for panel in spec.panels)
    assert float(np.max(first.levels)) > float(np.max(second.levels))
    assert all(panel.contour_layers[0].colorbar_label == "m.r.d." for panel in spec.panels)


def test_a_plate_refuses_a_ragged_grid_or_a_missing_label() -> None:
    figure = _pole_figure((0.0, 0.0, 0.0), 10.0, seed=4)
    with pytest.raises(ValueError, match="same, non-zero number"):
        PoleFigureSet(figures=((figure, figure), (figure,)), sample_labels=("a", "b"))
    with pytest.raises(ValueError, match="one label per sample"):
        PoleFigureSet(figures=((figure,), (figure,)), sample_labels=("a",))


def test_the_plate_and_the_single_figure_render(recwarn: pytest.WarningsRecorder) -> None:
    """Rendering must produce the grid asked for and warn about nothing."""

    strong = _pole_figure((0.0, 0.0, 0.0), 7.0, seed=1)
    weak = _pole_figure((0.0, 0.0, 0.0), 30.0, seed=2)
    style = PoleFigureStyle(resolution=41, halfwidth_deg=18.0)
    single = plot_pole_figure_contours(strong, style=style)
    plate = plot_pole_figure_comparison(
        PoleFigureSet(
            figures=((strong, weak), (weak, strong)),
            sample_labels=("first", "second"),
            style=style,
        ),
        suptitle="plate",
    )
    try:
        assert len(single.axes) >= 1
        # Four panels plus the shared colour bar.
        assert len(plate.axes) == 5
    finally:
        plt.close(single)
        plt.close(plate)
    assert [warning.message for warning in recwarn] == []
