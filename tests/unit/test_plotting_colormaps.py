from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pytex.plotting.colormaps import (
    OKABE_ITO_COLORS,
    PYTEX_COLORMAP_SPECS,
    ColormapSpec,
    categorical_colors,
    get_pytex_colormap,
    register_pytex_colormaps,
    srgb_to_lightness,
)


def test_register_is_idempotent_and_names_resolve() -> None:
    names = register_pytex_colormaps()
    assert names == register_pytex_colormaps()
    for name in names:
        assert name.startswith("pytex.")
        colormap = matplotlib.colormaps[name]
        assert colormap.N == 256


def test_get_pytex_colormap_rejects_unknown_name() -> None:
    with pytest.raises(ValueError, match="Unknown PyTex colormap"):
        get_pytex_colormap("pytex.bogus")


def test_sequential_ramps_have_monotone_decreasing_lightness() -> None:
    for spec in PYTEX_COLORMAP_SPECS:
        if spec.kind != "sequential":
            continue
        colormap = get_pytex_colormap(spec.name)
        samples = np.asarray(colormap(np.linspace(0.0, 1.0, 128)))[:, :3]
        lightness = srgb_to_lightness(samples)
        # magnitude must map to a perceptually ordered ramp (dark = more)
        assert np.all(np.diff(lightness) <= 1e-6), spec.name
        # white zero anchor: paper-background convention for m.r.d. plots
        assert lightness[0] == pytest.approx(100.0, abs=0.5)


def test_diverging_ramp_is_light_at_midpoint_dark_at_poles() -> None:
    colormap = get_pytex_colormap("pytex.diverging")
    samples = np.asarray(colormap(np.array([0.0, 0.5, 1.0])))[:, :3]
    lightness = srgb_to_lightness(samples)
    assert lightness[1] > lightness[0] + 25.0
    assert lightness[1] > lightness[2] + 25.0


def test_categorical_colors_fixed_order_and_cycling() -> None:
    three = categorical_colors(3)
    five = categorical_colors(5)
    # fixed assignment order: growing the class count never repaints earlier classes
    assert five[:3] == three
    assert three == OKABE_ITO_COLORS[:3]
    wrapped = categorical_colors(len(OKABE_ITO_COLORS) + 1)
    assert wrapped[-1] == OKABE_ITO_COLORS[0]
    assert categorical_colors(0) == ()
    with pytest.raises(ValueError, match="non-negative"):
        categorical_colors(-1)


def test_okabe_ito_palette_is_pairwise_distinct() -> None:
    rgb = np.array(
        [tuple(int(color[i : i + 2], 16) / 255.0 for i in (1, 3, 5)) for color in OKABE_ITO_COLORS]
    )
    for i in range(len(rgb)):
        for j in range(i + 1, len(rgb)):
            assert np.linalg.norm(rgb[i] - rgb[j]) > 0.15


def test_colormap_spec_validation() -> None:
    with pytest.raises(ValueError, match="kind"):
        ColormapSpec(name="x", colors=("#000000", "#ffffff"), kind="rainbow")
    with pytest.raises(ValueError, match="two color stops"):
        ColormapSpec(name="x", colors=("#000000",), kind="sequential")


def test_srgb_to_lightness_black_and_white_anchors() -> None:
    lightness = srgb_to_lightness(np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]))
    assert lightness[0] == pytest.approx(0.0, abs=1e-9)
    assert lightness[1] == pytest.approx(100.0, abs=1e-4)
    with pytest.raises(ValueError, match="shape"):
        srgb_to_lightness(np.array([0.0, 0.0, 0.0]))
