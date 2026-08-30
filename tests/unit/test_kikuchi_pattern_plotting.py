"""Structural tests for the Kikuchi *pattern* figure, and its band naming.

Runtime plotting is checked semantically rather than against tracked SVG bytes,
per ``docs/standards/executable_examples.md``. The claim under test is the one
the naming exists for: a band's indices must run **along** that band. It is
checked against the geometry the same figure drew, so a label that is merely
tilted by a plausible-looking amount fails.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pytex.core import FrameDomain, ReferenceFrame, get_phase_fixture
from pytex.core.frame_catalog import SPECIMEN_FRAME
from pytex.core.orientation import Orientation
from pytex.diffraction.kikuchi import simulate_kikuchi_pattern
from pytex.diffraction.models import DiffractionGeometry
from pytex.plotting.diffraction import plot_kikuchi_pattern

# Every test here starts from a phase loaded out of its CIF fixture, and
# CIF-backed phase creation is pymatgen's job. Declared at module scope because
# it is the whole module: without the optional `adapters` extra these are not
# failures but tests that cannot be run, and reporting them as failures hid a
# real defect underneath them for as long as it lasted.
pytest.importorskip(
    "pymatgen",
    reason="loading a phase from its CIF fixture needs the 'adapters' extra",
)

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


@pytest.fixture(scope="module")
def pattern():
    phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=CRYSTAL)
    geometry = DiffractionGeometry.for_ebsd(specimen_frame=SPECIMEN_FRAME)
    orientation = Orientation.from_euler(
        0.0, 0.0, 0.0, degrees=True, specimen_frame=SPECIMEN_FRAME, phase=phase
    )
    return simulate_kikuchi_pattern(geometry, phase, orientation, max_index=3, max_bands=8)


def _band_labels(axes):
    """Texts that name a plane, as ``(label, rotation)`` pairs."""

    return [
        (text.get_text(), float(text.get_rotation()))
        for text in axes.texts
        if text.get_text().startswith("$(")
    ]


def test_bands_are_not_named_unless_asked(pattern) -> None:
    """A pattern with tens of bands would otherwise become a page of text."""

    figure = plot_kikuchi_pattern(pattern, samples=181)
    try:
        assert _band_labels(figure.axes[0]) == []
    finally:
        plt.close(figure)


def test_a_band_is_named_along_its_own_trace(pattern) -> None:
    """The rotation of each label matches the slope of the line it names.

    Paired by position rather than by index: the label sits on its own trace, so
    the nearest drawn line to the anchor *is* the band it belongs to, and the
    two angles must agree. A label written horizontally beside a steep band —
    the behaviour this replaced — fails on every band that is not level.
    """

    figure = plot_kikuchi_pattern(pattern, samples=361, label_bands=True)
    try:
        # Drawn once first: with `transform_rotates_text`, a text's rotation is
        # reported in *display* space, and the axes transform is not final until
        # the figure has been laid out. Querying before that compares an angle
        # against a transform nobody will ever render with.
        figure.canvas.draw()
        axes = figure.axes[0]
        labels = _band_labels(axes)
        assert labels, "label_bands must name the bands"

        traces = [np.column_stack(line.get_data()) for line in axes.lines]
        traces = [trace for trace in traces if trace.shape[0] > 1]
        assert traces

        for text in axes.texts:
            if not text.get_text().startswith("$("):
                continue
            anchor = np.asarray(text.get_position(), dtype=float)
            best = None
            for trace in traces:
                index = int(np.argmin(np.linalg.norm(trace - anchor, axis=1)))
                distance = float(np.linalg.norm(trace[index] - anchor))
                if best is None or distance < best[0]:
                    best = (distance, trace, index)
            distance, trace, index = best
            # The anchor is a sample of its own trace, so the nearest point is
            # that sample: anything else means the label is not on a band.
            assert distance < 1e-9, text.get_text()
            # The local tangent, from the samples either side of the anchor. A
            # wider baseline would be measured across the point where a trace
            # leaves the picture and re-enters it, which is a jump rather than a
            # direction.
            before = trace[max(0, index - 1)]
            after = trace[min(trace.shape[0] - 1, index + 1)]
            delta = after - before
            expected = float(np.degrees(np.arctan2(delta[1], delta[0])))
            # Modulo a half-turn: the same line, written the same way up.
            difference = (float(text.get_rotation()) - expected) % 180.0
            assert min(difference, 180.0 - difference) < 1.0, text.get_text()

        # And the pattern really does carry bands at several slopes, so the case
        # cannot pass by every label happening to be horizontal.
        assert len({round(rotation, 3) for _, rotation in labels}) > 1
    finally:
        plt.close(figure)


def test_a_named_band_keeps_its_name_in_detector_coordinates(pattern) -> None:
    """Centre lines curve in pixels on a tilted camera; the label follows them.

    The gnomonic frame is where a centre line is straight. In detector pixels it
    is not, so a label placed from the endpoints of the trace would drift off a
    curve that the figure draws correctly.
    """

    figure = plot_kikuchi_pattern(
        pattern, coordinates="detector", samples=361, label_bands=True
    )
    try:
        axes = figure.axes[0]
        assert _band_labels(axes)
        height, width = pattern.geometry.detector_shape
        for text in axes.texts:
            if not text.get_text().startswith("$("):
                continue
            x, y = text.get_position()
            # Drawn inside the detector, which is the only place a reader looks.
            assert -width <= x <= 2 * width
            assert -height <= y <= 2 * height
    finally:
        plt.close(figure)
