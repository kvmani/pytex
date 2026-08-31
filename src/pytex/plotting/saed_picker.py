"""Click-to-pick front end for building a measured SAED pattern.

Picking spots off a micrograph is inherently interactive, but an interactive
tool that cannot be tested is a liability in a scientific library. So the
picking *logic* lives in :class:`SpotPickerState` — a plain object with no
Matplotlib dependency, fully exercised headlessly — and
:class:`SAEDSpotPicker` is a thin event adapter over it.

The output is a :class:`~pytex.diffraction.solving.MeasuredSAEDPattern`, which
serializes to YAML. That file, not the click session, is the reproducibility
boundary: a solved pattern must be reproducible from a committed text file,
whether the spots were clicked or typed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from pytex.diffraction.solving import (
    MeasuredSAEDPattern,
    MeasuredSpot,
    PatternCalibration,
)


def _matplotlib() -> Any:
    """pyplot, imported on demand. See `pytex.plotting.crystal3d._to_hex`."""

    import matplotlib.pyplot as plt

    return plt


@dataclass
class SpotPickerState:
    """The picking session's state machine, independent of any GUI.

    Holds the picked positions and the transmitted-beam centre, and supports
    adding, removing the nearest pick, undoing, moving the centre, and clearing.
    Every operation is a plain method so the behaviour can be tested without a
    display, which is the point: the science is in what gets picked, not in how
    the clicks arrive.

    Positions are stored in whatever units the picker was told to work in, and
    are handed to `PatternCalibration` unchanged.
    """

    centre: tuple[float, float] = (0.0, 0.0)
    positions: list[tuple[float, float]] = field(default_factory=list)
    labels: list[str | None] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.positions) != len(self.labels):
            raise ValueError("positions and labels must be the same length.")

    def __len__(self) -> int:
        return len(self.positions)

    def add(self, x: float, y: float, *, label: str | None = None) -> int:
        """Record a pick and return its index."""

        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError("A picked position must be finite.")
        self.positions.append((float(x), float(y)))
        self.labels.append(label)
        return len(self.positions) - 1

    def set_centre(self, x: float, y: float) -> None:
        """Move the transmitted-beam position.

        Every spot is measured relative to this, so moving it re-scales the whole
        pattern; it does not remove or renumber any pick.
        """

        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError("The centre must be finite.")
        self.centre = (float(x), float(y))

    def remove_nearest(self, x: float, y: float, *, radius: float | None = None) -> int | None:
        """Remove the pick nearest ``(x, y)`` and return its index, or ``None``.

        ``radius`` bounds how far a click may be from a pick and still remove it;
        without it the nearest pick is removed however far away the click was,
        which is rarely what a user means.
        """

        if not self.positions:
            return None
        distances = np.linalg.norm(
            np.asarray(self.positions, dtype=np.float64) - np.array([x, y]), axis=1
        )
        index = int(np.argmin(distances))
        if radius is not None and distances[index] > radius:
            return None
        del self.positions[index]
        del self.labels[index]
        return index

    def undo(self) -> int | None:
        """Remove the most recent pick and return its index, or ``None`` if empty."""

        if not self.positions:
            return None
        self.positions.pop()
        self.labels.pop()
        return len(self.positions)

    def clear(self) -> None:
        """Discard every pick, keeping the centre."""

        self.positions.clear()
        self.labels.clear()

    def to_pattern(
        self, *, name: str, calibration: PatternCalibration
    ) -> MeasuredSAEDPattern:
        """Build the measured pattern, taking the centre from this session.

        The session's own centre wins over the calibration's, because the centre
        is something the user set by clicking; the calibration supplies the units
        and the camera constant.
        """

        from dataclasses import replace

        return MeasuredSAEDPattern(
            name=name,
            spots=tuple(
                MeasuredSpot(position=position, label=label)
                for position, label in zip(self.positions, self.labels, strict=True)
            ),
            calibration=replace(calibration, centre=self.centre),
        )


class SAEDSpotPicker:
    """Matplotlib click-to-pick front end producing a `MeasuredSAEDPattern`.

    Purpose: turn a displayed diffraction pattern into a list of spot
    coordinates without leaving Python. Left-click adds a spot, right-click
    removes the nearest one, and middle-click moves the transmitted-beam centre;
    ``u`` undoes the last pick and ``c`` clears every pick.

    When to use: interactively, in a notebook or a script with an interactive
    Matplotlib backend. For reproducible work, save the result with
    `save_yaml` and solve from the file — the file is the contract, and a
    committed one gives the same answer on any machine.

    Inputs: ``image`` — an optional 2-D array to display beneath the picks;
    ``calibration`` — the `PatternCalibration` describing the coordinate units
    (defaults to reciprocal angstroms); ``extent`` — the Matplotlib extent to
    map image pixels onto those units; ``state`` — an existing
    `SpotPickerState` to continue.

    All the picking logic lives in `SpotPickerState`, which is testable without
    a display; this class only wires events to it.
    """

    def __init__(
        self,
        image: Any | None = None,
        *,
        calibration: PatternCalibration | None = None,
        extent: tuple[float, float, float, float] | None = None,
        state: SpotPickerState | None = None,
        removal_radius: float | None = None,
    ) -> None:
        self.calibration = calibration or PatternCalibration()
        self.state = state or SpotPickerState(centre=self.calibration.centre)
        self.removal_radius = removal_radius
        self._image = image
        self._extent = extent
        self.figure: Any | None = None
        self.axes: Any | None = None
        self._connections: list[int] = []
        self._scatter: Any | None = None
        self._centre_marker: Any | None = None

    def show(self, *, ax: Any | None = None) -> SAEDSpotPicker:
        """Draw the pattern and connect the mouse and key handlers."""

        plt = _matplotlib()
        if ax is None:
            self.figure, self.axes = plt.subplots()
        else:
            self.axes = ax
            self.figure = ax.figure
        if self._image is not None:
            self.axes.imshow(
                self._image, extent=self._extent, origin="lower", cmap="gray"
            )
        self.axes.set_aspect("equal")
        self.axes.set_xlabel(f"x ({self.calibration.units})")
        self.axes.set_ylabel(f"y ({self.calibration.units})")
        self.axes.set_title(
            "left click: add spot | right click: remove | middle click: set centre | "
            "u: undo | c: clear"
        )
        self._connections = [
            self.figure.canvas.mpl_connect("button_press_event", self._on_click),
            self.figure.canvas.mpl_connect("key_press_event", self._on_key),
        ]
        self._redraw()
        return self

    def disconnect(self) -> None:
        """Detach the event handlers, leaving the figure and the picks intact."""

        if self.figure is not None:
            for connection in self._connections:
                self.figure.canvas.mpl_disconnect(connection)
        self._connections = []

    def _on_click(self, event: Any) -> None:
        if event.inaxes is not self.axes or event.xdata is None or event.ydata is None:
            return
        if event.button == 1:
            self.state.add(float(event.xdata), float(event.ydata))
        elif event.button == 3:
            self.state.remove_nearest(
                float(event.xdata), float(event.ydata), radius=self.removal_radius
            )
        elif event.button == 2:
            self.state.set_centre(float(event.xdata), float(event.ydata))
        self._redraw()

    def _on_key(self, event: Any) -> None:
        if event.key == "u":
            self.state.undo()
        elif event.key == "c":
            self.state.clear()
        else:
            return
        self._redraw()

    def _redraw(self) -> None:
        if self.axes is None:
            return
        if self._scatter is not None:
            self._scatter.remove()
            self._scatter = None
        if self._centre_marker is not None:
            self._centre_marker.remove()
            self._centre_marker = None
        if self.state.positions:
            coordinates = np.asarray(self.state.positions, dtype=np.float64)
            self._scatter = self.axes.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                s=60,
                facecolors="none",
                edgecolors="#2563eb",
                linewidths=1.4,
                gid="pytex-picker:spots",
            )
        self._centre_marker = self.axes.scatter(
            [self.state.centre[0]],
            [self.state.centre[1]],
            marker="+",
            s=120,
            color="#7c3aed",
            gid="pytex-picker:centre",
        )
        if self.figure is not None:
            self.figure.canvas.draw_idle()

    def pattern(self, *, name: str = "picked_pattern") -> MeasuredSAEDPattern:
        """The picks as a `MeasuredSAEDPattern`, ready to solve or serialize."""

        return self.state.to_pattern(name=name, calibration=self.calibration)

    def save_yaml(self, path: str | Path, *, name: str = "picked_pattern") -> Path:
        """Write the picks to a measured-pattern YAML file and return the path."""

        return self.pattern(name=name).to_yaml(path)


__all__ = ["SAEDSpotPicker", "SpotPickerState"]
