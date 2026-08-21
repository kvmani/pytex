"""Turning a sampled curve into the pieces of it that are worth drawing.

Every Kikuchi figure in the application has the same problem. A Kossel-cone edge
is a conic whose far branch runs to the horizon, and a band centre line curves
off the picture as soon as the detector is tilted; sampled along its whole
extent and handed to a renderer as one polyline, either is closed by a chord
straight across the pattern — a line the crystal never produced, drawn in the
middle of the figure the reader is measuring.

The fix is the same wherever it appears, so it lives here rather than in each
panel's service: drop the samples that are far outside the frame, and break what
remains wherever the trace jumped.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["clipped_runs"]


def clipped_runs(points: np.ndarray, width: float, height: float) -> list[list[list[float]]]:
    """Split a sampled trace into the runs that are near the picture.

    Purpose
    -------
    Produce polylines a renderer can draw without inventing segments. Points
    beyond one frame-size of margin are dropped — they are off the picture and
    keeping them would only stretch the coordinate range — and the remaining
    samples are broken wherever consecutive points are further apart than half
    the frame diagonal, which is what a trace passing through infinity looks
    like after sampling.

    Parameters
    ----------
    points : np.ndarray
        ``(n, 2)`` samples along one trace, in the frame's own pixels.
    width, height : float
        Frame size in the same pixels.

    Returns
    -------
    list of list of list of float
        Zero or more runs, each at least two points long, as plain nested lists
        ready for JSON. A trace that misses the frame gives an empty list, which
        is the honest answer rather than an empty polyline.
    """

    if points.size == 0:
        return []
    margin = max(width, height)
    inside = (
        (points[:, 0] > -margin)
        & (points[:, 0] < width + margin)
        & (points[:, 1] > -margin)
        & (points[:, 1] < height + margin)
    )
    runs: list[list[list[float]]] = []
    current: list[list[float]] = []
    limit = 0.5 * math.hypot(width, height)
    previous: np.ndarray | None = None
    for index, keep in enumerate(inside):
        point = points[index]
        if not keep:
            if len(current) > 1:
                runs.append(current)
            current = []
            previous = None
            continue
        if previous is not None and float(np.linalg.norm(point - previous)) > limit:
            if len(current) > 1:
                runs.append(current)
            current = []
        current.append([float(point[0]), float(point[1])])
        previous = point
    if len(current) > 1:
        runs.append(current)
    return runs
