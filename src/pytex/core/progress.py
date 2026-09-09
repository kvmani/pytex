"""Reporting how far a long computation has got, without knowing who is watching.

Why this is in the core
-----------------------
The code that knows how much work remains is the code doing the work: the loop
over multislice slices, over candidate phases, over refinement iterations. The
code that knows how to *show* it is the workbench shell, a CLI progress bar, or
a notebook widget — and a scientific module must not have to know which of those
is listening, or import a user interface to say "one of eight done".

So this module holds the mechanism and none of the presentation. A sink is
installed for the duration of a call by whichever surface is running it, and the
scientific code calls :func:`report` or iterates through :func:`tracking`. With
no sink installed both are no-ops, so instrumented code stays callable from a
test, a script or a notebook with no guard at the call site and no import of
anything above it in the stack.

What a fraction here means
--------------------------
**A measured fraction of the work this computation planned**, never elapsed time
against a guessed total. Code that cannot count its own work reports nothing,
and the surface showing it must then show elapsed time without a percentage. A
fabricated percentage is worse than none: it claims a measurement nobody took.

Examples
--------
>>> from pytex.core.progress import report, reporting, tracking
>>> seen = []
>>> class Collect:
...     def report(self, fraction, stage=None):
...         seen.append((round(fraction, 3), stage))
...         return True
>>> with reporting(Collect()):
...     for _ in tracking(range(4), stage="Scoring"):
...         pass
>>> seen
[(0.0, 'Scoring'), (0.25, 'Scoring'), (0.5, 'Scoring'), (0.75, 'Scoring'), (1.0, 'Scoring')]

Outside a :func:`reporting` block nothing is recorded and nothing fails:

>>> report(0.5, stage="ignored")
False

See Also
--------
pytex.app.progress : the workbench sink, which adds elapsed time and an ETA.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Protocol, TypeVar, runtime_checkable

__all__ = [
    "ProgressSink",
    "active_sink",
    "report",
    "reporting",
    "tracking",
]

T = TypeVar("T")


@runtime_checkable
class ProgressSink(Protocol):
    """Something that can be told how far a computation has got.

    One method, because that is the whole contract a scientific module needs to
    know about. Everything a surface adds — elapsed time, an estimated
    remaining time, rate limiting, where the bar is drawn — belongs to the sink
    and never to the code being measured.
    """

    def report(self, fraction: float, stage: str | None = None) -> bool:
        """Accept a completed fraction in ``[0, 1]``; return whether it was recorded."""
        ...


_ACTIVE: ContextVar[ProgressSink | None] = ContextVar("pytex_progress_sink", default=None)


def active_sink() -> ProgressSink | None:
    """The sink installed on this thread's context, if any."""
    return _ACTIVE.get()


@contextmanager
def reporting(sink: ProgressSink | None) -> Iterator[ProgressSink | None]:
    """Install ``sink`` for the duration of the block.

    A :class:`~contextvars.ContextVar` rather than a global, so two computations
    running on two threads report to their own sink and not to each other's.
    Passing ``None`` deliberately silences reporting inside the block.
    """
    token = _ACTIVE.set(sink)
    try:
        yield sink
    finally:
        _ACTIVE.reset(token)


def report(fraction: float, *, stage: str | None = None) -> bool:
    """Report a completed fraction to the active sink.

    Parameters
    ----------
    fraction : float
        Fraction of this computation's planned work that is done.
    stage : str, optional
        What is being worked on now, in words the person waiting can read.

    Returns
    -------
    bool
        Whether a sink recorded it. ``False`` when none is installed, which is
        the normal case outside an interactive surface.
    """
    sink = _ACTIVE.get()
    if sink is None:
        return False
    return sink.report(fraction, stage)


def tracking(
    items: Iterable[T], *, total: int | None = None, stage: str | None = None
) -> Iterator[T]:
    """Iterate ``items``, reporting the fraction consumed.

    Parameters
    ----------
    items : iterable
        The work to do.
    total : int, optional
        Number of items. Taken from ``len(items)`` when the iterable is sized.
    stage : str, optional
        What this loop is doing.

    Yields
    ------
    Each item of ``items``, unchanged.

    Notes
    -----
    The fraction reported before item ``i`` is ``i / total`` — work *completed*,
    not work started — so the bar reaches 100% only once the loop has finished.
    An unsized iterable with no ``total`` is yielded without reporting: a
    fraction of an unknown quantity is not a measurement.
    """
    if total is None:
        try:
            total = len(items)  # type: ignore[arg-type]
        except TypeError:
            total = None
    if not total:
        yield from items
        return
    for index, item in enumerate(items):
        report(index / total, stage=stage)
        yield item
    report(1.0, stage=stage)
