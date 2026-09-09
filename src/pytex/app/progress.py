"""The workbench's progress sink: fraction, elapsed time, and an estimated ETA.

Why this exists
---------------
A workbench operation that runs for thirty seconds and says nothing is
indistinguishable, to the person waiting, from one that has hung. The console
already carries a progress record type — a fraction and an estimated remaining
time — and the shell already polls for records while a call is in flight. What
was missing was a way for the code doing the work to say how far along it is,
and something to turn that into a record.

:mod:`pytex.core.progress` is the first half: a one-method sink protocol that
scientific code reports into without knowing who is listening. This module is
the second: the sink the workbench installs, which adds what a person waiting
needs and a scientific module has no business computing — how long this has
taken, how long the rest is likely to take, and a rate limit so a tight loop
cannot flood a bounded log buffer.

What the two times mean
-----------------------
The elapsed time is measured. The remaining time is an **extrapolation of the
rate observed so far onto the work not yet done**, and is presented as an
estimate everywhere it is shown. It is deliberately not reported before any
progress arrives — a rate cannot be measured from no completed work — nor at
completion, where nothing remains.

See Also
--------
pytex.core.progress : the sink protocol and the ``report``/``tracking`` helpers.
pytex.app.contracts.execute : installs one of these for the duration of a call.
"""

from __future__ import annotations

import time

from pytex.app.logbook import APP_LOG, Logbook

__all__ = ["MINIMUM_TICK_INTERVAL_S", "ProgressReporter"]

#: Shortest interval between emitted ticks, in seconds.
#:
#: A tight loop can report hundreds of thousands of times. Emitting every one
#: would push every other message out of a bounded log buffer and cost more than
#: the work being measured. The first tick, the last, and every change of stage
#: always go out; the rest are thinned to this rate, which is well below the
#: interval the shell polls at, so nothing visible is lost.
MINIMUM_TICK_INTERVAL_S = 0.4


class ProgressReporter:
    """Turns reported fractions into console progress records.

    Purpose
    -------
    Implements :class:`pytex.core.progress.ProgressSink` for the workbench: each
    accepted fraction becomes one progress record carrying the fraction, the
    stage, the measured elapsed time and the estimated remaining time, which the
    shell renders as a bar.

    Parameters
    ----------
    task : str
        Identifier grouping the ticks of one operation so each replaces the last
        rather than accumulating. The operation id is used.
    title : str
        Human-readable operation name, used when no stage has been named.
    book : Logbook, optional
        Where ticks are emitted. Defaults to the application log.
    minimum_interval_s : float, optional
        Rate limit for emitted ticks. See :data:`MINIMUM_TICK_INTERVAL_S`.

    Notes
    -----
    The reported fraction is monotonic. A caller reporting a smaller fraction
    than it last reported is measuring a different quantity, and a bar that ran
    backwards reads as a fault, so the largest fraction seen is kept.
    """

    def __init__(
        self,
        task: str,
        title: str,
        *,
        book: Logbook | None = None,
        minimum_interval_s: float = MINIMUM_TICK_INTERVAL_S,
    ) -> None:
        self.task = task
        self.title = title
        self._book = book if book is not None else APP_LOG
        self._minimum_interval_s = float(minimum_interval_s)
        self._started = time.monotonic()
        self._last_emit = 0.0
        self._emitted = 0
        self._fraction = 0.0
        self._stage: str | None = None

    @property
    def elapsed_s(self) -> float:
        """Seconds since this reporter was created."""
        return time.monotonic() - self._started

    @property
    def fraction(self) -> float:
        """Largest fraction reported so far, in ``[0, 1]``."""
        return self._fraction

    @property
    def emitted(self) -> int:
        """How many ticks reached the log, after rate limiting."""
        return self._emitted

    def estimated_remaining_s(self) -> float | None:
        """Extrapolate the remaining seconds from the rate measured so far.

        ``None`` before any progress is reported, because a rate cannot be
        measured from no completed work, and at completion, because nothing
        remains. Both are cases where a surface must show elapsed time rather
        than invent a countdown.
        """
        if self._fraction <= 0.0 or self._fraction >= 1.0:
            return None
        elapsed = self.elapsed_s
        if elapsed <= 0.0:
            return None
        return elapsed * (1.0 - self._fraction) / self._fraction

    def report(self, fraction: float, stage: str | None = None, *, force: bool = False) -> bool:
        """Record a completed fraction; return whether a record was emitted.

        Parameters
        ----------
        fraction : float
            Completed fraction, clamped into ``[0, 1]``. A value below the
            largest already reported does not move the bar backwards.
        stage : str, optional
            What is being worked on now. A change of stage always emits,
            whatever the rate limit, because the stage is the part of the
            message that says where the work has got to.
        force : bool, optional
            Emit regardless of the rate limit. Used for the closing tick.
        """
        clamped = min(max(float(fraction), 0.0), 1.0)
        self._fraction = max(self._fraction, clamped)
        stage_changed = stage is not None and stage != self._stage
        if stage is not None:
            self._stage = stage

        now = time.monotonic()
        first = self._emitted == 0
        due = force or first or stage_changed or (now - self._last_emit) >= self._minimum_interval_s
        if not due:
            return False
        self._last_emit = now
        self._emitted += 1

        remaining = self.estimated_remaining_s()
        headline = self._stage or self.title
        self._book.progress(
            self.task,
            self._fraction,
            message=f"{headline}: {round(self._fraction * 100.0)}% complete.",
            eta_seconds=remaining,
            source=self.task,
            detail={
                "stage": self._stage,
                "elapsed_s": round(self.elapsed_s, 3),
                "title": self.title,
            },
        )
        return True

    def finish(self) -> None:
        """Close the bar, so one left at 97% is completed rather than abandoned.

        An operation that reported nothing emits nothing here: a lone 100% tick
        would claim a measurement of work that was never counted.
        """
        if self._emitted == 0:
            return
        self.report(1.0, force=True)
