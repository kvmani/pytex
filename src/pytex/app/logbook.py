"""The application's log: one graded, ordered stream of what happened.

Why a logbook rather than scattered messages
--------------------------------------------
Before this module the workbench told the user things in four unrelated places:
a toast that vanished after nine seconds, an "activity" strip that counted
service calls, an error line beside a control, and Python's own ``logging``
output in a terminal the desktop user never sees. Nothing joined them, so the
answer to "what did the application just do, and in what order" did not exist.

A logbook is that answer. Every notable event — a calculation starting, a spot
being picked, an input being rejected, a long simulation reaching 50% — becomes
a :class:`LogRecord` with a severity, a monotonic sequence number, a wall-clock
time, the surface it came from, and optional structured detail. The records go
into one bounded, thread-safe buffer that both shells read, and the frontend
renders them in a console pinned to the bottom of the window.

How records reach the user
--------------------------
Two paths, deliberately:

1. **Attached to the call envelope.** :func:`collecting` captures whatever a
   service emitted while handling one operation, and
   :func:`pytex.app.contracts.execute` puts those records in the response. They
   arrive with the result, so a message about the calculation cannot land before
   or after the thing it describes.
2. **Polled from the buffer.** Records emitted outside any call — the server
   binding a port, a background thread, a stdlib logger routed here through
   :class:`LogbookHandler` — are read over ``GET /api/log?since=N``. The
   sequence number is what makes that poll exact rather than approximate: the
   client asks for what it has not seen, not for "the last few".

The frontend merges the two by sequence number, so a record delivered on both
paths appears once.

Severity
--------
Seven levels, ordered by :attr:`LogLevel.severity`. ``NOTICE`` is the "important
but not wrong" level — a result worth noticing, a convention that was applied, a
fallback that was taken — and it exists because forcing such messages to be
either ``INFO`` (ignored) or ``WARNING`` (alarming) is what makes warnings stop
meaning anything.

Example
-------
>>> book = Logbook()
>>> _ = book.info("Spot 1 selected", source="tem", detail={"x": 512.0, "y": 480.0})
>>> _ = book.progress("cbed", 0.5, eta_seconds=150.0, source="cbed")
>>> [record.message for record in book.records()]
['Spot 1 selected', '50% progress. ETA: 2 min 30 sec']
>>> book.records()[-1].level.name
'PROGRESS'
"""

from __future__ import annotations

import logging
import threading
import time
from collections import deque
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

__all__ = [
    "APP_LOG",
    "LOG_RECORD_SCHEMA",
    "LogLevel",
    "LogRecord",
    "Logbook",
    "LogbookHandler",
    "ProgressReporter",
    "collecting",
    "format_duration",
    "install_logging_bridge",
]

#: Schema identifier of one wire-format log record.
LOG_RECORD_SCHEMA = "pytex.log_record/1"


class LogLevel(Enum):
    """How much a record matters, and what the console should call it.

    Purpose
    -------
    The severity vocabulary shared by Python, the wire format, and the console's
    filter chips. Declaring the label and the wire token here — rather than
    letting the frontend invent them — is what stops the console from showing a
    level Python never emits, or hiding one it does.

    Attributes
    ----------
    token : str
        Lower-case identifier used on the wire and as a CSS modifier.
    label : str
        What the console calls this level in its filter row.
    severity : int
        Ordering key. Higher is more serious; :data:`PROGRESS` deliberately
        sits below :data:`INFO` because a progress tick is the least important
        thing in the book and must be the first thing a filter drops.
    """

    PROGRESS = ("progress", "Progress", 5)
    DEBUG = ("debug", "Debug", 10)
    INFO = ("info", "Info", 20)
    NOTICE = ("notice", "Important", 25)
    SUCCESS = ("success", "Success", 27)
    WARNING = ("warning", "Warning", 30)
    ERROR = ("error", "Error", 40)
    CRITICAL = ("critical", "Critical", 50)

    def __init__(self, token: str, label: str, severity: int) -> None:
        self.token = token
        self.label = label
        self.severity = severity

    @classmethod
    def from_token(cls, token: str) -> LogLevel:
        """Return the level with this wire token.

        Raises
        ------
        ValueError
            If no level uses the token, which is a client sending a severity
            this application does not define.
        """

        for level in cls:
            if level.token == token:
                return level
        allowed = ", ".join(level.token for level in cls)
        raise ValueError(f"Unknown log level {token!r}; expected one of {allowed}.")


def format_duration(seconds: float) -> str:
    """Render a duration the way a person reads a wait.

    Purpose
    -------
    Produces the ``"2 min 30 sec"`` half of a progress message. Written here
    rather than in each caller so that every estimate in the application is
    phrased identically — an ETA that says "150 s" in one place and "2.5 min"
    in another reads as two different quantities.

    Parameters
    ----------
    seconds : float
        A non-negative duration. Negative and non-finite values are reported as
        ``"unknown"`` rather than rendered, because a nonsense ETA on screen is
        worse than an admitted absence of one.

    Returns
    -------
    str
        ``"12 sec"``, ``"2 min 30 sec"``, ``"1 hr 05 min"``, or ``"unknown"``.

    Examples
    --------
    >>> format_duration(150.0)
    '2 min 30 sec'
    >>> format_duration(9.4)
    '9 sec'
    >>> format_duration(3930.0)
    '1 hr 05 min'
    """

    if seconds != seconds or seconds in (float("inf"), float("-inf")) or seconds < 0.0:
        return "unknown"
    whole = round(seconds)
    if whole < 60:
        return f"{whole} sec"
    if whole < 3600:
        return f"{whole // 60} min {whole % 60:02d} sec"
    return f"{whole // 3600} hr {(whole % 3600) // 60:02d} min"


@dataclass(frozen=True)
class LogRecord:
    """One thing that happened, graded and timed.

    Attributes
    ----------
    sequence : int
        Monotonic per-:class:`Logbook` counter, starting at 1. The client polls
        with ``since=<last seen sequence>``, so this is what makes "everything I
        have not seen" an exact request rather than a guess at a time window.
    time : float
        Unix timestamp of the emission, for display.
    level : LogLevel
        Severity.
    message : str
        One sentence, already written for a person. Services do not send format
        strings; they send the finished sentence, because the console has no
        vocabulary with which to complete one.
    source : str
        Where it came from: a panel id (``"tem"``), an operation id
        (``"cbed.pattern"``), or a subsystem name (``"server"``). The console
        shows it and filters on it.
    detail : mapping
        Structured extras — picked coordinates, the offending field, the
        rejected value. Rendered as a key/value list under the message when the
        entry is expanded, and never required to understand it.
    task : str, optional
        Groups the ticks of one long operation. Records sharing a task replace
        one another in the console rather than stacking, which is what keeps a
        thousand progress ticks from burying everything else.
    progress : float, optional
        Completion fraction in ``[0, 1]``.
    eta_seconds : float, optional
        Estimated remaining time.
    """

    sequence: int
    time: float
    level: LogLevel
    message: str
    source: str = "app"
    detail: Mapping[str, Any] = field(default_factory=dict)
    task: str | None = None
    progress: float | None = None
    eta_seconds: float | None = None

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this record."""

        payload: dict[str, Any] = {
            "schema": LOG_RECORD_SCHEMA,
            "sequence": self.sequence,
            "time": self.time,
            "level": self.level.token,
            "message": self.message,
            "source": self.source,
        }
        if self.detail:
            payload["detail"] = dict(self.detail)
        if self.task is not None:
            payload["task"] = self.task
        if self.progress is not None:
            payload["progress"] = self.progress
        if self.eta_seconds is not None:
            payload["eta_seconds"] = self.eta_seconds
        return payload

    def describe(self) -> str:
        """Return the record as one line of plain text.

        This is what the console's "Copy" button writes, so a user reporting a
        problem can paste the log into an issue and have it stay readable.
        """

        stamp = time.strftime("%H:%M:%S", time.localtime(self.time))
        line = f"{stamp}  {self.level.label.upper():<9} [{self.source}] {self.message}"
        if self.detail:
            extras = ", ".join(f"{key}={value}" for key, value in self.detail.items())
            line = f"{line} ({extras})"
        return line


class Logbook:
    """A bounded, thread-safe buffer of :class:`LogRecord`.

    Purpose
    -------
    The single place application events are collected. Bounded because a session
    left open for a day must not grow without limit, and thread-safe because the
    HTTP server answers requests on several threads at once — two colleagues
    running calculations must not interleave into a corrupt deque.

    When to use it
    --------------
    Call the level methods from anywhere in the application layer that does
    something a user would want narrated. Do not use it for control flow, and do
    not use it in :mod:`pytex.core`, :mod:`pytex.ebsd` or the other scientific
    packages: those are a library, and a library that writes to an application's
    console is a library that cannot be used without one.

    Parameters
    ----------
    capacity : int
        How many records to retain. Older records are dropped from the front.

    Examples
    --------
    >>> book = Logbook(capacity=2)
    >>> for index in range(3):
    ...     _ = book.info(f"message {index}")
    >>> [record.message for record in book.records()]
    ['message 1', 'message 2']
    """

    def __init__(self, capacity: int = 500) -> None:
        if capacity < 1:
            raise ValueError("Logbook capacity must be at least 1.")
        self._capacity = capacity
        self._records: deque[LogRecord] = deque(maxlen=capacity)
        self._lock = threading.Lock()
        self._sequence = 0

    @property
    def capacity(self) -> int:
        """How many records this book retains."""

        return self._capacity

    def emit(
        self,
        level: LogLevel,
        message: str,
        *,
        source: str = "app",
        detail: Mapping[str, Any] | None = None,
        task: str | None = None,
        progress: float | None = None,
        eta_seconds: float | None = None,
    ) -> LogRecord:
        """Append one record and return it.

        The record is also handed to the active :func:`collecting` collector, if
        one is open, so a service's narration travels back with the result of
        the call that produced it.
        """

        with self._lock:
            self._sequence += 1
            record = LogRecord(
                sequence=self._sequence,
                time=time.time(),
                level=level,
                message=message,
                source=source,
                detail=dict(detail or {}),
                task=task,
                progress=progress,
                eta_seconds=eta_seconds,
            )
            self._records.append(record)
        collector = _COLLECTOR.get()
        if collector is not None:
            collector.append(record)
        return record

    def debug(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit a record only a developer wants."""

        return self.emit(LogLevel.DEBUG, message, **kwargs)

    def info(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit an ordinary narration record."""

        return self.emit(LogLevel.INFO, message, **kwargs)

    def notice(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit an important-but-not-wrong record."""

        return self.emit(LogLevel.NOTICE, message, **kwargs)

    def success(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit a record for something that completed as intended."""

        return self.emit(LogLevel.SUCCESS, message, **kwargs)

    def warning(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit a record for something suspect that did not stop the work."""

        return self.emit(LogLevel.WARNING, message, **kwargs)

    def error(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit a record for work that failed."""

        return self.emit(LogLevel.ERROR, message, **kwargs)

    def critical(self, message: str, **kwargs: Any) -> LogRecord:
        """Emit a record for a failure that leaves the application unusable."""

        return self.emit(LogLevel.CRITICAL, message, **kwargs)

    def progress(
        self,
        task: str,
        fraction: float,
        *,
        message: str | None = None,
        eta_seconds: float | None = None,
        source: str = "app",
        detail: Mapping[str, Any] | None = None,
    ) -> LogRecord:
        """Emit a progress tick for a long operation.

        Parameters
        ----------
        task : str
            Identifier grouping the ticks of one operation, so the console
            updates a single line rather than appending a hundred.
        fraction : float
            Completion in ``[0, 1]``; values outside are clamped, because a
            progress bar that exceeds its own end is a defect the user should
            not be shown.
        message : str, optional
            Overrides the generated sentence.
        eta_seconds : float, optional
            Remaining time. Appended to the generated sentence when present.
        """

        clamped = min(max(float(fraction), 0.0), 1.0)
        if message is None:
            message = f"{clamped * 100:.0f}% progress."
            if eta_seconds is not None:
                message = f"{message} ETA: {format_duration(eta_seconds)}"
        return self.emit(
            LogLevel.PROGRESS,
            message,
            source=source,
            detail=detail,
            task=task,
            progress=clamped,
            eta_seconds=eta_seconds,
        )

    def records(
        self, *, since: int = 0, minimum_level: LogLevel | None = None
    ) -> tuple[LogRecord, ...]:
        """Return retained records after ``since``, oldest first.

        Parameters
        ----------
        since : int
            Return only records whose sequence number exceeds this. ``0``
            returns everything retained.
        minimum_level : LogLevel, optional
            Drop records less severe than this.

        Notes
        -----
        A client that falls further behind than :attr:`capacity` silently misses
        records. That is the honest consequence of a bounded buffer, and the
        console reports it by noticing that the first sequence it received is
        more than one past the last it held.
        """

        with self._lock:
            snapshot = tuple(self._records)
        selected = tuple(record for record in snapshot if record.sequence > since)
        if minimum_level is not None:
            floor = minimum_level.severity
            selected = tuple(record for record in selected if record.level.severity >= floor)
        return selected

    def latest_sequence(self) -> int:
        """The sequence number of the most recent record, or ``0`` if empty."""

        with self._lock:
            return self._sequence

    def clear(self) -> None:
        """Drop every retained record without resetting the sequence counter.

        The counter is deliberately not reset: a client holding ``since=42``
        must never be handed a *different* record 42 after a clear.
        """

        with self._lock:
            self._records.clear()


#: The application-wide logbook. One per process, like the registry.
APP_LOG = Logbook()

#: Records emitted inside the innermost open :func:`collecting` block.
_COLLECTOR: ContextVar[list[LogRecord] | None] = ContextVar("pytex_logbook_collector", default=None)

#: Distinguishes one :class:`ProgressReporter` run from the next. Process-wide
#: and lock-guarded, because two calls served on two threads must not be handed
#: the same task id and start overwriting each other's console line.
_RUN_COUNTER = 0
_RUN_COUNTER_LOCK = threading.Lock()


@contextmanager
def collecting() -> Iterator[list[LogRecord]]:
    """Capture the records emitted inside this block.

    Purpose
    -------
    Lets :func:`pytex.app.contracts.execute` attach a call's own narration to
    that call's envelope. A :class:`~contextvars.ContextVar` rather than a
    global, so two calls served on two threads capture their own records and
    not each other's.

    Yields
    ------
    list of LogRecord
        Filled as records are emitted; ordered.

    Examples
    --------
    >>> with collecting() as records:
    ...     _ = APP_LOG.info("inside")
    >>> [record.message for record in records]
    ['inside']
    """

    collected: list[LogRecord] = []
    token = _COLLECTOR.set(collected)
    try:
        yield collected
    finally:
        _COLLECTOR.reset(token)


class ProgressReporter:
    """Turns "I am 3 steps into 40" into a throttled, ETA-bearing narration.

    Purpose
    -------
    Long simulations — a Bloch-wave CBED pattern, a grain segmentation over a
    million points — must say how far along they are. Doing that well needs
    three things every caller would otherwise re-derive: an estimate of the
    remaining time from the elapsed time and the fraction done, a throttle so a
    tight loop does not emit thousands of records, and a terminal record so the
    console shows a finished task as finished rather than as stalled at 97%.

    When to use it
    --------------
    Wrap any application-layer loop that can exceed roughly a second. Do not
    push it down into the scientific packages; pass a callback there instead if
    a library routine ever needs to report.

    Parameters
    ----------
    task : str
        Base name grouping this operation's ticks in the console. A per-instance
        counter is appended, so two runs of the same operation occupy two lines:
        a task id shared across invocations would let the second run overwrite
        the first run's entry *in place*, silently rewriting history the reader
        may still be looking at.
    total : int
        Number of steps. Must be positive.
    source : str
        Panel or operation the work belongs to.
    label : str
        Noun phrase naming the work, used in the opening and closing records —
        ``"CBED disc simulation"``.
    minimum_interval_s : float
        Shortest gap between emitted ticks.
    logbook : Logbook, optional
        Defaults to :data:`APP_LOG`.

    Examples
    --------
    >>> book = Logbook()
    >>> reporter = ProgressReporter(
    ...     "demo", total=2, source="demo", label="Demo work", minimum_interval_s=0.0, logbook=book
    ... )
    >>> reporter.advance()
    >>> reporter.finish()
    >>> book.records()[-1].level.name
    'SUCCESS'
    """

    def __init__(
        self,
        task: str,
        *,
        total: int,
        source: str = "app",
        label: str = "Work",
        minimum_interval_s: float = 0.25,
        logbook: Logbook | None = None,
    ) -> None:
        if total < 1:
            raise ValueError("ProgressReporter total must be at least 1.")
        with _RUN_COUNTER_LOCK:
            global _RUN_COUNTER
            _RUN_COUNTER += 1
            run = _RUN_COUNTER
        self._task = f"{task}#{run}"
        self._total = total
        self._source = source
        self._label = label
        self._interval = max(float(minimum_interval_s), 0.0)
        self._book = logbook if logbook is not None else APP_LOG
        self._completed = 0
        self._started = time.monotonic()
        # Seeded with the start rather than with zero: `time.monotonic()` counts
        # from an arbitrary epoch, so a zero here makes the very first step look
        # overdue and the throttle never applies to it.
        self._last_emit = self._started
        self._finished = False
        self._book.info(f"{label} started: {total} step{'' if total == 1 else 's'}.", source=source)

    @property
    def completed(self) -> int:
        """Steps advanced so far."""

        return self._completed

    def advance(self, steps: int = 1) -> None:
        """Record ``steps`` more completed steps, emitting a tick if due."""

        self._completed = min(self._completed + steps, self._total)
        now = time.monotonic()
        due = (now - self._last_emit) >= self._interval
        if not due and self._completed < self._total:
            return
        self._last_emit = now
        fraction = self._completed / self._total
        elapsed = now - self._started
        # Linear extrapolation from the work done so far. It is wrong whenever
        # the steps are not equally expensive, which is why the console labels
        # it an estimate rather than a countdown -- but a stated estimate that
        # shrinks is far more useful than a spinner that says nothing at all.
        remaining = (elapsed / fraction) - elapsed if fraction > 0.0 else None
        self._book.progress(
            self._task,
            fraction,
            eta_seconds=remaining,
            source=self._source,
            detail={"completed": self._completed, "total": self._total},
        )

    def finish(self, message: str | None = None) -> None:
        """Close the task with a success record. Idempotent."""

        if self._finished:
            return
        self._finished = True
        elapsed = time.monotonic() - self._started
        self._book.success(
            message or f"{self._label} complete in {format_duration(elapsed)}.",
            source=self._source,
            task=self._task,
            detail={"completed": self._completed, "total": self._total},
        )

    def fail(self, message: str) -> None:
        """Close the task with an error record. Idempotent."""

        if self._finished:
            return
        self._finished = True
        self._book.error(message, source=self._source, task=self._task)

    def __enter__(self) -> ProgressReporter:
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        if exc is None:
            self.finish()
        else:
            self.fail(f"{self._label} failed: {exc}")


#: stdlib severity to logbook level. ``logging`` has no notion of "important but
#: fine", so nothing maps onto NOTICE or SUCCESS; those stay reserved for code
#: that deliberately chose them.
_STDLIB_LEVELS = (
    (logging.CRITICAL, LogLevel.CRITICAL),
    (logging.ERROR, LogLevel.ERROR),
    (logging.WARNING, LogLevel.WARNING),
    (logging.INFO, LogLevel.INFO),
)


class LogbookHandler(logging.Handler):
    """Routes standard-library logging into the logbook.

    Purpose
    -------
    The server, the desktop shell and the export layer already log through
    :mod:`logging`, and a desktop user has no terminal to read it in. Bridging
    the two is what makes the console *central* rather than merely another
    channel: "port already in use" and "spot 1 selected" belong in one ordered
    stream, because to the person reading it they are one story.

    Parameters
    ----------
    logbook : Logbook, optional
        Defaults to :data:`APP_LOG`.
    exclude : sequence of str
        Logger-name prefixes to drop. The one that matters is the HTTP access
        log: every stylesheet and every icon is an ``INFO`` record, and routing
        those into a user-facing console would bury the science under a
        transcript of the transport.
    """

    def __init__(
        self, logbook: Logbook | None = None, *, exclude: Sequence[str] = ()
    ) -> None:
        super().__init__()
        self._book = logbook if logbook is not None else APP_LOG
        self._exclude = tuple(exclude)

    def emit(self, record: logging.LogRecord) -> None:
        """Translate one stdlib record. Never raises: a failing log is not news."""

        try:
            if any(record.name.startswith(prefix) for prefix in self._exclude):
                return
            level = LogLevel.DEBUG
            for threshold, mapped in _STDLIB_LEVELS:
                if record.levelno >= threshold:
                    level = mapped
                    break
            # The logger name minus the "pytex." prefix is already the subsystem
            # name the console wants to show, so nothing has to be passed twice.
            source = record.name.removeprefix("pytex.")
            self._book.emit(level, record.getMessage(), source=source)
        except Exception:  # pragma: no cover - defensive, per logging convention
            self.handleError(record)


def install_logging_bridge(
    *,
    logger_name: str = "pytex",
    level: int = logging.INFO,
    logbook: Logbook | None = None,
    exclude: Sequence[str] = (),
) -> LogbookHandler:
    """Attach a :class:`LogbookHandler` to a logger, once.

    Purpose
    -------
    Called by the shells at start-up. Idempotent by inspection rather than by a
    module-level flag, so a test that builds two servers does not end up with
    every message duplicated.

    Returns
    -------
    LogbookHandler
        The attached handler — the existing one if this logger already had one.
    """

    logger = logging.getLogger(logger_name)
    for existing in logger.handlers:
        if isinstance(existing, LogbookHandler):
            return existing
    handler = LogbookHandler(logbook, exclude=exclude)
    handler.setLevel(level)
    logger.addHandler(handler)
    if logger.level == logging.NOTSET or logger.level > level:
        logger.setLevel(level)
    return handler
