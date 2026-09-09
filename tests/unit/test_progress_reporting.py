"""Progress reporting: the mechanism, the arithmetic, and what it refuses to say.

The expectations here are properties of the contract rather than recorded
outputs: a bar must not run backwards, a rate must not be quoted before any work
is measured, a tight loop must not flood a bounded buffer, and an operation that
counted nothing must not be made to claim it finished a measurement.
"""

from __future__ import annotations

import time

import pytest

from pytex.app.contracts import execute
from pytex.app.logbook import Logbook
from pytex.app.progress import ProgressReporter
from pytex.core.progress import active_sink, report, reporting, tracking


class _Collect:
    """A sink that keeps what it was told, and nothing else."""

    def __init__(self) -> None:
        self.seen: list[tuple[float, str | None]] = []

    def report(self, fraction: float, stage: str | None = None) -> bool:
        self.seen.append((fraction, stage))
        return True


class TestTheCoreMechanism:
    def test_reporting_outside_a_block_is_a_no_op_rather_than_an_error(self) -> None:
        """Instrumented scientific code must stay callable from a test or a script."""
        assert active_sink() is None
        assert report(0.5, stage="ignored") is False
        assert list(tracking([1, 2, 3], stage="ignored")) == [1, 2, 3]

    def test_tracking_reports_work_completed_not_work_started(self) -> None:
        sink = _Collect()
        with reporting(sink):
            consumed = list(tracking("abcd", stage="Scoring"))

        assert consumed == ["a", "b", "c", "d"]
        assert [round(fraction, 6) for fraction, _ in sink.seen] == [0.0, 0.25, 0.5, 0.75, 1.0]
        assert {stage for _, stage in sink.seen} == {"Scoring"}

    def test_an_unsized_iterable_without_a_total_is_not_measured(self) -> None:
        """A fraction of an unknown quantity is not a measurement."""
        sink = _Collect()
        with reporting(sink):
            consumed = list(tracking((value for value in range(3)), stage="Unknown"))

        assert consumed == [0, 1, 2]
        assert sink.seen == []

    def test_a_total_makes_an_unsized_iterable_measurable(self) -> None:
        sink = _Collect()
        with reporting(sink):
            list(tracking((value for value in range(4)), total=4, stage="Known"))

        assert [round(fraction, 6) for fraction, _ in sink.seen] == [0.0, 0.25, 0.5, 0.75, 1.0]

    def test_the_sink_is_restored_after_the_block(self) -> None:
        outer, inner = _Collect(), _Collect()
        with reporting(outer):
            with reporting(inner):
                report(0.5)
            report(0.25)

        assert inner.seen == [(0.5, None)]
        assert outer.seen == [(0.25, None)]
        assert active_sink() is None

    def test_reporting_none_deliberately_silences_a_block(self) -> None:
        sink = _Collect()
        with reporting(sink), reporting(None):
            report(0.5)

        assert sink.seen == []


class TestTheWorkbenchReporter:
    def _reporter(self, **kwargs: object) -> tuple[ProgressReporter, Logbook]:
        book = Logbook(capacity=200)
        return ProgressReporter("op.id", "An operation", book=book, **kwargs), book  # type: ignore[arg-type]

    def test_the_bar_never_runs_backwards(self) -> None:
        reporter, _ = self._reporter(minimum_interval_s=0.0)
        reporter.report(0.6)
        reporter.report(0.2)

        assert reporter.fraction == pytest.approx(0.6)

    def test_a_fraction_outside_the_unit_interval_is_clamped(self) -> None:
        reporter, _ = self._reporter(minimum_interval_s=0.0)
        reporter.report(-1.0)
        assert reporter.fraction == 0.0
        reporter.report(4.0)
        assert reporter.fraction == 1.0

    def test_no_remaining_time_is_quoted_before_any_work_is_measured(self) -> None:
        """A rate cannot be measured from no completed work."""
        reporter, _ = self._reporter()
        assert reporter.estimated_remaining_s() is None

    def test_no_remaining_time_is_quoted_at_completion(self) -> None:
        reporter, _ = self._reporter(minimum_interval_s=0.0)
        reporter.report(1.0)
        assert reporter.estimated_remaining_s() is None

    def test_the_remaining_time_extrapolates_the_measured_rate(self) -> None:
        """At a quarter done, three quarters remain, so about three times elapsed."""
        reporter, _ = self._reporter(minimum_interval_s=0.0)
        time.sleep(0.05)
        reporter.report(0.25)
        remaining = reporter.estimated_remaining_s()

        assert remaining is not None
        assert remaining == pytest.approx(reporter.elapsed_s * 3.0, rel=0.2)

    def test_a_tight_loop_is_thinned_rather_than_flooding_the_buffer(self) -> None:
        reporter, book = self._reporter(minimum_interval_s=10.0)
        for index in range(1000):
            reporter.report(index / 1000.0)

        # The first tick always goes out; the rest are inside the interval.
        assert reporter.emitted == 1
        assert len(book.records()) == 1

    def test_a_change_of_stage_always_reports_whatever_the_rate_limit(self) -> None:
        """The stage is the part that tells the user where the work has got to."""
        reporter, _ = self._reporter(minimum_interval_s=10.0)
        reporter.report(0.1, stage="First")
        reporter.report(0.2, stage="First")
        reporter.report(0.3, stage="Second")

        assert reporter.emitted == 2

    def test_an_emitted_tick_carries_the_fraction_stage_and_elapsed_time(self) -> None:
        reporter, book = self._reporter(minimum_interval_s=0.0)
        reporter.report(0.5, stage="Halfway")
        record = book.records()[-1]

        assert record.level.token == "progress"
        assert record.task == "op.id"
        assert record.progress == pytest.approx(0.5)
        assert record.detail["stage"] == "Halfway"
        assert record.detail["title"] == "An operation"
        assert record.detail["elapsed_s"] >= 0.0
        assert "50% complete" in record.message

    def test_finishing_closes_a_bar_that_was_left_part_way(self) -> None:
        reporter, book = self._reporter(minimum_interval_s=0.0)
        reporter.report(0.97)
        reporter.finish()

        assert book.records()[-1].progress == pytest.approx(1.0)

    def test_finishing_claims_nothing_for_an_operation_that_counted_nothing(self) -> None:
        """A lone 100% tick would report a measurement that was never taken."""
        reporter, book = self._reporter()
        reporter.finish()

        assert reporter.emitted == 0
        assert book.records() == ()


class TestOperationsReportThroughTheDispatchPath:
    def test_a_long_operation_narrates_its_own_progress(self) -> None:
        """Phase identification scores candidates one at a time and says so."""
        envelope, status = execute(
            "xrd.phase_identification",
            {
                "candidates": {
                    "phases": [
                        {"phase": {"builtin": "ni_fcc"}},
                        {"phase": {"builtin": "cu_fcc"}},
                        {"phase": {"builtin": "fe_bcc"}},
                    ]
                },
            },
        )
        assert status == 200, envelope

        ticks = [record for record in envelope["log"] if record["level"] == "progress"]
        assert ticks, "the operation reported no progress at all"
        fractions = [tick["progress"] for tick in ticks]
        assert fractions == sorted(fractions), "the reported fraction went backwards"
        assert fractions[-1] == pytest.approx(1.0)
        assert all(tick["task"] == "xrd.phase_identification" for tick in ticks)

    def test_an_operation_that_counts_nothing_reports_nothing(self) -> None:
        """Silence is the honest answer, and the shell shows elapsed time instead."""
        envelope, status = execute("calc.catalog", {})
        assert status == 200

        assert [record for record in envelope["log"] if record["level"] == "progress"] == []
