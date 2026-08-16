"""The centralized logbook: ordering, severity, capture, and the call envelope.

What these tests protect is the property the console depends on: that the log is
*one* ordered stream with exact positions in it. A record delivered twice, a
sequence number reused after a clear, or a call whose narration is attached to a
different call would each make the console quietly wrong rather than visibly
broken — which is the failure mode worth a test.
"""

from __future__ import annotations

import logging
import threading

import pytest

from pytex.app.contracts import execute
from pytex.app.errors import InvalidInputError
from pytex.app.logbook import (
    APP_LOG,
    Logbook,
    LogbookHandler,
    LogLevel,
    ProgressReporter,
    collecting,
    format_duration,
    install_logging_bridge,
)
from pytex.app.registry import NumberParameter, ServiceRegistry


class TestLogLevel:
    """Severity is a shared vocabulary, not a per-module opinion."""

    def test_levels_are_ordered_from_progress_to_critical(self) -> None:
        ordered = sorted(LogLevel, key=lambda level: level.severity)
        assert [level.token for level in ordered] == [
            "progress",
            "debug",
            "info",
            "notice",
            "success",
            "warning",
            "error",
            "critical",
        ]

    def test_progress_is_the_first_thing_a_filter_drops(self) -> None:
        """A tight loop's ticks must never outrank a warning in the console."""

        assert LogLevel.PROGRESS.severity < LogLevel.INFO.severity

    def test_an_unknown_token_is_refused_by_name(self) -> None:
        with pytest.raises(ValueError, match="Unknown log level 'verbose'"):
            LogLevel.from_token("verbose")

    def test_every_token_round_trips(self) -> None:
        for level in LogLevel:
            assert LogLevel.from_token(level.token) is level


class TestFormatDuration:
    """The ETA half of a progress message, phrased one way everywhere."""

    @pytest.mark.parametrize(
        ("seconds", "expected"),
        [
            (0.0, "0 sec"),
            (9.4, "9 sec"),
            (59.0, "59 sec"),
            (60.0, "1 min 00 sec"),
            (150.0, "2 min 30 sec"),
            (3930.0, "1 hr 05 min"),
        ],
    )
    def test_durations_read_the_way_a_person_says_them(
        self, seconds: float, expected: str
    ) -> None:
        assert format_duration(seconds) == expected

    @pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
    def test_a_nonsense_duration_is_admitted_rather_than_rendered(self, value: float) -> None:
        """A wrong countdown on screen is worse than an absent one."""

        assert format_duration(value) == "unknown"


class TestLogbook:
    """The buffer itself."""

    def test_records_arrive_in_order_with_consecutive_sequences(self) -> None:
        book = Logbook()
        book.info("first")
        book.warning("second")
        sequences = [record.sequence for record in book.records()]
        assert sequences == [1, 2]
        assert [record.message for record in book.records()] == ["first", "second"]

    def test_the_buffer_is_bounded_and_drops_the_oldest(self) -> None:
        book = Logbook(capacity=3)
        for index in range(5):
            book.info(f"message {index}")
        assert [record.message for record in book.records()] == [
            "message 2",
            "message 3",
            "message 4",
        ]

    def test_a_zero_capacity_book_is_refused(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            Logbook(capacity=0)

    def test_since_returns_only_what_the_caller_has_not_seen(self) -> None:
        book = Logbook()
        book.info("one")
        book.info("two")
        cursor = book.latest_sequence()
        book.info("three")
        assert [record.message for record in book.records(since=cursor)] == ["three"]

    def test_clearing_does_not_reuse_a_sequence_number(self) -> None:
        """A client holding ``since=2`` must never be handed a different 3."""

        book = Logbook()
        book.info("one")
        book.info("two")
        book.clear()
        book.info("three")
        assert book.records()[0].sequence == 3

    def test_a_minimum_level_filters_by_severity(self) -> None:
        book = Logbook()
        book.debug("noise")
        book.progress("task", 0.5)
        book.warning("look at this")
        book.error("and this")
        selected = book.records(minimum_level=LogLevel.WARNING)
        assert [record.message for record in selected] == ["look at this", "and this"]

    def test_every_level_has_a_convenience_method(self) -> None:
        book = Logbook()
        for name in ("debug", "info", "notice", "success", "warning", "error", "critical"):
            getattr(book, name)(f"a {name} record")
        assert [record.level.token for record in book.records()] == [
            "debug",
            "info",
            "notice",
            "success",
            "warning",
            "error",
            "critical",
        ]

    def test_concurrent_writers_lose_nothing_and_reuse_no_sequence(self) -> None:
        """The HTTP server answers on many threads; the log is shared state."""

        book = Logbook(capacity=1000)

        def write() -> None:
            for index in range(50):
                book.info(f"message {index}")

        threads = [threading.Thread(target=write) for _ in range(8)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        records = book.records()
        assert len(records) == 400
        assert len({record.sequence for record in records}) == 400


class TestProgressRecords:
    """The "50% progress. ETA: 2 min 30 sec" surface."""

    def test_a_tick_writes_the_sentence_the_user_reads(self) -> None:
        book = Logbook()
        book.progress("cbed", 0.5, eta_seconds=150.0)
        record = book.records()[-1]
        assert record.message == "50% progress. ETA: 2 min 30 sec"
        assert record.level is LogLevel.PROGRESS
        assert record.progress == pytest.approx(0.5)
        assert record.task == "cbed"

    def test_a_tick_without_an_estimate_claims_none(self) -> None:
        book = Logbook()
        book.progress("cbed", 0.25)
        assert book.records()[-1].message == "25% progress."

    @pytest.mark.parametrize(("given", "expected"), [(-0.5, 0.0), (1.5, 1.0)])
    def test_the_fraction_is_clamped_to_its_own_ends(
        self, given: float, expected: float
    ) -> None:
        """A bar past 100% is a defect the user should not be shown."""

        book = Logbook()
        book.progress("task", given)
        assert book.records()[-1].progress == pytest.approx(expected)

    def test_ticks_share_a_task_so_the_console_can_replace_them(self) -> None:
        book = Logbook()
        book.progress("segmentation", 0.2)
        book.progress("segmentation", 0.9)
        assert {record.task for record in book.records()} == {"segmentation"}

    def test_two_runs_of_one_reporter_take_two_console_lines(self) -> None:
        """A shared task id would let the second run rewrite the first in place."""

        book = Logbook()
        first = ProgressReporter("raster", total=1, minimum_interval_s=0.0, logbook=book)
        second = ProgressReporter("raster", total=1, minimum_interval_s=0.0, logbook=book)
        first.advance()
        second.advance()
        tasks = {record.task for record in book.records() if record.task is not None}
        assert len(tasks) == 2
        assert all(task.startswith("raster#") for task in tasks)


class TestProgressReporter:
    """Throttled narration of a long loop."""

    def test_it_opens_with_a_count_and_closes_with_a_success(self) -> None:
        book = Logbook()
        reporter = ProgressReporter(
            "demo", total=4, source="demo", label="Demo work", minimum_interval_s=0.0, logbook=book
        )
        for _ in range(4):
            reporter.advance()
        reporter.finish()
        records = book.records()
        assert records[0].message == "Demo work started: 4 steps."
        assert records[0].level is LogLevel.INFO
        assert records[-1].level is LogLevel.SUCCESS
        assert reporter.completed == 4

    def test_the_final_step_always_emits_even_when_throttled(self) -> None:
        """A task that stops reporting at 97% looks stalled rather than done."""

        book = Logbook()
        reporter = ProgressReporter(
            "demo", total=2, source="demo", minimum_interval_s=3600.0, logbook=book
        )
        reporter.advance()
        reporter.advance()
        ticks = [record for record in book.records() if record.level is LogLevel.PROGRESS]
        assert len(ticks) == 1
        assert ticks[-1].progress == pytest.approx(1.0)

    def test_finishing_twice_reports_once(self) -> None:
        book = Logbook()
        reporter = ProgressReporter("demo", total=1, minimum_interval_s=0.0, logbook=book)
        reporter.finish()
        reporter.finish()
        assert sum(record.level is LogLevel.SUCCESS for record in book.records()) == 1

    def test_used_as_a_context_manager_a_failure_closes_the_task(self) -> None:
        book = Logbook()
        with pytest.raises(RuntimeError):
            with ProgressReporter(
                "demo", total=1, label="Demo work", minimum_interval_s=0.0, logbook=book
            ):
                raise RuntimeError("the solver diverged")
        last = book.records()[-1]
        assert last.level is LogLevel.ERROR
        assert "the solver diverged" in last.message

    def test_a_zero_step_task_is_refused(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            ProgressReporter("demo", total=0)


class TestWireForm:
    """What the console actually receives."""

    def test_a_record_carries_its_severity_and_source(self) -> None:
        book = Logbook()
        book.warning("Careful", source="ebsd", detail={"field": "step"})
        payload = book.records()[-1].to_json()
        assert payload["level"] == "warning"
        assert payload["source"] == "ebsd"
        assert payload["detail"] == {"field": "step"}
        assert payload["schema"] == "pytex.log_record/1"

    def test_absent_optional_fields_are_omitted_rather_than_null(self) -> None:
        book = Logbook()
        book.info("plain")
        payload = book.records()[-1].to_json()
        assert "progress" not in payload
        assert "task" not in payload
        assert "detail" not in payload

    def test_describe_renders_one_readable_line(self) -> None:
        book = Logbook()
        book.error("Only integers are allowed", source="controls", detail={"field": "zone"})
        line = book.records()[-1].describe()
        assert "ERROR" in line
        assert "[controls]" in line
        assert "Only integers are allowed" in line
        assert "field=zone" in line


class TestCollecting:
    """Capturing one call's narration without capturing another's."""

    def test_records_emitted_inside_the_block_are_captured(self) -> None:
        with collecting() as captured:
            APP_LOG.info("inside the block")
        assert [record.message for record in captured] == ["inside the block"]

    def test_records_emitted_outside_the_block_are_not(self) -> None:
        APP_LOG.info("before")
        with collecting() as captured:
            pass
        APP_LOG.info("after")
        assert captured == []

    def test_two_threads_capture_only_their_own_records(self) -> None:
        """Two colleagues calculating at once must not read each other's log."""

        captured: dict[str, list[str]] = {}
        barrier = threading.Barrier(2)

        def run(name: str) -> None:
            with collecting() as records:
                barrier.wait()
                APP_LOG.info(f"from {name}")
                barrier.wait()
            captured[name] = [record.message for record in records]

        threads = [threading.Thread(target=run, args=(name,)) for name in ("left", "right")]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        assert captured["left"] == ["from left"]
        assert captured["right"] == ["from right"]


class TestLoggingBridge:
    """Standard-library logging joins the same stream."""

    def test_stdlib_severities_map_onto_logbook_levels(self) -> None:
        book = Logbook()
        logger = logging.getLogger("pytex.test.bridge.levels")
        logger.propagate = False
        logger.setLevel(logging.DEBUG)
        handler = LogbookHandler(book)
        logger.addHandler(handler)
        try:
            logger.debug("d")
            logger.info("i")
            logger.warning("w")
            logger.error("e")
            logger.critical("c")
        finally:
            logger.removeHandler(handler)
        assert [record.level.token for record in book.records()] == [
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ]

    def test_the_logger_name_becomes_the_console_source(self) -> None:
        book = Logbook()
        logger = logging.getLogger("pytex.app.export")
        logger.propagate = False
        handler = LogbookHandler(book)
        logger.addHandler(handler)
        try:
            logger.warning("Nothing to export")
        finally:
            logger.removeHandler(handler)
        assert book.records()[-1].source == "app.export"

    def test_an_excluded_logger_is_dropped(self) -> None:
        """The HTTP access log must not bury the science in the console."""

        book = Logbook()
        logger = logging.getLogger("pytex.app.server.requests")
        logger.propagate = False
        handler = LogbookHandler(book, exclude=("pytex.app.server.requests",))
        logger.addHandler(handler)
        try:
            logger.info("GET /app.css")
        finally:
            logger.removeHandler(handler)
        assert book.records() == ()

    def test_installing_the_bridge_twice_does_not_duplicate_messages(self) -> None:
        logger_name = "pytex.test.bridge.idempotent"
        first = install_logging_bridge(logger_name=logger_name)
        second = install_logging_bridge(logger_name=logger_name)
        try:
            assert first is second
            handlers = [
                handler
                for handler in logging.getLogger(logger_name).handlers
                if isinstance(handler, LogbookHandler)
            ]
            assert len(handlers) == 1
        finally:
            logging.getLogger(logger_name).removeHandler(first)


class TestEnvelopeNarration:
    """Every operation of every module narrates itself, without opting in."""

    @staticmethod
    def _registry() -> ServiceRegistry:
        registry = ServiceRegistry()

        @registry.operation(
            "demo.double",
            title="Double a number",
            summary="Return twice the input.",
            help_text="Exists to prove that dispatch narrates itself.",
            parameters=(
                NumberParameter("value", label="Value", help_text="Any real number."),
            ),
            returns="An object with the doubled value.",
        )
        def _double(request: dict[str, object]) -> dict[str, object]:
            return {"value": float(request["value"]) * 2.0}  # type: ignore[arg-type]

        return registry

    def test_a_successful_call_reports_its_title_and_duration(self) -> None:
        envelope, status = execute("demo.double", {"value": 3}, registry=self._registry())
        assert status == 200
        messages = [record["message"] for record in envelope["log"]]
        assert messages[0] == "Double a number started."
        assert messages[-1].startswith("Double a number completed in ")
        assert envelope["log"][-1]["level"] == "success"
        assert envelope["log"][-1]["source"] == "demo.double"

    def test_a_rejected_input_is_logged_verbatim_beside_its_field(self) -> None:
        """The console and the message beside the control must not differ."""

        envelope, status = execute("demo.double", {"value": "abc"}, registry=self._registry())
        assert status == InvalidInputError("x").status
        failure = envelope["log"][-1]
        assert failure["level"] == "error"
        assert failure["message"] == envelope["error"]["message"]
        assert failure["detail"]["field"] == "value"

    def test_an_unknown_operation_is_named_by_its_id(self) -> None:
        """There is no registered title to use, and inventing one would lie."""

        envelope, _status = execute("demo.missing", {}, registry=self._registry())
        assert envelope["log"][0]["message"] == "demo.missing started."

    def test_the_envelope_always_carries_a_log_field(self) -> None:
        envelope, _status = execute("demo.double", {"value": 1}, registry=self._registry())
        assert isinstance(envelope["log"], list)
        assert all("sequence" in record for record in envelope["log"])
