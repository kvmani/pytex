"""The feedback surface: configuration, the JSON store, the relay, the route.

The contract these tests hold is the one the module docstring states: a
submission that reached the server is on disk *before* anything is said about
e-mail, and a relay that is unreachable costs a notification rather than a
note.
"""

from __future__ import annotations

import json
import smtplib
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

from pytex.app.config import (
    AppConfig,
    ConfigError,
    FeedbackConfig,
    RelayConfig,
    TourConfig,
    clear_config_cache,
    load_app_config,
)
from pytex.app.errors import InvalidInputError
from pytex.app.feedback import (
    CATEGORIES,
    FeedbackSubmission,
    append_to_store,
    normalise_submission,
    read_store,
    record_feedback,
    relay_feedback,
)


@pytest.fixture(autouse=True)
def _forget_configuration() -> Any:
    """Each test reads its own configuration, never one another's."""

    clear_config_cache()
    yield
    clear_config_cache()


def _config(tmp_path: Path, **relay: Any) -> AppConfig:
    return AppConfig(
        feedback=FeedbackConfig(store_path=tmp_path / "feedback.json"),
        relay=RelayConfig(**relay),
        tour=TourConfig(),
    )


class TestConfiguration:
    def test_defaults_are_a_working_deployment(self) -> None:
        """No file at all must leave a usable application, not a broken one."""

        config = AppConfig()
        assert config.feedback.enabled is True
        assert config.tour.enabled is True
        assert config.relay.enabled is False
        assert config.relay.recipients == ("kvmani@barc.gov.in",)

    def test_a_file_overlays_only_what_it_names(self, tmp_path: Path) -> None:
        path = tmp_path / "pytex_app.yml"
        path.write_text(
            "tour:\n  enabled: false\nrelay:\n  enabled: true\n  host: relay.invalid\n",
            encoding="utf-8",
        )
        config = load_app_config(path)
        assert config.tour.enabled is False
        assert config.relay.host == "relay.invalid"
        # Untouched sections keep their defaults rather than being emptied.
        assert config.feedback.enabled is True
        assert config.relay.subject_prefix == "[PyTex feedback]"
        assert config.source == path

    def test_the_shipped_example_is_a_file_this_module_can_read(self) -> None:
        """The template must stay in step with the loader that reads it.

        This is the failure the template exists to prevent: a documented key
        that the code renamed, discovered by an operator at deployment time
        rather than here.
        """

        example = Path("config/pytex_app.example.yml")
        assert example.is_file()
        config = load_app_config(example)
        assert config.relay.recipients == ("kvmani@barc.gov.in",)
        assert config.feedback.enabled is True
        assert config.tour.enabled is True
        # The example must not ship as a live relay, whatever else it says.
        assert config.relay.enabled is False

    def test_an_unknown_key_is_an_error_rather_than_a_shrug(self, tmp_path: Path) -> None:
        path = tmp_path / "pytex_app.yml"
        path.write_text("relay:\n  smtp_host: relay.invalid\n", encoding="utf-8")
        with pytest.raises(ConfigError, match="smtp_host"):
            load_app_config(path)

    def test_an_unknown_section_names_the_sections_that_exist(self, tmp_path: Path) -> None:
        path = tmp_path / "pytex_app.yml"
        path.write_text("email:\n  enabled: true\n", encoding="utf-8")
        with pytest.raises(ConfigError, match="feedback, relay and tour"):
            load_app_config(path)

    def test_credentials_come_from_the_environment_first(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        relay = RelayConfig(username="in-file", password="in-file")
        monkeypatch.setenv("PYTEX_SMTP_USERNAME", "from-env")
        monkeypatch.setenv("PYTEX_SMTP_PASSWORD", "secret")
        assert relay.resolved_credentials() == ("from-env", "secret")

    def test_no_credentials_at_all_is_the_unauthenticated_relay(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("PYTEX_SMTP_USERNAME", raising=False)
        monkeypatch.delenv("PYTEX_SMTP_PASSWORD", raising=False)
        assert RelayConfig().resolved_credentials() == ("", "")


class TestSubmission:
    def test_a_message_is_the_only_thing_required(self) -> None:
        submission = normalise_submission(
            {"message": "The IPF legend needs a key."},
            config=FeedbackConfig(),
            environment={"pytex_version": "0.0.0"},
        )
        assert submission.message == "The IPF legend needs a key."
        assert submission.category == "feedback"
        assert submission.name == ""

    def test_an_empty_message_is_refused_with_encouragement(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            normalise_submission({"message": "   "}, config=FeedbackConfig())
        assert excinfo.value.details["field"] == "message"
        assert excinfo.value.hint is not None

    def test_the_clock_is_the_servers(self) -> None:
        """A workstation with a wrong clock must not file a note under next year."""

        submission = normalise_submission(
            {"message": "hello", "received_at": "1999-01-01T00:00:00Z"},
            config=FeedbackConfig(),
            environment={},
            now=datetime(2026, 8, 21, 9, 0, tzinfo=UTC),
        )
        assert submission.received_at == "2026-08-21T09:00:00Z"

    def test_an_unknown_category_lists_the_known_ones(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            normalise_submission(
                {"message": "hello", "category": "complaint"}, config=FeedbackConfig()
            )
        assert "feature" in (excinfo.value.hint or "")

    @pytest.mark.parametrize("rating", [0, 6, "many"])
    def test_a_rating_outside_one_to_five_is_refused(self, rating: object) -> None:
        with pytest.raises(InvalidInputError):
            normalise_submission({"message": "hello", "rating": rating}, config=FeedbackConfig())

    def test_an_oversized_message_says_what_the_limit_is(self) -> None:
        config = FeedbackConfig(max_message_characters=20)
        with pytest.raises(InvalidInputError, match="limit is 20"):
            normalise_submission({"message": "x" * 21}, config=config)

    def test_the_prose_names_every_field_including_the_blank_ones(self) -> None:
        """So a reader can tell "no name given" from "the field was lost"."""

        submission = normalise_submission(
            {"message": "A note.", "email": "someone@example.invalid"},
            config=FeedbackConfig(),
            environment={"pytex_version": "0.0.0-under-test", "shell": "web"},
        )
        prose = submission.describe()
        assert "(not given)" in prose
        assert "someone@example.invalid" in prose
        assert "shell" in prose

    def test_the_subject_carries_the_category_and_the_first_line(self) -> None:
        submission = FeedbackSubmission(
            received_at="2026-08-21T09:00:00Z",
            category="feature",
            message="Add a Kearns factor panel.\nWith the three directions.",
            name="A Researcher",
        )
        subject = submission.subject("[PyTex feedback]")
        assert subject.startswith("[PyTex feedback] Feature request from A Researcher:")
        assert "Add a Kearns factor panel." in subject
        assert "three directions" not in subject


class TestStore:
    def test_the_store_is_a_json_array_that_grows(self, tmp_path: Path) -> None:
        path = tmp_path / "nested" / "feedback.json"
        for index in range(3):
            submission = normalise_submission(
                {"message": f"note {index}"}, config=FeedbackConfig(), environment={}
            )
            append_to_store(submission, path=path)
        document = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(document, list)
        assert [entry["message"] for entry in document] == ["note 0", "note 1", "note 2"]
        assert read_store(path) == document

    def test_a_missing_store_reads_as_empty_rather_than_failing(self, tmp_path: Path) -> None:
        assert read_store(tmp_path / "absent.json") == []

    def test_a_file_that_is_not_a_feedback_store_is_never_overwritten(self, tmp_path: Path) -> None:
        """Refusing is the only safe answer: the file belongs to someone."""

        path = tmp_path / "feedback.json"
        path.write_text('{"not": "an array"}', encoding="utf-8")
        with pytest.raises(ValueError, match="not a JSON array"):
            read_store(path)


class TestRelay:
    def test_a_deployment_with_no_relay_says_so_without_pretending(self, tmp_path: Path) -> None:
        submission = normalise_submission({"message": "hi"}, config=FeedbackConfig())
        delivered, detail = relay_feedback(submission, relay=RelayConfig(enabled=False))
        assert delivered is False
        assert "filed locally" in detail

    def test_an_enabled_relay_sends_one_message_to_the_configured_recipients(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        sent: list[Any] = []

        class _FakeSMTP:
            def __init__(self, host: str, port: int, timeout: float) -> None:
                sent.append(("connect", host, port, timeout))

            def __enter__(self) -> _FakeSMTP:
                return self

            def __exit__(self, *_: object) -> bool:
                return False

            def starttls(self) -> None:
                sent.append(("starttls",))

            def login(self, username: str, password: str) -> None:
                sent.append(("login", username, password))

            def send_message(self, message: Any) -> None:
                sent.append(("send", message))

        monkeypatch.setattr(smtplib, "SMTP", _FakeSMTP)
        monkeypatch.delenv("PYTEX_SMTP_USERNAME", raising=False)
        monkeypatch.delenv("PYTEX_SMTP_PASSWORD", raising=False)
        relay = RelayConfig(
            enabled=True,
            host="relay.invalid",
            port=25,
            recipients=("kvmani@barc.gov.in",),
            from_address="pytex-noreply@example.invalid",
        )
        submission = normalise_submission(
            {"message": "A note.", "email": "author@example.invalid"},
            config=FeedbackConfig(),
            environment={},
        )
        delivered, detail = relay_feedback(submission, relay=relay)
        assert delivered is True
        assert "e-mailed" in detail
        assert sent[0] == ("connect", "relay.invalid", 25, 10.0)
        message = next(entry[1] for entry in sent if entry[0] == "send")
        assert message["To"] == "kvmani@barc.gov.in"
        assert message["From"] == "pytex-noreply@example.invalid"
        # So a maintainer can answer with Reply rather than by copying an
        # address out of the body.
        assert message["Reply-To"] == "author@example.invalid"
        assert "A note." in message.get_content()
        # An unauthenticated relay must not be sent a login.
        assert not any(entry[0] == "login" for entry in sent)

    def test_an_unreachable_relay_is_reported_without_naming_the_host(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The submitter is not the administrator; the log gets the detail."""

        def _refuse(*_: object, **__: object) -> None:
            raise OSError("connection refused")

        monkeypatch.setattr(smtplib, "SMTP", _refuse)
        relay = RelayConfig(enabled=True, host="relay.invalid", recipients=("a@b.invalid",))
        submission = normalise_submission({"message": "hi"}, config=FeedbackConfig())
        delivered, detail = relay_feedback(submission, relay=relay)
        assert delivered is False
        assert "relay.invalid" not in detail
        assert "filed locally" in detail


class TestRecordFeedback:
    def test_a_note_is_stored_even_when_the_relay_is_down(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The ordering the whole module exists to guarantee."""

        def _refuse(*_: object, **__: object) -> None:
            raise OSError("connection refused")

        monkeypatch.setattr(smtplib, "SMTP", _refuse)
        config = _config(tmp_path, enabled=True, host="relay.invalid", recipients=("a@b.invalid",))
        _submission, receipt = record_feedback({"message": "Still worth saying."}, config=config)
        assert receipt.stored is True
        assert receipt.delivered is False
        stored = read_store(config.feedback.store_path)
        assert stored[0]["message"] == "Still worth saying."

    def test_the_server_metadata_travels_with_the_note(self, tmp_path: Path) -> None:
        config = _config(tmp_path)
        _submission, receipt = record_feedback(
            {"message": "hello", "context": {"panel": "texture"}},
            config=config,
            environment={"pytex_version": "0.0.0-under-test", "shell": "desktop"},
        )
        assert receipt.stored is True
        record = read_store(config.feedback.store_path)[0]
        assert record["environment"]["shell"] == "desktop"
        # Client-supplied context is kept, and kept separate from what the
        # server knows, because only one of the two can be trusted.
        assert record["context"] == {"panel": "texture"}

    def test_an_unwritable_store_is_reported_rather_than_swallowed(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Nothing holds the note, so the user must not be told it was received."""

        config = _config(tmp_path)

        def _fail(*_: object, **__: object) -> None:
            raise OSError("read-only file system")

        monkeypatch.setattr("pytex.app.feedback.append_to_store", _fail)
        _submission, receipt = record_feedback({"message": "hello"}, config=config)
        assert receipt.stored is False
        assert receipt.delivered is False
        assert "could not write" in receipt.delivery_detail

    def test_every_category_is_accepted(self, tmp_path: Path) -> None:
        config = _config(tmp_path)
        for value, _label in CATEGORIES:
            _submission, receipt = record_feedback(
                {"message": "a note", "category": value}, config=config
            )
            assert receipt.stored is True
        assert len(read_store(config.feedback.store_path)) == len(CATEGORIES)
