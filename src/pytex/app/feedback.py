"""Feedback and feature requests: recorded first, delivered second.

Purpose
-------
Somebody using the workbench notices something — a control that reads wrongly,
a quantity that is missing, a paper whose method they wish were here. This
module is what happens to that note.

Two things happen to every submission, in this order and never the other way
round:

1. **It is appended to a JSON file.** This is not optional and does not depend
   on any network. A submission that reached the server is a submission that is
   on disk before the user is told anything.
2. **It is relayed by e-mail**, if — and only if — the deployment configured an
   internal SMTP relay. Delivery is best-effort by construction: the store is
   already written, so a relay that is down costs a notification, not a note.

That ordering is the whole design. The alternative — mail it, and store it if
the mail failed — loses submissions exactly when the site is having a bad day,
which is when people have the most to say.

When and where to use it
------------------------
:func:`record_feedback` is called by the ``POST /api/feedback`` route in
:mod:`pytex.app.server`. Call it directly to replay a submission, or from a
script that migrates a store.

Expected inputs
---------------
A submission mapping as the browser sends it (see :func:`normalise_submission`
for the accepted fields) plus the metadata the *server* knows and the client
must not be trusted to state: the running version, which shell answered, and
when it arrived.

Expected outputs
----------------
:class:`FeedbackReceipt`, which says where it was stored and what became of the
delivery attempt. The receipt is what the browser shows the user, so it never
carries a relay hostname or a credential.
"""

from __future__ import annotations

import getpass
import json
import logging
import platform
import smtplib
import socket
import threading
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from email.message import EmailMessage
from pathlib import Path
from typing import Any

from pytex.app.config import AppConfig, FeedbackConfig, RelayConfig, load_app_config
from pytex.app.errors import InvalidInputError

__all__ = [
    "CATEGORIES",
    "FeedbackReceipt",
    "FeedbackSubmission",
    "append_to_store",
    "normalise_submission",
    "read_store",
    "record_feedback",
    "relay_feedback",
    "server_environment",
]

_LOGGER = logging.getLogger("pytex.app.feedback")

#: What a submission can be about. Kept short deliberately: a long taxonomy
#: makes the submitter classify their thought before they have written it down,
#: and the useful sorting is done by reading, not by the dropdown.
CATEGORIES: tuple[tuple[str, str], ...] = (
    ("feedback", "General feedback"),
    ("feature", "Feature request"),
    ("problem", "Something went wrong"),
    ("science", "A scientific question or correction"),
)

#: The store is appended to under this lock. One process serves both shells, so
#: a lock is enough; the file is rewritten atomically besides, so a second
#: process cannot leave a half-written array behind either.
_STORE_LOCK = threading.Lock()

#: Longest value accepted for a short free-text field (name, e-mail, group).
_SHORT_FIELD_LIMIT = 200


@dataclass(frozen=True)
class FeedbackSubmission:
    """One note, as it will be stored.

    Attributes
    ----------
    received_at : str
        UTC, ISO 8601 with a ``Z`` suffix. Set by the server, not the client:
        a workstation with a wrong clock must not be able to file a submission
        under next year.
    category : str
        One of :data:`CATEGORIES`.
    message : str
        What the user wrote. The only field that is mandatory.
    name, email, organisation : str
        Who they are, all optional. A note worth reading is worth reading
        anonymously, and demanding a name suppresses exactly the frank
        criticism that is most useful.
    rating : int or None
        1-5, optional. Coarse on purpose; it sorts submissions, it does not
        measure anything.
    contact_consent : bool
        Whether they are happy to be written back to.
    context : dict
        Where in the application they were: panel, operation, and the page's
        own reading of the shell. Client-supplied, so it is recorded as
        *claimed* context and never used for anything but reading.
    environment : dict
        What the server knows: PyTex version, shell, Python, platform, host,
        and the account the server runs as. Server-supplied, so this is the
        part that can be trusted.
    """

    received_at: str
    category: str
    message: str
    name: str = ""
    email: str = ""
    organisation: str = ""
    rating: int | None = None
    contact_consent: bool = False
    context: dict[str, Any] = field(default_factory=dict)
    environment: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        """The record as it is written to the store."""

        return asdict(self)

    def describe(self) -> str:
        """The submission as prose, which is also the body of the e-mail.

        Returns
        -------
        str
            A plain-text rendering with every field named, including the ones
            left blank, so that a reader of the mailbox can tell "no name
            given" from "the name field was dropped in transit".

        Examples
        --------
        >>> submission = FeedbackSubmission(
        ...     received_at="2026-08-21T09:00:00Z",
        ...     category="feature",
        ...     message="An inverse pole figure of the sample normal, please.",
        ...     name="A Researcher",
        ... )
        >>> print(submission.describe().splitlines()[0])
        PyTex feedback - Feature request
        """

        label = dict(CATEGORIES).get(self.category, self.category)
        lines = [
            f"PyTex feedback - {label}",
            "",
            self.message.strip(),
            "",
            "-- who --",
            f"Name         : {self.name or '(not given)'}",
            f"E-mail       : {self.email or '(not given)'}",
            f"Group        : {self.organisation or '(not given)'}",
            f"Rating       : {self.rating if self.rating is not None else '(not given)'}",
            f"May we reply : {'yes' if self.contact_consent else 'not stated'}",
            "",
            "-- where --",
            f"Received at  : {self.received_at}",
        ]
        for key, value in sorted(self.environment.items()):
            lines.append(f"{key:<13}: {value}")
        for key, value in sorted(self.context.items()):
            lines.append(f"{key:<13}: {value} (as reported by the page)")
        return "\n".join(lines)

    def subject(self, prefix: str) -> str:
        """The e-mail subject: the prefix, the category, and the first line."""

        label = dict(CATEGORIES).get(self.category, self.category)
        first = self.message.strip().splitlines()[0] if self.message.strip() else ""
        summary = first[:80] + ("…" if len(first) > 80 else "")
        who = f" from {self.name}" if self.name else ""
        return f"{prefix} {label}{who}: {summary}".strip()


@dataclass(frozen=True)
class FeedbackReceipt:
    """What became of one submission.

    Attributes
    ----------
    stored : bool
        Whether the JSON record was written. False only when the store itself
        could not be written, which is reported to the user as a failure —
        there is nothing else holding the note.
    store_path : str or None
        Where it went. Reported so an operator running the desktop shell can
        find their own file; it is a local path, not a secret.
    delivered : bool
        Whether the relay accepted the message.
    delivery_detail : str
        One sentence about the delivery attempt, written for the submitter
        rather than for an administrator: it never names the relay host.
    """

    stored: bool
    store_path: str | None
    delivered: bool
    delivery_detail: str

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


def normalise_submission(
    payload: Mapping[str, Any],
    *,
    config: FeedbackConfig,
    environment: Mapping[str, Any] | None = None,
    now: datetime | None = None,
) -> FeedbackSubmission:
    """Validate what the browser sent and stamp what the server knows.

    Parameters
    ----------
    payload : Mapping
        ``message`` (required), and optionally ``category``, ``name``,
        ``email``, ``organisation``, ``rating``, ``contact_consent`` and
        ``context``.
    config : FeedbackConfig
        Supplies the message length limit.
    environment : Mapping, optional
        Server-known facts. Defaults to :func:`server_environment`.
    now : datetime, optional
        Overrides the received-at stamp; the tests use it.

    Returns
    -------
    FeedbackSubmission

    Raises
    ------
    InvalidInputError
        The message is missing or empty, the category is not one of
        :data:`CATEGORIES`, or a field is longer than its limit. Raised as the
        application's own error type so the route answers with the same
        envelope every other rejected input gets.
    """

    message = _text(payload.get("message"), "message", limit=config.max_message_characters)
    if not message:
        raise InvalidInputError(
            "A feedback note needs something to say.",
            field="message",
            hint="Even one line is useful — say what you expected and what happened instead.",
        )
    category = str(payload.get("category") or CATEGORIES[0][0])
    if category not in dict(CATEGORIES):
        raise InvalidInputError(
            f"{category!r} is not one of the feedback categories.",
            field="category",
            hint="Choose one of: " + ", ".join(key for key, _ in CATEGORIES) + ".",
        )
    rating: Any = payload.get("rating")
    rating_value: int | None
    if rating is None or rating == "":
        rating_value = None
    else:
        try:
            rating_value = int(rating)
        except (TypeError, ValueError) as error:
            raise InvalidInputError(
                "The rating must be a whole number from 1 to 5.", field="rating"
            ) from error
        if not 1 <= rating_value <= 5:
            raise InvalidInputError(
                f"The rating must be from 1 to 5; got {rating_value}.", field="rating"
            )
    context = payload.get("context")
    stamp = (now or datetime.now(UTC)).astimezone(UTC)
    return FeedbackSubmission(
        received_at=stamp.isoformat(timespec="seconds").replace("+00:00", "Z"),
        category=category,
        message=message,
        name=_text(payload.get("name"), "name", limit=_SHORT_FIELD_LIMIT),
        email=_text(payload.get("email"), "email", limit=_SHORT_FIELD_LIMIT),
        organisation=_text(payload.get("organisation"), "organisation", limit=_SHORT_FIELD_LIMIT),
        rating=rating_value,
        contact_consent=bool(payload.get("contact_consent")),
        context=dict(context) if isinstance(context, Mapping) else {},
        environment=dict(environment) if environment is not None else server_environment(),
    )


def _text(value: Any, field_name: str, *, limit: int) -> str:
    if value is None:
        return ""
    if not isinstance(value, str):
        raise InvalidInputError(f"The {field_name} must be text.", field=field_name)
    stripped = value.strip()
    if len(stripped) > limit:
        raise InvalidInputError(
            f"The {field_name} is {len(stripped)} characters; the limit is {limit}.",
            field=field_name,
            hint="Attach the detail to a follow-up e-mail rather than pasting it here.",
        )
    return stripped


def server_environment() -> dict[str, Any]:
    """What the server knows about itself, for the stored record.

    Deliberately modest: the host name and the account the server runs as are
    what an administrator needs to find the machine a report came from, and
    nothing here is read from the request, so a client cannot forge it.
    """

    from pytex import __version__

    try:
        account = getpass.getuser()
    except Exception:  # pragma: no cover - getuser has no user on some daemons
        account = "unknown"
    return {
        "pytex_version": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "host": socket.gethostname(),
        "account": account,
    }


def read_store(path: str | Path) -> list[dict[str, Any]]:
    """Every submission the store holds, oldest first.

    Returns an empty list for a store that does not exist yet, which is the
    normal state of a fresh deployment rather than an error.

    Raises
    ------
    ValueError
        The file exists but is not a JSON array. Raised rather than replaced,
        because overwriting a file that turned out to hold something else would
        destroy it.
    """

    target = Path(path)
    if not target.is_file():
        return []
    text = target.read_text(encoding="utf-8").strip()
    if not text:
        return []
    document = json.loads(text)
    if not isinstance(document, list):
        raise ValueError(f"{target} is not a JSON array of feedback records.")
    return [entry for entry in document if isinstance(entry, dict)]


def append_to_store(submission: FeedbackSubmission, *, path: str | Path) -> Path:
    """Append one record to the JSON store and return the file written.

    How the append is done, and why
    -------------------------------
    The store is a JSON *array*, which cannot be appended to by writing at the
    end of the file — so it is read, extended and written whole, through a
    temporary file that is then moved into place. The move is atomic on both
    platform families, so a process that dies mid-write leaves the previous
    complete store rather than a truncated one.

    A line-delimited format would append without the rewrite, and was rejected:
    the file is meant to be opened and read by whoever maintains the
    deployment, and a JSON array is what every tool they have already opens.
    At the volume a feedback file grows — a note at a time, by hand — the
    rewrite costs nothing.
    """

    target = Path(path)
    with _STORE_LOCK:
        target.parent.mkdir(parents=True, exist_ok=True)
        records = read_store(target)
        records.append(submission.to_json())
        temporary = target.with_suffix(target.suffix + ".tmp")
        temporary.write_text(
            json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        temporary.replace(target)
    return target


def relay_feedback(submission: FeedbackSubmission, *, relay: RelayConfig) -> tuple[bool, str]:
    """Send one submission through the configured internal relay.

    Never raises. A relay is a network service on somebody else's schedule, and
    the submission is already on disk by the time this is called, so every
    failure mode here is reported rather than propagated.

    Returns
    -------
    tuple of (bool, str)
        Whether it was delivered, and one sentence saying what happened —
        written for the submitter, so it never names the relay host or a
        credential.
    """

    if not relay.enabled:
        return False, "This deployment does not forward feedback by e-mail; it was filed locally."
    if not relay.recipients:
        return False, "No feedback recipient is configured, so it was filed locally."
    if relay.dry_run:
        _LOGGER.info(
            "Feedback relay dry-run: would send %r to %s",
            submission.subject(relay.subject_prefix),
            ", ".join(relay.recipients),
        )
        return False, "E-mail delivery is in dry-run mode; the note was filed locally."
    if not relay.host:
        _LOGGER.warning("Feedback relay is enabled but relay.host is empty; nothing was sent.")
        return False, "The e-mail relay is not fully configured, so it was filed locally."

    message = EmailMessage()
    message["Subject"] = submission.subject(relay.subject_prefix)
    message["From"] = relay.from_address
    message["To"] = ", ".join(relay.recipients)
    if submission.email:
        # So that a maintainer can answer with Reply rather than by copying the
        # address out of the body — and only when the submitter gave one.
        message["Reply-To"] = submission.email
    message.set_content(submission.describe())

    username, password = relay.resolved_credentials()
    try:
        with smtplib.SMTP(relay.host, relay.port, timeout=relay.timeout_seconds) as server:
            if relay.use_tls:
                server.starttls()
            if username or password:
                server.login(username, password)
            server.send_message(message)
    except Exception as error:
        _LOGGER.warning(
            "Feedback relay failed (host=%s port=%s): %s", relay.host, relay.port, error
        )
        return False, "The e-mail relay could not be reached, so the note was filed locally."
    _LOGGER.info("Feedback relayed to %s", ", ".join(relay.recipients))
    return True, "Your note was filed and e-mailed to the maintainer."


def record_feedback(
    payload: Mapping[str, Any],
    *,
    config: AppConfig | None = None,
    environment: Mapping[str, Any] | None = None,
) -> tuple[FeedbackSubmission, FeedbackReceipt]:
    """Store one submission, then try to deliver it.

    Parameters
    ----------
    payload : Mapping
        What the browser sent; see :func:`normalise_submission`.
    config : AppConfig, optional
        Defaults to :func:`pytex.app.config.load_app_config`.
    environment : Mapping, optional
        Server-known facts to stamp, for a caller that knows more than this
        module does — the server adds which shell answered.

    Returns
    -------
    tuple of (FeedbackSubmission, FeedbackReceipt)

    Raises
    ------
    InvalidInputError
        The submission itself is not usable. Nothing is stored in that case,
        which is correct: there is nothing to store.
    """

    resolved = config if config is not None else load_app_config()
    submission = normalise_submission(payload, config=resolved.feedback, environment=environment)
    try:
        stored_at = append_to_store(submission, path=resolved.feedback.store_path)
    except (OSError, ValueError) as error:
        # The one failure the user must hear about plainly: nothing is holding
        # their note, so telling them it was received would be a lie.
        _LOGGER.error("Feedback could not be stored at %s: %s", resolved.feedback.store_path, error)
        return submission, FeedbackReceipt(
            stored=False,
            store_path=None,
            delivered=False,
            delivery_detail=(
                "PyTex could not write its feedback file, so nothing was kept. "
                "Please e-mail the maintainer directly."
            ),
        )
    delivered, detail = relay_feedback(submission, relay=resolved.relay)
    return submission, FeedbackReceipt(
        stored=True,
        store_path=str(stored_at),
        delivered=delivered,
        delivery_detail=detail,
    )
