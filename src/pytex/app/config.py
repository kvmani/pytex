"""Deployment configuration for the workbench: what the operator decides.

Purpose
-------
Almost everything the application does is decided by the code, and that is
deliberate — a scientific result must not depend on a file somebody edited. A
small number of things are genuinely site decisions, and they are collected
here rather than scattered as constants:

- whether the feedback form is offered, and where its JSON record is written;
- whether feedback is *also* relayed by e-mail, and through which internal
  SMTP server, to which address;
- whether a first-time visitor is greeted and shown the tour.

None of these touch a number the application computes. That is the test for
whether something belongs in this file at all.

When and where to use it
------------------------
Call :func:`load_app_config` once, at server start; the result is cached per
resolved path, so the routes that need it may call it freely. Point the loader
at a file with the ``PYTEX_APP_CONFIG`` environment variable, or place one at
``./pytex_app.yml`` or ``~/.pytex/pytex_app.yml``. With no file at all the
defaults in this module apply, which is a working application: feedback is
collected to a JSON file under the user's home directory, no e-mail is sent,
and the tour is offered.

``config/pytex_app.example.yml`` in the repository is the annotated template.
Every key it carries is a key this module reads.

Secrets
-------
A relay password is never written in the YAML file in a real deployment. The
file names an *environment variable* instead (``username_env_var``,
``password_env_var``), and this module reads it. That is the same split the
sibling ``project_management_software`` deployment uses, for the same reason:
the configuration file is reviewed, copied and backed up, and a password in it
is a password in all three places.

Returns
-------
:class:`AppConfig`, an immutable tree of plain values; see the class.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

__all__ = [
    "DEFAULT_FEEDBACK_RECIPIENTS",
    "AppConfig",
    "ConfigError",
    "FeedbackConfig",
    "RelayConfig",
    "TourConfig",
    "clear_config_cache",
    "config_search_paths",
    "load_app_config",
]

_LOGGER = logging.getLogger("pytex.app.config")

#: Where feedback is mailed when a deployment turns the relay on and does not
#: name a recipient of its own. This is the address the application's author
#: reads; it is the default rather than a hard-coded destination so that a
#: group running its own copy can send feedback to its own maintainer.
DEFAULT_FEEDBACK_RECIPIENTS: tuple[str, ...] = ("kvmani@barc.gov.in",)

#: The environment variable naming the configuration file, checked first.
CONFIG_ENV_VAR = "PYTEX_APP_CONFIG"

#: Filename looked for in the working directory and under ``~/.pytex``.
CONFIG_FILENAME = "pytex_app.yml"


class ConfigError(RuntimeError):
    """The configuration file exists but cannot be used.

    Raised rather than defaulted around: an operator who wrote a configuration
    file meant it to take effect, and silently ignoring a typo in the relay
    host would mean feedback quietly stops being delivered with nothing on
    screen to say so.
    """


@dataclass(frozen=True)
class RelayConfig:
    """The internal SMTP relay feedback is forwarded through.

    Off by default. A deployment with no relay still collects every submission
    — the JSON store is not optional — so turning this on adds delivery rather
    than enabling the feature.

    Attributes
    ----------
    enabled : bool
        Whether to attempt delivery at all.
    dry_run : bool
        Prepare and log the message without connecting to the relay. Useful for
        confirming the recipients and the body on a machine with no relay.
    host, port : str, int
        The relay. Port 25 is the usual internal-relay port, unauthenticated.
    use_tls : bool
        Issue ``STARTTLS`` after connecting.
    username, password : str or None
        Written here only for a test rig. Prefer the ``*_env_var`` fields.
    username_env_var, password_env_var : str or None
        Environment variables holding the credentials, read at load time.
    from_address : str
        The envelope sender. An internal relay usually requires this to be a
        local address it is willing to originate.
    recipients : tuple of str
        Who receives the feedback.
    timeout_seconds : float
        How long to wait on the relay before giving up and recording the
        failure. Short on purpose: a submission must never leave the user
        watching a spinner because a relay is down.
    """

    enabled: bool = False
    dry_run: bool = False
    host: str = ""
    port: int = 25
    use_tls: bool = False
    username: str | None = None
    password: str | None = None
    username_env_var: str | None = "PYTEX_SMTP_USERNAME"
    password_env_var: str | None = "PYTEX_SMTP_PASSWORD"
    from_address: str = "pytex-noreply@localhost"
    recipients: tuple[str, ...] = DEFAULT_FEEDBACK_RECIPIENTS
    subject_prefix: str = "[PyTex feedback]"
    timeout_seconds: float = 10.0

    def resolved_credentials(self) -> tuple[str, str]:
        """The username and password to authenticate with, environment first.

        Returns
        -------
        tuple of str
            Both empty when the relay is unauthenticated, which is the common
            case for an internal relay on port 25.
        """

        username = self.username or ""
        password = self.password or ""
        if self.username_env_var:
            username = os.environ.get(self.username_env_var, username) or ""
        if self.password_env_var:
            password = os.environ.get(self.password_env_var, password) or ""
        return username, password


@dataclass(frozen=True)
class FeedbackConfig:
    """The feedback and feature-request form.

    Attributes
    ----------
    enabled : bool
        Whether the workbench offers the form at all.
    store_path : Path
        The JSON file every submission is appended to. Written whether or not
        the relay is configured, and written *before* delivery is attempted,
        so a relay outage cannot lose a submission.
    invitation : str
        The warm sentence at the top of the form. Configurable because a group
        running its own copy has its own maintainer to name.
    acknowledgement : str
        What the user is told after a successful submission.
    max_message_characters : int
        Length limit on the message body, so the store cannot be filled by a
        pasted file.
    """

    enabled: bool = True
    store_path: Path = field(default_factory=lambda: Path.home() / ".pytex" / "feedback.json")
    invitation: str = (
        "PyTex is built by researchers, for researchers — and it gets better every time "
        "one of you tells us what is missing, what is confusing, or what would save you an "
        "afternoon. Nothing is too small to send. We read every note."
    )
    acknowledgement: str = (
        "Thank you — your note has been recorded, and it genuinely helps. "
        "If you left an e-mail address we may write back with a question or when it is done."
    )
    max_message_characters: int = 8000


@dataclass(frozen=True)
class TourConfig:
    """The welcome message and the guided tour.

    Attributes
    ----------
    enabled : bool
        Whether a visitor who has not seen the tour is offered it. On by
        default; a shared instrument PC where everyone has already seen it is
        the reason it can be turned off.
    show_every_visit : bool
        Off by default: the tour is remembered as seen in the browser's local
        storage, so a returning user is not greeted again. A demonstration
        machine may want it every time.
    """

    enabled: bool = True
    show_every_visit: bool = False


@dataclass(frozen=True)
class AppConfig:
    """Everything the operator decides, as one immutable object.

    Attributes
    ----------
    feedback, relay, tour
        The three sections.
    source : Path or None
        The file this was read from, or ``None`` when the defaults applied.
        Reported in the log at startup so an operator can see *which* file the
        running server actually used.
    """

    feedback: FeedbackConfig = field(default_factory=FeedbackConfig)
    relay: RelayConfig = field(default_factory=RelayConfig)
    tour: TourConfig = field(default_factory=TourConfig)
    source: Path | None = None


def config_search_paths() -> tuple[Path, ...]:
    """Where a configuration file is looked for, in order.

    Returns
    -------
    tuple of Path
        The environment variable's path first when it is set, then the working
        directory, then the user's ``~/.pytex`` directory.

    Examples
    --------
    >>> paths = config_search_paths()
    >>> paths[-1].name
    'pytex_app.yml'
    """

    paths: list[Path] = []
    from_env = os.environ.get(CONFIG_ENV_VAR)
    if from_env:
        paths.append(Path(from_env).expanduser())
    paths.append(Path.cwd() / CONFIG_FILENAME)
    paths.append(Path.home() / ".pytex" / CONFIG_FILENAME)
    return tuple(paths)


def load_app_config(path: str | Path | None = None) -> AppConfig:
    """Read the deployment configuration, or return the defaults.

    Parameters
    ----------
    path : str or Path, optional
        A specific file. When omitted, :func:`config_search_paths` decides, and
        the first file that exists wins.

    Returns
    -------
    AppConfig
        Never ``None``: a deployment with no configuration file at all gets a
        working application with the defaults documented on each dataclass.

    Raises
    ------
    ConfigError
        The named file is missing, is not a mapping, or carries a key this
        module does not read. An unknown key is an error rather than a warning
        because the common case is a misspelling, and a misspelled
        ``smtp_host`` is a relay that silently never connects.

    Examples
    --------
    >>> config = AppConfig()
    >>> config.feedback.enabled, config.relay.enabled, config.tour.enabled
    (True, False, True)
    """

    resolved = Path(path).expanduser() if path is not None else None
    if resolved is not None and not resolved.is_file():
        raise ConfigError(f"No configuration file at {resolved}.")
    if resolved is None:
        for candidate in config_search_paths():
            if candidate.is_file():
                resolved = candidate
                break
    return _load_cached(str(resolved) if resolved is not None else None)


@lru_cache(maxsize=8)
def _load_cached(resolved: str | None) -> AppConfig:
    if resolved is None:
        return AppConfig()
    path = Path(resolved)
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise ConfigError(f"{path} could not be read as YAML: {error}") from error
    if document is None:
        document = {}
    if not isinstance(document, dict):
        raise ConfigError(f"{path} must contain a mapping at the top level.")
    unknown = set(document) - {"feedback", "relay", "tour"}
    if unknown:
        raise ConfigError(
            f"{path} has unknown top-level section(s): {', '.join(sorted(unknown))}. "
            "The sections are feedback, relay and tour."
        )
    config = AppConfig(
        feedback=_section(FeedbackConfig(), document.get("feedback"), path, "feedback"),
        relay=_section(RelayConfig(), document.get("relay"), path, "relay"),
        tour=_section(TourConfig(), document.get("tour"), path, "tour"),
        source=path,
    )
    _LOGGER.info("Workbench configuration read from %s", path)
    return config


def _section(defaults: Any, values: Any, path: Path, name: str) -> Any:
    """Overlay one YAML section onto its dataclass of defaults."""

    if values is None:
        return defaults
    if not isinstance(values, dict):
        raise ConfigError(f"{path}: section '{name}' must be a mapping.")
    known = {entry.name for entry in defaults.__dataclass_fields__.values()}
    unknown = set(values) - known
    if unknown:
        raise ConfigError(
            f"{path}: section '{name}' has unknown key(s): {', '.join(sorted(unknown))}. "
            f"Accepted keys: {', '.join(sorted(known))}."
        )
    coerced: dict[str, Any] = {}
    for key, value in values.items():
        if value is None and key in {
            "username",
            "password",
            "username_env_var",
            "password_env_var",
        }:
            coerced[key] = None
        elif key == "store_path":
            coerced[key] = Path(str(value)).expanduser()
        elif key == "recipients":
            coerced[key] = tuple(str(entry) for entry in _as_list(value, path, name, key))
        else:
            coerced[key] = value
    return replace(defaults, **coerced)


def clear_config_cache() -> None:
    """Forget every file already read.

    The loader caches by path so that a route may call it per request. A test
    that writes a configuration file, and an operator who edits one and expects
    a restart to pick it up, both need the cache dropped; the restart does it
    by starting a process, and this does it in place.
    """

    _load_cached.cache_clear()


def _as_list(value: Any, path: Path, section: str, key: str) -> list[Any]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return value
    raise ConfigError(f"{path}: {section}.{key} must be a string or a list of strings.")
