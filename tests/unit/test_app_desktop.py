"""The desktop shell: the one thing it does differently, and whether it works.

The desktop shell is the web shell in a window, and there is exactly one place
they diverge — writing a file. That divergence is worth testing precisely,
because the failure it replaced was silent: the embedded web view accepts a
download click and does nothing at all with it, so every export in the
application appeared to succeed and produced no file.

pywebview is an optional extra, so the module is exercised against a stub. What
is under test here is PyTex's half of the contract: that the bytes the page
sends arrive on disk unchanged, that a dismissed dialog is reported as "not
saved" rather than as success, and that the shell tells the frontend what it can
do rather than leaving the page to guess.
"""

from __future__ import annotations

import base64
import sys
import types
from typing import Any

import pytest

from pytex.app import desktop
from pytex.app.desktop import SaveBridge, shell_capabilities


class _StubWindow:
    """A window whose save dialog answers with whatever the test wants."""

    def __init__(self, answer: Any) -> None:
        self.answer = answer
        self.calls: list[dict[str, Any]] = []

    def create_file_dialog(self, kind: Any, **kwargs: Any) -> Any:
        self.calls.append({"kind": kind, **kwargs})
        return self.answer


@pytest.fixture
def stub_webview(monkeypatch: pytest.MonkeyPatch):
    """Install a stub ``webview`` module and return a factory for windows."""

    module = types.ModuleType("webview")
    module.SAVE_DIALOG = 30  # type: ignore[attr-defined]
    module.windows = []  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "webview", module)

    def use(answer: Any) -> _StubWindow:
        window = _StubWindow(answer)
        module.windows = [window]  # type: ignore[attr-defined]
        return window

    return use


class TestSaveBridge:
    def test_the_bytes_the_page_sent_reach_the_chosen_path(self, stub_webview, tmp_path) -> None:
        target = tmp_path / "spots.csv"
        window = stub_webview(str(target))
        payload = b"h,k,l\n1,1,1\n"

        written = SaveBridge().save_file("spots.csv", base64.b64encode(payload).decode())

        assert written == str(target)
        assert target.read_bytes() == payload
        assert window.calls[0]["save_filename"] == "spots.csv"

    def test_binary_content_survives_the_round_trip(self, stub_webview, tmp_path) -> None:
        """An xlsx is a zip: one altered byte and the reader refuses it."""

        target = tmp_path / "table.xlsx"
        stub_webview(str(target))
        payload = bytes(range(256)) * 8

        SaveBridge().save_file("table.xlsx", base64.b64encode(payload).decode())

        assert target.read_bytes() == payload

    def test_a_dialog_answering_with_a_sequence_is_understood(
        self, stub_webview, tmp_path
    ) -> None:
        """pywebview has returned both a string and a tuple across versions."""

        target = tmp_path / "figure.png"
        stub_webview((str(target),))

        written = SaveBridge().save_file("figure.png", base64.b64encode(b"\x89PNG").decode())

        assert written == str(target)
        assert target.read_bytes() == b"\x89PNG"

    def test_a_dismissed_dialog_writes_nothing_and_says_so(self, stub_webview, tmp_path) -> None:
        stub_webview(None)
        bridge = SaveBridge()

        assert bridge.save_file("spots.csv", base64.b64encode(b"x").decode()) is None
        assert bridge.last_path is None
        assert list(tmp_path.iterdir()) == []

    def test_only_the_file_name_is_offered_to_the_dialog(self, stub_webview, tmp_path) -> None:
        """A result title becomes a filename; it must not become a directory."""

        target = tmp_path / "result.csv"
        window = stub_webview(str(target))

        SaveBridge().save_file("/etc/passwd", base64.b64encode(b"x").decode())

        assert window.calls[0]["save_filename"] == "passwd"


class TestShellCapabilities:
    def test_the_web_shell_cannot_write_local_files(self) -> None:
        capabilities = shell_capabilities(desktop=False)
        assert capabilities["shell"] == "web"
        assert capabilities["can_write_local_files"] is False

    def test_the_desktop_shell_can(self) -> None:
        capabilities = shell_capabilities(desktop=True)
        assert capabilities["shell"] == "desktop"
        assert capabilities["can_write_local_files"] is True


def test_native_window_opens_maximized(monkeypatch: pytest.MonkeyPatch) -> None:
    """The status/activity bar must not begin underneath the OS taskbar."""

    captured: dict[str, Any] = {}
    module = types.ModuleType("webview")

    def create_window(title: str, url: str, **kwargs: Any) -> object:
        captured.update(title=title, url=url, **kwargs)
        return object()

    module.create_window = create_window  # type: ignore[attr-defined]
    module.start = lambda: None  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "webview", module)

    assert desktop._open_native_window(object(), "http://127.0.0.1:1/") == 0
    assert captured["maximized"] is True
    assert captured["min_size"] == (960, 640)
