"""The desktop shell: the web shell in a window.

There is no second application here. This module starts the same server on a
loopback port and opens it in a window, so what a desktop user runs is what an
intranet user runs, down to the byte. That is the whole point of the split — a
feature cannot exist in one and not the other, because there is only one.

Two window strategies, in order of preference:

1. **pywebview**, when installed: a native, chrome-free window using the
   platform's own web view. Optional extra ``pytex[desktop]``.
2. **The default browser**, otherwise. Less polished, but it means
   ``python -m pytex.app desktop`` works on a bare install with no optional
   dependency at all — which matters for the same reason the server has no web
   framework.

The only genuine difference between the shells is file handling, and it is
reported to the frontend through :func:`shell_capabilities` rather than sniffed.
"""

from __future__ import annotations

import logging
import webbrowser
from typing import Any

from pytex.app.server import AppServer, create_server

__all__ = ["native_window_available", "open_desktop", "shell_capabilities"]

_LOGGER = logging.getLogger("pytex.app.desktop")

#: Window size on first open. Wide enough that the stage keeps a usable aspect
#: ratio once the fixed-width control rail is subtracted.
WINDOW_SIZE = (1440, 900)


def native_window_available() -> bool:
    """Whether a native window can be opened in this environment."""

    try:
        import webview  # noqa: F401  # type: ignore[import-not-found]
    except ImportError:
        return False
    return True


def shell_capabilities(*, desktop: bool) -> dict[str, Any]:
    """Report what this shell can do that the other cannot.

    Feature code never branches on this; only file handling does. Keeping the
    difference in one reported object, rather than in scattered checks, is what
    stops the two shells from drifting into two applications.
    """

    return {
        "shell": "desktop" if desktop else "web",
        # The desktop shell can be handed a path and write to it; the web shell
        # can only stream a download the browser then places wherever it places
        # downloads.
        "can_write_local_files": desktop,
        "can_read_local_paths": desktop,
        "native_window": desktop and native_window_available(),
    }


def open_desktop(*, port: int = 0, prefer_browser: bool = False) -> int:
    """Start the loopback server and open the application in a window.

    Parameters
    ----------
    port : int
        Loopback port. ``0`` picks a free one, which is the default because a
        fixed port would collide with a second instance.
    prefer_browser : bool
        Skip the native window even when pywebview is installed.

    Returns
    -------
    int
        Process exit code.
    """

    server = create_server("127.0.0.1", port)
    server.serve_in_background()
    url = server.url
    _LOGGER.info("Desktop shell serving at %s", url)
    try:
        if prefer_browser or not native_window_available():
            return _open_in_browser(server, url)
        return _open_native_window(server, url)
    finally:
        server.shutdown()
        server.server_close()


def _open_native_window(server: AppServer, url: str) -> int:
    import webview

    print(f"PyTex is open in a desktop window ({url}).")
    window = webview.create_window(
        "PyTex Workbench",
        url,
        width=WINDOW_SIZE[0],
        height=WINDOW_SIZE[1],
        min_size=(960, 640),
    )
    # webview.start() owns the main thread until the window closes, which is
    # why the server runs on a daemon thread rather than the other way round.
    webview.start()
    del window
    del server
    return 0


def _open_in_browser(server: AppServer, url: str) -> int:
    print(f"PyTex is open in your browser at {url}.")
    print("Install the optional extra 'pytex[desktop]' for a native window.")
    print("Press Ctrl+C here to close it.")
    webbrowser.open(url)
    try:
        # Nothing else to do: the browser is the window, and closing a browser
        # tab cannot be detected, so the process ends when the user says so.
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping.")
    return 0
