"""The intranet shell: a standard-library HTTP server over the service layer.

Why no web framework
--------------------
The stated deployment is one host on an internal network that colleagues use
without installing anything, and such hosts are routinely offline or behind a
proxy that blocks PyPI. ``python -m pytex.app serve`` has to work on a machine
that has PyTex and nothing else, so the transport is :mod:`http.server` and the
frontend is static files. The full argument, and the conditions under which this
decision should be revisited, are in ``docs/architecture/application_platform.md``.

What this module is *not*
-------------------------
It is not where anything is computed. Every route either returns the manifest,
returns a static file, or hands a decoded JSON body to
:func:`pytex.app.contracts.execute`. A route that did crystallography would be a
capability the desktop shell silently lacked.

Routes
------
``GET /``
    The application page.
``GET /api/health``
    Liveness plus the version, for a reverse proxy or a colleague's bookmark.
``GET /api/manifest``
    The self-describing operation manifest the frontend builds itself from.
``GET /api/shell``
    What this shell can do that the other cannot — currently, whether a result
    can be written to a chosen path instead of downloaded.
``POST /api/call``
    ``{"operation": "...", "params": {...}}`` in, the call envelope out.
``POST /api/export``
    ``{"result": {...}, "format": "csv"|"xlsx"|"json"}`` in, a file download out.
``GET /<path>``
    A file from the bundled static tree.
"""

from __future__ import annotations

import json
import logging
import socket
import sys
import threading
from collections.abc import Mapping
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from pytex.app.contracts import dumps, execute
from pytex.app.errors import ServiceError
from pytex.app.registry import REGISTRY, ServiceRegistry

__all__ = [
    "DEFAULT_HOST",
    "DEFAULT_PORT",
    "STATIC_ROOT",
    "AppServer",
    "PortInUseError",
    "create_server",
    "serve",
]


#: Loopback by default. Serving an intranet is an explicit act, not an accident:
#: the app has no authentication, so binding every interface must be something
#: the operator typed rather than something they inherited from a default.
DEFAULT_HOST = "127.0.0.1"

#: Default port. 8765 is outside the ranges the common development servers use,
#: so the app and a documentation build can run side by side.
DEFAULT_PORT = 8765

#: The bundled frontend.
STATIC_ROOT = Path(__file__).resolve().parent / "static"

_LOGGER = logging.getLogger("pytex.app.server")

#: Content types served from the static tree. Anything not listed is refused
#: rather than sent as ``application/octet-stream``: the tree holds a known set
#: of file kinds, and an unknown extension appearing in it is a packaging
#: mistake worth surfacing.
_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".webmanifest": "application/manifest+json",
    ".woff2": "font/woff2",
}

#: Largest request body accepted, in bytes. Uploaded diffraction patterns are
#: the reason it is measured in megabytes rather than kilobytes.
MAX_REQUEST_BYTES = 64 * 1024 * 1024


class PortInUseError(RuntimeError):
    """The requested port is already taken.

    Raised rather than papered over: a server that reports success while another
    process holds the port is worse than one that refuses to start, because the
    user then goes looking for their data in someone else's application.
    """

    def __init__(self, host: str, port: int, cause: OSError) -> None:
        super().__init__(
            f"Cannot bind {host}:{port} - something is already listening there ({cause}). "
            f"Choose another port with --port, or stop the other server."
        )
        self.host = host
        self.port = port


class _Handler(BaseHTTPRequestHandler):
    """One request. Holds no state; the registry is on the server object."""

    server_version = "PyTex"
    sys_version = ""
    protocol_version = "HTTP/1.1"

    @property
    def registry(self) -> ServiceRegistry:
        return getattr(self.server, "registry", REGISTRY)

    # -- responses ---------------------------------------------------------

    def _send(
        self,
        status: int,
        body: bytes,
        *,
        content_type: str,
        headers: Mapping[str, str] | None = None,
    ) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        # The page loads only what it ships with, so a strict policy costs
        # nothing and removes the whole class of injected-script problems from
        # a tool that will be pointed at unpublished data.
        self.send_header(
            "Content-Security-Policy",
            "default-src 'self'; img-src 'self' data: blob:; style-src 'self' 'unsafe-inline'; "
            "script-src 'self'; connect-src 'self'; base-uri 'none'; form-action 'none'",
        )
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        for key, value in (headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _send_json(self, status: int, payload: Any) -> None:
        self._send(
            status,
            dumps(payload).encode("utf-8"),
            content_type="application/json; charset=utf-8",
            headers={"Cache-Control": "no-store"},
        )

    def _send_error_json(
        self, status: int, code: str, message: str, *, hint: str | None = None
    ) -> None:
        error: dict[str, Any] = {"code": code, "message": message}
        if hint is not None:
            error["hint"] = hint
        self._send_json(status, {"ok": False, "error": error})

    # -- routing -----------------------------------------------------------

    def do_GET(self) -> None:
        path = unquote(urlparse(self.path).path)
        if path == "/api/health":
            from pytex import __version__

            self._send_json(HTTPStatus.OK, {"ok": True, "version": __version__})
            return
        if path == "/api/manifest":
            self._send_json(HTTPStatus.OK, self.registry.manifest())
            return
        if path == "/api/shell":
            from pytex.app.desktop import shell_capabilities

            self._send_json(
                HTTPStatus.OK,
                shell_capabilities(desktop=bool(getattr(self.server, "desktop", False))),
            )
            return
        if path.startswith("/api/"):
            self._send_error_json(
                HTTPStatus.NOT_FOUND,
                "route.unknown",
                f"No API route {path!r}.",
                hint=(
                    "The API routes are /api/health, /api/manifest, /api/shell, /api/call and "
                    "/api/export."
                ),
            )
            return
        self._serve_static(path)

    def do_HEAD(self) -> None:
        self.do_GET()

    def do_POST(self) -> None:
        path = unquote(urlparse(self.path).path)
        if path not in {"/api/call", "/api/export"}:
            self._send_error_json(
                HTTPStatus.NOT_FOUND,
                "route.unknown",
                f"No POST route {path!r}.",
                hint="Operations are invoked by POST to /api/call; files by POST to /api/export.",
            )
            return
        body = self._read_body()
        if body is None:
            return
        try:
            payload = json.loads(body.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            self._send_error_json(
                HTTPStatus.BAD_REQUEST,
                "request.malformed",
                f"The request body is not valid JSON: {error}.",
            )
            return
        if not isinstance(payload, Mapping):
            self._send_error_json(
                HTTPStatus.BAD_REQUEST,
                "request.malformed",
                "The request body must be a JSON object.",
                hint='Send {"operation": "calc.plane_angles", "params": {...}}.',
            )
            return
        if path == "/api/export":
            self._export(payload)
            return
        operation = payload.get("operation")
        if not isinstance(operation, str):
            self._send_error_json(
                HTTPStatus.BAD_REQUEST,
                "request.malformed",
                "The request must name an operation.",
                hint="Fetch /api/manifest for the list of operations.",
            )
            return
        params = payload.get("params") or {}
        if not isinstance(params, Mapping):
            self._send_error_json(
                HTTPStatus.BAD_REQUEST,
                "request.malformed",
                "The 'params' field must be a JSON object.",
            )
            return
        try:
            envelope, status = execute(operation, params, registry=self.registry)
        except Exception:
            # Anything reaching here is a PyTex defect rather than user error:
            # deliberate failures are ServiceErrors and were handled inside
            # execute(). Log it in full, tell the user only that it happened.
            _LOGGER.exception("Unhandled error in operation %s", operation)
            self._send_error_json(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "internal.error",
                f"{operation} failed with an internal error. This is a defect in PyTex, "
                "not a problem with your input.",
                hint="The server log holds the traceback; please report it with your inputs.",
            )
            return
        self._send_json(status, envelope)

    def _export(self, payload: Mapping[str, Any]) -> None:
        """Turn a result the client already has into a downloadable file.

        The result travels back rather than being recomputed from its inputs.
        Recomputing would be tidier to write and wrong to ship: the file must be
        the numbers the user is looking at, and a second evaluation is a second
        chance to differ from them.
        """

        from pytex.app.export import export_result

        result = payload.get("result")
        fmt = payload.get("format")
        if not isinstance(result, Mapping) or not isinstance(fmt, str):
            self._send_error_json(
                HTTPStatus.BAD_REQUEST,
                "request.malformed",
                "An export needs a result object and a format.",
                hint='Send {"result": {...}, "format": "xlsx"}.',
            )
            return
        try:
            data, mime, filename = export_result(result, fmt=fmt)
        except ServiceError as error:
            self._send_json(error.status, {"ok": False, "error": error.to_json()})
            return
        self._send(
            HTTPStatus.OK,
            data,
            content_type=mime,
            headers={
                "Content-Disposition": f'attachment; filename="{filename}"',
                "Cache-Control": "no-store",
            },
        )

    def _read_body(self) -> bytes | None:
        raw_length = self.headers.get("Content-Length")
        if raw_length is None:
            self._send_error_json(
                HTTPStatus.LENGTH_REQUIRED, "request.malformed", "The request needs a body."
            )
            return None
        try:
            length = int(raw_length)
        except ValueError:
            self._send_error_json(
                HTTPStatus.BAD_REQUEST, "request.malformed", "Content-Length is not a number."
            )
            return None
        if length > MAX_REQUEST_BYTES:
            self._send_error_json(
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                "request.too_large",
                f"The request body is {length} bytes; the limit is {MAX_REQUEST_BYTES}.",
                hint="Downsample the image before uploading it.",
            )
            return None
        return self.rfile.read(length)

    # -- static files ------------------------------------------------------

    def _serve_static(self, path: str) -> None:
        relative = "index.html" if path in {"", "/"} else path.lstrip("/")
        candidate = (STATIC_ROOT / relative).resolve()
        try:
            # resolve() plus this check is the whole traversal defence: a path
            # that escapes the static root cannot be relative_to it, whatever
            # combination of "..", encoding, or symlink produced it.
            candidate.relative_to(STATIC_ROOT.resolve())
        except ValueError:
            self._send_error_json(
                HTTPStatus.FORBIDDEN, "path.outside_root", "That path is outside the app."
            )
            return
        if candidate.is_dir():
            candidate = candidate / "index.html"
        if not candidate.is_file():
            self._send_error_json(HTTPStatus.NOT_FOUND, "file.not_found", f"No file at {path!r}.")
            return
        content_type = _CONTENT_TYPES.get(candidate.suffix.lower())
        if content_type is None:
            self._send_error_json(
                HTTPStatus.FORBIDDEN,
                "file.unsupported_type",
                f"The app does not serve {candidate.suffix!r} files.",
            )
            return
        # The frontend is versioned with the package, so a running server may
        # cache it; a colleague picking up a new release gets a new process.
        self._send(
            HTTPStatus.OK,
            candidate.read_bytes(),
            content_type=content_type,
            headers={"Cache-Control": "no-cache"},
        )

    # -- logging -----------------------------------------------------------

    def log_message(self, format: str, *args: Any) -> None:
        """Route request logging through :mod:`logging` rather than stderr."""

        _LOGGER.info("%s %s", self.address_string(), format % args)


class AppServer(ThreadingHTTPServer):
    """A threading HTTP server carrying the registry it serves.

    Threading matters for one reason: a long calculation must not block a
    colleague loading the page. Requests are otherwise independent — the
    application holds no per-user server-side state, which is what makes this
    safe without a session layer.
    """

    daemon_threads = True

    #: ``SO_REUSEADDR`` means two different things on the two platform families,
    #: and only one of them is wanted.
    #:
    #: On Unix it allows rebinding a port still in ``TIME_WAIT`` from a previous
    #: run, which is exactly right for a server a user restarts repeatedly. On
    #: Windows it allows binding a port **another process is actively listening
    #: on**: both sockets bind, and which one receives a connection is not
    #: defined. That is not port reuse, it is silently shadowing whatever else
    #: the user is running — observed in development against a documentation
    #: server on the same port, where PyTex started, reported success, and
    #: served nothing.
    #:
    #: So the option is enabled only where it means what it says. On Windows a
    #: busy port raises, which is the honest outcome.
    allow_reuse_address = sys.platform != "win32"

    def __init__(
        self,
        address: tuple[str, int],
        *,
        registry: ServiceRegistry | None = None,
        desktop: bool = False,
    ) -> None:
        super().__init__(address, _Handler)
        self.registry: ServiceRegistry = registry if registry is not None else REGISTRY
        #: Whether this server was started by the desktop shell. The frontend
        #: asks over ``/api/shell`` rather than sniffing for a window object,
        #: so the one difference between the shells is stated in one place.
        self.desktop: bool = desktop

    @property
    def url(self) -> str:
        """The URL this server is reachable at."""

        address = self.server_address
        host = address[0] if isinstance(address[0], str) else address[0].decode()
        port = int(address[1])
        if host in {"0.0.0.0", "::"}:
            # Reporting, not binding: a colleague needs a name they can type,
            # and the wildcard address is not one.
            host = socket.gethostname()
        return f"http://{host}:{port}"

    def serve_in_background(self) -> threading.Thread:
        """Start serving on a daemon thread and return it.

        Used by the desktop shell, which needs the server running while it owns
        the main thread for the window, and by the tests.
        """

        thread = threading.Thread(target=self.serve_forever, name="pytex-app-server", daemon=True)
        thread.start()
        return thread


def create_server(
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    *,
    registry: ServiceRegistry | None = None,
    desktop: bool = False,
) -> AppServer:
    """Bind a server without serving.

    Parameters
    ----------
    host : str
        Interface to bind. Loopback by default; pass ``"0.0.0.0"`` to serve the
        intranet, which is a deliberate act because the app has no
        authentication.
    port : int
        Port to bind. ``0`` picks a free one, which is what the desktop shell
        and the tests use.
    registry : ServiceRegistry, optional
        Defaults to the application-wide registry.
    desktop : bool
        Whether the desktop shell is starting this server. Reported to the
        frontend over ``/api/shell``, which is how the page learns that it can
        save through a native dialog rather than a browser download.

    Returns
    -------
    AppServer
        Bound and ready; call ``serve_forever`` or ``serve_in_background``.
    """

    if not STATIC_ROOT.is_dir():
        raise RuntimeError(
            f"The frontend is missing from {STATIC_ROOT}. This is a packaging fault: "
            "pytex.app.static must ship with the package."
        )
    try:
        return AppServer((host, port), registry=registry, desktop=desktop)
    except OSError as error:
        raise PortInUseError(host, port, error) from error


def serve(
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    *,
    registry: ServiceRegistry | None = None,
) -> None:
    """Serve until interrupted, announcing the URL."""

    server = create_server(host, port, registry=registry)
    _LOGGER.info("PyTex is serving at %s", server.url)
    print(f"PyTex is serving at {server.url}")
    if host not in {"127.0.0.1", "localhost", "::1"}:
        print("This server has no authentication. Serve it only on a trusted internal network.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping.")
    finally:
        server.shutdown()
        server.server_close()
