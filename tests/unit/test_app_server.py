"""The HTTP layer, tested against a real server on a real socket.

No browser and no mocking of the transport: a server is bound on an ephemeral
loopback port, served on a background thread, and driven with
:mod:`urllib`. That covers the parts a unit test of the handler class would
miss — status codes, content types, header behaviour, and the traversal defence
— while staying fast enough for the base lane.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from collections.abc import Iterator
from typing import Any

import pytest

from pytex.app.server import STATIC_ROOT, AppServer, create_server


@pytest.fixture(scope="module")
def server() -> Iterator[AppServer]:
    instance = create_server("127.0.0.1", 0)
    instance.serve_in_background()
    try:
        yield instance
    finally:
        instance.shutdown()
        instance.server_close()


def get(server: AppServer, path: str) -> tuple[int, dict[str, str], bytes]:
    request = urllib.request.Request(f"{server.url}{path}")
    try:
        with urllib.request.urlopen(request) as response:
            return response.status, dict(response.headers), response.read()
    except urllib.error.HTTPError as error:
        return error.code, dict(error.headers), error.read()


def post(server: AppServer, path: str, body: Any) -> tuple[int, dict[str, Any]]:
    payload = body if isinstance(body, bytes) else json.dumps(body).encode("utf-8")
    request = urllib.request.Request(
        f"{server.url}{path}", data=payload, headers={"Content-Type": "application/json"}
    )
    try:
        with urllib.request.urlopen(request) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


class TestApiRoutes:
    def test_health_reports_the_version(self, server: AppServer) -> None:
        from pytex import __version__

        status, _, body = get(server, "/api/health")
        assert status == 200
        assert json.loads(body) == {"ok": True, "version": __version__}

    def test_the_shell_reports_that_it_cannot_write_local_files(self, server: AppServer) -> None:
        """A page served over HTTP must not be told it can save through Python.

        The desktop shell writes files through a bridge only its own window
        provides. A web page that believed the same would call a bridge that is
        not there, and the export would fail where the browser's own download
        would have worked.
        """

        status, _, body = get(server, "/api/shell")
        assert status == 200
        assert json.loads(body)["shell"] == "web"
        assert json.loads(body)["can_write_local_files"] is False

    def test_the_desktop_shell_says_so_over_the_same_route(self) -> None:
        from pytex.app.server import create_server

        instance = create_server("127.0.0.1", 0, desktop=True)
        instance.serve_in_background()
        try:
            status, _, body = get(instance, "/api/shell")
        finally:
            instance.shutdown()
            instance.server_close()
        assert status == 200
        assert json.loads(body)["shell"] == "desktop"
        assert json.loads(body)["can_write_local_files"] is True

    def test_manifest_is_served(self, server: AppServer) -> None:
        status, headers, body = get(server, "/api/manifest")
        assert status == 200
        assert headers["Content-Type"].startswith("application/json")
        manifest = json.loads(body)
        assert manifest["schema"] == "pytex.app_manifest/1"
        assert manifest["operations"]
        assert manifest["examples"]

    def test_a_call_returns_a_result(self, server: AppServer) -> None:
        status, payload = post(
            server,
            "/api/call",
            {"operation": "calc.plane_angles", "params": {"phase": {"builtin": "ni_fcc"}}},
        )
        assert status == 200
        assert payload["ok"] is True
        assert payload["result"]["table"]["rows"]

    def test_a_service_error_becomes_a_400_envelope(self, server: AppServer) -> None:
        status, payload = post(
            server,
            "/api/call",
            {
                "operation": "calc.plane_angles",
                "params": {"phase": {"builtin": "ni_fcc"}, "planes": [[0, 0, 0]]},
            },
        )
        assert status == 400
        assert payload["ok"] is False
        assert payload["error"]["code"] == "input.invalid"
        assert payload["error"]["message"]

    def test_an_unknown_operation_is_404(self, server: AppServer) -> None:
        status, payload = post(server, "/api/call", {"operation": "nope", "params": {}})
        assert status == 404
        assert payload["error"]["code"] == "operation.unknown"

    def test_a_malformed_body_is_reported_not_crashed(self, server: AppServer) -> None:
        status, payload = post(server, "/api/call", b"{not json")
        assert status == 400
        assert payload["error"]["code"] == "request.malformed"

    def test_a_body_without_an_operation_is_refused(self, server: AppServer) -> None:
        status, payload = post(server, "/api/call", {"params": {}})
        assert status == 400
        assert payload["error"]["hint"]

    def test_an_unknown_api_route_is_404_json(self, server: AppServer) -> None:
        status, headers, body = get(server, "/api/nothing")
        assert status == 404
        assert headers["Content-Type"].startswith("application/json")
        assert json.loads(body)["error"]["code"] == "route.unknown"


class TestStaticFiles:
    def test_the_root_serves_the_application_page(self, server: AppServer) -> None:
        status, headers, body = get(server, "/")
        assert status == 200
        assert headers["Content-Type"].startswith("text/html")
        assert b"PyTex Workbench" in body

    def test_the_stylesheet_and_entry_script_are_reachable(self, server: AppServer) -> None:
        for path, kind in (("/app.css", "text/css"), ("/js/main.js", "text/javascript")):
            status, headers, body = get(server, path)
            assert status == 200, path
            assert headers["Content-Type"].startswith(kind)
            assert body

    def test_a_missing_file_is_404(self, server: AppServer) -> None:
        status, _, _ = get(server, "/js/not-a-file.js")
        assert status == 404

    @pytest.mark.parametrize(
        "path",
        [
            "/../pyproject.toml",
            "/../../../../etc/passwd",
            "/%2e%2e/%2e%2e/server.py",
            "/js/../../server.py",
        ],
    )
    def test_traversal_out_of_the_static_root_is_refused(
        self, server: AppServer, path: str
    ) -> None:
        status, _, body = get(server, path)
        assert status in {403, 404}
        assert b"def " not in body

    def test_only_known_file_kinds_are_served(self, server: AppServer, tmp_path: Any) -> None:
        # A stray file in the static tree must not be served as an opaque blob:
        # the tree holds a known set of kinds, and anything else is a packaging
        # mistake worth surfacing rather than shipping.
        stray = STATIC_ROOT / "stray.txt"
        stray.write_text("secret", encoding="utf-8")
        try:
            status, _, body = get(server, "/stray.txt")
            assert status == 403
            assert b"secret" not in body
        finally:
            stray.unlink()


class TestSecurityHeaders:
    def test_every_response_carries_a_strict_policy(self, server: AppServer) -> None:
        for path in ("/", "/api/manifest"):
            _, headers, _ = get(server, path)
            policy = headers["Content-Security-Policy"]
            assert "default-src 'self'" in policy
            assert "script-src 'self'" in policy
            assert headers["X-Content-Type-Options"] == "nosniff"


class TestBinding:
    def test_the_url_is_reported(self, server: AppServer) -> None:
        assert server.url.startswith("http://127.0.0.1:")

    def test_a_busy_port_raises_rather_than_shadowing(self, server: AppServer) -> None:
        """A second server on the same port must fail, not silently shadow.

        Windows treats ``SO_REUSEADDR`` as permission to bind a port another
        process is actively listening on, which produced exactly this bug in
        development: PyTex started, reported success, and served nothing while a
        different application answered the requests.
        """

        from pytex.app.server import PortInUseError

        port = int(server.server_address[1])
        with pytest.raises(PortInUseError):
            create_server("127.0.0.1", port)


class TestCommandLine:
    """``python -m pytex.app`` is how both shells are started."""

    def test_operations_prints_the_manifest(self, capsys: Any) -> None:
        from pytex.app.__main__ import main

        assert main(["operations"]) == 0
        manifest = json.loads(capsys.readouterr().out)
        assert manifest["schema"] == "pytex.app_manifest/1"
        assert manifest["operations"]

    def test_serve_defaults_to_loopback(self) -> None:
        # Binding every interface must be something an operator typed, because
        # the application has no authentication.
        from pytex.app.__main__ import build_parser

        args = build_parser().parse_args(["serve"])
        assert args.host == "127.0.0.1"

    def test_desktop_defaults_to_an_ephemeral_port(self) -> None:
        from pytex.app.__main__ import build_parser

        args = build_parser().parse_args(["desktop"])
        assert args.port == 0


class TestShellCapabilities:
    """The one place the two shells are allowed to differ."""

    def test_only_file_handling_differs(self) -> None:
        from pytex.app.desktop import shell_capabilities

        desktop = shell_capabilities(desktop=True)
        web = shell_capabilities(desktop=False)
        differing = {key for key in desktop if desktop[key] != web[key]}
        assert differing <= {
            "shell",
            "can_write_local_files",
            "can_read_local_paths",
            "native_window",
        }
        assert web["can_write_local_files"] is False


class TestFrontendIsSelfContained:
    """The page must load nothing from the network.

    An intranet host may be air-gapped, and a CDN reference that works on a
    developer's laptop is an application that shows a blank page in the lab.
    """

    def test_no_external_urls_in_the_static_tree(self) -> None:
        offenders = []
        for path in STATIC_ROOT.rglob("*"):
            if not path.is_file() or path.suffix not in {".html", ".css", ".js"}:
                continue
            text = path.read_text(encoding="utf-8")
            for marker in ("http://", "https://"):
                for line in text.splitlines():
                    if (
                        marker in line
                        and "www.w3.org" not in line
                        and not line.strip().startswith("*")
                    ):
                        offenders.append(f"{path.name}: {line.strip()[:80]}")
        assert not offenders, f"the frontend must not reference the network: {offenders}"

    def test_elements_marked_hidden_are_actually_hidden(self) -> None:
        """Every overlay in the page closes itself by setting ``hidden``.

        The palette and the help drawer are full-viewport ``display: flex``
        overlays, and a class selector outranks the user agent's
        ``[hidden] { display: none }``. Without an explicit rule in the
        stylesheet both sit open over the application from the moment the page
        loads: the workbench appears behind two dark scrims, and every click a
        user makes is swallowed by the topmost one instead of reaching the
        control it was aimed at. The page is unusable, and nothing in the
        JavaScript is wrong, so the stylesheet is where this is pinned.
        """

        import re

        html = (STATIC_ROOT / "index.html").read_text(encoding="utf-8")
        stylesheet = (STATIC_ROOT / "app.css").read_text(encoding="utf-8")
        css = re.sub(r"/\*.*?\*/", "", stylesheet, flags=re.S)
        assert " hidden" in html, "no element relies on the attribute; revisit this test"
        rule = re.search(r"\[hidden\]\s*\{([^}]*)\}", css)
        assert rule is not None, "app.css must give [hidden] a rule of its own"
        body = rule.group(1).replace(" ", "")
        assert "display:none!important" in body

    def test_the_crystal_drag_is_not_bound_to_the_drawing_it_replaces(self) -> None:
        """Turning the crystal must survive the redraw that turning causes.

        Each step of a drag rebuilds the SVG and discards the old one. A
        `pointermove` handler owned by that SVG therefore fires once and then
        belongs to a detached element, and the pointer capture goes with it: the
        crystal nudges once however far the user pulls, silently. The handler
        has to live on the frame, which is built once per mount.

        This is asserted against the source because there is no JavaScript test
        runner in this repository, and the alternative to a blunt check is no
        check at all for a defect that a passing Python suite cannot see.
        """

        import re

        source = (STATIC_ROOT / "js" / "panels" / "crystal.js").read_text(encoding="utf-8")
        attach = re.search(r"function attachPointer\(node\) \{(.*?)\n  \}", source, flags=re.S)
        assert attach is not None, "attachPointer has been renamed; revisit this invariant"
        assert "pointermove" not in attach.group(1)
        assert "pointerdown" not in attach.group(1)
        assert "frame.element.addEventListener('pointermove'" in source

    def test_the_entry_points_exist(self) -> None:
        assert (STATIC_ROOT / "index.html").is_file()
        assert (STATIC_ROOT / "app.css").is_file()
        assert (STATIC_ROOT / "js" / "main.js").is_file()
