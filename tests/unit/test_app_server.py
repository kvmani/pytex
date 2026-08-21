"""The HTTP layer, tested against a real server on a real socket.

No browser and no mocking of the transport: a server is bound on an ephemeral
loopback port, served on a background thread, and driven with
:mod:`urllib`. That covers the parts a unit test of the handler class would
miss — status codes, content types, header behaviour, and the traversal defence
— while staying fast enough for the base lane.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from collections.abc import Iterator
from http.server import BaseHTTPRequestHandler
from typing import Any
from unittest import mock
from unittest.mock import Mock

import pytest

from pytex.app.server import STATIC_ROOT, AppServer, _Handler, create_server


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


class TestExperienceAndFeedbackRoutes:
    """The two routes that read the deployment's configuration.

    They are the only routes whose answer a site can change without a build, so
    what they publish — and what they refuse to publish — is a contract.
    """

    def test_experience_publishes_what_the_page_needs_to_greet_a_user(
        self, server: AppServer
    ) -> None:
        status, _, body = get(server, "/api/experience")
        assert status == 200
        document = json.loads(body)
        assert document["ok"] is True
        assert document["feedback"]["enabled"] is True
        assert document["tour"]["enabled"] is True
        assert document["feedback"]["invitation"].strip()
        values = {entry["value"] for entry in document["feedback"]["categories"]}
        assert {"feedback", "feature"} <= values

    def test_experience_never_publishes_the_relay(self, server: AppServer) -> None:
        """The page is told *whether* a note is mailed, never where or as whom.

        The workbench is served without authentication on an internal network,
        so anything this route returns is readable by everyone who can reach
        the port. A relay host and an envelope sender are the administrator's
        business, and the page has no use for either.
        """

        status, _, body = get(server, "/api/experience")
        assert status == 200
        document = json.loads(body)
        assert set(document["feedback"]) == {
            "enabled",
            "invitation",
            "acknowledgement",
            "max_message_characters",
            "categories",
            "relayed",
        }
        serialised = body.decode("utf-8")
        for forbidden in ("smtp", "host", "password", "username", "from_address", "store_path"):
            assert forbidden not in serialised

    def test_a_submission_is_stored_and_acknowledged(
        self, server: AppServer, tmp_path: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from pytex.app.config import AppConfig, FeedbackConfig
        from pytex.app.feedback import read_store

        store = tmp_path / "feedback.json"
        config = AppConfig(feedback=FeedbackConfig(store_path=store))
        monkeypatch.setattr(_Handler, "_app_config", lambda self: config)

        status, document = post(
            server,
            "/api/feedback",
            {
                "message": "A Kearns factor panel would save me an afternoon.",
                "category": "feature",
                "name": "A Researcher",
                "email": "someone@example.invalid",
                "rating": 5,
                "contact_consent": True,
                "context": {"panel": "texture"},
            },
        )
        assert status == 200
        assert document["ok"] is True
        assert document["acknowledgement"].strip()
        assert document["receipt"]["stored"] is True
        assert document["receipt"]["delivered"] is False

        stored = read_store(store)
        assert len(stored) == 1
        assert stored[0]["category"] == "feature"
        assert stored[0]["rating"] == 5
        assert stored[0]["context"] == {"panel": "texture"}
        # Stamped by the server, not by the page: a submission cannot claim to
        # have come from a shell it did not.
        assert stored[0]["environment"]["shell"] == "web"
        assert stored[0]["environment"]["pytex_version"]

    def test_an_empty_note_is_rejected_with_the_ordinary_error_envelope(
        self, server: AppServer
    ) -> None:
        status, document = post(server, "/api/feedback", {"message": "  "})
        assert status == 400
        assert document["ok"] is False
        assert document["error"]["code"] == "input.invalid"
        assert document["error"]["details"]["field"] == "message"

    def test_a_deployment_that_turned_the_form_off_says_so(
        self, server: AppServer, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from pytex.app.config import AppConfig, FeedbackConfig

        config = AppConfig(feedback=FeedbackConfig(enabled=False))
        monkeypatch.setattr(_Handler, "_app_config", lambda self: config)
        status, document = post(server, "/api/feedback", {"message": "hello"})
        assert status == 403
        assert document["error"]["code"] == "feedback.disabled"

    def test_an_unusable_configuration_file_does_not_take_the_workbench_down(
        self, server: AppServer, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A misspelled relay host must not stop anyone drawing a pole figure."""

        from pytex.app.config import ConfigError

        def _broken() -> None:
            raise ConfigError("relay.smtp_host is not a key")

        monkeypatch.setattr("pytex.app.config.load_app_config", _broken)
        status, _, body = get(server, "/api/experience")
        assert status == 200
        assert json.loads(body)["feedback"]["enabled"] is True


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

    @pytest.mark.parametrize(
        "disconnect",
        [BrokenPipeError(), ConnectionAbortedError(), ConnectionResetError()],
    )
    def test_a_client_disconnect_does_not_escape_the_connection_handler(
        self, disconnect: OSError, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Closing a browser page mid-response is expected transport behaviour."""

        base_handle = Mock(side_effect=disconnect)
        monkeypatch.setattr(BaseHTTPRequestHandler, "handle", base_handle)
        handler = object.__new__(_Handler)
        handler.handle()

        base_handle.assert_called_once_with()


class TestConcurrency:
    """One operation runs at a time; everything else still answers.

    The defect behind this: the scientific stack under the services is not
    thread-safe — pyplot's state is global by construction — and two operations
    arriving together were executing on two handler threads. The failure was not
    a wrong number but a Windows access violation that killed the server
    mid-session, with no Python traceback. It became routine the moment a panel
    started drawing two figures on mount.
    """

    def test_operations_do_not_overlap(self, server: AppServer) -> None:
        import threading

        from pytex.app import registry as registry_module

        overlapping = False
        active = 0
        guard = threading.Lock()
        original = registry_module.REGISTRY.call

        def observed(operation: str, request: Any = None) -> Any:
            nonlocal overlapping, active
            with guard:
                active += 1
                overlapping = overlapping or active > 1
            try:
                time.sleep(0.05)
                return original(operation, request)
            finally:
                with guard:
                    active -= 1

        body = {"operation": "calc.plane_angles", "params": {"phase": {"builtin": "ni_fcc"}}}
        with mock.patch.object(registry_module.REGISTRY, "call", observed):
            threads = [
                threading.Thread(target=lambda: post(server, "/api/call", body)) for _ in range(6)
            ]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=60)

        assert not overlapping, "two operations executed at the same time"

    def test_a_page_still_loads_while_an_operation_runs(self, server: AppServer) -> None:
        """The lock must not have quietly turned the server single-threaded."""

        import threading

        from pytex.app import registry as registry_module

        release = threading.Event()
        original = registry_module.REGISTRY.call

        def blocking(operation: str, request: Any = None) -> Any:
            release.wait(timeout=30)
            return original(operation, request)

        body = {"operation": "calc.plane_angles", "params": {"phase": {"builtin": "ni_fcc"}}}
        with mock.patch.object(registry_module.REGISTRY, "call", blocking):
            worker = threading.Thread(target=lambda: post(server, "/api/call", body))
            worker.start()
            try:
                # The calculation is held open; the page must still be served.
                status, _, page = get(server, "/")
                assert status == 200
                assert b"PyTex Workbench" in page
            finally:
                release.set()
                worker.join(timeout=60)


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

    What belongs in this class, now that a browser suite exists
    ----------------------------------------------------------
    These tests read the frontend as *text*, and that is a blunt instrument:
    it pins how the code is written rather than what it does, and it breaks on
    a rename that changed nothing. It was once the only instrument available.
    It is not any more — `tests/browser/workbench.spec.js` drives the real page
    in a real browser — so a behaviour that can be observed at runtime is
    asserted there and not here. Several tests that used to live in this class
    moved: the table cap and its export, the legend keeping focus, the theme
    surviving a reload, the tabs and toolbars that must not overflow, and the
    viewport contract of the shared plot frame.

    What remains is what no runtime observation can make: statements about the
    *shape* of the source. That nothing reaches the network; that only one
    module saves a file; that the browser does no crystallography; that a
    gallery is read from the manifest rather than hard-coded; that a drawing
    exists once rather than in two copies. A passing page cannot demonstrate
    any of those, because each is a claim about what the code does *not*
    contain.
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

    def test_only_one_place_in_the_frontend_saves_a_file(self) -> None:
        """Every file leaves through `saveBlob`, or the desktop shell loses it.

        Handing a blob to an anchor with a `download` attribute works in a
        browser and does nothing at all in the embedded web view, so a second
        copy of those four lines anywhere in the tree is a feature that saves
        in one shell and silently fails in the other. That is exactly how the
        crystal viewer came to publish an SVG the desktop shell could save and
        a PNG it could not, from adjacent lines of the same function.
        """

        offenders = []
        for path in (STATIC_ROOT / "js").rglob("*.js"):
            text = path.read_text(encoding="utf-8")
            if "createObjectURL" in text and path.name != "result.js":
                offenders.append(path.name)
        assert not offenders, f"these files save a file without going through saveBlob: {offenders}"

    def test_the_shell_reserves_one_centralized_message_console(self) -> None:
        """Every surface reports to one log, and the page reserves its strip.

        The console's markup is built by ``js/core/logbook.js`` rather than
        written in the page, so what the static HTML must carry is the mount
        point and nothing else. Two definitions of an entry — one in markup, one
        in script — is exactly the drift this application avoids elsewhere.
        """

        html = (STATIC_ROOT / "index.html").read_text(encoding="utf-8")
        console = (STATIC_ROOT / "js" / "core" / "logbook.js").read_text(encoding="utf-8")
        api = (STATIC_ROOT / "js" / "core" / "api.js").read_text(encoding="utf-8")
        shell = (STATIC_ROOT / "js" / "main.js").read_text(encoding="utf-8")
        stylesheet = (STATIC_ROOT / "app.css").read_text(encoding="utf-8")

        assert 'id="console"' in html
        assert "log.mountConsole(dom.console)" in shell

        # The console grades severity with the same vocabulary Python emits;
        # a level defined on one side only cannot be filtered or coloured.
        for level in ("progress", "debug", "info", "notice", "success", "warning", "error"):
            assert f"{level}:" in console, f"the console must define the {level} level"
        assert "critical:" in console

        for element_id in ("console-toggle", "console-panel", "console-stream", "console-clear"):
            assert element_id in console

        # Records ride back on every envelope, so no panel has to remember to
        # report what it ran, and the poll catches what belongs to no call.
        assert "log.ingest(payload.log)" in api
        assert "/api/log?since=" in console
        assert "CAPACITY = 1000" in console, "a day-long session must stay bounded"

        assert ".console__indicator--busy" in stylesheet
        assert "@keyframes console-spin" in stylesheet

    def test_diffraction_marker_style_is_a_presentation_only_shared_control(self) -> None:
        """Spot styling redraws computed rows; it must not become fake science input."""

        panel = (STATIC_ROOT / "js" / "panels" / "diffraction.js").read_text(encoding="utf-8")
        shared = (STATIC_ROOT / "js" / "core" / "visualstyle.js").read_text(encoding="utf-8")
        assert "markerStyleControl" in panel
        assert "if (state.result) draw()" in panel
        assert "call(operation.id" not in shared
        for control in (
            "Spot shape",
            "Spot fill",
            "Variant encoding",
            "Spot-size scale",
            "Intensity sizing",
            "Parent colour",
            "Variant palette",
            "Product colour",
        ):
            assert control in shared
        for shape in ("circle", "square", "triangle", "diamond", "star", "cross"):
            assert f"['{shape}'" in shared
        for fill in ("filled", "outline"):
            assert f"['{fill}'" in shared
        for encoding in (
            "color",
            "shape",
            "size",
            "color_shape",
            "color_size",
            "shape_size",
            "color_shape_size",
        ):
            assert f"['{encoding}'" in shared
        assert "variantMarkerStyle" in shared
        assert "Transmitted beam (000)" in panel
        assert "(000) transmitted" in panel
        assert "appearance.fill === 'outline'" in panel
        assert "dashed: spot.double_diffraction" in panel

    def test_crystal_object_properties_redraw_geometry_without_recomputing_it(self) -> None:
        """Atoms and planes are view properties; changing them must not call Python again."""

        source = (STATIC_ROOT / "js" / "panels" / "crystal.js").read_text(encoding="utf-8")
        control = source.split("function appearanceControl")[1].split("export function mount")[0]
        assert "call(" not in control
        assert "onChange();" in control
        assert "publicationAppearance(state.appearance, camera)" in source
        for label in (
            "Object properties",
            "Atom size",
            "Atom opacity",
            "Surface finish",
            "Light azimuth",
            "Light elevation",
            "Ambient light",
            "Diffuse light",
            "Specular highlight",
            "Highlight sharpness",
            "Depth cue",
            "Bond colour",
            "Cell opacity",
            "Plane colour",
            "Plane opacity",
            "Direction colour",
            "Annotation size",
        ):
            assert label in source
        for object_name in (
            "showAtoms",
            "showBonds",
            "showCells",
            "showPlanes",
            "showDirections",
            "showGizmo",
        ):
            assert object_name in source
        for lighting_contract in (
            "function spherePaint",
            "function bondGlyph",
            "radialGradient",
            "data-surface",
            "data-depth-opacity",
            "lightingDirection(appearance)",
            "applyTranspose(camera.rotation",
        ):
            assert lighting_contract in source

    def test_texture_contours_have_real_isolines_and_adjustable_shared_levels(self) -> None:
        """A coloured point mosaic is not a contour plot, however smooth it looks."""

        source = (STATIC_ROOT / "js" / "panels" / "texture.js").read_text(encoding="utf-8")
        controls = source.split("function contourStyleControl")[1].split("export function mount")[0]
        assert "call(" not in controls, "contour appearance must not recompute the texture"
        for algorithm in (
            "function interpolatePoleFigure",
            "function crossingPoint",
            "function contourPath",
            "function drawContourGrid",
            "data-contour-level",
            "vector-effect",
        ):
            assert algorithm in source
        for control in (
            "Filled + lines",
            "Filled contours",
            "Contour lines",
            "Automatic levels",
            "Custom levels",
            "Upper colour limit",
            "Colour palette",
            "Line width",
            "Fill opacity",
            "Display grid",
        ):
            assert control in source
        assert "preserveViewport" in source
        assert "pytex-texture-figure.svg" in source
        density = source.split("function renderDensity")[1].split("function renderScatter")[0]
        assert "frame.hoverable(node, point, columns)" in density
        assert "fill: 'transparent'" in density

    def test_the_tem_gallery_is_read_from_the_manifest_not_hard_coded(self) -> None:
        """A fourth practice plate must appear by adding it in Python, not here."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "tem.gallery_pattern" in source
        assert "parameter.name === 'pattern'" in source
        for identifier in ("fcc_al_001", "bcc_fe_110", "hcp_zr_2-1-10"):
            assert identifier not in source

    def test_the_practice_pattern_carries_its_calibration_across(self) -> None:
        """A camera constant retyped by hand is a camera constant retyped wrongly."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "const calibration = result.data.calibration;" in source
        assert "camera_constant_mm_angstrom: calibration.camera_constant_mm_angstrom" in source
        assert "pixel_size_mm: calibration.pixel_size_mm" in source

    def test_the_answer_check_is_made_in_python_not_in_the_browser(self) -> None:
        """Symmetry lives in the library; [110] and [101] are the same bcc pattern."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "request.expected_zone_axis = state.gallery.data.pattern.zone_axis;" in source
        # The verdict is read from the payload, never recomputed here: nothing in
        # the panel touches symmetry operators or compares index triples.
        assert "check.correct" in source
        assert "check.deviation_deg" in source
        for arithmetic in ("operators", "symmetry.operators"):
            assert arithmetic not in source
        # `Math.acos` appears once, in the measurement overlay, and measures the
        # angle at the beam between two *clicked image points*. That is the same
        # class of image-plane arithmetic as the cursor readout's distance and
        # d-spacing, which this panel has always done: it reports what was picked
        # and carries no crystallography. The line that must never appear is one
        # that turns picks into an orientation or an index triple, which is what
        # the operator bans above are for.
        assert source.count("Math.acos") == 1
        assert "angleToReference" in source

    def test_the_simulated_plate_is_legible_as_a_diffraction_pattern(self) -> None:
        """A dark ground, an unmistakable direct beam, and a reciprocal scale bar.

        Asserted against the shared drawing rather than against the solver panel:
        the solver and the SAED simulator paint the same plate, because it is the
        same calculation, and one drawing is what keeps them from drifting apart.
        """

        source = (STATIC_ROOT / "js" / "core" / "saedplot.js").read_text(encoding="utf-8")
        assert "#05070d" in source
        # The beam is sized against the nearest reflection so a dense zone does
        # not have its innermost spots swallowed.
        assert "Math.min(0.5 * nearest" in source
        assert "drawScaleBar" in source
        assert "Å⁻¹" in source
        # And both panels draw through it rather than carrying a copy.
        for panel in ("tem.js", "saedsim.js"):
            panel_source = (STATIC_ROOT / "js" / "panels" / panel).read_text(encoding="utf-8")
            assert "drawSimulatedPattern" in panel_source
            assert "from '../core/saedplot.js'" in panel_source
            assert "function drawSimulatedPattern" not in panel_source

    def test_an_error_on_a_hidden_field_is_still_visible(self) -> None:
        """The likeliest failure in the TEM panel has no on-screen field to land on."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "hiddenFields.add(name)" in source
        assert "function reportError(form, error)" in source
        assert "hiddenFields.has(error.field)" in source
        # Every action routes failures through it; none calls showError directly.
        assert source.count("reportError(state.") == 3

    def test_a_failed_index_does_not_leave_a_correct_verdict_standing(self) -> None:
        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "card--verdict" in source
        assert "details.querySelectorAll('.card--verdict')" in source

    def test_loading_your_own_image_warns_that_the_calibration_is_not_yours(self) -> None:
        """A camera constant from another exposure fails silently and plausibly."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "Set the camera constant and pixel size for this image before" in source
        assert "without ever looking wrong" in source

    def test_the_lattice_overlay_is_refitted_whenever_the_picks_move(self) -> None:
        """An overlay that lags the picks is worse than none: it looks authoritative."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "tem.fit_lattice" in source
        assert "function scheduleFit()" in source
        # Every path that changes a pick or the centre has to re-fit.
        assert source.count("scheduleFit()") >= 6
        # Debounced, because picking is a burst of clicks.
        assert "setTimeout" in source and "clearTimeout" in source

    def test_the_centre_can_be_nudged_and_refined_and_put_back(self) -> None:
        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        # The pad moves whichever pick is selected, the beam included, so there
        # is one nudge path rather than one per kind of pick.
        assert "function nudgeSelected(dx, dy)" in source
        assert "function adoptRefinedCentre()" in source
        assert "function restorePickedCentre()" in source
        assert "state.nudgeStep" in source
        # A centre the fit did not solve for must not be offered for adoption.
        assert "if (!state.fit?.data.centre_refined) return;" in source

    def test_the_basis_vectors_are_drawn_at_the_picks_that_generate_them(self) -> None:
        """The arrows must follow the spots, which is what makes them worth drawing."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "function drawBasisVectors(" in source
        assert "state.fit?.data.basis_vectors" in source
        assert "polygon" in source

    def test_overlay_colours_are_fixed_and_haloed_not_theme_tokens(self) -> None:
        """The plate is dark in both themes; a theme token would vanish in one."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        for constant in (
            "LATTICE_COLOUR",
            "CALCULATED_COLOUR",
            "HALO_COLOUR",
            # Each overlay means a different thing, so each carries its own
            # colour: the fitted grid and the basis vectors shared one, and an
            # arrow into empty space then read as part of the grid.
            "BASIS_COLOUR",
            "PICK_COLOUR",
            "REFINED_COLOUR",
        ):
            assert constant in source
        drawing = source.split("function drawFittedLattice")[1].split("function drawCalculated")[0]
        assert "var(--teal)" not in drawing
        assert "var(--violet)" not in drawing

    def test_a_candidate_is_selected_to_look_at_and_accepted_separately(self) -> None:
        """Looking at a solution must cost nothing and commit to nothing."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "function solutionsCard(" in source
        assert "function acceptSolution(" in source
        assert "Accept this solution" in source
        # Selecting redraws; accepting is what carries the axis onward.
        assert "state.selected = index" in source
        assert "state.tiltForm.setValues({ phase, current_zone_axis: axis })" in source

    def test_the_scored_terms_are_shown_beside_the_fused_number(self) -> None:
        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "function scoreBar(" in source
        for term in ("length_agreement", "angle_agreement", "coverage_agreement"):
            assert term in source
        assert "function deviationCard(" in source

    def test_an_atlas_row_can_be_acted_on_rather_than_transcribed(self) -> None:
        """Retyping indices from a table into a form is how wrong indices get planned."""

        source = (STATIC_ROOT / "js" / "panels" / "tem.js").read_text(encoding="utf-8")
        assert "Choose a destination" in source
        assert "chooseTarget(row.indices, row.target)" in source
        assert "state.tiltForm.setValues({ target_zone_axis: indices })" in source
        # Reachability is never carried by colour alone.
        assert "out of the holder" in source
        assert "s range in one move" in source

    def test_composite_visibility_controls_support_toggle_bulk_and_focus(self) -> None:
        """A 24-variant pattern must not require 23 clicks to isolate one variant."""

        source = (STATIC_ROOT / "js" / "panels" / "diffraction.js").read_text(encoding="utf-8")
        assert "Click a chip to toggle one source." in source
        assert "Show all" in source
        assert "Parent only" in source
        assert "Focus a variant" in source
        assert "state.hidden.clear()" in source
        assert "sourceKeys.filter" in source
        assert "key !== focusKey" in source
        assert "key !== parentKey" in source

    def test_export_buttons_are_generated_from_the_manifest(self) -> None:
        """A format added in Python must reach every panel without an edit here."""

        source = (STATIC_ROOT / "js" / "core" / "result.js").read_text(encoding="utf-8")
        shell = (STATIC_ROOT / "js" / "main.js").read_text(encoding="utf-8")
        assert "export function setExportFormats(" in source
        assert "EXPORT_FORMATS.map((format) =>" in source
        assert "setExportFormats(app.manifest.export_formats)" in shell
        # The hard-coded trio is gone from the render path.
        card = source.split("function tableCard")[1].split("function buildTable")[0]
        for hard_coded in ("'csv'", "'xlsx'", "'json'"):
            assert hard_coded not in card

    def test_the_entry_points_exist(self) -> None:
        assert (STATIC_ROOT / "index.html").is_file()
        assert (STATIC_ROOT / "app.css").is_file()
        assert (STATIC_ROOT / "js" / "main.js").is_file()
