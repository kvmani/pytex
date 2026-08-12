"""Command line for the application: ``python -m pytex.app``.

Two verbs, matching the two shells:

``serve``
    Run the intranet server. Loopback unless ``--host`` says otherwise, because
    the application has no authentication and exposing it must be deliberate.
``desktop``
    Open the desktop window, which is the same application against a server on a
    loopback port nobody else can reach.

A third, ``operations``, prints the manifest — useful for scripting against the
service layer without a browser, and for confirming what a build exposes.
"""

from __future__ import annotations

import argparse
import logging
import sys

from pytex.app.contracts import dumps
from pytex.app.registry import REGISTRY
from pytex.app.server import DEFAULT_HOST, DEFAULT_PORT


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for ``python -m pytex.app``."""

    parser = argparse.ArgumentParser(
        prog="python -m pytex.app",
        description="The PyTex workbench: a desktop app and an intranet server over one codebase.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    serve = subparsers.add_parser("serve", help="Serve the application over HTTP.")
    serve.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help=(
            "Interface to bind. Defaults to loopback. Pass 0.0.0.0 to serve colleagues on the "
            "intranet — the application has no authentication, so do that only on a trusted "
            "network."
        ),
    )
    serve.add_argument(
        "--port", type=int, default=DEFAULT_PORT, help=f"Port to bind (default: {DEFAULT_PORT})."
    )

    desktop = subparsers.add_parser("desktop", help="Open the desktop window.")
    desktop.add_argument(
        "--port",
        type=int,
        default=0,
        help="Loopback port for the embedded server. 0 picks a free one (default).",
    )
    desktop.add_argument(
        "--browser",
        action="store_true",
        help="Skip the native window and open the default browser instead.",
    )

    subparsers.add_parser("operations", help="Print the operation manifest as JSON.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command line. Returns the process exit code."""

    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.command == "operations":
        print(dumps(REGISTRY.manifest(), indent=2))
        return 0

    if args.command == "serve":
        from pytex.app.server import serve

        serve(args.host, args.port)
        return 0

    from pytex.app.desktop import open_desktop

    return open_desktop(port=args.port, prefer_browser=args.browser)


if __name__ == "__main__":
    sys.exit(main())
