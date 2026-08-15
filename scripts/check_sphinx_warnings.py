from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

_SUMMARY_PATTERN = re.compile(r"build succeeded,\s+(\d+)\s+warnings?\.", re.IGNORECASE)


def warning_count_from_output(output: str) -> int:
    """Return Sphinx's summary count, with a line-count fallback for failed builds."""

    matches = _SUMMARY_PATTERN.findall(output)
    if matches:
        return int(matches[-1])
    return sum("WARNING:" in line for line in output.splitlines())


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the Sphinx site and fail when its warning baseline grows."
    )
    parser.add_argument("--max-warnings", type=int, required=True)
    parser.add_argument("--source", type=Path, default=Path("docs/site"))
    parser.add_argument("--output", type=Path, default=Path("docs/_build/html"))
    parser.add_argument("--builder", default="html")
    parser.add_argument("--nitpicky", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.max_warnings < 0:
        raise SystemExit("--max-warnings must be non-negative")

    command = [
        sys.executable,
        "-m",
        "sphinx",
        "-E",
        "-b",
        str(args.builder),
    ]
    if args.nitpicky:
        command.append("-n")
    command.extend((str(args.source), str(args.output)))
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    combined = completed.stdout + completed.stderr
    sys.stdout.write(combined)
    if completed.returncode != 0:
        print(f"Sphinx failed with exit code {completed.returncode}.", file=sys.stderr)
        return completed.returncode

    warning_count = warning_count_from_output(combined)
    if warning_count > args.max_warnings:
        print(
            f"Sphinx warning ratchet failed: {warning_count} warnings exceed "
            f"the baseline of {args.max_warnings}.",
            file=sys.stderr,
        )
        return 1
    print(f"Sphinx warning ratchet passed: {warning_count} warnings (maximum {args.max_warnings}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
