"""Render the Sphinx site into the package, so an installed PyTex carries its docs.

Why this script exists
----------------------
The workbench serves its documentation from ``/docs/``, and
:func:`pytex.app.server.docs_root` looks for it in four places. On a development
machine the third and fourth candidates -- ``docs/_build/html`` and
``docs/site/_build/html`` -- are usually populated, because somebody has run
Sphinx in the checkout. That is *why the feature appears to work in development
and does not work anywhere else*: both directories are build output, both are in
``.gitignore``, and neither travels in a wheel, an sdist, or a clone.

The deployment this repository is built for is an office intranet host, often
air-gapped, where PyTex is installed from a wheel and there is no checkout, no
network, and no Sphinx. For ``/docs/`` to work there, the built HTML has to be
*inside the distribution*. This script puts it there: it renders ``docs/site``
into ``src/pytex/app/static/docs``, which
``[tool.setuptools.package-data]`` ships and which ``docs_root`` finds as its
second candidate -- ahead of both checkout paths, so an installed copy prefers
its own bundled docs over whatever happens to be lying in a source tree.

Run it before building a distribution::

    python scripts/build_docs_bundle.py
    python -m build

The bundle is generated, so it is git-ignored and must never be committed; the
repository's cardinal rule is that nothing a command here can regenerate is
tracked unless documentation, a test, or a pinned baseline names it. Nothing
does: the *packaging* names it, and packaging happens on the machine that builds
the wheel.

MathJax
-------
Math must render with no network. ``docs/site/conf.py`` already vendors the
``tex-chtml-full`` bundle under ``_static/mathjax`` and points ``mathjax_path``
at it, and those files *are* tracked, so the bundle produced here is offline by
construction. This script checks that the rendered output actually contains it
rather than assuming, because a silently CDN-linked build looks perfect on the
machine that made it and shows raw TeX on the one that matters.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Where the wheel expects to find the docs. Must stay in step with
#: `docs_root()` candidate 2 and with the `package-data` glob in pyproject.toml.
BUNDLE_ROOT = REPO_ROOT / "src" / "pytex" / "app" / "static" / "docs"

#: Doctrees are Sphinx's incremental cache. Written outside the bundle so they
#: are not copied into the wheel, where they would add tens of megabytes of
#: pickles that nothing reads.
DOCTREE_ROOT = REPO_ROOT / "docs" / "_build" / "bundle-doctrees"

SOURCE_ROOT = REPO_ROOT / "docs" / "site"

#: Rendered subtrees a *served* copy does not need. `_sources` is the MyST
#: original of every page, which the reader reaches through the site itself, and
#: `.doctrees` is a build cache. Excluding them is the difference between a
#: wheel of tens of megabytes and one substantially larger, for no capability.
PRUNE = (".doctrees", "_sources")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--keep-sources",
        action="store_true",
        help="Keep the _sources tree in the bundle (larger; lets readers view page source).",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete any existing bundle before building, rather than building over it.",
    )
    return parser


def build(*, keep_sources: bool = False, clean: bool = False) -> int:
    """Render the site into the package and report what was produced."""

    if clean and BUNDLE_ROOT.exists():
        shutil.rmtree(BUNDLE_ROOT)

    command = [
        sys.executable,
        "-m",
        "sphinx",
        "-b",
        "html",
        "-d",
        str(DOCTREE_ROOT),
        str(SOURCE_ROOT),
        str(BUNDLE_ROOT),
    ]
    print("$", " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    if completed.returncode != 0:
        print("Sphinx failed; the bundle is incomplete and must not be shipped.", file=sys.stderr)
        return completed.returncode

    if not keep_sources:
        for name in PRUNE:
            victim = BUNDLE_ROOT / name
            if victim.exists():
                shutil.rmtree(victim)

    # myst-nb writes its executed notebooks *beside* the output directory, which
    # for this build means inside the package. Nothing serves them and nothing
    # should ship them, so they go before the wheel is built.
    stray = BUNDLE_ROOT.parent / "jupyter_execute"
    if stray.is_dir():
        shutil.rmtree(stray)

    problems = _verify()
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    if problems:
        return 1

    files = sum(1 for path in BUNDLE_ROOT.rglob("*") if path.is_file())
    size_mb = sum(p.stat().st_size for p in BUNDLE_ROOT.rglob("*") if p.is_file()) / 1e6
    print(f"bundled {files} files ({size_mb:.1f} MB) into {BUNDLE_ROOT.relative_to(REPO_ROOT)}")
    return 0


def _verify() -> list[str]:
    """The checks worth failing the build over, rather than discovering in the office."""

    problems: list[str] = []
    index = BUNDLE_ROOT / "index.html"
    if not index.is_file():
        problems.append(f"{index} was not produced; /docs/ would 404 on the server.")
        return problems

    mathjax = BUNDLE_ROOT / "_static" / "mathjax" / "tex-chtml-full.js"
    if not mathjax.is_file():
        problems.append(
            "the vendored MathJax bundle is missing from the build; every derivation "
            "would render as raw TeX on a host with no network."
        )

    # A page that still points at a CDN is the failure this whole arrangement
    # exists to prevent, and it is invisible on a machine that has a network.
    for page in list(BUNDLE_ROOT.glob("*.html"))[:20]:
        text = page.read_text(encoding="utf-8", errors="ignore")
        if "cdn.jsdelivr.net" in text or "cdnjs.cloudflare.com" in text:
            problems.append(f"{page.name} references a CDN; the build is not offline-safe.")
            break
    return problems


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    return build(keep_sources=args.keep_sources, clean=args.clean)


if __name__ == "__main__":
    raise SystemExit(main())
