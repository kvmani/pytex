from __future__ import annotations

import argparse
from pathlib import Path

from pytex import __version__
from pytex.core import (
    PYTEX_CANONICAL_CONVENTIONS,
    FrameDomain,
    FrameTransform,
    Handedness,
    Orientation,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _cmd_info(_: argparse.Namespace) -> int:
    print(f"PyTex {__version__}")
    print(f"Canonical convention: {PYTEX_CANONICAL_CONVENTIONS.name}")
    print(f"Repository root: {_repo_root()}")
    return 0


def _cmd_docs_inventory(_: argparse.Namespace) -> int:
    repo_root = _repo_root()
    tex_docs = sorted((repo_root / "docs" / "tex").rglob("*.tex"))
    figures = sorted((repo_root / "docs" / "figures").glob("*.svg"))
    print("LaTeX documents:")
    for path in tex_docs:
        print(f"  - {path.relative_to(repo_root)}")
    print("SVG figures:")
    for path in figures:
        print(f"  - {path.relative_to(repo_root)}")
    return 0


def _cmd_core_demo(_: argparse.Namespace) -> int:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
        provenance=ProvenanceRecord.minimal("pytex-demo"),
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
        provenance=ProvenanceRecord.minimal("pytex-demo"),
    )
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    orientation = Orientation(
        rotation=Rotation.from_bunge_euler(45.0, 35.0, 15.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    transform = FrameTransform.identity(specimen)
    print("Orientation matrix:")
    print(orientation.as_matrix())
    print("Identity transform matrix:")
    print(transform.rotation_matrix)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pytex")
    subparsers = parser.add_subparsers(dest="command", required=True)

    info_parser = subparsers.add_parser("info", help="Show repository and convention information.")
    info_parser.set_defaults(func=_cmd_info)

    docs_parser = subparsers.add_parser("docs", help="Inspect documentation assets.")
    docs_subparsers = docs_parser.add_subparsers(dest="docs_command", required=True)
    docs_inventory_parser = docs_subparsers.add_parser(
        "inventory",
        help="List LaTeX and SVG assets.",
    )
    docs_inventory_parser.set_defaults(func=_cmd_docs_inventory)

    core_parser = subparsers.add_parser("core", help="Run core-model demonstrations.")
    core_subparsers = core_parser.add_subparsers(dest="core_command", required=True)
    core_demo_parser = core_subparsers.add_parser("demo", help="Show a minimal core-model demo.")
    core_demo_parser.set_defaults(func=_cmd_core_demo)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
