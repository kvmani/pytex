from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from pytex import __version__
from pytex.adapters import (
    read_benchmark_manifest,
    read_ebsd_import_manifest,
    read_experiment_manifest,
    read_transformation_manifest,
    read_validation_manifest,
    read_workflow_result_manifest,
)
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


def _cmd_docs_build(args: argparse.Namespace) -> int:
    repo_root = _repo_root()
    output_dir = repo_root / "docs" / "_build" / "html"
    if getattr(args, "clean", False) and output_dir.exists():
        for path in sorted(output_dir.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    command = [
        sys.executable,
        "-m",
        "sphinx",
        "-b",
        "html",
        str(repo_root / "docs" / "site"),
        str(output_dir),
    ]
    completed = subprocess.run(command, check=False, cwd=repo_root)
    if completed.returncode == 0:
        print(output_dir)
    return int(completed.returncode)


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


def _cmd_validate_repo(_: argparse.Namespace) -> int:
    repo_root = _repo_root()
    return int(
        subprocess.run(
            [sys.executable, "scripts/check_repo_integrity.py"],
            check=False,
            cwd=repo_root,
        ).returncode
    )


def _cmd_validate_manifests(_: argparse.Namespace) -> int:
    repo_root = _repo_root()
    manifest_paths = sorted((repo_root / "benchmarks").glob("**/*.json"))
    readers = {
        "pytex.benchmark_manifest": read_benchmark_manifest,
        "pytex.ebsd_import_manifest": read_ebsd_import_manifest,
        "pytex.experiment_manifest": read_experiment_manifest,
        "pytex.transformation_manifest": read_transformation_manifest,
        "pytex.validation_manifest": read_validation_manifest,
        "pytex.workflow_result_manifest": read_workflow_result_manifest,
    }
    for manifest_path in manifest_paths:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        schema_id = payload.get("schema_id")
        reader = readers.get(schema_id)
        if reader is None:
            continue
        reader(manifest_path)
        print(manifest_path.relative_to(repo_root))
    return 0


def _cmd_bench_inventory(_: argparse.Namespace) -> int:
    repo_root = _repo_root()
    for manifest_path in sorted((repo_root / "benchmarks").glob("**/*.json")):
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        schema_id = payload.get("schema_id", "unknown")
        print(f"{manifest_path.relative_to(repo_root)} [{schema_id}]")
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
    docs_build_parser = docs_subparsers.add_parser("build", help="Build the Sphinx HTML site.")
    docs_build_parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove the existing HTML build directory before rebuilding.",
    )
    docs_build_parser.set_defaults(func=_cmd_docs_build)

    core_parser = subparsers.add_parser("core", help="Run core-model demonstrations.")
    core_subparsers = core_parser.add_subparsers(dest="core_command", required=True)
    core_demo_parser = core_subparsers.add_parser("demo", help="Show a minimal core-model demo.")
    core_demo_parser.set_defaults(func=_cmd_core_demo)

    validate_parser = subparsers.add_parser("validate", help="Run repository validation helpers.")
    validate_subparsers = validate_parser.add_subparsers(dest="validate_command", required=True)
    validate_repo_parser = validate_subparsers.add_parser(
        "repo",
        help="Run repository integrity checks.",
    )
    validate_repo_parser.set_defaults(func=_cmd_validate_repo)
    validate_manifests_parser = validate_subparsers.add_parser(
        "manifests",
        help="Read and validate benchmark and workflow manifests.",
    )
    validate_manifests_parser.set_defaults(func=_cmd_validate_manifests)

    benchmarks_parser = subparsers.add_parser(
        "benchmarks",
        help="Inspect benchmark and validation manifest assets.",
    )
    benchmarks_subparsers = benchmarks_parser.add_subparsers(
        dest="benchmarks_command",
        required=True,
    )
    benchmarks_inventory_parser = benchmarks_subparsers.add_parser(
        "inventory",
        help="List benchmark and validation manifest files.",
    )
    benchmarks_inventory_parser.set_defaults(func=_cmd_bench_inventory)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
