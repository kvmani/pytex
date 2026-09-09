from __future__ import annotations

import argparse
import base64
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
from pytex.contracts import JSON_CONTRACT_SCHEMA_VERSION
from pytex.core import (
    PYTEX_CANONICAL_CONVENTIONS,
    FrameTransform,
    Orientation,
    ProvenanceRecord,
    Rotation,
    SymmetrySpec,
    crystal_frame,
    list_phase_fixtures,
    sample_frame,
    specimen_frame,
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
    notes_root = repo_root / "docs" / "site" / "theory"
    notes = sorted(p for p in notes_root.glob("*.md") if p.stem != "index")
    figures = sorted((repo_root / "docs" / "figures").glob("*.svg"))
    print("Theory and algorithm notes:")
    for path in notes:
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
    provenance = ProvenanceRecord.minimal("pytex-demo")
    crystal = crystal_frame(provenance=provenance)
    specimen = specimen_frame(provenance=provenance)
    sample = sample_frame(provenance=provenance)
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    orientation = Orientation(
        rotation=Rotation.from_bunge_euler(45.0, 35.0, 15.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    transform = FrameTransform.from_axis_correspondence(
        specimen,
        sample,
        {"x": "RD", "y": "TD", "z": "ND"},
        provenance=provenance,
    )
    print("Orientation matrix:")
    print(orientation.as_matrix())
    print("Specimen frame:")
    print(specimen.describe())
    print("Sample frame:")
    print(sample.describe())
    print("Specimen-to-sample transform:")
    print(transform.describe())
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


def _cmd_phases_inventory(_: argparse.Namespace) -> int:
    for fixture in list_phase_fixtures():
        crystal_system = str(fixture.metadata.get("crystal_system", "unknown"))
        print(f"{fixture.fixture_id}: {fixture.phase_name} [{crystal_system}]")
    return 0


def _cmd_contracts_inventory(_: argparse.Namespace) -> int:
    print(f"JSON contract schema version: {JSON_CONTRACT_SCHEMA_VERSION}")
    for schema_id in (
        "pytex.core.miller_bravais_plane",
        "pytex.core.miller_bravais_direction",
        "pytex.core.pattern_center",
        "pytex.core.ebsd_detector_geometry",
        "pytex.core.ebsd_calibration_geometry",
        "pytex.diffraction.diffraction_pattern",
    ):
        print(f"  - {schema_id}")
    return 0


def _cmd_examples_inventory(_: argparse.Namespace) -> int:
    repo_root = _repo_root()
    examples = sorted((repo_root / "examples").glob("**/*"))
    notebooks = sorted((repo_root / "docs" / "site" / "tutorials" / "notebooks").glob("*.ipynb"))
    print("Examples:")
    for path in examples:
        if path.is_file():
            print(f"  - {path.relative_to(repo_root)}")
    print("Tutorial notebooks:")
    for path in notebooks:
        print(f"  - {path.relative_to(repo_root)}")
    return 0


def _cmd_hrem_ctf(args: argparse.Namespace) -> int:
    from pytex.diffraction.hrem import DoubleCorrectionMode, MicroscopeAberrations

    mode_map = {
        "conventional": DoubleCorrectionMode.UNCORRECTED,
        "cs_corrected": DoubleCorrectionMode.CS_CORRECTED,
        "double_corrected": DoubleCorrectionMode.DOUBLE_CORRECTED,
        "ncsi": DoubleCorrectionMode.NCSI,
    }
    mode = mode_map.get(args.mode, DoubleCorrectionMode.DOUBLE_CORRECTED)
    cs_mm = (
        (args.cs * 1e-3)
        if args.cs is not None
        else (0.0 if mode != DoubleCorrectionMode.UNCORRECTED else 1.0)
    )
    aberr = MicroscopeAberrations(
        energy_kev=args.voltage,
        defocus_angstrom=args.defocus,
        cs_mm=cs_mm,
        c5_mm=args.c5,
        astigmatism_angstrom=args.astigmatism,
        astigmatism_angle_deg=args.astigmatism_angle,
        coma_angstrom=args.coma,
        coma_angle_deg=args.coma_angle,
        trefoil_angstrom=args.trefoil,
        trefoil_angle_deg=args.trefoil_angle,
        focal_spread_angstrom=args.focal_spread,
        aperture_cutoff_mrad=args.aperture,
        mode=mode,
    )
    ctf = aberr.evaluate_ctf_1d(azimuth_deg=args.azimuth)
    print("Contrast Transfer Function (CTF) Diagnostics:")
    print(f"  Voltage: {aberr.energy_kev:.1f} kV (lambda = {aberr.wavelength_angstrom:.5f} A)")
    print(f"  Mode: {aberr.mode.value}")
    print(f"  Defocus: {aberr.defocus_angstrom:.1f} A")
    print(f"  Cs: {aberr.cs_um:.2f} um ({aberr.cs_mm:.4f} mm)")
    print(f"  Focal spread: {aberr.focal_spread_angstrom:.1f} A")
    first_z = ctf.first_zero_q_inv_angstrom
    first_d = ctf.first_zero_d_spacing_angstrom
    info_q = ctf.information_limit_q_inv_angstrom
    info_d = ctf.information_limit_d_spacing_angstrom
    print(f"  First zero crossing: {first_z:.3f} A^-1 (d = {first_d:.2f} A)")
    print(f"  Information limit (1/e^2): {info_q:.3f} A^-1 (d = {info_d:.2f} A)")
    if aberr.has_azimuthal_aberrations:
        band = aberr.evaluate_ctf_azimuthal()
        best, worst = band.point_resolution_range_angstrom
        print(f"  Cut azimuth: {ctf.azimuth_deg:.1f} deg (transfer is not isotropic)")
        print(
            f"  Point resolution over azimuth: {best:.2f} - {worst:.2f} A "
            f"(anisotropy {band.resolution_anisotropy_angstrom:.2f} A, "
            f"worst at {band.worst_azimuth_deg:.1f} deg)"
        )
    if getattr(args, "report", False):
        print("\nExplainable Diagnostics:\n" + ctf.describe())
        if aberr.has_azimuthal_aberrations:
            print("\n" + aberr.evaluate_ctf_azimuthal().describe())
    return 0


def _cmd_hrem_simulate(args: argparse.Namespace) -> int:
    from pytex.adapters.abtem import simulate_hrem
    from pytex.core import crystal_frame
    from pytex.core.fixtures import get_phase_fixture
    from pytex.diffraction.hrem import AtomicSnapshot, DoubleCorrectionMode, MicroscopeAberrations

    sample_type = getattr(args, "sample_type", "crystalline")
    if sample_type == "amorphous":
        snap = AtomicSnapshot.amorphous_sample(
            species=getattr(args, "species", "C"),
            density_g_cm3=1.8,
            dimensions_angstrom=(16.0, 16.0, 10.0),
            seed=42,
        )
    elif Path(args.phase).is_file():
        snap = AtomicSnapshot.from_xyz(Path(args.phase).read_text(encoding="utf-8"))
    else:
        fixture = get_phase_fixture(args.phase)
        cs_sys = fixture.metadata.get("crystal_system", "fcc")
        phase = fixture.load_phase(crystal_frame=crystal_frame(cs_sys))
        if sample_type == "vacancy":
            snap = AtomicSnapshot.crystalline_with_vacancy(
                phase, supercell=(2, 2, 2), vacancy_count=1
            )
        elif sample_type == "dislocation":
            snap = AtomicSnapshot.crystalline_with_dislocation(phase, supercell=(3, 3, 2))
        else:
            snap = AtomicSnapshot.from_phase(phase, supercell=(2, 2, 2), zone_axis=(0, 0, 1))

    mode_map = {
        "conventional": DoubleCorrectionMode.UNCORRECTED,
        "cs_corrected": DoubleCorrectionMode.CS_CORRECTED,
        "double_corrected": DoubleCorrectionMode.DOUBLE_CORRECTED,
        "ncsi": DoubleCorrectionMode.NCSI,
    }
    mode = mode_map.get(args.mode, DoubleCorrectionMode.DOUBLE_CORRECTED)
    defocus = (
        args.defocus
        if args.defocus is not None
        else (-30.0 if mode != DoubleCorrectionMode.UNCORRECTED else -500.0)
    )
    cs_mm = (
        (args.cs * 1e-3)
        if args.cs is not None
        else (0.0 if mode != DoubleCorrectionMode.UNCORRECTED else 1.0)
    )

    aberr = MicroscopeAberrations(
        energy_kev=args.voltage,
        defocus_angstrom=defocus,
        cs_mm=cs_mm,
        focal_spread_angstrom=getattr(args, "focal_spread", 10.0),
        mode=mode,
    )

    result = simulate_hrem(snap, aberr, sampling_angstrom=args.sampling)
    print("HREM Simulation Completed:")
    print(f"  Sample: {snap.label} ({snap.natoms} atoms)")
    print(
        f"  Microscope: {aberr.energy_kev:.1f} kV, {aberr.mode.value} "
        f"(Delta_f = {aberr.defocus_angstrom:.1f} A, Cs = {aberr.cs_um:.1f} um)"
    )
    print(
        f"  Image size: {result.image.shape[1]}x{result.image.shape[0]} px "
        f"({result.extent_angstrom[0]:.1f} x {result.extent_angstrom[1]:.1f} A)"
    )
    print(f"  Pixel size: {result.pixel_size_angstrom:.3f} A/px")
    print(f"  Michelson contrast: {result.michelson_contrast * 100.0:.1f}%")
    print(f"  Point resolution: {result.point_resolution_angstrom:.2f} A")
    print(f"  Information limit: {result.information_limit_angstrom:.2f} A")

    if args.output:
        out_p = Path(args.output)
        if out_p.suffix.lower() == ".json":
            data = {
                "sample": snap.label,
                "natoms": snap.natoms,
                "voltage_kv": aberr.energy_kev,
                "mode": aberr.mode.value,
                "defocus_angstrom": aberr.defocus_angstrom,
                "cs_um": aberr.cs_um,
                "michelson_contrast": result.michelson_contrast,
                "weber_contrast": result.weber_contrast,
                "point_resolution_angstrom": result.point_resolution_angstrom,
                "information_limit_angstrom": result.information_limit_angstrom,
            }
            out_p.write_text(json.dumps(data, indent=2), encoding="utf-8")
            print(f"  Wrote JSON result: {out_p}")
        else:
            b64 = result.to_png_base64()
            raw = base64.b64decode(b64.split(",", 1)[1])
            out_p.write_bytes(raw)
            print(f"  Wrote micrograph image: {out_p}")

    if getattr(args, "report", False):
        print("\n" + result.describe())

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
        help="List scientific-note and SVG assets.",
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

    phases_parser = subparsers.add_parser("phases", help="Inspect built-in phase fixtures.")
    phases_subparsers = phases_parser.add_subparsers(dest="phases_command", required=True)
    phases_inventory_parser = phases_subparsers.add_parser(
        "inventory",
        help="List built-in phase fixtures.",
    )
    phases_inventory_parser.set_defaults(func=_cmd_phases_inventory)

    contracts_parser = subparsers.add_parser("contracts", help="Inspect JSON contract surfaces.")
    contracts_subparsers = contracts_parser.add_subparsers(
        dest="contracts_command",
        required=True,
    )
    contracts_inventory_parser = contracts_subparsers.add_parser(
        "inventory",
        help="List high-value public JSON contract schema ids.",
    )
    contracts_inventory_parser.set_defaults(func=_cmd_contracts_inventory)

    examples_parser = subparsers.add_parser("examples", help="Inspect example assets.")
    examples_subparsers = examples_parser.add_subparsers(dest="examples_command", required=True)
    examples_inventory_parser = examples_subparsers.add_parser(
        "inventory",
        help="List examples and tutorial notebooks.",
    )
    examples_inventory_parser.set_defaults(func=_cmd_examples_inventory)

    hrem_parser = subparsers.add_parser(
        "hrem",
        help="Simulate High-Resolution TEM micrographs and evaluate CTF.",
    )
    hrem_subparsers = hrem_parser.add_subparsers(dest="hrem_command", required=True)

    hrem_sim_parser = hrem_subparsers.add_parser(
        "simulate",
        help="Run multislice or phase-object HRTEM simulation.",
    )
    hrem_sim_parser.add_argument(
        "--phase",
        type=str,
        default="ni_fcc",
        help="Built-in phase fixture ID or path to an XYZ structure file.",
    )
    hrem_sim_parser.add_argument(
        "--sample-type",
        choices=["crystalline", "vacancy", "dislocation", "amorphous"],
        default="crystalline",
        help="Specimen microstructure morphology.",
    )
    hrem_sim_parser.add_argument(
        "--voltage",
        type=float,
        default=200.0,
        help="Accelerating voltage in keV (kV).",
    )
    hrem_sim_parser.add_argument(
        "--mode",
        choices=["conventional", "cs_corrected", "double_corrected", "ncsi"],
        default="double_corrected",
        help="Electron optical aberration correction mode.",
    )
    hrem_sim_parser.add_argument(
        "--defocus",
        type=float,
        default=None,
        help="Objective lens defocus in Angstrom.",
    )
    hrem_sim_parser.add_argument(
        "--cs",
        type=float,
        default=None,
        help="Spherical aberration Cs in micrometers.",
    )
    hrem_sim_parser.add_argument(
        "--sampling",
        type=float,
        default=0.15,
        help="Real-space pixel sampling in Angstrom.",
    )
    hrem_sim_parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save simulation results (JSON or PNG).",
    )
    hrem_sim_parser.add_argument(
        "--report",
        action="store_true",
        help="Print detailed scientific explainable report.",
    )
    hrem_sim_parser.set_defaults(func=_cmd_hrem_simulate)

    hrem_ctf_parser = hrem_subparsers.add_parser(
        "ctf",
        help=(
            "Evaluate the contrast transfer function along one azimuth, with damping "
            "envelopes and the azimuthal resolution spread of a non-round lens."
        ),
    )
    hrem_ctf_parser.add_argument(
        "--voltage",
        type=float,
        default=200.0,
        help="Accelerating voltage in keV (kV).",
    )
    hrem_ctf_parser.add_argument(
        "--mode",
        choices=["conventional", "cs_corrected", "double_corrected", "ncsi"],
        default="double_corrected",
        help="Electron optical aberration correction mode.",
    )
    hrem_ctf_parser.add_argument(
        "--defocus",
        type=float,
        default=-50.0,
        help="Objective lens defocus in Angstrom.",
    )
    hrem_ctf_parser.add_argument(
        "--cs",
        type=float,
        default=0.0,
        help="Spherical aberration Cs in micrometers.",
    )
    hrem_ctf_parser.add_argument(
        "--focal-spread",
        type=float,
        default=10.0,
        help="Focal spread Delta in Angstrom.",
    )
    hrem_ctf_parser.add_argument(
        "--c5",
        type=float,
        default=0.0,
        help="Fifth-order spherical aberration C5 in millimetres.",
    )
    hrem_ctf_parser.add_argument(
        "--astigmatism",
        type=float,
        default=0.0,
        help="Two-fold astigmatism amplitude C12 in Angstrom.",
    )
    hrem_ctf_parser.add_argument(
        "--astigmatism-angle",
        type=float,
        default=0.0,
        help="Azimuth phi12 of two-fold astigmatism in degrees.",
    )
    hrem_ctf_parser.add_argument(
        "--coma",
        type=float,
        default=0.0,
        help="Axial coma amplitude C21 in Angstrom.",
    )
    hrem_ctf_parser.add_argument(
        "--coma-angle",
        type=float,
        default=0.0,
        help="Azimuth phi21 of axial coma in degrees.",
    )
    hrem_ctf_parser.add_argument(
        "--trefoil",
        type=float,
        default=0.0,
        help="Three-fold astigmatism (trefoil) amplitude C23 in Angstrom.",
    )
    hrem_ctf_parser.add_argument(
        "--trefoil-angle",
        type=float,
        default=0.0,
        help="Azimuth phi23 of trefoil in degrees.",
    )
    hrem_ctf_parser.add_argument(
        "--azimuth",
        type=float,
        default=0.0,
        help=(
            "Azimuth in degrees of the radial cut through the back focal plane. "
            "Only a lens with a non-round aberration transfers differently by azimuth."
        ),
    )
    hrem_ctf_parser.add_argument(
        "--aperture",
        type=float,
        default=None,
        help="Objective aperture cutoff semiangle in mrad.",
    )
    hrem_ctf_parser.add_argument(
        "--report",
        action="store_true",
        help="Print detailed scientific explainable report.",
    )
    hrem_ctf_parser.set_defaults(func=_cmd_hrem_ctf)

    return parser


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
