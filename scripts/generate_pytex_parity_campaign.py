from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from pytex import __version__
from pytex.core import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
)
from pytex.plotting import IPFColorKey
from pytex.texture import ODF, KernelSpec

REPO_ROOT = Path(__file__).resolve().parents[1]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable.")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _canonical_quaternion(rotation: Rotation) -> np.ndarray:
    return rotation.canonicalized().quaternion


def _axis_angle(rotation: Rotation) -> dict[str, Any]:
    return {
        "axis": rotation.axis,
        "angle_deg": rotation.angle_deg,
    }


def _reference_frames() -> tuple[ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    return crystal, specimen


def _phase_from_payload(payload: dict[str, Any]) -> tuple[Phase, ReferenceFrame, ReferenceFrame]:
    crystal, specimen = _reference_frames()
    lattice_payload = payload["lattice_parameters"]
    lattice = Lattice(
        float(lattice_payload["a"]),
        float(lattice_payload["b"]),
        float(lattice_payload["c"]),
        float(lattice_payload["alpha"]),
        float(lattice_payload["beta"]),
        float(lattice_payload["gamma"]),
        crystal_frame=crystal,
    )
    symmetry = SymmetrySpec.from_point_group(str(payload["point_group"]), reference_frame=crystal)
    phase = Phase(
        name=str(payload["phase_id"]),
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    return phase, crystal, specimen


def _orientation_payload(orientation: Orientation) -> dict[str, Any]:
    return {
        "quaternion_wxyz": _canonical_quaternion(orientation.rotation),
        "rotation_matrix": orientation.as_matrix(),
        "euler_bunge_deg": orientation.rotation.to_bunge_euler(degrees=True),
        "axis_angle": _axis_angle(orientation.rotation),
        "symmetry_equivalent_count": orientation.symmetry.order if orientation.symmetry else 1,
    }


def _orientation_from_case(case: dict[str, Any]) -> Orientation:
    phase, crystal, specimen = _phase_from_payload(case["phase"])
    input_payload = case["input"]
    operation = case["operation"]
    if operation == "orientation_from_euler":
        return Orientation.from_euler(
            *input_payload["euler_deg"],
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
            convention="bunge",
            degrees=True,
        )
    if operation == "orientation_from_axis_angle":
        return Orientation.from_axis_angle(
            input_payload["axis"],
            np.deg2rad(float(input_payload["angle_deg"])),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
        )
    if operation == "orientation_from_quaternion":
        return Orientation.from_quaternion(
            input_payload["quaternion_wxyz"],
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
        )
    if operation == "orientation_from_matrix":
        return Orientation.from_matrix(
            input_payload["rotation_matrix"],
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
        )
    if operation == "orientation_from_miller":
        plane = MillerPlane.from_hkl(input_payload["plane_hkl"], phase=phase)
        direction = MillerDirection.from_uvw(input_payload["direction_uvw"], phase=phase)
        return Orientation.from_miller(
            plane,
            direction,
            specimen_frame=specimen,
            phase=phase,
            specimen_plane_normal=input_payload["specimen_plane_normal"],
            specimen_direction=input_payload["specimen_direction"],
        )
    raise ValueError(f"Cannot construct orientation for operation {operation!r}.")


def _run_orientation_case(case: dict[str, Any]) -> dict[str, Any]:
    if case["operation"] == "orientation_operations":
        phase, crystal, specimen = _phase_from_payload(case["phase"])
        input_payload = case["input"]
        left = Orientation.from_euler(
            *input_payload["left_euler_deg"],
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
        )
        right = Orientation.from_euler(
            *input_payload["right_euler_deg"],
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
        )
        composed = left.rotation.compose(right.rotation)
        inverse_left = left.rotation.inverse()
        return {
            "left_quaternion_wxyz": _canonical_quaternion(left.rotation),
            "right_quaternion_wxyz": _canonical_quaternion(right.rotation),
            "composed_quaternion_wxyz": _canonical_quaternion(composed),
            "inverse_left_quaternion_wxyz": _canonical_quaternion(inverse_left),
            "misorientation_angle_deg": np.rad2deg(left.distance_to(right)),
            "mapped_test_vector": left.map_crystal_vector(input_payload["test_vector_crystal"]),
        }
    orientation = _orientation_from_case(case)
    result = _orientation_payload(orientation)
    input_payload = case["input"]
    if "map_crystal_vectors" in input_payload:
        vectors = np.asarray(input_payload["map_crystal_vectors"], dtype=np.float64)
        result["mapped_crystal_vectors"] = np.vstack(
            [orientation.map_crystal_vector(vector) for vector in vectors]
        )
    if case["operation"] == "orientation_from_miller":
        plane = MillerPlane.from_hkl(input_payload["plane_hkl"], phase=orientation.phase)
        direction = MillerDirection.from_uvw(
            input_payload["direction_uvw"],
            phase=orientation.phase,
        )
        result["mapped_plane_normal"] = orientation.map_crystal_vector(plane.normal_cartesian)
        result["mapped_direction"] = orientation.map_crystal_vector(direction.unit_vector_cartesian)
    return result


def _run_ipf_case(case: dict[str, Any]) -> dict[str, Any]:
    phase, crystal, specimen = _phase_from_payload(case["phase"])
    input_payload = case["input"]
    orientations = OrientationSet.from_euler_angles(
        input_payload["euler_deg"],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    key = IPFColorKey(
        crystal_symmetry=phase.symmetry,
        specimen_direction=input_payload["specimen_direction"],
        saturation_gamma=float(input_payload.get("saturation_gamma", 0.5)),
    )
    crystal_directions = orientations.map_sample_directions_to_crystal(key.specimen_direction)
    reduced = phase.symmetry.reduce_vectors_to_fundamental_sector(
        crystal_directions,
        antipodal=key.antipodal,
    )
    return {
        "rgb": key.colors_from_orientations(orientations),
        "color_space": "srgb_0_1",
        "specimen_direction": key.specimen_direction,
        "crystal_directions": crystal_directions,
        "sector_reduced_directions": reduced,
    }


def _miller_from_case(case: dict[str, Any]) -> dict[str, Any]:
    phase, _, _ = _phase_from_payload(case["phase"])
    input_payload = case["input"]
    if "plane_hkil_a" in input_payload:
        plane_a = MillerPlane.from_hkil(input_payload["plane_hkil_a"], phase=phase)
        plane_b = MillerPlane.from_hkil(input_payload["plane_hkil_b"], phase=phase)
        direction_a = MillerDirection.from_UVTW(input_payload["direction_uvtw_a"], phase=phase)
        direction_b = MillerDirection.from_UVTW(input_payload["direction_uvtw_b"], phase=phase)
    else:
        plane_a = MillerPlane.from_hkl(input_payload["plane_hkl_a"], phase=phase)
        plane_b = MillerPlane.from_hkl(input_payload["plane_hkl_b"], phase=phase)
        direction_a = MillerDirection.from_uvw(input_payload["direction_uvw_a"], phase=phase)
        direction_b = MillerDirection.from_uvw(input_payload["direction_uvw_b"], phase=phase)
    plane_family = plane_a.symmetry_equivalents().indices
    direction_family = direction_a.symmetry_equivalents().indices
    result = {
        "plane_plane_angle_deg": np.rad2deg(angle_plane_plane_rad(plane_a, plane_b)),
        "direction_direction_angle_deg": np.rad2deg(angle_dir_dir_rad(direction_a, direction_b)),
        "plane_a_d_spacing_angstrom": plane_a.d_spacing_angstrom,
        "plane_a_family_hkl": plane_family,
        "direction_a_family_uvw": direction_family,
    }
    if "plane_hkil_a" in input_payload:
        result["plane_a_hkl"] = plane_a.indices
        result["direction_a_uvw"] = direction_a.indices
    return result


def _run_odf_case(case: dict[str, Any]) -> dict[str, Any]:
    if case["operation"] == "odf_from_xrdml":
        raise ValueError("XRDML ODF parity cases are pending until fixture files are provided.")
    phase, crystal, specimen = _phase_from_payload(case["phase"])
    input_payload = case["input"]
    support = OrientationSet.from_euler_angles(
        input_payload["support_euler_deg"],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    query = OrientationSet.from_euler_angles(
        input_payload["query_euler_deg"],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    kernel_payload = input_payload["kernel"]
    kernel = KernelSpec(
        name=kernel_payload["type"],
        halfwidth_deg=float(kernel_payload["halfwidth_deg"]),
    )
    odf = ODF.from_orientations(support, weights=input_payload["weights"], kernel=kernel)
    density = np.asarray(odf.evaluate(query), dtype=np.float64)
    return {
        "odf_representation": "sampled_density",
        "support_euler_deg": input_payload["support_euler_deg"],
        "query_euler_deg": input_payload["query_euler_deg"],
        "normalized_weights": odf.normalized_weights,
        "query_density": density,
        "summary": {
            "min": float(np.min(density)),
            "max": float(np.max(density)),
            "mean": float(np.mean(density)),
        },
    }


def _run_pole_figure_case(case: dict[str, Any]) -> dict[str, Any]:
    raise ValueError(
        f"Pole figure operation {case['operation']!r} is pending until XRDML fixtures are provided."
    )


def run_case(case: dict[str, Any]) -> dict[str, Any]:
    family = case["operation_family"]
    if family == "orientation":
        return _run_orientation_case(case)
    if family == "ipf_color":
        return _run_ipf_case(case)
    if family == "miller_geometry":
        return _miller_from_case(case)
    if family == "odf":
        return _run_odf_case(case)
    if family == "pole_figure":
        return _run_pole_figure_case(case)
    raise ValueError(f"Unsupported operation family {family!r}.")


def result_for_case(
    campaign: dict[str, Any],
    case: dict[str, Any],
    *,
    case_file: Path,
    created_utc: str,
    include_pending: bool,
) -> dict[str, Any] | None:
    if case["status"] != "active":
        if not include_pending:
            return None
        status = "skipped"
        results: dict[str, Any] = {"reason_pending": case.get("reason_pending", "")}
        notes = ["Case is pending and was not computed by PyTex."]
    else:
        try:
            results = run_case(case)
        except ValueError:
            if not include_pending:
                raise
            status = "skipped"
            results = {"reason_pending": case.get("reason_pending", "Operation is pending.")}
            notes = ["Case operation is not active in the PyTex generator."]
        else:
            status = "active"
            notes = list(case.get("notes", ()))
    return {
        "schema_id": "pytex.parity_result",
        "schema_version": "0.1.0",
        "campaign_id": campaign["campaign_id"],
        "case_id": case["case_id"],
        "case_status": status,
        "producer": {
            "system": "pytex",
            "system_version": __version__,
            "runtime": f"Python {platform.python_version()}",
            "script": "scripts/generate_pytex_parity_campaign.py",
            "created_utc": created_utc,
        },
        "conventions": campaign["conventions"],
        "phase": case["phase"],
        "tolerances": case["tolerances"],
        "results": results,
        "artifacts": [],
        "provenance": {
            "input_sha256": _sha256(case_file),
            "case_file": str(case_file),
            "target_baseline": campaign["target_baseline"],
        },
        "notes": notes,
    }


def generate_campaign(
    case_file: Path,
    output_root: Path,
    *,
    created_utc: str | None = None,
    include_pending: bool = True,
) -> list[Path]:
    case_file = case_file.resolve()
    output_root = output_root.resolve()
    campaign = _read_json(case_file)
    if campaign["schema_id"] != "pytex.parity_case_campaign":
        raise ValueError(f"Unsupported campaign schema in {case_file}.")
    timestamp = created_utc or datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    written: list[Path] = []
    for case in campaign["cases"]:
        result = result_for_case(
            campaign,
            case,
            case_file=case_file,
            created_utc=timestamp,
            include_pending=include_pending,
        )
        if result is None:
            continue
        result_path = output_root / campaign["campaign_id"] / f"{case['case_id']}.json"
        _write_json(result_path, result)
        written.append(result_path)
    return written


def _campaign_files(path: Path) -> list[Path]:
    if path.is_dir():
        return sorted(path.glob("*.json"))
    return [path]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate PyTex parity result JSON from shared MTEX/PyTex campaign cases."
    )
    parser.add_argument("case_path", type=Path, help="Campaign JSON file or directory.")
    parser.add_argument("output_root", type=Path, help="Output root for pytex result JSON files.")
    parser.add_argument(
        "--created-utc",
        help="Override result timestamp, useful for deterministic tests.",
    )
    parser.add_argument(
        "--skip-pending",
        action="store_true",
        help="Do not write skipped result files for pending cases.",
    )
    args = parser.parse_args()
    count = 0
    for case_file in _campaign_files(args.case_path):
        written = generate_campaign(
            case_file,
            args.output_root,
            created_utc=args.created_utc,
            include_pending=not args.skip_pending,
        )
        count += len(written)
    print(f"Wrote {count} PyTex parity result files to {args.output_root}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
