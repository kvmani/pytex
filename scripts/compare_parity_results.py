from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ComparisonIssue:
    case_id: str
    field_path: str
    message: str


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _result_files(root: Path) -> dict[str, Path]:
    files: dict[str, Path] = {}
    for path in sorted(root.glob("**/*.json")):
        payload = _read_json(path)
        if payload.get("schema_id") == "pytex.parity_result":
            files[str(payload["case_id"])] = path
    return files


def _numeric_tolerance(tolerances: dict[str, float], field_path: str) -> tuple[float, float]:
    path = field_path.lower()
    if "rgb" in path:
        return float(tolerances.get("rgb_abs", 1e-8)), 0.0
    if "density" in path:
        return float(tolerances.get("density_abs", 1e-8)), float(
            tolerances.get("density_relative", 0.0)
        )
    if "angle" in path:
        return float(tolerances.get("angle_deg", 1e-8)), 0.0
    if "matrix" in path:
        return float(tolerances.get("matrix_abs", 1e-10)), 0.0
    if "vector" in path or "direction" in path or "axis" in path:
        return float(tolerances.get("vector_abs", 1e-10)), 0.0
    if "weight" in path:
        return float(tolerances.get("weight_abs", 1e-12)), 0.0
    if "spacing" in path or "length" in path:
        return float(tolerances.get("length_abs", 1e-10)), 0.0
    return 1e-10, 0.0


def _as_array(value: Any) -> np.ndarray | None:
    if isinstance(value, int | float):
        return np.asarray(value, dtype=np.float64)
    if isinstance(value, list):
        try:
            return np.asarray(value, dtype=np.float64)
        except (TypeError, ValueError):
            return None
    return None


def _quaternion_close(left: np.ndarray, right: np.ndarray, atol: float, rtol: float) -> bool:
    return bool(
        np.allclose(left, right, atol=atol, rtol=rtol)
        or np.allclose(left, -right, atol=atol, rtol=rtol)
    )


def _compare_values(
    *,
    case_id: str,
    field_path: str,
    expected: Any,
    actual: Any,
    tolerances: dict[str, float],
) -> list[ComparisonIssue]:
    issues: list[ComparisonIssue] = []
    if isinstance(expected, dict) and isinstance(actual, dict):
        expected_keys = set(expected)
        actual_keys = set(actual)
        for missing in sorted(expected_keys - actual_keys):
            issues.append(
                ComparisonIssue(case_id, f"{field_path}.{missing}", "missing actual field")
            )
        for extra in sorted(actual_keys - expected_keys):
            issues.append(
                ComparisonIssue(case_id, f"{field_path}.{extra}", "unexpected actual field")
            )
        for key in sorted(expected_keys & actual_keys):
            issues.extend(
                _compare_values(
                    case_id=case_id,
                    field_path=f"{field_path}.{key}",
                    expected=expected[key],
                    actual=actual[key],
                    tolerances=tolerances,
                )
            )
        return issues
    expected_array = _as_array(expected)
    actual_array = _as_array(actual)
    if expected_array is not None and actual_array is not None:
        if expected_array.shape != actual_array.shape:
            issues.append(
                ComparisonIssue(
                    case_id,
                    field_path,
                    f"shape mismatch: expected {expected_array.shape}, got {actual_array.shape}",
                )
            )
            return issues
        atol, rtol = _numeric_tolerance(tolerances, field_path)
        if "quaternion_wxyz" in field_path.lower():
            close = _quaternion_close(expected_array, actual_array, atol, rtol)
        else:
            close = bool(np.allclose(expected_array, actual_array, atol=atol, rtol=rtol))
        if not close:
            delta = float(np.max(np.abs(expected_array - actual_array)))
            issues.append(
                ComparisonIssue(
                    case_id,
                    field_path,
                    "numeric mismatch: "
                    f"max abs delta {delta:.6g}, atol {atol:.6g}, rtol {rtol:.6g}",
                )
            )
        return issues
    if expected != actual:
        issues.append(
            ComparisonIssue(
                case_id,
                field_path,
                f"value mismatch: expected {expected!r}, got {actual!r}",
            )
        )
    return issues


def compare_result_payloads(
    expected: dict[str, Any],
    actual: dict[str, Any],
) -> list[ComparisonIssue]:
    case_id = str(expected["case_id"])
    issues: list[ComparisonIssue] = []
    for key in ("schema_id", "schema_version", "campaign_id", "case_id", "case_status"):
        if expected.get(key) != actual.get(key):
            issues.append(
                ComparisonIssue(
                    case_id,
                    key,
                    f"expected {expected.get(key)!r}, got {actual.get(key)!r}",
                )
            )
    if expected.get("conventions") != actual.get("conventions"):
        issues.append(ComparisonIssue(case_id, "conventions", "convention metadata differs"))
    if expected.get("phase") != actual.get("phase"):
        issues.append(ComparisonIssue(case_id, "phase", "phase metadata differs"))
    if expected.get("case_status") != "active":
        return issues
    tolerances = {key: float(value) for key, value in expected.get("tolerances", {}).items()}
    issues.extend(
        _compare_values(
            case_id=case_id,
            field_path="results",
            expected=expected.get("results", {}),
            actual=actual.get("results", {}),
            tolerances=tolerances,
        )
    )
    return issues


def compare_result_roots(expected_root: Path, actual_root: Path) -> list[ComparisonIssue]:
    expected_files = _result_files(expected_root)
    actual_files = _result_files(actual_root)
    issues: list[ComparisonIssue] = []
    for case_id in sorted(set(expected_files) - set(actual_files)):
        issues.append(ComparisonIssue(case_id, "", "missing actual result file"))
    for case_id in sorted(set(actual_files) - set(expected_files)):
        issues.append(ComparisonIssue(case_id, "", "unexpected actual result file"))
    for case_id in sorted(set(expected_files) & set(actual_files)):
        issues.extend(
            compare_result_payloads(
                _read_json(expected_files[case_id]),
                _read_json(actual_files[case_id]),
            )
        )
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare MTEX and PyTex parity result JSON files field by field."
    )
    parser.add_argument("expected_root", type=Path, help="Expected/baseline result root.")
    parser.add_argument("actual_root", type=Path, help="Actual result root.")
    args = parser.parse_args()
    issues = compare_result_roots(args.expected_root, args.actual_root)
    if issues:
        for issue in issues:
            print(f"{issue.case_id} {issue.field_path}: {issue.message}")
        return 1
    print("Parity result comparison passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
