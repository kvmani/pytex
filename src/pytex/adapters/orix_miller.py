from __future__ import annotations

from typing import Any

import numpy as np

from pytex.core import (
    MillerDirection,
    MillerDirectionSet,
    MillerPlane,
    MillerPlaneSet,
    Phase,
    direction_uvtw_to_uvw_array,
    plane_hkil_to_hkl_array,
)


def _require_orix() -> tuple[Any, Any]:
    try:
        from orix.crystal_map import Phase as OrixPhase  # type: ignore[import-untyped]
        from orix.vector import Miller as OrixMiller  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "The PyTex Miller orix adapter requires the optional 'orix' dependency. "
            "Install PyTex with the 'adapters' extra."
        ) from exc
    return OrixPhase, OrixMiller


def _shape_rows(values: Any, *, columns: int, name: str) -> tuple[np.ndarray, bool]:
    array = np.asarray(values)
    if array.size == 0:
        raise ValueError(f"{name} must expose {columns}-index coordinates.")
    if array.ndim == 1:
        if array.shape[0] != columns:
            raise ValueError(f"{name} must expose shape ({columns},) or (n, {columns}).")
        return np.ascontiguousarray(array[None, :]), True
    if array.ndim != 2 or array.shape[1] != columns:
        raise ValueError(f"{name} must expose shape ({columns},) or (n, {columns}).")
    return np.ascontiguousarray(array), int(array.shape[0]) == 1


def _rounded_int_rows(values: Any, *, columns: int, name: str) -> tuple[np.ndarray, bool]:
    rows, scalar = _shape_rows(values, columns=columns, name=name)
    rounded = np.rint(rows).astype(np.int64)
    if not np.allclose(rows, rounded, atol=1e-8, rtol=0.0):
        raise ValueError(f"{name} must expose integer Miller indices.")
    return rounded, scalar


def to_orix_phase(pytex_phase: Phase) -> Any:
    orix_phase_cls, _ = _require_orix()
    return orix_phase_cls(
        name=pytex_phase.name,
        point_group=pytex_phase.symmetry.point_group,
        space_group=pytex_phase.space_group_number,
    )


def to_orix_miller_plane(planes: MillerPlane | MillerPlaneSet) -> Any:
    _, orix_miller_cls = _require_orix()
    if isinstance(planes, MillerPlane):
        return orix_miller_cls(hkl=planes.indices, phase=to_orix_phase(planes.phase))
    return orix_miller_cls(hkl=planes.indices, phase=to_orix_phase(planes.phase))


def to_orix_miller_direction(directions: MillerDirection | MillerDirectionSet) -> Any:
    _, orix_miller_cls = _require_orix()
    if isinstance(directions, MillerDirection):
        return orix_miller_cls(uvw=directions.indices, phase=to_orix_phase(directions.phase))
    return orix_miller_cls(uvw=directions.indices, phase=to_orix_phase(directions.phase))


def from_orix_miller(
    miller: Any,
    *,
    phase: Phase,
) -> MillerPlane | MillerPlaneSet | MillerDirection | MillerDirectionSet:
    preferred = getattr(miller, "coordinate_format", None)
    attribute_order: tuple[tuple[str, str], ...]
    if isinstance(preferred, str) and preferred in {"hkl", "hkil", "uvw", "UVTW"}:
        kind = "plane" if preferred in {"hkl", "hkil"} else "direction"
        attribute_order = ((preferred, kind),)
    else:
        attribute_order = (
            ("hkl", "plane"),
            ("hkil", "plane"),
            ("uvw", "direction"),
            ("UVTW", "direction"),
        )
    for attribute, kind in attribute_order:
        raw = getattr(miller, attribute, None)
        if raw is None:
            continue
        columns = 4 if attribute in {"hkil", "UVTW"} else 3
        rows, scalar = _rounded_int_rows(raw, columns=columns, name=f"orix Miller.{attribute}")
        if kind == "plane":
            indices = plane_hkil_to_hkl_array(rows) if attribute == "hkil" else rows
            if scalar:
                return MillerPlane(indices=indices[0], phase=phase)
            return MillerPlaneSet(indices=indices, phase=phase)
        indices = direction_uvtw_to_uvw_array(rows) if attribute == "UVTW" else rows
        if scalar:
            return MillerDirection(indices=indices[0], phase=phase)
        return MillerDirectionSet(indices=indices, phase=phase)
    raise ValueError("Expected an orix Miller object exposing hkl, hkil, uvw, or UVTW data.")


__all__ = [
    "from_orix_miller",
    "to_orix_miller_direction",
    "to_orix_miller_phase",
    "to_orix_miller_plane",
    "to_orix_phase",
]


to_orix_miller_phase = to_orix_phase
