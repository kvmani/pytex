from __future__ import annotations

from typing import Any

import numpy as np

from pytex.core.batches import RotationSet
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex, Phase
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec

_PROPER_POINT_GROUP_TO_ORIX = {
    "1": "C1",
    "2": "C2",
    "222": "D2",
    "3": "C3",
    "32": "D3",
    "4": "C4",
    "422": "D4",
    "6": "C6",
    "622": "D6",
    "23": "T",
    "432": "O",
}


def _require_orix() -> tuple[Any, Any, Any, Any, Any]:
    try:
        from orix.crystal_map import Phase as OrixPhase  # type: ignore[import-untyped]
        from orix.quaternion import Orientation as OrixOrientation  # type: ignore[import-untyped]
        from orix.quaternion import Rotation as OrixRotation
        from orix.quaternion import symmetry as orix_symmetry
        from orix.vector import Miller as OrixMiller  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "The PyTex orix adapter requires the optional 'orix' dependency. "
            "Install PyTex with the 'adapters' extra."
        ) from exc
    return OrixPhase, OrixOrientation, OrixRotation, orix_symmetry, OrixMiller


def _quaternion_array_from_orix(obj: Any, *, name: str) -> tuple[np.ndarray, bool]:
    source = getattr(obj, "data", obj)
    array = np.asarray(source, dtype=np.float64)
    if array.shape == (4,):
        return np.ascontiguousarray(array[None, :], dtype=np.float64), True
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError(f"{name} must expose quaternion data with shape (4,) or (n, 4).")
    return np.ascontiguousarray(array, dtype=np.float64), int(array.shape[0]) == 1


def to_orix_symmetry(symmetry: SymmetrySpec) -> Any:
    _, _, _, orix_symmetry, _ = _require_orix()
    attribute = _PROPER_POINT_GROUP_TO_ORIX[symmetry.proper_point_group]
    return getattr(orix_symmetry, attribute)


def from_orix_symmetry(orix_symmetry: Any, *, reference_frame: ReferenceFrame) -> SymmetrySpec:
    point_group = str(getattr(orix_symmetry, "name", "")).strip()
    if not point_group:
        raise ValueError("Unable to derive a point-group name from the orix symmetry object.")
    return SymmetrySpec.from_point_group(point_group, reference_frame=reference_frame)


def to_orix_rotation(rotations: Rotation | RotationSet) -> Any:
    _, _, orix_rotation_cls, _, _ = _require_orix()
    quaternions = (
        rotations.quaternion[None, :]
        if isinstance(rotations, Rotation)
        else rotations.quaternions
    )
    return orix_rotation_cls(quaternions)


def from_orix_rotation(
    rotation: Any,
    *,
    provenance: ProvenanceRecord | None = None,
) -> Rotation | RotationSet:
    quaternions, scalar = _quaternion_array_from_orix(rotation, name="orix rotation")
    if scalar:
        return Rotation(quaternion=quaternions[0], provenance=provenance)
    return RotationSet(quaternions=quaternions, provenance=provenance)


def to_orix_orientation(orientations: Orientation | OrientationSet) -> Any:
    _, orix_orientation_cls, _, _, _ = _require_orix()
    if isinstance(orientations, Orientation):
        quaternions = orientations.rotation.quaternion[None, :]
        symmetry = orientations.symmetry
    else:
        quaternions = orientations.quaternions
        symmetry = orientations.symmetry
    orix_orientation = orix_orientation_cls(quaternions)
    if symmetry is not None:
        orix_orientation.symmetry = to_orix_symmetry(symmetry)
    return orix_orientation


def from_orix_orientation(
    orientation: Any,
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    symmetry: SymmetrySpec | None = None,
    phase: Phase | None = None,
    provenance: ProvenanceRecord | None = None,
) -> Orientation | OrientationSet:
    quaternions, scalar = _quaternion_array_from_orix(orientation, name="orix orientation")
    resolved_symmetry = symmetry
    if resolved_symmetry is None and phase is not None:
        resolved_symmetry = phase.symmetry
    if resolved_symmetry is None and getattr(orientation, "symmetry", None) is not None:
        resolved_symmetry = from_orix_symmetry(
            orientation.symmetry,
            reference_frame=crystal_frame,
        )
    if scalar:
        return Orientation(
            rotation=Rotation(quaternion=quaternions[0], provenance=provenance),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=resolved_symmetry,
            phase=phase,
            provenance=provenance,
        )
    return OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        symmetry=resolved_symmetry,
        phase=phase,
        provenance=provenance,
    )


def _to_orix_phase(phase: Phase) -> Any:
    orix_phase_cls, _, _, _, _ = _require_orix()
    return orix_phase_cls(
        name=phase.name,
        point_group=phase.symmetry.point_group,
        space_group=phase.space_group_number,
    )


def to_orix_plane(plane: CrystalPlane) -> Any:
    _, _, _, _, orix_miller_cls = _require_orix()
    return orix_miller_cls(hkl=plane.miller.indices, phase=_to_orix_phase(plane.phase))


def to_orix_direction(direction: CrystalDirection) -> Any:
    _, _, _, _, orix_miller_cls = _require_orix()
    return orix_miller_cls(uvw=direction.coordinates, phase=_to_orix_phase(direction.phase))


def plane_from_orix_miller(miller: Any, *, phase: Phase) -> CrystalPlane:
    hkil = np.asarray(getattr(miller, "hkil", []), dtype=np.int64).reshape(-1)
    if hkil.size == 4:
        return CrystalPlane.from_miller_bravais(tuple(int(value) for value in hkil), phase=phase)
    hkl = np.asarray(getattr(miller, "hkl", []), dtype=np.int64).reshape(-1)
    if hkl.size != 3:
        raise ValueError("Expected an orix Miller plane with hkl or hkil indices.")
    return CrystalPlane(miller=MillerIndex(hkl, phase=phase), phase=phase)


def direction_from_orix_miller(miller: Any, *, phase: Phase) -> CrystalDirection:
    uvw = np.asarray(getattr(miller, "uvw", []), dtype=np.float64).reshape(-1)
    if uvw.size == 3 and not np.allclose(uvw, 0.0):
        return CrystalDirection(uvw, phase=phase)
    uvtw = np.asarray(getattr(miller, "UVTW", []), dtype=np.int64).reshape(-1)
    if uvtw.size == 4:
        return CrystalDirection.from_miller_bravais(
            tuple(int(value) for value in uvtw),
            phase=phase,
        )
    raise ValueError("Expected an orix Miller direction with uvw or UVTW indices.")
