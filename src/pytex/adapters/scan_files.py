"""Direct readers for vendor EBSD scan files (.ang, .ctf, .oh5/.h5).

These readers parse EDAX/TSL ``.ang`` and Oxford/HKL Channel 5 ``.ctf`` text
scans, and EDAX OIM ``.oh5``/``.h5`` HDF5 scans, into `NormalizedEBSDDataset`
objects (a `CrystalMap` plus an auto-generated `EBSDImportManifest`). The text
readers need nothing beyond NumPy; the HDF5 reader needs the optional ``h5py``
dependency, imported lazily so the module stays importable without it.
`read_scan` picks the right one from a path's extension.

Frame policy: Euler angles are imported as Bunge angles expressed in the
vendor specimen frame; no axis remapping between vendor and PyTex specimen
conventions is applied by the reader itself. The generated manifest records
the source system and reader policy so downstream normalization stays
auditable.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from itertools import pairwise
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from pytex.adapters.ebsd import NormalizedEBSDDataset, _normalize_vendor_payload
from pytex.core import frame_catalog
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase

# TSL/EDAX .ang "Symmetry" codes mapped onto Laue-class Hermann-Mauguin
# symbols. The codes identify the Laue group of the indexed phase.
_TSL_SYMMETRY_CODES = {
    "1": "-1",
    "2": "2/m",
    "20": "2/m",
    "22": "mmm",
    "3": "-3",
    "32": "-3m",
    "4": "4/m",
    "42": "4/mmm",
    "6": "6/m",
    "62": "6/mmm",
    "23": "m-3",
    "43": "m-3m",
}

# Oxford/HKL .ctf phase records carry the Laue group as an integer 1-11.
_HKL_LAUE_CODES = {
    1: "-1",
    2: "2/m",
    3: "mmm",
    4: "4/m",
    5: "4/mmm",
    6: "-3",
    7: "-3m",
    8: "6/m",
    9: "6/mmm",
    10: "m-3",
    11: "m-3m",
}

_ANG_PROPERTY_COLUMNS = (
    (5, "image_quality"),
    (6, "confidence_index"),
    (8, "detector_signal"),
    (9, "fit"),
)

_CTF_PROPERTY_NAMES = {
    "bands": "bands",
    "error": "error",
    "mad": "mean_angular_deviation",
    "bc": "band_contrast",
    "bs": "band_slope",
}


def default_ebsd_frames() -> tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]:
    """The crystal, specimen, and map frames a vendor scan is imported into.

    Built from `pytex.core.frame_catalog` so an imported scan's frames compare
    equal to the same frames anywhere else in the library. No axis remapping
    between vendor and PyTex specimen conventions is applied here; that policy
    is recorded in the generated import manifest.
    """

    return (
        frame_catalog.crystal_frame(),
        frame_catalog.specimen_frame(),
        frame_catalog.map_frame(),
    )


def _readonly_float(values: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class EBSDScanFileResult:
    dataset: NormalizedEBSDDataset
    properties: Mapping[str, np.ndarray] = field(default_factory=dict)
    header_metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        point_count = len(self.dataset.crystal_map.orientations)
        frozen: dict[str, np.ndarray] = {}
        for name, values in dict(self.properties).items():
            array = _readonly_float(values)
            if array.shape != (point_count,):
                raise ValueError(
                    f"EBSDScanFileResult property '{name}' must have one value per map point."
                )
            frozen[str(name)] = array
        object.__setattr__(self, "properties", MappingProxyType(frozen))
        object.__setattr__(
            self,
            "header_metadata",
            MappingProxyType({str(k): str(v) for k, v in dict(self.header_metadata).items()}),
        )
        # Attach the parsed per-point channels to the CrystalMap so downstream
        # `result.crystal_map` carries them directly (e.g. plot_property_map),
        # while `EBSDScanFileResult.properties` remains for back-compat access.
        if frozen:
            map_with_properties = self.dataset.crystal_map.with_properties(frozen)
            object.__setattr__(
                self,
                "dataset",
                replace(self.dataset, crystal_map=map_with_properties),
            )

    @property
    def crystal_map(self) -> Any:
        return self.dataset.crystal_map

    @property
    def manifest(self) -> Any:
        return self.dataset.manifest


@dataclass
class _AngPhaseBlock:
    phase_id: int
    material_name: str | None = None
    formula: str | None = None
    symmetry_code: str | None = None
    lattice_constants: str | None = None


def _parse_ang_header(
    lines: list[str],
) -> tuple[dict[str, str], list[_AngPhaseBlock]]:
    header: dict[str, str] = {}
    blocks: list[_AngPhaseBlock] = []
    current: _AngPhaseBlock | None = None
    block_keys = {"materialname", "formula", "symmetry", "latticeconstants"}
    for raw_line in lines:
        content = raw_line[1:].strip()
        if not content:
            continue
        if ":" in content:
            key, _, value = content.partition(":")
            key = key.strip()
            value = value.strip()
        else:
            parts = content.split(None, 1)
            key = parts[0].strip()
            value = parts[1].strip() if len(parts) > 1 else ""
        lowered = key.lower()
        if lowered == "phase":
            try:
                phase_id = int(value)
            except ValueError as error:
                raise ValueError(f"Invalid .ang phase block id '{value}'.") from error
            current = _AngPhaseBlock(phase_id=phase_id)
            blocks.append(current)
            continue
        if lowered in block_keys:
            if current is None:
                current = _AngPhaseBlock(phase_id=1)
                blocks.append(current)
            if lowered == "materialname":
                current.material_name = value
            elif lowered == "formula":
                current.formula = value
            elif lowered == "symmetry":
                current.symmetry_code = value
            elif lowered == "latticeconstants":
                current.lattice_constants = value
            continue
        if value:
            header[key] = value
    return header, blocks


def _point_group_for_tsl_code(code: str | None, *, phase_name: str, source: str = ".ang") -> str:
    """The Laue-class symbol a TSL symmetry code names.

    Both EDAX formats carry the same code set - the ``.ang`` ``# Symmetry``
    header line and the OIM HDF5 ``LGsymID`` dataset - so one table serves both
    and a code PyTex does not know fails the same way in either.
    """

    if code is None:
        raise ValueError(f"{source} phase '{phase_name}' does not declare a symmetry code.")
    normalized = code.strip()
    point_group = _TSL_SYMMETRY_CODES.get(normalized)
    if point_group is None:
        supported = ", ".join(sorted(_TSL_SYMMETRY_CODES))
        raise ValueError(
            f"Unsupported TSL symmetry code '{code}' for phase '{phase_name}'. "
            f"Supported codes: {supported}."
        )
    return point_group


def read_ang(
    path: str | Path,
    *,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    phases: dict[int | str, Phase] | tuple[Phase, ...] | list[Phase] | None = None,
) -> EBSDScanFileResult:
    """Read an EDAX/TSL ``.ang`` scan into a normalized EBSD dataset.

    Use this direct pure-Python reader when orientations, map coordinates,
    phase declarations, and scalar quality channels are needed without a live
    vendor or optional EBSD package. ``SqrGrid`` metadata becomes a rectangular
    ``grid_shape``. ``HexGrid`` metadata becomes an explicit staggered
    ``grid_kind`` with alternating row lengths and six-neighbour topology.

    Parameters
    ----------
    path : str or Path
        EDAX/TSL text ``.ang`` file. Euler angles are interpreted as Bunge
        radians, following the format convention.
    crystal_frame, specimen_frame, map_frame : ReferenceFrame, optional
        Explicit canonical frames. EBSD defaults are constructed when omitted.
    phase : Phase, optional
        Full phase semantics for a single-phase file.
    phases : mapping or sequence of Phase, optional
        Full phase semantics resolved against a multiphase file's ids or names.

    Returns
    -------
    EBSDScanFileResult
        A ``CrystalMap``, generated import manifest, immutable per-point IQ/CI/
        fit channels, and the parsed header metadata.

    Notes
    -----
    If multiphase filtering drops points, the reader deliberately drops the
    logical grid claim too: a row with holes is no longer the complete vendor
    topology. Hexagonal curvature/GND stencils remain unsupported even though
    graph-backed KAM and segmentation preserve the six-neighbour scan.
    """

    file_path = Path(path)
    text = file_path.read_text(encoding="utf-8", errors="replace")
    header_lines = [line for line in text.splitlines() if line.startswith("#")]
    data_lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not data_lines:
        raise ValueError(f".ang file '{file_path}' contains no data rows.")
    header, blocks = _parse_ang_header(header_lines)
    grid_kind = header.get("GRID", "SqrGrid").strip()
    normalized_grid_kind = grid_kind.lower()
    if normalized_grid_kind not in {"sqrgrid", "hexgrid"}:
        raise ValueError(
            f".ang grid type '{grid_kind}' is not supported; expected SqrGrid or HexGrid."
        )
    if not blocks:
        raise ValueError(f".ang file '{file_path}' declares no phase information.")

    data = np.loadtxt(data_lines, dtype=np.float64, ndmin=2)
    if data.shape[1] < 5:
        raise ValueError(".ang data rows require at least phi1, PHI, phi2, x, y columns.")
    euler_rad = data[:, 0:3]
    coordinates = data[:, 3:5]
    phase_column = (
        data[:, 7].astype(np.int64) if data.shape[1] > 7 else np.zeros(len(data), dtype=np.int64)
    )

    declared_ids = {block.phase_id for block in blocks}
    dropped_points = 0
    if len(blocks) == 1:
        keep_mask = np.ones(len(data), dtype=bool)
    else:
        keep_mask = np.isin(phase_column, sorted(declared_ids))
        dropped_points = int(np.count_nonzero(~keep_mask))
        if not np.any(keep_mask):
            raise ValueError(".ang data phase column does not reference any declared phase block.")

    metadata = {
        "reader": "pytex.adapters.scan_files.read_ang",
        "grid": grid_kind,
        "dropped_unindexed_points": str(dropped_points),
        "frame_policy": "vendor specimen frame, no axis remapping applied",
    }
    for key in ("XSTEP", "YSTEP", "NCOLS_ODD", "NCOLS_EVEN", "NROWS", "OPERATOR", "SAMPLEID"):
        if key in header:
            metadata[key.lower()] = header[key]

    payload: dict[str, Any] = {
        "coordinates": coordinates[keep_mask],
        "euler_angles": euler_rad[keep_mask],
        "orientation_convention": "bunge",
        "angle_unit": "radian",
        "source_file": str(file_path),
        "metadata": metadata,
    }

    step_x = float(header.get("XSTEP", "0") or 0.0)
    step_y = float(header.get("YSTEP", "0") or 0.0)
    if step_x > 0.0 and step_y > 0.0:
        payload["step_sizes"] = (step_x, step_y)
    n_rows = int(header.get("NROWS", "0") or 0)
    n_cols_odd = int(header.get("NCOLS_ODD", "0") or 0)
    n_cols_even = int(header.get("NCOLS_EVEN", "0") or 0)
    kept_count = int(np.count_nonzero(keep_mask))
    if normalized_grid_kind == "sqrgrid":
        if n_rows > 0 and n_cols_odd > 0 and n_rows * n_cols_odd == kept_count:
            payload["grid_shape"] = (n_rows, n_cols_odd)
    elif n_rows <= 0 or n_cols_odd <= 0 or n_cols_even <= 0:
        raise ValueError(
            ".ang HexGrid scans require positive NROWS, NCOLS_ODD, and NCOLS_EVEN headers."
        )
    else:
        row_lengths = tuple(n_cols_odd if row % 2 == 0 else n_cols_even for row in range(n_rows))
        expected_count = sum(row_lengths)
        if expected_count != len(data):
            raise ValueError(
                ".ang HexGrid row metadata does not match the number of data rows "
                f"(expected {expected_count}, found {len(data)})."
            )
        if kept_count == len(data):
            payload["grid_kind"] = "hexagonal"
            payload["row_lengths"] = row_lengths

    if len(blocks) == 1:
        block = blocks[0]
        name = block.material_name or block.formula or f"phase_{block.phase_id}"
        payload["phase_name"] = name
        payload["point_group"] = _point_group_for_tsl_code(block.symmetry_code, phase_name=name)
    else:
        payload["phases"] = [
            {
                "phase_id": str(block.phase_id),
                "name": block.material_name or block.formula or f"phase_{block.phase_id}",
                "point_group": _point_group_for_tsl_code(
                    block.symmetry_code,
                    phase_name=block.material_name or f"phase_{block.phase_id}",
                ),
            }
            for block in blocks
        ]
        payload["phase_ids"] = phase_column[keep_mask]

    resolved_crystal, resolved_specimen, resolved_map = default_ebsd_frames()
    dataset = _normalize_vendor_payload(
        payload,
        source_system="tsl_ang",
        crystal_frame=crystal_frame or resolved_crystal,
        specimen_frame=specimen_frame or resolved_specimen,
        map_frame=map_frame or resolved_map,
        angle_key_candidates=("euler_angles",),
        phase=phase,
        phases=phases,
    )
    properties = {
        name: data[keep_mask, column]
        for column, name in _ANG_PROPERTY_COLUMNS
        if data.shape[1] > column
    }
    return EBSDScanFileResult(
        dataset=dataset,
        properties=properties,
        header_metadata=header,
    )


def _parse_ctf_phase_line(line: str, *, phase_id: int) -> dict[str, str]:
    fields = [item.strip() for item in line.split("\t")]
    if len(fields) < 4:
        raise ValueError(
            f".ctf phase record {phase_id} requires lattice, angles, name, and Laue group."
        )
    try:
        laue_code = int(float(fields[3]))
    except ValueError as error:
        raise ValueError(
            f".ctf phase record {phase_id} has a non-numeric Laue group '{fields[3]}'."
        ) from error
    point_group = _HKL_LAUE_CODES.get(laue_code)
    if point_group is None:
        raise ValueError(
            f".ctf phase record {phase_id} declares unsupported Laue code {laue_code}."
        )
    return {
        "phase_id": str(phase_id),
        "name": fields[2] or f"phase_{phase_id}",
        "point_group": point_group,
    }


def read_ctf(
    path: str | Path,
    *,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    phases: dict[int | str, Phase] | tuple[Phase, ...] | list[Phase] | None = None,
) -> EBSDScanFileResult:
    """Read an Oxford/HKL Channel 5 .ctf scan into a normalized EBSD dataset.

    Non-indexed points (phase id 0) are dropped; the count is recorded in the
    generated manifest metadata.
    """

    file_path = Path(path)
    lines = file_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines or not lines[0].strip().lower().startswith("channel text file"):
        raise ValueError(f"'{file_path}' does not start with a Channel Text File header.")

    header: dict[str, str] = {}
    phase_records: list[dict[str, str]] = []
    data_start: int | None = None
    column_names: list[str] = []
    index = 1
    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if not stripped:
            index += 1
            continue
        first_token = stripped.split("\t")[0].split()[0] if stripped else ""
        if first_token.lower() == "phases":
            parts = stripped.split()
            try:
                declared = int(parts[-1])
            except ValueError as error:
                raise ValueError(".ctf 'Phases' line must end with a phase count.") from error
            for offset in range(1, declared + 1):
                if index + offset >= len(lines):
                    raise ValueError(".ctf file ends before all declared phases are listed.")
                phase_records.append(_parse_ctf_phase_line(lines[index + offset], phase_id=offset))
            index += declared + 1
            continue
        if first_token.lower() == "phase" and "euler1" in stripped.lower():
            column_names = [name.strip().lower() for name in stripped.split("\t") if name.strip()]
            data_start = index + 1
            break
        key, _, value = stripped.partition("\t")
        header[key.strip()] = value.strip()
        index += 1

    if data_start is None:
        raise ValueError(f"'{file_path}' contains no .ctf data header row.")
    if not phase_records:
        raise ValueError(f"'{file_path}' declares no phases.")

    data_lines = [line for line in lines[data_start:] if line.strip()]
    if not data_lines:
        raise ValueError(f"'{file_path}' contains no data rows.")
    data = np.loadtxt(data_lines, dtype=np.float64, ndmin=2)
    if data.shape[1] < len(column_names):
        raise ValueError(".ctf data rows have fewer columns than the data header declares.")

    def column(name: str) -> np.ndarray:
        try:
            position = column_names.index(name)
        except ValueError as error:
            raise ValueError(f".ctf data header is missing the '{name}' column.") from error
        return data[:, position]

    phase_column = column("phase").astype(np.int64)
    keep_mask = phase_column > 0
    dropped_points = int(np.count_nonzero(~keep_mask))
    if not np.any(keep_mask):
        raise ValueError(f"'{file_path}' contains only non-indexed points.")

    coordinates = np.column_stack([column("x"), column("y")])[keep_mask]
    euler_deg = np.column_stack([column("euler1"), column("euler2"), column("euler3")])[keep_mask]

    metadata = {
        "reader": "pytex.adapters.scan_files.read_ctf",
        "dropped_unindexed_points": str(dropped_points),
        "frame_policy": "vendor specimen frame, no axis remapping applied",
    }
    for key in ("JobMode", "XCells", "YCells", "XStep", "YStep", "Mag", "Prj"):
        if key in header:
            metadata[key.lower()] = header[key]

    payload: dict[str, Any] = {
        "coordinates": coordinates,
        "euler_angles_deg": euler_deg,
        "orientation_convention": "bunge",
        "angle_unit": "degree",
        "source_file": str(file_path),
        "metadata": metadata,
    }

    step_x = float(header.get("XStep", "0") or 0.0)
    step_y = float(header.get("YStep", "0") or 0.0)
    if step_x > 0.0 and step_y > 0.0:
        payload["step_sizes"] = (step_x, step_y)
    x_cells = int(header.get("XCells", "0") or 0)
    y_cells = int(header.get("YCells", "0") or 0)
    if x_cells > 0 and y_cells > 0 and x_cells * y_cells == int(np.count_nonzero(keep_mask)):
        payload["grid_shape"] = (y_cells, x_cells)

    if len(phase_records) == 1:
        payload["phase_name"] = phase_records[0]["name"]
        payload["point_group"] = phase_records[0]["point_group"]
    else:
        payload["phases"] = phase_records
        payload["phase_ids"] = phase_column[keep_mask]

    resolved_crystal, resolved_specimen, resolved_map = default_ebsd_frames()
    dataset = _normalize_vendor_payload(
        payload,
        source_system="hkl_ctf",
        crystal_frame=crystal_frame or resolved_crystal,
        specimen_frame=specimen_frame or resolved_specimen,
        map_frame=map_frame or resolved_map,
        angle_key_candidates=("euler_angles_deg",),
        phase=phase,
        phases=phases,
    )
    properties = {
        friendly: column(short)[keep_mask]
        for short, friendly in _CTF_PROPERTY_NAMES.items()
        if short in column_names
    }
    return EBSDScanFileResult(
        dataset=dataset,
        properties=properties,
        header_metadata=header,
    )


# EDAX OIM writes the same HDF5 container under two extensions: `.oh5` for a
# scan saved by OIM Analysis and `.h5` for the same data exported as h5ebsd.
# One reader serves both.
# Top-level members of an OIM HDF5 file that are file metadata rather than scans.
_OH5_NON_SCAN_MEMBERS = frozenset({"Manufacturer", "Version"})

# Per-point channels that carry the map's geometry and indexing rather than a
# measured quality value, so they are not exposed as property channels.
_OH5_STRUCTURAL_FIELDS = frozenset({"Phi1", "Phi", "Phi2", "X Position", "Y Position", "Phase"})

# EDAX channel names mapped onto the names `read_ang` gives the same
# quantities, so a workflow cannot tell which of the two exports fed it.
_OH5_PROPERTY_FIELDS = {
    "IQ": "image_quality",
    "CI": "confidence_index",
    "SEM Signal": "detector_signal",
    "Fit": "fit",
}


def _h5py() -> Any:
    """The HDF5 binding the OIM reader uses.

    Required, not optional; imported at call time so the .ang and .ctf readers,
    which are pure NumPy, do not pay for it.
    """

    import h5py

    return h5py


def _oh5_value(group: Any, name: str) -> Any:
    """The single value held by an OIM header dataset.

    OIM writes header scalars as length-1 datasets rather than as HDF5 scalars,
    so a plain ``dataset[()]`` gives back an array. This unwraps both spellings.
    """

    dataset = group.get(name)
    if dataset is None:
        return None
    value = dataset[()]
    if np.ndim(value) == 0:
        return value
    flat = np.ravel(value)
    return flat[0] if flat.size else None


def _oh5_number(group: Any, name: str) -> float | None:
    value = _oh5_value(group, name)
    if value is None or isinstance(value, bytes | str):
        return None
    return float(value)


def _oh5_text(group: Any, name: str) -> str | None:
    value = _oh5_value(group, name)
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    return str(value).strip()


def _oh5_scan_group(handle: Any, scan: str | None, *, file_path: Path) -> str:
    """Name the scan group to read, of possibly several in one file."""

    h5py = _h5py()
    candidates: list[str] = [
        str(name)
        for name in handle
        if name not in _OH5_NON_SCAN_MEMBERS
        and isinstance(handle[name], h5py.Group)
        and "EBSD" in handle[name]
    ]
    if scan is not None:
        if scan not in candidates:
            available = ", ".join(candidates) or "none"
            raise ValueError(
                f"'{file_path}' has no EBSD scan named '{scan}'. Available scans: {available}."
            )
        return scan
    if not candidates:
        raise ValueError(f"'{file_path}' contains no group holding an EBSD scan.")
    return candidates[0]


def _oh5_phase_blocks(header: Any, *, file_path: Path) -> list[dict[str, str]]:
    """The declared phases, in the header's own numbering.

    OIM names the phase subgroups ``1``, ``2``, ... and records the Laue class
    as ``LGsymID``, using the same code set the ``.ang`` ``# Symmetry`` line
    uses - so the two EDAX formats resolve symmetry through one table.
    """

    group = header.get("Phase")
    if group is None:
        raise ValueError(f"'{file_path}' declares no phases in its EBSD header.")

    def sort_key(name: str) -> tuple[int, str]:
        return (int(name), "") if name.strip().isdigit() else (2**31, name)

    blocks: list[dict[str, str]] = []
    for position, name in enumerate(sorted(group.keys(), key=sort_key), start=1):
        phase_group = group[name]
        phase_id = int(name) if name.strip().isdigit() else position
        material = _oh5_text(phase_group, "MaterialName")
        formula = _oh5_text(phase_group, "Formula")
        phase_name = material or formula or f"phase_{phase_id}"
        code = _oh5_number(phase_group, "LGsymID")
        symmetry_code = None if code is None else str(int(code))
        blocks.append(
            {
                "phase_id": str(phase_id),
                "name": phase_name,
                "point_group": _point_group_for_tsl_code(
                    symmetry_code, phase_name=phase_name, source="OIM HDF5"
                ),
            }
        )
    if not blocks:
        raise ValueError(f"'{file_path}' declares no phases in its EBSD header.")
    return blocks


def _oh5_phase_column_base(phase_index: np.ndarray, *, phase_count: int) -> int:
    """Whether the ``Phase`` channel counts phases from zero or from one.

    OIM Analysis 8 writes a zero-based index into the header phase list with
    ``-1`` for an unindexed point, which is what the ``.ang`` ``Phase index``
    column carries for the same scan. Older writers use a one-based index with
    ``0`` for unindexed. The sentinel settles it: a negative value can only be
    the zero-based spelling, and so can a zero that is inside the declared
    range.
    """

    if phase_index.size and int(phase_index.min()) < 0:
        return 0
    indexed = phase_index[phase_index >= 0]
    if indexed.size and int(indexed.min()) == 0 and int(indexed.max()) <= phase_count - 1:
        return 0
    return 1


def _oh5_property_name(field: str) -> str:
    mapped = _OH5_PROPERTY_FIELDS.get(field)
    if mapped is not None:
        return mapped
    return "_".join(field.strip().lower().replace("-", " ").replace("_", " ").split())


def _oh5_hex_row_lengths(
    y_positions: np.ndarray,
    *,
    step_y: float,
    rows: int,
    file_path: Path,
) -> tuple[int, ...]:
    """Recover a hexagonal scan's staggered row lengths from its Y coordinates.

    A ``.ang`` header states ``NCOLS_ODD`` and ``NCOLS_EVEN`` outright; the OIM
    HDF5 header does not, so the rows are counted from the scan itself. Points
    are written in scan order and a row holds one Y value, so a row boundary is
    a change in Y - with half a step as the tolerance, which is the largest
    value that can neither merge two rows nor split one.
    """

    if step_y <= 0.0:
        raise ValueError(
            f"'{file_path}' is a HexGrid scan without a positive Step Y, so its staggered "
            "row lengths cannot be recovered."
        )
    boundaries = np.flatnonzero(np.abs(np.diff(y_positions)) > 0.5 * step_y) + 1
    edges = [0, *(int(value) for value in boundaries), len(y_positions)]
    row_lengths = tuple(int(end - start) for start, end in pairwise(edges))
    if rows > 0 and len(row_lengths) != rows:
        raise ValueError(
            f"'{file_path}' is a HexGrid scan whose Y positions describe {len(row_lengths)} "
            f"rows while its header declares {rows}."
        )
    if len(set(row_lengths)) > 2:
        raise ValueError(
            f"'{file_path}' is a HexGrid scan whose rows have {len(set(row_lengths))} distinct "
            "lengths; a staggered grid alternates between two."
        )
    return row_lengths


def read_oh5(
    path: str | Path,
    *,
    scan: str | None = None,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    phases: dict[int | str, Phase] | tuple[Phase, ...] | list[Phase] | None = None,
) -> EBSDScanFileResult:
    """Read an EDAX OIM HDF5 (``.oh5`` or ``.h5``) scan into a normalized dataset.

    Use this when the scan was saved from OIM Analysis rather than exported as
    text. ``.oh5`` and ``.h5`` are the same container under two extensions, and
    the same scan read from either one - or from the ``.ang`` OIM exports beside
    it - gives the same orientations, phases, grid and quality channels: the
    Euler angles are Bunge radians in both formats, and the HDF5 ``LGsymID``
    holds the same TSL symmetry code as the ``.ang`` ``# Symmetry`` line.

    The reason to prefer the HDF5 export is the channels: every per-point scalar
    the file carries is read, not only the columns a ``.ang`` row has room for.
    The four EDAX channels that ``.ang`` also carries keep the names `read_ang`
    gives them (``image_quality``, ``confidence_index``, ``detector_signal``,
    ``fit``); every other channel - vendor PRIAS signals, the ``valid`` flag,
    anything a processing tool wrote back - is exposed under its own lower-cased
    name.

    Parameters
    ----------
    path : str or Path
        The ``.oh5`` or ``.h5`` file. Requires the optional ``h5py`` dependency.
    scan : str, optional
        Which scan group to read, for a file holding more than one. The first
        group carrying an ``EBSD`` record is read when omitted.
    crystal_frame, specimen_frame, map_frame : ReferenceFrame, optional
        Explicit canonical frames. EBSD defaults are constructed when omitted.
    phase : Phase, optional
        Full phase semantics for a single-phase file.
    phases : mapping or sequence of Phase, optional
        Full phase semantics resolved against a multiphase file's ids or names.

    Returns
    -------
    EBSDScanFileResult
        A ``CrystalMap``, generated import manifest, the immutable per-point
        channels, and the header values worth carrying as text metadata.

    Notes
    -----
    Point keeping follows `read_ang`, so that the two EDAX formats import the
    same scan identically: a single-phase file keeps every measured point,
    including the unindexed ones - read the ``confidence_index`` channel before
    believing an orientation - while a multiphase file keeps only the points
    referencing a declared phase and records how many it dropped.

    ``HexGrid`` scans state only ``nRows`` and ``nColumns`` in the header, with
    no odd/even row split, so the staggered row lengths are recovered from the
    ``Y Position`` channel and cross-checked against ``nRows``.

    See Also
    --------
    read_ang : The text export of the same EDAX scan.
    read_scan : Pick the reader for a path by its extension.
    """

    h5py = _h5py()
    file_path = Path(path)
    header_metadata: dict[str, str] = {}
    channels: dict[str, np.ndarray] = {}
    with h5py.File(file_path, "r") as handle:
        manufacturer = _oh5_text(handle, "Manufacturer")
        version = _oh5_text(handle, "Version")
        scan_name = _oh5_scan_group(handle, scan, file_path=file_path)
        ebsd = handle[f"{scan_name}/EBSD"]
        if "Header" not in ebsd or "Data" not in ebsd:
            raise ValueError(f"'{file_path}' scan '{scan_name}' has no EBSD Header and Data pair.")
        header = ebsd["Header"]
        data = ebsd["Data"]

        grid_kind = _oh5_text(header, "Grid Type") or "SqrGrid"
        normalized_grid_kind = grid_kind.lower()
        if normalized_grid_kind not in {"sqrgrid", "hexgrid"}:
            raise ValueError(
                f"OIM HDF5 grid type '{grid_kind}' is not supported; expected SqrGrid or HexGrid."
            )
        n_columns = _oh5_number(header, "nColumns")
        n_rows = _oh5_number(header, "nRows")
        step_x = _oh5_number(header, "Step X") or 0.0
        step_y = _oh5_number(header, "Step Y") or 0.0

        blocks = _oh5_phase_blocks(header, file_path=file_path)

        for field in ("Phi1", "Phi", "Phi2", "X Position", "Y Position"):
            if field not in data:
                raise ValueError(
                    f"'{file_path}' scan '{scan_name}' has no '{field}' channel; it is not an "
                    "indexed OIM EBSD scan."
                )
        euler_rad = np.column_stack(
            [
                np.asarray(data["Phi1"][()], dtype=np.float64).ravel(),
                np.asarray(data["Phi"][()], dtype=np.float64).ravel(),
                np.asarray(data["Phi2"][()], dtype=np.float64).ravel(),
            ]
        )
        coordinates = np.column_stack(
            [
                np.asarray(data["X Position"][()], dtype=np.float64).ravel(),
                np.asarray(data["Y Position"][()], dtype=np.float64).ravel(),
            ]
        )
        point_count = len(euler_rad)
        if point_count == 0:
            raise ValueError(f"'{file_path}' scan '{scan_name}' contains no measured points.")
        if len(coordinates) != point_count:
            raise ValueError(
                f"'{file_path}' scan '{scan_name}' has {len(coordinates)} positions for "
                f"{point_count} orientations."
            )

        if "Phase" in data:
            phase_index = np.asarray(data["Phase"][()], dtype=np.int64).ravel()
            if len(phase_index) != point_count:
                raise ValueError(
                    f"'{file_path}' scan '{scan_name}' has a Phase channel of "
                    f"{len(phase_index)} values for {point_count} points."
                )
        else:
            phase_index = np.zeros(point_count, dtype=np.int64)

        for field, dataset in data.items():
            if field in _OH5_STRUCTURAL_FIELDS or not isinstance(dataset, h5py.Dataset):
                continue
            if dataset.dtype.kind not in {"b", "i", "u", "f"}:
                continue
            if dataset.ndim != 1 or dataset.shape[0] != point_count:
                continue
            channels[_oh5_property_name(field)] = np.asarray(dataset[()], dtype=np.float64)

        for key in ("Sample Tilt", "Working Distance", "Voltage[kV]", "Camera Elevation Angle"):
            number = _oh5_number(header, key)
            if number is not None:
                header_metadata[key] = repr(number)
        notes = _oh5_text(header, "Notes")
        if notes:
            header_metadata["Notes"] = notes

    base = _oh5_phase_column_base(phase_index, phase_count=len(blocks))
    phase_column = phase_index + 1 if base == 0 else phase_index
    declared_ids = {int(block["phase_id"]) for block in blocks}
    if len(blocks) == 1:
        keep_mask = np.ones(point_count, dtype=bool)
    else:
        keep_mask = np.isin(phase_column, sorted(declared_ids))
        if not np.any(keep_mask):
            raise ValueError(f"'{file_path}' Phase channel references none of its declared phases.")
    dropped_points = int(np.count_nonzero(~keep_mask))
    kept_count = int(np.count_nonzero(keep_mask))

    metadata = {
        "reader": "pytex.adapters.scan_files.read_oh5",
        "scan": scan_name,
        "grid": grid_kind,
        "phase_column_base": str(base),
        "dropped_unindexed_points": str(dropped_points),
        "frame_policy": "vendor specimen frame, no axis remapping applied",
    }
    if manufacturer:
        metadata["manufacturer"] = manufacturer
    if version:
        metadata["version"] = version
    if n_rows is not None:
        metadata["nrows"] = str(int(n_rows))
    if n_columns is not None:
        metadata["ncolumns"] = str(int(n_columns))
    if step_x > 0.0:
        metadata["xstep"] = repr(step_x)
    if step_y > 0.0:
        metadata["ystep"] = repr(step_y)

    payload: dict[str, Any] = {
        "coordinates": coordinates[keep_mask],
        "euler_angles": euler_rad[keep_mask],
        "orientation_convention": "bunge",
        "angle_unit": "radian",
        "source_file": str(file_path),
        "metadata": metadata,
    }
    if step_x > 0.0 and step_y > 0.0:
        payload["step_sizes"] = (step_x, step_y)

    rows = int(n_rows or 0)
    columns = int(n_columns or 0)
    if normalized_grid_kind == "sqrgrid":
        if rows > 0 and columns > 0 and rows * columns == kept_count:
            payload["grid_shape"] = (rows, columns)
    elif kept_count == point_count:
        payload["grid_kind"] = "hexagonal"
        payload["row_lengths"] = _oh5_hex_row_lengths(
            coordinates[:, 1], step_y=step_y, rows=rows, file_path=file_path
        )

    if len(blocks) == 1:
        payload["phase_name"] = blocks[0]["name"]
        payload["point_group"] = blocks[0]["point_group"]
    else:
        payload["phases"] = blocks
        payload["phase_ids"] = phase_column[keep_mask]

    resolved_crystal, resolved_specimen, resolved_map = default_ebsd_frames()
    dataset = _normalize_vendor_payload(
        payload,
        source_system="edax_oh5",
        crystal_frame=crystal_frame or resolved_crystal,
        specimen_frame=specimen_frame or resolved_specimen,
        map_frame=map_frame or resolved_map,
        angle_key_candidates=("euler_angles",),
        phase=phase,
        phases=phases,
    )
    properties = {name: values[keep_mask] for name, values in sorted(channels.items())}
    return EBSDScanFileResult(
        dataset=dataset,
        properties=properties,
        header_metadata=header_metadata,
    )


#: The scan file extensions PyTex reads, and the reader each one dispatches to.
_SCAN_READERS: dict[str, str] = {
    ".ang": "read_ang",
    ".ctf": "read_ctf",
    ".oh5": "read_oh5",
    ".h5": "read_oh5",
}

#: Every scan file extension `read_scan` accepts, in dispatch order.
SCAN_FILE_SUFFIXES: tuple[str, ...] = tuple(_SCAN_READERS)


def scan_reader_for(path: str | Path) -> Any:
    """The reader that opens a scan file of this path's kind.

    Use it when the format is only known at run time - a file a user chose, a
    directory being walked - so extension dispatch is written once, here, rather
    than once per caller.

    Raises
    ------
    ValueError
        If the extension is not one of `SCAN_FILE_SUFFIXES`.
    """

    suffix = Path(path).suffix.lower()
    name = _SCAN_READERS.get(suffix)
    if name is None:
        supported = ", ".join(SCAN_FILE_SUFFIXES)
        raise ValueError(
            f"'{Path(path).name}' is not a scan file PyTex reads: its extension is "
            f"{suffix or 'missing'}. Supported extensions: {supported}."
        )
    return globals()[name]


def read_scan(path: str | Path, **kwargs: Any) -> EBSDScanFileResult:
    """Read an EBSD scan file of any supported vendor format.

    Dispatches on the extension - ``.ang`` to `read_ang`, ``.ctf`` to
    `read_ctf`, ``.oh5`` and ``.h5`` to `read_oh5` - and passes the keyword
    arguments straight through, so anything the chosen reader accepts (explicit
    frames, phase semantics) works here too.
    """

    reader = scan_reader_for(path)
    result: EBSDScanFileResult = reader(path, **kwargs)
    return result


__all__ = [
    "SCAN_FILE_SUFFIXES",
    "EBSDScanFileResult",
    "default_ebsd_frames",
    "read_ang",
    "read_ctf",
    "read_oh5",
    "read_scan",
    "scan_reader_for",
]
