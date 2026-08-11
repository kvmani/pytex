"""Tabular and file exports for composite SAED patterns.

A simulated pattern is only useful once it leaves the process: as a reflection
table that can be pasted into a paper, as a figure, and as a manifest recording
exactly how it was produced. This module is that boundary. The simulation
itself lives in :mod:`pytex.diffraction.composite` and the renderer in
:mod:`pytex.plotting.composite_saed`; nothing here recomputes crystallography.

Everything a table row reports is taken from the `SpotTable` the engine
produced, so the exported numbers are the rendered numbers by construction
rather than by a parallel calculation that could drift.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.composite import (
    CompositeSAEDPattern,
    SpotCoincidenceReport,
    find_spot_coincidences,
)

#: Schema identifier of the composite-SAED manifest payload.
COMPOSITE_SAED_MANIFEST_SCHEMA = "pytex.composite_saed_manifest/1"

#: Schema identifier of the reflection-table payload.
COMPOSITE_REFLECTION_TABLE_SCHEMA = "pytex.composite_reflection_table/1"

#: Column order of the exported reflection table.
#:
#: Fixed and public because downstream tools and pinned test fixtures read it
#: positionally as often as by name. Appending is a compatible change; removing
#: or reordering is not.
REFLECTION_TABLE_COLUMNS: tuple[str, ...] = (
    "source",
    "variant",
    "phase",
    "h",
    "k",
    "l",
    "hkl_label",
    "zone_axis",
    "d_angstrom",
    "g_inv_angstrom",
    "detector_x_mm",
    "detector_y_mm",
    "detector_r_mm",
    "excitation_error_inv_angstrom",
    "structure_factor_amplitude",
    "relative_intensity",
    "double_diffraction",
    "double_diffraction_origin",
)


@dataclass(frozen=True, slots=True)
class ReflectionTableRow:
    """One rendered reflection of a composite pattern, with its provenance.

    ``source`` is ``"parent"`` or ``"variant"``; ``variant_index`` is ``None``
    for the parent. ``relative_intensity`` is normalized **within its own
    sub-pattern**, because kinematic theory does not define an intensity ratio
    between two different phases — comparing intensities across sources in this
    table is not meaningful, and `ReflectionTable.describe` says so.

    ``double_diffraction_origin`` is non-empty exactly for a reflection that is
    kinematically **forbidden** and appears only because a diffracted beam
    re-diffracts; it names the path, as ``g1 + g2 = g``. Such a row reports an
    observability estimate, not a kinematic intensity — its structure-factor
    amplitude is (near) zero, which is why it is marked.
    """

    source: str
    variant_index: int | None
    phase_name: str
    hkl: tuple[int, int, int]
    hkl_label: str
    zone_axis_label: str
    d_angstrom: float
    g_inv_angstrom: float
    detector_mm: tuple[float, float]
    excitation_error_inv_angstrom: float
    structure_factor_amplitude: float
    relative_intensity: float
    double_diffraction_origin: str = ""

    def __post_init__(self) -> None:
        if self.source not in {"parent", "variant"}:
            raise ValueError("ReflectionTableRow.source must be 'parent' or 'variant'.")
        if self.source == "parent" and self.variant_index is not None:
            raise ValueError("A parent row must not carry a variant index.")
        if self.source == "variant" and self.variant_index is None:
            raise ValueError("A variant row must carry a variant index.")
        if not np.isfinite(self.d_angstrom) or self.d_angstrom <= 0.0:
            raise ValueError("d_angstrom must be finite and positive.")
        if not np.isfinite(self.g_inv_angstrom) or self.g_inv_angstrom <= 0.0:
            raise ValueError("g_inv_angstrom must be finite and positive.")

    @property
    def detector_radius_mm(self) -> float:
        return float(np.hypot(*self.detector_mm))

    @property
    def is_double_diffraction(self) -> bool:
        """Whether this row is a forbidden reflection revived by double diffraction."""

        return bool(self.double_diffraction_origin)

    def as_record(self) -> dict[str, Any]:
        """The row as a flat dictionary keyed by `REFLECTION_TABLE_COLUMNS`."""

        return {
            "source": self.source,
            "variant": "" if self.variant_index is None else self.variant_index,
            "phase": self.phase_name,
            "h": self.hkl[0],
            "k": self.hkl[1],
            "l": self.hkl[2],
            "hkl_label": self.hkl_label,
            "zone_axis": self.zone_axis_label,
            "d_angstrom": self.d_angstrom,
            "g_inv_angstrom": self.g_inv_angstrom,
            "detector_x_mm": self.detector_mm[0],
            "detector_y_mm": self.detector_mm[1],
            "detector_r_mm": self.detector_radius_mm,
            "excitation_error_inv_angstrom": self.excitation_error_inv_angstrom,
            "structure_factor_amplitude": self.structure_factor_amplitude,
            "relative_intensity": self.relative_intensity,
            "double_diffraction": self.is_double_diffraction,
            "double_diffraction_origin": self.double_diffraction_origin,
        }


@dataclass(frozen=True, slots=True)
class ReflectionTable:
    """Every rendered reflection of a composite pattern, as one table.

    The tabular counterpart of the rendered figure: parent rows first, then each
    variant in the composite's own order, each sub-pattern keeping the engine's
    deterministic sort (decreasing intensity, then detector radius, then
    lexicographic ``hkl``). Exports through `to_csv`, `to_markdown`,
    `to_records` and `to_json_dict`.
    """

    relationship_name: str
    parent_phase_name: str
    child_phase_name: str
    parent_zone_axis_label: str
    camera_constant_mm_angstrom: float
    beam_energy_kev: float
    wavelength_angstrom: float
    centering_audit: tuple[tuple[str, str, bool], ...]
    rows: tuple[ReflectionTableRow, ...]
    intensity_threshold: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", tuple(self.rows))
        object.__setattr__(self, "centering_audit", tuple(self.centering_audit))
        if not 0.0 <= self.intensity_threshold < 1.0:
            raise ValueError("intensity_threshold must lie in the interval [0, 1).")

    def __len__(self) -> int:
        return len(self.rows)

    @property
    def variant_indices(self) -> tuple[int, ...]:
        seen: list[int] = []
        for row in self.rows:
            if row.variant_index is not None and row.variant_index not in seen:
                seen.append(row.variant_index)
        return tuple(seen)

    def rows_for_source(self, source: str, *, variant_index: int | None = None) -> tuple[
        ReflectionTableRow, ...
    ]:
        """Rows of one sub-pattern: the parent, or one variant."""

        return tuple(
            row
            for row in self.rows
            if row.source == source
            and (variant_index is None or row.variant_index == variant_index)
        )

    def to_records(self) -> list[dict[str, Any]]:
        """One flat dictionary per row, in `REFLECTION_TABLE_COLUMNS` order."""

        return [row.as_record() for row in self.rows]

    def to_csv(self, path: str | Path) -> Path:
        """Write the table as UTF-8 CSV with a header row and return the path."""

        output = Path(path)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(REFLECTION_TABLE_COLUMNS))
            writer.writeheader()
            writer.writerows(self.to_records())
        return output

    def to_markdown(self, *, max_rows: int | None = None) -> str:
        """Render as a GitHub-flavoured Markdown table, optionally truncated."""

        header = (
            "| source | hkl | d (A) | r (mm) | s_g (1/A) | I_rel |\n"
            "| --- | --- | --- | --- | --- | --- |\n"
        )
        rows = self.rows if max_rows is None else self.rows[:max_rows]
        body = "".join(
            f"| {row.phase_name}"
            f"{'' if row.variant_index is None else f' V{row.variant_index}'} "
            f"| {row.hkl_label}{' (dd)' if row.is_double_diffraction else ''} "
            f"| {row.d_angstrom:.4f} | {row.detector_radius_mm:.3f} "
            f"| {row.excitation_error_inv_angstrom:+.5f} | {row.relative_intensity:.4f} |\n"
            for row in rows
        )
        omitted = (
            ""
            if max_rows is None or len(self.rows) <= max_rows
            else f"\n({len(self.rows) - max_rows} further row(s) omitted.)\n"
        )
        # Without the note a reader would take a forbidden reflection for a
        # genuine one, which is the whole hazard of showing these at all.
        footnote = (
            "\n(dd) kinematically forbidden; present only through double diffraction, at an "
            "indicative intensity.\n"
            if any(row.is_double_diffraction for row in rows)
            else ""
        )
        return header + body + footnote + omitted

    def to_json_dict(self) -> dict[str, Any]:
        """A JSON-ready payload carrying the table and the settings behind it."""

        return {
            "schema": COMPOSITE_REFLECTION_TABLE_SCHEMA,
            "relationship": self.relationship_name,
            "parent_phase": self.parent_phase_name,
            "child_phase": self.child_phase_name,
            "parent_zone_axis": self.parent_zone_axis_label,
            "camera_constant_mm_angstrom": self.camera_constant_mm_angstrom,
            "beam_energy_kev": self.beam_energy_kev,
            "wavelength_angstrom": self.wavelength_angstrom,
            "intensity_threshold": self.intensity_threshold,
            "centering_audit": [
                {"phase": name, "centering": centering, "declared": declared}
                for name, centering, declared in self.centering_audit
            ],
            "columns": list(REFLECTION_TABLE_COLUMNS),
            "rows": self.to_records(),
        }

    def describe(self) -> str:
        """Prose summary: what the table holds and what its numbers mean."""

        parent_rows = self.rows_for_source("parent")
        variant_count = len(self.variant_indices)
        assumed = [name for name, _, declared in self.centering_audit if not declared]
        lines = [
            f"Reflection table for composite pattern '{self.relationship_name}' along parent "
            f"zone axis {self.parent_zone_axis_label}: {len(self.rows)} reflection(s) — "
            f"{len(parent_rows)} from the parent '{self.parent_phase_name}' and "
            f"{len(self.rows) - len(parent_rows)} across {variant_count} variant(s) of "
            f"'{self.child_phase_name}'. Recorded at {self.beam_energy_kev:g} kV "
            f"(lambda = {self.wavelength_angstrom:.6f} A), camera constant "
            f"{self.camera_constant_mm_angstrom:g} mm*A. Detector radius is the camera "
            "constant times the *in-plane* part of g, so it is slightly smaller than "
            "camera constant x |g|: the difference is the out-of-plane component that the "
            "excitation error records, and d = 1/|g| uses the full vector.",
            "Intensities are kinematic and max-normalized within each sub-pattern "
            "separately, because kinematic theory defines no intensity ratio between two "
            "different phases: compare intensities within one source, never across sources.",
        ]
        forbidden = [row for row in self.rows if row.is_double_diffraction]
        if forbidden:
            lines.append(
                f"{len(forbidden)} row(s) are kinematically forbidden and listed only because "
                "double diffraction is enabled: their indices are the algebraic sum of two "
                "excited reflections, so a real pattern shows them even though the structure "
                "factor is zero. Their intensity is an observability estimate, not a kinematic "
                "intensity; the double_diffraction_origin column names the path."
            )
        if self.intensity_threshold > 0.0:
            lines.append(
                f"Reflections below a relative intensity of {self.intensity_threshold:g} "
                "were excluded from this table (they remain in the pattern)."
            )
        if assumed:
            lines.append(
                "WARNING: lattice centering was assumed primitive for "
                + ", ".join(assumed)
                + " because the phase declares no space group, so systematically absent "
                "reflections may be listed."
            )
        return " ".join(lines)


def composite_reflection_table(
    pattern: CompositeSAEDPattern,
    *,
    intensity_threshold: float = 0.0,
    provenance: ProvenanceRecord | None = None,
) -> ReflectionTable:
    """Tabulate every rendered reflection of a composite SAED pattern.

    Purpose: turns a simulated composite into the reflection list a paper or a
    lab notebook needs — one row per spot, carrying its source (parent or
    variant ``k``), phase, Miller indices and formatted label, d-spacing,
    ``|g|``, detector position and radius in millimetres, excitation error,
    structure-factor amplitude, and relative intensity.

    When to use: whenever the pattern must leave the process — for publication
    tables, for comparison against a measured pattern, or as the tabular
    companion to the rendered figure. `export_composite_saed` writes it to disk
    alongside the figure and a manifest.

    Inputs: the `CompositeSAEDPattern`; ``intensity_threshold``, a relative
    intensity below which reflections are omitted from the *table* (they remain
    in the pattern, and `describe()` states that the filter was applied).

    Output: a `ReflectionTable` — read its ``describe()``.

    Every value is read from the `SpotTable` objects the engine produced, so
    the table and the rendered figure cannot disagree.

    See also
    --------
    `pytex.diffraction.composite.simulate_composite_saed` : produces the pattern.
    `export_composite_saed` : writes table, figure and manifest together.
    """

    if not 0.0 <= intensity_threshold < 1.0:
        raise ValueError("intensity_threshold must lie in the interval [0, 1).")

    rows: list[ReflectionTableRow] = []

    def _append(spots: Any, *, source: str, variant_index: int | None) -> None:
        labels = spots.hkl_labels()
        zone_label = spots.zone_axis_label()
        for position in range(len(spots)):
            intensity = float(spots.intensity[position])
            if intensity < intensity_threshold:
                continue
            magnitude = float(np.linalg.norm(spots.g_crystal[position]))
            rows.append(
                ReflectionTableRow(
                    source=source,
                    variant_index=variant_index,
                    phase_name=spots.phase.name,
                    hkl=(
                        int(spots.hkl[position, 0]),
                        int(spots.hkl[position, 1]),
                        int(spots.hkl[position, 2]),
                    ),
                    hkl_label=labels[position],
                    zone_axis_label=zone_label,
                    d_angstrom=float(spots.d_spacing_angstrom[position]),
                    g_inv_angstrom=magnitude,
                    detector_mm=(
                        float(spots.detector_mm[position, 0]),
                        float(spots.detector_mm[position, 1]),
                    ),
                    excitation_error_inv_angstrom=float(
                        spots.excitation_error_inv_angstrom[position]
                    ),
                    structure_factor_amplitude=float(
                        spots.structure_factor_amplitude[position]
                    ),
                    relative_intensity=intensity,
                    double_diffraction_origin=spots.double_diffraction_origin_label(position),
                )
            )

    if pattern.parent_spots is not None:
        _append(pattern.parent_spots, source="parent", variant_index=None)
    for variant_pattern in pattern.variant_patterns:
        _append(
            variant_pattern.spots,
            source="variant",
            variant_index=variant_pattern.variant_index,
        )


    return ReflectionTable(
        relationship_name=pattern.relationship.name,
        parent_phase_name=pattern.relationship.parent_phase.name,
        child_phase_name=pattern.relationship.child_phase.name,
        parent_zone_axis_label=pattern.parent_zone_axis_label(),
        camera_constant_mm_angstrom=float(pattern.config.camera_constant_mm_angstrom),
        beam_energy_kev=float(pattern.config.beam_energy_kev),
        wavelength_angstrom=float(pattern.config.wavelength_angstrom),
        centering_audit=pattern.centering_audit(),
        rows=tuple(rows),
        intensity_threshold=float(intensity_threshold),
        provenance=provenance or pattern.provenance,
    )


def _parent_zone_axis_payload(pattern: CompositeSAEDPattern) -> list[float]:
    """The parent zone axis as JSON numbers, integer or not.

    A parent-anchored pattern has integer indices. A pattern anchored on a
    child variant's zone axis has the exact — generally irrational — parent
    direction that corresponds to it, so the manifest records real components
    and relies on ``parent_zone_axis_nearest`` for the readable label.
    """

    from pytex.core.lattice import ZoneAxis

    axis = pattern.parent_zone_axis
    if isinstance(axis, ZoneAxis):
        return [float(int(value)) for value in axis.indices]
    return [float(value) for value in axis.coordinates]


def composite_saed_manifest(
    pattern: CompositeSAEDPattern,
    *,
    table: ReflectionTable | None = None,
    coincidences: SpotCoincidenceReport | None = None,
    files: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Machine-readable record of how a composite pattern was produced.

    Purpose: the data contract for a pattern crossing a tool boundary. It
    records the relationship, both phases and their applied lattice centering,
    the parent zone axis and every variant's exact and nearest-rational child
    zone axis, the full simulation configuration for parent and children, the
    reflection counts, and the inventory of files written beside it.

    Inputs: the pattern; optionally the `ReflectionTable` and
    `SpotCoincidenceReport` computed from it, and a mapping of role to filename
    for the files written alongside.

    Output: a JSON-serializable dictionary carrying
    `COMPOSITE_SAED_MANIFEST_SCHEMA`.
    """

    def _config(config: Any) -> dict[str, Any]:
        return {
            "beam_energy_kev": float(config.beam_energy_kev),
            "wavelength_angstrom": float(config.wavelength_angstrom),
            "camera_constant_mm_angstrom": float(config.camera_constant_mm_angstrom),
            "max_index": int(config.max_index),
            "g_max_inv_angstrom": (
                None if config.g_max_inv_angstrom is None else float(config.g_max_inv_angstrom)
            ),
            "max_excitation_error_inv_angstrom": float(
                config.max_excitation_error_inv_angstrom
            ),
            "intensity_model": str(config.intensity_model),
            "relrod_sigma_inv_angstrom": (
                None
                if config.relrod_sigma_inv_angstrom is None
                else float(config.relrod_sigma_inv_angstrom)
            ),
            "apply_centering_absences": bool(config.apply_centering_absences),
            "min_relative_intensity": float(config.min_relative_intensity),
            "include_double_diffraction": bool(config.include_double_diffraction),
            "double_diffraction_coupling": float(config.double_diffraction_coupling),
        }


    manifest: dict[str, Any] = {
        "schema": COMPOSITE_SAED_MANIFEST_SCHEMA,
        "relationship": pattern.relationship.name,
        "parent_phase": pattern.relationship.parent_phase.name,
        "child_phase": pattern.relationship.child_phase.name,
        "parent_zone_axis": _parent_zone_axis_payload(pattern),
        "parent_zone_axis_label": pattern.parent_zone_axis_label(),
        "anchor_variant_index": pattern.anchor_variant_index,
        "parent_zone_axis_nearest": (
            None
            if pattern.nearest_parent_zone_axis is None
            else {
                "indices": [
                    int(value) for value in pattern.nearest_parent_zone_axis.indices
                ],
                "label": pattern.nearest_parent_zone_axis.label(),
                "deviation_deg": float(
                    pattern.nearest_parent_zone_axis.deviation_deg
                ),
            }
        ),
        "includes_parent": pattern.parent_spots is not None,
        "parent_reflection_count": (
            0 if pattern.parent_spots is None else len(pattern.parent_spots)
        ),
        "zone_basis_parent": [
            [float(value) for value in row] for row in pattern.zone_basis_parent
        ],
        "centering_audit": [
            {"phase": name, "centering": centering, "declared": declared}
            for name, centering, declared in pattern.centering_audit()
        ],
        "parent_config": _config(pattern.config),
        "child_config": _config(pattern.child_config),
        "variants": [
            {
                "variant_index": variant_pattern.variant_index,
                "child_zone_axis_exact": [
                    float(value) for value in variant_pattern.zone_axis_child.coordinates
                ],
                "child_zone_axis_nearest": [
                    int(value) for value in variant_pattern.nearest_zone_axis.indices
                ],
                "child_zone_axis_label": variant_pattern.nearest_zone_axis.label(),
                "child_zone_axis_deviation_deg": float(
                    variant_pattern.nearest_zone_axis.deviation_deg
                ),
                "reflection_count": len(variant_pattern.spots),
            }
            for variant_pattern in pattern.variant_patterns
        ],
        "total_reflection_count": pattern.spot_count(),
        "intensity_normalization": "per_sub_pattern_max",
        "theory_level": "kinematic",
        "files": dict(files or {}),
    }
    if table is not None:
        manifest["table_row_count"] = len(table)
        manifest["table_intensity_threshold"] = table.intensity_threshold
    if coincidences is not None:
        manifest["coincidence_count"] = len(coincidences.coincidences)
        manifest["coincidence_tolerance_mm"] = float(coincidences.tolerance_mm)
    return manifest


@dataclass(frozen=True, slots=True)
class CompositeSAEDExport:
    """Inventory of the files written for one composite pattern."""

    directory: Path
    reflection_table_path: Path
    manifest_path: Path
    figure_paths: tuple[Path, ...] = ()
    coincidence_table_path: Path | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "figure_paths", tuple(self.figure_paths))

    def paths(self) -> tuple[Path, ...]:
        """Every file written, manifest last."""

        written: list[Path] = [self.reflection_table_path, *self.figure_paths]
        if self.coincidence_table_path is not None:
            written.append(self.coincidence_table_path)
        written.append(self.manifest_path)
        return tuple(written)

    def describe(self) -> str:
        """Prose summary of what was written and where."""

        names = ", ".join(path.name for path in self.paths())
        return (
            f"Composite SAED export in '{self.directory}': {len(self.paths())} file(s) — "
            f"{names}. The manifest records the relationship, both phases and their applied "
            "lattice centering, every variant's exact and nearest-rational child zone axis, "
            "and the full simulation configuration, so the export is reproducible without "
            "the originating script."
        )


def export_composite_saed(
    pattern: CompositeSAEDPattern,
    directory: str | Path,
    *,
    stem: str = "composite_saed",
    figure_formats: tuple[str, ...] = ("svg",),
    include_coincidences: bool = True,
    coincidence_tolerance_mm: float = 0.05,
    intensity_threshold: float = 0.0,
    plot_config: Any | None = None,
) -> CompositeSAEDExport:
    """Write a composite pattern's table, figure(s) and manifest to a directory.

    Purpose: one call that produces everything needed to report a simulated
    composite pattern — the reflection table as CSV, the rendered figure in the
    requested vector or raster formats, the parent/child spot-coincidence table
    when requested, and a JSON manifest recording how all of it was produced.

    When to use: at the end of a simulation, when the results are to be shared,
    published, or compared against a measured pattern.

    Inputs: the pattern and an output ``directory`` (created if missing);
    ``stem``, the shared filename prefix; ``figure_formats``, matplotlib output
    formats (pass an empty tuple to skip rendering, which also avoids importing
    matplotlib); ``include_coincidences`` and ``coincidence_tolerance_mm`` for
    the parent/child overlap table; ``intensity_threshold`` forwarded to the
    reflection table; ``plot_config``, an optional
    `pytex.plotting.composite_saed.CompositeSAEDPlotConfig`.

    Output: a `CompositeSAEDExport` naming every file written — read its
    ``describe()``.

    Figures are closed after writing, so calling this in a loop or a test does
    not leak matplotlib figures.
    """

    output_directory = Path(directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    table = composite_reflection_table(pattern, intensity_threshold=intensity_threshold)
    table_path = table.to_csv(output_directory / f"{stem}_reflections.csv")

    coincidence_report: SpotCoincidenceReport | None = None
    coincidence_path: Path | None = None
    if include_coincidences and pattern.parent_spots is not None:
        coincidence_report = find_spot_coincidences(
            pattern, tolerance_mm=coincidence_tolerance_mm
        )
        coincidence_path = _write_coincidence_csv(
            coincidence_report, output_directory / f"{stem}_coincidences.csv"
        )

    figure_paths: list[Path] = []
    if figure_formats:
        from pytex.plotting.composite_saed import render_composite_saed

        figure = render_composite_saed(pattern, config=plot_config)
        try:
            for suffix in figure_formats:
                path = output_directory / f"{stem}.{suffix.lstrip('.')}"
                figure.savefig(path, bbox_inches="tight")
                figure_paths.append(path)
        finally:
            import matplotlib.pyplot as plt

            plt.close(figure)

    manifest = composite_saed_manifest(
        pattern,
        table=table,
        coincidences=coincidence_report,
        files={
            "reflection_table": table_path.name,
            **(
                {"coincidence_table": coincidence_path.name}
                if coincidence_path is not None
                else {}
            ),
            **{
                f"figure_{path.suffix.lstrip('.')}": path.name for path in figure_paths
            },
        },
    )
    manifest_path = output_directory / f"{stem}_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    return CompositeSAEDExport(
        directory=output_directory,
        reflection_table_path=table_path,
        manifest_path=manifest_path,
        figure_paths=tuple(figure_paths),
        coincidence_table_path=coincidence_path,
    )


def _write_coincidence_csv(report: SpotCoincidenceReport, path: Path) -> Path:
    """Write a parent/child spot-coincidence report as CSV."""

    fieldnames = (
        "variant",
        "parent_hkl",
        "child_hkl",
        "separation_mm",
        "parent_detector_x_mm",
        "parent_detector_y_mm",
        "child_detector_x_mm",
        "child_detector_y_mm",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for coincidence in report.coincidences:
            writer.writerow(
                {
                    "variant": coincidence.variant_index,
                    "parent_hkl": " ".join(
                        str(int(value)) for value in coincidence.parent_hkl
                    ),
                    "child_hkl": " ".join(
                        str(int(value)) for value in coincidence.child_hkl
                    ),
                    "separation_mm": float(coincidence.separation_mm),
                    "parent_detector_x_mm": float(coincidence.parent_detector_mm[0]),
                    "parent_detector_y_mm": float(coincidence.parent_detector_mm[1]),
                    "child_detector_x_mm": float(coincidence.child_detector_mm[0]),
                    "child_detector_y_mm": float(coincidence.child_detector_mm[1]),
                }
            )
    return path


__all__ = [
    "COMPOSITE_REFLECTION_TABLE_SCHEMA",
    "COMPOSITE_SAED_MANIFEST_SCHEMA",
    "REFLECTION_TABLE_COLUMNS",
    "CompositeSAEDExport",
    "ReflectionTable",
    "ReflectionTableRow",
    "composite_reflection_table",
    "composite_saed_manifest",
    "export_composite_saed",
]
