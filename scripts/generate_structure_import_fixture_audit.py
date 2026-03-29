from __future__ import annotations

import json
from pathlib import Path

from pytex import FrameDomain, Handedness, ReferenceFrame, list_phase_fixtures


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _make_crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def _summarize_fixture(record_fixture_id: str, *, primitive: bool) -> dict[str, object]:
    crystal = _make_crystal_frame()
    record = next(item for item in list_phase_fixtures() if item.fixture_id == record_fixture_id)
    phase = record.load_phase(
        crystal_frame=crystal,
        primitive=primitive,
        phase_name=record.phase_name,
    )
    return {
        "space_group_symbol": phase.space_group_symbol,
        "space_group_number": phase.space_group_number,
        "point_group": phase.symmetry.point_group,
        "chemical_formula": phase.chemical_formula,
        "site_count": 0 if phase.unit_cell is None else len(phase.unit_cell.sites),
        "lattice_parameters_angstrom": {
            "a": phase.lattice.a,
            "b": phase.lattice.b,
            "c": phase.lattice.c,
        },
        "lattice_angles_deg": {
            "alpha": phase.lattice.alpha_deg,
            "beta": phase.lattice.beta_deg,
            "gamma": phase.lattice.gamma_deg,
        },
    }


def build_summary() -> dict[str, object]:
    repo_root = _repo_root()
    fixtures: list[dict[str, object]] = []
    for record in list_phase_fixtures():
        fixtures.append(
            {
                "fixture_id": record.fixture_id,
                "display_name": record.display_name,
                "phase_name": record.phase_name,
                "chemical_formula": record.chemical_formula,
                "artifact_path": record.artifact_path.relative_to(repo_root).as_posix(),
                "metadata_path": record.metadata_path.relative_to(repo_root).as_posix(),
                "intended_uses": list(record.metadata["intended_uses"]),
                "conventional": _summarize_fixture(record.fixture_id, primitive=False),
                "primitive": _summarize_fixture(record.fixture_id, primitive=True),
            }
        )
    return {
        "schema_id": "pytex.structure_import_fixture_audit",
        "schema_version": "1.0.0",
        "fixture_catalog_id": "pytex.phase_fixture_catalog",
        "fixtures": fixtures,
    }


def main() -> int:
    repo_root = _repo_root()
    output_path = repo_root / "benchmarks/structure_import/phase_fixture_audit_summary.json"
    output_path.write_text(json.dumps(build_summary(), indent=2) + "\n", encoding="utf-8")
    print(f"Wrote structure-import fixture audit to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
