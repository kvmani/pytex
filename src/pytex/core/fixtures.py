from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def phase_fixture_catalog_path() -> Path:
    return _repo_root() / "fixtures/phases/catalog.json"


@lru_cache(maxsize=1)
def _phase_fixture_catalog_payload() -> dict[str, Any]:
    return json.loads(phase_fixture_catalog_path().read_text(encoding="utf-8"))


@dataclass(frozen=True, slots=True)
class PhaseFixtureRecord:
    fixture_id: str
    artifact_path: Path
    metadata_path: Path
    artifact_sha256: str
    metadata_sha256: str
    metadata: dict[str, Any]

    @property
    def display_name(self) -> str:
        return str(self.metadata["display_name"])

    @property
    def phase_name(self) -> str:
        return str(self.metadata["phase_name"])

    @property
    def chemical_formula(self) -> str:
        return str(self.metadata["chemical_formula"])

    def read_cif_text(self) -> str:
        return self.artifact_path.read_text(encoding="utf-8")

    def load_phase(
        self,
        *,
        crystal_frame: ReferenceFrame,
        primitive: bool = False,
        phase_name: str | None = None,
    ) -> Phase:
        return Phase.from_cif(
            self.artifact_path,
            crystal_frame=crystal_frame,
            primitive=primitive,
            phase_name=phase_name or self.phase_name,
        )


def list_phase_fixtures() -> tuple[PhaseFixtureRecord, ...]:
    repo_root = _repo_root()
    records: list[PhaseFixtureRecord] = []
    for entry in _phase_fixture_catalog_payload()["fixtures"]:
        metadata_path = repo_root / entry["metadata_path"]
        records.append(
            PhaseFixtureRecord(
                fixture_id=str(entry["fixture_id"]),
                artifact_path=repo_root / str(entry["artifact_path"]),
                metadata_path=metadata_path,
                artifact_sha256=str(entry["artifact_sha256"]),
                metadata_sha256=str(entry["metadata_sha256"]),
                metadata=json.loads(metadata_path.read_text(encoding="utf-8")),
            )
        )
    return tuple(records)


def get_phase_fixture(fixture_id: str) -> PhaseFixtureRecord:
    for record in list_phase_fixtures():
        if record.fixture_id == fixture_id:
            return record
    raise KeyError(f"Unknown phase fixture id: {fixture_id}")
