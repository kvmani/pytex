from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, cast

from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def phase_fixture_catalog_path() -> Path:
    """Path to the repository's phase-fixture catalog JSON.

    The corpus is a **repository asset**, not package data: it is checksum-pinned
    in ``fixtures/phases/`` and validated by the repository integrity script, so
    it is available in a source checkout but is not shipped inside an installed
    wheel. See :func:`list_phase_fixtures` for what to do in an installed
    environment.
    """

    return _repo_root() / "fixtures/phases/catalog.json"


def phase_fixtures_available() -> bool:
    """Whether the checksum-pinned phase-fixture corpus is reachable.

    Purpose
    -------
    The corpus ships with the repository, not with the wheel, so code that can
    run in either context should test for it rather than catch a file error.

    Returns
    -------
    bool
        ``True`` in a source checkout with the corpus present; ``False`` in an
        installed environment.
    """

    return phase_fixture_catalog_path().is_file()


@lru_cache(maxsize=1)
def _phase_fixture_catalog_payload() -> dict[str, Any]:
    catalog_path = phase_fixture_catalog_path()
    if not catalog_path.is_file():
        raise FileNotFoundError(
            "The PyTex phase-fixture corpus is not available. It is a checksum-pinned "
            "repository asset under 'fixtures/phases/', not package data, so it is "
            "present in a source checkout but not inside an installed wheel. Either "
            "work from a clone of the repository, or construct the Phase directly from "
            "its lattice, symmetry and crystal frame, or load it from your own CIF with "
            "Phase.from_cif(...). Use phase_fixtures_available() to branch on this. "
            f"Looked for: {catalog_path}"
        )
    payload = json.loads(catalog_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("fixtures/phases/catalog.json must decode to a JSON object.")
    return cast(dict[str, Any], payload)


@dataclass(frozen=True, slots=True)
class PhaseFixtureRecord:
    """A checksummed reference crystal structure from the repository catalog.

    Purpose
    -------
    Documentation, tests, and worked examples must rest on pinned structural
    data rather than on hand-typed lattice parameters. Each record pairs a
    CIF artifact with its metadata and the SHA-256 of both, so a silently
    changed fixture is detectable rather than quietly shifting every number
    that depends on it.

    Attributes
    ----------
    fixture_id : str
        Catalog identifier.
    artifact_path : Path
        The CIF file.
    metadata_path : Path
    artifact_sha256, metadata_sha256 : str
        Pinned checksums.
    metadata : dict
        Display name, phase name, chemical formula, and source information.
    """

    fixture_id: str
    artifact_path: Path
    metadata_path: Path
    artifact_sha256: str
    metadata_sha256: str
    metadata: dict[str, Any]

    @property
    def display_name(self) -> str:
        """Human-readable name of the fixture, for figures and prose.
        """

        return str(self.metadata["display_name"])

    @property
    def phase_name(self) -> str:
        """Canonical phase name recorded in the fixture metadata.
        """

        return str(self.metadata["phase_name"])

    @property
    def chemical_formula(self) -> str:
        """Chemical formula recorded in the fixture metadata.
        """

        return str(self.metadata["chemical_formula"])

    def read_cif_text(self) -> str:
        """The raw CIF text of the fixture artifact.
        """

        return self.artifact_path.read_text(encoding="utf-8")

    def load_phase(
        self,
        *,
        crystal_frame: ReferenceFrame,
        primitive: bool = False,
        phase_name: str | None = None,
    ) -> Phase:
        """Load this fixture as a fully specified :class:`~pytex.core.lattice.Phase`.

        Requires the optional pymatgen dependency, since the CIF must be parsed
        and its symmetry determined.

        Parameters
        ----------
        crystal_frame : ReferenceFrame
            The crystal-domain frame to attach.
        primitive : bool
            Reduce to the primitive cell. Off by default, so the conventional
            cell of the CIF is kept and reflection conditions remain those of the
            conventional setting.
        phase_name : str, optional
            Overrides the name derived from the structure.
        """

        return Phase.from_cif(
            self.artifact_path,
            crystal_frame=crystal_frame,
            primitive=primitive,
            phase_name=phase_name or self.phase_name,
        )


def list_phase_fixtures() -> tuple[PhaseFixtureRecord, ...]:
    """Every phase fixture in the repository catalog.

    Purpose
    -------
    The curated set of real crystal structures — with checksummed CIF
    artifacts and metadata — that documentation, tests, and worked examples
    build on, so scientific claims rest on pinned structural data rather than
    on hand-typed lattice parameters.

    Returns
    -------
    tuple of PhaseFixtureRecord
        In catalog order, each carrying artifact and metadata paths together
        with their SHA-256 checksums.

    Raises
    ------
    FileNotFoundError
        When the corpus is not reachable. It is a repository asset rather than
        package data, so it is present in a source checkout but not inside an
        installed wheel; test with :func:`phase_fixtures_available` and build
        the phase directly, or from your own CIF, in that case.
    """

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
    """One phase fixture by its catalog id.

    Raises ``KeyError`` for an unknown id, and ``FileNotFoundError`` when the
    corpus is unavailable; see :func:`list_phase_fixtures`.
    """

    for record in list_phase_fixtures():
        if record.fixture_id == fixture_id:
            return record
    raise KeyError(f"Unknown phase fixture id: {fixture_id}")
