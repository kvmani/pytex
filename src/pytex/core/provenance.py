from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
    """Where a value came from, carried alongside the value itself.

    Purpose
    -------
    Traceability is a stated engineering priority of this library, ranking
    above maintainability and speed. Provenance records travel through
    derived objects, so a number in a report can be traced back to the file,
    tool, and conventions that produced it.

    Attributes
    ----------
    source_system : str
        The producing tool or system.
    source_file : str, optional
    notes : tuple of str
        Free-text notes, including original-convention remarks from an
        imported dataset.
    """

    source_system: str
    source_identifier: str | None = None
    source_path: str | None = None
    source_version: str | None = None
    imported_at: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
        object.__setattr__(self, "notes", tuple(self.notes))

    @classmethod
    def minimal(cls, source_system: str, *, note: str | None = None) -> ProvenanceRecord:
        """A provenance record carrying just a source system and an optional note.

        For cases where full provenance is not available but the origin of a
        value should still travel with it.
        """

        notes = (note,) if note else ()
        return cls(source_system=source_system, notes=notes)
