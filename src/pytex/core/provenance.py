from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
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
        notes = (note,) if note else ()
        return cls(source_system=source_system, notes=notes)
