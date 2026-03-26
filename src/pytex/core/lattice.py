from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array, normalize_vector
from pytex.core.conventions import PYTEX_CANONICAL_CONVENTIONS, BasisKind, FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec


def _require_pymatgen() -> tuple[Any, Any]:
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:
        raise ImportError(
            "CIF-backed phase creation requires the optional 'pymatgen' dependency. "
            "Install PyTex with the 'adapters' extra."
        ) from exc
    return Structure, SpacegroupAnalyzer


def _sites_from_pymatgen_structure(structure: Any) -> tuple[AtomicSite, ...]:
    sites: list[AtomicSite] = []
    for site_index, site in enumerate(structure.sites, start=1):
        fractional_coordinates = np.asarray(site.frac_coords, dtype=np.float64)
        label_base = str(getattr(site, "label", "") or f"site_{site_index}")
        species_items = list(site.species.items())
        if not species_items:
            raise ValueError("Encountered a structure site with no species information.")
        b_iso_raw = dict(getattr(site, "properties", {})).get("B_iso")
        b_iso = None if b_iso_raw is None else float(b_iso_raw)
        for species_index, (species, occupancy) in enumerate(species_items, start=1):
            label = label_base if len(species_items) == 1 else f"{label_base}_{species_index}"
            sites.append(
                AtomicSite(
                    label=label,
                    species=str(species),
                    fractional_coordinates=fractional_coordinates,
                    occupancy=float(occupancy),
                    b_iso=b_iso,
                )
            )
    return tuple(sites)


@dataclass(frozen=True, slots=True)
class Basis:
    frame: ReferenceFrame
    kind: BasisKind
    matrix: np.ndarray
    unit: str = "angstrom"

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", as_float_array(self.matrix, shape=(3, 3)))
        if self.kind is BasisKind.DIRECT and self.frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("Direct bases must use a crystal-domain reference frame.")
        if self.kind is BasisKind.RECIPROCAL and self.frame.domain is not FrameDomain.RECIPROCAL:
            raise ValueError("Reciprocal bases must use a reciprocal-domain reference frame.")

    def vector(self, index: int) -> np.ndarray:
        return as_float_array(self.matrix[:, index], shape=(3,))


@dataclass(frozen=True, slots=True)
class Lattice:
    a: float
    b: float
    c: float
    alpha_deg: float
    beta_deg: float
    gamma_deg: float
    crystal_frame: ReferenceFrame
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.crystal_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("Lattice.crystal_frame must belong to the crystal domain.")
        for name, value in (("a", self.a), ("b", self.b), ("c", self.c)):
            if value <= 0.0:
                raise ValueError(f"Lattice parameter '{name}' must be strictly positive.")
        for name, value in (
            ("alpha_deg", self.alpha_deg),
            ("beta_deg", self.beta_deg),
            ("gamma_deg", self.gamma_deg),
        ):
            if not 0.0 < value < 180.0:
                raise ValueError(f"Lattice angle '{name}' must lie strictly between 0 and 180.")

    @classmethod
    def from_pymatgen_lattice(
        cls,
        pymatgen_lattice: Any,
        *,
        crystal_frame: ReferenceFrame,
        provenance: ProvenanceRecord | None = None,
    ) -> Lattice:
        lengths = tuple(float(value) for value in pymatgen_lattice.abc)
        angles = tuple(float(value) for value in pymatgen_lattice.angles)
        return cls(
            a=lengths[0],
            b=lengths[1],
            c=lengths[2],
            alpha_deg=angles[0],
            beta_deg=angles[1],
            gamma_deg=angles[2],
            crystal_frame=crystal_frame,
            provenance=provenance,
        )

    def direct_basis(self) -> Basis:
        alpha = np.deg2rad(self.alpha_deg)
        beta = np.deg2rad(self.beta_deg)
        gamma = np.deg2rad(self.gamma_deg)
        sin_gamma = np.sin(gamma)
        if np.isclose(sin_gamma, 0.0):
            raise ValueError("gamma must not yield a degenerate lattice basis.")
        a_vec = np.array([self.a, 0.0, 0.0])
        b_vec = np.array([self.b * np.cos(gamma), self.b * sin_gamma, 0.0])
        c_x = self.c * np.cos(beta)
        c_y = self.c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_gamma
        c_z_sq = self.c**2 - c_x**2 - c_y**2
        if c_z_sq <= 0:
            raise ValueError("Lattice parameters produce a non-physical unit cell.")
        c_vec = np.array([c_x, c_y, np.sqrt(c_z_sq)])
        matrix = np.column_stack([a_vec, b_vec, c_vec])
        return Basis(frame=self.crystal_frame, kind=BasisKind.DIRECT, matrix=matrix)

    def reciprocal_basis(self) -> Basis:
        direct = self.direct_basis().matrix
        reciprocal = np.linalg.inv(direct).T
        reciprocal_frame = ReferenceFrame(
            name=f"{self.crystal_frame.name}_reciprocal",
            domain=FrameDomain.RECIPROCAL,
            axes=self.crystal_frame.axes,
            handedness=self.crystal_frame.handedness,
            convention=PYTEX_CANONICAL_CONVENTIONS,
            description=f"Reciprocal basis associated with {self.crystal_frame.name}.",
            provenance=self.provenance,
        )
        return Basis(
            frame=reciprocal_frame,
            kind=BasisKind.RECIPROCAL,
            matrix=reciprocal,
            unit="1/angstrom",
        )


@dataclass(frozen=True, slots=True)
class AtomicSite:
    label: str
    species: str
    fractional_coordinates: np.ndarray
    occupancy: float = 1.0
    b_iso: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "fractional_coordinates", as_float_array(self.fractional_coordinates, shape=(3,))
        )
        if not 0.0 < self.occupancy <= 1.0:
            raise ValueError("AtomicSite.occupancy must lie in the interval (0, 1].")
        if self.b_iso is not None and self.b_iso < 0.0:
            raise ValueError("AtomicSite.b_iso must be non-negative when provided.")


@dataclass(frozen=True, slots=True)
class UnitCell:
    lattice: Lattice
    sites: tuple[AtomicSite, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "sites", tuple(self.sites))

    @classmethod
    def from_pymatgen_structure(
        cls,
        structure: Any,
        *,
        crystal_frame: ReferenceFrame,
        lattice: Lattice | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> UnitCell:
        unit_cell_lattice = lattice or Lattice.from_pymatgen_lattice(
            structure.lattice,
            crystal_frame=crystal_frame,
            provenance=provenance,
        )
        return cls(
            lattice=unit_cell_lattice,
            sites=_sites_from_pymatgen_structure(structure),
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class Phase:
    name: str
    lattice: Lattice
    symmetry: SymmetrySpec
    crystal_frame: ReferenceFrame
    unit_cell: UnitCell | None = None
    space_group_symbol: str | None = None
    space_group_number: int | None = None
    chemical_formula: str | None = None
    aliases: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "aliases", tuple(self.aliases))
        if self.crystal_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("Phase.crystal_frame must belong to the crystal domain.")
        if self.lattice.crystal_frame != self.crystal_frame:
            raise ValueError("Phase.lattice.crystal_frame must match Phase.crystal_frame.")
        if (
            self.symmetry.reference_frame is not None
            and self.symmetry.reference_frame != self.crystal_frame
        ):
            raise ValueError("Phase.symmetry.reference_frame must match Phase.crystal_frame.")
        if self.unit_cell is not None and self.unit_cell.lattice != self.lattice:
            raise ValueError("Phase.unit_cell.lattice must match Phase.lattice.")
        if self.space_group_symbol is not None and not self.space_group_symbol:
            raise ValueError("Phase.space_group_symbol must be non-empty when provided.")
        if self.space_group_number is not None and not 1 <= self.space_group_number <= 230:
            raise ValueError("Phase.space_group_number must lie in the interval [1, 230].")
        if self.chemical_formula is not None and not self.chemical_formula:
            raise ValueError("Phase.chemical_formula must be non-empty when provided.")

    @classmethod
    def from_pymatgen_structure(
        cls,
        structure: Any,
        *,
        crystal_frame: ReferenceFrame,
        phase_name: str | None = None,
        aliases: tuple[str, ...] = (),
        symprec: float = 1e-3,
        angle_tolerance: float = 5.0,
        provenance: ProvenanceRecord | None = None,
    ) -> Phase:
        _, spacegroup_analyzer_cls = _require_pymatgen()
        analyzer = spacegroup_analyzer_cls(
            structure,
            symprec=float(symprec),
            angle_tolerance=float(angle_tolerance),
        )
        point_group = str(analyzer.get_point_group_symbol())
        symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal_frame)
        lattice = Lattice.from_pymatgen_lattice(
            structure.lattice,
            crystal_frame=crystal_frame,
            provenance=provenance,
        )
        unit_cell = UnitCell.from_pymatgen_structure(
            structure,
            crystal_frame=crystal_frame,
            lattice=lattice,
            provenance=provenance,
        )
        formula = str(structure.composition.reduced_formula)
        return cls(
            name=phase_name or formula,
            lattice=lattice,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
            unit_cell=unit_cell,
            space_group_symbol=str(analyzer.get_space_group_symbol()),
            space_group_number=int(analyzer.get_space_group_number()),
            chemical_formula=formula,
            aliases=aliases,
            provenance=provenance,
        )

    @classmethod
    def from_cif(
        cls,
        path: str | Path,
        *,
        crystal_frame: ReferenceFrame,
        phase_name: str | None = None,
        aliases: tuple[str, ...] = (),
        primitive: bool = False,
        symprec: float = 1e-3,
        angle_tolerance: float = 5.0,
        provenance: ProvenanceRecord | None = None,
    ) -> Phase:
        structure_cls, _ = _require_pymatgen()
        cif_path = Path(path)
        structure = structure_cls.from_file(str(cif_path))
        if primitive:
            structure = structure.get_primitive_structure()
        record = provenance or ProvenanceRecord(
            source_system="cif",
            source_path=str(cif_path),
            metadata={"reader": "pymatgen.Structure.from_file"},
        )
        return cls.from_pymatgen_structure(
            structure,
            crystal_frame=crystal_frame,
            phase_name=phase_name,
            aliases=aliases,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            provenance=record,
        )

    @classmethod
    def from_cif_string(
        cls,
        cif_text: str,
        *,
        crystal_frame: ReferenceFrame,
        phase_name: str | None = None,
        aliases: tuple[str, ...] = (),
        primitive: bool = False,
        symprec: float = 1e-3,
        angle_tolerance: float = 5.0,
        provenance: ProvenanceRecord | None = None,
    ) -> Phase:
        structure_cls, _ = _require_pymatgen()
        structure = structure_cls.from_str(cif_text, fmt="cif")
        if primitive:
            structure = structure.get_primitive_structure()
        record = provenance or ProvenanceRecord(
            source_system="cif",
            metadata={"reader": "pymatgen.Structure.from_str", "format": "cif"},
        )
        return cls.from_pymatgen_structure(
            structure,
            crystal_frame=crystal_frame,
            phase_name=phase_name,
            aliases=aliases,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            provenance=record,
        )


@dataclass(frozen=True, slots=True)
class MillerIndex:
    indices: np.ndarray
    phase: Phase
    basis_kind: BasisKind = BasisKind.RECIPROCAL

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("Miller indices must not be the zero triplet.")

    def as_array(self) -> np.ndarray:
        return self.indices


@dataclass(frozen=True, slots=True)
class CrystalDirection:
    coordinates: np.ndarray
    phase: Phase
    basis_kind: BasisKind = BasisKind.DIRECT

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", as_float_array(self.coordinates, shape=(3,)))
        if np.allclose(self.coordinates, 0.0):
            raise ValueError("CrystalDirection coordinates must not define the zero direction.")

    @property
    def unit_vector(self) -> np.ndarray:
        if self.basis_kind is BasisKind.DIRECT:
            basis = self.phase.lattice.direct_basis()
        else:
            basis = self.phase.lattice.reciprocal_basis()
        cartesian = basis.matrix @ self.coordinates
        return normalize_vector(cartesian)


@dataclass(frozen=True, slots=True)
class ZoneAxis:
    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("ZoneAxis indices must not be the zero triplet.")

    @property
    def direction(self) -> CrystalDirection:
        return CrystalDirection(self.indices.astype(np.float64), phase=self.phase)

    @property
    def unit_vector(self) -> np.ndarray:
        return self.direction.unit_vector

    def zone_law_value(self, miller: MillerIndex) -> int:
        if miller.phase != self.phase:
            raise ValueError("MillerIndex.phase must match ZoneAxis.phase.")
        return int(np.dot(self.indices, miller.indices))

    def contains_miller_index(self, miller: MillerIndex) -> bool:
        return self.zone_law_value(miller) == 0


@dataclass(frozen=True, slots=True)
class ReciprocalLatticeVector:
    coordinates: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", as_float_array(self.coordinates, shape=(3,)))
        if np.allclose(self.coordinates, 0.0):
            raise ValueError("ReciprocalLatticeVector coordinates must not be the zero vector.")

    @classmethod
    def from_miller_index(cls, miller: MillerIndex) -> ReciprocalLatticeVector:
        return cls(coordinates=miller.indices.astype(np.float64), phase=miller.phase)

    @property
    def cartesian_vector(self) -> np.ndarray:
        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        return as_float_array(reciprocal @ self.coordinates, shape=(3,))

    @property
    def magnitude_inv_angstrom(self) -> float:
        return float(np.linalg.norm(self.cartesian_vector))

    @property
    def unit_vector(self) -> np.ndarray:
        return normalize_vector(self.cartesian_vector)


@dataclass(frozen=True, slots=True)
class CrystalPlane:
    miller: MillerIndex
    phase: Phase

    def __post_init__(self) -> None:
        if self.miller.phase != self.phase:
            raise ValueError("CrystalPlane.miller.phase must match CrystalPlane.phase.")

    @property
    def normal(self) -> np.ndarray:
        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        normal = reciprocal @ self.miller.indices.astype(np.float64)
        return normalize_vector(normal)

    @property
    def d_spacing_angstrom(self) -> float:
        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        reciprocal_vector = reciprocal @ self.miller.indices.astype(np.float64)
        magnitude = float(np.linalg.norm(reciprocal_vector))
        if np.isclose(magnitude, 0.0):
            raise ValueError("CrystalPlane reciprocal vector magnitude must be non-zero.")
        return 1.0 / magnitude

    @property
    def reciprocal_lattice_vector(self) -> ReciprocalLatticeVector:
        return ReciprocalLatticeVector.from_miller_index(self.miller)
