from __future__ import annotations

import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array, normalize_vector
from pytex.core.conventions import BasisKind, FrameDomain
from pytex.core.frame_catalog import reciprocal_frame_for
from pytex.core.frames import ReferenceFrame
from pytex.core.hexagonal import direction_uvtw_to_uvw, plane_hkil_to_hkl
from pytex.core.point_groups import normalize_point_group_symbol
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec


@contextmanager
def _spglib_dict_shim_silenced() -> Iterator[None]:
    """Silence spglib's force-enabled dict-interface DeprecationWarning.

    Older pymatgen accessors read spglib's symmetry dataset through its dict
    shim, which re-enables its own DeprecationWarning inside a private filter
    context, so neither caller filters nor pytest configuration can silence
    it. This adapter-boundary context records emissions and re-emits every
    warning except that one shim message.
    """

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        yield
    for item in captured:
        if not str(item.message).startswith("dict interface is deprecated"):
            warnings.warn_explicit(item.message, item.category, item.filename, item.lineno)


def _require_pymatgen() -> tuple[Any, Any]:
    try:
        from pymatgen.core.structure import Structure
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
class SpaceGroupSpec:
    """A space-group identification: symbol, number, and reference frame.

    Purpose
    -------
    Records the space group so that centring-based reflection conditions can
    be applied. Its leading letter drives
    :meth:`~pytex.diffraction.ReflectionCondition.from_phase`; without it,
    systematic absences cannot be determined and forbidden reflections will
    be listed as present.

    Attributes
    ----------
    symbol : str
        Hermann-Mauguin symbol.
    number : int
        International Tables number, in ``[1, 230]``.
    reference_frame : ReferenceFrame
    provenance : ProvenanceRecord, optional
    """

    symbol: str
    number: int
    reference_frame: ReferenceFrame
    setting: str | None = None
    crystal_system: str | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.reference_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("SpaceGroupSpec.reference_frame must belong to the crystal domain.")
        normalized_symbol = self.symbol.strip()
        if not normalized_symbol:
            raise ValueError("SpaceGroupSpec.symbol must be non-empty.")
        if not 1 <= self.number <= 230:
            raise ValueError("SpaceGroupSpec.number must lie in the interval [1, 230].")
        if self.setting is not None and not self.setting.strip():
            raise ValueError("SpaceGroupSpec.setting must be non-empty when provided.")
        if self.crystal_system is not None and not self.crystal_system.strip():
            raise ValueError("SpaceGroupSpec.crystal_system must be non-empty when provided.")
        object.__setattr__(self, "symbol", normalized_symbol)
        if self.setting is not None:
            object.__setattr__(self, "setting", self.setting.strip())
        if self.crystal_system is not None:
            object.__setattr__(self, "crystal_system", self.crystal_system.strip().lower())

    @classmethod
    def from_pymatgen_analyzer(
        cls,
        analyzer: Any,
        *,
        reference_frame: ReferenceFrame,
        provenance: ProvenanceRecord | None = None,
    ) -> SpaceGroupSpec:
        """Build a space-group specification from a pymatgen symmetry analyzer.

        Optional-dependency bridge: it reads the symbol and number that
        pymatgen's spglib-backed analysis determined, and attaches the PyTex
        reference frame, so the space group enters the data model with explicit
        frame meaning rather than as a bare string.
        """

        with _spglib_dict_shim_silenced():
            return cls(
                symbol=str(analyzer.get_space_group_symbol()),
                number=int(analyzer.get_space_group_number()),
                reference_frame=reference_frame,
                crystal_system=str(analyzer.get_crystal_system()),
                provenance=provenance,
            )


@dataclass(frozen=True, slots=True)
class Basis:
    """A set of three basis vectors with its frame, kind, and unit declared.

    Purpose
    -------
    Distinguishes a direct basis ``(a, b, c)`` from a reciprocal basis
    ``(a*, b*, c*)`` at the type level. Outside the cubic system the two are
    not interchangeable, and a bare matrix cannot say which it is — so the
    kind and unit travel with the numbers.

    Attributes
    ----------
    frame : ReferenceFrame
    kind : BasisKind
        ``DIRECT`` or ``RECIPROCAL``.
    matrix : np.ndarray
        ``(3, 3)`` with the basis vectors as columns.
    unit : str
        ``"angstrom"`` for direct, ``"1/angstrom"`` for reciprocal.
    """

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
        """One basis vector, as a Cartesian 3-vector in this basis's frame.

        Column ``index`` of the basis matrix: ``0`` is ``a``, ``1`` is ``b``,
        ``2`` is ``c`` for a direct basis, and ``a*``, ``b*``, ``c*`` for a
        reciprocal one. Units follow the basis: angstroms for direct, inverse
        angstroms for reciprocal.
        """

        return as_float_array(self.matrix[:, index], shape=(3,))


@dataclass(frozen=True, slots=True)
class Lattice:
    """A crystal lattice: six cell parameters plus the frame they are stated in.

    Purpose
    -------
    The geometric foundation of every crystallographic calculation in the
    library. It fixes the crystal-to-Cartesian convention once, in
    :meth:`direct_basis`, so that directions, plane normals, d-spacings,
    structure factors, and orientation matrices all derive from one
    definition rather than from per-site formulae.

    Attributes
    ----------
    a, b, c : float
        Cell edge lengths in angstroms.
    alpha_deg, beta_deg, gamma_deg : float
        Cell angles in degrees.
    crystal_frame : ReferenceFrame
        The crystal-domain frame the basis is expressed in.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Reciprocal quantities follow the crystallographic convention *without*
    ``2*pi``, so ``|g_hkl| = 1/d_hkl`` directly.
    """

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
        """Build a lattice from a pymatgen ``Lattice``, attaching a crystal frame.

        Optional-dependency bridge. Only cell lengths and angles cross the
        boundary; pymatgen's own basis-vector convention is not adopted, because
        PyTex fixes its own crystal-to-Cartesian convention in
        :meth:`direct_basis`.
        """

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
        """The direct-space basis ``(a, b, c)`` as Cartesian columns, in angstroms.

        Purpose
        -------
        The single definition of the crystal-to-Cartesian convention in PyTex.
        Everything that turns indices into geometry — directions, plane normals,
        structure factors, orientation matrices — goes through this matrix, so
        the convention is fixed in exactly one place.

        Convention
        ----------
        ``a`` along ``x``; ``b`` in the ``x-y`` plane, making angle ``gamma``
        with ``a``; ``c`` completing the cell with the specified ``alpha`` and
        ``beta``. This is the standard crystallographic setting used by the
        IUCr-facing literature.

        Returns
        -------
        Basis
            Carrying the crystal frame, the basis kind, and the unit, so the
            matrix cannot be reused in the wrong frame by accident.
        """

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
        """The reciprocal basis ``(a*, b*, c*)`` as Cartesian columns.

        Computed as the inverse transpose of the direct basis, which realizes
        the duality ``a_i . a*_j = delta_ij``. Note the normalization: PyTex
        uses the crystallographic convention *without* a factor of ``2*pi``, so
        ``|g_hkl| = 1/d_hkl`` directly. Units are inverse angstroms, and the
        returned basis carries the reciprocal frame.
        """

        direct = self.direct_basis().matrix
        reciprocal = np.linalg.inv(direct).T
        reciprocal_frame = reciprocal_frame_for(self.crystal_frame, provenance=self.provenance)
        return Basis(
            frame=reciprocal_frame,
            kind=BasisKind.RECIPROCAL,
            matrix=reciprocal,
            unit="1/angstrom",
        )

    def metric_tensor(self) -> np.ndarray:
        """The direct metric tensor of this lattice; see :func:`metric_tensor`.
        """

        return metric_tensor(self)

    def reciprocal_metric_tensor(self) -> np.ndarray:
        """The reciprocal metric tensor of this lattice; see
        :func:`reciprocal_metric_tensor`.
        """

        return reciprocal_metric_tensor(self)

    def direct_to_reciprocal_components(self, components: Any) -> np.ndarray:
        """Convert direct-basis components to reciprocal-basis components in this
        lattice; see :func:`direct_to_reciprocal_components`.
        """

        return direct_to_reciprocal_components(components, self)

    def reciprocal_to_direct_components(self, components: Any) -> np.ndarray:
        """Convert reciprocal-basis components to direct-basis components in this
        lattice; see :func:`reciprocal_to_direct_components`.
        """

        return reciprocal_to_direct_components(components, self)


def metric_tensor(lattice: Lattice) -> np.ndarray:
    """The direct metric tensor ``G = A^T A`` of a lattice.

    Purpose
    -------
    ``G`` encodes the whole geometry of a lattice: the dot product of two
    direct-basis index vectors is ``u^T G v``, so lengths, angles, and
    volumes in a triclinic cell follow from one matrix instead of a
    per-system formula. Its diagonal holds ``a^2, b^2, c^2`` and its
    off-diagonal entries the products ``ab cos(gamma)`` and so on.

    Returns
    -------
    np.ndarray
        ``(3, 3)`` symmetric positive-definite tensor in angstrom squared,
        read-only.
    """

    basis = lattice.direct_basis().matrix
    tensor = basis.T @ basis
    tensor = np.ascontiguousarray(tensor, dtype=np.float64)
    tensor.setflags(write=False)
    return tensor


def reciprocal_metric_tensor(lattice: Lattice) -> np.ndarray:
    """The reciprocal metric tensor ``G* = B^T B`` of a lattice.

    The reciprocal-space counterpart of :func:`metric_tensor`: the dot
    product of two reciprocal-basis index vectors is ``h^T G* k``, which is
    where the general interplanar-spacing formula
    ``1/d^2 = h^T G* h`` comes from. Units are inverse angstrom squared;
    returned read-only.
    """

    basis = lattice.reciprocal_basis().matrix
    tensor = basis.T @ basis
    tensor = np.ascontiguousarray(tensor, dtype=np.float64)
    tensor.setflags(write=False)
    return tensor


def direct_to_reciprocal_components(components: Any, lattice: Lattice) -> np.ndarray:
    """Convert direct-basis components to reciprocal-basis components.

    Purpose
    -------
    Index-lowering with the metric tensor. The direction ``[uvw]`` and the
    plane ``(hkl)`` with the same numbers are *not* the same geometric
    object outside the cubic system; this is the correct conversion between
    the two component sets, and the reason PyTex never equates them.

    Parameters
    ----------
    components : ArrayLike
        Any array ending in dimension 3.
    lattice : Lattice

    Returns
    -------
    np.ndarray
        Same shape, in reciprocal-basis components; read-only.
    """

    array = np.asarray(components, dtype=np.float64)
    if array.shape[-1] != 3:
        raise ValueError("direct components must end with dimension 3.")
    transformed = np.asarray(array, dtype=np.float64) @ metric_tensor(lattice)
    transformed = np.ascontiguousarray(transformed)
    transformed.setflags(write=False)
    return transformed


def reciprocal_to_direct_components(components: Any, lattice: Lattice) -> np.ndarray:
    """Convert reciprocal-basis components to direct-basis components.

    Index-raising with the reciprocal metric tensor; the inverse operation
    of :func:`direct_to_reciprocal_components`.
    """

    array = np.asarray(components, dtype=np.float64)
    if array.shape[-1] != 3:
        raise ValueError("reciprocal components must end with dimension 3.")
    transformed = np.asarray(array, dtype=np.float64) @ reciprocal_metric_tensor(lattice)
    transformed = np.ascontiguousarray(transformed)
    transformed.setflags(write=False)
    return transformed


@dataclass(frozen=True, slots=True)
class AtomicSite:
    """One atom in a unit cell: species, fractional position, and occupancy.

    Attributes
    ----------
    label : str
        Site label, as used in CIF files.
    species : str
        Element symbol; sets the scattering factor.
    fractional_coordinates : np.ndarray
        Position in fractional cell coordinates.
    occupancy : float
        Fractional site occupancy, for disordered or partially filled sites.
    """

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
    """A lattice together with the atomic sites in one cell.

    Purpose
    -------
    The atomic basis that structure-factor calculations sum over, and the
    input to crystal-structure visualization. A phase without a unit cell
    supports lattice geometry but not intensity calculation.

    Attributes
    ----------
    lattice : Lattice
    sites : tuple of AtomicSite
        Atoms in the cell, in fractional coordinates.
    provenance : ProvenanceRecord, optional
    """

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
        """Build a unit cell with atomic sites from a pymatgen ``Structure``.

        Optional-dependency bridge that carries fractional coordinates, species,
        and site labels across. Pass ``lattice`` to reuse an already-constructed
        PyTex lattice instead of deriving one from the structure.
        """

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
    """A crystallographic phase: lattice, symmetry, frame, and optional structure.

    Purpose
    -------
    The unit of material identity in PyTex. Anything that depends on *which
    material* is being described — d-spacings, structure factors, symmetry
    orbits, index conversions, transformation relationships — takes a phase
    rather than loose parameters, so those pieces cannot become inconsistent
    with one another.

    Construction enforces internal consistency: the lattice, symmetry, unit
    cell, and space group must all agree on the crystal frame and on each
    other, and a mismatch raises rather than propagating.

    Attributes
    ----------
    name : str
        Canonical phase name.
    lattice : Lattice
    symmetry : SymmetrySpec
        Crystal symmetry, in ``crystal_frame``.
    crystal_frame : ReferenceFrame
        Must belong to the crystal domain.
    unit_cell : UnitCell, optional
        Atomic basis. Required for structure factors; without it, only
        lattice-level reasoning is available.
    space_group : SpaceGroupSpec, optional
    space_group_symbol : str, optional
    space_group_number : int, optional
        Kept consistent with ``space_group`` when both are given. The space
        group determines centring absences, so omitting it means forbidden
        reflections will be reported as present.
    chemical_formula : str, optional
    aliases : tuple of str
        Alternative names, for matching against vendor phase names.
    provenance : ProvenanceRecord, optional
    """

    name: str
    lattice: Lattice
    symmetry: SymmetrySpec
    crystal_frame: ReferenceFrame
    unit_cell: UnitCell | None = None
    space_group: SpaceGroupSpec | None = None
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
        if self.space_group is not None and self.space_group.reference_frame != self.crystal_frame:
            raise ValueError("Phase.space_group.reference_frame must match Phase.crystal_frame.")
        if self.space_group is not None:
            if self.space_group_symbol is None:
                object.__setattr__(self, "space_group_symbol", self.space_group.symbol)
            elif self.space_group_symbol != self.space_group.symbol:
                raise ValueError("Phase.space_group_symbol must match Phase.space_group.symbol.")
            if self.space_group_number is None:
                object.__setattr__(self, "space_group_number", self.space_group.number)
            elif self.space_group_number != self.space_group.number:
                raise ValueError("Phase.space_group_number must match Phase.space_group.number.")
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
        """Build a phase from a pymatgen ``Structure``.

        Purpose
        -------
        The route from an external structural model to a fully specified PyTex
        phase: lattice, unit cell, symmetry, and space group are all derived
        together, so the phase is internally consistent by construction.

        Parameters
        ----------
        structure : pymatgen Structure
        crystal_frame : ReferenceFrame
            The crystal-domain frame to attach.
        phase_name : str, optional
            Defaults to the structure's reduced formula.
        aliases : tuple of str
            Additional names this phase should answer to, for matching against
            vendor phase names.
        symprec, angle_tolerance : float
            Symmetry-detection tolerances passed to the spglib-backed analysis.
            Loose values can promote a distorted cell to a higher symmetry than
            the data supports; they are exposed so that choice is explicit.

        Returns
        -------
        Phase
        """

        _, spacegroup_analyzer_cls = _require_pymatgen()
        analyzer = spacegroup_analyzer_cls(
            structure,
            symprec=float(symprec),
            angle_tolerance=float(angle_tolerance),
        )
        with _spglib_dict_shim_silenced():
            point_group = str(analyzer.get_point_group_symbol())
        symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal_frame)
        space_group = SpaceGroupSpec.from_pymatgen_analyzer(
            analyzer,
            reference_frame=crystal_frame,
            provenance=provenance,
        )
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
            space_group=space_group,
            space_group_symbol=space_group.symbol,
            space_group_number=space_group.number,
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
        """Build a phase from a CIF file.

        Purpose
        -------
        The most common entry point for real structural data. Requires the
        optional pymatgen dependency; the CIF is parsed, its symmetry
        determined, and the result attached to the given crystal frame.

        Parameters
        ----------
        path : str or Path
        crystal_frame : ReferenceFrame
        phase_name : str, optional
            Defaults to the reduced formula.
        aliases : tuple of str
            Additional names for phase matching.
        primitive : bool
            Reduce to the primitive cell. Off by default, so the conventional
            cell of the CIF is kept and reflection conditions stay those of the
            conventional setting.
        symprec, angle_tolerance : float
            Symmetry-detection tolerances.

        Returns
        -------
        Phase

        See Also
        --------
        from_cif_string : The same construction from CIF text already in memory.
        """

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
        """Build a phase from CIF text held in memory.

        Identical in contract to :meth:`from_cif`, for CIF content that came
        from a database query, an archive, or a test fixture rather than a file
        on disk.
        """

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


def phases_semantically_match(left: Phase | None, right: Phase | None) -> bool:
    """Whether two phases carry the same crystallographic identity.

    Purpose: the single, shared definition of phase-identity used by
    transformation, reconstruction, and workflow validation code when deciding
    whether two ``Phase`` objects describe the same material phase.

    When to use: prefer this over ``Phase == Phase`` whenever the two objects
    may have been constructed independently. Identity is defined as equal phase
    name, crystal frame, lattice parameters, and normalized point-group symbol;
    provenance, aliases, and unit-cell contents are deliberately excluded, and
    symmetry is compared through its normalized point-group symbol rather than
    operator arrays so equivalent specs constructed separately still match.

    Inputs: two ``Phase`` objects or ``None``. Two ``None`` values match; a
    ``None`` never matches a real phase.

    Output: ``bool``.
    """

    if left is None or right is None:
        return left is right
    return (
        left.name == right.name
        and left.crystal_frame == right.crystal_frame
        and left.lattice == right.lattice
        and normalize_point_group_symbol(left.symmetry.point_group)
        == normalize_point_group_symbol(right.symmetry.point_group)
    )


@dataclass(frozen=True, slots=True)
class MillerIndex:
    """An ``(hkl)`` index triple bound to a phase and a basis kind.

    Purpose
    -------
    The minimal typed index. It defaults to the reciprocal basis, since
    Miller indices are reciprocal-basis components — the distinction that
    makes ``[hkl]`` and ``(hkl)`` different objects outside the cubic system.

    For index algebra — families, symmetry orbits, reductions, angle
    tables — use :class:`~pytex.core.miller.MillerPlane` and its set form,
    which build on this.

    Attributes
    ----------
    indices : np.ndarray
        Integer index triple.
    phase : Phase
    basis_kind : BasisKind
        ``RECIPROCAL`` by default.
    """

    indices: np.ndarray
    phase: Phase
    basis_kind: BasisKind = BasisKind.RECIPROCAL

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("Miller indices must not be the zero triplet.")

    def as_array(self) -> np.ndarray:
        """The stored index triple as an integer array.
        """

        return self.indices

    @classmethod
    def from_miller_bravais(cls, indices: Any, *, phase: Phase) -> MillerIndex:
        """Miller index from a hexagonal four-index ``(hkil)`` quadruple.

        The constraint ``i = -(h + k)`` is checked on conversion.
        """

        return cls(indices=plane_hkil_to_hkl(indices), phase=phase)


@dataclass(frozen=True, slots=True)
class CrystalDirection:
    """A crystallographic direction ``[uvw]`` on a phase.

    Purpose
    -------
    A direction in the lattice, resolved to Cartesian coordinates through the
    *direct* basis. Keeping the phase attached is what makes the resolution
    correct for non-cubic lattices, where the Cartesian direction of
    ``[uvw]`` depends on the cell parameters.

    Attributes
    ----------
    coordinates : np.ndarray
        Direction components; must not be the zero vector.
    phase : Phase
    basis_kind : BasisKind
        ``DIRECT`` by default. A reciprocal-basis direction is resolved
        through the reciprocal basis instead, which is why the kind is
        stored rather than assumed.
    """

    coordinates: np.ndarray
    phase: Phase
    basis_kind: BasisKind = BasisKind.DIRECT

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", as_float_array(self.coordinates, shape=(3,)))
        if np.allclose(self.coordinates, 0.0):
            raise ValueError("CrystalDirection coordinates must not define the zero direction.")

    @property
    def unit_vector(self) -> np.ndarray:
        """Unit vector of this direction in the Cartesian crystal frame.

        Resolved through the direct basis for a direct-basis direction and
        through the reciprocal basis for a reciprocal-basis one, following the
        direction's declared ``basis_kind`` — which is why the basis is stored
        with the coordinates rather than assumed.
        """

        if self.basis_kind is BasisKind.DIRECT:
            basis = self.phase.lattice.direct_basis()
        else:
            basis = self.phase.lattice.reciprocal_basis()
        cartesian = basis.matrix @ self.coordinates
        return normalize_vector(cartesian)

    @classmethod
    def from_miller_bravais(cls, indices: Any, *, phase: Phase) -> CrystalDirection:
        """Direction from a hexagonal four-index ``[UVTW]`` quadruple.

        The constraint ``U + V + T = 0`` is checked on conversion.
        """

        return cls(coordinates=direction_uvtw_to_uvw(indices).astype(np.float64), phase=phase)

    @classmethod
    def from_cartesian(cls, vector: Any, *, phase: Phase) -> CrystalDirection:
        """Direction from a Cartesian vector expressed in the phase crystal frame.

        Purpose: the inverse of ``unit_vector`` — converts a crystal-frame
        Cartesian vector into direct-basis ``[uvw]`` coordinates so the stored
        direction keeps index meaning.

        Inputs: a nonzero 3-vector in the crystal Cartesian frame of ``phase``.

        Output: a ``CrystalDirection`` whose ``unit_vector`` reproduces the
        normalized input vector.
        """

        cartesian = as_float_array(vector, shape=(3,))
        if np.allclose(cartesian, 0.0):
            raise ValueError("CrystalDirection.from_cartesian requires a nonzero vector.")
        coordinates = np.linalg.solve(phase.lattice.direct_basis().matrix, cartesian)
        return cls(coordinates=coordinates, phase=phase)


@dataclass(frozen=True, slots=True)
class ZoneAxis:
    """A zone axis ``[uvw]``: the beam direction of a diffraction pattern.

    Purpose
    -------
    A lattice direction seen from the diffraction side. Reflections belong to
    its zone exactly when the zone law ``hu + kv + lw = 0`` holds, which is
    the selection rule determining what appears in a zone-axis pattern.

    Attributes
    ----------
    indices : np.ndarray
        Integer direction indices.
    phase : Phase
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("ZoneAxis indices must not be the zero triplet.")

    @property
    def direction(self) -> CrystalDirection:
        """This zone axis as a :class:`CrystalDirection`.

        A zone axis is a lattice direction; this is the same vector in the
        general direction type.
        """

        return CrystalDirection(self.indices.astype(np.float64), phase=self.phase)

    @property
    def unit_vector(self) -> np.ndarray:
        """Unit vector of the zone axis in the Cartesian crystal frame.

        For electron diffraction this is the direction that must be brought
        parallel to the beam for the corresponding zone-axis pattern.
        """

        return self.direction.unit_vector

    def zone_law_value(self, miller: MillerIndex) -> int:
        """The zone-law value ``hu + kv + lw`` for a reflection on this zone.

        Zero means the reflection belongs to the zone — that is, it appears in
        this zone-axis pattern. The phases of the operands must match.
        """

        if miller.phase != self.phase:
            raise ValueError("MillerIndex.phase must match ZoneAxis.phase.")
        return int(np.dot(self.indices, miller.indices))

    def contains_miller_index(self, miller: MillerIndex) -> bool:
        """Whether a reflection belongs to this zone (zone-law value zero).

        The selection rule that determines which reflections appear in a
        zone-axis diffraction pattern.
        """

        return self.zone_law_value(miller) == 0

    @classmethod
    def from_miller_bravais(cls, indices: Any, *, phase: Phase) -> ZoneAxis:
        """Zone axis from a hexagonal four-index ``[UVTW]`` quadruple.
        """

        return cls(indices=direction_uvtw_to_uvw(indices), phase=phase)


@dataclass(frozen=True, slots=True)
class ReciprocalLatticeVector:
    """A reciprocal-lattice vector ``g = h a* + k b* + l c*`` on a phase.

    Purpose
    -------
    The reciprocal-space view of a plane: its direction is the plane normal
    and its magnitude is ``1/d``. This is the natural object for diffraction
    reasoning, where the Ewald construction and the scattering condition are
    both statements about ``g``.

    Attributes
    ----------
    coordinates : np.ndarray
        Reciprocal-basis components; must not be the zero vector.
    phase : Phase

    Notes
    -----
    Uses the crystallographic convention without ``2*pi``.
    """

    coordinates: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", as_float_array(self.coordinates, shape=(3,)))
        if np.allclose(self.coordinates, 0.0):
            raise ValueError("ReciprocalLatticeVector coordinates must not be the zero vector.")

    @classmethod
    def from_miller_index(cls, miller: MillerIndex) -> ReciprocalLatticeVector:
        """The reciprocal-lattice vector ``g`` of a Miller index.

        The ``(hkl)`` components are reinterpreted as reciprocal-basis
        coordinates, which is what they already are — a Miller index *is* a
        reciprocal-lattice vector in index form.
        """

        return cls(coordinates=miller.indices.astype(np.float64), phase=miller.phase)

    @property
    def cartesian_vector(self) -> np.ndarray:
        """The vector ``g = h a* + k b* + l c*`` in Cartesian crystal-frame
        coordinates, in inverse angstroms.

        Not normalized: its magnitude carries the d-spacing information.
        """

        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        return as_float_array(reciprocal @ self.coordinates, shape=(3,))

    @property
    def magnitude_inv_angstrom(self) -> float:
        """``|g| = 1 / d``, in inverse angstroms.

        The reciprocal of the interplanar spacing, under the crystallographic
        (no ``2*pi``) reciprocal-lattice convention.
        """

        return float(np.linalg.norm(self.cartesian_vector))

    @property
    def unit_vector(self) -> np.ndarray:
        """Unit vector along ``g`` in the Cartesian crystal frame.

        The plane-normal direction of the corresponding ``(hkl)``.
        """

        return normalize_vector(self.cartesian_vector)


@dataclass(frozen=True, slots=True)
class CrystalPlane:
    """A crystallographic plane ``(hkl)`` on a phase.

    Purpose
    -------
    A lattice plane, resolved through the *reciprocal* basis — the only
    correct route for a normal outside the cubic system, where the direction
    ``[hkl]`` is not parallel to the normal of ``(hkl)``. Carries the
    d-spacing that Bragg's law consumes.

    Attributes
    ----------
    miller : MillerIndex
        The index triple; its phase must match ``phase``.
    phase : Phase
    """

    miller: MillerIndex
    phase: Phase

    def __post_init__(self) -> None:
        if self.miller.phase != self.phase:
            raise ValueError("CrystalPlane.miller.phase must match CrystalPlane.phase.")

    @property
    def normal(self) -> np.ndarray:
        """Unit normal of the plane in the Cartesian crystal frame.

        Resolved through the reciprocal basis, the only correct route outside
        the cubic system: for a non-cubic lattice the direction ``[hkl]`` is not
        parallel to the normal of the plane ``(hkl)``.
        """

        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        normal = reciprocal @ self.miller.indices.astype(np.float64)
        return normalize_vector(normal)

    @property
    def d_spacing_angstrom(self) -> float:
        """Interplanar spacing ``d`` in angstroms.

        The reciprocal of ``|g|``, correct for every crystal system without a
        per-system formula, and the quantity Bragg's law consumes.
        """

        reciprocal = self.phase.lattice.reciprocal_basis().matrix
        reciprocal_vector = reciprocal @ self.miller.indices.astype(np.float64)
        magnitude = float(np.linalg.norm(reciprocal_vector))
        if np.isclose(magnitude, 0.0):
            raise ValueError("CrystalPlane reciprocal vector magnitude must be non-zero.")
        return 1.0 / magnitude

    @property
    def reciprocal_lattice_vector(self) -> ReciprocalLatticeVector:
        """This plane as a :class:`ReciprocalLatticeVector`.

        The same object seen from reciprocal space, where diffraction reasoning
        happens.
        """

        return ReciprocalLatticeVector.from_miller_index(self.miller)

    @classmethod
    def from_miller_bravais(cls, indices: Any, *, phase: Phase) -> CrystalPlane:
        """Plane from a hexagonal four-index ``(hkil)`` quadruple.

        The constraint ``i = -(h + k)`` is checked on conversion.
        """

        return cls(miller=MillerIndex.from_miller_bravais(indices, phase=phase), phase=phase)
