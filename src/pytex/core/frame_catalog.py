"""The catalog of standard reference frames every PyTex workflow can reach for by name.

Most workflows need the same handful of frames: the canonical Cartesian triad
``X, Y, Z``; the specimen frame ``x, y, z``; the rolling-geometry sample frame
``RD, TD, ND``; a crystal frame ``a, b, c``; an EBSD scan (map) frame; a
detector frame; a laboratory frame; and the reciprocal frame attached to a
crystal frame. Building those by hand in each module is how convention drift
starts, so this module builds them **once**.

Two ways to use it
------------------

- **Constants** — `CARTESIAN_FRAME`, `SPECIMEN_FRAME`, `SAMPLE_RD_TD_ND_FRAME`,
  `CRYSTAL_FRAME`, `MAP_FRAME`, `DETECTOR_FRAME`, `LABORATORY_FRAME` are ready
  to use directly and compare equal wherever they appear.
- **Builders** — `cartesian_frame`, `specimen_frame`, `sample_frame`,
  `crystal_frame`, `map_frame`, `detector_frame`, `laboratory_frame`, and
  `reciprocal_frame_for` accept a name, axis geometry, and provenance for the
  cases where one workflow holds several frames of the same kind (two phases,
  two detectors, a parent and a child crystal).

Identity preservation
---------------------

Frame equality is load-bearing across PyTex: `pytex.core.batches.VectorSet`,
`FrameTransform`, `Orientation`, and `SymmetrySpec` all compare frames. The
builder defaults here are pinned to exactly the field values that the modules
of this repository already used before the catalog existed (``crystal`` with
``a, b, c``; ``specimen`` with ``x, y, z``; ``map`` with ``x, y, z``), so
adopting the catalog never changes a frame's identity. That invariant is
asserted directly in ``tests/unit/test_frame_catalog.py``.

Domain vocabulary
-----------------

Every frame here uses a member of the fixed `FrameDomain` vocabulary
(``crystal, specimen, map, detector, laboratory, reciprocal``). New domains may
not be invented (`docs/standards/notation_and_conventions.md`).

See also
--------
`pytex.core.frames` : the underlying `ReferenceFrame`/`FrameTransform`/`FrameGraph` model.
`pytex.plotting.frames` : rendering any of these frames as a triad, gizmo, or SVG.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import (
    IDENTITY_AXIS_VECTORS,
    FrameGraph,
    FrameTransform,
    ReferenceFrame,
    as_axis_vectors,
)
from pytex.core.notation import format_reciprocal_axis_labels
from pytex.core.provenance import ProvenanceRecord

__all__ = [
    "CARTESIAN_FRAME",
    "CRYSTAL_FRAME",
    "DETECTOR_FRAME",
    "LABORATORY_FRAME",
    "MAP_FRAME",
    "SAMPLE_RD_TD_ND_FRAME",
    "SPECIMEN_FRAME",
    "STANDARD_FRAMES",
    "cartesian_frame",
    "crystal_frame",
    "detector_frame",
    "get_standard_frame",
    "laboratory_frame",
    "list_standard_frames",
    "map_frame",
    "reciprocal_frame_for",
    "rolling_frame_graph",
    "sample_frame",
    "specimen_frame",
]

_AxisLabels = Sequence[str]


def _axis_labels(labels: _AxisLabels, default: tuple[str, str, str]) -> tuple[str, str, str]:
    """Validate a three-label axis tuple, falling back to ``default``."""

    resolved = tuple(str(label) for label in labels) if labels else default
    if len(resolved) != 3:
        raise ValueError("A reference frame requires exactly three axis labels.")
    return (resolved[0], resolved[1], resolved[2])


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #


def cartesian_frame(
    name: str = "cartesian",
    *,
    axes: _AxisLabels = ("X", "Y", "Z"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = (
        "Canonical right-handed Cartesian reference in which every PyTex frame's axis "
        "vectors are expressed."
    ),
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the canonical Cartesian frame ``X, Y, Z``.

    When to use it
        As the neutral reference every other frame's `ReferenceFrame.axis_vectors`
        are quoted in, and as the world frame of 3D scenes. It is placed in the
        ``laboratory`` domain because it is instrument-fixed rather than attached
        to the crystal or the specimen.

    Returns
    -------
    ReferenceFrame
        A right-handed Cartesian frame with unit axes.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.LABORATORY,
        axes=_axis_labels(axes, ("X", "Y", "Z")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
        axis_descriptions=("Cartesian X axis", "Cartesian Y axis", "Cartesian Z axis"),
    )


def specimen_frame(
    name: str = "specimen",
    *,
    axes: _AxisLabels = ("x", "y", "z"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = "",
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the specimen frame ``x, y, z``.

    When to use it
        As the target frame of an `pytex.core.orientation.Orientation` (PyTex
        orientations map crystal to specimen), and as the frame texture and EBSD
        results are interpreted in. Use `sample_frame` instead when the
        rolling-geometry labels ``RD, TD, ND`` carry the meaning.

    Notes
    -----
    The defaults reproduce the specimen frame this repository used before the
    catalog existed, so this frame compares equal to those constructions.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.SPECIMEN,
        axes=_axis_labels(axes, ("x", "y", "z")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
    )


def sample_frame(
    name: str = "sample_rd_td_nd",
    *,
    axes: _AxisLabels = ("RD", "TD", "ND"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = (
        "Rolling-geometry sample frame: RD along the rolling direction, TD transverse to it "
        "in the sheet plane, ND along the sheet normal."
    ),
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the rolling-geometry sample frame ``RD, TD, ND``.

    What it does
        Names the specimen axes the way sheet-texture work names them: the
        rolling direction, the transverse direction, and the sheet normal
        direction, in that right-handed order.

    When to use it
        For rolled-sheet texture work — pole figures quoted as ``RD`` at the top
        and ``TD`` to the right, fibre and component definitions, and any figure
        whose axes a metallurgist expects to read as RD/TD/ND rather than x/y/z.

    Parameters
    ----------
    name:
        Frame name; give distinct names when several sample frames coexist.
    axes:
        Axis labels, defaulting to ``("RD", "TD", "ND")``.
    axis_vectors:
        Where the three axes point, in canonical Cartesian components. The
        default identity triad means ``RD = X``, ``TD = Y``, ``ND = Z``. Pass
        explicit vectors (or use `rolling_frame_graph`) to record a sample
        mounted at an angle.

    Returns
    -------
    ReferenceFrame
        A specimen-domain frame carrying the long axis names as
        `ReferenceFrame.axis_descriptions`.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.SPECIMEN,
        axes=_axis_labels(axes, ("RD", "TD", "ND")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
        axis_descriptions=(
            "rolling direction",
            "transverse direction",
            "normal direction",
        ),
    )


def crystal_frame(
    name: str = "crystal",
    *,
    axes: _AxisLabels = ("a", "b", "c"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = "",
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build a crystal-attached frame ``a, b, c``.

    When to use it
        As the source frame of an orientation, and as the home of crystal
        directions, planes, and symmetry operators. Give distinct names when a
        workflow holds more than one phase (``"parent"``/``"child"``, or the
        phase name).

    Notes
    -----
    The axis *lengths and angles* of a real unit cell live in
    `pytex.core.lattice.Lattice.direct_basis`, not here: this frame carries axis
    orientation only, and defaults to the Cartesian triad.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.CRYSTAL,
        axes=_axis_labels(axes, ("a", "b", "c")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
    )


def map_frame(
    name: str = "map",
    *,
    axes: _AxisLabels = ("x", "y", "z"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = "",
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the EBSD scan (map) frame.

    When to use it
        For scan-grid layout and neighbour topology. PyTex deliberately keeps the
        map frame distinct from the specimen frame: they coincide only when a
        workflow declares that relationship
        (`docs/standards/notation_and_conventions.md`).
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.MAP,
        axes=_axis_labels(axes, ("x", "y", "z")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
    )


def detector_frame(
    name: str = "detector",
    *,
    axes: _AxisLabels = ("u", "v", "n"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = "",
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build a detector frame ``u, v, n``.

    What it does
        Names the two in-plane detector axes ``u`` and ``v`` and the
        detector-plane normal ``n``, matching how PyTex projects SAED spots and
        EBSD patterns onto an image plane.

    When to use it
        For anything expressed in image or pattern coordinates: simulated
        diffractograms, pattern centres, projected traces. Keep it separate from
        specimen and crystal semantics — that separation is exactly what makes a
        diffraction figure's coordinates unambiguous.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.DETECTOR,
        axes=_axis_labels(axes, ("u", "v", "n")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
        axis_descriptions=(
            "in-plane detector axis",
            "in-plane detector axis",
            "detector-plane normal",
        ),
    )


def laboratory_frame(
    name: str = "laboratory",
    *,
    axes: _AxisLabels = ("x_lab", "y_lab", "z_lab"),
    axis_vectors: ArrayLike = IDENTITY_AXIS_VECTORS,
    description: str = (
        "Instrument-fixed laboratory frame. The beam or optical axis is not assumed: a "
        "workflow that needs it must declare which axis carries it."
    ),
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the instrument-fixed laboratory frame.

    When to use it
        As the common ground between a specimen stage and a detector in a
        diffraction geometry. PyTex does not silently assign the beam direction
        to one of its axes; declare that relationship explicitly so the geometry
        stays auditable.
    """

    return ReferenceFrame(
        name=str(name),
        domain=FrameDomain.LABORATORY,
        axes=_axis_labels(axes, ("x_lab", "y_lab", "z_lab")),
        handedness=Handedness.RIGHT,
        description=description,
        provenance=provenance,
        axis_vectors=as_axis_vectors(axis_vectors),
    )


def reciprocal_frame_for(
    frame: ReferenceFrame,
    *,
    provenance: ProvenanceRecord | None = None,
) -> ReferenceFrame:
    """Build the reciprocal frame attached to a crystal frame.

    What it does
        Returns a ``reciprocal``-domain frame named ``"<crystal>_reciprocal"``
        whose axis labels carry the IUCr reciprocal star (``a -> a*``), applied
        through `pytex.core.notation.format_reciprocal_axis_labels` so the
        repository has one starring rule, and inheriting the crystal frame's
        handedness and convention set.

    When to use it
        Wherever reciprocal-space quantities are produced —
        `pytex.core.lattice.Lattice.reciprocal_basis` uses it — so that a
        reciprocal-space vector can never be mistaken for a direct-space one.

    Parameters
    ----------
    frame:
        The crystal frame the reciprocal frame is dual to.
    provenance:
        Optional record, normally inherited from the lattice.

    Returns
    -------
    ReferenceFrame
        The dual frame. Its numerical reciprocal basis (with units of inverse
        angstrom) lives in `pytex.core.lattice.Basis`, not on the frame.
    """

    starred = format_reciprocal_axis_labels(frame.axes)
    return ReferenceFrame(
        name=f"{frame.name}_reciprocal",
        domain=FrameDomain.RECIPROCAL,
        axes=(starred[0], starred[1], starred[2]),
        handedness=frame.handedness,
        convention=frame.convention,
        description=f"Reciprocal basis associated with {frame.name}.",
        provenance=provenance,
        axis_descriptions=tuple(
            f"reciprocal axis dual to {label}" for label in frame.axes
        ),
    )


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

#: The canonical Cartesian reference ``X, Y, Z``.
CARTESIAN_FRAME: ReferenceFrame = cartesian_frame()

#: The default specimen frame ``x, y, z``.
SPECIMEN_FRAME: ReferenceFrame = specimen_frame()

#: The rolling-geometry sample frame ``RD, TD, ND``.
SAMPLE_RD_TD_ND_FRAME: ReferenceFrame = sample_frame()

#: The default crystal frame ``a, b, c``.
CRYSTAL_FRAME: ReferenceFrame = crystal_frame()

#: The default EBSD scan frame.
MAP_FRAME: ReferenceFrame = map_frame()

#: The default detector frame ``u, v, n``.
DETECTOR_FRAME: ReferenceFrame = detector_frame()

#: The default instrument-fixed laboratory frame.
LABORATORY_FRAME: ReferenceFrame = laboratory_frame()

#: Every standard frame, keyed by a short slug.
STANDARD_FRAMES: Mapping[str, ReferenceFrame] = MappingProxyType(
    {
        "cartesian": CARTESIAN_FRAME,
        "specimen": SPECIMEN_FRAME,
        "sample": SAMPLE_RD_TD_ND_FRAME,
        "crystal": CRYSTAL_FRAME,
        "map": MAP_FRAME,
        "detector": DETECTOR_FRAME,
        "laboratory": LABORATORY_FRAME,
    }
)


def list_standard_frames() -> tuple[str, ...]:
    """The slugs accepted by `get_standard_frame`, in a stable order."""

    return tuple(STANDARD_FRAMES)


def get_standard_frame(key: str) -> ReferenceFrame:
    """Look up a standard frame by slug.

    Parameters
    ----------
    key:
        One of the slugs from `list_standard_frames`: ``"cartesian"``,
        ``"specimen"``, ``"sample"``, ``"crystal"``, ``"map"``, ``"detector"``,
        ``"laboratory"``. Matching is case-insensitive.

    Returns
    -------
    ReferenceFrame
        The shared catalog frame for that slug.

    Raises
    ------
    KeyError
        If the slug is not in the catalog.
    """

    slug = str(key).strip().lower()
    try:
        return STANDARD_FRAMES[slug]
    except KeyError:
        available = ", ".join(list_standard_frames())
        raise KeyError(
            f"'{key}' is not a standard PyTex frame. Available frames: {available}."
        ) from None


# --------------------------------------------------------------------------- #
# Ready-made graphs
# --------------------------------------------------------------------------- #


def rolling_frame_graph(
    *,
    rd_offset_deg: float = 0.0,
    name: str = "rolling_geometry",
) -> FrameGraph:
    """A `FrameGraph` for sheet-texture geometry: Cartesian, specimen, and RD/TD/ND.

    What it does
        Registers the canonical Cartesian frame, the specimen frame, and a
        rolling-geometry sample frame, and declares the two relationships between
        them. The sample frame is rotated about ``ND`` by ``rd_offset_deg``
        relative to the specimen ``x`` axis, which is how a sheet mounted at an
        angle to the stage is recorded honestly rather than assumed away.

    When to use it
        As a worked starting point for rolled-sheet workflows, and as the example
        the documentation and tests use to demonstrate multi-hop resolution:
        asking for ``cartesian -> sample_rd_td_nd`` composes both declared edges
        automatically.

    Parameters
    ----------
    rd_offset_deg:
        Angle in degrees from the specimen ``x`` axis to ``RD``, measured about
        ``ND`` in the right-handed sense. Zero (the default) means RD is aligned
        with specimen ``x``.
    name:
        Graph name used in `FrameGraph.describe`.

    Returns
    -------
    FrameGraph
        A graph with three frames and two declared relationships.

    Examples
    --------
    A sample mounted 30 degrees off the stage axis, asked for the full chain::

        graph = rolling_frame_graph(rd_offset_deg=30.0)
        transform = graph.transform_between("cartesian", "sample_rd_td_nd")
        transform.rotation_angle_deg  # 30.0
    """

    angle_rad = float(np.deg2rad(rd_offset_deg))
    cosine = float(np.cos(angle_rad))
    sine = float(np.sin(angle_rad))
    # RD/TD rotated about ND (= Z) by rd_offset_deg, expressed in canonical Cartesian.
    sample = sample_frame(
        axis_vectors=(
            (cosine, sine, 0.0),
            (-sine, cosine, 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    graph = FrameGraph(name=name)
    graph.add_transform(
        FrameTransform.between_frames(CARTESIAN_FRAME, SPECIMEN_FRAME),
    )
    graph.add_transform(
        FrameTransform.between_frames(SPECIMEN_FRAME, sample),
    )
    return graph
