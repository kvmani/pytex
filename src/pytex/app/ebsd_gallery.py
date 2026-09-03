"""Practice orientation maps with known answers.

Why a gallery
-------------
An EBSD panel that needs a scan file before it will show anything is unusable
for the two people most likely to open it: someone evaluating whether the tool
does what they need, and someone learning what a GROD map *is*. Worse, a real
scan comes with no answer key — a mistake in the workflow and a genuine feature
of the microstructure look identical on screen.

Each entry here is built from a stated construction, so the answer is known
before the calculation runs:

``bicrystal_gradient``
    Two grains meeting at a straight vertical boundary, misoriented by exactly
    40 degrees about ``[001]``. The right-hand grain carries a *linear*
    orientation gradient of a stated degrees-per-step about ``[010]``. Three
    quantities are therefore known in advance:

    - the boundary misorientation, exactly 40 degrees;
    - the KAM inside the gradient, which is **half** the per-step rotation at
      every interior point. Half, not all of it: the four-neighbour kernel
      averages two neighbours along the gradient, each a full step away, with
      two neighbours across it that are identical. A KAM is an average over a
      kernel, not a gradient magnitude, and this entry is where that stops being
      a technicality;
    - the GROD, which is the *absolute* deviation from the grain's own
      reference orientation. It therefore falls to zero at the reference point
      and rises linearly on both sides of it, rather than ramping across the
      grain.

    This is the entry to check a KAM or GROD pipeline against.

``sigma3_twin``
    A matrix grain crossed by two twin lamellae related to it by exactly 60
    degrees about ``[111]`` — the coherent twin of a cubic metal, the Sigma 3
    boundary. Every boundary segment in this map has the same misorientation,
    and it is a number quoted in every textbook.

``equiaxed_polycrystal``
    Twelve grains grown from seed points, each with its own orientation and a
    small orientation spread. The realistic case: many boundaries at many
    misorientations, which is what a boundary-character distribution or a grain
    size needs. Its answer key is the grain count and the seed orientations.

Every entry carries the three scalar channels a vendor file carries — confidence
index, fit and image quality — because the panel's greyscale modulation exists
to show them, and because they behave the way real ones do: confidence falls and
fit worsens at the boundaries, where a pattern overlaps two lattices.

These are constructions, not measurements. They contain no detector geometry, no
indexing, and no noise model beyond the stated spread; nothing here should be
read as a statement about how a real scan behaves.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import OrientationSet
from pytex.core.symmetry import SymmetrySpec
from pytex.ebsd.models import CrystalMap

__all__ = [
    "DEFAULT_ENTRY_ID",
    "GALLERY",
    "GalleryEntry",
    "build_map",
    "default_entry",
    "entry_ids",
    "get_entry",
]

#: Side of a gallery map when the caller does not choose one. Large enough that
#: grains have interiors and boundaries have length. The workbench asks for a
#: smaller map than this, because it draws one the moment the panel opens and
#: the IPF symmetry reduction dominates the cost; a script calling
#: :func:`build_map` directly is under no such deadline.
_GRID = 80

#: Step size in micrometres. Reported, and used for every length in the result,
#: so a grain diameter is in specimen units rather than in pixels.
_STEP_UM = 0.5


def _frames() -> tuple[ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    return crystal, specimen


def _nickel(crystal: ReferenceFrame) -> Phase:
    """Face-centred cubic nickel, the canonical cubic reference of this repo."""

    return Phase(
        name="Nickel (fcc)",
        lattice=Lattice(3.524, 3.524, 3.524, 90.0, 90.0, 90.0, crystal_frame=crystal),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
    )


@dataclass(frozen=True)
class GalleryEntry:
    """One practice map, and the answer it is built to have.

    Attributes
    ----------
    id : str
        Stable identifier used as the request value.
    title : str
        Name shown in the picker.
    summary : str
        One line: what this microstructure is.
    teaches : str
        The thing to notice once it is on screen.
    known_answer : str
        What is true of this map by construction, phrased so a reader can check
        the displayed numbers against it. This is what makes the entry a test
        rather than a picture.
    """

    id: str
    title: str
    summary: str
    teaches: str
    known_answer: str

    def describe(self) -> dict[str, str]:
        """Return the manifest-facing description of this entry."""

        return {
            "id": self.id,
            "title": self.title,
            "summary": self.summary,
            "teaches": self.teaches,
            "known_answer": self.known_answer,
        }


#: Degrees of lattice rotation per measurement step in ``bicrystal_gradient``.
#: Chosen at 0.35 so the KAM lands in the range a deformed metal actually shows
#: — a few tenths of a degree — rather than at a round number that would look
#: like an artefact of the construction.
GRADIENT_DEG_PER_STEP = 0.35

#: Misorientation across the bicrystal boundary, about [001].
BICRYSTAL_BOUNDARY_DEG = 40.0

#: The coherent twin of a cubic metal: 60 degrees about <111>.
SIGMA3_ANGLE_DEG = 60.0

GALLERY: tuple[GalleryEntry, ...] = (
    GalleryEntry(
        id="bicrystal_gradient",
        title="Bicrystal with a deformation gradient",
        summary=(
            "Two grains across a straight boundary; the right-hand one carries a linear "
            "orientation gradient."
        ),
        teaches=(
            "What KAM and GROD each measure. KAM is flat across the whole gradient, because a "
            "linear gradient presents the same step from every point to its neighbours; GROD "
            "rises away from the grain's reference orientation in both directions. The same "
            "microstructure therefore looks uniform in one map and strongly graded in the other, "
            "which is the distinction between local and grain-referenced misorientation."
        ),
        known_answer=(
            f"The boundary misorientation is exactly {BICRYSTAL_BOUNDARY_DEG:g} degrees about "
            f"[001]. Inside the gradient the KAM is {GRADIENT_DEG_PER_STEP / 2:g} degrees at "
            f"every interior point — half the {GRADIENT_DEG_PER_STEP:g} degree per-step rotation, "
            "because two of the four kernel neighbours lie across the gradient and are identical. "
            "The GROD falls to zero at the grain's reference orientation and rises linearly away "
            f"from it, reaching roughly a quarter of the map width in steps times "
            f"{GRADIENT_DEG_PER_STEP:g} degrees — so it scales with the map size, where the two "
            "misorientation figures above do not."
        ),
    ),
    GalleryEntry(
        id="sigma3_twin",
        title="Annealing twins in a cubic metal",
        summary="A matrix grain crossed by two coherent twin lamellae.",
        teaches=(
            "That a boundary map is a map of *misorientation*, not of contrast: the twin "
            "boundaries are as high-angle as any grain boundary, at 60 degrees, while the IPF "
            "colours on either side are related by a symmetry the eye does not read as similar."
        ),
        known_answer=(
            f"Every boundary segment is a Sigma 3 twin at exactly {SIGMA3_ANGLE_DEG:g} degrees "
            "about [111], so the boundary misorientation histogram is a single spike."
        ),
    ),
    GalleryEntry(
        id="equiaxed_polycrystal",
        title="Equiaxed polycrystal",
        summary="Twelve grains grown from seeds, each with its own orientation and a small spread.",
        teaches=(
            "The realistic case, and what the scalar channels are for: confidence index falls at "
            "every boundary, where a diffraction pattern overlaps two lattices. Modulating an IPF "
            "map by confidence draws that indexing quality straight onto the orientation map."
        ),
        known_answer=(
            "Twelve grains by construction, so a segmentation at any threshold below the smallest "
            "boundary misorientation must find twelve."
        ),
    ),
)


#: The dataset every EBSD view opens on unless told otherwise.
#:
#: Named rather than taken as ``GALLERY[0]``, because display order and default
#: are two different decisions and coupling them means reordering the gallery
#: silently changes what every map shows. The equiaxed polycrystal is the
#: default because it is the realistic case: twelve grains, a spread within each,
#: and quality channels that fall at the boundaries. A user arriving at the IPF
#: map, the KAM map or the grain-size distribution for the first time should see
#: a microstructure, not a two-grain construction built to isolate one teaching
#: point. The constructions stay one choice away, and each still states the known
#: answer that makes it checkable.
DEFAULT_ENTRY_ID = "equiaxed_polycrystal"


def default_entry() -> GalleryEntry:
    """The entry every EBSD operation defaults to."""

    return get_entry(DEFAULT_ENTRY_ID)


def entry_ids() -> tuple[str, ...]:
    """Identifiers of every gallery entry, in display order."""

    return tuple(entry.id for entry in GALLERY)


def get_entry(entry_id: str) -> GalleryEntry:
    """Return one entry.

    Raises
    ------
    KeyError
        If no entry uses that identifier.
    """

    for entry in GALLERY:
        if entry.id == entry_id:
            return entry
    raise KeyError(f"Unknown EBSD gallery entry {entry_id!r}; expected one of {entry_ids()}.")


def _grid_coordinates(rows: int, cols: int, step: float) -> np.ndarray:
    """Row-major point positions of a regular grid, in specimen units."""

    row_index, col_index = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    return np.column_stack(
        [col_index.ravel().astype(float) * step, row_index.ravel().astype(float) * step]
    )


def _axis_angle_quaternions(axes: np.ndarray, angles_deg: np.ndarray) -> np.ndarray:
    """Quaternions for a batch of axis-angle rotations, vectorised."""

    axes = np.asarray(axes, dtype=np.float64).reshape(-1, 3)
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
    half = np.deg2rad(np.asarray(angles_deg, dtype=np.float64).ravel()) * 0.5
    return np.column_stack([np.cos(half), axes * np.sin(half)[:, None]])


def _compose(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Hamilton product of two quaternion arrays, broadcast over the batch."""

    left = np.atleast_2d(left)
    right = np.atleast_2d(right)
    w1, x1, y1, z1 = (left[:, index] for index in range(4))
    w2, x2, y2, z2 = (right[:, index] for index in range(4))
    return np.column_stack(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ]
    )


def _bicrystal(rows: int, cols: int) -> tuple[np.ndarray, np.ndarray]:
    """Two grains, the right one linearly graded. Returns quaternions and labels."""

    labels = np.zeros((rows, cols), dtype=np.int64)
    labels[:, cols // 2 :] = 1

    left = _axis_angle_quaternions(np.array([[0.0, 0.0, 1.0]]), np.array([0.0]))
    right = _axis_angle_quaternions(
        np.array([[0.0, 0.0, 1.0]]), np.array([BICRYSTAL_BOUNDARY_DEG])
    )

    quaternions = np.repeat(left, rows * cols, axis=0)
    flat_labels = labels.ravel()
    # The gradient is a rotation about [010] whose angle grows with the column
    # index measured from the boundary. Linear in position, so the step between
    # neighbouring columns is constant -- which is exactly why KAM is flat here
    # and GROD is not.
    _, col_index = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    distance = np.clip(col_index.ravel() - cols // 2, 0, None).astype(float)
    gradient = _axis_angle_quaternions(
        np.tile(np.array([0.0, 1.0, 0.0]), (rows * cols, 1)),
        distance * GRADIENT_DEG_PER_STEP,
    )
    right_side = _compose(np.repeat(right, rows * cols, axis=0), gradient)
    quaternions[flat_labels == 1] = right_side[flat_labels == 1]
    return quaternions, flat_labels


def _sigma3(rows: int, cols: int) -> tuple[np.ndarray, np.ndarray]:
    """A matrix grain with two twin lamellae, all related by 60 degrees about [111]."""

    row_index, _ = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    flat_rows = row_index.ravel()
    # Two lamellae, each a sixth of the map tall, so the matrix survives between
    # and around them and every lamella has two boundaries rather than one.
    in_twin = ((flat_rows >= rows // 6) & (flat_rows < rows // 3)) | (
        (flat_rows >= 3 * rows // 5) & (flat_rows < 3 * rows // 5 + rows // 6)
    )
    labels = in_twin.astype(np.int64)

    matrix = _axis_angle_quaternions(np.array([[0.0, 0.0, 1.0]]), np.array([0.0]))
    twin = _axis_angle_quaternions(np.array([[1.0, 1.0, 1.0]]), np.array([SIGMA3_ANGLE_DEG]))
    quaternions = np.repeat(matrix, rows * cols, axis=0)
    quaternions[in_twin] = np.repeat(twin, int(in_twin.sum()), axis=0)
    return quaternions, labels


def _equiaxed(rows: int, cols: int, *, seed: int = 20260816) -> tuple[np.ndarray, np.ndarray]:
    """Twelve grains grown from seeds, each with a small orientation spread."""

    generator = np.random.default_rng(seed)
    grain_count = 12
    seeds = generator.integers(0, [rows, cols], size=(grain_count, 2))

    row_index, col_index = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    points = np.column_stack([row_index.ravel(), col_index.ravel()]).astype(float)
    # Nearest-seed assignment: the grains are Voronoi cells, which is the
    # standard idealisation of an equiaxed structure and gives straight
    # boundaries meeting near 120 degrees.
    distances = np.linalg.norm(points[:, None, :] - seeds[None, :, :], axis=2)
    labels = np.argmin(distances, axis=1).astype(np.int64)

    # Well separated grain orientations: evenly spaced angles about scattered
    # axes, so no two grains are accidentally within the segmentation threshold
    # of each other and the answer "twelve grains" holds.
    axes = generator.normal(size=(grain_count, 3))
    angles = np.linspace(12.0, 165.0, grain_count)
    grain_quaternions = _axis_angle_quaternions(axes, angles)

    quaternions = grain_quaternions[labels]
    # A small intragranular spread, so GROD and KAM are non-zero but far below
    # any boundary: real grains are never perfectly uniform.
    spread = _axis_angle_quaternions(
        generator.normal(size=(rows * cols, 3)),
        generator.uniform(0.0, 0.6, size=rows * cols),
    )
    return _compose(quaternions, spread), labels


def _quality_channels(labels: np.ndarray, rows: int, cols: int) -> dict[str, np.ndarray]:
    """Confidence, fit and image quality that behave the way real ones do.

    At a boundary the diffraction pattern is a superposition of two lattices, so
    the indexing confidence falls and the fit worsens. That is the whole reason
    the panel can modulate a colour map by a scalar channel, so the channels
    have to carry the effect rather than being decorative noise.
    """

    grid = labels.reshape(rows, cols)
    # A point is "near a boundary" when any of its four neighbours belongs to
    # another grain. Computed by comparison with shifted copies, which is exact
    # and needs no loop.
    edge = np.zeros((rows, cols), dtype=bool)
    edge[:-1, :] |= grid[:-1, :] != grid[1:, :]
    edge[1:, :] |= grid[:-1, :] != grid[1:, :]
    edge[:, :-1] |= grid[:, :-1] != grid[:, 1:]
    edge[:, 1:] |= grid[:, :-1] != grid[:, 1:]
    edge_flat = edge.ravel()

    generator = np.random.default_rng(7)
    confidence = np.where(edge_flat, 0.22, 0.86) + generator.normal(0.0, 0.03, labels.size)
    fit = np.where(edge_flat, 2.1, 0.55) + generator.normal(0.0, 0.08, labels.size)
    quality = np.where(edge_flat, 320.0, 980.0) + generator.normal(0.0, 25.0, labels.size)
    return {
        "confidence_index": np.clip(confidence, 0.0, 1.0),
        "fit": np.clip(fit, 0.0, None),
        "image_quality": np.clip(quality, 0.0, None),
    }


_BUILDERS: dict[str, Any] = {
    "bicrystal_gradient": _bicrystal,
    "sigma3_twin": _sigma3,
    "equiaxed_polycrystal": _equiaxed,
}


def build_map(entry_id: str, *, grid: int = _GRID, step_um: float = _STEP_UM) -> CrystalMap:
    """Build one gallery map.

    Purpose
    -------
    Returns a fully formed :class:`~pytex.ebsd.models.CrystalMap` — orientations,
    a regular square grid with a physical step, a declared fcc phase, and the
    three scalar channels — so every downstream calculation the panel offers
    (segmentation, boundaries, KAM, GROD, IPF colouring) runs on it exactly as
    it would on an imported scan.

    Parameters
    ----------
    entry_id : str
        One of :func:`entry_ids`.
    grid : int
        Side of the square map, in measurement points.
    step_um : float
        Step size in micrometres, used for coordinates and hence for every
        reported length.

    Returns
    -------
    CrystalMap

    Raises
    ------
    KeyError
        If ``entry_id`` names no entry.

    Examples
    --------
    >>> crystal_map = build_map("sigma3_twin", grid=16)
    >>> crystal_map.grid_shape
    (16, 16)
    >>> crystal_map.property_names
    ('confidence_index', 'fit', 'image_quality')
    """

    get_entry(entry_id)
    crystal, specimen = _frames()
    phase = _nickel(crystal)
    quaternions, labels = _BUILDERS[entry_id](grid, grid)
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
        symmetry=phase.symmetry,
    )
    return CrystalMap(
        coordinates=_grid_coordinates(grid, grid, step_um),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(grid, grid),
        grid_kind="square",
        step_sizes=(step_um, step_um),
        properties=_quality_channels(labels, grid, grid),
    )
