# ruff: noqa: RUF001
"""EBSD orientation maps for the shared web and desktop workbench.

One operation, four orthogonal choices
--------------------------------------
An EBSD plotter could be a dozen operations — an IPF map, a KAM map, a boundary
map, a CI map — and it would be the wrong shape. What a user actually varies is
four independent things:

1. **What the colour means.** IPF along a chosen specimen direction, grain
   identity, or a scalar field rendered through a colour map: GROD, KAM, or a
   measured channel such as confidence index or fit.
2. **Whether a scalar modulates it.** Any coloured map can be darkened by any
   scalar channel, which is how an IPF map is made to show where the indexing
   was poor without giving up the orientation information.
3. **Whether the boundaries are drawn on top.** Any map at all, with the
   boundary network superimposed and classified into low- and high-angle.
4. **What counts as a grain**, which is one misorientation threshold and
   changes every derived quantity beneath it.

Declaring them as four parameters of :func:`ebsd.map` rather than as a dozen
operations is what makes "grain boundaries superimposed on a GROD map, greyed by
confidence index" reachable — a combination nobody would have thought to
enumerate, and the one a real analysis asks for.

How the map is delivered
------------------------
As a raster, because it *is* one: a crystal map on a regular grid has one
measurement per pixel, so the image is the data at its native resolution with no
interpolation anywhere. The boundary network travels separately as line segments
in map coordinates, since boundaries are geometry rather than pixels and must
stay sharp when the figure is zoomed.
"""

from __future__ import annotations

import base64
from typing import Any

import numpy as np

from pytex.adapters.scan_files import SCAN_FILE_SUFFIXES
from pytex.app.ebsd_gallery import (
    DEFAULT_ENTRY_ID,
    GALLERY,
    build_map,
    entry_ids,
    get_entry,
)
from pytex.app.errors import InvalidInputError, UnsupportedRequestError
from pytex.app.logbook import APP_LOG
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    DocumentationLink,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    Parameter,
)
from pytex.app.results import AppResult, Column, ResultTable

__all__: tuple[str, ...] = ()

_CITATION_RANDLE_ENGLER = (
    "Randle & Engler, Introduction to Texture Analysis, 2nd ed., CRC Press 2010, Chs. 6-7."
)
_CITATION_WRIGHT_KAM = (
    "Wright, Nowell & Field, Microsc. Microanal. 17 (2011) 316, doi:10.1017/S1431927611000055."
)
_CITATION_NOLZE_IPF = (
    "Nolze & Hielscher, J. Appl. Crystallogr. 49 (2016) 1786, doi:10.1107/S1600576716012942."
)

#: Scalar channels a map can be coloured or modulated by. The key is the request
#: value; the tuple is (channel name on the CrystalMap, label, unit).
_SCALAR_CHANNELS = {
    "confidence_index": ("confidence_index", "Confidence index", ""),
    "fit": ("fit", "Fit", "°"),
    "image_quality": ("image_quality", "Image quality", ""),
}

_GRAIN_COLUMNS = (
    Column("grain_id", "Grain", numeric=True),
    Column("size", "Points", numeric=True),
    Column(
        "equivalent_diameter_um",
        "Equivalent diameter",
        units="µm",
        numeric=True,
        digits=3,
        help_text="Diameter of the circle with the same area as the grain.",
    ),
    Column("area_um2", "Area", units="µm²", numeric=True, digits=3),
    Column(
        "grain_orientation_spread_deg",
        "GOS",
        units="°",
        numeric=True,
        digits=3,
        help_text=(
            "Grain orientation spread: the mean deviation of the grain's points from its own "
            "reference orientation. One number per grain, where GROD is one per pixel."
        ),
    ),
    Column("mean_kam_deg", "Mean KAM", units="°", numeric=True, digits=4),
    Column(
        "mean_phi1_deg",
        "phi1",
        units="°",
        numeric=True,
        digits=3,
        help_text=(
            "Bunge phi1 of the grain's symmetry-aware mean orientation — the average over its "
            "points, not the single reference point GROD is measured against. Together with Phi "
            "and phi2 this is what the measured-pair views of the Variants workspace take, so a "
            "relationship between two grains needs no retyping. Read the GOS column beside it: a "
            "mean has no scatter of its own, and a relationship computed from two means reports a "
            "zero residual however noisy the grains are."
        ),
    ),
    Column("mean_Phi_deg", "Phi", units="°", numeric=True, digits=3),
    Column("mean_phi2_deg", "phi2", units="°", numeric=True, digits=3),
    Column(
        "phase_name",
        "Phase",
        help_text="Which declared phase the grain's points carry; blank if the map declares none.",
    ),
)

_BOUNDARY_COLUMNS = (
    Column("character", "Character"),
    Column("count", "Segments", numeric=True),
    Column("length_um", "Total length", units="µm", numeric=True, digits=3),
    Column("fraction", "Fraction of length", numeric=True, digits=4),
    Column("mean_misorientation_deg", "Mean misorientation", units="°", numeric=True, digits=3),
)


def _colouring_parameter() -> ChoiceParameter:
    return ChoiceParameter(
        name="colouring",
        label="Colour by",
        help_text=(
            "What the colour of a pixel means. IPF colouring shows which crystal direction is "
            "along the chosen specimen axis, and is the standard orientation map. The scalar "
            "choices render a computed or measured field through a sequential colour map."
        ),
        options=(
            (
                "ipf",
                "IPF (inverse pole figure)",
                "Crystal direction along the specimen axis chosen below, in the standard "
                "fundamental-sector colours.",
            ),
            (
                "grain",
                "Grain identity",
                "One arbitrary colour per grain: the map to read a segmentation by.",
            ),
            (
                "grod",
                "GROD (grain reference orientation deviation)",
                "Per pixel, the misorientation from its own grain's reference orientation — "
                "intragranular gradients and stored plastic strain.",
            ),
            (
                "kam",
                "KAM (kernel average misorientation)",
                "Per pixel, the mean misorientation to its neighbours — local deformation, and "
                "the geometrically necessary dislocation content behind it.",
            ),
            (
                "confidence_index",
                "Confidence index",
                "The measured indexing confidence channel.",
            ),
            ("fit", "Fit", "The measured indexing fit, in degrees. Lower is better."),
            ("image_quality", "Image quality", "The measured pattern-quality channel."),
        ),
        default="ipf",
    )


def _imported_map(payload: Any) -> tuple[Any, Any]:
    """Read a user's own scan file into the same objects the gallery produces.

    The file goes through the library's own importer — `read_scan`, the same
    call a script would make, dispatching to `read_ang`, `read_ctf` or
    `read_oh5` on the extension — so phases, symmetry, grid topology and the
    quality channels all come from the file's header rather than from anything
    assumed here. What comes back is a `CrystalMap` and an entry describing it,
    which is exactly the pair a practice dataset produces, so every calculation
    below this point cannot tell the two apart.

    The "known answer" of an imported map is the honest one: there isn't one.
    Saying so is the point — the practice datasets are checkable because they
    were constructed, and a measurement is not, so a result built from a file
    must not be presented as if it carried a construction's guarantee.
    """

    from pytex.adapters.scan_files import scan_reader_for
    from pytex.app.ebsd_gallery import GalleryEntry
    from pytex.app.uploads import uploaded_file

    with uploaded_file(payload, field="scan_file", suffixes=SCAN_FILE_SUFFIXES) as (path, name):
        reader = scan_reader_for(path)
        APP_LOG.info(
            f"Reading the scan file {name}.",
            source="ebsd.map",
            detail={"file": name, "reader": reader.__name__},
        )
        try:
            scan = reader(path)
        except ImportError as error:
            raise InvalidInputError(
                f"{name} is an HDF5 scan, and this PyTex was installed without the library that "
                "reads HDF5.",
                details={"field": "scan_file"},
                hint=(
                    "Install PyTex with the 'hdf5' extra (`pip install pytex[hdf5]`), or open "
                    "the `.ang` export of the same scan instead."
                ),
            ) from error
        except (ValueError, KeyError, IndexError, OSError) as error:
            raise InvalidInputError(
                f"{name} could not be read as a {path.suffix.lower()} scan: {error}",
                details={"field": "scan_file"},
                hint=(
                    "The header is what fails first: a phase without a symmetry declaration, a "
                    "column layout the format does not describe, or — for an HDF5 scan — a file "
                    "holding patterns but no indexed orientations."
                ),
            ) from error

    crystal_map = scan.dataset.crystal_map
    channels = ", ".join(sorted(crystal_map.property_names)) or "none"
    phases = ", ".join(sorted(_phase_names(crystal_map))) or "unnamed"
    entry = GalleryEntry(
        id=f"file:{name}",
        title=name,
        summary=(
            f"{len(crystal_map.orientations)} measurements read from {name}; "
            f"phases: {phases}; scalar channels: {channels}."
        ),
        teaches=(
            "This is your own measurement, so nothing here is guaranteed in advance. Read the "
            "confidence-index or fit channel alongside the orientation map before believing a "
            "feature in it: a low-confidence region is where the indexing, not the "
            "microstructure, is producing the pattern you are looking at."
        ),
        known_answer=(
            "None — this is a measurement, not a construction. The practice datasets carry an "
            "answer fixed before the calculation ran; a file cannot, and a number from one is "
            "only as good as the scan behind it."
        ),
    )
    APP_LOG.success(
        f"Read {len(crystal_map.orientations)} measurements from {name}.",
        source="ebsd.map",
        detail={"file": name, "points": len(crystal_map.orientations), "channels": channels},
    )
    return crystal_map, entry


def _phase_names(crystal_map: Any) -> list[str]:
    """Names of the phases a crystal map declares, for the import summary."""

    candidates = list(getattr(crystal_map, "phases", None) or ())
    single = getattr(crystal_map, "phase", None)
    if single is not None:
        candidates.append(single)
    names: list[str] = []
    for phase in candidates:
        name = getattr(phase, "name", None)
        if name and str(name) not in names:
            names.append(str(name))
    return names


def _request_map(request: dict[str, Any]) -> tuple[Any, Any]:
    """Return the crystal map named by the request, and its gallery entry.

    A user's own scan file, when one is open, in preference to the practice
    dataset: someone who has opened their data is looking at their data, and
    silently analysing the example beside it would be the worst possible answer.
    """

    scan_file = request.get("scan_file")
    if scan_file:
        return _imported_map(scan_file)

    entry_id = str(request["dataset"])
    try:
        entry = get_entry(entry_id)
    except KeyError as error:
        raise InvalidInputError(
            f"Unknown dataset {entry_id!r}.",
            details={"field": "dataset"},
            hint=f"Choose one of: {', '.join(entry_ids())}.",
        ) from error
    APP_LOG.info(
        f"Building the '{entry.title}' practice map.", source="ebsd.map", detail={"id": entry_id}
    )
    return build_map(entry_id, grid=int(request["grid_points"])), entry


def _colour_lookup(name: str) -> np.ndarray:
    """A 256-entry RGB lookup table, as a ``(256, 3)`` uint8 array.

    ``viridis`` is the default for scalar fields because it is perceptually
    uniform and survives greyscale printing and colour-vision deficiency; a
    rainbow map would invent boundaries in a smooth GROD field where none exist.
    """

    positions = np.linspace(0.0, 1.0, 256)
    if name == "viridis":
        # Four anchors of the viridis ramp, linearly interpolated. Enough to
        # reproduce its perceptual behaviour without depending on matplotlib,
        # which this application deliberately does not import at request time.
        anchors = np.array(
            [
                [0.267, 0.005, 0.329],
                [0.229, 0.322, 0.545],
                [0.127, 0.567, 0.551],
                [0.369, 0.789, 0.383],
                [0.993, 0.906, 0.144],
            ]
        )
    elif name == "magma":
        anchors = np.array(
            [
                [0.001, 0.000, 0.014],
                [0.316, 0.072, 0.485],
                [0.716, 0.215, 0.475],
                [0.983, 0.529, 0.380],
                [0.987, 0.991, 0.750],
            ]
        )
    else:
        anchors = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    anchor_positions = np.linspace(0.0, 1.0, len(anchors))
    table = np.stack(
        [np.interp(positions, anchor_positions, anchors[:, channel]) for channel in range(3)],
        axis=1,
    )
    return np.clip(np.rint(table * 255.0), 0, 255).astype(np.uint8)


def _scalar_to_rgb(values: np.ndarray, *, colour_map: str) -> tuple[np.ndarray, float, float]:
    """Map a scalar field onto RGB, returning the range the colours stand for.

    The range is returned rather than kept private because a colour bar without
    numbers on it is decoration: a GROD map is unreadable unless the reader
    knows whether its brightest pixel is 2 degrees or 20.
    """

    finite = values[np.isfinite(values)]
    low = float(finite.min()) if finite.size else 0.0
    high = float(finite.max()) if finite.size else 1.0
    span = high - low
    normalised = np.zeros_like(values) if span <= 0.0 else (values - low) / span
    index = np.clip(np.rint(np.nan_to_num(normalised) * 255.0), 0, 255).astype(np.int64)
    return _colour_lookup(colour_map)[index], low, high


def _grain_colours(labels: np.ndarray) -> np.ndarray:
    """A distinct, stable colour per grain label.

    The hue is taken from the label through a golden-ratio step, which spreads
    consecutive labels as far apart in hue as possible — neighbouring grains
    usually carry consecutive labels, and a sequential palette would give them
    near-identical colours and hide the very boundary the map exists to show.
    """

    unique = np.unique(labels)
    hues = (np.arange(unique.size) * 0.61803398875) % 1.0
    # HSV to RGB at full saturation and value, written out because this is the
    # only place the application needs it.
    sector = np.floor(hues * 6.0).astype(np.int64) % 6
    fraction = hues * 6.0 - np.floor(hues * 6.0)
    rising = fraction
    falling = 1.0 - fraction
    zeros = np.zeros_like(hues)
    ones = np.ones_like(hues)
    table = np.stack(
        [
            np.choose(sector, [ones, falling, zeros, zeros, rising, ones]),
            np.choose(sector, [rising, ones, ones, falling, zeros, zeros]),
            np.choose(sector, [zeros, zeros, rising, ones, ones, falling]),
        ],
        axis=1,
    )
    palette = np.clip(np.rint(table * 235.0 + 10.0), 0, 255).astype(np.uint8)
    remap = np.zeros(int(unique.max()) + 1, dtype=np.int64)
    remap[unique] = np.arange(unique.size)
    coloured: np.ndarray = palette[remap[labels]]
    return coloured


def _boundary_lines(network: Any, crystal_map: Any) -> list[dict[str, Any]]:
    """Turn boundary segments into drawable line segments in map coordinates.

    A segment is a *pixel face*: the shared edge of two neighbouring points. Its
    endpoints are not stored, because they are implied — the face is
    perpendicular to the line joining the two point centres, centred on the
    stored midpoint, and one step long. Reconstructing it here keeps the drawn
    boundary exactly on the face rather than on a smoothed contour through the
    midpoints, which would round off the corners where three grains meet.
    """

    coordinates = np.asarray(crystal_map.coordinates, dtype=float)
    steps = crystal_map.step_sizes or (1.0, 1.0)
    step_x = float(steps[0])
    step_y = float(steps[-1])
    lines: list[dict[str, Any]] = []
    for segment in network.segments:
        left = coordinates[segment.left_index]
        right = coordinates[segment.right_index]
        joining = right[:2] - left[:2]
        norm = float(np.linalg.norm(joining))
        if norm <= 0.0:
            continue
        # Perpendicular to the join, in the plane of the map.
        normal = np.array([-joining[1], joining[0]]) / norm
        # The face runs along the normal, so its length is the step of *that*
        # axis: a face between two points separated in x is one y-step long.
        # Equal on a square grid, and not on a rectangular one.
        half = 0.5 * (step_y if abs(normal[1]) >= abs(normal[0]) else step_x)
        midpoint = np.asarray(segment.midpoint, dtype=float)[:2]
        start = midpoint - normal * half
        end = midpoint + normal * half
        lines.append(
            {
                "x1": float(start[0]),
                "y1": float(start[1]),
                "x2": float(end[0]),
                "y2": float(end[1]),
                "misorientation_deg": float(segment.misorientation_deg),
            }
        )
    return lines


#: Colour of a raster cell no measurement lands in. Only a staggered scan
#: produces any: the raster is on a square pitch and a hexagonal scan does not
#: fill it. Mid grey rather than black or white, so it is legible as "no
#: measurement here" against both a dark IPF colour and a pale one.
_EMPTY_CELL = np.array([44, 48, 56], dtype=np.uint8)


class _Raster:
    """Where each measurement point lands in the image the panel draws.

    A square scan is already a raster: the points are in row-major order and the
    grid shape says so. A **hexagonal** scan — the EDAX default, and therefore
    most of the ``.ang`` files anyone opens — is not. Its rows are offset by half
    a step and hold alternating counts, so it has no rectangular shape at all and
    the reader reports none.

    Rather than refuse those files, the points are placed into a square raster of
    half-step pitch, which is the pitch the stagger actually lives on: every
    measurement lands on its own cell, alternate cells between them stay empty,
    and the drawn map has the offset rows a hexagonal scan really has. Nothing is
    interpolated and no measurement moves — the empty cells are drawn as empty,
    because inventing a value there would be inventing data.

    Attributes
    ----------
    rows, cols : int
        Raster size.
    placement : np.ndarray or None
        For each measurement point, its flat index in the raster. ``None`` when
        the points already fill the raster in row-major order.
    step_x, step_y : float
        Physical size of one raster cell, which is the scan step for a square
        grid and half of it across a staggered one.
    """

    def __init__(self, crystal_map: Any) -> None:
        steps = crystal_map.step_sizes or (1.0, 1.0)
        shape = crystal_map.grid_shape
        if shape:
            self.rows, self.cols = (int(shape[0]), int(shape[1]))
            self.placement: np.ndarray | None = None
            self.step_x = float(steps[0])
            self.step_y = float(steps[-1])
            return

        coordinates = np.asarray(crystal_map.coordinates, dtype=float)[:, :2]
        if coordinates.size == 0:
            raise InvalidInputError(
                "This scan holds no measurement points.",
                details={"field": "scan_file"},
                hint="Check that the file is a complete scan and not only a header.",
            )
        self.step_x = float(steps[0]) / 2.0
        self.step_y = float(steps[-1])
        origin = coordinates.min(axis=0)
        columns = np.rint((coordinates[:, 0] - origin[0]) / self.step_x).astype(np.int64)
        row_index = np.rint((coordinates[:, 1] - origin[1]) / self.step_y).astype(np.int64)
        self.cols = int(columns.max()) + 1
        self.rows = int(row_index.max()) + 1
        self.placement = row_index * self.cols + columns

    @property
    def size(self) -> int:
        return self.rows * self.cols

    def extent_um(self) -> list[float]:
        """The map's physical extent, as ``[x0, y0, x1, y1]``."""

        return [0.0, 0.0, (self.cols - 1) * self.step_x, (self.rows - 1) * self.step_y]


def _encode_rgb(rgb: np.ndarray, raster: _Raster) -> dict[str, Any]:
    """Base64 the RGB raster, row-major from the top."""

    values = rgb.astype(np.uint8)
    if raster.placement is None:
        grid = values.reshape(raster.rows, raster.cols, 3)
    else:
        grid = np.repeat(_EMPTY_CELL[None, :], raster.size, axis=0)
        grid[raster.placement] = values
        grid = grid.reshape(raster.rows, raster.cols, 3)
    return {
        "width": raster.cols,
        "height": raster.rows,
        "encoding": "base64-rgb8",
        "data": base64.b64encode(np.ascontiguousarray(grid).tobytes()).decode("ascii"),
    }


def _encode_grain_ids(labels: np.ndarray, raster: _Raster) -> dict[str, Any]:
    """Base64 the per-pixel grain labels onto the same raster as the image.

    Purpose
    -------
    A click on the map has to answer *which grain is that*, and the only
    honest answer is the label the segmentation gave that measurement point.
    Sending the labels beside the pixels lets the browser answer it by lookup
    rather than by guessing from colour — two grains of one orientation share a
    colour exactly, and a colour-matching pick would join them silently.

    The array is row-major from the top, one ``int32`` per raster cell, aligned
    cell for cell with :func:`_encode_rgb`. Cells with no measurement — the gaps
    a hexagonal scan leaves in a square raster — carry ``-1``, the same value
    :attr:`GrainSegmentation.label_grid` uses for them, so an empty cell is
    picked as nothing rather than as grain zero.
    """

    values = np.asarray(labels, dtype=np.int32).reshape(-1)
    if raster.placement is None:
        grid = values.copy()
    else:
        grid = np.full(raster.size, -1, dtype=np.int32)
        grid[raster.placement] = values
    return {
        "width": raster.cols,
        "height": raster.rows,
        "encoding": "base64-int32",
        "data": base64.b64encode(
            np.ascontiguousarray(grid.reshape(raster.rows, raster.cols)).tobytes()
        ).decode("ascii"),
    }


def _builtin_phase_ids(grain_rows: list[dict[str, Any]]) -> dict[str, str | None]:
    """Which built-in phase, if any, each phase the scan names corresponds to.

    Purpose
    -------
    A grain handed to a relationship calculation needs a *phase* — a lattice and
    a point group — and a scan carries only a name. Where that name is one this
    application already knows, the correspondence can be made once here and
    checked by a test; where it is not, the answer is ``None`` and the surface
    that uses it must ask the user rather than choose.

    The names are taken from the grain rows rather than from the map header,
    so what is offered is exactly what the pickable grains are labelled with.

    Matching is by name, case- and punctuation-insensitively, against
    :data:`~pytex.app.phases.BUILTIN_PHASES`. Deliberately nothing cleverer: a
    lattice-parameter match would name a phase the scan never claimed, and a
    fuzzy name match would confuse alpha and beta zirconium, which differ in
    symmetry and therefore in every number computed from them.
    """

    from pytex.app.phases import BUILTIN_PHASES

    def key(name: str) -> str:
        return "".join(character for character in name.lower() if character.isalnum())

    catalogue = {key(spec.name): identifier for identifier, spec in BUILTIN_PHASES.items()}
    names = {str(row["phase_name"]) for row in grain_rows if row.get("phase_name")}
    return {name: catalogue.get(key(name)) for name in sorted(names)}


def _source_parameters() -> tuple[Parameter, ...]:
    """Which scan an operation works on: a practice dataset, or the user's file.

    Declared once and shared by every EBSD operation. They are the same question
    for all of them — *which map* — and a second copy of the wording would let
    the map and the summary beside it describe their inputs differently while
    analysing the same scan.
    """

    return (
        ChoiceParameter(
            name="dataset",
            label="Dataset",
            help_text=(
                "Which practice map to analyse. Each is a construction with a known answer, so "
                "the numbers on screen can be checked rather than trusted. Ignored when a scan "
                "file of your own is open."
            ),
            options=tuple((entry.id, entry.title, entry.summary) for entry in GALLERY),
            default=DEFAULT_ENTRY_ID,
        ),
        ObjectParameter(
            name="scan_file",
            label="Your own scan",
            help_text=(
                "An EDAX/TSL `.ang`, Oxford/HKL `.ctf`, or EDAX OIM HDF5 `.oh5`/`.h5` file, "
                "opened with **Open a scan** above. When one is open it replaces the practice "
                "dataset, and every choice "
                "below — the colouring, the modulation, the boundaries, the grain threshold — "
                "means exactly what it means for a practice map.\n\n"
                "The file is read by the same library importer a script would call, so the "
                "phases, the symmetry and the scalar channels come from the file's own header "
                "rather than from anything assumed here."
            ),
            required=False,
        ),
        IntegerParameter(
            name="grid_points",
            label="Map size",
            help_text=(
                "Side of the square practice map, in measurement points. Most derived "
                "quantities cost linearly in the point count; the grain reference orientations "
                "cost the square of each grain's size, so a map of a few large grains is the "
                "slow case rather than a map of many small ones. 200 points a side draws in "
                "well under a second on the practice datasets. The known answers hold at every "
                "size — only the GROD magnitude scales with the map, since it is a deviation "
                "accumulated across it."
            ),
            units="points",
            default=56,
            minimum=16,
            maximum=200,
            advanced=True,
        ),
    )


@REGISTRY.operation(
    "ebsd.map",
    title="EBSD orientation map",
    summary=(
        "IPF, grain, GROD, KAM and scalar maps, greyed by any channel, with boundaries on top."
    ),
    help_text=(
        "Draws one orientation map from four independent choices, so combinations nobody would "
        "enumerate as separate buttons are still reachable.\n\n"
        "**Colour by** decides what a pixel's colour means. *IPF* is the standard orientation "
        "map: each pixel is coloured by which crystal direction lies along the specimen axis you "
        "choose, folded into the symmetry fundamental sector, so symmetrically equivalent "
        "orientations get the same colour. Two maps of different point groups are not "
        "colour-comparable. *Grain identity* colours each segmented grain arbitrarily — the map "
        "to check a segmentation by. *GROD* is the deviation of each pixel from its own grain's "
        "reference orientation, which is intragranular gradient and stored plastic strain; *KAM* "
        "is the mean misorientation to its neighbours, which is local deformation. The remaining "
        "choices render the measured channels.\n\n"
        "**Modulate by** darkens any of those colourings by a scalar channel. This is how an IPF "
        "map is made to show *where the indexing was poor* without giving up the orientation "
        "information: pixels of low confidence go dark while their hue stays.\n\n"
        "**Grain boundaries** superimposes the boundary network on whatever is underneath, "
        "classified into low- and high-angle at the threshold you set. Boundaries are drawn as "
        "line geometry rather than as pixels, so they stay sharp when the figure is zoomed.\n\n"
        "**Grain threshold** is what counts as one grain, and it changes everything computed "
        "below it — the grain table, GROD, and the boundary network. Conventionally 5 to 15 "
        "degrees; the value decides whether subgrains are resolved as separate grains.\n\n"
        "The KAM threshold is separate and does a different job: it excludes neighbour pairs "
        "above it from the average, which is the standard way to keep grain boundaries out of an "
        "intragranular KAM. Without it, boundary pixels report the boundary misorientation "
        "instead of the local gradient."
    ),
    parameters=(
        *_source_parameters(),
        _colouring_parameter(),
        ChoiceParameter(
            name="ipf_direction",
            label="IPF specimen direction",
            help_text=(
                "Which specimen axis the IPF colour refers to. Three maps of the same scan along "
                "X, Y and Z together fix the orientation; any one alone does not."
            ),
            options=(
                ("X", "X (rolling direction)", "The first specimen axis."),
                ("Y", "Y (transverse direction)", "The second specimen axis."),
                ("Z", "Z (normal direction)", "The map normal; the usual default."),
            ),
            default="Z",
            group="Colour",
        ),
        ChoiceParameter(
            name="modulate_by",
            label="Modulate by",
            help_text=(
                "Darken the colour by a measured scalar channel, keeping its hue. The standard "
                "use is an IPF map greyed by confidence index, so poorly indexed pixels recede "
                "without being removed."
            ),
            options=(
                ("none", "Nothing", "Full brightness everywhere."),
                ("confidence_index", "Confidence index", "Darken where indexing confidence fell."),
                ("fit", "Fit", "Darken where the fit worsened; the scale is inverted for this."),
                ("image_quality", "Image quality", "Darken where the pattern was weak."),
            ),
            default="none",
            group="Colour",
        ),
        NumberParameter(
            name="modulation_floor",
            label="Darkest modulation",
            help_text=(
                "How dark the worst pixel becomes, as a fraction of full brightness. Zero makes "
                "the worst pixels black, which hides them; 0.25 keeps them legible."
            ),
            default=0.25,
            minimum=0.0,
            maximum=1.0,
            group="Colour",
            advanced=True,
        ),
        ChoiceParameter(
            name="colour_map",
            label="Scalar colour map",
            help_text=(
                "Used for the GROD, KAM and measured-channel colourings. Both choices are "
                "perceptually uniform; a rainbow map is deliberately not offered, because it "
                "invents visible boundaries in a smooth field where none exist."
            ),
            options=(
                ("viridis", "Viridis", "Perceptually uniform, colour-vision safe."),
                ("magma", "Magma", "Perceptually uniform, dark to bright."),
            ),
            default="viridis",
            group="Colour",
            advanced=True,
        ),
        BooleanParameter(
            name="show_boundaries",
            label="Superimpose grain boundaries",
            help_text=(
                "Draw the boundary network over whichever map is beneath. Works on every "
                "colouring, which is the point: a GROD map with boundaries on it answers "
                "'is this gradient inside a grain or across one'."
            ),
            default=True,
            required=False,
            group="Boundaries",
        ),
        NumberParameter(
            name="grain_threshold_deg",
            label="Grain threshold",
            help_text=(
                "Neighbouring points join the same grain below this misorientation. "
                "Conventionally 5 to 15 degrees; it decides whether subgrains count as grains, "
                "and changes the grain table, GROD and the boundary network together."
            ),
            units="°",
            default=5.0,
            minimum=0.1,
            maximum=62.8,
            group="Boundaries",
            row="Thresholds",
        ),
        NumberParameter(
            name="high_angle_threshold_deg",
            label="High-angle threshold",
            help_text=(
                "The Read-Shockley dividing line between low-angle (dislocation-wall) and "
                "high-angle boundaries. A convention, not a physical constant."
            ),
            units="°",
            default=15.0,
            minimum=1.0,
            maximum=62.8,
            group="Boundaries",
            row="Thresholds",
        ),
        NumberParameter(
            name="kam_threshold_deg",
            label="KAM threshold",
            help_text=(
                "Neighbour pairs above this misorientation are excluded from the KAM average. "
                "This is what keeps grain boundaries out of an intragranular KAM; without it the "
                "boundary pixels report the boundary misorientation."
            ),
            units="°",
            default=5.0,
            minimum=0.1,
            maximum=62.8,
            group="Misorientation",
        ),
        IntegerParameter(
            name="kam_order",
            label="KAM neighbour shell",
            help_text="Higher orders average over a wider kernel, smoothing the local field.",
            default=1,
            minimum=1,
            maximum=4,
            group="Misorientation",
            advanced=True,
        ),
    ),
    returns=(
        "The map under `data.image` as base64 RGB at native resolution, the boundary network as "
        "line segments in map coordinates, a grain table, and the scalar range the colours stand "
        "for."
    ),
    panel="ebsd",
    documentation=DocumentationLink(
        "EBSD KAM and dislocation density", "theory/ebsd_kam_parameterization"
    ),
    tags=(
        "ebsd",
        "ipf",
        "grain boundaries",
        "grod",
        "kam",
        "orientation map",
        "confidence index",
        "microstructure",
    ),
    citations=(_CITATION_RANDLE_ENGLER, _CITATION_NOLZE_IPF, _CITATION_WRIGHT_KAM),
)
def _map(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.plotting.ipf import ipf_colors

    crystal_map, entry = _request_map(request)
    raster = _Raster(crystal_map)
    points = len(crystal_map.orientations)
    threshold = float(request["grain_threshold_deg"])
    high_angle = float(request["high_angle_threshold_deg"])
    colouring = str(request["colouring"])

    APP_LOG.info(
        f"Segmenting {points} points at a {threshold:g}° grain threshold.",
        source="ebsd.map",
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=threshold)
    APP_LOG.notice(
        f"{len(segmentation.grains)} grains found at {threshold:g}°.",
        source="ebsd.map",
        detail={"grains": len(segmentation.grains), "threshold_deg": threshold},
    )

    kam = np.asarray(
        crystal_map.kernel_average_misorientation_deg(
            threshold_deg=float(request["kam_threshold_deg"]),
            order=int(request["kam_order"]),
        ),
        dtype=float,
    ).reshape(points)
    grod = np.asarray(segmentation.grod_map_deg(), dtype=float).reshape(points)

    rgb, scale = _colour_field(request, crystal_map, segmentation, kam, grod, ipf_colors)
    rgb = _modulate(request, crystal_map, rgb)

    network = segmentation.boundary_network(high_angle_threshold_deg=high_angle)
    lines = _boundary_lines(network, crystal_map) if request.get("show_boundaries") else []
    boundary_rows = _boundary_rows(network, high_angle)
    grain_rows = _grain_rows(crystal_map, segmentation, kam)

    APP_LOG.notice(
        f"{network.count} boundary segments, {network.high_angle_count} of them high-angle above "
        f"{high_angle:g}°; mean misorientation {network.mean_misorientation_deg:.2f}°.",
        source="ebsd.map",
        detail={"segments": network.count, "high_angle": network.high_angle_count},
    )

    step = float((crystal_map.step_sizes or (1.0, 1.0))[0])
    result = AppResult(
        title=f"{_colouring_title(request)} of {entry.title}",
        summary=_summary(request, entry, crystal_map, segmentation, network, scale, kam, grod),
        table=ResultTable(
            columns=_GRAIN_COLUMNS,
            rows=tuple(grain_rows),
            caption=f"Grains at a {threshold:g}° threshold, largest first.",
        ),
        data={
            "image": _encode_rgb(rgb, raster),
            "extent_um": raster.extent_um(),
            "step_um": step,
            "grid_shape": [raster.rows, raster.cols],
            "colouring": colouring,
            "colour_scale": scale,
            "boundaries": lines,
            "boundary_summary": boundary_rows,
            "boundary_columns": [column.to_json() for column in _BOUNDARY_COLUMNS],
            "high_angle_threshold_deg": high_angle,
            "grains": grain_rows,
            "columns": [column.to_json() for column in _GRAIN_COLUMNS],
            "grain_count": len(segmentation.grains),
            # The labels beside the pixels, so a click on the map resolves to a
            # grain by lookup rather than by colour. See `_encode_grain_ids`.
            "grain_ids": _encode_grain_ids(segmentation.labels, raster),
            "phase_builtins": _builtin_phase_ids(grain_rows),
            "dataset": entry.describe(),
            "ipf_direction": str(request["ipf_direction"]),
            "modulate_by": str(request["modulate_by"]),
        },
        inputs={
            "dataset": entry.id,
            "colouring": colouring,
            "ipf_direction": str(request["ipf_direction"]),
            "modulate_by": str(request["modulate_by"]),
            "grain_threshold_deg": threshold,
            "high_angle_threshold_deg": high_angle,
            "kam_threshold_deg": float(request["kam_threshold_deg"]),
            "kam_order": int(request["kam_order"]),
            "grid_points": int(request["grid_points"]),
            "show_boundaries": bool(request.get("show_boundaries")),
        },
        notes=_notes(request, entry),
        citations=(_CITATION_RANDLE_ENGLER, _CITATION_NOLZE_IPF, _CITATION_WRIGHT_KAM),
    )
    return result.to_json()


_SUMMARY_COLUMNS = (
    Column("group", "Section"),
    Column("metric", "Quantity"),
    Column("value", "Value"),
    Column("note", "What it says"),
)


def _channel_statistics(values: np.ndarray) -> dict[str, float] | None:
    """Mean, spread and extremes of one measured channel, ignoring absent points."""

    finite = np.asarray(values, dtype=float).reshape(-1)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    return {
        "count": int(finite.size),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "std": float(np.std(finite)),
        "minimum": float(np.min(finite)),
        "maximum": float(np.max(finite)),
        "p05": float(np.percentile(finite, 5.0)),
        "p95": float(np.percentile(finite, 95.0)),
    }


@REGISTRY.operation(
    "ebsd.scan_summary",
    title="Scan summary",
    summary="What this scan is: points, grid, phases, indexing quality, microstructure.",
    help_text=(
        "The first thing to look at after a scan is opened, and the last thing anyone builds by "
        "hand. Four sections, in the order the questions are asked.\n\n"
        "**Acquisition** is the geometry of the measurement: how many points, on what grid, at "
        "what step, over what area. A step much larger than the microstructure means the grains "
        "are undersampled whatever else the map shows.\n\n"
        "**Indexing quality** is the confidence index, the fit and the image quality, each with "
        "its mean, median, spread and 5th and 95th percentiles. The mean alone hides the shape "
        "that matters: a scan with CI 0.6 everywhere and one with half its points at 0.9 and "
        "half at 0.3 have the same mean and are not the same scan. The *indexed fraction* above "
        "a stated confidence threshold is the number most often quoted, so it is computed here "
        "with the threshold visible rather than assumed.\n\n"
        "**Phases** are the point counts and area fractions of each phase the file declares, "
        "with its point group — the check that the file was read as the material it is.\n\n"
        "**Microstructure** is what segmentation makes of it: how many grains, how large they "
        "are, how much orientation spread they carry, and how much of the boundary length is "
        "high-angle. Each depends on the grain threshold, which is stated with them, because a "
        "grain count without the threshold that produced it is not a measurement.\n\n"
        "Nothing here is a substitute for looking at the map. It is the page of numbers that "
        "says whether the map is worth trusting."
    ),
    parameters=(
        *_source_parameters(),
        NumberParameter(
            name="confidence_threshold",
            label="Counts as indexed above",
            help_text=(
                "The confidence index above which a point is counted as indexed. 0.1 is the "
                "conventional cut-off for EDAX/TSL data; the number is stated with the answer "
                "because the answer changes with it."
            ),
            default=0.1,
            minimum=0.0,
            maximum=1.0,
        ),
        NumberParameter(
            name="grain_threshold_deg",
            label="Grain threshold",
            help_text=(
                "Neighbouring points join the same grain below this misorientation. Every "
                "microstructural number below it moves when this moves."
            ),
            units="deg",
            default=5.0,
            minimum=0.1,
            maximum=62.8,
            row="Thresholds",
        ),
        NumberParameter(
            name="high_angle_threshold_deg",
            label="High-angle threshold",
            help_text="Boundaries at or above this misorientation are counted as high-angle.",
            units="deg",
            default=15.0,
            minimum=0.1,
            maximum=62.8,
            row="Thresholds",
        ),
    ),
    returns="One row per reported quantity; the whole summary, sectioned, under `data`.",
    panel="ebsd",
    documentation=DocumentationLink(
        "EBSD foundation and crystal map model", "concepts/ebsd_foundation"
    ),
    citations=(_CITATION_RANDLE_ENGLER,),
    tags=("EBSD", "summary", "metadata", "statistics", "quality", "OIM"),
)
def _scan_summary(request: dict[str, Any]) -> dict[str, Any]:
    crystal_map, entry = _request_map(request)
    points = len(crystal_map.orientations)
    summary = crystal_map.summary()
    threshold = float(request["grain_threshold_deg"])
    high_angle = float(request["high_angle_threshold_deg"])
    confidence_threshold = float(request["confidence_threshold"])

    rows: list[dict[str, Any]] = []

    def add(group: str, metric: str, value: Any, note: str) -> None:
        rows.append({"group": group, "metric": metric, "value": value, "note": note})

    step = summary.get("step_sizes")
    grid_shape = summary.get("grid_shape")
    add("Acquisition", "Measurement points", f"{points:,}", "How many orientations were recorded.")
    add(
        "Acquisition",
        "Grid",
        (
            f"{summary['grid_kind']} {grid_shape[0]} x {grid_shape[1]}"
            if grid_shape
            else str(summary["grid_kind"])
        ),
        "The sampling topology. A hexagonal grid has six first neighbours, a square four.",
    )
    if step:
        add(
            "Acquisition",
            "Step size",
            " x ".join(f"{float(value):.4g}" for value in step) + " um",
            "The distance between measurements: the finest microstructure that can be resolved.",
        )
        if grid_shape and len(grid_shape) == 2 and len(step) >= 2:
            width = float(step[0]) * float(grid_shape[1])
            height = float(step[1]) * float(grid_shape[0])
            add(
                "Acquisition",
                "Scanned area",
                f"{width:.4g} x {height:.4g} um ({width * height:.4g} um2)",
                "How much material the statistics below are drawn from.",
            )
    add(
        "Acquisition",
        "Specimen frame",
        str(summary["specimen_frame"]),
        "The frame every orientation and every specimen direction is expressed in.",
    )

    channels: dict[str, dict[str, Any]] = {}
    for key, (channel, label, units) in _SCALAR_CHANNELS.items():
        if channel not in crystal_map.property_names:
            continue
        statistics = _channel_statistics(crystal_map.get_property(channel))
        if statistics is None:
            continue
        channels[key] = {"label": label, "units": units, **statistics}
        unit_text = f" {units}" if units else ""
        add(
            "Indexing quality",
            f"{label}, mean",
            f"{statistics['mean']:.4g}{unit_text}",
            f"Median {statistics['median']:.4g}, 5th to 95th percentile "
            f"{statistics['p05']:.4g} to {statistics['p95']:.4g}.",
        )

    indexed_fraction: float | None = None
    if "confidence_index" in crystal_map.property_names:
        confidence = np.asarray(crystal_map.get_property("confidence_index"), dtype=float)
        finite = confidence[np.isfinite(confidence)]
        if finite.size:
            indexed_fraction = float(np.mean(finite >= confidence_threshold))
            add(
                "Indexing quality",
                f"Indexed above CI {confidence_threshold:g}",
                f"{100.0 * indexed_fraction:.2f} %",
                "The fraction of points a conventional cut-off would keep.",
            )

    phases = []
    for name, count in crystal_map.phase_summary().items():
        fraction = float(count) / float(points) if points else 0.0
        phases.append({"name": name, "points": int(count), "fraction": fraction})
        add(
            "Phases",
            name,
            f"{count:,} points ({100.0 * fraction:.2f} %)",
            "Point count and area fraction, as the file declares the phase.",
        )
    for phase_entry in crystal_map.resolved_phase_entries:
        add(
            "Phases",
            f"{phase_entry.name}, point group",
            str(phase_entry.point_group),
            "The symmetry every misorientation and every IPF colour is reduced by.",
        )

    segmentation = crystal_map.segment_grains(max_misorientation_deg=threshold)
    diameters = np.asarray(list(segmentation.grain_equivalent_diameters().values()), dtype=float)
    gos = np.asarray(segmentation.gos_map_deg(), dtype=float).reshape(-1)
    network = segmentation.boundary_network(high_angle_threshold_deg=high_angle)
    lengths = np.asarray([segment.length for segment in network.segments], dtype=float)
    misorientations = np.asarray(
        [segment.misorientation_deg for segment in network.segments], dtype=float
    )
    high_angle_length = float(lengths[misorientations >= high_angle].sum()) if lengths.size else 0.0
    total_length = float(lengths.sum()) if lengths.size else 0.0

    add(
        "Microstructure",
        "Grains",
        f"{len(segmentation.grains):,}",
        f"At a {threshold:g} deg grain threshold. Lower it and subgrains become grains.",
    )
    if diameters.size:
        add(
            "Microstructure",
            "Equivalent diameter, mean",
            f"{float(np.mean(diameters)):.4g} um",
            f"Median {float(np.median(diameters)):.4g} um, largest "
            f"{float(np.max(diameters)):.4g} um. The circle of equal area.",
        )
    finite_gos = gos[np.isfinite(gos)]
    if finite_gos.size:
        add(
            "Microstructure",
            "Grain orientation spread, mean",
            f"{float(np.mean(finite_gos)):.4g} deg",
            "Mean deviation of a grain's points from its own reference orientation: how "
            "deformed the grains are.",
        )
    if total_length > 0.0:
        add(
            "Microstructure",
            "Boundary length",
            f"{total_length:.4g} um",
            f"{100.0 * high_angle_length / total_length:.2f} % of it at or above "
            f"{high_angle:g} deg.",
        )

    result = AppResult(
        title=f"Scan summary: {entry.title}",
        summary=(
            f"{points:,} points on a {summary['grid_kind']} grid"
            + (f" of {grid_shape[0]} x {grid_shape[1]}" if grid_shape else "")
            + f", {len(crystal_map.phase_summary())} phase(s), "
            f"{len(segmentation.grains):,} grains at a {threshold:g} deg threshold"
            + (
                f", {100.0 * indexed_fraction:.1f} % indexed above CI {confidence_threshold:g}."
                if indexed_fraction is not None
                else "."
            )
        ),
        table=ResultTable(
            columns=_SUMMARY_COLUMNS,
            rows=tuple(rows),
            caption=f"Summary of {entry.title}, at a {threshold:g} deg grain threshold.",
        ),
        data={
            "acquisition": {
                "point_count": points,
                "grid_kind": summary["grid_kind"],
                # Lists rather than the tuples the library reports: this payload
                # is a JSON contract, and a caller comparing it with parsed JSON
                # must get the same answer either side of the wire.
                "grid_shape": None if grid_shape is None else [int(v) for v in grid_shape],
                "step_sizes_um": None if step is None else [float(v) for v in step],
                "map_frame": summary["map_frame"],
                "specimen_frame": summary["specimen_frame"],
                "is_multiphase": bool(summary["is_multiphase"]),
            },
            "channels": channels,
            "indexed_fraction": indexed_fraction,
            "confidence_threshold": confidence_threshold,
            "phases": phases,
            "microstructure": {
                "grain_count": len(segmentation.grains),
                "grain_threshold_deg": threshold,
                "high_angle_threshold_deg": high_angle,
                "mean_equivalent_diameter_um": (
                    float(np.mean(diameters)) if diameters.size else None
                ),
                "median_equivalent_diameter_um": (
                    float(np.median(diameters)) if diameters.size else None
                ),
                "mean_gos_deg": float(np.mean(finite_gos)) if finite_gos.size else None,
                "boundary_length_um": total_length,
                "high_angle_length_fraction": (
                    high_angle_length / total_length if total_length > 0.0 else None
                ),
            },
            "dataset": entry.describe(),
            "describe": crystal_map.describe(),
        },
        inputs={
            "dataset": entry.id,
            "grain_threshold_deg": threshold,
            "high_angle_threshold_deg": high_angle,
            "confidence_threshold": confidence_threshold,
            "grid_points": int(request["grid_points"]),
        },
        notes=(
            "Every microstructural number depends on the grain threshold, which is reported "
            "with them; a grain count without its threshold is not a measurement.",
            "The mean of a quality channel hides its shape, so the median and the 5th and 95th "
            "percentiles are given beside it.",
        ),
        citations=(_CITATION_RANDLE_ENGLER,),
    )
    return result.to_json()


def _histogram(
    values: np.ndarray, *, bins: int, weights: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Bin values, tolerating a quantity that has no spread at all.

    A coherent twin's boundaries are every one of them at 60 degrees, and a
    two-orientation map's random pairs likewise. That is a distribution — a very
    informative one — but numpy refuses a zero-width range. Widening it about the
    single value keeps the answer readable and honest: the bins are visible on
    the axis, and the whole population sits in the middle one.
    """

    lowest = float(np.min(values))
    highest = float(np.max(values))
    span = highest - lowest
    # The test is relative, not against zero. Sixty degrees computed two ways
    # differ in the last bit, which is a span of 1e-14 — not a distribution with
    # structure in it, but enough to get past `span > 0` and then produce
    # duplicate bin edges inside numpy.
    scale = max(abs(lowest), abs(highest), 1.0)
    if span > scale * 1e-9:
        return np.histogram(values, bins=bins, weights=weights)
    # One bin, not `bins` of them. Splitting a single value across the middle
    # two bins of a padded range would report half the boundary length at 59.998
    # degrees and half at 60.002, which is an artefact of the padding rather than
    # anything measured. One bin says what is true: all of it, at one value.
    centre = 0.5 * (lowest + highest)
    pad = max(abs(centre) * 1e-3, 1e-6)
    return np.histogram(values, bins=1, weights=weights, range=(centre - pad, centre + pad))


_DISTRIBUTION_QUANTITIES = {
    "grain_diameter": (
        "Grain equivalent diameter",
        "um",
        "The diameter of the circle with each grain's area. One entry per grain, so it is a "
        "number-weighted distribution: the many small grains dominate it, which is what makes "
        "it different from the area-weighted one a micrograph suggests.",
    ),
    "grain_area": (
        "Grain area",
        "um2",
        "One entry per grain. Its long tail is the usual reason a grain-size distribution is "
        "plotted against the logarithm of size rather than against size.",
    ),
    "grain_aspect_ratio": (
        "Grain aspect ratio",
        "",
        "Major over minor axis of the best-fit ellipse, one per grain: how elongated the "
        "microstructure is, which a diameter alone cannot say.",
    ),
    "misorientation_angle": (
        "Boundary misorientation angle",
        "deg",
        "One entry per boundary segment, weighted by its length, which is what makes it a "
        "distribution of *boundary* rather than of *pixels*. A peak near 60 degrees in a cubic "
        "material is the signature of annealing twins.",
    ),
    "kam": (
        "Kernel average misorientation",
        "deg",
        "Per point: the mean misorientation to its neighbours, excluding pairs above the KAM "
        "threshold so grain boundaries do not enter it. Local plastic gradient.",
    ),
    "grod": (
        "Grain reference orientation deviation",
        "deg",
        "Per point: the deviation from its own grain's reference orientation. Intragranular "
        "gradient, and therefore stored deformation, at a scale KAM cannot see.",
    ),
    "confidence_index": (
        "Confidence index",
        "",
        "Per point, as measured. Its shape says whether a scan is uniformly good or two "
        "populations averaged together.",
    ),
    "fit": (
        "Fit",
        "deg",
        "Per point, as measured: how well the indexed solution matched the pattern.",
    ),
    "image_quality": (
        "Image quality",
        "",
        "Per point, as measured. It falls with deformation and with surface damage alike, so it "
        "is read beside the orientation data rather than alone.",
    ),
}


@REGISTRY.operation(
    "ebsd.distribution",
    title="Distribution of a scan quantity",
    summary="Histogram of grain size, misorientation angle, KAM, GROD or a measured channel.",
    help_text=(
        "A map answers *where*; a distribution answers *how much of what*, and the two are read "
        "together. This histograms one quantity of the open scan and reports the statistics of "
        "it beside the bins.\n\n"
        "**Grain quantities** — equivalent diameter, area, aspect ratio — have one entry per "
        "grain, so they are *number-weighted*: a microstructure of a few large grains and many "
        "small ones looks small-grained here and large-grained in a micrograph, and neither is "
        "wrong. The grain threshold decides what a grain is and therefore what is being "
        "counted.\n\n"
        "**The misorientation-angle distribution** is per boundary segment, weighted by segment "
        "length. Its shape is diagnostic: a peak at 60 degrees in a cubic material is annealing "
        "twins, a concentration below the high-angle threshold is substructure, and a "
        "featureless spread matching the random distribution is a texture-free material. The "
        "random reference is computed too, from randomly paired points of this scan rather than "
        "from a formula, so it carries this material's own symmetry.\n\n"
        "**Point quantities** — KAM, GROD, and the measured channels — have one entry per "
        "measurement point, and are the same numbers the corresponding map is coloured by. "
        "Reading the map and this histogram together is how a colour scale is judged: a map "
        "whose range is set by a handful of outlying pixels shows a nearly uniform field, and "
        "the histogram says so at a glance."
    ),
    parameters=(
        *_source_parameters(),
        ChoiceParameter(
            name="quantity",
            label="Distribution of",
            help_text="Which quantity to histogram.",
            options=tuple(
                (key, label if not units else f"{label} ({units})", note)
                for key, (label, units, note) in _DISTRIBUTION_QUANTITIES.items()
            ),
            default="grain_diameter",
        ),
        IntegerParameter(
            name="bins",
            label="Bins",
            help_text=(
                "How many bins the range is divided into. Too few hides structure; too many "
                "turns a distribution into its own noise."
            ),
            default=24,
            minimum=4,
            maximum=120,
        ),
        NumberParameter(
            name="grain_threshold_deg",
            label="Grain threshold",
            help_text="What counts as one grain, for the grain and GROD quantities.",
            units="deg",
            default=5.0,
            minimum=0.1,
            maximum=62.8,
        ),
        NumberParameter(
            name="kam_threshold_deg",
            label="KAM threshold",
            help_text=(
                "Neighbour pairs above this are excluded from KAM, which is how grain "
                "boundaries are kept out of an intragranular measurement."
            ),
            units="deg",
            default=5.0,
            minimum=0.1,
            maximum=62.8,
            advanced=True,
        ),
        IntegerParameter(
            name="kam_order",
            label="KAM neighbour order",
            help_text="How many shells of neighbours enter the kernel average.",
            default=1,
            minimum=1,
            maximum=4,
            advanced=True,
        ),
    ),
    returns="One row per bin; the values' statistics and any reference series under `data`.",
    panel="ebsd",
    citations=(_CITATION_RANDLE_ENGLER, _CITATION_WRIGHT_KAM),
    documentation=DocumentationLink(
        "EBSD Kernel Average Misorientation and Disorientation",
        "theory/ebsd_kam_parameterization",
    ),
    tags=("EBSD", "distribution", "histogram", "grain size", "misorientation", "KAM", "GROD"),
)
def _distribution(request: dict[str, Any]) -> dict[str, Any]:
    crystal_map, entry = _request_map(request)
    quantity = str(request["quantity"])
    label, units, note = _DISTRIBUTION_QUANTITIES[quantity]
    threshold = float(request["grain_threshold_deg"])
    bins = int(request["bins"])

    weights: np.ndarray | None = None
    reference: dict[str, Any] | None = None

    if quantity in {"grain_diameter", "grain_area", "grain_aspect_ratio"}:
        segmentation = crystal_map.segment_grains(max_misorientation_deg=threshold)
        source = {
            "grain_diameter": segmentation.grain_equivalent_diameters,
            "grain_area": segmentation.grain_areas,
            "grain_aspect_ratio": segmentation.grain_aspect_ratios,
        }[quantity]()
        values = np.asarray(list(source.values()), dtype=float)
        population = "grains"
    elif quantity == "misorientation_angle":
        segmentation = crystal_map.segment_grains(max_misorientation_deg=threshold)
        network = segmentation.boundary_network()
        values = np.asarray(
            [segment.misorientation_deg for segment in network.segments], dtype=float
        )
        weights = np.asarray([segment.length for segment in network.segments], dtype=float)
        population = "boundary segments"
        reference = _random_pair_reference(crystal_map, bins=bins)
    elif quantity == "kam":
        values = np.asarray(
            crystal_map.kernel_average_misorientation_deg(
                threshold_deg=float(request["kam_threshold_deg"]),
                order=int(request["kam_order"]),
            ),
            dtype=float,
        ).reshape(-1)
        population = "points"
    elif quantity == "grod":
        segmentation = crystal_map.segment_grains(max_misorientation_deg=threshold)
        values = np.asarray(segmentation.grod_map_deg(), dtype=float).reshape(-1)
        population = "points"
    else:
        channel = _SCALAR_CHANNELS[quantity][0]
        if channel not in crystal_map.property_names:
            raise InvalidInputError(
                f"This scan carries no {label.lower()} channel.",
                details={"field": "quantity"},
                hint=(
                    "The channels a scan has come from its own file. Choose a quantity the "
                    "scan actually measured, or a computed one such as KAM or GROD."
                ),
            )
        values = np.asarray(crystal_map.get_property(channel), dtype=float).reshape(-1)
        population = "points"

    finite = np.isfinite(values)
    values = values[finite]
    if weights is not None:
        weights = weights[finite]
    if values.size == 0:
        raise InvalidInputError(
            f"No {label.lower()} values to histogram in this scan.",
            details={"field": "quantity"},
            hint="A segmentation that finds one grain leaves nothing to distribute.",
        )

    counts, edges = _histogram(values, bins=bins, weights=weights)
    total = float(counts.sum())
    centres = 0.5 * (edges[:-1] + edges[1:])
    cumulative = np.cumsum(counts)
    rows = [
        {
            "centre": float(centre),
            "lower": float(edges[index]),
            "upper": float(edges[index + 1]),
            "count": float(counts[index]),
            "fraction": float(counts[index] / total) if total > 0.0 else 0.0,
            "cumulative": float(cumulative[index] / total) if total > 0.0 else 0.0,
        }
        for index, centre in enumerate(centres)
    ]

    statistics = {
        "count": int(values.size),
        "mean": float(np.average(values, weights=weights)),
        "median": float(np.median(values)),
        "std": float(np.std(values)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "p10": float(np.percentile(values, 10.0)),
        "p90": float(np.percentile(values, 90.0)),
    }

    unit_text = f" {units}" if units else ""
    result = AppResult(
        title=f"{label} distribution: {entry.title}",
        summary=(
            f"{statistics['count']:,} {population}, mean {statistics['mean']:.4g}{unit_text}, "
            f"median {statistics['median']:.4g}{unit_text}, spanning "
            f"{statistics['minimum']:.4g} to {statistics['maximum']:.4g}{unit_text} in "
            f"{bins} bins."
            + (
                " The random reference is the same measurement on randomly paired points of "
                "this scan, so it carries this material's own symmetry rather than a formula's."
                if reference
                else ""
            )
        ),
        table=ResultTable(
            columns=(
                Column("centre", "Bin centre", units=units, numeric=True, digits=4),
                Column("lower", "From", units=units, numeric=True, digits=4),
                Column("upper", "To", units=units, numeric=True, digits=4),
                Column(
                    "count",
                    "Length" if weights is not None else "Count",
                    units="um" if weights is not None else "",
                    numeric=True,
                    digits=3,
                ),
                Column("fraction", "Fraction", numeric=True, digits=4),
                Column("cumulative", "Cumulative", numeric=True, digits=4),
            ),
            rows=tuple(rows),
            caption=f"{label} over {statistics['count']:,} {population} of {entry.title}.",
        ),
        data={
            "quantity": quantity,
            "label": label,
            "units": units,
            "population": population,
            "weighted_by_length": weights is not None,
            "edges": [float(value) for value in edges],
            "counts": [float(value) for value in counts],
            "statistics": statistics,
            "reference": reference,
            "grain_threshold_deg": threshold,
            "dataset": entry.describe(),
        },
        inputs={
            "dataset": entry.id,
            "quantity": quantity,
            "bins": bins,
            "grain_threshold_deg": threshold,
            "grid_points": int(request["grid_points"]),
        },
        notes=(note,),
        citations=(_CITATION_RANDLE_ENGLER, _CITATION_WRIGHT_KAM),
    )
    return result.to_json()


def _random_pair_reference(crystal_map: Any, *, bins: int, sample: int = 160) -> dict[str, Any]:
    """The misorientation distribution of randomly paired points of this scan.

    The reference a measured misorientation distribution is read against. Taken
    from the scan's own orientations rather than from an analytic random
    distribution, so it carries this material's symmetry and its texture: a
    strongly textured scan has a random-*pair* distribution unlike the
    texture-free one, and comparing against the texture-free curve would call
    the texture a boundary preference.
    """

    count = len(crystal_map.orientations)
    if count < 4:
        return {}
    generator = np.random.default_rng(0)
    size = min(sample, count // 2)
    if size < 2:
        return {}
    # Two disjoint random subsets, taken through the map's own point selection so
    # the phase and symmetry travel with them.
    chosen = generator.choice(count, size=2 * size, replace=False)
    left_mask = np.zeros(count, dtype=bool)
    right_mask = np.zeros(count, dtype=bool)
    left_mask[chosen[:size]] = True
    right_mask[chosen[size:]] = True
    left = crystal_map.select_points(left_mask).orientations
    right = crystal_map.select_points(right_mask).orientations
    angles = np.degrees(np.asarray(left.misorientation_angles_to(right), dtype=float)).reshape(-1)
    angles = angles[np.isfinite(angles) & (angles > 1e-9)]
    if angles.size == 0:
        return {}
    counts, edges = _histogram(angles, bins=bins)
    total = float(counts.sum())
    return {
        "label": "Randomly paired points",
        "edges": [float(value) for value in edges],
        "fractions": [float(value / total) if total else 0.0 for value in counts],
        "pair_count": int(angles.size),
        "mean": float(np.mean(angles)),
    }


@REGISTRY.operation(
    "ebsd.discrete_figure",
    title="Discrete pole and inverse pole figure",
    summary="Every measured orientation as a point, rather than as a contoured density.",
    help_text=(
        "The scatter a contoured pole figure is made from. A density map is an estimate with a "
        "kernel width in it; the discrete figure is the measurement, and the two answer "
        "different questions.\n\n"
        "**Where a discrete figure is the right one.** When the count is small — a few hundred "
        "grains rather than a million pixels — a contour is mostly kernel and the scatter is "
        "the honest picture. When looking for *structure* rather than intensity: variant "
        "clusters, a fibre that is a line of points rather than a smear, a handful of outliers "
        "a contour would smooth away entirely.\n\n"
        "**Pole figure** plots where a chosen crystal plane normal points in the specimen "
        "frame, one point per symmetry-related pole per measurement. **Inverse pole figure** "
        "plots which crystal direction lies along a chosen specimen axis, folded into the "
        "fundamental sector.\n\n"
        "**Subsampling is stated, not hidden.** A scan of a hundred thousand points cannot be "
        "drawn as a hundred thousand markers, so a random subset is taken with a fixed seed and "
        "its size is reported. It is a random subset of the *points*, so it is unbiased with "
        "respect to orientation; it is not unbiased with respect to grain size, because a large "
        "grain contributes more points, exactly as it does in the map."
    ),
    parameters=(
        *_source_parameters(),
        ChoiceParameter(
            name="kind",
            label="Figure",
            help_text="Which of the two projections to draw.",
            options=(
                (
                    "pole",
                    "Pole figure",
                    "Where a crystal plane normal points in the specimen frame.",
                ),
                (
                    "inverse",
                    "Inverse pole figure",
                    "Which crystal direction lies along a specimen axis.",
                ),
            ),
            default="pole",
        ),
        IndicesParameter(
            name="pole",
            label="Plane (hkl)",
            help_text=(
                "The plane whose normals are plotted, for a pole figure. Its whole symmetry "
                "family is drawn, as a measured pole figure inevitably contains."
            ),
            default=(1, 1, 1),
        ),
        ChoiceParameter(
            name="sample_direction",
            label="Specimen direction",
            help_text="Which specimen axis the inverse pole figure refers to.",
            options=(
                ("x", "X", "The first specimen axis."),
                ("y", "Y", "The second specimen axis."),
                ("z", "Z", "The map normal; the usual default."),
            ),
            default="z",
        ),
        IntegerParameter(
            name="max_points",
            label="Points drawn",
            help_text=(
                "How many measurement points to draw, at most. Beyond a few thousand markers a "
                "figure is a solid disc and says less than a contour would."
            ),
            default=1000,
            minimum=50,
            maximum=20000,
        ),
    ),
    returns="One row per drawn point; the projected coordinates and the figure's frame under "
    "`data`.",
    panel="ebsd",
    citations=(_CITATION_RANDLE_ENGLER, _CITATION_NOLZE_IPF),
    tags=("EBSD", "pole figure", "inverse pole figure", "discrete", "scatter", "texture"),
)
def _discrete_figure(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.sphere import project_directions

    crystal_map, entry = _request_map(request)
    kind = str(request["kind"])
    limit = int(request["max_points"])
    count = len(crystal_map.orientations)

    # A fixed seed, so the same scan gives the same figure twice: a scatter that
    # reshuffles on every redraw cannot be compared with the one beside it.
    generator = np.random.default_rng(0)
    if count > limit:
        mask = np.zeros(count, dtype=bool)
        mask[generator.choice(count, size=limit, replace=False)] = True
        view = crystal_map.select_points(mask)
        drawn_from = limit
    else:
        view = crystal_map
        drawn_from = count

    if kind == "pole":
        indices = tuple(int(value) for value in request["pole"])
        figure = view.pole_figure(indices)
        directions = np.asarray(figure.sample_directions, dtype=float)
        label = f"{{{' '.join(str(value) for value in indices)}}} pole figure"
        frame_text = "specimen frame: X right, Y up"
    else:
        direction = str(request["sample_direction"])
        figure = view.inverse_pole_figure(direction)
        directions = np.asarray(figure.crystal_directions, dtype=float)
        label = f"Inverse pole figure of specimen {direction.upper()}"
        frame_text = "crystal frame, folded into the fundamental sector"

    projected = np.asarray(
        project_directions(directions, method="stereographic", antipodal=True), dtype=float
    )
    points = [
        {"x": round(float(point[0]), 4), "y": round(float(point[1]), 4)} for point in projected
    ]

    rows = [
        {"point": index + 1, "x": entry_point["x"], "y": entry_point["y"]}
        for index, entry_point in enumerate(points[:200])
    ]
    result = AppResult(
        title=f"{label}: {entry.title}",
        summary=(
            f"{len(points):,} projected points from {drawn_from:,} measurement point(s)"
            + (f", subsampled from {count:,}" if count > limit else "")
            + f". Stereographic, upper hemisphere, {frame_text}."
            + (
                " Every symmetry-related member of the family is plotted, which is why there "
                "are more points than measurements."
                if kind == "pole"
                else ""
            )
        ),
        table=ResultTable(
            columns=(
                Column("point", "Point", numeric=True),
                Column("x", "x", numeric=True, digits=4),
                Column("y", "y", numeric=True, digits=4),
            ),
            rows=tuple(rows),
            caption=f"The first {len(rows):,} projected points of the {label.lower()}.",
        ),
        data={
            "kind": kind,
            "label": label,
            "points": points,
            "drawn_points": len(points),
            "measurement_points": int(drawn_from),
            "scan_points": count,
            "subsampled": bool(count > limit),
            "projection": "stereographic",
            "dataset": entry.describe(),
        },
        inputs={
            "dataset": entry.id,
            "kind": kind,
            "pole": [int(value) for value in request["pole"]],
            "sample_direction": str(request["sample_direction"]),
            "max_points": limit,
            "grid_points": int(request["grid_points"]),
        },
        notes=(
            "A discrete figure is the measurement; a contoured one is an estimate of a density "
            "from it, with a kernel width that must be stated. Neither replaces the other.",
            "The subset is random over measurement points and therefore unbiased in "
            "orientation, but not in grain size: a large grain contributes more points, as it "
            "does in the map.",
        ),
        citations=(_CITATION_RANDLE_ENGLER, _CITATION_NOLZE_IPF),
    )
    return result.to_json()


def _colouring_title(request: dict[str, Any]) -> str:
    titles = {
        "ipf": f"IPF-{request['ipf_direction']} map",
        "grain": "Grain map",
        "grod": "GROD map",
        "kam": "KAM map",
        "confidence_index": "Confidence-index map",
        "fit": "Fit map",
        "image_quality": "Image-quality map",
    }
    return titles.get(str(request["colouring"]), "Orientation map")


def _colour_field(
    request: dict[str, Any],
    crystal_map: Any,
    segmentation: Any,
    kam: np.ndarray,
    grod: np.ndarray,
    ipf_colors: Any,
) -> tuple[np.ndarray, dict[str, Any] | None]:
    """Colour every point, and report the range a scalar colouring stands for."""

    colouring = str(request["colouring"])
    colour_map = str(request["colour_map"])

    if colouring == "ipf":
        count = len(crystal_map.orientations)
        # Narrated because a wait nobody explained reads as a hang. The cost is
        # the symmetry reduction: every direction is folded into the
        # fundamental sector, which pytex.core.symmetry does for the whole map
        # in one blocked pass over the symmetry orbit.
        APP_LOG.info(
            f"Folding {count} orientations into the fundamental sector for IPF-"
            f"{request['ipf_direction']} colouring.",
            source="ebsd.map",
        )
        colours = np.asarray(
            ipf_colors(crystal_map.orientations, direction=str(request["ipf_direction"])),
            dtype=float,
        )
        return np.clip(np.rint(colours * 255.0), 0, 255).astype(np.uint8), None
    if colouring == "grain":
        labels = np.asarray(segmentation.labels, dtype=np.int64).reshape(-1)
        return _grain_colours(labels), None

    if colouring == "grod":
        values, label, units = grod, "GROD", "°"
    elif colouring == "kam":
        values, label, units = kam, "KAM", "°"
    else:
        channel, label, units = _SCALAR_CHANNELS[colouring]
        try:
            values = np.asarray(crystal_map.get_property(channel), dtype=float).reshape(-1)
        except KeyError as error:
            raise UnsupportedRequestError(
                f"This dataset carries no {label.lower()} channel.",
                details={"field": "colouring"},
                hint=f"Available channels: {', '.join(crystal_map.property_names) or 'none'}.",
            ) from error

    rgb, low, high = _scalar_to_rgb(values, colour_map=colour_map)
    return rgb, {
        "label": label,
        "units": units,
        "minimum": low,
        "maximum": high,
        "colour_map": colour_map,
    }


def _modulate(request: dict[str, Any], crystal_map: Any, rgb: np.ndarray) -> np.ndarray:
    """Darken the colours by a scalar channel, keeping their hue."""

    choice = str(request["modulate_by"])
    if choice == "none":
        return rgb
    channel, label, _units = _SCALAR_CHANNELS[choice]
    try:
        values = np.asarray(crystal_map.get_property(channel), dtype=float).reshape(-1)
    except KeyError as error:
        raise UnsupportedRequestError(
            f"This dataset carries no {label.lower()} channel to modulate by.",
            details={"field": "modulate_by"},
            hint=f"Available channels: {', '.join(crystal_map.property_names) or 'none'}.",
        ) from error

    finite = values[np.isfinite(values)]
    low = float(finite.min()) if finite.size else 0.0
    high = float(finite.max()) if finite.size else 1.0
    span = high - low
    quality = np.full_like(values, 1.0) if span <= 0.0 else (values - low) / span
    if choice == "fit":
        # Fit is an error: a large value is a *worse* measurement, so the scale
        # inverts. Modulating by it directly would darken exactly the pixels
        # that were indexed best.
        quality = 1.0 - quality
    floor = float(request["modulation_floor"])
    weight = floor + (1.0 - floor) * np.nan_to_num(quality, nan=0.0)
    APP_LOG.info(
        f"Modulating brightness by {label.lower()} over {low:.3g} to {high:.3g}.",
        source="ebsd.map",
    )
    return np.clip(np.rint(rgb.astype(float) * weight[:, None]), 0, 255).astype(np.uint8)


def _grain_orientation_rows(crystal_map: Any, segmentation: Any) -> dict[int, dict[str, Any]]:
    """Each grain's mean orientation as Bunge angles, and the phase it belongs to.

    Purpose
    -------
    Closing the loop from the map. The measured-pair operations take two
    orientations as Euler angles and two phases; without this a user reads six
    numbers off one screen and types them into another, which is exactly the
    kind of hand transcription this repository refuses everywhere else.

    What the numbers are, precisely
    -------------------------------
    The mean is the symmetry-aware average over the grain's member points
    (`GrainSegmentation.grain_mean_orientation`), **not** the reference point
    GROD is measured against — that one is a single measured orientation chosen
    as representative, and averaging is the right operation for handing a grain
    to a relationship calculation. The angles are Bunge (ZXZ), which is what
    every EBSD vendor exports and what the measured-pair operations default to.

    The grain-orientation spread travels with them in the same table row. That
    is deliberate: a mean has no scatter of its own, so a relationship computed
    from two means would report a residual of zero however noisy the two grains
    are, and the spread is the only honest measure of what that zero conceals.
    """

    from pytex.core.representations import quaternions_to_euler_angles

    entries = tuple(crystal_map.phase_entries or ())
    phase_ids = crystal_map.phase_ids
    rows: dict[int, dict[str, Any]] = {}
    for grain in segmentation.grains:
        mean = segmentation.grain_mean_orientation(grain)
        angles = quaternions_to_euler_angles(
            np.asarray(mean.rotation.quaternion, dtype=float)[None, :], degrees=True
        )[0]
        name: str | None = None
        if phase_ids is not None and entries:
            # A multiphase map carries a phase per point; a grain does not cross
            # a phase boundary, so its first member names the whole grain.
            members = np.asarray(grain.member_indices, dtype=np.int64)
            identifier = int(np.asarray(phase_ids, dtype=np.int64)[members][0])
            match = next(
                (entry for entry in entries if int(entry.phase_id) == identifier), None
            )
            name = None if match is None else str(match.name)
        elif entries:
            name = str(entries[0].name)
        elif crystal_map.orientations.phase is not None:
            name = str(crystal_map.orientations.phase.name)
        rows[int(grain.grain_id)] = {
            "mean_phi1_deg": float(angles[0]),
            "mean_Phi_deg": float(angles[1]),
            "mean_phi2_deg": float(angles[2]),
            "phase_name": name,
        }
    return rows


def _grain_rows(crystal_map: Any, segmentation: Any, kam: np.ndarray) -> list[dict[str, Any]]:
    """One row per grain, largest first."""

    spreads = segmentation.grain_orientation_spread_deg()
    orientations = _grain_orientation_rows(crystal_map, segmentation)
    step = float((crystal_map.step_sizes or (1.0, 1.0))[0])
    point_area = step * step
    rows = []
    for grain in segmentation.grains:
        members = np.asarray(grain.member_indices, dtype=np.int64)
        area = float(grain.size) * point_area
        rows.append(
            {
                "grain_id": int(grain.grain_id),
                "size": int(grain.size),
                "area_um2": area,
                "equivalent_diameter_um": float(2.0 * np.sqrt(area / np.pi)),
                "grain_orientation_spread_deg": float(spreads.get(int(grain.grain_id), 0.0)),
                "mean_kam_deg": float(kam[members].mean()) if members.size else 0.0,
                **orientations[int(grain.grain_id)],
            }
        )
    rows.sort(key=lambda row: row["size"], reverse=True)
    return rows


def _boundary_rows(network: Any, high_angle: float) -> list[dict[str, Any]]:
    """Low- and high-angle boundary totals, by length rather than by count.

    By length, because that is the quantity a boundary-character distribution is
    defined on: counting segments weights a boundary by how finely the grid
    happened to sample it.
    """

    if not network.segments:
        return []
    angles = np.array([segment.misorientation_deg for segment in network.segments], dtype=float)
    lengths = np.array([segment.length for segment in network.segments], dtype=float)
    total = float(lengths.sum())
    rows = []
    for character, mask in (
        ("Low-angle", angles < high_angle),
        ("High-angle", angles >= high_angle),
    ):
        selected = lengths[mask]
        rows.append(
            {
                "character": f"{character} (<{high_angle:g}°)"
                if character == "Low-angle"
                else f"{character} (≥{high_angle:g}°)",
                "count": int(mask.sum()),
                "length_um": float(selected.sum()),
                "fraction": float(selected.sum() / total) if total > 0.0 else 0.0,
                "mean_misorientation_deg": float(angles[mask].mean()) if mask.any() else 0.0,
            }
        )
    return rows


def _summary(
    request: dict[str, Any],
    entry: Any,
    crystal_map: Any,
    segmentation: Any,
    network: Any,
    scale: dict[str, Any] | None,
    kam: np.ndarray,
    grod: np.ndarray,
) -> str:
    raster = _Raster(crystal_map)
    step = float((crystal_map.step_sizes or (1.0, 1.0))[0])
    threshold = float(request["grain_threshold_deg"])
    extent = raster.extent_um()
    shape = (
        f"{raster.rows}×{raster.cols} points"
        if raster.placement is None
        else f"{len(crystal_map.orientations)} points on a staggered scan"
    )
    parts = [
        f"{shape} at a {step:g} µm step, covering "
        f"{extent[2]:g} × {extent[3]:g} µm. "
        f"{len(segmentation.grains)} grains at a {threshold:g}° threshold, with "
        f"{network.count} boundary segments whose mean misorientation is "
        f"{network.mean_misorientation_deg:.2f}°."
    ]
    if scale is not None:
        parts.append(
            f"The colour runs from {scale['minimum']:.4g} to {scale['maximum']:.4g}"
            f"{scale['units']} of {scale['label']}."
        )
    elif str(request["colouring"]) == "ipf":
        parts.append(
            f"Colour is IPF along {request['ipf_direction']}, folded into the cubic fundamental "
            "sector so symmetrically equivalent orientations share a colour."
        )
    else:
        parts.append("Colour is grain identity, one arbitrary hue per grain.")
    parts.append(
        f"KAM averages {kam.mean():.3f}° and reaches {kam.max():.3f}°; GROD reaches "
        f"{grod.max():.3f}°."
    )
    modulation = str(request["modulate_by"])
    if modulation != "none":
        parts.append(f"Brightness is modulated by {_SCALAR_CHANNELS[modulation][1].lower()}.")
    return " ".join(parts)


def _notes(request: dict[str, Any], entry: Any) -> tuple[str, ...]:
    notes = [
        f"Known answer for this dataset: {entry.known_answer}",
        "This is a constructed microstructure, not a measurement. It carries no detector "
        "geometry and no indexing step; its value is that the answer is known before the "
        "calculation runs.",
    ]
    if str(request["colouring"]) == "ipf":
        notes.append(
            "IPF colour is defined only up to the crystal symmetry, and the colour key belongs to "
            "the point group. Maps of two different point groups are not colour-comparable, and "
            "one direction alone does not fix an orientation — X, Y and Z together do."
        )
    if str(request["colouring"]) == "kam":
        notes.append(
            f"Neighbour pairs above {float(request['kam_threshold_deg']):g}° are excluded from the "
            "average, so grain boundaries do not leak into the intragranular field. Raising that "
            "threshold above the grain threshold makes boundary pixels report the boundary."
        )
    if str(request["colouring"]) == "grod":
        notes.append(
            "GROD is measured from each grain's own reference orientation, so it falls to zero at "
            "that point and rises away from it in every direction — it is a deviation, not a ramp."
        )
    if str(request["modulate_by"]) == "fit":
        notes.append(
            "Fit is an error, so the modulation scale is inverted: a large fit darkens the pixel."
        )
    return tuple(notes)


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="ebsd.example.scan_summary",
            title="What is in this scan",
            panel="ebsd",
            summary="Points, grid, phases, indexing quality and microstructure, in one page.",
            teaches=(
                "The page of numbers that decides whether a map is worth reading, and the "
                "habit of looking at it first. Three questions are answered together and are "
                "usually asked apart. Is the *step* small enough for the microstructure - a "
                "scan whose grains are three pixels across cannot support a grain-size "
                "distribution whatever the map looks like. Is the *indexing* good, and good "
                "uniformly - the mean confidence index hides the difference between a scan "
                "that is uniformly fair and one that is half excellent and half unindexed, "
                "which is why the median and the percentiles are beside it. And is the "
                "*segmentation* sensible - the grain count is quoted with the threshold that "
                "produced it, because without that it is not a measurement of anything."
            ),
            operation="ebsd.scan_summary",
            request={"dataset": "equiaxed_polycrystal", "grain_threshold_deg": 5.0},
        ),
        ExampleScenario(
            id="ebsd.example.grain_size_distribution",
            title="Grain-size distribution",
            panel="ebsd",
            summary="Equivalent circular diameter over every segmented grain.",
            teaches=(
                "That a grain-size distribution is a statement about a *segmentation*, not "
                "only about a material. Every grain here is one entry, so the distribution is "
                "number-weighted: many small grains dominate it, while a micrograph - where "
                "each grain occupies its own area - reads as though the large ones do. Neither "
                "weighting is wrong and they answer different questions, which is why the "
                "weighting is named on the axis. Change the grain threshold and the whole "
                "distribution moves, because subgrains become grains."
            ),
            operation="ebsd.distribution",
            request={
                "dataset": "equiaxed_polycrystal",
                "quantity": "grain_diameter",
                "bins": 20,
            },
        ),
        ExampleScenario(
            id="ebsd.example.misorientation_distribution",
            title="Misorientation angles, against random pairs",
            panel="ebsd",
            summary="The boundary misorientation distribution, with its random reference.",
            teaches=(
                "How a misorientation distribution is read, which is never on its own. The "
                "measured curve is per boundary segment and weighted by length, so it is a "
                "distribution of boundary rather than of pixels. Beside it is the same "
                "measurement made on randomly paired points of this same scan - the reference "
                "a departure is judged against, and deliberately computed from this material "
                "rather than from the texture-free formula, because a textured material has a "
                "random-pair distribution of its own and comparing against the formula would "
                "report its texture as a boundary preference."
            ),
            operation="ebsd.distribution",
            request={
                "dataset": "equiaxed_polycrystal",
                "quantity": "misorientation_angle",
                "bins": 24,
            },
        ),
        ExampleScenario(
            id="ebsd.example.discrete_pole_figure",
            title="Discrete pole figure of the scan",
            panel="ebsd",
            summary="Every measured orientation as a point, rather than as a contour.",
            teaches=(
                "The difference between a measurement and an estimate of its density. A "
                "contoured pole figure has a kernel width in it, chosen by whoever drew it; "
                "the scatter has none, and shows the structure a kernel smooths away - variant "
                "clusters, a fibre that is a line of points rather than a smear, and the "
                "handful of outliers that a contour removes entirely. The subset drawn is "
                "random over measurement points and its size is stated, because a figure that "
                "silently drew a tenth of the data would be a different figure."
            ),
            operation="ebsd.discrete_figure",
            request={
                "dataset": "equiaxed_polycrystal",
                "kind": "pole",
                "pole": [1, 1, 1],
                "max_points": 600,
            },
        ),
        ExampleScenario(
            id="ebsd.example.ipf_with_boundaries",
            title="IPF map with grain boundaries",
            panel="ebsd",
            summary="The standard orientation map of an equiaxed polycrystal, boundaries on top.",
            teaches=(
                "What an EBSD map is before anything is derived from it: colour is crystal "
                "direction along Z, and the boundary network drawn over it is where that "
                "direction changes by more than the grain threshold."
            ),
            operation="ebsd.map",
            request={
                "dataset": "equiaxed_polycrystal",
                "colouring": "ipf",
                "ipf_direction": "Z",
                "modulate_by": "none",
                "modulation_floor": 0.25,
                "colour_map": "viridis",
                "show_boundaries": True,
                "grain_threshold_deg": 5.0,
                "high_angle_threshold_deg": 15.0,
                "kam_threshold_deg": 5.0,
                "kam_order": 1,
                "grid_points": 56,
            },
        ),
        ExampleScenario(
            id="ebsd.example.ipf_greyed_by_confidence",
            title="IPF map greyed by confidence index",
            panel="ebsd",
            summary="The same map, darkened where the indexing confidence fell.",
            teaches=(
                "Why scalar modulation exists: the boundaries go dark because a pattern collected "
                "there overlaps two lattices, so the map shows *where it should be believed* "
                "without giving up the orientation it is showing."
            ),
            operation="ebsd.map",
            request={
                "dataset": "equiaxed_polycrystal",
                "colouring": "ipf",
                "ipf_direction": "Z",
                "modulate_by": "confidence_index",
                "modulation_floor": 0.1,
                "colour_map": "viridis",
                "show_boundaries": False,
                "grain_threshold_deg": 5.0,
                "high_angle_threshold_deg": 15.0,
                "kam_threshold_deg": 5.0,
                "kam_order": 1,
                "grid_points": 56,
            },
        ),
        ExampleScenario(
            id="ebsd.example.grod_gradient",
            title="GROD across a deformation gradient",
            panel="ebsd",
            summary="A bicrystal whose right-hand grain carries a linear orientation gradient.",
            teaches=(
                "That GROD is a deviation from the grain's own reference orientation, not a ramp "
                "across it: the field falls to zero at the reference point and rises on both "
                "sides. Switch the colouring to KAM on the same map and the gradient becomes "
                "uniform, because a linear gradient presents the same step everywhere."
            ),
            operation="ebsd.map",
            request={
                "dataset": "bicrystal_gradient",
                "colouring": "grod",
                "ipf_direction": "Z",
                "modulate_by": "none",
                "modulation_floor": 0.25,
                "colour_map": "viridis",
                "show_boundaries": True,
                "grain_threshold_deg": 5.0,
                "high_angle_threshold_deg": 15.0,
                "kam_threshold_deg": 5.0,
                "kam_order": 1,
                "grid_points": 56,
            },
        ),
        ExampleScenario(
            id="ebsd.example.kam_gradient",
            title="KAM across the same gradient",
            panel="ebsd",
            summary="The bicrystal again, coloured by kernel average misorientation.",
            teaches=(
                "The complement of the GROD example. KAM is flat across the whole gradient at "
                "half the per-step rotation — half, because two of the four kernel neighbours lie "
                "across the gradient and are identical to the centre point."
            ),
            operation="ebsd.map",
            request={
                "dataset": "bicrystal_gradient",
                "colouring": "kam",
                "ipf_direction": "Z",
                "modulate_by": "none",
                "modulation_floor": 0.25,
                "colour_map": "viridis",
                "show_boundaries": True,
                "grain_threshold_deg": 5.0,
                "high_angle_threshold_deg": 15.0,
                "kam_threshold_deg": 5.0,
                "kam_order": 1,
                "grid_points": 56,
            },
        ),
        ExampleScenario(
            id="ebsd.example.sigma3_boundaries",
            title="Twin boundaries on a grain map",
            panel="ebsd",
            summary="Coherent annealing twins, every boundary at 60 degrees about <111>.",
            teaches=(
                "That a boundary map maps misorientation rather than contrast: every segment here "
                "is a Sigma 3 twin at exactly 60 degrees, so the boundary misorientation "
                "histogram is a single spike and every segment classifies as high-angle."
            ),
            operation="ebsd.map",
            request={
                "dataset": "sigma3_twin",
                "colouring": "grain",
                "ipf_direction": "Z",
                "modulate_by": "none",
                "modulation_floor": 0.25,
                "colour_map": "viridis",
                "show_boundaries": True,
                "grain_threshold_deg": 5.0,
                "high_angle_threshold_deg": 15.0,
                "kam_threshold_deg": 5.0,
                "kam_order": 1,
                "grid_points": 56,
            },
        ),
    )
)
