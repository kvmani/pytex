"""The crystal viewer: a structure you can turn in your hands.

What it does
    Builds the three-dimensional scene of a structure — atoms, bonds, the unit
    cell, and any number of superimposed planes and directions — and sends it to
    the browser as geometry in Cartesian angstrom. The browser holds a camera and
    redraws on drag; nothing round-trips while the mouse is down.

When to use it
    To see what a plane or a direction actually *is* in a structure: which atoms
    a slip plane cuts, how a Burgers direction lies in it, how a family of planes
    is arranged around the cell. The calculator gives the numbers; this gives the
    picture the numbers describe.

The division of labour
    Every crystallographic decision is made here, in Python: which atoms exist,
    where the plane polygon is clipped, what a direction arrow's endpoints are,
    what colour a species takes. The browser receives finished vertices and
    applies a rotation and an orthographic divide — under fifty lines of viewing
    arithmetic, no crystallography. When a figure is exported, the camera comes
    back and :mod:`pytex.plotting.crystal3d` renders the same scene at
    publication quality, so the picture on screen and the figure in the paper are
    the same geometry seen the same way.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesListParameter,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import (
    direction_label,
    family_label,
    phase_parameter,
    plane_label,
)
from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex, Phase
from pytex.core.sphere import project_directions

__all__ = [
    "camera_angles_from_matrix",
    "camera_matrix_from_euler",
    "euler_from_camera_matrix",
    "orientation_overlay",
    "scene_payload",
]

_CITATION_VESTA = (
    "Momma & Izumi, VESTA 3, J. Appl. Crystallogr. 44 (2011) 1272 (visual conventions)."
)

_RENDER_STYLES = (
    (
        "ball_and_stick",
        "Ball and stick",
        "Spheres at covalent-radius scale joined by bonds. The default, and the clearest for "
        "seeing how planes cut a structure.",
    ),
    (
        "space_filling",
        "Space filling",
        "Atomic-radius spheres, bonds suppressed. Shows packing and free volume.",
    ),
    ("stick", "Stick", "Uniform thin bonds with small atoms. Best for open framework structures."),
    ("wireframe", "Wireframe", "Line bonds and marker atoms. Lightest to read behind overlays."),
    (
        "polyhedral",
        "Polyhedral",
        "Coordination polyhedra instead of individual bonds. The mineralogist's view.",
    ),
)

_HEX_COLOR = re.compile(r"^#[0-9a-fA-F]{6}$")
_APPEARANCE_DEFAULTS: dict[str, Any] = {
    "show_atoms": True,
    "show_bonds": True,
    "show_cells": True,
    "show_planes": True,
    "show_directions": True,
    "show_labels": True,
    "show_gizmo": True,
    "atom_scale": 1.0,
    "atom_opacity": 1.0,
    "surface_finish": "glossy",
    "light_direction": [-0.5, 0.5, 0.7071067811865476],
    "light_ambient": 0.42,
    "light_diffuse": 0.78,
    "light_specular": 0.38,
    "atom_shininess": 26.0,
    "depth_cue_strength": 0.18,
    "species_colors": {},
    "bond_color": "#64748b",
    "bond_width": 1.0,
    "bond_opacity": 0.85,
    "cell_color": "#64748b",
    "cell_width": 1.0,
    "cell_opacity": 0.5,
    "plane_color": "#0f766e",
    "plane_opacity": 0.28,
    "direction_color": "#2563eb",
    "direction_width": 1.0,
    "direction_opacity": 0.96,
    "annotation_scale": 1.0,
}


def _appearance(request_value: Any) -> dict[str, Any]:
    """Validate presentation-only crystal properties from the shared frontend."""

    raw = dict(request_value or {})
    unknown = sorted(set(raw) - set(_APPEARANCE_DEFAULTS))
    if unknown:
        raise InvalidInputError(
            f"Unknown crystal appearance setting(s): {', '.join(unknown)}.",
            field="appearance",
            hint="Reset object properties in the Crystal Viewer and try again.",
        )
    result = {**_APPEARANCE_DEFAULTS, **raw}
    for key in (
        "show_atoms",
        "show_bonds",
        "show_cells",
        "show_planes",
        "show_directions",
        "show_labels",
        "show_gizmo",
    ):
        if not isinstance(result[key], bool):
            raise InvalidInputError(f"appearance.{key} must be true or false.", field="appearance")
    if result["surface_finish"] not in {"flat", "matte", "glossy"}:
        raise InvalidInputError(
            "appearance.surface_finish must be 'flat', 'matte', or 'glossy'.",
            field="appearance",
        )
    bounds = {
        "atom_scale": (0.2, 2.5),
        "atom_opacity": (0.1, 1.0),
        "light_ambient": (0.05, 1.0),
        "light_diffuse": (0.0, 1.25),
        "light_specular": (0.0, 1.0),
        "atom_shininess": (2.0, 96.0),
        "depth_cue_strength": (0.0, 0.75),
        "bond_width": (0.2, 3.0),
        "bond_opacity": (0.05, 1.0),
        "cell_width": (0.2, 3.0),
        "cell_opacity": (0.05, 1.0),
        "plane_opacity": (0.02, 0.85),
        "direction_width": (0.2, 3.0),
        "direction_opacity": (0.05, 1.0),
        "annotation_scale": (0.5, 2.5),
    }
    for key, (minimum, maximum) in bounds.items():
        try:
            value = float(result[key])
        except (TypeError, ValueError) as error:
            raise InvalidInputError(
                f"appearance.{key} must be a number.", field="appearance"
            ) from error
        if not np.isfinite(value) or not minimum <= value <= maximum:
            raise InvalidInputError(
                f"appearance.{key} must be between {minimum:g} and {maximum:g}.",
                field="appearance",
            )
        result[key] = value
    try:
        light_direction = np.asarray(result["light_direction"], dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise InvalidInputError(
            "appearance.light_direction must contain three numbers.", field="appearance"
        ) from error
    if (
        light_direction.shape != (3,)
        or not np.all(np.isfinite(light_direction))
        or np.linalg.norm(light_direction) <= 1e-12
    ):
        raise InvalidInputError(
            "appearance.light_direction must be a finite, non-zero three-vector.",
            field="appearance",
        )
    result["light_direction"] = (light_direction / np.linalg.norm(light_direction)).tolist()
    for key in ("bond_color", "cell_color", "plane_color", "direction_color"):
        if not isinstance(result[key], str) or not _HEX_COLOR.fullmatch(result[key]):
            raise InvalidInputError(
                f"appearance.{key} must be a #RRGGBB colour.", field="appearance"
            )
        result[key] = result[key].lower()
    species_colors = result["species_colors"]
    if not isinstance(species_colors, dict) or len(species_colors) > 64:
        raise InvalidInputError(
            "appearance.species_colors must be an object with at most 64 entries.",
            field="appearance",
        )
    clean_species: dict[str, str] = {}
    for species, color in species_colors.items():
        if not isinstance(species, str) or not species.strip() or len(species) > 16:
            raise InvalidInputError(
                "Invalid species name in appearance colours.", field="appearance"
            )
        if not isinstance(color, str) or not _HEX_COLOR.fullmatch(color):
            raise InvalidInputError(
                f"The appearance colour for {species!r} must be #RRGGBB.", field="appearance"
            )
        clean_species[species] = color.lower()
    result["species_colors"] = clean_species
    return result


def _appearance_style(appearance: dict[str, Any], render_style: str) -> dict[str, Any]:
    atom_base = {
        "space_filling": 1.0,
        "stick": 0.24,
        "wireframe": 0.3,
        "polyhedral": 0.4,
    }.get(render_style, 0.55)
    bond_base = 1.0 if render_style == "stick" else 0.22
    finish = str(appearance["surface_finish"])
    finish_specular = {"flat": 0.0, "matte": 0.12, "glossy": 1.0}[finish]
    light_ambient = 1.0 if finish == "flat" else appearance["light_ambient"]
    light_diffuse = 0.0 if finish == "flat" else appearance["light_diffuse"]
    light_specular = 0.0 if finish == "flat" else appearance["light_specular"]
    crystal: dict[str, Any] = {
        "atom_radius_scale": atom_base * appearance["atom_scale"],
        "atom_alpha": appearance["atom_opacity"],
        "species_colors": appearance["species_colors"],
        "light_direction": appearance["light_direction"],
        "light_ambient": light_ambient,
        "light_diffuse": light_diffuse,
        "light_specular": light_specular,
        "atom_specular_strength": finish_specular,
        "bond_specular_strength": 0.28 * finish_specular,
        "atom_shininess": appearance["atom_shininess"],
        "bond_shininess": max(4.0, appearance["atom_shininess"] * 0.55),
        "depth_cue_strength": appearance["depth_cue_strength"],
        "bond_color_mode": "uniform",
        "bond_color": appearance["bond_color"],
        "bond_alpha": appearance["bond_opacity"],
        "bond_radius_scale": bond_base * appearance["bond_width"],
        "bond_radius": 1.2 * appearance["bond_width"],
        "cell_color": appearance["cell_color"],
        "cell_alpha": appearance["cell_opacity"],
        "cell_linewidth": appearance["cell_width"],
        "lattice_color": appearance["cell_color"],
        "lattice_linewidth": 1.5 * appearance["cell_width"],
        "plane_color": appearance["plane_color"],
        "plane_alpha": appearance["plane_opacity"],
        "direction_color": appearance["direction_color"],
        "direction_alpha": appearance["direction_opacity"],
        "direction_linewidth": 2.2 * appearance["direction_width"],
        "atom_label_fontsize": 10.0 * appearance["annotation_scale"],
        "plane_label_fontsize": 11.0 * appearance["annotation_scale"],
        "direction_label_fontsize": 11.0 * appearance["annotation_scale"],
    }
    if not appearance["show_atoms"]:
        crystal["atom_render_mode"] = "none"
    if not appearance["show_cells"]:
        crystal["cell_alpha"] = 0.0
        crystal["cell_linewidth"] = 0.0
        crystal["lattice_linewidth"] = 0.0
    if not appearance["show_labels"]:
        crystal["atom_label_fontsize"] = 0.01
        crystal["plane_label_fontsize"] = 0.01
        crystal["direction_label_fontsize"] = 0.01
    return {"crystal": crystal}


def _plane_overlay(indices: tuple[int, ...], phase: Phase) -> CrystalPlane:
    return CrystalPlane(
        miller=MillerIndex(indices=np.asarray(indices, dtype=int), phase=phase), phase=phase
    )


def _direction_overlay(indices: tuple[int, ...], phase: Phase) -> CrystalDirection:
    return CrystalDirection(coordinates=np.asarray(indices, dtype=float), phase=phase)


def _overlay_label(default: str, replacements: Sequence[str] | None, index: int) -> str:
    """The caller's label for overlay ``index``, or the scene's own."""

    if replacements is None or index >= len(replacements):
        return str(default)
    return str(replacements[index])


def camera_angles_from_matrix(matrix: Any) -> tuple[float, float]:
    """Convert a browser camera matrix into the renderer's elevation and azimuth.

    Purpose
    -------
    The interactive view holds its camera as a 3x3 rotation, because that is
    what a drag composes cleanly; the publication renderer is driven by an
    elevation and an azimuth. This is the one conversion between them, and it
    lives here rather than in JavaScript so that the exported figure and the
    on-screen view cannot disagree through two slightly different derivations.

    The third row of the camera matrix is the direction the viewer looks along,
    expressed in crystal coordinates. Elevation is its angle above the xy-plane
    and azimuth is its bearing within that plane, which is precisely how
    matplotlib's 3-D axes define ``elev`` and ``azim``.

    Parameters
    ----------
    matrix : array-like
        Nine numbers in row-major order, or a ``(3, 3)`` array.

    Returns
    -------
    tuple of (float, float)
        Elevation and azimuth in degrees.

    Raises
    ------
    InvalidInputError
        If the matrix is not nine finite numbers, or its view row is degenerate.
    """

    values = np.asarray(matrix, dtype=float).reshape(-1)
    if values.size != 9 or not np.all(np.isfinite(values)):
        raise InvalidInputError(
            "The camera matrix must be nine finite numbers in row-major order.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    forward = values.reshape(3, 3)[2]
    norm = float(np.linalg.norm(forward))
    if norm < 1e-9:
        raise InvalidInputError(
            "The camera matrix has a degenerate viewing direction.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    forward = forward / norm
    elevation = float(np.degrees(np.arcsin(np.clip(forward[2], -1.0, 1.0))))
    azimuth = float(np.degrees(np.arctan2(forward[1], forward[0])))
    return elevation, azimuth



# ---------------------------------------------------------------- orientation

#: The pole families a viewer offers, by crystal system.
#:
#: Three-index ``(hkl)`` throughout, because that is what
#: :class:`~pytex.core.lattice.CrystalPlane` takes; the labels the user sees are
#: produced by :func:`~pytex.app.services.calculator.family_label`, which writes
#: four-index Miller-Bravais for hexagonal and trigonal phases. The lists are
#: the families a stereographic projection of that system is conventionally read
#: against, not an exhaustive catalogue: a pole figure of six families is
#: unreadable, and any plane the user actually cares about can be added as a
#: scene overlay, which appears in the figure in its own right.
_POLE_FAMILIES: dict[str, tuple[tuple[int, int, int], ...]] = {
    "cubic": ((1, 0, 0), (1, 1, 0), (1, 1, 1)),
    "hexagonal": ((0, 0, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0)),
    "trigonal": ((0, 0, 1), (1, 0, 0), (1, 0, 1)),
    "tetragonal": ((0, 0, 1), (1, 0, 0), (1, 1, 0)),
    "orthorhombic": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    "monoclinic": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    "triclinic": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
}
_DEFAULT_POLE_FAMILIES: tuple[tuple[int, int, int], ...] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

#: The specimen frame the viewer works in: the screen itself.
#:
#: ``RD`` is screen right, ``TD`` is screen up, ``ND`` points out of the screen
#: towards the viewer. That is a right-handed triad, and it is
#: :data:`~pytex.core.frame_catalog.SAMPLE_RD_TD_ND_FRAME` with its default
#: identity axis vectors, so no new frame is invented here. It also matches how
#: the texture panel already draws a pole figure -- RD to the right, TD up --
#: which is what lets a reader move between the two panels without relearning
#: the picture.
_VIEWER_SPECIMEN_AXES: tuple[tuple[str, str, tuple[float, float, float], str], ...] = (
    ("rd", "RD", (1.0, 0.0, 0.0), "Screen right."),
    ("td", "TD", (0.0, 1.0, 0.0), "Screen up."),
    ("nd", "ND", (0.0, 0.0, 1.0), "Out of the screen, towards you."),
)

#: Points per sector edge in the projected standard-triangle outline. Enough
#: that the arc stays smooth at any size the dock reaches, and small enough that
#: the whole boundary is a few hundred numbers on the wire.
_SECTOR_ARC_SAMPLES = 48

#: Index bound for the direction *labels* in the orientation readout.
#:
#: Three, not the eight :func:`~pytex.core.miller.nearest_low_index_direction`
#: defaults to. A general orientation puts no specimen axis on a low-index
#: direction, so the search's job here is to name the nearest signpost, not to
#: match the direction: on a hexagonal phase a bound of eight buys three degrees
#: of accuracy and pays for it with ``[11 -1 -10 12]``, which nobody can read
#: and nobody wanted. The label comes with its angular residual attached, so the
#: approximation is stated rather than implied.
_READOUT_INDEX_BOUND = 3

#: Named ideal orientations offered as one-press presets, by crystal system.
#:
#: Only cubic, because the catalogue in :mod:`pytex.texture.components` is the
#: rolling-texture catalogue of cubic metals: "Goss" names a specific
#: relationship between ``{011}``, ``<100>`` and the rolling geometry, and
#: pressing it on a hexagonal phase would set an orientation the name does not
#: describe. Systems without a catalogue get none, rather than a plausible-
#: looking wrong one.
_PRESET_SYSTEMS = frozenset({"cubic"})


def _preset_components(crystal_system: str) -> list[dict[str, Any]]:
    """Named ideal orientations for one crystal system, as Bunge angle triples."""

    if crystal_system.lower() not in _PRESET_SYSTEMS:
        return []
    from pytex.texture.components import (
        BRASS,
        COPPER,
        CUBE,
        GOSS,
        ROTATED_CUBE,
        ROTATED_GOSS,
        S_COMPONENT,
    )

    return [
        {
            "name": component.name.replace("_", " "),
            "label": component.miller_label,
            "angles_deg": [float(value) for value in component.bunge_euler_deg],
        }
        for component in (CUBE, ROTATED_CUBE, GOSS, ROTATED_GOSS, BRASS, COPPER, S_COMPONENT)
    ]

_EULER_CONVENTIONS = (
    (
        "bunge",
        "Bunge (phi1, Phi, phi2)",
        "The ZXZ convention of Bunge, used by essentially all texture software and EBSD "
        "vendors. Choose this unless you have a specific reason not to.",
    ),
    (
        "matthies",
        "Matthies (alpha, beta, gamma)",
        "The ZYZ convention of Matthies. The same rotation, a different angle triple.",
    ),
)

_CITATION_BUNGE = (
    "Bunge, Texture Analysis in Materials Science, Butterworths (1982) "
    "(Euler-angle convention and pole-figure definitions)."
)
_CITATION_RANDLE = (
    "Randle & Engler, Introduction to Texture Analysis, 2nd ed., CRC Press (2009) "
    "(stereographic projection and the standard triangle)."
)


def _symmetry_orbit(operators: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Every symmetry image of one crystal direction, with its antipodes.

    Both signs are kept because a pole figure is antipodal: the browser draws
    whichever member of each pair currently points out of the screen, so no
    hemisphere folding has to happen while the mouse is down.
    """

    orbit = np.einsum("nij,j->ni", np.asarray(operators, dtype=float), vector)
    orbit = np.vstack([orbit, -orbit])
    norms = np.linalg.norm(orbit, axis=1)
    keep = norms > 1e-12
    orbit = orbit[keep] / norms[keep][:, None]
    _, unique = np.unique(np.round(orbit, 6), axis=0, return_index=True)
    return np.ascontiguousarray(orbit[np.sort(unique)])


def _sector_outline(vertices: np.ndarray, *, method: str) -> list[list[float]]:
    """The standard triangle's boundary, projected, as one closed polyline.

    Each edge of the fundamental sector is a great-circle arc between adjacent
    corner directions, so the boundary is walked on the sphere and the samples
    projected -- not drawn by joining the projected corners with straight lines,
    which cuts visibly inside the true boundary of a cubic triangle.
    """

    corners = np.asarray(vertices, dtype=float)
    if corners.shape[0] < 3:
        return []
    fractions = np.linspace(0.0, 1.0, _SECTOR_ARC_SAMPLES, endpoint=False)
    samples: list[np.ndarray] = []
    for index in range(corners.shape[0]):
        start = corners[index]
        end = corners[(index + 1) % corners.shape[0]]
        angle = float(np.arccos(np.clip(float(np.dot(start, end)), -1.0, 1.0)))
        if angle < 1e-9:
            continue
        arc = (
            np.sin((1.0 - fractions)[:, None] * angle) * start
            + np.sin(fractions[:, None] * angle) * end
        ) / np.sin(angle)
        samples.append(arc)
    if not samples:
        return []
    points = np.vstack(samples)
    points = points / np.linalg.norm(points, axis=1)[:, None]
    projected = np.asarray(project_directions(points, method=method), dtype=float)
    return [[float(x), float(y)] for x, y in projected]


def orientation_overlay(
    spec: PhaseSpec,
    *,
    plane_rows: Sequence[Sequence[int]] = (),
    direction_rows: Sequence[Sequence[int]] = (),
    method: str = "stereographic",
) -> dict[str, Any]:
    """Everything the browser needs to draw orientation figures while dragging.

    Purpose
    -------
    A pole figure that only updates when the mouse is released is not a pole
    figure of the view; it is a pole figure of the last view. To keep it live,
    the browser must be able to redraw without asking Python anything -- and
    without deciding any crystallography either. This function is that split
    made concrete: it sends the *crystallographic content* of the figures once,
    in the Cartesian crystal frame, and the browser thereafter only multiplies
    by the camera and projects, exactly as it already does for atoms.

    When and where to use it
    ------------------------
    Called by the ``crystal.scene`` operation, whose payload it joins under
    ``orientation``. It is not itself a user-facing calculation: every number a
    user reads comes from the ``crystal.orientation`` operation, which computes
    it in Python for the camera the browser reports.

    Parameters
    ----------
    spec : PhaseSpec
        The phase being viewed; supplies the lattice and the point group.
    plane_rows, direction_rows : sequence of index sequences
        The overlays already drawn in the 3-D scene. Their poles and directions
        are carried through so the figures show the same features as the
        structure, labelled the same way.
    method : str
        ``"stereographic"`` (the default, and what the dock draws) or
        ``"equal_area"``; passed to
        :func:`~pytex.core.sphere.project_directions` for the sector outline.

    Returns
    -------
    dict
        ``point_group``, ``crystal_system``, the proper-rotation ``operators``
        as flat row-major triples, the ``pole_families`` with their full
        symmetry orbits, the ``overlay_poles`` and ``overlay_directions``, the
        fundamental ``sector`` (corner directions, inward edge normals, labelled
        corners, and the projected outline), the ``specimen_axes``, the offered
        ``euler_conventions``, and the named ``components`` that can be set with
        one press.

    See Also
    --------
    camera_matrix_from_euler, euler_from_camera_matrix :
        The conversions between the browser's camera and an orientation.
    """

    from pytex.core.miller import nearest_low_index_direction

    phase = spec.to_phase()
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    families = _POLE_FAMILIES.get(spec.crystal_system.lower(), _DEFAULT_POLE_FAMILIES)

    pole_families: list[dict[str, Any]] = []
    for indices in families:
        normal = np.asarray(_plane_overlay(tuple(indices), phase).normal, dtype=float)
        orbit = _symmetry_orbit(operators, normal)
        pole_families.append(
            {
                "key": ",".join(str(value) for value in indices),
                "label": family_label(indices, spec=spec, family="plane"),
                "indices": [int(value) for value in indices],
                "vectors": [[float(value) for value in row] for row in orbit],
            }
        )

    overlay_poles = [
        {
            "label": plane_label(row, spec=spec),
            "vector": [float(value) for value in _plane_overlay(tuple(row), phase).normal],
        }
        for row in plane_rows
    ]
    overlay_directions = [
        {
            "label": direction_label(row, spec=spec),
            "vector": [
                float(value) for value in _direction_overlay(tuple(row), phase).unit_vector
            ],
        }
        for row in direction_rows
    ]

    sector = phase.symmetry.fundamental_sector(antipodal=True)
    sector_vertices = np.asarray(sector.vertices, dtype=float)
    corners = []
    for vertex in sector_vertices:
        corner_indices, residual_deg = nearest_low_index_direction(vertex, phase=phase)
        corners.append(
            {
                "label": direction_label(
                    tuple(int(value) for value in corner_indices), spec=spec
                ),
                "vector": [float(value) for value in vertex],
                "residual_deg": float(residual_deg),
            }
        )

    return {
        "point_group": spec.point_group,
        "crystal_system": spec.crystal_system,
        "projection": method,
        "operators": [[float(value) for value in row.reshape(-1)] for row in operators],
        "pole_families": pole_families,
        "overlay_poles": overlay_poles,
        "overlay_directions": overlay_directions,
        "sector": {
            "vertices": [[float(value) for value in row] for row in sector_vertices],
            "edge_normals": [
                [float(value) for value in row]
                for row in np.asarray(sector.edge_normals, dtype=float)
            ],
            "corners": corners,
            "outline": _sector_outline(sector_vertices, method=method),
        },
        "specimen_axes": [
            {"key": key, "label": label, "vector": list(vector), "help_text": help_text}
            for key, label, vector, help_text in _VIEWER_SPECIMEN_AXES
        ],
        "euler_conventions": [
            {"key": key, "label": label, "help_text": help_text}
            for key, label, help_text in _EULER_CONVENTIONS
        ],
        "components": _preset_components(spec.crystal_system),
    }


def _euler_convention(value: Any) -> str:
    """Validate a convention name against the two the viewer offers."""

    name = str(value).strip().lower()
    known = {key for key, _label, _help in _EULER_CONVENTIONS}
    if name not in known:
        raise InvalidInputError(
            f"Unknown Euler-angle convention {value!r}.",
            field="euler_convention",
            hint="Available: " + ", ".join(sorted(known)) + ".",
        )
    return name


def camera_matrix_from_euler(
    angle1: float,
    angle2: float,
    angle3: float,
    *,
    convention: str = "bunge",
) -> list[float]:
    """Build the viewer's camera matrix from a triple of Euler angles.

    Purpose
    -------
    "Set the view to the cube component" is an orientation statement, and the
    viewer's camera *is* an orientation: with the screen taken as the specimen
    frame -- RD right, TD up, ND out of the screen -- the camera matrix ``C``
    satisfies ``v_specimen = C v_crystal``, which is exactly PyTex's
    crystal-to-specimen convention for
    :meth:`~pytex.core.orientation.Orientation.as_matrix`. So this is not a new
    derivation; it is :class:`~pytex.core.orientation.Rotation` written out in
    the nine numbers the browser holds.

    When and where to use it
    ------------------------
    Behind the crystal viewer's Euler-angle entry. Use
    :func:`euler_from_camera_matrix` for the return trip, and derive neither in
    JavaScript: one convention, defined once, is what keeps the on-screen view,
    the reported angles, and the exported figure describing the same
    orientation.

    Parameters
    ----------
    angle1, angle2, angle3 : float
        The angle triple in degrees, in the order the convention names them --
        ``(phi1, Phi, phi2)`` for Bunge.
    convention : str
        ``"bunge"`` (ZXZ, the default) or ``"matthies"`` (ZYZ).

    Returns
    -------
    list of float
        Nine numbers in row-major order, ready to be assigned to the camera.

    Raises
    ------
    InvalidInputError
        If an angle is not finite, or the convention is not one of the two
        offered.
    """

    from pytex.core.orientation import Rotation

    angles = np.asarray([angle1, angle2, angle3], dtype=float)
    if not np.all(np.isfinite(angles)):
        raise InvalidInputError(
            "Euler angles must be finite numbers in degrees.",
            field="phi1",
            hint="Enter three angles in degrees.",
        )
    rotation = Rotation.from_euler(
        float(angles[0]),
        float(angles[1]),
        float(angles[2]),
        convention=_euler_convention(convention),
        degrees=True,
    )
    return [float(value) for value in np.asarray(rotation.as_matrix(), dtype=float).reshape(-1)]


def euler_from_camera_matrix(
    matrix: Any,
    *,
    convention: str = "bunge",
) -> tuple[float, float, float]:
    """Read the Euler angles of the viewer's camera.

    Purpose
    -------
    The inverse of :func:`camera_matrix_from_euler`, and the only place the
    viewer's rotation becomes an angle triple.

    Method
    ------
    The matrix is orthonormalized before it is read. A camera accumulated from
    thousands of drag increments drifts from orthogonality by a part in
    ``1e-12`` or so -- harmless on screen, and enough to put a spurious final
    digit on a reported angle. The nearest rotation in the Frobenius sense is
    the polar factor, obtained here from the singular-value decomposition.

    Parameters
    ----------
    matrix : array-like
        Nine finite numbers in row-major order, or a ``(3, 3)`` array.
    convention : str
        ``"bunge"`` (ZXZ, the default) or ``"matthies"`` (ZYZ).

    Returns
    -------
    tuple of (float, float, float)
        The angle triple in degrees, each wrapped into ``[0, 360)``.

    Raises
    ------
    InvalidInputError
        If the matrix is not nine finite numbers, or is too far from a rotation
        to be a drifted camera.
    """

    from pytex.core.orientation import Rotation

    values = np.asarray(matrix, dtype=float).reshape(-1)
    if values.size != 9 or not np.all(np.isfinite(values)):
        raise InvalidInputError(
            "The camera matrix must be nine finite numbers in row-major order.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    name = _euler_convention(convention)
    square = values.reshape(3, 3)
    left, _singular_values, right = np.linalg.svd(square)
    nearest = left @ right
    if float(np.linalg.det(nearest)) < 0.0 or float(np.max(np.abs(nearest - square))) > 1e-3:
        raise InvalidInputError(
            "The camera matrix is not a rotation.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    angles = Rotation.from_matrix(nearest).to_euler(convention=name, degrees=True)
    # A wrap into [0, 360) turns a decomposition that lands a hair below zero
    # into 359.999999999, which the viewer then shows as "360.00 deg" -- an angle
    # that is not in the half-open range it claims to be in, next to two angles
    # that are exact. Round to a picodegree first, so a tiny negative becomes a
    # clean zero rather than a full turn.
    return tuple(float(np.round(value, 9) % 360.0) for value in angles)  # type: ignore[return-value]

def scene_payload(
    scene: Any,
    *,
    spec: PhaseSpec,
    plane_labels: Sequence[str] | None = None,
    direction_labels: Sequence[str] | None = None,
    plane_rows: Sequence[Sequence[int]] = (),
    direction_rows: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    """Convert a :class:`CrystalScene` into the JSON the browser draws.

    Only what a renderer needs crosses the wire: positions, radii, colours,
    labels, and the polygons already clipped to the cell. Everything is in
    Cartesian angstrom in the crystal frame, so the browser's camera is a pure
    rotation with no crystallography hidden inside it.

    Each atom and each bond carries enough of its own description to answer a
    hover — species, occupancy, bond length — because a scene the user cannot
    interrogate is a picture rather than an instrument.

    Overlay labels are replaced by ``plane_labels`` and ``direction_labels``
    when given, in scene order. The scene builder labels for matplotlib, whose
    mathtext ``$(1\\bar{1}0)$`` is markup rather than text; the browser draws
    the string literally, so it must be handed the plain form the rest of the
    application already uses in tables and prose.
    """

    bounds = np.asarray(scene.bounds(), dtype=float)
    centre = (bounds[0] + bounds[1]) / 2.0
    radius = float(np.linalg.norm(bounds[1] - bounds[0]) / 2.0) or 1.0

    atoms = [
        {
            "position": [float(value) for value in atom.position_angstrom],
            "radius": float(atom.radius_angstrom),
            "color": str(atom.color),
            "species": str(atom.species),
            "occupancy": float(atom.occupancy),
            "label": atom.label,
        }
        for atom in scene.atoms
    ]
    bonds = [
        {
            "start": [float(value) for value in bond.start_angstrom],
            "end": [float(value) for value in bond.end_angstrom],
            "color": str(bond.color),
            "start_color": bond.start_color,
            "end_color": bond.end_color,
            "length": float(bond.length_angstrom),
            "species": f"{bond.start_species or '?'}-{bond.end_species or '?'}",
        }
        for bond in scene.bonds
    ]
    cell_edges = [
        [[float(value) for value in point] for point in edge] for edge in scene.lattice_edges
    ]
    for cell in scene.cells:
        cell_edges.extend(
            [[float(value) for value in point] for point in edge] for edge in cell.edges_angstrom
        )
    planes = [
        {
            "vertices": [[float(value) for value in point] for point in plane.vertices_angstrom],
            "normal": [float(value) for value in plane.normal_angstrom],
            "color": str(plane.color),
            "alpha": float(plane.alpha),
            "label": _overlay_label(plane.label, plane_labels, index),
        }
        for index, plane in enumerate(scene.planes)
    ]
    directions = [
        {
            "start": [float(value) for value in direction.start_angstrom],
            "end": [float(value) for value in direction.end_angstrom],
            "color": str(direction.color),
            "label": _overlay_label(direction.label, direction_labels, index),
        }
        for index, direction in enumerate(scene.directions)
    ]
    return {
        "atoms": atoms,
        "bonds": bonds,
        "cell_edges": cell_edges,
        "planes": planes,
        "directions": directions,
        "centre": [float(value) for value in centre],
        "radius": radius,
        "bounds": [[float(value) for value in row] for row in bounds],
        "phase": spec.to_json(),
        "axes": _axis_arrows(spec),
        "orientation": orientation_overlay(
            spec, plane_rows=plane_rows, direction_rows=direction_rows
        ),
    }


def _axis_arrows(spec: PhaseSpec) -> list[dict[str, Any]]:
    """The a, b, c axes as arrows, so the view always says which way is which.

    Without them a rotated structure is unreadable: every cubic cell looks the
    same from every direction, and the whole point of turning it is to know
    where you now are.
    """

    phase = spec.to_phase()
    basis = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    return [
        {
            "label": name,
            "vector": [float(value) for value in basis[:, index]],
        }
        for index, name in enumerate(("a", "b", "c"))
    ]


@REGISTRY.operation(
    "crystal.scene",
    title="Crystal structure",
    summary="The 3D structure, with any number of planes and directions superimposed.",
    help_text=(
        "Builds the structure as a scene you can rotate: drag to turn it, scroll to zoom, and "
        "hover any atom or bond for its identity and its bond length.\n\n"
        "Enter as many planes and directions as you like — each plane is clipped to the cell and "
        "drawn as a translucent polygon with its indices, and each direction is drawn as a "
        "labelled arrow from the cell origin. Superimposing a slip plane and its Burgers "
        "direction is the canonical use: the arrow must lie *in* the polygon, and seeing that it "
        "does is worth more than checking that the zone-law integer is zero.\n\n"
        "The render style changes the whole visual system at once. Ball-and-stick reads best "
        "behind plane overlays; space-filling shows packing; polyhedral is the mineralogist's "
        "view of a framework structure."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="render_style",
            label="Style",
            help_text="The visual system for atoms and bonds.",
            options=_RENDER_STYLES,
            default="ball_and_stick",
        ),
        IntegerParameter(
            name="repeat_a",
            label="Repeats along a",
            help_text="How many cells to build along a. More cells show longer-range structure.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="a",
            row="Repeats",
        ),
        IntegerParameter(
            name="repeat_b",
            label="Repeats along b",
            help_text="How many cells to build along b.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="b",
            row="Repeats",
        ),
        IntegerParameter(
            name="repeat_c",
            label="Repeats along c",
            help_text="How many cells to build along c.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="c",
            row="Repeats",
        ),
        IndicesListParameter(
            name="planes",
            label="Planes to superimpose (hkl)",
            help_text=(
                "One plane per row. Each is drawn as a translucent polygon clipped to the cell "
                "and labelled with its indices. Leave empty for none."
            ),
            required=False,
            group="Overlays",
        ),
        IndicesListParameter(
            name="directions",
            label="Directions to superimpose [uvw]",
            help_text="One direction per row, drawn as a labelled arrow from the cell origin.",
            required=False,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_bonds",
            label="Draw bonds",
            help_text=(
                "Bonds are inferred from covalent radii plus a tolerance, not read from a file, "
                "so they are an aid to reading the picture rather than a statement about "
                "chemistry."
            ),
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_unit_cells",
            label="Outline every cell",
            help_text="Draw the edges of each repeated cell, not only the outer box.",
            default=False,
            group="Overlays",
        ),
        BooleanParameter(
            name="hexagonal_prism",
            label="Draw hexagonal phases as the prism",
            help_text=(
                "Draw a hexagonal phase as its hexagonal prism — three cells about one "
                "atomic column — rather than as the 120-degree rhombus the lattice is "
                "written on. The prism is the figure the sixfold symmetry is visible in; the "
                "rhombus shows a third of it. Plane overlays follow, so a basal plane is "
                "clipped to the hexagon. Ignored for phases that are not on hexagonal axes."
            ),
            default=True,
            group="Extent",
        ),
        ChoiceParameter(
            name="atom_labels",
            label="Atom labels",
            help_text="Label every atom with its element or its site name.",
            options=(
                ("none", "None", "No labels; cleanest for dense structures."),
                ("species", "Element", "The element symbol on each atom."),
                ("site", "Site", "The crystallographic site label, as in the CIF."),
            ),
            default="none",
            advanced=True,
            group="Overlays",
        ),
        NumberParameter(
            name="bond_tolerance_angstrom",
            label="Bond tolerance",
            help_text=(
                "Added to the sum of covalent radii when deciding whether two atoms are bonded. "
                "Raise it to catch long bonds; raise it too far and everything bonds to "
                "everything."
            ),
            units="Å",
            default=0.45,
            minimum=0.0,
            maximum=2.0,
            advanced=True,
            group="Overlays",
        ),
    ),
    returns="The scene under `data.scene`; the atom list as the table.",
    panel="crystal",
    citations=(_CITATION_VESTA,),
    tags=("crystal", "structure", "3D", "viewer", "plane", "direction", "unit cell", "bonds"),
)
def _crystal_scene(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.plotting.crystal3d import build_crystal_scene

    spec, phase = phase_from_request(request["phase"])
    if not spec.has_structure:
        raise InvalidInputError(
            f"{spec.name} carries no atomic basis, so there is nothing to draw.",
            field="phase",
            hint=(
                "Choose a built-in phase, or add atomic sites to the phase description. Lattice "
                "geometry alone supports the calculator but not the viewer."
            ),
        )
    repeats = (int(request["repeat_a"]), int(request["repeat_b"]), int(request["repeat_c"]))
    plane_rows = tuple(request.get("planes") or ())
    direction_rows = tuple(request.get("directions") or ())

    scene = build_crystal_scene(
        phase,
        repeats=repeats,
        render_style=str(request["render_style"]),
        show_bonds=bool(request["show_bonds"]),
        bond_tolerance_angstrom=float(request["bond_tolerance_angstrom"]),
        show_unit_cells=bool(request["show_unit_cells"]),
        hexagonal_prism=bool(request["hexagonal_prism"]),
        atom_label_mode=str(request["atom_labels"]),
        plane_overlays=tuple(_plane_overlay(row, phase) for row in plane_rows),
        direction_overlays=tuple(_direction_overlay(row, phase) for row in direction_rows),
    )
    plane_texts = tuple(plane_label(row, spec=spec) for row in plane_rows)
    direction_texts = tuple(direction_label(row, spec=spec) for row in direction_rows)
    payload = scene_payload(
        scene,
        spec=spec,
        plane_labels=plane_texts,
        direction_labels=direction_texts,
        plane_rows=plane_rows,
        direction_rows=direction_rows,
    )

    rows = tuple(
        {
            "index": index + 1,
            "species": atom["species"],
            "x": atom["position"][0],
            "y": atom["position"][1],
            "z": atom["position"][2],
            "occupancy": atom["occupancy"],
        }
        for index, atom in enumerate(payload["atoms"])
    )
    bond_summary = {
        f"{pair[0]}-{pair[1]}": stats for pair, stats in scene.bond_length_summary().items()
    }
    overlay_text = ", ".join([*plane_texts, *direction_texts])
    result = AppResult(
        title=f"{spec.name}: crystal structure",
        summary=(
            f"{len(rows)} atoms in a {repeats[0]}x{repeats[1]}x{repeats[2]} block of "
            f"{spec.name} ({spec.crystal_system}, {spec.point_group})"
            + (f", with {overlay_text} superimposed" if overlay_text else "")
            + ". Positions are Cartesian angstrom in the crystal frame. Bonds are inferred from "
            "covalent radii plus the stated tolerance, so they aid reading rather than assert "
            "chemistry."
        ),
        table=ResultTable(
            columns=(
                Column("index", "#", numeric=True),
                Column("species", "Element"),
                Column("x", "x", units="Å", numeric=True, digits=5),
                Column("y", "y", units="Å", numeric=True, digits=5),
                Column("z", "z", units="Å", numeric=True, digits=5),
                Column("occupancy", "Occupancy", numeric=True, digits=3),
            ),
            rows=rows,
            caption=f"Atom positions in the displayed block of {spec.name}.",
        ),
        data={
            "scene": payload,
            "bond_summary": bond_summary,
            "repeats": list(repeats),
        },
        inputs={
            "phase": spec.to_json(),
            "render_style": request["render_style"],
            "repeat_a": repeats[0],
            "repeat_b": repeats[1],
            "repeat_c": repeats[2],
            "planes": [list(row) for row in plane_rows],
            "directions": [list(row) for row in direction_rows],
            "show_bonds": bool(request["show_bonds"]),
            "show_unit_cells": bool(request["show_unit_cells"]),
            "hexagonal_prism": bool(request["hexagonal_prism"]),
            "atom_labels": request["atom_labels"],
            "bond_tolerance_angstrom": float(request["bond_tolerance_angstrom"]),
        },
        citations=(_CITATION_VESTA,),
    )
    return result.to_json()



@REGISTRY.operation(
    "crystal.orientation",
    title="Orientation of the view",
    summary="Euler angles, poles, and the crystal direction along each specimen axis.",
    help_text=(
        "The crystal viewer's camera is an orientation. Turning the structure by hand and "
        "setting it by Euler angles are the same act, and this operation is the conversion "
        "between them — in both directions, in one place, so the picture and the numbers "
        "cannot describe different orientations.\n\n"
        "**The frame.** The screen is the specimen frame: RD points right, TD points up, and "
        "ND points out of the screen towards you. That is a right-handed triad, and the "
        "returned camera matrix `C` satisfies `v_specimen = C v_crystal`, PyTex's "
        "crystal-to-specimen convention throughout.\n\n"
        "**A deviation, stated.** Pole figures are conventionally drawn with RD at the top. "
        "Here RD is to the right, because the pole figure beside the structure is the *same "
        "view* of the *same crystal* — turning the structure turns the pole figure with it — "
        "and rotating the projection by ninety degrees would break exactly the correspondence "
        "the figure exists to show. The texture panel draws its pole figures the same way.\n\n"
        "**What comes back.** The angle triple in both offered conventions, the misorientation "
        "from the identity as an axis and angle, where each requested pole family sits in the "
        "specimen frame, and — the table — which crystal direction lies along RD, TD and ND. "
        "That last is the inverse pole figure written out: the direction is folded into the "
        "fundamental sector and labelled with the nearest low-index direction, with the "
        "angular residual of that label stated rather than hidden."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="euler_convention",
            label="Euler convention",
            help_text="Which axis sequence the three angles name.",
            options=_EULER_CONVENTIONS,
            default="bunge",
        ),
        NumberParameter(
            name="angle1",
            label="First Euler angle",
            help_text="First Euler angle, in degrees.",
            units="deg",
            default=0.0,
            minimum=-360.0,
            maximum=720.0,
            row="Euler angles",
        ),
        NumberParameter(
            name="angle2",
            label="Second Euler angle",
            help_text="Second Euler angle, in degrees.",
            units="deg",
            default=0.0,
            minimum=-360.0,
            maximum=720.0,
            row="Euler angles",
        ),
        NumberParameter(
            name="angle3",
            label="Third Euler angle",
            help_text="Third Euler angle, in degrees.",
            units="deg",
            default=0.0,
            minimum=-360.0,
            maximum=720.0,
            row="Euler angles",
        ),
        TextParameter(
            name="camera_matrix",
            label="Camera matrix",
            help_text=(
                "Nine numbers, row-major, for the rotation the viewer currently holds. Sent by "
                "the viewer as you turn the structure; when it is given it wins, and the three "
                "angles are outputs rather than inputs. Leave empty to build the orientation "
                "from the angles."
            ),
            required=False,
            default="",
            advanced=True,
        ),
        IndicesListParameter(
            name="poles",
            label="Poles to locate (hkl)",
            help_text=(
                "One plane per row. Each is expanded over the point group and reported in the "
                "specimen frame. Leave empty to use the conventional families for the crystal "
                "system."
            ),
            required=False,
            group="Poles",
        ),
        ChoiceParameter(
            name="projection",
            label="Projection",
            help_text="How the sphere is flattened onto the page.",
            options=(
                (
                    "stereographic",
                    "Stereographic (Wulff)",
                    "Conformal: angles and circles are preserved. What the viewer draws.",
                ),
                (
                    "equal_area",
                    "Equal area (Schmidt)",
                    "Area is proportional to solid angle. The right choice for densities.",
                ),
            ),
            default="stereographic",
        ),
    ),
    returns=(
        "The crystal direction along each specimen axis as the table; the camera matrix, both "
        "angle triples, and the projected poles under `data`."
    ),
    panel="crystal",
    citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    tags=(
        "crystal",
        "orientation",
        "Euler angles",
        "Bunge",
        "pole figure",
        "inverse pole figure",
        "stereographic",
    ),
)
def _crystal_orientation(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.miller import nearest_low_index_direction

    spec, phase = phase_from_request(request["phase"])
    convention = _euler_convention(request["euler_convention"])
    method = str(request["projection"])
    matrix_text = str(request.get("camera_matrix") or "").strip()

    if matrix_text:
        values = _camera_matrix_values(matrix_text)
        angles = euler_from_camera_matrix(values, convention=convention)
        camera = [float(value) for value in np.asarray(values, dtype=float).reshape(-1)]
        # Round-trip through the conversion just applied, so the reported matrix
        # is the rotation the reported angles name rather than the drifted one
        # the drag accumulated. The two differ by parts in 1e-12; reporting both
        # would invite a reader to wonder which is authoritative.
        camera = camera_matrix_from_euler(*angles, convention=convention)
        source = "the camera"
    else:
        angles = (
            float(request["angle1"]),
            float(request["angle2"]),
            float(request["angle3"]),
        )
        camera = camera_matrix_from_euler(*angles, convention=convention)
        angles = euler_from_camera_matrix(camera, convention=convention)
        source = "the angles entered"

    rotation_matrix = np.asarray(camera, dtype=float).reshape(3, 3)
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    sector = phase.symmetry.fundamental_sector(antipodal=True)

    pole_rows = tuple(request.get("poles") or ()) or _POLE_FAMILIES.get(
        spec.crystal_system.lower(), _DEFAULT_POLE_FAMILIES
    )
    pole_points: list[dict[str, Any]] = []
    for indices in pole_rows:
        normal = np.asarray(_plane_overlay(tuple(indices), phase).normal, dtype=float)
        orbit = _symmetry_orbit(operators, normal)
        specimen = orbit @ rotation_matrix.T
        upper = specimen[specimen[:, 2] >= -1e-12]
        projected = (
            np.asarray(project_directions(upper, method=method), dtype=float)
            if upper.shape[0]
            else np.zeros((0, 2))
        )
        pole_points.append(
            {
                "label": family_label(indices, spec=spec, family="plane"),
                "indices": [int(value) for value in indices],
                "points": [
                    {
                        "x": float(projected[index, 0]),
                        "y": float(projected[index, 1]),
                        "specimen": [float(value) for value in upper[index]],
                    }
                    for index in range(upper.shape[0])
                ],
            }
        )

    rows: list[dict[str, Any]] = []
    ipf_points: list[dict[str, Any]] = []
    for key, label, vector, _help_text in _VIEWER_SPECIMEN_AXES:
        axis = np.asarray(vector, dtype=float)
        # The specimen axis in crystal coordinates: C maps crystal to specimen,
        # so its transpose brings the specimen axis home to the crystal.
        crystal = rotation_matrix.T @ axis
        reduced = np.asarray(
            phase.symmetry.reduce_vectors_to_fundamental_sector(crystal[None, :], antipodal=True),
            dtype=float,
        ).reshape(3)
        axis_indices, residual_deg = nearest_low_index_direction(
            reduced, phase=phase, max_index=_READOUT_INDEX_BOUND
        )
        projected = np.asarray(
            project_directions(reduced[None, :], method=method), dtype=float
        ).reshape(2)
        polar_deg = float(np.degrees(np.arccos(np.clip(reduced[2], -1.0, 1.0))))
        # On the pole itself the azimuth is undefined, and arctan2 of two
        # rounding errors reports it with a straight face. Say zero instead.
        azimuth_deg = (
            0.0
            if min(polar_deg, 180.0 - polar_deg) < 1e-6
            else float(np.degrees(np.arctan2(reduced[1], reduced[0]))) % 360.0
        )
        rows.append(
            {
                "axis": label,
                "direction": direction_label(
                    tuple(int(value) for value in axis_indices), spec=spec
                ),
                "residual_deg": float(residual_deg),
                "u": float(reduced[0]),
                "v": float(reduced[1]),
                "w": float(reduced[2]),
                "polar_deg": polar_deg,
                "azimuth_deg": azimuth_deg,
            }
        )
        ipf_points.append(
            {
                "key": key,
                "label": label,
                "direction": rows[-1]["direction"],
                "residual_deg": float(residual_deg),
                "crystal": [float(value) for value in reduced],
                "x": float(projected[0]),
                "y": float(projected[1]),
            }
        )

    from pytex.core.orientation import Rotation

    rotation = Rotation.from_matrix(rotation_matrix)
    axis_vector = rotation.axis
    angle_deg = rotation.angle_deg
    other = "matthies" if convention == "bunge" else "bunge"
    other_angles = euler_from_camera_matrix(camera, convention=other)
    angle_names = (
        ("phi1", "Phi", "phi2") if convention == "bunge" else ("alpha", "beta", "gamma")
    )
    angle_text = ", ".join(
        f"{name} = {value:.2f} deg" for name, value in zip(angle_names, angles, strict=True)
    )

    columns: tuple[Column, ...] = (
        Column("axis", "Specimen axis"),
        Column(
            "direction",
            "Crystal direction",
            help_text="Nearest low-index direction to the specimen axis, in the standard triangle.",
        ),
        Column(
            "residual_deg",
            "Label residual",
            units="deg",
            numeric=True,
            digits=2,
            help_text=(
                "Angle between the true direction and the label; zero means the label is exact."
            ),
        ),
        Column("u", "u", numeric=True, digits=5),
        Column("v", "v", numeric=True, digits=5),
        Column("w", "w", numeric=True, digits=5),
        Column(
            "polar_deg",
            "From c",
            units="deg",
            numeric=True,
            digits=2,
            help_text="Angle between the specimen axis and the crystal c axis.",
        ),
        Column("azimuth_deg", "Azimuth", units="deg", numeric=True, digits=2),
    )

    result = AppResult(
        title=f"{spec.name}: orientation of the view",
        summary=(
            f"From {source}: {angle_text} in the "
            f"{'Bunge' if convention == 'bunge' else 'Matthies'} convention, equivalently a "
            f"rotation of {angle_deg:.2f} deg about "
            f"[{axis_vector[0]:.3f} {axis_vector[1]:.3f} {axis_vector[2]:.3f}] in the specimen "
            "frame. The screen is the specimen frame — RD right, TD up, ND out of the screen — "
            f"so the camera matrix C satisfies v_specimen = C v_crystal. The table gives the "
            f"crystal direction along each specimen axis, folded into the {spec.point_group} "
            "fundamental sector; that is the inverse pole figure, one point per axis."
        ),
        table=ResultTable(
            columns=columns,
            rows=tuple(rows),
            caption=f"Crystal direction along each specimen axis for this view of {spec.name}.",
        ),
        data={
            "camera_matrix": camera,
            "euler": {
                "convention": convention,
                "angles_deg": [float(value) for value in angles],
                "names": list(angle_names),
            },
            "euler_other": {
                "convention": other,
                "angles_deg": [float(value) for value in other_angles],
            },
            # Bunge is stated unconditionally, whichever convention was asked
            # for. It is the convention every EBSD file, every ODF section and
            # every published orientation is written in, so a reader who has
            # switched the picker to Matthies still needs it in front of them —
            # and having it only in whichever of the two slots above happens to
            # hold it makes every consumer decode the pair to find it.
            "euler_bunge": {
                "convention": "bunge",
                "angles_deg": [
                    float(value)
                    for value in (
                        angles if convention == "bunge" else other_angles
                    )
                ],
                "names": ["phi1", "Phi", "phi2"],
            },
            "axis_angle": {
                "axis": [float(value) for value in axis_vector],
                "angle_deg": float(angle_deg),
            },
            "poles": pole_points,
            "ipf_points": ipf_points,
            "projection": method,
            "sector_outline": _sector_outline(
                np.asarray(sector.vertices, dtype=float), method=method
            ),
            "columns": [column.to_json() for column in columns],
        },
        inputs={
            "phase": spec.to_json(),
            "euler_convention": convention,
            "angle1": float(angles[0]),
            "angle2": float(angles[1]),
            "angle3": float(angles[2]),
            "camera_matrix": " ".join(f"{value!r}" for value in camera),
            "poles": [list(row) for row in pole_rows],
            "projection": method,
        },
        notes=(
            "RD is screen right and TD is screen up, so the pole figure is the same view as the "
            "structure beside it. Pole figures elsewhere in the literature are usually drawn "
            "with RD at the top; the ninety-degree difference is deliberate and is stated here "
            "rather than left for the reader to discover.",
            "The crystal direction along each axis is folded into the fundamental sector, so it "
            "is one of a symmetry-equivalent set rather than a unique answer. The label is the "
            "nearest low-index direction and its angular residual is given beside it.",
        ),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    )
    return result.to_json()


def _camera_matrix_values(text: str) -> list[float]:
    """Parse the nine numbers of a camera matrix out of one text field.

    Accepts whitespace, commas, or both, because the viewer sends
    space-separated numbers and a person typing a matrix by hand will not.
    """

    parts = [part for part in re.split(r"[,\s]+", text.strip()) if part]
    try:
        values = [float(part) for part in parts]
    except ValueError as error:
        raise InvalidInputError(
            "The camera matrix must be nine numbers separated by spaces or commas.",
            field="camera_matrix",
            hint="Leave it empty to use the Euler angles instead.",
        ) from error
    if len(values) != 9:
        raise InvalidInputError(
            f"The camera matrix needs nine numbers; got {len(values)}.",
            field="camera_matrix",
            hint="Leave it empty to use the Euler angles instead.",
        )
    return values

_CITATION_KIKUCHI_MAP = (
    "Kikuchi, Japanese Journal of Physics 5 (1928) 83; Williams & Carter, Transmission "
    "Electron Microscopy, 2nd ed., chapter 19, on Kikuchi maps and their use in tilting."
)

#: Points sampled along each band trace before it is sent to the browser.
#:
#: A deliberate compromise. A band edge is a small circle of angular radius near
#: ninety degrees, so its projection is a long curve; the plotting layer uses 721
#: points for a publication figure, which here would be most of a megabyte of
#: JSON for a picture two hundred pixels across. At 181 the arcs are smooth at
#: that size, and the whole map is tens of kilobytes.
_MAP_TRACE_SAMPLES = 181


def _trace_runs(directions: Any) -> list[list[list[float]]]:
    """Projected polylines of a sampled curve, rounded for transport.

    Four decimals on a coordinate that spans [-1, 1] is a ten-thousandth of the
    map: far below a pixel on any figure this will be drawn at, and it halves
    the size of the payload.
    """

    from pytex.diffraction.kikuchi_map import projected_trace_runs

    return [
        [[round(float(point[0]), 4), round(float(point[1]), 4)] for point in run]
        for run in projected_trace_runs(directions)
    ]


@REGISTRY.operation(
    "crystal.kikuchi_map",
    title="Kikuchi map about a zone axis",
    summary="The band network of this phase, on a stereogram centred on a chosen axis.",
    help_text=(
        "The operator's road atlas. Where a simulated pattern shows one orientation on one "
        "detector, this shows the *whole* band network of the crystal on the sphere, projected "
        "stereographically about a chosen zone axis - which is the map a microscopist tilts "
        "by.\n\n"
        "**What is drawn.** Each visible lattice plane contributes a band: a centre line, which "
        "is the great-circle trace of the plane, between two edges, which are the Kossel cones "
        "at the Bragg angle either side of it. A band's angular width is twice its Bragg angle, "
        "so it grows as the spacing falls - the widest bands on a map are its finest-spaced "
        "planes. Where bands cross, a zone axis sits; those are the destinations, and the "
        "number of bands meeting at one is the n-fold symmetry the pattern shows on "
        "arrival.\n\n"
        "**What the centre means.** The projection is centred on the direction given, so that "
        "direction is the one on the beam. Changing it re-centres the same crystal, exactly as "
        "turning to a different standard projection in a textbook does; nothing about the "
        "crystal changes, only which part of its atlas is in view.\n\n"
        "**Limits.** The geometry is exact and the contrast is not modelled at all: excess and "
        "deficient sides, relative darkness and HOLZ lines are dynamical. Band ordering uses "
        "kinematic structure factors, so a phase carrying only a lattice still gives the "
        "correct traces, widths and zone axes - only the ranking is lost, and the result says "
        "so."
    ),
    parameters=(
        phase_parameter(help_text="The phase whose band network to map."),
        IndicesParameter(
            name="centre_direction",
            label="Centre on [uvw]",
            help_text=(
                "The zone axis at the middle of the map - the direction on the beam. The "
                "classical standard projection of a cubic crystal is centred on [001]."
            ),
            default=(0, 0, 1),
        ),
        IndicesParameter(
            name="horizontal_direction",
            label="Along +x",
            help_text=(
                "The direction drawn to the right. Orthogonalized against the centre, so it "
                "need not be exactly perpendicular to it."
            ),
            default=(1, 0, 0),
            advanced=True,
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text=(
                "Enters only through the wavelength, so it scales every band width and changes "
                "nothing about which bands or zone axes exist."
            ),
            units="kV",
            default=200.0,
            minimum=1.0,
            maximum=1000.0,
        ),
        NumberParameter(
            name="max_polar_angle_deg",
            label="Map radius",
            help_text=(
                "How far from the centre the map extends. Sixty degrees covers the useful tilt "
                "range of most holders; ninety is the full hemisphere."
            ),
            units="deg",
            default=60.0,
            minimum=5.0,
            maximum=90.0,
        ),
        IntegerParameter(
            name="max_bands",
            label="Bands drawn",
            help_text=(
                "The strongest this many bands. A map of every band of a crystal is black, and "
                "the point of a map is that it can be navigated."
            ),
            default=14,
            minimum=1,
            maximum=60,
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest |h|, |k| or |l| enumerated for bands.",
            default=3,
            minimum=1,
            maximum=6,
            advanced=True,
        ),
        IntegerParameter(
            name="zone_axis_max_index",
            label="Zone-axis index limit",
            help_text="Largest |u|, |v| or |w| a labelled crossing may carry.",
            default=3,
            minimum=1,
            maximum=6,
            advanced=True,
        ),
        IntegerParameter(
            name="min_zone_axis_order",
            label="Bands to count as a crossing",
            help_text=(
                "How many bands must meet before a crossing is labelled a zone axis. Two is "
                "every intersection; three or more keeps the ones worth tilting to."
            ),
            default=3,
            minimum=2,
            maximum=8,
            advanced=True,
        ),
    ),
    returns=(
        "One row per band with its plane, spacing, width and prominence; the traces, the zone "
        "axes and the map's frame under `data`."
    ),
    panel="crystal",
    citations=(_CITATION_KIKUCHI_MAP,),
    tags=("Kikuchi", "map", "stereogram", "zone axis", "TEM", "navigation", "crystal"),
)
def _kikuchi_map(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.diffraction.kikuchi_map import compute_kikuchi_map

    spec, phase = phase_from_request(request["phase"])
    centre = tuple(int(value) for value in request["centre_direction"])
    horizontal = tuple(int(value) for value in request["horizontal_direction"])
    max_polar = float(request["max_polar_angle_deg"])
    try:
        kikuchi_map = compute_kikuchi_map(
            phase,
            beam_energy_kev=float(request["beam_energy_kev"]),
            centre_direction=centre,
            horizontal_direction=horizontal,
            max_index=int(request["max_index"]),
            max_bands=int(request["max_bands"]),
            zone_axis_max_index=int(request["zone_axis_max_index"]),
            min_zone_axis_order=int(request["min_zone_axis_order"]),
            max_polar_angle_deg=max_polar,
        )
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            field="centre_direction",
            hint=(
                "The centre and the +x direction must both be non-zero and must not be "
                "parallel: the map frame is built from the two of them."
            ),
        ) from error

    if not kikuchi_map.bands:
        raise InvalidInputError(
            "No band of this phase survives the cut-offs, so there is no map to draw.",
            field="max_index",
            hint="Raise the index limit, or lower the number of bands required at a crossing.",
        )

    # The traces are sampled and split here rather than in the browser: the same
    # projection, the same equator splitting and the same folding as every other
    # figure in the application, done once in the place that owns them.
    bands = []
    for band in kikuchi_map.bands:
        indices = tuple(int(value) for value in band.indices)
        narrow, far = band.edge_directions(samples=_MAP_TRACE_SAMPLES)
        bands.append(
            {
                "hkl": list(indices),
                "label": plane_label(indices, spec=spec),
                "d_angstrom": float(band.d_spacing_angstrom),
                "width_deg": float(band.angular_width_deg),
                "bragg_angle_deg": float(band.bragg_angle_deg),
                "intensity": float(band.relative_intensity),
                "multiplicity": int(band.family_multiplicity),
                "centre": _trace_runs(band.centre_directions(samples=_MAP_TRACE_SAMPLES)),
                "edges": [_trace_runs(narrow), _trace_runs(far)],
            }
        )

    axes = []
    for axis in kikuchi_map.zone_axes:
        indices = tuple(int(value) for value in axis.indices)
        point = project_directions(
            np.asarray(axis.direction_map, dtype=float).reshape(1, 3),
            method="stereographic",
            antipodal=True,
        )[0]
        axes.append(
            {
                "uvw": list(indices),
                "label": direction_label(indices, spec=spec),
                "x": round(float(point[0]), 5),
                "y": round(float(point[1]), 5),
                "order": int(axis.order),
                "polar_angle_deg": float(axis.polar_angle_deg),
            }
        )

    centre_text = direction_label(centre, spec=spec)
    rows = [
        {
            "plane": band["label"],
            "d": band["d_angstrom"],
            "width": band["width_deg"],
            "intensity": band["intensity"],
            "multiplicity": band["multiplicity"],
        }
        for band in bands
    ]
    result = AppResult(
        title=f"{spec.name}: Kikuchi map about {centre_text}",
        summary=(
            f"{len(bands)} band(s) and {len(axes)} labelled zone axis / axes within "
            f"{max_polar:g} deg of {centre_text}, at "
            f"{float(request['beam_energy_kev']):g} kV. A band's angular width is twice its "
            "Bragg angle, so the widest bands belong to the finest-spaced planes; where bands "
            "cross, a zone axis sits, and the number meeting there is the symmetry of the "
            "pattern seen on arriving."
            + (
                ""
                if kikuchi_map.has_intensity_model
                else " This phase carries no atomic basis, so the geometry is exact but the "
                "band ordering is not: every reflection is treated as equally strong."
            )
        ),
        table=ResultTable(
            columns=(
                Column("plane", "Band"),
                Column("d", "d", units="A", numeric=True, digits=4),
                Column(
                    "width",
                    "Angular width",
                    units="deg",
                    numeric=True,
                    digits=3,
                    help_text="Twice the Bragg angle. It grows as the spacing falls.",
                ),
                Column(
                    "intensity",
                    "Prominence",
                    numeric=True,
                    digits=3,
                    help_text=(
                        "Kinematic, relative to the strongest band. Indicative only: band "
                        "contrast is dynamical and is not modelled."
                    ),
                ),
                Column("multiplicity", "Family size", numeric=True),
            ),
            rows=tuple(rows),
            caption=f"Kikuchi bands of {spec.name} in the region mapped about {centre_text}.",
        ),
        data={
            "bands": bands,
            "zone_axes": axes,
            "centre": list(centre),
            "centre_label": centre_text,
            "horizontal": list(horizontal),
            "beam_energy_kev": float(request["beam_energy_kev"]),
            "wavelength_angstrom": float(kikuchi_map.wavelength_angstrom),
            "max_polar_angle_deg": max_polar,
            # Where the map's rim falls in projected coordinates, so the browser
            # can draw the boundary without repeating the projection: the
            # stereographic radius of a direction at the polar angle is
            # tan(rho / 2).
            "boundary_radius": round(float(np.tan(np.radians(max_polar) / 2.0)), 6),
            "has_intensity_model": bool(kikuchi_map.has_intensity_model),
            # Crystal Cartesian to map frame, row-major. With it a caller can put
            # any crystal direction on this map -- in particular the direction a
            # viewer currently has on the beam, which is what makes the map
            # something to navigate by rather than a static atlas.
            "view_matrix": [
                float(value)
                for value in np.asarray(kikuchi_map.view_matrix, dtype=float).reshape(-1)
            ],
            "projection": "stereographic",
            "describe": kikuchi_map.describe(),
        },
        inputs={
            "phase": spec.to_json(),
            "centre_direction": list(centre),
            "horizontal_direction": list(horizontal),
            "beam_energy_kev": float(request["beam_energy_kev"]),
            "max_polar_angle_deg": max_polar,
            "max_bands": int(request["max_bands"]),
            "max_index": int(request["max_index"]),
        },
        notes=(
            "Geometry only: band positions and widths are exact, and contrast - excess and "
            "deficient sides, relative darkness, HOLZ lines - is dynamical and not modelled.",
            "The map is of the crystal rather than of a specimen: it says where the bands lie "
            "relative to each other, and a measured pattern is matched to it by finding the "
            "orientation that brings the two into register.",
        ),
        citations=(_CITATION_KIKUCHI_MAP,),
    )
    return result.to_json()


@REGISTRY.operation(
    "crystal.render",
    title="Publication figure of the structure",
    summary="Render the current view through the publication renderer.",
    help_text=(
        "Takes the camera you set by dragging and renders the same scene through "
        "`pytex.plotting.crystal3d` in the journal style, with lit spheres and correct "
        "depth ordering.\n\n"
        "The browser view and this figure are the same geometry seen from the same direction, "
        "because the camera angles travel with the request. What you framed on screen is what "
        "the figure shows.\n\n"
        "**Choosing a format.** A lit sphere is a mesh, and a vector format must write every "
        "facet of it: an SVG of a two-cell rock-salt block runs to megabytes and will defeat "
        "most journal submission systems. PNG at 600 dpi is therefore the default, and is the "
        "better publication artefact for ball-and-stick and space-filling views. Choose SVG "
        "when the figure must stay editable, or for the stick and wireframe styles, which are "
        "line art and stay small; the sphere mesh is coarsened automatically for SVG so the "
        "file remains openable."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="render_style",
            label="Style",
            help_text="The visual system for atoms and bonds.",
            options=_RENDER_STYLES,
            default="ball_and_stick",
        ),
        IntegerParameter(
            name="repeat_a",
            label="Repeats along a",
            help_text="Cells along a.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="a",
            row="Repeats",
        ),
        IntegerParameter(
            name="repeat_b",
            label="Repeats along b",
            help_text="Cells along b.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="b",
            row="Repeats",
        ),
        IntegerParameter(
            name="repeat_c",
            label="Repeats along c",
            help_text="Cells along c.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
            symbol="c",
            row="Repeats",
        ),
        IndicesListParameter(
            name="planes",
            label="Planes (hkl)",
            help_text="Planes to superimpose, one per row.",
            required=False,
            group="Overlays",
        ),
        IndicesListParameter(
            name="directions",
            label="Directions [uvw]",
            help_text="Directions to superimpose, one per row.",
            required=False,
            group="Overlays",
        ),
        NumberParameter(
            name="elevation_deg",
            label="Camera elevation",
            help_text="Angle above the horizon, matching the on-screen view.",
            units="°",
            default=22.0,
            minimum=-90.0,
            maximum=90.0,
            group="Camera",
            row="Camera",
        ),
        NumberParameter(
            name="azimuth_deg",
            label="Camera azimuth",
            help_text="Rotation about the vertical, matching the on-screen view.",
            units="°",
            default=34.0,
            group="Camera",
            row="Camera",
            field_width="short",
        ),
        BooleanParameter(
            name="show_bonds",
            label="Draw bonds",
            help_text="Include inferred bonds.",
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_unit_cells",
            label="Outline every cell",
            help_text="Draw the edges of each repeated cell, not only the outer box.",
            default=False,
            group="Overlays",
        ),
        BooleanParameter(
            name="hexagonal_prism",
            label="Draw hexagonal phases as the prism",
            help_text=(
                "Draw a hexagonal phase as its hexagonal prism — three cells about one "
                "atomic column — rather than as the 120-degree rhombus the lattice is "
                "written on. The prism is the figure the sixfold symmetry is visible in; the "
                "rhombus shows a third of it. Plane overlays follow, so a basal plane is "
                "clipped to the hexagon. Ignored for phases that are not on hexagonal axes."
            ),
            default=True,
            group="Extent",
        ),
        ChoiceParameter(
            name="atom_labels",
            label="Atom labels",
            help_text="Label every atom with its element or its site name.",
            options=(
                ("none", "None", "No labels; cleanest for dense structures."),
                ("species", "Element", "The element symbol on each atom."),
                ("site", "Site", "The crystallographic site label, as in the CIF."),
            ),
            default="none",
            advanced=True,
            group="Overlays",
        ),
        NumberParameter(
            name="bond_tolerance_angstrom",
            label="Bond tolerance",
            help_text=(
                "Added to the sum of covalent radii when deciding whether two atoms are bonded. "
                "Must match the viewer's value, or the figure will show a different set of "
                "bonds from the one on screen."
            ),
            units="Å",
            default=0.45,
            minimum=0.0,
            maximum=2.0,
            advanced=True,
            group="Overlays",
        ),
        ObjectParameter(
            name="appearance",
            label="Object properties",
            help_text=(
                "Presentation-only atom, bond, cell, plane, direction and annotation settings. "
                "The interactive viewer supplies this object automatically when publishing."
            ),
            editor="json",
            required=False,
            default={},
            advanced=True,
            group="Appearance",
        ),
        BooleanParameter(
            name="show_legend",
            label="Species legend",
            help_text="Add a colour key for the elements present.",
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_frame_indicator",
            label="Axis gizmo",
            help_text=(
                "Add a small indicator showing where a, b and c point in the rendered view, so "
                "the figure states its own orientation."
            ),
            default=True,
            group="Overlays",
        ),
        TextParameter(
            name="camera_matrix",
            label="Camera matrix",
            help_text=(
                "Nine numbers, row-major, describing the on-screen camera. Sent automatically by "
                "the viewer so the figure shows exactly the view you framed; it overrides the "
                "elevation and azimuth above when present."
            ),
            required=False,
            advanced=True,
            group="Camera",
        ),
        ChoiceParameter(
            name="format",
            label="Format",
            help_text="PNG for lit spheres; SVG when the figure must stay editable.",
            options=(
                ("png", "PNG", "Raster at the chosen resolution. Compact, and exact for spheres."),
                ("svg", "SVG", "Vector and editable. Large for sphere-based styles."),
            ),
            default="png",
            group="Output",
        ),
        IntegerParameter(
            name="dpi",
            label="Resolution",
            help_text=(
                "Dots per inch for PNG. 300 is the usual journal minimum; 600 is safe for "
                "anything. Ignored for SVG."
            ),
            units="dpi",
            default=600,
            minimum=72,
            maximum=1200,
            group="Output",
        ),
    ),
    returns="The image under `data.image`: markup for SVG, base64 for PNG.",
    panel="crystal",
    citations=(_CITATION_VESTA,),
    tags=("export", "figure", "SVG", "publication", "crystal"),
)
def _crystal_render(request: dict[str, Any]) -> dict[str, Any]:
    import matplotlib

    # A server process has no display, and the renderer is called from request
    # handling. `force=False` leaves an already-chosen backend alone.
    matplotlib.use("Agg", force=False)
    import base64
    import io

    import matplotlib.pyplot as plt

    from pytex.plotting.crystal3d import plot_crystal_structure_3d

    spec, phase = phase_from_request(request["phase"])
    repeats = (int(request["repeat_a"]), int(request["repeat_b"]), int(request["repeat_c"]))
    plane_rows = tuple(request.get("planes") or ())
    direction_rows = tuple(request.get("directions") or ())
    appearance = _appearance(request.get("appearance"))
    image_format = str(request["format"])
    dpi = int(request["dpi"])
    elevation = float(request["elevation_deg"])
    azimuth = float(request["azimuth_deg"])
    if request.get("camera_matrix"):
        elevation, azimuth = camera_angles_from_matrix(
            [float(value) for value in str(request["camera_matrix"]).replace(",", " ").split()]
        )

    # A sphere is drawn as a quad mesh, and a vector format writes every facet:
    # at the renderer's default resolution a two-cell block becomes tens of
    # thousands of polygons and an SVG of several megabytes. Coarsening the mesh
    # for SVG keeps the file openable, at a cost invisible at figure size. PNG
    # keeps the full mesh, because there the cost is a fixed pixel count either
    # way.
    style_overrides = _appearance_style(appearance, str(request["render_style"]))
    if image_format == "svg":
        style_overrides["crystal"].update(
            {"atom_surface_resolution": 18, "bond_surface_resolution": 16}
        )

    axes = plot_crystal_structure_3d(
        phase,
        repeats=repeats,
        render_style=str(request["render_style"]),
        show_bonds=bool(request["show_bonds"]) and appearance["show_bonds"],
        bond_tolerance_angstrom=float(request["bond_tolerance_angstrom"]),
        show_unit_cells=bool(request["show_unit_cells"]) and appearance["show_cells"],
        hexagonal_prism=bool(request["hexagonal_prism"]),
        atom_label_mode=str(request["atom_labels"]) if appearance["show_labels"] else "none",
        plane_overlays=(
            tuple(_plane_overlay(row, phase) for row in plane_rows)
            if appearance["show_planes"]
            else ()
        ),
        direction_overlays=(
            tuple(_direction_overlay(row, phase) for row in direction_rows)
            if appearance["show_directions"]
            else ()
        ),
        elev_deg=elevation,
        azim_deg=azimuth,
        show_legend=bool(request["show_legend"]),
        show_frame_indicator=bool(request["show_frame_indicator"]) and appearance["show_gizmo"],
        style_overrides=style_overrides,
    )
    figure = axes.get_figure()
    buffer = io.BytesIO()
    try:
        figure.savefig(buffer, format=image_format, dpi=dpi, bbox_inches="tight")
    finally:
        # Tests treat a leaked figure as a defect, and a server that leaked one
        # per export would exhaust memory in an afternoon.
        plt.close(figure)
    payload = buffer.getvalue()
    image = (
        payload.decode("utf-8")
        if image_format == "svg"
        else base64.b64encode(payload).decode("ascii")
    )

    result = AppResult(
        title=f"{spec.name}: publication figure",
        summary=(
            f"{image_format.upper()} figure of {spec.name} at elevation {elevation:.1f}° and "
            f"azimuth {azimuth:.1f}°, rendered through the journal style"
            + (f" at {dpi} dpi" if image_format == "png" else "")
            + f" ({len(payload) / 1024:.0f} kB). The geometry is identical to the interactive "
            "view; only the renderer differs."
        ),
        notes=(
            (
                "The sphere mesh is coarsened for SVG so the file stays openable. For a "
                "sphere-based style, PNG at 600 dpi is both smaller and closer to what the "
                "renderer intends.",
            )
            if image_format == "svg"
            else ()
        ),
        data={
            "image": image,
            "format": image_format,
            "encoding": "text" if image_format == "svg" else "base64",
            "bytes": len(payload),
            "dpi": dpi,
            "appearance": appearance,
            "elevation_deg": elevation,
            "azimuth_deg": azimuth,
        },
        inputs={
            "phase": spec.to_json(),
            "render_style": request["render_style"],
            "repeat_a": repeats[0],
            "repeat_b": repeats[1],
            "repeat_c": repeats[2],
            "planes": [list(row) for row in plane_rows],
            "directions": [list(row) for row in direction_rows],
            "elevation_deg": elevation,
            "azimuth_deg": azimuth,
            "format": image_format,
            "dpi": dpi,
        },
        citations=(_CITATION_VESTA,),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Canonical examples
# --------------------------------------------------------------------------

REGISTRY.add_examples(
    (
        ExampleScenario(
            id="crystal.example.nacl",
            title="Rock salt: two interpenetrating fcc lattices",
            panel="crystal",
            summary="NaCl with the (100) cube face and the (111) close-packed plane drawn in.",
            teaches=(
                "Rotate until you look down [100]: sodium and chlorine alternate along every "
                "cube edge. The (111) plane, by contrast, cuts through atoms of only one species "
                "at a time — which is why the 111 reflection of NaCl is weak while 200 is strong."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "nacl"},
                "planes": [[1, 0, 0], [1, 1, 1]],
                "atom_labels": "species",
            },
        ),
        ExampleScenario(
            id="crystal.example.fcc_slip",
            title="The fcc slip system, drawn",
            panel="crystal",
            summary="Austenite with (111) and the [1-10] slip direction superimposed.",
            teaches=(
                "The arrow lies *in* the polygon — that is what h·u + k·v + l·w = 0 looks like. "
                "Turn the structure until the (111) plane is edge-on and the atoms in it line up "
                "into the close-packed rows the dislocation actually moves along."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "planes": [[1, 1, 1]],
                "directions": [[1, -1, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
        ExampleScenario(
            id="crystal.example.bcc",
            title="Ferrite and its {110} planes",
            panel="crystal",
            summary="bcc iron with two members of the {110} family drawn in.",
            teaches=(
                "The body-centring atom sits in every (110) plane, not between them — which is "
                "why {110} is the most densely packed plane of a bcc metal and carries most of "
                "its slip, despite bcc having no truly close-packed plane at all."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "fe_bcc"},
                "planes": [[1, 1, 0], [1, -1, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
        ExampleScenario(
            id="crystal.example.zr_basal",
            title="Zirconium: basal and prism planes together",
            panel="crystal",
            summary="hcp zirconium with (0001) and a prism plane, plus the a-direction.",
            teaches=(
                "The basal plane is close-packed and the prism plane is not, and the picture "
                "shows why: turn it so the basal plane is edge-on and the ABAB stacking is "
                "visible as two offset layers. That stacking, not symmetry alone, is what makes "
                "prism slip compete with basal slip in zirconium."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "zr_hcp"},
                "planes": [[0, 0, 1], [1, 0, 0]],
                "directions": [[1, 0, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
    )
)
