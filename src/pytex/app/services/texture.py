"""Texture: pole figures, inverse pole figures, and ODF sections of a model texture.

What it does
    Builds a polycrystal from named texture components, then shows it the three
    ways texture is read: as a **pole figure** (where a crystal plane points, in
    specimen axes), as an **inverse pole figure** (which crystal direction lies
    along a specimen axis), and as **ODF sections** (the orientation density
    itself, sliced).

Why a model texture and not a data file
    The application has no measurement to load, and a texture panel that can
    only say "import an EBSD map first" teaches nothing. A model built from the
    components the literature names — cube, Goss, brass, copper, S — is the
    thing every textbook figure is of, and it has the decisive property that a
    file does not: the answer is known before the calculation runs. A random
    texture must give 1 m.r.d. everywhere; a single sharp component must put its
    poles exactly where its Miller label says. Those are the checks in
    `test_app_texture.py`, and neither is available from a data set.

The one number that matters
    Everything here is reported in **multiples of a random distribution**. An
    intensity of 4 m.r.d. means four times as many poles point that way as would
    in a texture-free material, and the scale is what makes two figures from
    different instruments comparable at all. It is not a convenience: an
    unnormalised pole figure is a picture, not a measurement. See
    `docs/site/theory/pole_figure_arithmetic_and_mrd.md`.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.logbook import APP_LOG
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

__all__: tuple[str, ...] = ()

_CITATION_BUNGE = "Bunge, Texture Analysis in Materials Science (1982), chapters 2 and 4."
_CITATION_RANDLE = "Randle & Engler, Introduction to Texture Analysis, 2nd ed., chapters 2 and 5."
_CITATION_HIRSCH = "Hirsch & Lucke, Acta Metall. 36 (1988) 2863 (fcc rolling texture components)."

#: The specimen frame every texture in this panel is expressed in.
#:
#: Named RD/TD/ND rather than x/y/z because that is what a rolling texture is
#: read in, and because a pole figure with unlabelled axes is unreadable: the
#: whole content of "brass" versus "copper" is where the poles sit relative to
#: the rolling direction.
_SPECIMEN_AXES = ("RD", "TD", "ND")


def _specimen_frame() -> Any:
    from pytex.core import FrameDomain, Handedness, ReferenceFrame

    return ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=_SPECIMEN_AXES,
        handedness=Handedness.RIGHT,
    )


#: The texture models the panel offers, as choice options.
_TEXTURE_MODELS: tuple[tuple[str, str, str], ...] = (
    (
        "fcc_rolling",
        "fcc rolling texture (cube, Goss, brass, copper, S)",
        "The five components of a rolled fcc metal, in equal proportions.",
    ),
    (
        "bcc_rolling",
        "bcc rolling texture (rotated cube, Goss, rotated Goss)",
        "The components of a rolled bcc metal.",
    ),
    ("cube", "Cube {001}<100>", "The recrystallisation component of fcc metals."),
    ("goss", "Goss {011}<100>", "The component of grain-oriented electrical steel."),
    ("brass", "Brass {011}<211>", "The dominant component of low stacking-fault fcc."),
    ("copper", "Copper {112}<111>", "The dominant component of high stacking-fault fcc."),
    ("s", "S {123}<634>", "The third fcc rolling component."),
    (
        "random",
        "Random",
        "No texture at all: the 1 m.r.d. baseline everything is measured against.",
    ),
)

#: Component identifiers to their catalogue entries, resolved lazily.
_SINGLE_COMPONENTS = {
    "cube": "CUBE",
    "goss": "GOSS",
    "brass": "BRASS",
    "copper": "COPPER",
    "s": "S_COMPONENT",
    "rotated_cube": "ROTATED_CUBE",
    "rotated_goss": "ROTATED_GOSS",
}


def _components(model: str) -> tuple[Any, ...]:
    """Return the catalogue components a model is built from."""

    from pytex.texture import components as catalogue

    if model == "fcc_rolling":
        return tuple(catalogue.STANDARD_FCC_ROLLING_COMPONENTS)
    if model == "bcc_rolling":
        return tuple(catalogue.STANDARD_BCC_ROLLING_COMPONENTS)
    if model == "random":
        return ()
    name = _SINGLE_COMPONENTS.get(model)
    if name is None:  # pragma: no cover - the choice list is closed
        raise InvalidInputError(f"Unknown texture model {model!r}.", field="model")
    return (getattr(catalogue, name),)


def _model_parameters() -> tuple[Any, ...]:
    """The controls that define the texture, shared by all three views."""

    return (
        phase_parameter(
            help_text=(
                "The crystal the texture belongs to. Its symmetry is what folds the pole "
                "figure — a cubic texture repeats every 90° where a hexagonal one does not."
            ),
            builtin="ni_fcc",
        ),
        ChoiceParameter(
            name="model",
            label="Texture",
            help_text=(
                "Which model texture to build. The rolling textures are mixtures of the named "
                "components in equal proportions; the single components are what those mixtures "
                "are made of; random is the baseline."
            ),
            options=_TEXTURE_MODELS,
            default="fcc_rolling",
        ),
        NumberParameter(
            name="spread_deg",
            label="Component spread",
            help_text=(
                "How far grains scatter about each ideal orientation, as a standard deviation "
                "in each Euler angle. Real textures are spread: 0° gives the textbook ideal "
                "orientation and an unphysically sharp figure, 10-15° is a realistic rolled "
                "metal, and beyond about 30° the components merge into each other."
            ),
            units="°",
            default=10.0,
            minimum=0.0,
            maximum=40.0,
        ),
        IntegerParameter(
            name="grain_count",
            label="Grains per component",
            help_text=(
                "How many orientations are drawn around each component. More is smoother and "
                "slower; the shape of the figure is settled by a few hundred."
            ),
            default=200,
            minimum=10,
            maximum=2000,
        ),
        NumberParameter(
            name="halfwidth_deg",
            label="Kernel halfwidth",
            help_text=(
                "The width of the bell placed on each orientation before the density is "
                "evaluated. This is a *smoothing* choice, not a property of the material: too "
                "small and the figure shows the individual grains, too large and real detail is "
                "washed out. It is the single setting most often reported wrongly, so it is here "
                "rather than buried."
            ),
            units="°",
            default=10.0,
            minimum=2.0,
            maximum=30.0,
            advanced=True,
        ),
        IntegerParameter(
            name="seed",
            label="Random seed",
            help_text=(
                "Fixes the scatter, so the same settings give the same figure. Change it to see "
                "how much of a feature is the texture and how much is the sample of grains."
            ),
            default=7,
            minimum=0,
            maximum=1_000_000,
            advanced=True,
            field_width="short",
        ),
    )


def _build_texture(request: dict[str, Any]) -> tuple[Any, Any, PhaseSpec, str]:
    """Build the model orientation set and its ODF.

    Returns the orientation set, the ODF, the phase specification, and the
    display name of the model.
    """

    from pytex.core.orientation import OrientationSet
    from pytex.texture.models import ODF, KernelSpec

    spec, phase = phase_from_request(request["phase"])
    model = str(request["model"])
    spread = float(request["spread_deg"])
    per_component = int(request["grain_count"])
    generator = np.random.default_rng(int(request["seed"]))
    specimen = _specimen_frame()

    components = _components(model)
    if components:
        # Scatter is applied to the Euler angles rather than as a rotation about
        # a random axis. That is the cruder of the two and it is the honest one
        # here: this is a teaching model, not a measurement, and stating the
        # spread in the same coordinates the component is quoted in keeps the
        # control's meaning obvious.
        blocks = [
            np.asarray(component.bunge_euler_deg, dtype=float)
            + generator.normal(0.0, spread, size=(per_component, 3))
            for component in components
        ]
        angles = np.concatenate(blocks, axis=0)
    else:
        # A random texture must be uniform on SO(3), which uniform Euler angles
        # are *not*: they crowd the poles. Sampling Phi through its sine measure
        # is what makes the resulting pole figure come out flat at 1 m.r.d.,
        # which is the property the whole panel is calibrated against.
        count = per_component * 5
        phi1 = generator.uniform(0.0, 360.0, size=count)
        big_phi = np.degrees(np.arccos(generator.uniform(-1.0, 1.0, size=count)))
        phi2 = generator.uniform(0.0, 360.0, size=count)
        angles = np.column_stack([phi1, big_phi, phi2])

    orientations = OrientationSet.from_euler_angles(
        angles,
        specimen_frame=specimen,
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF.from_orientations(
        orientations, kernel=KernelSpec(halfwidth_deg=float(request["halfwidth_deg"]))
    )
    label = next(
        (title for key, title, _ in _TEXTURE_MODELS if key == model),
        model,
    )
    return orientations, odf, spec, label


def _model_inputs(request: dict[str, Any], spec: PhaseSpec) -> dict[str, Any]:
    """The texture-defining inputs, echoed so a result can be reproduced."""

    return {
        "phase": spec.to_json(),
        "model": str(request["model"]),
        "spread_deg": float(request["spread_deg"]),
        "grain_count": int(request["grain_count"]),
        "halfwidth_deg": float(request["halfwidth_deg"]),
        "seed": int(request["seed"]),
    }


def _odf_random_reference(phase: Any) -> float:
    r"""What a random texture reads as from ``ODF.evaluate(normalized=True)``.

    Not 1, which is the trap this function exists to close. The kernel is
    normalised to integrate to one over the whole of SO(3), but a
    symmetry-aware evaluation folds every query into the fundamental zone,
    which is :math:`1/|G|` of SO(3) for a proper rotation group :math:`G`. A
    uniform distribution therefore reads :math:`|G|`, not 1 — and dividing by
    :math:`|G|` is what puts an ODF section on the same multiples-of-random
    scale as the pole figures beside it.

    Measured before it was derived: a random texture reads 23.9 for m-3m and
    11.9 for 6/mmm, against operator counts of 24 and 12.
    ``test_app_texture.py::test_a_random_texture_reads_one_mrd_in_every_view``
    checks both, because a factor of 24 in a quantity labelled "m.r.d." is the
    kind of error that survives every test that only compares a figure with
    itself.
    """

    symmetry = getattr(phase, "symmetry", None)
    if symmetry is None:
        return 1.0
    return float(len(np.asarray(symmetry.operators)))


def _crystal_plane(phase: Any, indices: tuple[int, ...]) -> Any:
    from pytex.core.lattice import CrystalPlane, MillerIndex

    return CrystalPlane(
        miller=MillerIndex(np.asarray(indices, dtype=int), phase=phase), phase=phase
    )


def _project(vectors: np.ndarray, method: str) -> np.ndarray:
    """Project unit vectors onto the unit disc, rim at radius 1 in both methods.

    The same normalisation the variants panel applies, and for the same reason:
    the library returns equal-area at its natural radius of √2, and two figures
    that are meant to be compared must share a rim.
    """

    from pytex.texture.projections import project_directions

    projected = np.asarray(project_directions(vectors, method=method), dtype=float)
    return projected / (math.sqrt(2.0) if method == "equal_area" else 1.0)


_PROJECTION_PARAMETER = ChoiceParameter(
    name="projection",
    label="Projection",
    help_text=(
        "How the hemisphere is flattened. Equal-area is the right default for a pole figure, "
        "because the figure is about how much intensity sits where, and a stereographic "
        "projection exaggerates area toward the rim."
    ),
    options=(
        ("equal_area", "Equal area (Schmidt)", "Preserves area; the pole-figure convention."),
        ("stereographic", "Stereographic (Wulff)", "Preserves angles; what a Wulff net measures."),
    ),
    default="equal_area",
)

#: Columns of the pole-figure grid, shared by the hover card and every export.
_DENSITY_COLUMNS: tuple[Column, ...] = (
    Column(
        "mrd",
        "Intensity",
        units="m.r.d.",
        numeric=True,
        digits=4,
        help_text=(
            "Multiples of a random distribution. 1 is what a texture-free material gives "
            "everywhere; 4 means four times as many poles point this way."
        ),
    ),
    Column("x", "x", numeric=True, digits=5, help_text="Position on the projection."),
    Column("y", "y", numeric=True, digits=5),
    Column(
        "polar_deg",
        "Polar angle",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle from ND, the centre of the figure.",
    ),
    Column(
        "azimuth_deg",
        "Azimuth",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle from RD, measured toward TD.",
    ),
)


@REGISTRY.operation(
    "texture.pole_figure",
    title="Pole figure",
    summary="Where a crystal plane points, over a whole polycrystal, in specimen axes.",
    help_text=(
        "Builds a model texture from named components and evaluates the density of a chosen "
        "crystal plane's poles over the specimen hemisphere. This is the figure texture is read "
        "in: the centre is ND, the rim is the rolling plane, and RD is to the right.\n\n"
        "**Everything is in multiples of a random distribution.** An intensity of 4 m.r.d. means "
        "four times as many poles point that way as would in a texture-free material. The scale "
        "is what makes two figures comparable at all, and it has an exact consequence worth "
        "knowing: the area-weighted mean over the hemisphere is 1 m.r.d. by construction, "
        "whatever the texture. A figure whose mean is not 1 has not been normalised, and its "
        "numbers mean nothing outside itself.\n\n"
        "**Start with random.** It gives a flat figure at 1 m.r.d. everywhere, which is both the "
        "baseline and the check that the normalisation is right. Then switch to a single "
        "component and watch the poles land exactly where its Miller label says they should — "
        "Goss is {011}<100>, so its (011) poles sit at the centre of the figure.\n\n"
        "**The kernel halfwidth is a choice, not a measurement.** It sets how much the "
        "individual grains are smoothed before the density is read. Too small and the figure "
        "shows the sample rather than the texture; too large and real detail is washed out. It "
        "is the setting most often left unreported in the literature, which is why it is a "
        "control here rather than a constant."
    ),
    parameters=(
        *_model_parameters(),
        IndicesParameter(
            name="pole",
            label="Plane to plot {hkl}",
            help_text=(
                "The crystal plane whose poles are counted. The whole symmetry family is "
                "included, which is what a measured pole figure contains too. (111) and (200) "
                "are the usual pair for an fcc metal."
            ),
            default=(1, 1, 1),
        ),
        _PROJECTION_PARAMETER,
        NumberParameter(
            name="resolution_deg",
            label="Grid resolution",
            help_text=(
                "Spacing of the evaluation grid on the hemisphere. Finer is smoother and "
                "slower, and does not add information beyond the kernel halfwidth."
            ),
            units="°",
            default=5.0,
            minimum=2.0,
            maximum=15.0,
            advanced=True,
        ),
    ),
    returns="One row per grid point, with intensity in m.r.d.; the texture summary under `data`.",
    panel="texture",
    citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_HIRSCH),
    tags=("texture", "pole figure", "m.r.d.", "ODF", "rolling", "component", "specimen"),
)
def _pole_figure(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.sphere import S2Grid
    from pytex.texture.models import KernelSpec, random_pole_density

    orientations, odf, spec, model_label = _build_texture(request)
    _spec, phase = phase_from_request(request["phase"])
    pole = _crystal_plane(phase, tuple(request["pole"]))
    method = str(request["projection"])

    grid = S2Grid.equispaced(
        float(request["resolution_deg"]), reference_frame=_specimen_frame(), hemisphere="upper"
    )
    directions = np.asarray(grid.vectors, dtype=float)
    weights = np.asarray(grid.weights, dtype=float)

    # `evaluate_pole_density` returns a density in the kernel's own units. The
    # m.r.d. scale divides by what a *random* texture gives under the same
    # kernel, which is a closed form rather than a second simulation — so the
    # normalisation carries no sampling noise of its own.
    raw = np.asarray(odf.evaluate_pole_density(pole, directions), dtype=float)
    kernel = KernelSpec(halfwidth_deg=float(request["halfwidth_deg"]))
    reference = float(random_pole_density(kernel))
    mrd = raw / reference

    projected = _project(directions, method)
    polar_deg = np.degrees(np.arccos(np.clip(directions[:, 2], -1.0, 1.0)))
    azimuth_deg = np.degrees(np.arctan2(directions[:, 1], directions[:, 0])) % 360.0

    rows = [
        {
            "mrd": float(mrd[index]),
            "x": float(projected[index, 0]),
            "y": float(projected[index, 1]),
            "polar_deg": float(polar_deg[index]),
            "azimuth_deg": float(azimuth_deg[index]),
        }
        for index in range(directions.shape[0])
    ]

    mean_mrd = float(np.average(mrd, weights=weights))
    peak = int(np.argmax(mrd))
    pole_text = plane_label(tuple(int(value) for value in request["pole"]), spec=spec)
    fractions = _volume_fractions(orientations, str(request["model"]))

    result = AppResult(
        title=f"{pole_text} pole figure: {model_label} in {spec.name}",
        summary=(
            f"{pole_text} poles of {len(orientations.quaternions)} grains, on a "
            f"{float(request['resolution_deg']):g}° hemisphere grid smoothed with a "
            f"{float(request['halfwidth_deg']):g}° kernel. The texture reaches "
            f"{float(mrd.max()):.2f} m.r.d. at {polar_deg[peak]:.0f}° from ND, and falls to "
            f"{float(mrd.min()):.2f} m.r.d. at its weakest. The area-weighted mean is "
            f"{mean_mrd:.3f} m.r.d., which is 1 by construction and is the check that the "
            "figure is normalised."
        ),
        table=ResultTable(
            columns=_DENSITY_COLUMNS,
            rows=tuple(rows),
            caption=f"{pole_text} pole density over the specimen hemisphere, in m.r.d.",
        ),
        data={
            "points": rows,
            "projection": method,
            "max_mrd": float(mrd.max()),
            "min_mrd": float(mrd.min()),
            "mean_mrd": mean_mrd,
            "pole_label": pole_text,
            "model_label": model_label,
            "grain_count": len(orientations.quaternions),
            "specimen_axes": list(_SPECIMEN_AXES),
            "component_fractions": fractions,
            "columns": [column.to_json() for column in _DENSITY_COLUMNS],
        },
        inputs=dict(
            _model_inputs(request, spec),
            pole=[int(value) for value in request["pole"]],
            projection=method,
            resolution_deg=float(request["resolution_deg"]),
        ),
        notes=(
            "The area-weighted mean of any correctly normalised pole figure is exactly 1 m.r.d. "
            "An unweighted mean over a polar raster is not: it over-counts the crowded grid near "
            "the centre. See the m.r.d. theory note.",
            "This is a model texture built from named components, not a measurement. Its value "
            "is that the answer is known in advance — random gives 1 m.r.d. everywhere, and a "
            "single component puts its poles where its Miller label says.",
        ),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_HIRSCH),
    )
    return result.to_json()


def _volume_fractions(orientations: Any, model: str) -> list[dict[str, Any]]:
    """Fraction of grains within 15° of each named component of this model.

    Reported beside every figure because it is the quantitative reading of what
    the figure shows: "the brass component is 30% of this texture" is a claim a
    pole figure supports only qualitatively.
    """

    from pytex.texture.components import component_volume_fractions

    components = _components(model)
    if not components:
        return []
    fractions = component_volume_fractions(orientations, list(components), tolerance_deg=15.0)
    return [
        {
            "component": component.name,
            "miller": component.miller_label,
            "fraction": float(fractions.get(component.name, 0.0)),
        }
        for component in components
    ]


@REGISTRY.operation(
    "texture.inverse_pole_figure",
    title="Inverse pole figure",
    summary="Which crystal direction lies along a chosen specimen axis.",
    help_text=(
        "The pole figure asks where a crystal plane points in the specimen. The inverse pole "
        "figure asks the opposite and more directly useful question: given a specimen "
        "direction — the rolling direction, the normal direction, the tensile axis — which "
        "crystal direction is aligned with it?\n\n"
        "That is the form in which texture enters a property calculation. A Schmid factor, a "
        "Young's modulus along a loading axis, a Taylor factor: each depends on where the "
        "loading direction sits in the crystal, which is exactly what this figure shows.\n\n"
        "**One point per grain, folded into the fundamental sector.** Symmetry makes many "
        "crystal directions equivalent, so each grain is reduced to the one representative "
        "inside the standard triangle. For a cubic crystal that triangle has [001], [101] and "
        "[111] at its corners, and a point near a corner means that grain has that direction "
        "along the chosen specimen axis.\n\n"
        "**Read it against random.** A random texture fills the triangle evenly by area — which "
        "means more points near the middle than at the corners, because there is more area "
        "there. Clustering at a corner is texture; an even scatter is not."
    ),
    parameters=(
        *_model_parameters(),
        ChoiceParameter(
            name="sample_direction",
            label="Specimen direction",
            help_text="Which specimen axis to look along.",
            options=(
                ("nd", "ND (normal direction)", "The sheet normal: what a rolling plane means."),
                ("rd", "RD (rolling direction)", "The rolling direction."),
                ("td", "TD (transverse direction)", "The transverse direction."),
            ),
            default="nd",
        ),
        _PROJECTION_PARAMETER,
    ),
    returns="One row per grain, with its crystal direction in the fundamental sector.",
    panel="texture",
    citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    tags=("texture", "inverse pole figure", "IPF", "fundamental sector", "specimen direction"),
)
def _inverse_pole_figure(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.texture.models import InversePoleFigure

    orientations, _odf, spec, model_label = _build_texture(request)
    axis_key = str(request["sample_direction"])
    axis_vector = {"rd": (1.0, 0.0, 0.0), "td": (0.0, 1.0, 0.0), "nd": (0.0, 0.0, 1.0)}[axis_key]
    method = str(request["projection"])

    figure = InversePoleFigure.from_orientations(orientations, np.asarray(axis_vector))
    directions = np.asarray(figure.crystal_directions, dtype=float)
    projected = _project(directions, method)
    polar_deg = np.degrees(np.arccos(np.clip(directions[:, 2], -1.0, 1.0)))
    azimuth_deg = np.degrees(np.arctan2(directions[:, 1], directions[:, 0])) % 360.0

    columns: tuple[Column, ...] = (
        Column("grain", "Grain", numeric=True),
        Column(
            "direction",
            "Nearest direction",
            help_text="The low-index crystal direction this grain sits closest to.",
        ),
        Column("x", "x", numeric=True, digits=5),
        Column("y", "y", numeric=True, digits=5),
        Column(
            "polar_deg",
            "From [001]",
            units="°",
            numeric=True,
            digits=3,
            help_text="Angle between the specimen axis and the crystal [001].",
        ),
        Column("azimuth_deg", "Azimuth", units="°", numeric=True, digits=3),
    )
    rows = [
        {
            "grain": index + 1,
            "direction": direction_label(_nearest_low_index(directions[index]), spec=spec),
            "x": float(projected[index, 0]),
            "y": float(projected[index, 1]),
            "polar_deg": float(polar_deg[index]),
            "azimuth_deg": float(azimuth_deg[index]),
        }
        for index in range(directions.shape[0])
    ]

    vertices = np.asarray(figure.project_sector_vertices(method=method), dtype=float)
    vertices = vertices / (math.sqrt(2.0) if method == "equal_area" else 1.0)

    axis_name = _SPECIMEN_AXES[{"rd": 0, "td": 1, "nd": 2}[axis_key]]
    common = _most_common(row["direction"] for row in rows)
    result = AppResult(
        title=f"Inverse pole figure along {axis_name}: {model_label} in {spec.name}",
        summary=(
            f"{len(rows)} grains of the {model_label} texture, each reduced to the crystal "
            "direction "
            f"lying along {axis_name} and folded into the fundamental sector. The most common "
            f"nearest direction is {common[0]}, which {common[1]} of the grains sit closest to. "
            "A random texture scatters evenly by area rather than evenly by eye, so clustering "
            "at a corner is the signature of texture."
        ),
        table=ResultTable(
            columns=columns,
            rows=tuple(rows),
            caption=f"Crystal direction along {axis_name} for every grain.",
        ),
        data={
            "points": rows,
            "sector_vertices": [[float(x), float(y)] for x, y in vertices],
            "projection": method,
            "axis_label": axis_name,
            "model_label": model_label,
            "grain_count": len(rows),
            "columns": [column.to_json() for column in columns],
        },
        inputs=dict(
            _model_inputs(request, spec),
            sample_direction=axis_key,
            projection=method,
        ),
        notes=(
            "One point per grain, not a density: the scatter is the data. A density map would "
            "need a kernel choice, and the point count here is small enough that the scatter is "
            "the more honest picture.",
            "A random texture fills the sector evenly by *area*, so an even-looking scatter is "
            "the absence of texture rather than a uniform distribution of directions.",
        ),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    )
    return result.to_json()


def _nearest_low_index(vector: np.ndarray) -> tuple[int, int, int]:
    """The low-index cubic direction a unit vector sits closest to.

    Deliberately a short list rather than a search: the corners and edges of the
    cubic fundamental sector are what an IPF is read against, and a label of
    ⟨7 3 1⟩ would be precise and useless. A grain in the middle of the triangle
    is genuinely between the named directions, and gets the nearest of them.
    """

    candidates = ((0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 1, 2), (1, 0, 3), (1, 1, 3))
    unit = np.asarray(vector, dtype=float)
    unit = unit / max(float(np.linalg.norm(unit)), 1e-12)
    best = max(
        candidates,
        key=lambda triple: abs(
            float(np.dot(np.asarray(triple, dtype=float) / np.linalg.norm(triple), unit))
        ),
    )
    return best


def _most_common(values: Any) -> tuple[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    if not counts:
        return ("—", 0)
    best = max(counts.items(), key=lambda item: item[1])
    return best


@REGISTRY.operation(
    "texture.odf_sections",
    title="ODF sections",
    summary="The orientation distribution itself, sliced at constant phi-2.",
    help_text=(
        "A pole figure is a projection of the orientation distribution, and projections lose "
        "information — the ghost problem is exactly that loss. The ODF is the distribution "
        "itself, a density over the three Euler angles, and the way it is read is in sections at "
        "constant φ₂.\n\n"
        "**Why φ₂ = 0°, 45° and 65°.** For cubic crystals those three sections carry every "
        "standard rolling and recrystallisation component between them, which is why they are "
        "the sections every fcc texture paper prints. Cube sits at the origin of the φ₂ = 0° "
        "section; brass, copper and S all appear on φ₂ = 45° and 65°.\n\n"
        "**The scale is m.r.d. again**, and means the same thing: how many times more "
        "orientations sit near this point of Euler space than would in a random material. "
        "Unlike a pole figure, an ODF section is not area-preserving in any useful sense — Euler "
        "space is not a metric space for orientations — so read the peaks, not the areas.\n\n"
        "**A caution the sections cannot show.** Euler space distorts badly near Φ = 0, where "
        "φ₁ and φ₂ become degenerate and a single orientation smears along a line. A feature "
        "there is a coordinate artefact as often as it is a texture."
    ),
    parameters=(
        *_model_parameters(),
        NumberParameter(
            name="section_resolution_deg",
            label="Section resolution",
            help_text="Grid spacing in phi-1 and Phi within each section.",
            units="°",
            default=5.0,
            minimum=2.5,
            maximum=15.0,
            advanced=True,
        ),
    ),
    returns="One row per section grid point; the sections themselves under `data.sections`.",
    panel="texture",
    citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_HIRSCH),
    tags=("texture", "ODF", "Euler", "sections", "phi2", "Bunge", "orientation distribution"),
)
def _odf_sections(request: dict[str, Any]) -> dict[str, Any]:
    orientations, odf, spec, model_label = _build_texture(request)
    _spec, phase = phase_from_request(request["phase"])
    resolution = float(request["section_resolution_deg"])
    sections = odf.phi2_sections(
        phi2_deg=(0.0, 45.0, 65.0),
        resolution_deg=resolution,
        normalized=True,
    )

    phi1 = np.asarray(sections.phi1_deg, dtype=float)
    big_phi = np.asarray(sections.big_phi_deg, dtype=float)
    # `normalized=True` normalises the kernel over SO(3), not over the
    # fundamental zone, so a random texture reads |G| rather than 1. Dividing by
    # the reference is what makes this figure's scale the same m.r.d. scale as
    # the pole figures beside it, instead of one 24 times larger wearing the
    # same unit.
    densities = np.asarray(sections.densities, dtype=float) / _odf_random_reference(phase)
    phi2_values = np.asarray(sections.phi2_deg, dtype=float)

    columns: tuple[Column, ...] = (
        Column("phi2_deg", "phi-2", units="°", numeric=True, digits=1),
        Column("phi1_deg", "phi-1", units="°", numeric=True, digits=1),
        Column("big_phi_deg", "Phi", units="°", numeric=True, digits=1),
        Column(
            "mrd",
            "Intensity",
            units="m.r.d.",
            numeric=True,
            digits=4,
            help_text="Multiples of a random distribution at this point of Euler space.",
        ),
    )
    rows: list[dict[str, Any]] = []
    payload: list[dict[str, Any]] = []
    for index, phi2 in enumerate(phi2_values):
        plane = densities[index]
        payload.append(
            {
                "phi2_deg": float(phi2),
                "phi1_deg": [float(value) for value in phi1],
                "big_phi_deg": [float(value) for value in big_phi],
                "densities": [[float(value) for value in row] for row in plane],
                "max_mrd": float(plane.max()),
            }
        )
        for row_index, phi1_value in enumerate(phi1):
            for column_index, big_phi_value in enumerate(big_phi):
                rows.append(
                    {
                        "phi2_deg": float(phi2),
                        "phi1_deg": float(phi1_value),
                        "big_phi_deg": float(big_phi_value),
                        "mrd": float(plane[row_index, column_index]),
                    }
                )

    peak = max(payload, key=lambda entry: float(entry["max_mrd"]))
    fractions = _volume_fractions(orientations, str(request["model"]))
    result = AppResult(
        title=f"ODF sections: {model_label} in {spec.name}",
        summary=(
            f"Three φ₂ sections at 0°, 45° and 65° through the orientation distribution of "
            f"{len(orientations.quaternions)} grains, on a {resolution:g}° grid. The strongest "
            f"feature is {float(peak['max_mrd']):.2f} m.r.d. in the "
            f"φ₂ = {float(peak['phi2_deg']):g}° section. These three sections are the ones every "
            "fcc texture paper prints, because between them they carry every standard rolling "
            "and recrystallisation component."
        ),
        table=ResultTable(
            columns=columns,
            rows=tuple(rows),
            caption="Orientation density on the three sections, in m.r.d.",
        ),
        data={
            "sections": payload,
            "model_label": model_label,
            "grain_count": len(orientations.quaternions),
            "max_mrd": float(max(float(entry["max_mrd"]) for entry in payload)),
            "component_fractions": fractions,
            "resolution_deg": resolution,
            "columns": [column.to_json() for column in columns],
        },
        inputs=dict(_model_inputs(request, spec), section_resolution_deg=resolution),
        notes=(
            "Euler space is not a metric space for orientations: equal volumes of it do not hold "
            "equal ranges of orientation. Read the peaks and where they sit, not the areas.",
            "Near Φ = 0 the coordinates degenerate — φ₁ and φ₂ become the same rotation — so a "
            "single orientation smears along a line there. A feature at the top of a section is "
            "a coordinate artefact as often as it is texture.",
        ),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_HIRSCH),
    )
    return result.to_json()


#: File kinds the measured pole-figure reader accepts.
POLE_FIGURE_SUFFIXES = (".xrdml",)


def _contour_levels(
    text: str,
    *,
    count: int,
    minimum: float,
    maximum: float,
) -> list[float]:
    """The contour levels to draw, from what the user asked for.

    Purpose
    -------
    Contour levels are a reading decision, not a property of the data. The
    conventional set for texture is a chosen sequence — 1, 2, 4, 7, 10 m.r.d. is
    the one most papers use — and picking those is often the difference between
    a figure that shows the texture and one that shows a smooth blob. So an
    explicit list is accepted, and the automatic fallback is only a fallback.

    Parameters
    ----------
    text : str
        Levels as a comma- or space-separated list, in m.r.d. Empty for
        automatic levels.
    count : int
        How many automatic levels to place, when none are given.
    minimum, maximum : float
        The range the automatic levels span.

    Returns
    -------
    list of float
        Ascending, deduplicated, and never empty.

    Raises
    ------
    InvalidInputError
        If the text contains something that is not a number, or every level is
        negative. Silently dropping an unparsable level would draw a figure with
        fewer contours than the user asked for and no indication why.
    """

    cleaned = text.replace(",", " ").split()
    if cleaned:
        levels: list[float] = []
        for token in cleaned:
            try:
                value = float(token)
            except ValueError as error:
                raise InvalidInputError(
                    f"{token!r} is not a number, so it cannot be a contour level.",
                    field="contour_levels",
                    hint="Give levels in m.r.d., separated by commas: 1, 2, 4, 7, 10.",
                ) from error
            if value >= 0.0:
                levels.append(value)
        if not levels:
            raise InvalidInputError(
                "Contour levels must be non-negative intensities in m.r.d.",
                field="contour_levels",
                hint="For example: 1, 2, 4, 7, 10.",
            )
        return sorted(set(levels))

    if not math.isfinite(minimum) or not math.isfinite(maximum) or maximum <= minimum:
        return [max(minimum, 0.0)]
    step = (maximum - minimum) / (count + 1)
    return [minimum + step * (index + 1) for index in range(count)]


def _harmonic_sections(odf: Any, *, resolution_deg: float) -> list[dict[str, Any]]:
    """Slice a harmonic ODF at the three phi-2 sections texture papers print.

    A harmonic ODF has no section method of its own: it is a series, evaluated
    wherever it is asked. It also needs no ``|G|`` correction of the kind
    :func:`_odf_random_reference` applies to a discrete ODF, because its
    densities are already in multiples of a random distribution by
    construction — its mean density over SO(3) is 1.
    """

    from pytex.core.orientation import OrientationSet

    phi1 = np.arange(0.0, 360.0 + 1e-9, resolution_deg)
    big_phi = np.arange(0.0, 90.0 + 1e-9, resolution_deg)
    phi1_mesh, big_phi_mesh = np.meshgrid(phi1, big_phi, indexing="xy")
    payload: list[dict[str, Any]] = []
    for phi2_value in (0.0, 45.0, 65.0):
        angles = np.column_stack(
            [
                phi1_mesh.reshape(-1),
                big_phi_mesh.reshape(-1),
                np.full(phi1_mesh.size, float(phi2_value)),
            ]
        )
        orientations = OrientationSet.from_euler_angles(
            angles,
            crystal_frame=odf.crystal_frame,
            specimen_frame=odf.specimen_frame,
            symmetry=odf.crystal_symmetry,
            phase=odf.phase,
            convention="bunge",
            degrees=True,
        )
        densities = np.asarray(odf.evaluate(orientations), dtype=float).reshape(
            big_phi.size, phi1.size
        )
        payload.append(
            {
                "phi2_deg": float(phi2_value),
                "phi1_deg": [float(value) for value in phi1],
                "big_phi_deg": [float(value) for value in big_phi],
                "densities": [[float(value) for value in row] for row in densities],
                "max_mrd": float(densities.max()),
            }
        )
    return payload


def _measured_harmonic_odf(
    pole_figures: list[Any],
    *,
    degree_bandlimit: int,
    quadrature_step_deg: float,
    halfwidth_deg: float,
    regularization: float,
    ghost_correction: str,
    resolution_deg: float,
) -> dict[str, Any]:
    """Reconstruct an ODF by the series-expansion route, optionally de-ghosted.

    Purpose
    -------
    The classical Bunge method, and the only route in this application on which
    ghost correction is defined. A pole figure obeys Friedel's law, so it
    determines the even-degree half of the expansion and nothing at all of the
    odd half; the even-only solution that results carries false maxima where the
    specimen is empty and depresses the true maxima to pay for them.

    Method and honesty
    ------------------
    Regularized least squares for the symmetry-projected harmonic coefficients,
    then — when asked — an odd part recovered from positivity by
    :func:`pytex.texture.correct_ghosts`. The correction's own cost travels back
    with the sections, because the odd part it supplies is an *inference*, not a
    measurement: no pole-figure experiment can confirm or refute it.

    Returns
    -------
    dict
        The same shape :func:`_measured_odf` returns, plus ``ghost`` describing
        the correction when one was applied.
    """

    from pytex.texture.ghosts import GhostCorrectionSpec
    from pytex.texture.harmonics import HarmonicODF
    from pytex.texture.models import KernelSpec

    spec = (
        None
        if ghost_correction == "none"
        else GhostCorrectionSpec(method=ghost_correction)
    )
    report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=degree_bandlimit,
        regularization=regularization,
        pole_kernel=KernelSpec(halfwidth_deg=halfwidth_deg),
        phi1_step_deg=quadrature_step_deg,
        big_phi_step_deg=quadrature_step_deg,
        phi2_step_deg=quadrature_step_deg,
        ghost_correction=spec,
    )
    odf = report.final_odf
    sections = _harmonic_sections(odf, resolution_deg=resolution_deg)
    ghost: dict[str, Any] | None = None
    if report.ghost_correction is not None:
        correction = report.ghost_correction
        ghost = {
            "method": correction.method,
            "odd_basis_size": correction.odd_basis_size,
            "amplitude_ratio": correction.ghost_amplitude_ratio,
            "negative_before": correction.negative_density_fraction_before,
            "negative_after": correction.negative_density_fraction_after,
            "minimum_before": correction.minimum_density_before,
            "minimum_after": correction.minimum_density_after,
            "maximum_before": correction.maximum_density_before,
            "maximum_after": correction.maximum_density_after,
            "texture_index_before": correction.texture_index_before,
            "texture_index_after": correction.texture_index_after,
            "pole_figure_max_change": correction.pole_figure_max_change,
            "describe": correction.describe(),
        }
    return {
        "sections": sections,
        "max_mrd": float(max(section["max_mrd"] for section in sections)),
        "residual": float(report.relative_residual_norm),
        "dictionary_count": int(report.basis_size),
        "coefficient_count": int(report.basis_size),
        "observation_count": int(report.observation_count),
        "pole_figure_count": len(pole_figures),
        "method": "harmonic",
        "method_label": (
            f"harmonic series to degree {degree_bandlimit}, {report.basis_size} coefficients"
        ),
        "ghost": ghost,
    }


def _measured_odf(
    pole_figures: list[Any],
    *,
    phase: Any,
    dictionary_count: int,
    halfwidth_deg: float,
    resolution_deg: float,
    seed: int = 12345,
) -> dict[str, Any]:
    """Reconstruct an ODF from measured pole figures, and slice it.

    Purpose
    -------
    A pole figure is a projection, and projections lose information: the ghost
    problem is exactly that loss. What physical models need is the orientation
    distribution, and getting it from measurements is the classical inverse
    problem of quantitative texture analysis.

    Method and honesty
    ------------------
    The inversion is :meth:`pytex.texture.ODF.invert_pole_figures` — a
    non-negative regularized least squares over a dictionary of orientations —
    and it is **ill-posed**. The answer depends on the dictionary, the kernel and
    the regularization, and one pole figure cannot constrain it at all. So the
    residual travels back with the sections rather than being discarded, and the
    caller reports it: a reconstruction whose residual is large is a picture of
    the regularization, not of the specimen.

    The dictionary is sampled uniformly on SO(3) — Phi through its sine measure,
    because uniform Euler angles crowd the poles and would weight the
    reconstruction towards them before any data was seen.

    Returns
    -------
    dict
        ``sections`` in the same shape :func:`_odf_sections` produces, plus
        ``residual``, ``dictionary_count`` and ``pole_figure_count``.
    """

    from pytex.core.orientation import OrientationSet
    from pytex.texture.models import ODF, KernelSpec

    generator = np.random.default_rng(seed)
    phi1 = generator.uniform(0.0, 360.0, size=dictionary_count)
    big_phi = np.degrees(np.arccos(generator.uniform(-1.0, 1.0, size=dictionary_count)))
    phi2 = generator.uniform(0.0, 360.0, size=dictionary_count)
    dictionary = OrientationSet.from_euler_angles(
        np.column_stack([phi1, big_phi, phi2]),
        specimen_frame=_specimen_frame(),
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    report = ODF.invert_pole_figures(
        pole_figures,
        orientation_dictionary=dictionary,
        kernel=KernelSpec(halfwidth_deg=halfwidth_deg),
    )
    sections = report.odf.phi2_sections(
        phi2_deg=(0.0, 45.0, 65.0),
        resolution_deg=resolution_deg,
        normalized=True,
    )
    reference = _odf_random_reference(phase)
    densities = np.asarray(sections.densities, dtype=float) / reference
    phi1_axis = [float(value) for value in np.asarray(sections.phi1_deg, dtype=float)]
    big_phi_axis = [float(value) for value in np.asarray(sections.big_phi_deg, dtype=float)]
    payload = [
        {
            "phi2_deg": float(phi2_value),
            "phi1_deg": phi1_axis,
            "big_phi_deg": big_phi_axis,
            "densities": [[float(value) for value in row] for row in densities[index]],
            "max_mrd": float(densities[index].max()),
        }
        for index, phi2_value in enumerate(np.asarray(sections.phi2_deg, dtype=float))
    ]
    return {
        "sections": payload,
        "max_mrd": float(densities.max()),
        "residual": float(report.residual_norm),
        "dictionary_count": int(dictionary_count),
        "pole_figure_count": len(pole_figures),
        "method": "dictionary",
        "method_label": f"non-negative dictionary of {int(dictionary_count)} orientations",
        # Ghost correction is defined on the harmonic expansion, where the odd
        # part is an explicit set of coefficients. The dictionary route has no
        # such split: its non-negativity constraint acts on the weights.
        "ghost": None,
    }


@REGISTRY.operation(
    "texture.measured_pole_figures",
    title="Measured pole figures",
    summary="Open XRDML pole-figure files and draw them on one shared intensity scale.",
    help_text=(
        "Reads Panalytical **XRDML** pole-figure files — one file per measured reflection — and "
        "draws every one of them, in tabs, so a set of {111}, {200} and {220} figures is one "
        "result rather than three.\n\n"
        "**Say which reflection each file is.** The file records the diffraction angle, not the "
        "plane, so the poles are given here in the order the files were opened. Getting the "
        "order wrong is the mistake this makes easiest, and it is visible: the figures will not "
        "be consistent with any single texture.\n\n"
        "**One scale across all of them, by default.** Two pole figures of the same specimen "
        "drawn on separate scales cannot be compared, and comparing them is the entire reason "
        "for measuring more than one. Turn the shared scale off only to look at a weak figure "
        "in isolation.\n\n"
        "**Contour levels are a reading decision.** Give them explicitly — `1, 2, 4, 7, 10` is "
        "the sequence most of the texture literature uses — or leave the field empty for evenly "
        "spaced ones. The levels apply to every figure, which is what makes the set "
        "comparable.\n\n"
        "**Normalisation decides what the numbers mean.** A measured figure arrives in detector "
        "counts, which are a property of the instrument and the counting time. *m.r.d.* rescales "
        "it so that a texture-free specimen would read 1 everywhere, which is the only form in "
        "which two instruments' figures mean the same thing."
    ),
    parameters=(
        ObjectParameter(
            name="files",
            label="XRDML files",
            help_text=(
                'The opened files, as `{"items": [{"name": ..., "text": ...}]}`. Supplied '
                "by the **Open pole figures** control rather than typed."
            ),
            required=True,
        ),
        phase_parameter(help_text="The phase the measured reflections belong to."),
        IndicesListParameter(
            name="poles",
            label="Plane of each file {hkl}",
            help_text=(
                "One plane per file, in the order the files were opened — 1 1 1 on one line, "
                "2 0 0 on the next. With fewer lines than files, the last plane is reused."
            ),
            default=((1, 1, 1),),
        ),
        ChoiceParameter(
            name="intensity_normalization",
            label="Normalisation",
            help_text=(
                "What the intensities are rescaled to. m.r.d. is the only one that makes two "
                "instruments comparable; the others are for inspecting a file as recorded."
            ),
            options=(
                (
                    "mrd",
                    "m.r.d.",
                    "Multiples of a random distribution: 1 everywhere for a texture-free specimen.",
                ),
                ("max", "Peak = 1", "Scaled so the strongest point is 1."),
                ("none", "As recorded", "Raw detector counts, exactly as the file holds them."),
            ),
            default="mrd",
        ),
        _PROJECTION_PARAMETER,
        TextParameter(
            name="contour_levels",
            label="Contour levels",
            help_text=(
                "Levels to draw, in the normalised intensity unit, separated by commas — "
                "`1, 2, 4, 7, 10` is the conventional texture sequence. Leave empty for evenly "
                "spaced levels."
            ),
            required=False,
            placeholder="1, 2, 4, 7, 10",
        ),
        IntegerParameter(
            name="contour_count",
            label="Automatic levels",
            help_text="How many evenly spaced levels to place when none are given.",
            default=6,
            minimum=1,
            maximum=20,
        ),
        BooleanParameter(
            name="reconstruct_odf",
            label="Reconstruct the ODF",
            help_text=(
                "Invert the opened figures into an orientation distribution and add it as a "
                "further tab, sliced at the three phi-2 sections texture papers print.\n\n"
                "**This is an ill-posed inversion.** Pole figures are projections and lose the "
                "odd-order information, so the answer depends on the dictionary, the kernel and "
                "the regularization. One pole figure cannot constrain it at all; three from "
                "different planes is the usual minimum. The residual is reported beside the "
                "sections, and a large one means you are looking at the regularization rather "
                "than at the specimen."
            ),
            default=False,
        ),
        IntegerParameter(
            name="dictionary_count",
            label="Dictionary orientations",
            help_text=(
                "How many orientations the inversion solves over. More resolves a sharper "
                "texture and costs time roughly linearly."
            ),
            default=800,
            minimum=100,
            maximum=5000,
            advanced=True,
        ),
        NumberParameter(
            name="odf_halfwidth_deg",
            label="ODF kernel halfwidth",
            help_text=(
                "Width of the bell on each dictionary orientation. As for a model texture this "
                "is a smoothing choice, and here it is also the regularization: too small and "
                "the inversion fits noise."
            ),
            units="°",
            default=10.0,
            minimum=2.0,
            maximum=30.0,
            advanced=True,
        ),
        ChoiceParameter(
            name="odf_method",
            label="Inversion route",
            help_text=(
                "How the orientation distribution is solved for.\n\n"
                "**Dictionary** fits non-negative weights on a cloud of orientations. It cannot "
                "produce a negative density because it is not allowed to, and it has no explicit "
                "odd part to correct.\n\n"
                "**Harmonic series** is the classical Bunge expansion. It is the only route on "
                "which ghost correction is defined, because it is the only one that separates "
                "the even part a pole figure determines from the odd part it cannot see."
            ),
            options=(
                (
                    "dictionary",
                    "Non-negative dictionary",
                    "Weights on an orientation cloud; never negative, no explicit odd part.",
                ),
                (
                    "harmonic",
                    "Harmonic series (Bunge)",
                    "Symmetry-projected coefficients; the route ghost correction applies to.",
                ),
            ),
            default="dictionary",
            advanced=True,
        ),
        ChoiceParameter(
            name="ghost_correction",
            label="Ghost correction",
            help_text=(
                "Recover the odd part of the distribution, which a pole figure cannot measure.\n\n"
                "A diffraction pole figure obeys Friedel's law and cannot tell a plane normal "
                "from its opposite, so it fixes only the even-degree half of the expansion. "
                "Setting the odd half to zero is not neutral: it puts false maxima where the "
                "specimen is empty and depresses the true maxima to pay for them. Correction "
                "adds the smallest odd part that makes the density non-negative.\n\n"
                "**Zero range** asks for more — that the density stay at zero wherever the "
                "measurement says the specimen has no such orientations at all.\n\n"
                "**What it costs.** The odd part is an inference from positivity, not a "
                "measurement, and no pole-figure experiment can confirm or refute it. The size "
                "of the inference is reported beside the sections. Applies to the harmonic "
                "route only."
            ),
            options=(
                (
                    "none",
                    "None (even part only)",
                    "Report what the data determine, and nothing else.",
                ),
                (
                    "positivity",
                    "Positivity",
                    "The smallest odd part that makes the density non-negative.",
                ),
                (
                    "zero_range",
                    "Zero range",
                    "Positivity, plus holding the measured empty range empty.",
                ),
            ),
            default="none",
            advanced=True,
        ),
        IntegerParameter(
            name="odf_bandlimit",
            label="Harmonic bandlimit",
            help_text=(
                "Highest harmonic degree retained by the series route. Higher resolves a sharper "
                "texture, costs steeply more, and is more sensitive to noise. Ignored by the "
                "dictionary route.\n\n"
                "**It also decides whether ghost correction can do anything.** The odd part is "
                "expanded to the same degree, and a symmetry admits odd terms only where it has "
                "an odd-degree invariant: for a cubic material the first one is at **degree 9**, "
                "so a cubic ODF expanded to degree 6 or 8 has no ghost part to correct and the "
                "correction will say so. Lower symmetries admit odd terms much earlier."
            ),
            default=6,
            minimum=2,
            maximum=16,
            advanced=True,
        ),
        NumberParameter(
            name="odf_regularization",
            label="Harmonic regularization",
            help_text=(
                "Tikhonov weight on the harmonic coefficients. This is the whole defence against "
                "an under-determined fit: with fewer measured points than coefficients, the "
                "unregularized solution is a picture of the null space rather than of the "
                "specimen. Larger is smoother and more stable."
            ),
            default=0.01,
            minimum=1e-8,
            maximum=1.0,
            advanced=True,
        ),
        BooleanParameter(
            name="shared_scale",
            label="One scale for every figure",
            help_text=(
                "Put every figure on the same intensity range, so they can be compared. Off "
                "gives each its own range, which makes a weak figure legible and makes it "
                "impossible to read against the others."
            ),
            default=True,
        ),
    ),
    returns=(
        "One row per measured point of every figure; the figures, the shared scale and the "
        "contour levels under `data`."
    ),
    panel="texture",
    citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    tags=("texture", "pole figure", "XRDML", "measurement", "import", "contour", "m.r.d."),
)
def _measured_pole_figures(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.adapters.xrdml import read_xrdml_pole_figure
    from pytex.app.uploads import uploaded_file

    spec, phase = phase_from_request(request["phase"])
    payload = request["files"]
    items = payload.get("items") if isinstance(payload, dict) else None
    if not isinstance(items, list) or not items:
        raise InvalidInputError(
            "No pole-figure file has been opened.",
            field="files",
            hint="Open one or more .xrdml files with the control above.",
        )
    poles = list(request["poles"]) or [(1, 1, 1)]
    normalization = str(request["intensity_normalization"])
    method = str(request["projection"])
    specimen = _specimen_frame()

    figures: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    reconstructed: list[Any] = []
    for index, item in enumerate(items):
        indices = tuple(int(value) for value in poles[min(index, len(poles) - 1)])
        with uploaded_file(item, field="files", suffixes=POLE_FIGURE_SUFFIXES) as (path, name):
            try:
                measurement = read_xrdml_pole_figure(path)
            # Broad on purpose: the XML reader raises ParseError, ValueError,
            # KeyError and TypeError depending on which part of the document
            # is wrong, and every one of them means the same thing to a user.
            except Exception as error:
                raise InvalidInputError(
                    f"{name} could not be read as an XRDML pole figure: {error}",
                    field="files",
                    hint=(
                        "The file must be a pole-figure measurement rather than a line scan: "
                        "it needs both a Phi and a Psi (or Chi) axis."
                    ),
                ) from error
        pole_figure = measurement.to_pole_figure(
            _crystal_plane(phase, indices),
            specimen_frame=specimen,
            intensity_normalization=normalization,
        )
        directions = np.asarray(pole_figure.sample_directions, dtype=float)
        # A pole figure is antipodal, so a point measured below the equator is
        # the same pole as its reflection above it. Folding here rather than
        # discarding keeps every measured count in the figure.
        directions = np.where(directions[:, 2:3] < 0.0, -directions, directions)
        values = np.asarray(pole_figure.intensities, dtype=float).reshape(-1)
        projected = _project(directions, method)
        # A measured pole figure collects the whole symmetry family: the
        # specimen is rotated through every orientation that puts any member of
        # {hkl} into the diffraction condition. So it is written {hkl}, not
        # (hkl) — the notation registry treats that distinction as meaning, not
        # as style.
        label = family_label(indices, spec=spec, family="plane")
        points = [
            {
                "mrd": float(values[point]),
                "x": float(projected[point, 0]),
                "y": float(projected[point, 1]),
                "polar_deg": float(
                    math.degrees(math.acos(min(1.0, max(-1.0, float(directions[point, 2])))))
                ),
                "azimuth_deg": float(
                    math.degrees(math.atan2(directions[point, 1], directions[point, 0])) % 360.0
                ),
            }
            for point in range(directions.shape[0])
        ]
        rows.extend({**point, "figure": label, "file": name} for point in points)
        reconstructed.append(pole_figure)
        # The identifier drawn on the figure itself. A plate of six discs is
        # unreadable without one, and the file name is the only identifier a
        # measurement always has; the sample name in the file is better when it
        # is there, which for a lab instrument is most of the time.
        sample_label = (measurement.sample_name or "").strip() or name.rsplit(".", 1)[0]
        figures.append(
            {
                "file": name,
                "sample_label": sample_label,
                "label": label,
                "indices": list(indices),
                "points": points,
                "minimum": float(values.min()),
                "maximum": float(values.max()),
                "mean": float(values.mean()),
                "count": int(values.size),
                "two_theta_deg": (
                    float(np.mean(measurement.two_theta_deg))
                    if measurement.two_theta_deg is not None
                    else None
                ),
                "sample_name": measurement.sample_name,
            }
        )

    # A label that is the same on every panel identifies nothing. The sample
    # name in the file is the better identifier when it distinguishes the
    # figures, and the file name is the only one that always does.
    if len({figure["sample_label"] for figure in figures}) < len(figures):
        for figure in figures:
            figure["sample_label"] = str(figure["file"]).rsplit(".", 1)[0]

    overall_min = min(figure["minimum"] for figure in figures)
    overall_max = max(figure["maximum"] for figure in figures)
    shared = bool(request["shared_scale"])
    levels = _contour_levels(
        str(request["contour_levels"] or ""),
        count=int(request["contour_count"]),
        minimum=overall_min,
        maximum=overall_max,
    )
    for figure in figures:
        figure["scale"] = (
            {"minimum": overall_min, "maximum": overall_max}
            if shared
            else {"minimum": figure["minimum"], "maximum": figure["maximum"]}
        )

    odf: dict[str, Any] | None = None
    if bool(request["reconstruct_odf"]):
        APP_LOG.info(
            f"Inverting {len(reconstructed)} pole figure(s) over "
            f"{int(request['dictionary_count'])} orientations.",
            source="texture.measured_pole_figures",
        )
        if str(request["odf_method"]) == "harmonic":
            odf = _measured_harmonic_odf(
                reconstructed,
                degree_bandlimit=int(request["odf_bandlimit"]),
                # Coarser than the section grid on purpose: the quadrature is
                # three-dimensional, so its cost is the cube of the step, and it
                # only has to resolve the response kernel.
                quadrature_step_deg=20.0,
                halfwidth_deg=float(request["odf_halfwidth_deg"]),
                regularization=float(request["odf_regularization"]),
                ghost_correction=str(request["ghost_correction"]),
                resolution_deg=5.0,
            )
        else:
            odf = _measured_odf(
                reconstructed,
                phase=phase,
                dictionary_count=int(request["dictionary_count"]),
                halfwidth_deg=float(request["odf_halfwidth_deg"]),
                resolution_deg=5.0,
            )
        APP_LOG.notice(
            f"ODF reconstructed: peak {odf['max_mrd']:.2f} m.r.d., residual {odf['residual']:.4g}.",
            source="texture.measured_pole_figures",
            detail={"residual": odf["residual"], "peak_mrd": odf["max_mrd"]},
        )

    unit = {"mrd": "m.r.d.", "max": "peak = 1", "none": "counts"}[normalization]
    names = ", ".join(figure["label"] for figure in figures)
    notes = [
        "These are measurements, so no answer is known in advance. A figure whose intensities "
        "are not normalised to m.r.d. means nothing outside its own file: detector counts "
        "depend on the counting time and the instrument.",
        "Defocusing and absorption corrections are not applied here. A reflection measured to "
        "high tilt falls off for instrumental reasons as well as textural ones, and the "
        "outer rim of an uncorrected figure is the part to distrust.",
        "The plane assigned to each file comes from the order the files were opened, not from "
        "the file: XRDML records the diffraction angle rather than the reflection.",
        "ODF reconstruction from pole figures is ill-posed: projections lose the odd-order "
        "information, so the result depends on the dictionary, the kernel and the "
        "regularization. Three figures from different planes is the usual minimum, and the "
        "reported residual is what says whether the estimate is worth anything.",
    ]
    if odf is not None and odf.get("coefficient_count", 0) > int(odf.get("observation_count", 0)):
        # Stated rather than blocked: an under-determined fit is legitimate when
        # the regularization is doing the work knowingly, and misleading when it
        # is not. The reader is the one who can tell the difference.
        notes.append(
            f"The harmonic fit solved for {odf['coefficient_count']} coefficients from "
            f"{odf['observation_count']} measured intensities. There are fewer data than "
            "unknowns, so the regularization, not the specimen, is deciding the part of the "
            "answer the data leave free: lower the bandlimit, measure more directions, or read "
            "the result as a smoothed lower bound on the texture."
        )
    if odf is not None and odf.get("ghost") is not None:
        ghost = odf["ghost"]
        notes.append(
            "Ghost correction was requested but did nothing: this symmetry admits no odd-degree "
            "harmonic term at or below the chosen bandlimit, so the even-only solution already "
            "spans every function the symmetry allows. For a cubic material the first odd "
            "invariant is at degree 9."
            if ghost["odd_basis_size"] == 0
            else (
                "The odd part the ghost correction supplied is an inference from positivity, not "
                "a measurement: no pole-figure experiment can confirm or refute it. It changed "
                "the pole densities the ODF predicts at the measured directions by at most "
                f"{ghost['pole_figure_max_change']:.3g} m.r.d., which is the check that it did "
                "not buy positivity with data agreement."
            )
        )
    result = AppResult(
        title=f"Measured pole figures of {spec.name}: {names}",
        summary=(
            f"{len(figures)} pole figure(s) read from XRDML, {sum(f['count'] for f in figures)} "
            f"measured points in total, in {unit}. "
            + (
                f"All of them are drawn on one scale, {overall_min:.3g} to {overall_max:.3g}, "
                "so they can be compared."
                if shared
                else "Each is drawn on its own scale, so they cannot be compared with each other."
            )
            + f" Contours at {', '.join(f'{level:g}' for level in levels)}."
            + (
                ""
                if odf is None
                else (
                    f" The ODF reconstructed from them peaks at {odf['max_mrd']:.2f} m.r.d. with "
                    f"a residual of {odf['residual']:.4g}, over a "
                    f"{odf.get('method_label', 'dictionary')}."
                )
            )
            + (
                ""
                if odf is None or odf.get("ghost") is None
                else (
                    " Ghost correction did nothing, and could not: the crystal and specimen "
                    f"symmetries admit no odd-degree harmonic term at degree "
                    f"{int(request['odf_bandlimit'])} or below, so the even-only solution "
                    "already spans every function the symmetry allows. For a cubic material the "
                    "first odd invariant is at degree 9."
                    if odf["ghost"]["odd_basis_size"] == 0
                    else f" Ghost correction by {odf['ghost']['method']} added an odd part "
                    f"{odf['ghost']['amplitude_ratio']:.3f} times the even one, taking the "
                    f"negative density from {odf['ghost']['negative_before']:.1%} of orientation "
                    f"space to {odf['ghost']['negative_after']:.1%} and the peak from "
                    f"{odf['ghost']['maximum_before']:.2f} to "
                    f"{odf['ghost']['maximum_after']:.2f} m.r.d. That odd part is an inference "
                    "from positivity, not a measurement."
                )
            )
        ),
        table=ResultTable(
            columns=(
                Column("figure", "Figure"),
                *_DENSITY_COLUMNS,
                Column("file", "Source file"),
            ),
            rows=tuple(rows),
            caption=f"Every measured point of every opened figure, in {unit}.",
        ),
        data={
            "figures": figures,
            "levels": levels,
            "shared_scale": shared,
            "scale": {"minimum": overall_min, "maximum": overall_max},
            "odf": odf,
            "projection": method,
            "unit": unit,
            "specimen_axes": list(_SPECIMEN_AXES),
            "columns": [column.to_json() for column in _DENSITY_COLUMNS],
        },
        inputs={
            "phase": spec.to_json(),
            "poles": [list(pole) for pole in poles],
            "intensity_normalization": normalization,
            "projection": method,
            "contour_levels": str(request["contour_levels"] or ""),
            "contour_count": int(request["contour_count"]),
            "shared_scale": shared,
            "reconstruct_odf": bool(request["reconstruct_odf"]),
            "dictionary_count": int(request["dictionary_count"]),
            "odf_halfwidth_deg": float(request["odf_halfwidth_deg"]),
            "odf_method": str(request["odf_method"]),
            "odf_bandlimit": int(request["odf_bandlimit"]),
            "odf_regularization": float(request["odf_regularization"]),
            "ghost_correction": str(request["ghost_correction"]),
            "files": [str(item.get("name", "")) for item in items if isinstance(item, dict)],
        },
        notes=tuple(notes),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="texture.example.random_baseline",
            title="What no texture looks like",
            panel="texture",
            summary="A random polycrystal: 1 m.r.d. everywhere, by construction.",
            teaches=(
                "Run this first and read the number under the figure: the mean is 1.000 m.r.d. "
                "and so is the maximum, to within the sampling noise of a few thousand grains. "
                "That flatness is what every other figure in this panel is measured against, and "
                "checking it is how you know a pole figure has been normalised at all."
            ),
            operation="texture.pole_figure",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "random",
                "pole": [1, 1, 1],
                "spread_deg": 10.0,
            },
        ),
        ExampleScenario(
            id="texture.example.goss",
            title="Goss: the poles land where the label says",
            panel="texture",
            summary="{011}<100> — so the (011) poles sit at the centre of the figure.",
            teaches=(
                "The Goss component is written {011}<100>, which means the {011} plane lies in "
                "the sheet plane. Plot (011) and the strongest pole is at the centre of the "
                "figure, which is ND — exactly what the notation asserts. This is the check that "
                "the whole pipeline puts poles where crystallography says, and it needs no "
                "reference figure to verify."
            ),
            operation="texture.pole_figure",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "goss",
                "pole": [0, 1, 1],
                "spread_deg": 8.0,
            },
        ),
        ExampleScenario(
            id="texture.example.fcc_rolling",
            title="The fcc rolling texture, on (111)",
            panel="texture",
            summary="Cube, Goss, brass, copper and S together, as a rolled sheet shows them.",
            teaches=(
                "This is the figure in every rolling-texture paper. Five components put five "
                "sets of poles on one figure, and they overlap into the arcs that make a rolled "
                "fcc texture recognisable at a glance. Compare it against the single-component "
                "figures to see which arc belongs to which component."
            ),
            operation="texture.pole_figure",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "fcc_rolling",
                "pole": [1, 1, 1],
                "spread_deg": 10.0,
            },
        ),
        ExampleScenario(
            id="texture.example.ipf_nd",
            title="Which direction lies along the sheet normal",
            panel="texture",
            summary="The inverse pole figure of the fcc rolling texture along ND.",
            teaches=(
                "The same texture asked the other way round: not where a plane points, but which "
                "crystal direction is along ND. This is the form texture enters a property "
                "calculation in — a modulus or a Schmid factor along a loading axis depends on "
                "exactly this. Note the clustering: an untextured material would scatter evenly "
                "over the triangle by area."
            ),
            operation="texture.inverse_pole_figure",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "fcc_rolling",
                "sample_direction": "nd",
                "spread_deg": 10.0,
            },
        ),
        ExampleScenario(
            id="texture.example.odf_cube",
            title="Cube sits at the origin of Euler space",
            panel="texture",
            summary="ODF sections of the cube component, whose ideal orientation is (0, 0, 0).",
            teaches=(
                "The cube component is Bunge (0°, 0°, 0°), so its density peaks at the corner of "
                "the φ₂ = 0° section — and, because of cubic symmetry, at the equivalent corners "
                "of the others. It is the clearest demonstration that an ODF section is a slice "
                "through a density over Euler angles, not a projection of anything."
            ),
            operation="texture.odf_sections",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "cube",
                "spread_deg": 8.0,
            },
        ),
        ExampleScenario(
            id="texture.example.odf_rolling",
            title="The rolling texture, as an ODF",
            panel="texture",
            summary="The three canonical sections of the five-component fcc rolling texture.",
            teaches=(
                "The φ₂ = 45° and 65° sections are where brass, copper and S live, and this is "
                "the pair of sections a rolling-texture paper prints to argue about which "
                "component dominates. Compare the peaks against the volume fractions listed "
                "beside the figure — the ODF shows where, the fractions say how much."
            ),
            operation="texture.odf_sections",
            request={
                "phase": {"builtin": "ni_fcc"},
                "model": "fcc_rolling",
                "spread_deg": 10.0,
            },
        ),
    )
)
