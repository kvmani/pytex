"""Measured texture analysis: pole figures loaded once, read every way.

What it does
    One request takes a set of measured pole figures, the crystal (from the phase)
    and an assumed sample symmetry, and returns everything quantitative texture
    analysis reads from them:

    * the **measured** figures, as recorded and normalised to m.r.d.;
    * the figures with the **sample symmetry imposed**, which is what is inverted;
    * the **ODF** reconstructed from them, sliced at constant phi2 (the default),
      phi1 or sigma, including a LaboTex-style plate of every section;
    * the figures **recalculated** from that ODF and their **difference** from the
      measurement, with the residual numbers that say whether the ODF is usable;
    * the **volume fractions** of the ideal orientations of the crystal system.

Why one operation rather than one per view
    The inputs are the same for every view, and asking a user to restate the
    file, the crystal and the sample symmetry for each is how two views of "the
    same" texture end up disagreeing. The panel therefore sends one request and
    shows its parts in tabs. The inversion is cached on its inputs, so choosing a
    different section or tolerance re-reads the ODF instead of re-solving it.

Honesty about the inversion
    Pole-figure inversion is ill-posed. The recalculated and difference figures
    are here precisely so the reader can see how well the ODF reproduces the data
    it came from - an ODF that cannot is a picture of the regularisation.

See ``docs/site/workflows/texture_analysis_workbench.md``.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import OrderedDict
from dataclasses import replace
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.logbook import APP_LOG
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    ExampleScenario,
    IndicesListParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultMetric, ResultStage, ResultTable
from pytex.app.services.calculator import family_label, phase_parameter
from pytex.app.services.texture import (
    _CITATION_BUNGE,
    _CITATION_RANDLE,
    _PROJECTION_PARAMETER,
    _SPECIMEN_AXES,
    POLE_FIGURE_SUFFIXES,
    _contour_levels,
    _crystal_plane,
    _project,
    _specimen_frame,
)
from pytex.app.services.texture_figures import (
    fraction_figure,
    parity_figure,
    tilt_misfit_figure,
)

__all__: tuple[str, ...] = ()

_CITATION_MATTHIES = (
    "Matthies, Wenk & Vinel, J. Appl. Cryst. 21 (1988) 285 (the RP factor of recalculated "
    "pole figures)."
)
_CITATION_LABOTEX = (
    "Pawlik & Ozga, LaboTex: The Texture Analysis Software, Goettinger Arbeiten zur Geologie und "
    "Palaeontologie SB4 (1999) (section plates and recalculated-figure diagnostics)."
)

_PANEL = "texture_analysis"

_SAMPLE_SYMMETRY_OPTIONS = (
    ("triclinic", "None (triclinic)", "No symmetry assumed: every direction is independent."),
    (
        "monoclinic",
        "Monoclinic",
        "One two-fold axis along ND, as in asymmetric or cross rolling.",
    ),
    (
        "orthorhombic",
        "Orthorhombic (rolling)",
        "Two-fold axes along RD, TD and ND: the conventional assumption for rolled sheet.",
    ),
    (
        "axial",
        "Axial (fibre)",
        "Every rotation about ND: drawn wire, extruded rod, a tube read along its axis.",
    ),
)

_SECTION_KIND_OPTIONS = (
    ("phi2", "Constant phi-2", "phi-1 across, Phi down: the convention of most texture papers."),
    ("phi1", "Constant phi-1", "phi-2 across, Phi down."),
    ("sigma", "Constant sigma", "sigma = phi-1 + phi-2: the view of the bcc gamma fibre."),
)

_SECTION_PRESET_OPTIONS = (
    (
        "standard",
        "Standard sections",
        "The sections the literature prints: phi-2 = 0, 45 and 65 degrees for cubic, 0 and 30 "
        "for hexagonal.",
    ),
    (
        "labotex",
        "Every section (LaboTex plate)",
        "Every section at the chosen step over the whole symmetry range, on one scale.",
    ),
    ("custom", "Chosen values", "The values typed in Section values."),
)

_INVERSION_GROUP = "ODF inversion"

#: Display names of the catalogue components, so a table never shows a key.
_COMPONENT_TITLES = {
    "cube": "Cube",
    "goss": "Goss",
    "brass": "Brass",
    "copper": "Copper",
    "s": "S",
    "rotated_cube": "Rotated cube",
    "rotated_goss": "Rotated Goss",
    "basal": "Basal",
    "basal_30_td": "Basal, 30° to TD",
    "basal_30_rd": "Basal, 30° to RD",
    "c_along_td": "c along TD",
    "c_along_rd": "c along RD",
}

#: Cached inversions, keyed on every input that changes the ODF.
_INVERSION_CACHE: OrderedDict[str, dict[str, Any]] = OrderedDict()
_INVERSION_CACHE_SIZE = 4

#: The demonstration raster: what a laboratory goniometer measures.
_DEMO_PSI_DEG = np.arange(0.0, 75.0 + 1e-9, 5.0)
_DEMO_PHI_DEG = np.arange(0.0, 360.0 - 1e-9, 5.0)

#: Measured intensities below this many m.r.d. are left out of the RP factor,
#: whose relative deviation is meaningless where the measurement is near zero.
_RP_THRESHOLD_MRD = 0.5


def _file_items(payload: Any) -> list[dict[str, Any]]:
    if payload is None:
        return []
    items = payload.get("items") if isinstance(payload, dict) else None
    if items is None:
        return []
    if not isinstance(items, list):
        raise InvalidInputError(
            "The opened files arrived in an unexpected shape.",
            field="files",
            hint="Open the files again with the control above.",
        )
    return [item for item in items if isinstance(item, dict)]


def _item_digest(item: dict[str, Any]) -> str:
    body = str(item.get("text") or item.get("data_base64") or "")
    return hashlib.sha256(f"{item.get('name', '')}\0{body}".encode()).hexdigest()


def _is_cubic(phase: Any) -> bool:
    symmetry = getattr(phase, "symmetry", None)
    return symmetry is not None and int(np.asarray(symmetry.operators).shape[0]) == 24


def _components_for(phase: Any) -> tuple[Any, ...]:
    from pytex.core.hexagonal import is_hexagonal_phase
    from pytex.texture import components as catalogue

    if _is_cubic(phase):
        seen: dict[str, Any] = {}
        for component in (
            *catalogue.STANDARD_FCC_ROLLING_COMPONENTS,
            *catalogue.STANDARD_BCC_ROLLING_COMPONENTS,
        ):
            seen.setdefault(component.name, component)
        return tuple(seen.values())
    if is_hexagonal_phase(phase):
        return tuple(catalogue.STANDARD_HCP_COMPONENTS)
    return ()


def _demonstration_figures(
    phase: Any, spec: PhaseSpec, poles: list[tuple[int, ...]], seed: int
) -> list[dict[str, Any]]:
    """Pole figures a goniometer would record from a known model texture.

    The truth is known, so the demonstration can be read as a test: the ODF
    recalculates the figures to within the 3 percent counting noise added here,
    and the volume fractions name the components the model was built from.
    """

    from pytex.core.hexagonal import is_hexagonal_phase
    from pytex.core.orientation import OrientationSet
    from pytex.core.symmetry import SymmetrySpec
    from pytex.texture import components as catalogue
    from pytex.texture.models import ODF, KernelSpec, PoleFigure, random_pole_density

    if _is_cubic(phase):
        centres = [c.bunge_euler_deg for c in catalogue.STANDARD_FCC_ROLLING_COMPONENTS]
        texture_name = "an fcc rolling texture (cube, Goss, brass, copper and S)"
    elif is_hexagonal_phase(phase):
        centres = [(180.0, 30.0, 0.0), (0.0, 30.0, 0.0)]
        texture_name = "a split-basal texture (basal poles 30° either side of ND towards TD)"
    else:
        raise InvalidInputError(
            f"The demonstration figures are defined for cubic and hexagonal phases; "
            f"{spec.name} is neither.",
            field="files",
            hint="Open your own pole-figure files, or choose a cubic or hexagonal phase.",
        )
    generator = np.random.default_rng(seed)
    angles = np.concatenate(
        [np.asarray(c, dtype=float) + generator.normal(0.0, 8.0, size=(150, 3)) for c in centres]
    )
    scattered = OrientationSet.from_euler_angles(
        angles,
        specimen_frame=_specimen_frame(),
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    # The model is given the sample symmetry of the process it imitates - rolling
    # is orthorhombic - so an analysis assuming that symmetry is tested against a
    # specimen that actually has it, and one assuming axial symmetry is not.
    # Specimen operations act on the left of a crystal-to-specimen matrix.
    operators = np.asarray(SymmetrySpec.specimen("orthorhombic").operators, dtype=float)
    matrices = np.einsum("oij,njk->onik", operators, scattered.as_matrices()).reshape(-1, 3, 3)
    orientations = OrientationSet.from_matrices(
        matrices,
        specimen_frame=_specimen_frame(),
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    kernel = KernelSpec(halfwidth_deg=8.0)
    model = ODF.from_orientations(orientations, kernel=kernel)
    psi, phi = np.meshgrid(_DEMO_PSI_DEG, _DEMO_PHI_DEG, indexing="ij")
    theta = np.radians(psi.ravel())
    azimuth = np.radians(phi.ravel())
    directions = np.column_stack(
        [np.sin(theta) * np.cos(azimuth), np.sin(theta) * np.sin(azimuth), np.cos(theta)]
    )
    reference = float(random_pole_density(kernel))
    figures = []
    for indices in poles:
        plane = _crystal_plane(phase, indices)
        mrd = np.asarray(model.evaluate_pole_density(plane, directions), dtype=float) / reference
        noisy = np.clip(mrd * generator.normal(1.0, 0.03, size=mrd.shape), 0.0, None)
        label = family_label(indices, spec=spec, family="plane")
        figures.append(
            {
                "file": f"demonstration {label}",
                "sample_label": "demonstration",
                "figure": PoleFigure(
                    pole=plane,
                    sample_directions=directions,
                    intensities=noisy,
                    specimen_frame=_specimen_frame(),
                    sampling="sampled_density",
                ),
                "texture_name": texture_name,
            }
        )
    return figures


def _read_figures(
    items: list[dict[str, Any]], phase: Any, poles: list[tuple[int, ...]]
) -> list[dict[str, Any]]:
    from pytex.adapters.xrdml import read_xrdml_pole_figure
    from pytex.app.uploads import uploaded_file

    figures = []
    for index, item in enumerate(items):
        indices = poles[min(index, len(poles) - 1)]
        with uploaded_file(item, field="files", suffixes=POLE_FIGURE_SUFFIXES) as (path, name):
            try:
                measurement = read_xrdml_pole_figure(path)
            # Broad on purpose, as in the measured pole-figure view: every parser
            # failure means the same thing to the person who opened the file.
            except Exception as error:
                raise InvalidInputError(
                    f"{name} could not be read as an XRDML pole figure: {error}",
                    field="files",
                    hint="The file must be a pole-figure measurement with Phi and Psi axes.",
                ) from error
        figure = measurement.to_pole_figure(
            _crystal_plane(phase, indices),
            specimen_frame=_specimen_frame(),
            intensity_normalization="mrd",
        )
        sample_label = (measurement.sample_name or "").strip() or name.rsplit(".", 1)[0]
        figures.append(
            {"file": name, "sample_label": sample_label, "figure": figure, "texture_name": None}
        )
    return figures


def _invert(
    figures: list[Any], phase: Any, request: dict[str, Any], symmetry_name: str
) -> tuple[Any, dict[str, Any]]:
    """Solve for the ODF by the chosen route, returning it with its evidence."""

    from pytex.core.orientation import OrientationSet
    from pytex.core.symmetry import SymmetrySpec
    from pytex.texture.ghosts import GhostCorrectionSpec
    from pytex.texture.harmonics import HarmonicODF
    from pytex.texture.models import ODF, KernelSpec

    halfwidth = float(request["odf_halfwidth_deg"])
    if str(request["odf_method"]) == "harmonic":
        ghost = str(request["ghost_correction"])
        # The axial group has no finite operator list the harmonic projection
        # could use exactly; the symmetrized data already carry the assumption.
        specimen = (
            SymmetrySpec.specimen(symmetry_name, reference_frame=_specimen_frame())
            if symmetry_name in {"monoclinic", "orthorhombic"}
            else None
        )
        bandlimit = int(request["odf_bandlimit"])
        harmonic = HarmonicODF.invert_pole_figures(
            figures,
            degree_bandlimit=bandlimit,
            regularization=float(request["odf_regularization"]),
            pole_kernel=KernelSpec(halfwidth_deg=halfwidth),
            specimen_symmetry=specimen,
            phi1_step_deg=20.0,
            big_phi_step_deg=20.0,
            phi2_step_deg=20.0,
            ghost_correction=None if ghost == "none" else GhostCorrectionSpec(method=ghost),
        )
        odf = harmonic.final_odf
        ghost_payload = None
        if harmonic.ghost_correction is not None:
            correction = harmonic.ghost_correction
            ghost_payload = {
                "method": correction.method,
                "odd_basis_size": correction.odd_basis_size,
                "amplitude_ratio": correction.ghost_amplitude_ratio,
                "describe": correction.describe(),
            }
        return odf, {
            "method": "harmonic",
            "method_label": (
                f"harmonic series to degree {bandlimit}, {harmonic.basis_size} coefficients"
            ),
            "inversion_residual": float(harmonic.relative_residual_norm),
            "unknowns": int(harmonic.basis_size),
            "observations": int(harmonic.observation_count),
            "texture_index": float(odf.texture_index),
            "ghost": ghost_payload,
        }

    count = int(request["dictionary_count"])
    generator = np.random.default_rng(12345)
    dictionary = OrientationSet.from_euler_angles(
        np.column_stack(
            [
                generator.uniform(0.0, 360.0, size=count),
                np.degrees(np.arccos(generator.uniform(-1.0, 1.0, size=count))),
                generator.uniform(0.0, 360.0, size=count),
            ]
        ),
        specimen_frame=_specimen_frame(),
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    try:
        report = ODF.invert_pole_figures(
            figures, orientation_dictionary=dictionary, kernel=KernelSpec(halfwidth_deg=halfwidth)
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The pole figures could not be inverted: {error}",
            field="files",
            hint="Every figure must share one specimen frame; three planes is the usual minimum.",
        ) from error
    return report.odf, {
        "method": "dictionary",
        "method_label": f"non-negative dictionary of {count} orientations",
        "inversion_residual": float(report.relative_residual_norm),
        "unknowns": count,
        "observations": int(report.observation_count),
        "texture_index": None,
        "ghost": None,
    }


def _section_values(
    request: dict[str, Any], *, kind: str, phase: Any, ranges: dict[str, float]
) -> tuple[list[float] | None, str]:
    """Which sections to draw, and a phrase naming them."""

    preset = str(request["section_preset"])
    maximum = {
        "phi2": ranges["phi2_max_deg"],
        "phi1": ranges["phi1_max_deg"],
        "sigma": ranges["phi1_max_deg"] + ranges["phi2_max_deg"],
    }[kind]
    if preset == "labotex":
        return None, f"every section at {float(request['section_step_deg']):g}° steps"
    if preset == "custom":
        tokens = str(request["section_values"] or "").replace(",", " ").split()
        if not tokens:
            raise InvalidInputError(
                "Chosen section values were asked for, but none were given.",
                field="section_values",
                hint="Type the constant angle of each section, separated by commas: 0, 45, 65.",
            )
        values: list[float] = []
        for token in tokens:
            try:
                value = float(token)
            except ValueError as error:
                raise InvalidInputError(
                    f"{token!r} is not an angle.",
                    field="section_values",
                    hint="Give section values in degrees, separated by commas: 0, 45, 65.",
                ) from error
            if not 0.0 <= value <= maximum + 1e-9:
                raise InvalidInputError(
                    f"A section at {value:g}° lies outside the {maximum:g}° range this crystal "
                    "and sample symmetry need.",
                    field="section_values",
                    hint=f"Every distinct section lies between 0 and {maximum:g} degrees.",
                )
            values.append(value)
        return sorted(set(values)), "the chosen sections"
    if kind == "phi2" and _is_cubic(phase):
        return [0.0, 45.0, 65.0], "the standard cubic sections"
    if kind == "phi2" and math.isclose(maximum, 60.0):
        return [0.0, 30.0], "the standard hexagonal sections"
    step = maximum / 3.0
    return [round(step * index, 6) for index in range(3)], "three evenly spaced sections"


def _residual_summary(measured: np.ndarray, recalculated: np.ndarray) -> dict[str, float]:
    difference = recalculated - measured
    mask = measured >= _RP_THRESHOLD_MRD
    rp = (
        float(np.mean(np.abs(difference[mask]) / measured[mask]) * 100.0)
        if np.any(mask)
        else float("nan")
    )
    return {
        "rp_percent": rp,
        "relative_residual": float(
            np.linalg.norm(difference) / max(float(np.linalg.norm(measured)), 1e-12)
        ),
        "mean_absolute": float(np.mean(np.abs(difference))),
        "max_absolute": float(np.max(np.abs(difference))),
        "rms": float(np.sqrt(np.mean(difference * difference))),
    }


def _analysis_parameters() -> tuple[Any, ...]:
    from pytex.texture.sample_symmetry import SAMPLE_SYMMETRY_NAMES

    assert tuple(option[0] for option in _SAMPLE_SYMMETRY_OPTIONS) == SAMPLE_SYMMETRY_NAMES
    return (
        ObjectParameter(
            name="files",
            label="Pole-figure files",
            help_text=(
                'The opened XRDML files, as `{"items": [{"name": ..., "text": ...}]}`. Supplied by '
                "the **Open pole figures** control. With no file open, a demonstration set "
                "measured from a known model texture is analysed instead."
            ),
            required=False,
        ),
        phase_parameter(
            help_text=(
                "The crystal the figures were measured on. Its symmetry is the crystal symmetry "
                "of the analysis: it folds the ODF, fixes the Euler range of the sections and "
                "chooses the ideal orientations whose volume fractions are reported."
            ),
            builtin="ni_fcc",
        ),
        IndicesListParameter(
            name="poles",
            label="Plane of each file {hkl}",
            help_text=(
                "One plane per file, in the order the files were opened: the file records the "
                "diffraction angle, not the reflection. With fewer lines than files, the last "
                "plane is reused. Three-index Miller indices for every crystal system - "
                "`0 0 2` is the basal plane of a hexagonal metal."
            ),
            default=((1, 1, 1), (2, 0, 0), (2, 2, 0)),
        ),
        ChoiceParameter(
            name="sample_symmetry",
            label="Sample symmetry",
            help_text=(
                "The statistical symmetry the process imposed on the specimen. It is applied to "
                "the measured figures before inversion, and the Euler range of the sections "
                "follows from it.\n\n"
                "**Axial** is the fibre symmetry of drawn, extruded or axially read tubular "
                "product: the texture is the same after any rotation about ND. It is imposed "
                "exactly, as an azimuthal average on each ring of constant tilt.\n\n"
                "Imposing a symmetry the specimen does not have fabricates it. Compare the "
                "measured and symmetrized figures: a large difference between them is the "
                "specimen telling you the assumption is wrong."
            ),
            options=_SAMPLE_SYMMETRY_OPTIONS,
            default="orthorhombic",
        ),
        _PROJECTION_PARAMETER,
        ChoiceParameter(
            name="section_kind",
            label="ODF sections",
            help_text=(
                "Which Euler angle is held constant in each ODF section. Constant phi-2 is the "
                "default of the texture literature and of LaboTex; constant sigma is the view "
                "the gamma fibre of rolled steel is read in."
            ),
            options=_SECTION_KIND_OPTIONS,
            default="phi2",
            group="ODF sections",
        ),
        ChoiceParameter(
            name="section_preset",
            label="Which sections",
            help_text=(
                "The standard sections, every section at a fixed step (the LaboTex plate), or "
                "values you type."
            ),
            options=_SECTION_PRESET_OPTIONS,
            default="standard",
            group="ODF sections",
        ),
        TextParameter(
            name="section_values",
            label="Section values",
            help_text=(
                "The constant angle of each section, in degrees, separated by commas. Used "
                "when Which sections is set to chosen values."
            ),
            required=False,
            default="",
            placeholder="0, 45, 65",
            group="ODF sections",
            field_width="medium",
        ),
        NumberParameter(
            name="section_step_deg",
            label="Section step",
            help_text="Spacing of the sections in the LaboTex plate.",
            units="°",
            default=5.0,
            minimum=2.5,
            maximum=30.0,
            group="ODF sections",
        ),
        NumberParameter(
            name="section_resolution_deg",
            label="Grid in a section",
            help_text="Spacing of the density grid inside each section.",
            units="°",
            default=5.0,
            minimum=2.5,
            maximum=15.0,
            group="ODF sections",
        ),
        NumberParameter(
            name="component_tolerance_deg",
            label="Component tolerance",
            help_text=(
                "Misorientation radius about each ideal orientation that counts as belonging to "
                "it. A fraction means nothing without its tolerance, so it is always reported "
                "with it. Keep it under 30° for hexagonal crystals, whose symmetry-equivalent "
                "balls touch at 30°."
            ),
            units="°",
            default=15.0,
            minimum=5.0,
            maximum=30.0,
            group="Volume fractions",
        ),
        TextParameter(
            name="contour_levels",
            label="Contour levels",
            help_text=(
                "Pole-figure contour levels in m.r.d., separated by commas; `1, 2, 4, 7, 10` is "
                "the conventional sequence. Empty for evenly spaced levels."
            ),
            required=False,
            default="",
            placeholder="1, 2, 4, 7",
            group="Pole figures",
            field_width="medium",
        ),
        ChoiceParameter(
            name="odf_method",
            label="Inversion route",
            help_text=(
                "**Non-negative dictionary** fits weights on a cloud of orientations and cannot "
                "go negative. **Harmonic series** is the Bunge expansion, the route ghost "
                "correction is defined on."
            ),
            options=(
                ("dictionary", "Non-negative dictionary", "Weights on an orientation cloud."),
                ("harmonic", "Harmonic series (Bunge)", "Symmetry-projected coefficients."),
            ),
            default="dictionary",
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        IntegerParameter(
            name="dictionary_count",
            label="Dictionary orientations",
            help_text=(
                "Orientations the dictionary route solves over; more resolves sharper texture."
            ),
            default=800,
            minimum=100,
            maximum=5000,
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        NumberParameter(
            name="odf_halfwidth_deg",
            label="ODF kernel halfwidth",
            help_text=(
                "Width of the bell on each orientation. It is the smoothing and, for an "
                "inversion, the regularisation: too small and the ODF fits the counting noise."
            ),
            units="°",
            default=10.0,
            minimum=2.0,
            maximum=30.0,
            symbol="halfwidth",
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        IntegerParameter(
            name="odf_bandlimit",
            label="Harmonic bandlimit",
            help_text="Highest harmonic degree of the series route. Ignored by the dictionary.",
            default=8,
            minimum=2,
            maximum=16,
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        NumberParameter(
            name="odf_regularization",
            label="Harmonic regularisation",
            help_text="Tikhonov weight on the harmonic coefficients; larger is smoother.",
            default=0.01,
            minimum=1e-8,
            maximum=1.0,
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        ChoiceParameter(
            name="ghost_correction",
            label="Ghost correction",
            help_text=(
                "Recover the odd part a pole figure cannot measure, from positivity. Harmonic "
                "route only; the odd part is an inference, not a measurement."
            ),
            options=(
                ("none", "None (even part only)", "Report what the data determine."),
                ("positivity", "Positivity", "The smallest odd part making the density positive."),
            ),
            default="none",
            group=_INVERSION_GROUP,
            group_collapsed=True,
        ),
        IntegerParameter(
            name="demonstration_seed",
            label="Demonstration seed",
            help_text="Seeds the model texture and counting noise of the demonstration figures.",
            default=11,
            minimum=0,
            maximum=1_000_000,
            group=_INVERSION_GROUP,
            group_collapsed=True,
            field_width="short",
        ),
    )


@REGISTRY.operation(
    "texture.analysis",
    title="Texture analysis of measured pole figures",
    summary=(
        "Load pole figures, the crystal and the sample symmetry once; get measured, recalculated "
        "and difference figures, ODF sections and component volume fractions together."
    ),
    help_text=(
        "The whole of quantitative pole-figure texture analysis from one set of inputs.\n\n"
        "**1. Measured figures.** Each XRDML file is read, assigned its plane, and normalised to "
        "multiples of a random distribution.\n\n"
        "**2. Sample symmetry.** The chosen symmetry is imposed on the measured figures. Axial "
        "symmetry is exact: an azimuthal average on each tilt ring.\n\n"
        "**3. ODF.** The symmetrized figures are inverted into an orientation distribution and "
        "sliced into sections - constant phi-2 by default, or phi-1, or sigma, or every section "
        "at once as a LaboTex plate - over the Euler range the two symmetries require.\n\n"
        "**4. Recalculated and difference figures.** The ODF is projected back onto every "
        "measured direction. The difference, recalculated minus measured, is the test of the "
        "whole analysis: an ODF that cannot reproduce its own data is not a result. Its size is "
        "given as the RP factor, the mean relative deviation where the measurement exceeds "
        "0.5 m.r.d.\n\n"
        "**5. Volume fractions.** The ODF is integrated over a ball of the chosen tolerance about "
        "each ideal orientation of the crystal system, and the fraction is set beside what a "
        "random texture holds in the same ball.\n\n"
        "**Changing a view is cheap.** The inversion is remembered for the same files, crystal, "
        "sample symmetry and inversion settings, so a different section or tolerance returns "
        "without solving again.\n\n"
        "**With no file open** a demonstration set is analysed: figures measured from a known "
        "model texture with 3% counting noise, so every number has an answer to be read against."
    ),
    parameters=_analysis_parameters(),
    returns=(
        "One row per figure with its residual numbers; the figures (measured, symmetrized, "
        "recalculated, difference), the ODF sections and the volume fractions under `data`."
    ),
    panel=_PANEL,
    citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_MATTHIES, _CITATION_LABOTEX),
    tags=(
        "texture", "pole figure", "ODF", "recalculated", "difference", "sample symmetry", "axial",
        "volume fraction", "LaboTex", "sections",
    ),
)
def _texture_analysis(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.texture.components import odf_component_volume_fractions
    from pytex.texture.reconstruction import PoleFigureResidualReport
    from pytex.texture.sample_symmetry import impose_sample_symmetry
    from pytex.texture.sections import euler_section_ranges, odf_sections

    spec, phase = phase_from_request(request["phase"])
    poles = [tuple(int(value) for value in pole) for pole in request["poles"]] or [(1, 1, 1)]
    symmetry_name = str(request["sample_symmetry"])
    items = _file_items(request.get("files"))
    demonstration = not items
    seed = int(request["demonstration_seed"])

    key = hashlib.sha256(
        json.dumps(
            {
                "files": [_item_digest(item) for item in items],
                "phase": spec.to_json(),
                "poles": poles,
                "symmetry": symmetry_name,
                "seed": seed if demonstration else None,
                "inversion": [
                    str(request["odf_method"]),
                    int(request["dictionary_count"]),
                    float(request["odf_halfwidth_deg"]),
                    int(request["odf_bandlimit"]),
                    float(request["odf_regularization"]),
                    str(request["ghost_correction"]),
                ],
            },
            sort_keys=True,
            default=str,
        ).encode()
    ).hexdigest()
    cached = _INVERSION_CACHE.get(key)
    reused = cached is not None
    if cached is None:
        entries = (
            _demonstration_figures(phase, spec, poles, seed)
            if demonstration
            else _read_figures(items, phase, poles)
        )
        symmetrized = [impose_sample_symmetry(entry["figure"], symmetry_name) for entry in entries]
        APP_LOG.info(
            f"Inverting {len(entries)} pole figure(s) with {symmetry_name} sample symmetry.",
            source="texture.analysis",
        )
        odf, evidence = _invert(symmetrized, phase, request, symmetry_name)
        recalculated = [
            np.asarray(
                PoleFigureResidualReport.from_odf(odf, entry["figure"]).predicted_intensities,
                dtype=float,
            )
            for entry in entries
        ]
        cached = {
            "entries": entries,
            "symmetrized": symmetrized,
            "odf": odf,
            "evidence": evidence,
            "recalculated": recalculated,
            "fractions": {},
            "sections": {},
        }
        _INVERSION_CACHE[key] = cached
        while len(_INVERSION_CACHE) > _INVERSION_CACHE_SIZE:
            _INVERSION_CACHE.popitem(last=False)
    else:
        _INVERSION_CACHE.move_to_end(key)

    entries = cached["entries"]
    odf = cached["odf"]
    evidence = cached["evidence"]
    method = str(request["projection"])

    # ---- pole figures: measured, symmetrized, recalculated, difference ----
    figures: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for index, entry in enumerate(entries):
        figure = entry["figure"]
        directions = np.asarray(figure.sample_directions, dtype=float)
        upper = np.where(directions[:, 2:3] < 0.0, -directions, directions)
        projected = _project(upper, method)
        measured = np.asarray(figure.intensities, dtype=float)
        symmetric = np.asarray(cached["symmetrized"][index].intensities, dtype=float)
        recalculated = cached["recalculated"][index]
        difference = recalculated - measured
        polar = np.degrees(np.arccos(np.clip(upper[:, 2], -1.0, 1.0)))
        azimuth = np.degrees(np.arctan2(upper[:, 1], upper[:, 0])) % 360.0
        points = [
            {
                "x": float(projected[point, 0]),
                "y": float(projected[point, 1]),
                "polar_deg": float(polar[point]),
                "azimuth_deg": float(azimuth[point]),
                "measured": float(measured[point]),
                "symmetrized": float(symmetric[point]),
                "recalculated": float(recalculated[point]),
                "difference": float(difference[point]),
                "fit_difference": float(recalculated[point] - symmetric[point]),
            }
            for point in range(directions.shape[0])
        ]
        indices = poles[min(index, len(poles) - 1)]
        label = family_label(indices, spec=spec, family="plane")
        residual = _residual_summary(measured, recalculated)
        # Against the figures that were inverted, the residual measures the fit;
        # against the measurement it measures the fit plus the sample-symmetry
        # assumption. The gap between the two is what the assumption costs.
        fit = _residual_summary(symmetric, recalculated)
        figures.append(
            {
                "label": label,
                "file": entry["file"],
                "sample_label": entry["sample_label"],
                "indices": list(indices),
                "points": points,
                "count": len(points),
                "max_polar_deg": float(polar.max()),
                "ranges": {
                    name: {"minimum": float(values.min()), "maximum": float(values.max())}
                    for name, values in (
                        ("measured", measured),
                        ("symmetrized", symmetric),
                        ("recalculated", recalculated),
                        ("difference", difference),
                    )
                },
                "symmetry_change_rms": float(np.sqrt(np.mean((symmetric - measured) ** 2))),
                "residual": residual,
                "fit_residual": fit,
            }
        )
        summary_rows.append(
            {
                "figure": label,
                "file": entry["file"],
                "points": len(points),
                "measured_max": float(measured.max()),
                "recalculated_max": float(recalculated.max()),
                "rp_percent": residual["rp_percent"],
                "rp_fit_percent": fit["rp_percent"],
                "rms": residual["rms"],
                "symmetry_rms": float(np.sqrt(np.mean((symmetric - measured) ** 2))),
            }
        )

    intensity_min = min(
        min(f["ranges"][name]["minimum"] for name in ("measured", "symmetrized", "recalculated"))
        for f in figures
    )
    intensity_max = max(
        max(f["ranges"][name]["maximum"] for name in ("measured", "symmetrized", "recalculated"))
        for f in figures
    )
    difference_max = max(
        max(abs(f["ranges"]["difference"]["minimum"]), abs(f["ranges"]["difference"]["maximum"]))
        for f in figures
    )
    levels = _contour_levels(
        str(request["contour_levels"] or ""),
        count=6,
        minimum=intensity_min,
        maximum=intensity_max,
    )

    # ---- ODF sections -------------------------------------------------------
    kind = str(request["section_kind"])
    ranges = euler_section_ranges(phase.symmetry, symmetry_name)
    values, section_phrase = _section_values(request, kind=kind, phase=phase, ranges=ranges)
    step = float(request["section_step_deg"])
    resolution = float(request["section_resolution_deg"])
    section_key = json.dumps([kind, values, step, resolution])
    sections = cached["sections"].get(section_key)
    if sections is None:
        sections = odf_sections(
            odf,
            kind=kind,
            values_deg=values,
            ranges=ranges,
            step_deg=step,
            resolution_deg=resolution,
        )
        cached["sections"][section_key] = sections
    densities = np.asarray(sections.densities, dtype=float)
    across = [float(value) for value in sections.phi1_deg]
    big_phi = [float(value) for value in sections.big_phi_deg]
    section_payload = [
        {
            "value_deg": float(value),
            "across_deg": across,
            "big_phi_deg": big_phi,
            "densities": [[float(cell) for cell in row] for row in densities[index]],
            "max_mrd": float(densities[index].max()),
        }
        for index, value in enumerate(np.asarray(sections.phi2_deg, dtype=float))
    ]

    # ---- volume fractions ---------------------------------------------------
    tolerance = float(request["component_tolerance_deg"])
    components = _components_for(phase)
    fractions = cached["fractions"].get(tolerance)
    if fractions is None:
        fractions = (
            odf_component_volume_fractions(
                odf, components, tolerance_deg=tolerance, sample_count=600, seed=0
            )
            if components
            else []
        )
        cached["fractions"][tolerance] = fractions
    fraction_rows = [
        {
            "component": _COMPONENT_TITLES.get(entry["component"], entry["component"]),
            "miller": entry["miller"],
            "fraction": entry["fraction"],
            "percent": 100.0 * entry["fraction"],
            "random_percent": 100.0 * entry["random_fraction"],
            "times_random": entry["times_random"],
        }
        for entry in fractions
    ]

    # ---- prose, stages, result ---------------------------------------------
    symmetry_label = next(
        label for option, label, _ in _SAMPLE_SYMMETRY_OPTIONS if option == symmetry_name
    )
    kind_label = next(label for option, label, _ in _SECTION_KIND_OPTIONS if option == kind)
    names = ", ".join(figure["label"] for figure in figures)
    mean_rp = float(np.nanmean([row["rp_percent"] for row in summary_rows]))
    mean_fit_rp = float(np.nanmean([row["rp_fit_percent"] for row in summary_rows]))
    strongest = max(fraction_rows, key=lambda row: row["fraction"]) if fraction_rows else None
    source = (
        f"demonstration figures measured from {entries[0]['texture_name']}"
        if demonstration
        else f"{len(entries)} XRDML file(s)"
    )
    quality = (
        "reproduces its data well"
        if mean_fit_rp <= 10.0
        else "reproduces its data only approximately"
        if mean_fit_rp <= 25.0
        else "does not reproduce its data"
    )
    summary = (
        f"{len(figures)} pole figure(s) of {spec.name} ({names}) from {source}, analysed with "
        f"{symmetry_label.lower()} sample symmetry. The ODF, a "
        f"{evidence['method_label']}, {quality}: recalculated, it matches the figures it was "
        f"inverted from to an RP factor of {mean_fit_rp:.1f}% and the measured figures to "
        f"{mean_rp:.1f}%, and the gap between the two is what the sample-symmetry assumption "
        f"costs. Its sections are "
        f"{kind_label.lower()} ({section_phrase}), peaking at {float(densities.max()):.2f} m.r.d."
        + (
            f" The strongest ideal orientation within {tolerance:g}° is "
            f"{strongest['component']} at {strongest['percent']:.1f}% of the volume, "
            f"{strongest['times_random']:.2f} times what a random texture holds there."
            if strongest is not None
            else " No ideal-orientation catalogue is defined for this crystal system, so no "
            "volume fractions are reported."
        )
        + (" The inversion was reused from the previous request." if reused else "")
    )

    figure_columns = (
        Column("figure", "Figure"),
        Column("file", "Source"),
        Column("points", "Points", numeric=True),
        Column("measured_max", "Measured max", units="m.r.d.", numeric=True, digits=3),
        Column("recalculated_max", "Recalculated max", units="m.r.d.", numeric=True, digits=3),
        Column(
            "rp_percent",
            "RP",
            units="%",
            numeric=True,
            digits=2,
            help_text=(
                "Mean of |recalculated - measured| / measured over points measured above "
                f"{_RP_THRESHOLD_MRD:g} m.r.d. Below about 10% the ODF reproduces the figure well."
            ),
        ),
        Column(
            "rp_fit_percent",
            "RP to inverted",
            units="%",
            numeric=True,
            digits=2,
            help_text=(
                "The same RP factor against the symmetrized figures the ODF was inverted from: "
                "the quality of the fit alone, without the sample-symmetry assumption."
            ),
        ),
        Column("rms", "RMS difference", units="m.r.d.", numeric=True, digits=4),
        Column(
            "symmetry_rms",
            "Symmetry change",
            units="m.r.d.",
            numeric=True,
            digits=4,
            help_text="RMS change the imposed sample symmetry made to the measured figure.",
        ),
    )
    fraction_columns = (
        Column("component", "Component"),
        Column("miller", "Ideal orientation"),
        Column("percent", "Volume", units="%", numeric=True, digits=2),
        Column(
            "random_percent",
            "Random",
            units="%",
            numeric=True,
            digits=2,
            help_text="What a texture-free specimen holds in the same ball: |G|(w - sin w)/pi.",
        ),
        Column("times_random", "Times random", numeric=True, digits=2),
    )
    stages = (
        ResultStage(
            key="measured",
            title="1. Measured pole figures",
            summary=(
                f"{len(figures)} figure(s), {sum(f['count'] for f in figures)} measured points, "
                f"normalised to m.r.d.; intensities span {intensity_min:.3g} to "
                f"{intensity_max:.3g}."
            ),
            metrics=tuple(
                ResultMetric(
                    f"{figure['label']} tilt reached", figure["max_polar_deg"], units="°"
                )
                for figure in figures
            ),
            explanation=(
                "A figure measured only to a limited tilt constrains the ODF only over that cap; "
                "the rim of an uncorrected figure is where defocusing, not texture, lowers the "
                "counts."
            ),
        ),
        ResultStage(
            key="symmetry",
            title="2. Sample symmetry imposed",
            summary=(
                f"{symmetry_label} sample symmetry changed the measured figures by an RMS of "
                + ", ".join(f"{row['symmetry_rms']:.3g}" for row in summary_rows)
                + " m.r.d."
            ),
            explanation=(
                "A change comparable to the counting noise means the specimen has the assumed "
                "symmetry; a change comparable to the texture means it does not, and the "
                "assumption is fabricating what it averages away."
            ),
        ),
        ResultStage(
            key="inversion",
            title="3. ODF inversion",
            summary=(
                f"A {evidence['method_label']} fitted {evidence['observations']} intensities "
                f"with a relative residual of {evidence['inversion_residual']:.4g}."
            ),
            metrics=(
                ResultMetric("Relative residual", evidence["inversion_residual"]),
                ResultMetric("Unknowns", evidence["unknowns"]),
                ResultMetric("Observations", evidence["observations"]),
                *(
                    (ResultMetric("Texture index", evidence["texture_index"]),)
                    if evidence["texture_index"] is not None
                    else ()
                ),
            ),
            explanation=(
                "Pole-figure inversion is ill-posed: the dictionary or bandlimit, the kernel and "
                "the regularisation all shape the answer. Read the next stage before the ODF."
            ),
            status="ok" if evidence["inversion_residual"] <= 0.25 else "warning",
        ),
        ResultStage(
            key="recalculated",
            title="4. Recalculated and difference figures",
            summary=(
                f"Mean RP factor {mean_fit_rp:.2f}% against the inverted figures and "
                f"{mean_rp:.2f}% against the measured ones; the ODF {quality}."
            ),
            table=ResultTable(columns=figure_columns, rows=tuple(summary_rows)),
            explanation=(
                "The difference figure is recalculated minus measured on the measured directions. "
                "Noise spread over the whole figure is expected; a coherent lobe is a component "
                "the ODF missed, or a systematic error in the measurement."
            ),
            status="ok" if mean_fit_rp <= 25.0 else "warning",
        ),
        ResultStage(
            key="sections",
            title="5. ODF sections",
            summary=(
                f"{len(section_payload)} {kind_label.lower()} section(s) over phi-1 0-"
                f"{ranges['phi1_max_deg']:g}°, Phi 0-{ranges['big_phi_max_deg']:g}°, phi-2 0-"
                f"{ranges['phi2_max_deg']:g}°, the range {spec.name} and {symmetry_label.lower()} "
                "sample symmetry require."
            ),
            metrics=(ResultMetric("Peak density", float(densities.max()), units="m.r.d."),),
        ),
        ResultStage(
            key="fractions",
            title="6. Volume fractions of ideal orientations",
            summary=(
                f"{len(fraction_rows)} ideal orientation(s) within {tolerance:g}°."
                if fraction_rows
                else "No catalogue for this crystal system."
            ),
            table=ResultTable(columns=fraction_columns, rows=tuple(fraction_rows))
            if fraction_rows
            else None,
            explanation=(
                "Each fraction integrates the ODF over a ball about the ideal orientation and its "
                "symmetry equivalents. Balls of neighbouring components can overlap, so the "
                "fractions need not sum to one."
            ),
        ),
    )
    notes = [
        "Pole-figure inversion is ill-posed; the recalculated and difference figures are the "
        "acceptance test of the ODF, not decoration.",
        "The sample symmetry is imposed on the measured figures before inversion. Axial symmetry "
        "is applied exactly, as a trapezoid-weighted azimuthal average on each ring of tilt.",
        "Volume fractions integrate the smoothed ODF over a misorientation ball, estimated from "
        "600 orientations sampled uniformly inside it; they carry a sampling uncertainty of a "
        "few percent of their value.",
        "Defocusing and absorption corrections are not applied to the measured figures.",
    ]
    if len(figures) < 3:
        notes.append(
            f"Only {len(figures)} figure(s) constrain this inversion; three from different planes "
            "is the usual minimum, and with fewer the ODF is largely the regularisation."
        )
    if evidence["ghost"] is not None:
        notes.append(str(evidence["ghost"]["describe"]))

    # The figures that test the ODF, attached to the stages they are evidence for.
    stage_figures: dict[str, tuple[Any, ...]] = {
        "recalculated": (parity_figure(figures), tilt_misfit_figure(figures)),
    }
    if fraction_rows:
        stage_figures["fractions"] = (fraction_figure(fraction_rows, tolerance_deg=tolerance),)
    figured_stages = tuple(
        replace(stage, figures=stage_figures[stage.key]) if stage.key in stage_figures else stage
        for stage in stages
    )

    result = AppResult(
        title=f"Texture analysis of {spec.name}: {names}",
        summary=summary,
        table=ResultTable(
            columns=figure_columns,
            rows=tuple(summary_rows),
            caption="How well the ODF reproduces each measured figure.",
        ),
        data={
            "demonstration": demonstration,
            "cached_inversion": reused,
            "figures": figures,
            "projection": method,
            "specimen_axes": list(_SPECIMEN_AXES),
            "levels": levels,
            "intensity_scale": {"minimum": intensity_min, "maximum": intensity_max},
            "difference_scale": {"maximum": difference_max},
            "sample_symmetry": symmetry_name,
            "sample_symmetry_label": symmetry_label,
            "odf": {
                **{k: v for k, v in evidence.items() if k != "ghost"},
                "ghost": evidence["ghost"],
                "section_kind": kind,
                "section_kind_label": kind_label,
                "horizontal_coordinate": sections.horizontal_coordinate,
                "ranges": ranges,
                "sections": section_payload,
                "max_mrd": float(densities.max()),
                "min_mrd": float(densities.min()),
            },
            "volume_fractions": {"tolerance_deg": tolerance, "rows": fraction_rows},
            "mean_rp_percent": mean_rp,
            "mean_fit_rp_percent": mean_fit_rp,
        },
        inputs={
            "phase": spec.to_json(),
            "poles": [list(pole) for pole in poles],
            "sample_symmetry": symmetry_name,
            "files": [str(item.get("name", "")) for item in items],
            "section_kind": kind,
            "section_preset": str(request["section_preset"]),
            "section_values": values,
            "section_step_deg": step,
            "section_resolution_deg": resolution,
            "component_tolerance_deg": tolerance,
            "odf_method": str(request["odf_method"]),
            "dictionary_count": int(request["dictionary_count"]),
            "odf_halfwidth_deg": float(request["odf_halfwidth_deg"]),
            "odf_bandlimit": int(request["odf_bandlimit"]),
            "odf_regularization": float(request["odf_regularization"]),
            "ghost_correction": str(request["ghost_correction"]),
        },
        notes=tuple(notes),
        citations=(_CITATION_BUNGE, _CITATION_RANDLE, _CITATION_MATTHIES, _CITATION_LABOTEX),
        stages=figured_stages,
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="texture_analysis.example.fcc_rolling",
            title="A rolled fcc sheet, read every way",
            panel=_PANEL,
            summary="Three pole figures of a known rolling texture, orthorhombic sample symmetry.",
            teaches=(
                "Every tab is the same analysis. Start with the difference figures: the RP factor "
                "is a few percent and the difference is speckle, which is what an ODF that "
                "reproduces its data looks like. Then read the phi-2 = 45 and 65 degree sections, "
                "where brass, copper and S sit, and check the volume fractions: the five "
                "components the figures were measured from are the five well above random."
            ),
            operation="texture.analysis",
            request={
                "phase": {"builtin": "ni_fcc"},
                "poles": [[1, 1, 1], [2, 0, 0], [2, 2, 0]],
                "sample_symmetry": "orthorhombic",
                "dictionary_count": 600,
            },
        ),
        ExampleScenario(
            id="texture_analysis.example.axial_zirconium",
            title="Axial symmetry on a zirconium tube",
            panel=_PANEL,
            summary="A split-basal texture with axial sample symmetry imposed about ND.",
            teaches=(
                "The figures were measured from basal poles tilted 30 degrees towards TD, which "
                "is not axial. Imposing axial symmetry turns the two lobes into a ring: compare "
                "the measured and symmetrized figures, then the difference, which now shows the "
                "two lobes the assumption averaged away. The symmetry change in stage 2 is as "
                "large as the texture, which is how the analysis says the assumption is wrong for "
                "this specimen."
            ),
            operation="texture.analysis",
            request={
                "phase": {"builtin": "zr_hcp"},
                "poles": [[0, 0, 2], [1, 0, 0], [1, 0, 1]],
                "sample_symmetry": "axial",
                "dictionary_count": 600,
            },
        ),
        ExampleScenario(
            id="texture_analysis.example.labotex_plate",
            title="Every phi-2 section, LaboTex style",
            panel=_PANEL,
            summary="The rolled fcc demonstration as a plate of every 5 degree phi-2 section.",
            teaches=(
                "The standard three sections are a choice made for cubic rolling textures. The "
                "plate shows all nineteen, 0 to 90 degrees, on one scale, so a component that sits "
                "between the standard sections is not missed. Watch the brass and S maxima move "
                "through the sections: an ODF is a density in three dimensions, and a section is "
                "one slice of it."
            ),
            operation="texture.analysis",
            request={
                "phase": {"builtin": "ni_fcc"},
                "poles": [[1, 1, 1], [2, 0, 0], [2, 2, 0]],
                "sample_symmetry": "orthorhombic",
                "section_preset": "labotex",
                "section_resolution_deg": 7.5,
                "dictionary_count": 600,
            },
        ),
    )
)
