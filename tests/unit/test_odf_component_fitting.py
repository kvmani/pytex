"""Named-component ODF fitting against an exactly constructed mixture."""

from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    ODF,
    FrameDomain,
    Handedness,
    KernelSpec,
    OrientationSet,
    ReferenceFrame,
    SymmetrySpec,
    TextureComponent,
    fit_odf_components,
    from_json_contract,
    to_json_contract,
)


def _exact_cube_goss_odf() -> ODF:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame(
        "specimen", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    # The first two support points carry an exactly known 70/30 mixture. The
    # remaining zero-weight points are an explicit evaluation dictionary that
    # makes the component-plus-random design identifiable.
    euler = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 45.0, 0.0],
            [15.0, 20.0, 10.0],
            [30.0, 55.0, 15.0],
            [70.0, 35.0, 40.0],
            [10.0, 75.0, 80.0],
        ]
    )
    orientations = OrientationSet.from_euler_angles(
        euler,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        degrees=True,
    )
    return ODF.from_orientations(
        orientations,
        weights=np.array([0.7, 0.3, 0.0, 0.0, 0.0, 0.0]),
        kernel=KernelSpec(halfwidth_deg=12.0),
    )


def _components() -> tuple[TextureComponent, TextureComponent]:
    return (
        TextureComponent("cube", (0.0, 0.0, 0.0), "{001}<100>"),
        TextureComponent("goss", (0.0, 45.0, 0.0), "{011}<100>"),
    )


def test_exact_named_mixture_is_recovered_with_zero_residual() -> None:
    fit = fit_odf_components(_exact_cube_goss_odf(), _components())

    assert fit.fraction_for("cube") == pytest.approx(0.7, abs=1e-10)
    assert fit.fraction_for("goss") == pytest.approx(0.3, abs=1e-10)
    assert fit.random_fraction == pytest.approx(0.0, abs=1e-10)
    assert fit.rms_residual == pytest.approx(0.0, abs=1e-10)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-12)
    assert "Fractions describe only this declared component basis" in fit.describe()


def test_component_fit_json_contract_round_trip() -> None:
    fit = fit_odf_components(_exact_cube_goss_odf(), _components())
    payload = to_json_contract(fit)
    restored = from_json_contract(payload)

    assert payload["schema_id"] == "pytex.texture.odf_component_fit"
    assert restored.fractions == pytest.approx(fit.fractions)
    assert restored.predicted_density == pytest.approx(fit.predicted_density)
    assert restored.describe() == fit.describe()


def test_rank_deficient_evaluation_support_raises() -> None:
    odf = _exact_cube_goss_odf()
    sparse = odf.orientations[:2]

    with pytest.raises(ValueError, match="rank deficient"):
        fit_odf_components(odf, _components(), evaluation_orientations=sparse)


def test_duplicate_component_names_raise() -> None:
    cube = TextureComponent("cube", (0.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="unique"):
        fit_odf_components(_exact_cube_goss_odf(), (cube, cube))
