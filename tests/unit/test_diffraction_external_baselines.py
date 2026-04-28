from __future__ import annotations

import json
import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pytest

from pytex import (
    FrameDomain,
    Handedness,
    RadiationSpec,
    ReferenceFrame,
    ZoneAxis,
    generate_powder_reflections,
    generate_saed_pattern,
    get_phase_fixture,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
BASELINE_ROOT = REPO_ROOT / "fixtures" / "diffraction"


def _canonical_family(indices: tuple[int, int, int] | np.ndarray) -> tuple[int, int, int]:
    return tuple(sorted((abs(int(value)) for value in indices), reverse=True))


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


@pytest.mark.parametrize("fixture_id", ["ni_fcc", "fe_bcc"])
def test_xrd_external_baseline_matches_pytex_peak_families(fixture_id: str) -> None:
    pytest.importorskip(
        "pymatgen.core",
        reason="CIF-backed diffraction baseline tests require the optional pymatgen dependency.",
    )
    payload = json.loads(
        (BASELINE_ROOT / f"{fixture_id}_pymatgen_xrd_cuka.json").read_text(encoding="utf-8")
    )
    assert payload["schema_id"] == "pytex.diffraction_external_baseline.xrd"
    assert payload["fixture_id"] == fixture_id

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        phase = get_phase_fixture(fixture_id).load_phase(crystal_frame=_crystal_frame())

    families: OrderedDict[tuple[int, int, int], dict[str, float | int]] = OrderedDict()
    for reflection in generate_powder_reflections(
        phase,
        radiation=RadiationSpec.cu_ka(),
        two_theta_range_deg=(20.0, 120.0),
        max_index=6,
    ):
        family = _canonical_family(reflection.miller_indices)
        families.setdefault(
            family,
            {
                "two_theta_deg": float(reflection.two_theta_deg),
                "multiplicity": int(reflection.multiplicity),
            },
        )

    expected_families = [
        tuple(row["representative_hkl"]) for row in payload["reference_peaks"]
    ]
    assert list(families)[: len(expected_families)] == expected_families
    for expected_row in payload["reference_peaks"]:
        family = tuple(expected_row["representative_hkl"])
        assert family in families
        candidate = families[family]
        assert candidate["multiplicity"] == expected_row["multiplicity"]
        assert abs(candidate["two_theta_deg"] - expected_row["two_theta_deg"]) <= 0.15

    assert payload["pytex_reference_families"] == [
        {
            "representative_hkl": list(family),
            "two_theta_deg": round(float(values["two_theta_deg"]), 6),
            "multiplicity": int(values["multiplicity"]),
        }
        for family, values in list(families.items())[: len(payload["pytex_reference_families"])]
    ]


@pytest.mark.parametrize("fixture_id", ["ni_fcc", "fe_bcc"])
def test_saed_external_baseline_matches_pytex_shell_geometry(fixture_id: str) -> None:
    pytest.importorskip(
        "pymatgen.core",
        reason="CIF-backed diffraction baseline tests require the optional pymatgen dependency.",
    )
    payload = json.loads(
        (BASELINE_ROOT / f"{fixture_id}_diffsims_saed_001_200kev.json").read_text(encoding="utf-8")
    )
    assert payload["schema_id"] == "pytex.diffraction_external_baseline.saed"
    assert payload["fixture_id"] == fixture_id
    assert payload["zone_axis"] == [0, 0, 1]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        phase = get_phase_fixture(fixture_id).load_phase(crystal_frame=_crystal_frame())

    zone_axis = ZoneAxis(np.array(payload["zone_axis"], dtype=np.int64), phase=phase)
    pattern = generate_saed_pattern(
        phase,
        zone_axis,
        camera_constant_mm_angstrom=1.0,
        max_index=4,
        max_g_inv_angstrom=float(payload["reciprocal_radius_inv_angstrom"]),
    )
    shell_rows: OrderedDict[tuple[int, int, int], dict[str, float | int]] = OrderedDict()
    for spot in pattern.spots:
        if spot.intensity <= 1.0:
            continue
        family = _canonical_family(spot.miller_indices)
        row = shell_rows.setdefault(
            family,
            {
                "spot_count": 0,
                "radius_inv_angstrom": float(np.linalg.norm(spot.detector_coordinates)),
            },
        )
        row["spot_count"] = int(row["spot_count"]) + 1

    expected_shells = payload["reference_shells"]
    expected_families = [tuple(row["representative_hkl"]) for row in expected_shells]
    assert list(shell_rows)[: len(expected_families)] == expected_families
    for expected_row in expected_shells:
        family = tuple(expected_row["representative_hkl"])
        assert family in shell_rows
        candidate = shell_rows[family]
        assert candidate["spot_count"] == expected_row["spot_count"]
        assert abs(candidate["radius_inv_angstrom"] - expected_row["radius_inv_angstrom"]) <= 1e-6

    assert payload["pytex_reference_shells"] == [
        {
            "representative_hkl": list(family),
            "spot_count": int(values["spot_count"]),
            "radius_inv_angstrom": round(float(values["radius_inv_angstrom"]), 6),
        }
        for family, values in shell_rows.items()
    ]
