from __future__ import annotations

import json
import warnings
from collections import OrderedDict, defaultdict
from importlib.metadata import version
from pathlib import Path

import numpy as np
from diffpy.structure import loadStructure
from diffsims.generators.diffraction_generator import DiffractionGenerator
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure

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


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _canonical_family(indices: tuple[int, int, int] | np.ndarray) -> tuple[int, int, int]:
    return tuple(sorted((abs(int(value)) for value in indices), reverse=True))


def _make_crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def _xrd_baseline(repo_root: Path) -> dict[str, object]:
    fixture = get_phase_fixture("ni_fcc")
    fixture_path = repo_root / fixture.artifact_path.relative_to(repo_root)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        structure = Structure.from_file(fixture_path)
        pattern = XRDCalculator(wavelength="CuKa").get_pattern(
            structure,
            scaled=True,
            two_theta_range=(20.0, 120.0),
        )
    peaks: list[dict[str, object]] = []
    for two_theta_deg, relative_intensity, hkl_rows in zip(
        pattern.x,
        pattern.y,
        pattern.hkls,
        strict=True,
    ):
        representative = tuple(int(value) for value in hkl_rows[0]["hkl"])
        peaks.append(
            {
                "representative_hkl": list(_canonical_family(representative)),
                "multiplicity": int(hkl_rows[0]["multiplicity"]),
                "two_theta_deg": round(float(two_theta_deg), 6),
                "relative_intensity": round(float(relative_intensity), 6),
            }
        )

    crystal = _make_crystal_frame()
    phase = fixture.load_phase(crystal_frame=crystal)
    families: OrderedDict[tuple[int, int, int], dict[str, object]] = OrderedDict()
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
                "representative_hkl": list(family),
                "two_theta_deg": round(float(reflection.two_theta_deg), 6),
                "multiplicity": int(reflection.multiplicity),
            },
        )

    return {
        "schema_id": "pytex.diffraction_external_baseline.xrd",
        "schema_version": "1.0.0",
        "fixture_id": fixture.fixture_id,
        "phase_name": fixture.phase_name,
        "source_kind": "open_source_reference_implementation",
        "source_package": {
            "name": "pymatgen",
            "version": version("pymatgen"),
            "calculator": "XRDCalculator",
        },
        "radiation": "CuKa",
        "two_theta_range_deg": [20.0, 120.0],
        "reference_peaks": peaks,
        "pytex_reference_families": list(families.values()),
        "notes": [
            "Peak positions and multiplicities come from pymatgen's "
            "XRDCalculator using the pinned ni_fcc fixture.",
            "PyTex comparisons should treat peak-position agreement as the hard "
            "external baseline and intensity differences as informative only.",
        ],
    }


def _saed_baseline(repo_root: Path) -> dict[str, object]:
    fixture = get_phase_fixture("ni_fcc")
    fixture_path = repo_root / fixture.artifact_path.relative_to(repo_root)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="CUDA Toolkit .* unsupported.*")
        simulation = DiffractionGenerator(accelerating_voltage=200).calculate_ed_data(
            loadStructure(str(fixture_path)),
            reciprocal_radius=1.0,
            with_direct_beam=False,
        )
    spot_shells: defaultdict[tuple[int, int, int], list[dict[str, object]]] = defaultdict(list)
    for indices, coordinates, intensity in zip(
        simulation.indices,
        simulation.coordinates,
        simulation.intensities,
        strict=True,
    ):
        family = _canonical_family(tuple(int(value) for value in indices))
        spot_shells[family].append(
            {
                "miller_indices": [int(value) for value in indices],
                "reciprocal_coordinates_inv_angstrom": [
                    round(float(coordinates[0]), 6),
                    round(float(coordinates[1]), 6),
                ],
                "relative_intensity": round(float(intensity), 6),
            }
        )

    crystal = _make_crystal_frame()
    phase = fixture.load_phase(crystal_frame=crystal)
    pytex_pattern = generate_saed_pattern(
        phase,
        ZoneAxis(np.array([0, 0, 1]), phase=phase),
        camera_constant_mm_angstrom=1.0,
        max_index=4,
        max_g_inv_angstrom=1.0,
    )
    pytex_shells: dict[tuple[int, int, int], dict[str, object]] = {}
    for spot in pytex_pattern.spots:
        if spot.intensity <= 1.0:
            continue
        family = _canonical_family(spot.miller_indices)
        shell = pytex_shells.setdefault(
            family,
            {
                "representative_hkl": list(family),
                "spot_count": 0,
                "radius_inv_angstrom": round(float(np.linalg.norm(spot.detector_coordinates)), 6),
            },
        )
        shell["spot_count"] = int(shell["spot_count"]) + 1

    baseline_shells: list[dict[str, object]] = []
    for family, spots in sorted(spot_shells.items()):
        radii = [
            float(
                np.linalg.norm(np.asarray(spot["reciprocal_coordinates_inv_angstrom"], dtype=float))
            )
            for spot in spots
        ]
        baseline_shells.append(
            {
                "representative_hkl": list(family),
                "spot_count": len(spots),
                "radius_inv_angstrom": round(float(np.mean(radii)), 6),
                "reference_spots": spots,
            }
        )

    return {
        "schema_id": "pytex.diffraction_external_baseline.saed",
        "schema_version": "1.0.0",
        "fixture_id": fixture.fixture_id,
        "phase_name": fixture.phase_name,
        "zone_axis": [0, 0, 1],
        "accelerating_voltage_kev": 200.0,
        "reciprocal_radius_inv_angstrom": 1.0,
        "source_kind": "open_source_reference_implementation",
        "source_package": {
            "name": "diffsims",
            "version": version("diffsims"),
            "calculator": "DiffractionGenerator",
        },
        "reference_shells": baseline_shells,
        "pytex_reference_shells": list(pytex_shells.values()),
        "notes": [
            "The diffsims baseline records shell radii and indexed spots for "
            "the ni_fcc [001] SAED case.",
            "PyTex comparisons should treat shell geometry and indexed-family "
            "coverage as the hard external baseline while allowing in-plane "
            "detector-basis rotations.",
        ],
    }


def main() -> int:
    repo_root = _repo_root()
    output_root = repo_root / "fixtures/diffraction"
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "ni_fcc_pymatgen_xrd_cuka.json").write_text(
        json.dumps(_xrd_baseline(repo_root), indent=2) + "\n",
        encoding="utf-8",
    )
    (output_root / "ni_fcc_diffsims_saed_001_200kev.json").write_text(
        json.dumps(_saed_baseline(repo_root), indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote diffraction baseline artifacts to {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
