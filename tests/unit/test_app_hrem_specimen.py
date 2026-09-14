"""The HRTEM specimen: thickness, imported .xyz structures, and the viewer payload.

The expectations come from geometry, not from recorded outputs: a slab of n unit
cells along [001] is n·c thick, a file's atoms are the atoms simulated, and a PNG
written at native resolution has one pixel per simulated pixel.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

import numpy as np
import pytest
from matplotlib.image import imread

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot

_SIMULATE = "tem.simulate_hrem"

_WATER_LIKE = "3\nsmall cluster\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"


def _simulate(**overrides: object) -> dict:
    request: dict[str, object] = {
        "phase": {"builtin": "si_diamond"},
        "sample_type": "crystalline",
        "zone_axis": [0, 0, 1],
        "supercell_xy": 1,
        "thickness_angstrom": 10.0,
        "mode": "double_corrected",
        "sampling_angstrom": 0.25,
    }
    request.update(overrides)
    return REGISTRY.call(_SIMULATE, request)


def _silicon():  # type: ignore[no-untyped-def]
    return phase_from_request({"builtin": "si_diamond"})[1]


def _png_shape(data_url: str) -> tuple[int, int]:
    raw = base64.b64decode(data_url.split(",", 1)[1])
    image = imread(io.BytesIO(raw), format="png")
    return int(image.shape[0]), int(image.shape[1])


class TestThickness:
    def test_a_zone_axis_slab_is_a_whole_number_of_periods_at_or_above_the_request(self) -> None:
        phase = _silicon()
        c = float(phase.lattice.c)
        repeats, slab = AtomicSnapshot.repeats_for_thickness(phase, (0, 0, 1), 30.0, 2)
        assert repeats == (2, 2, int(np.ceil(30.0 / c)))
        assert slab == pytest.approx(repeats[2] * c)
        assert slab >= 30.0

    def test_the_thickness_control_reaches_the_specimen(self) -> None:
        thin = _simulate(thickness_angstrom=6.0)
        thick = _simulate(thickness_angstrom=40.0)
        thickness = "specimen_thickness_angstrom"
        assert thick["data"][thickness] > thin["data"][thickness]
        assert thick["data"]["natoms"] > thin["data"]["natoms"]
        assert thick["data"]["unit_cell_repeats"][2] > thin["data"]["unit_cell_repeats"][2]

    def test_the_delivered_thickness_is_reported_beside_the_request(self) -> None:
        result = _simulate(thickness_angstrom=12.0)
        assert result["data"]["requested_thickness_angstrom"] == 12.0
        assert result["data"]["specimen_thickness_angstrom"] >= 12.0
        assert "thick along the beam" in result["summary"]

    def test_the_amorphous_foil_takes_the_thickness_as_its_depth(self) -> None:
        result = _simulate(sample_type="amorphous", thickness_angstrom=12.0)
        assert result["data"]["specimen_thickness_angstrom"] == pytest.approx(12.0)

    def test_an_oversized_specimen_is_refused_beside_the_thickness(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            _simulate(supercell_xy=10, thickness_angstrom=500.0)
        assert caught.value.details["field"] == "thickness_angstrom"


class TestImportedStructure:
    def test_an_xyz_file_is_the_specimen_that_is_simulated(self) -> None:
        result = _simulate(
            sample_type="imported",
            structure_file={"name": "cluster.xyz", "text": _WATER_LIKE},
        )
        assert result["data"]["natoms"] == 3
        assert result["data"]["structure_file"]["name"] == "cluster.xyz"
        assert "cluster.xyz" in result["summary"]
        # The provenance travels with the result, not a second copy of the text.
        assert "text" not in result["inputs"]["structure_file"]

    def test_extended_xyz_lattice_and_column_layout_are_read(self) -> None:
        text = (
            "2\n"
            'Lattice="8 0 0 0 9 0 0 0 10" Properties=id:I:1:species:S:1:pos:R:3 relaxed\n'
            "1 Au 1.0 2.0 3.0\n"
            "2 Au 5.0 6.0 7.0\n"
        )
        snapshot = AtomicSnapshot.from_xyz(text)
        assert snapshot.species == ("Au", "Au")
        np.testing.assert_allclose(snapshot.cell, np.diag([8.0, 9.0, 10.0]))
        np.testing.assert_allclose(snapshot.positions[1], [5.0, 6.0, 7.0])
        assert snapshot.label == "relaxed"

    def test_a_periodic_cell_wraps_atoms_and_a_cluster_is_boxed_with_a_margin(self) -> None:
        periodic = AtomicSnapshot.from_xyz(
            '1\nLattice="10 0 0 0 10 0 0 0 10"\nC -1.0 11.0 4.0\n'
        ).prepared_for_imaging(periodic_xy=True)
        np.testing.assert_allclose(periodic.positions[0, :2], [9.0, 1.0])
        cluster = AtomicSnapshot.from_xyz(_WATER_LIKE).prepared_for_imaging(periodic_xy=False)
        assert np.all(cluster.positions >= 1.0 - 1e-12)
        assert np.all(cluster.positions <= np.diag(cluster.cell) - 1.0 + 1e-12)

    def test_a_sheared_cell_is_refused_rather_than_squared(self) -> None:
        sheared = AtomicSnapshot.from_xyz('1\nLattice="10 0 0 3 10 0 0 0 10"\nC 1 1 1\n')
        with pytest.raises(ValueError, match="orthogonal"):
            sheared.prepared_for_imaging(periodic_xy=True)

    def test_a_single_line_upload_is_never_read_as_a_server_path(self, tmp_path: Path) -> None:
        secret = tmp_path / "secret.xyz"
        secret.write_text("1\nsecret\nC 0 0 0\n", encoding="utf-8")
        with pytest.raises(InvalidInputError) as caught:
            _simulate(sample_type="imported", structure_file={"name": "a.xyz", "text": str(secret)})
        assert caught.value.details["field"] == "structure_file"

    def test_species_that_are_not_elements_are_named(self) -> None:
        with pytest.raises(InvalidInputError, match="not elements: Zz"):
            _simulate(
                sample_type="imported",
                structure_file={"name": "typed.xyz", "text": "1\n\nZz 0 0 0\n"},
            )

    def test_the_imported_specimen_needs_a_file(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            _simulate(sample_type="imported")
        assert caught.value.details["field"] == "structure_file"

    def test_the_wrong_file_kind_is_refused(self) -> None:
        with pytest.raises(InvalidInputError, match="extension"):
            _simulate(
                sample_type="imported",
                structure_file={"name": "cell.cif", "text": _WATER_LIKE},
            )


class TestViewerPayload:
    def test_the_images_hold_one_pixel_per_simulated_pixel(self) -> None:
        result = _simulate()
        rows, columns = result["data"]["image_shape_px"]
        assert _png_shape(result["data"]["image_png"]) == (rows, columns)
        assert _png_shape(result["data"]["power_spectrum_png"]) == (rows, columns)

    def test_the_nyquist_extent_follows_the_sampling(self) -> None:
        data = _simulate()["data"]
        rows, columns = data["image_shape_px"]
        lx, ly = data["extent_angstrom"]
        assert data["nyquist_inv_angstrom"] == pytest.approx([0.5 * columns / lx, 0.5 * rows / ly])

    def test_the_panel_draws_the_figures_as_zoomable_svg(self) -> None:
        source = (
            Path(__file__).resolve().parents[2] / "src/pytex/app/static/js/panels/hrem.js"
        ).read_text(encoding="utf-8")
        assert "img.hrem-image" not in source
        assert "svg('image'" in source
        assert "structure_file" in source
