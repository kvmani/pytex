"""The workbench's multislice HRTEM: engines, periodic specimens, and focal/thickness series.

Expected values come from geometry and from the engine's own contracts (a periodic
slab is a whole number of lattice periods; a series of R thicknesses by D defoci
has R x D images), never from recorded outputs.
"""

from __future__ import annotations

import pytest

from pytex.adapters.abtem import is_abtem_available
from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

_SIMULATE = "tem.simulate_hrem"
_SERIES = "tem.hrtem_series"
_A_SI = 5.43102


def _request(**overrides: object) -> dict[str, object]:
    request: dict[str, object] = {
        "phase": {"builtin": "si_diamond"},
        "sample_type": "crystalline",
        "zone_axis": [0, 0, 1],
        "supercell_xy": 2,
        "thickness_angstrom": 20.0,
        "mode": "double_corrected",
        "sampling_angstrom": 0.2,
    }
    request.update(overrides)
    return request


def _rows(result: dict) -> dict[str, str]:
    return {row["metric"]: row["value"] for row in result["table"]["rows"]}


class TestEngines:
    def test_the_default_engine_is_the_pytex_multislice(self) -> None:
        result = REGISTRY.call(_SIMULATE, _request())
        rows = _rows(result)
        assert result["data"]["engine"] == "multislice"
        assert rows["Simulation engine"] == "PyTex multislice"
        assert result["data"]["multislice"]["num_slices"] >= 1
        assert "PyTex multislice" in result["notes"][0]

    def test_the_phase_object_says_what_it_leaves_out(self) -> None:
        result = REGISTRY.call(_SIMULATE, _request(engine="phase_object"))
        assert result["data"]["engine"] == "phase_object"
        assert result["data"]["multislice"] is None
        assert "no propagation" in result["notes"][0]

    @pytest.mark.skipif(is_abtem_available(), reason="checks the refusal without abTEM")
    def test_abtem_is_refused_where_it_is_not_installed(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            REGISTRY.call(_SIMULATE, _request(engine="abtem"))
        assert caught.value.details["field"] == "engine"

    def test_the_multislice_diagnostics_reach_the_table(self) -> None:
        rows = _rows(REGISTRY.call(_SIMULATE, _request(slice_thickness_angstrom=_A_SI / 4)))
        assert rows["Slices (distinct potentials)"].startswith("16 (4)")
        assert float(rows["Retained intensity at the exit surface"]) <= 1.0

    def test_tilt_and_frozen_phonons_are_reported(self) -> None:
        result = REGISTRY.call(
            _SIMULATE,
            _request(tilt_x_mrad=5.0, frozen_phonons=3, thermal_rms_angstrom=0.08),
        )
        rows = _rows(result)
        assert rows["Beam tilt (x, y)"] == "5.00, 0.00"
        assert rows["Frozen-phonon configurations"] == "3"

    def test_an_oversized_grid_is_refused_beside_the_sampling(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            REGISTRY.call(_SIMULATE, _request(supercell_xy=10, sampling_angstrom=0.05))
        assert caught.value.details["field"] == "sampling_angstrom"


class TestPeriodicSpecimen:
    def test_a_crystal_box_is_a_whole_number_of_lattice_periods(self) -> None:
        result = REGISTRY.call(_SIMULATE, _request(zone_axis=[1, 1, 0], supercell_xy=2))
        lx, ly = result["data"]["extent_angstrom"]
        assert lx == pytest.approx(2 * _A_SI)
        assert ly == pytest.approx(2 * _A_SI * 2**0.5)
        assert any("Orthogonal periodic cell" in note for note in result["notes"])

    def test_the_vacancy_specimen_lacks_exactly_one_atom(self) -> None:
        perfect = REGISTRY.call(_SIMULATE, _request())
        vacancy = REGISTRY.call(_SIMULATE, _request(sample_type="vacancy"))
        assert vacancy["data"]["natoms"] == perfect["data"]["natoms"] - 1


class TestSeries:
    def test_a_tableau_has_one_image_per_thickness_and_defocus(self) -> None:
        result = REGISTRY.call(
            _SERIES,
            _request(
                zone_axis=[1, 1, 0],
                supercell_xy=1,
                thickness_angstrom=60.0,
                thickness_steps=3,
                defocus_start_angstrom=-100.0,
                defocus_stop_angstrom=100.0,
                defocus_step_angstrom=50.0,
            ),
        )
        data = result["data"]
        assert data["defoci_angstrom"] == [-100.0, -50.0, 0.0, 50.0, 100.0]
        assert len(data["thicknesses_angstrom"]) == 3
        assert len(data["tiles"]) == 15
        assert len(data["contrasts"]) == 3 and len(data["contrasts"][0]) == 5
        assert data["thicknesses_angstrom"][-1] == pytest.approx(
            data["beam_thickness_angstrom"][-1]
        )
        assert [figure["key"] for figure in result["figures"]] == [
            "hrtem_defocus_thickness_tableau",
            "hrtem_focal_contrast",
            "hrtem_beam_intensities_with_thickness",
        ]

    def test_the_beams_are_the_transmitted_one_and_the_strongest_diffracted(self) -> None:
        result = REGISTRY.call(
            _SERIES, _request(zone_axis=[1, 1, 0], supercell_xy=1, thickness_angstrom=60.0)
        )
        labels = [beam["label"] for beam in result["data"]["beams"]]
        assert labels[0] == "(000)"
        assert len(labels) == 5
        # The [110] zone's strongest reflections are {111}, {220} and {004}.
        assert "(1 -1 1)" in labels or "(-1 1 1)" in labels
        for beam in result["data"]["beams"]:
            assert len(beam["intensity"]) == len(result["data"]["beam_thickness_angstrom"])

    def test_an_amorphous_series_has_no_bragg_beams(self) -> None:
        result = REGISTRY.call(_SERIES, _request(sample_type="amorphous", thickness_steps=1))
        assert result["data"]["beams"] == []
        assert len(result["figures"]) == 2

    def test_a_reversed_focal_range_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            REGISTRY.call(
                _SERIES, _request(defocus_start_angstrom=100.0, defocus_stop_angstrom=-100.0)
            )
        assert caught.value.details["field"] == "defocus_stop_angstrom"

    def test_too_many_images_are_refused(self) -> None:
        with pytest.raises(InvalidInputError) as caught:
            REGISTRY.call(
                _SERIES,
                _request(
                    thickness_steps=8,
                    defocus_start_angstrom=-600.0,
                    defocus_stop_angstrom=600.0,
                    defocus_step_angstrom=50.0,
                ),
            )
        assert caught.value.details["field"] == "thickness_steps"


def _decode(data_url: str) -> bytes:
    import base64

    return base64.b64decode(data_url.split(",", 1)[1])


class TestFullResolutionDownloads:
    def test_the_micrograph_tiff_is_the_computed_intensity(self) -> None:
        import io

        import numpy as np
        from PIL import Image

        result = REGISTRY.call(_SIMULATE, _request(sampling_angstrom=0.2))
        tiff = result["data"]["image_tiff"]
        values = np.asarray(Image.open(io.BytesIO(_decode(tiff["data"]))))
        rows, columns = result["data"]["image_shape_px"]
        assert values.dtype == np.float32
        assert values.shape == (rows, columns) == (tiff["height"], tiff["width"])
        assert tiff["filename"].endswith("-float32.tif")
        # The file holds the image intensity itself, the same way up as the
        # displayed PNG: that PNG is the same array through an 8-bit grey map, so
        # the TIFF rescaled to 0-255 must match it everywhere within the two levels
        # that matplotlib's binned colour lookup can differ from rounding.
        png = np.asarray(Image.open(io.BytesIO(_decode(result["data"]["image_png"]))))
        grey = (png[..., 0] if png.ndim == 3 else png).astype(float)
        scaled = 255.0 * (values - values.min()) / (values.max() - values.min())
        assert float(np.max(np.abs(scaled - grey))) <= 2.0
        assert np.isfinite(values).all() and float(values.min()) >= 0.0

    def test_the_spectrum_tiff_matches_the_spectrum_png_size(self) -> None:
        result = REGISTRY.call(_SIMULATE, _request())
        tiff = result["data"]["power_spectrum_tiff"]
        rows, columns = result["data"]["image_shape_px"]
        assert (tiff["height"], tiff["width"]) == (rows, columns)
        assert "log10" in tiff["label"]

    def test_the_series_zip_holds_every_image_and_its_index(self) -> None:
        import csv
        import io
        import zipfile

        import numpy as np
        from PIL import Image

        result = REGISTRY.call(
            _SERIES,
            _request(
                zone_axis=[1, 1, 0], supercell_xy=1, thickness_angstrom=40.0, thickness_steps=2,
                defocus_start_angstrom=-50.0, defocus_stop_angstrom=50.0,
                defocus_step_angstrom=50.0,
            ),
        )
        tableau, bundle = result["data"]["downloads"]
        rows, columns = result["data"]["image_shape_px"]
        montage = Image.open(io.BytesIO(_decode(tableau["data"])))
        assert montage.size == (3 * columns + 2 * 4, 2 * rows + 4)
        archive = zipfile.ZipFile(io.BytesIO(_decode(bundle["data"])))
        names = archive.namelist()
        assert {"tableau.png", "series.csv", "README.txt"} <= set(names)
        assert sum(name.endswith(".tif") for name in names) == 6
        assert sum(name.endswith(".png") for name in names) == 7
        index = list(csv.DictReader(io.StringIO(archive.read("series.csv").decode())))
        assert len(index) == 6
        first = np.asarray(Image.open(io.BytesIO(archive.read(index[0]["file"]))))
        assert first.shape == (rows, columns) and first.dtype == np.float32
        assert float(first.mean()) == pytest.approx(float(index[0]["mean_intensity"]), rel=1e-5)
