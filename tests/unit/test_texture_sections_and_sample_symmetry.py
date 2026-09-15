"""Sample symmetry, ODF sections and ODF volume fractions, against fixed answers.

Every expectation here is fixed before the code runs:

* Imposing a symmetry on a figure that already has it changes nothing, and
  imposing it on one that lacks it removes exactly the terms the group forbids
  (an orthorhombic average kills anything odd in x or y; an axial average leaves
  the polar dependence alone).
* The Euler box a section must cover is set by the operators: 90 degrees of
  phi2 for cubic, 60 for hexagonal, and the specimen symmetry folds phi1.
* A random texture is 1 m.r.d. on every section and holds the Haar-measure
  fraction ``|G| (w - sin w) / pi`` of orientation space within ``w`` of any
  ideal orientation.
* A constant-phi1 section and a constant-phi2 section cross along a line, and
  must agree on it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app.phases import phase_from_request
from pytex.core import FrameDomain, Handedness, ReferenceFrame
from pytex.core.lattice import CrystalPlane, MillerIndex
from pytex.core.orientation import OrientationSet
from pytex.core.sphere import S2Grid
from pytex.core.symmetry import SymmetrySpec
from pytex.texture import (
    ODF,
    KernelSpec,
    PoleFigure,
    euler_section_ranges,
    impose_sample_symmetry,
    odf_sections,
)
from pytex.texture.components import (
    HCP_BASAL,
    HCP_BASAL_RD_SPLIT,
    HCP_BASAL_TD_SPLIT,
    HCP_C_ALONG_RD,
    HCP_C_ALONG_TD,
    odf_component_volume_fractions,
    random_component_fraction,
)

SPECIMEN = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
_NI_SPEC, NICKEL = phase_from_request({"builtin": "ni_fcc"})
_ZR_SPEC, ZIRCONIUM = phase_from_request({"builtin": "zr_hcp"})


def _odf(phase, angles: np.ndarray, halfwidth: float = 10.0) -> ODF:
    orientations = OrientationSet.from_euler_angles(
        angles,
        specimen_frame=SPECIMEN,
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    return ODF.from_orientations(orientations, kernel=KernelSpec(halfwidth_deg=halfwidth))


def _random_angles(count: int, seed: int = 3) -> np.ndarray:
    generator = np.random.default_rng(seed)
    return np.column_stack(
        [
            generator.uniform(0.0, 360.0, count),
            np.degrees(np.arccos(generator.uniform(-1.0, 1.0, count))),
            generator.uniform(0.0, 360.0, count),
        ]
    )


def _cube_angles(count: int = 300, spread: float = 6.0, seed: int = 5) -> np.ndarray:
    generator = np.random.default_rng(seed)
    return generator.normal(0.0, spread, size=(count, 3))


def _raster_figure(values_of) -> PoleFigure:
    """A goniometer-like raster: 5 degree tilt rings to 70 degrees, 10 degree azimuth."""

    psi, phi = np.meshgrid(np.arange(0.0, 71.0, 5.0), np.arange(0.0, 360.0, 10.0), indexing="ij")
    theta = np.radians(psi.ravel())
    azimuth = np.radians(phi.ravel())
    directions = np.column_stack(
        [np.sin(theta) * np.cos(azimuth), np.sin(theta) * np.sin(azimuth), np.cos(theta)]
    )
    plane = CrystalPlane(miller=MillerIndex(np.array([1, 1, 1]), phase=NICKEL), phase=NICKEL)
    return PoleFigure(
        pole=plane,
        sample_directions=directions,
        intensities=values_of(directions),
        specimen_frame=SPECIMEN,
        sampling="sampled_density",
    )


class TestAxialSpecimenSymmetry:
    def test_axial_is_a_supported_specimen_symmetry_with_aliases(self) -> None:
        for name in ("axial", "fibre", "fiber", "cylindrical"):
            spec = SymmetrySpec.specimen(name, reference_frame=SPECIMEN)
            assert spec.specimen_symmetry == "axial"
            assert spec.point_group == "∞/mm"

    def test_the_operator_list_is_a_closed_group_of_proper_rotations(self) -> None:
        operators = SymmetrySpec.specimen("axial").operators
        assert operators.shape == (144, 3, 3)
        assert np.allclose(np.linalg.det(operators), 1.0)
        products = np.einsum("aij,bjk->abik", operators[::7], operators[::11]).reshape(-1, 9)
        flat = operators.reshape(-1, 9)
        distance = np.min(np.linalg.norm(products[:, None, :] - flat[None, :, :], axis=2), axis=1)
        assert np.max(distance) < 1e-9

    def test_an_unknown_name_lists_axial_among_the_supported(self) -> None:
        with pytest.raises(ValueError, match="axial"):
            SymmetrySpec.specimen("hexagonal")


class TestImposeSampleSymmetry:
    @staticmethod
    def asymmetric(directions: np.ndarray) -> np.ndarray:
        x, y, z = directions.T
        return 2.0 + 0.6 * x + 0.3 * x * y + 0.5 * z * z

    def test_triclinic_changes_nothing(self) -> None:
        figure = _raster_figure(self.asymmetric)
        result = impose_sample_symmetry(figure, "triclinic")
        assert np.allclose(result.intensities, figure.intensities)
        assert result.sample_symmetry.specimen_symmetry == "triclinic"

    def test_orthorhombic_removes_exactly_the_terms_its_two_folds_forbid(self) -> None:
        figure = _raster_figure(self.asymmetric)
        result = impose_sample_symmetry(figure, "orthorhombic")
        z = figure.sample_directions[:, 2]
        assert np.allclose(result.intensities, 2.0 + 0.5 * z * z, atol=1e-9)

    def test_axial_keeps_only_the_polar_dependence(self) -> None:
        figure = _raster_figure(self.asymmetric)
        result = impose_sample_symmetry(figure, "axial")
        z = figure.sample_directions[:, 2]
        assert np.allclose(result.intensities, 2.0 + 0.5 * z * z, atol=1e-9)
        assert result.sample_symmetry.specimen_symmetry == "axial"

    def test_an_already_axial_figure_is_left_alone(self) -> None:
        figure = _raster_figure(lambda d: 1.0 + np.cos(3.0 * np.arccos(d[:, 2])) ** 2)
        result = impose_sample_symmetry(figure, "axial")
        assert np.allclose(result.intensities, figure.intensities)

    def test_an_unevenly_sampled_ring_is_integrated_not_averaged(self) -> None:
        """Four points at 0, 10, 20 and 180 degrees: a plain mean over-counts the
        clustered side. The arcs between neighbours are 10, 10, 160 and 180
        degrees, so each point owns half of its two arcs: 95, 10, 85 and 170."""

        azimuth = np.radians([0.0, 10.0, 20.0, 180.0])
        theta = math.radians(40.0)
        directions = np.column_stack(
            [
                np.sin(theta) * np.cos(azimuth),
                np.sin(theta) * np.sin(azimuth),
                np.full(4, np.cos(theta)),
            ]
        )
        values = np.array([4.0, 4.0, 4.0, 1.0])
        plane = CrystalPlane(miller=MillerIndex(np.array([1, 1, 1]), phase=NICKEL), phase=NICKEL)
        figure = PoleFigure(
            pole=plane, sample_directions=directions, intensities=values,
            specimen_frame=SPECIMEN, sampling="sampled_density",
        )
        weights = np.array([95.0, 10.0, 85.0, 170.0])
        expected = float(np.dot(weights, values) / weights.sum())
        result = impose_sample_symmetry(figure, "axial")
        assert np.allclose(result.intensities, expected)
        assert expected != pytest.approx(float(values.mean()))

    def test_axial_averaging_preserves_the_area_weighted_mean(self) -> None:
        grid = S2Grid.equispaced(5.0, reference_frame=SPECIMEN, hemisphere="upper")
        directions = np.asarray(grid.vectors, dtype=float)
        weights = np.asarray(grid.weights, dtype=float)
        plane = CrystalPlane(miller=MillerIndex(np.array([1, 1, 1]), phase=NICKEL), phase=NICKEL)
        values = 1.0 + 0.7 * directions[:, 0] ** 2 + 0.2 * directions[:, 1]
        figure = PoleFigure(
            pole=plane, sample_directions=directions, intensities=values,
            specimen_frame=SPECIMEN, sampling="sampled_density",
        )
        result = impose_sample_symmetry(figure, "axial")
        assert np.average(result.intensities, weights=weights) == pytest.approx(
            np.average(values, weights=weights), rel=2e-3
        )


class TestEulerSectionRanges:
    def test_cubic_crystal_orthorhombic_specimen_is_the_ninety_degree_cube(self) -> None:
        assert euler_section_ranges(NICKEL.symmetry, "orthorhombic") == {
            "phi1_max_deg": 90.0, "big_phi_max_deg": 90.0, "phi2_max_deg": 90.0,
        }

    def test_hexagonal_crystal_repeats_every_sixty_degrees_of_phi2(self) -> None:
        ranges = euler_section_ranges(ZIRCONIUM.symmetry, None)
        assert ranges == {"phi1_max_deg": 360.0, "big_phi_max_deg": 90.0, "phi2_max_deg": 60.0}

    def test_the_specimen_symmetry_folds_phi1(self) -> None:
        folded = {
            name: euler_section_ranges(ZIRCONIUM.symmetry, name)["phi1_max_deg"]
            for name in ("triclinic", "monoclinic", "orthorhombic", "axial")
        }
        assert folded == {"triclinic": 360.0, "monoclinic": 180.0, "orthorhombic": 90.0,
                          "axial": 90.0}

    def test_a_triclinic_crystal_and_specimen_need_the_whole_box(self) -> None:
        triclinic = SymmetrySpec.from_point_group("1")
        assert euler_section_ranges(triclinic, None) == {
            "phi1_max_deg": 360.0, "big_phi_max_deg": 180.0, "phi2_max_deg": 360.0,
        }


class TestOdfSections:
    def test_a_random_texture_reads_one_mrd_on_every_kind_of_section(self) -> None:
        odf = _odf(NICKEL, _random_angles(3000), halfwidth=15.0)
        for kind in ("phi2", "phi1", "sigma"):
            sections = odf_sections(
                odf, kind=kind, values_deg=(0.0, 45.0), specimen_symmetry="orthorhombic",
                resolution_deg=15.0,
            )
            assert float(np.mean(sections.densities)) == pytest.approx(1.0, abs=0.08)
            assert sections.section_kind == kind

    def test_the_default_is_the_labotex_plate_over_the_symmetry_range(self) -> None:
        cubic = odf_sections(_odf(NICKEL, _cube_angles(60)), specimen_symmetry="orthorhombic",
                             resolution_deg=15.0)
        assert list(cubic.phi2_deg) == [float(v) for v in range(0, 91, 5)]
        hexagonal = odf_sections(
            _odf(ZIRCONIUM, _cube_angles(60)), specimen_symmetry="orthorhombic",
            resolution_deg=15.0,
        )
        assert list(hexagonal.phi2_deg) == [float(v) for v in range(0, 61, 5)]
        assert float(hexagonal.phi1_deg[-1]) == 90.0

    def test_the_cube_component_peaks_at_the_origin_of_the_phi2_zero_section(self) -> None:
        sections = odf_sections(
            _odf(NICKEL, _cube_angles()), values_deg=(0.0,), specimen_symmetry="orthorhombic"
        )
        plane = sections.densities[0]
        assert plane[0, 0] == pytest.approx(plane.max(), rel=0.1)
        assert plane.max() > 10.0

    def test_phi1_and_phi2_sections_agree_where_they_cross(self) -> None:
        odf = _odf(NICKEL, _cube_angles(spread=12.0))
        by_phi2 = odf_sections(odf, kind="phi2", values_deg=(45.0,), resolution_deg=15.0,
                               specimen_symmetry="orthorhombic")
        by_phi1 = odf_sections(odf, kind="phi1", values_deg=(0.0,), resolution_deg=15.0,
                               specimen_symmetry="orthorhombic")
        assert by_phi1.horizontal_coordinate == "phi2"
        column_phi2_45 = by_phi1.densities[0][:, list(by_phi1.phi1_deg).index(45.0)]
        column_phi1_0 = by_phi2.densities[0][:, 0]
        assert np.allclose(column_phi2_45, column_phi1_0, rtol=1e-9, atol=1e-9)

    def test_sigma_zero_meets_phi2_zero_at_phi1_zero(self) -> None:
        odf = _odf(NICKEL, _cube_angles(spread=12.0))
        sigma = odf_sections(odf, kind="sigma", values_deg=(0.0,), resolution_deg=15.0,
                             specimen_symmetry="orthorhombic")
        phi2 = odf_sections(odf, kind="phi2", values_deg=(0.0,), resolution_deg=15.0,
                            specimen_symmetry="orthorhombic")
        assert np.allclose(sigma.densities[0][:, 0], phi2.densities[0][:, 0])

    def test_an_unknown_kind_is_refused(self) -> None:
        with pytest.raises(ValueError, match="kind"):
            odf_sections(_odf(NICKEL, _cube_angles(20)), kind="omega")


class TestHexagonalComponentsAndVolumeFractions:
    @pytest.mark.parametrize(
        ("component", "expected"),
        [
            (HCP_BASAL, (0.0, 0.0, 1.0)),
            (HCP_BASAL_TD_SPLIT, (0.0, 0.5, math.sqrt(3.0) / 2.0)),
            (HCP_BASAL_RD_SPLIT, (0.5, 0.0, math.sqrt(3.0) / 2.0)),
            (HCP_C_ALONG_TD, (0.0, 1.0, 0.0)),
            (HCP_C_ALONG_RD, (1.0, 0.0, 0.0)),
        ],
    )
    def test_the_basal_pole_lies_where_the_component_name_says(self, component, expected) -> None:
        orientations = OrientationSet.from_euler_angles(
            np.array([component.bunge_euler_deg]),
            specimen_frame=SPECIMEN,
            crystal_frame=ZIRCONIUM.crystal_frame,
            symmetry=ZIRCONIUM.symmetry,
            phase=ZIRCONIUM,
        )
        mapped = orientations.map_crystal_directions(np.array([0.0, 0.0, 1.0]))
        pole = np.asarray(getattr(mapped, "values", mapped), dtype=float)[0]
        # A pole is an axis: +c and -c are the same pole.
        assert min(np.linalg.norm(pole - expected), np.linalg.norm(pole + expected)) < 1e-9

    def test_the_random_fraction_is_the_haar_ball_volume(self) -> None:
        omega = math.radians(15.0)
        assert random_component_fraction(15.0, 24) == pytest.approx(
            24 * (omega - math.sin(omega)) / math.pi
        )
        assert random_component_fraction(15.0, 24) == pytest.approx(
            24 * omega**3 / (6.0 * math.pi), rel=0.02
        )
        assert random_component_fraction(180.0, 24) == 1.0

    def test_a_random_texture_is_one_times_random_in_every_component(self) -> None:
        odf = _odf(ZIRCONIUM, _random_angles(3000), halfwidth=15.0)
        fractions = odf_component_volume_fractions(
            odf, [HCP_BASAL, HCP_C_ALONG_RD], tolerance_deg=20.0, sample_count=600
        )
        for entry in fractions:
            assert entry["times_random"] == pytest.approx(1.0, abs=0.12)
            assert entry["fraction"] == pytest.approx(entry["random_fraction"], rel=0.12)

    def test_a_basal_texture_is_dominated_by_the_basal_component(self) -> None:
        odf = _odf(ZIRCONIUM, _cube_angles(400, spread=5.0), halfwidth=8.0)
        fractions = {
            entry["component"]: entry
            for entry in odf_component_volume_fractions(
                odf, [HCP_BASAL, HCP_C_ALONG_RD], tolerance_deg=15.0, sample_count=800
            )
        }
        assert fractions["basal"]["times_random"] > 5.0
        assert fractions["basal"]["fraction"] > 0.4
        assert fractions["c_along_rd"]["fraction"] < 0.02

    def test_the_sampling_is_reproducible(self) -> None:
        odf = _odf(ZIRCONIUM, _cube_angles(100, spread=8.0))
        first = odf_component_volume_fractions(odf, [HCP_BASAL], sample_count=200, seed=4)
        second = odf_component_volume_fractions(odf, [HCP_BASAL], sample_count=200, seed=4)
        assert first == second
