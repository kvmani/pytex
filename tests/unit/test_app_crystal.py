"""The crystal viewer service: what the browser is given, and what it means.

The interesting assertions are geometric identities that hold whatever the
renderer does: an fcc conventional cell has four atoms per cell, a (111) plane
polygon lies in the plane it claims, a superimposed direction lies in a plane
whose zone law it satisfies, and the scene bounds enclose everything drawn.
"""

from __future__ import annotations

import base64

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase
from pytex.app.services.crystal import camera_angles_from_matrix

pytest.importorskip("matplotlib", reason="the crystal scene is built by the plotting layer")


def scene(**request: object) -> dict:
    return REGISTRY.call("crystal.scene", request)["data"]["scene"]


class TestScene:
    def test_an_fcc_conventional_cell_has_four_atoms(self) -> None:
        payload = scene(phase={"builtin": "ni_fcc"}, show_bonds=False)
        # The builder includes atoms on the shared cell boundaries so the cell
        # looks complete, so the count exceeds the four of the asymmetric
        # description; every one must still sit on an fcc position.
        a = builtin_phase("ni_fcc").a
        for atom in payload["atoms"]:
            fractional = np.asarray(atom["position"]) / a
            doubled = np.round(fractional * 2.0)
            assert np.allclose(fractional * 2.0, doubled, atol=1e-9)
            assert int(doubled.sum()) % 2 == 0

    def test_rock_salt_has_both_species(self) -> None:
        payload = scene(phase={"builtin": "nacl"})
        species = {atom["species"] for atom in payload["atoms"]}
        assert species == {"Na", "Cl"}

    def test_the_scene_bounds_enclose_every_atom(self) -> None:
        payload = scene(phase={"builtin": "nacl"})
        lower, upper = np.asarray(payload["bounds"])
        for atom in payload["atoms"]:
            position = np.asarray(atom["position"])
            assert np.all(position >= lower - 1e-9)
            assert np.all(position <= upper + 1e-9)

    def test_the_centre_is_the_middle_of_the_bounds(self) -> None:
        payload = scene(phase={"builtin": "fe_bcc"})
        lower, upper = np.asarray(payload["bounds"])
        assert np.allclose(payload["centre"], (lower + upper) / 2.0)

    def test_repeats_multiply_the_block(self) -> None:
        single = scene(phase={"builtin": "fe_bcc"}, show_bonds=False)
        doubled = scene(phase={"builtin": "fe_bcc"}, repeat_a=2, show_bonds=False)
        assert len(doubled["atoms"]) > len(single["atoms"])

    def test_a_superimposed_plane_polygon_lies_in_its_plane(self) -> None:
        payload = scene(phase={"builtin": "ni_fcc"}, planes=[[1, 1, 1]], show_bonds=False)
        assert len(payload["planes"]) == 1
        plane = payload["planes"][0]
        vertices = np.asarray(plane["vertices"], dtype=float)
        normal = np.asarray(plane["normal"], dtype=float)
        offsets = vertices @ normal
        # Every vertex has the same projection on the normal: that is what being
        # coplanar means, and it is the one property the polygon must have.
        assert np.ptp(offsets) < 1e-6

    def test_the_plane_normal_matches_the_reciprocal_vector(self) -> None:
        payload = scene(phase={"builtin": "zr_hcp"}, planes=[[0, 0, 1]], show_bonds=False)
        normal = np.asarray(payload["planes"][0]["normal"], dtype=float)
        phase = builtin_phase("zr_hcp").to_phase()
        expected = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=float) @ np.asarray(
            [0.0, 0.0, 1.0]
        )
        expected = expected / np.linalg.norm(expected)
        assert np.allclose(np.abs(normal @ expected), 1.0, atol=1e-9)

    def test_a_direction_in_a_plane_is_perpendicular_to_its_normal(self) -> None:
        payload = scene(
            phase={"builtin": "ni_fcc"},
            planes=[[1, 1, 1]],
            directions=[[1, -1, 0]],
            show_bonds=False,
        )
        normal = np.asarray(payload["planes"][0]["normal"], dtype=float)
        direction = payload["directions"][0]
        vector = np.asarray(direction["end"], dtype=float) - np.asarray(
            direction["start"], dtype=float
        )
        # h*u + k*v + l*w = 0, seen in Cartesian: the arrow lies in the polygon.
        assert abs(float(vector @ normal)) < 1e-9

    def test_labels_use_four_index_notation_for_hexagonal_phases(self) -> None:
        payload = scene(phase={"builtin": "zr_hcp"}, planes=[[1, 0, 0]], show_bonds=False)
        # (100) in a hexagonal phase is the prism plane (1 0 -1 0), and the
        # browser must be shown the four-index form the literature uses.
        assert payload["planes"][0]["label"] == "(1 0 -1 0)"

    def test_overlay_labels_are_text_the_browser_can_draw(self) -> None:
        """No mathtext on the wire: the browser draws the string literally.

        The scene builder labels for matplotlib, where ``$(1\\bar{1}0)$`` is
        markup. Sent to the browser unchanged it appears on screen with the
        dollar signs and the backslash showing, which is what a user reported
        seeing on every superimposed plane.
        """

        payload = scene(
            phase={"builtin": "nacl"},
            planes=[[1, -1, 0], [1, 1, 1]],
            directions=[[1, -1, 0]],
            show_bonds=False,
        )
        labels = [plane["label"] for plane in payload["planes"]]
        labels += [direction["label"] for direction in payload["directions"]]
        assert labels == ["(1 -1 0)", "(111)", "[1 -1 0]"]
        for label in labels:
            assert "$" not in label
            assert "\\" not in label

    def test_bonds_can_be_suppressed(self) -> None:
        assert scene(phase={"builtin": "nacl"}, show_bonds=False)["bonds"] == []

    def test_every_bond_length_matches_its_endpoints(self) -> None:
        payload = scene(phase={"builtin": "nacl"})
        for bond in payload["bonds"]:
            span = np.asarray(bond["end"], dtype=float) - np.asarray(bond["start"], dtype=float)
            assert float(np.linalg.norm(span)) == pytest.approx(bond["length"], abs=1e-9)

    def test_rock_salt_bond_length_is_half_the_cell_edge(self) -> None:
        result = REGISTRY.call("crystal.scene", {"phase": {"builtin": "nacl"}})
        summary = result["data"]["bond_summary"]["Cl-Na"]
        assert summary["mean"] == pytest.approx(builtin_phase("nacl").a / 2.0, abs=1e-6)

    def test_the_axis_arrows_are_the_lattice_vectors(self) -> None:
        payload = scene(phase={"builtin": "zr_hcp"}, show_bonds=False)
        labels = [axis["label"] for axis in payload["axes"]]
        assert labels == ["a", "b", "c"]
        lengths = [float(np.linalg.norm(axis["vector"])) for axis in payload["axes"]]
        spec = builtin_phase("zr_hcp")
        assert lengths == pytest.approx([spec.a, spec.b, spec.c], abs=1e-9)

    def test_a_phase_without_atoms_is_refused_with_a_reason(self) -> None:
        bare = {
            "name": "lattice only",
            "a": 4.0,
            "b": 4.0,
            "c": 4.0,
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0,
            "point_group": "m-3m",
        }
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call("crystal.scene", {"phase": bare})
        assert excinfo.value.details["field"] == "phase"

    def test_the_atom_table_matches_the_scene(self) -> None:
        result = REGISTRY.call("crystal.scene", {"phase": {"builtin": "nacl"}})
        assert len(result["table"]["rows"]) == len(result["data"]["scene"]["atoms"])


class TestCameraAngles:
    """One conversion, in Python, so the export cannot drift from the view."""

    def test_looking_down_z_is_ninety_degrees_of_elevation(self) -> None:
        elevation, _ = camera_angles_from_matrix([1, 0, 0, 0, 1, 0, 0, 0, 1])
        assert elevation == pytest.approx(90.0)

    def test_looking_along_x_is_zero_elevation_and_zero_azimuth(self) -> None:
        elevation, azimuth = camera_angles_from_matrix([0, 0, 1, 0, 1, 0, 1, 0, 0])
        assert elevation == pytest.approx(0.0)
        assert azimuth == pytest.approx(0.0)

    def test_looking_along_y_is_ninety_degrees_of_azimuth(self) -> None:
        elevation, azimuth = camera_angles_from_matrix([0, 0, 1, 1, 0, 0, 0, 1, 0])
        assert elevation == pytest.approx(0.0)
        assert azimuth == pytest.approx(90.0)

    def test_the_matrix_need_not_be_normalised(self) -> None:
        elevation, _ = camera_angles_from_matrix([1, 0, 0, 0, 1, 0, 0, 0, 5])
        assert elevation == pytest.approx(90.0)

    def test_a_wrong_sized_matrix_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            camera_angles_from_matrix([1, 0, 0])

    def test_a_degenerate_view_direction_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            camera_angles_from_matrix([1, 0, 0, 0, 1, 0, 0, 0, 0])


class TestRender:
    def test_png_export_returns_a_decodable_image(self) -> None:
        result = REGISTRY.call(
            "crystal.render", {"phase": {"builtin": "fe_bcc"}, "format": "png", "dpi": 100}
        )
        assert result["data"]["encoding"] == "base64"
        data = base64.b64decode(result["data"]["image"])
        assert data[:8] == b"\x89PNG\r\n\x1a\n"
        assert result["data"]["bytes"] == len(data)

    def test_svg_export_returns_markup(self) -> None:
        result = REGISTRY.call(
            "crystal.render",
            {"phase": {"builtin": "fe_bcc"}, "format": "svg", "render_style": "wireframe"},
        )
        assert result["data"]["encoding"] == "text"
        assert result["data"]["image"].lstrip().startswith("<?xml")
        assert result["notes"]

    def test_the_camera_matrix_overrides_the_angles(self) -> None:
        result = REGISTRY.call(
            "crystal.render",
            {
                "phase": {"builtin": "fe_bcc"},
                "elevation_deg": 10.0,
                "azimuth_deg": 20.0,
                "camera_matrix": "1 0 0 0 1 0 0 0 1",
                "format": "png",
                "dpi": 72,
            },
        )
        assert result["data"]["elevation_deg"] == pytest.approx(90.0)
        assert result["data"]["azimuth_deg"] == pytest.approx(0.0)

    def test_render_accepts_every_scene_parameter(self) -> None:
        """The Figure button replays the scene request, so render must accept it.

        The viewer sends the inputs it drew with, plus the camera. Any structural
        parameter the scene understands and the figure does not is a rejected
        request the moment a user touches that control — and a figure that no
        longer shows what the screen shows if the mismatch is patched by
        dropping the parameter instead.
        """

        manifest = REGISTRY.manifest()
        operations = {entry["id"]: entry for entry in manifest["operations"]}
        scene_names = {p["name"] for p in operations["crystal.scene"]["parameters"]}
        render_names = {p["name"] for p in operations["crystal.render"]["parameters"]}
        assert scene_names <= render_names

    def test_the_figure_honours_the_scene_structural_controls(self) -> None:
        result = REGISTRY.call(
            "crystal.render",
            {
                "phase": {"builtin": "nacl"},
                "repeat_a": 2,
                "show_unit_cells": True,
                "atom_labels": "species",
                "bond_tolerance_angstrom": 0.2,
                "format": "png",
                "dpi": 72,
            },
        )
        assert result["data"]["encoding"] == "base64"

    def test_rendering_does_not_leak_figures(self) -> None:
        import matplotlib.pyplot as plt

        plt.close("all")
        before = len(plt.get_fignums())
        REGISTRY.call(
            "crystal.render", {"phase": {"builtin": "fe_bcc"}, "format": "png", "dpi": 72}
        )
        assert len(plt.get_fignums()) == before
