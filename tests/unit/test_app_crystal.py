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
from pytex.app.services.crystal import (
    _appearance,
    _appearance_style,
    camera_angles_from_matrix,
    camera_matrix_from_euler,
    euler_from_camera_matrix,
    orientation_overlay,
)


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
    def test_object_properties_are_validated_and_mapped_to_the_renderer(self) -> None:
        appearance = _appearance(
            {
                "atom_scale": 1.4,
                "surface_finish": "matte",
                "light_direction": [0.0, 3.0, 4.0],
                "light_ambient": 0.3,
                "light_diffuse": 0.9,
                "light_specular": 0.6,
                "atom_shininess": 48,
                "depth_cue_strength": 0.25,
                "plane_color": "#123456",
                "plane_opacity": 0.42,
                "species_colors": {"Fe": "#abcdef"},
                "show_cells": False,
            }
        )
        style = _appearance_style(appearance, "ball_and_stick")["crystal"]
        assert style["atom_radius_scale"] == pytest.approx(0.55 * 1.4)
        assert style["plane_color"] == "#123456"
        assert style["plane_alpha"] == pytest.approx(0.42)
        assert style["species_colors"] == {"Fe": "#abcdef"}
        assert style["lattice_linewidth"] == 0.0
        assert style["light_direction"] == pytest.approx([0.0, 0.6, 0.8])
        assert style["light_ambient"] == pytest.approx(0.3)
        assert style["light_diffuse"] == pytest.approx(0.9)
        assert style["light_specular"] == pytest.approx(0.6)
        assert style["atom_specular_strength"] == pytest.approx(0.12)
        assert style["atom_shininess"] == pytest.approx(48)
        assert style["depth_cue_strength"] == pytest.approx(0.25)

    def test_invalid_object_properties_are_rejected_at_the_boundary(self) -> None:
        with pytest.raises(InvalidInputError, match="plane_opacity"):
            _appearance({"plane_opacity": 4.2})
        with pytest.raises(InvalidInputError, match="RRGGBB"):
            _appearance({"direction_color": "blue"})
        with pytest.raises(InvalidInputError, match="surface_finish"):
            _appearance({"surface_finish": "chrome"})
        with pytest.raises(InvalidInputError, match="light_direction"):
            _appearance({"light_direction": [0, 0, 0]})

    def test_per_species_colours_reach_the_scene_glyphs(self) -> None:
        from pytex.plotting.crystal3d import build_crystal_scene

        phase = builtin_phase("fe_bcc").to_phase()
        rendered = build_crystal_scene(
            phase,
            show_bonds=False,
            style_overrides={"crystal": {"species_colors": {"Fe": "#123456"}}},
        )
        assert {atom.color for atom in rendered.atoms} == {"#123456"}

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


class TestOrientationOverlay:
    """What the browser is handed so it can draw a live pole figure.

    The assertions are the invariants a pole figure depends on and a renderer
    cannot repair: the poles are unit vectors, each family is closed under the
    point group, the standard triangle's corners are exact low-index directions,
    and its outline lies on the sector rather than inside it.
    """

    def test_every_pole_is_a_unit_vector(self) -> None:
        overlay = orientation_overlay(builtin_phase("cu_fcc"))
        for family in overlay["pole_families"]:
            vectors = np.asarray(family["vectors"], dtype=float)
            assert np.allclose(np.linalg.norm(vectors, axis=1), 1.0)

    def test_the_cubic_families_have_their_textbook_multiplicities(self) -> None:
        overlay = orientation_overlay(builtin_phase("cu_fcc"))
        counts = {family["label"]: len(family["vectors"]) for family in overlay["pole_families"]}
        # Signed poles, so twice the usual family size: {100} is six faces,
        # {110} twelve, {111} eight.
        assert counts == {"{100}": 6, "{110}": 12, "{111}": 8}

    def test_each_family_is_closed_under_the_point_group(self) -> None:
        overlay = orientation_overlay(builtin_phase("ti_hcp"))
        operators = np.asarray(overlay["operators"], dtype=float).reshape(-1, 3, 3)
        for family in overlay["pole_families"]:
            vectors = np.asarray(family["vectors"], dtype=float)
            images = np.einsum("nij,mj->nmi", operators, vectors).reshape(-1, 3)
            for image in images:
                distances = np.linalg.norm(vectors - image, axis=1)
                assert float(distances.min()) < 1e-9

    def test_the_standard_triangle_corners_are_exact_low_index_directions(self) -> None:
        for key, expected in (
            ("cu_fcc", {"[001]", "[101]", "[111]"}),
            ("ti_hcp", {"[0001]", "[2 -1 -1 0]", "[1 0 -1 0]"}),
        ):
            corners = orientation_overlay(builtin_phase(key))["sector"]["corners"]
            assert {corner["label"] for corner in corners} == expected
            assert all(corner["residual_deg"] < 1e-6 for corner in corners)

    def test_the_sector_outline_lies_on_the_sector_boundary(self) -> None:
        overlay = orientation_overlay(builtin_phase("cu_fcc"))
        normals = np.asarray(overlay["sector"]["edge_normals"], dtype=float)
        outline = np.asarray(overlay["sector"]["outline"], dtype=float)
        assert outline.shape[0] > 3
        # Invert the stereographic projection and check that every sampled point
        # is inside the sector and on one of its bounding planes. A straight line
        # between projected corners would fail the second half.
        radius_squared = np.sum(outline**2, axis=1)
        z = (1.0 - radius_squared) / (1.0 + radius_squared)
        scale = 1.0 + z
        directions = np.column_stack([outline[:, 0] * scale, outline[:, 1] * scale, z])
        products = directions @ normals.T
        assert np.all(products > -1e-9)
        assert np.all(np.min(np.abs(products), axis=1) < 1e-9)

    def test_the_overlays_of_the_scene_travel_with_it(self) -> None:
        overlay = orientation_overlay(
            builtin_phase("cu_fcc"), plane_rows=((1, 1, 1),), direction_rows=((1, -1, 0),)
        )
        assert [entry["label"] for entry in overlay["overlay_poles"]] == ["(111)"]
        assert [entry["label"] for entry in overlay["overlay_directions"]] == ["[1 -1 0]"]
        pole = np.asarray(overlay["overlay_poles"][0]["vector"], dtype=float)
        direction = np.asarray(overlay["overlay_directions"][0]["vector"], dtype=float)
        # The Burgers direction lies in its slip plane, which is why the viewer
        # draws them together; the pole figure must say the same thing.
        assert float(np.dot(pole, direction)) == pytest.approx(0.0, abs=1e-12)

    def test_the_specimen_frame_is_the_screen(self) -> None:
        overlay = orientation_overlay(builtin_phase("fe_bcc"))
        axes = {entry["label"]: entry["vector"] for entry in overlay["specimen_axes"]}
        assert axes == {"RD": [1.0, 0.0, 0.0], "TD": [0.0, 1.0, 0.0], "ND": [0.0, 0.0, 1.0]}

    def test_the_scene_carries_the_overlay(self) -> None:
        payload = scene(phase={"builtin": "cu_fcc"}, planes=[[1, 1, 1]])
        assert payload["orientation"]["point_group"] == "m-3m"
        assert [entry["label"] for entry in payload["orientation"]["overlay_poles"]] == ["(111)"]


class TestEulerAngles:
    """The camera and an orientation are the same object, so this is a bijection."""

    def test_the_identity_camera_is_the_identity_orientation(self) -> None:
        assert euler_from_camera_matrix([1, 0, 0, 0, 1, 0, 0, 0, 1]) == pytest.approx(
            (0.0, 0.0, 0.0)
        )
        assert camera_matrix_from_euler(0.0, 0.0, 0.0) == pytest.approx(
            [1, 0, 0, 0, 1, 0, 0, 0, 1]
        )

    def test_the_camera_matrix_is_the_crystal_to_specimen_orientation(self) -> None:
        from pytex.core.orientation import Rotation

        camera = np.asarray(camera_matrix_from_euler(30.0, 45.0, 60.0), dtype=float).reshape(3, 3)
        expected = Rotation.from_bunge_euler(30.0, 45.0, 60.0).as_matrix()
        assert np.allclose(camera, expected)

    @pytest.mark.parametrize(
        "angles",
        [(0.0, 45.0, 0.0), (30.0, 45.0, 60.0), (359.0, 179.0, 1.0), (12.5, 0.0, 0.0)],
    )
    def test_angles_survive_the_round_trip(self, angles: tuple[float, float, float]) -> None:
        camera = camera_matrix_from_euler(*angles)
        recovered = euler_from_camera_matrix(camera)
        assert np.allclose(camera_matrix_from_euler(*recovered), camera, atol=1e-12)

    def test_both_conventions_name_the_same_rotation(self) -> None:
        camera = camera_matrix_from_euler(30.0, 45.0, 60.0, convention="bunge")
        matthies = euler_from_camera_matrix(camera, convention="matthies")
        assert np.allclose(
            camera_matrix_from_euler(*matthies, convention="matthies"), camera, atol=1e-12
        )

    def test_a_drifted_camera_is_repaired_rather_than_refused(self) -> None:
        camera = np.asarray(camera_matrix_from_euler(30.0, 45.0, 60.0), dtype=float)
        drifted = camera + 1e-9
        assert euler_from_camera_matrix(drifted) == pytest.approx((30.0, 45.0, 60.0), abs=1e-5)

    def test_a_matrix_that_is_not_a_rotation_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            euler_from_camera_matrix([1, 0, 0, 0, 1, 0, 0, 0, 2])
        with pytest.raises(InvalidInputError):
            euler_from_camera_matrix([1, 0, 0, 0, 1, 0, 0, 0, -1])
        with pytest.raises(InvalidInputError):
            euler_from_camera_matrix([1, 0, 0])

    def test_an_unknown_convention_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            camera_matrix_from_euler(0.0, 0.0, 0.0, convention="kocks")
        assert excinfo.value.details["field"] == "euler_convention"

    def test_a_non_finite_angle_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            camera_matrix_from_euler(float("nan"), 0.0, 0.0)


class TestOrientationOperation:
    """The numbers the dock reports, and the conventions they hold under."""

    def orientation(self, **request: object) -> dict:
        return REGISTRY.call("crystal.orientation", {"phase": {"builtin": "cu_fcc"}, **request})

    def test_the_identity_view_looks_down_the_c_axis(self) -> None:
        result = self.orientation()
        rows = {row["axis"]: row for row in result["table"]["rows"]}
        # With the identity camera the screen axes are the crystal axes, so ND
        # -- out of the screen -- is [001], and RD and TD are equivalent to it
        # under cubic symmetry.
        assert rows["ND"]["direction"] == "[001]"
        assert rows["ND"]["polar_deg"] == pytest.approx(0.0, abs=1e-9)

    def test_a_camera_matrix_overrides_the_angles(self) -> None:
        result = self.orientation(
            angle1=10.0, angle2=20.0, angle3=30.0, camera_matrix="1 0 0 0 1 0 0 0 1"
        )
        assert result["data"]["euler"]["angles_deg"] == pytest.approx([0.0, 0.0, 0.0])

    def test_the_reported_matrix_is_the_orientation_the_angles_name(self) -> None:
        result = self.orientation(angle1=30.0, angle2=45.0, angle3=60.0)
        assert result["data"]["camera_matrix"] == pytest.approx(
            camera_matrix_from_euler(30.0, 45.0, 60.0)
        )

    def test_the_poles_are_in_the_upper_hemisphere_and_on_the_disc(self) -> None:
        result = self.orientation(angle1=17.0, angle2=41.0, angle3=63.0)
        for family in result["data"]["poles"]:
            for point in family["points"]:
                assert point["specimen"][2] >= -1e-9
                assert float(np.hypot(point["x"], point["y"])) <= 1.0 + 1e-9

    def test_the_specimen_axes_land_inside_the_fundamental_sector(self) -> None:
        overlay = orientation_overlay(builtin_phase("cu_fcc"))
        normals = np.asarray(overlay["sector"]["edge_normals"], dtype=float)
        result = self.orientation(angle1=17.0, angle2=41.0, angle3=63.0)
        for point in result["data"]["ipf_points"]:
            crystal = np.asarray(point["crystal"], dtype=float)
            assert np.all(crystal @ normals.T > -1e-9)

    def test_the_specimen_axis_is_recovered_from_its_crystal_direction(self) -> None:
        result = self.orientation(angle1=17.0, angle2=41.0, angle3=63.0)
        camera = np.asarray(result["data"]["camera_matrix"], dtype=float).reshape(3, 3)
        overlay = orientation_overlay(builtin_phase("cu_fcc"))
        operators = np.asarray(overlay["operators"], dtype=float).reshape(-1, 3, 3)
        axes = {"rd": [1.0, 0.0, 0.0], "td": [0.0, 1.0, 0.0], "nd": [0.0, 0.0, 1.0]}
        for point in result["data"]["ipf_points"]:
            crystal = np.asarray(point["crystal"], dtype=float)
            # Some symmetry image of the reported direction maps back onto the
            # specimen axis it was measured along; that is the whole content of
            # "folded into the fundamental sector".
            images = camera @ np.einsum("nij,j->ni", operators, crystal).T
            target = np.asarray(axes[point["key"]], dtype=float)
            best = float(np.max(np.abs(images.T @ target)))
            assert best == pytest.approx(1.0, abs=1e-9)

    def test_a_hexagonal_phase_is_labelled_with_four_indices(self) -> None:
        result = REGISTRY.call(
            "crystal.orientation", {"phase": {"builtin": "ti_hcp"}, "angle2": 30.0}
        )
        assert all(row["direction"].count(" ") == 3 for row in result["table"]["rows"])

    def test_a_malformed_camera_matrix_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            self.orientation(camera_matrix="1 0 0 0 1")
        assert excinfo.value.details["field"] == "camera_matrix"
        with pytest.raises(InvalidInputError):
            self.orientation(camera_matrix="one two three four five six seven eight nine")

    def test_the_result_explains_its_conventions(self) -> None:
        result = self.orientation(angle2=45.0)
        assert "v_specimen = C v_crystal" in result["summary"]
        assert any("RD is screen right" in note for note in result["notes"])
