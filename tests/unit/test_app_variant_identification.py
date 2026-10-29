"""The step-by-step variant identification, checked against constructed truth.

Every pair here is built to have a known answer: the child is ``C = P V_kᵀ S_c``
for a chosen variant ``k`` and a chosen child operator ``S_c``. The operation
must name ``k``, write out the operators that produce it, and print numbers in
its equations that are the numbers of the calculation -- so the displayed
matrices are read back out of the MathML and compared with the library's.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.export import result_to_html, result_to_markdown
from pytex.app.mathml import ALLOWED_TAGS, MINUS
from pytex.app.phases import phase_from_request
from pytex.core.orientation import Rotation
from pytex.core.transformation import OrientationRelationship

BETA = {"builtin": "zr_bcc_beta"}
ALPHA = {"builtin": "zr_hcp"}
AUSTENITE = {"builtin": "austenite_fcc"}
FERRITE = {"builtin": "fe_bcc"}
PARENT_EULER = (30.0, 40.0, 10.0)

#: Kurdjumov-Sachs and Nishiyama-Wassermann differ by 5.26 degrees (published,
#: e.g. Nolze, Cryst. Res. Technol. 41 (2006) 72); a KS child tested against NW
#: therefore misses every NW variant by about that much.
KS_TO_NW_DEG = 5.26


def _constructed(
    relationship: OrientationRelationship, position: int, child_operator: int
) -> tuple[float, ...]:
    parent_matrix = Rotation.from_bunge_euler(*PARENT_EULER).as_matrix()
    variant = relationship.generate_variants()[position]
    operator = relationship.child_phase.symmetry.operators[child_operator]
    child = parent_matrix @ variant.parent_to_child_rotation.as_matrix().T @ operator
    return tuple(float(value) for value in Rotation.from_matrix(child).to_bunge_euler(degrees=True))


def _request(parent: dict, child: dict, relationship: str, child_euler, **extra) -> dict:
    request = {
        "phase": parent,
        "child_phase": child,
        "relationship": relationship,
        "euler_convention": "bunge",
        "parent_angle1": PARENT_EULER[0],
        "parent_angle2": PARENT_EULER[1],
        "parent_angle3": PARENT_EULER[2],
        "child_angle1": child_euler[0],
        "child_angle2": child_euler[1],
        "child_angle3": child_euler[2],
        "digits": 6,
    }
    request.update(extra)
    return request


@pytest.fixture(scope="module")
def burgers() -> OrientationRelationship:
    _, parent = phase_from_request(BETA)
    _, child = phase_from_request(ALPHA)
    return OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=child
    )


@pytest.fixture(scope="module")
def ks() -> OrientationRelationship:
    _, parent = phase_from_request(AUSTENITE)
    _, child = phase_from_request(FERRITE)
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )


@pytest.fixture(scope="module")
def burgers_v6(burgers) -> dict:
    euler = _constructed(burgers, 5, 7)
    return REGISTRY.call("variants.identify_variant", _request(BETA, ALPHA, "burgers", euler))


def _stage(result: dict, key: str) -> dict:
    return next(stage for stage in result["stages"] if stage["key"] == key)


def _matrices(node: dict) -> list[np.ndarray]:
    """Every numeric matrix typeset in a tree, in reading order."""

    found: list[np.ndarray] = []

    def walk(item: dict) -> None:
        if item.get("tag") == "mtable":
            rows = []
            for tr in item["children"]:
                cells = []
                for td in tr["children"]:
                    leaf = td["children"][0]
                    if leaf.get("tag") == "mrow":
                        # A phantom minus pads a positive entry in a signed column.
                        leaf = leaf["children"][-1]
                    if leaf.get("tag") != "mn":
                        return
                    cells.append(float(leaf["text"].replace(MINUS, "-")))
                rows.append(cells)
            found.append(np.array(rows))
            return
        for child in item.get("children", ()):
            walk(child)

    walk(node)
    return found


def _tags(node: dict) -> set[str]:
    tags = {node["tag"]}
    for child in node.get("children", ()):
        tags |= _tags(child)
    return tags


class TestIdentification:
    def test_a_constructed_burgers_child_is_its_variant(self, burgers_v6) -> None:
        assert burgers_v6["data"]["variant_index"] == 6
        assert burgers_v6["data"]["misfit_deg"] < 1e-3
        assert burgers_v6["data"]["verdict"] == "assigned"

    @pytest.mark.parametrize("position", [0, 7, 19])
    def test_a_constructed_ks_child_is_its_variant(self, ks, position) -> None:
        euler = _constructed(ks, position, 11)
        result = REGISTRY.call(
            "variants.identify_variant", _request(AUSTENITE, FERRITE, "kurdjumov_sachs", euler)
        )
        assert result["data"]["variant_index"] == position + 1

    def test_a_ks_child_tested_against_nw_does_not_fit(self, ks) -> None:
        euler = _constructed(ks, 0, 0)
        result = REGISTRY.call(
            "variants.identify_variant",
            _request(AUSTENITE, FERRITE, "nishiyama_wassermann", euler),
        )
        assert result["data"]["verdict"] == "no_fit"
        assert result["data"]["misfit_deg"] == pytest.approx(KS_TO_NW_DEG, abs=0.01)
        assert result["warnings"]

    def test_the_matthies_convention_gives_the_same_answer(self, burgers) -> None:
        bunge = _constructed(burgers, 3, 0)

        def to_matthies(angles):
            matrix = Rotation.from_bunge_euler(*angles).as_matrix()
            return tuple(Rotation.from_matrix(matrix).to_euler(convention="matthies", degrees=True))

        parent = to_matthies(PARENT_EULER)
        child = to_matthies(bunge)
        request = _request(BETA, ALPHA, "burgers", child, euler_convention="matthies")
        request.update(parent_angle1=parent[0], parent_angle2=parent[1], parent_angle3=parent[2])
        result = REGISTRY.call("variants.identify_variant", request)
        assert result["data"]["variant_index"] == 4
        # The ZYZ product is written with R_y in the middle, not R_x.
        first = _stage(result, "parent_orientation")["equations"][0]["math"]
        scripts = [
            child["children"][1]["text"]
            for child in first["children"][0]["children"]
            if child.get("tag") == "mrow"
            for child in child["children"]
            if child.get("tag") == "msub"
        ]
        assert scripts == ["z", "y", "z"]

    def test_a_stated_relationship_identifies_like_the_catalogued_one(self, ks) -> None:
        euler = _constructed(ks, 9, 2)
        result = REGISTRY.call(
            "variants.identify_variant",
            _request(
                AUSTENITE,
                FERRITE,
                "custom",
                euler,
                custom_parent_plane=[1, 1, 1],
                custom_child_plane=[0, 1, 1],
                custom_parent_direction=[-1, 0, 1],
                custom_child_direction=[-1, -1, 1],
            ),
        )
        assert result["data"]["misfit_deg"] < 1e-3


class TestDerivation:
    def test_nine_steps_in_the_derivation_section(self, burgers_v6) -> None:
        stages = burgers_v6["stages"]
        assert [stage["title"].split(".")[0] for stage in stages] == [
            f"Step {number}" for number in range(1, 10)
        ]
        assert {stage["section"] for stage in stages} == {"derivation"}

    def test_equations_are_numbered_once_and_in_order(self, burgers_v6) -> None:
        tags = [
            int(equation["tag"])
            for stage in burgers_v6["stages"]
            for equation in stage.get("equations", ())
        ]
        assert tags == list(range(1, len(tags) + 1))
        assert len(tags) >= 20

    def test_every_tree_uses_only_whitelisted_tags(self, burgers_v6) -> None:
        for stage in burgers_v6["stages"]:
            for equation in stage.get("equations", ()):
                assert _tags(equation["math"]) <= ALLOWED_TAGS

    def test_the_typeset_orientation_matrices_are_the_measured_ones(self, burgers_v6) -> None:
        identification = burgers_v6["data"]["identification"]
        parent_eq = _stage(burgers_v6, "parent_orientation")["equations"][-1]["math"]
        child_eq = _stage(burgers_v6, "child_orientation")["equations"][-1]["math"]
        assert np.allclose(_matrices(parent_eq)[-1], identification["parent_matrix"], atol=1e-6)
        assert np.allclose(_matrices(child_eq)[-1], identification["child_matrix"], atol=1e-6)

    def test_the_relationship_matrix_is_the_product_of_the_triads(
        self, burgers_v6, burgers
    ) -> None:
        stage = _stage(burgers_v6, "relationship")
        rotation = _matrices(stage["equations"][-1]["math"])[-1]
        assert np.allclose(rotation, burgers.parent_to_child_rotation.as_matrix(), atol=1e-6)
        triads = [
            _matrices(equation["math"])[0]
            for equation in stage["equations"]
            if len(_matrices(equation["math"])) == 1
            and _matrices(equation["math"])[0].shape == (3, 3)
        ]
        f_parent, f_child = triads[0], triads[1]
        assert np.allclose(f_child @ f_parent.T, rotation, atol=1e-5)

    def test_the_winning_operators_are_written_out(self, burgers_v6, burgers) -> None:
        stage = _stage(burgers_v6, "winner")
        matrices = [_matrices(equation["math"])[-1] for equation in stage["equations"][:5]]
        parent_operator, variant, predicted, child_operator, residual = matrices
        variant_rotation = burgers.generate_variants()[5]
        assert np.allclose(parent_operator, variant_rotation.parent_symmetry_operator, atol=1e-6)
        assert np.allclose(
            variant, variant_rotation.parent_to_child_rotation.as_matrix(), atol=1e-6
        )
        # The child was written with operator 7; the match is its inverse.
        expected = burgers.child_phase.symmetry.operators[7].T
        assert np.allclose(child_operator, expected, atol=1e-6)
        assert np.allclose(residual, np.eye(3), atol=1e-5)
        identification = burgers_v6["data"]["identification"]
        assert np.allclose(
            predicted,
            np.asarray(identification["parent_matrix"]) @ variant.T,
            atol=1e-5,
        )

    def test_the_misfit_table_ranks_the_identified_variant_first(self, burgers_v6) -> None:
        rows = _stage(burgers_v6, "misfits")["table"]["rows"]
        first = next(row for row in rows if row["rank"] == 1)
        assert first["variant"] == 6
        assert len(rows) == 12

    def test_the_parallelism_check_is_exact_for_an_exact_child(self, burgers_v6) -> None:
        rows = _stage(burgers_v6, "parallelism")["table"]["rows"]
        assert {row["kind"] for row in rows} == {"plane", "direction"}
        assert all(row["angle_deg"] < 1e-3 for row in rows)

    def test_a_negative_index_is_overbarred_in_the_equations(self, burgers_v6) -> None:
        stage = _stage(burgers_v6, "relationship")
        assert any("mover" in _tags(equation["math"]) for equation in stage["equations"])


class TestReports:
    def test_the_markdown_report_carries_the_numbered_derivation(self, burgers_v6) -> None:
        text = result_to_markdown(burgers_v6).decode()
        assert "## Step-by-step derivation" in text
        assert r"\tag{1}" in text
        assert r"\begin{matrix}" in text

    def test_the_html_report_carries_mathml(self, burgers_v6) -> None:
        page = result_to_html(burgers_v6).decode()
        assert page.count("<math ") >= 20
        assert "Step 7. Variant V6, substituted in full" in page
