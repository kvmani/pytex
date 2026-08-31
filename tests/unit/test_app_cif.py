"""CIF files enter every Workbench calculation through the shared phase contract."""

from __future__ import annotations

from dataclasses import replace

import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.phases import PhaseSpec, builtin_phase, phase_from_request
from pytex.core.fixtures import get_phase_fixture, phase_fixtures_available
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord


def _cif_payload(*, name: str = "nickel.cif", text: str = "data_nickel") -> dict:
    return {"cif": {"name": name, "text": text}}


def _fake_cif_phase(
    cif_text: str,
    *,
    crystal_frame: ReferenceFrame,
    provenance: ProvenanceRecord,
    **_: object,
) -> Phase:
    assert cif_text == "data_nickel"
    phase = builtin_phase("ni_fcc").to_phase()
    assert phase.crystal_frame == crystal_frame
    return replace(phase, provenance=provenance)


def test_cif_upload_becomes_the_same_phase_spec_services_already_use(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(Phase, "from_cif_string", staticmethod(_fake_cif_phase))

    spec, phase = phase_from_request(_cif_payload())

    assert spec.name == "Nickel (fcc)"
    assert spec.point_group == "m-3m"
    assert spec.space_group_symbol == "Fm-3m"
    assert spec.space_group_number == 225
    assert len(spec.sites) == 4
    assert spec.source is not None and "nickel.cif" in spec.source
    assert phase.provenance is not None
    assert phase.provenance.source_system == "cif"
    assert phase.provenance.source_identifier == "nickel.cif"
    assert phase.provenance.metadata["reader"] == "pymatgen.Structure.from_str"


def test_phase_spec_from_phase_preserves_metric_symmetry_and_atomic_basis() -> None:
    original = builtin_phase("zr_hcp")
    restored = PhaseSpec.from_phase(original.to_phase(), source="test source")

    assert restored.name == original.name
    assert restored.a == original.a
    assert restored.c == original.c
    assert restored.point_group == original.point_group
    assert restored.space_group_number == original.space_group_number
    assert [site.species for site in restored.sites] == [site.species for site in original.sites]
    assert [(site.x, site.y, site.z) for site in restored.sites] == [
        (site.x, site.y, site.z) for site in original.sites
    ]
    assert restored.source == "test source"


def test_cif_upload_rejects_another_extension_before_calling_the_parser(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser_called = False

    def parser(*_: object, **__: object) -> Phase:
        nonlocal parser_called
        parser_called = True
        return builtin_phase("ni_fcc").to_phase()

    monkeypatch.setattr(Phase, "from_cif_string", staticmethod(parser))
    with pytest.raises(InvalidInputError, match=r"extension is \.txt") as excinfo:
        phase_from_request(_cif_payload(name="nickel.txt"))

    assert not parser_called
    assert excinfo.value.details["field"] == "phase"


def test_cif_upload_cannot_silently_mix_file_and_manual_fields() -> None:
    with pytest.raises(InvalidInputError, match="cannot also override") as excinfo:
        phase_from_request(_cif_payload() | {"a": 4.0})
    assert excinfo.value.details["field"] == "phase"


def test_a_real_cif_upload_is_parsed_into_a_phase() -> None:
    """The CIF path end to end, with the parser actually present.

    pymatgen is a required dependency as of 0.5.0, so there is no longer a lane
    in which this route answers with an install hint instead of a phase. This
    asserts the outcome the workbench user gets rather than the absence of the
    error that used to stand in for it, and it reads the checksum-pinned
    nickel fixture rather than a CIF typed into the test.
    """

    if not phase_fixtures_available():
        pytest.skip("the checksum-pinned phase-fixture corpus is a source-checkout asset")


    fixture = get_phase_fixture("ni_fcc")
    spec, phase = phase_from_request(_cif_payload(text=fixture.read_cif_text()))

    assert phase.lattice.a == pytest.approx(3.52387)
    assert len(phase.unit_cell.sites) == 4
    assert {site.species for site in phase.unit_cell.sites} == {"Ni"}
    assert phase.space_group_number == 225
    assert "pymatgen" in spec.source


def test_malformed_cif_is_reported_against_the_shared_phase_control(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def malformed(*_: object, **__: object) -> Phase:
        raise ValueError("missing cell length")

    monkeypatch.setattr(Phase, "from_cif_string", staticmethod(malformed))
    with pytest.raises(InvalidInputError, match="missing cell length") as excinfo:
        phase_from_request(_cif_payload())

    assert excinfo.value.details["field"] == "phase"
    assert "valid text CIF" in (excinfo.value.hint or "")
