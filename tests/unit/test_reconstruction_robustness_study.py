from __future__ import annotations

from scripts.study_reconstruction_robustness import run


def test_robustness_study_quick_mode_runs_and_reports_clean_partitions() -> None:
    """The study lane must stay runnable, and its noise-free cell must be perfect.

    The quick sweep is a smoke size, so the assertions pin only what must hold
    by construction: every cell reports the bookkeeping fields, and with zero
    orientation noise the planted partition is recovered exactly with zero
    parent error. The full sweep's quantitative envelope lives in
    `docs/testing/reconstruction_robustness_study.md`.
    """

    payload = run(quick=True)
    assert payload["study"] == "parent_grain_reconstruction_robustness"
    assert payload["relationship"] == "kurdjumov_sachs"
    cells = payload["cells"]
    assert isinstance(cells, list)
    assert cells

    for cell in cells:
        assert cell["trials_judged"] > 0, cell
        # The fingerprint edge test must never merge separable parents. This is
        # the property the study exists to defend; an angle-only edge test
        # drives it well above zero.
        assert cell["false_link_rate"] == 0.0, cell

    noise_free = [cell for cell in cells if cell["noise_deg"] == 0.0]
    assert noise_free, "the sweep must include a noise-free cell"
    for cell in noise_free:
        assert cell["partition_exact"] == 1.0, cell
        assert cell["missed_link_rate"] == 0.0, cell
        assert cell["parent_error_max_deg"] < 1e-6, cell

    # Map-scale sweep: tiled parents with a dense four-connected grain graph.
    map_cells = payload["map_scale_cells"]
    assert map_cells, "the sweep must include map-scale cells"
    assert {str(cell["relationship"]) for cell in map_cells} == {
        "kurdjumov_sachs",
        "burgers",
    }, "map scale must cover a cubic-cubic and a cubic-hexagonal relationship"
    for cell in map_cells:
        # Same defended property, now on a graph dense enough that every parent
        # pair shares several boundaries. Merges here must come only from
        # genuinely ambiguous boundaries, which are reported separately.
        assert cell["false_link_rate"] == 0.0, cell
        assert cell["clusters_mean"] <= cell["true_parent_count"], cell
        assert 0.0 <= cell["ambiguous_cross_fraction"] < 1.0, cell
