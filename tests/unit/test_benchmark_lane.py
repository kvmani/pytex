from __future__ import annotations

from scripts.benchmark_transformation_performance import run


def test_transformation_benchmark_quick_mode_runs_all_pinned_cases() -> None:
    payload = run(quick=True)
    assert payload["schema_id"] == "pytex.benchmarks.transformation_performance"
    results = payload["results"]
    assert isinstance(results, dict)
    assert set(results) == {
        "intervariant_misorientations_ks_276_pairs",
        "or_deviation_200_pairs",
        "fit_orientation_relationship_200_pairs",
        "weighted_or_fit_200_pairs",
        "reconstruct_parent_grains_30_grains",
    }
    assert results["weighted_or_fit_200_pairs"]["peak_traced_bytes"] > 0
    for entry in results.values():
        assert entry["best_seconds"] > 0.0
