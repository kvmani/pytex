from __future__ import annotations

from scripts.benchmark_transformation_performance import run


def test_transformation_benchmark_quick_mode_runs_all_pinned_cases() -> None:
    payload = run(quick=True)
    assert payload["schema_id"] == "pytex.benchmarks.transformation_performance"
    results = payload["results"]
    assert isinstance(results, dict)
    assert len(results) == 4
    for entry in results.values():
        assert entry["best_seconds"] > 0.0
