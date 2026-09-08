from __future__ import annotations

from pytex.cli import build_parser


def test_cli_parser_includes_new_docs_validate_and_benchmark_commands() -> None:
    parser = build_parser()

    docs_args = parser.parse_args(["docs", "build"])
    assert callable(docs_args.func)
    assert docs_args.clean is False

    validate_repo_args = parser.parse_args(["validate", "repo"])
    assert callable(validate_repo_args.func)

    validate_manifest_args = parser.parse_args(["validate", "manifests"])
    assert callable(validate_manifest_args.func)

    benchmark_args = parser.parse_args(["benchmarks", "inventory"])
    assert callable(benchmark_args.func)


def test_cli_docs_build_supports_clean_flag() -> None:
    parser = build_parser()
    args = parser.parse_args(["docs", "build", "--clean"])
    assert args.clean is True


def test_cli_hrem_commands() -> None:
    parser = build_parser()

    ctf_args = parser.parse_args(["hrem", "ctf", "--voltage", "200", "--defocus", "-50"])
    assert callable(ctf_args.func)
    assert ctf_args.voltage == 200.0
    assert ctf_args.defocus == -50.0

    sim_args = parser.parse_args([
        "hrem",
        "simulate",
        "--phase",
        "ni_fcc",
        "--sample-type",
        "crystalline",
        "--voltage",
        "200",
        "--mode",
        "double_corrected",
    ])
    assert callable(sim_args.func)
    assert sim_args.phase == "ni_fcc"
    assert sim_args.sample_type == "crystalline"
    assert sim_args.voltage == 200.0
    assert sim_args.mode == "double_corrected"

