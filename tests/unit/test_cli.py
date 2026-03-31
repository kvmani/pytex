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
