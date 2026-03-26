from __future__ import annotations

from scripts.check_repo_integrity import main


def test_repo_integrity_script_passes() -> None:
    assert main() == 0
