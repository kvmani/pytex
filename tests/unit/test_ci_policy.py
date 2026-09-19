"""CI platform policy: Linux and Windows only, never macOS.

The maintainer removed macOS from CI on 2026-09-19 and asked that no future
change reintroduce it (AGENTS.md, "Continuous Integration Platforms"). A runner
label is easy to add back by copying a matrix from elsewhere, so the rule is
checked here against every workflow file rather than left to review.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

WORKFLOWS = Path(__file__).resolve().parents[2] / ".github" / "workflows"


def _workflow_files() -> list[Path]:
    return sorted([*WORKFLOWS.glob("*.yml"), *WORKFLOWS.glob("*.yaml")])


def test_there_are_workflows_to_check() -> None:
    assert _workflow_files(), "no workflow files found; the policy check would be vacuous"


@pytest.mark.parametrize("workflow", _workflow_files(), ids=lambda path: path.name)
def test_no_workflow_schedules_a_macos_runner(workflow: Path) -> None:
    offending = [
        f"{workflow.name}:{number}: {line.strip()}"
        for number, line in enumerate(workflow.read_text(encoding="utf-8").splitlines(), 1)
        if not line.lstrip().startswith("#") and re.search(r"mac[\s_-]*os", line, re.IGNORECASE)
    ]
    assert not offending, (
        "macOS is not a CI platform (AGENTS.md, 'Continuous Integration Platforms'): "
        + "; ".join(offending)
    )


def test_the_policy_is_documented_where_contributors_look() -> None:
    root = WORKFLOWS.parents[1]
    agents = (root / "AGENTS.md").read_text(encoding="utf-8")
    strategy = (root / "docs" / "testing" / "strategy.md").read_text(encoding="utf-8")
    assert "## Continuous Integration Platforms" in agents
    assert "macOS is\nnot a CI platform" in strategy or "macOS is not a CI platform" in strategy
