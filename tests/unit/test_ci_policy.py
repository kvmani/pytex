"""CI policy: Linux and Windows only, never macOS; and started by hand only.

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


def _triggers(workflow: Path) -> set[str]:
    """The event names under a workflow's top-level ``on:`` key."""

    import yaml

    document = yaml.safe_load(workflow.read_text(encoding="utf-8"))
    # PyYAML reads the bare key ``on`` as the boolean True (YAML 1.1).
    triggers = document.get("on", document.get(True))
    if isinstance(triggers, str):
        return {triggers}
    if isinstance(triggers, list):
        return set(triggers)
    return set(triggers or {})


def test_ci_is_started_by_hand_only() -> None:
    """CI runs only when the maintainer starts it (decision of 2026-10-29).

    No push, pull request or schedule may start the CI workflow. The release
    workflow is the one exception, on a pushed version tag, because pushing a
    tag is itself the deliberate act of releasing; it may also be started by
    hand with the version to release, and then creates the tag itself.
    """

    assert _triggers(WORKFLOWS / "ci.yml") == {"workflow_dispatch"}


def test_the_release_workflow_runs_only_on_a_pushed_version_tag() -> None:
    workflow = WORKFLOWS / "release.yml"
    assert _triggers(workflow) <= {"push", "workflow_dispatch"}
    import yaml

    document = yaml.safe_load(workflow.read_text(encoding="utf-8"))
    push = (document.get("on", document.get(True)) or {}).get("push") or {}
    assert set(push) == {"tags"}, "release.yml may react to pushed tags only, never to branches"


@pytest.mark.parametrize("workflow", _workflow_files(), ids=lambda path: path.name)
def test_no_workflow_is_scheduled_or_runs_on_pull_requests(workflow: Path) -> None:
    triggers = _triggers(workflow)
    assert not triggers & {"schedule", "pull_request", "pull_request_target"}, (
        f"{workflow.name} would start without the maintainer: {sorted(triggers)}"
    )


def test_the_policy_is_documented_where_contributors_look() -> None:
    root = WORKFLOWS.parents[1]
    agents = (root / "AGENTS.md").read_text(encoding="utf-8")
    strategy = (root / "docs" / "testing" / "strategy.md").read_text(encoding="utf-8")
    assert "## Continuous Integration Platforms" in agents
    assert "started by hand only" in agents
    assert "manual only" in strategy
    assert "macOS is\nnot a CI platform" in strategy or "macOS is not a CI platform" in strategy
