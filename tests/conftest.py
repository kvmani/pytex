"""Repository-wide test fixtures enforcing the hygiene rules of the development guide."""

from __future__ import annotations

from collections.abc import Iterator

import pytest


@pytest.fixture(autouse=True)
def _close_matplotlib_figures() -> Iterator[None]:
    """Close every matplotlib figure a test leaves open.

    Plotting tests create real figures; without this, figures accumulate
    across the session and matplotlib emits "More than 20 figures" resource
    warnings (development-guide finding 18). Closing after each test keeps
    memory flat and lets the warnings-as-errors policy stay strict.
    """

    yield
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plt.close("all")
