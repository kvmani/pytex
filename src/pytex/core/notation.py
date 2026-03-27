from __future__ import annotations

from collections.abc import Sequence


def _normalize_indices(indices: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(int(value) for value in indices)
    if len(normalized) not in {3, 4}:
        raise ValueError("Miller notation helpers require a 3-index or 4-index tuple.")
    return normalized


def _mathtext_component(value: int) -> str:
    if value < 0:
        return rf"\bar{{{abs(value)}}}"
    return str(value)


def _plain_component(value: int) -> str:
    return str(value)


def format_miller_indices(
    indices: Sequence[int],
    *,
    family: str,
    style: str = "mathtext",
) -> str:
    normalized = _normalize_indices(indices)
    if family not in {"direction", "plane"}:
        raise ValueError("family must be either 'direction' or 'plane'.")
    if style not in {"mathtext", "plain"}:
        raise ValueError("style must be either 'mathtext' or 'plain'.")
    opener, closer = ("[", "]") if family == "direction" else ("(", ")")
    component_formatter = _mathtext_component if style == "mathtext" else _plain_component
    payload = "".join(component_formatter(value) for value in normalized)
    if style == "mathtext":
        return f"${opener}{payload}{closer}$"
    return f"{opener}{payload}{closer}"


def format_direction_indices(indices: Sequence[int], *, style: str = "mathtext") -> str:
    return format_miller_indices(indices, family="direction", style=style)


def format_plane_indices(indices: Sequence[int], *, style: str = "mathtext") -> str:
    return format_miller_indices(indices, family="plane", style=style)
