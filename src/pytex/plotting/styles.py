from __future__ import annotations

import importlib.resources as resources
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

import yaml

_THEME_PACKAGE = "pytex.plotting.themes"
_DEFAULT_THEME = "journal"
_TOP_LEVEL_SECTIONS = {"common", "xrd", "saed", "crystal"}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _validate_style_mapping(style: dict[str, Any]) -> dict[str, Any]:
    unknown = set(style) - _TOP_LEVEL_SECTIONS
    if unknown:
        raise ValueError(f"Unknown top-level style sections: {sorted(unknown)!r}")
    return style


def list_style_themes() -> tuple[str, ...]:
    theme_dir = resources.files(_THEME_PACKAGE)
    paths = [Path(str(path)) for path in theme_dir.iterdir()]
    return tuple(
        sorted(
            path.stem for path in paths if path.is_file() and path.suffix == ".yaml"
        )
    )


def read_style_yaml(path: str | Path) -> dict[str, Any]:
    content = Path(path).read_text(encoding="utf-8")
    payload = cast(dict[str, Any], yaml.safe_load(content) or {})
    if not isinstance(payload, dict):
        raise ValueError("YAML style payload must be a mapping.")
    return _validate_style_mapping(payload)


def load_style_theme(name: str = _DEFAULT_THEME) -> dict[str, Any]:
    if name not in list_style_themes():
        raise ValueError(f"Unknown style theme {name!r}.")
    theme_dir = resources.files(_THEME_PACKAGE)
    with resources.as_file(theme_dir / f"{name}.yaml") as path:
        return read_style_yaml(path)


def resolve_style(
    *,
    theme: str = _DEFAULT_THEME,
    style_path: str | Path | None = None,
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    base = load_style_theme("base")
    if theme != "base":
        base = _deep_merge(base, load_style_theme(theme))
    if style_path is not None:
        base = _deep_merge(base, read_style_yaml(style_path))
    if overrides is not None:
        base = _deep_merge(base, _validate_style_mapping(overrides))
    return base
