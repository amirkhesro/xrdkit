"""Tests for what the package exports."""

import ast
import importlib
from pathlib import Path

import pytest

import xrdkit

INIT = Path(xrdkit.__file__)


def _reexports() -> dict[str, set[str]]:
    """Every name ``xrdkit/__init__.py`` imports, by the module it comes from."""
    tree = ast.parse(INIT.read_text(encoding="utf-8"))
    imported: dict[str, set[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith(
            "xrdkit."
        ):
            names = {alias.name for alias in node.names}
            imported.setdefault(node.module, set()).update(names)
    return imported


REEXPORTS = _reexports()


def test_the_package_reexports_from_every_module() -> None:
    assert len(REEXPORTS) >= 19


@pytest.mark.parametrize("module", sorted(REEXPORTS))
def test_a_reexported_name_is_in_its_modules_all(module: str) -> None:
    """A name the package re-exports is public, so the module that defines it
    has to say so: ``from xrdkit import x`` and ``from xrdkit.m import x``
    must agree, and Section 31 of the user guide counts a name once, under
    the module whose ``__all__`` holds it."""
    declared = getattr(importlib.import_module(module), "__all__", None)
    assert declared is not None, f"{module} declares no __all__"

    missing = sorted(REEXPORTS[module] - set(declared))
    assert missing == []
