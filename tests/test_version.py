"""Tests for the version xrdkit reports."""

import tomllib
from pathlib import Path

import xrdkit

PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


def test_version_matches_pyproject() -> None:
    with PYPROJECT.open("rb") as handle:
        declared = tomllib.load(handle)["project"]["version"]

    assert xrdkit.__version__ == declared


def test_version_is_exported() -> None:
    assert "__version__" in xrdkit.__all__
