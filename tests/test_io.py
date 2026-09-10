"""Tests for xrdkit.io."""

from pathlib import Path

import numpy as np
import pytest

from xrdkit import XRDScan, read_xrdml

RAW_DIR = Path(r"C:\Users\amirk\Source\repos\XRD-Analysis\data\raw")


def _first_xrdml() -> Path:
    if not RAW_DIR.is_dir():
        pytest.skip(f"raw data folder not present: {RAW_DIR}")
    files = sorted(RAW_DIR.glob("*.xrdml"))
    if not files:
        pytest.skip(f"no .xrdml files in {RAW_DIR}")
    return files[0]


def test_read_xrdml_real_file() -> None:
    scan = read_xrdml(_first_xrdml())

    assert isinstance(scan, XRDScan)
    assert scan.two_theta.size == scan.intensity.size
    assert scan.two_theta.size > 0
    assert np.all(np.diff(scan.two_theta) > 0), "two_theta must be strictly increasing"
    assert np.all(scan.intensity >= 0), "intensities must be non-negative"
    assert 1.5 < scan.wavelength < 1.6
