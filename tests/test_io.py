"""Tests for xrdkit.io."""

import os
from pathlib import Path

import numpy as np
import pytest

from xrdkit import XRDScan, read_xrdml

# A folder of measured .xrdml scans to read. No data lives in this repository,
# so the tests that need a real scan skip unless this names one.
RAW_DIR_VARIABLE = "XRDKIT_TEST_RAW_DIR"


def _first_xrdml() -> Path:
    setting = os.environ.get(RAW_DIR_VARIABLE)
    if not setting:
        pytest.skip(f"{RAW_DIR_VARIABLE} is not set to a folder of .xrdml scans")
    raw_dir = Path(setting)
    if not raw_dir.is_dir():
        pytest.skip(f"raw data folder not present: {raw_dir}")
    files = sorted(raw_dir.glob("*.xrdml"))
    if not files:
        pytest.skip(f"no .xrdml files in {raw_dir}")
    return files[0]


def test_read_xrdml_real_file() -> None:
    scan = read_xrdml(_first_xrdml())

    assert isinstance(scan, XRDScan)
    assert scan.two_theta.size == scan.intensity.size
    assert scan.two_theta.size > 0
    assert np.all(np.diff(scan.two_theta) > 0), "two_theta must be strictly increasing"
    assert np.all(scan.intensity >= 0), "intensities must be non-negative"
    assert 1.5 < scan.wavelength < 1.6
