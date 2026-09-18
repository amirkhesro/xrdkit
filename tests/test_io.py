"""Tests for xrdkit.io."""

import os
from pathlib import Path

import numpy as np
import pytest

from xrdkit import XRDScan, read_scan, read_xrdml, read_xy

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


# A minimal PANalytical scan, written by the tests that need read_scan to
# dispatch on the suffix; three points to match the .xy patterns below.
MINIMAL_XRDML = """<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">
  <sample><id>synthetic</id></sample>
  <xrdMeasurement>
    <usedWavelength><kAlpha1>1.540598</kAlpha1></usedWavelength>
    <scan>
      <dataPoints>
        <positions axis="2Theta" unit="deg">
          <startPosition>10.0</startPosition>
          <endPosition>10.04</endPosition>
        </positions>
        <commonCountingTime>1.5</commonCountingTime>
        <counts>300 400 500</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
"""


def test_read_xrdml_real_file() -> None:
    scan = read_xrdml(_first_xrdml())

    assert isinstance(scan, XRDScan)
    assert scan.two_theta.size == scan.intensity.size
    assert scan.two_theta.size > 0
    assert np.all(np.diff(scan.two_theta) > 0), "two_theta must be strictly increasing"
    assert np.all(scan.intensity >= 0), "intensities must be non-negative"
    assert 1.5 < scan.wavelength < 1.6


# A two column pattern, as most diffractometer software and every converter
# writes one. The tests write their own files: no pattern lives in this
# repository, and .gitignore ignores *.xy and *.xye in any case.
PLAIN_XY = "10.00000 300.000\n10.02000 301.500\n10.04000 9000.000\n"


def _write(path: Path, text: str, newline: str = "\n") -> Path:
    """``text`` written to ``path`` with the line ending asked for."""
    path.write_bytes(text.replace("\n", newline).encode("utf-8"))
    return path


def test_read_xy_plain_two_columns_with_crlf(tmp_path) -> None:
    path = _write(tmp_path / "plain.xy", PLAIN_XY, newline="\r\n")

    scan = read_xy(path, wavelength=1.540598)

    assert scan.two_theta == pytest.approx([10.0, 10.02, 10.04])
    assert scan.intensity == pytest.approx([300.0, 301.5, 9000.0])
    assert scan.wavelength == 1.540598
    assert scan.start_angle == pytest.approx(10.0)
    assert scan.end_angle == pytest.approx(10.04)
    assert scan.step_size == pytest.approx(0.02)
    assert scan.time_per_step is None
    assert scan.esd is None
    assert scan.sample_id == "plain"
    assert scan.source_path == str(path)


def test_read_xy_without_a_wavelength(tmp_path) -> None:
    # The reader allows it; the command line is what insists on one.
    scan = read_xy(_write(tmp_path / "bare.xy", PLAIN_XY))

    assert scan.wavelength is None


def test_read_xy_skips_a_header_and_comment_lines(tmp_path) -> None:
    text = (
        "Angle Intensity\n"
        "# exported 2026-09-18\n"
        "; a second comment\n"
        "\n"
        "10.0 300\n"
        "10.02 400\n"
    )
    path = _write(tmp_path / "header.xy", text)

    scan = read_xy(path)

    assert scan.two_theta == pytest.approx([10.0, 10.02])
    assert scan.intensity == pytest.approx([300.0, 400.0])


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("tabs", "10.0\t300\n10.02\t400\n10.04\t500\n"),
        ("commas", "10.0,300\n10.02,400\n10.04,500\n"),
        ("comma and space", "10.0, 300\n10.02, 400\n10.04, 500\n"),
        ("runs of spaces", "10.0    300\n10.02   400\n10.04   500\n"),
    ],
    ids=["tabs", "commas", "comma and space", "runs of spaces"],
)
def test_read_xy_separators(tmp_path, name, text) -> None:
    scan = read_xy(_write(tmp_path / "sep.xy", text))

    assert scan.two_theta == pytest.approx([10.0, 10.02, 10.04]), name
    assert scan.intensity == pytest.approx([300.0, 400.0, 500.0]), name


def test_read_xy_integer_intensities(tmp_path) -> None:
    scan = read_xy(_write(tmp_path / "ints.xy", "10 300\n11 400\n12 500\n"))

    assert scan.intensity.dtype == np.dtype(float)
    assert scan.intensity == pytest.approx([300.0, 400.0, 500.0])
    assert scan.two_theta == pytest.approx([10.0, 11.0, 12.0])


def test_read_xye_third_column_is_the_esd(tmp_path) -> None:
    text = "10.0 300 17.3\n10.02 400 20.0\n10.04 500 22.4\n"

    scan = read_xy(_write(tmp_path / "three.xye", text))

    assert scan.esd == pytest.approx([17.3, 20.0, 22.4])
    assert scan.intensity == pytest.approx([300.0, 400.0, 500.0])


def test_read_xy_blank_lines_inside_the_data(tmp_path) -> None:
    text = "10.0 300\n\n10.02 400\n   \n10.04 500\n\n"

    scan = read_xy(_write(tmp_path / "gaps.xy", text))

    assert scan.two_theta == pytest.approx([10.0, 10.02, 10.04])


def test_read_xy_refuses_a_text_field_once_the_data_has_begun(tmp_path) -> None:
    path = _write(tmp_path / "bad.xy", "10.0 300\n10.02 rubbish\n10.04 500\n")

    with pytest.raises(ValueError) as raised:
        read_xy(path)

    assert str(raised.value) == (f"{path}, line 2: 'rubbish' is not a number")


def test_read_xy_refuses_a_two_theta_that_does_not_increase(tmp_path) -> None:
    path = _write(tmp_path / "back.xy", "10.0 300\n10.02 400\n10.02 500\n")

    with pytest.raises(ValueError) as raised:
        read_xy(path)

    assert str(raised.value) == (
        f"{path}, line 3: two theta 10.02 does not increase on 10.02 before it"
    )


def test_read_xy_refuses_a_line_of_one_number(tmp_path) -> None:
    path = _write(tmp_path / "short.xy", "10.0 300\n10.02\n")

    with pytest.raises(ValueError) as raised:
        read_xy(path)

    assert str(raised.value) == (
        f"{path}, line 2: needs two numbers, a two theta and an intensity"
    )


def test_read_xy_refuses_a_file_with_no_data(tmp_path) -> None:
    path = _write(tmp_path / "empty.xy", "# header only\n\n")

    with pytest.raises(ValueError, match="holds no data"):
        read_xy(path)


def test_read_scan_dispatches_on_the_suffix(tmp_path) -> None:
    xy = _write(tmp_path / "pattern.xy", PLAIN_XY)
    xrdml = tmp_path / "pattern.xrdml"
    xrdml.write_text(MINIMAL_XRDML, encoding="utf-8")

    from_xy = read_scan(xy, wavelength=1.789)
    from_xrdml = read_scan(xrdml)
    upper = read_scan(_write(tmp_path / "shouted.XY", PLAIN_XY))

    assert from_xy.two_theta.size == 3
    assert from_xy.wavelength == 1.789
    assert from_xrdml.wavelength == pytest.approx(1.540598)
    assert from_xrdml.intensity == pytest.approx([300.0, 400.0, 500.0])
    assert upper.two_theta.size == 3


def test_read_scan_overrides_the_wavelength_of_an_xrdml(tmp_path) -> None:
    xrdml = tmp_path / "pattern.xrdml"
    xrdml.write_text(MINIMAL_XRDML, encoding="utf-8")

    assert read_scan(xrdml, wavelength=0.7093).wavelength == 0.7093


def test_read_scan_refuses_an_unknown_suffix(tmp_path) -> None:
    path = tmp_path / "pattern.dat"
    path.write_text(PLAIN_XY, encoding="utf-8")

    with pytest.raises(ValueError) as raised:
        read_scan(path)

    assert str(raised.value) == (
        f"{path}: cannot read a '.dat' scan; xrdkit reads .xrdml, .xy and .xye"
    )
