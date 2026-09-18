"""Tests for xrdkit.instrument.

Every scan here is drawn from a Caglioti relation chosen in the test, so the
width fit can be judged against the numbers it was drawn with. No measured
data is read.
"""

import math
from pathlib import Path

import numpy as np
import pytest

from xrdkit import Cell, generate_reflections
from xrdkit.broadening import KALPHA2_INTENSITY_RATIO
from xrdkit.instrument import (
    REFINED_KEYS,
    WIDTH_WINDOW,
    fit_instrument_widths,
    kalpha2_wavelength,
)
from xrdkit.io import XRDScan

# A cubic standard near LaB6, measured with a Cu K alpha doublet.
STANDARD_A = 4.15683
LAMBDA1 = 1.540598
LAMBDA2 = 1.544426

# The relation the patterns below are drawn from, in degrees squared.
DRAWN_U = 0.0060
DRAWN_V = -0.0035
DRAWN_W = 0.0045


def drawn_fwhm(two_theta: float) -> float:
    """The FWHM the Caglioti relation above gives at ``two_theta``, in degrees."""
    tan_theta = math.tan(math.radians(two_theta / 2.0))
    return math.sqrt(DRAWN_U * tan_theta**2 + DRAWN_V * tan_theta + DRAWN_W)


def synthetic_standard(
    two_theta_max: float = 100.0, step: float = 0.013, background: float = 200.0
) -> XRDScan:
    """A cubic Pm-3m standard drawn as K alpha doublets of the drawn widths."""
    two_theta = np.arange(8.0, two_theta_max, step)
    intensity = np.full(two_theta.size, background)
    reflections = generate_reflections(
        Cell.cubic(STANDARD_A), LAMBDA1, two_theta_max - 1.0, space_group="Pm-3m"
    )
    for reflection in reflections:
        centre = reflection.two_theta
        if centre < 10.0:
            continue
        fwhm = drawn_fwhm(centre)
        sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        height = 6000.0 * reflection.multiplicity / 8.0
        intensity += height * np.exp(-0.5 * ((two_theta - centre) / sigma) ** 2)
        # The satellite, where its parent's d spacing puts the longer wavelength.
        theta = math.radians(centre / 2.0)
        sin_two = min(1.0, math.sin(theta) * LAMBDA2 / LAMBDA1)
        satellite = 2.0 * math.degrees(math.asin(sin_two))
        intensity += (
            KALPHA2_INTENSITY_RATIO
            * height
            * np.exp(-0.5 * ((two_theta - satellite) / sigma) ** 2)
        )
    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=LAMBDA1,
        start_angle=float(two_theta[0]),
        end_angle=float(two_theta[-1]),
        step_size=step,
        time_per_step=1.0,
        sample_id="standard",
        source_path="standard.xrdml",
    )


def test_width_fit_recovers_the_relation_it_was_drawn_with() -> None:
    fit = fit_instrument_widths(synthetic_standard())

    # W sets the width at low angle and is the best determined of the three;
    # U and V trade against each other, so they are judged less tightly.
    assert fit.n_peaks >= 6
    assert fit.w == pytest.approx(DRAWN_W, abs=0.0005)
    assert fit.u == pytest.approx(DRAWN_U, abs=0.002)
    assert fit.v == pytest.approx(DRAWN_V, abs=0.003)
    assert fit.rms < 0.01
    assert fit.wavelength == LAMBDA1
    assert fit.window == WIDTH_WINDOW
    assert len(fit.two_theta) == len(fit.fwhm) == fit.n_peaks


def test_width_fit_widths_follow_the_drawn_curve(tmp_path) -> None:
    fit = fit_instrument_widths(synthetic_standard())

    for position, width in zip(fit.two_theta, fit.fwhm):
        assert width == pytest.approx(drawn_fwhm(position), abs=0.004), position


def test_width_fit_honours_the_window() -> None:
    whole = fit_instrument_widths(synthetic_standard())
    narrow = fit_instrument_widths(synthetic_standard(), window=(10.0, 60.0))

    assert narrow.window == (10.0, 60.0)
    assert narrow.n_peaks < whole.n_peaks
    assert max(narrow.two_theta) < 60.0


def test_width_fit_as_a_caglioti_for_write_instprm() -> None:
    fit = fit_instrument_widths(synthetic_standard())

    caglioti = fit.caglioti

    assert (caglioti.u, caglioti.v, caglioti.w) == (fit.u, fit.v, fit.w)
    assert caglioti.n_peaks == fit.n_peaks


def test_width_fit_refuses_a_scan_with_no_wavelength() -> None:
    scan = synthetic_standard()
    bare = XRDScan(**{**scan.__dict__, "wavelength": None})

    with pytest.raises(ValueError, match="carries no wavelength"):
        fit_instrument_widths(bare)


def test_width_fit_refuses_too_few_reflections() -> None:
    with pytest.raises(ValueError, match="the width fit needs three"):
        fit_instrument_widths(synthetic_standard(), window=(10.0, 21.0))


def test_refined_keys_are_the_profile_terms() -> None:
    assert REFINED_KEYS == ("Zero", "U", "V", "W", "X", "Y", "SH/L")


XRDML_WITH_KA2 = """<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">
  <xrdMeasurement>
    <usedWavelength>
      <kAlpha1>1.540598</kAlpha1>
      <kAlpha2>1.544426</kAlpha2>
    </usedWavelength>
  </xrdMeasurement>
</xrdMeasurements>
"""


def test_kalpha2_wavelength_of_an_xrdml(tmp_path) -> None:
    path = tmp_path / "scan.xrdml"
    path.write_text(XRDML_WITH_KA2, encoding="utf-8")

    assert kalpha2_wavelength(path) == pytest.approx(1.544426)


def test_kalpha2_wavelength_of_a_text_pattern(tmp_path) -> None:
    # A two column pattern carries no wavelength at all, let alone two.
    path = tmp_path / "scan.xy"
    path.write_text("10.0 300\n10.02 400\n", encoding="utf-8")

    assert kalpha2_wavelength(path) is None


def test_kalpha2_wavelength_of_an_xrdml_without_one(tmp_path) -> None:
    path = tmp_path / "scan.xrdml"
    path.write_text(XRDML_WITH_KA2.replace("1.544426", ""), encoding="utf-8")

    assert kalpha2_wavelength(path) is None


def test_kalpha2_wavelength_of_an_unreadable_file(tmp_path: Path) -> None:
    path = tmp_path / "broken.xrdml"
    path.write_text("not xml at all", encoding="utf-8")

    assert kalpha2_wavelength(path) is None
