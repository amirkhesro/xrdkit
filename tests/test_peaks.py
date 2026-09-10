"""Tests for xrdkit.peaks."""

import csv

import numpy as np
import pytest

from xrdkit import Peak, XRDScan, find_peaks, peaks_to_csv

WAVELENGTH = 1.5406

# position (degrees), height (counts), fwhm (degrees)
SYNTHETIC_PEAKS = (
    (20.000, 8000.0, 0.180),
    (35.500, 5000.0, 0.220),
    (60.250, 3000.0, 0.260),
)

STEP = 0.01
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def make_scan() -> XRDScan:
    """Three Gaussians of known position, height and width on a sloping background."""
    two_theta = np.arange(10.0, 80.0 + STEP / 2, STEP)
    # Sloping background, so the peak finder cannot rely on a flat baseline.
    intensity = 300.0 - 1.2 * (two_theta - two_theta[0])
    for position, height, fwhm in SYNTHETIC_PEAKS:
        sigma = fwhm * FWHM_TO_SIGMA
        intensity += height * np.exp(-((two_theta - position) ** 2) / (2 * sigma**2))

    rng = np.random.default_rng(12345)
    intensity += rng.normal(0.0, 3.0, two_theta.size)

    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=WAVELENGTH,
        start_angle=float(two_theta[0]),
        end_angle=float(two_theta[-1]),
        step_size=STEP,
        time_per_step=1.0,
        sample_id="synthetic",
        source_path="synthetic://three-gaussians",
    )


def test_finds_exactly_the_synthetic_peaks() -> None:
    peaks = find_peaks(make_scan())

    assert len(peaks) == len(SYNTHETIC_PEAKS)
    assert all(isinstance(peak, Peak) for peak in peaks)
    # find_peaks promises results ordered by 2theta.
    assert [peak.two_theta for peak in peaks] == sorted(
        peak.two_theta for peak in peaks
    )


def test_peak_positions_within_hundredth_of_a_degree() -> None:
    peaks = find_peaks(make_scan())
    expected = [position for position, _, _ in SYNTHETIC_PEAKS]

    for peak, position in zip(peaks, expected):
        assert peak.two_theta == pytest.approx(position, abs=0.01)


def test_fwhm_within_fifteen_percent() -> None:
    peaks = find_peaks(make_scan())
    expected = [fwhm for _, _, fwhm in SYNTHETIC_PEAKS]

    for peak, fwhm in zip(peaks, expected):
        assert peak.fwhm == pytest.approx(fwhm, rel=0.15)


def test_d_spacing_matches_bragg() -> None:
    peaks = find_peaks(make_scan())

    for peak in peaks:
        expected = WAVELENGTH / (2.0 * np.sin(np.radians(peak.two_theta / 2.0)))
        assert peak.d_spacing == pytest.approx(expected)

    # d spacing must fall as 2theta rises.
    spacings = [peak.d_spacing for peak in peaks]
    assert spacings == sorted(spacings, reverse=True)


def test_relative_intensity_is_percent_of_strongest() -> None:
    peaks = find_peaks(make_scan())
    strongest = max(peaks, key=lambda peak: peak.intensity)

    assert strongest.relative_intensity == pytest.approx(100.0)
    # The 8000 count peak is the tallest of the three.
    assert strongest.two_theta == pytest.approx(20.0, abs=0.01)
    assert all(peak.relative_intensity <= 100.0 for peak in peaks)

    weakest = min(peaks, key=lambda peak: peak.intensity)
    assert weakest.relative_intensity == pytest.approx(
        100.0 * weakest.intensity / strongest.intensity
    )


def test_two_theta_range_excludes_outside_peaks() -> None:
    scan = make_scan()
    peaks = find_peaks(scan, two_theta_range=(15.0, 45.0))

    assert len(peaks) == 2
    assert [peak.two_theta for peak in peaks] == pytest.approx([20.0, 35.5], abs=0.01)
    # The strongest peak inside the window sets the 100 percent reference.
    assert max(peak.relative_intensity for peak in peaks) == pytest.approx(100.0)


def test_two_theta_range_outside_the_scan_is_rejected() -> None:
    with pytest.raises(ValueError, match="selects no points"):
        find_peaks(make_scan(), two_theta_range=(150.0, 160.0))


def test_min_prominence_filters_weak_peaks() -> None:
    scan = make_scan()
    # The threshold is a fraction of the strongest intensity, near 8300 counts,
    # so 0.5 keeps the 8000 and 5000 peaks and drops the 3000 one.
    assert len(find_peaks(scan, min_prominence=0.02)) == 3
    assert len(find_peaks(scan, min_prominence=0.5)) == 2
    assert len(find_peaks(scan, min_prominence=0.8)) == 1
    assert find_peaks(scan, min_prominence=2.0) == []


def test_peaks_to_csv_writes_a_row_per_peak(tmp_path) -> None:
    peaks = find_peaks(make_scan())
    path = peaks_to_csv(peaks, tmp_path / "tables" / "peaks.csv")

    assert path.is_file()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    header, data = rows[0], rows[1:]
    assert header == [
        "two_theta",
        "intensity",
        "prominence",
        "fwhm",
        "d_spacing",
        "relative_intensity",
    ]
    assert len(data) == len(peaks)
    # Positions to 3 decimals, d spacings to 4.
    assert all(len(row[0].split(".")[1]) == 3 for row in data)
    assert all(len(row[3].split(".")[1]) == 3 for row in data)
    assert all(len(row[4].split(".")[1]) == 4 for row in data)
    assert float(data[0][0]) == pytest.approx(20.0, abs=0.01)


def test_peaks_to_csv_writes_header_only_for_no_peaks(tmp_path) -> None:
    path = peaks_to_csv([], tmp_path / "empty.csv")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    assert len(rows) == 1


def test_refinement_improves_on_the_grid(tmp_path) -> None:
    # A peak deliberately placed off-grid: the parabola vertex must beat the
    # nearest grid point, which can only ever be within half a step.
    two_theta = np.arange(20.0, 30.0 + STEP / 2, STEP)
    centre = 25.0 + 0.004
    sigma = 0.2 * FWHM_TO_SIGMA
    intensity = 100.0 + 5000.0 * np.exp(-((two_theta - centre) ** 2) / (2 * sigma**2))
    scan = XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=WAVELENGTH,
        start_angle=20.0,
        end_angle=30.0,
        step_size=STEP,
        time_per_step=1.0,
        sample_id="offgrid",
        source_path="synthetic://offgrid",
    )

    peak = find_peaks(scan)[0]
    nearest_grid = float(two_theta[np.argmin(np.abs(two_theta - centre))])

    assert abs(peak.two_theta - centre) < abs(nearest_grid - centre)
    assert peak.two_theta == pytest.approx(centre, abs=0.001)
