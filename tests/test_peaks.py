"""Tests for xrdkit.peaks."""

import csv

import numpy as np
import pytest

from xrdkit import (
    Peak,
    XRDScan,
    exclude_kalpha2,
    find_peaks,
    flag_kalpha2,
    peaks_to_csv,
)
from xrdkit.peaks import (
    KALPHA2_INTENSITY_BAND,
    KALPHA2_POSITION_TOLERANCE,
    KALPHA2_RATIO,
)

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
        "kalpha2_of",
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


# A K alpha 2 doublet needs narrower peaks than the widely spaced scan above,
# or the pair merges into one maximum instead of being resolved.
DOUBLET_FWHM = {60.0: 0.08, 85.0: 0.12}

# Height of the isolated peak, then of the two K alpha 1 parents.
ISOLATED_HEIGHT = 9000.0
PARENT_HEIGHTS = {60.0: 6000.0, 85.0: 4000.0}

# A satellite carries about half its parent, which is what flag_kalpha2 looks
# for, so the synthetic pair is built at exactly half.
SATELLITE_FRACTION = 0.5


def satellite_position(two_theta: float) -> float:
    """Where the K alpha 2 satellite of a K alpha 1 peak at ``two_theta`` falls."""
    return 2.0 * float(
        np.degrees(np.arcsin(KALPHA2_RATIO * np.sin(np.radians(two_theta / 2.0))))
    )


def make_doublet_scan() -> XRDScan:
    """One clean peak at 30 degrees plus resolved K alpha 1/2 pairs at 60 and 85."""
    two_theta = np.arange(25.0, 90.0 + STEP / 2, STEP)
    intensity = 300.0 - 1.2 * (two_theta - two_theta[0])

    components = [(30.0, ISOLATED_HEIGHT, 0.20)]
    for parent, height in PARENT_HEIGHTS.items():
        fwhm = DOUBLET_FWHM[parent]
        components.append((parent, height, fwhm))
        components.append(
            (satellite_position(parent), height * SATELLITE_FRACTION, fwhm)
        )

    for position, height, fwhm in components:
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
        sample_id="doublets",
        source_path="synthetic://kalpha2-doublets",
    )


def test_the_doublet_scan_resolves_into_five_peaks() -> None:
    peaks = find_peaks(make_doublet_scan())

    # 30, then the 60 pair, then the 85 pair, in 2theta order.
    assert len(peaks) == 5
    assert [peak.two_theta for peak in peaks] == pytest.approx(
        [
            30.0,
            60.0,
            satellite_position(60.0),
            85.0,
            satellite_position(85.0),
        ],
        abs=0.01,
    )


def test_flags_exactly_the_two_satellites_with_their_parents() -> None:
    peaks = find_peaks(make_doublet_scan())

    # The satellites sit directly above their parents in the list.
    assert peaks[2].kalpha2_of == 1
    assert peaks[4].kalpha2_of == 3


def test_the_parents_and_the_isolated_peak_are_not_flagged() -> None:
    peaks = find_peaks(make_doublet_scan())

    assert peaks[0].kalpha2_of is None
    assert peaks[1].kalpha2_of is None
    assert peaks[3].kalpha2_of is None
    assert sum(peak.kalpha2_of is not None for peak in peaks) == 2


def test_widely_spaced_peaks_keep_no_flags() -> None:
    """The three Gaussian scan has no satellites, and default flagging finds none."""
    peaks = find_peaks(make_scan())

    assert all(peak.kalpha2_of is None for peak in peaks)


def test_flag_satellites_can_be_turned_off() -> None:
    peaks = find_peaks(make_doublet_scan(), flag_satellites=False)

    assert all(peak.kalpha2_of is None for peak in peaks)
    # The same peaks are still found; only the flagging is skipped.
    assert len(peaks) == 5


def test_flag_kalpha2_returns_the_same_list() -> None:
    peaks = find_peaks(make_doublet_scan(), flag_satellites=False)

    assert flag_kalpha2(peaks) is peaks
    assert [peak.kalpha2_of for peak in peaks] == [None, None, 1, None, 3]


def test_flag_kalpha2_rejects_a_bad_intensity_ratio() -> None:
    peaks = find_peaks(make_doublet_scan())

    with pytest.raises(ValueError, match="intensity_ratio"):
        flag_kalpha2(peaks, intensity_ratio=(0.75, 0.25))


def test_a_tight_tolerance_flags_nothing() -> None:
    peaks = find_peaks(make_doublet_scan(), flag_satellites=False)
    # A thousandth of the 0.08 and 0.12 degree FWHM, about the old 0.0001 degrees.
    flag_kalpha2(peaks, position_tolerance=0.001)

    assert all(peak.kalpha2_of is None for peak in peaks)


# A pair at 70 degrees, where the K alpha 2 split is about 0.20 degrees. The
# lines are narrow and finely sampled, so the found positions are good to a few
# thousandths of a degree and the half FWHM position window is 0.02 degrees.
PAIR_PARENT = 70.0
PAIR_HEIGHT = 6000.0
PAIR_FWHM = 0.04
PAIR_STEP = 0.002
PAIR_BACKGROUND = 50.0


def make_pair_scan(
    fraction: float, offset: float = 0.0, background: float = PAIR_BACKGROUND
) -> XRDScan:
    """A peak at 70 degrees and a second peak ``offset`` degrees above its
    calculated K alpha 2 position, ``fraction`` of its height, on a flat
    ``background``."""
    two_theta = np.arange(69.0, 71.0 + PAIR_STEP / 2, PAIR_STEP)
    intensity = np.full(two_theta.size, background)
    sigma = PAIR_FWHM * FWHM_TO_SIGMA
    second = satellite_position(PAIR_PARENT) + offset
    for position, height in (
        (PAIR_PARENT, PAIR_HEIGHT),
        (second, PAIR_HEIGHT * fraction),
    ):
        intensity += height * np.exp(-((two_theta - position) ** 2) / (2 * sigma**2))

    rng = np.random.default_rng(12345)
    intensity += rng.normal(0.0, 3.0, two_theta.size)

    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=WAVELENGTH,
        start_angle=float(two_theta[0]),
        end_angle=float(two_theta[-1]),
        step_size=PAIR_STEP,
        time_per_step=1.0,
        sample_id="pair",
        source_path="synthetic://pair",
    )


def test_a_true_doublet_at_half_intensity_is_flagged() -> None:
    peaks = find_peaks(make_pair_scan(0.5))

    assert len(peaks) == 2
    ratio = peaks[1].intensity / peaks[0].intensity
    # (3000 + 50) / (6000 + 50), inside the band.
    assert ratio == pytest.approx(0.504, abs=0.005)
    assert KALPHA2_INTENSITY_BAND[0] <= ratio <= KALPHA2_INTENSITY_BAND[1]
    gap = peaks[1].two_theta - satellite_position(peaks[0].two_theta)
    assert abs(gap) <= 0.002
    assert abs(gap) <= KALPHA2_POSITION_TOLERANCE * peaks[0].fwhm
    assert [peak.kalpha2_of for peak in peaks] == [None, 0]


def test_a_doublet_on_a_high_background_is_flagged_on_net_heights() -> None:
    # 6000 and 3000 counts on 20000: the raw heights give 23000 / 26000.
    peaks = find_peaks(make_pair_scan(0.5, background=20000.0))

    assert len(peaks) == 2
    assert [peak.kalpha2_of for peak in peaks] == [None, 0]
    assert peaks[1].intensity / peaks[0].intensity == pytest.approx(0.885, abs=0.005)
    assert peaks[0].background == pytest.approx(20000.0, abs=20.0)
    assert peaks[1].background == pytest.approx(20000.0, abs=20.0)
    net = (peaks[1].intensity - peaks[1].background) / (
        peaks[0].intensity - peaks[0].background
    )
    assert net == pytest.approx(0.50, abs=0.01)


def test_two_reflections_of_comparable_intensity_are_not_flagged() -> None:
    # The 4 2 2 beside 5 5 1 case: the right spacing, but 0.85 of the height.
    peaks = find_peaks(make_pair_scan(0.85))

    assert len(peaks) == 2
    # (5100 + 50) / (6000 + 50), above the 0.8 ceiling.
    assert peaks[1].intensity / peaks[0].intensity == pytest.approx(0.851, abs=0.005)
    assert abs(peaks[1].two_theta - satellite_position(peaks[0].two_theta)) <= 0.002
    assert [peak.kalpha2_of for peak in peaks] == [None, None]


def test_the_right_ratio_outside_the_position_window_is_not_flagged() -> None:
    # 0.026 degrees high, against a window of half the 0.04 degree FWHM.
    peaks = find_peaks(make_pair_scan(0.5, offset=0.026))

    assert len(peaks) == 2
    assert peaks[0].fwhm == pytest.approx(PAIR_FWHM, rel=0.05)
    assert peaks[1].intensity / peaks[0].intensity == pytest.approx(0.504, abs=0.005)
    gap = peaks[1].two_theta - satellite_position(peaks[0].two_theta)
    assert gap == pytest.approx(0.026, abs=0.002)
    assert gap > KALPHA2_POSITION_TOLERANCE * peaks[0].fwhm
    assert [peak.kalpha2_of for peak in peaks] == [None, None]


def test_an_instrument_without_kalpha2_flags_nothing() -> None:
    peaks = find_peaks(make_pair_scan(0.5))
    assert [peak.kalpha2_of for peak in peaks] == [None, 0]

    flag_kalpha2(peaks, wavelength_ratio=None)

    assert len(peaks) == 2
    assert [peak.kalpha2_of for peak in peaks] == [None, None]


def make_peak(two_theta: float, intensity: float, fwhm: float) -> Peak:
    return Peak(
        two_theta=two_theta,
        intensity=intensity,
        prominence=intensity,
        fwhm=fwhm,
        d_spacing=WAVELENGTH / (2.0 * np.sin(np.radians(two_theta / 2.0))),
        relative_intensity=100.0,
    )


def test_the_parent_is_the_nearest_peak_below_in_the_window() -> None:
    # The third peak sits exactly on the first's satellite position and 0.030
    # degrees below the second's, inside the second's 0.10 degree window.
    candidate = satellite_position(70.0)
    peaks = [
        make_peak(70.0, 6000.0, 0.10),
        make_peak(70.03, 6000.0, 0.20),
        make_peak(candidate, 3000.0, 0.10),
    ]
    assert satellite_position(70.03) - candidate == pytest.approx(0.030, abs=0.001)

    flag_kalpha2(peaks)
    assert [peak.kalpha2_of for peak in peaks] == [None, None, 1]

    # When that nearest parent fails the intensity band, at 3000 / 3200, the
    # peak is not handed on to the farther one.
    peaks[1].intensity = 3200.0
    flag_kalpha2(peaks)
    assert [peak.kalpha2_of for peak in peaks] == [None, None, None]


def test_exclude_kalpha2_removes_exactly_the_satellites() -> None:
    peaks = find_peaks(make_doublet_scan())
    kept = exclude_kalpha2(peaks)

    assert len(peaks) == 5
    assert len(kept) == 3
    assert [peak.two_theta for peak in kept] == pytest.approx(
        [30.0, 60.0, 85.0], abs=0.01
    )
    assert all(peak.kalpha2_of is None for peak in kept)

    # The input is left alone, flags and all.
    assert [peak.kalpha2_of for peak in peaks] == [None, None, 1, None, 3]
    assert all(kept_peak is not peak for kept_peak in kept for peak in peaks)


def test_exclude_kalpha2_renormalises_relative_intensity() -> None:
    kept = exclude_kalpha2(find_peaks(make_doublet_scan()))

    strongest = max(peak.intensity for peak in kept)
    for peak in kept:
        assert peak.relative_intensity == pytest.approx(
            100.0 * peak.intensity / strongest
        )
    assert max(peak.relative_intensity for peak in kept) == pytest.approx(100.0)


def test_exclude_kalpha2_of_nothing() -> None:
    assert exclude_kalpha2([]) == []


def test_peaks_to_csv_writes_the_kalpha2_column(tmp_path) -> None:
    peaks = find_peaks(make_doublet_scan())
    path = peaks_to_csv(peaks, tmp_path / "doublets.csv")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    header, data = rows[0], rows[1:]

    assert header[-1] == "kalpha2_of"
    # Blank for a peak that is not a satellite, the parent's index otherwise.
    assert [row[-1] for row in data] == ["", "", "1", "", "3"]
