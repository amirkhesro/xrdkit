"""Peak finding for X-ray diffraction patterns."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from scipy import signal

from xrdkit.io import XRDScan

__all__ = [
    "Peak",
    "exclude_kalpha2",
    "find_peaks",
    "flag_kalpha2",
    "peaks_to_csv",
]

# Minimum peak prominence, as a fraction of the strongest intensity in range.
DEFAULT_MIN_PROMINENCE = 0.02

# Minimum separation between peaks, in degrees.
DEFAULT_MIN_DISTANCE = 0.15

# scipy measures the width at this fraction of the prominence below the peak,
# so 0.5 gives the full width at half maximum.
HALF_PROMINENCE = 0.5

# The background under a peak is the lowest intensity within this many degrees
# either side of it. That is wide enough to reach past the peak and a resolved
# K alpha 2 satellite to the floor between reflections, and narrow enough to
# follow a sloping background; a satellite and its parent, a fraction of a
# degree apart, share nearly the same window and so the same background.
BACKGROUND_HALF_WIDTH = 1.0

# Copper K alpha 1 and K alpha 2 wavelengths, in angstroms. A satellite
# diffracts at the same d spacing as its parent, so Bragg's law puts it at
# 2theta_2 = 2 arcsin(KALPHA2_RATIO sin(theta_1)), always to high angle.
KALPHA1_WAVELENGTH = 1.54056
KALPHA2_WAVELENGTH = 1.54439
KALPHA2_RATIO = KALPHA2_WAVELENGTH / KALPHA1_WAVELENGTH

# Largest gap between a peak and the satellite position predicted for its
# parent, as a fraction of the parent's FWHM.
KALPHA2_POSITION_TOLERANCE = 0.5

# Allowed range of a satellite's intensity over its parent's. The theoretical
# ratio is 0.5; a genuine reflection that happens to sit at the K alpha 2
# spacing above a neighbour of comparable height falls above the ceiling.
KALPHA2_INTENSITY_BAND = (0.2, 0.8)

CSV_COLUMNS = (
    "two_theta",
    "intensity",
    "prominence",
    "fwhm",
    "d_spacing",
    "relative_intensity",
    "kalpha2_of",
)

# Decimal places per column when writing CSV; counts keep one place.
CSV_DECIMALS = {
    "two_theta": 3,
    "intensity": 1,
    "prominence": 1,
    "fwhm": 3,
    "d_spacing": 4,
    "relative_intensity": 2,
}


@dataclass
class Peak:
    """A single reflection located in a diffraction pattern."""

    two_theta: float
    intensity: float
    prominence: float
    fwhm: float
    d_spacing: float
    relative_intensity: float
    # Index in the peak list of the K alpha 1 parent this peak is a satellite
    # of, or None if it is not judged to be one.
    kalpha2_of: int | None = None
    # Background counts under the peak: the lowest intensity within
    # BACKGROUND_HALF_WIDTH either side of it, over the points find_peaks
    # searched, which is what _local_background returns. It stays 0.0 on a
    # Peak built by hand, since only find_peaks estimates it.
    background: float = 0.0


def _step_size(scan: XRDScan) -> float:
    """Return the 2theta step of ``scan``, falling back to the axis spacing."""
    if scan.step_size > 0.0:
        return float(scan.step_size)
    if scan.two_theta.size > 1:
        return float(np.median(np.diff(scan.two_theta)))
    raise ValueError("Cannot determine the 2theta step of a single point scan")


def _d_spacing(two_theta: float, wavelength: float) -> float:
    """Return the d spacing in angstroms from Bragg's law."""
    theta = np.radians(two_theta / 2.0)
    sin_theta = float(np.sin(theta))
    if sin_theta <= 0.0:
        return float("nan")
    return wavelength / (2.0 * sin_theta)


def _refine_position(
    two_theta: np.ndarray, intensity: np.ndarray, index: int, step: float
) -> float:
    """Refine a peak position by the vertex of a parabola through three points.

    Falls back to the grid position when the peak sits on the first or last
    point, when the three points are collinear, or when the vertex lands outside
    the neighbouring points, since none of those give a meaningful vertex.
    """
    grid_position = float(two_theta[index])
    if index <= 0 or index >= intensity.size - 1:
        return grid_position

    left, centre, right = (float(value) for value in intensity[index - 1 : index + 2])
    curvature = left - 2.0 * centre + right
    if curvature == 0.0:
        return grid_position

    shift = 0.5 * (left - right) / curvature
    if not -1.0 < shift < 1.0:
        return grid_position
    return grid_position + shift * step


def _local_background(
    two_theta: np.ndarray, intensity: np.ndarray, indices: np.ndarray
) -> list[float]:
    """The lowest intensity within BACKGROUND_HALF_WIDTH of each peak index,
    on a rising 2theta axis."""
    centres = two_theta[indices]
    starts = np.searchsorted(two_theta, centres - BACKGROUND_HALF_WIDTH, "left")
    ends = np.searchsorted(two_theta, centres + BACKGROUND_HALF_WIDTH, "right")
    return [float(np.min(intensity[s:e])) for s, e in zip(starts, ends)]


def find_peaks(
    scan: XRDScan,
    min_prominence: float = DEFAULT_MIN_PROMINENCE,
    min_distance: float = DEFAULT_MIN_DISTANCE,
    two_theta_range: tuple[float, float] | None = None,
    flag_satellites: bool = True,
) -> list[Peak]:
    """Locate reflections in ``scan``, sorted by 2theta.

    Parameters
    ----------
    scan
        The scan to search.
    min_prominence
        Minimum prominence a peak must have, as a fraction of the strongest
        intensity within the searched range.
    min_distance
        Minimum separation between peaks, in degrees.
    two_theta_range
        Optional ``(low, high)`` window, in degrees, to restrict the search
        to. The window is applied to the points before anything is measured,
        so both the prominence threshold and ``relative_intensity`` are taken
        from the strongest intensity inside it: narrowing the window changes
        which peaks clear ``min_prominence`` and changes what 100 per cent
        means.
    flag_satellites
        Whether to run :func:`flag_kalpha2` over the result, which sets
        ``kalpha2_of`` on the peaks that look like K alpha 2 satellites.

    Returns
    -------
    list[Peak]
        Peaks ordered by 2theta. ``relative_intensity`` is a percentage of the
        strongest peak found, which is therefore always 100. Satellites are
        still returned; they are flagged rather than dropped, so that
        :func:`exclude_kalpha2` can be applied when it suits.

    Raises
    ------
    ValueError
        If ``two_theta_range`` selects no data, or the step cannot be determined.
    """
    two_theta = np.asarray(scan.two_theta, dtype=float)
    intensity = np.asarray(scan.intensity, dtype=float)

    if two_theta_range is not None:
        low, high = two_theta_range
        selected = (two_theta >= low) & (two_theta <= high)
        if not np.any(selected):
            raise ValueError(
                f"two_theta_range {two_theta_range} selects no points of "
                f"{two_theta[0]:.3f}-{two_theta[-1]:.3f} degrees"
            )
        two_theta = two_theta[selected]
        intensity = intensity[selected]

    if intensity.size == 0:
        return []

    step = _step_size(scan)
    prominence = min_prominence * float(np.max(intensity))
    # scipy counts separation in points, and needs at least one.
    distance = max(1, round(min_distance / step))

    indices, properties = signal.find_peaks(
        intensity, prominence=prominence, distance=distance
    )
    if indices.size == 0:
        return []

    widths = signal.peak_widths(intensity, indices, rel_height=HALF_PROMINENCE)[0]

    heights = intensity[indices]
    strongest = float(np.max(heights))
    backgrounds = _local_background(two_theta, intensity, indices)

    peaks = []
    for index, height, peak_prominence, width, background in zip(
        indices, heights, properties["prominences"], widths, backgrounds
    ):
        position = _refine_position(two_theta, intensity, int(index), step)
        peaks.append(
            Peak(
                two_theta=position,
                intensity=float(height),
                prominence=float(peak_prominence),
                fwhm=float(width) * step,
                d_spacing=_d_spacing(position, scan.wavelength),
                relative_intensity=100.0 * float(height) / strongest,
                background=background,
            )
        )
    peaks.sort(key=lambda peak: peak.two_theta)
    if flag_satellites:
        flag_kalpha2(peaks)
    return peaks


def _satellite_position(two_theta: float, wavelength_ratio: float) -> float:
    """Return where the K alpha 2 satellite of a peak at ``two_theta`` falls.

    ``nan`` when the scaled sine passes 1, which is beyond the reach of the
    longer wavelength and so has no satellite.
    """
    sin_theta = wavelength_ratio * np.sin(np.radians(two_theta / 2.0))
    if sin_theta > 1.0:
        return float("nan")
    return 2.0 * float(np.degrees(np.arcsin(sin_theta)))


def flag_kalpha2(
    peaks: list[Peak],
    position_tolerance: float = KALPHA2_POSITION_TOLERANCE,
    intensity_ratio: tuple[float, float] = KALPHA2_INTENSITY_BAND,
    wavelength_ratio: float | None = KALPHA2_RATIO,
) -> list[Peak]:
    """Flag the peaks that look like K alpha 2 satellites, in place.

    A satellite sits at the position the parent's d spacing would give at the
    longer K alpha 2 wavelength, and carries roughly half the parent's
    intensity. A peak is flagged only when both hold. Its parent is the nearest
    peak below it in 2theta whose predicted satellite position lies within
    ``position_tolerance`` times that peak's FWHM of it; the peak is then
    flagged only if its intensity over that parent's lies in
    ``intensity_ratio``. A peak whose nearest parent fails the intensity test
    is not handed on to a farther one. Position alone is not enough: a genuine
    reflection can sit at the K alpha 2 spacing above a neighbour, and is told
    apart by being of comparable height.

    Intensities are the heights above the background, ``Peak.intensity -
    Peak.background``, for both peaks. Raw heights on a high background give a
    ratio pulled towards one. Prominences are no substitute: a satellite on its
    parent's tail has its prominence measured down to the saddle between the
    two, which pulls the ratio the other way.

    Only resolved satellites can be caught this way. Below about 50 degrees the
    pair is not separated enough for the peak finder to report two peaks, so
    there is nothing to flag.

    Parameters
    ----------
    peaks
        Peaks to flag, as returned by :func:`find_peaks`. Modified in place.
    position_tolerance
        Largest allowed gap between a peak and its parent's predicted
        satellite position, as a fraction of the parent's FWHM.
    intensity_ratio
        Allowed ``(low, high)`` range of the peak's height above the background
        over its parent's.
    wavelength_ratio
        K alpha 2 over K alpha 1, 1.002486 for copper, or None for radiation
        without K alpha 2, which flags nothing.

    Returns
    -------
    list[Peak]
        The same list, for chaining.

    Raises
    ------
    ValueError
        If ``intensity_ratio`` is not a rising pair of positive numbers.
    """
    low, high = intensity_ratio
    if not 0.0 < low < high:
        raise ValueError(
            f"Need 0 < low < high for intensity_ratio, got {intensity_ratio}"
        )

    if wavelength_ratio is None:
        for peak in peaks:
            peak.kalpha2_of = None
        return peaks

    # Predicting once per peak keeps this a single pass over the pairs.
    predicted = [
        _satellite_position(peak.two_theta, wavelength_ratio) for peak in peaks
    ]

    for peak in peaks:
        parent_index = None
        for index, other in enumerate(peaks):
            if other.two_theta >= peak.two_theta:
                continue
            # A nan prediction, beyond the reach of K alpha 2, fails this.
            if not abs(peak.two_theta - predicted[index]) <= (
                position_tolerance * other.fwhm
            ):
                continue
            # A tie in 2theta keeps the earlier peak.
            if parent_index is None or other.two_theta > peaks[parent_index].two_theta:
                parent_index = index
        peak.kalpha2_of = None
        if parent_index is not None:
            parent = peaks[parent_index]
            parent_net = parent.intensity - parent.background
            net = peak.intensity - peak.background
            if parent_net > 0.0 and low <= net / parent_net <= high:
                peak.kalpha2_of = parent_index

    return peaks


def exclude_kalpha2(peaks: list[Peak]) -> list[Peak]:
    """Return a new list without the flagged satellites.

    The survivors are copied, not shared, and their ``relative_intensity`` is
    recomputed against the strongest of them, so it again reaches 100.
    """
    kept = [replace(peak) for peak in peaks if peak.kalpha2_of is None]
    if not kept:
        return []

    strongest = max(peak.intensity for peak in kept)
    for peak in kept:
        peak.relative_intensity = 100.0 * peak.intensity / strongest
    return kept


def peaks_to_csv(peaks: list[Peak], path: str | Path) -> Path:
    """Write ``peaks`` to ``path`` as CSV with a header, returning the path.

    Positions and widths are written to 3 decimal places and d spacings to 4;
    counts keep one place and relative intensities two. ``kalpha2_of`` is the
    parent's index, or blank for a peak that is not a satellite.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        for peak in peaks:
            values = asdict(peak)
            writer.writerow(
                [
                    f"{values[column]:.{CSV_DECIMALS[column]}f}"
                    if column in CSV_DECIMALS
                    else ("" if values[column] is None else values[column])
                    for column in CSV_COLUMNS
                ]
            )
    return path
