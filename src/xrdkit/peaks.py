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

# Copper K alpha 1 and K alpha 2 wavelengths, in angstroms. A satellite
# diffracts at the same d spacing as its parent, so Bragg's law puts it at
# 2theta_2 = 2 arcsin(KALPHA2_RATIO sin(theta_1)), always to high angle.
KALPHA1_WAVELENGTH = 1.54056
KALPHA2_WAVELENGTH = 1.54439
KALPHA2_RATIO = KALPHA2_WAVELENGTH / KALPHA1_WAVELENGTH

# Largest gap between a peak and the satellite position predicted for it, in
# degrees.
DEFAULT_KALPHA2_TOLERANCE = 0.03

# A K alpha 2 satellite carries about half the intensity of its parent, but a
# satellite resolved on the tail of its parent measures well above that, so the
# ceiling is set high enough to admit the partly overlapped ones.
DEFAULT_KALPHA2_INTENSITY_RATIO = (0.25, 0.90)

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
        Optional ``(low, high)`` window, in degrees, to restrict the search to.
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

    peaks = []
    for index, height, peak_prominence, width in zip(
        indices, heights, properties["prominences"], widths
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
    tolerance: float = DEFAULT_KALPHA2_TOLERANCE,
    intensity_ratio: tuple[float, float] = DEFAULT_KALPHA2_INTENSITY_RATIO,
    wavelength_ratio: float = KALPHA2_RATIO,
) -> list[Peak]:
    """Flag the peaks that look like K alpha 2 satellites, in place.

    A satellite sits at the position the parent's d spacing would give at the
    longer K alpha 2 wavelength, and carries roughly half the parent's
    intensity. Every stronger peak is tested as a possible parent and the
    closest one that satisfies both conditions wins, so a peak that happens to
    fall near two predicted positions is attributed to the better match.

    Only resolved satellites can be caught this way. Below about 50 degrees the
    pair is not separated enough for the peak finder to report two peaks, so
    there is nothing to flag. A satellite that is only just resolved, still
    sitting on the tail of its parent, measures high because its height is taken
    above a baseline the parent has raised, which is why the default ceiling on
    the intensity ratio is well above the half that clean separation would give.

    Parameters
    ----------
    peaks
        Peaks to flag, as returned by :func:`find_peaks`. Modified in place.
    tolerance
        Largest allowed gap between a peak and a predicted satellite position,
        in degrees.
    intensity_ratio
        Allowed ``(low, high)`` range of the peak's intensity over its parent's.
    wavelength_ratio
        K alpha 2 over K alpha 1, 1.002486 for copper.

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

    # Predicting once per peak keeps this a single pass over the pairs.
    predicted = [
        _satellite_position(peak.two_theta, wavelength_ratio) for peak in peaks
    ]

    for index, peak in enumerate(peaks):
        best_parent = None
        best_gap = float("inf")
        for parent_index, parent in enumerate(peaks):
            if parent_index == index or parent.intensity <= peak.intensity:
                continue
            if not low <= peak.intensity / parent.intensity <= high:
                continue
            gap = abs(peak.two_theta - predicted[parent_index])
            # A tie keeps the earlier parent, which is the lower angle one.
            if gap <= tolerance and gap < best_gap:
                best_parent = parent_index
                best_gap = gap
        peak.kalpha2_of = best_parent

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
