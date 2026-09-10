"""Peak finding for X-ray diffraction patterns."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import signal

from xrdkit.io import XRDScan

__all__ = ["Peak", "find_peaks", "peaks_to_csv"]

# Minimum peak prominence, as a fraction of the strongest intensity in range.
DEFAULT_MIN_PROMINENCE = 0.02

# Minimum separation between peaks, in degrees.
DEFAULT_MIN_DISTANCE = 0.15

# scipy measures the width at this fraction of the prominence below the peak,
# so 0.5 gives the full width at half maximum.
HALF_PROMINENCE = 0.5

CSV_COLUMNS = (
    "two_theta",
    "intensity",
    "prominence",
    "fwhm",
    "d_spacing",
    "relative_intensity",
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

    Returns
    -------
    list[Peak]
        Peaks ordered by 2theta. ``relative_intensity`` is a percentage of the
        strongest peak found, which is therefore always 100.

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
    return sorted(peaks, key=lambda peak: peak.two_theta)


def peaks_to_csv(peaks: list[Peak], path: str | Path) -> Path:
    """Write ``peaks`` to ``path`` as CSV with a header, returning the path.

    Positions and widths are written to 3 decimal places and d spacings to 4;
    counts keep one place and relative intensities two.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        for peak in peaks:
            values = asdict(peak)
            writer.writerow(
                [f"{values[column]:.{CSV_DECIMALS[column]}f}" for column in CSV_COLUMNS]
            )
    return path
