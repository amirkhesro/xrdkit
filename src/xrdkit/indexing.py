"""Indexing of diffraction peaks against a tetragonal cell."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from xrdkit.peaks import Peak

__all__ = [
    "TTB_CELL",
    "IndexedPeak",
    "Reflection",
    "TetragonalCell",
    "generate_reflections",
    "index_peaks",
    "indexed_to_csv",
    "indexing_summary",
]

# Maximum difference between corrected and calculated 2theta, in degrees.
DEFAULT_TOLERANCE = 0.05

# Zero point correction subtracted from every observed position, in degrees.
DEFAULT_ZERO_OFFSET = 0.0

DEFAULT_SPACE_GROUP = "P4bm"

# Space groups whose reflection conditions this module knows about.
SUPPORTED_SPACE_GROUPS = ("P4bm",)

CSV_COLUMNS = (
    "two_theta",
    "corrected_two_theta",
    "d_spacing",
    "relative_intensity",
    "h",
    "k",
    "l",
    "calculated_two_theta",
    "difference",
    "n_candidates",
)

# Decimal places per column when writing CSV; hkl and counts stay integers.
CSV_DECIMALS = {
    "two_theta": 3,
    "corrected_two_theta": 3,
    "d_spacing": 4,
    "relative_intensity": 2,
    "calculated_two_theta": 3,
    "difference": 3,
}


@dataclass
class TetragonalCell:
    """A tetragonal unit cell, in angstroms."""

    a: float
    c: float

    def d_spacing(self, h: int, k: int, l: int) -> float:
        """Return the d spacing of ``(h k l)`` from 1/d^2 = (h^2+k^2)/a^2 + l^2/c^2.

        Raises
        ------
        ValueError
            If ``(h k l)`` is ``(0 0 0)``, which has no d spacing.
        """
        if h == 0 and k == 0 and l == 0:
            raise ValueError("(000) has no d spacing")
        inverse_squared = (h**2 + k**2) / self.a**2 + l**2 / self.c**2
        return float(1.0 / np.sqrt(inverse_squared))


# Provisional starting cell for the tungsten bronze phase; a later refinement
# stage is expected to replace these values.
TTB_CELL = TetragonalCell(a=12.45, c=3.94)


@dataclass
class Reflection:
    """A calculated reflection of a cell at a given wavelength."""

    h: int
    k: int
    l: int
    d_spacing: float
    two_theta: float

    @property
    def hkl(self) -> tuple[int, int, int]:
        """The Miller indices as a tuple."""
        return (self.h, self.k, self.l)


@dataclass
class IndexedPeak:
    """An observed peak with its assignment against a calculated reflection."""

    peak: Peak
    corrected_two_theta: float
    reflection: Reflection | None
    difference: float
    candidates: list[Reflection] = field(default_factory=list)

    @property
    def is_indexed(self) -> bool:
        """Whether a reflection was assigned."""
        return self.reflection is not None


def _is_allowed(h: int, k: int, l: int, space_group: str | None) -> bool:
    """Whether ``(h k l)`` survives the reflection conditions of ``space_group``.

    ``None`` applies no conditions. For P4bm the zonal conditions 0kl with k even
    and h0l with h even also cover the axial 0k0 and h00 cases, which are just
    those zones with the third index zero.

    Raises
    ------
    ValueError
        If ``space_group`` is not one this module knows.
    """
    if space_group is None:
        return True
    if space_group not in SUPPORTED_SPACE_GROUPS:
        raise ValueError(
            f"Unknown space group {space_group!r}; "
            f"expected None or one of {SUPPORTED_SPACE_GROUPS}"
        )

    # P4bm: 0kl present only for k even, h0l only for h even. 00l satisfies both
    # trivially, and general reflections meet neither.
    if h == 0:
        return k % 2 == 0
    if k == 0:
        return h % 2 == 0
    return True


def generate_reflections(
    cell: TetragonalCell,
    wavelength: float,
    two_theta_max: float,
    two_theta_min: float = 0.0,
    space_group: str | None = DEFAULT_SPACE_GROUP,
) -> list[Reflection]:
    """Enumerate the allowed reflections of ``cell`` in a 2theta window.

    Only ``h >= k >= 0`` and ``l >= 0`` are enumerated, so each family that is
    equivalent under 4/mmm Laue symmetry appears exactly once.

    Parameters
    ----------
    cell
        The cell to calculate d spacings from.
    wavelength
        Radiation wavelength in angstroms.
    two_theta_max
        Upper limit of the window, in degrees. Must be above ``two_theta_min``
        and below 180.
    two_theta_min
        Lower limit of the window, in degrees.
    space_group
        Space group whose reflection conditions to apply, or ``None`` for none.

    Returns
    -------
    list[Reflection]
        Reflections inside the window, sorted by 2theta.

    Raises
    ------
    ValueError
        If the window is empty or outside 0 to 180 degrees, or the space group
        is not known.
    """
    if not 0.0 <= two_theta_min < two_theta_max < 180.0:
        raise ValueError(
            f"Need 0 <= two_theta_min < two_theta_max < 180, got "
            f"{two_theta_min} and {two_theta_max}"
        )

    # The smallest d spacing that can diffract within the window fixes how far
    # the indices need to run, so the limits follow from the cell rather than
    # from an arbitrary cut-off.
    d_min = wavelength / (2.0 * np.sin(np.radians(two_theta_max / 2.0)))
    # Along a*, 1/d = h/a, so h can reach a/d_min before falling below d_min.
    h_max = int(np.floor(cell.a / d_min))
    l_max = int(np.floor(cell.c / d_min))

    reflections = []
    for h in range(h_max + 1):
        for k in range(h + 1):
            for l in range(l_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                if not _is_allowed(h, k, l, space_group):
                    continue

                d = cell.d_spacing(h, k, l)
                sin_theta = wavelength / (2.0 * d)
                if sin_theta > 1.0:
                    continue

                two_theta = 2.0 * float(np.degrees(np.arcsin(sin_theta)))
                if two_theta_min <= two_theta <= two_theta_max:
                    reflections.append(
                        Reflection(h=h, k=k, l=l, d_spacing=d, two_theta=two_theta)
                    )

    return sorted(reflections, key=lambda reflection: reflection.two_theta)


def index_peaks(
    peaks: list[Peak],
    cell: TetragonalCell,
    wavelength: float,
    tolerance: float = DEFAULT_TOLERANCE,
    zero_offset: float = DEFAULT_ZERO_OFFSET,
    space_group: str | None = DEFAULT_SPACE_GROUP,
) -> list[IndexedPeak]:
    """Assign reflections of ``cell`` to ``peaks``, in peak order.

    Each observed position has ``zero_offset`` subtracted from it before it is
    compared with the calculated positions.

    Parameters
    ----------
    peaks
        Observed peaks to index.
    cell
        The cell to index against.
    wavelength
        Radiation wavelength in angstroms.
    tolerance
        Largest allowed difference between corrected and calculated 2theta,
        in degrees.
    zero_offset
        Zero point correction, in degrees, subtracted from each observed peak.
    space_group
        Space group whose reflection conditions to apply, or ``None`` for none.

    Returns
    -------
    list[IndexedPeak]
        One entry per peak, in the order given. Peaks with no reflection within
        ``tolerance`` keep ``reflection=None`` and a difference of ``nan``.
    """
    if not peaks:
        return []

    corrected = [peak.two_theta - zero_offset for peak in peaks]
    # Widen the calculated window by the tolerance so peaks at either end can
    # still find a partner just outside the observed range.
    window_min = max(0.0, min(corrected) - tolerance)
    window_max = min(179.0, max(corrected) + tolerance)
    reflections = generate_reflections(
        cell, wavelength, window_max, window_min, space_group
    )

    indexed = []
    for peak, position in zip(peaks, corrected):
        candidates = [
            reflection
            for reflection in reflections
            if abs(position - reflection.two_theta) <= tolerance
        ]
        # Closest first; ties break on hkl so the choice is reproducible.
        candidates.sort(
            key=lambda reflection: (
                abs(position - reflection.two_theta),
                reflection.hkl,
            )
        )
        best = candidates[0] if candidates else None
        indexed.append(
            IndexedPeak(
                peak=peak,
                corrected_two_theta=position,
                reflection=best,
                difference=(
                    position - best.two_theta if best is not None else float("nan")
                ),
                candidates=candidates,
            )
        )
    return indexed


def indexing_summary(indexed: list[IndexedPeak]) -> dict[str, float]:
    """Summarise an indexing run.

    Returns
    -------
    dict
        ``n_peaks``, ``n_indexed``, ``n_unindexed``, ``n_ambiguous`` (peaks with
        more than one candidate within tolerance) and ``rms_difference``, the
        root mean square of the differences of the indexed peaks in degrees.
        The root mean square is ``nan`` when nothing was indexed.
    """
    differences = [entry.difference for entry in indexed if entry.is_indexed]
    rms = (
        float(np.sqrt(np.mean(np.square(differences)))) if differences else float("nan")
    )
    return {
        "n_peaks": len(indexed),
        "n_indexed": len(differences),
        "n_unindexed": len(indexed) - len(differences),
        "n_ambiguous": sum(1 for entry in indexed if len(entry.candidates) > 1),
        "rms_difference": rms,
    }


def indexed_to_csv(indexed: list[IndexedPeak], path: str | Path) -> Path:
    """Write ``indexed`` to ``path`` as CSV with a header, returning the path.

    Unindexed peaks leave the hkl, calculated position and difference blank.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        for entry in indexed:
            reflection = entry.reflection
            values = {
                "two_theta": entry.peak.two_theta,
                "corrected_two_theta": entry.corrected_two_theta,
                "d_spacing": entry.peak.d_spacing,
                "relative_intensity": entry.peak.relative_intensity,
                "h": reflection.h if reflection else "",
                "k": reflection.k if reflection else "",
                "l": reflection.l if reflection else "",
                "calculated_two_theta": reflection.two_theta if reflection else "",
                "difference": entry.difference if reflection else "",
                "n_candidates": len(entry.candidates),
            }
            writer.writerow(
                [
                    f"{values[column]:.{CSV_DECIMALS[column]}f}"
                    if column in CSV_DECIMALS and values[column] != ""
                    else values[column]
                    for column in CSV_COLUMNS
                ]
            )
    return path
