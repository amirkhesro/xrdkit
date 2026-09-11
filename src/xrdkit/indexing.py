"""Indexing of diffraction peaks against a tetragonal cell."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from xrdkit.peaks import Peak

__all__ = [
    "TTB_CELL",
    "CellFit",
    "IndexedPeak",
    "Reflection",
    "TetragonalCell",
    "ZeroSearch",
    "estimate_zero_offset",
    "generate_reflections",
    "index_and_refine",
    "index_peaks",
    "indexed_to_csv",
    "indexing_summary",
    "refine_cell",
]

# Maximum difference between corrected and calculated 2theta, in degrees.
DEFAULT_TOLERANCE = 0.05

# Zero point correction subtracted from every observed position, in degrees.
DEFAULT_ZERO_OFFSET = 0.0

DEFAULT_SPACE_GROUP = "P4bm"

# Range and step of the automatic zero offset search, in degrees. Wide enough
# for a pellet standing proud of its holder, fine enough to land inside the
# fine tolerance.
DEFAULT_ZERO_SEARCH = (-0.4, 0.4)
DEFAULT_ZERO_SEARCH_STEP = 0.01

# Every trial offset is given a cell of its own, seeded from the peaks below
# this position at this tolerance. Holding one cell for the whole search does
# not work: a cell that is a per cent out shifts the low angle peaks by about
# as much as a zero offset does, so the trial that wins is whichever one best
# papers over the cell error rather than the one that is right.
ZERO_SEED_TOLERANCE = 0.4
ZERO_SEED_TWO_THETA_MAX = 35.0

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


def estimate_zero_offset(
    peaks: list[Peak],
    cell: TetragonalCell,
    wavelength: float,
    search: tuple[float, float] = DEFAULT_ZERO_SEARCH,
    step: float = DEFAULT_ZERO_SEARCH_STEP,
    tolerance: float = DEFAULT_TOLERANCE,
    two_theta_max: float | None = None,
    space_group: str | None = DEFAULT_SPACE_GROUP,
) -> ZeroSearch:
    """Find the zero offset that indexes the most peaks unambiguously.

    A specimen that sits proud of its holder moves every reflection by close to
    a constant, which no cell can absorb: indexed against an uncorrected cell
    such a pattern loses most of its peaks. Trying offsets across ``search`` and
    keeping the one that leaves the most peaks with exactly one candidate
    recovers the shift without anyone having to measure it by hand. Ties go to
    the offset whose matched peaks sit closest to their calculated positions.

    Each trial gets a cell refined for it, from the low angle peaks, rather than
    sharing one. This matters: a constant offset and a cell that is a per cent
    out shift the low angle peaks by much the same amount, so a search that
    holds the cell fixed returns whichever offset best hides the error in that
    cell. Letting the cell move with the offset asks the question that was meant
    instead, which is how much of the pattern each offset can account for.

    Parameters
    ----------
    peaks
        Observed peaks to search on.
    cell
        Cell each trial's own refinement starts from.
    wavelength
        Radiation wavelength in angstroms.
    search
        ``(low, high)`` bounds of the offsets tried, in degrees.
    step
        Spacing of the offsets tried, in degrees.
    tolerance
        Tolerance the trials are scored at, in degrees.
    two_theta_max
        Score on the peaks below this position only. All of them by default,
        since the whole range is what tells an offset from a cell.
    space_group
        Space group whose reflection conditions to apply, or ``None`` for none.

    Returns
    -------
    ZeroSearch
        The best offset, how many peaks it indexed unambiguously, their root
        mean square difference, and the whole profile that was tried.

    Raises
    ------
    ValueError
        If ``search`` is not a rising pair, ``step`` is not positive, or there
        are no peaks to search on.
    """
    low, high = search
    if not low < high:
        raise ValueError(f"Need search low < high, got {search}")
    if step <= 0.0:
        raise ValueError(f"Need a positive step, got {step}")

    candidates = [
        peak
        for peak in peaks
        if two_theta_max is None or peak.two_theta <= two_theta_max
    ]
    if not candidates:
        raise ValueError(
            f"No peak at or below {two_theta_max} degrees to search a zero offset on"
        )
    seed_peaks = [
        peak for peak in candidates if peak.two_theta <= ZERO_SEED_TWO_THETA_MAX
    ]

    trials = np.arange(low, high + step / 2.0, step)
    counts = np.zeros(trials.size, dtype=int)
    deviations = np.full(trials.size, np.inf)

    for index, offset in enumerate(trials):
        trial_cell = _seed_cell(
            seed_peaks, cell, wavelength, float(offset), space_group
        )
        indexed = index_peaks(
            candidates, trial_cell, wavelength, tolerance, float(offset), space_group
        )
        differences = [
            entry.difference
            for entry in indexed
            if entry.is_indexed and len(entry.candidates) == 1
        ]
        counts[index] = len(differences)
        if differences:
            deviations[index] = float(np.sqrt(np.mean(np.square(differences))))

    # Most peaks wins; a tie is broken on how well those peaks actually sit.
    best = int(np.lexsort((deviations, -counts))[0])
    return ZeroSearch(
        offset=float(trials[best]),
        n_indexed=int(counts[best]),
        rms=float(deviations[best]),
        offsets=trials,
        counts=counts,
    )


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


@dataclass
class CellFit:
    """The outcome of a least squares cell refinement."""

    cell: TetragonalCell
    n_peaks: int
    rms_two_theta: float
    c_fitted: bool
    # Zero offset the indexing was run with, in degrees. Zero unless
    # index_and_refine searched for one.
    zero_offset: float = 0.0


@dataclass
class ZeroSearch:
    """The outcome of a scan over trial zero offsets."""

    offset: float
    n_indexed: int
    rms: float
    offsets: np.ndarray
    counts: np.ndarray


def _two_theta(d_spacing: float, wavelength: float) -> float:
    """Return the 2theta of a d spacing, or ``nan`` if it cannot diffract."""
    sin_theta = wavelength / (2.0 * d_spacing)
    if sin_theta > 1.0:
        return float("nan")
    return 2.0 * float(np.degrees(np.arcsin(sin_theta)))


def refine_cell(
    indexed: list[IndexedPeak],
    wavelength: float | None = None,
    start_cell: TetragonalCell | None = None,
) -> CellFit:
    """Refine a tetragonal cell from indexed peaks by linear least squares.

    For a tetragonal cell 1/d^2 = (h^2 + k^2) A + l^2 C with A = 1/a^2 and
    C = 1/c^2, so A and C are the coefficients of an ordinary linear least
    squares problem in the observed 1/d^2. The observed d spacings are
    recomputed from ``corrected_two_theta`` rather than taken from the peaks, so
    that any zero point correction already applied is carried through.

    Parameters
    ----------
    indexed
        Result of :func:`index_peaks`; only entries with an assignment are used.
    wavelength
        Radiation wavelength in angstroms. :class:`IndexedPeak` does not carry
        one, so it must be given here.
    start_cell
        Cell the indexing started from. Only needed when no peak has ``l != 0``,
        in which case its ``c`` is kept unchanged.

    Returns
    -------
    CellFit
        The fitted cell, the number of peaks used, the root mean square
        difference in degrees between their corrected and recalculated 2theta,
        and whether ``c`` was fitted.

    Raises
    ------
    ValueError
        If ``wavelength`` is None, if fewer than three indexed peaks are
        available, if ``c`` cannot be fitted and no ``start_cell`` is given, or
        if the fit returns a coefficient that is not a positive number.
    """
    if wavelength is None:
        raise ValueError(
            "refine_cell needs a wavelength; IndexedPeak does not carry one"
        )

    used = [entry for entry in indexed if entry.is_indexed]
    if len(used) < 3:
        raise ValueError(
            f"Need at least 3 indexed peaks to refine a cell, got {len(used)}"
        )

    positions = np.array([entry.corrected_two_theta for entry in used], dtype=float)
    d_observed = wavelength / (2.0 * np.sin(np.radians(positions / 2.0)))
    inverse_squared = 1.0 / d_observed**2

    hk = np.array(
        [entry.reflection.h**2 + entry.reflection.k**2 for entry in used], dtype=float
    )
    ll = np.array([entry.reflection.l**2 for entry in used], dtype=float)

    c_fitted = bool(np.any(ll > 0.0))
    if c_fitted:
        design = np.column_stack((hk, ll))
    else:
        if start_cell is None:
            raise ValueError(
                "No peak with l != 0, so c cannot be fitted; pass start_cell to "
                "keep its c"
            )
        design = hk.reshape(-1, 1)

    solution, *_ = np.linalg.lstsq(design, inverse_squared, rcond=None)

    a_coefficient = float(solution[0])
    if not a_coefficient > 0.0:
        raise ValueError(
            f"Refinement gave a non-positive 1/a^2 of {a_coefficient}; the "
            "assignments are probably wrong"
        )
    a = 1.0 / np.sqrt(a_coefficient)

    if c_fitted:
        c_coefficient = float(solution[1])
        if not c_coefficient > 0.0:
            raise ValueError(
                f"Refinement gave a non-positive 1/c^2 of {c_coefficient}; the "
                "assignments are probably wrong"
            )
        c = 1.0 / np.sqrt(c_coefficient)
    else:
        c = start_cell.c

    cell = TetragonalCell(a=float(a), c=float(c))
    recalculated = np.array(
        [
            _two_theta(cell.d_spacing(*entry.reflection.hkl), wavelength)
            for entry in used
        ],
        dtype=float,
    )
    rms = float(np.sqrt(np.mean(np.square(positions - recalculated))))

    return CellFit(cell=cell, n_peaks=len(used), rms_two_theta=rms, c_fitted=c_fitted)


def _seed_cell(
    peaks: list[Peak],
    cell: TetragonalCell,
    wavelength: float,
    zero_offset: float,
    space_group: str | None,
) -> TetragonalCell:
    """Refine ``cell`` on the low angle peaks at one trial offset.

    Falls back to the cell it was given whenever there is too little to refine
    on, which leaves that trial to be judged on the starting cell rather than
    dropped.
    """
    indexed = index_peaks(
        peaks, cell, wavelength, ZERO_SEED_TOLERANCE, zero_offset, space_group
    )
    unambiguous = [entry for entry in indexed if len(entry.candidates) == 1]
    try:
        return refine_cell(unambiguous, wavelength, cell).cell
    except ValueError:
        return cell


def index_and_refine(
    peaks: list[Peak],
    start_cell: TetragonalCell,
    wavelength: float,
    zero_offset: float = DEFAULT_ZERO_OFFSET,
    coarse_tolerance: float = 0.4,
    coarse_two_theta_max: float = 35.0,
    fine_tolerance: float = DEFAULT_TOLERANCE,
    space_group: str | None = DEFAULT_SPACE_GROUP,
    n_cycles: int = 2,
    search_zero: bool = True,
) -> tuple[list[IndexedPeak], CellFit]:
    """Index and refine in cycles, starting coarse and low angle.

    The first cycle indexes only the peaks below ``coarse_two_theta_max`` with
    ``coarse_tolerance``, where a cell that is still some way off can be trusted
    to put reflections near the right peaks. Every later cycle indexes the whole
    list with ``fine_tolerance`` against the cell of the previous cycle. Each
    cycle refines on the peaks that matched exactly one reflection, so an
    ambiguous peak never chooses between two candidates on the strength of a
    cell that is not yet converged.

    Parameters
    ----------
    peaks
        Observed peaks to index.
    start_cell
        Cell the first cycle indexes against.
    wavelength
        Radiation wavelength in angstroms.
    zero_offset
        Zero point correction, in degrees, subtracted from each observed peak.
        Giving one turns the search off, since it is then already known.
    coarse_tolerance
        Tolerance of the first cycle, in degrees.
    coarse_two_theta_max
        Upper limit of the first cycle, in corrected degrees.
    fine_tolerance
        Tolerance of the later cycles and of the returned indexing, in degrees.
    space_group
        Space group whose reflection conditions to apply, or ``None`` for none.
    n_cycles
        Number of refinement cycles, at least one.
    search_zero
        Look for the zero offset with :func:`estimate_zero_offset` before
        indexing, starting from ``start_cell`` and scoring at the fine
        tolerance. Only done when ``zero_offset`` is left at zero.

    Returns
    -------
    tuple[list[IndexedPeak], CellFit]
        The whole peak list indexed against the final cell with
        ``fine_tolerance``, and the fit of the last cycle. The fit carries the
        offset everything was indexed with in its ``zero_offset``.

    Raises
    ------
    ValueError
        If ``n_cycles`` is below one, if no peak falls below
        ``coarse_two_theta_max``, or if a cycle has too few peaks to refine on.
    """
    if n_cycles < 1:
        raise ValueError(f"Need at least one cycle, got {n_cycles}")

    if search_zero and zero_offset == 0.0 and peaks:
        # Searched on the same window and tolerance the first cycle will use,
        # so the offset that wins is the one that cycle can act on.
        # Scored over the whole pattern at the fine tolerance: the spread in
        # 2theta is exactly what separates a constant offset from a cell error,
        # so confining the search to the coarse window throws that away.
        zero_offset = estimate_zero_offset(
            peaks,
            start_cell,
            wavelength,
            tolerance=fine_tolerance,
            space_group=space_group,
        ).offset

    coarse_peaks = [
        peak for peak in peaks if peak.two_theta - zero_offset <= coarse_two_theta_max
    ]
    if not coarse_peaks:
        raise ValueError(
            f"No peak below {coarse_two_theta_max} degrees to start the refinement from"
        )

    cell = start_cell
    fit = None
    for cycle in range(n_cycles):
        if cycle == 0:
            candidates = index_peaks(
                coarse_peaks,
                cell,
                wavelength,
                coarse_tolerance,
                zero_offset,
                space_group,
            )
        else:
            candidates = index_peaks(
                peaks, cell, wavelength, fine_tolerance, zero_offset, space_group
            )
        unambiguous = [entry for entry in candidates if len(entry.candidates) == 1]
        fit = refine_cell(unambiguous, wavelength, cell)
        cell = fit.cell

    final = index_peaks(
        peaks, cell, wavelength, fine_tolerance, zero_offset, space_group
    )
    return final, replace(fit, zero_offset=zero_offset)
