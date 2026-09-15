"""Indexing of diffraction peaks against a unit cell of any crystal system.

Reflections come from :class:`~xrdkit.cell.Cell` and :mod:`xrdkit.symmetry`:
each family of reflections equivalent under the Laue group appears once,
labelled by :func:`~xrdkit.symmetry.representative` and carrying its
multiplicity, and a space group, when given, removes the systematic absences.
"""

from __future__ import annotations

import csv
import itertools
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from xrdkit.cell import Cell
from xrdkit.library import CELL_PARAMETERS
from xrdkit.peaks import Peak
from xrdkit.symmetry import (
    SUPPORTED_SPACE_GROUPS,
    holohedry,
    is_absent,
    laue_group,
    laue_orbit,
    representative,
    space_group_operations,
)

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

# The space group the command line assumes for a tetragonal start cell given
# without one. The functions here default to no space group.
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


# Each component of the reciprocal metric G* a crystal system leaves free:
# its name, the cell parameter it stands for, and the symmetric matrix M it adds
# to G* per unit value, so that its column in 1/d^2 = h . G* h is h . M h. The
# names follow 1/d^2 = h^2 A + k^2 B + l^2 C + 2kl D + 2hl E + 2hk F.
_A = ((1, 0, 0), (0, 0, 0), (0, 0, 0))
_B = ((0, 0, 0), (0, 1, 0), (0, 0, 0))
_C = ((0, 0, 0), (0, 0, 0), (0, 0, 1))
_D = ((0, 0, 0), (0, 0, 1), (0, 1, 0))
_E = ((0, 0, 1), (0, 0, 0), (1, 0, 0))
_F = ((0, 1, 0), (1, 0, 0), (0, 0, 0))
_HEXAGONAL_A = ((1, 0.5, 0), (0.5, 1, 0), (0, 0, 0))
RECIPROCAL_COMPONENTS = {
    "cubic": (("A", "a", ((1, 0, 0), (0, 1, 0), (0, 0, 1))),),
    "tetragonal": (("A", "a", ((1, 0, 0), (0, 1, 0), (0, 0, 0))), ("C", "c", _C)),
    "orthorhombic": (("A", "a", _A), ("B", "b", _B), ("C", "c", _C)),
    "hexagonal": (("A", "a", _HEXAGONAL_A), ("C", "c", _C)),
    "trigonal": (("A", "a", _HEXAGONAL_A), ("C", "c", _C)),
    "monoclinic": (("A", "a", _A), ("B", "b", _B), ("C", "c", _C), ("E", "beta", _E)),
    "triclinic": (
        ("A", "a", _A),
        ("B", "b", _B),
        ("C", "c", _C),
        ("D", "alpha", _D),
        ("E", "beta", _E),
        ("F", "gamma", _F),
    ),
}

# The crystal system refine_cell assumes when it is given no start cell.
DEFAULT_CRYSTAL_SYSTEM = "tetragonal"


def TetragonalCell(a: float, c: float) -> Cell:
    """A tetragonal :class:`~xrdkit.cell.Cell`, in angstroms.

    Deprecated: kept so that existing code keeps working. Use
    :meth:`Cell.tetragonal <xrdkit.cell.Cell.tetragonal>` instead.
    """
    return Cell.tetragonal(a, c)


# Provisional starting cell for the tungsten bronze phase; a later refinement
# stage is expected to replace these values.
TTB_CELL = Cell.tetragonal(12.45, 3.94)


@dataclass
class Reflection:
    """A calculated reflection of a cell at a given wavelength."""

    h: int
    k: int
    l: int
    d_spacing: float
    two_theta: float
    # Reflections in the family under the Laue group; 1 for one built by hand.
    multiplicity: int = 1

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


@lru_cache(maxsize=64)
def _families(
    limits: tuple[int, int, int], space_group: str | None, crystal_system: str
) -> tuple[np.ndarray, np.ndarray]:
    """The representative and multiplicity of every family of reflections with
    a member in the signed box ``limits``, absences removed.

    It depends on the box and the symmetry but not on the cell, so it is
    cached: a zero offset search asks for the same box many times over.
    """
    if space_group is None:
        operations = None
        laue = holohedry(crystal_system)
    else:
        operations = space_group_operations(space_group)
        laue = laue_group(operations)

    seen = set()
    representatives = []
    multiplicities = []
    for hkl in itertools.product(*(range(-limit, limit + 1) for limit in limits)):
        if hkl in seen or not any(hkl):
            continue
        orbit = laue_orbit(hkl, laue)
        seen.update(orbit)
        # Equivalent reflections share their absence, so one member decides.
        if operations is not None and is_absent(hkl, operations):
            continue
        representatives.append(representative(hkl, laue))
        multiplicities.append(len(orbit))

    hkl_array = np.array(representatives, dtype=int).reshape(-1, 3)
    multiplicity_array = np.array(multiplicities, dtype=int)
    hkl_array.flags.writeable = False
    multiplicity_array.flags.writeable = False
    return hkl_array, multiplicity_array


def generate_reflections(
    cell: Cell,
    wavelength: float,
    two_theta_max: float,
    two_theta_min: float = 0.0,
    space_group: str | None = None,
) -> list[Reflection]:
    """Enumerate the allowed reflections of ``cell`` in a 2theta window.

    Every h, k and l in the signed box that can reach the window is tried.
    Reflections equivalent under the Laue group of ``space_group``, or of the
    holohedry of ``cell.crystal_system`` without one, are merged into one
    :class:`Reflection`, labelled by :func:`~xrdkit.symmetry.representative`
    (h >= k >= 0 and l >= 0 for a tetragonal cell) and carrying the family's
    multiplicity.

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
        Space group whose systematic absences to remove, one of
        :data:`~xrdkit.symmetry.SUPPORTED_SPACE_GROUPS`, or ``None`` for none.

    Returns
    -------
    list[Reflection]
        Reflections inside the window, sorted by 2theta and then hkl.

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
    if space_group is not None and space_group not in SUPPORTED_SPACE_GROUPS:
        raise ValueError(
            f"Unknown space group {space_group!r}; "
            f"expected None or one of {SUPPORTED_SPACE_GROUPS}"
        )

    # The smallest d spacing that can diffract within the window fixes how far
    # the indices need to run, so the limits follow from the cell rather than
    # from an arbitrary cut-off. |h| = |d* . a| <= a / d, so h can reach
    # a / d_min before falling below d_min, and likewise k with b and l with c.
    d_min = wavelength / (2.0 * np.sin(np.radians(two_theta_max / 2.0)))
    limits = tuple(int(np.floor(length / d_min)) for length in (cell.a, cell.b, cell.c))
    hkl, multiplicities = _families(limits, space_group, cell.crystal_system)
    if not len(hkl):
        return []

    d_spacings = cell.d_spacings(hkl)
    sin_theta = wavelength / (2.0 * d_spacings)
    reachable = sin_theta <= 1.0
    two_theta = np.full(len(hkl), np.nan)
    two_theta[reachable] = 2.0 * np.degrees(np.arcsin(sin_theta[reachable]))
    inside = reachable & (two_theta >= two_theta_min) & (two_theta <= two_theta_max)

    reflections = [
        Reflection(
            h=int(h),
            k=int(k),
            l=int(l),
            d_spacing=float(d),
            two_theta=float(position),
            multiplicity=int(multiplicity),
        )
        for (h, k, l), d, position, multiplicity in zip(
            hkl[inside], d_spacings[inside], two_theta[inside], multiplicities[inside]
        )
    ]
    return sorted(
        reflections, key=lambda reflection: (reflection.two_theta, reflection.hkl)
    )


def estimate_zero_offset(
    peaks: list[Peak],
    cell: Cell,
    wavelength: float,
    search: tuple[float, float] = DEFAULT_ZERO_SEARCH,
    step: float = DEFAULT_ZERO_SEARCH_STEP,
    tolerance: float = DEFAULT_TOLERANCE,
    two_theta_max: float | None = None,
    space_group: str | None = None,
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
    cell: Cell,
    wavelength: float,
    tolerance: float = DEFAULT_TOLERANCE,
    zero_offset: float = DEFAULT_ZERO_OFFSET,
    space_group: str | None = None,
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

    cell: Cell
    n_peaks: int
    rms_two_theta: float
    # The cell parameters whose reciprocal metric component no peak carried
    # information on, held at the start cell's value; empty when all were fitted.
    held: tuple[str, ...]
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
    start_cell: Cell | None = None,
) -> CellFit:
    """Refine a cell of any crystal system from indexed peaks by linear least
    squares on its reciprocal metric.

    1/d^2 = h^2 A + k^2 B + l^2 C + 2kl D + 2hl E + 2hk F, linear in the
    components of the reciprocal metric G*. The crystal system of
    ``start_cell`` decides which are free: A alone for a cubic cell, with the
    column h^2 + k^2 + l^2; A with h^2 + k^2 and C for a tetragonal one; A, B
    and C for orthorhombic; A with h^2 + hk + k^2 and C for hexagonal and
    trigonal; A, B, C and E for monoclinic (b unique); all six for triclinic.
    A component no peak carries information on (its column is all zero) is
    held at the value ``start_cell`` gives it. The fitted G* is inverted to the
    direct metric, from which the cell parameters are read. The observed d
    spacings are recomputed from ``corrected_two_theta`` rather than taken
    from the peaks, so that any zero point correction already applied is
    carried through.

    Parameters
    ----------
    indexed
        Result of :func:`index_peaks`; only entries with an assignment are used.
    wavelength
        Radiation wavelength in angstroms. :class:`IndexedPeak` does not carry
        one, so it must be given here.
    start_cell
        Cell the indexing started from, which gives the crystal system
        (tetragonal without one) and the value of any component that has to
        be held.

    Returns
    -------
    CellFit
        The fitted cell, the number of peaks used, the root mean square
        difference in degrees between their corrected and recalculated 2theta,
        and the names of the parameters held.

    Raises
    ------
    ValueError
        If ``wavelength`` is None; if there are not more indexed peaks than
        free components; if a component has to be held and no ``start_cell``
        is given; if the peaks cannot separate the free components, naming
        them; or if the fitted reciprocal metric is not positive definite.
    """
    if wavelength is None:
        raise ValueError(
            "refine_cell needs a wavelength; IndexedPeak does not carry one"
        )

    system = start_cell.crystal_system if start_cell else DEFAULT_CRYSTAL_SYSTEM
    components = RECIPROCAL_COMPONENTS[system]
    used = [entry for entry in indexed if entry.is_indexed]
    needed = len(components) + 1
    if len(used) < needed:
        raise ValueError(
            f"Need at least {needed} indexed peaks to refine a {system} cell, "
            f"got {len(used)}"
        )

    positions = np.array([entry.corrected_two_theta for entry in used], dtype=float)
    d_observed = wavelength / (2.0 * np.sin(np.radians(positions / 2.0)))
    inverse_squared = 1.0 / d_observed**2
    hkl = np.array([entry.reflection.hkl for entry in used], dtype=float)

    fitted = []
    held = []
    for component in components:
        _, _, unit = component
        column = np.einsum("ni,ij,nj->n", hkl, np.array(unit, dtype=float), hkl)
        if np.any(column != 0.0):
            fitted.append((component, column))
        else:
            held.append(component)

    held_names = tuple(parameter for _, parameter, _ in held)
    if held and start_cell is None:
        raise ValueError(
            f"No indexed peak carries information on {', '.join(held_names)}, so "
            "it cannot be fitted; pass start_cell to hold it at its value"
        )
    fitted_names = ", ".join(
        f"{name} ({parameter})" for (name, parameter, _), _ in fitted
    )
    design = np.column_stack([column for _, column in fitted]) if fitted else None
    if design is None or np.linalg.matrix_rank(design) < len(fitted):
        raise ValueError(
            "The indexed peaks cannot separate the reciprocal metric components "
            f"{fitted_names or 'none'}; the fit is rank deficient"
        )

    solution, *_ = np.linalg.lstsq(design, inverse_squared, rcond=None)

    reciprocal = np.zeros((3, 3))
    for ((_, _, unit), _), value in zip(fitted, solution):
        reciprocal += float(value) * np.array(unit, dtype=float)
    if held:
        start = start_cell.reciprocal_metric
        for _, _, unit in held:
            matrix = np.array(unit, dtype=float)
            i, j = np.argwhere(matrix)[0]
            reciprocal += start[i, j] / matrix[i, j] * matrix

    if not np.all(np.linalg.eigvalsh(reciprocal) > 0.0):
        raise ValueError(
            "Refinement gave a reciprocal metric that is not positive definite; "
            "the assignments are probably wrong"
        )
    direct = np.linalg.inv(reciprocal)
    a, b, c = np.sqrt(np.diag(direct))
    cosines = (direct[1, 2] / (b * c), direct[0, 2] / (a * c), direct[0, 1] / (a * b))
    alpha, beta, gamma = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    six = dict(
        zip(("a", "b", "c", "alpha", "beta", "gamma"), (a, b, c, alpha, beta, gamma))
    )
    cell = Cell.from_parameters(
        system, {name: float(six[name]) for name in CELL_PARAMETERS[system]}
    )

    recalculated = np.array(
        [
            _two_theta(cell.d_spacing(*entry.reflection.hkl), wavelength)
            for entry in used
        ],
        dtype=float,
    )
    rms = float(np.sqrt(np.mean(np.square(positions - recalculated))))

    return CellFit(cell=cell, n_peaks=len(used), rms_two_theta=rms, held=held_names)


def _seed_cell(
    peaks: list[Peak],
    cell: Cell,
    wavelength: float,
    zero_offset: float,
    space_group: str | None,
) -> Cell:
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
    start_cell: Cell,
    wavelength: float,
    zero_offset: float = DEFAULT_ZERO_OFFSET,
    coarse_tolerance: float = 0.4,
    coarse_two_theta_max: float = 35.0,
    fine_tolerance: float = DEFAULT_TOLERANCE,
    space_group: str | None = None,
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
