"""Tests for xrdkit.indexing."""

import csv

import numpy as np
import pytest

from xrdkit import (
    TTB_CELL,
    IndexedPeak,
    Peak,
    TetragonalCell,
    estimate_zero_offset,
    generate_reflections,
    index_and_refine,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
    refine_cell,
)
from xrdkit.indexing import _is_allowed

WAVELENGTH = 1.5406

ZERO_OFFSET = 0.15

# Reflections further than this from their neighbours are unambiguous targets
# for the synthetic indexing test.
ISOLATION = 0.25


def make_peak(two_theta: float, relative_intensity: float = 50.0) -> Peak:
    """Build a Peak at a given position; only position and intensity matter here."""
    return Peak(
        two_theta=two_theta,
        intensity=1000.0,
        prominence=900.0,
        fwhm=0.15,
        d_spacing=WAVELENGTH / (2.0 * np.sin(np.radians(two_theta / 2.0))),
        relative_intensity=relative_intensity,
    )


def isolated_reflections(
    two_theta_max: float = 40.0,
    cell: TetragonalCell = TTB_CELL,
    separation: float = ISOLATION,
) -> list:
    """Reflections that are well separated from their neighbours."""
    reflections = generate_reflections(cell, WAVELENGTH, two_theta_max)
    positions = [reflection.two_theta for reflection in reflections]
    return [
        reflection
        for index, reflection in enumerate(reflections)
        if (index == 0 or positions[index] - positions[index - 1] > separation)
        and (
            index == len(reflections) - 1
            or positions[index + 1] - positions[index] > separation
        )
    ]


def test_d_spacing_matches_hand_calculation() -> None:
    cell = TetragonalCell(a=12.45, c=3.94)

    # (001) is just c; (110) is a/sqrt(2); (310) is a/sqrt(10).
    assert cell.d_spacing(0, 0, 1) == pytest.approx(3.94)
    assert cell.d_spacing(1, 1, 0) == pytest.approx(12.45 / np.sqrt(2.0))
    assert cell.d_spacing(3, 1, 0) == pytest.approx(12.45 / np.sqrt(10.0))

    # Written out in full, to check the formula and not just its factorisation.
    assert cell.d_spacing(3, 1, 0) == pytest.approx(
        1.0 / np.sqrt((3**2 + 1**2) / 12.45**2 + 0**2 / 3.94**2)
    )
    assert cell.d_spacing(1, 1, 1) == pytest.approx(
        1.0 / np.sqrt((1**2 + 1**2) / 12.45**2 + 1**2 / 3.94**2)
    )


def test_d_spacing_of_000_is_rejected() -> None:
    with pytest.raises(ValueError, match="no d spacing"):
        TTB_CELL.d_spacing(0, 0, 0)


def test_p4bm_conditions_forbid_odd_axial_reflections() -> None:
    # h00 needs h even and 0k0 needs k even, so both (100) and (010) are absent.
    assert not _is_allowed(1, 0, 0, "P4bm")
    assert not _is_allowed(0, 1, 0, "P4bm")
    assert _is_allowed(2, 0, 0, "P4bm")
    assert _is_allowed(0, 2, 0, "P4bm")

    # The same conditions applied to the zones they come from.
    assert not _is_allowed(0, 1, 1, "P4bm")
    assert _is_allowed(0, 2, 1, "P4bm")
    assert not _is_allowed(1, 0, 1, "P4bm")
    assert _is_allowed(2, 0, 1, "P4bm")

    # 00l carries no condition, and general reflections are untouched.
    assert _is_allowed(0, 0, 1, "P4bm")
    assert _is_allowed(3, 1, 0, "P4bm")

    # Without a space group nothing is filtered.
    for indices in ((1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1)):
        assert _is_allowed(*indices, None)


def test_unknown_space_group_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown space group"):
        generate_reflections(TTB_CELL, WAVELENGTH, 40.0, space_group="Pnma")


def test_generate_reflections_excludes_forbidden_and_000() -> None:
    allowed = generate_reflections(TTB_CELL, WAVELENGTH, 40.0)
    unfiltered = generate_reflections(TTB_CELL, WAVELENGTH, 40.0, space_group=None)

    assert (1, 0, 0) not in [reflection.hkl for reflection in allowed]
    assert (1, 0, 0) in [reflection.hkl for reflection in unfiltered]
    # (010) is equivalent to (100) under 4/mmm, so only the h >= k form is ever
    # enumerated, in either mode.
    assert (0, 1, 0) not in [reflection.hkl for reflection in unfiltered]

    assert len(allowed) < len(unfiltered)
    for reflections in (allowed, unfiltered):
        assert all(reflection.hkl != (0, 0, 0) for reflection in reflections)


def test_generate_reflections_is_sorted_and_within_range() -> None:
    low, high = 15.0, 60.0
    reflections = generate_reflections(TTB_CELL, WAVELENGTH, high, low)

    assert reflections
    positions = [reflection.two_theta for reflection in reflections]
    assert positions == sorted(positions)
    assert all(low <= position <= high for position in positions)

    # Only the h >= k >= 0, l >= 0 octant is enumerated.
    assert all(
        reflection.h >= reflection.k >= 0 and reflection.l >= 0
        for reflection in reflections
    )

    # Each position must agree with Bragg's law for its own d spacing.
    for reflection in reflections:
        expected = 2.0 * np.degrees(
            np.arcsin(WAVELENGTH / (2.0 * reflection.d_spacing))
        )
        assert reflection.two_theta == pytest.approx(expected)


def test_generate_reflections_rejects_an_empty_window() -> None:
    with pytest.raises(ValueError, match="two_theta_min < two_theta_max"):
        generate_reflections(TTB_CELL, WAVELENGTH, 20.0, 30.0)


def test_index_peaks_recovers_a_known_zero_offset() -> None:
    reflections = isolated_reflections()
    assert len(reflections) > 8

    rng = np.random.default_rng(20260910)
    shifts = rng.uniform(-0.01, 0.01, len(reflections))
    peaks = [
        make_peak(reflection.two_theta + ZERO_OFFSET + shift)
        for reflection, shift in zip(reflections, shifts)
    ]

    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, zero_offset=ZERO_OFFSET)

    assert len(indexed) == len(peaks)
    assert all(isinstance(entry, IndexedPeak) for entry in indexed)
    # Every peak is indexed, and with the hkl it was built from.
    assert [entry.reflection.hkl for entry in indexed] == [
        reflection.hkl for reflection in reflections
    ]
    # Results stay in the order the peaks were given.
    assert [entry.peak for entry in indexed] == peaks

    for entry, shift in zip(indexed, shifts):
        assert entry.corrected_two_theta == pytest.approx(
            entry.peak.two_theta - ZERO_OFFSET
        )
        assert entry.difference == pytest.approx(shift, abs=1e-9)

    summary = indexing_summary(indexed)
    assert summary["n_peaks"] == len(peaks)
    assert summary["n_indexed"] == len(peaks)
    assert summary["n_unindexed"] == 0
    assert summary["rms_difference"] < 0.01


def test_index_peaks_without_the_offset_leaves_peaks_unindexed() -> None:
    # The same peaks indexed with no zero point correction: a 0.15 degree shift
    # is three times the default tolerance, so the assignment should fail.
    reflections = isolated_reflections()
    peaks = [
        make_peak(reflection.two_theta + ZERO_OFFSET) for reflection in reflections
    ]

    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH)
    summary = indexing_summary(indexed)

    assert summary["n_unindexed"] > 0
    unindexed = [entry for entry in indexed if not entry.is_indexed]
    assert all(entry.reflection is None for entry in unindexed)
    assert all(np.isnan(entry.difference) for entry in unindexed)
    assert all(entry.candidates == [] for entry in unindexed)


def test_index_peaks_records_every_candidate() -> None:
    # A deliberately wide tolerance must expose the ambiguity rather than hide it.
    reflections = generate_reflections(TTB_CELL, WAVELENGTH, 40.0)
    peaks = [make_peak(reflections[0].two_theta)]

    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, tolerance=5.0)
    entry = indexed[0]

    assert len(entry.candidates) > 1
    # The assigned reflection is the closest of the candidates, listed first.
    assert entry.reflection is entry.candidates[0]
    distances = [
        abs(entry.corrected_two_theta - candidate.two_theta)
        for candidate in entry.candidates
    ]
    assert distances == sorted(distances)
    assert indexing_summary(indexed)["n_ambiguous"] == 1


def test_index_peaks_of_nothing() -> None:
    assert index_peaks([], TTB_CELL, WAVELENGTH) == []
    summary = indexing_summary([])
    assert summary["n_peaks"] == 0
    assert np.isnan(summary["rms_difference"])


def test_indexed_to_csv_writes_a_row_per_peak(tmp_path) -> None:
    reflections = isolated_reflections()
    peaks = [
        make_peak(reflection.two_theta + ZERO_OFFSET) for reflection in reflections
    ]
    # One peak far from any reflection, so the unindexed columns are exercised.
    peaks.append(make_peak(41.234))
    peaks.sort(key=lambda peak: peak.two_theta)

    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, zero_offset=ZERO_OFFSET)
    path = indexed_to_csv(indexed, tmp_path / "tables" / "indexed.csv")

    assert path.is_file()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    header, data = rows[0], rows[1:]
    assert header == [
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
    ]
    assert len(data) == len(peaks)
    assert all(len(row[0].split(".")[1]) == 3 for row in data)
    assert all(len(row[2].split(".")[1]) == 4 for row in data)

    # The unindexed peak leaves hkl, position and difference blank.
    blank = [row for row in data if row[4] == ""]
    assert len(blank) == 1
    assert blank[0][5] == "" and blank[0][6] == ""
    assert blank[0][7] == "" and blank[0][8] == ""


# The cell the synthetic refinement peaks come from, about 1% larger in a than
# TTB_CELL so the refinement has something real to recover.
REFINED_CELL = TetragonalCell(a=12.58, c=3.96)

# Largest synthetic position error, in degrees. Well inside the fine tolerance,
# so every peak still indexes once the cell is close.
NOISE = 0.008

# Reflections of REFINED_CELL closer together than this overlap badly enough to
# be genuinely ambiguous, so they make poor synthetic peaks.
REFINEMENT_ISOLATION = 0.30

# The 0.13 angstrom gap between TTB_CELL and REFINED_CELL moves peaks by up to
# 0.3 degrees by 30 degrees 2theta, which the 0.15 degree default cannot bridge.
BRIDGING_TOLERANCE = 0.5


def synthetic_peaks(
    cell: TetragonalCell = REFINED_CELL,
    two_theta_max: float = 45.0,
    seed: int = 20240501,
) -> tuple[list, list[Peak]]:
    """Isolated reflections of ``cell`` and the noisy peaks they would produce."""
    reflections = isolated_reflections(two_theta_max, cell, REFINEMENT_ISOLATION)
    rng = np.random.default_rng(seed)
    peaks = [
        make_peak(round(reflection.two_theta + rng.uniform(-NOISE, NOISE), 3))
        for reflection in reflections
    ]
    return reflections, peaks


def test_index_and_refine_recovers_the_cell_the_peaks_came_from() -> None:
    reflections, peaks = synthetic_peaks()

    indexed, fit = index_and_refine(
        peaks, TTB_CELL, WAVELENGTH, coarse_tolerance=BRIDGING_TOLERANCE
    )

    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert fit.cell.c == pytest.approx(REFINED_CELL.c, abs=0.005)
    assert fit.c_fitted
    assert fit.n_peaks == len(peaks)
    assert fit.rms_two_theta < 0.01

    assert all(entry.is_indexed for entry in indexed)
    assert [entry.reflection.hkl for entry in indexed] == [
        reflection.hkl for reflection in reflections
    ]


def test_index_and_refine_cannot_start_from_too_tight_a_coarse_tolerance() -> None:
    """A tenth of a degree is narrower than the shift TTB_CELL produces.

    With the zero search off, so that this is about the coarse tolerance alone
    and not about an offset standing in for the gap between the cells.
    """
    _, peaks = synthetic_peaks()

    with pytest.raises(ValueError, match="at least 3 indexed peaks"):
        index_and_refine(
            peaks, TTB_CELL, WAVELENGTH, coarse_tolerance=0.15, search_zero=False
        )


def test_index_and_refine_rejects_a_cycle_count_below_one() -> None:
    _, peaks = synthetic_peaks()

    with pytest.raises(ValueError, match="at least one cycle"):
        index_and_refine(peaks, TTB_CELL, WAVELENGTH, n_cycles=0)


def test_refine_cell_needs_three_indexed_peaks() -> None:
    reflections = isolated_reflections()[:2]
    peaks = [make_peak(reflection.two_theta) for reflection in reflections]
    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH)
    assert all(entry.is_indexed for entry in indexed)

    with pytest.raises(ValueError, match="at least 3 indexed peaks"):
        refine_cell(indexed, WAVELENGTH)


def test_refine_cell_needs_a_wavelength() -> None:
    reflections = isolated_reflections()[:5]
    peaks = [make_peak(reflection.two_theta) for reflection in reflections]
    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH)

    with pytest.raises(ValueError, match="needs a wavelength"):
        refine_cell(indexed)


def hk0_indexed() -> list[IndexedPeak]:
    """Indexed peaks of REFINED_CELL that all have l == 0."""
    reflections, peaks = synthetic_peaks()
    hk0 = [
        (reflection, peak)
        for reflection, peak in zip(reflections, peaks)
        if reflection.l == 0
    ]
    indexed = index_peaks([peak for _, peak in hk0], REFINED_CELL, WAVELENGTH)
    assert all(entry.is_indexed and entry.reflection.l == 0 for entry in indexed)
    return indexed


def test_refine_cell_keeps_c_when_no_reflection_has_l() -> None:
    indexed = hk0_indexed()

    fit = refine_cell(indexed, WAVELENGTH, start_cell=TTB_CELL)

    assert not fit.c_fitted
    # hk0 fixes a on its own, and c is carried over untouched.
    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert fit.cell.c == TTB_CELL.c
    assert fit.n_peaks == len(indexed)
    assert fit.rms_two_theta < 0.01


def test_refine_cell_without_l_needs_a_start_cell() -> None:
    indexed = hk0_indexed()

    with pytest.raises(ValueError, match="pass start_cell"):
        refine_cell(indexed, WAVELENGTH)


# A pellet standing proud of its holder shifts every peak by about this much.
APPLIED_ZERO_OFFSET = 0.17


def shifted_peaks(offset: float = APPLIED_ZERO_OFFSET, seed: int = 11) -> tuple:
    """Isolated reflections of REFINED_CELL, every one moved by ``offset``."""
    reflections = isolated_reflections(80.0, REFINED_CELL, REFINEMENT_ISOLATION)
    rng = np.random.default_rng(seed)
    peaks = [
        make_peak(round(reflection.two_theta + offset + rng.uniform(-NOISE, NOISE), 3))
        for reflection in reflections
    ]
    return reflections, peaks


def test_estimate_zero_offset_finds_an_applied_shift() -> None:
    _, peaks = shifted_peaks()

    search = estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH)

    assert search.offset == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.02)
    # It wins on the count, not on a tie-break.
    assert search.n_indexed > int(search.counts[int(np.argmin(np.abs(search.offsets)))])


def test_estimate_zero_offset_returns_near_zero_for_unshifted_peaks() -> None:
    _, peaks = synthetic_peaks()

    search = estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH)

    assert search.offset == pytest.approx(0.0, abs=0.02)


def test_estimate_zero_offset_reports_the_profile_it_tried() -> None:
    _, peaks = shifted_peaks()

    search = estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, step=0.02)

    assert search.offsets.shape == search.counts.shape
    assert search.offsets[0] == pytest.approx(-0.4)
    assert search.offsets[-1] == pytest.approx(0.4)
    assert int(search.counts.max()) == search.n_indexed
    assert search.rms > 0.0


def test_estimate_zero_offset_rejects_a_bad_search_range() -> None:
    _, peaks = shifted_peaks()

    with pytest.raises(ValueError, match="search low < high"):
        estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, search=(0.4, -0.4))
    with pytest.raises(ValueError, match="positive step"):
        estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, step=0.0)
    with pytest.raises(ValueError, match="No peak at or below"):
        estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, two_theta_max=5.0)


def test_index_and_refine_indexes_a_shifted_pattern_by_default() -> None:
    reflections, peaks = shifted_peaks()

    indexed, fit = index_and_refine(peaks, TTB_CELL, WAVELENGTH)

    assert fit.zero_offset == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.02)
    # Every peak indexed, and to the reflection it was made from.
    assert all(entry.is_indexed for entry in indexed)
    assert [entry.reflection.hkl for entry in indexed] == [
        reflection.hkl for reflection in reflections
    ]
    # The cell is recovered too, rather than distorted to absorb the shift.
    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert fit.cell.c == pytest.approx(REFINED_CELL.c, abs=0.005)


def test_the_search_is_skipped_when_a_zero_offset_is_given() -> None:
    _, peaks = shifted_peaks()

    indexed, fit = index_and_refine(
        peaks, TTB_CELL, WAVELENGTH, zero_offset=APPLIED_ZERO_OFFSET
    )

    assert fit.zero_offset == pytest.approx(APPLIED_ZERO_OFFSET)
    assert all(entry.is_indexed for entry in indexed)


def test_the_search_can_be_turned_off() -> None:
    _, peaks = shifted_peaks()

    _, fit = index_and_refine(
        peaks, TTB_CELL, WAVELENGTH, coarse_tolerance=0.5, search_zero=False
    )

    # Without it the shift has to go somewhere, and the cell absorbs it.
    assert fit.zero_offset == 0.0
    assert abs(fit.cell.a - REFINED_CELL.a) > 0.005


def test_an_unshifted_pattern_is_unharmed_by_the_search() -> None:
    reflections, peaks = synthetic_peaks()

    indexed, fit = index_and_refine(peaks, TTB_CELL, WAVELENGTH)

    assert fit.zero_offset == pytest.approx(0.0, abs=0.02)
    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert [entry.reflection.hkl for entry in indexed] == [
        reflection.hkl for reflection in reflections
    ]
