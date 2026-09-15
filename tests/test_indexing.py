"""Tests for xrdkit.indexing."""

import csv

import numpy as np
import pytest

from xrdkit import (
    TTB_CELL,
    Cell,
    IndexedPeak,
    Peak,
    Reflection,
    TetragonalCell,
    estimate_zero_offset,
    generate_reflections,
    index_and_refine,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
    is_absent,
    laue_group,
    multiplicity,
    refine_cell,
    space_group_operations,
)

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
    cell: Cell = TTB_CELL,
    separation: float = ISOLATION,
    space_group: str | None = "P4bm",
) -> list:
    """Reflections that are well separated from their neighbours."""
    reflections = generate_reflections(
        cell, WAVELENGTH, two_theta_max, space_group=space_group
    )
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
    operations = space_group_operations("P4bm")

    def allowed(*hkl):
        return not is_absent(hkl, operations)

    # h00 needs h even and 0k0 needs k even, so both (100) and (010) are absent.
    assert not allowed(1, 0, 0)
    assert not allowed(0, 1, 0)
    assert allowed(2, 0, 0)
    assert allowed(0, 2, 0)

    # The same conditions applied to the zones they come from.
    assert not allowed(0, 1, 1)
    assert allowed(0, 2, 1)
    assert not allowed(1, 0, 1)
    assert allowed(2, 0, 1)

    # 00l carries no condition, and general reflections are untouched.
    assert allowed(0, 0, 1)
    assert allowed(3, 1, 0)

    # Without a space group nothing is filtered: (100) and (101) are generated,
    # and stand for (010) and (011).
    unfiltered = {
        reflection.hkl
        for reflection in generate_reflections(
            TTB_CELL, WAVELENGTH, 40.0, space_group=None
        )
    }
    assert {(1, 0, 0), (1, 0, 1)} <= unfiltered


def test_unknown_space_group_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown space group"):
        generate_reflections(TTB_CELL, WAVELENGTH, 40.0, space_group="Pnma")


def test_generate_reflections_excludes_forbidden_and_000() -> None:
    allowed = generate_reflections(TTB_CELL, WAVELENGTH, 40.0, space_group="P4bm")
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
    reflections = generate_reflections(
        TTB_CELL, WAVELENGTH, high, low, space_group="P4bm"
    )

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
        generate_reflections(TTB_CELL, WAVELENGTH, 20.0, 30.0, space_group="P4bm")


def test_index_peaks_recovers_a_known_zero_offset() -> None:
    reflections = isolated_reflections()
    assert len(reflections) > 8

    rng = np.random.default_rng(20260910)
    shifts = rng.uniform(-0.01, 0.01, len(reflections))
    peaks = [
        make_peak(reflection.two_theta + ZERO_OFFSET + shift)
        for reflection, shift in zip(reflections, shifts)
    ]

    indexed = index_peaks(
        peaks, TTB_CELL, WAVELENGTH, zero_offset=ZERO_OFFSET, space_group="P4bm"
    )

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

    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")
    summary = indexing_summary(indexed)

    assert summary["n_unindexed"] > 0
    unindexed = [entry for entry in indexed if not entry.is_indexed]
    assert all(entry.reflection is None for entry in unindexed)
    assert all(np.isnan(entry.difference) for entry in unindexed)
    assert all(entry.candidates == [] for entry in unindexed)


def test_index_peaks_records_every_candidate() -> None:
    # A deliberately wide tolerance must expose the ambiguity rather than hide it.
    reflections = generate_reflections(TTB_CELL, WAVELENGTH, 40.0, space_group="P4bm")
    peaks = [make_peak(reflections[0].two_theta)]

    indexed = index_peaks(
        peaks, TTB_CELL, WAVELENGTH, tolerance=5.0, space_group="P4bm"
    )
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
    assert index_peaks([], TTB_CELL, WAVELENGTH, space_group="P4bm") == []
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

    indexed = index_peaks(
        peaks, TTB_CELL, WAVELENGTH, zero_offset=ZERO_OFFSET, space_group="P4bm"
    )
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
        peaks,
        TTB_CELL,
        WAVELENGTH,
        coarse_tolerance=BRIDGING_TOLERANCE,
        space_group="P4bm",
    )

    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert fit.cell.c == pytest.approx(REFINED_CELL.c, abs=0.005)
    assert fit.held == ()
    assert fit.n_peaks == len(peaks)
    # The bronze has lone peaks enough below 35 degrees, so the window stays.
    assert fit.coarse_two_theta_max == 35.0
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
            peaks,
            TTB_CELL,
            WAVELENGTH,
            coarse_tolerance=0.15,
            search_zero=False,
            space_group="P4bm",
        )


def test_index_and_refine_rejects_a_cycle_count_below_one() -> None:
    _, peaks = synthetic_peaks()

    with pytest.raises(ValueError, match="at least one cycle"):
        index_and_refine(peaks, TTB_CELL, WAVELENGTH, n_cycles=0, space_group="P4bm")


def test_refine_cell_needs_three_indexed_peaks() -> None:
    reflections = isolated_reflections()[:2]
    peaks = [make_peak(reflection.two_theta) for reflection in reflections]
    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")
    assert all(entry.is_indexed for entry in indexed)

    with pytest.raises(ValueError, match="at least 3 indexed peaks"):
        refine_cell(indexed, WAVELENGTH)


def test_refine_cell_needs_a_wavelength() -> None:
    reflections = isolated_reflections()[:5]
    peaks = [make_peak(reflection.two_theta) for reflection in reflections]
    indexed = index_peaks(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")

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
    indexed = index_peaks(
        [peak for _, peak in hk0], REFINED_CELL, WAVELENGTH, space_group="P4bm"
    )
    assert all(entry.is_indexed and entry.reflection.l == 0 for entry in indexed)
    return indexed


def test_refine_cell_keeps_c_when_no_reflection_has_l() -> None:
    indexed = hk0_indexed()

    fit = refine_cell(indexed, WAVELENGTH, start_cell=TTB_CELL)

    assert fit.held == ("c",)
    assert fit.coarse_two_theta_max is None
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

    search = estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")

    assert search.offset == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.02)
    # It wins on the count, not on a tie-break.
    assert search.n_indexed > int(search.counts[int(np.argmin(np.abs(search.offsets)))])


def test_estimate_zero_offset_returns_near_zero_for_unshifted_peaks() -> None:
    _, peaks = synthetic_peaks()

    search = estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")

    assert search.offset == pytest.approx(0.0, abs=0.02)


def test_estimate_zero_offset_reports_the_profile_it_tried() -> None:
    _, peaks = shifted_peaks()

    search = estimate_zero_offset(
        peaks, TTB_CELL, WAVELENGTH, step=0.02, space_group="P4bm"
    )

    assert search.offsets.shape == search.counts.shape
    assert search.offsets[0] == pytest.approx(-0.4)
    assert search.offsets[-1] == pytest.approx(0.4)
    assert int(search.counts.max()) == search.n_indexed
    assert search.rms > 0.0


def test_estimate_zero_offset_rejects_a_bad_search_range() -> None:
    _, peaks = shifted_peaks()

    with pytest.raises(ValueError, match="search low < high"):
        estimate_zero_offset(
            peaks, TTB_CELL, WAVELENGTH, search=(0.4, -0.4), space_group="P4bm"
        )
    with pytest.raises(ValueError, match="positive step"):
        estimate_zero_offset(peaks, TTB_CELL, WAVELENGTH, step=0.0, space_group="P4bm")
    with pytest.raises(ValueError, match="No peak at or below"):
        estimate_zero_offset(
            peaks, TTB_CELL, WAVELENGTH, two_theta_max=5.0, space_group="P4bm"
        )


def test_index_and_refine_indexes_a_shifted_pattern_by_default() -> None:
    reflections, peaks = shifted_peaks()

    indexed, fit = index_and_refine(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")

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
        peaks, TTB_CELL, WAVELENGTH, zero_offset=APPLIED_ZERO_OFFSET, space_group="P4bm"
    )

    assert fit.zero_offset == pytest.approx(APPLIED_ZERO_OFFSET)
    assert all(entry.is_indexed for entry in indexed)


def test_the_search_can_be_turned_off() -> None:
    _, peaks = shifted_peaks()

    _, fit = index_and_refine(
        peaks,
        TTB_CELL,
        WAVELENGTH,
        coarse_tolerance=0.5,
        search_zero=False,
        space_group="P4bm",
    )

    # Without it the shift has to go somewhere, and the cell absorbs it.
    assert fit.zero_offset == 0.0
    assert abs(fit.cell.a - REFINED_CELL.a) > 0.005


def test_an_unshifted_pattern_is_unharmed_by_the_search() -> None:
    reflections, peaks = synthetic_peaks()

    indexed, fit = index_and_refine(peaks, TTB_CELL, WAVELENGTH, space_group="P4bm")

    assert fit.zero_offset == pytest.approx(0.0, abs=0.02)
    assert fit.cell.a == pytest.approx(REFINED_CELL.a, abs=0.005)
    assert [entry.reflection.hkl for entry in indexed] == [
        reflection.hkl for reflection in reflections
    ]


# Cells of every crystal system, the space group to generate them in, and a
# start cell a little way off, for the general refinement tests.
SYSTEM_CASES = {
    "cubic": (Cell.cubic(3.905), "Pm-3m", Cell.cubic(3.95)),
    "tetragonal": (Cell.tetragonal(3.994, 4.034), "P4mm", Cell.tetragonal(4.0, 4.0)),
    "orthorhombic": (
        Cell.orthorhombic(5.38, 5.44, 7.64),
        "Pbnm",
        Cell.orthorhombic(5.4, 5.4, 7.6),
    ),
    "trigonal": (Cell.trigonal(5.148, 13.863), "R3c", Cell.trigonal(5.2, 13.9)),
    "monoclinic": (
        Cell.monoclinic(5.0, 6.0, 7.0, 100.0),
        None,
        Cell.monoclinic(5.05, 5.95, 7.1, 99.0),
    ),
    "triclinic": (
        Cell.triclinic(5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
        None,
        Cell.triclinic(5.1, 6.1, 6.9, 81.0, 91.0, 99.0),
    ),
}


def exact_indexed(reflections) -> list[IndexedPeak]:
    """Peaks sitting exactly on ``reflections``, each assigned to its own."""
    return [
        IndexedPeak(
            peak=make_peak(reflection.two_theta),
            corrected_two_theta=reflection.two_theta,
            reflection=reflection,
            difference=0.0,
            candidates=[reflection],
        )
        for reflection in reflections
    ]


def test_cubic_generation_merges_equivalents_with_multiplicity() -> None:
    reflections = generate_reflections(
        Cell.cubic(3.905), WAVELENGTH, 90.0, space_group="Pm-3m"
    )
    hkl = [reflection.hkl for reflection in reflections]

    assert hkl.count((1, 0, 0)) == 1
    assert (0, 0, 1) not in hkl and (0, 1, 0) not in hkl
    by_hkl = {reflection.hkl: reflection for reflection in reflections}
    assert by_hkl[(1, 0, 0)].multiplicity == 6
    assert by_hkl[(1, 1, 0)].multiplicity == 12
    assert by_hkl[(1, 1, 1)].multiplicity == 8
    assert len(set(hkl)) == len(hkl)


def test_generation_without_a_space_group_still_merges_equivalents() -> None:
    cubic = generate_reflections(Cell.cubic(3.905), WAVELENGTH, 90.0)
    hkl = [reflection.hkl for reflection in cubic]
    assert hkl.count((1, 0, 0)) == 1
    assert (0, 0, 1) not in hkl
    assert hkl == [
        reflection.hkl
        for reflection in generate_reflections(
            Cell.cubic(3.905), WAVELENGTH, 90.0, space_group="Pm-3m"
        )
    ]

    # No absences: the (100) P4bm forbids is there, merged with (010).
    tetragonal = {
        reflection.hkl: reflection.multiplicity
        for reflection in generate_reflections(TTB_CELL, WAVELENGTH, 30.0)
    }
    assert tetragonal[(1, 0, 0)] == 4
    assert tetragonal[(1, 1, 0)] == 4
    assert tetragonal[(2, 1, 0)] == 8
    assert tetragonal[(0, 0, 1)] == 2
    assert (0, 1, 0) not in tetragonal


def test_multiplicity_and_d_spacing_follow_the_cell_and_laue_group() -> None:
    laue = laue_group(space_group_operations("P4bm"))
    for reflection in generate_reflections(
        TTB_CELL, WAVELENGTH, 60.0, space_group="P4bm"
    ):
        assert reflection.multiplicity == multiplicity(reflection.hkl, laue)
        assert reflection.d_spacing == pytest.approx(
            TTB_CELL.d_spacing(*reflection.hkl), rel=1e-12
        )


def test_r3c_generation_leaves_out_every_absent_reflection() -> None:
    cell, space_group, _ = SYSTEM_CASES["trigonal"]
    reflections = generate_reflections(cell, WAVELENGTH, 90.0, space_group=space_group)
    operations = space_group_operations("R3c")

    assert len(reflections) > 20
    assert not any(is_absent(reflection.hkl, operations) for reflection in reflections)
    hkl = {reflection.hkl for reflection in reflections}
    assert {(0, 1, 2), (1, 0, 4), (1, 1, 0), (0, 0, 6)} <= hkl
    assert not {(0, 0, 3), (1, 0, 1), (0, 1, 4)} & hkl


def test_monoclinic_generation_keeps_101_and_10_1_apart() -> None:
    cell, _, _ = SYSTEM_CASES["monoclinic"]
    by_hkl = {
        reflection.hkl: reflection
        for reflection in generate_reflections(cell, WAVELENGTH, 60.0)
    }

    assert (1, 0, 1) in by_hkl and (1, 0, -1) in by_hkl
    assert by_hkl[(1, 0, 1)].two_theta != pytest.approx(by_hkl[(1, 0, -1)].two_theta)


@pytest.mark.parametrize("system", list(SYSTEM_CASES))
def test_refine_cell_recovers_every_free_parameter(system) -> None:
    cell, space_group, start = SYSTEM_CASES[system]
    reflections = generate_reflections(cell, WAVELENGTH, 90.0, space_group=space_group)

    fit = refine_cell(exact_indexed(reflections), WAVELENGTH, start_cell=start)

    assert fit.cell.crystal_system == cell.crystal_system
    assert fit.held == ()
    assert fit.n_peaks == len(reflections)
    assert fit.rms_two_theta < 1e-6
    assert tuple(fit.cell.parameters) == cell.parameter_names
    for name, value in cell.parameters.items():
        assert fit.cell.parameters[name] == pytest.approx(value, abs=1e-6), name


def test_index_and_refine_recovers_a_cubic_cell_from_a_distant_start() -> None:
    cell = Cell.cubic(3.905)
    reflections = generate_reflections(cell, WAVELENGTH, 90.0, space_group="Pm-3m")
    peaks = [make_peak(reflection.two_theta) for reflection in reflections]

    indexed, fit = index_and_refine(
        peaks,
        Cell.cubic(3.95),
        WAVELENGTH,
        coarse_tolerance=0.5,
        space_group="Pm-3m",
        search_zero=False,
    )

    assert fit.cell.crystal_system == "cubic"
    assert fit.cell.a == pytest.approx(3.905, abs=1e-6)
    # Every peak lands on its own reflection; (221) and (300), which share a d
    # spacing, are one peak with two candidates.
    assert all(entry.is_indexed for entry in indexed)
    for entry, reflection in zip(indexed, reflections):
        assert entry.reflection.two_theta == pytest.approx(reflection.two_theta)


@pytest.mark.parametrize(("system", "needed"), [("cubic", 2), ("triclinic", 7)])
def test_refine_cell_needs_one_more_peak_than_free_components(system, needed) -> None:
    cell, space_group, start = SYSTEM_CASES[system]
    reflections = generate_reflections(cell, WAVELENGTH, 90.0, space_group=space_group)

    with pytest.raises(ValueError, match=f"at least {needed} indexed peaks"):
        refine_cell(
            exact_indexed(reflections[: needed - 1]), WAVELENGTH, start_cell=start
        )


@pytest.mark.parametrize(
    ("system", "held"), [("trigonal", ("c",)), ("monoclinic", ("c", "beta"))]
)
def test_refine_cell_holds_what_no_peak_carries(system, held) -> None:
    cell, space_group, start = SYSTEM_CASES[system]
    reflections = [
        reflection
        for reflection in generate_reflections(
            cell, WAVELENGTH, 90.0, space_group=space_group
        )
        if reflection.l == 0
    ]

    fit = refine_cell(exact_indexed(reflections), WAVELENGTH, start_cell=start)

    assert fit.held == held
    if system == "trigonal":
        # On hexagonal axes the hk0 reflections fix a alone, and c is kept.
        assert fit.cell.a == pytest.approx(cell.a, abs=1e-6)
        assert fit.cell.c == pytest.approx(start.c, abs=1e-9)


def test_refine_cell_rank_deficient_names_the_components() -> None:
    cell = Cell.tetragonal(4.0, 4.0)
    # h0h reflections move h^2 + k^2 and l^2 together, so a and c cannot part.
    reflections = []
    for h in (1, 2, 3):
        d = cell.d_spacing(h, 0, h)
        two_theta = 2.0 * float(np.degrees(np.arcsin(WAVELENGTH / (2.0 * d))))
        reflections.append(Reflection(h, 0, h, d, two_theta))

    with pytest.raises(ValueError, match=r"A \(a\), C \(c\)"):
        refine_cell(exact_indexed(reflections), WAVELENGTH, start_cell=cell)


# The pseudo-cubic Pbnm perovskite and the start it is indexed from, a and b
# split the wrong way round from it by about a per cent.
PBNM_CELL = Cell.orthorhombic(5.38, 5.44, 7.64)
PBNM_START = Cell.orthorhombic(5.30, 5.50, 7.60)


def lone_peaks(cell: Cell, space_group: str) -> list[Peak]:
    """Peaks on the reflections of ``cell`` to 80 degrees that stand clear of
    their neighbours, as a scan would resolve them."""
    reflections = isolated_reflections(80.0, cell, REFINEMENT_ISOLATION, space_group)
    return [make_peak(round(reflection.two_theta, 3)) for reflection in reflections]


def test_the_coarse_window_widens_for_a_pseudo_cubic_perovskite() -> None:
    peaks = lone_peaks(PBNM_CELL, "Pbnm")

    indexed, fit = index_and_refine(peaks, PBNM_START, WAVELENGTH, space_group="Pbnm")

    # Below 35 degrees it has three lone peaks; an orthorhombic cell needs five.
    assert fit.coarse_two_theta_max > 35.0
    assert fit.coarse_two_theta_max == 50.0
    for name, value in PBNM_CELL.parameters.items():
        assert fit.cell.parameters[name] == pytest.approx(value, abs=1e-3), name
    assert all(entry.is_indexed for entry in indexed)


def test_a_fixed_35_degree_window_is_too_narrow_for_the_pbnm_perovskite() -> None:
    peaks = lone_peaks(PBNM_CELL, "Pbnm")

    with pytest.raises(ValueError, match="at least 4 indexed peaks"):
        index_and_refine(
            peaks,
            PBNM_START,
            WAVELENGTH,
            space_group="Pbnm",
            coarse_two_theta_max=35.0,
        )


def test_a_cubic_cell_with_few_low_angle_peaks_is_indexed() -> None:
    peaks = lone_peaks(Cell.cubic(3.905), "Pm-3m")

    _, fit = index_and_refine(peaks, Cell.cubic(3.95), WAVELENGTH, space_group="Pm-3m")

    assert fit.cell.a == pytest.approx(3.905, abs=1e-3)
    # From 3.95 only (100) and (110) fall inside the coarse tolerance at any
    # window below the top of the pattern, so the whole pattern is taken.
    assert fit.coarse_two_theta_max == pytest.approx(
        max(peak.two_theta for peak in peaks)
    )


def test_too_few_peaks_in_the_whole_pattern_still_raises() -> None:
    peaks = lone_peaks(Cell.cubic(3.905), "Pm-3m")[:1]

    with pytest.raises(ValueError, match="at least 2 indexed peaks"):
        index_and_refine(peaks, Cell.cubic(3.95), WAVELENGTH, space_group="Pm-3m")
