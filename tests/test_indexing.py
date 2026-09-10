"""Tests for xrdkit.indexing."""

import csv

import numpy as np
import pytest

from xrdkit import (
    TTB_CELL,
    IndexedPeak,
    Peak,
    TetragonalCell,
    generate_reflections,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
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


def isolated_reflections(two_theta_max: float = 40.0) -> list:
    """Reflections that are well separated from their neighbours."""
    reflections = generate_reflections(TTB_CELL, WAVELENGTH, two_theta_max)
    positions = [reflection.two_theta for reflection in reflections]
    return [
        reflection
        for index, reflection in enumerate(reflections)
        if (index == 0 or positions[index] - positions[index - 1] > ISOLATION)
        and (
            index == len(reflections) - 1
            or positions[index + 1] - positions[index] > ISOLATION
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
