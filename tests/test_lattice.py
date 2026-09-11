"""Tests for xrdkit.lattice."""

import numpy as np
import pytest

from xrdkit import (
    TTB_CELL,
    LatticeFit,
    Peak,
    TetragonalCell,
    generate_reflections,
    index_and_refine,
    lattice_fit_to_dict,
    refine_lattice,
)

WAVELENGTH = 1.5406

# The cell the synthetic peaks come from, and the systematic errors applied to
# them before they are handed to the indexer.
TRUE_CELL = TetragonalCell(a=12.48, c=3.93)
ZERO_OFFSET = 0.17
DISPLACEMENT = 0.1
RADIUS_MM = 145.0

# Largest random error added to a synthetic position, in degrees.
NOISE = 0.005

# Reflections closer together than this overlap badly enough to be ambiguous.
ISOLATION = 0.40

TWO_THETA_MAX = 80.0

# A zero point error of 0.17 degrees is more than the 0.05 degree default lets
# a peak stray, so the indexing that feeds the refinement is run wider. The
# cell it settles on is wrong either way, which is the point: the refinement
# has to take the offset out again.
FINE_TOLERANCE = 0.15

# What the fit has to get back to, in angstroms, degrees and millimetres.
CELL_TOLERANCE = 0.003
ZERO_TOLERANCE = 0.01
DISPLACEMENT_TOLERANCE = 0.02


def make_peak(two_theta: float) -> Peak:
    """A Peak at a position; only the position matters to the refinement."""
    return Peak(
        two_theta=two_theta,
        intensity=1000.0,
        prominence=900.0,
        fwhm=0.15,
        d_spacing=WAVELENGTH / (2.0 * np.sin(np.radians(two_theta / 2.0))),
        relative_intensity=50.0,
    )


def synthetic_peaks(
    zero: float = ZERO_OFFSET,
    displacement: float = 0.0,
    seed: int = 7,
) -> list[Peak]:
    """Isolated reflections of TRUE_CELL, displaced and offset as a scan would be."""
    reflections = generate_reflections(TRUE_CELL, WAVELENGTH, TWO_THETA_MAX)
    positions = [reflection.two_theta for reflection in reflections]
    isolated = [
        reflection
        for index, reflection in enumerate(reflections)
        if (index == 0 or positions[index] - positions[index - 1] > ISOLATION)
        and (
            index == len(reflections) - 1
            or positions[index + 1] - positions[index] > ISOLATION
        )
    ]

    rng = np.random.default_rng(seed)
    peaks = []
    for reflection in isolated:
        theta = np.radians(reflection.two_theta / 2.0)
        shift = zero + np.degrees(-2.0 * displacement * np.cos(theta) / RADIUS_MM)
        peaks.append(
            make_peak(
                round(reflection.two_theta + shift + rng.uniform(-NOISE, NOISE), 3)
            )
        )
    return peaks


def indexed_synthetic(
    zero: float = ZERO_OFFSET, displacement: float = 0.0, seed: int = 7
):
    """Index synthetic peaks without correcting the offset, so the cell is wrong.

    The automatic search is off on purpose: it would find the offset itself and
    hand refine_lattice a cell that is already right, leaving these tests with
    nothing to recover. Finding it is index_and_refine's job and is tested
    there; taking it back out of a distorted cell is what is tested here.
    """
    peaks = synthetic_peaks(zero=zero, displacement=displacement, seed=seed)
    return index_and_refine(
        peaks,
        TTB_CELL,
        WAVELENGTH,
        zero_offset=0.0,
        fine_tolerance=FINE_TOLERANCE,
        search_zero=False,
    )


def test_indexing_without_the_offset_distorts_the_cell() -> None:
    """The premise of the refinement tests: the starting cell really is wrong."""
    _, cell_fit = indexed_synthetic()

    # A constant offset cannot be absorbed by a cell, so the linear fit splits
    # the difference and lands well away from the truth.
    assert abs(cell_fit.cell.a - TRUE_CELL.a) > CELL_TOLERANCE


def test_refine_lattice_recovers_the_cell_and_the_zero_offset() -> None:
    indexed, cell_fit = indexed_synthetic()

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    assert fit.converged
    assert fit.a == pytest.approx(TRUE_CELL.a, abs=CELL_TOLERANCE)
    assert fit.c == pytest.approx(TRUE_CELL.c, abs=CELL_TOLERANCE)
    assert fit.zero == pytest.approx(ZERO_OFFSET, abs=ZERO_TOLERANCE)
    # The cell and the offset are separated, not traded off against each other.
    assert abs(fit.a - TRUE_CELL.a) < abs(cell_fit.cell.a - TRUE_CELL.a)


def test_the_esds_are_positive_and_well_inside_the_tolerances() -> None:
    indexed, cell_fit = indexed_synthetic()

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    assert fit.esd_a > 0.0
    assert fit.esd_c > 0.0
    assert fit.esd_zero > 0.0
    assert fit.esd_a < CELL_TOLERANCE
    assert fit.esd_c < CELL_TOLERANCE
    assert fit.esd_zero < ZERO_TOLERANCE
    # An honest error bar covers the error actually made, within a factor of a
    # few; one far smaller than that would be claiming too much.
    assert abs(fit.a - TRUE_CELL.a) < 5.0 * fit.esd_a
    assert abs(fit.zero - ZERO_OFFSET) < 5.0 * fit.esd_zero


def test_refine_lattice_recovers_a_specimen_displacement() -> None:
    indexed, cell_fit = indexed_synthetic(zero=0.0, displacement=DISPLACEMENT)

    # No zero point error was applied, and it is left fixed: over this range a
    # constant and a cos(theta) term are too alike to refine together.
    fit = refine_lattice(
        indexed,
        WAVELENGTH,
        cell_fit.cell,
        fit_zero=False,
        fit_displacement=True,
        radius_mm=RADIUS_MM,
    )

    assert fit.converged
    assert fit.displacement == pytest.approx(DISPLACEMENT, abs=DISPLACEMENT_TOLERANCE)
    assert fit.esd_displacement > 0.0
    assert fit.a == pytest.approx(TRUE_CELL.a, abs=CELL_TOLERANCE)
    assert fit.radius_mm == RADIUS_MM


def test_a_zero_and_a_displacement_together_are_poorly_determined() -> None:
    """Both at once is allowed, but the error bar says how little it is worth."""
    indexed, cell_fit = indexed_synthetic(zero=0.0, displacement=DISPLACEMENT)

    alone = refine_lattice(
        indexed,
        WAVELENGTH,
        cell_fit.cell,
        fit_zero=False,
        fit_displacement=True,
        radius_mm=RADIUS_MM,
    )
    together = refine_lattice(
        indexed,
        WAVELENGTH,
        cell_fit.cell,
        fit_zero=True,
        fit_displacement=True,
        radius_mm=RADIUS_MM,
    )

    # Refining the pair inflates the displacement's error bar by an order of
    # magnitude, which is the correlation between them showing up honestly.
    assert together.esd_displacement > 10.0 * alone.esd_displacement


def test_fit_displacement_without_a_radius_is_rejected() -> None:
    indexed, cell_fit = indexed_synthetic()

    with pytest.raises(ValueError, match="radius_mm"):
        refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_displacement=True)


def test_a_fixed_parameter_reports_no_esd() -> None:
    indexed, cell_fit = indexed_synthetic()

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=False)

    assert fit.zero == 0.0
    assert fit.esd_zero is None
    assert fit.displacement is None
    assert fit.esd_displacement is None
    assert fit.radius_mm is None


def test_refine_lattice_needs_more_peaks_than_parameters() -> None:
    indexed, cell_fit = indexed_synthetic()
    single = [
        entry for entry in indexed if entry.is_indexed and len(entry.candidates) == 1
    ]

    # Three parameters, so six peaks are the fewest that mean anything.
    with pytest.raises(ValueError, match="at least 6 peaks"):
        refine_lattice(single[:5], WAVELENGTH, cell_fit.cell, fit_zero=True)


def test_only_unambiguous_peaks_are_used() -> None:
    indexed, cell_fit = indexed_synthetic()
    single = sum(
        1 for entry in indexed if entry.is_indexed and len(entry.candidates) == 1
    )

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    assert fit.n_peaks == single
    assert len(fit.hkl) == single
    assert fit.residuals.shape == (single,)
    # rms is the root mean square of the residuals it reports.
    assert fit.rms_two_theta == pytest.approx(float(np.sqrt(np.mean(fit.residuals**2))))


def test_the_residuals_are_small_once_the_offset_is_taken_out() -> None:
    indexed, cell_fit = indexed_synthetic()

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    # Only the noise should be left, which was drawn inside NOISE.
    assert fit.rms_two_theta < NOISE
    assert np.max(np.abs(fit.residuals)) < 3.0 * NOISE


def test_the_fitted_cell_is_available_as_a_cell() -> None:
    indexed, cell_fit = indexed_synthetic()

    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    assert isinstance(fit.cell, TetragonalCell)
    assert fit.cell.a == fit.a
    assert fit.cell.c == fit.c


def test_lattice_fit_to_dict_is_flat_and_scalar() -> None:
    indexed, cell_fit = indexed_synthetic()
    fit = refine_lattice(indexed, WAVELENGTH, cell_fit.cell, fit_zero=True)

    row = lattice_fit_to_dict(fit)

    assert isinstance(fit, LatticeFit)
    assert set(row) == {
        "a",
        "esd_a",
        "c",
        "esd_c",
        "c_over_a",
        "zero",
        "esd_zero",
        "displacement",
        "esd_displacement",
        "radius_mm",
        "n_peaks",
        "rms_two_theta",
        "converged",
    }
    # Every value is something a CSV writer can take as it stands.
    assert all(
        value is None or isinstance(value, (int, float, bool)) for value in row.values()
    )
    assert row["a"] == fit.a
    assert row["c_over_a"] == pytest.approx(fit.c / fit.a)
