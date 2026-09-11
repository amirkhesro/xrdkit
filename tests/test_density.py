"""Tests for xrdkit.density."""

import numpy as np
import pytest

from xrdkit import (
    ATOMIC_MASSES,
    TetragonalCell,
    cell_volume,
    formula_mass,
    theoretical_density,
)

SRTIO3 = {"Sr": 1, "Ti": 1, "O": 3}

# 87.62 + 47.867 + 3 x 15.999, worked by hand.
SRTIO3_MASS = 183.484

# Cubic SrTiO3, a = 3.905 A, one formula unit per cell.
SRTIO3_A = 3.905


def test_formula_mass_srtio3_matches_hand_value() -> None:
    assert formula_mass(SRTIO3) == pytest.approx(SRTIO3_MASS, abs=1e-9)


def test_formula_mass_takes_fractional_coefficients() -> None:
    sbn = {"Sr": 0.5, "Ba": 0.5, "Nb": 2, "O": 6}
    expected = 0.5 * 87.62 + 0.5 * 137.327 + 2 * 92.90637 + 6 * 15.999

    assert formula_mass(sbn) == pytest.approx(expected)


@pytest.mark.parametrize(
    "composition",
    [{}, {"Xx": 1}, {"Sr": 1, "O": -1}],
    ids=["empty", "unknown element", "negative"],
)
def test_formula_mass_rejects_bad_compositions(composition: dict) -> None:
    with pytest.raises(ValueError):
        formula_mass(composition)


def test_atomic_masses_cover_the_elements_asked_for() -> None:
    wanted = [
        "Sr",
        "Ba",
        "La",
        "Nb",
        "Ti",
        "O",
        "Ca",
        "Na",
        "K",
        "Bi",
        "Pb",
        "Zr",
        "Fe",
        "Mn",
        "Mg",
        "Al",
        "Si",
        "Li",
        "Ta",
        "W",
        "Sn",
        "Ce",
        "Nd",
        "Sm",
        "Gd",
        "Y",
        "Zn",
        "Cu",
        "Ni",
        "Co",
        "Cr",
        "Hf",
    ]

    assert set(wanted) <= set(ATOMIC_MASSES)


def test_density_of_cubic_srtio3() -> None:
    volume, _ = cell_volume(TetragonalCell(a=SRTIO3_A, c=SRTIO3_A))
    density, esd = theoretical_density(SRTIO3, 1, volume)

    assert density == pytest.approx(5.11, abs=0.01)
    assert esd is None


def test_cell_volume_without_esds() -> None:
    volume, esd = cell_volume(TetragonalCell(a=12.5, c=3.9))

    assert volume == pytest.approx(12.5**2 * 3.9)
    assert esd is None


@pytest.mark.parametrize(
    ("esd_a", "esd_c", "expected"),
    [
        (0.001, None, 2 * 12.5 * 3.9 * 0.001),
        (None, 0.002, 12.5**2 * 0.002),
        (0.001, 0.002, np.hypot(2 * 12.5 * 3.9 * 0.001, 12.5**2 * 0.002)),
    ],
    ids=["a only", "c only", "both"],
)
def test_cell_volume_esd_propagation(
    esd_a: float | None, esd_c: float | None, expected: float
) -> None:
    _, esd = cell_volume(TetragonalCell(a=12.5, c=3.9), esd_a=esd_a, esd_c=esd_c)

    assert esd == pytest.approx(expected)


def test_cell_volume_esd_agrees_with_monte_carlo() -> None:
    cell = TetragonalCell(a=12.47, c=3.93)
    esd_a, esd_c = 0.002, 0.001
    _, esd = cell_volume(cell, esd_a=esd_a, esd_c=esd_c)

    rng = np.random.default_rng(1)
    a = rng.normal(cell.a, esd_a, 200_000)
    c = rng.normal(cell.c, esd_c, 200_000)

    assert np.std(a**2 * c) == pytest.approx(esd, rel=0.01)


def test_density_esd_has_the_relative_esd_of_the_volume() -> None:
    volume, esd_volume = 600.0, 0.6
    density, esd = theoretical_density(SRTIO3, 5, volume, esd_volume)

    assert esd / density == pytest.approx(esd_volume / volume)


@pytest.mark.parametrize(("z", "volume"), [(0, 60.0), (1, 0.0), (1, -60.0)])
def test_density_rejects_non_positive_inputs(z: float, volume: float) -> None:
    with pytest.raises(ValueError):
        theoretical_density(SRTIO3, z, volume)
