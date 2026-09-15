"""Tests for xrdkit.density."""

import itertools
import re

import numpy as np
import pytest

from xrdkit import (
    ATOMIC_MASSES,
    Cell,
    cell_volume,
    formula_mass,
    parse_formula,
    relative_density,
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
    volume, _ = cell_volume(Cell.cubic(SRTIO3_A))
    density, esd = theoretical_density(SRTIO3, 1, volume)

    assert density == pytest.approx(5.11, abs=0.01)
    assert esd is None


def test_cell_volume_without_esds() -> None:
    volume, esd = cell_volume(Cell.tetragonal(12.5, 3.9))

    assert volume == pytest.approx(12.5**2 * 3.9)
    assert esd is None


@pytest.mark.parametrize(
    ("esd", "expected"),
    [
        ({"a": 0.001}, 2 * 12.5 * 3.9 * 0.001),
        ({"c": 0.002}, 12.5**2 * 0.002),
        ({"a": 0.001, "c": 0.002}, np.hypot(2 * 12.5 * 3.9 * 0.001, 12.5**2 * 0.002)),
    ],
    ids=["a only", "c only", "both"],
)
def test_cell_volume_esd_propagation(esd: dict[str, float], expected: float) -> None:
    # A free parameter left out of the mapping counts as zero.
    _, propagated = cell_volume(Cell.tetragonal(12.5, 3.9), esd)

    assert propagated == pytest.approx(expected, rel=1e-6)


def test_cell_volume_esd_agrees_with_monte_carlo() -> None:
    cell = Cell.tetragonal(12.47, 3.93)
    esd_a, esd_c = 0.002, 0.001
    _, esd = cell_volume(cell, {"a": esd_a, "c": esd_c})

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


# parse_formula


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("LaB6", {"La": 1.0, "B": 6.0}),
        (
            "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6",
            {"Sr": 0.4, "Ba": 0.5, "La": 0.1, "Nb": 1.9, "Ti": 0.1, "O": 6.0},
        ),
        ("Ca(OH)2", {"Ca": 1.0, "O": 2.0, "H": 2.0}),
        ("Mg3(PO4)2", {"Mg": 3.0, "P": 2.0, "O": 8.0}),
        ("K4(Fe(CN)6)", {"K": 4.0, "Fe": 1.0, "C": 6.0, "N": 6.0}),
        ("CH3COOH", {"C": 2.0, "H": 4.0, "O": 2.0}),
        ("CO", {"C": 1.0, "O": 1.0}),
        (" SrTiO3 ", {"Sr": 1.0, "Ti": 1.0, "O": 3.0}),
    ],
    ids=[
        "LaB6",
        "TTB",
        "parentheses",
        "group of groups",
        "nested",
        "repeat",
        "CO",
        "spaces",
    ],
)
def test_parse_formula(formula: str, expected: dict) -> None:
    parsed = parse_formula(formula)

    assert parsed == pytest.approx(expected)
    assert list(parsed) == list(expected)


@pytest.mark.parametrize(
    ("formula", "offending"),
    [
        ("", "formula is empty"),
        ("Sr0.4 Ba0.5", "' Ba0.5'"),
        ("srTiO3", "'srTiO3'"),
        ("Xx2O3", "'Xx'"),
        ("Ca(OH2", "'(OH2'"),
        ("CaOH)2", "')2'"),
        ("Ca()2", "'()'"),
        ("O6.", "'.'"),
        ("Ti-O2", "'-O2'"),
    ],
    ids=[
        "empty",
        "space",
        "lower case",
        "unknown",
        "open",
        "close",
        "empty group",
        "stray dot",
        "stray sign",
    ],
)
def test_parse_formula_rejects_malformed_text(formula: str, offending: str) -> None:
    with pytest.raises(ValueError, match=re.escape(offending)):
        parse_formula(formula)


def test_formula_mass_takes_a_formula_string() -> None:
    assert formula_mass("SrTiO3") == pytest.approx(SRTIO3_MASS, abs=1e-9)
    assert formula_mass("Ca(OH)2") == pytest.approx(40.078 + 2 * (15.999 + 1.008))


def test_theoretical_density_takes_a_formula_string() -> None:
    volume, esd_volume = 610.0, 0.3
    composition = {"Sr": 0.4, "Ba": 0.5, "La": 0.1, "Nb": 1.9, "Ti": 0.1, "O": 6}

    from_string = theoretical_density(
        "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6", 5, volume, esd_volume
    )

    assert from_string == pytest.approx(
        theoretical_density(composition, 5, volume, esd_volume)
    )


def test_atomic_masses_run_from_hydrogen_to_uranium() -> None:
    # fmt: off
    symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
        "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
        "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",
        "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    ]
    # fmt: on

    assert len(symbols) == 92
    assert list(ATOMIC_MASSES) == symbols
    assert all(isinstance(mass, float) for mass in ATOMIC_MASSES.values())
    # Masses rise with atomic number except at the known inversions.
    inversions = {("Ar", "K"), ("Co", "Ni"), ("Te", "I"), ("Th", "Pa")}
    for lighter, heavier in itertools.pairwise(symbols):
        if (lighter, heavier) not in inversions:
            assert ATOMIC_MASSES[lighter] < ATOMIC_MASSES[heavier], (lighter, heavier)


def test_atomic_masses_keep_their_earlier_values() -> None:
    assert ATOMIC_MASSES["Na"] == 22.98976928
    assert ATOMIC_MASSES["Gd"] == 157.25
    assert ATOMIC_MASSES["Zr"] == 91.224
    assert ATOMIC_MASSES["Bi"] == 208.98040
    assert ATOMIC_MASSES["O"] == 15.999


# relative_density


def test_relative_density_without_esds() -> None:
    ratio, esd = relative_density(5.15, 5.40)

    assert ratio == pytest.approx(100 * 5.15 / 5.40)
    assert esd is None


@pytest.mark.parametrize(
    ("esd_measured", "esd_theoretical", "relative"),
    [
        (0.02, None, 0.02 / 5.15),
        (None, 0.004, 0.004 / 5.40),
        (0.02, 0.004, np.hypot(0.02 / 5.15, 0.004 / 5.40)),
    ],
    ids=["measured only", "theoretical only", "both"],
)
def test_relative_density_adds_relative_errors_in_quadrature(
    esd_measured: float | None, esd_theoretical: float | None, relative: float
) -> None:
    ratio, esd = relative_density(5.15, 5.40, esd_measured, esd_theoretical)

    assert esd == pytest.approx(ratio * relative)


@pytest.mark.parametrize(
    ("measured", "theoretical"), [(0.0, 5.4), (5.1, 0.0), (-1, 5.4)]
)
def test_relative_density_rejects_non_positive_densities(
    measured: float, theoretical: float
) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        relative_density(measured, theoretical)


def test_cell_volume_esd_for_a_hexagonal_cell() -> None:
    a, c, esd_a, esd_c = 5.148, 13.863, 0.0004, 0.002
    volume, esd = cell_volume(Cell.hexagonal(a, c), {"a": esd_a, "c": esd_c})

    # V = (sqrt 3 / 2) a^2 c, so dV/da = sqrt 3 a c and dV/dc = (sqrt 3 / 2) a^2.
    assert volume == pytest.approx(np.sqrt(3) / 2 * a**2 * c)
    expected = np.hypot(np.sqrt(3) * a * c * esd_a, np.sqrt(3) / 2 * a**2 * esd_c)
    assert esd == pytest.approx(expected, rel=1e-6)
    _, c_only = cell_volume(Cell.trigonal(a, c), {"c": esd_c})
    assert c_only == pytest.approx(np.sqrt(3) / 2 * a**2 * esd_c, rel=1e-6)


def test_cell_volume_esd_with_nothing_to_propagate_is_none() -> None:
    cell = Cell.tetragonal(12.5, 3.9)

    assert cell_volume(cell)[1] is None
    assert cell_volume(cell, {})[1] is None
    assert cell_volume(cell, {"a": 0.0, "c": 0.0})[1] is None


def test_cell_volume_rejects_an_esd_for_a_parameter_that_is_not_free() -> None:
    with pytest.raises(ValueError, match="'b'"):
        cell_volume(Cell.tetragonal(12.5, 3.9), {"a": 0.001, "b": 0.001})
