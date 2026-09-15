"""Tests for xrdkit.cell."""

import math

import numpy as np
import pytest

from xrdkit import Cell, TetragonalCell
from xrdkit.library import CELL_PARAMETERS, list_entries, load_entry

REL = 1e-10

# Synthetic values for any cell parameter, used where only the names matter.
VALUES = {"a": 5.0, "b": 6.0, "c": 7.0, "alpha": 80.0, "beta": 100.0, "gamma": 95.0}

# One cell of each crystal system.
CELLS = {
    "cubic": Cell.cubic(3.905),
    "tetragonal": Cell.tetragonal(12.45, 3.94),
    "orthorhombic": Cell.orthorhombic(5.38, 5.44, 7.64),
    "hexagonal": Cell.hexagonal(5.148, 13.863),
    "trigonal": Cell.trigonal(5.148, 13.863),
    "monoclinic": Cell.monoclinic(5.0, 6.0, 7.0, 100.0),
    "triclinic": Cell.triclinic(5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
}


def _triclinic_inverse_squared(a, b, c, alpha, beta, gamma, h, k, l):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    sa, sb, sg = (math.sin(math.radians(x)) for x in (alpha, beta, gamma))
    volume = _triclinic_volume(a, b, c, alpha, beta, gamma)
    numerator = (
        h * h * b * b * c * c * sa * sa
        + k * k * a * a * c * c * sb * sb
        + l * l * a * a * b * b * sg * sg
        + 2 * h * k * a * b * c * c * (ca * cb - cg)
        + 2 * k * l * a * a * b * c * (cb * cg - ca)
        + 2 * h * l * a * b * b * c * (ca * cg - cb)
    )
    return numerator / volume**2


def _triclinic_volume(a, b, c, alpha, beta, gamma):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    return a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)


@pytest.mark.parametrize("hkl", [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0)])
def test_cubic_d_spacing(hkl):
    a = 3.905
    expected = a / math.sqrt(sum(i * i for i in hkl))
    assert Cell.cubic(a).d_spacing(*hkl) == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize("hkl", [(0, 0, 1), (3, 1, 0), (3, 1, 1)])
def test_tetragonal_d_spacing(hkl):
    a, c = 12.45, 3.94
    h, k, l = hkl
    expected = 1 / math.sqrt((h * h + k * k) / a**2 + l * l / c**2)
    d = Cell.tetragonal(a, c).d_spacing(*hkl)
    assert d == pytest.approx(expected, rel=REL)
    assert abs(d - TetragonalCell(a=a, c=c).d_spacing(*hkl)) < 1e-12


@pytest.mark.parametrize("hkl", [(1, 1, 0), (0, 2, 0), (1, 1, 2)])
def test_orthorhombic_d_spacing(hkl):
    a, b, c = 5.38, 5.44, 7.64
    h, k, l = hkl
    expected = 1 / math.sqrt(h * h / a**2 + k * k / b**2 + l * l / c**2)
    d = Cell.orthorhombic(a, b, c).d_spacing(*hkl)
    assert d == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize("system", ["hexagonal", "trigonal"])
@pytest.mark.parametrize("hkl", [(0, 1, 2), (1, 0, 4), (1, 1, 0), (0, 0, 6)])
def test_hexagonal_and_trigonal_d_spacing(system, hkl):
    a, c = 5.148, 13.863
    h, k, l = hkl
    expected = 1 / math.sqrt(4 * (h * h + h * k + k * k) / (3 * a * a) + l * l / c**2)
    cell = getattr(Cell, system)(a, c)
    assert cell.crystal_system == system
    assert cell.d_spacing(*hkl) == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize("hkl", [(1, 0, 0), (0, 0, 1), (1, 0, 1), (1, 0, -1)])
def test_monoclinic_d_spacing(hkl):
    a, b, c, beta = 5.0, 6.0, 7.0, 100.0
    h, k, l = hkl
    sin_beta = math.sin(math.radians(beta))
    cos_beta = math.cos(math.radians(beta))
    inverse_squared = (
        h * h / a**2
        + k * k * sin_beta**2 / b**2
        + l * l / c**2
        - 2 * h * l * cos_beta / (a * c)
    ) / sin_beta**2
    cell = Cell.monoclinic(a, b, c, beta)
    assert cell.d_spacing(*hkl) == pytest.approx(
        1 / math.sqrt(inverse_squared), rel=REL
    )


def test_monoclinic_101_and_10_1_differ():
    cell = Cell.monoclinic(5.0, 6.0, 7.0, 100.0)
    assert abs(cell.d_spacing(1, 0, 1) - cell.d_spacing(1, 0, -1)) > 0.1


def test_triclinic_d_spacing():
    parameters = (5.0, 6.0, 7.0, 80.0, 90.0, 100.0)
    expected = 1 / math.sqrt(_triclinic_inverse_squared(*parameters, 1, 1, 1))
    cell = Cell.triclinic(*parameters)
    assert cell.d_spacing(1, 1, 1) == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize(
    ("cell", "expected"),
    [
        (Cell.cubic(3.905), 3.905**3),
        (Cell.tetragonal(12.45, 3.94), 12.45**2 * 3.94),
        (Cell.hexagonal(5.148, 13.863), math.sqrt(3) / 2 * 5.148**2 * 13.863),
        (
            Cell.monoclinic(5.0, 6.0, 7.0, 100.0),
            5.0 * 6.0 * 7.0 * math.sin(math.radians(100.0)),
        ),
        (
            Cell.triclinic(5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
            _triclinic_volume(5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
        ),
    ],
    ids=["cubic", "tetragonal", "hexagonal", "monoclinic", "triclinic"],
)
def test_volume(cell, expected):
    assert cell.volume == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize("system", list(CELLS))
def test_d_spacings_matches_d_spacing(system):
    cell = CELLS[system]
    hkl = np.array([[1, 0, 0], [0, 1, 1], [1, 1, 1], [2, -1, 3], [1, 0, -1]])
    d = cell.d_spacings(hkl)
    assert d.shape == (len(hkl),)
    np.testing.assert_allclose(d, [cell.d_spacing(*row) for row in hkl], rtol=1e-14)


def test_reciprocal_metric_is_inverse_of_metric():
    cell = CELLS["triclinic"]
    np.testing.assert_allclose(
        cell.metric_tensor @ cell.reciprocal_metric, np.eye(3), atol=1e-12
    )
    assert cell.reciprocal_metric is cell.reciprocal_metric


def test_d_spacing_of_000_raises():
    with pytest.raises(ValueError, match="000"):
        Cell.cubic(3.905).d_spacing(0, 0, 0)
    with pytest.raises(ValueError, match="000"):
        Cell.cubic(3.905).d_spacings([[1, 0, 0], [0, 0, 0]])


@pytest.mark.parametrize("a", [0.0, -3.905])
def test_non_positive_length_raises(a):
    with pytest.raises(ValueError, match="parameter a "):
        Cell.cubic(a)


def test_angle_of_180_raises():
    with pytest.raises(ValueError, match="parameter beta "):
        Cell.monoclinic(5.0, 6.0, 7.0, 180.0)


@pytest.mark.parametrize("gamma", [120.0, 130.0])
def test_angles_that_make_no_cell_raise(gamma):
    with pytest.raises(ValueError, match="positive definite"):
        Cell.triclinic(5.0, 6.0, 7.0, 120.0, 120.0, gamma)


def test_tetragonal_six_keys_with_b_not_a_raises():
    six = {"a": 12.45, "b": 12.5, "c": 3.94, "alpha": 90, "beta": 90, "gamma": 90}
    with pytest.raises(ValueError, match="parameter b "):
        Cell.from_parameters("tetragonal", six)


def test_missing_key_raises():
    with pytest.raises(ValueError, match="missing cell parameter 'c'"):
        Cell.from_parameters("tetragonal", {"a": 12.45})


def test_unknown_key_raises():
    with pytest.raises(ValueError, match="unknown cell parameter 'd'"):
        Cell.from_parameters("tetragonal", {"a": 12.45, "c": 3.94, "d": 1.0})


def test_key_not_free_raises():
    with pytest.raises(ValueError, match="'b' is not free"):
        Cell.from_parameters("tetragonal", {"a": 12.45, "b": 12.45, "c": 3.94})


def test_unknown_crystal_system_raises():
    with pytest.raises(ValueError, match="unknown crystal_system 'rhombic'"):
        Cell.from_parameters("rhombic", {"a": 5.0})
    with pytest.raises(ValueError, match="unknown crystal_system 'rhombic'"):
        Cell(5.0, 5.0, 5.0, 90.0, 90.0, 90.0, "rhombic")


def test_replace_with_name_not_free_raises():
    with pytest.raises(ValueError, match="'b' is not free"):
        CELLS["tetragonal"].replace(b=12.5)


@pytest.mark.parametrize("system", list(CELLS))
def test_own_keys_and_six_keys_give_equal_cells(system):
    cell = CELLS[system]
    six = {name: getattr(cell, name) for name in VALUES}
    own = Cell.from_parameters(system, cell.parameters)
    assert own == Cell.from_parameters(system, six) == cell


@pytest.mark.parametrize("name", list_entries())
def test_from_entry_loads_every_library_entry(name):
    entry = load_entry(name)
    values = {parameter: VALUES[parameter] for parameter in entry.cell_parameters}
    cell = Cell.from_entry(entry, values)
    assert cell.crystal_system == entry.crystal_system
    assert cell.parameters == values


@pytest.mark.parametrize("system", list(CELLS))
def test_parameters_and_parameter_names(system):
    cell = CELLS[system]
    assert cell.parameter_names == CELL_PARAMETERS[system]
    assert tuple(cell.parameters) == CELL_PARAMETERS[system]
    assert cell.parameters == {
        name: getattr(cell, name) for name in CELL_PARAMETERS[system]
    }


@pytest.mark.parametrize("system", list(CELLS))
def test_replace_changes_only_the_named_parameter(system):
    cell = CELLS[system]
    for name in cell.parameter_names:
        changed = cell.replace(**{name: getattr(cell, name) + 1.0})
        assert changed.crystal_system == system
        for other, value in cell.parameters.items():
            expected = value + 1.0 if other == name else value
            assert changed.parameters[other] == expected


@pytest.mark.parametrize("system", list(CELLS))
def test_to_dict_round_trips(system):
    cell = CELLS[system]
    data = cell.to_dict()
    assert set(data) == set(VALUES) | {"crystal_system"}
    crystal_system = data.pop("crystal_system")
    assert Cell.from_parameters(crystal_system, data) == cell


def test_cells_are_hashable_and_compare_equal():
    first = Cell.tetragonal(12.45, 3.94)
    second = Cell.from_parameters("tetragonal", {"a": 12.45, "c": 3.94})
    _ = first.reciprocal_metric  # caching must not change equality or hash
    assert first == second
    assert hash(first) == hash(second)
    assert len({first, second, Cell.tetragonal(12.45, 3.95)}) == 2
