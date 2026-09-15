"""Tests for xrdkit.symmetry."""

import itertools
import math
from fractions import Fraction

import pytest

from xrdkit import (
    Cell,
    TetragonalCell,
    generate_reflections,
    holohedry,
    is_absent,
    laue_group,
    laue_orbit,
    multiplicity,
    parse_xyz,
    representative,
    space_group_operations,
    symmetry,
)
from xrdkit.symmetry import SUPPORTED_SPACE_GROUPS

ORDERS = {
    "Pm-3m": 48,
    "P4mm": 8,
    "P4bm": 8,
    "P4/mbm": 16,
    "R3c": 18,
    "R3m": 18,
    "Pbnm": 8,
    "Amm2": 8,
}

HALF = Fraction(1, 2)

# Allowed and forbidden reflections of each group, from the reflection
# conditions of International Tables Volume A.
P4BM_ALLOWED = [(0, 2, 1), (1, 1, 0), (3, 1, 0), (0, 0, 1), (2, 0, 3), (1, 2, 3)]
P4BM_FORBIDDEN = [(0, 1, 1), (1, 0, 1), (0, 3, 0), (3, 0, 2), (0, -1, 4)]
ABSENCES = {
    "Pm-3m": (
        [(1, 0, 0), (1, 1, 0), (1, 1, 1), (0, 0, 1), (2, 1, 0), (1, 2, 3)],
        [],
    ),
    "P4mm": (
        [(1, 0, 0), (0, 1, 1), (1, 0, 1), (0, 0, 1), (3, 1, 0), (1, 2, 3)],
        [],
    ),
    "P4bm": (P4BM_ALLOWED, P4BM_FORBIDDEN),
    "P4/mbm": (P4BM_ALLOWED, P4BM_FORBIDDEN),
    "R3c": (
        [(0, 1, 2), (1, 0, 4), (1, 1, 0), (0, 0, 6), (1, -1, 2), (1, 1, 3)],
        [(0, 0, 3), (0, 0, 9), (1, 0, -1), (1, 1, 1), (1, -1, 4), (2, -2, 1)],
    ),
    "R3m": (
        [(0, 0, 3), (0, 0, 6), (1, -1, 2), (1, 0, 1), (0, 1, 2), (2, -2, 1)],
        [(0, 0, 1), (0, 0, 2), (1, -1, 1), (0, 1, 1), (1, 1, 1)],
    ),
    "Pbnm": (
        [(0, 2, 1), (1, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 1, 1)],
        [(0, 1, 1), (1, 0, 2), (1, 0, 0), (0, 1, 0), (0, 0, 1), (3, 0, 0)],
    ),
    "Amm2": (
        [(1, 1, 1), (0, 1, 1), (1, 0, 2), (1, 2, 0), (0, 0, 2), (2, 0, 0)],
        [(1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (2, 1, 2), (0, 3, 2)],
    ),
}


@pytest.mark.parametrize(
    ("text", "rotation", "translation"),
    [
        ("x,y,z", ((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0)),
        ("-y,x,z", ((0, -1, 0), (1, 0, 0), (0, 0, 1)), (0, 0, 0)),
        (
            "x+1/2,-y+1/2,z",
            ((1, 0, 0), (0, -1, 0), (0, 0, 1)),
            (HALF, HALF, 0),
        ),
        ("-y,x-y,z", ((0, -1, 0), (1, -1, 0), (0, 0, 1)), (0, 0, 0)),
        (
            "-x+y, 1/2-y ,Z-1/3",
            ((-1, 1, 0), (0, -1, 0), (0, 0, 1)),
            (0, HALF, Fraction(2, 3)),
        ),
        (
            "x+2/3,y+1/3,z+4/3",
            ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            (Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
        ),
    ],
)
def test_parse_xyz(text, rotation, translation):
    parsed_rotation, parsed_translation = parse_xyz(text)
    assert parsed_rotation == rotation
    assert parsed_translation == tuple(Fraction(value) for value in translation)
    assert all(isinstance(value, Fraction) for value in parsed_translation)


@pytest.mark.parametrize("text", ["x,y", "x,y,1/2", "x,y,z+1/7", "x,yz,z", "x,y,w"])
def test_parse_xyz_rejects_bad_text(text):
    with pytest.raises(ValueError):
        parse_xyz(text)


def test_supported_space_groups():
    assert SUPPORTED_SPACE_GROUPS == tuple(ORDERS)


@pytest.mark.parametrize("symbol", SUPPORTED_SPACE_GROUPS)
def test_group_order_identity_and_closure(symbol):
    operations = space_group_operations(symbol)
    assert len(operations) == ORDERS[symbol]
    assert len(set(operations)) == len(operations)
    assert parse_xyz("x,y,z") in operations
    group = set(operations)
    for first, second in itertools.product(operations, repeat=2):
        assert symmetry._compose(first, second) in group


def test_operations_are_cached():
    assert space_group_operations("P4bm") is space_group_operations("P4bm")


def test_unknown_space_group_raises():
    with pytest.raises(ValueError, match="Pnma.*P4bm"):
        space_group_operations("Pnma")


def test_wrong_order_raises(monkeypatch):
    monkeypatch.setitem(symmetry.SPACE_GROUPS, "P4", (("-y,x,z",), 8))
    with pytest.raises(ValueError, match="give 4 operations, not 8"):
        space_group_operations("P4")
    assert "P4" not in symmetry._OPERATIONS


@pytest.mark.parametrize("symbol", SUPPORTED_SPACE_GROUPS)
def test_absences(symbol):
    operations = space_group_operations(symbol)
    allowed, forbidden = ABSENCES[symbol]
    assert len(allowed) >= 4
    assert len(forbidden) >= 4 or symbol in ("Pm-3m", "P4mm")
    for hkl in allowed:
        assert not is_absent(hkl, operations), hkl
    for hkl in forbidden:
        assert is_absent(hkl, operations), hkl


@pytest.mark.parametrize("symbol", ["Pm-3m", "P4mm"])
def test_primitive_groups_without_glides_forbid_nothing(symbol):
    operations = space_group_operations(symbol)
    box = itertools.product(range(-3, 4), repeat=3)
    assert not any(is_absent(hkl, operations) for hkl in box if any(hkl))


def test_r3c_integral_and_glide_conditions_over_a_box():
    operations = space_group_operations("R3c")
    for h, k, l in itertools.product(range(-4, 5), repeat=3):
        integral = (-h + k + l) % 3 == 0
        # The c glide acts on the h -h l zone and its equivalents h0l and 0kl.
        glide = (h != 0 and k != 0 and k != -h) or l % 2 == 0
        assert is_absent((h, k, l), operations) == (not (integral and glide))


@pytest.mark.parametrize("symbol", ["Pbnm", "Amm2"])
def test_orthorhombic_conditions_over_a_box(symbol):
    operations = space_group_operations(symbol)
    for h, k, l in itertools.product(range(-4, 5), repeat=3):
        if symbol == "Pbnm":
            present = (
                (h != 0 or k % 2 == 0)
                and (k != 0 or (h + l) % 2 == 0)
                and not (k == l == 0 and h % 2)
                and not (h == l == 0 and k % 2)
                and not (h == k == 0 and l % 2)
            )
        else:
            present = (k + l) % 2 == 0
        assert is_absent((h, k, l), operations) == (not present), (h, k, l)


@pytest.mark.parametrize(
    ("system", "hkl", "expected"),
    [
        ("cubic", (1, 0, 0), 6),
        ("cubic", (1, 1, 0), 12),
        ("cubic", (1, 1, 1), 8),
        ("cubic", (1, 2, 3), 48),
        ("tetragonal", (1, 0, 0), 4),
        ("tetragonal", (0, 0, 1), 2),
        ("tetragonal", (1, 1, 0), 4),
        ("tetragonal", (2, 1, 0), 8),
        ("tetragonal", (3, 1, 0), 8),
        ("orthorhombic", (1, 0, 0), 2),
        ("orthorhombic", (1, 1, 1), 8),
        ("hexagonal", (1, 0, 0), 6),
        ("hexagonal", (0, 0, 1), 2),
        ("hexagonal", (1, 0, 1), 12),
    ],
)
def test_multiplicity_under_holohedry(system, hkl, expected):
    laue = holohedry(system)
    assert multiplicity(hkl, laue) == expected
    assert len(laue_orbit(hkl, laue)) == expected
    assert hkl in laue_orbit(hkl, laue)


def test_r3c_laue_orbit_of_101():
    laue = laue_group(space_group_operations("R3c"))
    assert len(laue) == 12
    orbit = laue_orbit((1, 0, 1), laue)
    assert len(orbit) == 6
    assert (0, 1, 1) not in orbit


def test_monoclinic_101_and_10_1_are_in_different_orbits():
    laue = holohedry("monoclinic")
    assert (1, 0, -1) not in laue_orbit((1, 0, 1), laue)
    assert multiplicity((1, 0, 1), laue) == 2


@pytest.mark.parametrize(
    ("system", "orbit", "label"),
    [
        (
            "tetragonal",
            itertools.product([3, -3, 1, -1], [3, -3, 1, -1], [2, -2]),
            (3, 1, 2),
        ),
        ("cubic", [(0, 0, -1), (0, 1, 0), (-1, 0, 0)], (1, 0, 0)),
        ("cubic", [(-1, 2, -3), (3, 1, 2)], (3, 2, 1)),
        ("hexagonal", [(2, -1, 0), (-1, -1, 0), (1, -2, 0)], (1, 1, 0)),
        ("monoclinic", [(1, 0, -1), (-1, 0, 1)], (1, 0, -1)),
        ("monoclinic", [(-1, 0, -1), (1, 0, 1)], (1, 0, 1)),
    ],
)
def test_representative_under_holohedry(system, orbit, label):
    laue = holohedry(system)
    for hkl in orbit:
        if abs(hkl[0]) == abs(hkl[1]) and system == "tetragonal":
            continue
        assert representative(hkl, laue) == label, hkl


@pytest.mark.parametrize(
    ("symbol", "label", "absent"),
    [
        ("R3m", (1, 0, 1), (0, 1, 1)),
        ("R3c", (1, 0, 4), (0, 1, 4)),
        ("R3c", (0, 1, 2), (1, 0, 2)),
    ],
)
def test_representative_on_obverse_axes_avoids_absent_label(symbol, label, absent):
    operations = space_group_operations(symbol)
    laue = laue_group(operations)
    for hkl in laue_orbit(label, laue):
        assert representative(hkl, laue) == label
    assert absent not in laue_orbit(label, laue)
    assert is_absent(absent, operations)
    assert not is_absent(label, operations)


def test_representative_matches_indexing_convention_for_p4bm():
    a, c, wavelength, two_theta_max = 12.45, 3.94, 1.5406, 60.0
    expected = {
        reflection.hkl
        for reflection in generate_reflections(
            TetragonalCell(a=a, c=c), wavelength, two_theta_max, space_group="P4bm"
        )
    }
    operations = space_group_operations("P4bm")
    laue = laue_group(operations)
    cell = Cell.tetragonal(a, c)
    d_min = wavelength / (2 * math.sin(math.radians(two_theta_max / 2)))
    h_max, l_max = int(a / d_min), int(c / d_min)
    found = set()
    for hkl in itertools.product(
        range(-h_max, h_max + 1), range(-h_max, h_max + 1), range(-l_max, l_max + 1)
    ):
        if (
            any(hkl)
            and not is_absent(hkl, operations)
            and cell.d_spacing(*hkl) >= d_min
        ):
            found.add(representative(hkl, laue))
    assert len(expected) > 50
    assert found == expected


@pytest.mark.parametrize(
    ("system", "order"),
    [
        ("cubic", 48),
        ("tetragonal", 16),
        ("orthorhombic", 8),
        ("hexagonal", 24),
        ("trigonal", 24),
        ("monoclinic", 4),
        ("triclinic", 2),
    ],
)
def test_holohedry_order(system, order):
    laue = holohedry(system)
    assert len(laue) == order
    assert len(set(laue)) == order
    group = set(laue)
    for first, second in itertools.product(laue, repeat=2):
        assert symmetry._multiply(first, second) in group


def test_holohedry_unknown_system_raises():
    with pytest.raises(ValueError, match="unknown crystal_system"):
        holohedry("rhombic")


def test_laue_group_of_p4bm_is_4_mmm():
    assert set(laue_group(space_group_operations("P4bm"))) == set(
        holohedry("tetragonal")
    )
