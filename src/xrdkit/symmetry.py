"""Reflection conditions, Laue equivalence and multiplicity from the operations
of a space group.

An operation is a pair (R, t): R a 3 by 3 integer rotation matrix, as a tuple
of rows, and t a translation of three :class:`~fractions.Fraction` reduced
modulo 1, so that x' = R x + t for a fractional position x and comparisons
are exact. Operations are written in the xyz notation of the CIF
``_symmetry_equiv_pos_as_xyz`` field, such as ``"-y,x,z"`` or
``"x+1/2,-y+1/2,z"``.

A reflection h, a row of three integers, transforms as h' = h R, which is R
transposed acting on h as a column. h is systematically absent when some
operation leaves it fixed, h R = h, with h . t not an integer; that gives the
integral, zonal and serial conditions of the group at once.

Each supported space group is stored as a few generators from International
Tables Volume A in the setting its symbol names, R3c and R3m on hexagonal
axes (obverse), with the centring translations among the generators, and is
expanded by closure once.
"""

from __future__ import annotations

import re
from fractions import Fraction

__all__ = [
    "SUPPORTED_SPACE_GROUPS",
    "holohedry",
    "is_absent",
    "laue_group",
    "laue_orbit",
    "multiplicity",
    "parse_xyz",
    "representative",
    "space_group_operations",
]

Matrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
Translation = tuple[Fraction, Fraction, Fraction]
Operation = tuple[Matrix, Translation]

# Centring translations of an obverse rhombohedral cell on hexagonal axes.
OBVERSE = ("x+2/3,y+1/3,z+1/3", "x+1/3,y+2/3,z+2/3")

# Generators of each supported space group, keyed by the symbol the library
# entries use, and the order of the full group counting centring.
SPACE_GROUPS: dict[str, tuple[tuple[str, ...], int]] = {
    "Pm-3m": (("z,x,y", "-y,x,z", "-x,-y,-z"), 48),
    "P4mm": (("-y,x,z", "x,-y,z"), 8),
    "P4bm": (("-y,x,z", "x+1/2,-y+1/2,z"), 8),
    "P4/mbm": (("-y,x,z", "-x+1/2,y+1/2,-z", "-x,-y,-z"), 16),
    "R3c": (("-y,x-y,z", "-y,-x,z+1/2", *OBVERSE), 18),
    "R3m": (("-y,x-y,z", "-y,-x,z", *OBVERSE), 18),
    # The cab setting of Pnma (No. 62).
    "Pbnm": (("-x,-y,z+1/2", "x+1/2,-y+1/2,-z", "-x,-y,-z"), 8),
    "Amm2": (("-x,-y,z", "x,-y,z", "x,y+1/2,z+1/2"), 8),
}
SUPPORTED_SPACE_GROUPS = tuple(SPACE_GROUPS)

# Generators of the Laue group of each crystal system's holohedry: m-3m,
# 4/mmm, mmm, 6/mmm on hexagonal axes, 2/m with b unique, and -1.
HOLOHEDRY_GENERATORS = {
    "cubic": ("z,x,y", "-y,x,z"),
    "tetragonal": ("-y,x,z", "x,-y,z"),
    "orthorhombic": ("-x,-y,z", "x,-y,-z"),
    "hexagonal": ("x-y,x,z", "-y,-x,z"),
    "trigonal": ("x-y,x,z", "-y,-x,z"),
    "monoclinic": ("-x,y,-z",),
    "triclinic": (),
}

IDENTITY: Matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
INVERSION: Matrix = ((-1, 0, 0), (0, -1, 0), (0, 0, -1))
ZERO: Translation = (Fraction(0), Fraction(0), Fraction(0))

# Translations must be multiples of 1/12, which covers every space group.
DENOMINATOR = 12

_TERM = re.compile(r"([+-]?)(?:(\d+)/(\d+)|(\d+)|([xyz]))")

_OPERATIONS: dict[str, tuple[Operation, ...]] = {}


def parse_xyz(text: str) -> Operation:
    """The operation (R, t) of ``text`` in xyz notation, such as
    ``"-x+y,y,z+1/2"``, with t reduced modulo 1.

    Raises
    ------
    ValueError
        If ``text`` does not have three components of x, y, z terms and
        fractions, a component has no x, y or z, or a translation is not a
        multiple of 1/12.
    """
    components = text.replace(" ", "").lower().split(",")
    if len(components) != 3:
        raise ValueError(f"xyz operation {text!r} must have three components")
    rows = []
    translation = []
    for component in components:
        row = [0, 0, 0]
        shift = Fraction(0)
        position = 0
        while position < len(component):
            match = _TERM.match(component, position)
            if match is None or match.end() == position:
                raise ValueError(
                    f"cannot read {component[position:]!r} in xyz operation {text!r}"
                )
            if position and not match.group(1):
                raise ValueError(f"missing sign before a term in {text!r}")
            sign = -1 if match.group(1) == "-" else 1
            if match.group(5):
                row["xyz".index(match.group(5))] += sign
            elif match.group(2):
                shift += sign * Fraction(int(match.group(2)), int(match.group(3)))
            else:
                shift += sign * int(match.group(4))
            position = match.end()
        if row == [0, 0, 0]:
            raise ValueError(f"component {component!r} of {text!r} has no x, y or z")
        if (shift * DENOMINATOR).denominator != 1:
            raise ValueError(
                f"translation {shift} in {text!r} is not a multiple of 1/{DENOMINATOR}"
            )
        rows.append(tuple(row))
        translation.append(shift % 1)
    return tuple(rows), tuple(translation)


def _multiply(first: Matrix, second: Matrix) -> Matrix:
    return tuple(
        tuple(sum(first[i][m] * second[m][j] for m in range(3)) for j in range(3))
        for i in range(3)
    )


def _compose(first: Operation, second: Operation) -> Operation:
    """``first`` after ``second``: x -> R1 (R2 x + t2) + t1."""
    (r1, t1), (r2, t2) = first, second
    shift = tuple(
        (sum(r1[i][m] * t2[m] for m in range(3)) + t1[i]) % 1 for i in range(3)
    )
    return _multiply(r1, r2), shift


def _closure(generators: tuple[Operation, ...]) -> tuple[Operation, ...]:
    """Every operation the ``generators`` make, the identity first."""
    group = [(IDENTITY, ZERO)]
    seen = set(group)
    frontier = list(group)
    while frontier:
        found = []
        for operation in frontier:
            for generator in generators:
                product = _compose(operation, generator)
                if product not in seen:
                    seen.add(product)
                    found.append(product)
        group.extend(found)
        frontier = found
    return tuple(group)


def space_group_operations(symbol: str) -> tuple[Operation, ...]:
    """Every operation of space group ``symbol``, centring included, the
    identity first, expanded from its generators once and cached.

    Raises
    ------
    ValueError
        If ``symbol`` is not one of :data:`SUPPORTED_SPACE_GROUPS`, or its
        generators do not give the group's order.
    """
    if symbol in _OPERATIONS:
        return _OPERATIONS[symbol]
    if symbol not in SPACE_GROUPS:
        # Other symbols could take their operations from pymatgen's
        # SpaceGroup, an optional dependency, here.
        raise ValueError(
            f"unsupported space group {symbol!r}; the supported space groups are "
            + ", ".join(SUPPORTED_SPACE_GROUPS)
        )
    generators, order = SPACE_GROUPS[symbol]
    operations = _closure(tuple(parse_xyz(text) for text in generators))
    if len(operations) != order:
        raise ValueError(
            f"the generators of {symbol} give {len(operations)} operations, not {order}"
        )
    _OPERATIONS[symbol] = operations
    return operations


def _transform(hkl: tuple[int, int, int], rotation: Matrix) -> tuple[int, int, int]:
    """h R, the reflection ``hkl`` as a row times ``rotation``."""
    return tuple(sum(hkl[i] * rotation[i][j] for i in range(3)) for j in range(3))


def is_absent(hkl: tuple[int, int, int], operations: tuple[Operation, ...]) -> bool:
    """Whether reflection ``hkl`` is systematically absent under
    ``operations``: some operation has h R = h with h . t not an integer."""
    h = tuple(int(index) for index in hkl)
    for rotation, translation in operations:
        if _transform(h, rotation) == h:
            phase = sum(index * shift for index, shift in zip(h, translation))
            if Fraction(phase).denominator != 1:
                return True
    return False


def laue_group(operations: tuple[Operation, ...]) -> tuple[Matrix, ...]:
    """The Laue group of ``operations``: their distinct rotations together
    with the inversion, closed under multiplication, sorted."""
    generators = {rotation for rotation, _ in operations} | {INVERSION}
    closed = _closure(tuple((rotation, ZERO) for rotation in generators))
    return tuple(sorted(rotation for rotation, _ in closed))


def holohedry(crystal_system: str) -> tuple[Matrix, ...]:
    """The Laue group of the holohedral class of ``crystal_system``, for when
    no space group is given: m-3m, 4/mmm, mmm, 6/mmm on hexagonal axes for
    hexagonal and trigonal, 2/m with b unique, or -1.

    Raises
    ------
    ValueError
        If ``crystal_system`` is not one of the seven.
    """
    if crystal_system not in HOLOHEDRY_GENERATORS:
        raise ValueError(
            f"unknown crystal_system {crystal_system!r}; the crystal systems are "
            + ", ".join(HOLOHEDRY_GENERATORS)
        )
    generators = HOLOHEDRY_GENERATORS[crystal_system]
    return laue_group(tuple(parse_xyz(text) for text in generators))


def laue_orbit(
    hkl: tuple[int, int, int], laue: tuple[Matrix, ...]
) -> tuple[tuple[int, int, int], ...]:
    """The distinct reflections equivalent to ``hkl`` under the Laue group
    ``laue``, largest first."""
    h = tuple(int(index) for index in hkl)
    return tuple(sorted({_transform(h, rotation) for rotation in laue}, reverse=True))


def multiplicity(hkl: tuple[int, int, int], laue: tuple[Matrix, ...]) -> int:
    """The number of reflections equivalent to ``hkl`` under ``laue``."""
    return len(laue_orbit(hkl, laue))


def representative(
    hkl: tuple[int, int, int], laue: tuple[Matrix, ...]
) -> tuple[int, int, int]:
    """The label of the orbit of ``hkl`` under ``laue``: of its members with
    no negative index, or of all members if none, the lexicographically
    largest, h first, then k, then l.

    That gives h >= k >= 0 and l >= 0 for tetragonal groups, (100) for the
    cubic {100}, (110) rather than (2 -1 0) on hexagonal axes, and keeps the
    sign of l in a monoclinic (h 0 -l). On obverse hexagonal axes -3m keeps
    (h k l) and (k h l), l not 0, in separate orbits, so a label never swaps
    h and k: (101) rather than (011), which is absent in R3m, and (104) and
    (012) in R3c, where (014) and (102) are absent.
    """
    orbit = laue_orbit(hkl, laue)
    non_negative = [member for member in orbit if min(member) >= 0]
    return max(non_negative or orbit)
