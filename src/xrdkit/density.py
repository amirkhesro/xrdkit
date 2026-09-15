"""Theoretical density from a refined cell and a composition."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping

import numpy as np
from scipy.constants import Avogadro

from xrdkit.indexing import TetragonalCell

__all__ = [
    "ATOMIC_MASSES",
    "cell_volume",
    "formula_mass",
    "parse_formula",
    "relative_density",
    "theoretical_density",
]

# Standard atomic weights in g/mol, IUPAC 2021 (Prohaska et al., Pure Appl.
# Chem. 94 (2022) 573), for every element from H to U in order of atomic number,
# to four decimal places where the table gives that many. Elements whose
# standard atomic weight is an interval (H, Li, B, C, N, O, Mg, Si, S, Cl, Ar,
# Br, Tl, Pb) take the IUPAC conventional value. The 2024 table revised Gd to
# 157.249 and Zr to 91.222; the 2021 values are kept here, as are the longer
# values some entries were first given with.
#
# Tc, Pm, Po, At, Rn, Fr, Ra and Ac have no stable isotope and no standard
# atomic weight; each takes the mass number of its longest lived isotope
# (98Tc, 145Pm, 209Po, 210At, 222Rn, 223Fr, 226Ra, 227Ac). The half-lives of
# 97Tc and 98Tc agree within their uncertainties; 98 is the value tables give.
ATOMIC_MASSES: dict[str, float] = {
    "H": 1.008,
    "He": 4.0026,
    "Li": 6.94,
    "Be": 9.0122,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.9984,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815384,
    "Si": 28.085,
    "P": 30.9738,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.95,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.9559,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938043,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.9216,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.905838,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Mo": 95.95,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.9045,
    "Xe": 131.293,
    "Cs": 132.9055,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.9077,
    "Nd": 144.242,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.9254,
    "Dy": 162.500,
    "Ho": 164.9303,
    "Er": 167.259,
    "Tm": 168.9342,
    "Yb": 173.045,
    "Lu": 174.9668,
    "Hf": 178.486,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.9666,
    "Hg": 200.592,
    "Tl": 204.38,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209.0,
    "At": 210.0,
    "Rn": 222.0,
    "Fr": 223.0,
    "Ra": 226.0,
    "Ac": 227.0,
    "Th": 232.0377,
    "Pa": 231.0359,
    "U": 238.0289,
}

# Cubic angstroms in a cubic centimetre.
A3_PER_CM3 = 1e24


_ELEMENT = re.compile(r"[A-Z][a-z]?")
_COUNT = re.compile(r"\d+(?:\.\d+)?|\.\d+")


def _count(text: str, position: int) -> tuple[float, int]:
    """Read the count at ``position``, or 1 if there is none."""
    match = _COUNT.match(text, position)
    if match is None:
        return 1.0, position
    return float(match.group()), match.end()


def parse_formula(text: str) -> dict[str, float]:
    """Return the composition a chemical formula describes.

    A formula is a run of element symbols, each followed by an optional count
    that may be fractional, such as ``"Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"`` or
    ``"LaB6"``. Parentheses group symbols under one multiplier, as in
    ``"Ca(OH)2"``, and may nest. An element named more than once has its
    counts added.

    Returns
    -------
    dict[str, float]
        Element symbol to atoms per formula unit, in order of first appearance.

    Raises
    ------
    ValueError
        If the formula is empty, holds anything other than element symbols,
        counts and balanced parentheses, names an element not in
        :data:`ATOMIC_MASSES`, or has an empty pair of parentheses. The message
        quotes the offending text.
    """
    formula = text.strip()
    if not formula:
        raise ValueError("formula is empty")
    groups: list[dict[str, float]] = [{}]
    opened: list[int] = []
    position = 0
    while position < len(formula):
        element = _ELEMENT.match(formula, position)
        if element is not None:
            symbol = element.group()
            if symbol not in ATOMIC_MASSES:
                raise ValueError(f"unknown element {symbol!r} in formula {formula!r}")
            count, position = _count(formula, element.end())
            groups[-1][symbol] = groups[-1].get(symbol, 0.0) + count
        elif formula[position] == "(":
            groups.append({})
            opened.append(position)
            position += 1
        elif formula[position] == ")":
            if not opened:
                raise ValueError(
                    f"unmatched ')' at {formula[position:]!r} in formula {formula!r}"
                )
            start = opened.pop()
            group = groups.pop()
            if not group:
                raise ValueError(
                    f"empty group {formula[start : position + 1]!r} "
                    f"in formula {formula!r}"
                )
            count, position = _count(formula, position + 1)
            for symbol, n in group.items():
                groups[-1][symbol] = groups[-1].get(symbol, 0.0) + n * count
        else:
            raise ValueError(
                f"cannot read {formula[position:]!r} in formula {formula!r}"
            )
    if opened:
        raise ValueError(
            f"unmatched '(' at {formula[opened[-1] :]!r} in formula {formula!r}"
        )
    return groups[0]


def formula_mass(composition: Mapping[str, float] | str) -> float:
    """Return the mass of one formula unit in g/mol.

    Parameters
    ----------
    composition
        Element symbol to stoichiometric coefficient, such as
        ``{"Sr": 1, "Ti": 1, "O": 3}``, or a formula that :func:`parse_formula`
        reads, such as ``"SrTiO3"``.

    Raises
    ------
    ValueError
        If the composition is empty, names an element not in
        :data:`ATOMIC_MASSES`, or has a negative coefficient, or if the formula
        cannot be read.
    """
    if isinstance(composition, str):
        composition = parse_formula(composition)
    if not composition:
        raise ValueError("composition is empty")
    unknown = sorted(set(composition) - set(ATOMIC_MASSES))
    if unknown:
        raise ValueError(f"No atomic mass for {', '.join(unknown)}")
    negative = sorted(element for element, n in composition.items() if n < 0)
    if negative:
        raise ValueError(f"Negative coefficient for {', '.join(negative)}")
    return float(sum(ATOMIC_MASSES[element] * n for element, n in composition.items()))


def cell_volume(
    cell: TetragonalCell,
    esd_a: float | None = None,
    esd_c: float | None = None,
) -> tuple[float, float | None]:
    """Return the volume a^2 c of a tetragonal cell and its esd, in cubic angstroms.

    The esd is propagated to first order, sqrt((2ac esd_a)^2 + (a^2 esd_c)^2),
    treating a and c as uncorrelated since the fits here report no covariance.
    An esd left as ``None`` counts as zero, and the esd returned is ``None``
    only when neither is given.
    """
    volume = cell.a**2 * cell.c
    if esd_a is None and esd_c is None:
        return float(volume), None
    from_a = 2.0 * cell.a * cell.c * (esd_a or 0.0)
    from_c = cell.a**2 * (esd_c or 0.0)
    return float(volume), float(np.hypot(from_a, from_c))


def theoretical_density(
    composition: Mapping[str, float] | str,
    z: float,
    volume_a3: float,
    esd_volume_a3: float | None = None,
) -> tuple[float, float | None]:
    """Return the X-ray density in g/cm^3 and its esd.

    The density is Z M / (N_A V), with M the formula mass and V the cell
    volume. Only the volume carries an error here, so the density has the same
    relative esd as the volume; the esd returned is ``None`` when no volume esd
    is given.

    Parameters
    ----------
    composition
        Element symbol to stoichiometric coefficient for one formula unit, or
        a formula such as ``"Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"``.
    z
        Formula units per cell.
    volume_a3
        Cell volume in cubic angstroms.
    esd_volume_a3
        Estimated standard deviation of the volume, in cubic angstroms.

    Raises
    ------
    ValueError
        If ``z`` or ``volume_a3`` is not positive, or from :func:`formula_mass`.
    """
    if z <= 0:
        raise ValueError(f"z must be positive, got {z}")
    if volume_a3 <= 0:
        raise ValueError(f"volume_a3 must be positive, got {volume_a3}")
    density = z * formula_mass(composition) / (Avogadro * volume_a3 / A3_PER_CM3)
    if esd_volume_a3 is None:
        return float(density), None
    return float(density), float(density * esd_volume_a3 / volume_a3)


def relative_density(
    measured: float,
    theoretical: float,
    esd_measured: float | None = None,
    esd_theoretical: float | None = None,
) -> tuple[float, float | None]:
    """Return a measured density as a percentage of the theoretical one, and its esd.

    The esd adds the relative errors of the two densities in quadrature,
    ratio sqrt((esd_measured / measured)^2 + (esd_theoretical / theoretical)^2).
    An esd left as ``None`` counts as zero, and the esd returned is ``None``
    only when neither is given.

    Raises
    ------
    ValueError
        If either density is not positive.
    """
    if measured <= 0:
        raise ValueError(f"measured density must be positive, got {measured}")
    if theoretical <= 0:
        raise ValueError(f"theoretical density must be positive, got {theoretical}")
    ratio = 100.0 * measured / theoretical
    if esd_measured is None and esd_theoretical is None:
        return float(ratio), None
    relative = math.hypot(
        (esd_measured or 0.0) / measured, (esd_theoretical or 0.0) / theoretical
    )
    return float(ratio), float(ratio * relative)
