"""Theoretical density from a refined cell and a composition."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.constants import Avogadro

from xrdkit.indexing import TetragonalCell

__all__ = [
    "ATOMIC_MASSES",
    "cell_volume",
    "formula_mass",
    "theoretical_density",
]

# Standard atomic weights in g/mol, IUPAC 2021 (Prohaska et al., Pure Appl.
# Chem. 94 (2022) 573). Elements whose standard atomic weight is an interval
# (Li, O, Mg, Si, Pb) take the IUPAC conventional value. The 2024 table revised
# Gd to 157.249 and Zr to 91.222; the 2021 values are kept here.
ATOMIC_MASSES: dict[str, float] = {
    "Li": 6.94,
    "O": 15.999,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815384,
    "Si": 28.085,
    "K": 39.0983,
    "Ca": 40.078,
    "Ti": 47.867,
    "Cr": 51.9961,
    "Mn": 54.938043,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Sr": 87.62,
    "Y": 88.905838,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Sn": 118.710,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Nd": 144.242,
    "Sm": 150.36,
    "Gd": 157.25,
    "Hf": 178.486,
    "Ta": 180.94788,
    "W": 183.84,
    "Pb": 207.2,
    "Bi": 208.98040,
}

# Cubic angstroms in a cubic centimetre.
A3_PER_CM3 = 1e24


def formula_mass(composition: Mapping[str, float]) -> float:
    """Return the mass of one formula unit in g/mol.

    Parameters
    ----------
    composition
        Element symbol to stoichiometric coefficient, such as
        ``{"Sr": 1, "Ti": 1, "O": 3}``.

    Raises
    ------
    ValueError
        If the composition is empty, names an element not in
        :data:`ATOMIC_MASSES`, or has a negative coefficient.
    """
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
    composition: Mapping[str, float],
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
        Element symbol to stoichiometric coefficient for one formula unit.
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
