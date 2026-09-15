"""A general unit cell: six parameters and a crystal system.

Lengths are in angstroms and angles in degrees. A :class:`Cell` always holds
all six parameters, and its crystal system decides which of them are free, in
the order :data:`xrdkit.library.CELL_PARAMETERS` lists them; the rest follow
from the free ones (b = a in a tetragonal cell, say). Monoclinic cells take b
as the unique axis. Hexagonal and trigonal cells are on hexagonal axes, so a
rhombohedral setting has to be converted to hexagonal axes first. d spacings
come from the reciprocal metric, 1/d^2 = h . G* h, with G* the inverse of the
direct metric tensor G, and the volume is sqrt(det G).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from xrdkit.library import CELL_PARAMETERS, StructureEntry
from xrdkit.structure import metric_tensor

__all__ = ["Cell"]

LENGTHS = ("a", "b", "c")
ANGLES = ("alpha", "beta", "gamma")
PARAMETERS = LENGTHS + ANGLES

# How far a dependent parameter may stray from the value its crystal system
# fixes, in angstroms or degrees.
TOLERANCE = 1e-6

# A cell whose squared volume, over (abc)^2, is no more than this is flat.
FLAT = 1e-12


def _dependent(system: str, free: Mapping[str, float]) -> dict[str, float]:
    """All six parameters of a ``system`` cell from its free ones."""
    if system == "triclinic":
        return {name: free[name] for name in PARAMETERS}
    a = free["a"]
    six = {"a": a, "b": a, "c": a, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    if system in ("tetragonal", "hexagonal", "trigonal"):
        six["c"] = free["c"]
    if system in ("orthorhombic", "monoclinic"):
        six["b"] = free["b"]
        six["c"] = free["c"]
    if system in ("hexagonal", "trigonal"):
        six["gamma"] = 120.0
    if system == "monoclinic":
        six["beta"] = free["beta"]
    return six


@dataclass(frozen=True)
class Cell:
    """A unit cell of any crystal system, in angstroms and degrees.

    Build one with a named constructor (:meth:`cubic`, :meth:`tetragonal` and
    so on), :meth:`from_parameters` or :meth:`from_entry`; the six parameters
    are checked against ``crystal_system`` either way.

    Raises
    ------
    ValueError
        Naming the parameter, if a length or angle is not positive, an angle
        is not below 180, a parameter disagrees with the crystal system by
        more than 1e-6, the crystal system is unknown, or the angles make no
        cell (the metric is not positive definite).
    """

    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    crystal_system: str

    def __post_init__(self) -> None:
        for name in PARAMETERS:
            value = getattr(self, name)
            try:
                number = float(value)
            except (TypeError, ValueError):
                raise ValueError(
                    f"cell parameter {name} must be a number, not {value!r}"
                ) from None
            if not number > 0:
                raise ValueError(
                    f"cell parameter {name} must be greater than 0, not {value!r}"
                )
            if name in ANGLES and not number < 180:
                raise ValueError(
                    f"cell parameter {name} must be below 180 degrees, not {value!r}"
                )
            object.__setattr__(self, name, number)

        system = self.crystal_system
        if system not in CELL_PARAMETERS:
            raise ValueError(
                f"unknown crystal_system {system!r}; the crystal systems are "
                + ", ".join(CELL_PARAMETERS)
            )
        free = {name: getattr(self, name) for name in CELL_PARAMETERS[system]}
        for name, expected in _dependent(system, free).items():
            actual = getattr(self, name)
            if abs(actual - expected) > TOLERANCE:
                raise ValueError(
                    f"cell parameter {name} must be {expected:g} in a {system} "
                    f"cell, not {actual:g}"
                )

        cosines = np.cos(np.radians([self.alpha, self.beta, self.gamma]))
        product = 2 * cosines.prod()
        if 1 - (cosines**2).sum() + product <= FLAT:
            raise ValueError(
                f"cell parameters alpha, beta and gamma ({self.alpha:g}, "
                f"{self.beta:g}, {self.gamma:g}) make no cell: the metric is not "
                "positive definite"
            )

    @classmethod
    def cubic(cls, a: float) -> Cell:
        """A cubic cell of edge ``a``."""
        return cls._build("cubic", {"a": a})

    @classmethod
    def tetragonal(cls, a: float, c: float) -> Cell:
        """A tetragonal cell, b = a."""
        return cls._build("tetragonal", {"a": a, "c": c})

    @classmethod
    def orthorhombic(cls, a: float, b: float, c: float) -> Cell:
        """An orthorhombic cell, every angle 90."""
        return cls._build("orthorhombic", {"a": a, "b": b, "c": c})

    @classmethod
    def hexagonal(cls, a: float, c: float) -> Cell:
        """A hexagonal cell, b = a and gamma = 120."""
        return cls._build("hexagonal", {"a": a, "c": c})

    @classmethod
    def trigonal(cls, a: float, c: float) -> Cell:
        """A trigonal cell on hexagonal axes, b = a and gamma = 120.

        A cell in the rhombohedral setting has to be converted to hexagonal
        axes first.
        """
        return cls._build("trigonal", {"a": a, "c": c})

    @classmethod
    def monoclinic(cls, a: float, b: float, c: float, beta: float) -> Cell:
        """A monoclinic cell with b the unique axis, alpha = gamma = 90."""
        return cls._build("monoclinic", {"a": a, "b": b, "c": c, "beta": beta})

    @classmethod
    def triclinic(
        cls, a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ) -> Cell:
        """A triclinic cell, every parameter free."""
        return cls(a, b, c, alpha, beta, gamma, "triclinic")

    @classmethod
    def _build(cls, system: str, free: Mapping[str, float]) -> Cell:
        return cls(**_dependent(system, free), crystal_system=system)

    @classmethod
    def from_parameters(
        cls, crystal_system: str, parameters: Mapping[str, float]
    ) -> Cell:
        """The ``crystal_system`` cell of ``parameters``, which gives either
        exactly the parameters the system leaves free, as
        :data:`~xrdkit.library.CELL_PARAMETERS` lists them, or all six of a,
        b, c, alpha, beta and gamma, which must agree with the system.

        Raises
        ------
        ValueError
            If the crystal system is unknown, or naming the key, if a key is
            unknown, missing or not free for the system, or a value is out of
            range (see :class:`Cell`).
        """
        if crystal_system not in CELL_PARAMETERS:
            raise ValueError(
                f"unknown crystal_system {crystal_system!r}; the crystal systems "
                "are " + ", ".join(CELL_PARAMETERS)
            )
        allowed = CELL_PARAMETERS[crystal_system]
        unknown = [name for name in parameters if name not in PARAMETERS]
        if unknown:
            raise ValueError(
                f"unknown cell parameter {', '.join(map(repr, unknown))}; the cell "
                "parameters are " + ", ".join(PARAMETERS)
            )
        if set(parameters) == set(PARAMETERS):
            return cls(**parameters, crystal_system=crystal_system)
        given = ", ".join(allowed)
        extra = [name for name in parameters if name not in allowed]
        if extra:
            raise ValueError(
                f"cell parameter {', '.join(map(repr, extra))} is not free in a "
                f"{crystal_system} cell; give {given} or all six"
            )
        missing = [name for name in allowed if name not in parameters]
        if missing:
            raise ValueError(
                f"missing cell parameter {', '.join(map(repr, missing))}; a "
                f"{crystal_system} cell needs {given}"
            )
        return cls._build(crystal_system, parameters)

    @classmethod
    def from_entry(cls, entry: StructureEntry, parameters: Mapping[str, float]) -> Cell:
        """The cell of library ``entry``, whose crystal system it takes, with
        the values ``parameters`` gives for the entry's cell parameters (or
        all six); see :meth:`from_parameters`.

        An entry names the parameters its cell leaves free but holds no
        values for them, so they come separately.
        """
        return cls.from_parameters(entry.crystal_system, parameters)

    @property
    def metric_tensor(self) -> np.ndarray:
        """The direct metric tensor G, so a fractional vector v has length
        sqrt(v . G v)."""
        return metric_tensor(
            (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        )

    @cached_property
    def reciprocal_metric(self) -> np.ndarray:
        """The reciprocal metric tensor G*, the inverse of G, read only."""
        inverse = np.linalg.inv(self.metric_tensor)
        inverse.flags.writeable = False
        return inverse

    @property
    def volume(self) -> float:
        """The cell volume in cubic angstroms, sqrt(det G)."""
        return float(np.sqrt(np.linalg.det(self.metric_tensor)))

    def d_spacing(self, h: int, k: int, l: int) -> float:
        """The d spacing of ``(h k l)``, from 1/d^2 = h . G* h.

        Raises
        ------
        ValueError
            If ``(h k l)`` is ``(0 0 0)``, which has no d spacing.
        """
        if h == 0 and k == 0 and l == 0:
            raise ValueError("(000) has no d spacing")
        vector = np.array([h, k, l], dtype=float)
        return float(1.0 / np.sqrt(vector @ self.reciprocal_metric @ vector))

    def d_spacings(self, hkl: np.ndarray) -> np.ndarray:
        """The d spacings of the rows of ``hkl``, an array of shape (n, 3).

        Raises
        ------
        ValueError
            If ``hkl`` is not of shape (n, 3), or a row is ``(0 0 0)``.
        """
        indices = np.asarray(hkl, dtype=float)
        if indices.ndim != 2 or indices.shape[1] != 3:
            raise ValueError(f"hkl must have shape (n, 3), not {indices.shape}")
        if np.any(np.all(indices == 0, axis=1)):
            raise ValueError("(000) has no d spacing")
        inverse_squared = np.einsum(
            "ij,jk,ik->i", indices, self.reciprocal_metric, indices
        )
        return 1.0 / np.sqrt(inverse_squared)

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """The names of the parameters the crystal system leaves free."""
        return CELL_PARAMETERS[self.crystal_system]

    @property
    def parameters(self) -> dict[str, float]:
        """The free parameters by name, in the order of :attr:`parameter_names`."""
        return {name: getattr(self, name) for name in self.parameter_names}

    def replace(self, **changes: float) -> Cell:
        """A cell of the same crystal system with the free parameters named
        in ``changes`` set to new values.

        Raises
        ------
        ValueError
            Naming the parameter, if one is not free for the crystal system
            or its new value is out of range.
        """
        fixed = [name for name in changes if name not in self.parameter_names]
        if fixed:
            raise ValueError(
                f"cell parameter {', '.join(map(repr, fixed))} is not free in a "
                f"{self.crystal_system} cell; the free parameters are "
                + ", ".join(self.parameter_names)
            )
        return self.from_parameters(self.crystal_system, {**self.parameters, **changes})

    def to_dict(self) -> dict[str, float | str]:
        """All six parameters by name, and ``crystal_system``."""
        return {name: getattr(self, name) for name in PARAMETERS} | {
            "crystal_system": self.crystal_system
        }
