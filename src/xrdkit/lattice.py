"""Lattice parameter refinement against indexed peaks, with systematic errors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize

from xrdkit.indexing import IndexedPeak, TetragonalCell

__all__ = ["LatticeFit", "lattice_fit_to_dict", "refine_lattice"]

# A peak needs this many more than the parameters before a fit means anything.
MIN_EXTRA_PEAKS = 3

# Keeps the arcsine in range while the optimiser tries a cell that cannot
# diffract at all; the residual there is large either way.
MAX_SIN_THETA = 1.0 - 1e-9

# Order of the parameter vector the model is written in.
PARAMETER_NAMES = ("a", "c", "zero", "displacement")


@dataclass
class LatticeFit:
    """The outcome of a lattice parameter refinement.

    ``zero`` is the zero point error in degrees, added to every calculated
    position. ``displacement`` is the specimen displacement in millimetres,
    positive below the focusing circle, or ``None`` when it was not refined.
    An estimated standard deviation of a parameter that was held fixed is
    ``None``.
    """

    a: float
    c: float
    esd_a: float
    esd_c: float
    zero: float
    esd_zero: float | None
    displacement: float | None
    esd_displacement: float | None
    radius_mm: float | None
    n_peaks: int
    rms_two_theta: float
    residuals: np.ndarray
    hkl: list[tuple[int, int, int]]
    converged: bool

    @property
    def cell(self) -> TetragonalCell:
        """The refined cell."""
        return TetragonalCell(a=self.a, c=self.c)


def _two_theta(
    a: float, c: float, hkl: np.ndarray, wavelength: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the Bragg 2theta of ``hkl`` for a tetragonal cell, and its theta.

    Both in degrees and radians respectively, from 1/d^2 = (h^2+k^2)/a^2 +
    l^2/c^2 and sin(theta) = lambda / 2d.
    """
    hk, ll = hkl[:, 0] ** 2 + hkl[:, 1] ** 2, hkl[:, 2] ** 2
    inverse_squared = hk / a**2 + ll / c**2
    d_spacing = 1.0 / np.sqrt(inverse_squared)
    sin_theta = np.clip(wavelength / (2.0 * d_spacing), -MAX_SIN_THETA, MAX_SIN_THETA)
    theta = np.arcsin(sin_theta)
    return 2.0 * np.degrees(theta), theta


def _model(
    parameters: np.ndarray, hkl: np.ndarray, wavelength: float, radius_mm: float | None
) -> np.ndarray:
    """Calculated 2theta for ``parameters`` = (a, c, zero, displacement).

    The specimen displacement shifts a reflection by -2 s cos(theta) / R
    radians, which is the usual flat plate term: a specimen sitting below the
    focusing circle moves every peak to low angle, most of all at low angle.
    """
    a, c, zero, displacement = parameters
    calculated, theta = _two_theta(a, c, hkl, wavelength)
    shift = zero
    if radius_mm is not None:
        shift = shift + np.degrees(-2.0 * displacement * np.cos(theta) / radius_mm)
    return calculated + shift


def refine_lattice(
    indexed: list[IndexedPeak],
    wavelength: float,
    start_cell: TetragonalCell,
    fit_zero: bool = True,
    fit_displacement: bool = False,
    radius_mm: float | None = None,
) -> LatticeFit:
    """Refine a tetragonal cell, and its systematic errors, by least squares.

    Each peak is modelled as its Bragg position for the cell plus a zero point
    error plus a specimen displacement term, and the observed positions are
    taken as measured: ``peak.two_theta``, not ``corrected_two_theta``, since
    the zero point is what is being refined rather than assumed.

    Only peaks that matched exactly one reflection are used, so an assignment
    that was a choice between rivals cannot pull the cell.

    A zero point error and a specimen displacement are not independent over a
    short angular range, since one is constant and the other follows
    cos(theta); refining both needs peaks spread widely enough in 2theta to
    tell them apart.

    Parameters
    ----------
    indexed
        Result of :func:`~xrdkit.indexing.index_peaks`, or the first element of
        :func:`~xrdkit.indexing.index_and_refine`.
    wavelength
        Radiation wavelength in angstroms.
    start_cell
        Cell to start from.
    fit_zero
        Refine the zero point error. Held at zero when ``False``.
    fit_displacement
        Refine the specimen displacement. Held at zero when ``False``.
    radius_mm
        Goniometer radius in millimetres. Needed only to refine a displacement.

    Returns
    -------
    LatticeFit
        The refined parameters, their estimated standard deviations, and the
        residuals of the peaks used.

    Raises
    ------
    ValueError
        If ``fit_displacement`` is asked for without a ``radius_mm``, or fewer
        than three peaks more than there are parameters are available.
    """
    if fit_displacement and radius_mm is None:
        raise ValueError(
            "refine_lattice needs radius_mm to refine a specimen displacement"
        )

    used = [
        entry for entry in indexed if entry.is_indexed and len(entry.candidates) == 1
    ]

    free = np.array([True, True, bool(fit_zero), bool(fit_displacement)])
    n_parameters = int(free.sum())
    if len(used) < n_parameters + MIN_EXTRA_PEAKS:
        raise ValueError(
            f"Need at least {n_parameters + MIN_EXTRA_PEAKS} peaks with a single "
            f"candidate to refine {n_parameters} parameters, got {len(used)}"
        )

    observed = np.array([entry.peak.two_theta for entry in used], dtype=float)
    hkl = np.array([entry.reflection.hkl for entry in used], dtype=float)
    hkl_tuples = [entry.reflection.hkl for entry in used]

    # The whole vector is carried about so the model stays readable; only the
    # free entries are handed to the optimiser.
    start = np.array([start_cell.a, start_cell.c, 0.0, 0.0], dtype=float)
    # A displacement term of zero leaves the model unchanged, so the radius is
    # only consulted when one is being refined.
    radius = radius_mm if fit_displacement else None

    def residual(varying: np.ndarray) -> np.ndarray:
        parameters = start.copy()
        parameters[free] = varying
        return observed - _model(parameters, hkl, wavelength, radius)

    result = optimize.least_squares(residual, start[free])

    parameters = start.copy()
    parameters[free] = result.x
    a, c, zero, displacement = parameters

    residuals = np.asarray(result.fun, dtype=float)
    degrees_of_freedom = len(used) - n_parameters
    # Scaling by the reduced chi squared turns the curvature of the surface
    # into an error bar that reflects how well the model actually fits.
    reduced_chi_squared = float(np.sum(residuals**2) / degrees_of_freedom)
    jacobian = np.asarray(result.jac, dtype=float)
    covariance = np.linalg.pinv(jacobian.T @ jacobian) * reduced_chi_squared
    deviations = np.sqrt(np.abs(np.diag(covariance)))

    esd = dict(zip(np.array(PARAMETER_NAMES)[free], deviations))

    return LatticeFit(
        a=float(a),
        c=float(c),
        esd_a=float(esd["a"]),
        esd_c=float(esd["c"]),
        zero=float(zero),
        esd_zero=float(esd["zero"]) if fit_zero else None,
        displacement=float(displacement) if fit_displacement else None,
        esd_displacement=float(esd["displacement"]) if fit_displacement else None,
        radius_mm=radius_mm,
        n_peaks=len(used),
        rms_two_theta=float(np.sqrt(np.mean(residuals**2))),
        residuals=residuals,
        hkl=hkl_tuples,
        converged=bool(result.success),
    )


def lattice_fit_to_dict(fit: LatticeFit) -> dict[str, float | int | bool | None]:
    """Flatten ``fit`` to the scalar fields, ready for a CSV row.

    The residuals and the hkl are left out, being one row each rather than one
    value.
    """
    return {
        "a": fit.a,
        "esd_a": fit.esd_a,
        "c": fit.c,
        "esd_c": fit.esd_c,
        "c_over_a": fit.c / fit.a,
        "zero": fit.zero,
        "esd_zero": fit.esd_zero,
        "displacement": fit.displacement,
        "esd_displacement": fit.esd_displacement,
        "radius_mm": fit.radius_mm,
        "n_peaks": fit.n_peaks,
        "rms_two_theta": fit.rms_two_theta,
        "converged": fit.converged,
    }
