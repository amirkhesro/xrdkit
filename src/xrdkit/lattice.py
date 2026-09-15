"""Lattice parameter refinement against indexed peaks, with systematic errors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize

from xrdkit.cell import Cell
from xrdkit.density import volume_gradient
from xrdkit.indexing import IndexedPeak

__all__ = ["LatticeFit", "lattice_fit_to_dict", "refine_lattice"]

# A peak needs this many more than the parameters before a fit means anything.
MIN_EXTRA_PEAKS = 3

# The residual, in degrees, of every peak for a trial cell that Cell rejects or
# that cannot diffract a reflection: far larger than any real misfit, so the
# optimiser backs out of the step, and finite, so it can.
REJECTED_RESIDUAL = 180.0

# The systematic errors refined after the free cell parameters.
ERROR_NAMES = ("zero", "displacement")

# Crystal systems whose c over a is written out.
C_OVER_A_SYSTEMS = ("tetragonal", "hexagonal", "trigonal")


@dataclass
class LatticeFit:
    """The outcome of a lattice parameter refinement.

    All six cell parameters are given, in angstroms and degrees; one the
    crystal system fixes has the value the cell implies. ``zero`` is the zero
    point error in degrees, added to every calculated position.
    ``displacement`` is the specimen displacement in millimetres, positive
    below the focusing circle, or ``None`` when it was not refined. An
    estimated standard deviation of a parameter that was not refined, whether
    held or fixed by the crystal system, is ``None``.

    ``covariance`` is the covariance matrix of the refined parameters, scaled
    by the reduced chi squared: the free cell parameters in the order of
    ``parameter_names``, then the zero and the displacement where refined.
    ``esd_volume`` is propagated through its cell block, correlations
    included.
    """

    cell: Cell
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    esd_a: float | None
    esd_b: float | None
    esd_c: float | None
    esd_alpha: float | None
    esd_beta: float | None
    esd_gamma: float | None
    volume: float
    esd_volume: float
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
    parameter_names: tuple[str, ...]
    covariance: np.ndarray


def _model(
    cell: Cell,
    zero: float,
    displacement: float,
    hkl: np.ndarray,
    wavelength: float,
    radius_mm: float | None,
) -> np.ndarray | None:
    """Calculated 2theta of ``hkl`` for ``cell`` and the systematic errors, or
    ``None`` if a reflection cannot diffract.

    The specimen displacement shifts a reflection by -2 s cos(theta) / R
    radians, which is the usual flat plate term: a specimen sitting below the
    focusing circle moves every peak to low angle, most of all at low angle.
    """
    sin_theta = wavelength / (2.0 * cell.d_spacings(hkl))
    if np.any(sin_theta > 1.0):
        return None
    theta = np.arcsin(sin_theta)
    shift = zero
    if radius_mm is not None:
        shift = shift + np.degrees(-2.0 * displacement * np.cos(theta) / radius_mm)
    return 2.0 * np.degrees(theta) + shift


def refine_lattice(
    indexed: list[IndexedPeak],
    wavelength: float,
    start_cell: Cell,
    fit_zero: bool = True,
    fit_displacement: bool = False,
    radius_mm: float | None = None,
    start_zero: float = 0.0,
) -> LatticeFit:
    """Refine a cell of any crystal system, and its systematic errors, by
    least squares.

    The free parameters of ``start_cell``'s crystal system are refined, in the
    order of :attr:`~xrdkit.cell.Cell.parameter_names`. Each peak is modelled
    as its Bragg position for the cell plus a zero point error plus a specimen
    displacement term, and the observed positions are taken as measured:
    ``peak.two_theta``, not ``corrected_two_theta``, since the zero point is
    what is being refined rather than assumed.

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
        Cell to start from, which sets the crystal system.
    fit_zero
        Refine the zero point error. When ``False`` it is held at
        ``start_zero``, which is the thing to do when the instrument zero is
        known from a standard.
    fit_displacement
        Refine the specimen displacement. Held at zero when ``False``.
    radius_mm
        Goniometer radius in millimetres. Needed only to refine a displacement.
    start_zero
        Zero point error to start from, in degrees: refined from there when
        ``fit_zero``, and held there when not.

    Returns
    -------
    LatticeFit
        The refined parameters, their estimated standard deviations and
        covariance, the cell volume with its esd, and the residuals of the
        peaks used.

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

    system = start_cell.crystal_system
    parameter_names = start_cell.parameter_names
    n_cell = len(parameter_names)
    names = np.array(parameter_names + ERROR_NAMES)
    free = np.array([True] * n_cell + [bool(fit_zero), bool(fit_displacement)])
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
    start = np.array([*start_cell.parameters.values(), start_zero, 0.0], dtype=float)
    # A displacement term of zero leaves the model unchanged, so the radius is
    # only consulted when one is being refined.
    radius = radius_mm if fit_displacement else None

    def unpack(parameters: np.ndarray) -> tuple[Cell, float, float]:
        cell = Cell.from_parameters(
            system, dict(zip(parameter_names, map(float, parameters[:n_cell])))
        )
        return cell, float(parameters[n_cell]), float(parameters[n_cell + 1])

    def residual(varying: np.ndarray) -> np.ndarray:
        parameters = start.copy()
        parameters[free] = varying
        try:
            cell, zero, displacement = unpack(parameters)
        except ValueError:
            return np.full(len(used), REJECTED_RESIDUAL)
        calculated = _model(cell, zero, displacement, hkl, wavelength, radius)
        if calculated is None:
            return np.full(len(used), REJECTED_RESIDUAL)
        return observed - calculated

    result = optimize.least_squares(residual, start[free])

    parameters = start.copy()
    parameters[free] = result.x
    cell, zero, displacement = unpack(parameters)

    residuals = np.asarray(result.fun, dtype=float)
    degrees_of_freedom = len(used) - n_parameters
    # Scaling by the reduced chi squared turns the curvature of the surface
    # into an error bar that reflects how well the model actually fits.
    reduced_chi_squared = float(np.sum(residuals**2) / degrees_of_freedom)
    jacobian = np.asarray(result.jac, dtype=float)
    covariance = np.linalg.pinv(jacobian.T @ jacobian) * reduced_chi_squared
    deviations = np.sqrt(np.abs(np.diag(covariance)))
    esd = {str(name): float(value) for name, value in zip(names[free], deviations)}

    # The cell parameters lead the free vector, so their block is the corner.
    gradient = volume_gradient(cell)
    cell_covariance = covariance[:n_cell, :n_cell]
    esd_volume = float(np.sqrt(max(float(gradient @ cell_covariance @ gradient), 0.0)))

    return LatticeFit(
        cell=cell,
        a=cell.a,
        b=cell.b,
        c=cell.c,
        alpha=cell.alpha,
        beta=cell.beta,
        gamma=cell.gamma,
        esd_a=esd.get("a"),
        esd_b=esd.get("b"),
        esd_c=esd.get("c"),
        esd_alpha=esd.get("alpha"),
        esd_beta=esd.get("beta"),
        esd_gamma=esd.get("gamma"),
        volume=cell.volume,
        esd_volume=esd_volume,
        zero=zero,
        esd_zero=esd.get("zero"),
        displacement=displacement if fit_displacement else None,
        esd_displacement=esd.get("displacement"),
        radius_mm=radius_mm,
        n_peaks=len(used),
        rms_two_theta=float(np.sqrt(np.mean(residuals**2))),
        residuals=residuals,
        hkl=hkl_tuples,
        converged=bool(result.success),
        parameter_names=parameter_names,
        covariance=covariance,
    )


def lattice_fit_to_dict(fit: LatticeFit) -> dict[str, float | int | bool | None]:
    """Flatten ``fit`` to the scalar fields, ready for a CSV row.

    All six cell parameters with their esds, the volume and its esd, and
    ``c_over_a`` for a tetragonal, hexagonal or trigonal cell only. The
    residuals, the hkl and the covariance are left out, being more than one
    value each.
    """
    row: dict[str, float | int | bool | None] = {}
    for name in ("a", "b", "c", "alpha", "beta", "gamma"):
        row[name] = getattr(fit, name)
        row[f"esd_{name}"] = getattr(fit, f"esd_{name}")
    if fit.cell.crystal_system in C_OVER_A_SYSTEMS:
        row["c_over_a"] = fit.c / fit.a
    row.update(
        {
            "volume": fit.volume,
            "esd_volume": fit.esd_volume,
            "zero": fit.zero,
            "esd_zero": fit.esd_zero,
            "displacement": fit.displacement,
            "esd_displacement": fit.esd_displacement,
            "radius_mm": fit.radius_mm,
            "n_peaks": fit.n_peaks,
            "rms_two_theta": fit.rms_two_theta,
            "converged": fit.converged,
        }
    )
    return row
