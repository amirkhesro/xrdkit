"""Peak widths: profile fitting of reflections and the instrumental resolution.

The widths reported by :func:`~xrdkit.peaks.find_peaks` are read straight off
the counts at half the prominence. They carry no error bar, and they include
whatever of the K alpha 2 satellite lies inside the half maximum, which is all
of it at low angle and a changing fraction of it once the pair starts to
separate. That is good enough to find peaks, but not to measure broadening.

:func:`fit_profile` fits a reflection instead as a K alpha 1 and K alpha 2
doublet of split pseudo-Voigt lines on a linear background, and reports the
width of the K alpha 1 line alone with its esd. The split lets the low angle
half be wider than the high angle half, which is what axial divergence does to
the low angle reflections; a symmetric line misses their tops and so misreads
their widths.

:func:`fit_caglioti` fits the Caglioti relation to the widths from a standard,
which gives the instrumental FWHM at any angle. Sample widths must be measured
with the same :func:`fit_profile` before the instrumental width is taken out of
them, so that the satellite is handled the same way on both sides.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize

from xrdkit.peaks import KALPHA2_RATIO

__all__ = [
    "Caglioti",
    "ProfileFit",
    "fit_caglioti",
    "fit_profile",
    "pseudo_voigt",
    "split_pseudo_voigt",
]

# Integrated intensity of K alpha 2 over K alpha 1, as the Aeris records it.
KALPHA2_INTENSITY_RATIO = 0.5

# Half the fitting window beyond the doublet, as a multiple of the starting
# width, so that the background either side is well sampled.
WINDOW_FWHM = 10.0

# The data must run at least this many fitted widths below K alpha 1 and beyond
# K alpha 2 for the fit to count as whole; less and the tails, and with them the
# background, are not seen.
MIN_MARGIN_FWHM = 3.0

# Each end of the window, as a fraction of its points, gives the starting
# background.
BACKGROUND_FRACTION = 0.1

# Parameters of a split pseudo-Voigt doublet on a linear background, in the
# order fit_profile hands them to the optimiser.
PROFILE_PARAMETERS = (
    "two_theta",
    "fwhm",
    "eta",
    "asymmetry",
    "area",
    "background",
    "slope",
)

# Bounds on the asymmetry, short of the +-1 at which one half width vanishes.
MAX_ASYMMETRY = 0.95

# The Caglioti relation has three parameters; one more point is the least that
# leaves a residual to scale the esds by.
CAGLIOTI_PARAMETERS = 3
MIN_CAGLIOTI_POINTS = CAGLIOTI_PARAMETERS + 1

FOUR_LN2 = 4.0 * np.log(2.0)


def pseudo_voigt(
    two_theta: np.ndarray, centre: float, fwhm: float, eta: float
) -> np.ndarray:
    """A pseudo-Voigt line of unit area.

    The Lorentzian and Gaussian components share one FWHM and are mixed as
    ``eta`` L + (1 - ``eta``) G, so ``eta`` runs from 0 for a pure Gaussian to
    1 for a pure Lorentzian.
    """
    x = (np.asarray(two_theta, dtype=float) - centre) / fwhm
    gaussian = np.sqrt(FOUR_LN2 / np.pi) / fwhm * np.exp(-FOUR_LN2 * x**2)
    lorentzian = 2.0 / (np.pi * fwhm) / (1.0 + 4.0 * x**2)
    return eta * lorentzian + (1.0 - eta) * gaussian


def split_pseudo_voigt(
    two_theta: np.ndarray, centre: float, fwhm: float, eta: float, asymmetry: float
) -> np.ndarray:
    """A pseudo-Voigt line of unit area with different half widths either side.

    The half width below ``centre`` is ``fwhm`` (1 - ``asymmetry``) / 2 and the
    half width above it ``fwhm`` (1 + ``asymmetry``) / 2, so ``fwhm`` is still
    the full width at half maximum and a negative ``asymmetry`` gives the low
    angle tail that axial divergence puts on low angle reflections. The two
    halves meet at the same height, and zero asymmetry is :func:`pseudo_voigt`.
    """
    two_theta = np.asarray(two_theta, dtype=float)
    below = fwhm * (1.0 - asymmetry)
    above = fwhm * (1.0 + asymmetry)
    # Each half is half of a symmetric line of its own width, scaled so the
    # halves join at the peak and the whole integrates to one.
    return np.where(
        two_theta < centre,
        below / fwhm * pseudo_voigt(two_theta, centre, below, eta),
        above / fwhm * pseudo_voigt(two_theta, centre, above, eta),
    )


def _satellite(two_theta: float | np.ndarray, wavelength_ratio: float):
    """The K alpha 2 position of a K alpha 1 line at ``two_theta``, in degrees."""
    sin_theta = wavelength_ratio * np.sin(np.radians(np.asarray(two_theta) / 2.0))
    return 2.0 * np.degrees(np.arcsin(np.clip(sin_theta, -1.0, 1.0)))


@dataclass
class ProfileFit:
    """A reflection fitted as a split pseudo-Voigt K alpha doublet.

    ``two_theta`` and ``fwhm`` belong to the K alpha 1 line; the satellite sits
    at ``kalpha2_two_theta`` with the same width and shape. ``area`` is the
    integrated intensity of K alpha 1 above the background, in counts times
    degrees. Esds are scaled by the reduced chi squared; the esd of a parameter
    held fixed is ``None``.
    """

    two_theta: float
    esd_two_theta: float
    fwhm: float
    esd_fwhm: float
    eta: float
    esd_eta: float
    asymmetry: float
    esd_asymmetry: float | None
    area: float
    esd_area: float
    background: float
    slope: float
    kalpha2_two_theta: float
    window: tuple[float, float]
    n_points: int
    reduced_chi_squared: float
    # Weighted profile R factor, as a fraction.
    r_wp: float
    converged: bool
    # Whether the data stop short of MIN_MARGIN_FWHM widths beyond the doublet
    # at either end, so that a tail and the background there go unseen.
    truncated: bool
    # The fitted curve over the window, for plotting and inspection.
    fitted: np.ndarray


def fit_profile(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    centre: float,
    fwhm_guess: float,
    window: tuple[float, float] | None = None,
    fit_asymmetry: bool = True,
    wavelength_ratio: float = KALPHA2_RATIO,
    intensity_ratio: float = KALPHA2_INTENSITY_RATIO,
) -> ProfileFit:
    """Fit one reflection as a K alpha 1 and K alpha 2 split pseudo-Voigt doublet.

    The satellite is tied to its parent: it sits where the parent's d spacing
    puts the longer wavelength, carries ``intensity_ratio`` of its area, and
    shares its width, mixing parameter and asymmetry. Both lines stand on one
    linear background. The free parameters are therefore the K alpha 1
    position, FWHM, mixing parameter, asymmetry and area, and the background
    level and slope. Counts are weighted as Poisson, 1 / max(counts, 1).

    The same function, with the same settings, has to measure the standard and
    the samples, so that an instrumental width is only ever compared with a
    sample width that treats the doublet and the asymmetry the same way.

    Parameters
    ----------
    two_theta, intensity
        The pattern, in degrees and counts.
    centre
        Starting K alpha 1 position, in degrees.
    fwhm_guess
        Starting FWHM, in degrees. A width read off the raw doublet, such as
        ``Peak.fwhm``, is close enough.
    window
        ``(low, high)`` range to fit, in degrees. By default it runs
        ``WINDOW_FWHM`` starting widths below ``centre`` and the same beyond
        the satellite, cut at the ends of the data.
    fit_asymmetry
        Refine the asymmetry of :func:`split_pseudo_voigt`. When ``False`` it
        is held at zero and the lines are symmetric pseudo-Voigts.
    wavelength_ratio
        K alpha 2 over K alpha 1 wavelength.
    intensity_ratio
        K alpha 2 over K alpha 1 integrated intensity. Zero fits a single line.

    Returns
    -------
    ProfileFit
        The fitted parameters, their esds and the quality of the fit.

    Raises
    ------
    ValueError
        If ``fwhm_guess`` is not positive, or the window holds no more points
        than there are parameters.
    """
    if fwhm_guess <= 0.0:
        raise ValueError(f"fwhm_guess must be positive, got {fwhm_guess}")

    two_theta = np.asarray(two_theta, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if window is None:
        margin = WINDOW_FWHM * fwhm_guess
        window = (
            centre - margin,
            float(_satellite(centre, wavelength_ratio)) + margin,
        )
    low, high = window
    selected = (two_theta >= low) & (two_theta <= high)
    x = two_theta[selected]
    y = intensity[selected]

    free = np.ones(len(PROFILE_PARAMETERS), dtype=bool)
    free[PROFILE_PARAMETERS.index("asymmetry")] = fit_asymmetry
    n_parameters = int(free.sum())
    if x.size <= n_parameters:
        raise ValueError(
            f"window {low:.3f}-{high:.3f} holds {x.size} points, "
            f"need more than {n_parameters}"
        )
    middle = 0.5 * (x[0] + x[-1])
    sigma = np.sqrt(np.maximum(y, 1.0))

    def model(parameters: np.ndarray) -> np.ndarray:
        position, fwhm, eta, asymmetry, area, background, slope = parameters
        satellite = _satellite(position, wavelength_ratio)
        peak = split_pseudo_voigt(x, position, fwhm, eta, asymmetry)
        peak += intensity_ratio * split_pseudo_voigt(x, satellite, fwhm, eta, asymmetry)
        return background + slope * (x - middle) + area * peak

    # Starting background: a straight line through the two ends of the window.
    edge = max(2, int(BACKGROUND_FRACTION * x.size))
    left = float(np.mean(y[:edge]))
    right = float(np.mean(y[-edge:]))
    slope0 = (right - left) / (np.mean(x[-edge:]) - np.mean(x[:edge]))
    background0 = 0.5 * (left + right)
    height = max(float(np.max(y)) - background0, 1.0)
    # A pseudo-Voigt of eta one half is about 1.3 times its height by its FWHM.
    area0 = 1.3 * height * fwhm_guess / (1.0 + intensity_ratio)

    step = float(np.median(np.diff(x)))
    lower = np.array([x[0], step / 2.0, 0.0, -MAX_ASYMMETRY, 0.0, -np.inf, -np.inf])
    upper = np.array([x[-1], x[-1] - x[0], 1.0, MAX_ASYMMETRY, np.inf, np.inf, np.inf])
    start = np.array([centre, fwhm_guess, 0.5, 0.0, area0, background0, slope0])
    start = np.clip(start, lower, upper)

    def residual(varying: np.ndarray) -> np.ndarray:
        parameters = start.copy()
        parameters[free] = varying
        return (y - model(parameters)) / sigma

    result = optimize.least_squares(
        residual, start[free], bounds=(lower[free], upper[free]), x_scale="jac"
    )

    parameters = start.copy()
    parameters[free] = result.x
    residuals = np.asarray(result.fun, dtype=float)
    reduced_chi_squared = float(np.sum(residuals**2) / (x.size - n_parameters))
    jacobian = np.asarray(result.jac, dtype=float)
    covariance = np.linalg.pinv(jacobian.T @ jacobian) * reduced_chi_squared
    deviations = np.sqrt(np.abs(np.diag(covariance)))
    esd = {
        name: float(value)
        for name, value in zip(np.array(PROFILE_PARAMETERS)[free], deviations)
    }
    position, fwhm, eta, asymmetry, area, background, slope = (
        float(value) for value in parameters
    )
    satellite = float(_satellite(position, wavelength_ratio))
    margin = MIN_MARGIN_FWHM * fwhm
    r_wp = float(np.sqrt(np.sum(residuals**2) / np.sum((y / sigma) ** 2)))

    return ProfileFit(
        two_theta=position,
        esd_two_theta=esd["two_theta"],
        fwhm=fwhm,
        esd_fwhm=esd["fwhm"],
        eta=eta,
        esd_eta=esd["eta"],
        asymmetry=asymmetry,
        esd_asymmetry=esd.get("asymmetry"),
        area=area,
        esd_area=esd["area"],
        background=background,
        slope=slope,
        kalpha2_two_theta=satellite,
        window=(float(x[0]), float(x[-1])),
        n_points=int(x.size),
        reduced_chi_squared=reduced_chi_squared,
        r_wp=r_wp,
        converged=bool(result.success),
        truncated=bool(position - x[0] < margin or x[-1] - satellite < margin),
        fitted=model(parameters),
    )


def _tan_theta(two_theta: float | np.ndarray) -> np.ndarray:
    """tan(theta) of ``two_theta`` in degrees, which must lie in (0, 180)."""
    two_theta = np.asarray(two_theta, dtype=float)
    if np.any((two_theta <= 0.0) | (two_theta >= 180.0)):
        raise ValueError("two_theta must lie strictly between 0 and 180 degrees")
    return np.tan(np.radians(two_theta / 2.0))


@dataclass
class Caglioti:
    """The Caglioti resolution function, FWHM^2 = U tan^2 + V tan + W.

    U, V and W are in degrees squared of 2theta, theta being half the
    diffraction angle. ``covariance`` is the 3 by 3 covariance of (U, V, W),
    scaled by the reduced chi squared of the fit; ``rms`` is the root mean
    square of the FWHM residuals, in degrees.
    """

    u: float
    v: float
    w: float
    esd_u: float
    esd_v: float
    esd_w: float
    n_peaks: int
    rms: float
    covariance: np.ndarray

    def _squared(self, two_theta: float | np.ndarray) -> np.ndarray:
        tan_theta = _tan_theta(two_theta)
        return self.u * tan_theta**2 + self.v * tan_theta + self.w

    def fwhm(self, two_theta: float | np.ndarray) -> float | np.ndarray:
        """The FWHM at ``two_theta``, in degrees.

        Where the quadratic dips below zero, which a fitted V < 0 can do outside
        the range of the data, the width is taken as zero rather than the
        square root of a negative number.

        Raises
        ------
        ValueError
            If any ``two_theta`` lies outside (0, 180) degrees.
        """
        width = np.sqrt(np.clip(self._squared(two_theta), 0.0, None))
        return float(width) if width.ndim == 0 else width

    def fwhm_esd(self, two_theta: float | np.ndarray) -> float | np.ndarray:
        """The esd of :meth:`fwhm` at ``two_theta``, from the covariance.

        ``nan`` where the width is zero, since the esd of a square root is not
        defined there.
        """
        tan_theta = _tan_theta(two_theta)
        gradient = np.stack([tan_theta**2, tan_theta, np.ones_like(tan_theta)])
        variance = np.einsum("i...,ij,j...->...", gradient, self.covariance, gradient)
        width = np.asarray(self.fwhm(two_theta), dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            esd = np.where(
                width > 0.0, np.sqrt(np.abs(variance)) / (2.0 * width), np.nan
            )
        return float(esd) if esd.ndim == 0 else esd


def fit_caglioti(
    two_theta: np.ndarray,
    fwhm: np.ndarray,
    weights: np.ndarray | None = None,
) -> Caglioti:
    """Fit FWHM^2 = U tan^2(theta) + V tan(theta) + W by weighted least squares.

    The relation is linear in U, V and W, so this is a single linear solve.
    The covariance is scaled by the reduced chi squared, as
    :func:`~xrdkit.lattice.refine_lattice` scales its own, so the esds reflect
    how well the relation actually describes the widths.

    Parameters
    ----------
    two_theta
        Peak positions, in degrees.
    fwhm
        Peak widths, in degrees.
    weights
        Weight of each point in the sum of squared FWHM^2 residuals. For a
        width with esd s this is 1 / (2 FWHM s)^2. Equal weights by default.

    Raises
    ------
    ValueError
        If the arrays differ in length, there are fewer than four points, a
        position lies outside (0, 180) degrees, a width is not positive, or a
        weight is negative or not finite.
    """
    two_theta = np.asarray(two_theta, dtype=float).ravel()
    fwhm = np.asarray(fwhm, dtype=float).ravel()
    weights = (
        np.ones_like(fwhm) if weights is None else np.asarray(weights, float).ravel()
    )
    if not two_theta.size == fwhm.size == weights.size:
        raise ValueError(
            f"two_theta, fwhm and weights differ in length: {two_theta.size}, "
            f"{fwhm.size}, {weights.size}"
        )
    if two_theta.size < MIN_CAGLIOTI_POINTS:
        raise ValueError(
            f"Need at least {MIN_CAGLIOTI_POINTS} widths to fit U, V and W, "
            f"got {two_theta.size}"
        )
    if np.any(~np.isfinite(fwhm) | (fwhm <= 0.0)):
        raise ValueError("Every fwhm must be positive")
    if np.any(~np.isfinite(weights) | (weights < 0.0)):
        raise ValueError("Every weight must be finite and not negative")

    tan_theta = _tan_theta(two_theta)
    design = np.column_stack([tan_theta**2, tan_theta, np.ones_like(tan_theta)])
    root = np.sqrt(weights)
    parameters, *_ = np.linalg.lstsq(design * root[:, None], fwhm**2 * root, rcond=None)

    weighted = (fwhm**2 - design @ parameters) * root
    degrees_of_freedom = two_theta.size - CAGLIOTI_PARAMETERS
    reduced_chi_squared = float(np.sum(weighted**2) / degrees_of_freedom)
    normal = design.T @ (design * weights[:, None])
    covariance = np.linalg.pinv(normal) * reduced_chi_squared
    esd = np.sqrt(np.abs(np.diag(covariance)))

    u, v, w = (float(p) for p in parameters)
    fit = Caglioti(
        u=u,
        v=v,
        w=w,
        esd_u=float(esd[0]),
        esd_v=float(esd[1]),
        esd_w=float(esd[2]),
        n_peaks=int(two_theta.size),
        rms=0.0,
        covariance=covariance,
    )
    fit.rms = float(np.sqrt(np.mean((fwhm - fit.fwhm(two_theta)) ** 2)))
    return fit
