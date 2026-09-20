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

:func:`fit_breadth_models` asks what the corrected breadths follow with angle:
crystallite size, strain, or a spread of specimen heights across the surface,
:func:`height_spread_breadth`, which the LaB6 standard would not share.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special

from xrdkit.peaks import KALPHA2_RATIO

__all__ = [
    "BreadthModelFit",
    "BreadthModels",
    "BroadeningCorrection",
    "Caglioti",
    "ProfileFit",
    "correct_broadening",
    "doublet_gaps",
    "fit_breadth_models",
    "fit_caglioti",
    "fit_profile",
    "height_spread_breadth",
    "integral_breadth",
    "kalpha2_position",
    "pseudo_voigt",
    "pseudo_voigt_components",
    "pseudo_voigt_from_components",
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

# Thompson, Cox and Hastings, J. Appl. Cryst. 20 (1987) 79: the FWHM of a
# Voigt as the fifth root of a quintic in its Gaussian and Lorentzian FWHM,
# coefficients of G^5, G^4 L, ... L^5, and the pseudo-Voigt mixing parameter
# that matches it as a cubic in L / FWHM, coefficients of the first to third
# powers.
TCH_FWHM = (1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0)
TCH_ETA = (1.36603, -0.47719, 0.11116)

# A sample width must exceed the instrumental width by this many combined
# esds to count as resolved.
DEFAULT_SIGNIFICANCE = 2.0

# Relative step of the numerical derivatives used to propagate esds through
# the correction.
DERIVATIVE_STEP = 1e-6

# fit_breadth_models: the fewest breadths that leave the two parameter model a
# residual, and the Scherrer constant of an integral breadth.
MIN_BREADTH_MODEL_POINTS = 3
BREADTH_MODEL_K = 1.0

# The squared terms of the size plus height fit start at least this far above
# zero, as a fraction of the largest squared breadth. A free fit that beats the
# better one term fit by less than QUADRATURE_BOUND_CHI_SQUARED in chi squared,
# far below the one or so a real improvement needs, has only crept towards a
# bound from inside, where the esd of the vanishing term runs away; the fit is
# then taken as on that bound.
QUADRATURE_START_FLOOR = 1e-6
QUADRATURE_BOUND_CHI_SQUARED = 1e-3


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


def kalpha2_position(
    two_theta: float, wavelength_ratio: float = KALPHA2_RATIO
) -> float:
    """Where the K alpha 2 line of a K alpha 1 line at ``two_theta`` falls.

    The satellite diffracts at the same d spacing at the longer wavelength, so
    it always lies to high angle, by more the higher the angle.
    """
    return float(_satellite(two_theta, wavelength_ratio))


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
        Weight of each point in the sum of squared FWHM^2 residuals, used as
        given: the residual of every point is scaled by the root of its
        weight, and the reduced chi squared the covariance is scaled by
        follows from the same weights. Equal weights by default.

        Which weights to pass is the caller's choice. The kit's own caller,
        :func:`~xrdkit.instrument.fit_instrument_widths`, passes 1 / s^2 with
        s the esd of the fitted FWHM, which weights each width by its own
        precision. A caller wanting the esd of a width carried through to the
        square of that width should pass 1 / (2 FWHM s)^2 instead, since the
        variance of FWHM^2 is (2 FWHM s)^2. The two are not the same
        weighting and give slightly different U, V and W.

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


def doublet_gaps(
    two_theta: float,
    others: list[float] | np.ndarray,
    wavelength_ratio: float = KALPHA2_RATIO,
) -> tuple[float, float]:
    """The clear space either side of a K alpha doublet, in degrees.

    ``others`` are the K alpha 1 positions of every other reflection or peak
    that might lie nearby; each brings its own K alpha 2 line too. The first
    gap runs down from the K alpha 1 line at ``two_theta`` to the nearest line
    below it, the second up from its K alpha 2 line to the nearest line above.
    A line falling between the two members of the doublet closes both gaps to
    zero, and a side with no line at all is infinitely clear.
    """
    satellite = float(_satellite(two_theta, wavelength_ratio))
    others = np.asarray(others, dtype=float).ravel()
    lines = np.concatenate([others, _satellite(others, wavelength_ratio)])
    if np.any((lines >= two_theta) & (lines <= satellite)):
        return 0.0, 0.0
    below = lines[lines < two_theta]
    above = lines[lines > satellite]
    return (
        float(two_theta - below.max()) if below.size else float("inf"),
        float(above.min() - satellite) if above.size else float("inf"),
    )


def pseudo_voigt_from_components(
    fwhm_gaussian: float, fwhm_lorentzian: float
) -> tuple[float, float]:
    """The FWHM and mixing parameter of the pseudo-Voigt matching a Voigt.

    Uses the Thompson, Cox and Hastings approximation, from the FWHM of the
    Gaussian and Lorentzian the Voigt convolves.

    Raises
    ------
    ValueError
        If either width is negative, or both are zero.
    """
    if fwhm_gaussian < 0.0 or fwhm_lorentzian < 0.0:
        raise ValueError(
            f"Component widths cannot be negative, got {fwhm_gaussian} "
            f"and {fwhm_lorentzian}"
        )
    if fwhm_gaussian == 0.0 and fwhm_lorentzian == 0.0:
        raise ValueError("At least one component width must be positive")
    fwhm = _tch_fwhm(fwhm_gaussian, fwhm_lorentzian)
    # The cubic reaches 1 at a pure Lorentzian only to rounding.
    return fwhm, min(_tch_eta(fwhm_lorentzian / fwhm), 1.0)


def pseudo_voigt_components(fwhm: float, eta: float) -> tuple[float, float]:
    """The Gaussian and Lorentzian FWHM of the Voigt a pseudo-Voigt matches.

    The inverse of :func:`pseudo_voigt_from_components`: the Thompson, Cox and
    Hastings cubic is solved for the Lorentzian fraction of the width, and the
    quintic then for the Gaussian width that makes up the rest.

    Returns
    -------
    tuple[float, float]
        The Gaussian and Lorentzian FWHM, in the units of ``fwhm``.

    Raises
    ------
    ValueError
        If ``fwhm`` is not positive or ``eta`` lies outside 0 to 1.
    """
    if not fwhm > 0.0:
        raise ValueError(f"fwhm must be positive, got {fwhm}")
    if not 0.0 <= eta <= 1.0:
        raise ValueError(f"eta must lie between 0 and 1, got {eta}")

    # The cubic rises steadily from 0 to 1 across the unit interval, so it has
    # exactly one root there; the ends are settled directly so that rounding in
    # the coefficients cannot push the root out of the bracket.
    if eta <= 0.0:
        fraction = 0.0
    elif eta >= _tch_eta(1.0):
        fraction = 1.0
    else:
        fraction = optimize.brentq(
            lambda ratio: _tch_eta(ratio) - eta, 0.0, 1.0, xtol=1e-15
        )
    lorentzian = fraction * fwhm
    if fraction >= 1.0:
        return 0.0, fwhm

    # The quintic grows with the Gaussian width, from the Lorentzian width
    # alone at zero up past fwhm, so again there is one root in the bracket.
    # Where the Lorentzian is so small that rounding leaves the top of the
    # bracket short of fwhm, the profile is Gaussian to working precision.
    def shortfall(width: float) -> float:
        return _tch_fwhm(width, lorentzian) - fwhm

    if shortfall(fwhm) <= 0.0:
        return fwhm, lorentzian
    gaussian = optimize.brentq(shortfall, 0.0, fwhm, xtol=1e-15)
    return float(gaussian), float(lorentzian)


def integral_breadth(fwhm_gaussian: float, fwhm_lorentzian: float) -> float:
    """The integral breadth of a Voigt, from its component FWHM.

    Area over height, exact for the Voigt: beta = beta_G / erfcx(k), with
    k = beta_L / (sqrt(pi) beta_G), beta_G = (FWHM_G / 2) sqrt(pi / ln 2) and
    beta_L = (pi / 2) FWHM_L.
    """
    beta_gaussian = 0.5 * fwhm_gaussian * np.sqrt(np.pi / np.log(2.0))
    beta_lorentzian = 0.5 * np.pi * fwhm_lorentzian
    if beta_gaussian == 0.0:
        return float(beta_lorentzian)
    k = beta_lorentzian / (np.sqrt(np.pi) * beta_gaussian)
    return float(beta_gaussian / special.erfcx(k))


def _tch_fwhm(gaussian: float, lorentzian: float) -> float:
    powers = [
        gaussian ** (5 - order) * lorentzian**order for order in range(len(TCH_FWHM))
    ]
    return float(np.dot(TCH_FWHM, powers) ** 0.2)


def _tch_eta(ratio: float) -> float:
    return float(sum(c * ratio ** (order + 1) for order, c in enumerate(TCH_ETA)))


@dataclass
class BroadeningCorrection:
    """The sample breadth of a reflection, with the instrument taken out.

    Widths are in the units given, degrees of 2theta in practice. The profile
    is split into Gaussian and Lorentzian parts by Thompson, Cox and Hastings,
    the instrumental Lorentzian taken off the observed one linearly and the
    Gaussian in quadrature, and the two recombined into ``fwhm`` and ``eta``.
    ``integral_breadth`` is that of the corrected Voigt.

    ``fwhm_linear`` and ``fwhm_quadrature`` are the simple estimates for
    comparison: the observed minus the instrumental FWHM, as if both were
    Lorentzian, and the difference of their squares, as if both were Gaussian.

    A sample component that comes out below zero, which noise can do to the
    Gaussian part of a size broadened peak for instance, is set to zero and
    flagged; its esd is then that of the difference it was clipped from.

    When ``unresolved``, the observed width is within ``significance`` combined
    esds of the instrumental one, so the sample breadth is not measurable and
    every sample quantity is ``None``. ``excess`` is the observed minus the
    instrumental FWHM over their combined esd either way.
    """

    fwhm: float | None
    esd_fwhm: float | None
    eta: float | None
    esd_eta: float | None
    fwhm_gaussian: float | None
    esd_fwhm_gaussian: float | None
    fwhm_lorentzian: float | None
    esd_fwhm_lorentzian: float | None
    integral_breadth: float | None
    esd_integral_breadth: float | None
    fwhm_linear: float | None
    esd_fwhm_linear: float | None
    fwhm_quadrature: float | None
    esd_fwhm_quadrature: float | None
    gaussian_clipped: bool
    lorentzian_clipped: bool
    unresolved: bool
    excess: float


def _quadrature_esd(value: float, esd_squared: float) -> float:
    """The esd of sqrt(D) from the esd of D.

    Linear propagation gives esd(D) / (2 sqrt(D)), which runs away as D goes to
    zero; it is capped at sqrt(esd(D)), the width at which D would stand one
    esd above zero, so a width at or near zero keeps a finite, honest esd.
    """
    if value <= 0.0:
        return float(np.sqrt(esd_squared))
    return float(min(esd_squared / (2.0 * value), np.sqrt(esd_squared)))


def _corrected(inputs: np.ndarray) -> np.ndarray:
    """Sample (FWHM, eta, G, L, beta, raw L difference, raw G^2 difference)."""
    fwhm_obs, eta_obs, fwhm_inst, eta_inst = inputs
    gaussian_obs, lorentzian_obs = pseudo_voigt_components(fwhm_obs, eta_obs)
    gaussian_inst, lorentzian_inst = pseudo_voigt_components(fwhm_inst, eta_inst)
    lorentzian_difference = lorentzian_obs - lorentzian_inst
    squared_difference = gaussian_obs**2 - gaussian_inst**2
    gaussian = np.sqrt(max(squared_difference, 0.0))
    lorentzian = max(lorentzian_difference, 0.0)
    fwhm, eta = pseudo_voigt_from_components(gaussian, lorentzian)
    return np.array(
        [
            fwhm,
            eta,
            gaussian,
            lorentzian,
            integral_breadth(gaussian, lorentzian),
            lorentzian_difference,
            squared_difference,
        ]
    )


def _jacobian(function, inputs: np.ndarray) -> np.ndarray:
    """Numerical derivatives of ``function`` at ``inputs``, one column each.

    Central differences, except that a mixing parameter at 0 or 1 is stepped
    inwards only, since the components are not defined beyond it.
    """
    centre = function(inputs)
    columns = []
    for index, value in enumerate(inputs):
        step = DERIVATIVE_STEP * max(abs(value), 1.0e-3)
        is_eta = index in (1, 3)
        up = inputs.copy()
        down = inputs.copy()
        if is_eta and value + step > 1.0:
            down[index] -= step
            columns.append((centre - function(down)) / step)
        elif is_eta and value - step < 0.0:
            up[index] += step
            columns.append((function(up) - centre) / step)
        else:
            up[index] += step
            down[index] -= step
            columns.append((function(up) - function(down)) / (2.0 * step))
    return np.column_stack(columns)


def correct_broadening(
    fwhm_obs: float,
    eta_obs: float,
    fwhm_inst: float,
    eta_inst: float,
    esd_fwhm_obs: float = 0.0,
    esd_eta_obs: float = 0.0,
    esd_fwhm_inst: float = 0.0,
    esd_eta_inst: float = 0.0,
    significance: float = DEFAULT_SIGNIFICANCE,
) -> BroadeningCorrection:
    """Take the instrumental broadening out of an observed pseudo-Voigt.

    Both profiles are split into their Gaussian and Lorentzian FWHM with
    :func:`pseudo_voigt_components`. Lorentzians convolve by adding widths and
    Gaussians by adding squares, so the sample Lorentzian is the observed one
    less the instrumental one, and the sample Gaussian the root of the
    difference of squares. The two are recombined with
    :func:`pseudo_voigt_from_components`.

    Esds are propagated linearly through numerical derivatives, treating the
    four inputs as independent. The observed width and mixing parameter of a
    fit are in fact correlated, so the esds are indicative rather than exact.

    Parameters
    ----------
    fwhm_obs, eta_obs
        The observed FWHM and mixing parameter.
    fwhm_inst, eta_inst
        The instrumental FWHM and mixing parameter at the same angle. The
        caller has to supply ``eta_inst`` itself: the Caglioti relation is
        about widths only, so neither :class:`Caglioti` nor
        :class:`~xrdkit.instrument.WidthFit` carries a mixing parameter. It
        comes from the profile fits of the standard the instrument was
        measured on, interpolated to the angle wanted.
    esd_fwhm_obs, esd_eta_obs, esd_fwhm_inst, esd_eta_inst
        Their esds; zero by default.
    significance
        How many combined esds the observed width must clear the instrumental
        one by for the sample breadth to be resolved.

    Raises
    ------
    ValueError
        If a width is not positive, a mixing parameter lies outside 0 to 1, or
        an esd is negative.
    """
    for name, width in (("fwhm_obs", fwhm_obs), ("fwhm_inst", fwhm_inst)):
        if not width > 0.0:
            raise ValueError(f"{name} must be positive, got {width}")
    for name, eta in (("eta_obs", eta_obs), ("eta_inst", eta_inst)):
        if not 0.0 <= eta <= 1.0:
            raise ValueError(f"{name} must lie between 0 and 1, got {eta}")
    esds = np.array([esd_fwhm_obs, esd_eta_obs, esd_fwhm_inst, esd_eta_inst])
    if np.any(~np.isfinite(esds) | (esds < 0.0)):
        raise ValueError("Every esd must be finite and not negative")

    combined = float(np.hypot(esd_fwhm_obs, esd_fwhm_inst))
    difference = fwhm_obs - fwhm_inst
    # With no esds at all, any positive difference counts as resolved.
    excess = (
        difference / combined
        if combined > 0.0
        else float(np.copysign(np.inf, difference))
    )
    if difference <= 0.0 or difference < significance * combined:
        return BroadeningCorrection(
            *([None] * 14),
            gaussian_clipped=False,
            lorentzian_clipped=False,
            unresolved=True,
            excess=float(excess),
        )

    inputs = np.array([fwhm_obs, eta_obs, fwhm_inst, eta_inst], dtype=float)
    values = _corrected(inputs)
    jacobian = _jacobian(_corrected, inputs)
    propagated = np.sqrt((jacobian**2) @ esds**2)
    fwhm, eta, gaussian, lorentzian, breadth, lorentzian_difference, squared = values

    gaussian_clipped = bool(squared < 0.0)
    lorentzian_clipped = bool(lorentzian_difference < 0.0)
    esd_lorentzian = float(propagated[5])
    esd_gaussian = _quadrature_esd(float(gaussian), float(propagated[6]))

    linear = difference
    esd_linear = combined
    quadrature_squared = fwhm_obs**2 - fwhm_inst**2
    quadrature = float(np.sqrt(quadrature_squared))
    esd_quadrature_squared = 2.0 * float(
        np.hypot(fwhm_obs * esd_fwhm_obs, fwhm_inst * esd_fwhm_inst)
    )

    return BroadeningCorrection(
        fwhm=float(fwhm),
        esd_fwhm=float(propagated[0]),
        eta=float(eta),
        esd_eta=float(propagated[1]),
        fwhm_gaussian=float(gaussian),
        esd_fwhm_gaussian=esd_gaussian,
        fwhm_lorentzian=float(lorentzian),
        esd_fwhm_lorentzian=esd_lorentzian,
        integral_breadth=float(breadth),
        esd_integral_breadth=float(propagated[4]),
        fwhm_linear=float(linear),
        esd_fwhm_linear=esd_linear,
        fwhm_quadrature=quadrature,
        esd_fwhm_quadrature=_quadrature_esd(quadrature, esd_quadrature_squared),
        gaussian_clipped=gaussian_clipped,
        lorentzian_clipped=lorentzian_clipped,
        unresolved=False,
        excess=float(excess),
    )


def _cos_theta(two_theta: float | np.ndarray) -> np.ndarray:
    """cos(theta) of ``two_theta`` in degrees, which must lie in (0, 180)."""
    two_theta = np.asarray(two_theta, dtype=float)
    if np.any(~np.isfinite(two_theta) | (two_theta <= 0.0) | (two_theta >= 180.0)):
        raise ValueError("two_theta must lie strictly between 0 and 180 degrees")
    return np.cos(np.radians(two_theta / 2.0))


def height_spread_breadth(
    two_theta: float | np.ndarray, delta_s_mm: float, radius_mm: float
) -> float | np.ndarray:
    """The breadth a spread of specimen heights gives a line, in degrees.

    A specimen displaced by s shifts a line by -2 s cos(theta) / R radians of
    2theta. A surface whose heights spread uniformly over ``delta_s_mm``
    spreads the shift over 2 delta_s cos(theta) / R radians, and that is both
    the FWHM and the integral breadth of the broadening it adds. Unlike size
    or strain broadening it narrows as the angle rises. The result is in
    degrees, like every other width here.

    Parameters
    ----------
    two_theta
        Position, in degrees.
    delta_s_mm
        Full spread of heights across the irradiated surface, in mm.
    radius_mm
        Goniometer radius, in mm.

    Raises
    ------
    ValueError
        If ``delta_s_mm`` is negative, ``radius_mm`` is not positive or a
        position lies outside (0, 180) degrees.
    """
    if not delta_s_mm >= 0.0:
        raise ValueError(f"delta_s_mm cannot be negative, got {delta_s_mm}")
    if not radius_mm > 0.0:
        raise ValueError(f"radius_mm must be positive, got {radius_mm}")
    breadth = np.degrees(2.0 * delta_s_mm * _cos_theta(two_theta) / radius_mm)
    return float(breadth) if breadth.ndim == 0 else breadth


@dataclass
class BreadthModelFit:
    """One model of sample breadth against angle, fitted by weighted least squares.

    The model is sqrt((K lambda / (D cos))^2 + (2 delta_s cos / R)^2)
    + 4 epsilon tan, with whichever of the three terms the model has.
    ``size`` is in the units of the wavelength and ``delta_s`` in mm; each is
    ``None`` when the model lacks it or it fits at zero, a size that fits at
    zero breadth being infinite. Esds are scaled by the reduced chi squared.
    ``at_bound`` marks a two parameter fit with one term held at zero, where
    the fit has in effect fallen back to a one parameter model.
    ``normalised_residuals`` are the observed less the fitted breadths over
    their esds.
    """

    model: str
    size: float | None
    esd_size: float | None
    strain: float | None
    esd_strain: float | None
    delta_s: float | None
    esd_delta_s: float | None
    chi_squared: float
    degrees_of_freedom: int
    reduced_chi_squared: float
    at_bound: bool
    normalised_residuals: np.ndarray
    # The terms in radians: K lambda / D, epsilon and 2 delta_s / R.
    size_term: float
    strain_term: float
    height_term: float

    def breadth(self, two_theta: float | np.ndarray) -> float | np.ndarray:
        """The fitted breadth at ``two_theta``, in degrees."""
        cos_theta = _cos_theta(two_theta)
        tan_theta = np.sqrt(1.0 - cos_theta**2) / cos_theta
        breadth = np.degrees(
            np.hypot(self.size_term / cos_theta, self.height_term * cos_theta)
            + 4.0 * self.strain_term * tan_theta
        )
        return float(breadth) if breadth.ndim == 0 else breadth


@dataclass
class BreadthModels:
    """The four fits of :func:`fit_breadth_models`, by model name."""

    fits: dict[str, BreadthModelFit]

    @property
    def preferred(self) -> BreadthModelFit:
        """The fit with the lowest reduced chi squared."""
        return min(self.fits.values(), key=lambda fit: fit.reduced_chi_squared)


def _breadth_model_fit(
    model: str,
    terms: dict[str, tuple[float, float]],
    residuals: np.ndarray,
    n_parameters: int,
    at_bound: bool,
    k_lambda: float,
    radius_mm: float,
) -> BreadthModelFit:
    """Turn fitted terms in radians into a size, a strain and a height spread.

    ``terms`` maps "size", "strain" and "height" to (term, esd) for the terms
    the model has; ``residuals`` are normalised by the breadth esds.
    """
    chi_squared = float(np.sum(residuals**2))
    dof = residuals.size - n_parameters
    size = esd_size = strain = esd_strain = delta_s = esd_delta_s = None
    size_term, esd_size_term = terms.get("size", (0.0, 0.0))
    strain_term, esd_strain_term = terms.get("strain", (0.0, 0.0))
    height_term, esd_height_term = terms.get("height", (0.0, 0.0))
    if size_term > 0.0:
        size = k_lambda / size_term
        esd_size = size * esd_size_term / size_term
    if "strain" in terms:
        strain, esd_strain = strain_term, esd_strain_term
    if height_term > 0.0:
        delta_s = height_term * radius_mm / 2.0
        esd_delta_s = esd_height_term * radius_mm / 2.0
    return BreadthModelFit(
        model=model,
        size=size,
        esd_size=esd_size,
        strain=strain,
        esd_strain=esd_strain,
        delta_s=delta_s,
        esd_delta_s=esd_delta_s,
        chi_squared=chi_squared,
        degrees_of_freedom=dof,
        reduced_chi_squared=chi_squared / dof,
        at_bound=at_bound,
        normalised_residuals=residuals,
        size_term=max(size_term, 0.0),
        strain_term=strain_term,
        height_term=max(height_term, 0.0),
    )


def fit_breadth_models(
    two_theta: np.ndarray,
    breadth: np.ndarray,
    esd: np.ndarray,
    wavelength: float,
    radius_mm: float,
    k: float = BREADTH_MODEL_K,
) -> BreadthModels:
    """Fit sample breadths with size, strain and height spread models.

    Three one parameter models, each a breadth proportional to its own
    angular dependence and fitted linearly:

    ``size``
        beta = K lambda / (D cos(theta)).
    ``strain``
        beta = 4 epsilon tan(theta).
    ``height``
        beta = 2 delta_s cos(theta) / R, as :func:`height_spread_breadth`.

    and one with two, ``size+height``, the size and height breadths added in
    quadrature. That one is fitted in the squares of the two terms,
    P = (K lambda / D)^2 and Q = (2 delta_s / R)^2, each held at zero or
    above, which keeps the derivatives finite when either term vanishes.
    Every fit weights by the breadth esds and is judged by its reduced chi
    squared on the breadths themselves, so the four compare directly.

    Size and height breadths go as 1 / cos(theta) and cos(theta), both close
    to flat over a narrow range of angle, so the two parameter fit is strongly
    correlated and its esds are large unless the angles span widely.

    Parameters
    ----------
    two_theta
        Positions, in degrees.
    breadth, esd
        Sample breadths with the instrument taken out, and their esds, in
        degrees of 2theta; every esd must be positive.
    wavelength
        Wavelength, which sets the units of the size.
    radius_mm
        Goniometer radius, in mm.
    k
        Scherrer constant, 1 for integral breadths by default and 0.9 for a
        FWHM.

    Raises
    ------
    ValueError
        If the arrays differ in length, there are fewer than three breadths, a
        position lies outside (0, 180) degrees, a breadth is not finite, an
        esd is not positive, or the wavelength, radius or ``k`` is not
        positive.
    """
    two_theta = np.asarray(two_theta, dtype=float).ravel()
    beta = np.radians(np.asarray(breadth, dtype=float).ravel())
    sigma = np.radians(np.asarray(esd, dtype=float).ravel())
    if not two_theta.size == beta.size == sigma.size:
        raise ValueError(
            f"two_theta, breadth and esd differ in length: {two_theta.size}, "
            f"{beta.size}, {sigma.size}"
        )
    if two_theta.size < MIN_BREADTH_MODEL_POINTS:
        raise ValueError(
            f"Need at least {MIN_BREADTH_MODEL_POINTS} breadths, got {two_theta.size}"
        )
    if not np.all(np.isfinite(beta)):
        raise ValueError("Every breadth must be finite")
    if np.any(~np.isfinite(sigma) | (sigma <= 0.0)):
        raise ValueError("Every esd must be positive")
    for name, value in (("wavelength", wavelength), ("radius_mm", radius_mm), ("k", k)):
        if not value > 0.0:
            raise ValueError(f"{name} must be positive, got {value}")

    cos_theta = _cos_theta(two_theta)
    tan_theta = np.sqrt(1.0 - cos_theta**2) / cos_theta
    weights = 1.0 / sigma**2
    k_lambda = k * wavelength

    fits: dict[str, BreadthModelFit] = {}
    single_terms: dict[str, tuple[float, float]] = {}
    for name, shape in (
        ("size", 1.0 / cos_theta),
        ("strain", 4.0 * tan_theta),
        ("height", cos_theta),
    ):
        normal = float(np.sum(weights * shape**2))
        term = float(np.sum(weights * shape * beta)) / normal
        residuals = (beta - term * shape) / sigma
        reduced = float(np.sum(residuals**2)) / (beta.size - 1)
        esd_term = float(np.sqrt(reduced / normal))
        single_terms[name] = (term, esd_term)
        fits[name] = _breadth_model_fit(
            name, {name: (term, esd_term)}, residuals, 1, False, k_lambda, radius_mm
        )

    # Size and height in quadrature, in P and Q, started from a linear fit of
    # the squared breadths and polished on the breadths themselves.
    design = np.column_stack([1.0 / cos_theta**2, cos_theta**2])
    root_weights = 1.0 / (2.0 * np.maximum(np.abs(beta), sigma) * sigma)
    start, *_ = np.linalg.lstsq(
        design * root_weights[:, None], beta**2 * root_weights, rcond=None
    )
    scale = float(np.max(beta**2))
    start = np.clip(start, 0.0, None) + QUADRATURE_START_FLOOR * scale

    def residual(parameters: np.ndarray) -> np.ndarray:
        return (beta - np.sqrt(design @ parameters)) / sigma

    result = optimize.least_squares(
        residual, start, bounds=(0.0, np.inf), x_scale=scale, xtol=1e-15, ftol=1e-15
    )
    residuals = np.asarray(result.fun, dtype=float)
    # With one term at zero the fit is the other term alone, which the one
    # parameter fits have already found exactly. The optimiser only approaches
    # a bound from inside, so when either of those does as well, to within
    # QUADRATURE_BOUND_CHI_SQUARED, the best fit lies on that bound and is
    # taken from it.
    boundary = min((fits["size"], fits["height"]), key=lambda fit: fit.chi_squared)
    free_chi_squared = float(np.sum(residuals**2))
    if boundary.chi_squared <= free_chi_squared + QUADRATURE_BOUND_CHI_SQUARED:
        name = boundary.model
        term, esd_term = single_terms[name]
        # Its esd rescaled to the one fewer degree of freedom of this model.
        dof_ratio = (beta.size - 1) / (beta.size - 2)
        fits["size+height"] = _breadth_model_fit(
            "size+height",
            {name: (term, esd_term * np.sqrt(dof_ratio))},
            boundary.normalised_residuals,
            2,
            True,
            k_lambda,
            radius_mm,
        )
        return BreadthModels(fits)

    reduced = float(np.sum(residuals**2)) / (beta.size - 2)
    jacobian = np.asarray(result.jac, dtype=float)
    covariance = np.linalg.pinv(jacobian.T @ jacobian) * reduced
    # The esd of a square root, from the esd of its square.
    roots = np.sqrt(result.x)
    esd_roots = np.sqrt(np.abs(np.diag(covariance))) / (2.0 * roots)
    fits["size+height"] = _breadth_model_fit(
        "size+height",
        {
            "size": (float(roots[0]), float(esd_roots[0])),
            "height": (float(roots[1]), float(esd_roots[1])),
        },
        residuals,
        2,
        False,
        k_lambda,
        radius_mm,
    )
    return BreadthModels(fits)
