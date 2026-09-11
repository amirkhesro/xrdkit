"""Crystallite size and microstrain from sample breadths.

Every function here takes breadths that already have the instrument taken out,
such as those :func:`~xrdkit.broadening.correct_broadening` gives, in degrees
of 2theta, and works in radians internally. Sizes come out in the units of the
wavelength, angstroms in practice.

:func:`scherrer_size` gives an apparent size per reflection from its breadth
alone, putting all of it down to size. :func:`williamson_hall` separates size
from strain across reflections by their different angular dependence,

    beta cos(theta) = K lambda / D + 4 epsilon sin(theta),

fitted as a straight line with a free intercept. :func:`component_size_strain`
is the cross check of de Keijser et al., J. Appl. Cryst. 15 (1982) 308, which
puts the Lorentzian part of each profile down to size and the Gaussian part to
strain. :func:`resolution_limit` gives the size above which a reflection could
not have been told apart from the instrument, which is the lower bound to
quote when a reflection shows no measurable broadening.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from xrdkit.broadening import DEFAULT_SIGNIFICANCE

__all__ = [
    "ComponentSizeStrain",
    "ResolutionLimit",
    "WilliamsonHall",
    "component_size_strain",
    "resolution_limit",
    "scherrer_size",
    "scherrer_size_integral",
    "williamson_hall",
]

# Scherrer constant for the FWHM of roughly spherical crystallites, and for the
# integral breadth, which gives the volume weighted size with K = 1.
SCHERRER_K = 0.9
INTEGRAL_BREADTH_K = 1.0

# A straight line has two parameters; one more point is the least that leaves
# a residual to scale the esds by.
WILLIAMSON_HALL_PARAMETERS = 2
MIN_WILLIAMSON_HALL_POINTS = WILLIAMSON_HALL_PARAMETERS + 1

# Each component fit has one parameter.
MIN_COMPONENT_POINTS = 2

# How the smallest detectable sample breadth is worked out from the smallest
# detectable excess of observed over instrumental FWHM: as Lorentzians, whose
# widths add, or as Gaussians, whose squares add.
CONVOLUTIONS = ("quadrature", "linear")


def _theta(two_theta: float | np.ndarray) -> np.ndarray:
    """theta in radians of ``two_theta`` in degrees, which must lie in (0, 180)."""
    two_theta = np.asarray(two_theta, dtype=float)
    if np.any(~np.isfinite(two_theta) | (two_theta <= 0.0) | (two_theta >= 180.0)):
        raise ValueError("two_theta must lie strictly between 0 and 180 degrees")
    return np.radians(two_theta / 2.0)


def _scalar(value: np.ndarray) -> float | np.ndarray:
    return float(value) if np.ndim(value) == 0 else value


def _check_esd(esd: np.ndarray, name: str, positive: bool = False) -> None:
    bad = esd <= 0.0 if positive else esd < 0.0
    if np.any(~np.isfinite(esd) | bad):
        rule = "positive" if positive else "finite and not negative"
        raise ValueError(f"Every {name} must be {rule}")


def scherrer_size(
    fwhm_deg: float | np.ndarray,
    two_theta: float | np.ndarray,
    wavelength: float,
    k: float = SCHERRER_K,
    esd_fwhm: float | np.ndarray = 0.0,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    """The Scherrer size D = K lambda / (beta cos(theta)) and its esd.

    The whole breadth is put down to size, so D is an apparent size along the
    normal to the reflecting planes, and a lower bound on it if the line is
    also strain broadened. The esd comes from that of the breadth alone, as
    esd(D) = D esd(beta) / beta; the position is taken as exact.

    Parameters
    ----------
    fwhm_deg
        Sample breadth, instrument taken out, in degrees of 2theta.
    two_theta
        Position of the reflection, in degrees.
    wavelength
        Wavelength, which sets the units of the size.
    k
        Scherrer constant, 0.9 for the FWHM by default.
    esd_fwhm
        Esd of the breadth, in degrees. Zero by default.

    Returns
    -------
    tuple
        The size and its esd, in the units of ``wavelength``, as floats for
        scalar input and arrays otherwise.

    Raises
    ------
    ValueError
        If a breadth, the wavelength or ``k`` is not positive, a position lies
        outside (0, 180) degrees, or an esd is negative.
    """
    breadth, esd = np.broadcast_arrays(
        np.asarray(fwhm_deg, dtype=float), np.asarray(esd_fwhm, dtype=float)
    )
    if np.any(~np.isfinite(breadth) | (breadth <= 0.0)):
        raise ValueError("Every breadth must be positive")
    if not wavelength > 0.0 or not k > 0.0:
        raise ValueError(f"wavelength and k must be positive, got {wavelength}, {k}")
    _check_esd(esd, "esd")
    size = k * wavelength / (np.radians(breadth) * np.cos(_theta(two_theta)))
    return _scalar(size), _scalar(size * esd / breadth)


def scherrer_size_integral(
    breadth_deg: float | np.ndarray,
    two_theta: float | np.ndarray,
    wavelength: float,
    k: float = INTEGRAL_BREADTH_K,
    esd_breadth: float | np.ndarray = 0.0,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    """The Scherrer size from an integral breadth, with K = 1 by default.

    Area over height is the breadth for which K = 1 gives the volume weighted
    column length whatever the crystallite shape, so this is
    :func:`scherrer_size` with that constant.
    """
    return scherrer_size(breadth_deg, two_theta, wavelength, k, esd_breadth)


def _size_from(
    intercept: float, esd_intercept: float, k: float, wavelength: float
) -> tuple[float | None, float | None]:
    """D = K lambda / intercept and its esd, or ``None`` for a non-positive one."""
    if not intercept > 0.0:
        return None, None
    size = k * wavelength / intercept
    return float(size), float(size * esd_intercept / intercept)


@dataclass
class WilliamsonHall:
    """A weighted straight line through beta cos(theta) against 4 sin(theta).

    ``x`` is 4 sin(theta), ``y`` the breadth times cos(theta) and ``esd_y`` its
    esd, all per reflection and ``y`` in radians. The intercept is K lambda / D
    and the slope the strain epsilon, the apparent strain of beta = 4 epsilon
    tan(theta), which is an upper limit on the local lattice strain. Esds are
    scaled by the reduced chi squared, as every other fit in xrdkit is.

    ``size`` is ``None`` when the intercept is not positive, since no size then
    fits. ``residuals`` are y less the line, in radians, and
    ``normalised_residuals`` the same over ``esd_y``.
    """

    size: float | None
    esd_size: float | None
    strain: float
    esd_strain: float
    intercept: float
    esd_intercept: float
    slope: float
    esd_slope: float
    covariance: np.ndarray
    reduced_chi_squared: float
    n_reflections: int
    k: float
    wavelength: float
    x: np.ndarray
    y: np.ndarray
    esd_y: np.ndarray
    residuals: np.ndarray
    normalised_residuals: np.ndarray

    @property
    def slope_significance(self) -> float:
        """The slope over its esd."""
        return self.slope / self.esd_slope if self.esd_slope > 0.0 else np.inf

    def line(self, x: float | np.ndarray) -> float | np.ndarray:
        """The fitted beta cos(theta) at ``x`` = 4 sin(theta), in radians."""
        return _scalar(self.intercept + self.slope * np.asarray(x, dtype=float))

    def line_esd(self, x: float | np.ndarray) -> float | np.ndarray:
        """The esd of :meth:`line` at ``x``, from the covariance."""
        x = np.asarray(x, dtype=float)
        gradient = np.stack([np.ones_like(x), x])
        variance = np.einsum("i...,ij,j...->...", gradient, self.covariance, gradient)
        return _scalar(np.sqrt(np.abs(variance)))


def williamson_hall(
    two_theta: np.ndarray,
    breadth_deg: np.ndarray,
    esd: np.ndarray,
    wavelength: float,
    k: float = SCHERRER_K,
) -> WilliamsonHall:
    """Fit beta cos(theta) = K lambda / D + epsilon 4 sin(theta).

    A weighted linear least squares fit with both the intercept and the slope
    free; the line is never forced through the origin, so a size too large to
    broaden anything shows up as an intercept consistent with zero rather than
    being assumed. Each point is weighted by 1 / (esd(beta) cos(theta))^2.

    Parameters
    ----------
    two_theta
        Positions, in degrees.
    breadth_deg
        Sample breadths, instrument taken out, in degrees of 2theta. The
        corrected FWHM with K = 0.9, or the integral breadth with K = 1.
    esd
        Their esds, in degrees; every one must be positive.
    wavelength
        Wavelength, which sets the units of the size.
    k
        Scherrer constant matching the kind of breadth.

    Raises
    ------
    ValueError
        If the arrays differ in length, there are fewer than three points, a
        position lies outside (0, 180) degrees, a breadth is not finite, an
        esd is not positive, or the wavelength or ``k`` is not positive.
    """
    two_theta = np.asarray(two_theta, dtype=float).ravel()
    breadth = np.asarray(breadth_deg, dtype=float).ravel()
    esd = np.asarray(esd, dtype=float).ravel()
    if not two_theta.size == breadth.size == esd.size:
        raise ValueError(
            f"two_theta, breadth_deg and esd differ in length: {two_theta.size}, "
            f"{breadth.size}, {esd.size}"
        )
    if two_theta.size < MIN_WILLIAMSON_HALL_POINTS:
        raise ValueError(
            f"Need at least {MIN_WILLIAMSON_HALL_POINTS} breadths to fit a line "
            f"with a residual, got {two_theta.size}"
        )
    if not np.all(np.isfinite(breadth)):
        raise ValueError("Every breadth must be finite")
    _check_esd(esd, "esd", positive=True)
    if not wavelength > 0.0 or not k > 0.0:
        raise ValueError(f"wavelength and k must be positive, got {wavelength}, {k}")

    theta = _theta(two_theta)
    x = 4.0 * np.sin(theta)
    y = np.radians(breadth) * np.cos(theta)
    esd_y = np.radians(esd) * np.cos(theta)
    weights = 1.0 / esd_y**2

    design = np.column_stack([np.ones_like(x), x])
    normal = design.T @ (design * weights[:, None])
    parameters = np.linalg.solve(normal, design.T @ (weights * y))
    residuals = y - design @ parameters
    degrees_of_freedom = x.size - WILLIAMSON_HALL_PARAMETERS
    reduced_chi_squared = float(np.sum(weights * residuals**2) / degrees_of_freedom)
    covariance = np.linalg.inv(normal) * reduced_chi_squared
    esd_intercept, esd_slope = (float(v) for v in np.sqrt(np.diag(covariance)))
    intercept, slope = (float(p) for p in parameters)
    size, esd_size = _size_from(intercept, esd_intercept, k, wavelength)

    return WilliamsonHall(
        size=size,
        esd_size=esd_size,
        strain=slope,
        esd_strain=esd_slope,
        intercept=intercept,
        esd_intercept=esd_intercept,
        slope=slope,
        esd_slope=esd_slope,
        covariance=covariance,
        reduced_chi_squared=reduced_chi_squared,
        n_reflections=int(x.size),
        k=float(k),
        wavelength=float(wavelength),
        x=x,
        y=y,
        esd_y=esd_y,
        residuals=residuals,
        normalised_residuals=residuals / esd_y,
    )


@dataclass
class ComponentSizeStrain:
    """Size from the Lorentzian breadths and strain from the Gaussian ones.

    ``size_term`` is the weighted mean of beta_L cos(theta), which is K lambda
    / D, and ``strain`` the weighted slope of beta_G against 4 tan(theta). Each
    is a one parameter fit of a component against its own angular dependence,
    and each reduced chi squared says how well that dependence holds; esds are
    scaled by them. ``size`` is ``None`` when ``size_term`` is not positive.
    Residuals are in radians, of beta_L cos(theta) and of beta_G.
    """

    size: float | None
    esd_size: float | None
    size_term: float
    esd_size_term: float
    strain: float
    esd_strain: float
    size_reduced_chi_squared: float
    strain_reduced_chi_squared: float
    n_reflections: int
    k: float
    wavelength: float
    size_residuals: np.ndarray
    strain_residuals: np.ndarray


def _proportional(
    x: np.ndarray, y: np.ndarray, esd_y: np.ndarray
) -> tuple[float, float, float, np.ndarray]:
    """Weighted y = p x: p, its esd scaled by the reduced chi squared, that, and
    the residuals."""
    weights = 1.0 / esd_y**2
    normal = float(np.sum(weights * x**2))
    value = float(np.sum(weights * x * y)) / normal
    residuals = y - value * x
    reduced_chi_squared = float(np.sum(weights * residuals**2) / (x.size - 1))
    return (
        value,
        float(np.sqrt(reduced_chi_squared / normal)),
        reduced_chi_squared,
        residuals,
    )


def component_size_strain(
    two_theta: np.ndarray,
    lorentzian_deg: np.ndarray,
    esd_lorentzian: np.ndarray,
    gaussian_deg: np.ndarray,
    esd_gaussian: np.ndarray,
    wavelength: float,
    k: float = INTEGRAL_BREADTH_K,
) -> ComponentSizeStrain:
    """Size from the Lorentzian part of each breadth, strain from the Gaussian.

    The cross check of de Keijser, Langford, Mittemeijer and Vogels, J. Appl.
    Cryst. 15 (1982) 308, on :func:`williamson_hall`: size broadening is taken
    to be Lorentzian, beta_L = K lambda / (D cos(theta)), and strain broadening
    Gaussian, beta_G = 4 epsilon tan(theta). Each component is fitted against
    its own dependence alone, which by that model has no constant term, so
    the size comes from the weighted mean of beta_L cos(theta) and the strain
    from the weighted slope of beta_G through 4 tan(theta).

    A component clipped to zero by the correction is a measurement of zero,
    with its esd, and is used as such.

    Parameters
    ----------
    two_theta
        Positions, in degrees.
    lorentzian_deg, esd_lorentzian
        Integral breadths of the sample Lorentzian and their esds, in degrees;
        (pi / 2) times the Lorentzian FWHM.
    gaussian_deg, esd_gaussian
        Integral breadths of the sample Gaussian and their esds, in degrees;
        sqrt(pi / ln 2) / 2 times the Gaussian FWHM.
    wavelength
        Wavelength, which sets the units of the size.
    k
        Scherrer constant, 1 for integral breadths by default.

    Raises
    ------
    ValueError
        If the arrays differ in length, there are fewer than two points, a
        position lies outside (0, 180) degrees, a breadth is negative or not
        finite, an esd is not positive, or the wavelength or ``k`` is not
        positive.
    """
    arrays = [
        np.asarray(a, dtype=float).ravel()
        for a in (two_theta, lorentzian_deg, esd_lorentzian, gaussian_deg, esd_gaussian)
    ]
    two_theta, lorentzian, esd_l, gaussian, esd_g = arrays
    if len({a.size for a in arrays}) != 1:
        raise ValueError("two_theta and the breadths and esds differ in length")
    if two_theta.size < MIN_COMPONENT_POINTS:
        raise ValueError(
            f"Need at least {MIN_COMPONENT_POINTS} reflections, got {two_theta.size}"
        )
    for name, breadth in (("Lorentzian", lorentzian), ("Gaussian", gaussian)):
        if np.any(~np.isfinite(breadth) | (breadth < 0.0)):
            raise ValueError(f"Every {name} breadth must be finite and not negative")
    _check_esd(esd_l, "Lorentzian esd", positive=True)
    _check_esd(esd_g, "Gaussian esd", positive=True)
    if not wavelength > 0.0 or not k > 0.0:
        raise ValueError(f"wavelength and k must be positive, got {wavelength}, {k}")

    theta = _theta(two_theta)
    cos_theta = np.cos(theta)
    size_term, esd_size_term, size_chi2, size_residuals = _proportional(
        np.ones_like(theta),
        np.radians(lorentzian) * cos_theta,
        np.radians(esd_l) * cos_theta,
    )
    strain, esd_strain, strain_chi2, strain_residuals = _proportional(
        4.0 * np.tan(theta), np.radians(gaussian), np.radians(esd_g)
    )
    size, esd_size = _size_from(size_term, esd_size_term, k, wavelength)
    return ComponentSizeStrain(
        size=size,
        esd_size=esd_size,
        size_term=size_term,
        esd_size_term=esd_size_term,
        strain=strain,
        esd_strain=esd_strain,
        size_reduced_chi_squared=size_chi2,
        strain_reduced_chi_squared=strain_chi2,
        n_reflections=int(two_theta.size),
        k=float(k),
        wavelength=float(wavelength),
        size_residuals=size_residuals,
        strain_residuals=strain_residuals,
    )


@dataclass
class ResolutionLimit:
    """The smallest sample breadth a reflection would show, and its size.

    ``threshold`` is the smallest excess of observed over instrumental FWHM
    that counts as resolved, ``significance`` combined esds, in degrees.
    ``breadth`` is the sample FWHM that would produce that excess under
    ``convolution``, in degrees, and ``size`` the Scherrer size of that
    breadth: the largest size the reflection could tell from the instrument,
    and so the lower bound on the size when the reflection is unresolved.
    """

    threshold: float | np.ndarray
    breadth: float | np.ndarray
    size: float | np.ndarray
    convolution: str
    significance: float
    k: float


def resolution_limit(
    two_theta: float | np.ndarray,
    fwhm_inst: float | np.ndarray,
    esd_inst: float | np.ndarray,
    esd_obs: float | np.ndarray,
    wavelength: float,
    k: float = SCHERRER_K,
    significance: float = DEFAULT_SIGNIFICANCE,
    convolution: str = "quadrature",
) -> ResolutionLimit:
    """The largest crystallite size a reflection can tell from the instrument.

    :func:`~xrdkit.broadening.correct_broadening` counts a reflection as
    resolved when its FWHM exceeds the instrumental FWHM by ``significance``
    times the combined esd, sqrt(esd_obs^2 + esd_inst^2). The sample breadth
    that just reaches that threshold depends on how the two profiles combine:

    ``"quadrature"``
        As Gaussians, sqrt((FWHM_inst + threshold)^2 - FWHM_inst^2). This is
        the larger breadth, so the smaller size, and the conservative bound:
        a sample breadth of any shape below it could have gone unseen.
    ``"linear"``
        As Lorentzians, the threshold itself. The smaller breadth and the
        larger size, the bound if the sample broadening is known to be
        Lorentzian.

    The breadth is then turned into a size by :func:`scherrer_size`.

    Parameters
    ----------
    two_theta
        Position of the reflection, in degrees.
    fwhm_inst, esd_inst
        Instrumental FWHM at that position and its esd, in degrees.
    esd_obs
        Esd of the observed FWHM, in degrees.
    wavelength
        Wavelength, which sets the units of the size.
    k
        Scherrer constant, 0.9 for the FWHM by default.
    significance
        Combined esds an excess must reach to count, as for
        :func:`~xrdkit.broadening.correct_broadening`.
    convolution
        ``"quadrature"`` or ``"linear"``, as above.

    Raises
    ------
    ValueError
        If ``convolution`` is unknown, an instrumental width or
        ``significance`` is not positive, both esds of a reflection are zero,
        or an esd is negative.
    """
    if convolution not in CONVOLUTIONS:
        raise ValueError(
            f"convolution must be one of {CONVOLUTIONS}, got {convolution!r}"
        )
    if not significance > 0.0:
        raise ValueError(f"significance must be positive, got {significance}")
    fwhm_inst = np.asarray(fwhm_inst, dtype=float)
    esd_inst = np.asarray(esd_inst, dtype=float)
    esd_obs = np.asarray(esd_obs, dtype=float)
    if np.any(~np.isfinite(fwhm_inst) | (fwhm_inst <= 0.0)):
        raise ValueError("Every instrumental width must be positive")
    _check_esd(esd_inst, "instrumental esd")
    _check_esd(esd_obs, "observed esd")

    threshold = significance * np.hypot(esd_obs, esd_inst)
    if np.any(threshold <= 0.0):
        raise ValueError("A reflection with no esd at all has no resolution limit")
    if convolution == "linear":
        breadth = threshold
    else:
        breadth = np.sqrt((fwhm_inst + threshold) ** 2 - fwhm_inst**2)
    size, _ = scherrer_size(breadth, two_theta, wavelength, k)
    return ResolutionLimit(
        threshold=_scalar(threshold),
        breadth=_scalar(breadth),
        size=size,
        convolution=convolution,
        significance=float(significance),
        k=float(k),
    )
