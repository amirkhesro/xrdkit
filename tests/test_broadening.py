"""Tests for xrdkit.broadening."""

import numpy as np
import pytest

from xrdkit import (
    Caglioti,
    fit_caglioti,
    fit_profile,
    pseudo_voigt,
    split_pseudo_voigt,
)
from xrdkit.broadening import KALPHA2_INTENSITY_RATIO, MIN_CAGLIOTI_POINTS
from xrdkit.peaks import KALPHA2_RATIO

# A Caglioti function of the size the Aeris gives, in degrees squared, with the
# negative V that puts the narrowest peaks in the middle of the range.
U, V, W = 0.006, -0.006, 0.0067

# LaB6 positions for Cu K alpha 1 across 20-100 degrees.
LAB6_TWO_THETA = np.array(
    [21.36, 30.38, 37.44, 43.51, 48.96, 53.99, 63.22, 67.55, 71.75, 75.85,
     79.88, 83.85, 87.79, 95.68, 99.66]
)  # fmt: skip

# Relative scatter put on the synthetic widths.
WIDTH_NOISE = 0.01

# How closely U, V and W must come back, in degrees squared.
CAGLIOTI_EXACT = 1e-10
CAGLIOTI_NOISY = 0.002

# The profile the synthetic reflection is drawn with.
STEP = 0.0217
CENTRE = 43.48
FWHM = 0.072
ETA = 0.6
ASYMMETRY = -0.25
AREA = 400.0
BACKGROUND = 550.0

# How closely the profile fit must recover it.
POSITION_TOLERANCE = 0.002
FWHM_TOLERANCE = 0.003
ETA_TOLERANCE = 0.08
ASYMMETRY_TOLERANCE = 0.05


def caglioti_widths(two_theta: np.ndarray) -> np.ndarray:
    tan_theta = np.tan(np.radians(two_theta / 2.0))
    return np.sqrt(U * tan_theta**2 + V * tan_theta + W)


def make_fit(u: float, v: float, w: float) -> Caglioti:
    return Caglioti(
        u=u,
        v=v,
        w=w,
        esd_u=0.0,
        esd_v=0.0,
        esd_w=0.0,
        n_peaks=0,
        rms=0.0,
        covariance=np.zeros((3, 3)),
    )


def synthetic_doublet(
    centre: float = CENTRE,
    asymmetry: float = ASYMMETRY,
    end: float | None = None,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Counts of a split pseudo-Voigt K alpha doublet on a flat background."""
    two_theta = np.arange(centre - 2.0, end or centre + 2.0, STEP)
    satellite = 2.0 * np.degrees(
        np.arcsin(KALPHA2_RATIO * np.sin(np.radians(centre / 2.0)))
    )
    profile = split_pseudo_voigt(two_theta, centre, FWHM, ETA, asymmetry)
    profile += KALPHA2_INTENSITY_RATIO * split_pseudo_voigt(
        two_theta, satellite, FWHM, ETA, asymmetry
    )
    expected = BACKGROUND + AREA * profile
    counts = np.random.default_rng(seed).poisson(expected).astype(float)
    return two_theta, counts


class TestLineShapes:
    def test_pseudo_voigt_has_unit_area(self):
        x = np.linspace(-200.0, 200.0, 400001)
        area = np.trapezoid(pseudo_voigt(x, 0.0, 0.1, 0.5), x)
        assert area == pytest.approx(1.0, abs=1e-3)

    @pytest.mark.parametrize("asymmetry", [-0.5, 0.0, 0.3])
    def test_split_pseudo_voigt_has_unit_area(self, asymmetry):
        x = np.linspace(-200.0, 200.0, 400001)
        area = np.trapezoid(split_pseudo_voigt(x, 0.0, 0.1, 0.5, asymmetry), x)
        assert area == pytest.approx(1.0, abs=1e-3)

    def test_split_pseudo_voigt_without_asymmetry_is_symmetric(self):
        x = np.linspace(-1.0, 1.0, 201)
        np.testing.assert_allclose(
            split_pseudo_voigt(x, 0.0, 0.1, 0.4, 0.0), pseudo_voigt(x, 0.0, 0.1, 0.4)
        )

    def test_split_pseudo_voigt_keeps_its_fwhm(self):
        fwhm, asymmetry = 0.1, -0.4
        peak = float(split_pseudo_voigt(0.0, 0.0, fwhm, 0.5, asymmetry))
        below = -fwhm * (1.0 - asymmetry) / 2.0
        above = fwhm * (1.0 + asymmetry) / 2.0
        for x in (below, above):
            assert float(split_pseudo_voigt(x, 0.0, fwhm, 0.5, asymmetry)) == (
                pytest.approx(peak / 2.0)
            )


class TestFitProfile:
    def test_recovers_the_doublet(self):
        two_theta, counts = synthetic_doublet()
        fit = fit_profile(two_theta, counts, CENTRE + 0.01, 0.09)

        assert fit.converged
        assert not fit.truncated
        assert fit.two_theta == pytest.approx(CENTRE, abs=POSITION_TOLERANCE)
        assert fit.fwhm == pytest.approx(FWHM, abs=FWHM_TOLERANCE)
        assert fit.eta == pytest.approx(ETA, abs=ETA_TOLERANCE)
        assert fit.asymmetry == pytest.approx(ASYMMETRY, abs=ASYMMETRY_TOLERANCE)
        assert fit.reduced_chi_squared == pytest.approx(1.0, abs=0.5)
        # The esds describe the scatter: the truth is within a few of them.
        assert abs(fit.fwhm - FWHM) < 4.0 * fit.esd_fwhm
        assert fit.esd_asymmetry is not None and fit.esd_asymmetry > 0.0

    def test_satellite_sits_at_the_longer_wavelength(self):
        two_theta, counts = synthetic_doublet()
        fit = fit_profile(two_theta, counts, CENTRE, 0.09)
        sin_theta = np.sin(np.radians(fit.two_theta / 2.0))
        expected = 2.0 * np.degrees(np.arcsin(KALPHA2_RATIO * sin_theta))
        assert fit.kalpha2_two_theta == pytest.approx(expected)
        assert fit.kalpha2_two_theta > fit.two_theta

    def test_symmetric_fit_holds_asymmetry_at_zero(self):
        two_theta, counts = synthetic_doublet(asymmetry=0.0)
        fit = fit_profile(two_theta, counts, CENTRE, 0.09, fit_asymmetry=False)
        assert fit.asymmetry == 0.0
        assert fit.esd_asymmetry is None
        assert fit.fwhm == pytest.approx(FWHM, abs=FWHM_TOLERANCE)

    def test_window_limits_the_points(self):
        two_theta, counts = synthetic_doublet()
        fit = fit_profile(two_theta, counts, CENTRE, 0.09, window=(43.0, 44.0))
        assert 43.0 <= fit.window[0] and fit.window[1] <= 44.0
        assert fit.n_points == fit.fitted.size

    def test_flags_a_doublet_cut_by_the_end_of_the_scan(self):
        # The data stop just past the satellite, as they do for the last LaB6
        # reflection of a scan to 100 degrees.
        two_theta, counts = synthetic_doublet(end=CENTRE + 0.18)
        fit = fit_profile(two_theta, counts, CENTRE, 0.09)
        assert fit.truncated

    def test_rejects_a_width_that_is_not_positive(self):
        two_theta, counts = synthetic_doublet()
        with pytest.raises(ValueError, match="fwhm_guess"):
            fit_profile(two_theta, counts, CENTRE, 0.0)

    def test_rejects_a_window_with_too_few_points(self):
        two_theta, counts = synthetic_doublet()
        with pytest.raises(ValueError, match="points"):
            fit_profile(two_theta, counts, CENTRE, 0.09, window=(43.47, 43.52))


class TestFitCaglioti:
    def test_recovers_exact_widths(self):
        fit = fit_caglioti(LAB6_TWO_THETA, caglioti_widths(LAB6_TWO_THETA))
        assert fit.u == pytest.approx(U, abs=CAGLIOTI_EXACT)
        assert fit.v == pytest.approx(V, abs=CAGLIOTI_EXACT)
        assert fit.w == pytest.approx(W, abs=CAGLIOTI_EXACT)
        assert fit.rms == pytest.approx(0.0, abs=1e-10)
        assert fit.n_peaks == LAB6_TWO_THETA.size

    def test_recovers_noisy_widths_within_their_esds(self):
        rng = np.random.default_rng(1)
        widths = caglioti_widths(LAB6_TWO_THETA)
        widths *= 1.0 + WIDTH_NOISE * rng.standard_normal(widths.size)
        esd = WIDTH_NOISE * widths
        fit = fit_caglioti(LAB6_TWO_THETA, widths, 1.0 / (2.0 * widths * esd) ** 2)

        for value, true, esd_value in (
            (fit.u, U, fit.esd_u),
            (fit.v, V, fit.esd_v),
            (fit.w, W, fit.esd_w),
        ):
            assert value == pytest.approx(true, abs=CAGLIOTI_NOISY)
            assert abs(value - true) < 4.0 * esd_value
        assert fit.rms == pytest.approx(WIDTH_NOISE * np.mean(widths), rel=0.6)

    def test_esds_are_the_root_of_the_covariance_diagonal(self):
        rng = np.random.default_rng(2)
        widths = caglioti_widths(LAB6_TWO_THETA) + 0.001 * rng.standard_normal(
            LAB6_TWO_THETA.size
        )
        fit = fit_caglioti(LAB6_TWO_THETA, widths)
        assert fit.covariance.shape == (3, 3)
        np.testing.assert_allclose(fit.covariance, fit.covariance.T)
        np.testing.assert_allclose(
            np.sqrt(np.diag(fit.covariance)), [fit.esd_u, fit.esd_v, fit.esd_w]
        )

    def test_weights_pull_the_fit_towards_the_trusted_points(self):
        widths = caglioti_widths(LAB6_TWO_THETA)
        widths[0] += 0.02
        weights = np.ones_like(widths)
        weights[0] = 1e-6
        plain = fit_caglioti(LAB6_TWO_THETA, widths)
        weighted = fit_caglioti(LAB6_TWO_THETA, widths, weights)
        assert abs(weighted.w - W) < abs(plain.w - W)
        assert weighted.w == pytest.approx(W, abs=1e-4)

    def test_fwhm_matches_the_widths_it_was_fitted_to(self):
        widths = caglioti_widths(LAB6_TWO_THETA)
        fit = fit_caglioti(LAB6_TWO_THETA, widths)
        np.testing.assert_allclose(fit.fwhm(LAB6_TWO_THETA), widths)
        assert isinstance(fit.fwhm(45.0), float)

    @pytest.mark.parametrize(
        ("two_theta", "fwhm", "weights", "match"),
        [
            (LAB6_TWO_THETA[:3], [0.07] * 3, None, "at least"),
            (LAB6_TWO_THETA[:5], [0.07] * 4, None, "differ in length"),
            ([-20.0, 30.0, 40.0, 50.0], [0.07] * 4, None, "between 0 and 180"),
            ([20.0, 30.0, 40.0, 180.0], [0.07] * 4, None, "between 0 and 180"),
            (LAB6_TWO_THETA[:4], [0.07, -0.07, 0.07, 0.07], None, "positive"),
            (LAB6_TWO_THETA[:4], [0.07, 0.0, 0.07, 0.07], None, "positive"),
            (LAB6_TWO_THETA[:4], [0.07] * 4, [1.0, -1.0, 1.0, 1.0], "weight"),
        ],
    )
    def test_rejects_bad_input(self, two_theta, fwhm, weights, match):
        assert MIN_CAGLIOTI_POINTS == 4
        with pytest.raises(ValueError, match=match):
            fit_caglioti(two_theta, fwhm, weights)


class TestCagliotiFwhm:
    def test_negative_argument_under_the_root_gives_zero(self):
        # U tan^2 + V tan + W is negative for tan(theta) between 0.4 and 0.6.
        fit = make_fit(u=1.0, v=-1.0, w=0.24)
        widths = fit.fwhm(np.array([20.0, 57.0, 120.0]))
        assert np.all(np.isfinite(widths))
        assert widths[1] == 0.0
        assert widths[0] > 0.0 and widths[2] > 0.0

    def test_negative_argument_esd_is_nan(self):
        fit = make_fit(u=1.0, v=-1.0, w=0.24)
        assert np.isnan(fit.fwhm_esd(57.0))

    def test_rejects_angles_outside_the_range(self):
        fit = make_fit(U, V, W)
        for two_theta in (-10.0, 0.0, 180.0):
            with pytest.raises(ValueError, match="between 0 and 180"):
                fit.fwhm(two_theta)

    def test_fwhm_esd_propagates_the_covariance(self):
        rng = np.random.default_rng(3)
        widths = caglioti_widths(LAB6_TWO_THETA) + 0.001 * rng.standard_normal(
            LAB6_TWO_THETA.size
        )
        fit = fit_caglioti(LAB6_TWO_THETA, widths)
        tan_theta = np.tan(np.radians(30.0 / 2.0))
        gradient = np.array([tan_theta**2, tan_theta, 1.0])
        expected = np.sqrt(gradient @ fit.covariance @ gradient) / (
            2.0 * fit.fwhm(30.0)
        )
        assert fit.fwhm_esd(30.0) == pytest.approx(expected)
        assert fit.fwhm_esd(np.array([30.0, 60.0])).shape == (2,)
