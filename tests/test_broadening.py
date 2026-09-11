"""Tests for xrdkit.broadening."""

import numpy as np
import pytest
from scipy import special

from xrdkit import (
    Caglioti,
    ProfileFit,
    correct_broadening,
    doublet_gaps,
    fit_caglioti,
    fit_profile,
    integral_breadth,
    kalpha2_position,
    pseudo_voigt,
    pseudo_voigt_components,
    pseudo_voigt_from_components,
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


def voigt_counts(
    x: np.ndarray, fwhm_gaussian: float, fwhm_lorentzian: float, area: float = 1e5
) -> np.ndarray:
    """A noise free Voigt of the given component widths, centred on zero."""
    sigma = fwhm_gaussian / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return area * special.voigt_profile(x, sigma, fwhm_lorentzian / 2.0)


def measured_fwhm(x: np.ndarray, y: np.ndarray) -> float:
    """The FWHM of a finely sampled single peak, by interpolation."""
    half = y.max() / 2.0
    above = np.nonzero(y >= half)[0]
    left, right = above[0], above[-1]
    low = np.interp(half, y[left - 1 : left + 1], x[left - 1 : left + 1])
    high = np.interp(half, y[right : right + 2][::-1], x[right : right + 2][::-1])
    return float(high - low)


def fitted_pseudo_voigt(fwhm_gaussian: float, fwhm_lorentzian: float) -> ProfileFit:
    """A symmetric single line pseudo-Voigt fitted to an exact Voigt."""
    x = np.arange(-3.0, 3.0, 0.002)
    y = voigt_counts(x, fwhm_gaussian, fwhm_lorentzian) + 100.0
    return fit_profile(
        x, y, 0.0, 0.1, window=(-3.0, 3.0), fit_asymmetry=False, intensity_ratio=0.0
    )


def kalpha2_of(two_theta: float) -> float:
    sin_theta = KALPHA2_RATIO * np.sin(np.radians(two_theta / 2.0))
    return 2.0 * np.degrees(np.arcsin(sin_theta))


class TestPseudoVoigtComponents:
    def test_pure_gaussian(self):
        assert pseudo_voigt_components(0.1, 0.0) == pytest.approx((0.1, 0.0))
        assert pseudo_voigt_from_components(0.1, 0.0) == pytest.approx((0.1, 0.0))

    def test_pure_lorentzian(self):
        assert pseudo_voigt_components(0.1, 1.0) == pytest.approx((0.0, 0.1))
        assert pseudo_voigt_from_components(0.0, 0.1) == pytest.approx((0.1, 1.0))

    @pytest.mark.parametrize("eta", [0.05, 0.3, 0.6, 0.9, 0.99])
    def test_round_trip(self, eta):
        gaussian, lorentzian = pseudo_voigt_components(0.08, eta)
        assert pseudo_voigt_from_components(gaussian, lorentzian) == pytest.approx(
            (0.08, eta), rel=1e-10
        )

    def test_combined_width_matches_the_true_voigt(self):
        x = np.linspace(-5.0, 5.0, 200001)
        true = measured_fwhm(x, voigt_counts(x, 1.0, 1.0))
        fwhm, _ = pseudo_voigt_from_components(1.0, 1.0)
        assert fwhm == pytest.approx(true, rel=0.005)

    @pytest.mark.parametrize(
        ("gaussian", "lorentzian"), [(0.06, 0.04), (0.05, 0.05), (0.03, 0.07)]
    )
    def test_components_of_a_fitted_voigt_are_recovered(self, gaussian, lorentzian):
        fit = fitted_pseudo_voigt(gaussian, lorentzian)
        recovered = pseudo_voigt_components(fit.fwhm, fit.eta)
        assert recovered == pytest.approx((gaussian, lorentzian), abs=0.004)

    @pytest.mark.parametrize(
        ("fwhm", "eta"), [(0.0, 0.5), (-0.1, 0.5), (0.1, -0.1), (0.1, 1.1)]
    )
    def test_rejects_bad_input(self, fwhm, eta):
        with pytest.raises(ValueError):
            pseudo_voigt_components(fwhm, eta)

    def test_rejects_bad_components(self):
        with pytest.raises(ValueError, match="negative"):
            pseudo_voigt_from_components(-0.1, 0.1)
        with pytest.raises(ValueError, match="positive"):
            pseudo_voigt_from_components(0.0, 0.0)


class TestIntegralBreadth:
    def test_limits(self):
        assert integral_breadth(0.1, 0.0) == pytest.approx(
            0.05 * np.sqrt(np.pi / np.log(2.0))
        )
        assert integral_breadth(0.0, 0.1) == pytest.approx(0.05 * np.pi)

    def test_matches_area_over_height_of_the_voigt(self):
        x = np.linspace(-200.0, 200.0, 2000001)
        y = voigt_counts(x, 0.06, 0.04, area=1.0)
        assert integral_breadth(0.06, 0.04) == pytest.approx(
            np.trapezoid(y, x) / y.max(), rel=1e-3
        )


class TestCorrectBroadening:
    # Sample and instrumental Voigts, as Gaussian and Lorentzian FWHM.
    SAMPLE = (0.05, 0.04)
    INSTRUMENT = (0.055, 0.03)

    def observed(self) -> tuple[float, float]:
        """The pseudo-Voigt of the two Voigts convolved, exactly."""
        return pseudo_voigt_from_components(
            np.hypot(self.SAMPLE[0], self.INSTRUMENT[0]),
            self.SAMPLE[1] + self.INSTRUMENT[1],
        )

    def test_recovers_the_sample_components(self):
        fwhm_obs, eta_obs = self.observed()
        fwhm_inst, eta_inst = pseudo_voigt_from_components(*self.INSTRUMENT)
        result = correct_broadening(fwhm_obs, eta_obs, fwhm_inst, eta_inst)

        assert not result.unresolved
        assert result.fwhm_gaussian == pytest.approx(self.SAMPLE[0], rel=1e-8)
        assert result.fwhm_lorentzian == pytest.approx(self.SAMPLE[1], rel=1e-8)
        fwhm, eta = pseudo_voigt_from_components(*self.SAMPLE)
        assert result.fwhm == pytest.approx(fwhm, rel=1e-8)
        assert result.eta == pytest.approx(eta, rel=1e-8)
        assert result.integral_breadth == pytest.approx(
            integral_breadth(*self.SAMPLE), rel=1e-8
        )
        # The simple estimates bracket the right answer.
        assert result.fwhm_linear < result.fwhm < result.fwhm_quadrature

    def test_fitted_broadened_peak_comes_back_to_the_sample_width(self):
        # Voigts convolve into a Voigt, so the observed profile is exact; it and
        # the instrument are then measured by fitting, as a pattern would be.
        observed = fitted_pseudo_voigt(
            np.hypot(self.SAMPLE[0], self.INSTRUMENT[0]),
            self.SAMPLE[1] + self.INSTRUMENT[1],
        )
        instrument = fitted_pseudo_voigt(*self.INSTRUMENT)
        result = correct_broadening(
            observed.fwhm, observed.eta, instrument.fwhm, instrument.eta
        )

        x = np.linspace(-2.0, 2.0, 400001)
        true = measured_fwhm(x, voigt_counts(x, *self.SAMPLE))
        assert result.fwhm == pytest.approx(true, rel=0.03)

    def test_simple_estimates_and_their_esds(self):
        result = correct_broadening(0.15, 0.6, 0.08, 0.6, 0.004, 0.0, 0.003, 0.0)
        assert result.fwhm_linear == pytest.approx(0.07)
        assert result.esd_fwhm_linear == pytest.approx(0.005)
        assert result.fwhm_quadrature == pytest.approx(np.sqrt(0.15**2 - 0.08**2))
        expected = np.hypot(0.15 * 0.004, 0.08 * 0.003) / result.fwhm_quadrature
        assert result.esd_fwhm_quadrature == pytest.approx(expected)

    def test_lorentzian_esd_is_the_combined_width_esd(self):
        # With both profiles Lorentzian the correction is a plain difference.
        result = correct_broadening(0.15, 1.0, 0.08, 1.0, 0.004, 0.0, 0.003, 0.0)
        assert result.fwhm == pytest.approx(0.07)
        assert result.esd_fwhm == pytest.approx(0.005, rel=1e-4)
        assert result.esd_fwhm_lorentzian == pytest.approx(0.005, rel=1e-4)

    def test_esds_grow_with_the_input_esds(self):
        small = correct_broadening(0.15, 0.7, 0.075, 0.62, 0.002, 0.02, 0.001, 0.01)
        large = correct_broadening(0.15, 0.7, 0.075, 0.62, 0.004, 0.04, 0.002, 0.02)
        assert large.esd_fwhm == pytest.approx(2.0 * small.esd_fwhm, rel=1e-3)
        assert small.esd_integral_breadth > 0.0

    def test_unresolved_within_two_combined_esds(self):
        # 0.004 above the instrument against a combined esd of 0.0025.
        result = correct_broadening(0.084, 0.6, 0.08, 0.6, 0.0015, 0.03, 0.002, 0.02)
        assert result.unresolved
        assert result.excess == pytest.approx(0.004 / 0.0025)
        assert result.fwhm is None
        assert result.fwhm_gaussian is None
        assert result.fwhm_linear is None
        assert result.integral_breadth is None

    def test_narrower_than_the_instrument_is_unresolved(self):
        result = correct_broadening(0.07, 0.6, 0.08, 0.6)
        assert result.unresolved
        assert result.excess == -np.inf

    def test_resolved_beyond_two_combined_esds(self):
        result = correct_broadening(0.086, 0.6, 0.08, 0.6, 0.0015, 0.03, 0.002, 0.02)
        assert not result.unresolved
        assert result.fwhm is not None and result.fwhm > 0.0

    def test_significance_is_a_parameter(self):
        arguments = (0.086, 0.6, 0.08, 0.6, 0.0015, 0.03, 0.002, 0.02)
        assert correct_broadening(*arguments, significance=3.0).unresolved

    def test_negative_gaussian_part_is_clipped_and_flagged(self):
        # A Lorentzian observed peak on a mostly Gaussian instrument.
        result = correct_broadening(0.10, 0.95, 0.075, 0.3, 0.001, 0.03, 0.0008, 0.02)
        assert result.gaussian_clipped
        assert not result.lorentzian_clipped
        assert result.fwhm_gaussian == 0.0
        assert result.esd_fwhm_gaussian > 0.0
        assert result.eta == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "arguments",
        [
            (0.0, 0.5, 0.08, 0.5),
            (0.1, 0.5, -0.08, 0.5),
            (0.1, 1.5, 0.08, 0.5),
            (0.1, 0.5, 0.08, -0.5),
            (0.1, 0.5, 0.08, 0.5, -0.001),
        ],
    )
    def test_rejects_bad_input(self, arguments):
        with pytest.raises(ValueError):
            correct_broadening(*arguments)


class TestDoubletGaps:
    def test_gaps_either_side(self):
        below, above = doublet_gaps(30.0, [29.0, 31.0])
        # Below, the nearest line is the K alpha 2 of 29.0, not 29.0 itself.
        assert below == pytest.approx(30.0 - kalpha2_of(29.0))
        assert above == pytest.approx(31.0 - kalpha2_of(30.0))

    def test_the_satellite_of_a_neighbour_counts(self):
        # The K alpha 2 line of 29.8 falls nearer to 30 than 29.8 itself.
        below, _ = doublet_gaps(30.0, [29.8])
        assert below == pytest.approx(30.0 - kalpha2_of(29.8))

    def test_line_inside_the_doublet_closes_both_gaps(self):
        assert doublet_gaps(60.0, [60.1]) == (0.0, 0.0)

    def test_open_sides_are_infinite(self):
        assert doublet_gaps(30.0, []) == (float("inf"), float("inf"))


def test_kalpha2_position():
    assert kalpha2_position(30.0) == pytest.approx(kalpha2_of(30.0))
    # The split widens with angle.
    assert kalpha2_position(90.0) - 90.0 > kalpha2_position(30.0) - 30.0
