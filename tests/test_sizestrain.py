"""Tests for xrdkit.sizestrain."""

import numpy as np
import pytest

from xrdkit import (
    component_size_strain,
    correct_broadening,
    integral_breadth,
    resolution_limit,
    scherrer_size,
    scherrer_size_integral,
    williamson_hall,
)
from xrdkit.peaks import KALPHA1_WAVELENGTH

# Positions spread over the range the samples give isolated reflections in.
TWO_THETA = np.array([25.7, 27.7, 29.5, 34.5, 42.2, 49.5, 54.5, 55.7, 62.7, 68.5])

# The size and strain the synthetic breadths are built from, in angstroms and
# as a fraction, with the Scherrer constant of a FWHM.
SIZE = 800.0
STRAIN = 1.5e-3
K = 0.9

# Esd of each synthetic breadth, in degrees.
BREADTH_ESD = 0.004

# Relative tolerance on anything recovered from exact synthetic data.
EXACT = 1e-9

# Noisy copies of the synthetic breadths fitted to check the esds.
MONTE_CARLO_RUNS = 1000


def wh_breadths(
    two_theta: np.ndarray, size: float, strain: float, k: float = K
) -> np.ndarray:
    """Breadths in degrees obeying beta cos = K lambda / D + 4 epsilon sin."""
    theta = np.radians(two_theta / 2.0)
    size_term = 0.0 if np.isinf(size) else k * KALPHA1_WAVELENGTH / size
    beta = (size_term + 4.0 * strain * np.sin(theta)) / np.cos(theta)
    return np.degrees(beta)


class TestScherrer:
    def test_hand_calculated_value(self):
        # 0.2 degrees is 0.00349066 rad and cos(30 degrees) is 0.866025, so
        # D = 0.9 x 1.54056 / 0.00302300 = 458.652 angstroms.
        size, esd = scherrer_size(0.2, 60.0, KALPHA1_WAVELENGTH, esd_fwhm=0.02)
        assert size == pytest.approx(458.652, abs=1e-3)
        # A 10% esd on the breadth is a 10% esd on the size.
        assert esd == pytest.approx(45.8652, abs=1e-4)

    def test_integral_breadth_uses_k_of_one(self):
        # 1.54056 / 0.00302300 = 509.613 angstroms.
        size, _ = scherrer_size_integral(0.2, 60.0, KALPHA1_WAVELENGTH)
        assert size == pytest.approx(509.613, abs=1e-3)

    def test_arrays(self):
        breadth = wh_breadths(TWO_THETA, SIZE, 0.0)
        size, esd = scherrer_size(breadth, TWO_THETA, KALPHA1_WAVELENGTH)
        np.testing.assert_allclose(size, SIZE, rtol=EXACT)
        np.testing.assert_array_equal(esd, 0.0)

    @pytest.mark.parametrize(
        ("breadth", "two_theta"), [(0.0, 30.0), (-0.1, 30.0), (0.1, 0.0), (0.1, 180.0)]
    )
    def test_rejects_bad_input(self, breadth, two_theta):
        with pytest.raises(ValueError):
            scherrer_size(breadth, two_theta, KALPHA1_WAVELENGTH)

    def test_rejects_negative_esd(self):
        with pytest.raises(ValueError, match="esd"):
            scherrer_size(0.1, 30.0, KALPHA1_WAVELENGTH, esd_fwhm=-0.01)


class TestWilliamsonHall:
    def test_recovers_size_and_strain_exactly(self):
        breadth = wh_breadths(TWO_THETA, SIZE, STRAIN)
        esd = np.full_like(breadth, BREADTH_ESD)
        fit = williamson_hall(TWO_THETA, breadth, esd, KALPHA1_WAVELENGTH, K)
        assert fit.size == pytest.approx(SIZE, rel=EXACT)
        assert fit.strain == pytest.approx(STRAIN, rel=EXACT)
        assert fit.intercept == pytest.approx(K * KALPHA1_WAVELENGTH / SIZE, rel=EXACT)
        assert fit.slope == fit.strain
        assert fit.reduced_chi_squared == pytest.approx(0.0, abs=1e-12)
        np.testing.assert_allclose(fit.residuals, 0.0, atol=1e-15)
        assert fit.n_reflections == TWO_THETA.size

    def test_noisy_breadths_are_unbiased_with_calibrated_esds(self):
        # Over many noisy copies the intercept and slope average to the truth,
        # scatter as their esds say, and the reduced chi squared averages one.
        rng = np.random.default_rng(8)
        exact = wh_breadths(TWO_THETA, SIZE, STRAIN)
        esd = np.full_like(exact, BREADTH_ESD)
        fits = [
            williamson_hall(
                TWO_THETA,
                exact + rng.normal(0.0, BREADTH_ESD, exact.size),
                esd,
                KALPHA1_WAVELENGTH,
                K,
            )
            for _ in range(MONTE_CARLO_RUNS)
        ]
        intercept = np.array([fit.intercept for fit in fits])
        slope = np.array([fit.slope for fit in fits])
        chi2 = np.array([fit.reduced_chi_squared for fit in fits])
        # Unscaled esds, since each run's own is scaled by its own chi squared.
        esd_intercept = np.mean(
            [f.esd_intercept / np.sqrt(f.reduced_chi_squared) for f in fits]
        )
        esd_slope = np.mean(
            [f.esd_slope / np.sqrt(f.reduced_chi_squared) for f in fits]
        )
        standard_error = 1.0 / np.sqrt(MONTE_CARLO_RUNS)
        truth = K * KALPHA1_WAVELENGTH / SIZE
        assert abs(intercept.mean() - truth) < 4.0 * esd_intercept * standard_error
        assert abs(slope.mean() - STRAIN) < 4.0 * esd_slope * standard_error
        assert intercept.std() == pytest.approx(esd_intercept, rel=0.1)
        assert slope.std() == pytest.approx(esd_slope, rel=0.1)
        assert chi2.mean() == pytest.approx(1.0, rel=0.1)

        fit = fits[0]
        # Normalised residuals are the residuals over the esd of beta cos.
        np.testing.assert_allclose(
            fit.normalised_residuals, fit.residuals / fit.esd_y, rtol=EXACT
        )

    def test_intercept_is_not_forced_through_the_origin(self):
        # Strain alone: the fitted intercept must come out at zero on its own,
        # and a strain free line must keep its intercept.
        breadth = wh_breadths(TWO_THETA, np.inf, STRAIN)
        rng = np.random.default_rng(3)
        noisy = breadth + rng.normal(0.0, BREADTH_ESD, breadth.size)
        esd = np.full_like(breadth, BREADTH_ESD)
        fit = williamson_hall(TWO_THETA, noisy, esd, KALPHA1_WAVELENGTH, K)
        assert abs(fit.intercept) < 3.0 * fit.esd_intercept
        if fit.intercept <= 0.0:
            assert fit.size is None and fit.esd_size is None

        size_only = williamson_hall(
            TWO_THETA, wh_breadths(TWO_THETA, SIZE, 0.0), esd, KALPHA1_WAVELENGTH, K
        )
        assert size_only.slope == pytest.approx(0.0, abs=1e-12)
        assert size_only.size == pytest.approx(SIZE, rel=EXACT)

    def test_line_esd_follows_the_covariance(self):
        rng = np.random.default_rng(5)
        breadth = wh_breadths(TWO_THETA, SIZE, STRAIN)
        breadth = breadth + rng.normal(0.0, BREADTH_ESD, breadth.size)
        esd = np.full_like(breadth, BREADTH_ESD)
        fit = williamson_hall(TWO_THETA, breadth, esd, KALPHA1_WAVELENGTH, K)
        assert fit.line_esd(0.0) == pytest.approx(fit.esd_intercept, rel=EXACT)
        assert fit.line(0.0) == pytest.approx(fit.intercept, rel=EXACT)
        assert fit.slope_significance == pytest.approx(fit.slope / fit.esd_slope)

    def test_weights_follow_the_esds(self):
        # A wild point with a huge esd barely moves the line.
        breadth = wh_breadths(TWO_THETA, SIZE, STRAIN)
        esd = np.full_like(breadth, BREADTH_ESD)
        breadth[0] *= 3.0
        esd[0] = 1e3
        fit = williamson_hall(TWO_THETA, breadth, esd, KALPHA1_WAVELENGTH, K)
        assert fit.size == pytest.approx(SIZE, rel=1e-6)

    def test_needs_three_points(self):
        with pytest.raises(ValueError, match="at least 3"):
            williamson_hall(TWO_THETA[:2], [0.1, 0.1], [0.01, 0.01], 1.54)

    def test_needs_positive_esds(self):
        breadth = wh_breadths(TWO_THETA, SIZE, STRAIN)
        with pytest.raises(ValueError, match="esd"):
            williamson_hall(TWO_THETA, breadth, np.zeros_like(breadth), 1.54)


class TestComponents:
    def test_recovers_size_from_lorentzian_and_strain_from_gaussian(self):
        theta = np.radians(TWO_THETA / 2.0)
        lorentzian = np.degrees(KALPHA1_WAVELENGTH / (SIZE * np.cos(theta)))
        gaussian = np.degrees(4.0 * STRAIN * np.tan(theta))
        esd = np.full_like(lorentzian, BREADTH_ESD)
        fit = component_size_strain(
            TWO_THETA, lorentzian, esd, gaussian, esd, KALPHA1_WAVELENGTH
        )
        assert fit.size == pytest.approx(SIZE, rel=EXACT)
        assert fit.strain == pytest.approx(STRAIN, rel=EXACT)
        assert fit.size_reduced_chi_squared == pytest.approx(0.0, abs=1e-12)
        assert fit.strain_reduced_chi_squared == pytest.approx(0.0, abs=1e-12)
        assert fit.n_reflections == TWO_THETA.size

    def test_agrees_with_williamson_hall_on_lorentzian_breadths(self):
        # With Lorentzian size and strain broadening the breadths add, so the
        # single line and the components must find the same size and strain.
        theta = np.radians(TWO_THETA / 2.0)
        size_breadth = np.degrees(KALPHA1_WAVELENGTH / (SIZE * np.cos(theta)))
        strain_breadth = np.degrees(4.0 * STRAIN * np.tan(theta))
        esd = np.full_like(size_breadth, BREADTH_ESD)
        wh = williamson_hall(
            TWO_THETA, size_breadth + strain_breadth, esd, KALPHA1_WAVELENGTH, 1.0
        )
        components = component_size_strain(
            TWO_THETA, size_breadth, esd, strain_breadth, esd, KALPHA1_WAVELENGTH
        )
        assert wh.size == pytest.approx(components.size, rel=EXACT)
        assert wh.strain == pytest.approx(components.strain, rel=EXACT)

    def test_rejects_negative_breadths(self):
        esd = np.full(3, BREADTH_ESD)
        with pytest.raises(ValueError, match="Lorentzian"):
            component_size_strain(
                TWO_THETA[:3], [-0.1, 0.1, 0.1], esd, [0.0] * 3, esd, 1.54
            )


class TestResolutionLimit:
    # An instrumental FWHM of 0.075 +- 0.001 degrees and an observed FWHM esd of
    # 0.004 degrees at 40 degrees: a threshold of 2 sqrt(0.004^2 + 0.001^2) =
    # 0.0082462 degrees.
    TWO_THETA = 40.0
    FWHM_INST = 0.075
    ESD_INST = 0.001
    ESD_OBS = 0.004

    def limit(self, convolution="quadrature"):
        return resolution_limit(
            self.TWO_THETA,
            self.FWHM_INST,
            self.ESD_INST,
            self.ESD_OBS,
            KALPHA1_WAVELENGTH,
            convolution=convolution,
        )

    def test_quadrature_by_hand(self):
        limit = self.limit()
        assert limit.threshold == pytest.approx(0.0082462, abs=1e-7)
        # sqrt(0.0832462^2 - 0.075^2) = 0.0361238 degrees, which at cos(20)
        # = 0.939693 is 0.9 x 1.54056 / (0.000630485 x 0.939693) = 2340.26 A.
        assert limit.breadth == pytest.approx(0.0361238, abs=1e-7)
        assert limit.size == pytest.approx(2340.26, abs=0.01)

    def test_linear_is_the_threshold_itself(self):
        limit = self.limit("linear")
        assert limit.breadth == pytest.approx(limit.threshold)
        assert limit.size == pytest.approx(10251.88, abs=0.01)
        # Lorentzian widths add, so the same threshold hides a smaller breadth
        # and the bound on the size is larger.
        assert limit.size > self.limit().size

    def test_matches_the_resolution_test_of_correct_broadening(self):
        limit = self.limit()
        eta = 0.6
        just_above = self.FWHM_INST + limit.threshold * (1.0 + 1e-6)
        just_below = self.FWHM_INST + limit.threshold * (1.0 - 1e-6)
        for observed, resolved in ((just_above, True), (just_below, False)):
            correction = correct_broadening(
                observed, eta, self.FWHM_INST, eta, self.ESD_OBS, 0.0, self.ESD_INST
            )
            assert correction.unresolved is not resolved

    def test_arrays(self):
        limit = resolution_limit(
            TWO_THETA,
            np.full(TWO_THETA.size, self.FWHM_INST),
            self.ESD_INST,
            np.full(TWO_THETA.size, self.ESD_OBS),
            KALPHA1_WAVELENGTH,
        )
        assert limit.size.shape == TWO_THETA.shape
        # The same breadth stands for a larger size at higher angle.
        assert np.all(np.diff(limit.size) > 0.0)

    def test_rejects_unknown_convolution(self):
        with pytest.raises(ValueError, match="convolution"):
            self.limit("voigt")

    def test_rejects_zero_esds(self):
        with pytest.raises(ValueError, match="no esd"):
            resolution_limit(40.0, 0.075, 0.0, 0.0, KALPHA1_WAVELENGTH)


def test_component_breadths_are_the_integral_breadths_of_pure_lines():
    # The script converts component FWHM to integral breadths with these.
    assert integral_breadth(0.0, 1.0) == pytest.approx(np.pi / 2.0)
    assert integral_breadth(1.0, 0.0) == pytest.approx(
        0.5 * np.sqrt(np.pi / np.log(2.0))
    )
