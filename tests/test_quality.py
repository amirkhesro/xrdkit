"""Tests for xrdkit.quality."""

import numpy as np
import pytest

from xrdkit import ScanQuality, XRDScan, assess_scan, format_report
from xrdkit.quality import CRITERIA, WORKFLOWS


def _scan(
    start: float,
    end: float,
    step: float,
    background: float,
    peak: float,
    high_background: float,
    high_peak: float,
    time_per_step: float = 2.0,
) -> XRDScan:
    """A flat scan with one peak point low and one in the high angle third."""
    n_points = round((end - start) / step) + 1
    two_theta = np.linspace(start, end, n_points)
    intensity = np.full(n_points, background)
    high = two_theta >= start + 2.0 * (end - start) / 3
    intensity[high] = high_background
    intensity[n_points // 4] = peak
    intensity[n_points - 5] = high_peak
    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=1.5406,
        start_angle=start,
        end_angle=end,
        step_size=(end - start) / (n_points - 1),
        time_per_step=time_per_step,
        sample_id="synthetic",
        source_path="synthetic.xrdml",
    )


def test_criteria_hold_the_guide_numbers() -> None:
    assert WORKFLOWS == ("plotting", "phase_identification", "le_bail", "rietveld")
    plotting = CRITERIA["plotting"]
    assert plotting["angular_range"] == (10.0, 80.0)
    assert plotting["step_size"] == (0.01, 0.03)
    assert plotting["maximum"] == 2000.0
    assert "peak_over_median" not in plotting
    assert CRITERIA["phase_identification"] == {
        **plotting,
        "median": 100.0,
        "peak_over_median": 20.0,
    }
    assert CRITERIA["le_bail"] == {
        "angular_range": (10.0, 120.0),
        "step_size": (0.013, 0.026),
        "maximum": 10000.0,
        "high_angle_peak_over_median": 10.0,
    }
    assert CRITERIA["rietveld"] == {
        "angular_range": (5.0, 130.0),
        "step_size": (0.01, 0.02),
        "maximum": 20000.0,
        "high_angle_median": 200.0,
    }


def test_a_good_scan_passes_every_workflow() -> None:
    scan = _scan(4.0, 132.0, 0.016, 500.0, 30000.0, 400.0, 8000.0)
    quality = assess_scan(scan)

    assert isinstance(quality, ScanQuality)
    assert quality.start_angle == 4.0
    assert quality.end_angle == 132.0
    assert quality.step_size == pytest.approx(0.016)
    assert quality.points == 8001
    assert quality.time_per_step == 2.0
    assert quality.maximum == 30000.0
    assert quality.median == 500.0
    assert quality.peak_over_median == pytest.approx(60.0)
    assert quality.high_angle_start == pytest.approx(4.0 + 128.0 * 2 / 3)
    assert quality.high_angle_maximum == 8000.0
    assert quality.high_angle_median == 400.0

    assert list(quality.verdicts) == list(WORKFLOWS)
    for workflow, verdict in quality.verdicts.items():
        assert verdict.suitable, (workflow, verdict.reasons)
        assert verdict.reasons == []


def test_weak_short_scan_fails_le_bail_on_counts_and_range() -> None:
    scan = _scan(10.0, 100.0, 0.02, 200.0, 6100.0, 200.0, 3000.0)
    quality = assess_scan(scan)

    assert quality.maximum == 6100.0
    assert quality.median == 200.0
    assert quality.peak_over_median == pytest.approx(30.5)
    assert quality.high_angle_median == 200.0
    assert quality.high_angle_maximum == 3000.0

    assert quality.verdicts["plotting"].suitable
    assert quality.verdicts["phase_identification"].suitable

    le_bail = quality.verdicts["le_bail"]
    assert not le_bail.suitable
    assert le_bail.reasons == [
        "range 10.00 to 100.00 degrees, short of 10 to 120",
        "strongest peak 6100 counts, below 10000",
    ]

    rietveld = quality.verdicts["rietveld"]
    assert not rietveld.suitable
    assert rietveld.reasons == [
        "range 10.00 to 100.00 degrees, short of 5 to 130",
        "strongest peak 6100 counts, below 20000",
    ]


def test_low_high_angle_background_fails_rietveld_only() -> None:
    scan = _scan(4.0, 132.0, 0.016, 500.0, 30000.0, 120.0, 2400.0)
    quality = assess_scan(scan)

    assert quality.median == 500.0
    assert quality.high_angle_median == 120.0
    assert quality.high_angle_maximum == 2400.0

    for workflow in ("plotting", "phase_identification", "le_bail"):
        assert quality.verdicts[workflow].suitable, workflow
    rietveld = quality.verdicts["rietveld"]
    assert not rietveld.suitable
    assert rietveld.reasons == ["high angle median 120 counts, below 200"]


def test_low_background_and_contrast_fail_phase_identification() -> None:
    # Background 80 counts, peak 1200: short of 2000, of 20 times, and of 100.
    scan = _scan(10.0, 90.0, 0.02, 80.0, 1200.0, 80.0, 400.0)
    quality = assess_scan(scan)

    assert quality.verdicts["plotting"].reasons == [
        "strongest peak 1200 counts, below 2000",
    ]
    assert quality.verdicts["phase_identification"].reasons == [
        "strongest peak 1200 counts, below 2000",
        "median 80 counts, below 100",
        "peak over median 15.0, below 20, weak phases may not be visible",
    ]
    assert "high angle peak over median 5.0, below 10" in (
        quality.verdicts["le_bail"].reasons
    )


def test_low_contrast_is_fine_for_plotting_but_not_phase_identification() -> None:
    # Peak over median 14, as in the guide's 10s.xrdml, with every other number
    # good for both workflows.
    scan = _scan(10.0, 90.0, 0.02, 1000.0, 14000.0, 1000.0, 3000.0)
    quality = assess_scan(scan)

    assert quality.peak_over_median == pytest.approx(14.0)
    assert quality.verdicts["plotting"].suitable
    assert quality.verdicts["plotting"].reasons == []
    phase = quality.verdicts["phase_identification"]
    assert not phase.suitable
    assert len(phase.reasons) == 1
    assert phase.reasons[0].startswith("peak over median 14.0, below 20")
    assert phase.reasons[0].endswith("weak phases may not be visible")


def test_step_and_range_limits_are_inclusive() -> None:
    # A step of exactly 0.02 lies inside Rietveld's 0.01 to 0.02, and a scan
    # starting at 10.01 covers plotting's 10 degrees to within a step.
    passes = assess_scan(_scan(5.0, 130.0, 0.02, 500.0, 30000.0, 400.0, 8000.0))
    assert passes.verdicts["rietveld"].suitable
    close = assess_scan(_scan(10.01, 80.01, 0.02, 500.0, 30000.0, 400.0, 8000.0))
    assert close.verdicts["plotting"].suitable

    coarse = assess_scan(_scan(5.0, 130.0, 0.025, 500.0, 30000.0, 400.0, 8000.0))
    assert coarse.verdicts["rietveld"].reasons == [
        "step 0.0250 degrees, outside 0.01 to 0.02"
    ]
    late = assess_scan(_scan(10.5, 80.0, 0.02, 500.0, 30000.0, 400.0, 8000.0))
    assert late.verdicts["plotting"].reasons == [
        "range 10.50 to 80.00 degrees, short of 10 to 80"
    ]


def test_zero_background_does_not_divide_by_zero() -> None:
    quality = assess_scan(_scan(10.0, 90.0, 0.02, 0.0, 5000.0, 0.0, 100.0))
    assert quality.peak_over_median == float("inf")
    assert quality.verdicts["plotting"].suitable


def test_empty_scan_is_rejected() -> None:
    scan = _scan(10.0, 90.0, 0.02, 500.0, 30000.0, 400.0, 8000.0)
    scan.intensity = np.array([])
    scan.two_theta = np.array([])
    with pytest.raises(ValueError, match="no intensities"):
        assess_scan(scan)


def test_format_report() -> None:
    scan = _scan(10.0, 100.0, 0.02, 200.0, 6100.0, 200.0, 3000.0)
    lines = format_report(assess_scan(scan)).splitlines()

    assert lines[:7] == [
        "range   10.00 to 100.00 degrees",
        "step    0.0200 degrees",
        "points  4501",
        "time    2.0 s per step",
        "maximum 6100 counts",
        "median  200 counts",
        "peak over median 30",
    ]
    assert lines[7] == "high angle maximum 3000 counts, from 70.00 degrees"
    assert lines[8] == "high angle median  200 counts"
    assert lines[-4:] == [
        "plotting: suitable",
        "phase_identification: suitable",
        (
            "le_bail: not suitable (range 10.00 to 100.00 degrees, short of 10 to 120; "
            "strongest peak 6100 counts, below 10000)"
        ),
        (
            "rietveld: not suitable (range 10.00 to 100.00 degrees, short of 5 to 130; "
            "strongest peak 6100 counts, below 20000)"
        ),
    ]
