"""Data quality of a scan, judged against what each workflow needs.

The thresholds are those of the criteria table in Section 4.2 of
``docs/USER_GUIDE.md``. Only the criteria that can be read off the scan itself
are applied: angular range, step size and counting statistics. Sample
preparation, standards and radiation are left to the user.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from xrdkit.io import XRDScan

__all__ = [
    "CRITERIA",
    "WORKFLOWS",
    "ScanQuality",
    "Verdict",
    "assess_scan",
    "format_report",
]

# Thresholds per workflow, from the criteria table of USER_GUIDE.md Section
# 4.2. CRITERIA is public, and is the programmatic route to those numbers: a
# dict of workflow name to a dict of criterion name to threshold, for a caller
# that would rather read them than copy them out of the guide. The criteria
# are:
#   angular_range                (low, high) degrees two theta the scan must cover
#   step_size                    (low, high) degrees the step must lie within
#   maximum                      least strongest intensity, counts
#   median                       least median intensity (the background), counts
#   peak_over_median             least ratio of maximum to median
#   high_angle_median            least median in the high angle third, counts
#   high_angle_peak_over_median  least ratio of maximum to median in the high
#                                angle third
_PLOTTING: dict[str, float | tuple[float, float]] = {
    "angular_range": (10.0, 80.0),
    "step_size": (0.01, 0.03),
    "maximum": 2000.0,
}

CRITERIA: dict[str, dict[str, float | tuple[float, float]]] = {
    "plotting": _PLOTTING,
    "phase_identification": {**_PLOTTING, "median": 100.0, "peak_over_median": 20.0},
    "le_bail": {
        "angular_range": (10.0, 120.0),
        "step_size": (0.013, 0.026),
        "maximum": 10000.0,
        "high_angle_peak_over_median": 10.0,
    },
    "rietveld": {
        "angular_range": (5.0, 130.0),
        "step_size": (0.01, 0.02),
        "maximum": 20000.0,
        "high_angle_median": 200.0,
    },
}

# The workflow names, in the order assess_scan judges them and format_report
# prints them. Public beside CRITERIA, so that a caller can walk the
# thresholds in the report's own order without hard coding the four names.
WORKFLOWS: tuple[str, ...] = tuple(CRITERIA)

# Slack on the step size bounds, in degrees, so that a step computed from the
# scan limits as 0.020000000001 still counts as 0.02.
STEP_TOLERANCE = 1e-6


@dataclass
class Verdict:
    """Whether a scan is good enough for one workflow, and if not, why not."""

    suitable: bool
    reasons: list[str] = field(default_factory=list)


@dataclass
class ScanQuality:
    """The data quality numbers of a scan, and a verdict per workflow.

    The first seven numbers are the ones ``scan_quality.py`` in the user guide
    prints. The high angle numbers are taken over the last third of the
    scanned range.
    """

    start_angle: float
    end_angle: float
    step_size: float
    points: int
    time_per_step: float | None
    maximum: float
    median: float
    peak_over_median: float
    high_angle_start: float
    high_angle_maximum: float
    high_angle_median: float
    verdicts: dict[str, Verdict] = field(default_factory=dict)


def _ratio(peak: float, background: float) -> float:
    if background > 0:
        return peak / background
    return math.inf if peak > 0 else 0.0


def _judge(quality: ScanQuality, criteria: dict) -> Verdict:
    reasons = []

    if "angular_range" in criteria:
        low, high = criteria["angular_range"]
        # A scan that starts or ends within one step of a limit covers it.
        slack = quality.step_size
        if quality.start_angle > low + slack or quality.end_angle < high - slack:
            reasons.append(
                f"range {quality.start_angle:.2f} to {quality.end_angle:.2f} "
                f"degrees, short of {low:g} to {high:g}"
            )

    if "step_size" in criteria:
        low, high = criteria["step_size"]
        step = quality.step_size
        if not low - STEP_TOLERANCE <= step <= high + STEP_TOLERANCE:
            reasons.append(f"step {step:.4f} degrees, outside {low:g} to {high:g}")

    if "maximum" in criteria and quality.maximum < criteria["maximum"]:
        reasons.append(
            f"strongest peak {quality.maximum:.0f} counts, "
            f"below {criteria['maximum']:g}"
        )

    if "median" in criteria and quality.median < criteria["median"]:
        reasons.append(
            f"median {quality.median:.0f} counts, below {criteria['median']:g}"
        )

    if (
        "peak_over_median" in criteria
        and quality.peak_over_median < criteria["peak_over_median"]
    ):
        reasons.append(
            f"peak over median {quality.peak_over_median:.1f}, "
            f"below {criteria['peak_over_median']:g}, weak phases may not be visible"
        )

    if (
        "high_angle_median" in criteria
        and quality.high_angle_median < criteria["high_angle_median"]
    ):
        reasons.append(
            f"high angle median {quality.high_angle_median:.0f} counts, "
            f"below {criteria['high_angle_median']:g}"
        )

    if "high_angle_peak_over_median" in criteria:
        ratio = _ratio(quality.high_angle_maximum, quality.high_angle_median)
        if ratio < criteria["high_angle_peak_over_median"]:
            reasons.append(
                f"high angle peak over median {ratio:.1f}, "
                f"below {criteria['high_angle_peak_over_median']:g}"
            )

    return Verdict(suitable=not reasons, reasons=reasons)


def assess_scan(scan: XRDScan) -> ScanQuality:
    """Measure the data quality of ``scan`` and judge it for each workflow.

    The median intensity stands in for the background, since most of the
    points in a powder pattern are background. The high angle third is the last
    third of the scanned range: the points at or above
    ``start_angle + 2/3 * (end_angle - start_angle)``. Its median is the high
    angle background and its maximum the strongest high angle peak.

    A scan covers an angular range when it starts and ends within one step of
    its limits. A step size passes when it lies inside the range the table
    gives, limits included.

    Raises
    ------
    ValueError
        If the scan holds no intensities.
    """
    intensity = np.asarray(scan.intensity, dtype=float)
    two_theta = np.asarray(scan.two_theta, dtype=float)
    if intensity.size == 0:
        raise ValueError("the scan holds no intensities")

    high_angle_start = scan.start_angle + 2.0 * (scan.end_angle - scan.start_angle) / 3
    high = intensity[two_theta >= high_angle_start]
    if high.size == 0:
        high = intensity[-1:]

    maximum = float(intensity.max())
    median = float(np.median(intensity))
    quality = ScanQuality(
        start_angle=float(scan.start_angle),
        end_angle=float(scan.end_angle),
        step_size=float(scan.step_size),
        points=int(intensity.size),
        time_per_step=(
            None if scan.time_per_step is None else float(scan.time_per_step)
        ),
        maximum=maximum,
        median=median,
        peak_over_median=_ratio(maximum, median),
        high_angle_start=float(high_angle_start),
        high_angle_maximum=float(high.max()),
        high_angle_median=float(np.median(high)),
    )
    quality.verdicts = {
        workflow: _judge(quality, criteria) for workflow, criteria in CRITERIA.items()
    }
    return quality


def format_report(quality: ScanQuality) -> str:
    """Return a plain text report of ``quality``.

    The numbers come first, in the order and wording of ``scan_quality.py`` in
    the user guide, followed by the high angle numbers, then one line per
    workflow: ``plotting: suitable`` or ``le_bail: not suitable (reason;
    reason)``.
    """
    q = quality
    lines = [
        f"range   {q.start_angle:.2f} to {q.end_angle:.2f} degrees",
        f"step    {q.step_size:.4f} degrees",
        f"points  {q.points}",
        "time    "
        + (
            "unknown"
            if q.time_per_step is None
            else f"{q.time_per_step:.1f} s per step"
        ),
        f"maximum {q.maximum:.0f} counts",
        f"median  {q.median:.0f} counts",
        f"peak over median {q.peak_over_median:.0f}",
        (
            f"high angle maximum {q.high_angle_maximum:.0f} counts, "
            f"from {q.high_angle_start:.2f} degrees"
        ),
        f"high angle median  {q.high_angle_median:.0f} counts",
        "",
    ]
    for workflow, verdict in q.verdicts.items():
        if verdict.suitable:
            lines.append(f"{workflow}: suitable")
        else:
            lines.append(f"{workflow}: not suitable ({'; '.join(verdict.reasons)})")
    return "\n".join(lines)
