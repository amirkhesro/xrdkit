"""Tests for xrdkit.plotting."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from xrdkit import XRDScan, plot_pattern, plot_stacked, save_figure
from xrdkit.plotting import (
    LABEL_HEIGHT,
    OFFSET_FACTOR,
    X_LABEL,
    X_MAJOR_TICK,
    X_MINOR_TICK,
    Y_LABEL,
    apply_style,
)


def make_scan(
    sample_id: str = "test", peaks: tuple[float, ...] = (22.0, 32.0)
) -> XRDScan:
    """Build a synthetic scan: Gaussian peaks on a flat background."""
    two_theta = np.linspace(10.0, 80.0, 500)
    intensity = np.full_like(two_theta, 100.0)
    for position in peaks:
        intensity += 5000.0 * np.exp(-((two_theta - position) ** 2) / (2 * 0.2**2))
    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=1.540598,
        start_angle=float(two_theta[0]),
        end_angle=float(two_theta[-1]),
        step_size=float(two_theta[1] - two_theta[0]),
        time_per_step=1.0,
        sample_id=sample_id,
        source_path=f"synthetic://{sample_id}",
    )


def test_plot_pattern_returns_figure_and_axes() -> None:
    scan = make_scan()
    fig, ax = plot_pattern(scan)

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_xlabel() == X_LABEL
    assert ax.get_ylabel() == Y_LABEL
    assert ax.get_xlim() == pytest.approx((scan.two_theta[0], scan.two_theta[-1]))
    assert len(ax.get_yticks()) == 0
    assert len(ax.lines) == 1


def test_plot_pattern_draws_on_supplied_axes() -> None:
    scan = make_scan()
    fig, ax = plot_pattern(scan)
    same_fig, same_ax = plot_pattern(make_scan("second"), ax=ax)

    assert same_ax is ax
    assert same_fig is fig
    assert len(ax.lines) == 2


@pytest.mark.parametrize("scale", ["linear", "sqrt", "log"])
@pytest.mark.parametrize("normalise", [False, True])
def test_plot_pattern_scales(scale: str, normalise: bool) -> None:
    _, ax = plot_pattern(make_scan(), scale=scale, normalise=normalise)
    values = ax.lines[0].get_ydata()

    assert np.all(np.isfinite(values))
    if normalise:
        assert np.max(values) == pytest.approx(1.0)


def test_plot_pattern_rejects_unknown_scale() -> None:
    with pytest.raises(ValueError, match="Unknown scale"):
        plot_pattern(make_scan(), scale="cubic")


def test_plot_stacked_rejects_unknown_scale() -> None:
    with pytest.raises(ValueError, match="Unknown scale"):
        plot_stacked([make_scan()], scale="cubic")


def test_plot_stacked_one_line_and_label_per_scan() -> None:
    scans = [make_scan("10"), make_scan("12", peaks=(25.0, 45.0))]
    fig, ax = plot_stacked(scans)

    assert isinstance(fig, Figure)
    assert len(ax.lines) == len(scans)
    assert len(ax.texts) == len(scans)
    assert [text.get_text() for text in ax.texts] == ["10", "12"]
    assert ax.get_xlabel() == X_LABEL
    assert ax.get_ylabel() == Y_LABEL


def test_plot_stacked_offsets_traces() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax = plot_stacked(scans, offset=2.0)
    lower, upper = (line.get_ydata() for line in ax.lines)

    assert np.min(upper) - np.min(lower) == pytest.approx(2.0)


def test_plot_stacked_custom_labels_and_length_check() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax = plot_stacked(scans, labels=["x = 0.10", "x = 0.12"])
    assert [text.get_text() for text in ax.texts] == ["x = 0.10", "x = 0.12"]

    with pytest.raises(ValueError, match="lengths must match"):
        plot_stacked(scans, labels=["only one"])


def test_plot_stacked_rejects_empty_list() -> None:
    with pytest.raises(ValueError, match="at least one scan"):
        plot_stacked([])


def test_save_figure_writes_png_and_pdf(tmp_path) -> None:
    fig, _ = plot_pattern(make_scan())
    written = save_figure(fig, tmp_path / "figures" / "pattern")

    assert [path.suffix for path in written] == [".png", ".pdf"]
    for path in written:
        assert path.is_file()
        assert path.stat().st_size > 0


def test_save_figure_keeps_dots_in_stem(tmp_path) -> None:
    fig, _ = plot_pattern(make_scan())
    written = save_figure(fig, tmp_path / "x0.10_calcined", formats=("png",))

    assert written[0].name == "x0.10_calcined.png"
    assert written[0].is_file()


def test_x_axis_uses_fixed_tick_spacing() -> None:
    # A scan running to nearly 100 degrees must still be labelled at its far end.
    scan = make_scan()
    scan.two_theta = np.linspace(10.0, 99.98, 500)
    _, ax = plot_pattern(scan)

    major = [tick for tick in ax.get_xticks() if 10.0 <= tick <= 99.98]
    minor = [tick for tick in ax.get_xticks(minor=True) if 10.0 <= tick <= 99.98]

    assert np.allclose(np.diff(major), X_MAJOR_TICK)
    assert max(major) == pytest.approx(90.0)

    # Minor ticks skip the positions already taken by a major tick, so check the
    # combined ladder rather than the minor ticks alone.
    combined = sorted(major + minor)
    assert np.allclose(np.diff(combined), X_MINOR_TICK)


def test_plot_stacked_labels_anchored_to_slot() -> None:
    scans = [make_scan("a"), make_scan("b"), make_scan("c")]
    offset = 2.0
    _, ax = plot_stacked(scans, offset=offset)

    expected = [index * offset + LABEL_HEIGHT * offset for index in range(len(scans))]
    assert [text.get_position()[1] for text in ax.texts] == pytest.approx(expected)

    # Labels sit inside the data range, just in from the right-hand end.
    x_max = float(scans[0].two_theta[-1])
    for text in ax.texts:
        assert text.get_position()[0] < x_max
        assert text.get_horizontalalignment() == "right"
        assert text.get_verticalalignment() == "top"


def test_plot_stacked_label_clears_a_tall_high_angle_peak() -> None:
    # The old anchoring used the trace maximum, so a strong high angle peak
    # dragged the label down onto the trace. The slot anchor must not move.
    low = make_scan("low", peaks=(22.0,))
    high = make_scan("high", peaks=(78.0,))
    offset = 2.0
    _, ax = plot_stacked([low, high], offset=offset)

    assert [text.get_position()[1] for text in ax.texts] == pytest.approx(
        [LABEL_HEIGHT * offset, offset + LABEL_HEIGHT * offset]
    )


def test_plot_stacked_y_limits_leave_room_for_top_label() -> None:
    scans = [make_scan("a"), make_scan("b")]
    offset = 2.0
    _, ax = plot_stacked(scans, offset=offset)
    bottom, top = ax.get_ylim()

    assert bottom <= 0.0
    assert top == pytest.approx(len(scans) * offset)
    # The topmost label must fall inside the axes.
    assert max(text.get_position()[1] for text in ax.texts) < top


def test_plot_stacked_default_offset_clears_the_trace_above() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax = plot_stacked(scans)
    lower, upper = (line.get_ydata() for line in ax.lines)
    offset = float(np.min(upper) - np.min(lower))

    assert offset == pytest.approx(OFFSET_FACTOR)
    # Tallest peak of the lower trace stays below the baseline of the upper one.
    assert np.max(lower) < np.min(upper)


def test_apply_style_sets_rcparams() -> None:
    apply_style()

    assert matplotlib.rcParams["font.size"] == 9
    assert matplotlib.rcParams["axes.linewidth"] == 0.8
    assert matplotlib.rcParams["xtick.direction"] == "in"
    assert matplotlib.rcParams["ytick.direction"] == "in"
    assert matplotlib.rcParams["xtick.top"] is True
    assert matplotlib.rcParams["ytick.right"] is True
    assert matplotlib.rcParams["xtick.minor.visible"] is True
    assert matplotlib.rcParams["axes.grid"] is False
    assert matplotlib.rcParams["axes.spines.top"] is True
    assert matplotlib.rcParams["axes.spines.right"] is True
    assert matplotlib.rcParams["savefig.dpi"] == 300
    assert matplotlib.rcParams["savefig.bbox"] == "tight"
