"""Tests for xrdkit.plotting."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.text import Text

from xrdkit import (
    IndexedPeak,
    Peak,
    Reflection,
    XRDScan,
    annotate_hkl,
    mark_peaks,
    plot_pattern,
    plot_stacked,
    save_figure,
)
from xrdkit.plotting import (
    AMBIGUOUS_SEPARATOR,
    CHARACTER_WIDTH,
    HKL_AMBIGUOUS,
    HKL_FONTSIZE,
    HKL_LABEL_HEIGHT,
    HKL_MAX_LEVELS,
    HKL_MIN_RELATIVE_INTENSITY,
    HKL_MIN_SEPARATION,
    HKL_THIN_SPACE,
    LABEL_GAP_FACTOR,
    LABEL_HEIGHT,
    OFFSET_FACTOR,
    PEAK_MARKER,
    UPRIGHT_LABEL_WIDTH,
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
    fig, ax, _ = plot_stacked(scans)

    assert isinstance(fig, Figure)
    assert len(ax.lines) == len(scans)
    assert len(ax.texts) == len(scans)
    assert [text.get_text() for text in ax.texts] == ["10", "12"]
    assert ax.get_xlabel() == X_LABEL
    assert ax.get_ylabel() == Y_LABEL


def test_plot_stacked_offsets_traces() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax, _ = plot_stacked(scans, offset=2.0)
    lower, upper = (line.get_ydata() for line in ax.lines)

    assert np.min(upper) - np.min(lower) == pytest.approx(2.0)


def test_plot_stacked_custom_labels_and_length_check() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax, _ = plot_stacked(scans, labels=["x = 0.10", "x = 0.12"])
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
    _, ax, _ = plot_stacked(scans, offset=offset)

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
    _, ax, _ = plot_stacked([low, high], offset=offset)

    assert [text.get_position()[1] for text in ax.texts] == pytest.approx(
        [LABEL_HEIGHT * offset, offset + LABEL_HEIGHT * offset]
    )


def test_plot_stacked_y_limits_leave_room_for_top_label() -> None:
    scans = [make_scan("a"), make_scan("b")]
    offset = 2.0
    _, ax, _ = plot_stacked(scans, offset=offset)
    bottom, top = ax.get_ylim()

    assert bottom <= 0.0
    assert top == pytest.approx(len(scans) * offset)
    # The topmost label must fall inside the axes.
    assert max(text.get_position()[1] for text in ax.texts) < top


def test_plot_stacked_default_offset_clears_the_trace_above() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax, _ = plot_stacked(scans)
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


# Two peaks inside HKL_MIN_SEPARATION of each other, then one well clear.
CLOSE_PAIR = (40.000, 40.300)
FAR_PEAK = 60.000

# Base the labels are written at in the synthetic annotation tests.
LABEL_BASE = 1.0


def make_indexed_peak(
    two_theta: float,
    hkl: tuple[int, int, int],
    relative_intensity: float = 50.0,
    extra_candidates: tuple[tuple[int, int, int], ...] = (),
) -> IndexedPeak:
    """An IndexedPeak at a position, assigned to ``hkl``, with optional rivals."""
    peak = Peak(
        two_theta=two_theta,
        intensity=1000.0,
        prominence=900.0,
        fwhm=0.15,
        d_spacing=1.5406 / (2.0 * np.sin(np.radians(two_theta / 2.0))),
        relative_intensity=relative_intensity,
    )

    def reflection(indices: tuple[int, int, int]) -> Reflection:
        h, k, l = indices
        return Reflection(h=h, k=k, l=l, d_spacing=peak.d_spacing, two_theta=two_theta)

    assigned = reflection(hkl)
    candidates = [assigned] + [reflection(other) for other in extra_candidates]
    return IndexedPeak(
        peak=peak,
        corrected_two_theta=two_theta,
        reflection=assigned,
        difference=0.0,
        candidates=candidates,
    )


def make_indexed() -> list[IndexedPeak]:
    """Two peaks 0.3 degrees apart and one far away, all singly indexed.

    The close pair carry different intensities, so which of them annotate_hkl
    drops when they cannot both be placed is well defined.
    """
    return [
        make_indexed_peak(CLOSE_PAIR[0], (3, 1, 1), relative_intensity=60.0),
        make_indexed_peak(CLOSE_PAIR[1], (4, 2, 0), relative_intensity=40.0),
        make_indexed_peak(FAR_PEAK, (5, 5, 0), relative_intensity=50.0),
    ]


def annotated_axes() -> Axes:
    """Axes carrying one trace, ready to be annotated."""
    _, ax = plot_pattern(make_scan())
    return ax


def test_annotate_hkl_returns_one_text_per_indexed_peak() -> None:
    texts = annotate_hkl(annotated_axes(), make_indexed(), LABEL_BASE, max_levels=2)

    assert len(texts) == 3
    assert all(isinstance(text, Text) for text in texts)
    assert [text.get_text() for text in texts] == ["311", "420", "550"]


def test_annotate_hkl_places_labels_at_the_peak_positions() -> None:
    texts = annotate_hkl(annotated_axes(), make_indexed(), LABEL_BASE, max_levels=2)

    assert [text.get_position()[0] for text in texts] == pytest.approx(
        [CLOSE_PAIR[0], CLOSE_PAIR[1], FAR_PEAK]
    )


def test_annotate_hkl_drops_the_weaker_of_a_crowded_pair_in_one_row() -> None:
    texts = annotate_hkl(
        annotated_axes(),
        make_indexed(),
        LABEL_BASE,
        min_separation=HKL_MIN_SEPARATION,
    )

    # The pair is 0.3 degrees apart, inside the separation given, and one row
    # has nowhere to put the second of them.
    assert CLOSE_PAIR[1] - CLOSE_PAIR[0] < HKL_MIN_SEPARATION
    # 420 is the weaker of the pair, so 311 keeps the room.
    assert [text.get_text() for text in texts] == ["311", "550"]
    assert all(text.get_position()[1] == pytest.approx(LABEL_BASE) for text in texts)


def test_annotate_hkl_puts_a_crowded_pair_on_two_levels() -> None:
    ax = annotated_axes()
    texts = annotate_hkl(ax, make_indexed(), LABEL_BASE, max_levels=2)
    heights = [text.get_position()[1] for text in texts]

    assert [text.get_text() for text in texts] == ["311", "420", "550"]
    low, high = ax.get_ylim()
    step = HKL_LABEL_HEIGHT * (high - low)
    # The stronger 311 takes level 0 and the weaker 420 is raised one level.
    assert heights[0] == pytest.approx(LABEL_BASE)
    assert heights[1] == pytest.approx(LABEL_BASE + step)
    # The far peak has room on level 0 of its own.
    assert heights[2] == pytest.approx(LABEL_BASE)


def test_annotate_hkl_rejects_a_level_count_below_one() -> None:
    with pytest.raises(ValueError, match="at least one level"):
        annotate_hkl(annotated_axes(), make_indexed(), LABEL_BASE, max_levels=0)


def test_annotate_hkl_honours_an_explicit_label_height() -> None:
    texts = annotate_hkl(
        annotated_axes(), make_indexed(), LABEL_BASE, label_height=0.5, max_levels=2
    )
    heights = [text.get_position()[1] for text in texts]

    assert heights == pytest.approx([LABEL_BASE, LABEL_BASE + 0.5, LABEL_BASE])


def test_annotate_hkl_fills_levels_strongest_first() -> None:
    indexed = [
        make_indexed_peak(40.0, (3, 1, 1), relative_intensity=30.0),
        make_indexed_peak(40.2, (4, 2, 0), relative_intensity=90.0),
        make_indexed_peak(40.4, (5, 5, 0), relative_intensity=60.0),
    ]
    texts = annotate_hkl(
        annotated_axes(), indexed, LABEL_BASE, label_height=0.5, max_levels=3
    )
    heights = {text.get_text(): text.get_position()[1] for text in texts}

    # 420 is strongest and takes level 0, then 550, then the weakest 311.
    assert heights == {
        "420": pytest.approx(LABEL_BASE),
        "550": pytest.approx(LABEL_BASE + 0.5),
        "311": pytest.approx(LABEL_BASE + 1.0),
    }


def test_annotate_hkl_skips_what_will_not_fit_in_the_levels_given() -> None:
    indexed = [
        make_indexed_peak(40.0, (3, 1, 1), relative_intensity=30.0),
        make_indexed_peak(40.2, (4, 2, 0), relative_intensity=90.0),
        make_indexed_peak(40.4, (5, 5, 0), relative_intensity=60.0),
    ]
    texts = annotate_hkl(
        annotated_axes(), indexed, LABEL_BASE, label_height=0.5, max_levels=2
    )

    # Only the two strongest of the run find a level; 311 loses its label.
    assert [text.get_text() for text in texts] == ["420", "550"]


def test_annotate_hkl_min_relative_intensity_filters() -> None:
    indexed = make_indexed()
    indexed[1].peak.relative_intensity = 2.0

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE)

    assert [text.get_text() for text in texts] == ["311", "550"]
    # Below the default cut-off but above a lower one.
    assert (
        len(annotate_hkl(annotated_axes(), indexed, LABEL_BASE, 1.0, max_levels=2)) == 3
    )


def test_annotate_hkl_skips_unindexed_peaks() -> None:
    indexed = make_indexed()
    indexed[1].reflection = None
    indexed[1].candidates = []

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE)

    assert [text.get_text() for text in texts] == ["311", "550"]


def test_annotate_hkl_ambiguous_first_labels_the_assignment() -> None:
    indexed = make_indexed()
    indexed[0] = make_indexed_peak(
        CLOSE_PAIR[0], (3, 1, 1), extra_candidates=((4, 2, 0),)
    )

    texts = annotate_hkl(
        annotated_axes(), indexed, LABEL_BASE, ambiguous="first", max_levels=2
    )

    assert [text.get_text() for text in texts] == ["311", "420", "550"]


def test_annotate_hkl_ambiguous_all_joins_candidates_with_a_solidus() -> None:
    indexed = make_indexed()
    indexed[0] = make_indexed_peak(
        CLOSE_PAIR[0], (3, 1, 1), extra_candidates=((4, 2, 0),)
    )

    texts = annotate_hkl(
        annotated_axes(), indexed, LABEL_BASE, ambiguous="all", max_levels=2
    )

    assert texts[0].get_text() == f"311{AMBIGUOUS_SEPARATOR}420"
    assert AMBIGUOUS_SEPARATOR in texts[0].get_text()
    # A peak with a single candidate is untouched.
    assert texts[2].get_text() == "550"


def test_annotate_hkl_ambiguous_skip_drops_the_peak() -> None:
    indexed = make_indexed()
    indexed[0] = make_indexed_peak(
        CLOSE_PAIR[0], (3, 1, 1), extra_candidates=((4, 2, 0),)
    )

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE, ambiguous="skip")

    assert len(texts) == 2
    assert [text.get_text() for text in texts] == ["420", "550"]


def test_annotate_hkl_rejects_an_unknown_ambiguous_mode() -> None:
    with pytest.raises(ValueError, match="Unknown ambiguous mode"):
        annotate_hkl(annotated_axes(), make_indexed(), LABEL_BASE, ambiguous="both")


def test_annotate_hkl_of_nothing() -> None:
    assert annotate_hkl(annotated_axes(), [], LABEL_BASE) == []


def test_mark_peaks_returns_one_text_per_position() -> None:
    positions = [22.0, 45.5, 71.25]
    texts = mark_peaks(annotated_axes(), positions, LABEL_BASE)

    assert len(texts) == len(positions)
    assert all(isinstance(text, Text) for text in texts)
    assert [text.get_position()[0] for text in texts] == pytest.approx(positions)
    assert all(text.get_position()[1] == pytest.approx(LABEL_BASE) for text in texts)
    assert all(text.get_text() == PEAK_MARKER for text in texts)


def test_mark_peaks_takes_a_custom_marker() -> None:
    texts = mark_peaks(annotated_axes(), [30.0], LABEL_BASE, marker="x", fontsize=11)

    assert texts[0].get_text() == "x"
    assert texts[0].get_fontsize() == 11


def test_mark_peaks_of_nothing() -> None:
    assert mark_peaks(annotated_axes(), [], LABEL_BASE) == []


def test_plot_stacked_returns_a_base_per_slot() -> None:
    scans = [make_scan("a"), make_scan("b"), make_scan("c")]
    fig, ax, bases = plot_stacked(scans, offset=2.0)

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert bases == pytest.approx([0.0, 2.0, 4.0])


def test_annotate_hkl_onto_a_chosen_stack_slot() -> None:
    scans = [make_scan("a"), make_scan("b")]
    _, ax, bases = plot_stacked(scans, offset=2.0)

    texts = annotate_hkl(ax, make_indexed(), bases[1] + 1.2, label_height=0.1)

    assert texts[0].get_position()[1] == pytest.approx(bases[1] + 1.2)


def test_hkl_label_runs_single_digit_indices_together() -> None:
    indexed = [make_indexed_peak(FAR_PEAK, (3, 1, 1))]

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE)

    assert texts[0].get_text() == "311"
    assert HKL_THIN_SPACE not in texts[0].get_text()


def test_hkl_label_spaces_a_two_digit_index() -> None:
    indexed = [make_indexed_peak(FAR_PEAK, (10, 0, 2))]

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE)

    # Every index is spaced, not just the wide one, so the label reads evenly.
    assert texts[0].get_text() == f"10{HKL_THIN_SPACE}0{HKL_THIN_SPACE}2"
    # The compact form would be indistinguishable from (1 0 0) with a spare 2.
    assert texts[0].get_text() != "1002"


def test_hkl_label_switches_at_ten() -> None:
    nine = [make_indexed_peak(FAR_PEAK, (9, 5, 1))]
    ten = [make_indexed_peak(FAR_PEAK, (9, 5, 10))]

    assert annotate_hkl(annotated_axes(), nine, LABEL_BASE)[0].get_text() == "951"
    assert (
        annotate_hkl(annotated_axes(), ten, LABEL_BASE)[0].get_text()
        == f"9{HKL_THIN_SPACE}5{HKL_THIN_SPACE}10"
    )


# A reference peak with three weaker ones at these offsets, in degrees, to sit
# either side of whatever separation the geometry works out to.
REFERENCE_PEAK = 40.0
OFFSETS = (0.5, 1.0, 3.0)

LABEL_FONT = 7


def sized_axes(width_inches: float) -> Axes:
    """Bare axes of a known figure width and x range, ready to annotate."""
    fig = Figure(figsize=(width_inches, 3.0))
    ax = fig.add_subplot()
    ax.set_xlim(10.0, 100.0)
    return ax


def upright_separation(ax: Axes, fontsize: float = LABEL_FONT) -> float:
    """The separation annotate_hkl should work out for a 90 degree rotation."""
    width_points = ax.get_position().width * ax.figure.get_figwidth() * 72.0
    low, high = ax.get_xlim()
    return (
        LABEL_GAP_FACTOR
        * UPRIGHT_LABEL_WIDTH
        * fontsize
        * abs(high - low)
        / width_points
    )


def make_spaced_indexed() -> list[IndexedPeak]:
    """One strong peak with three weaker ones at OFFSETS above it."""
    peaks = [make_indexed_peak(REFERENCE_PEAK, (3, 1, 1), relative_intensity=100.0)]
    hkls = ((4, 2, 0), (5, 5, 0), (6, 3, 0))
    for offset, hkl, intensity in zip(OFFSETS, hkls, (50.0, 40.0, 30.0)):
        peaks.append(
            make_indexed_peak(
                REFERENCE_PEAK + offset, hkl, relative_intensity=intensity
            )
        )
    return peaks


def labelled_offsets(texts: list[Text]) -> list[float]:
    """The offsets from REFERENCE_PEAK that kept a label, nearest first."""
    return sorted(
        round(text.get_position()[0] - REFERENCE_PEAK, 3)
        for text in texts
        if text.get_position()[0] != REFERENCE_PEAK
    )


def test_min_separation_is_computed_from_the_figure_geometry() -> None:
    ax = sized_axes(3.5)
    separation = upright_separation(ax)

    texts = annotate_hkl(ax, make_spaced_indexed(), LABEL_BASE, fontsize=LABEL_FONT)

    # The strongest peak always claims its room first.
    assert REFERENCE_PEAK in [text.get_position()[0] for text in texts]
    # Exactly those far enough from it to clear a label of this size.
    assert labelled_offsets(texts) == [
        offset for offset in OFFSETS if offset > separation
    ]


def test_a_wider_figure_computes_a_smaller_separation_and_labels_more() -> None:
    narrow, wide = sized_axes(3.5), sized_axes(7.0)

    assert upright_separation(wide) < upright_separation(narrow)

    narrow_texts = annotate_hkl(
        narrow, make_spaced_indexed(), LABEL_BASE, fontsize=LABEL_FONT
    )
    wide_texts = annotate_hkl(
        wide, make_spaced_indexed(), LABEL_BASE, fontsize=LABEL_FONT
    )

    assert len(wide_texts) > len(narrow_texts)
    assert labelled_offsets(wide_texts) == [
        offset for offset in OFFSETS if offset > upright_separation(wide)
    ]


def test_a_narrower_x_range_computes_a_smaller_separation() -> None:
    ax = sized_axes(3.5)
    wide_range = upright_separation(ax)
    ax.set_xlim(30.0, 50.0)

    assert upright_separation(ax) < wide_range


def test_an_explicit_min_separation_overrides_the_computed_one() -> None:
    ax = sized_axes(3.5)
    # The computed value would drop all three; this keeps every one of them.
    assert upright_separation(ax) > max(OFFSETS)

    texts = annotate_hkl(
        ax,
        make_spaced_indexed(),
        LABEL_BASE,
        fontsize=LABEL_FONT,
        min_separation=0.4,
    )

    assert labelled_offsets(texts) == list(OFFSETS)


def test_a_label_that_is_not_upright_is_spaced_by_its_length() -> None:
    upright = annotate_hkl(
        sized_axes(7.0), make_spaced_indexed(), LABEL_BASE, fontsize=LABEL_FONT
    )
    flat = annotate_hkl(
        sized_axes(7.0),
        make_spaced_indexed(),
        LABEL_BASE,
        fontsize=LABEL_FONT,
        rotation=0,
    )

    # A three character label lying flat is wider than the same label upright,
    # so it needs more room and fewer of them fit.
    assert CHARACTER_WIDTH * 3 > UPRIGHT_LABEL_WIDTH
    assert len(flat) < len(upright)


def test_the_defaults_are_the_standing_annotation_rule() -> None:
    """A single row of major peaks, each labelled with its assignment."""
    assert HKL_MAX_LEVELS == 1
    assert HKL_MIN_RELATIVE_INTENSITY == 10.0
    assert HKL_AMBIGUOUS == "first"
    assert HKL_FONTSIZE == 7


def test_the_defaults_put_every_label_in_one_row() -> None:
    ax = sized_axes(7.0)
    separation = upright_separation(ax, HKL_FONTSIZE)

    texts = annotate_hkl(ax, make_spaced_indexed(), LABEL_BASE)

    # Crowded peaks are dropped, never raised into a second row.
    assert len(texts) < len(make_spaced_indexed())
    assert all(text.get_position()[1] == pytest.approx(LABEL_BASE) for text in texts)
    assert labelled_offsets(texts) == [
        offset for offset in OFFSETS if offset > separation
    ]


def test_the_default_intensity_cut_off_drops_a_minor_peak() -> None:
    indexed = make_indexed()
    # Above the old 5 per cent default, below the 10 per cent one.
    indexed[2].peak.relative_intensity = 7.0

    texts = annotate_hkl(annotated_axes(), indexed, LABEL_BASE, max_levels=2)

    assert [text.get_text() for text in texts] == ["311", "420"]
