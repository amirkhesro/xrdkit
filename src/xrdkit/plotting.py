"""Publication quality plotting for X-ray diffraction patterns."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import Bbox

from xrdkit.broadening import Caglioti
from xrdkit.indexing import IndexedPeak, Reflection
from xrdkit.io import XRDScan

__all__ = [
    "annotate_hkl",
    "apply_style",
    "mark_peaks",
    "plot_caglioti",
    "plot_pattern",
    "plot_rietveld",
    "plot_stacked",
    "save_figure",
]

SCALES = ("linear", "sqrt", "log")

X_LABEL = "2θ (degrees)"
Y_LABEL = "Intensity (arb. units)"
FWHM_LABEL = "FWHM (degrees)"

# Legend entries of plot_caglioti: the curve, the widths it was fitted to and
# those left out.
CAGLIOTI_LABELS = ("Caglioti fit", "Used in fit", "Excluded")

SINGLE_COLUMN = (3.5, 2.6)

X_MAJOR_TICK = 10.0
X_MINOR_TICK = 2.0

# Automatic stack spacing, as a multiple of the tallest scaled trace.
OFFSET_FACTOR = 1.2

# Stack labels sit this far below the top of their slot, and this fraction of
# the x range in from the right-hand end of the data.
LABEL_HEIGHT = 0.92
LABEL_MARGIN = 0.02

# The defaults of annotate_hkl, which together set the standing rule for hkl
# annotation: the labels go in a single row, only the major peaks are labelled,
# a weaker peak that would collide with a stronger one loses its label, and
# labels are never stacked into a column.
#
# Every peak of any visible size competes for room; which of them actually get
# a label is settled by what fits, not by this, so the cut-off only has to keep
# noise out of the running.
HKL_MIN_RELATIVE_INTENSITY = 5.0

HKL_FONTSIZE = 7

# Upright labels take the least horizontal room, which is what a crowded 2theta
# axis is short of.
HKL_ROTATION = 90

# One row. Raising this stacks labels into a column, which is off by default
# because a column drifts away from the peaks it names.
HKL_MAX_LEVELS = 1

# A peak that matched more than one reflection is labelled with the one it was
# assigned; see AMBIGUOUS_MODES for the alternatives.
HKL_AMBIGUOUS = "first"

# Fallback separation, in degrees, for a caller that wants a fixed one rather
# than the value annotate_hkl works out from the rendered label size.
HKL_MIN_SEPARATION = 0.6

# Clear space left between two neighbouring labels, in points.
LABEL_GAP_POINTS = 1.0

# Clear space left between a label and the top of its own peak, in points.
LABEL_PAD_POINTS = 2

# How far either side of a peak to look along a trace for its top, in degrees.
LABEL_PEAK_WINDOW = 0.15

# One level up, as a fraction of the y range of the axes.
HKL_LABEL_HEIGHT = 0.04

# How annotate_hkl treats a peak that matched more than one reflection.
AMBIGUOUS_MODES = ("first", "all", "skip")

# Joins the candidate labels of an ambiguous peak under ambiguous="all".
AMBIGUOUS_SEPARATOR = "/"

# plot_rietveld: the pattern columns the GSAS-II driver exports, the legend
# entries of the four curves, and how each is drawn. The package draws in
# black and white, so the observed points are grey to let the black
# calculated line show through them.
RIETVELD_COLUMNS = ("two_theta", "observed", "calculated", "background", "difference")
RIETVELD_LABELS = ("Observed", "Calculated", "Background", "Difference")
RIETVELD_FIGSIZE = (3.5, 3.0)
OBSERVED_STYLE = {
    "linestyle": "none",
    "marker": "o",
    "markersize": 1.6,
    "markerfacecolor": "none",
    "markeredgecolor": "0.45",
    "markeredgewidth": 0.4,
}
CALCULATED_STYLE = {"color": "black", "linewidth": 0.7}
BACKGROUND_STYLE = {"color": "0.55", "linewidth": 0.5}
DIFFERENCE_STYLE = {"color": "black", "linewidth": 0.5}

# Reflection ticks: one row per phase, the first black and the rest in
# lighter greys, each row this fraction of the pattern's height below the
# one above, with the difference curve a further row below the last.
TICK_COLOURS = ("black", "0.5", "0.7", "0.3")
TICK_MARKER_SIZE = 4.0
TICK_MARKER_WIDTH = 0.6
TICK_ROW_FRACTION = 0.06

# Room left above the pattern and below the difference curve, as a fraction
# of the whole height plotted.
RIETVELD_MARGIN = 0.03

# The block of fit statistics goes this many points below the legend.
STATISTICS_GAP_POINTS = 4.0
STATISTICS_FONTSIZE = 7

# The top of the y axis is raised until the pattern under the legend and the
# statistics stops this fraction of the axes height below them.
LEGEND_CLEARANCE = 0.02

# Two lattice parameters closer than this are taken as equal, and an angle
# this close to 90 degrees as a right angle, when the cell is written out.
CELL_TOLERANCE = 1e-6

# Written above a peak that the cell does not account for.
PEAK_MARKER = "*"
PEAK_MARKER_FONTSIZE = 9

# An index of 10 or more makes the compact form ambiguous, since (10 0 2) and
# (1 0 0) both run together as 1002. Those labels are spaced with a thin space,
# which separates the indices without opening the gap a full space would.
HKL_THIN_SPACE = "\u2009"
HKL_WIDE_INDEX = 10


def apply_style() -> None:
    """Set global rcParams to a clean single-column journal style.

    The change is to ``matplotlib.rcParams`` itself, so it holds for every
    figure the process draws afterwards and nothing here puts the old values
    back. Call it once, before anything is drawn; a script that also draws
    figures of its own should call this first and then override what it
    wants, or make its own figures inside ``matplotlib.rc_context``.
    """
    families = {font.name for font in mpl.font_manager.fontManager.ttflist}
    sans = ["Arial", "DejaVu Sans"] if "Arial" in families else ["DejaVu Sans"]

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": sans,
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            # Keep the full box: all four spines stay on.
            "axes.linewidth": 0.8,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "axes.spines.left": True,
            "axes.spines.bottom": True,
            "axes.grid": False,
            # Inward ticks on all four sides, minor ticks on x only.
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": False,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.6,
            "lines.linewidth": 0.7,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def _scaled(intensity: np.ndarray, scale: str, normalise: bool) -> np.ndarray:
    """Apply ``scale`` to ``intensity``, normalised to a maximum of 1 if asked.

    Normalisation is applied to the transformed values so the result always
    peaks at 1. For ``"linear"`` and ``"sqrt"`` that is identical to normalising
    the counts first; for ``"log"`` it is the only meaningful order, since
    normalised counts would all clip to the floor of 1 and vanish.
    """
    if scale not in SCALES:
        raise ValueError(f"Unknown scale {scale!r}; expected one of {SCALES}")

    values = np.asarray(intensity, dtype=float)
    if scale == "sqrt":
        values = np.sqrt(np.clip(values, 0.0, None))
    elif scale == "log":
        values = np.log10(np.clip(values, 1.0, None))

    if normalise:
        peak = float(np.max(values)) if values.size else 0.0
        if peak > 0.0:
            values = values / peak
    return values


def _format_axes(ax: Axes, x_min: float, x_max: float) -> None:
    """Apply the shared pattern axis labelling, ticks and limits."""
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_xlim(x_min, x_max)
    # Fixed 10 degree major ticks so a scan running to nearly 100 degrees is
    # still labelled at its far end, which the automatic locator does not do.
    ax.xaxis.set_major_locator(MultipleLocator(X_MAJOR_TICK))
    ax.xaxis.set_minor_locator(MultipleLocator(X_MINOR_TICK))
    # Intensities are arbitrary units, so the y scale carries no ticks.
    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)


def plot_pattern(
    scan: XRDScan,
    ax: Axes | None = None,
    scale: str = "linear",
    normalise: bool = False,
    label: str | None = None,
    colour: str = "black",
    linewidth: float = 0.7,
) -> tuple[Figure, Axes]:
    """Plot a single diffraction pattern.

    Parameters
    ----------
    scan
        The scan to plot.
    ax
        Axes to draw on. A new single-column figure is created if omitted.
    scale
        Intensity scale: ``"linear"``, ``"sqrt"`` or ``"log"``.
    normalise
        Scale the trace to a maximum of 1.
    label
        Legend label for the line.
    colour, linewidth
        Line appearance.

    Raises
    ------
    ValueError
        If ``scale`` is not a known scale.
    """
    values = _scaled(scan.intensity, scale, normalise)

    if ax is None:
        fig = Figure(figsize=SINGLE_COLUMN)
        ax = fig.add_subplot()
    else:
        fig = ax.figure

    ax.plot(scan.two_theta, values, color=colour, linewidth=linewidth, label=label)
    _format_axes(ax, float(scan.two_theta[0]), float(scan.two_theta[-1]))
    return fig, ax


def plot_stacked(
    scans: list[XRDScan],
    labels: list[str] | None = None,
    offset: float | None = None,
    scale: str = "linear",
    normalise: bool = True,
    colour: str = "black",
    linewidth: float = 0.7,
    figsize: tuple[float, float] = (3.5, 5.0),
) -> tuple[Figure, Axes, list[float], list[Line2D]]:
    """Plot several patterns stacked vertically with a constant offset.

    Parameters
    ----------
    scans
        Scans to stack, drawn bottom to top in the order given.
    labels
        One label per scan, written at the upper right of each trace. Defaults
        to each scan's ``sample_id``.
    offset
        Vertical spacing between traces. Defaults to 1.2 times the tallest
        scaled trace, which keeps the tallest peak clear of the trace above.
    scale, normalise, colour, linewidth
        As for :func:`plot_pattern`.
    figsize
        Figure size in inches.

    Returns
    -------
    tuple[Figure, Axes, list[float], list[Line2D]]
        The figure, the axes, the vertical base of each slot, and the line
        drawn for each scan, both in the order the scans were given. The bases
        and the lines are what :func:`annotate_hkl` and :func:`mark_peaks` need
        to place labels against a chosen trace.

    Raises
    ------
    ValueError
        If ``scans`` is empty, ``scale`` is unknown, or ``labels`` has a
        different length to ``scans``.
    """
    if not scans:
        raise ValueError("plot_stacked needs at least one scan")
    if labels is not None and len(labels) != len(scans):
        raise ValueError(
            f"Got {len(labels)} labels for {len(scans)} scans; lengths must match"
        )

    traces = [_scaled(scan.intensity, scale, normalise) for scan in scans]
    if labels is None:
        labels = [scan.sample_id for scan in scans]

    if offset is None:
        offset = OFFSET_FACTOR * max(float(np.max(trace)) for trace in traces)

    x_min = min(float(scan.two_theta[0]) for scan in scans)
    x_max = max(float(scan.two_theta[-1]) for scan in scans)
    label_x = x_max - LABEL_MARGIN * (x_max - x_min)

    fig = Figure(figsize=figsize)
    ax = fig.add_subplot()

    bases: list[float] = []
    lines: list[Line2D] = []
    for index, (scan, trace, text) in enumerate(zip(scans, traces, labels)):
        base = index * offset
        bases.append(base)
        drawn = ax.plot(
            scan.two_theta, trace + base, color=colour, linewidth=linewidth, label=text
        )
        lines.append(drawn[0])
        # Anchored to the slot, not to the trace's peak, so a strong high angle
        # reflection cannot push a label into the trace above it.
        ax.text(
            label_x,
            base + LABEL_HEIGHT * offset,
            text,
            ha="right",
            va="top",
            fontsize=8,
        )

    _format_axes(ax, x_min, x_max)
    # Every slot gets the same height, so the topmost label always has room.
    ax.set_ylim(min(0.0, float(np.min(traces[0]))), (len(scans) - 1) * offset + offset)
    return fig, ax, bases, lines


def _hkl_label(reflection: Reflection) -> str:
    """Return the Miller indices as a label, so ``(3 1 1)`` becomes ``311``.

    Indices below 10 run together. As soon as one reaches ``HKL_WIDE_INDEX`` the
    compact form stops being readable, so every index of that label is separated
    by a thin space instead: ``(10 0 2)`` becomes ``10 0 2`` rather than
    ``1002``, which would otherwise be indistinguishable from ``(1 0 0)`` with a
    trailing 2.
    """
    indices = reflection.hkl
    separator = HKL_THIN_SPACE if max(indices) >= HKL_WIDE_INDEX else ""
    return separator.join(str(index) for index in indices)


def _renderer_for(figure: Figure) -> object:
    """Return a renderer able to measure text on ``figure``.

    A figure built directly, rather than through a backend, carries a plain
    ``FigureCanvasBase``, which cannot rasterise and so offers no renderer at
    all. Attaching an Agg canvas gives one without going near pyplot.
    """
    canvas = figure.canvas
    if not hasattr(canvas, "get_renderer"):
        canvas = FigureCanvasAgg(figure)
    try:
        return canvas.get_renderer()
    except AttributeError:
        canvas.draw()
        return canvas.get_renderer()


def _write_label(
    ax: Axes,
    position: float,
    height: float,
    label: str,
    fontsize: float,
    rotation: float,
) -> Text:
    """Write one hkl label with its bottom centre exactly at ``(position, height)``.

    The alignment is applied to the text as rotated, which is what
    ``rotation_mode="default"`` means and what puts the middle of the drawn
    label over its peak at any angle. Aligning before rotating instead, with
    ``rotation_mode="anchor"``, swings a label upright about its own baseline
    and leaves it a third of a degree to the left of the peak it names.
    """
    return ax.text(
        position,
        height,
        label,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        rotation=rotation,
        rotation_mode="default",
    )


def _peak_top(
    line: Line2D, position: float, window: float = LABEL_PEAK_WINDOW
) -> float:
    """Return the highest point of ``line`` within ``window`` of ``position``.

    Falls back to the nearest point it has when nothing lies inside the window,
    so a label is never stranded at the foot of the axes.
    """
    x = np.asarray(line.get_xdata(), dtype=float)
    values = np.asarray(line.get_ydata(), dtype=float)
    near = np.abs(x - position) <= window
    if not np.any(near):
        return float(values[np.argmin(np.abs(x - position))])
    return float(np.max(values[near]))


def _points_in_data(ax: Axes, points: float) -> float:
    """Return ``points`` as a vertical distance in the data units of ``ax``."""
    pixels = points * ax.figure.dpi / 72.0
    inverse = ax.transData.inverted()
    return inverse.transform((0.0, pixels))[1] - inverse.transform((0.0, 0.0))[1]


def annotate_hkl(
    ax: Axes,
    indexed: list[IndexedPeak],
    y: float = 0.0,
    min_relative_intensity: float = HKL_MIN_RELATIVE_INTENSITY,
    fontsize: float = HKL_FONTSIZE,
    rotation: float = HKL_ROTATION,
    min_separation: float | None = None,
    label_height: float | None = None,
    ambiguous: str = HKL_AMBIGUOUS,
    max_levels: int = HKL_MAX_LEVELS,
    line: Line2D | None = None,
) -> list[Text]:
    """Write an hkl label above the indexed peaks of one trace, strongest first.

    By default the labels go in a single row, only the major peaks are
    labelled, a weaker peak that would collide with a stronger one loses its
    label, and labels are never stacked into a column.

    Given a ``line``, each label rides on top of its own peak, a couple of
    points above the highest the trace reaches nearby, which keeps a label and
    the reflection it names together however the pattern rises and falls.
    Without one they share a fixed row with their base at ``y``, normally the
    base of the slot the trace occupies plus enough room to clear its tallest
    peak. Either way the bottom centre of the label sits exactly over the
    observed peak position. Peaks with no assignment are skipped; use
    :func:`mark_peaks` for those.

    Room is allotted by priority rather than by position. The peaks are taken
    strongest first, and each is given the lowest level on which no label has
    yet been placed within ``min_separation`` degrees of it. Level 0 has its
    base at ``y`` and every level above it is raised by another
    ``label_height``. **A peak that finds no free level keeps no label**, so in
    a crowded stretch it is the weaker peaks that lose theirs and the strong
    reflections a reader is looking for stay labelled. With the default
    ``max_levels`` of 1 the annotation is a single row.

    By default nothing is estimated: each label is written, its rendered
    bounding box measured, and kept only if that box, widened by
    ``LABEL_GAP_POINTS``, clears every box already placed. Riding the peaks
    puts the labels at all sorts of heights, so with a ``line`` a box is
    compared with every other and two labels a tenth of a degree apart both
    stay as long as their peaks differ enough in height; without one the
    comparison is per level. This packs
    the labels as tightly as the text really allows, whatever they say and at
    whatever angle. **The figure size and the x limits are read as they
    stand**, so set both before calling this; annotating and then resizing the
    figure or changing the limits will leave the labels where the old geometry
    put them.

    Parameters
    ----------
    ax
        Axes to write on.
    indexed
        Result of :func:`~xrdkit.indexing.index_peaks` for the trace.
    y
        Base of the level 0 labels, in data coordinates.
    min_relative_intensity
        Skip peaks below this percentage of the strongest peak of the scan.
    fontsize, rotation
        Label appearance. The rotation is in degrees, anticlockwise.
    min_separation
        Two labels on the same level must be at least this far apart, in
        degrees. ``None``, the default, measures the rendered labels instead
        and keeps whichever ones do not overlap.
    label_height
        Height of one level, in data coordinates. Defaults to
        ``HKL_LABEL_HEIGHT`` times the y range of ``ax``.
    ambiguous
        What to do with a peak that matched more than one reflection:
        ``"first"`` labels the assigned one, ``"all"`` joins every candidate
        with a solidus, ``"skip"`` leaves the peak unlabelled.
    max_levels
        How many levels to try. One keeps every label in a single row. Ignored
        when ``line`` is given, since a label then sits on its own peak rather
        than on a level.
    line
        The trace the peaks belong to, as returned by :func:`plot_stacked`. Its
        y data is what each label is stood on. ``None`` puts them all in a row
        at ``y`` instead. Giving one always measures, so ``min_separation`` has
        no effect.

    Returns
    -------
    list[Text]
        The labels written, in 2theta order. Peaks that found no free level are
        not represented, so this can be shorter than the peaks given.

    Raises
    ------
    ValueError
        If ``ambiguous`` is not one of ``AMBIGUOUS_MODES``, or ``max_levels`` is
        below one.
    """
    if ambiguous not in AMBIGUOUS_MODES:
        raise ValueError(
            f"Unknown ambiguous mode {ambiguous!r}; expected one of {AMBIGUOUS_MODES}"
        )
    if max_levels < 1:
        raise ValueError(f"Need at least one level, got {max_levels}")

    if label_height is None:
        low, high = ax.get_ylim()
        label_height = HKL_LABEL_HEIGHT * (high - low)

    entries = [
        entry
        for entry in indexed
        if entry.is_indexed
        and entry.peak.relative_intensity >= min_relative_intensity
        and not (ambiguous == "skip" and len(entry.candidates) > 1)
    ]
    # Strongest first, so the peaks that matter claim their room before the
    # weak ones. Ties fall back on position to keep the choice reproducible.
    entries.sort(
        key=lambda entry: (-entry.peak.relative_intensity, entry.peak.two_theta)
    )

    labelled = [
        (
            entry,
            AMBIGUOUS_SEPARATOR.join(
                _hkl_label(reflection) for reflection in entry.candidates
            )
            if ambiguous == "all" and len(entry.candidates) > 1
            else _hkl_label(entry.reflection),
        )
        for entry in entries
    ]

    # Standing the labels on their peaks always measures, since their heights
    # are then all different and only the boxes can say what really clashes.
    measuring = line is not None or min_separation is None
    renderer = _renderer_for(ax.figure) if measuring else None
    # A point is this many pixels, and the gap is left on each side of the box.
    gap = LABEL_GAP_POINTS * ax.figure.dpi / 72.0 if measuring else 0.0
    pad = _points_in_data(ax, LABEL_PAD_POINTS) if line is not None else 0.0

    # One list of levels to try, or a single pass when each label has its own
    # height and there are no levels to speak of.
    levels = 1 if line is not None else max_levels
    boxes: list[list[Bbox]] = [[] for _ in range(levels)]
    taken: list[list[float]] = [[] for _ in range(levels)]
    placed: list[tuple[float, Text]] = []

    for entry, label in labelled:
        position = entry.peak.two_theta
        for level in range(levels):
            if line is not None:
                height = _peak_top(line, position) + pad
            else:
                height = y + level * label_height

            if not measuring:
                if any(
                    abs(position - other) < min_separation for other in taken[level]
                ):
                    continue
                taken[level].append(position)
                placed.append(
                    (
                        position,
                        _write_label(ax, position, height, label, fontsize, rotation),
                    )
                )
                break

            # Written before it is measured, since only a real Text knows how
            # much room it takes, and taken away again if it does not fit.
            text = _write_label(ax, position, height, label, fontsize, rotation)
            box = text.get_window_extent(renderer)
            widened = Bbox.from_extents(box.x0 - gap, box.y0, box.x1 + gap, box.y1)
            if any(widened.overlaps(other) for other in boxes[level]):
                text.remove()
                continue
            boxes[level].append(box)
            placed.append((position, text))
            break

    placed.sort(key=lambda item: item[0])
    return [text for _, text in placed]


def mark_peaks(
    ax: Axes,
    positions: list[float],
    y: float,
    marker: str = PEAK_MARKER,
    fontsize: float = PEAK_MARKER_FONTSIZE,
) -> list[Text]:
    """Write ``marker`` above each 2theta in ``positions``, with its base at ``y``.

    For the peaks a cell does not account for, whether unindexed or from a
    second phase, which are worth pointing at even though they carry no hkl.

    Nothing is checked and nothing is placed: every position given gets a
    marker wherever it falls, inside the axes limits or outside them, and
    markers that land on top of one another stay. Unlike
    :func:`annotate_hkl`, which measures and drops what will not fit, this
    writes what it is told.

    Returns
    -------
    list[Text]
        One marker per position, in the order given.
    """
    return [
        ax.text(position, y, marker, ha="center", va="bottom", fontsize=fontsize)
        for position in positions
    ]


def plot_caglioti(
    fit: Caglioti,
    two_theta: np.ndarray,
    fwhm: np.ndarray,
    esd: np.ndarray | None = None,
    included: np.ndarray | None = None,
    ax: Axes | None = None,
    two_theta_range: tuple[float, float] | None = None,
    labels: tuple[str, str, str] = CAGLIOTI_LABELS,
) -> tuple[Figure, Axes]:
    """Plot measured peak widths against 2theta with a fitted Caglioti curve.

    Widths used in the fit are drawn as filled circles and the rest as open
    ones, so the two can be told apart without colour. With other ``labels``
    the same figure serves for sample widths against the instrumental curve,
    the open circles then marking whatever the caller singles out.

    Parameters
    ----------
    fit
        The Caglioti function to draw.
    two_theta, fwhm
        The measured positions and widths, in degrees.
    esd
        Esds of the widths, drawn as error bars when given.
    included
        Which widths the fit used. All of them by default.
    ax
        Axes to draw on. A new single-column figure is created if omitted.
    two_theta_range
        ``(low, high)`` x limits. By default the data span, widened out to
        whole major ticks.
    labels
        Legend entries for the curve, the filled points and the open points.
    """
    two_theta = np.asarray(two_theta, dtype=float)
    fwhm = np.asarray(fwhm, dtype=float)
    curve_label, filled_label, open_label = labels
    used = (
        np.ones(two_theta.size, dtype=bool)
        if included is None
        else np.asarray(included, dtype=bool)
    )

    if ax is None:
        fig = Figure(figsize=SINGLE_COLUMN)
        ax = fig.add_subplot()
    else:
        fig = ax.figure

    if two_theta_range is None:
        two_theta_range = (
            X_MAJOR_TICK * np.floor(two_theta.min() / X_MAJOR_TICK),
            X_MAJOR_TICK * np.ceil(two_theta.max() / X_MAJOR_TICK),
        )
    low, high = two_theta_range
    # The curve is undefined at 0 degrees, so it starts just above it.
    curve = np.linspace(max(low, 0.5), high, 400)
    ax.plot(curve, fit.fwhm(curve), color="black", label=curve_label)

    for mask, face, label in (
        (used, "black", filled_label),
        (~used, "white", open_label),
    ):
        if not np.any(mask):
            continue
        ax.errorbar(
            two_theta[mask],
            fwhm[mask],
            yerr=None if esd is None else np.asarray(esd, dtype=float)[mask],
            fmt="o",
            markersize=4,
            markerfacecolor=face,
            markeredgecolor="black",
            markeredgewidth=0.8,
            ecolor="black",
            elinewidth=0.6,
            capsize=1.5,
            label=label,
        )

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(FWHM_LABEL)
    ax.set_xlim(low, high)
    ax.xaxis.set_major_locator(MultipleLocator(X_MAJOR_TICK))
    ax.xaxis.set_minor_locator(MultipleLocator(X_MINOR_TICK))
    ax.legend(frameon=False)
    return fig, ax


def _read_table(table: object, columns: Sequence[str], what: str) -> list[np.ndarray]:
    """The named columns of a CSV file or of anything indexed by column name.

    A path is read as a CSV file with a header row. Anything else is indexed
    by the column names, which covers a dict of arrays, a numpy structured
    array and a pandas DataFrame alike without importing pandas.
    """
    if isinstance(table, (str, Path)):
        table = np.atleast_1d(np.genfromtxt(table, delimiter=",", names=True))
        names = table.dtype.names or ()
    else:
        names = None
    values = []
    for column in columns:
        if names is not None and column not in names:
            raise ValueError(f"{what} has no {column!r} column")
        try:
            values.append(np.asarray(table[column], dtype=float).ravel())
        except (KeyError, ValueError, IndexError) as error:
            raise ValueError(f"{what} has no {column!r} column") from error
    return values


def _format_esd(value: float, esd: float | None) -> str:
    """``value`` with its esd in brackets, as 4.15683(2).

    The esd keeps one significant figure, or two when it begins with a 1, and
    the value is rounded to match. With no esd, or a zero one for a value
    held fixed, the value is written as it stands.
    """
    if esd is None or not esd > 0.0 or not math.isfinite(esd):
        return f"{value:.10g}"
    decimals = -math.floor(math.log10(esd))
    if round(esd * 10.0**decimals) == 1:
        decimals += 1
    decimals = max(decimals, 0)
    digits = round(esd * 10.0**decimals)
    return f"{value:.{decimals}f}({digits})"


def _cell_text(cell: Mapping[str, float], esd: Mapping[str, float] | None) -> str:
    """The lattice parameters a reader needs, with the cell's symmetry used.

    Lengths that are equal are given once, as a for a cubic cell and a and c
    for a tetragonal or hexagonal one, and only angles other than 90 degrees
    are written.
    """
    esd = esd or {}
    a, b, c = (cell[f"length_{axis}"] for axis in "abc")
    if abs(a - b) < CELL_TOLERANCE and abs(b - c) < CELL_TOLERANCE:
        lengths = "a"
    elif abs(a - b) < CELL_TOLERANCE:
        lengths = "ac"
    else:
        lengths = "abc"
    parts = [
        f"{axis} = {_format_esd(cell[f'length_{axis}'], esd.get(f'length_{axis}'))} Å"
        for axis in lengths
    ]
    for key, symbol in (("alpha", "α"), ("beta", "β"), ("gamma", "γ")):
        angle = cell[f"angle_{key}"]
        if abs(angle - 90.0) > CELL_TOLERANCE:
            parts.append(f"{symbol} = {_format_esd(angle, esd.get(f'angle_{key}'))}°")
    return ", ".join(parts)


def _statistics_text(result: Mapping | str | Path) -> str:
    """Rwp, GOF and the refined cells of a GSAS-II refine result.

    The figures of merit are those of the last stage that succeeded, and the
    cells those of the result's final values.
    """
    if isinstance(result, (str, Path)):
        result = json.loads(Path(result).read_text(encoding="utf-8"))
    lines = []
    good = [stage for stage in result.get("stages", []) if "error" not in stage]
    if good:
        last = good[-1]
        lines.append(
            rf"$\mathregular{{R_{{wp}}}}$ = {last['rwp']:.2f}%,  GOF = {last['gof']:.2f}"
        )
    for phase in result.get("final", {}).get("phases", []):
        lines.append(
            f"{phase['name']}: {_cell_text(phase['cell'], phase.get('cell_esd'))}"
        )
    return "\n".join(lines)


def _below_legend(ax: Axes, legend: object, text: str) -> Text:
    """Write ``text`` right-aligned just below ``legend``, in axes coordinates."""
    renderer = _renderer_for(ax.figure)
    box = legend.get_window_extent(renderer)
    gap = STATISTICS_GAP_POINTS * ax.figure.dpi / 72.0
    x, y = ax.transAxes.inverted().transform((box.x1, box.y0 - gap))
    return ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=STATISTICS_FONTSIZE,
        linespacing=1.3,
    )


def _clear_of(
    ax: Axes, artists: list, two_theta: np.ndarray, curves: list[np.ndarray]
) -> None:
    """Raise the top of the y axis until ``curves`` pass under ``artists``.

    The legend and the statistics are placed in axes coordinates, so they
    keep their place in the frame as the y limits change while the pattern
    shrinks beneath them. Only the stretch of 2theta they span is considered,
    so a pattern whose tall peaks lie elsewhere is not squashed for nothing.
    """
    renderer = _renderer_for(ax.figure)
    box = Bbox.union([artist.get_window_extent(renderer) for artist in artists])
    low, high = sorted(
        ax.transData.inverted().transform([(box.x0, 0.0), (box.x1, 0.0)])[:, 0]
    )
    bottom = ax.transAxes.inverted().transform((box.x0, box.y0))[1] - LEGEND_CLEARANCE
    under = (two_theta >= low) & (two_theta <= high)
    if not np.any(under) or bottom <= 0.0:
        return
    peak = max(float(np.max(curve[under])) for curve in curves)
    y_min, y_max = ax.get_ylim()
    needed = y_min + (peak - y_min) / bottom
    if needed > y_max:
        ax.set_ylim(y_min, needed)


def plot_rietveld(
    pattern: object,
    reflections: Mapping[str, object] | Sequence[object] | None = None,
    title: str | None = None,
    phase_labels: Sequence[str] | None = None,
    ax: Axes | None = None,
    sqrt_scale: bool = False,
    difference_offset: float | None = None,
    result: Mapping | str | Path | None = None,
) -> tuple[Figure, Axes]:
    """Plot a Rietveld fit: observed, calculated, background and difference.

    The observed points are small open grey circles, the calculated pattern a
    black line through them and the background a thin grey line. Below the
    pattern each phase has a row of tick marks at its reflections, the first
    phase at the top, and below those the observed less the calculated
    pattern is drawn about a zero of its own. Every curve and tick row is
    named in the legend.

    Parameters
    ----------
    pattern
        The pattern as the GSAS-II driver exports it, with columns
        ``two_theta``, ``observed``, ``calculated``, ``background`` and
        ``difference``: the path of its CSV file, or a table indexed by
        column name such as a pandas DataFrame or a dict of arrays.
    reflections
        Reflection lists, each with a ``two_theta`` column, as paths or
        tables: a mapping of phase name to list, such as the driver result's
        ``exports["reflections"]``, or a sequence of lists. Reflections
        outside the pattern's range are left out.
    title
        Axes title.
    phase_labels
        Legend names of the phases, in the order of ``reflections``. By
        default the mapping's keys, or "Phase 1", "Phase 2" and so on for a
        sequence.
    ax
        Axes to draw on. A new figure is created if omitted.
    sqrt_scale
        Plot the square root of the counts, which brings up weak reflections.
        The difference is then that of the square roots, so that it is on the
        scale of the curves above it.
    difference_offset
        Height of the zero of the difference curve. By default it is put
        below zero and below the lowest tick row, clear of both.
    result
        The result of the driver's refine action, as a dict or the path of
        its JSON file. Given, Rwp and GOF of the last stage that succeeded and
        each phase's refined cell are written in a small block below the
        legend.

    Returns
    -------
    tuple[Figure, Axes]
        The figure and the axes drawn on. The whole fit shares one Axes: the
        curves, the tick rows and the difference are all drawn on it, so there
        is no further panel to return.

    Raises
    ------
    ValueError
        If a table lacks a column, or ``phase_labels`` does not have one name
        per reflection list.
    """
    two_theta, observed, calculated, background, difference = _read_table(
        pattern, RIETVELD_COLUMNS, "pattern"
    )
    if sqrt_scale:
        observed, calculated, background = (
            np.sqrt(np.clip(values, 0.0, None))
            for values in (observed, calculated, background)
        )
        difference = observed - calculated

    if reflections is None:
        tables, names = [], []
    elif isinstance(reflections, Mapping):
        tables, names = list(reflections.values()), [str(name) for name in reflections]
    else:
        tables = list(reflections)
        names = [f"Phase {index}" for index in range(1, len(tables) + 1)]
    if phase_labels is not None:
        if len(phase_labels) != len(tables):
            raise ValueError(
                f"Got {len(phase_labels)} phase labels for {len(tables)} "
                "reflection lists; lengths must match"
            )
        names = list(phase_labels)
    positions = [
        _read_table(table, ("two_theta",), f"reflections of {name}")[0]
        for table, name in zip(tables, names)
    ]

    if ax is None:
        fig = Figure(figsize=RIETVELD_FIGSIZE)
        ax = fig.add_subplot()
    else:
        fig = ax.figure

    low, high = float(np.min(two_theta)), float(np.max(two_theta))
    top = float(max(np.max(observed), np.max(calculated)))
    # The tick rows start below zero, or below the pattern should it dip
    # under zero, and step down by a fixed share of the pattern's height.
    base = min(0.0, float(min(np.min(observed), np.min(calculated))))
    step = TICK_ROW_FRACTION * (top - base)
    rows = [base - (index + 1) * step for index in range(len(positions))]
    if difference_offset is None:
        lowest = rows[-1] if rows else base
        difference_offset = lowest - step - float(np.max(difference))

    observed_label, calculated_label, background_label, difference_label = (
        RIETVELD_LABELS
    )
    # The background goes down before the calculated line, which then lies
    # over it where the two meet between the peaks.
    ax.plot(two_theta, observed, label=observed_label, **OBSERVED_STYLE)
    ax.plot(two_theta, background, label=background_label, **BACKGROUND_STYLE)
    ax.plot(two_theta, calculated, label=calculated_label, **CALCULATED_STYLE)
    ax.plot(
        two_theta,
        difference + difference_offset,
        label=difference_label,
        **DIFFERENCE_STYLE,
    )
    for index, (ticks, name, row) in enumerate(zip(positions, names, rows)):
        inside = ticks[(ticks >= low) & (ticks <= high)]
        ax.plot(
            inside,
            np.full(inside.shape, row),
            linestyle="none",
            marker="|",
            markersize=TICK_MARKER_SIZE,
            markeredgewidth=TICK_MARKER_WIDTH,
            color=TICK_COLOURS[index % len(TICK_COLOURS)],
            label=name,
        )

    _format_axes(ax, low, high)
    if sqrt_scale:
        ax.set_ylabel("√" + Y_LABEL)
    bottom = float(np.min(difference)) + difference_offset
    margin = RIETVELD_MARGIN * (top - bottom)
    ax.set_ylim(bottom - margin, top + margin)
    if title:
        ax.set_title(title)

    # Observed first in the legend, then the curves drawn over it.
    handles, texts = ax.get_legend_handles_labels()
    order = [texts.index(label) for label in RIETVELD_LABELS] + list(
        range(len(RIETVELD_LABELS), len(texts))
    )
    legend = ax.legend(
        [handles[i] for i in order],
        [texts[i] for i in order],
        loc="upper right",
        frameon=False,
        markerscale=2.0,
        handlelength=1.5,
        labelspacing=0.3,
        borderaxespad=0.4,
    )
    artists = [legend]
    if result is not None:
        text = _statistics_text(result)
        if text:
            artists.append(_below_legend(ax, legend, text))
    _clear_of(ax, artists, two_theta, [observed, calculated])
    return fig, ax


def save_figure(
    fig: Figure,
    path: str | Path,
    formats: tuple[str, ...] = ("png", "pdf"),
    dpi: int = 300,
) -> list[Path]:
    """Save ``fig`` once per format next to ``path``, returning the paths written.

    ``path`` is treated as a stem: an extension already naming one of ``formats``
    is replaced, so ``"pattern.png"`` and ``"pattern"`` behave the same. Other
    dots are kept, which matters for names like ``"x0.10_calcined"``.
    """
    path = Path(path)
    known = {fmt.lower().lstrip(".") for fmt in formats}
    stem = path.with_suffix("") if path.suffix.lower().lstrip(".") in known else path
    stem.parent.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for fmt in formats:
        suffix = fmt.lower().lstrip(".")
        out = stem.parent / f"{stem.name}.{suffix}"
        fig.savefig(out, format=suffix, dpi=dpi)
        written.append(out)
    return written
