"""Publication quality plotting for X-ray diffraction patterns."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator

from xrdkit.indexing import IndexedPeak, Reflection
from xrdkit.io import XRDScan

__all__ = [
    "annotate_hkl",
    "apply_style",
    "mark_peaks",
    "plot_pattern",
    "plot_stacked",
    "save_figure",
]

SCALES = ("linear", "sqrt", "log")

X_LABEL = "2θ (degrees)"
Y_LABEL = "Intensity (arb. units)"

SINGLE_COLUMN = (3.5, 2.6)

X_MAJOR_TICK = 10.0
X_MINOR_TICK = 2.0

# Automatic stack spacing, as a multiple of the tallest scaled trace.
OFFSET_FACTOR = 1.2

# Stack labels sit this far below the top of their slot, and this fraction of
# the x range in from the right-hand end of the data.
LABEL_HEIGHT = 0.92
LABEL_MARGIN = 0.02

# Peaks weaker than this percentage of the strongest carry no hkl label, since
# a crowded pattern is unreadable if every shoulder is annotated.
HKL_MIN_RELATIVE_INTENSITY = 5.0

HKL_FONTSIZE = 7

# Upright labels take the least horizontal room, which is what a crowded 2theta
# axis is short of.
HKL_ROTATION = 90

# Fallback separation, in degrees, for a caller that wants a fixed one rather
# than the value annotate_hkl works out from the rendered label size.
HKL_MIN_SEPARATION = 0.6

# Clear space left between two labels, as a multiple of a label's own width.
LABEL_GAP_FACTOR = 1.25

# A label rotated 90 degrees is about this many font sizes wide across the x
# direction, whatever it says, since only its height runs along the axis.
UPRIGHT_LABEL_WIDTH = 1.2

# An unrotated character is about this fraction of the font size wide, so a
# label lying flat is this times its length times the font size.
CHARACTER_WIDTH = 0.6

# One level up, as a fraction of the y range of the axes.
HKL_LABEL_HEIGHT = 0.04

# How many levels of labels to try before giving up on a peak. One keeps the
# annotation to a single row, which reads best when the labels are meant to be
# scanned along the axis rather than picked out one at a time.
HKL_MAX_LEVELS = 1

# How annotate_hkl treats a peak that matched more than one reflection.
AMBIGUOUS_MODES = ("first", "all", "skip")

# Joins the candidate labels of an ambiguous peak under ambiguous="all".
AMBIGUOUS_SEPARATOR = "/"

# Written above a peak that the cell does not account for.
PEAK_MARKER = "*"
PEAK_MARKER_FONTSIZE = 9

# An index of 10 or more makes the compact form ambiguous, since (10 0 2) and
# (1 0 0) both run together as 1002. Those labels are spaced with a thin space,
# which separates the indices without opening the gap a full space would.
HKL_THIN_SPACE = "\u2009"
HKL_WIDE_INDEX = 10


def apply_style() -> None:
    """Set global rcParams to a clean single-column journal style."""
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
) -> tuple[Figure, Axes, list[float]]:
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
    tuple[Figure, Axes, list[float]]
        The figure, the axes, and the vertical base of each slot in the order
        the scans were given. The bases are what :func:`annotate_hkl` and
        :func:`mark_peaks` need to place labels against a chosen trace.

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
    for index, (scan, trace, text) in enumerate(zip(scans, traces, labels)):
        base = index * offset
        bases.append(base)
        ax.plot(
            scan.two_theta, trace + base, color=colour, linewidth=linewidth, label=text
        )
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
    return fig, ax, bases


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


def _label_separation(
    ax: Axes, fontsize: float, rotation: float, labels: list[str]
) -> float:
    """Return the 2theta two labels need between them so they cannot overlap.

    The axes are measured as they stand, so the figure size and the x limits
    have to be settled before this is called.
    """
    # get_position is a fraction of the figure, so this is the drawn width of
    # the axes in points.
    width_points = ax.get_position().width * ax.figure.get_figwidth() * 72.0
    low, high = ax.get_xlim()
    degrees_per_point = abs(high - low) / width_points

    if rotation == 90:
        # Upright, so the label's width is its font size whatever it says.
        label_points = UPRIGHT_LABEL_WIDTH * fontsize
    else:
        # At any other angle the longest label sets the spacing for all of them.
        longest = max((len(label) for label in labels), default=1)
        label_points = longest * CHARACTER_WIDTH * fontsize

    return LABEL_GAP_FACTOR * label_points * degrees_per_point


def annotate_hkl(
    ax: Axes,
    indexed: list[IndexedPeak],
    y: float,
    min_relative_intensity: float = HKL_MIN_RELATIVE_INTENSITY,
    fontsize: float = HKL_FONTSIZE,
    rotation: float = HKL_ROTATION,
    min_separation: float | None = None,
    label_height: float | None = None,
    ambiguous: str = "first",
    max_levels: int = HKL_MAX_LEVELS,
) -> list[Text]:
    """Write an hkl label above the indexed peaks of one trace, strongest first.

    Labels are placed at the observed peak position with their base at ``y``,
    so ``y`` is normally the base of the slot the trace occupies plus enough
    room to clear its tallest peak. Peaks with no assignment are skipped; use
    :func:`mark_peaks` for those.

    Room is allotted by priority rather than by position. The peaks are taken
    strongest first, and each is given the lowest level on which no label has
    yet been placed within ``min_separation`` degrees of it. Level 0 has its
    base at ``y`` and every level above it is raised by another
    ``label_height``. **A peak that finds no free level keeps no label**, so in
    a crowded stretch it is the weaker peaks that lose theirs and the strong
    reflections a reader is looking for stay labelled. With the default
    ``max_levels`` of 1 the annotation is a single row.

    By default ``min_separation`` is worked out from how much room a label
    actually takes: the width of the axes in points gives the degrees of 2theta
    per point, and a label rotated 90 degrees is about
    ``UPRIGHT_LABEL_WIDTH`` font sizes across, or ``CHARACTER_WIDTH`` times its
    length times the font size when it lies at any other angle. **The figure
    size and the x limits are read as they stand**, so set both before calling
    this; annotating and then resizing the figure or changing the limits will
    leave the labels spaced for the old geometry.

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
        degrees. ``None``, the default, works it out from the font size and the
        width of the axes so that labels cannot overlap.
    label_height
        Height of one level, in data coordinates. Defaults to
        ``HKL_LABEL_HEIGHT`` times the y range of ``ax``.
    ambiguous
        What to do with a peak that matched more than one reflection:
        ``"first"`` labels the assigned one, ``"all"`` joins every candidate
        with a solidus, ``"skip"`` leaves the peak unlabelled.
    max_levels
        How many levels to try. One keeps every label in a single row.

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

    if min_separation is None:
        min_separation = _label_separation(
            ax, fontsize, rotation, [label for _, label in labelled]
        )

    taken: list[list[float]] = [[] for _ in range(max_levels)]
    placed: list[tuple[float, Text]] = []
    for entry, label in labelled:
        position = entry.peak.two_theta
        level = next(
            (
                candidate
                for candidate in range(max_levels)
                if all(
                    abs(position - other) >= min_separation
                    for other in taken[candidate]
                )
            ),
            None,
        )
        if level is None:
            continue

        taken[level].append(position)
        placed.append(
            (
                position,
                ax.text(
                    position,
                    y + level * label_height,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=fontsize,
                    rotation=rotation,
                    rotation_mode="anchor",
                ),
            )
        )

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

    Returns
    -------
    list[Text]
        One marker per position, in the order given.
    """
    return [
        ax.text(position, y, marker, ha="center", va="bottom", fontsize=fontsize)
        for position in positions
    ]


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
