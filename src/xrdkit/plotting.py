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

# Two hkl labels closer than this in degrees would overlap, so the second is
# stepped up above the first.
HKL_MIN_SEPARATION = 0.6

# One step up, as a fraction of the y range of the axes.
HKL_LABEL_HEIGHT = 0.04

# How annotate_hkl treats a peak that matched more than one reflection.
AMBIGUOUS_MODES = ("first", "all", "skip")

# Joins the candidate labels of an ambiguous peak under ambiguous="all".
AMBIGUOUS_SEPARATOR = "/"

# Written above a peak that the cell does not account for.
PEAK_MARKER = "*"
PEAK_MARKER_FONTSIZE = 9


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
    """Return the Miller indices run together, so ``(3 1 1)`` becomes ``311``."""
    return f"{reflection.h}{reflection.k}{reflection.l}"


def annotate_hkl(
    ax: Axes,
    indexed: list[IndexedPeak],
    y: float,
    min_relative_intensity: float = HKL_MIN_RELATIVE_INTENSITY,
    fontsize: float = HKL_FONTSIZE,
    rotation: float = HKL_ROTATION,
    min_separation: float = HKL_MIN_SEPARATION,
    label_height: float | None = None,
    ambiguous: str = "first",
) -> list[Text]:
    """Write an hkl label above each indexed peak of one trace.

    Labels are placed at the observed peak position with their base at ``y``,
    so ``y`` is normally the base of the slot the trace occupies plus enough
    room to clear its tallest peak. Peaks with no assignment are skipped; use
    :func:`mark_peaks` for those.

    A label within ``min_separation`` of the one before is stepped up by
    ``label_height``, each further label in a crowded run rising another step,
    and the first label with room of its own returning to ``y``. This is
    deliberately simple: it works from 2theta alone and never measures the
    rendered text, so a long label in a small figure can still collide.

    Parameters
    ----------
    ax
        Axes to write on.
    indexed
        Result of :func:`~xrdkit.indexing.index_peaks` for the trace.
    y
        Base of the unstepped labels, in data coordinates.
    min_relative_intensity
        Skip peaks below this percentage of the strongest peak of the scan.
    fontsize, rotation
        Label appearance. The rotation is in degrees, anticlockwise.
    min_separation
        Labels closer together than this, in degrees, are stepped up.
    label_height
        Height of one step, in data coordinates. Defaults to
        ``HKL_LABEL_HEIGHT`` times the y range of ``ax``.
    ambiguous
        What to do with a peak that matched more than one reflection:
        ``"first"`` labels the assigned one, ``"all"`` joins every candidate
        with a solidus, ``"skip"`` leaves the peak unlabelled.

    Returns
    -------
    list[Text]
        The labels written, in 2theta order.

    Raises
    ------
    ValueError
        If ``ambiguous`` is not one of ``AMBIGUOUS_MODES``.
    """
    if ambiguous not in AMBIGUOUS_MODES:
        raise ValueError(
            f"Unknown ambiguous mode {ambiguous!r}; expected one of {AMBIGUOUS_MODES}"
        )

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
    entries.sort(key=lambda entry: entry.peak.two_theta)

    texts: list[Text] = []
    previous_x: float | None = None
    previous_y = y
    for entry in entries:
        position = entry.peak.two_theta
        if ambiguous == "all" and len(entry.candidates) > 1:
            label = AMBIGUOUS_SEPARATOR.join(
                _hkl_label(reflection) for reflection in entry.candidates
            )
        else:
            label = _hkl_label(entry.reflection)

        # A crowded run climbs one step at a time; the first label with room of
        # its own drops back to the baseline.
        if previous_x is not None and position - previous_x < min_separation:
            height = previous_y + label_height
        else:
            height = y

        texts.append(
            ax.text(
                position,
                height,
                label,
                ha="center",
                va="bottom",
                fontsize=fontsize,
                rotation=rotation,
                rotation_mode="anchor",
            )
        )
        previous_x = position
        previous_y = height

    return texts


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
