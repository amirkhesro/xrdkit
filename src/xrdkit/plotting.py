"""Publication quality plotting for X-ray diffraction patterns."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.transforms import blended_transform_factory

from xrdkit.io import XRDScan

__all__ = ["apply_style", "plot_pattern", "plot_stacked", "save_figure"]

SCALES = ("linear", "sqrt", "log")

X_LABEL = "2θ (degrees)"
Y_LABEL = "Intensity (arb. units)"

SINGLE_COLUMN = (3.5, 2.6)


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
    """Apply the shared pattern axis labelling and limits."""
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_xlim(x_min, x_max)
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
) -> tuple[Figure, Axes]:
    """Plot several patterns stacked vertically with a constant offset.

    Parameters
    ----------
    scans
        Scans to stack, drawn bottom to top in the order given.
    labels
        One label per scan, written at the upper right of each trace. Defaults
        to each scan's ``sample_id``.
    offset
        Vertical spacing between traces. Defaults to 1.1 times the tallest
        scaled trace.
    scale, normalise, colour, linewidth
        As for :func:`plot_pattern`.
    figsize
        Figure size in inches.

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
        offset = 1.1 * max(float(np.max(trace)) for trace in traces)

    fig = Figure(figsize=figsize)
    ax = fig.add_subplot()
    # Labels sit a fixed inset from the right edge, at each trace's own height.
    label_transform = blended_transform_factory(ax.transAxes, ax.transData)

    for index, (scan, trace, text) in enumerate(zip(scans, traces, labels)):
        base = index * offset
        ax.plot(
            scan.two_theta, trace + base, color=colour, linewidth=linewidth, label=text
        )
        ax.text(
            0.98,
            base + float(np.max(trace)),
            text,
            transform=label_transform,
            ha="right",
            va="top",
            fontsize=8,
        )

    x_min = min(float(scan.two_theta[0]) for scan in scans)
    x_max = max(float(scan.two_theta[-1]) for scan in scans)
    _format_axes(ax, x_min, x_max)
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
