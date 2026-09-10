"""Reusable X-ray diffraction analysis toolkit for electroceramics research."""

from xrdkit.io import XRDScan, read_xrdml
from xrdkit.peaks import Peak, find_peaks, peaks_to_csv
from xrdkit.plotting import apply_style, plot_pattern, plot_stacked, save_figure

__all__ = [
    "Peak",
    "XRDScan",
    "apply_style",
    "find_peaks",
    "peaks_to_csv",
    "plot_pattern",
    "plot_stacked",
    "read_xrdml",
    "save_figure",
]
