"""Reusable X-ray diffraction analysis toolkit for electroceramics research."""

from xrdkit.density import (
    ATOMIC_MASSES,
    cell_volume,
    formula_mass,
    theoretical_density,
)
from xrdkit.indexing import (
    TTB_CELL,
    CellFit,
    IndexedPeak,
    Reflection,
    TetragonalCell,
    ZeroSearch,
    estimate_zero_offset,
    generate_reflections,
    index_and_refine,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
    refine_cell,
)
from xrdkit.io import XRDScan, read_xrdml
from xrdkit.lattice import LatticeFit, lattice_fit_to_dict, refine_lattice
from xrdkit.peaks import (
    Peak,
    exclude_kalpha2,
    find_peaks,
    flag_kalpha2,
    peaks_to_csv,
)
from xrdkit.plotting import (
    annotate_hkl,
    apply_style,
    mark_peaks,
    plot_pattern,
    plot_stacked,
    save_figure,
)

__all__ = [
    "ATOMIC_MASSES",
    "TTB_CELL",
    "CellFit",
    "IndexedPeak",
    "LatticeFit",
    "Peak",
    "Reflection",
    "TetragonalCell",
    "XRDScan",
    "ZeroSearch",
    "annotate_hkl",
    "apply_style",
    "cell_volume",
    "estimate_zero_offset",
    "exclude_kalpha2",
    "find_peaks",
    "flag_kalpha2",
    "formula_mass",
    "generate_reflections",
    "index_and_refine",
    "index_peaks",
    "indexed_to_csv",
    "indexing_summary",
    "lattice_fit_to_dict",
    "mark_peaks",
    "peaks_to_csv",
    "plot_pattern",
    "plot_stacked",
    "read_xrdml",
    "refine_cell",
    "refine_lattice",
    "save_figure",
    "theoretical_density",
]
