"""The ``xrdkit`` command.

Each command is a subparser with a ``handler`` default that takes the parsed
arguments and returns the exit status. A new command adds an ``_add_<name>``
function and a call to it in :func:`build_parser`. A handler that cannot go on
raises :class:`CommandError`, which :func:`main` prints as one line, prefixed
with the command name, and turns into exit status 1.

A pattern named on the command line of check, plot or stack is an existing
file, used as it is, or else a sample key of the project file found from the
current folder (:func:`_resolve`); density takes a sample key. A sample brings
its scan, its instrument's wavelength and its first structure from the project,
and its outputs go to ``results/<command>/<key>`` under the project root unless
``--out`` is given. An option given on the command line wins over the project's
value.

lebail and rietveld take a sample key and run :func:`xrdkit.pipeline.run_sequence`
on it: lebail the Le Bail mode, rietveld the Rietveld modes from ``--from``
through ``--through``, each from the saved result of the mode before. Their
files go to ``results/lebail/<key>`` and ``results/rietveld/<key>``, or to the
folder ``--out`` names itself; a line is printed as each mode starts and for
each stage's outcome, and a mode that fails prints the path of its
``failure.md`` and returns 1.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import re
import sys
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from urllib.error import URLError

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from xrdkit import __version__
from xrdkit.broadening import KALPHA2_INTENSITY_RATIO, fit_profile
from xrdkit.cell import Cell
from xrdkit.density import (
    cell_volume,
    formula_mass,
    parse_formula,
    relative_density,
    theoretical_density,
)
from xrdkit.gsas2 import Gsas2Error, find_gsas2
from xrdkit.indexing import (
    DEFAULT_TOLERANCE,
    DEFAULT_ZERO_OFFSET,
    SUPPORTED_SPACE_GROUPS,
    CellFit,
    IndexedPeak,
    index_and_refine,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
)
from xrdkit.instrument import (
    REFINED_KEYS,
    WIDTH_WINDOW,
    fit_instrument_widths,
    kalpha2_wavelength,
    refine_instrument,
)
from xrdkit.io import XRDScan, read_scan
from xrdkit.lattice import LatticeFit, refine_lattice
from xrdkit.library import CELL_PARAMETERS, CRYSTAL_SYSTEMS, load_entry
from xrdkit.peaks import (
    KALPHA2_RATIO,
    Peak,
    exclude_kalpha2,
    find_peaks,
    flag_kalpha2,
    peaks_to_csv,
)
from xrdkit.phases import (
    PHASES_TOLERANCE,
    PHASES_WINDOW,
    MissingPhasesExtra,
    attribute_unexplained,
    cod_fetch,
    cod_search,
    fetch_candidates,
    observed_peaks,
    rank_candidates,
    require_phases_extra,
    write_cif_index,
)
from xrdkit.pipeline import MODES, Options, mode_paths, run_sequence
from xrdkit.plotting import (
    annotate_hkl,
    apply_style,
    plot_pattern,
    plot_stacked,
    save_figure,
)
from xrdkit.project import (
    FORMS,
    PROJECT_FILE,
    Project,
    Sample,
    StructureSpec,
    find_project,
    load_project,
    load_project_text,
    project_template,
    resolved_cell,
    resolved_z,
    results_dir,
    toml_string,
)
from xrdkit.quality import assess_scan, format_report

__all__ = ["CommandError", "build_parser", "main"]

SCALES = ("linear", "sqrt", "log")

# Room above the tallest point of a trace, as a fraction of the y axis span,
# for the hkl label that stands on it.
HKL_HEADROOM = 0.20

# The crystal systems a --cell of each count can be, the one it is taken as
# without --system first; None where --system has to say.
CELL_COUNTS = {
    1: ("cubic",),
    2: ("tetragonal", "hexagonal", "trigonal"),
    3: ("orthorhombic",),
    4: (None, "monoclinic"),
    6: ("triclinic",),
}


class CommandError(Exception):
    """A command could not go on; the message is printed as one line."""


def _jsonable(value):
    """Return ``value`` with numpy types as Python ones, non-finite floats as None."""
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _read_scans(paths: Sequence[str]) -> list[XRDScan]:
    """Read every scan in ``paths``, having first checked that each one exists."""
    for path in paths:
        if not Path(path).is_file():
            raise CommandError(f"no such file: {path}")
    scans = []
    for path in paths:
        try:
            scans.append(read_scan(path))
        except (ET.ParseError, ValueError) as error:
            # A reader that already names the file says it once.
            if str(error).startswith(str(Path(path))):
                raise CommandError(str(error)) from error
            raise CommandError(f"cannot read {path}: {error}") from error
    return scans


@dataclass(frozen=True)
class _Input:
    """A pattern named on the command line: a file, or a sample of a project."""

    path: Path
    project: Project | None = None
    sample: Sample | None = None

    @property
    def name(self) -> str:
        """The sample key, or the file's stem."""
        return self.sample.key if self.sample else self.path.stem

    @property
    def structure(self) -> StructureSpec | None:
        """The sample's first structure."""
        if self.sample is None:
            return None
        return self.project.structures[self.sample.structures[0]]


def _project(cache: dict, argument: str) -> Project:
    """The project of the current folder, read once per command."""
    if "project" not in cache:
        try:
            path = find_project()
        except FileNotFoundError:
            raise CommandError(
                f"no such file: {argument}, and no {PROJECT_FILE} in the current "
                "folder or above it to look it up in as a sample key"
            ) from None
        try:
            cache["project"] = load_project(path)
        except ValueError as error:
            raise CommandError(str(error)) from None
    return cache["project"]


def _sample(cache: dict, key: str) -> _Input:
    """The sample ``key`` of the project of the current folder."""
    project = _project(cache, key)
    if key not in project.samples:
        raise CommandError(
            f"no such file or sample: {key}; the samples in "
            f"{project.root / PROJECT_FILE} are "
            + (", ".join(project.samples) or "none")
        )
    sample = project.samples[key]
    return _Input(path=sample.file, project=project, sample=sample)


def _resolve(argument: str, cache: dict) -> _Input:
    """The pattern ``argument`` names: an existing file, as it is; otherwise,
    when it has no folder separator in it, the sample of that key in the
    project of the current folder; otherwise a missing file."""
    if Path(argument).is_file():
        return _Input(path=Path(argument))
    if "/" in argument or "\\" in argument:
        raise CommandError(f"no such file: {argument}")
    return _sample(cache, argument)


def _read_inputs(
    inputs: Sequence[_Input], wavelength: float | None = None
) -> list[XRDScan]:
    """The scans of ``inputs``: each with ``wavelength`` when it is given,
    and otherwise a sample's with the K alpha 1 wavelength of its instrument
    and a file's with its own."""
    scans = _read_scans([str(item.path) for item in inputs])
    for index, item in enumerate(inputs):
        if wavelength is not None:
            scans[index] = replace(scans[index], wavelength=wavelength)
        elif item.sample is not None:
            instrument = item.project.instruments[item.sample.instrument]
            scans[index] = replace(scans[index], wavelength=instrument.wavelength[0])
    return scans


def _require_wavelength(scan: XRDScan, item: _Input) -> float:
    """The scan's wavelength, refusing a pattern that carries none and was
    given none: a two column file says nothing about the radiation."""
    if scan.wavelength is None:
        raise CommandError(
            f"{item.path} carries no wavelength; give --wavelength ANGSTROM "
            "or a sample key"
        )
    return scan.wavelength


def _output_folders(
    args: argparse.Namespace, project: Project | None, name: str
) -> tuple[Path, Path]:
    """The folders results and figures go in: ``results/`` and ``figures/``
    under ``--out``, the current folder by default; or, for a sample of a
    project and no ``--out``, both ``results/<command>/<name>`` under the
    project root."""
    if args.out is None and project is not None:
        folder = results_dir(project, args.command, name)
        return folder, folder
    out = Path(args.out or ".")
    return out / "results", out / "figures"


def _finish(args: argparse.Namespace, written: list[Path], **values) -> int:
    """Print the files written, one per line, or them and ``values`` as JSON."""
    if args.json:
        print(json.dumps(_jsonable({"files": written, **values}), indent=2))
    else:
        for path in written:
            print(path)
    return 0


def _run_check(args: argparse.Namespace) -> int:
    item = _resolve(args.scan, {})
    (scan,) = _read_inputs([item])
    quality = assess_scan(scan)
    if item.sample is None:
        if args.json:
            print(json.dumps(_jsonable(asdict(quality)), indent=2))
        else:
            print(format_report(quality))
        return 0

    # A sample's report is headed by its key and composition, and kept.
    folder = results_dir(item.project, "check", item.name)
    folder.mkdir(parents=True, exist_ok=True)
    composition = item.structure.composition
    if args.json:
        path = folder / f"check_{item.name}.json"
        result = {
            "sample": item.name,
            "composition": composition,
            **_jsonable(asdict(quality)),
        }
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({**result, "files": [str(path)]}, indent=2))
    else:
        path = folder / f"check_{item.name}.txt"
        report = f"sample  {item.name}, {composition}\n{format_report(quality)}"
        path.write_text(report + "\n", encoding="utf-8")
        print(report)
        print(path)
    return 0


def _add_check(subparsers) -> None:
    parser = subparsers.add_parser(
        "check",
        help="report the data quality of a scan, with a verdict per workflow",
        description=(
            "Read a scan, print its range, step, counting time and "
            "intensities, and say whether it is good enough for plotting, phase "
            "identification, a Le Bail fit and a Rietveld refinement, by the "
            "criteria of the user guide, Section 2. For a sample of the project "
            "file the report is headed by its key and composition and is also "
            "written to results/check/KEY under the project root."
        ),
    )
    parser.add_argument(
        "scan",
        metavar="SCAN",
        help="path to a .xrdml, .xy or .xye file, or a sample key",
    )
    parser.add_argument(
        "--json", action="store_true", help="print the result as JSON instead"
    )
    parser.set_defaults(handler=_run_check)


def _add_output_options(parser: argparse.ArgumentParser, stem_default: str) -> None:
    """Add the options every command that writes files shares."""
    parser.add_argument(
        "--stem", help=f"name the output files end in (default: {stem_default})"
    )
    parser.add_argument(
        "--out",
        metavar="DIR",
        help=(
            "output root; results/ and figures/ are created under it as needed "
            "(default: the current folder, or results/COMMAND/KEY under the "
            "project root for a sample of the project file)"
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="print the files written, and any results, as JSON",
    )


def _add_cell_options(
    parser: argparse.ArgumentParser, purpose: str, esds: bool = False
) -> None:
    """Add --cell and --system, and with ``esds`` --esd-cell."""
    parser.add_argument(
        "--cell",
        nargs="+",
        type=float,
        metavar="X",
        help=(
            f"{purpose}, as the free parameters of its crystal system in "
            "angstroms and degrees: a (cubic); a c (tetragonal, or hexagonal or "
            "trigonal with --system); a b c (orthorhombic); a b c beta "
            "(monoclinic, with --system monoclinic); or a b c alpha beta gamma "
            "(triclinic)"
        ),
    )
    parser.add_argument(
        "--system",
        choices=CRYSTAL_SYSTEMS,
        help="crystal system of --cell, where its count of numbers leaves a choice",
    )
    if esds:
        parser.add_argument(
            "--esd-cell",
            nargs="+",
            type=float,
            metavar="E",
            help="esds of the --cell numbers, as many and in the same order",
        )


def _cell_option(args: argparse.Namespace) -> Cell | None:
    """The cell ``--cell`` and ``--system`` give, or ``None`` without --cell.

    Raises
    ------
    CommandError
        If --system is given without --cell, the count of numbers fits no
        crystal system or not the one --system names, --system is needed and
        missing, or a value is not greater than 0 or makes no cell.
    """
    if args.cell is None:
        if args.system is not None:
            raise CommandError("--system needs --cell")
        return None
    count = len(args.cell)
    if count not in CELL_COUNTS:
        raise CommandError(f"--cell takes 1, 2, 3, 4 or 6 numbers, not {count}")
    system = args.system or CELL_COUNTS[count][0]
    if system is None:
        raise CommandError(
            f"--cell with {count} numbers needs --system monoclinic (a b c beta)"
        )
    if system not in CELL_COUNTS[count]:
        names = CELL_PARAMETERS[system]
        raise CommandError(
            f"--system {system} takes {len(names)} "
            f"number{'s' if len(names) > 1 else ''} in --cell "
            f"({' '.join(names)}), not {count}"
        )
    for value in args.cell:
        if not value > 0:
            raise CommandError(f"--cell values must be greater than 0, not {value:g}")
    try:
        return Cell.from_parameters(
            system, dict(zip(CELL_PARAMETERS[system], args.cell))
        )
    except ValueError as error:
        raise CommandError(f"--cell: {error}") from None


def _cell_esds(args: argparse.Namespace, cell: Cell) -> dict[str, float] | None:
    """The esds ``--esd-cell`` gives each free parameter of ``cell``, or
    ``None`` without it."""
    if args.esd_cell is None:
        return None
    if len(args.esd_cell) != len(args.cell):
        raise CommandError(
            f"--esd-cell takes as many numbers as --cell, {len(args.cell)}, "
            f"not {len(args.esd_cell)}"
        )
    return dict(zip(cell.parameter_names, args.esd_cell))


def _cell_text(parameters: dict[str, float], esds: dict | None = None) -> str:
    """The free parameters of a cell on one line, lengths then angles, each
    with its esd where there is one."""
    esds = esds or {}
    lengths = [name for name in parameters if name in ("a", "b", "c")]
    angles = [name for name in parameters if name not in lengths]
    parts = [
        ", ".join(
            f"{name} = {_plus_minus(parameters[name], esds.get(name), digits)}"
            for name in names
        )
        + f" {unit}"
        for names, digits, unit in ((lengths, 4, "angstrom"), (angles, 3, "degrees"))
        if names
    ]
    return ", ".join(parts)


def _add_figure_options(parser: argparse.ArgumentParser, stem_default: str) -> None:
    """Add the options every figure command shares."""
    _add_output_options(parser, stem_default)
    parser.add_argument(
        "--scale",
        choices=SCALES,
        default="linear",
        help="intensity scale (default: linear)",
    )


def _hkl_figure(
    scan: XRDScan, indexed: list[IndexedPeak], scale: str, label: str
) -> tuple[Figure, Axes]:
    """Plot ``scan`` with hkl labels standing on its peaks.

    The top of the y axis is raised to the highest point of the trace plus
    ``HKL_HEADROOM`` of the axis span, on the chosen scale, before the labels
    are placed, so that the label on the tallest peak stays inside the axes.
    """
    fig, ax = plot_pattern(scan, scale=scale, label=label)
    bottom, top = ax.get_ylim()
    highest = float(np.max(ax.lines[0].get_ydata()))
    ax.set_ylim(bottom, highest + HKL_HEADROOM * (top - bottom))
    ax.legend(frameon=False)
    annotate_hkl(ax, indexed, line=ax.lines[0])
    return fig, ax


def _peaks(item: _Input, scan: XRDScan) -> list:
    """The peaks of ``scan`` with the K alpha 2 satellites flagged; for a
    sample, by the ratio of its instrument's two wavelengths, or not at all
    for an instrument without K alpha 2."""
    if item.sample is None:
        return find_peaks(scan)
    instrument = item.project.instruments[item.sample.instrument]
    peaks = find_peaks(scan, flag_satellites=False)
    if instrument.ka2:
        ratio = instrument.wavelength[1] / instrument.wavelength[0]
        flag_kalpha2(peaks, wavelength_ratio=ratio)
    return peaks


def _crystal_system_of(cell: dict[str, float]) -> str:
    """The highest crystal system six cell parameters fit, for a structure
    whose own crystal system is not known."""
    a, b, c = cell["a"], cell["b"], cell["c"]
    angles = (cell["alpha"], cell["beta"], cell["gamma"])
    if angles == (90.0, 90.0, 90.0):
        if a == b == c:
            return "cubic"
        return "tetragonal" if a == b else "orthorhombic"
    if angles == (90.0, 90.0, 120.0) and a == b:
        return "hexagonal"
    return "monoclinic" if angles[0] == angles[2] == 90.0 else "triclinic"


def _structure_cell(item: _Input) -> tuple[Cell | None, str | None, str | None]:
    """The start cell and space group to index a sample from, taken from its
    first structure, or no cell and the reason.

    A library structure gives its entry's crystal system and space group. A
    CIF structure's crystal system is the highest its six parameters fit, and
    its space group is not known here.
    """
    spec = item.structure
    try:
        parameters = resolved_cell(spec)
        if spec.library is not None:
            entry = load_entry(spec.library)
            system, space_group = entry.crystal_system, entry.space_group
        else:
            system, space_group = _crystal_system_of(parameters), None
        cell = Cell.from_parameters(system, parameters)
    except ValueError as error:
        return None, None, str(error)
    return cell, space_group, None


def _run_plot(args: argparse.Namespace) -> int:
    cell = _cell_option(args)
    item = _resolve(args.scan, {})
    (scan,) = _read_inputs([item], args.wavelength)
    stem = args.stem or item.name
    results, figures = _output_folders(args, item.project, item.name)
    written: list[Path] = []
    values: dict = {}

    apply_style()
    peaks = _peaks(item, scan)
    if args.no_satellites:
        peaks = exclude_kalpha2(peaks)
    written.append(peaks_to_csv(peaks, results / f"peaks_{stem}.csv"))

    if args.label is not None:
        label = args.label
    elif item.sample is not None:
        label = f"{item.name}, {item.structure.composition}"
    else:
        label = scan.sample_id
    fig, ax = plot_pattern(scan, scale=args.scale, label=label)
    ax.legend(frameon=False)
    written += save_figure(fig, figures / f"pattern_{stem}")

    # A cell given on the command line asks for labels, so failing to index
    # from it is an error; one taken from a sample's structure is only tried.
    space_group, automatic = args.space_group, False
    if cell is None and item.sample is not None:
        cell, structure_group, reason = _structure_cell(item)
        if reason is not None:
            print(
                f"xrdkit plot: {reason}, so the peaks are not labelled with hkl",
                file=sys.stderr,
            )
        space_group, automatic = space_group or structure_group, True
        if cell is not None and space_group is None:
            print(
                "xrdkit plot: the space group of a CIF structure is not known "
                "here, so the peaks are indexed without reflection conditions and "
                "some labels may name absent reflections",
                file=sys.stderr,
            )
    if cell is not None:
        if space_group is not None and space_group not in SUPPORTED_SPACE_GROUPS:
            print(
                f"xrdkit plot: the reflection conditions of {space_group} are not "
                "known, so the peaks are indexed without them and some labels may "
                "name absent reflections",
                file=sys.stderr,
            )
            space_group = None
        try:
            indexed, fit = index_and_refine(
                peaks,
                start_cell=cell,
                wavelength=_require_wavelength(scan, item),
                zero_offset=args.zero,
                space_group=space_group,
            )
        except ValueError as error:
            if not automatic:
                raise CommandError(f"indexing failed: {error}") from error
            print(
                f"xrdkit plot: indexing failed: {error}; plotted without hkl labels",
                file=sys.stderr,
            )
            return _finish(args, written)
        written.append(indexed_to_csv(indexed, results / f"indexed_{stem}.csv"))

        fig, _ = _hkl_figure(scan, indexed, args.scale, label)
        written += save_figure(fig, figures / f"pattern_hkl_{stem}")

        summary = indexing_summary(indexed)
        values = {
            "crystal_system": fit.cell.crystal_system,
            **fit.cell.parameters,
            "zero": fit.zero_offset,
            "n_indexed": summary["n_indexed"],
        }
        if not args.json:
            print(
                f"{_cell_text(fit.cell.parameters)}, "
                f"zero {fit.zero_offset:.3f} degrees, "
                f"{summary['n_indexed']} of {summary['n_peaks']} peaks indexed, "
                f"rms {summary['rms_difference']:.4f} degrees"
            )

    return _finish(args, written, **values)


def _add_plot(subparsers) -> None:
    parser = subparsers.add_parser(
        "plot",
        help="plot one scan, list its peaks, and label them with hkl given a cell",
        description=(
            "Plot one scan and write its peak list. With --cell the peaks "
            "are indexed from that start cell, of any crystal system, the "
            "indexing is written out, and a second figure is labelled with hkl; "
            "if they cannot be indexed the command fails. Without --cell no "
            "indexing is done. A sample of the project file is read at its "
            "instrument's wavelength and labelled with its key and composition, "
            "and without --cell is indexed from the cell, crystal system and "
            "space group of its first structure; if that fails the pattern is "
            "plotted without hkl labels, with a note. A space group whose "
            "reflection conditions are not known is indexed without them, with a "
            "note."
        ),
    )
    parser.add_argument(
        "scan",
        metavar="SCAN",
        help="path to a .xrdml, .xy or .xye file, or a sample key",
    )
    parser.add_argument(
        "--label",
        metavar="TEXT",
        help=(
            "legend label for the trace (default: the scan's sample id, or a "
            "sample's key and composition)"
        ),
    )
    parser.add_argument(
        "--wavelength",
        type=float,
        metavar="ANGSTROM",
        help=(
            "K alpha 1 wavelength for the d spacings and the indexing (default: "
            "the scan's own, or a sample's instrument's)"
        ),
    )
    _add_cell_options(parser, "start cell to index from; hkl labels need one")
    parser.add_argument(
        "--space-group",
        metavar="SG",
        help=(
            "space group whose reflection conditions apply (default: none, or a "
            "sample's structure's)"
        ),
    )
    parser.add_argument(
        "--zero",
        type=float,
        default=DEFAULT_ZERO_OFFSET,
        metavar="DEG",
        help=(
            f"zero offset in degrees (default: {DEFAULT_ZERO_OFFSET:g}, which "
            "searches for one)"
        ),
    )
    parser.add_argument(
        "--no-satellites",
        action="store_true",
        help=(
            "drop the peaks flagged as K alpha 2 satellites from the peak list "
            "and the indexing"
        ),
    )
    _add_figure_options(parser, "the scan's file stem")
    parser.set_defaults(handler=_run_plot)


def _run_stack(args: argparse.Namespace) -> int:
    cache: dict = {}
    inputs = [_resolve(argument, cache) for argument in args.scans]
    scans = _read_inputs(inputs)
    if len(args.labels) != len(scans):
        raise CommandError(
            f"{len(args.labels)} labels for {len(scans)} scans; give one label per scan"
        )
    stem = args.stem or "_".join(item.name for item in inputs)
    project = next((item.project for item in inputs if item.sample), None)
    _, figures = _output_folders(args, project, stem)

    apply_style()
    fig, _, _, _ = plot_stacked(
        scans,
        labels=args.labels,
        offset=args.offset,
        scale=args.scale,
        normalise=args.normalise,
    )
    written = save_figure(fig, figures / f"stack_{stem}")
    return _finish(args, written)


def _add_stack(subparsers) -> None:
    parser = subparsers.add_parser(
        "stack",
        help="plot several scans stacked one above the other",
        description=(
            "Plot several scans stacked vertically, drawn bottom to top in "
            "the order given, each with its label. Files and sample keys of the "
            "project file may be mixed; with a sample among them and no --out "
            "the figure goes to results/stack/STEM under the project root."
        ),
    )
    parser.add_argument(
        "scans",
        metavar="SCAN",
        nargs="+",
        help="paths to .xrdml, .xy or .xye files, or sample keys",
    )
    parser.add_argument(
        "--labels",
        metavar="TEXT",
        nargs="+",
        required=True,
        help="one label per scan, in the same order",
    )
    parser.add_argument(
        "--offset",
        type=float,
        metavar="F",
        help="vertical spacing between traces (default: 1.2 times the tallest)",
    )
    parser.add_argument(
        "--no-normalise",
        dest="normalise",
        action="store_false",
        help="keep the measured intensities rather than scaling each trace to 1",
    )
    _add_figure_options(parser, "the scan stems joined by an underscore")
    parser.set_defaults(handler=_run_stack)


DENSITY_METHOD = (
    "theoretical density from cell volume and formula mass; "
    "relative density = measured / theoretical"
)


def _plus_minus(value: float, esd: float | None, digits: int) -> str:
    """Return ``value``, and ``+/- esd`` when there is one, to ``digits`` places."""
    if esd is None:
        return f"{value:.{digits}f}"
    return f"{value:.{digits}f} +/- {esd:.{digits}f}"


def _density_from_project(args: argparse.Namespace) -> tuple[_Input, Cell | None]:
    """Fill the options not given from the sample ``args.sample`` and its
    first structure: the formula, z and the Archimedes density; and return the
    structure's cell when neither --cell nor --volume is given, else None."""
    item = _sample({}, args.sample)
    spec = item.structure
    try:
        if args.formula is None:
            args.formula = spec.composition
        if args.z is None:
            args.z = resolved_z(spec)
        cell = None
        if args.cell is None and args.volume is None:
            parameters = resolved_cell(spec)
            if spec.library is not None:
                system = load_entry(spec.library).crystal_system
            else:
                system = _crystal_system_of(parameters)
            cell = Cell.from_parameters(system, parameters)
    except ValueError as error:
        raise CommandError(f"{error}; or give it on the command line") from None
    if args.archimedes is None and item.sample.archimedes is not None:
        args.archimedes = [item.sample.archimedes]
    return item, cell


def _cell_columns(cell: Cell | None, esds: dict | None) -> dict:
    """The six cell parameters of ``cell`` and their esds as density columns,
    each ``None`` without a cell, and an esd ``None`` where none was given."""
    esds = esds or {}
    columns = {}
    for name in CELL_PARAMETERS["triclinic"]:
        unit = "angstrom" if name in ("a", "b", "c") else "deg"
        columns[f"{name}_{unit}"] = getattr(cell, name) if cell else None
        columns[f"esd_{name}_{unit}"] = esds.get(name)
    return columns


def _run_density(args: argparse.Namespace) -> int:
    item, cell = (None, None) if args.sample is None else _density_from_project(args)

    # Checked here rather than by argparse, so that a missing or conflicting
    # option returns 1 with one line, as every other failure of a command does.
    if args.formula is None:
        raise CommandError("--formula is required")
    if args.z is None:
        raise CommandError("--z is required")
    if cell is None and (args.cell is None) == (args.volume is None):
        raise CommandError("give one of --cell and --volume, not both or neither")
    if args.esd_cell is not None and args.cell is None:
        raise CommandError("--esd-cell needs --cell")
    if args.esd_volume is not None and args.volume is None:
        raise CommandError("--esd-volume needs --volume")
    if args.archimedes is not None and len(args.archimedes) > 2:
        raise CommandError("--archimedes takes a density and at most one esd")

    esds = None
    if args.cell is not None:
        cell = _cell_option(args)
        esds = _cell_esds(args, cell)
    elif args.system is not None:
        raise CommandError("--system needs --cell")

    measured = esd_measured = relative = esd_relative = None
    try:
        mass = formula_mass(args.formula)
        if cell is not None:
            volume, esd_volume = cell_volume(cell, esds)
        else:
            volume, esd_volume = args.volume, args.esd_volume
        density, esd_density = theoretical_density(
            args.formula, args.z, volume, esd_volume
        )
        if args.archimedes is not None:
            measured = args.archimedes[0]
            esd_measured = args.archimedes[1] if len(args.archimedes) == 2 else None
            relative, esd_relative = relative_density(
                measured, density, esd_measured, esd_density
            )
    except ValueError as error:
        raise CommandError(str(error)) from error

    row = {
        **({"sample": item.name, "structure": item.structure.key} if item else {}),
        "formula": args.formula,
        "z": args.z,
        "crystal_system": cell.crystal_system if cell else None,
        **_cell_columns(cell, esds),
        "volume_a3": volume,
        "esd_volume_a3": esd_volume,
        "formula_mass_g_mol": mass,
        "theoretical_density_g_cm3": density,
        "esd_theoretical_density_g_cm3": esd_density,
        "archimedes_density_g_cm3": measured,
        "esd_archimedes_density_g_cm3": esd_measured,
        "relative_density_percent": relative,
        "esd_relative_density_percent": esd_relative,
        "method": DENSITY_METHOD,
        "date": datetime.datetime.now(datetime.UTC).astimezone().date().isoformat(),
        "xrdkit_version": __version__,
    }
    if args.stem:
        name = f"density_{args.stem}.csv"
    else:
        name = f"density_{item.name}.csv" if item else "density.csv"
    if item is None:
        results, _ = _output_folders(args, None, "")
    else:
        results, _ = _output_folders(args, item.project, item.name)
    path = results / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(row)
        writer.writerow(["" if value is None else value for value in row.values()])

    if not args.json:
        print(f"M = {mass:.3f} g/mol per formula unit")
        if cell is not None:
            print(f"{cell.crystal_system} cell {_cell_text(cell.parameters, esds)}")
        print(f"V = {_plus_minus(volume, esd_volume, 3)} cubic angstrom")
        print(f"theoretical density {_plus_minus(density, esd_density, 4)} g/cm3")
        if measured is not None:
            print(f"Archimedes density {_plus_minus(measured, esd_measured, 4)} g/cm3")
            print(f"relative density {_plus_minus(relative, esd_relative, 2)} per cent")
    return _finish(args, [path], **row)


def _add_density(subparsers) -> None:
    parser = subparsers.add_parser(
        "density",
        help="theoretical density from a formula and a cell, and relative density",
        description=(
            "Work out the formula mass, the cell volume and the theoretical "
            "density, Z M / (N_A V), each with its esd where one can be given, "
            "and with --archimedes the relative density. Give the cell of any "
            "crystal system with --cell, or its volume with --volume. One row "
            "with every input and result, the crystal system and cell used, the "
            "method and "
            "the date is written to results/density.csv, or to "
            "results/density_STEM.csv with --stem. Given a sample key of the "
            "project file, the formula, z, cell and Archimedes density not given "
            "as options come from the sample and its first structure, and the "
            "row goes to results/density/KEY/density_KEY.csv under the project "
            "root."
        ),
    )
    parser.add_argument(
        "sample",
        metavar="SAMPLE",
        nargs="?",
        help="a sample key of the project file to take the values from",
    )
    parser.add_argument(
        "--formula",
        metavar="TEXT",
        help='formula of one formula unit, such as "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"',
    )
    parser.add_argument(
        "--z", type=float, metavar="N", help="formula units per cell (required)"
    )
    _add_cell_options(parser, "cell", esds=True)
    parser.add_argument(
        "--volume",
        type=float,
        metavar="V",
        help="cell volume in cubic angstroms, for a cell of any symmetry",
    )
    parser.add_argument(
        "--esd-volume",
        type=float,
        metavar="EV",
        help="esd of the volume, in cubic angstroms",
    )
    parser.add_argument(
        "--archimedes",
        nargs="+",
        type=float,
        metavar=("RHO", "ESD"),
        help="measured density in g/cm3, and optionally its esd",
    )
    _add_output_options(parser, "none, which writes results/density.csv")
    parser.set_defaults(handler=_run_density)


# A refitted peak position is kept only if the fit converged and moved the
# position by no more than this many of the found peak's FWHM.
REFIT_MAX_SHIFT_FWHM = 1.0

LATTICE_CAUTION = (
    "a low indexed fraction or a high rms means the cell should not be trusted"
)


def _refit_peaks(
    scan: XRDScan, peaks: list[Peak], ka2: bool, wavelength_ratio: float
) -> tuple[list[Peak], list[dict]]:
    """The peaks with their positions refitted by :func:`fit_profile`, and a
    record of each fit.

    Each peak is fitted as a K alpha doublet, or a single line without K alpha
    2, from its found position and width, so that the position is that of K
    alpha 1 rather than the vertex of the unresolved doublet. A fit that fails,
    does not converge or moves the position by more than the peak's FWHM is
    rejected, and the found position kept. A peak flagged as a K alpha 2
    satellite is part of its parent's doublet and is not fitted.
    """
    refitted, records = [], []
    for peak in peaks:
        record = {"fitted": None, "esd": None, "rejected": None, "recovered": False}
        if peak.kalpha2_of is None:
            try:
                fit = fit_profile(
                    scan.two_theta,
                    scan.intensity,
                    peak.two_theta,
                    peak.fwhm,
                    wavelength_ratio=wavelength_ratio,
                    intensity_ratio=KALPHA2_INTENSITY_RATIO if ka2 else 0.0,
                )
            except ValueError:
                fit = None
            accepted = (
                fit is not None
                and fit.converged
                and math.isfinite(fit.two_theta)
                and abs(fit.two_theta - peak.two_theta)
                <= REFIT_MAX_SHIFT_FWHM * peak.fwhm
            )
            record["rejected"] = not accepted
            if accepted:
                record["fitted"], record["esd"] = fit.two_theta, fit.esd_two_theta
                theta = math.radians(fit.two_theta / 2.0)
                peak = replace(
                    peak,
                    two_theta=fit.two_theta,
                    d_spacing=scan.wavelength / (2.0 * math.sin(theta)),
                )
        refitted.append(peak)
        records.append(record)
    return refitted, records


def _index_and_fit(
    peaks: list[Peak],
    cell: Cell,
    wavelength: float,
    zero: float | None,
    space_group: str | None,
    zero_free: bool,
    displacement_free: bool,
    radius: float | None,
) -> tuple[list[IndexedPeak | None], CellFit, LatticeFit]:
    """Index and refine on the peaks not flagged as K alpha 2 satellites.

    Returns one entry per peak, None for a flagged peak, which takes no part,
    with the fit of the indexing and the refined lattice.

    Raises
    ------
    CommandError
        If the indexing or the refinement fails.
    """
    taking_part = [index for index, peak in enumerate(peaks) if peak.kalpha2_of is None]
    try:
        indexed, cell_fit = index_and_refine(
            [peaks[index] for index in taking_part],
            start_cell=cell,
            wavelength=wavelength,
            zero_offset=zero or 0.0,
            fine_tolerance=DEFAULT_TOLERANCE,
            space_group=space_group,
        )
    except ValueError as error:
        raise CommandError(f"indexing failed: {error}") from error
    try:
        fit = refine_lattice(
            indexed,
            wavelength,
            cell_fit.cell,
            fit_zero=zero_free,
            fit_displacement=displacement_free,
            radius_mm=radius,
            start_zero=cell_fit.zero_offset if zero_free else (zero or 0.0),
        )
    except ValueError as error:
        raise CommandError(f"refinement failed: {error}") from error
    entries: list[IndexedPeak | None] = [None] * len(peaks)
    for index, entry in zip(taking_part, indexed):
        entries[index] = entry
    return entries, cell_fit, fit


def _recover_satellites(
    peaks: list[Peak], cell_fit: CellFit, wavelength: float, space_group: str | None
) -> list[int]:
    """The indices of the flagged peaks that a reflection may need: those
    within the indexing's own tolerance of a reflection of ``cell_fit``'s cell,
    at the zero offset the indexing used. They are candidates for recovery,
    which :func:`_refit_recovered` settles."""
    flagged = [index for index, peak in enumerate(peaks) if peak.kalpha2_of is not None]
    if not flagged:
        return []
    matches = index_peaks(
        [peaks[index] for index in flagged],
        cell_fit.cell,
        wavelength,
        DEFAULT_TOLERANCE,
        cell_fit.zero_offset,
        space_group,
    )
    return [index for index, match in zip(flagged, matches) if match.is_indexed]


def _refit_recovered(
    scan: XRDScan,
    found: list[Peak],
    peaks: list[Peak],
    records: list[dict],
    candidates: list[int],
    ka2: bool,
    wavelength_ratio: float,
) -> list[int]:
    """Refit the flagged peaks at ``candidates`` as K alpha 1 lines, in place,
    and return the indices of those recovered.

    The doublet refit is the test of the claim that such a peak is the K alpha
    1 line of a reflection, so a candidate whose refit is rejected is not
    recovered: it keeps its flag and takes no part in the refinement. An
    ordinary peak whose refit is rejected is still used at its found position.

    An ordinary peak is already known to be a reflection, so a failed fit costs only precision, whereas a recovered peak has no evidence for being a reflection other than the fit itself, and admitting it at its raw position admits a satellite at a satellite's position.
    """
    recovered = []
    for index in candidates:
        cleared = replace(found[index], kalpha2_of=None)
        (peak,), (record,) = _refit_peaks(scan, [cleared], ka2, wavelength_ratio)
        if record["rejected"]:
            continue
        found[index], peaks[index] = cleared, peak
        records[index] = {**record, "recovered": True}
        recovered.append(index)
    return recovered


def _peak_counts(
    found: list[Peak],
    records: list[dict],
    indexed: list[IndexedPeak | None],
    fit: LatticeFit,
) -> dict[str, int]:
    """The peak counts of a lattice run, for the results row and the report.

    found = refitted + rejected + satellites; used <= indexed <= refitted + rejected; recovered <= refitted.

    A recovered peak is no longer a satellite and always counts as refitted,
    since a candidate whose refit is rejected is not recovered. The indexed
    fraction is indexed over refitted plus rejected, the peaks that took part.
    """
    return {
        "n_peaks_found": len(found),
        "n_peaks_refitted": sum(record["rejected"] is False for record in records),
        "n_fits_rejected": sum(record["rejected"] is True for record in records),
        "n_satellites": sum(peak.kalpha2_of is not None for peak in found),
        "n_peaks_recovered": sum(record["recovered"] for record in records),
        "n_peaks_indexed": sum(
            entry is not None and entry.is_indexed for entry in indexed
        ),
        "n_peaks_refined": fit.n_peaks,
    }


def _calculated_two_theta(
    fit: LatticeFit, hkl: tuple[int, int, int], wavelength: float
) -> float:
    """Where the refined cell, zero and displacement put reflection ``hkl``."""
    sin_theta = wavelength / (2.0 * fit.cell.d_spacing(*hkl))
    if sin_theta > 1.0:
        return float("nan")
    theta = math.asin(sin_theta)
    shift = fit.zero
    if fit.displacement is not None:
        shift += math.degrees(-2.0 * fit.displacement * math.cos(theta) / fit.radius_mm)
    return 2.0 * math.degrees(theta) + shift


def _corrected_two_theta(fit: LatticeFit, observed: float) -> float:
    """``observed`` with the refined zero and displacement taken off: the
    position whose Bragg angle, moved as :func:`_calculated_two_theta` moves a
    reflection, lands on ``observed``. The displacement shift depends on that
    angle, so it is found by iteration, which the shift's small size makes
    converge at once."""
    corrected = observed - fit.zero
    if fit.displacement is None:
        return corrected
    for _ in range(5):
        theta = math.radians(corrected / 2.0)
        shift = math.degrees(-2.0 * fit.displacement * math.cos(theta) / fit.radius_mm)
        corrected = observed - fit.zero - shift
    return corrected


def _number(value: float | None, digits: int) -> str:
    """``value`` to ``digits`` places, or an empty cell for None or nan."""
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def _write_lattice_peaks(
    path: Path,
    found: list[Peak],
    records: list[dict],
    indexed: list[IndexedPeak | None],
    fit: LatticeFit,
    wavelength: float,
) -> Path:
    """Write one row per peak: found and fitted positions, the fit, the
    position after the refined zero and displacement are taken off and its d
    spacing, whether it is a satellite or was recovered from one, its
    intensity absolute and relative to the strongest peak, and the reflection
    assigned with its difference from the refined position and the number of
    reflections within the indexing tolerance.

    The corrected position is of the position that took part in the indexing,
    fitted or found, or of the found position for a peak that took no part;
    its d spacing is from ``wavelength`` and the position as written. The
    candidate count is left empty for a peak that took no part."""
    strongest = max(peak.intensity for peak in found)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(LATTICE_PEAK_COLUMNS)
        for peak, record, entry in zip(found, records, indexed):
            reflection = entry.reflection if entry is not None else None
            calculated = difference = None
            if reflection is not None:
                calculated = _calculated_two_theta(fit, reflection.hkl, wavelength)
                difference = entry.peak.two_theta - calculated
            position = entry.peak.two_theta if entry is not None else peak.two_theta
            corrected = _number(_corrected_two_theta(fit, position), 4)
            d_spacing = None
            if corrected:
                d_spacing = wavelength / (
                    2.0 * math.sin(math.radians(float(corrected) / 2.0))
                )
            writer.writerow(
                [
                    _number(peak.two_theta, 4),
                    _number(record["fitted"], 4),
                    _number(record["esd"], 5),
                    corrected,
                    _number(d_spacing, 5),
                    "" if record["rejected"] is None else record["rejected"],
                    peak.kalpha2_of is not None,
                    record["recovered"],
                    _number(peak.intensity, 1),
                    _number(100.0 * peak.intensity / strongest, 1),
                    _number(peak.fwhm, 4),
                    *(reflection.hkl if reflection else ("", "", "")),
                    _number(calculated, 4),
                    _number(difference, 4),
                    "" if entry is None else len(entry.candidates),
                ]
            )
    return path


LATTICE_PEAK_COLUMNS = (
    "found_two_theta",
    "fitted_two_theta",
    "esd_fitted_two_theta",
    "corrected_two_theta",
    "d_spacing",
    "fit_rejected",
    "kalpha2_satellite",
    "recovered",
    "intensity",
    "relative_intensity",
    "fwhm",
    "h",
    "k",
    "l",
    "calculated_two_theta",
    "difference",
    "n_candidates",
)


def _append_row(path: Path, row: dict) -> Path:
    """Append ``row`` to the CSV at ``path``, writing the header first when the
    file is new.

    Raises
    ------
    CommandError
        If the file exists with other columns.
    """
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle), [])
        if header != list(row):
            raise CommandError(
                f"{path} has other columns than this version writes; move it aside"
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.is_file()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        if new:
            writer.writerow(row)
        writer.writerow(["" if value is None else value for value in row.values()])
    return path


def _lattice_method(
    zero_free: bool, displacement_free: bool, density: bool, recover: bool
) -> str:
    """How the lattice command got its numbers, for the results row."""
    zero = "refined" if zero_free else "held"
    displacement = "refined" if displacement_free else "held at zero"
    satellites = (
        "peaks flagged as K alpha 2 satellites excluded, except those within the "
        "indexing tolerance of a reflection of the refined cell whose doublet "
        "refit is accepted, which are then used and the indexing and refinement "
        "run once more with them"
        if recover
        else "peaks flagged as K alpha 2 satellites excluded (none recovered, "
        "--no-satellites)"
    )
    method = (
        "peak positions fitted as K alpha 1 of a split pseudo-Voigt doublet "
        f"(fit_profile); {satellites}; indexed by index_and_refine with an "
        "adaptive coarse window; cell refined by least squares on the peak "
        f"positions (refine_lattice), zero {zero}, specimen displacement "
        f"{displacement}; volume esd propagated through the covariance"
    )
    if density:
        method += (
            "; theoretical density from the refined volume and formula mass; "
            "relative density = measured / theoretical"
        )
    return method


def _run_lattice(args: argparse.Namespace) -> int:
    # Every option is checked before a scan is read or a file written.
    cell = _cell_option(args)
    if args.archimedes is not None and len(args.archimedes) > 2:
        raise CommandError("--archimedes takes a density and at most one esd")
    if args.radius is not None and not args.radius > 0:
        raise CommandError(f"--radius must be greater than 0, not {args.radius:g}")
    item = _resolve(args.scan, {})
    sample, spec = item.sample, item.structure

    instrument = item.project.instruments[sample.instrument] if sample else None
    form = sample.form if sample else "powder"
    zero_free = form == "powder" and args.zero is None
    displacement_free = form == "pellet" or args.displacement
    radius = (
        args.radius if args.radius is not None else getattr(instrument, "radius", None)
    )
    if displacement_free and radius is None:
        raise CommandError(
            "refining a specimen displacement needs the goniometer radius: give "
            "--radius, or radius for the instrument in the project file"
        )

    space_group = args.space_group
    if cell is None:
        if sample is None:
            raise CommandError("--cell is required for a scan that is not a sample")
        cell, structure_group, reason = _structure_cell(item)
        if cell is None:
            raise CommandError(f"{reason}; or give --cell")
        space_group = space_group or structure_group
    if space_group is not None and space_group not in SUPPORTED_SPACE_GROUPS:
        print(
            f"xrdkit lattice: the reflection conditions of {space_group} are not "
            "known, so the peaks are indexed without them",
            file=sys.stderr,
        )
    conditions = space_group if space_group in SUPPORTED_SPACE_GROUPS else None

    formula, z = args.formula, args.z
    if sample is None and (formula is None) != (z is None):
        raise CommandError("--formula and --z go together, for a theoretical density")
    if spec is not None:
        formula = formula or spec.composition
        if z is None:
            try:
                z = resolved_z(spec)
            except ValueError:
                z = None
    has_density = formula is not None and z is not None
    archimedes = args.archimedes
    if archimedes is None and sample is not None and sample.archimedes is not None:
        archimedes = [sample.archimedes]
    if args.archimedes is not None and not has_density:
        raise CommandError(
            "--archimedes needs a theoretical density: give --formula and --z"
        )
    mass = None
    if has_density:
        try:
            mass = formula_mass(formula)
        except ValueError as error:
            raise CommandError(str(error)) from error

    (scan,) = _read_inputs([item], args.wavelength)
    _require_wavelength(scan, item)
    found = _peaks(item, scan)
    if not found:
        raise CommandError(f"no peaks found in {item.path}")
    ka2 = instrument.ka2 if instrument else True
    ratio = (
        instrument.wavelength[1] / instrument.wavelength[0]
        if instrument and instrument.ka2
        else KALPHA2_RATIO
    )
    peaks, records = _refit_peaks(scan, found, ka2, ratio)

    def index_and_fit() -> tuple[list[IndexedPeak | None], CellFit, LatticeFit]:
        return _index_and_fit(
            peaks,
            cell,
            scan.wavelength,
            args.zero,
            conditions,
            zero_free,
            displacement_free,
            radius,
        )

    indexed, cell_fit, fit = index_and_fit()
    candidates = (
        []
        if args.no_satellites
        else _recover_satellites(peaks, cell_fit, scan.wavelength, conditions)
    )
    if _refit_recovered(scan, found, peaks, records, candidates, ka2, ratio):
        indexed, cell_fit, fit = index_and_fit()
    counts = _peak_counts(found, records, indexed, fit)
    taking_part = counts["n_peaks_refitted"] + counts["n_fits_rejected"]

    density = esd_density = relative = esd_relative = measured = esd_measured = None
    if has_density:
        density, esd_density = theoretical_density(
            formula, z, fit.volume, fit.esd_volume
        )
        if archimedes is not None:
            measured = archimedes[0]
            esd_measured = archimedes[1] if len(archimedes) == 2 else None
            relative, esd_relative = relative_density(
                measured, density, esd_measured, esd_density
            )

    esds = {name: getattr(fit, f"esd_{name}") for name in CELL_PARAMETERS["triclinic"]}
    row = {
        "sample": item.name if sample else None,
        "structure": spec.key if spec else None,
        "scan_file": str(item.path),
        "wavelength_angstrom": scan.wavelength,
        "form": form,
        "space_group": space_group,
        "crystal_system": fit.cell.crystal_system,
        **{
            f"start_{name}": value
            for name, value in _cell_columns(cell, None).items()
            if not name.startswith("esd_")
        },
        **_cell_columns(fit.cell, esds),
        "volume_a3": fit.volume,
        "esd_volume_a3": fit.esd_volume,
        "zero_deg": fit.zero,
        "esd_zero_deg": fit.esd_zero,
        "zero_refined": zero_free,
        "displacement_mm": fit.displacement,
        "esd_displacement_mm": fit.esd_displacement,
        "displacement_refined": bool(displacement_free),
        "radius_mm": radius,
        **counts,
        "rms_two_theta_deg": fit.rms_two_theta,
        "coarse_two_theta_max_deg": cell_fit.coarse_two_theta_max,
        "formula": formula if has_density else None,
        "z": z if has_density else None,
        "formula_mass_g_mol": mass,
        "theoretical_density_g_cm3": density,
        "esd_theoretical_density_g_cm3": esd_density,
        "archimedes_density_g_cm3": measured,
        "esd_archimedes_density_g_cm3": esd_measured,
        "relative_density_percent": relative,
        "esd_relative_density_percent": esd_relative,
        "method": _lattice_method(
            zero_free, displacement_free, has_density, not args.no_satellites
        ),
        "date": datetime.datetime.now(datetime.UTC).astimezone().date().isoformat(),
        "xrdkit_version": __version__,
    }

    if args.out is None and item.project is not None:
        folder = results_dir(item.project, "lattice", item.name)
    else:
        folder = Path(args.out or ".") / "results" / "lattice"
    results = folder / (f"lattice_{item.name}.csv" if sample else "lattice.csv")
    stem = args.stem or item.name
    peaks_path = _write_lattice_peaks(
        folder / f"peaks_{stem}.csv", found, records, indexed, fit, scan.wavelength
    )
    written = [peaks_path, _append_row(results, row)]

    if not args.json:
        if sample is not None:
            print(f"sample  {item.name}, {spec.composition}")
        indexed_percent = 100.0 * counts["n_peaks_indexed"] / taking_part
        print(
            f"peaks   {counts['n_peaks_found']} found, "
            f"{counts['n_peaks_refitted']} refitted, "
            f"{counts['n_fits_rejected']} fits rejected, "
            f"{counts['n_satellites']} satellites excluded, "
            f"{counts['n_peaks_recovered']} recovered; "
            f"{counts['n_peaks_indexed']} of {taking_part} indexed "
            f"({indexed_percent:.1f} per cent), "
            f"{counts['n_peaks_refined']} used in the refinement"
        )
        print(f"coarse window to {cell_fit.coarse_two_theta_max:.2f} degrees")
        print(f"{fit.cell.crystal_system} cell {_cell_text(fit.cell.parameters, esds)}")
        print(f"V = {_plus_minus(fit.volume, fit.esd_volume, 3)} cubic angstrom")
        held = "" if zero_free else " (held)"
        print(f"zero {_plus_minus(fit.zero, fit.esd_zero, 4)} degrees{held}")
        if fit.displacement is None:
            print("displacement not refined")
        else:
            print(
                f"displacement {_plus_minus(fit.displacement, fit.esd_displacement, 4)} "
                f"mm, radius {radius:g} mm"
            )
        print(f"rms {fit.rms_two_theta:.4f} degrees; {LATTICE_CAUTION}")
        if has_density:
            print(f"M = {mass:.3f} g/mol per formula unit, Z = {z:g}")
            print(f"theoretical density {_plus_minus(density, esd_density, 4)} g/cm3")
            if measured is not None:
                print(
                    f"Archimedes density {_plus_minus(measured, esd_measured, 4)} g/cm3"
                )
                print(
                    f"relative density {_plus_minus(relative, esd_relative, 2)} per cent"
                )
    return _finish(args, written, **row)


def _add_lattice(subparsers) -> None:
    parser = subparsers.add_parser(
        "lattice",
        help="refine the cell of a scan, with its volume and density",
        description=(
            "Find the peaks of a scan, refit each position as the K alpha 1 line "
            "of its doublet, index them against a start cell of any crystal "
            "system, and refine the cell with the zero or the specimen "
            "displacement by least squares. Peaks flagged as K alpha 2 "
            "satellites are left out, except those the refined cell puts within "
            "the indexing tolerance of a reflection and whose doublet refit is "
            "accepted, which are used and the indexing and refinement run once "
            "more with them. A pellet refines the displacement "
            "with the zero held (at --zero, the instrument zero from a standard, "
            "or 0); a powder, or a scan that is not a sample, refines the zero "
            "unless --zero is given, and the displacement only with "
            "--displacement. Prints the cell, volume, zero and displacement with "
            "their esds, and with a formula and Z the theoretical density, and "
            "with an Archimedes density the relative density. Writes the peaks to "
            "results/lattice/peaks_STEM.csv and appends one row to "
            "results/lattice/lattice.csv; for a sample of the project file both "
            "go to results/lattice/KEY under the project root, the row to "
            "lattice_KEY.csv, and the cell, space group, formula, Z, form, "
            "radius and Archimedes density not given as options come from the "
            "sample, its instrument and its first structure."
        ),
    )
    parser.add_argument(
        "scan",
        metavar="SCAN",
        help="path to a .xrdml, .xy or .xye file, or a sample key",
    )
    _add_cell_options(parser, "start cell to index from")
    parser.add_argument(
        "--space-group",
        metavar="SG",
        help=(
            "space group whose reflection conditions apply (default: none, or a "
            "sample's structure's)"
        ),
    )
    parser.add_argument(
        "--wavelength",
        type=float,
        metavar="ANGSTROM",
        help="K alpha 1 wavelength (default: the scan's own, or the instrument's)",
    )
    parser.add_argument(
        "--formula", metavar="TEXT", help="formula of one formula unit, for density"
    )
    parser.add_argument("--z", type=float, metavar="N", help="formula units per cell")
    parser.add_argument(
        "--archimedes",
        nargs="+",
        type=float,
        metavar=("RHO", "ESD"),
        help="measured density in g/cm3, and optionally its esd",
    )
    parser.add_argument(
        "--zero",
        type=float,
        metavar="DEG",
        help="hold the zero at this value, in degrees, instead of refining it",
    )
    parser.add_argument(
        "--displacement",
        action="store_true",
        help="refine the specimen displacement (always refined for a pellet)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        metavar="MM",
        help="goniometer radius in mm (default: the instrument's)",
    )
    parser.add_argument(
        "--no-satellites",
        action="store_true",
        help=(
            "exclude every peak flagged as a K alpha 2 satellite, recovering "
            "none that the refined cell puts on a reflection"
        ),
    )
    _add_output_options(parser, "the scan's file stem, or the sample key")
    parser.set_defaults(handler=_run_lattice)


# Le Bail and Rietveld refinement, run by xrdkit.pipeline

# The Rietveld modes, in order, and the parameters of a cell by name, as
# GSAS-II reports them.
RIETVELD_MODES = MODES[1:]
GSAS_CELL = (
    ("a", "length_a"),
    ("b", "length_b"),
    ("c", "length_c"),
    ("alpha", "angle_alpha"),
    ("beta", "angle_beta"),
    ("gamma", "angle_gamma"),
)


def _refinement_sample(key: str) -> tuple[Project, Sample]:
    """The project of the current folder and its sample ``key``, with the
    scan, the instrument parameter file and every CIF the refinement reads
    checked to exist."""
    try:
        path = find_project()
    except FileNotFoundError:
        raise CommandError(
            f"no {PROJECT_FILE} in the current folder or above it to look up the "
            f"sample {key} in; run xrdkit init"
        ) from None
    try:
        project = load_project(path)
    except ValueError as error:
        message = str(error)
        marker = "no such file: "
        if marker in message:
            raise CommandError(marker + message.split(marker, 1)[1]) from None
        raise CommandError(message) from None
    if key not in project.samples:
        raise CommandError(
            f"no sample {key} in {project.root / PROJECT_FILE}; the samples are "
            + (", ".join(project.samples) or "none")
        )
    sample = project.samples[key]
    instrument = project.instruments[sample.instrument]
    if instrument.instprm is None:
        raise CommandError(
            f"instrument {instrument.key} of sample {key} gives no instprm, which a "
            "GSAS-II refinement needs"
        )
    needed = [sample.file, instrument.instprm]
    needed += [
        project.structures[name].cif
        for name in sample.structures
        if project.structures[name].cif is not None
    ]
    for file in needed:
        if not Path(file).is_file():
            raise CommandError(f"no such file: {file}")
    return project, sample


def _printer(args: argparse.Namespace):
    """A reporter printing each line: to stdout, or to stderr under --json so
    that stdout holds the JSON alone."""
    stream = sys.stderr if args.json else sys.stdout

    def report(line: str) -> None:
        print(line, file=stream)

    return report


def _outcome_files(outcome) -> list[str]:
    """The files a finished mode wrote, the result first."""
    paths = outcome.paths
    files = [paths.get("result"), paths.get("gpx"), paths.get("markdown")]
    files += list(paths.get("figures") or [])
    files += [paths.get("histogram"), *(paths.get("reflections") or {}).values()]
    files.append(paths.get("instprm"))
    return [str(path) for path in files if path]


def _no_zero(esd):
    """An esd, or None for none or 0, as _plus_minus takes it."""
    return esd if esd else None


def _lebail_lines(project: Project, outcome) -> list[str]:
    """The cell of each phase with esds, its size and microstrain, the zero or
    displacement, Rwp and chi squared, and where the start cell came from."""
    result = json.loads(Path(outcome.paths["result"]).read_text(encoding="utf-8"))
    inputs = result.get("inputs") or {}
    final = outcome.final or {}
    lines = []
    for phase in final.get("phases") or []:
        six = {short: phase["cell"][gsas] for short, gsas in GSAS_CELL}
        esds = {
            short: _no_zero((phase.get("cell_esd") or {}).get(gsas))
            for short, gsas in GSAS_CELL
        }
        spec = project.structures.get(phase["name"])
        system = (
            load_entry(spec.library).crystal_system
            if spec is not None and spec.library
            else _crystal_system_of(six)
        )
        free = {name: six[name] for name in CELL_PARAMETERS[system]}
        volume = phase["cell"].get("volume")
        if volume is None:
            volume = Cell.from_parameters(system, free).volume
        volume_text = _plus_minus(
            volume, _no_zero((phase.get("cell_esd") or {}).get("volume")), 3
        )
        size, strain = phase.get("size") or {}, phase.get("mustrain") or {}
        size_text = _plus_minus(size.get("value"), _no_zero(size.get("esd")), 4)
        strain_text = _plus_minus(strain.get("value"), _no_zero(strain.get("esd")), 0)
        lines += [
            f"{phase['name']}: {system} cell {_cell_text(free, esds)}",
            f"V = {volume_text} cubic angstrom",
            f"size {size_text} micron, microstrain {strain_text}",
        ]
    if inputs.get("displacement"):
        shift = (final.get("sample") or {}).get("Shift") or {}
        lines.append(
            f"displacement {_plus_minus(shift.get('value'), _no_zero(shift.get('esd')), 4)} "
            "(GSAS-II Shift, micron), zero held"
        )
    else:
        zero = (final.get("instrument") or {}).get("Zero") or {}
        lines.append(
            f"zero {_plus_minus(zero.get('value'), _no_zero(zero.get('esd')), 4)} degrees"
        )
    lines.append(_residual_text(outcome))
    sources = inputs.get("start_cell_source") or {}
    for name, source in sources.items():
        lines.append(f"start cell of {name} from {source}")
    return lines


def _residual_text(outcome) -> str:
    residuals = outcome.residuals
    rwp, chi = residuals.get("rwp"), residuals.get("reduced_chi_squared")
    text = (
        f"Rwp {'—' if rwp is None else f'{rwp:.3f}'} per cent, reduced chi squared "
        f"{'—' if chi is None else f'{chi:.3f}'}"
    )
    fractions = residuals.get("weight_fractions")
    if fractions:
        text += "; weight fractions " + ", ".join(
            f"{name} {_plus_minus(entry.get('value'), _no_zero(entry.get('esd')), 3)}"
            for name, entry in fractions.items()
            if entry and entry.get("value") is not None
        )
    return text


def _run_refinement(
    args: argparse.Namespace,
    project: Project,
    sample: Sample,
    modes: Sequence[str],
    options: Options,
) -> int:
    """Run ``modes`` and print the files written and the outcome; 1 when a
    mode failed, with the path of its failure record."""
    outcomes = run_sequence(project, sample, modes, options, _printer(args))
    written = []
    for outcome in outcomes:
        if outcome.error is None:
            written += _outcome_files(outcome)
    failed = next((outcome for outcome in outcomes if outcome.error), None)
    if failed is not None:
        written.append(failed.paths["failure"])
    if outcomes and outcomes[-1].paths.get("summary"):
        written.append(outcomes[-1].paths["summary"])
    if args.json:
        print(
            json.dumps(
                _jsonable(
                    {
                        "files": written,
                        "sample": sample.key,
                        "outcomes": [asdict(outcome) for outcome in outcomes],
                    }
                ),
                indent=2,
            )
        )
    else:
        for path in written:
            print(path)
        for outcome in outcomes:
            if outcome.error is not None:
                continue
            if outcome.mode == "lebail":
                for line in _lebail_lines(project, outcome):
                    print(line)
            else:
                print(
                    f"{outcome.mode}: {', '.join(outcome.accepted)}; "
                    f"{_residual_text(outcome)}"
                )
    if failed is not None:
        print(
            f"xrdkit {args.command}: {failed.mode} failed: {failed.error}; see "
            f"{failed.paths['failure']}",
            file=sys.stderr,
        )
        return 1
    return 0


# What the instrument command records as its method, for the CSV row.
INSTRUMENT_METHOD = (
    "Caglioti width fit of a standard scan, then a GSAS-II refinement of the "
    "zero and the profile terms with the cell held at the certified one"
)


def _instrument_column(key: str) -> str:
    """A refined parameter's name as a column name: ``SH/L`` as ``sh_l``."""
    return key.lower().replace("/", "_")


def _instrument_table(
    key: str,
    wavelength: Sequence[float],
    ka2: bool,
    radius: float | None,
    instprm: str,
) -> list[str]:
    """The lines of the ``[instruments.<key>]`` table the command appends."""
    header = key if BARE_KEY.fullmatch(key) else toml_string(key)
    lines = [
        f"[instruments.{header}]",
        "wavelength = [" + ", ".join(_toml_number(v) for v in wavelength) + "]",
        f"ka2 = {str(bool(ka2)).lower()}",
    ]
    if radius is not None:
        lines.append(f"radius = {_toml_number(radius)}")
    lines.append(f"instprm = {toml_string(instprm)}")
    return lines


def _instrument_project(args: argparse.Namespace, item: _Input):
    """The project the command writes into, and its file, or ``(None, None)``.

    A bad project file is only fatal where ``--name`` has to be appended to it;
    otherwise the command falls back to the working folder.
    """
    if item.project is not None:
        return item.project, item.project.root / PROJECT_FILE
    try:
        path = find_project()
    except FileNotFoundError:
        return None, None
    try:
        return load_project(path), path
    except ValueError as error:
        if args.name is not None:
            raise CommandError(str(error)) from None
        return None, path


def _run_instrument(args: argparse.Namespace) -> int:
    # Every input and option is checked, and GSAS-II located, before anything
    # is written: a standard is refined once and its file read for months.
    cell = _cell_option(args)
    if cell is None:
        raise CommandError(
            "--cell is required: the certified cell of the standard, one number "
            "for a cubic one"
        )
    low, high = args.window
    if not low < high:
        raise CommandError(
            f"--window takes MIN MAX with MIN below MAX, not {low:g} {high:g}"
        )
    cif = Path(args.cif)
    if not cif.is_file():
        raise CommandError(f"no such file: {args.cif}")
    item = _resolve(args.scan, {})
    (scan,) = _read_inputs([item], args.wavelength)
    _require_wavelength(scan, item)

    project, project_path = _instrument_project(args, item)
    stem = args.stem or item.path.stem
    if args.out is not None:
        folder = Path(args.out)
    elif project is not None:
        folder = project.root / "data" / "standards"
    else:
        folder = Path(".")
    instprm_path = (folder / f"{stem}.instprm").resolve()

    relative = None
    if args.name is not None:
        if project is None:
            raise CommandError(
                f"no {PROJECT_FILE} in the current folder or above it to add "
                f"[instruments.{args.name}] to; run xrdkit init"
            )
        if args.name in project.instruments:
            raise CommandError(
                f"instrument {args.name!r} is in {project_path} already; give "
                "another --name"
            )
        try:
            relative = instprm_path.relative_to(project.root)
        except ValueError:
            raise CommandError(
                f"{instprm_path} is outside the project folder {project.root}; "
                "--name needs the instrument file inside it"
            ) from None
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        raise CommandError(str(error)) from None

    report = _printer(args)
    try:
        fit = fit_instrument_widths(scan, (low, high))
    except ValueError as error:
        raise CommandError(str(error)) from None
    report(f"peaks   {fit.n_peaks} fitted from {low:g} to {high:g} degrees")
    report(f"U = {fit.u:.4f}, V = {fit.v:.4f}, W = {fit.w:.4f} degrees squared")
    report(f"rms of the width fit {fit.rms:.5f} degrees")

    try:
        refined = refine_instrument(
            fit,
            item.path,
            cif,
            args.phase or cif.stem,
            (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma),
            folder,
            stem,
            folder / "work" / stem,
            install,
            args.radius,
        )
    except (Gsas2Error, KeyError, OSError, ValueError) as error:
        raise CommandError(
            f"the GSAS-II refinement of the standard failed: {error}"
        ) from None
    report("stages  " + ", ".join(refined.stages))
    report(f"Rwp {refined.rwp:.3f} per cent, GOF {refined.gof:.3f}")
    for key, entry in refined.parameters.items():
        esd = "held" if entry["esd"] is None else f"esd {entry['esd']:.3g}"
        report(f"  {key:5s} {entry['value']:11.6g}  {esd}")

    row = {
        "scan": str(item.path),
        "cif": str(cif),
        "phase": args.phase or cif.stem,
        "wavelength_angstrom": fit.wavelength,
        **{
            name: value
            for name, value in _cell_columns(cell, None).items()
            if not name.startswith("esd_")
        },
        "window_min_deg": fit.window[0],
        "window_max_deg": fit.window[1],
        "n_peaks_fitted": fit.n_peaks,
        "fit_u_deg2": fit.u,
        "fit_v_deg2": fit.v,
        "fit_w_deg2": fit.w,
        "fit_rms_deg": fit.rms,
    }
    for key in REFINED_KEYS:
        entry = refined.parameters.get(key) or {}
        name = _instrument_column(key)
        row[f"refined_{name}"] = entry.get("value")
        row[f"esd_{name}"] = entry.get("esd")
    row.update(
        {
            "rwp_percent": refined.rwp,
            "gof": refined.gof,
            "radius_mm": args.radius,
            "instprm": str(refined.instprm),
            "method": INSTRUMENT_METHOD,
            "date": datetime.datetime.now(datetime.UTC).astimezone().date().isoformat(),
            "xrdkit_version": __version__,
        }
    )
    written = [refined.instprm, _append_row(folder / f"{stem}_instrument.csv", row)]

    if args.name is not None:
        kalpha2 = kalpha2_wavelength(item.path)
        wavelengths = [fit.wavelength] if kalpha2 is None else [fit.wavelength, kalpha2]
        lines = _instrument_table(
            args.name,
            wavelengths,
            kalpha2 is not None,
            args.radius,
            relative.as_posix(),
        )
        # As add-sample does: the file as it is, byte for byte, with the table
        # appended in its own line endings, and only once it still reads.
        original = project_path.read_bytes()
        newline = "\r\n" if b"\r\n" in original else "\n"
        text = original.decode("utf-8")
        ending = "" if not text or text.endswith("\n") else newline
        addition = ending + newline + newline.join(lines) + newline
        try:
            load_project_text(text + addition, project_path)
        except ValueError as error:
            raise CommandError(str(error)) from None
        project_path.write_bytes(original + addition.encode("utf-8"))
        report("\n".join(lines))
        written.append(project_path)

    return _finish(
        args,
        written,
        n_peaks_fitted=fit.n_peaks,
        fit_u_deg2=fit.u,
        fit_v_deg2=fit.v,
        fit_w_deg2=fit.w,
        fit_rms_deg=fit.rms,
        stages=list(refined.stages),
        rwp_percent=refined.rwp,
        gof=refined.gof,
        parameters=refined.parameters,
        instrument=args.name,
    )


def _add_instrument(subparsers) -> None:
    parser = subparsers.add_parser(
        "instrument",
        help="instrument parameter file from a standard scan, refined in GSAS-II",
        description=(
            "Fit the widths of a standard's reflections to the Caglioti "
            "relation, refine the zero and the profile terms in GSAS-II with "
            "the cell held at --cell, and write the instrument parameter file "
            "every later Le Bail and Rietveld refinement reads, with a record "
            "of how it was made beside it. The starting file keeps the same "
            "stem with _start on the end, and the GSAS-II project and its logs "
            "go under work/ in the same folder. Every input and option is "
            "checked, and GSAS-II located, before anything is written."
        ),
    )
    parser.add_argument(
        "scan",
        metavar="SCAN",
        help="the standard's scan: a .xrdml, .xy or .xye file, or a sample key",
    )
    parser.add_argument(
        "--cif", metavar="CIF", required=True, help="the standard's CIF"
    )
    _add_cell_options(parser, "the standard's certified cell")
    parser.add_argument(
        "--phase",
        metavar="NAME",
        help="phase name for GSAS-II (default: the CIF's file stem)",
    )
    parser.add_argument(
        "--wavelength",
        type=float,
        metavar="ANGSTROM",
        help="K alpha 1 wavelength (default: the scan's own, or the instrument's)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        metavar="MM",
        help=(
            "goniometer radius in mm, written into the instrument file so that "
            "GSAS-II refines a specimen displacement in the right geometry"
        ),
    )
    parser.add_argument(
        "--window",
        nargs=2,
        type=float,
        default=list(WIDTH_WINDOW),
        metavar=("MIN", "MAX"),
        help=(
            "two theta range the widths are fitted over (default: "
            f"{WIDTH_WINDOW[0]:g} {WIDTH_WINDOW[1]:g})"
        ),
    )
    parser.add_argument(
        "--name",
        metavar="KEY",
        help=(
            f"also append an [instruments.KEY] table to {PROJECT_FILE}, with the "
            "wavelengths, ka2, the radius where given and the path of the file "
            "written; refused if the key is there already"
        ),
    )
    parser.add_argument(
        "--stem", help="name the output files take (default: the scan's file stem)"
    )
    parser.add_argument(
        "--out",
        metavar="DIR",
        help=(
            "folder for the instrument file and its record (default: "
            f"data/standards under the project root when {PROJECT_FILE} is "
            "found, else the current folder)"
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="print the files written and the numbers as JSON, progress on stderr",
    )
    parser.set_defaults(handler=_run_instrument)


# What the phases command records as its method, for its record row.
PHASES_METHOD = (
    "COD search by element set, each candidate's pattern simulated with "
    "pymatgen and weighed against the observed peaks"
)


def _phases_elements(args: argparse.Namespace, item: _Input) -> list[str]:
    """The elements to search on: ``--elements``, or a sample's structures.

    A sample's are the union of the compositions of the structures it lists,
    in the order they first appear, which is the order the formulae give.
    """
    if args.elements:
        return list(dict.fromkeys(args.elements))
    if item.sample is None:
        raise CommandError(
            "--elements is required for a scan named directly; a sample key takes "
            "them from the compositions of its structures"
        )
    elements: dict[str, None] = {}
    for key in item.sample.structures:
        composition = item.project.structures[key].composition
        try:
            elements.update(dict.fromkeys(parse_formula(composition)))
        except ValueError as error:
            raise CommandError(f"structures.{key}: {error}") from None
    if not elements:
        raise CommandError(
            f"the structures of sample {item.name} give no elements; give --elements"
        )
    return list(elements)


@contextmanager
def _quiet_pymatgen() -> Iterator[None]:
    """Hold pymatgen's own warnings back for the length of a call into it.

    pymatgen reports a doubtful CIF through the ``warnings`` module, as a
    plain ``UserWarning`` with no category of its own: an entry whose formula
    cannot be checked against its sites, or one with no symmetry operations,
    which it falls back to P1 for. Several COD entries draw one, and the
    lines land on stderr in the middle of the command's own output, where
    they read as errors of xrdkit's.

    The filter goes on only around a call that reaches pymatgen, so every
    other line the command writes to stderr is untouched, and a library
    caller of :mod:`xrdkit.phases` still sees the warnings.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        yield


def _run_phases(args: argparse.Namespace) -> int:
    # Every option, the scan, the extra and the folders are settled before any
    # request is made: a search takes a while and the COD is someone else's.
    low, high = args.window
    if not low < high:
        raise CommandError(
            f"--window takes MIN MAX with MIN below MAX, not {low:g} {high:g}"
        )
    if args.tolerance <= 0:
        raise CommandError(
            f"--tolerance must be greater than 0, not {args.tolerance:g}"
        )
    if args.max_candidates is not None and args.max_candidates < 1:
        raise CommandError(
            f"--max-candidates must be 1 or more, not {args.max_candidates}"
        )
    item = _resolve(args.scan, {})
    (scan,) = _read_inputs([item], args.wavelength)
    wavelength = _require_wavelength(scan, item)
    elements = _phases_elements(args, item)
    try:
        with _quiet_pymatgen():
            require_phases_extra()
    except MissingPhasesExtra:
        raise CommandError(
            "phase identification needs the phases extra; install with "
            "pip install xrdkit[phases]"
        ) from None

    project = item.project
    if project is None:
        try:
            project = load_project(find_project())
        except (FileNotFoundError, ValueError):
            project = None
    root = project.root if (project is not None and args.out is None) else Path(".")
    cifs = (Path(args.out) if args.out is not None else root) / "cifs" / "cod"
    stem = args.stem or item.name
    if args.out is not None:
        folder = Path(args.out) / "results" / "phases" / item.name
    elif project is not None:
        folder = results_dir(project, "phases", item.name)
    else:
        folder = Path("results") / "phases" / item.name

    report = _printer(args)
    try:
        observed = observed_peaks(scan, (low, high), args.zero)
    except ValueError as error:
        raise CommandError(str(error)) from None
    if not observed:
        raise CommandError(
            f"no peaks found in {item.path} between {low:g} and {high:g} degrees"
        )
    report(f"{len(observed)} peaks observed from {low:g} to {high:g} degrees")

    try:
        records = sorted(
            cod_search(elements, exact=True, space_group=args.space_group),
            key=lambda record: record.cod_id,
        )
    except (URLError, ValueError) as error:
        raise CommandError(f"the COD search failed: {error}") from None
    report(f"{len(records)} COD entries made of exactly {', '.join(elements)}")
    if not records:
        raise CommandError(
            "the COD holds no entry of exactly those elements; widen --elements "
            "or drop --space-group"
        )
    if args.max_candidates is not None:
        records = records[: args.max_candidates]
    for record in records:
        report(f"  {record.cod_id}  {record.formula:22s} {record.space_group}")

    try:
        fetched = fetch_candidates(records, cifs)
    except (URLError, ValueError) as error:
        raise CommandError(f"fetching a CIF from the COD failed: {error}") from None
    index = write_cif_index(records, cifs / "index.csv", notes=stem)
    report(f"index written to {index}")

    try:
        with _quiet_pymatgen():
            candidates = rank_candidates(
                records, cifs, observed, wavelength, (low, high), args.tolerance
            )
    except MissingPhasesExtra as error:
        raise CommandError(str(error)) from None
    except ValueError as error:
        raise CommandError(
            f"a candidate's pattern could not be simulated: {error}"
        ) from None
    for candidate in candidates:
        report(
            f"  {candidate.cod_id}  explained {candidate.explained:2d}/"
            f"{len(observed)}  missing {candidate.missing:2d}  "
            f"score {candidate.score:3d}  "
            f"{'rejected, ' + candidate.rejected if candidate.rejected else ''}".rstrip()
        )

    kept = [c for c in candidates if not c.rejected] or list(candidates)
    by_id = {candidate.cod_id: candidate for candidate in candidates}
    main = args.main or kept[0].cod_id
    # A peak no candidate accounts for is worth naming, so --main may point at
    # a phase outside the set searched: it is fetched like any other.
    if main in by_id:
        main_cif = by_id[main].cif
    else:
        if not re.fullmatch(r"\d{7}", main):
            raise CommandError(f"--main takes a seven digit COD id, not {main!r}")
        main_cif = cifs / f"{main}.cif"
        if not main_cif.is_file():
            try:
                main_cif = cod_fetch(main, cifs)
            except (URLError, ValueError) as error:
                raise CommandError(
                    f"fetching the main phase {main} from the COD failed: {error}"
                ) from None
            fetched.append(main_cif)
    explained_positions = sorted(
        {
            position
            for candidate in candidates
            for position in candidate.explained_positions
        }
    )
    try:
        with _quiet_pymatgen():
            unexplained = attribute_unexplained(
                observed,
                explained_positions,
                main_cif,
                wavelength,
                main,
                (low, high),
                args.tolerance,
            )
    except MissingPhasesExtra as error:
        raise CommandError(str(error)) from None
    for peak in unexplained:
        if peak.identified:
            report(
                f"{peak.two_theta:.3f} degrees, d = {peak.d_spacing:.4f} angstrom: "
                f"{peak.phase} {peak.hkl} at {peak.reflection_two_theta:.3f}, "
                f"{peak.intensity:.1f} per cent"
            )
        else:
            report(
                f"{peak.two_theta:.3f} degrees, d = {peak.d_spacing:.4f} angstrom: "
                "unidentified"
            )

    folder.mkdir(parents=True, exist_ok=True)
    candidates_path = folder / f"phases_{stem}.csv"
    _write_rows(
        candidates_path,
        [
            {
                "rank": candidate.rank,
                "cod_id": candidate.cod_id,
                "formula": candidate.formula,
                "space_group": candidate.space_group,
                "n_explained": candidate.explained,
                "n_missing": candidate.missing,
                "score": candidate.score,
                "rejected": candidate.rejected,
                "cif": str(candidate.cif),
            }
            for candidate in candidates
        ],
    )
    unexplained_path = folder / f"phases_{stem}_unexplained.csv"
    _write_rows(
        unexplained_path,
        [
            {
                "two_theta_deg": peak.two_theta,
                "d_angstrom": peak.d_spacing,
                "phase": peak.phase,
                "hkl": "" if peak.hkl is None else " ".join(map(str, peak.hkl)),
                "reflection_two_theta_deg": peak.reflection_two_theta,
                "intensity_percent": peak.intensity,
            }
            for peak in unexplained
        ],
        columns=[
            "two_theta_deg",
            "d_angstrom",
            "phase",
            "hkl",
            "reflection_two_theta_deg",
            "intensity_percent",
        ],
    )
    record_path = _append_row(
        folder / f"phases_{stem}_record.csv",
        {
            "scan": str(item.path),
            "sample": item.sample.key if item.sample else None,
            "wavelength_angstrom": wavelength,
            "elements": " ".join(elements),
            "space_group": args.space_group,
            "zero_deg": args.zero,
            "window_min_deg": low,
            "window_max_deg": high,
            "tolerance_deg": args.tolerance,
            "n_peaks_observed": len(observed),
            "n_candidates_searched": len(records),
            "n_candidates_fetched": len(fetched),
            "n_unexplained": len(unexplained),
            "main": main,
            "index": str(index),
            "method": PHASES_METHOD,
            "date": datetime.datetime.now(datetime.UTC).astimezone().date().isoformat(),
            "xrdkit_version": __version__,
        },
    )
    written = [*fetched, index, candidates_path, unexplained_path, record_path]
    return _finish(
        args,
        written,
        n_peaks_observed=len(observed),
        elements=elements,
        candidates=[
            asdict(candidate) | {"cif": str(candidate.cif)} for candidate in candidates
        ],
        main=main,
        unexplained=[asdict(peak) for peak in unexplained],
    )


def _write_rows(path: Path, rows: list[dict], columns: list[str] | None = None) -> Path:
    """Write ``rows`` to ``path`` as a CSV, header first, replacing the file."""
    columns = columns or (list(rows[0]) if rows else [])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "" if value is None else value for key, value in row.items()}
            )
    return path


def _add_phases(subparsers) -> None:
    parser = subparsers.add_parser(
        "phases",
        help="identify the phases of a scan against the COD",
        description=(
            "Search the Crystallography Open Database for entries made of "
            "exactly the elements given, fetch their CIFs, simulate each "
            "pattern and weigh it against the peaks of the scan, then say what "
            "explains every peak left over. Needs the phases extra, "
            "pip install xrdkit[phases], and an internet connection. The CIFs "
            "and their index are kept under cifs/cod, and the ranking, the "
            "peaks left over and a record of the run go to "
            "results/phases/KEY."
        ),
    )
    parser.add_argument(
        "scan",
        metavar="SCAN",
        help="path to a .xrdml, .xy or .xye file, or a sample key",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        metavar="EL",
        help=(
            "element symbols the entries are made of, at most eight (default: "
            "the elements of a sample's structures' compositions)"
        ),
    )
    parser.add_argument(
        "--space-group",
        metavar="SYMBOL",
        help="only entries of this space group symbol (default: any)",
    )
    parser.add_argument(
        "--zero",
        type=float,
        default=0.0,
        metavar="DEG",
        help="zero offset taken off every observed peak (default: 0)",
    )
    parser.add_argument(
        "--window",
        nargs=2,
        type=float,
        default=list(PHASES_WINDOW),
        metavar=("MIN", "MAX"),
        help=(
            "two theta range the peaks are taken and the patterns simulated "
            f"over (default: {PHASES_WINDOW[0]:g} {PHASES_WINDOW[1]:g})"
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=PHASES_TOLERANCE,
        metavar="DEG",
        help=(
            "how far an observed peak may lie from a simulated reflection "
            f"(default: {PHASES_TOLERANCE:g})"
        ),
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        metavar="N",
        help="fetch and weigh only the first N entries found (default: all)",
    )
    parser.add_argument(
        "--main",
        metavar="COD_ID",
        help=(
            "the phase the peaks left over are attributed to (default: the top "
            "ranked candidate that was not rejected)"
        ),
    )
    parser.add_argument(
        "--wavelength",
        type=float,
        metavar="ANGSTROM",
        help="K alpha 1 wavelength (default: the scan's own, or the instrument's)",
    )
    _add_output_options(parser, "the scan's file stem, or the sample key")
    parser.set_defaults(handler=_run_phases)


def _two_theta_option(args: argparse.Namespace, sample: Sample) -> tuple | None:
    """``--two-theta`` as a pair, or None for the project file's range.

    Refused here, before anything is written: a range the wrong way round,
    and one that does not reach the scan at all. A range that overhangs the
    scan at either end is allowed and is clipped to it by
    :func:`~xrdkit.pipeline.resolve_inputs`, which is how the whole scan is
    fitted where the sample's refine table narrows it."""
    if args.two_theta is None:
        return None
    low, high = (float(value) for value in args.two_theta)
    if not low < high:
        raise CommandError(
            f"--two-theta takes MIN MAX with MIN below MAX, not {low:g} {high:g}"
        )
    try:
        scan = read_scan(sample.file)
    except (OSError, ValueError) as error:
        raise CommandError(f"{sample.file}: cannot be read: {error}") from None
    start, end = float(np.min(scan.two_theta)), float(np.max(scan.two_theta))
    if high <= start or low >= end:
        raise CommandError(
            f"--two-theta {low:g} {high:g} lies outside the scan's {start:.2f} to "
            f"{end:.2f} degrees"
        )
    return (low, high)


def _add_two_theta_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--two-theta",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        help=(
            "two theta range to fit, in degrees, clipped to the scan "
            "(default: the project file's refine two_theta)"
        ),
    )


def _run_lebail(args: argparse.Namespace) -> int:
    # Every option and input is checked before anything is written.
    cell = _cell_option(args)
    if args.zero is not None and args.displacement:
        raise CommandError(
            "--zero and --displacement conflict: the displacement is refined with "
            "the zero held at the instrument's"
        )
    project, sample = _refinement_sample(args.sample)
    options = Options(
        cell=None
        if cell is None
        else {name: getattr(cell, name) for name, _ in GSAS_CELL},
        zero=args.zero,
        displacement=True if args.displacement else None,
        two_theta=_two_theta_option(args, sample),
        out=_out_folder(args),
    )
    return _run_refinement(args, project, sample, ["lebail"], options)


def _add_lebail(subparsers) -> None:
    parser = subparsers.add_parser(
        "lebail",
        help="Le Bail extraction of a sample in GSAS-II, for its cell",
        description=(
            "Refine a sample of the project file by Le Bail extraction in GSAS-II: "
            "background and scale, then the zero (or for a pellet, or with "
            "--displacement, the specimen displacement), the cell, the size and, "
            "always, the microstrain as a test. Every phase of the sample is "
            "extracted. The first phase's cell starts from --cell, else the "
            "sample's lattice results, else its structure's cell. Writes the result "
            "JSON, the GSAS-II project, the fitted pattern and reflections, the "
            "instrument parameters as used, a figure and lebail.md to "
            "results/lebail/KEY under the project root, or to --out; prints each "
            "file, then the cell, size, microstrain, zero or displacement, Rwp and "
            "chi squared."
        ),
    )
    parser.add_argument(
        "sample", metavar="SAMPLE", help="a sample key of the project file"
    )
    _add_cell_options(parser, "start cell of the first phase")
    parser.add_argument(
        "--zero",
        type=float,
        metavar="DEG",
        help="start the zero at this value, in degrees (default: the instrument's)",
    )
    parser.add_argument(
        "--displacement",
        action="store_true",
        help="refine the specimen displacement, the zero held (always for a pellet)",
    )
    _add_two_theta_option(parser)
    _add_refinement_output_options(parser, "lebail")
    parser.set_defaults(handler=_run_lebail)


def _out_folder(args: argparse.Namespace) -> Path | None:
    """``--out`` as an absolute path, resolved against the current folder:
    the GSAS-II driver runs in a work folder of its own, where a relative
    path would name somewhere else."""
    return None if args.out is None else Path(args.out).resolve()


def _add_refinement_output_options(
    parser: argparse.ArgumentParser, command: str
) -> None:
    parser.add_argument(
        "--out",
        metavar="DIR",
        help=(
            "folder every file is written to (default: results/"
            f"{command}/KEY under the project root)"
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="print the files written and the outcome as JSON",
    )


def _run_rietveld(args: argparse.Namespace) -> int:
    first, last = RIETVELD_MODES.index(args.first), RIETVELD_MODES.index(args.through)
    if first > last:
        raise CommandError(f"--from {args.first} comes after --through {args.through}")
    modes = list(RIETVELD_MODES[first : last + 1])
    axis = args.preferred_orientation
    if axis is not None:
        if "fixed_atoms" not in modes:
            raise CommandError(
                "--preferred-orientation applies to the fixed_atoms mode, which "
                f"--from {args.first} leaves out"
            )
        if not any(axis):
            raise CommandError("--preferred-orientation takes an axis H K L not all 0")
    project, sample = _refinement_sample(args.sample)
    options = Options(
        mustrain=args.mustrain,
        preferred_orientation=None if axis is None else tuple(axis),
        two_theta=_two_theta_option(args, sample),
        out=_out_folder(args),
    )
    before = MODES[MODES.index(modes[0]) - 1]
    path = mode_paths(project, sample, before, options)["result"]
    if not path.is_file():
        out = f" --out {args.out}" if args.out else ""
        command = (
            f"xrdkit lebail {sample.key}{out}"
            if before == "lebail"
            else f"xrdkit rietveld {sample.key} --through {before}{out}"
        )
        raise CommandError(f"no such file: {path}; {command} writes it")
    existing = [
        str(result)
        for mode in modes
        if (result := mode_paths(project, sample, mode, options)["result"]).is_file()
    ]
    if existing and not args.overwrite:
        raise CommandError(
            f"{', '.join(existing)} {'exists' if len(existing) == 1 else 'exist'} "
            "already; give --overwrite to replace the results of the modes run, "
            "or --out DIR to write them elsewhere"
        )
    return _run_refinement(args, project, sample, modes, options)


def _add_rietveld(subparsers) -> None:
    parser = subparsers.add_parser(
        "rietveld",
        help="Rietveld refinement of a sample in GSAS-II, mode by mode",
        description=(
            "Refine a sample of the project file in GSAS-II from its Le Bail "
            "result, in the modes fixed_atoms (the structure in at its nominal "
            "composition, the atoms fixed but for one Uiso), coordinates (the "
            "coordinates kind of site by kind) and occupancies (the exchanged "
            "elements traded between their sites), each from the result of the "
            "one before, from --from through --through. A mode that fails is "
            "written up in failure.md and the modes after it are not run. Writes "
            "each mode's result JSON, GSAS-II project, fitted pattern and "
            "reflections, instrument parameters, figure and MODE.md, and "
            "summary.md, to results/rietveld/KEY under the project root, or to "
            "--out, and refuses to replace a mode's existing result unless "
            "--overwrite is given; prints each file, then a line per mode with "
            "its accepted stages, Rwp and chi squared."
        ),
    )
    parser.add_argument(
        "sample", metavar="SAMPLE", help="a sample key of the project file"
    )
    parser.add_argument(
        "--from",
        dest="first",
        metavar="MODE",
        choices=RIETVELD_MODES,
        default=RIETVELD_MODES[0],
        help=f"first mode to run, one of {', '.join(RIETVELD_MODES)} (default: fixed_atoms)",
    )
    parser.add_argument(
        "--through",
        metavar="MODE",
        choices=RIETVELD_MODES,
        default=RIETVELD_MODES[-1],
        help="last mode to run (default: occupancies)",
    )
    parser.add_argument(
        "--mustrain",
        action="store_true",
        help="refine the microstrain with the size (held at zero by default)",
    )
    parser.add_argument(
        "--preferred-orientation",
        nargs=3,
        type=int,
        metavar=("H", "K", "L"),
        help="refine a March-Dollase ratio about this axis in the fixed_atoms mode",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "replace the results of the modes run where they exist already "
            "(refused by default)"
        ),
    )
    _add_two_theta_option(parser)
    _add_refinement_output_options(parser, "rietveld")
    parser.set_defaults(handler=_run_rietveld)


# The folders xrdkit init makes beside the project file.
INIT_FOLDERS = ("data/raw", "cifs", "results")


def _run_init(args: argparse.Namespace) -> int:
    root = Path.cwd()
    path = root / PROJECT_FILE
    if path.exists():
        raise CommandError(f"{path} already exists; it is left as it is")
    name = (root.name or "project") if args.name is None else args.name
    if not name.strip():
        raise CommandError("the project name must not be blank")
    path.write_text(project_template(name), encoding="utf-8")
    for folder in INIT_FOLDERS:
        (root / folder).mkdir(parents=True, exist_ok=True)
    print(path)
    return 0


def _add_init(subparsers) -> None:
    parser = subparsers.add_parser(
        "init",
        help=f"start a project: write {PROJECT_FILE} and make its folders",
        description=(
            f"Write {PROJECT_FILE} in the current folder, with [project] filled "
            "in and a commented example of an instrument, a structure and a "
            "sample, and make data/raw, cifs and results if they are missing. "
            f"An existing {PROJECT_FILE} is never overwritten."
        ),
    )
    parser.add_argument(
        "--name",
        metavar="TEXT",
        help="the project name (default: the name of the current folder)",
    )
    parser.set_defaults(handler=_run_init)


# A key TOML takes bare in a table header; any other is quoted.
BARE_KEY = re.compile(r"[A-Za-z0-9_-]+")


def _toml_number(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else repr(float(value))


def _sample_table(args: argparse.Namespace, key: str, file: str, instrument: str):
    """The lines of the ``[samples.<key>]`` table ``add-sample`` writes."""
    header = key if BARE_KEY.fullmatch(key) else toml_string(key)
    lines = [
        f"[samples.{header}]",
        f"file = {toml_string(file)}",
        f"instrument = {toml_string(instrument)}",
        "structures = [" + ", ".join(map(toml_string, args.structures)) + "]",
    ]
    if args.stage is not None:
        lines.append(f"stage = {toml_string(args.stage)}")
    lines.append(f"form = {toml_string(args.form)}")
    if args.temperature_c is not None:
        lines.append(f"temperature_c = {_toml_number(args.temperature_c)}")
    if args.archimedes is not None:
        lines.append(f"archimedes = {_toml_number(args.archimedes)}")
    if args.notes is not None:
        lines.append(f"notes = {toml_string(args.notes)}")
    return lines


def _run_add_sample(args: argparse.Namespace) -> int:
    try:
        path = find_project()
        project = load_project(path)
    except (FileNotFoundError, ValueError) as error:
        raise CommandError(str(error)) from None

    file = Path(args.file).resolve()
    try:
        relative = file.relative_to(project.root)
    except ValueError:
        raise CommandError(
            f"{args.file} is outside the project folder {project.root}; move it "
            "inside, under data/raw say"
        ) from None
    if not file.is_file():
        raise CommandError(f"no such file: {args.file}")

    key = args.name if args.name is not None else file.stem
    if key in project.samples:
        raise CommandError(f"sample {key!r} is in {path} already; give another --name")
    instrument = args.instrument
    if instrument is None:
        if len(project.instruments) != 1:
            known = ", ".join(project.instruments) or "none"
            raise CommandError(
                f"give --instrument; the project has {len(project.instruments)} "
                f"instruments, not one ({known})"
            )
        (instrument,) = project.instruments
    lines = _sample_table(args, key, relative.as_posix(), instrument)

    # The table goes after the file as it is, byte for byte and in its own line
    # endings, and only once the whole text has been checked.
    original = path.read_bytes()
    newline = "\r\n" if b"\r\n" in original else "\n"
    text = original.decode("utf-8")
    ending = "" if not text or text.endswith("\n") else newline
    addition = ending + newline + newline.join(lines) + newline
    try:
        load_project_text(text + addition, path)
    except ValueError as error:
        raise CommandError(str(error)) from None
    path.write_bytes(original + addition.encode("utf-8"))
    print("\n".join(lines))
    return 0


def _add_add_sample(subparsers) -> None:
    parser = subparsers.add_parser(
        "add-sample",
        help=f"add a sample to {PROJECT_FILE}",
        description=(
            f"Add a [samples.KEY] table to the end of the {PROJECT_FILE} found from "
            "the current folder, leaving the rest of the file as it is, and print "
            "the table. FILE must lie inside the project folder and is stored "
            "relative to it. The file is checked with the table added before "
            "anything is written, and an existing key is refused."
        ),
    )
    parser.add_argument("file", metavar="FILE", help="the sample's scan")
    parser.add_argument(
        "--structure",
        dest="structures",
        action="append",
        required=True,
        metavar="KEY",
        help="a structure key of the project; repeat it for each phase",
    )
    parser.add_argument(
        "--name", metavar="KEY", help="the sample key (default: FILE's stem)"
    )
    parser.add_argument(
        "--instrument",
        metavar="KEY",
        help="an instrument key (default: the only instrument, if there is one)",
    )
    parser.add_argument("--stage", metavar="TEXT", help="the processing stage")
    parser.add_argument(
        "--form", choices=FORMS, default="powder", help="(default: powder)"
    )
    parser.add_argument(
        "--temperature-c", type=float, metavar="N", help="temperature in degrees C"
    )
    parser.add_argument(
        "--archimedes", type=float, metavar="N", help="measured density in g/cm3"
    )
    parser.add_argument("--notes", metavar="TEXT", help="free text")
    parser.set_defaults(handler=_run_add_sample)


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the ``xrdkit`` command."""
    parser = argparse.ArgumentParser(
        prog="xrdkit", description="X-ray diffraction analysis toolkit."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
    _add_init(subparsers)
    _add_add_sample(subparsers)
    _add_check(subparsers)
    _add_plot(subparsers)
    _add_stack(subparsers)
    _add_density(subparsers)
    _add_lattice(subparsers)
    _add_instrument(subparsers)
    _add_phases(subparsers)
    _add_lebail(subparsers)
    _add_rietveld(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ``xrdkit`` command with ``argv``, or the process arguments."""
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except CommandError as error:
        print(f"xrdkit {args.command}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
