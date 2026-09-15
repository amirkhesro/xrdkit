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
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from xrdkit import __version__
from xrdkit.cell import Cell
from xrdkit.density import (
    cell_volume,
    formula_mass,
    relative_density,
    theoretical_density,
)
from xrdkit.indexing import (
    DEFAULT_SPACE_GROUP,
    DEFAULT_ZERO_OFFSET,
    SUPPORTED_SPACE_GROUPS,
    IndexedPeak,
    index_and_refine,
    indexed_to_csv,
    indexing_summary,
)
from xrdkit.io import XRDScan, read_xrdml
from xrdkit.library import load_entry
from xrdkit.peaks import exclude_kalpha2, find_peaks, flag_kalpha2, peaks_to_csv
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
            scans.append(read_xrdml(path))
        except (ET.ParseError, ValueError) as error:
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
            "Read a .xrdml scan, print its range, step, counting time and "
            "intensities, and say whether it is good enough for plotting, phase "
            "identification, a Le Bail fit and a Rietveld refinement, by the "
            "criteria of the user guide, Section 2. For a sample of the project "
            "file the report is headed by its key and composition and is also "
            "written to results/check/KEY under the project root."
        ),
    )
    parser.add_argument(
        "scan", metavar="SCAN", help="path to a .xrdml file, or a sample key"
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


def _structure_cell(
    item: _Input, space_group: str | None
) -> tuple[list[float] | None, str | None, str | None]:
    """The start cell and space group to index a sample from, taken from its
    first structure, with ``space_group`` given on the command line winning;
    or no cell, and the note saying why, when that structure's cell cannot be
    indexed in this version."""
    spec = item.structure
    try:
        cell = resolved_cell(spec)
    except ValueError as error:
        return None, None, f"{error}, so the peaks are not labelled with hkl"
    if spec.library is not None:
        entry = load_entry(spec.library)
        system, structure_group = entry.crystal_system, entry.space_group
    else:
        system, structure_group = _crystal_system_of(cell), None
    if system != "tetragonal":
        return (
            None,
            None,
            (
                f"hkl labelling for a {system} cell arrives in the next release; "
                "plotted without hkl labels"
            ),
        )
    chosen = space_group or structure_group
    if chosen not in SUPPORTED_SPACE_GROUPS:
        return (
            None,
            None,
            (
                f"hkl labelling in {chosen or 'the space group of a CIF structure'} "
                "arrives in the next release; plotted without hkl labels"
            ),
        )
    return [cell["a"], cell["c"]], chosen, None


def _run_plot(args: argparse.Namespace) -> int:
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

    cell, space_group = args.cell, args.space_group or DEFAULT_SPACE_GROUP
    if cell is None and item.sample is not None:
        cell, space_group, note = _structure_cell(item, args.space_group)
        if note is not None:
            print(f"xrdkit plot: {note}", file=sys.stderr)
    if cell is not None:
        a, c = cell
        try:
            indexed, fit = index_and_refine(
                peaks,
                start_cell=Cell.tetragonal(a, c),
                wavelength=scan.wavelength,
                zero_offset=args.zero,
                space_group=space_group,
            )
        except ValueError as error:
            raise CommandError(f"indexing failed: {error}") from error
        written.append(indexed_to_csv(indexed, results / f"indexed_{stem}.csv"))

        fig, _ = _hkl_figure(scan, indexed, args.scale, label)
        written += save_figure(fig, figures / f"pattern_hkl_{stem}")

        summary = indexing_summary(indexed)
        values = {
            "a": fit.cell.a,
            "c": fit.cell.c,
            "zero": fit.zero_offset,
            "n_indexed": summary["n_indexed"],
        }
        if not args.json:
            print(
                f"a = {fit.cell.a:.4f}, c = {fit.cell.c:.4f} angstrom, "
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
            "Plot one .xrdml scan and write its peak list. With --cell the peaks "
            "are indexed from that start cell, the indexing is written out, and "
            "a second figure is labelled with hkl. Without --cell no indexing "
            "is done: hkl labels need a start cell; tetragonal only in this "
            "version. A sample of the project file is read at its instrument's "
            "wavelength and labelled with its key and composition, and without "
            "--cell is indexed from the cell and space group of its first "
            "structure; a structure that is not tetragonal P4bm is plotted "
            "without hkl labels, with a note, until the next release."
        ),
    )
    parser.add_argument(
        "scan", metavar="SCAN", help="path to a .xrdml file, or a sample key"
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
    parser.add_argument(
        "--cell",
        nargs=2,
        type=float,
        metavar=("A", "C"),
        help=(
            "tetragonal start cell in angstroms, to index from; hkl labels need "
            "a start cell; tetragonal only in this version"
        ),
    )
    parser.add_argument(
        "--space-group",
        metavar="SG",
        help=(
            "space group whose reflection conditions apply (default: "
            f"{DEFAULT_SPACE_GROUP}, or a sample's structure's)"
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
            "Plot several .xrdml scans stacked vertically, drawn bottom to top in "
            "the order given, each with its label. Files and sample keys of the "
            "project file may be mixed; with a sample among them and no --out "
            "the figure goes to results/stack/STEM under the project root."
        ),
    )
    parser.add_argument(
        "scans", metavar="SCAN", nargs="+", help="paths to .xrdml files, or sample keys"
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


def _volume(cell: dict[str, float]) -> float:
    """The volume of a cell of any symmetry from its six parameters."""
    ca, cb, cg = (
        math.cos(math.radians(cell[angle])) for angle in ("alpha", "beta", "gamma")
    )
    root = math.sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg)
    return cell["a"] * cell["b"] * cell["c"] * root


def _density_from_project(args: argparse.Namespace) -> _Input:
    """Fill the options not given from the sample ``args.sample`` and its
    first structure: the formula, z, the cell (as ``--cell`` for a tetragonal
    a and c, as ``--volume`` for any other) and the Archimedes density."""
    item = _sample({}, args.sample)
    spec = item.structure
    try:
        if args.formula is None:
            args.formula = spec.composition
        if args.z is None:
            args.z = resolved_z(spec)
        if args.cell is None and args.volume is None:
            cell = resolved_cell(spec)
            angles = (cell["alpha"], cell["beta"], cell["gamma"])
            if (
                cell["a"] == cell["b"]
                and angles == (90.0, 90.0, 90.0)
                and (spec.cell is not None and set(spec.cell) == {"a", "c"})
            ):
                args.cell = [cell["a"], cell["c"]]
            else:
                args.volume = _volume(cell)
    except ValueError as error:
        raise CommandError(f"{error}; or give it on the command line") from None
    if args.archimedes is None and item.sample.archimedes is not None:
        args.archimedes = [item.sample.archimedes]
    return item


def _run_density(args: argparse.Namespace) -> int:
    item = None if args.sample is None else _density_from_project(args)

    # Checked here rather than by argparse, so that a missing or conflicting
    # option returns 1 with one line, as every other failure of a command does.
    if args.formula is None:
        raise CommandError("--formula is required")
    if args.z is None:
        raise CommandError("--z is required")
    if (args.cell is None) == (args.volume is None):
        raise CommandError("give one of --cell and --volume, not both or neither")
    if args.esd_cell is not None and args.cell is None:
        raise CommandError("--esd-cell needs --cell")
    if args.esd_volume is not None and args.volume is None:
        raise CommandError("--esd-volume needs --volume")
    if args.archimedes is not None and len(args.archimedes) > 2:
        raise CommandError("--archimedes takes a density and at most one esd")

    a = c = esd_a = esd_c = None
    measured = esd_measured = relative = esd_relative = None
    try:
        mass = formula_mass(args.formula)
        if args.cell is not None:
            a, c = args.cell
            if args.esd_cell is not None:
                esd_a, esd_c = args.esd_cell
            volume, esd_volume = cell_volume(
                Cell.tetragonal(a, c),
                None if args.esd_cell is None else {"a": esd_a, "c": esd_c},
            )
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
        "a_angstrom": a,
        "esd_a_angstrom": esd_a,
        "c_angstrom": c,
        "esd_c_angstrom": esd_c,
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
            "and with --archimedes the relative density. Give the cell as a and "
            "c of a tetragonal cell with --cell, or as a volume of any symmetry "
            "with --volume. One row with every input and result, the method and "
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
    parser.add_argument(
        "--cell",
        nargs=2,
        type=float,
        metavar=("A", "C"),
        help="tetragonal cell in angstroms; tetragonal only in this version",
    )
    parser.add_argument(
        "--esd-cell",
        nargs=2,
        type=float,
        metavar=("EA", "EC"),
        help="esds of a and c, in angstroms",
    )
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
