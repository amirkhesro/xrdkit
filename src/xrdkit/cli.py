"""The ``xrdkit`` command.

Each command is a subparser with a ``handler`` default that takes the parsed
arguments and returns the exit status. A new command adds an ``_add_<name>``
function and a call to it in :func:`build_parser`. A handler that cannot go on
raises :class:`CommandError`, which :func:`main` prints as one line, prefixed
with the command name, and turns into exit status 1.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from xrdkit import __version__
from xrdkit.density import (
    cell_volume,
    formula_mass,
    relative_density,
    theoretical_density,
)
from xrdkit.indexing import (
    DEFAULT_SPACE_GROUP,
    DEFAULT_ZERO_OFFSET,
    IndexedPeak,
    TetragonalCell,
    index_and_refine,
    indexed_to_csv,
    indexing_summary,
)
from xrdkit.io import XRDScan, read_xrdml
from xrdkit.peaks import exclude_kalpha2, find_peaks, peaks_to_csv
from xrdkit.plotting import (
    annotate_hkl,
    apply_style,
    plot_pattern,
    plot_stacked,
    save_figure,
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


def _finish(args: argparse.Namespace, written: list[Path], **values) -> int:
    """Print the files written, one per line, or them and ``values`` as JSON."""
    if args.json:
        print(json.dumps(_jsonable({"files": written, **values}), indent=2))
    else:
        for path in written:
            print(path)
    return 0


def _run_check(args: argparse.Namespace) -> int:
    (scan,) = _read_scans([args.scan])
    quality = assess_scan(scan)
    if args.json:
        print(json.dumps(_jsonable(asdict(quality)), indent=2))
    else:
        print(format_report(quality))
    return 0


def _add_check(subparsers) -> None:
    parser = subparsers.add_parser(
        "check",
        help="report the data quality of a scan, with a verdict per workflow",
        description=(
            "Read a .xrdml scan, print its range, step, counting time and "
            "intensities, and say whether it is good enough for plotting, phase "
            "identification, a Le Bail fit and a Rietveld refinement, by the "
            "criteria of the user guide, Section 2."
        ),
    )
    parser.add_argument("scan", metavar="SCAN", help="path to a .xrdml file")
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
        default=".",
        help=(
            "output root; results/ and figures/ are created under it as needed "
            "(default: the current folder)"
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


def _run_plot(args: argparse.Namespace) -> int:
    (scan,) = _read_scans([args.scan])
    stem = args.stem or Path(args.scan).stem
    out = Path(args.out)
    written: list[Path] = []
    values: dict = {}

    apply_style()
    peaks = find_peaks(scan)
    if args.no_satellites:
        peaks = exclude_kalpha2(peaks)
    written.append(peaks_to_csv(peaks, out / "results" / f"peaks_{stem}.csv"))

    label = args.label if args.label is not None else scan.sample_id
    fig, ax = plot_pattern(scan, scale=args.scale, label=label)
    ax.legend(frameon=False)
    written += save_figure(fig, out / "figures" / f"pattern_{stem}")

    if args.cell is not None:
        a, c = args.cell
        try:
            indexed, fit = index_and_refine(
                peaks,
                start_cell=TetragonalCell(a=a, c=c),
                wavelength=scan.wavelength,
                zero_offset=args.zero,
                space_group=args.space_group,
            )
        except ValueError as error:
            raise CommandError(f"indexing failed: {error}") from error
        written.append(indexed_to_csv(indexed, out / "results" / f"indexed_{stem}.csv"))

        fig, _ = _hkl_figure(scan, indexed, args.scale, label)
        written += save_figure(fig, out / "figures" / f"pattern_hkl_{stem}")

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
            "version."
        ),
    )
    parser.add_argument("scan", metavar="SCAN", help="path to a .xrdml file")
    parser.add_argument(
        "--label",
        metavar="TEXT",
        help="legend label for the trace (default: the scan's sample id)",
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
        default=DEFAULT_SPACE_GROUP,
        metavar="SG",
        help=f"space group whose reflection conditions apply (default: {DEFAULT_SPACE_GROUP})",
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
    scans = _read_scans(args.scans)
    if len(args.labels) != len(scans):
        raise CommandError(
            f"{len(args.labels)} labels for {len(scans)} scans; give one label per scan"
        )
    stem = args.stem or "_".join(Path(path).stem for path in args.scans)
    out = Path(args.out)

    apply_style()
    fig, _, _, _ = plot_stacked(
        scans,
        labels=args.labels,
        offset=args.offset,
        scale=args.scale,
        normalise=args.normalise,
    )
    written = save_figure(fig, out / "figures" / f"stack_{stem}")
    return _finish(args, written)


def _add_stack(subparsers) -> None:
    parser = subparsers.add_parser(
        "stack",
        help="plot several scans stacked one above the other",
        description=(
            "Plot several .xrdml scans stacked vertically, drawn bottom to top in "
            "the order given, each with its label."
        ),
    )
    parser.add_argument(
        "scans", metavar="SCAN", nargs="+", help="paths to .xrdml files"
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


def _run_density(args: argparse.Namespace) -> int:
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
            volume, esd_volume = cell_volume(TetragonalCell(a=a, c=c), esd_a, esd_c)
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
    path = Path(args.out) / "results" / f"density_{args.stem or 'density'}.csv"
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
            "the date is written to results/density_STEM.csv."
        ),
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
    _add_output_options(parser, "density")
    parser.set_defaults(handler=_run_density)


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the ``xrdkit`` command."""
    parser = argparse.ArgumentParser(
        prog="xrdkit", description="X-ray diffraction analysis toolkit."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
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
