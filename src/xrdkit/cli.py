"""The ``xrdkit`` command.

Each command is a subparser with a ``handler`` default that takes the parsed
arguments and returns the exit status. A new command adds an ``_add_<name>``
function and a call to it in :func:`build_parser`.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path

import numpy as np

from xrdkit import __version__
from xrdkit.io import read_xrdml
from xrdkit.quality import assess_scan, format_report

__all__ = ["build_parser", "main"]


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
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _error(command: str, message: str) -> int:
    print(f"xrdkit {command}: {message}", file=sys.stderr)
    return 1


def _run_check(args: argparse.Namespace) -> int:
    path = Path(args.scan)
    if not path.is_file():
        return _error("check", f"no such file: {path}")
    try:
        scan = read_xrdml(path)
    except (ET.ParseError, ValueError) as error:
        return _error("check", f"cannot read {path}: {error}")

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
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ``xrdkit`` command with ``argv``, or the process arguments."""
    args = build_parser().parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
