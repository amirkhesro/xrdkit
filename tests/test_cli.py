"""Tests for xrdkit.cli."""

import copy
import csv
import datetime
import json
import os
import sys
import warnings
from pathlib import Path
from urllib.error import URLError

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg

import xrdkit
from xrdkit import (
    Cell,
    CodRecord,
    Gsas2Install,
    SimulatedReflection,
    TetragonalCell,
    cli,
    find_gsas2,
    find_peaks,
    generate_reflections,
    index_and_refine,
    pipeline,
    read_xrdml,
)
from xrdkit.cli import HKL_HEADROOM, main
from xrdkit.indexing import DEFAULT_TOLERANCE
from xrdkit.instrument import REFINED_KEYS
from xrdkit.phases import MissingPhasesExtra
from xrdkit.project import PROJECT_FILE, Sample, load_project
from xrdkit.quality import WORKFLOWS

# A minimal PANalytical XRDMeasurement 2.x scan, with the byte order mark the
# instrument software writes.
XRDML = """﻿<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">
  <sample><id>synthetic</id></sample>
  <xrdMeasurement>
    <usedWavelength>
      <kAlpha1>1.5405980</kAlpha1>
      <kAlpha2>1.5444260</kAlpha2>
      <ratioKAlpha2KAlpha1>0.5</ratioKAlpha2KAlpha1>
    </usedWavelength>
    <scan>
      <dataPoints>
        <positions axis="2Theta" unit="deg">
          <startPosition>{start}</startPosition>
          <endPosition>{end}</endPosition>
        </positions>
        <commonCountingTime>1.5</commonCountingTime>
        <counts>{counts}</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
"""


def _write_xrdml(path: Path) -> Path:
    """A 10 to 100 degree scan at 0.02 degrees: background 300, one peak 9000."""
    counts = np.full(4501, 300, dtype=int)
    counts[1000] = 9000
    text = XRDML.format(start=10.0, end=100.0, counts=" ".join(map(str, counts)))
    path.write_text(text, encoding="utf-8")
    return path


# The pattern test_indexing indexes from TTB_CELL: the reflections of this cell
# that are at least 0.30 degrees from their neighbours, every one moved by a
# 0.17 degree zero offset. Here they are drawn as Gaussians into a scan.
REFINED_CELL = TetragonalCell(a=12.58, c=3.96)
APPLIED_ZERO_OFFSET = 0.17
ISOLATION = 0.30


def _write_indexable_xrdml(
    path: Path, cell: Cell = REFINED_CELL, space_group: str | None = "P4bm"
) -> tuple[Path, int]:
    """Write the shifted pattern of ``cell`` as a scan, returning it and its
    peak count."""
    reflections = generate_reflections(cell, 1.540598, 80.0, space_group=space_group)
    positions = [reflection.two_theta for reflection in reflections]
    isolated = [
        position
        for index, position in enumerate(positions)
        if position > 10.5
        and (index == 0 or position - positions[index - 1] > ISOLATION)
        and (index == len(positions) - 1 or positions[index + 1] - position > ISOLATION)
    ]
    two_theta = np.linspace(10.0, 80.0, 7001)
    heights = np.random.default_rng(11).uniform(1000.0, 5000.0, len(isolated))
    sigma = 0.08 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    counts = np.full(two_theta.size, 100.0)
    for position, height in zip(isolated, heights):
        centre = position + APPLIED_ZERO_OFFSET
        counts += height * np.exp(-0.5 * ((two_theta - centre) / sigma) ** 2)
    text = XRDML.format(
        start=10.0, end=80.0, counts=" ".join(str(round(v)) for v in counts)
    )
    path.write_text(text, encoding="utf-8")
    return path, len(isolated)


# A cubic standard for xrdkit instrument: Pm-3m near LaB6, drawn as K alpha
# doublets on a flat background, wide enough for the width fit to have peaks
# across the window.
STANDARD_A = 4.15683
STANDARD_CIF = """data_LaB6
_cell_length_a  4.15683
_cell_length_b  4.15683
_cell_length_c  4.15683
_cell_angle_alpha  90
_cell_angle_beta   90
_cell_angle_gamma  90
_symmetry_space_group_name_H-M  "P m -3 m"
loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_occupancy
   _atom_site_U_iso_or_equiv
La  La  0.0  0.0  0.0     1.0  0.0086
B   B   0.5  0.5  0.2021  1.0  0.0090
"""


def _write_standard_xrdml(path: Path) -> Path:
    """The standard drawn from a Caglioti relation, as an .xrdml."""
    step = 0.013
    two_theta = np.arange(8.0, 100.0, step)
    counts = np.full(two_theta.size, 200.0)
    reflections = generate_reflections(
        Cell.cubic(STANDARD_A), 1.540598, 99.0, space_group="Pm-3m"
    )
    for reflection in reflections:
        centre = reflection.two_theta
        if centre < 10.0:
            continue
        tan_theta = np.tan(np.radians(centre / 2.0))
        fwhm = np.sqrt(0.0060 * tan_theta**2 - 0.0035 * tan_theta + 0.0045)
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        height = 6000.0 * reflection.multiplicity / 8.0
        counts += height * np.exp(-0.5 * ((two_theta - centre) / sigma) ** 2)
        theta = np.radians(centre / 2.0)
        satellite = 2.0 * np.degrees(np.arcsin(np.sin(theta) * 1.544426 / 1.540598))
        counts += 0.5 * height * np.exp(-0.5 * ((two_theta - satellite) / sigma) ** 2)
    text = XRDML.format(
        start=two_theta[0],
        end=two_theta[-1],
        counts=" ".join(str(round(v)) for v in counts),
    )
    path.write_text(text, encoding="utf-8")
    return path


def test_check_prints_the_report(tmp_path, capsys) -> None:
    path = _write_xrdml(tmp_path / "scan.xrdml")

    assert main(["check", str(path)]) == 0

    out = capsys.readouterr().out
    lines = out.splitlines()
    assert lines[0] == "range   10.00 to 100.00 degrees"
    assert "points  4501" in lines
    assert "maximum 9000 counts" in lines
    assert "median  300 counts" in lines
    assert "plotting: suitable" in lines
    assert "phase_identification: suitable" in lines
    le_bail = [line for line in lines if line.startswith("le_bail: ")]
    assert le_bail == [
        (
            "le_bail: not suitable (range 10.00 to 100.00 degrees, short of 10 to 120; "
            "strongest peak 9000 counts, below 10000; "
            "high angle peak over median 1.0, below 10)"
        )
    ]
    assert any(line.startswith("rietveld: not suitable (") for line in lines)


def test_check_json(tmp_path, capsys) -> None:
    path = _write_xrdml(tmp_path / "scan.xrdml")

    assert main(["check", str(path), "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result["points"] == 4501
    assert result["time_per_step"] == 1.5
    assert result["maximum"] == 9000.0
    assert result["peak_over_median"] == 30.0
    assert list(result["verdicts"]) == list(WORKFLOWS)
    assert result["verdicts"]["plotting"] == {"suitable": True, "reasons": []}
    assert result["verdicts"]["le_bail"]["suitable"] is False


def test_check_missing_file(tmp_path, capsys) -> None:
    path = tmp_path / "absent.xrdml"

    assert main(["check", str(path)]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == f"xrdkit check: no such file: {path}"


def test_check_unreadable_file(tmp_path, capsys) -> None:
    path = tmp_path / "broken.xrdml"
    path.write_text("not xml", encoding="utf-8")

    assert main(["check", str(path)]) == 1
    assert capsys.readouterr().err.startswith(f"xrdkit check: cannot read {path}")


def test_a_command_is_required(capsys) -> None:
    with pytest.raises(SystemExit) as raised:
        main([])
    assert raised.value.code == 2
    assert "COMMAND" in capsys.readouterr().err


def test_plot_writes_the_peaks_and_the_figure(tmp_path, capsys) -> None:
    scan = _write_xrdml(tmp_path / "scan.xrdml")
    out = tmp_path / "out"

    assert main(["plot", str(scan), "--out", str(out)]) == 0

    written = [
        out / "results" / "peaks_scan.csv",
        out / "figures" / "pattern_scan.png",
        out / "figures" / "pattern_scan.pdf",
    ]
    for path in written:
        assert path.is_file(), path
    assert not (out / "results" / "indexed_scan.csv").exists()
    assert not (out / "figures" / "pattern_hkl_scan.png").exists()
    assert capsys.readouterr().out.splitlines() == [str(path) for path in written]


def test_plot_with_a_cell_indexes_and_labels(tmp_path, capsys) -> None:
    scan, n_peaks = _write_indexable_xrdml(tmp_path / "ttb.xrdml")
    out = tmp_path / "out"

    argv = ["plot", str(scan), "--cell", "12.45", "3.94", "--out", str(out)]
    assert main([*argv, "--stem", "x10", "--scale", "sqrt", "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result["files"] == [
        str(out / "results" / "peaks_x10.csv"),
        str(out / "figures" / "pattern_x10.png"),
        str(out / "figures" / "pattern_x10.pdf"),
        str(out / "results" / "indexed_x10.csv"),
        str(out / "figures" / "pattern_hkl_x10.png"),
        str(out / "figures" / "pattern_hkl_x10.pdf"),
    ]
    for path in result["files"]:
        assert Path(path).is_file(), path
    assert result["a"] == pytest.approx(REFINED_CELL.a, abs=0.01)
    assert result["c"] == pytest.approx(REFINED_CELL.c, abs=0.01)
    assert result["zero"] == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.03)
    assert result["n_indexed"] == n_peaks

    # Without --json the cell goes on one line before the files.
    assert main([*argv, "--stem", "x10"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("a = 12.58")
    assert f"{n_peaks} of {n_peaks} peaks indexed" in lines[0]
    assert lines[1:] == result["files"]


def test_plot_missing_file(tmp_path, capsys) -> None:
    path = tmp_path / "absent.xrdml"

    assert main(["plot", str(path), "--out", str(tmp_path)]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == f"xrdkit plot: no such file: {path}"
    assert not (tmp_path / "figures").exists()


def test_stack_writes_the_figure(tmp_path, capsys) -> None:
    first = _write_xrdml(tmp_path / "10s.xrdml")
    second = _write_xrdml(tmp_path / "12s.xrdml")
    out = tmp_path / "out"

    argv = ["stack", str(first), str(second), "--labels", "x = 0.10", "x = 0.12"]
    assert main([*argv, "--out", str(out), "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result == {
        "files": [
            str(out / "figures" / "stack_10s_12s.png"),
            str(out / "figures" / "stack_10s_12s.pdf"),
        ]
    }
    for path in result["files"]:
        assert Path(path).is_file(), path

    assert main([*argv, "--out", str(out), "--offset", "2", "--no-normalise"]) == 0
    assert capsys.readouterr().out.splitlines() == result["files"]


def test_stack_needs_one_label_per_scan(tmp_path, capsys) -> None:
    first = _write_xrdml(tmp_path / "a.xrdml")
    second = _write_xrdml(tmp_path / "b.xrdml")

    argv = ["stack", str(first), str(second), "--labels", "only one"]
    assert main([*argv, "--out", str(tmp_path)]) == 1

    assert capsys.readouterr().err.strip() == (
        "xrdkit stack: 1 labels for 2 scans; give one label per scan"
    )
    assert not (tmp_path / "figures").exists()


def test_stack_missing_file(tmp_path, capsys) -> None:
    first = _write_xrdml(tmp_path / "a.xrdml")
    missing = tmp_path / "b.xrdml"

    argv = ["stack", str(first), str(missing), "--labels", "a", "b"]
    assert main(argv) == 1
    assert capsys.readouterr().err.strip() == f"xrdkit stack: no such file: {missing}"


def test_hkl_labels_stay_inside_the_axes(tmp_path, monkeypatch) -> None:
    path, _ = _write_indexable_xrdml(tmp_path / "ttb.xrdml")
    scan = read_xrdml(path)
    peaks = find_peaks(scan)
    indexed, _ = index_and_refine(
        peaks, TetragonalCell(a=12.45, c=3.94), scan.wavelength, space_group="P4bm"
    )

    def label_tops(headroom: float) -> tuple[float, float]:
        monkeypatch.setattr(cli, "HKL_HEADROOM", headroom)
        fig, ax = cli._hkl_figure(scan, indexed, "linear", "TTB")
        renderer = FigureCanvasAgg(fig).get_renderer()
        highest = max(text.get_window_extent(renderer).y1 for text in ax.texts)
        return highest, ax.get_window_extent(renderer).y1

    assert len(indexed) > 0
    # Without headroom the label on the tallest peak runs past the top.
    highest, top = label_tops(0.0)
    assert highest > top
    highest, top = label_tops(HKL_HEADROOM)
    assert highest <= top


DENSITY = ["density", "--formula", "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6", "--z", "5"]


def test_density_from_a_cell(tmp_path, capsys) -> None:
    argv = [*DENSITY, "--cell", "12.45", "3.94", "--esd-cell", "0.001", "0.0005"]
    argv += ["--archimedes", "5.15", "0.02", "--out", str(tmp_path), "--stem", "x10"]

    assert main(argv) == 0

    path = tmp_path / "results" / "density_x10.csv"
    lines = capsys.readouterr().out.splitlines()
    assert lines == [
        "M = 394.905 g/mol per formula unit",
        "tetragonal cell a = 12.4500 +/- 0.0010, c = 3.9400 +/- 0.0005 angstrom",
        "V = 610.710 +/- 0.125 cubic angstrom",
        "theoretical density 5.3688 +/- 0.0011 g/cm3",
        "Archimedes density 5.1500 +/- 0.0200 g/cm3",
        "relative density 95.92 +/- 0.37 per cent",
        str(path),
    ]
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    row = rows[0]
    assert row["formula"] == "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
    assert float(row["z"]) == 5
    assert float(row["a_angstrom"]) == 12.45
    assert float(row["esd_c_angstrom"]) == 0.0005
    assert row["crystal_system"] == "tetragonal"
    assert float(row["b_angstrom"]) == 12.45
    assert row["esd_b_angstrom"] == ""
    assert float(row["gamma_deg"]) == 90.0
    assert float(row["theoretical_density_g_cm3"]) == pytest.approx(5.3688, abs=1e-4)
    assert float(row["relative_density_percent"]) == pytest.approx(95.92, abs=0.01)
    assert row["method"] == (
        "theoretical density from cell volume and formula mass; "
        "relative density = measured / theoretical"
    )
    today = datetime.datetime.now(datetime.UTC).astimezone().date()
    assert row["date"] == today.isoformat()


def test_density_from_a_volume_json(tmp_path, capsys) -> None:
    argv = [*DENSITY, "--volume", "610.0", "--esd-volume", "0.3", "--json"]

    assert main([*argv, "--out", str(tmp_path)]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result["files"] == [str(tmp_path / "results" / "density.csv")]
    assert Path(result["files"][0]).is_file()
    assert result["volume_a3"] == 610.0
    assert result["a_angstrom"] is None
    assert result["archimedes_density_g_cm3"] is None
    assert result["relative_density_percent"] is None
    assert result["esd_theoretical_density_g_cm3"] == pytest.approx(
        result["theoretical_density_g_cm3"] * 0.3 / 610.0
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--formula", "LaB6", "--volume", "71.8"], "--z is required"),
        (["--z", "1", "--volume", "71.8"], "--formula is required"),
        (DENSITY[1:] + ["--volume", "610", "--cell", "12.45", "3.94"], "one of --cell"),
        (DENSITY[1:], "one of --cell"),
        (DENSITY[1:] + ["--volume", "610", "--esd-cell", "0.1", "0.1"], "--esd-cell"),
        (DENSITY[1:] + ["--volume", "610", "--archimedes", "5", "0.1", "1"], "at most"),
        (["--formula", "La(B6", "--z", "1", "--volume", "71.8"], "unmatched '('"),
    ],
    ids=["no z", "no formula", "both", "neither", "stray esd", "archimedes", "formula"],
)
def test_density_rejects_bad_options(tmp_path, capsys, argv, message) -> None:
    assert main(["density", *argv, "--out", str(tmp_path)]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("xrdkit density: ")
    assert message in captured.err
    assert not (tmp_path / "results").exists()


# The project file: add-sample, and the commands given a sample key.

COMPOSITION = "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"

# Appended to the file xrdkit init writes.
PROJECT_TABLES = f"""
[instruments.lab_diffractometer]
wavelength = [1.540598, 1.544426]
ka2 = true

[structures.ttb_x010]
library = "ttb/P4bm"
composition = "{COMPOSITION}"
cell = {{ a = 12.45, c = 3.94 }}

[structures.bto]
library = "perovskite/P4mm"
composition = "BaTiO3"
cell = {{ a = 3.99, c = 4.03 }}
"""

BFO_STRUCTURE = """
[structures.bfo]
library = "perovskite/R3c"
composition = "BiFeO3"
cell = { a = 5.5876, c = 13.867 }
"""

MO_INSTRUMENT = """
[instruments.mo]
wavelength = [0.709300, 0.713590]
ka2 = true
"""


@pytest.fixture
def project(tmp_path, monkeypatch, capsys) -> Path:
    """A project made by xrdkit init, with an instrument, two structures and
    a scan under data/raw, the current folder its root."""
    root = (tmp_path / "project").resolve()
    root.mkdir()
    monkeypatch.chdir(root)
    assert main(["init", "--name", "demo"]) == 0
    with (root / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(PROJECT_TABLES)
    _write_xrdml(root / "data" / "raw" / "10c.xrdml")
    capsys.readouterr()
    return root


@pytest.fixture
def sample(project, capsys) -> Path:
    """The project with the scan added as sample 10c."""
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", "ttb_x010"]
    assert main([*argv, "--archimedes", "5.15"]) == 0
    capsys.readouterr()
    return project


def test_add_sample_appends_the_table(project, capsys) -> None:
    path = project / PROJECT_FILE
    before = path.read_bytes()
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", "ttb_x010"]
    argv += ["--structure", "bto", "--stage", "calcined", "--temperature-c", "1200"]
    argv += ["--archimedes", "5.15", "--notes", 'first "batch"']

    assert main(argv) == 0

    table = [
        "[samples.10c]",
        'file = "data/raw/10c.xrdml"',
        'instrument = "lab_diffractometer"',
        'structures = ["ttb_x010", "bto"]',
        'stage = "calcined"',
        'form = "powder"',
        "temperature_c = 1200",
        "archimedes = 5.15",
        'notes = "first \\"batch\\""',
    ]
    assert capsys.readouterr().out.splitlines() == table
    newline = "\r\n" if b"\r\n" in before else "\n"
    after = path.read_bytes()
    assert after[: len(before)] == before
    assert after[len(before) :] == (newline + newline.join(table) + newline).encode()
    assert load_project(project).samples["10c"] == Sample(
        key="10c",
        file=project / "data" / "raw" / "10c.xrdml",
        instrument="lab_diffractometer",
        structures=("ttb_x010", "bto"),
        form="powder",
        stage="calcined",
        temperature_c=1200.0,
        archimedes=5.15,
        notes='first "batch"',
    )


def test_add_sample_name_and_form(project, capsys) -> None:
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", "ttb_x010"]

    assert main([*argv, "--name", "x0.10.pellet", "--form", "pellet"]) == 0

    assert capsys.readouterr().out.splitlines()[0] == '[samples."x0.10.pellet"]'
    assert load_project(project).samples["x0.10.pellet"].form == "pellet"


def refused(project: Path, argv: list[str], message: str, capsys) -> None:
    """Run add-sample with ``argv``, expecting it to fail with ``message``
    and leave the project file alone."""
    path = project / PROJECT_FILE
    before = path.read_bytes()

    assert main(["add-sample", *argv]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("xrdkit add-sample: ")
    assert message in captured.err
    assert path.read_bytes() == before


def test_add_sample_refuses_a_duplicate_key(sample, capsys) -> None:
    argv = ["data/raw/10c.xrdml", "--structure", "bto"]
    refused(sample, argv, "sample '10c' is in", capsys)


def test_add_sample_refuses_a_file_outside_the_project(
    project, tmp_path, capsys
) -> None:
    outside = _write_xrdml(tmp_path / "elsewhere.xrdml")
    argv = [str(outside), "--structure", "ttb_x010"]
    refused(project, argv, "is outside the project folder", capsys)


def test_add_sample_refuses_an_unknown_structure(project, capsys) -> None:
    argv = ["data/raw/10c.xrdml", "--structure", "ttb_x012"]
    refused(project, argv, "no structure 'ttb_x012' under structures", capsys)


def test_add_sample_needs_an_instrument_when_there_are_two(project, capsys) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(MO_INSTRUMENT)
    argv = ["data/raw/10c.xrdml", "--structure", "ttb_x010"]

    refused(project, argv, "give --instrument; the project has 2 instruments", capsys)
    assert main(["add-sample", *argv, "--instrument", "mo"]) == 0
    assert load_project(project).samples["10c"].instrument == "mo"


def test_check_by_key(sample, capsys) -> None:
    assert main(["check", "10c"]) == 0

    path = sample / "results" / "check" / "10c" / "check_10c.txt"
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == f"sample  10c, {COMPOSITION}"
    assert lines[1] == "range   10.00 to 100.00 degrees"
    assert lines[-1] == str(path)
    assert path.read_text(encoding="utf-8").splitlines() == lines[:-1]

    assert main(["check", "10c", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["sample"] == "10c"
    assert result["composition"] == COMPOSITION
    assert result["points"] == 4501
    assert result["files"] == [str(path.with_suffix(".json"))]


def test_plot_by_key_indexes_from_the_structure_cell(project, capsys) -> None:
    # A tetragonal P4bm sample: the start cell 12.45, 3.94 of its structure.
    _, n_peaks = _write_indexable_xrdml(project / "data" / "raw" / "ttb.xrdml")
    argv = ["add-sample", "data/raw/ttb.xrdml", "--structure", "ttb_x010"]
    assert main(argv) == 0
    capsys.readouterr()

    assert main(["plot", "ttb", "--json"]) == 0

    captured = capsys.readouterr()
    assert captured.err == ""
    result = json.loads(captured.out)
    folder = project / "results" / "plot" / "ttb"
    assert result["files"] == [
        str(folder / name)
        for name in (
            "peaks_ttb.csv",
            "pattern_ttb.png",
            "pattern_ttb.pdf",
            "indexed_ttb.csv",
            "pattern_hkl_ttb.png",
            "pattern_hkl_ttb.pdf",
        )
    ]
    for path in result["files"]:
        assert Path(path).is_file(), path
    assert result["a"] == pytest.approx(REFINED_CELL.a, abs=0.01)
    assert result["c"] == pytest.approx(REFINED_CELL.c, abs=0.01)
    assert result["n_indexed"] == n_peaks

    # --out wins over the project's results folder, and --cell over its cell.
    out = project / "elsewhere"
    argv = ["plot", "ttb", "--out", str(out), "--stem", "x", "--json"]
    assert main([*argv, "--cell", "12.5", "3.95"]) == 0
    assert json.loads(capsys.readouterr().out)["n_indexed"] == n_peaks
    assert (out / "results" / "indexed_x.csv").is_file()
    assert (out / "figures" / "pattern_hkl_x.png").is_file()


@pytest.mark.parametrize("structure", ["bfo", "bto"])
def test_plot_by_key_that_cannot_be_indexed_is_plotted_with_a_note(
    project, capsys, structure
) -> None:
    # The scan has one peak, which no structure's cell can be refined from.
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(BFO_STRUCTURE)
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", structure]
    assert main(argv) == 0
    capsys.readouterr()

    assert main(["plot", "10c"]) == 0

    captured = capsys.readouterr()
    assert len(captured.err.splitlines()) == 1
    assert captured.err.startswith("xrdkit plot: indexing failed: ")
    assert captured.err.endswith("; plotted without hkl labels\n")
    folder = project / "results" / "plot" / "10c"
    written = [folder / name for name in ("peaks_10c.csv", "pattern_10c.png")]
    written.append(folder / "pattern_10c.pdf")
    assert captured.out.splitlines() == [str(path) for path in written]
    assert not (folder / "indexed_10c.csv").exists()


def test_plot_by_key_uses_the_instrument_wavelength(project, capsys) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(MO_INSTRUMENT + BFO_STRUCTURE)
    # A trigonal structure, which the plot command does not index in this version.
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", "bfo"]
    assert main([*argv, "--instrument", "mo"]) == 0

    def d_spacing(*options: str) -> tuple[float, float]:
        assert main(["plot", "10c", "--stem", "d", *options]) == 0
        path = project / "results" / "plot" / "10c" / "peaks_d.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            (peak,) = list(csv.DictReader(handle))
        theta = np.radians(float(peak["two_theta"]) / 2.0)
        return float(peak["d_spacing"]), 2.0 * np.sin(theta)

    # The scan file says 1.540598; the sample's instrument is Mo.
    d, two_sin_theta = d_spacing()
    assert d == pytest.approx(0.709300 / two_sin_theta, abs=1e-3)
    d, two_sin_theta = d_spacing("--wavelength", "1.0")
    assert d == pytest.approx(1.0 / two_sin_theta, abs=1e-3)


def test_density_by_key(sample, capsys) -> None:
    assert main(["density", "10c"]) == 0

    path = sample / "results" / "density" / "10c" / "density_10c.csv"
    assert capsys.readouterr().out.splitlines() == [
        "M = 394.905 g/mol per formula unit",
        "tetragonal cell a = 12.4500, c = 3.9400 angstrom",
        "V = 610.710 cubic angstrom",
        "theoretical density 5.3688 g/cm3",
        "Archimedes density 5.1500 g/cm3",
        "relative density 95.92 per cent",
        str(path),
    ]
    with path.open(newline="", encoding="utf-8") as handle:
        (row,) = list(csv.DictReader(handle))
    assert row["sample"] == "10c"
    assert row["structure"] == "ttb_x010"
    assert row["formula"] == COMPOSITION
    assert float(row["z"]) == 5
    assert float(row["a_angstrom"]) == 12.45

    # Options given on the command line win.
    assert main(["density", "10c", "--z", "4", "--volume", "600", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["z"] == 4
    assert result["volume_a3"] == 600
    assert result["a_angstrom"] is None


def test_density_by_key_of_a_cubic_or_hexagonal_cell(project, capsys) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(
            '\n[structures.bfo]\nlibrary = "perovskite/R3c"\ncomposition = "BiFeO3"\n'
            "cell = { a = 5.5876, c = 13.867 }\n"
        )
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", "bfo"]
    assert main([*argv, "--name", "bfo"]) == 0
    capsys.readouterr()

    assert main(["density", "bfo", "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    expected = 5.5876**2 * 13.867 * np.sqrt(3.0) / 2.0
    assert result["volume_a3"] == pytest.approx(expected)
    assert result["z"] == 6


def test_stack_mixes_keys_and_paths(sample, tmp_path, capsys) -> None:
    other = _write_xrdml(tmp_path / "other.xrdml")

    argv = ["stack", "10c", str(other), "--labels", "x = 0.10", "other"]
    assert main(argv) == 0

    folder = sample / "results" / "stack" / "10c_other"
    written = [folder / "stack_10c_other.png", folder / "stack_10c_other.pdf"]
    assert capsys.readouterr().out.splitlines() == [str(path) for path in written]
    for path in written:
        assert path.is_file(), path


def test_a_file_path_inside_a_project_is_used_as_before(sample, capsys) -> None:
    assert main(["plot", "data/raw/10c.xrdml"]) == 0

    assert capsys.readouterr().out.splitlines() == [
        str(Path("results") / "peaks_10c.csv"),
        str(Path("figures") / "pattern_10c.png"),
        str(Path("figures") / "pattern_10c.pdf"),
    ]
    assert not (sample / "results" / "plot").exists()

    assert main(["check", "data/raw/10c.xrdml"]) == 0
    assert capsys.readouterr().out.splitlines()[0].startswith("range ")


def test_an_unknown_key_lists_the_samples(sample, capsys) -> None:
    assert main(["check", "10d"]) == 1

    error = capsys.readouterr().err.strip()
    assert error == (
        f"xrdkit check: no such file or sample: 10d; the samples in "
        f"{sample / PROJECT_FILE} are 10c"
    )
    assert not (sample / "results" / "check").exists()


def test_a_key_outside_a_project(tmp_path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)

    assert main(["plot", "10c"]) == 1
    assert capsys.readouterr().err.strip() == (
        f"xrdkit plot: no such file: 10c, and no {PROJECT_FILE} in the current "
        "folder or above it to look it up in as a sample key"
    )


# Cells of other crystal systems for plot --cell: the numbers given, the true
# cell the scan is drawn from, and its space group.
PLOT_CELLS = {
    "cubic": (["3.905"], Cell.cubic(3.905), "Pm-3m"),
    # Not the pseudo-cubic Pbnm perovskite: below 35 degrees its reflections
    # lie too close together for the first, coarse cycle to find four lone
    # peaks, which an orthorhombic cell needs.
    "orthorhombic": (["3.8", "4.2", "5.0"], Cell.orthorhombic(3.8, 4.2, 5.0), None),
    "hexagonal": (
        ["5.148", "13.863", "--system", "hexagonal"],
        Cell.hexagonal(5.148, 13.863),
        "R3c",
    ),
    "triclinic": (
        ["5.0", "6.0", "7.0", "80", "90", "100"],
        Cell.triclinic(5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
        None,
    ),
}


@pytest.mark.parametrize("system", list(PLOT_CELLS))
def test_plot_with_a_cell_of_any_system_indexes_and_labels(
    tmp_path, capsys, system
) -> None:
    numbers, cell, space_group = PLOT_CELLS[system]
    scan, n_peaks = _write_indexable_xrdml(tmp_path / "s.xrdml", cell, space_group)
    out = tmp_path / "out"
    argv = ["plot", str(scan), "--cell", *numbers, "--out", str(out), "--json"]
    if space_group is not None:
        argv += ["--space-group", space_group]

    assert main(argv) == 0

    captured = capsys.readouterr()
    assert captured.err == ""
    result = json.loads(captured.out)
    assert result["files"][3:] == [
        str(out / "results" / "indexed_s.csv"),
        str(out / "figures" / "pattern_hkl_s.png"),
        str(out / "figures" / "pattern_hkl_s.pdf"),
    ]
    for path in result["files"]:
        assert Path(path).is_file(), path
    assert result["crystal_system"] == system
    for name, value in cell.parameters.items():
        assert result[name] == pytest.approx(value, abs=0.01), name
    assert result["zero"] == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.03)
    assert result["n_indexed"] > 0.8 * n_peaks


def test_plot_names_a_space_group_whose_conditions_are_not_known(
    tmp_path, capsys
) -> None:
    scan, _ = _write_indexable_xrdml(tmp_path / "s.xrdml", Cell.cubic(3.905), "Pm-3m")
    argv = ["plot", str(scan), "--cell", "3.905", "--space-group", "Fm-3m"]

    assert main([*argv, "--out", str(tmp_path / "out")]) == 0

    err = capsys.readouterr().err
    assert err == (
        "xrdkit plot: the reflection conditions of Fm-3m are not known, so the "
        "peaks are indexed without them and some labels may name absent "
        "reflections\n"
    )
    assert (tmp_path / "out" / "results" / "indexed_s.csv").is_file()


@pytest.mark.parametrize(
    ("library", "cell", "true_cell", "space_group"),
    [
        ("perovskite/Pm-3m", "{ a = 3.905 }", Cell.cubic(3.905), "Pm-3m"),
        (
            "perovskite/R3c",
            "{ a = 5.148, c = 13.863 }",
            Cell.trigonal(5.148, 13.863),
            "R3c",
        ),
    ],
    ids=["cubic", "hexagonal axes"],
)
def test_plot_by_key_labels_a_cubic_or_hexagonal_structure(
    project, capsys, library, cell, true_cell, space_group
) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(
            f'\n[structures.other]\nlibrary = "{library}"\n'
            f'composition = "SrTiO3"\ncell = {cell}\n'
        )
    path = project / "data" / "raw" / "other.xrdml"
    _, n_peaks = _write_indexable_xrdml(path, true_cell, space_group)
    assert main(["add-sample", "data/raw/other.xrdml", "--structure", "other"]) == 0
    capsys.readouterr()

    assert main(["plot", "other", "--json"]) == 0

    captured = capsys.readouterr()
    assert captured.err == ""
    result = json.loads(captured.out)
    folder = project / "results" / "plot" / "other"
    assert str(folder / "indexed_other.csv") in result["files"]
    assert str(folder / "pattern_hkl_other.png") in result["files"]
    assert result["crystal_system"] == true_cell.crystal_system
    assert result["a"] == pytest.approx(true_cell.a, abs=0.01)
    assert result["n_indexed"] > 0.8 * n_peaks


def test_plot_with_a_cell_it_cannot_index_fails(tmp_path, capsys) -> None:
    scan = _write_xrdml(tmp_path / "one.xrdml")

    argv = ["plot", str(scan), "--cell", "3.905", "--out", str(tmp_path / "out")]
    assert main(argv) == 1

    err = capsys.readouterr().err
    assert len(err.splitlines()) == 1
    assert err.startswith("xrdkit plot: indexing failed: ")


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (
            ["--cell", "5", "6", "7", "8", "9"],
            "--cell takes 1, 2, 3, 4 or 6 numbers, not 5",
        ),
        (
            ["--cell", "5", "6", "7", "100"],
            "--cell with 4 numbers needs --system monoclinic (a b c beta)",
        ),
        (
            ["--cell", "5", "6", "7", "--system", "cubic"],
            "--system cubic takes 1 number in --cell (a), not 3",
        ),
        (
            ["--cell", "5", "6", "--system", "monoclinic"],
            "--system monoclinic takes 4 numbers in --cell (a b c beta), not 2",
        ),
        (["--cell", "3.905", "-1"], "--cell values must be greater than 0, not -1"),
        (["--system", "cubic"], "--system needs --cell"),
    ],
    ids=["count", "missing system", "contradiction", "too few", "negative", "stray"],
)
def test_plot_rejects_bad_cell_options(tmp_path, capsys, options, message) -> None:
    scan = _write_xrdml(tmp_path / "one.xrdml")
    out = tmp_path / "out"

    assert main(["plot", str(scan), *options, "--out", str(out)]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == f"xrdkit plot: {message}\n"
    assert not out.exists()


@pytest.mark.parametrize(
    ("options", "system", "volume", "esd_volume"),
    [
        (
            ["--cell", "3.905", "--esd-cell", "0.001"],
            "cubic",
            3.905**3,
            3 * 3.905**2 * 0.001,
        ),
        (
            ["--cell", "5.38", "5.44", "7.64", "--esd-cell", "0.001", "0.002", "0.003"],
            "orthorhombic",
            5.38 * 5.44 * 7.64,
            float(
                np.sqrt(
                    (5.44 * 7.64 * 0.001) ** 2
                    + (5.38 * 7.64 * 0.002) ** 2
                    + (5.38 * 5.44 * 0.003) ** 2
                )
            ),
        ),
        (
            ["--cell", "5.148", "13.863", "--system", "hexagonal"]
            + ["--esd-cell", "0.0004", "0.002"],
            "hexagonal",
            np.sqrt(3) / 2 * 5.148**2 * 13.863,
            float(
                np.hypot(
                    np.sqrt(3) * 5.148 * 13.863 * 0.0004,
                    np.sqrt(3) / 2 * 5.148**2 * 0.002,
                )
            ),
        ),
    ],
    ids=["cubic", "orthorhombic", "hexagonal"],
)
def test_density_from_a_cell_of_any_system(
    tmp_path, capsys, options, system, volume, esd_volume
) -> None:
    argv = ["density", "--formula", "SrTiO3", "--z", "1", *options]

    assert main([*argv, "--out", str(tmp_path), "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result["volume_a3"] == pytest.approx(volume, rel=1e-9)
    assert result["esd_volume_a3"] == pytest.approx(esd_volume, rel=1e-6)
    with (tmp_path / "results" / "density.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        (row,) = list(csv.DictReader(handle))
    assert row["crystal_system"] == system
    cell = (
        Cell.from_parameters(
            system, dict(zip(cli.CELL_PARAMETERS[system], map(float, options[1:])))
        )
        if system != "hexagonal"
        else Cell.hexagonal(5.148, 13.863)
    )
    for name in ("a", "b", "c"):
        assert float(row[f"{name}_angstrom"]) == pytest.approx(getattr(cell, name))
    for name in ("alpha", "beta", "gamma"):
        assert float(row[f"{name}_deg"]) == pytest.approx(getattr(cell, name))
    assert float(row["esd_a_angstrom"]) == float(
        options[options.index("--esd-cell") + 1]
    )
    assert float(row["volume_a3"]) == pytest.approx(volume, rel=1e-9)


def test_density_prints_the_cell_it_used(tmp_path, capsys) -> None:
    argv = ["density", "--formula", "LaTiO3", "--z", "4"]
    argv += ["--cell", "5.60", "5.62", "7.91", "--out", str(tmp_path)]

    assert main(argv) == 0

    lines = capsys.readouterr().out.splitlines()
    assert lines[1] == "orthorhombic cell a = 5.6000, b = 5.6200, c = 7.9100 angstrom"


def test_density_rejects_an_esd_count_that_does_not_match(tmp_path, capsys) -> None:
    argv = ["density", "--formula", "SrTiO3", "--z", "1", "--cell", "3.905"]
    argv += ["--esd-cell", "0.001", "0.001", "--out", str(tmp_path)]

    assert main(argv) == 1

    assert capsys.readouterr().err == (
        "xrdkit density: --esd-cell takes as many numbers as --cell, 1, not 2\n"
    )


# xrdkit lattice.

LATTICE_WAVELENGTH = 1.540598
TTB_TRUE = Cell.tetragonal(12.58, 3.96)

LAB_INSTRUMENT = """
[instruments.lab]
wavelength = [1.540598, 1.544426]
ka2 = true
radius = 240
"""

PEAK_COLUMNS = [
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
]

# The columns of a lattice results row, which say what was done and how.
RESULT_COLUMNS = [
    "sample",
    "structure",
    "scan_file",
    "wavelength_angstrom",
    "form",
    "space_group",
    "crystal_system",
    "start_a_angstrom",
    "start_b_angstrom",
    "start_c_angstrom",
    "start_alpha_deg",
    "start_beta_deg",
    "start_gamma_deg",
    "a_angstrom",
    "esd_a_angstrom",
    "b_angstrom",
    "esd_b_angstrom",
    "c_angstrom",
    "esd_c_angstrom",
    "alpha_deg",
    "esd_alpha_deg",
    "beta_deg",
    "esd_beta_deg",
    "gamma_deg",
    "esd_gamma_deg",
    "volume_a3",
    "esd_volume_a3",
    "zero_deg",
    "esd_zero_deg",
    "zero_refined",
    "displacement_mm",
    "esd_displacement_mm",
    "displacement_refined",
    "radius_mm",
    "n_peaks_found",
    "n_peaks_refitted",
    "n_fits_rejected",
    "n_satellites",
    "n_peaks_recovered",
    "n_peaks_indexed",
    "n_peaks_refined",
    "rms_two_theta_deg",
    "coarse_two_theta_max_deg",
    "formula",
    "z",
    "formula_mass_g_mol",
    "theoretical_density_g_cm3",
    "esd_theoretical_density_g_cm3",
    "archimedes_density_g_cm3",
    "esd_archimedes_density_g_cm3",
    "relative_density_percent",
    "esd_relative_density_percent",
    "method",
    "date",
    "xrdkit_version",
]


def _write_doublet_xrdml(
    path: Path,
    cell: Cell,
    space_group: str | None,
    zero: float = 0.0,
    displacement: float = 0.0,
    radius: float = 240.0,
    extra: tuple[tuple[float, float], ...] = (),
) -> tuple[Path, int]:
    """A 10 to 80 degree scan of ``cell`` with every lone reflection drawn as
    a K alpha 1 and K alpha 2 doublet, moved by ``zero`` and by a specimen
    ``displacement`` in mm on a goniometer of ``radius``, and a doublet for
    each ``(position, height)`` in ``extra``. Returns the path and the number
    of lone reflections drawn."""
    reflections = generate_reflections(
        cell, LATTICE_WAVELENGTH, 79.0, space_group=space_group
    )
    positions = [reflection.two_theta for reflection in reflections]
    drawn = [
        position
        for index, position in enumerate(positions)
        if position > 11.0
        and (index == 0 or position - positions[index - 1] > 0.4)
        and (index == len(positions) - 1 or positions[index + 1] - position > 0.4)
    ]
    two_theta = np.round(np.arange(10.0, 80.0 + 1e-9, 0.01), 4)
    counts = np.full(two_theta.size, 100.0)
    heights = np.random.default_rng(3).uniform(2000.0, 8000.0, len(drawn))
    sigma = 0.10 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    ratio = 1.544426 / LATTICE_WAVELENGTH
    for position, height in [*zip(drawn, heights), *extra]:
        theta = np.radians(position / 2.0)
        centre = (
            position + zero + np.degrees(-2.0 * displacement * np.cos(theta) / radius)
        )
        satellite = 2.0 * np.degrees(
            np.arcsin(ratio * np.sin(np.radians(centre / 2.0)))
        )
        for line, share in ((centre, 1.0), (satellite, 0.5)):
            counts += share * height * np.exp(-0.5 * ((two_theta - line) / sigma) ** 2)
    text = XRDML.format(
        start=10.0, end=80.0, counts=" ".join(str(round(v)) for v in counts)
    )
    path.write_text(text, encoding="utf-8")
    return path, len(drawn)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _ttb_lattice_argv(tmp_path: Path) -> tuple[list[str], int]:
    scan, drawn = _write_doublet_xrdml(
        tmp_path / "ttb.xrdml", TTB_TRUE, "P4bm", zero=APPLIED_ZERO_OFFSET
    )
    argv = ["lattice", str(scan), "--cell", "12.45", "3.94", "--space-group", "P4bm"]
    argv += ["--formula", COMPOSITION, "--z", "5", "--archimedes", "5.15", "0.02"]
    return [*argv, "--out", str(tmp_path / "out")], drawn


def test_lattice_refines_a_tetragonal_scan_and_writes_both_files(
    tmp_path, capsys
) -> None:
    argv, drawn = _ttb_lattice_argv(tmp_path)

    assert main(argv) == 0

    folder = tmp_path / "out" / "results" / "lattice"
    peaks_path, results_path = folder / "peaks_ttb.csv", folder / "lattice.csv"
    lines = capsys.readouterr().out.splitlines()
    assert lines[-2:] == [str(peaks_path), str(results_path)]
    assert lines[0].startswith("peaks   ")
    assert f"{drawn} refitted, 0 fits rejected, " in lines[0]
    assert f"{drawn} of {drawn} indexed (100.0 per cent)" in lines[0]
    assert lines[1] == "coarse window to 35.00 degrees"
    assert lines[2].startswith("tetragonal cell a = 12.5")
    assert lines[3].startswith("V = ")
    assert lines[4].startswith("zero 0.1")
    assert lines[5] == "displacement not refined"
    assert lines[6].endswith(
        "a low indexed fraction or a high rms means the cell should not be trusted"
    )
    assert lines[7] == "M = 394.905 g/mol per formula unit, Z = 5"
    assert lines[8].startswith("theoretical density ")
    assert lines[9] == "Archimedes density 5.1500 +/- 0.0200 g/cm3"
    assert lines[10].startswith("relative density ")

    with peaks_path.open(newline="", encoding="utf-8") as handle:
        assert next(csv.reader(handle)) == PEAK_COLUMNS
    peaks = _rows(peaks_path)
    fitted = [row for row in peaks if row["fit_rejected"] == "False"]
    assert len(fitted) == drawn
    # Below 40 degrees the doublet is unresolved, and its vertex lies above K
    # alpha 1; the refit moves each of those positions down to K alpha 1.
    unresolved = [row for row in fitted if float(row["found_two_theta"]) < 40.0]
    assert unresolved
    assert all(
        float(row["fitted_two_theta"]) < float(row["found_two_theta"])
        for row in unresolved
    )
    assert all(abs(float(row["difference"])) < 0.01 for row in fitted)

    (row,) = _rows(results_path)
    assert list(row) == RESULT_COLUMNS
    assert row["crystal_system"] == "tetragonal"
    assert row["space_group"] == "P4bm"
    assert float(row["a_angstrom"]) == pytest.approx(12.58, abs=0.002)
    assert float(row["c_angstrom"]) == pytest.approx(3.96, abs=0.002)
    assert float(row["b_angstrom"]) == float(row["a_angstrom"])
    assert row["esd_b_angstrom"] == ""
    assert float(row["zero_deg"]) == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.01)
    assert row["zero_refined"] == "True" and row["displacement_refined"] == "False"
    assert float(row["start_a_angstrom"]) == 12.45
    assert float(row["coarse_two_theta_max_deg"]) == 35.0
    assert float(row["esd_volume_a3"]) > 0.0
    relative = float(row["relative_density_percent"])
    expected = 100.0 * 5.15 / float(row["theoretical_density_g_cm3"])
    assert relative == pytest.approx(expected, rel=1e-6)
    assert float(row["esd_relative_density_percent"]) > 0.0
    assert "refine_lattice" in row["method"]


def test_lattice_peaks_file_is_a_complete_indexing_table(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)

    assert main(argv) == 0

    folder = tmp_path / "out" / "results" / "lattice"
    with (folder / "peaks_ttb.csv").open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))
    assert header.index("corrected_two_theta") == (
        header.index("esd_fitted_two_theta") + 1
    )
    assert header.index("d_spacing") == header.index("corrected_two_theta") + 1
    assert header.index("relative_intensity") == header.index("intensity") + 1
    assert header[-1] == "n_candidates"
    # The results file keeps its columns, so rows written before still append.
    with (folder / "lattice.csv").open(newline="", encoding="utf-8") as handle:
        assert next(csv.reader(handle)) == RESULT_COLUMNS

    peaks = _rows(folder / "peaks_ttb.csv")
    (row,) = _rows(folder / "lattice.csv")
    zero = float(row["zero_deg"])
    for peak in peaks:
        corrected = float(peak["corrected_two_theta"])
        theta = np.radians(corrected / 2.0)
        assert float(peak["d_spacing"]) == pytest.approx(
            LATTICE_WAVELENGTH / (2.0 * np.sin(theta)), abs=1e-5
        )
        used = peak["fitted_two_theta"] or peak["found_two_theta"]
        if peak["kalpha2_satellite"] == "True":
            used = peak["found_two_theta"]
        assert corrected == pytest.approx(float(used) - zero, abs=2e-4)
    assert max(float(peak["relative_intensity"]) for peak in peaks) == 100.0
    indexed = [peak for peak in peaks if peak["h"] != ""]
    assert indexed
    assert all(int(peak["n_candidates"]) >= 1 for peak in indexed)


def test_lattice_appends_a_row_per_run(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)

    assert main(argv) == 0
    assert main(argv) == 0

    rows = _rows(tmp_path / "out" / "results" / "lattice" / "lattice.csv")
    assert len(rows) == 2
    assert rows[0]["a_angstrom"] == rows[1]["a_angstrom"]


def test_lattice_json_carries_the_numbers_of_the_row(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)

    assert main([*argv, "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    (row,) = _rows(tmp_path / "out" / "results" / "lattice" / "lattice.csv")
    assert result["files"][1] == str(
        tmp_path / "out" / "results" / "lattice" / "lattice.csv"
    )
    for key in (
        "a_angstrom",
        "c_angstrom",
        "volume_a3",
        "zero_deg",
        "rms_two_theta_deg",
    ):
        assert result[key] == pytest.approx(float(row[key]), rel=1e-12), key
    assert result["esd_b_angstrom"] is None
    assert result["n_peaks_indexed"] == int(row["n_peaks_indexed"])
    assert result["relative_density_percent"] == pytest.approx(
        float(row["relative_density_percent"])
    )


def test_lattice_by_key_of_a_cubic_sample(project, capsys) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(
            '\n[structures.sto]\nlibrary = "perovskite/Pm-3m"\n'
            'composition = "SrTiO3"\ncell = { a = 3.95 }\n'
        )
    _write_doublet_xrdml(
        project / "data" / "raw" / "sto.xrdml", Cell.cubic(3.905), "Pm-3m"
    )
    assert main(["add-sample", "data/raw/sto.xrdml", "--structure", "sto"]) == 0
    capsys.readouterr()

    assert main(["lattice", "sto"]) == 0

    folder = project / "results" / "lattice" / "sto"
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "sample  sto, SrTiO3"
    assert lines[3].startswith("cubic cell a = 3.905")
    assert lines[-2:] == [
        str(folder / "peaks_sto.csv"),
        str(folder / "lattice_sto.csv"),
    ]
    (row,) = _rows(folder / "lattice_sto.csv")
    assert row["sample"] == "sto" and row["structure"] == "sto"
    assert row["space_group"] == "Pm-3m"
    assert float(row["a_angstrom"]) == pytest.approx(3.905, abs=0.001)
    assert float(row["c_angstrom"]) == float(row["a_angstrom"])
    assert row["esd_c_angstrom"] == ""
    assert float(row["z"]) == 1
    assert float(row["theoretical_density_g_cm3"]) == pytest.approx(5.12, abs=0.02)


def test_lattice_results_give_the_rietveld_start_cell(project, capsys) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(LAB_INSTRUMENT)
    raw = project / "data" / "raw" / "ttb.xrdml"
    _write_doublet_xrdml(raw, TTB_TRUE, "P4bm", zero=APPLIED_ZERO_OFFSET)
    argv = ["add-sample", "data/raw/ttb.xrdml", "--structure", "ttb_x010"]
    assert main([*argv, "--instrument", "lab"]) == 0
    assert main(["lattice", "ttb"]) == 0
    capsys.readouterr()

    loaded = load_project(project)
    path = project / "results" / "lattice" / "ttb" / "lattice_ttb.csv"
    (row,) = _rows(path)

    cell, source = pipeline._start_cell(
        loaded,
        loaded.samples["ttb"],
        loaded.structures["ttb_x010"],
        pipeline.Options(),
    )

    # The cell is read from the columns the lattice command writes.
    assert source == f"lattice results {path}"
    assert cell == {
        "a": float(row["a_angstrom"]),
        "b": float(row["b_angstrom"]),
        "c": float(row["c_angstrom"]),
        "alpha": float(row["alpha_deg"]),
        "beta": float(row["beta_deg"]),
        "gamma": float(row["gamma_deg"]),
    }


@pytest.mark.parametrize("form", ["pellet", "powder"])
def test_lattice_refines_by_the_form_of_the_sample(project, capsys, form) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(LAB_INSTRUMENT)
    raw = project / "data" / "raw" / "ttb.xrdml"
    if form == "pellet":
        _write_doublet_xrdml(raw, TTB_TRUE, "P4bm", displacement=0.1)
    else:
        _write_doublet_xrdml(raw, TTB_TRUE, "P4bm", zero=APPLIED_ZERO_OFFSET)
    argv = ["add-sample", "data/raw/ttb.xrdml", "--structure", "ttb_x010"]
    assert main([*argv, "--instrument", "lab", "--form", form]) == 0
    capsys.readouterr()

    assert main(["lattice", "ttb"]) == 0

    lines = capsys.readouterr().out.splitlines()
    (row,) = _rows(project / "results" / "lattice" / "ttb" / "lattice_ttb.csv")
    assert row["form"] == form
    # With the zero or the displacement taken off, an indexed peak sits at the
    # Bragg position of its reflection in the refined cell.
    cell = Cell.tetragonal(float(row["a_angstrom"]), float(row["c_angstrom"]))
    peaks = _rows(project / "results" / "lattice" / "ttb" / "peaks_ttb.csv")
    indexed = [peak for peak in peaks if peak["h"] != ""]
    assert indexed
    for peak in indexed:
        d = cell.d_spacing(int(peak["h"]), int(peak["k"]), int(peak["l"]))
        bragg = 2.0 * np.degrees(np.arcsin(LATTICE_WAVELENGTH / (2.0 * d)))
        assert float(peak["corrected_two_theta"]) == pytest.approx(bragg, abs=0.01)
    assert float(row["a_angstrom"]) == pytest.approx(12.58, abs=0.003)
    assert float(row["c_angstrom"]) == pytest.approx(3.96, abs=0.003)
    if form == "pellet":
        assert "zero 0.0000 degrees (held)" in lines
        assert any(line.startswith("displacement 0.") for line in lines)
        assert row["zero_refined"] == "False" and row["esd_zero_deg"] == ""
        assert row["displacement_refined"] == "True"
        assert float(row["displacement_mm"]) == pytest.approx(0.1, abs=0.02)
        assert float(row["radius_mm"]) == 240.0
    else:
        assert "displacement not refined" in lines
        assert row["zero_refined"] == "True"
        assert float(row["zero_deg"]) == pytest.approx(APPLIED_ZERO_OFFSET, abs=0.01)
        assert row["displacement_mm"] == "" and row["displacement_refined"] == "False"


def test_lattice_displacement_without_a_radius_is_refused(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)

    assert main([*argv, "--displacement"]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        "xrdkit lattice: refining a specimen displacement needs the goniometer "
        "radius: give --radius, or radius for the instrument in the project file\n"
    )
    assert not (tmp_path / "out").exists()


def test_lattice_missing_file(tmp_path, capsys) -> None:
    path = tmp_path / "absent" / "scan.xrdml"

    assert main(["lattice", str(path), "--cell", "3.905"]) == 1

    assert capsys.readouterr().err == f"xrdkit lattice: no such file: {path}\n"


# Satellites excluded from the refinement, and recovered when a reflection
# needs them.


def _assert_counts_add_up(row: dict[str, str], peaks: list[dict[str, str]]) -> None:
    """The counts of a results row obey the stated relation and agree with the
    peaks file."""
    found, refitted, rejected, satellites, recovered, indexed, used = (
        int(row[name])
        for name in (
            "n_peaks_found",
            "n_peaks_refitted",
            "n_fits_rejected",
            "n_satellites",
            "n_peaks_recovered",
            "n_peaks_indexed",
            "n_peaks_refined",
        )
    )
    assert found == refitted + rejected + satellites
    assert recovered <= refitted
    assert used <= indexed <= refitted + rejected
    assert found == len(peaks)
    assert refitted == sum(peak["fit_rejected"] == "False" for peak in peaks)
    assert rejected == sum(peak["fit_rejected"] == "True" for peak in peaks)
    assert satellites == sum(peak["kalpha2_satellite"] == "True" for peak in peaks)
    assert recovered == sum(peak["recovered"] == "True" for peak in peaks)
    assert indexed == sum(peak["h"] != "" for peak in peaks)


def _lattice_run(argv: list[str], out: Path) -> tuple[dict, list[dict]]:
    """Run ``argv`` into ``out`` and return the results row and the peaks rows."""
    assert main([*argv, "--out", str(out)]) == 0
    folder = out / "results" / "lattice"
    (row,) = _rows(folder / "lattice.csv")
    (peaks_file,) = folder.glob("peaks_*.csv")
    return row, _rows(peaks_file)


# 2 2 2 of this cell falls exactly on the K alpha 2 position of 3 1 0, and
# every other reflection lies at least 0.6 degrees from both; the satellites of
# the lone reflections lie at least 0.15 degrees from any reflection.
COINCIDENT_CELL = Cell.tetragonal(4.05, 5.65762)


def _coincident_argv(tmp_path: Path, height: float = 1600.0) -> tuple[list[str], float]:
    """A lattice run on a scan with 3 1 0 at 8000 counts and 2 2 2 at
    ``height`` on its K alpha 2 line; at 1600 the two give a peak of 0.7 of
    3 1 0's height at the satellite position, and at 0 the peak there is 3 1
    0's satellite alone. Returns the argv and the 2 2 2 position."""
    positions = {
        reflection.hkl: reflection.two_theta
        for reflection in generate_reflections(
            COINCIDENT_CELL, LATTICE_WAVELENGTH, 79.0
        )
    }
    coincident = positions[(2, 2, 2)]
    lines = [(positions[(3, 1, 0)], 8000.0)]
    if height > 0.0:
        lines.append((coincident, height))
    scan, _ = _write_doublet_xrdml(
        tmp_path / "coincident.xrdml", COINCIDENT_CELL, None, extra=tuple(lines)
    )
    return ["lattice", str(scan), "--cell", "4.04", "5.67"], coincident


def test_lattice_excludes_a_satellite_that_matches_no_reflection(
    tmp_path, capsys
) -> None:
    argv, drawn = _ttb_lattice_argv(tmp_path)
    argv = argv[: argv.index("--out")]

    row, peaks = _lattice_run(argv, tmp_path / "run")
    line = capsys.readouterr().out.splitlines()[0]

    satellites = [peak for peak in peaks if peak["kalpha2_satellite"] == "True"]
    assert satellites
    for peak in satellites:
        assert peak["recovered"] == "False"
        assert peak["fit_rejected"] == "" and peak["fitted_two_theta"] == ""
        assert peak["h"] == "" and peak["difference"] == ""
    assert int(row["n_satellites"]) == len(satellites)
    assert int(row["n_peaks_recovered"]) == 0
    # Only the drawn K alpha 1 lines are indexed and refined.
    assert int(row["n_peaks_indexed"]) == drawn
    assert int(row["n_peaks_refined"]) == drawn
    assert f"{len(satellites)} satellites excluded, 0 recovered; " in line
    _assert_counts_add_up(row, peaks)


def test_lattice_recovers_a_satellite_that_matches_a_reflection(
    tmp_path, capsys
) -> None:
    argv, coincident = _coincident_argv(tmp_path)

    row, peaks = _lattice_run(argv, tmp_path / "run")
    line = capsys.readouterr().out.splitlines()[0]

    (recovered,) = [peak for peak in peaks if peak["recovered"] == "True"]
    assert float(recovered["found_two_theta"]) == pytest.approx(coincident, abs=0.01)
    assert recovered["kalpha2_satellite"] == "False"
    assert recovered["fit_rejected"] == "False"
    fitted = float(recovered["fitted_two_theta"])
    assert fitted != float(recovered["found_two_theta"])
    # Most of this peak is 3 1 0's K alpha 2 line, which a doublet fitted to
    # 2 2 2 alone does not model, so the refit is pulled towards 2 2 2's own
    # weaker satellite; it stays within the indexing tolerance.
    assert fitted == pytest.approx(coincident, abs=DEFAULT_TOLERANCE)
    assert (recovered["h"], recovered["k"], recovered["l"]) == ("2", "2", "2")
    assert abs(float(recovered["difference"])) <= DEFAULT_TOLERANCE
    # Every peak indexed has at least one candidate; an excluded satellite took
    # no part in the indexing and has no count.
    assert all(int(peak["n_candidates"]) >= 1 for peak in peaks if peak["h"] != "")
    # Every other flagged peak is a satellite of a lone reflection, and stays
    # out.
    others = [peak for peak in peaks if peak["kalpha2_satellite"] == "True"]
    assert others
    assert all(peak["h"] == "" for peak in others)
    assert all(peak["n_candidates"] == "" for peak in others)
    assert int(row["n_peaks_recovered"]) == 1
    assert int(row["n_satellites"]) == len(others)
    assert f"{len(others)} satellites excluded, 1 recovered; " in line
    _assert_counts_add_up(row, peaks)


def test_lattice_indexed_fraction_is_of_the_peaks_taking_part(tmp_path, capsys) -> None:
    argv, _ = _coincident_argv(tmp_path)

    row, _ = _lattice_run(argv, tmp_path / "run")
    line = capsys.readouterr().out.splitlines()[0]

    indexed = int(row["n_peaks_indexed"])
    taking_part = int(row["n_peaks_refitted"]) + int(row["n_fits_rejected"])
    percent = 100.0 * indexed / taking_part
    assert f"{indexed} of {taking_part} indexed ({percent:.1f} per cent)" in line


def test_no_satellites_leaves_an_unmatched_satellite_as_it_was(
    tmp_path, capsys
) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)
    argv = argv[: argv.index("--out")]

    row, peaks = _lattice_run(argv, tmp_path / "default")
    bare_row, bare_peaks = _lattice_run([*argv, "--no-satellites"], tmp_path / "bare")

    for name in (
        "n_peaks_found",
        "n_peaks_refitted",
        "n_fits_rejected",
        "n_satellites",
        "n_peaks_recovered",
        "n_peaks_indexed",
        "n_peaks_refined",
        "a_angstrom",
        "c_angstrom",
    ):
        assert bare_row[name] == row[name], name
    assert bare_peaks == peaks
    _assert_counts_add_up(bare_row, bare_peaks)


def test_no_satellites_excludes_a_satellite_that_matches_a_reflection(
    tmp_path, capsys
) -> None:
    argv, coincident = _coincident_argv(tmp_path)

    row, _ = _lattice_run(argv, tmp_path / "default")
    bare_row, bare_peaks = _lattice_run([*argv, "--no-satellites"], tmp_path / "bare")

    (peak,) = [
        peak
        for peak in bare_peaks
        if abs(float(peak["found_two_theta"]) - coincident) < 0.01
    ]
    assert peak["kalpha2_satellite"] == "True" and peak["recovered"] == "False"
    assert peak["fitted_two_theta"] == "" and peak["h"] == ""
    assert int(bare_row["n_peaks_recovered"]) == 0
    assert int(bare_row["n_satellites"]) == int(row["n_satellites"]) + 1
    assert int(bare_row["n_peaks_indexed"]) == int(row["n_peaks_indexed"]) - 1
    assert int(bare_row["n_peaks_refined"]) == int(row["n_peaks_refined"]) - 1
    _assert_counts_add_up(bare_row, bare_peaks)


def test_lattice_does_not_recover_a_satellite_whose_refit_is_rejected(
    tmp_path, capsys
) -> None:
    # 3 1 0's satellite alone, lying on the 2 2 2 position: the indexing would
    # recover it, but a doublet refit at 2 2 2 finds no K alpha 1 line there.
    argv, coincident = _coincident_argv(tmp_path, height=0.0)

    row, peaks = _lattice_run(argv, tmp_path / "default")
    bare_row, _ = _lattice_run([*argv, "--no-satellites"], tmp_path / "bare")
    line = capsys.readouterr().out.splitlines()[0]

    (peak,) = [
        peak
        for peak in peaks
        if abs(float(peak["found_two_theta"]) - coincident) < DEFAULT_TOLERANCE
    ]
    assert peak["kalpha2_satellite"] == "True" and peak["recovered"] == "False"
    assert peak["fit_rejected"] == "" and peak["fitted_two_theta"] == ""
    assert peak["h"] == "" and peak["difference"] == ""
    assert int(row["n_peaks_recovered"]) == 0
    assert int(row["n_fits_rejected"]) == 0
    assert f"{row['n_satellites']} satellites excluded, 0 recovered; " in line
    # Absent from the refinement: the run is the one that recovers nothing.
    for name in (
        "n_peaks_found",
        "n_peaks_refitted",
        "n_fits_rejected",
        "n_satellites",
        "n_peaks_recovered",
        "n_peaks_indexed",
        "n_peaks_refined",
        "a_angstrom",
        "c_angstrom",
        "rms_two_theta_deg",
    ):
        assert row[name] == bare_row[name], name
    _assert_counts_add_up(row, peaks)


# lebail and rietveld, with GSAS-II played by the fake of test_pipeline


@pytest.fixture
def refinement(tmp_path, monkeypatch):
    """A synthetic project in tmp_path, the working folder, with GSAS-II faked."""
    from test_pipeline import FakeGsas2, fake_project

    from xrdkit import pipeline

    project = fake_project(tmp_path)
    monkeypatch.chdir(project.root)
    gsas2 = FakeGsas2()
    monkeypatch.setattr(pipeline, "run_job", gsas2)
    return project


@pytest.mark.parametrize(
    "argv",
    [
        ["lebail"],
        ["lebail", "chain", "--mustrain"],
        ["rietveld", "chain", "--from", "lebail"],
        ["rietveld", "chain", "--through", "everything"],
        ["rietveld", "chain", "--preferred-orientation", "0", "1"],
    ],
    ids=[
        "no sample",
        "mustrain on lebail",
        "from lebail",
        "through unknown",
        "two indices",
    ],
)
def test_refinement_usage_errors(refinement, capsys, argv) -> None:
    with pytest.raises(SystemExit) as raised:
        main(argv)

    assert raised.value.code == 2
    assert "usage: xrdkit" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["lebail", "chain", "--zero", "0.1", "--displacement"],
            "xrdkit lebail: --zero and --displacement conflict",
        ),
        (
            ["lebail", "chain", "--system", "cubic"],
            "xrdkit lebail: --system needs --cell",
        ),
        (
            ["rietveld", "chain", "--from", "occupancies", "--through", "fixed_atoms"],
            "xrdkit rietveld: --from occupancies comes after --through fixed_atoms",
        ),
        (
            [
                "rietveld",
                "chain",
                "--from",
                "coordinates",
                "--preferred-orientation",
                "0",
                "0",
                "1",
            ],
            (
                "xrdkit rietveld: --preferred-orientation applies to the fixed_atoms "
                "mode, which --from coordinates leaves out"
            ),
        ),
        (
            ["rietveld", "chain", "--preferred-orientation", "0", "0", "0"],
            "xrdkit rietveld: --preferred-orientation takes an axis H K L not all 0",
        ),
        (
            ["lebail", "nope"],
            "xrdkit lebail: no sample nope in ",
        ),
        (
            ["lebail", "chain", "--two-theta", "60", "20"],
            "xrdkit lebail: --two-theta takes MIN MAX with MIN below MAX, not 60 20",
        ),
        (
            ["rietveld", "chain", "--two-theta", "100", "140"],
            (
                "xrdkit rietveld: --two-theta 100 140 lies outside the scan's "
                "20.00 to 60.00 degrees"
            ),
        ),
    ],
    ids=[
        "zero and displacement",
        "system without cell",
        "from after through",
        "orientation without fixed atoms",
        "orientation of zeros",
        "missing sample",
        "two theta the wrong way round",
        "two theta off the scan",
    ],
)
def test_refinement_conflicts_return_1(refinement, capsys, argv, message) -> None:
    assert main(argv) == 1

    error = capsys.readouterr().err
    assert error.startswith(message)
    assert len(error.strip().splitlines()) == 1
    assert not (refinement.root / "results").exists()


def test_lebail_missing_scan(refinement, capsys) -> None:
    (refinement.root / "toy.xrdml").unlink()

    assert main(["lebail", "chain"]) == 1

    scan = refinement.root / "toy.xrdml"
    assert capsys.readouterr().err.strip() == f"xrdkit lebail: no such file: {scan}"
    assert not (refinement.root / "results").exists()


def test_rietveld_needs_the_lebail_result(refinement, capsys) -> None:
    assert main(["rietveld", "chain"]) == 1

    result = (
        refinement.root / "results" / "lebail" / "chain" / "chain_lebail_result.json"
    )
    assert capsys.readouterr().err.strip() == (
        f"xrdkit rietveld: no such file: {result}; xrdkit lebail chain writes it"
    )
    assert not (refinement.root / "results").exists()


def test_lebail_prints_its_files_and_the_fit(refinement, capsys) -> None:
    assert main(["lebail", "x0.10.powder"]) == 0

    lines = capsys.readouterr().out.splitlines()
    folder = refinement.root / "results" / "lebail" / "x0.10.powder"
    assert lines[0] == "x0.10.powder: lebail started, 5 stages, at most 60 passes each"
    assert "x0.10.powder: lebail: background and scale clean" in lines
    for name in (
        "x0.10.powder_lebail_result.json",
        "x0.10.powder_lebail.gpx",
        "lebail.md",
        "x0.10.powder_lebail.png",
        "x0.10.powder_lebail.pdf",
        "summary.md",
    ):
        assert str(folder / name) in lines
        if name != "x0.10.powder_lebail.gpx":
            assert (folder / name).is_file()
    assert "bronze: tetragonal cell a = 6.0000, c = 4.0000 angstrom" in lines
    assert "size 0.4000 micron, microstrain 0" in lines
    assert "zero 0.0100 +/- 0.0010 degrees" in lines
    assert "Rwp 10.000 per cent, reduced chi squared 2.250" in lines
    assert "start cell of bronze from structures.bronze.cell" in lines


def _recorded_options(monkeypatch) -> list:
    """Every Options the handlers hand to run_sequence, which still runs."""
    seen: list = []
    real = cli.run_sequence

    def recorder(project, sample, modes, options=None, reporter=None):
        seen.append(options)
        return real(project, sample, modes, options, reporter)

    monkeypatch.setattr(cli, "run_sequence", recorder)
    return seen


def test_lebail_two_theta_reaches_the_options(refinement, capsys, monkeypatch) -> None:
    seen = _recorded_options(monkeypatch)

    assert main(["lebail", "chain", "--two-theta", "25", "55"]) == 0

    assert [options.two_theta for options in seen] == [(25.0, 55.0)]
    result = json.loads(
        (
            refinement.root
            / "results"
            / "lebail"
            / "chain"
            / "chain_lebail_result.json"
        ).read_text(encoding="utf-8")
    )
    assert result["inputs"]["limits"] == [25.0, 55.0]


def test_rietveld_two_theta_reaches_the_options(
    refinement, capsys, monkeypatch
) -> None:
    assert main(["lebail", "chain"]) == 0
    seen = _recorded_options(monkeypatch)

    assert main(["rietveld", "chain", "--two-theta", "25", "55"]) == 0

    assert [options.two_theta for options in seen] == [(25.0, 55.0)]
    folder = refinement.root / "results" / "rietveld" / "chain"
    # Every mode of the one call took the range given, not the Le Bail's.
    for mode in ("fixed_atoms", "coordinates", "occupancies"):
        result = json.loads(
            (folder / f"chain_{mode}_result.json").read_text(encoding="utf-8")
        )
        assert result["inputs"]["limits"] == [25.0, 55.0]


def test_two_theta_wider_than_the_scan_fits_all_of_it(refinement, capsys) -> None:
    """The sample table narrows the range; the scan's own limits widen it back."""
    assert main(["lebail", "chain", "--two-theta", "20", "60"]) == 0

    result = json.loads(
        (
            refinement.root
            / "results"
            / "lebail"
            / "chain"
            / "chain_lebail_result.json"
        ).read_text(encoding="utf-8")
    )
    assert result["inputs"]["limits"] == [20.0, 60.0]


def test_rietveld_from_coordinates(refinement, capsys) -> None:
    assert main(["lebail", "chain"]) == 0
    assert main(["rietveld", "chain", "--through", "fixed_atoms"]) == 0
    capsys.readouterr()

    assert main(["rietveld", "chain", "--from", "coordinates"]) == 0

    lines = capsys.readouterr().out.splitlines()
    assert (
        "coordinates: profile, Uiso groups, A sites, O sites; Rwp 10.000 per cent, reduced chi squared 2.250"
        in lines
    )
    assert (
        "occupancies: profile and Uiso, A site occupancies; Rwp 10.000 per cent, reduced chi squared 2.250"
        in lines
    )
    folder = refinement.root / "results" / "rietveld" / "chain"
    assert str(folder / "occupancies.md") in lines
    assert (folder / "summary.md").is_file()


def test_rietveld_prints_the_weight_fractions(refinement, capsys) -> None:
    assert main(["lebail", "two"]) == 0
    capsys.readouterr()

    assert main(["rietveld", "two", "--through", "fixed_atoms"]) == 0

    lines = capsys.readouterr().out.splitlines()
    assert any(
        line.startswith(
            "fixed_atoms: scale and background, zero and cell, size, overall Uiso;"
        )
        and line.endswith(
            "weight fractions toy 0.500 +/- 0.020, second 0.500 +/- 0.020"
        )
        for line in lines
    )


def test_rietveld_with_a_relative_out(refinement, monkeypatch, capsys) -> None:
    from test_pipeline import FakeGsas2

    from xrdkit import pipeline

    fake = FakeGsas2()

    def in_its_work_folder(job, workdir, install=None):
        # As run_job does, the work folder is resolved and the GSAS-II driver
        # runs in it, where a relative path names somewhere else.
        workdir = Path(workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        here = Path.cwd()
        os.chdir(workdir)
        try:
            return fake(job, workdir, install)
        finally:
            os.chdir(here)

    monkeypatch.setattr(pipeline, "run_job", in_its_work_folder)

    assert main(["lebail", "chain", "--out", "elsewhere"]) == 0
    assert (
        main(["rietveld", "chain", "--through", "fixed_atoms", "--out", "elsewhere"])
        == 0
    )

    folder = refinement.root / "elsewhere"
    assert (folder / "chain_lebail_result.json").is_file()
    assert (folder / "chain_fixed_atoms_result.json").is_file()
    assert str(folder / "chain_fixed_atoms_result.json") in capsys.readouterr().out
    assert all(
        Path(job[key]).is_absolute()
        for job in fake.jobs
        if job["action"] == "refine"
        for key in ("result", "export_prefix")
    )


def test_rietveld_refuses_to_overwrite_a_result(refinement, capsys) -> None:
    assert main(["lebail", "chain"]) == 0
    assert main(["rietveld", "chain", "--through", "fixed_atoms"]) == 0
    folder = refinement.root / "results" / "rietveld" / "chain"
    result = folder / "chain_fixed_atoms_result.json"
    before = result.read_bytes()
    result.write_bytes(before + b"\n")
    capsys.readouterr()

    assert main(["rietveld", "chain", "--through", "fixed_atoms"]) == 1

    captured = capsys.readouterr()
    assert captured.err.strip() == (
        f"xrdkit rietveld: {result} exists already; give --overwrite to replace "
        "the results of the modes run, or --out DIR to write them elsewhere"
    )
    assert captured.out == ""
    assert result.read_bytes() == before + b"\n"

    assert main(["rietveld", "chain", "--through", "fixed_atoms", "--overwrite"]) == 0
    assert result.read_bytes() != before + b"\n"


def test_a_failed_mode_returns_1_with_its_failure_record(
    refinement, monkeypatch, capsys
) -> None:
    from test_pipeline import FakeGsas2

    from xrdkit import pipeline

    monkeypatch.setattr(pipeline, "run_job", FakeGsas2(status="rejected"))

    assert main(["lebail", "chain"]) == 1

    captured = capsys.readouterr()
    failure = refinement.root / "results" / "lebail" / "chain" / "failure.md"
    assert failure.is_file()
    assert str(failure) in captured.out.splitlines()
    assert captured.err.strip().startswith(
        "xrdkit lebail: lebail failed: PipelineError:"
    )
    assert captured.err.strip().endswith(f"see {failure}")


def test_lebail_json(refinement, capsys) -> None:
    assert main(["lebail", "chain", "--json"]) == 0

    captured = capsys.readouterr()
    output = json.loads(captured.out)
    assert output["sample"] == "chain"
    (outcome,) = output["outcomes"]
    assert outcome["mode"] == "lebail" and outcome["error"] is None
    assert outcome["residuals"]["rwp"] == 10.0
    assert any(path.endswith("chain_lebail_result.json") for path in output["files"])
    # The progress lines go to stderr, so that stdout is the JSON alone.
    assert "chain: lebail started" in captured.err


def _as_xy(xrdml: Path, path: Path) -> Path:
    """The same pattern as ``xrdml`` written as two columns, the layout a
    converter produces. The tests write their own: .gitignore ignores *.xy."""
    scan = read_xrdml(xrdml)
    rows = (f"{t:.5f} {i:.3f}" for t, i in zip(scan.two_theta, scan.intensity))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _on_xy(argv: list[str], xrdml: Path, xy: Path) -> list[str]:
    """``argv`` with the scan swapped for the .xy copy."""
    return [str(xy) if part == str(xrdml) else part for part in argv]


def test_lattice_on_an_xy_matches_the_equivalent_xrdml(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)
    assert main([*argv, "--json"]) == 0
    from_xrdml = json.loads(capsys.readouterr().out)
    xrdml = tmp_path / "ttb.xrdml"
    xy = _as_xy(xrdml, tmp_path / "ttb.xy")

    argv = _on_xy(argv, xrdml, xy) + ["--wavelength", str(LATTICE_WAVELENGTH)]
    assert main([*argv, "--json"]) == 0

    result = json.loads(capsys.readouterr().out)
    assert result["wavelength_angstrom"] == LATTICE_WAVELENGTH
    assert result["n_peaks_indexed"] == from_xrdml["n_peaks_indexed"]
    for key in ("a_angstrom", "c_angstrom", "volume_a3", "zero_deg"):
        assert result[key] == pytest.approx(from_xrdml[key], rel=1e-5), key


def test_lattice_refuses_an_xy_with_no_wavelength(tmp_path, capsys) -> None:
    argv, _ = _ttb_lattice_argv(tmp_path)
    xrdml = tmp_path / "ttb.xrdml"
    xy = _as_xy(xrdml, tmp_path / "ttb.xy")

    assert main(_on_xy(argv, xrdml, xy)) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        f"xrdkit lattice: {xy} carries no wavelength; give --wavelength ANGSTROM "
        "or a sample key"
    )
    assert not (tmp_path / "out").exists()


def test_check_reads_an_xy_with_no_wavelength(tmp_path, capsys) -> None:
    # check never uses the wavelength, so a bare .xy needs nothing.
    xy = _as_xy(_write_xrdml(tmp_path / "scan.xrdml"), tmp_path / "scan.xy")

    assert main(["check", str(xy)]) == 0

    report = capsys.readouterr().out
    assert "range   10.00 to 100.00 degrees" in report
    assert "time    unknown" in report


def test_lattice_of_an_xy_sample_takes_the_instrument_wavelength(
    project, capsys
) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(
            '\n[structures.sto]\nlibrary = "perovskite/Pm-3m"\n'
            'composition = "SrTiO3"\ncell = { a = 3.95 }\n'
        )
    raw = project / "data" / "raw"
    _write_doublet_xrdml(raw / "sto.xrdml", Cell.cubic(3.905), "Pm-3m")
    _as_xy(raw / "sto.xrdml", raw / "sto.xy")
    assert main(["add-sample", "data/raw/sto.xy", "--structure", "sto"]) == 0
    capsys.readouterr()

    assert main(["lattice", "sto"]) == 0

    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "sample  sto, SrTiO3"
    assert lines[3].startswith("cubic cell a = 3.905")
    (row,) = _rows(project / "results" / "lattice" / "sto" / "lattice_sto.csv")
    assert float(row["a_angstrom"]) == pytest.approx(3.905, abs=0.001)


def test_lebail_reads_an_xy_sample(refinement, capsys) -> None:
    root = refinement.root
    _as_xy(root / "toy.xrdml", root / "toy.xy")
    path = root / PROJECT_FILE
    path.write_text(
        path.read_text(encoding="utf-8").replace('"toy.xrdml"', '"toy.xy"'),
        encoding="utf-8",
    )

    assert main(["lebail", "chain"]) == 0

    result = root / "results" / "lebail" / "chain" / "chain_lebail_result.json"
    assert result.is_file()
    assert str(result) in capsys.readouterr().out


# xrdkit instrument


class FakeInstrumentGsas2:
    """Plays run_job for the instrument command: it writes the exported
    instrument parameter file the command copies into place, and a result of
    the stages the job asked for, keeping every job it was given."""

    def __init__(self, rwp: float = 8.125, gof: float = 1.406):
        self.jobs = []
        self.rwp = rwp
        self.gof = gof

    def __call__(self, job, workdir, install=None):
        self.jobs.append(copy.deepcopy(job))
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "refine.log").write_text("fake GSAS-II log\n", encoding="utf-8")
        prefix = Path(job["export_prefix"])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        exported = prefix.with_name(prefix.name + ".instprm")
        # GSAS-II exports the file with the refined values in it; the starting
        # file the job points at is copied and its Zero replaced, which is
        # enough for the command to have something real to copy.
        start = Path(job["instprm"]).read_text(encoding="utf-8")
        exported.write_text(start.replace("Zero:0.0", "Zero:0.0123"), encoding="utf-8")
        values = {
            "Zero": 0.0123,
            "U": 6.01,
            "V": -3.52,
            "W": 4.49,
            "X": 0.5,
            "Y": -1.6,
            "SH/L": 0.0221,
        }
        return {
            "stages": [
                {
                    "name": stage["name"],
                    "status": "clean",
                    "rwp": self.rwp,
                    "gof": self.gof,
                }
                for stage in job["stages"]
            ],
            "completed": True,
            "final": {
                "instrument": {
                    key: {"value": value, "esd": abs(value) * 0.01 + 0.001}
                    for key, value in values.items()
                }
            },
            "exports": {"instprm": str(exported)},
        }


@pytest.fixture
def standard(tmp_path, monkeypatch):
    """A folder with a synthetic standard scan and its CIF, GSAS-II faked."""
    from xrdkit import instrument as instrument_module

    monkeypatch.chdir(tmp_path)
    (tmp_path / "cifs").mkdir()
    (tmp_path / "data" / "standards").mkdir(parents=True)
    _write_standard_xrdml(tmp_path / "data" / "standards" / "lab6.xrdml")
    (tmp_path / "cifs" / "lab6.cif").write_text(STANDARD_CIF, encoding="utf-8")
    monkeypatch.setattr(
        cli, "find_gsas2", lambda: Gsas2Install(Path("py"), Path("home"))
    )
    fake = FakeInstrumentGsas2()
    monkeypatch.setattr(instrument_module, "run_job", fake)
    return fake


def _instrument_argv(**extra) -> list[str]:
    argv = [
        "instrument",
        "data/standards/lab6.xrdml",
        "--cif",
        "cifs/lab6.cif",
        "--cell",
        "4.15683",
    ]
    for key, value in extra.items():
        argv.append(f"--{key.replace('_', '-')}")
        if value is not True:
            argv.extend(value if isinstance(value, list) else [str(value)])
    return argv


def test_instrument_writes_the_file_and_its_record(tmp_path, standard, capsys) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"))) == 0

    folder = tmp_path / "out"
    instprm, record = folder / "lab6.instprm", folder / "lab6_instrument.csv"
    lines = capsys.readouterr().out.splitlines()
    assert lines[-2:] == [str(instprm), str(record)]
    assert instprm.is_file() and record.is_file()
    # The starting file and the GSAS-II project are kept, under their own names.
    assert (folder / "lab6_start.instprm").is_file()
    assert (folder / "work" / "lab6").is_dir()
    assert "Zero:0.0123" in instprm.read_text(encoding="utf-8")


def test_instrument_record_columns(tmp_path, standard) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"), radius=145.0)) == 0

    (row,) = _rows(tmp_path / "out" / "lab6_instrument.csv")
    assert list(row) == [
        "scan",
        "cif",
        "phase",
        "wavelength_angstrom",
        "a_angstrom",
        "b_angstrom",
        "c_angstrom",
        "alpha_deg",
        "beta_deg",
        "gamma_deg",
        "window_min_deg",
        "window_max_deg",
        "n_peaks_fitted",
        "fit_u_deg2",
        "fit_v_deg2",
        "fit_w_deg2",
        "fit_rms_deg",
        "refined_zero",
        "esd_zero",
        "refined_u",
        "esd_u",
        "refined_v",
        "esd_v",
        "refined_w",
        "esd_w",
        "refined_x",
        "esd_x",
        "refined_y",
        "esd_y",
        "refined_sh_l",
        "esd_sh_l",
        "rwp_percent",
        "gof",
        "radius_mm",
        "instprm",
        "method",
        "date",
        "xrdkit_version",
    ]
    assert row["phase"] == "lab6"
    assert float(row["a_angstrom"]) == pytest.approx(4.15683)
    assert float(row["window_min_deg"]) == 10.0
    assert float(row["window_max_deg"]) == 98.0
    assert int(row["n_peaks_fitted"]) >= 6
    assert float(row["refined_zero"]) == pytest.approx(0.0123)
    assert float(row["rwp_percent"]) == pytest.approx(8.125)
    assert float(row["radius_mm"]) == 145.0
    assert row["method"].startswith("Caglioti width fit")
    assert row["xrdkit_version"] == xrdkit.__version__


def test_instrument_radius_reaches_the_instprm(tmp_path, standard) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"), radius=145.0)) == 0

    start = (tmp_path / "out" / "lab6_start.instprm").read_text(encoding="utf-8")
    assert "Gonio. radius:145.0" in start.splitlines()


def test_instrument_without_a_radius_writes_no_radius_line(tmp_path, standard) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"))) == 0

    start = (tmp_path / "out" / "lab6_start.instprm").read_text(encoding="utf-8")
    assert not [line for line in start.splitlines() if line.startswith("Gonio.")]


def test_instrument_prints_the_fit_and_the_parameters(
    tmp_path, standard, capsys
) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"), window=["10", "80"])) == 0

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("peaks   ") and "from 10 to 80 degrees" in lines[0]
    assert lines[1].startswith("U = ") and "degrees squared" in lines[1]
    assert lines[2].startswith("rms of the width fit ")
    assert lines[3] == "stages  background and scale, zero, U V W, X Y, SH/L"
    assert lines[4] == "Rwp 8.125 per cent, GOF 1.406"
    assert [line.split()[0] for line in lines[5:12]] == list(REFINED_KEYS)
    assert "esd" in lines[5]


def test_instrument_json_puts_progress_on_stderr(tmp_path, standard, capsys) -> None:
    assert main(_instrument_argv(out=str(tmp_path / "out"), json=True)) == 0

    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["files"] == [
        str(tmp_path / "out" / "lab6.instprm"),
        str(tmp_path / "out" / "lab6_instrument.csv"),
    ]
    assert result["rwp_percent"] == pytest.approx(8.125)
    assert result["gof"] == pytest.approx(1.406)
    assert result["stages"][0] == "background and scale"
    assert result["parameters"]["Zero"]["value"] == pytest.approx(0.0123)
    assert result["n_peaks_fitted"] >= 6
    assert "Rwp 8.125 per cent" in captured.err


def test_instrument_stem_and_phase(tmp_path, standard) -> None:
    argv = _instrument_argv(out=str(tmp_path / "out"), stem="lab", phase="LaB6")

    assert main(argv) == 0

    assert (tmp_path / "out" / "lab.instprm").is_file()
    (row,) = _rows(tmp_path / "out" / "lab_instrument.csv")
    assert row["phase"] == "LaB6"
    (job,) = [j for j in standard.jobs if j.get("phases")]
    assert job["phases"][0]["name"] == "LaB6"
    assert job["phases"][0]["cell"] == [4.15683] * 3 + [90.0] * 3
    assert [stage["name"] for stage in job["stages"]] == [
        "background and scale",
        "zero",
        "U V W",
        "X Y",
        "SH/L",
    ]


# --name, and the project file


def _instrument_project(tmp_path, monkeypatch) -> Path:
    root = (tmp_path / "project").resolve()
    root.mkdir()
    (root / "cifs").mkdir()
    (root / "data" / "standards").mkdir(parents=True)
    _write_standard_xrdml(root / "data" / "standards" / "lab6.xrdml")
    (root / "cifs" / "lab6.cif").write_text(STANDARD_CIF, encoding="utf-8")
    (root / PROJECT_FILE).write_text(
        '[project]\nname = "demo"\nversion = 1\n', encoding="utf-8"
    )
    monkeypatch.chdir(root)
    return root


def test_instrument_appends_the_instrument_table(
    tmp_path, standard, monkeypatch, capsys
) -> None:
    root = _instrument_project(tmp_path, monkeypatch)
    path = root / PROJECT_FILE
    before = path.read_bytes()

    assert main(_instrument_argv(name="lab_diffractometer", radius=145.0)) == 0

    table = [
        "[instruments.lab_diffractometer]",
        "wavelength = [1.540598, 1.544426]",
        "ka2 = true",
        "radius = 145",
        'instprm = "data/standards/lab6.instprm"',
    ]
    out = capsys.readouterr().out.splitlines()
    start = out.index("[instruments.lab_diffractometer]")
    assert out[start : start + len(table)] == table
    after = path.read_bytes()
    assert after[: len(before)] == before
    # Appended in the file's own line endings, as add-sample does.
    newline = "\r\n" if b"\r\n" in before else "\n"
    assert newline.join(table) in after.decode("utf-8")
    instrument = load_project(root).instruments["lab_diffractometer"]
    assert instrument.wavelength == (1.540598, 1.544426)
    assert instrument.ka2 is True
    assert instrument.radius == 145.0
    assert instrument.instprm == root / "data" / "standards" / "lab6.instprm"
    # The default folder of a project is data/standards under its root.
    assert (root / "data" / "standards" / "lab6.instprm").is_file()


def test_instrument_refuses_an_instrument_key_that_exists(
    tmp_path, standard, monkeypatch, capsys
) -> None:
    root = _instrument_project(tmp_path, monkeypatch)
    with (root / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(
            "\n[instruments.lab_diffractometer]\nwavelength = [1.540598]\nka2 = false\n"
        )
    before = (root / PROJECT_FILE).read_bytes()

    assert main(_instrument_argv(name="lab_diffractometer")) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        f"xrdkit instrument: instrument 'lab_diffractometer' is in {root / PROJECT_FILE} "
        "already; give another --name"
    )
    assert (root / PROJECT_FILE).read_bytes() == before
    assert not (root / "data" / "standards" / "lab6.instprm").exists()


def test_instrument_refuses_name_with_no_project_file(
    tmp_path, standard, capsys
) -> None:
    assert main(_instrument_argv(name="lab_diffractometer")) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        f"xrdkit instrument: no {PROJECT_FILE} in the current folder or above it "
        "to add [instruments.lab_diffractometer] to; run xrdkit init"
    )
    assert not (tmp_path / "data" / "standards" / "lab6.instprm").exists()


def test_instrument_refuses_an_xy_with_no_wavelength(
    tmp_path, standard, capsys
) -> None:
    scan = tmp_path / "data" / "standards" / "lab6.xy"
    source = read_xrdml(tmp_path / "data" / "standards" / "lab6.xrdml")
    scan.write_text(
        "\n".join(
            f"{t:.5f} {i:.3f}" for t, i in zip(source.two_theta, source.intensity)
        )
        + "\n",
        encoding="utf-8",
    )
    argv = ["instrument", str(scan), "--cif", "cifs/lab6.cif", "--cell", "4.15683"]

    assert main(argv) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        f"xrdkit instrument: {scan} carries no wavelength; give --wavelength "
        "ANGSTROM or a sample key"
    )


def test_instrument_refuses_when_gsas2_is_absent(
    tmp_path, standard, monkeypatch, capsys
) -> None:
    def absent():
        raise FileNotFoundError("no GSAS-II: set XRDKIT_GSAS2_PYTHON")

    monkeypatch.setattr(cli, "find_gsas2", absent)

    assert main(_instrument_argv(out=str(tmp_path / "out"))) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        "xrdkit instrument: no GSAS-II: set XRDKIT_GSAS2_PYTHON"
    )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["instrument", "data/standards/lab6.xrdml", "--cif", "cifs/lab6.cif"],
            "--cell is required",
        ),
        (_instrument_argv() + ["--window", "98", "10"], "--window takes MIN MAX"),
        (
            [
                "instrument",
                "data/standards/lab6.xrdml",
                "--cif",
                "cifs/absent.cif",
                "--cell",
                "4.15683",
            ],
            "no such file: cifs/absent.cif",
        ),
    ],
    ids=["no cell", "backwards window", "no cif"],
)
def test_instrument_option_errors(tmp_path, standard, capsys, argv, message) -> None:
    assert main([*argv, "--out", str(tmp_path / "out")]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("xrdkit instrument: ")
    assert message in captured.err
    assert not (tmp_path / "out").exists()


def test_instrument_with_gsas2(tmp_path, monkeypatch) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    (tmp_path / "cifs").mkdir()
    (tmp_path / "data" / "standards").mkdir(parents=True)
    _write_standard_xrdml(tmp_path / "data" / "standards" / "lab6.xrdml")
    (tmp_path / "cifs" / "lab6.cif").write_text(STANDARD_CIF, encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    assert main(_instrument_argv(out=str(tmp_path / "out"), phase="LaB6")) == 0

    instprm = tmp_path / "out" / "lab6.instprm"
    assert instprm.is_file()
    assert install.home.is_dir()
    (row,) = _rows(tmp_path / "out" / "lab6_instrument.csv")
    # The pattern is drawn by multiplicity rather than from structure
    # factors, so the residual is large; that it refined at all is the point.
    assert 0.0 < float(row["rwp_percent"]) < 100.0
    assert float(row["gof"]) > 0.0
    assert float(row["refined_w"]) != 0.0
    values = dict(line.split(":", 1) for line in instprm.read_text().splitlines()[1:])
    assert float(values["Lam1"]) == pytest.approx(1.540598, abs=1e-4)


# xrdkit phases

COD_CIF = """data_candidate
_cell_length_a  5.640
_symmetry_space_group_name_H-M  "F m -3 m"
"""


def _cod_record(cod_id: str, formula: str = "Na Cl") -> CodRecord:
    return CodRecord(
        cod_id=cod_id,
        formula=formula,
        space_group="F m -3 m",
        space_group_number=225,
        a=5.64,
        b=5.64,
        c=5.64,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        volume=179.4,
        authors="A Person",
        journal="Acta",
        year=2006,
    )


# One line per candidate, keyed by COD id, standing in for pymatgen.
COD_PATTERNS = {
    "2100720": [(20.0, 100.0, (1, 1, 1)), (30.0, 60.0, (2, 0, 0))],
    "2100721": [(55.0, 100.0, (2, 2, 0)), (20.0, 40.0, (1, 1, 1))],
    "2100722": [
        (20.0, 100.0, (1, 1, 1)),
        (26.8, 8.0, (2, 0, 1)),
        (30.0, 60.0, (2, 0, 0)),
    ],
}


@pytest.fixture
def cod(tmp_path, monkeypatch):
    """The COD faked: a fixed search result, CIFs written on fetch, and
    simulate_pattern answering from COD_PATTERNS, so no network and no
    pymatgen are needed."""
    from xrdkit import phases as phases_module

    fetched: list[str] = []

    def fake_search(elements, exact=True, space_group=None, extra=None):
        fake_search.asked = (list(elements), exact, space_group)
        return [_cod_record(key) for key in sorted(COD_PATTERNS)]

    def fake_fetch(cod_id, folder):
        fetched.append(str(cod_id))
        path = Path(folder) / f"{cod_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(COD_CIF, encoding="utf-8")
        return path

    def fake_simulate(cif_path, wavelength=1.0, two_theta_range=(10, 100)):
        lines = COD_PATTERNS[Path(cif_path).stem]
        return [SimulatedReflection(*line) for line in lines]

    monkeypatch.setattr(cli, "cod_search", fake_search)
    monkeypatch.setattr(cli, "cod_fetch", fake_fetch)
    monkeypatch.setattr(phases_module, "cod_fetch", fake_fetch)
    monkeypatch.setattr(phases_module, "simulate_pattern", fake_simulate)
    monkeypatch.setattr(cli, "require_phases_extra", lambda: None)
    fake_search.asked = None
    fake_search.fetched = fetched
    return fake_search


def _phases_scan(path: Path) -> Path:
    """A scan whose peaks sit where COD_PATTERNS puts its lines, plus one more."""
    two_theta = np.arange(10.0, 80.0, 0.02)
    counts = np.full(two_theta.size, 100.0)
    for centre in (20.0, 30.0, 26.8):
        counts += 5000.0 * np.exp(-0.5 * ((two_theta - centre) / 0.05) ** 2)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        XRDML.format(
            start=two_theta[0],
            end=two_theta[-1],
            counts=" ".join(str(round(v)) for v in counts),
        ),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def phases_scan(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return _phases_scan(tmp_path / "data" / "raw" / "sample.xrdml")


def _phases_argv(scan: Path, *extra: str) -> list[str]:
    return ["phases", str(scan), "--elements", "Na", "Cl", *extra]


def _warn_to_stderr(message, category, filename, lineno, file=None, line=None):
    """The warnings module's own handler, which pytest replaces while it
    records warnings; the test puts it back so that a warning reaching
    stderr can be seen."""
    sys.stderr.write(warnings.formatwarning(message, category, filename, lineno, line))


def test_phases_keeps_pymatgens_warnings_off_stderr(
    tmp_path, cod, phases_scan, monkeypatch, capsys
) -> None:
    """A UserWarning raised inside the simulation, which is what pymatgen's
    CIF parser does for a doubtful COD entry, does not reach stderr."""
    from xrdkit import phases as phases_module

    quiet = phases_module.simulate_pattern

    def noisy(cif_path, wavelength=1.0, two_theta_range=(10, 100)):
        warnings.warn(
            "Issues encountered while parsing CIF: Skipping relative "
            "stoichiometry check because CIF does not contain formula keys.",
            stacklevel=2,
        )
        return quiet(cif_path, wavelength=wavelength, two_theta_range=two_theta_range)

    monkeypatch.setattr(phases_module, "simulate_pattern", noisy)
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        warnings.showwarning = _warn_to_stderr
        assert main(_phases_argv(phases_scan)) == 0

    assert capsys.readouterr().err == ""


def test_phases_writes_its_files(tmp_path, cod, phases_scan, capsys) -> None:
    assert main(_phases_argv(phases_scan)) == 0

    # Without a project file the paths are relative to the working folder.
    cifs, folder = Path("cifs") / "cod", Path("results") / "phases" / "sample"
    lines = capsys.readouterr().out.splitlines()
    assert lines[-5:] == [
        str(cifs / "2100722.cif"),
        str(cifs / "index.csv"),
        str(folder / "phases_sample.csv"),
        str(folder / "phases_sample_unexplained.csv"),
        str(folder / "phases_sample_record.csv"),
    ]
    cifs = tmp_path / "cifs" / "cod"
    assert sorted(p.name for p in cifs.glob("*.cif")) == [
        "2100720.cif",
        "2100721.cif",
        "2100722.cif",
    ]


def test_phases_candidate_columns_and_ranking(tmp_path, cod, phases_scan) -> None:
    assert main(_phases_argv(phases_scan)) == 0

    folder = tmp_path / "results" / "phases" / "sample"
    rows = _rows(folder / "phases_sample.csv")
    assert list(rows[0]) == [
        "rank",
        "cod_id",
        "formula",
        "space_group",
        "n_explained",
        "n_missing",
        "score",
        "rejected",
        "cif",
    ]
    # 2100722 explains all three peaks, 2100720 two of them, and 2100721's
    # strongest line is nowhere in the scan.
    assert [row["cod_id"] for row in rows] == ["2100722", "2100720", "2100721"]
    assert [row["rank"] for row in rows] == ["1", "2", "3"]
    assert rows[0]["n_explained"] == "3" and rows[0]["rejected"] == ""
    assert rows[2]["rejected"] == "strongest line absent"
    assert rows[0]["cif"].endswith("2100722.cif")


def test_phases_record_and_unexplained_columns(tmp_path, cod, phases_scan) -> None:
    assert main(_phases_argv(phases_scan, "--zero", "0.0", "--tolerance", "0.2")) == 0

    folder = tmp_path / "results" / "phases" / "sample"
    (record,) = _rows(folder / "phases_sample_record.csv")
    assert list(record) == [
        "scan",
        "sample",
        "wavelength_angstrom",
        "elements",
        "space_group",
        "zero_deg",
        "window_min_deg",
        "window_max_deg",
        "tolerance_deg",
        "n_peaks_observed",
        "n_candidates_searched",
        "n_candidates_fetched",
        "n_unexplained",
        "main",
        "index",
        "method",
        "date",
        "xrdkit_version",
    ]
    assert record["elements"] == "Na Cl"
    assert float(record["window_min_deg"]) == 10.0
    assert float(record["window_max_deg"]) == 80.0
    assert float(record["tolerance_deg"]) == 0.2
    assert int(record["n_peaks_observed"]) == 3
    assert int(record["n_candidates_searched"]) == 3
    assert record["main"] == "2100722"
    assert record["method"].startswith("COD search by element set")
    assert record["xrdkit_version"] == xrdkit.__version__
    unexplained = _rows(folder / "phases_sample_unexplained.csv")
    assert list(unexplained[0] if unexplained else {}) == [] or list(
        unexplained[0]
    ) == [
        "two_theta_deg",
        "d_angstrom",
        "phase",
        "hkl",
        "reflection_two_theta_deg",
        "intensity_percent",
    ]


def test_phases_defaults_come_from_the_guide(tmp_path, cod, phases_scan) -> None:
    assert main(_phases_argv(phases_scan)) == 0

    (record,) = _rows(
        tmp_path / "results" / "phases" / "sample" / "phases_sample_record.csv"
    )
    assert (float(record["window_min_deg"]), float(record["window_max_deg"])) == (
        10.0,
        80.0,
    )
    assert float(record["tolerance_deg"]) == 0.15
    assert float(record["zero_deg"]) == 0.0


def test_phases_elements_from_a_sample_key(tmp_path, cod, project, capsys) -> None:
    _phases_scan(project / "data" / "raw" / "phase_sample.xrdml")
    assert (
        main(["add-sample", "data/raw/phase_sample.xrdml", "--structure", "bto"]) == 0
    )
    capsys.readouterr()

    assert main(["phases", "phase_sample"]) == 0

    elements, exact, space_group = cod.asked
    # BaTiO3 of the bto structure, in the order the formula gives.
    assert elements == ["Ba", "Ti", "O"]
    assert exact is True and space_group is None
    folder = project / "results" / "phases" / "phase_sample"
    assert (folder / "phases_phase_sample.csv").is_file()
    (record,) = _rows(folder / "phases_phase_sample_record.csv")
    assert record["sample"] == "phase_sample"
    assert record["elements"] == "Ba Ti O"


def test_phases_refuses_a_bare_file_without_elements(
    tmp_path, cod, phases_scan, capsys
) -> None:
    assert main(["phases", str(phases_scan)]) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        "xrdkit phases: --elements is required for a scan named directly; a "
        "sample key takes them from the compositions of its structures"
    )
    assert not (tmp_path / "cifs").exists()


def test_phases_refuses_without_the_extra(
    tmp_path, cod, phases_scan, monkeypatch, capsys
) -> None:
    def absent():
        raise MissingPhasesExtra("simulate_pattern needs pymatgen")

    monkeypatch.setattr(cli, "require_phases_extra", absent)

    assert main(_phases_argv(phases_scan)) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        "xrdkit phases: phase identification needs the phases extra; install "
        "with pip install xrdkit[phases]"
    )
    # Refused before the COD was asked anything.
    assert cod.asked is None
    assert not (tmp_path / "cifs").exists()


def test_phases_reports_a_search_that_could_not_reach_the_cod(
    tmp_path, cod, phases_scan, monkeypatch, capsys
) -> None:
    def unreachable(elements, exact=True, space_group=None, extra=None):
        raise URLError("getaddrinfo failed")

    monkeypatch.setattr(cli, "cod_search", unreachable)

    assert main(_phases_argv(phases_scan)) == 1

    captured = capsys.readouterr()
    assert captured.err.strip().startswith("xrdkit phases: the COD search failed: ")
    assert "getaddrinfo failed" in captured.err
    assert not (tmp_path / "cifs").exists()


def test_phases_reports_a_fetch_that_could_not_reach_the_cod(
    tmp_path, cod, phases_scan, monkeypatch, capsys
) -> None:
    def unreachable(cod_id, folder):
        raise URLError("connection reset")

    monkeypatch.setattr(cli, "cod_fetch", unreachable)
    from xrdkit import phases as phases_module

    monkeypatch.setattr(phases_module, "cod_fetch", unreachable)

    assert main(_phases_argv(phases_scan)) == 1

    captured = capsys.readouterr()
    assert captured.err.strip().startswith(
        "xrdkit phases: fetching a CIF from the COD failed: "
    )
    assert "connection reset" in captured.err


def test_phases_does_not_fetch_a_cif_already_there(tmp_path, cod, phases_scan) -> None:
    cifs = tmp_path / "cifs" / "cod"
    cifs.mkdir(parents=True)
    (cifs / "2100720.cif").write_text(COD_CIF, encoding="utf-8")

    assert main(_phases_argv(phases_scan)) == 0

    assert cod.fetched == ["2100721", "2100722"]


def test_phases_max_candidates(tmp_path, cod, phases_scan) -> None:
    assert main(_phases_argv(phases_scan, "--max-candidates", "2")) == 0

    rows = _rows(tmp_path / "results" / "phases" / "sample" / "phases_sample.csv")
    assert {row["cod_id"] for row in rows} == {"2100720", "2100721"}
    assert cod.fetched == ["2100720", "2100721"]


def test_phases_main_names_the_phase_of_the_unexplained(
    tmp_path, cod, phases_scan, capsys
) -> None:
    # 2100720 explains 20 and 30 but not 26.8, which 2100722 accounts for.
    assert (
        main(_phases_argv(phases_scan, "--max-candidates", "1", "--main", "2100720"))
        == 0
    )

    folder = tmp_path / "results" / "phases" / "sample"
    (record,) = _rows(folder / "phases_sample_record.csv")
    assert record["main"] == "2100720"
    rows = _rows(folder / "phases_sample_unexplained.csv")
    assert [float(row["two_theta_deg"]) for row in rows] == [
        pytest.approx(26.8, abs=0.05)
    ]
    assert rows[0]["phase"] == "" and rows[0]["hkl"] == ""
    assert "unidentified" in capsys.readouterr().out


def test_phases_json_puts_progress_on_stderr(
    tmp_path, cod, phases_scan, capsys
) -> None:
    assert main(_phases_argv(phases_scan, "--json")) == 0

    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["n_peaks_observed"] == 3
    assert result["elements"] == ["Na", "Cl"]
    assert result["main"] == "2100722"
    assert [c["cod_id"] for c in result["candidates"]] == [
        "2100722",
        "2100720",
        "2100721",
    ]
    assert result["files"][-1].endswith("phases_sample_record.csv")
    assert "peaks observed from 10 to 80 degrees" in captured.err


def test_phases_refuses_an_xy_with_no_wavelength(
    tmp_path, cod, phases_scan, capsys
) -> None:
    scan = tmp_path / "data" / "raw" / "sample.xy"
    source = read_xrdml(phases_scan)
    scan.write_text(
        "\n".join(
            f"{t:.5f} {i:.3f}" for t, i in zip(source.two_theta, source.intensity)
        )
        + "\n",
        encoding="utf-8",
    )

    assert main(_phases_argv(scan)) == 1

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() == (
        f"xrdkit phases: {scan} carries no wavelength; give --wavelength "
        "ANGSTROM or a sample key"
    )
