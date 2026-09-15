"""Tests for xrdkit.cli."""

import csv
import datetime
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg

from xrdkit import (
    TetragonalCell,
    cli,
    find_peaks,
    generate_reflections,
    index_and_refine,
    read_xrdml,
)
from xrdkit.cli import HKL_HEADROOM, main
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


def _write_indexable_xrdml(path: Path) -> tuple[Path, int]:
    """Write the shifted pattern as a scan, returning it and its peak count."""
    reflections = generate_reflections(REFINED_CELL, 1.540598, 80.0, space_group="P4bm")
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
[instruments.aeris]
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
        'instrument = "aeris"',
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
        instrument="aeris",
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


@pytest.mark.parametrize(
    ("structure", "note"),
    [
        ("bfo", "hkl labelling for a trigonal cell arrives in the next release"),
    ],
)
def test_plot_by_key_notes_a_cell_it_cannot_label_yet(
    project, capsys, structure, note
) -> None:
    with (project / PROJECT_FILE).open("a", encoding="utf-8") as handle:
        handle.write(BFO_STRUCTURE)
    argv = ["add-sample", "data/raw/10c.xrdml", "--structure", structure]
    assert main(argv) == 0
    capsys.readouterr()

    assert main(["plot", "10c"]) == 0

    captured = capsys.readouterr()
    assert captured.err == (f"xrdkit plot: {note}; plotted without hkl labels\n")
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
