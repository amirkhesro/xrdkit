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
    reflections = generate_reflections(REFINED_CELL, 1.540598, 80.0)
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
        peaks, TetragonalCell(a=12.45, c=3.94), scan.wavelength
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
