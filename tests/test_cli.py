"""Tests for xrdkit.cli."""

import json
from pathlib import Path

import numpy as np
import pytest

from xrdkit.cli import main
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
