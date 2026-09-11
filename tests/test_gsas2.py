"""Tests for xrdkit.gsas2.

None of these import GSASII: the job runner is exercised with a fake Python
that stands in for GSAS-II's, and the one test that runs GSAS-II itself is
skipped where it is not installed.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from xrdkit.broadening import Caglioti
from xrdkit.gsas2 import (
    GSAS2_HOME_VARIABLE,
    GSAS2_PYTHON_VARIABLE,
    Gsas2Error,
    Gsas2Install,
    find_gsas2,
    run_job,
    write_instprm,
)

CAGLIOTI = Caglioti(
    u=1.047e-02,
    v=-1.184e-02,
    w=8.359e-03,
    esd_u=1.714e-03,
    esd_v=2.049e-03,
    esd_w=5.663e-04,
    n_peaks=14,
    rms=0.00255,
    covariance=np.zeros((3, 3)),
)

# A minimal PANalytical scan of a few points, with the byte order mark the
# instrument software writes.
XRDML = """\ufeff<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.0">
  <sample><id>tiny</id></sample>
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
        <commonCountingTime>1.0</commonCountingTime>
        <counts>{counts}</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
"""

# Primitive cubic LaB6, for the one test that runs GSAS-II.
LAB6_CIF = """data_LaB6
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


def _write_xrdml(path: Path, start: float, end: float, counts: list[int]) -> Path:
    text = XRDML.format(start=start, end=end, counts=" ".join(map(str, counts)))
    path.write_text(text, encoding="utf-8")
    return path


def _fake_install(tmp_path: Path) -> Gsas2Install:
    """An install whose "python" is a fake that plays the driver's part.

    The fake is handed the same arguments as the GSAS-II Python, checks them,
    and writes the job back as its result. A job with ``"fail"`` set makes it
    report an error and exit with status 3.
    """
    script = tmp_path / "fake_driver.py"
    script.write_text(
        "import json, sys\n"
        "flag, driver, job_path = sys.argv[1:]\n"
        "job = json.load(open(job_path))\n"
        "if job.get('fail'):\n"
        "    sys.stderr.write('GSAS-II went wrong\\n')\n"
        "    sys.exit(3)\n"
        "print('fake GSAS-II output')\n"
        "result = {'flag': flag, 'driver': driver, 'job': job}\n"
        "json.dump(result, open(job['result'], 'w'))\n",
        encoding="utf-8",
    )
    if sys.platform == "win32":
        python = tmp_path / "python.cmd"
        python.write_text(f'@"{sys.executable}" "{script}" %*\n', encoding="utf-8")
    else:
        python = tmp_path / "python"
        python.write_text(
            f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n', encoding="utf-8"
        )
        python.chmod(0o755)
    home = tmp_path / "home"
    (home / "GSASII").mkdir(parents=True)
    return Gsas2Install(python=python, home=home)


def _parse_instprm(path: Path) -> tuple[str, dict[str, str]]:
    """The first line and the key:value pairs, split as GSAS-II splits them."""
    first, *lines = path.read_text(encoding="utf-8").splitlines()
    return first, dict(line.split(":", 1) for line in lines)


# find_gsas2


def _clear_environment(monkeypatch, home: Path) -> None:
    monkeypatch.delenv(GSAS2_PYTHON_VARIABLE, raising=False)
    monkeypatch.delenv(GSAS2_HOME_VARIABLE, raising=False)
    monkeypatch.setattr(Path, "home", lambda: home)


def test_find_gsas2_from_environment(tmp_path, monkeypatch) -> None:
    python = tmp_path / "somewhere" / "python.exe"
    python.parent.mkdir()
    python.touch()
    home = tmp_path / "elsewhere"
    (home / "GSASII").mkdir(parents=True)
    _clear_environment(monkeypatch, tmp_path / "nobody")
    monkeypatch.setenv(GSAS2_PYTHON_VARIABLE, str(python))
    monkeypatch.setenv(GSAS2_HOME_VARIABLE, str(home))

    assert find_gsas2() == Gsas2Install(python=python, home=home)


def test_find_gsas2_default_location(tmp_path, monkeypatch) -> None:
    folder = tmp_path / "gsas2main"
    python = folder / ("python.exe" if sys.platform == "win32" else "bin/python")
    python.parent.mkdir(parents=True, exist_ok=True)
    python.touch()
    (folder / "GSAS-II" / "GSASII").mkdir(parents=True)
    _clear_environment(monkeypatch, tmp_path)

    assert find_gsas2() == Gsas2Install(python=python, home=folder / "GSAS-II")


def test_find_gsas2_not_found(tmp_path, monkeypatch) -> None:
    _clear_environment(monkeypatch, tmp_path)

    with pytest.raises(FileNotFoundError) as caught:
        find_gsas2()
    assert GSAS2_PYTHON_VARIABLE in str(caught.value)
    assert GSAS2_HOME_VARIABLE in str(caught.value)


def test_find_gsas2_environment_pointing_nowhere(tmp_path, monkeypatch) -> None:
    home = tmp_path / "home"
    (home / "GSASII").mkdir(parents=True)
    _clear_environment(monkeypatch, tmp_path)
    monkeypatch.setenv(GSAS2_PYTHON_VARIABLE, str(tmp_path / "missing.exe"))
    monkeypatch.setenv(GSAS2_HOME_VARIABLE, str(home))

    with pytest.raises(FileNotFoundError, match="missing.exe"):
        find_gsas2()


# write_instprm


def test_write_instprm_content(tmp_path) -> None:
    path = write_instprm(tmp_path / "cu.instprm", CAGLIOTI, zero=-0.01, shl=0.003)
    first, values = _parse_instprm(path)

    assert first == "#GSAS-II instrument parameter file; do not add/delete items!"
    assert "GSAS-II" in first
    assert list(values) == [
        "Type",
        "Bank",
        "Lam1",
        "Lam2",
        "Zero",
        "I(L2)/I(L1)",
        "Polariz.",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z",
        "SH/L",
        "Azimuth",
    ]
    assert values["Type"] == "PXC"
    assert values["Bank"] == "1"
    assert float(values["Lam1"]) == 1.54056
    assert float(values["Lam2"]) == 1.54439
    assert float(values["I(L2)/I(L1)"]) == 0.5
    assert float(values["Polariz."]) == 0.7
    assert float(values["Zero"]) == -0.01
    assert float(values["SH/L"]) == 0.003
    assert float(values["X"]) == float(values["Y"]) == 0.0


def test_write_instprm_converts_fwhm_squared_to_gaussian_variance(tmp_path) -> None:
    _, values = _parse_instprm(write_instprm(tmp_path / "cu.instprm", CAGLIOTI))
    scale = 1.0e4 / (8.0 * math.log(2.0))

    for key, value in (("U", CAGLIOTI.u), ("V", CAGLIOTI.v), ("W", CAGLIOTI.w)):
        assert float(values[key]) == pytest.approx(value * scale, rel=1e-12)

    # The Gaussian FWHM GSAS-II computes from them, sqrt(8 ln 2) sigma in
    # centidegrees, is the Caglioti FWHM in degrees.
    tan_theta = math.tan(math.radians(30.0))
    variance = (
        float(values["U"]) * tan_theta**2
        + float(values["V"]) * tan_theta
        + float(values["W"])
    )
    fwhm = math.sqrt(8.0 * math.log(2.0) * variance) / 100.0
    assert fwhm == pytest.approx(CAGLIOTI.fwhm(60.0), rel=1e-12)


# run_job


def test_run_job_runs_driver_and_returns_result(tmp_path) -> None:
    install = _fake_install(tmp_path)
    job = {"action": "create", "gpx": "project.gpx"}
    workdir = tmp_path / "work"

    result = run_job(job, workdir, install)

    assert result["flag"] == "-B"
    assert Path(result["driver"]).name == "gsas2_driver.py"
    assert result["job"]["gsas2_home"] == str(install.home)
    assert result["job"]["gpx"] == "project.gpx"
    assert result["job"]["result"] == str(workdir.resolve() / "create_result.json")
    assert "gsas2_home" not in job
    written = json.loads((workdir / "create_job.json").read_text(encoding="utf-8"))
    assert written == result["job"]
    assert "fake GSAS-II output" in (workdir / "create.log").read_text("utf-8")


def test_run_job_writes_xy_fallback_for_xrdml(tmp_path) -> None:
    install = _fake_install(tmp_path)
    counts = [10, 20, 400, 30, 15]
    data_file = _write_xrdml(tmp_path / "scan.xrdml", 20.0, 20.08, counts)

    result = run_job(
        {"action": "create", "data_file": str(data_file)}, tmp_path / "w", install
    )

    fallback = Path(result["job"]["data_fallback"])
    assert fallback == (tmp_path / "w").resolve() / "scan.xy"
    columns = np.loadtxt(fallback)
    np.testing.assert_allclose(columns[:, 0], np.linspace(20.0, 20.08, 5))
    np.testing.assert_array_equal(columns[:, 1], counts)


def test_run_job_raises_with_stderr(tmp_path) -> None:
    install = _fake_install(tmp_path)

    with pytest.raises(Gsas2Error, match="GSAS-II went wrong") as caught:
        run_job({"action": "create", "fail": True}, tmp_path / "work", install)
    assert "exit code 3" in str(caught.value)
    assert "GSAS-II went wrong" in caught.value.stderr


# GSAS-II itself


def test_create_project_with_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 50.0, 1501)
    # LaB6 100, 110 and 111 on a flat background.
    counts = 50.0 + sum(
        height * np.exp(-0.5 * ((two_theta - centre) / 0.03) ** 2)
        for centre, height in ((21.36, 600.0), (30.39, 1000.0), (37.44, 400.0))
    )
    data_file = _write_xrdml(
        tmp_path / "lab6.xrdml", 20.0, 50.0, [round(c) for c in counts]
    )
    cif = tmp_path / "lab6.cif"
    cif.write_text(LAB6_CIF, encoding="utf-8")
    instprm = write_instprm(tmp_path / "cu.instprm", CAGLIOTI)
    gpx = tmp_path / "lab6.gpx"

    result = run_job(
        {
            "action": "create",
            "gpx": str(gpx),
            "data_file": str(data_file),
            "instprm": str(instprm),
            "phases": [{"cif": str(cif), "name": "LaB6"}],
        },
        tmp_path / "work",
        install,
    )

    assert gpx.is_file()
    assert result["source"]["importer"] == "file"
    histogram = result["histogram"]
    assert histogram["n_points"] == 1501
    assert histogram["two_theta_range"] == pytest.approx([20.0, 50.0])
    assert histogram["instrument"]["Type"] == "PXC"
    assert histogram["instrument"]["U"] == pytest.approx(
        CAGLIOTI.u * 1.0e4 / (8.0 * math.log(2.0))
    )
    assert histogram["sample"]["type"] == "Bragg-Brentano"
    (phase,) = result["phases"]
    assert phase["name"] == "LaB6"
    assert phase["cell"]["length_a"] == pytest.approx(4.15683)
    assert phase["histograms"] == [histogram["name"]]
