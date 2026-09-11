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
    build_refine_job,
    find_gsas2,
    gsas2_fwhm,
    run_job,
    standard_stages,
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


# standard_stages and build_refine_job


def test_standard_stages() -> None:
    stages = standard_stages()

    assert [stage["name"] for stage in stages] == [
        "background and scale",
        "zero",
        "cell",
        "U V W",
        "X Y",
        "SH/L",
    ]
    assert stages[0] == {
        "name": "background and scale",
        "background": {"type": "chebyschev-1", "terms": 6},
        "scale": True,
    }
    assert stages[1] == {"name": "zero", "zero": True}
    assert stages[2] == {"name": "cell", "cell": True}
    assert stages[3]["instrument"] == ["U", "V", "W"]
    assert stages[4]["instrument"] == ["X", "Y"]
    assert stages[5]["instrument"] == ["SH/L"]
    assert standard_stages(background_terms=4)[0]["background"]["terms"] == 4

    stages[0]["scale"] = False
    stages.pop()
    assert standard_stages()[0]["scale"] is True
    assert len(standard_stages()) == 6


def test_build_refine_job_creating_the_project(tmp_path) -> None:
    stages = [stage for stage in standard_stages() if stage["name"] != "cell"]
    job = build_refine_job(
        tmp_path / "lab6.gpx",
        stages,
        data_file=tmp_path / "scan.xrdml",
        instprm=tmp_path / "cu.instprm",
        phases=[
            {
                "cif": tmp_path / "lab6.cif",
                "name": "LaB6",
                "cell": [4.156826] * 3 + [90] * 3,
            }
        ],
        limits=(10, 98),
    )

    assert job == {
        "action": "refine",
        "gpx": str(tmp_path / "lab6.gpx"),
        "stages": stages,
        "cycles": 10,
        "export_prefix": str(tmp_path / "lab6"),
        "limits": [10.0, 98.0],
        "data_file": str(tmp_path / "scan.xrdml"),
        "instprm": str(tmp_path / "cu.instprm"),
        "phases": [
            {
                "cif": str(tmp_path / "lab6.cif"),
                "name": "LaB6",
                "cell": [4.156826, 4.156826, 4.156826, 90.0, 90.0, 90.0],
            }
        ],
    }
    job["stages"][0]["scale"] = False
    assert stages[0]["scale"] is True, "the stages are copied"
    json.dumps(job)


def test_build_refine_job_on_an_existing_project(tmp_path) -> None:
    job = build_refine_job(
        tmp_path / "lab6.gpx",
        [{"zero": True}],
        cycles=4,
        broadening={"*": {"size": 10, "mustrain": 0, "lgmix": 0}},
        export_prefix=tmp_path / "out" / "run1",
    )

    assert job == {
        "action": "refine",
        "gpx": str(tmp_path / "lab6.gpx"),
        "stages": [{"zero": True}],
        "cycles": 4,
        "broadening": {"*": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
        "export_prefix": str(tmp_path / "out" / "run1"),
    }


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"data_file": "scan.xrdml"}, "all of data_file"),
        ({"limits": (98, 10)}, "limits must increase"),
        ({"cycles": 0}, "cycles"),
        ({"broadening": {"LaB6": {"size": 100.0}}}, "size must lie between"),
        ({"stages": [{"cell": True, "celll": True}]}, "unknown keys"),
        (
            {"data_file": "s", "instprm": "i", "phases": [{"name": "LaB6"}]},
            "cif and a name",
        ),
    ],
)
def test_build_refine_job_rejects(arguments, message) -> None:
    arguments = {"gpx": "lab6.gpx", "stages": [{"zero": True}], **arguments}
    with pytest.raises(ValueError, match=message):
        build_refine_job(**arguments)


# gsas2_fwhm


def test_gsas2_fwhm_hand_calculated() -> None:
    # 2theta = 40: tan 20 = 0.3639702, cos 20 = 0.9396926.
    # sigma^2 = 10 tan^2 - 5 tan + 12 = 11.504892 cdeg^2, G = 2.35482 sigma
    # = 7.987282 cdeg; gamma = 2 / cos + 3 tan = 3.220266 cdeg.
    # (G^5 + 2.69269 G^4 L + 2.42843 G^3 L^2 + 4.47163 G^2 L^3
    #  + 0.07842 G L^4 + L^5)^(1/5) = 9.803916 cdeg.
    assert gsas2_fwhm(40.0, 10.0, -5.0, 12.0, 2.0, 3.0) == pytest.approx(
        0.09803916, abs=1e-8
    )
    assert gsas2_fwhm(40.0, 10.0, -5.0, 12.0, 2.0, 3.0, shl=0.02) == gsas2_fwhm(
        40.0, 10.0, -5.0, 12.0, 2.0, 3.0
    )


def test_gsas2_fwhm_gaussian_matches_caglioti(tmp_path) -> None:
    _, values = _parse_instprm(write_instprm(tmp_path / "cu.instprm", CAGLIOTI))
    u, v, w = (float(values[key]) for key in "UVW")
    angles = np.array([20.0, 45.0, 90.0, 120.0])

    widths = gsas2_fwhm(angles, u, v, w, 0.0, 0.0)

    # 2.35482 against sqrt(8 ln 2) = 2.3548200 leaves a part in 1e7.
    np.testing.assert_allclose(widths, CAGLIOTI.fwhm(angles), rtol=1e-6)


def test_gsas2_fwhm_lorentzian_and_floor() -> None:
    # A pure Lorentzian: FWHM = X / cos + Y tan + Z, in centidegrees.
    theta = math.radians(30.0)
    expected = (4.0 / math.cos(theta) + 2.0 * math.tan(theta) + 1.0) / 100.0
    lorentzian = gsas2_fwhm(60.0, 0.0, 0.0, -1.0, 4.0, 2.0, z=1.0)
    # The Gaussian variance is held at 0.001, so it barely shows.
    assert lorentzian == pytest.approx(expected, rel=1e-3)
    assert lorentzian > expected


# GSAS-II itself


def test_create_and_refine_with_gsas2(tmp_path) -> None:
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

    # One short stage on the project just made.
    job = build_refine_job(
        gpx,
        [{"name": "background and scale", "background": {"terms": 3}, "scale": True}],
        limits=(20.5, 49.5),
        cycles=3,
        broadening={"LaB6": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
    )
    refined = run_job(job, tmp_path / "work", install)

    assert refined["completed"], refined["stages"]
    (stage,) = refined["stages"]
    assert stage["n_variables"] == 4
    assert set(stage["parameters"]) == {
        ":0:Scale",
        ":0:Back;0",
        ":0:Back;1",
        ":0:Back;2",
    }
    assert all(p["esd"] > 0.0 for p in stage["parameters"].values())
    assert 0.0 < stage["rwp"] < 100.0
    assert stage["gof"] > 0.0
    assert refined["limits"] == pytest.approx([20.5, 49.5])
    final = refined["final"]
    assert final["phases"][0]["cell"]["length_a"] == pytest.approx(4.15683)
    assert final["instrument"]["U"]["esd"] is None
    assert final["phases"][0]["size"] == {
        "type": "isotropic",
        "value": 10.0,
        "lorentzian_fraction": 0.0,
    }
    assert final["phases"][0]["mustrain"]["value"] == 0.0
    exports = refined["exports"]
    table = np.loadtxt(exports["histogram"], delimiter=",", skiprows=1)
    assert table[0, 0] >= 20.5 and table[-1, 0] <= 49.5
    np.testing.assert_allclose(table[:, 4], table[:, 1] - table[:, 2], atol=1e-3)
    reflections = np.loadtxt(exports["reflections"]["LaB6"], delimiter=",", skiprows=1)
    np.testing.assert_array_equal(
        reflections[:3, :3], [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
    )
    assert reflections[0, 5] == pytest.approx(21.357, abs=0.01)
    first, values = _parse_instprm(Path(exports["instprm"]))
    assert "GSAS-II" in first
    assert float(values["U"]) == pytest.approx(
        CAGLIOTI.u * 1.0e4 / (8.0 * math.log(2.0))
    )
