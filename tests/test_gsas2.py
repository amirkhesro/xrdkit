"""Tests for xrdkit.gsas2.

None of these import GSASII: the job runner is exercised with a fake Python
that stands in for GSAS-II's, and the one test that runs GSAS-II itself is
skipped where it is not installed.
"""

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from xrdkit.broadening import Caglioti
from xrdkit.gsas2 import (
    GSAS2_HOME_VARIABLE,
    GSAS2_PYTHON_VARIABLE,
    LOG_TAIL_LINES,
    Gsas2Error,
    Gsas2Install,
    accepted_stages,
    build_refine_job,
    failure_markdown,
    find_gsas2,
    gsas2_fwhm,
    log_tail,
    run_job,
    stage_status,
    stage_status_table,
    stage_statuses,
    standard_stages,
    structure_edits,
    summary_markdown,
    write_instprm,
)
from xrdkit.pipeline import start_from_result

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
    and writes the job back as its result, or the job's ``fake_result`` if it
    has one. A job with ``"fail"`` set makes it report an error and exit with
    status 3.
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
        "result = job.get('fake_result') or {'flag': flag, 'driver': driver, 'job': job}\n"
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


# A two atom phase as the driver reports it: Ti at the centre of a 4 A cube
# written out in P1, O at the middle of a face.
FINAL_PHASE = {
    "name": "P",
    "cell": {
        "length_a": 4.0,
        "length_b": 4.0,
        "length_c": 4.0,
        "angle_alpha": 90.0,
        "angle_beta": 90.0,
        "angle_gamma": 90.0,
    },
    "atoms": [
        {"label": "Ti", "type": "Ti", "xyz": [0.5, 0.5, 0.5]},
        {"label": "O1", "type": "O", "xyz": [0.5, 0.5, 0.0]},
    ],
    "operators": [{"rotation": np.eye(3).tolist(), "translation": [0.0, 0.0, 0.0]}],
}


def test_run_job_adds_the_bond_lengths_a_job_asks_for(tmp_path) -> None:
    install = _fake_install(tmp_path)
    job = {
        "action": "refine",
        "bonds": {"dmax": 2.5},
        "fake_result": {"final": {"phases": [FINAL_PHASE]}},
    }

    result = run_job(job, tmp_path / "work", install)

    (bond,) = result["bonds"]["P"]
    assert bond["centre"] == "Ti" and bond["target"] == "O1"
    assert bond["distance"] == pytest.approx(2.0) and bond["count"] == 2
    written = json.loads((tmp_path / "work" / "refine_result.json").read_text())
    assert written["bonds"] == result["bonds"]

    del job["bonds"]
    assert "bonds" not in run_job(job, tmp_path / "plain", install)


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
    tmp_path = tmp_path.resolve()
    _touch(tmp_path, "scan.xrdml", "cu.instprm", "lab6.cif")
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


def test_build_refine_job_le_bail_cycles() -> None:
    job = build_refine_job(
        "ttb.gpx",
        [{"le_bail": True, "size": True}],
        le_bail_cycles=3,
        max_passes=20,
        pass_tolerance=0.05,
    )

    assert job["le_bail_cycles"] == 3
    assert job["max_passes"] == 20
    assert job["pass_tolerance"] == 0.05
    plain = build_refine_job("ttb.gpx", [{"zero": True}])
    assert not {"le_bail_cycles", "max_passes", "pass_tolerance"} & set(plain)


def test_build_refine_job_structure_background_and_uiso(tmp_path) -> None:
    tmp_path = tmp_path.resolve()
    _touch(tmp_path, "scan.xrdml", "cu.instprm", "ttb.cif")
    job = build_refine_job(
        tmp_path / "ttb.gpx",
        [{"scale": True}, {"overall_uiso": True}],
        data_file=tmp_path / "scan.xrdml",
        instprm=tmp_path / "cu.instprm",
        phases=[
            {
                "cif": tmp_path / "ttb.cif",
                "name": "TTB",
                "atoms": [
                    {"label": "Sr1", "occupancy": 0.6},
                    {"label": "La1", "type": "La", "copy": "Sr1", "occupancy": "0.15"},
                ],
            }
        ],
        background_start={"coefficients": [1000, -50]},
        overall_uiso_start={"TTB": 0.01},
    )

    assert job["phases"][0]["atoms"] == [
        {"label": "Sr1", "occupancy": 0.6},
        {"label": "La1", "type": "La", "copy": "Sr1", "occupancy": 0.15},
    ]
    assert job["background_start"] == {
        "type": "chebyschev-1",
        "coefficients": [1000.0, -50.0],
    }
    assert job["overall_uiso_start"] == {"TTB": 0.01}
    json.dumps(job)


def test_build_refine_job_scale_sanity_and_bonds() -> None:
    reference = {"TTB": [{"label": "O1", "xyz": [0.3, 0.1, 0.1]}]}
    job = build_refine_job(
        "ttb.gpx",
        [{"scale": True}],
        scale_start=14,
        sanity={"max_shift": 0.02, "reference": reference},
        bonds=True,
    )

    assert job["scale_start"] == 14.0
    assert job["sanity"] == {"max_shift": 0.02, "reference": reference}
    assert job["bonds"] == {}
    job = build_refine_job(
        "ttb.gpx",
        [{"scale": True}],
        sanity={"max_shift": 0.1},
        bonds={"anions": ("O", "F"), "dmax": 2.8},
    )
    assert job["sanity"] == {"max_shift": 0.1}
    assert job["bonds"] == {"anions": ["O", "F"], "dmax": 2.8}
    plain = build_refine_job("ttb.gpx", [{"scale": True}], bonds=False)
    assert not {"scale_start", "sanity", "bonds"} & set(plain)


def test_structure_edits_set_up_a_refined_structure() -> None:
    base = [
        {"label": "Sr1", "type": "Sr", "xyz": [0.0, 0.0, 0.49]},
        {"label": "O1", "type": "O", "xyz": [0.3, 0.1, 0.1]},
    ]
    refined = [
        {
            "label": "Sr1",
            "type": "Sr",
            "xyz": [0.0, 0.0, 0.47],
            "occupancy": 0.5,
            "adp": "I",
            "uiso": 0.02,
        },
        {
            "label": "La1",
            "type": "La",
            "xyz": [0.0, 0.0, 1.47],
            "occupancy": 0.1,
            "adp": "I",
            "uiso": 0.02,
        },
        {
            "label": "O1",
            "type": "O",
            "xyz": [0.31, 0.1, 0.1],
            "occupancy": 1.0,
            "adp": "A",
            "uiso": None,
        },
    ]

    edits = structure_edits(refined, base)

    # La1, which the CIF lacks, is added last, beside Sr1, whose site it
    # shares a cell along; O1, anisotropic, keeps its Uij.
    assert edits == [
        {"label": "Sr1", "occupancy": 0.5, "xyz": [0.0, 0.0, 0.47], "uiso": 0.02},
        {"label": "O1", "occupancy": 1.0, "xyz": [0.31, 0.1, 0.1]},
        {
            "label": "La1",
            "occupancy": 0.1,
            "xyz": [0.0, 0.0, 1.47],
            "uiso": 0.02,
            "type": "La",
            "copy": "Sr1",
        },
    ]
    lone = {**refined[1], "label": "Ce1", "xyz": [0.5, 0.5, 0.5]}
    with pytest.raises(ValueError, match="shares no site"):
        structure_edits([*refined, lone], base)


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"data_file": "scan.xrdml"}, "all of data_file"),
        ({"background_start": {"coefficients": []}}, "at least one"),
        ({"overall_uiso_start": {"*": -0.01}}, "must be positive"),
        (
            {
                "data_file": "s",
                "instprm": "i",
                "phases": [{"cif": "c", "name": "P", "atoms": [{"label": "A"}]}],
            },
            "needs a label",
        ),
        ({"limits": (98, 10)}, "limits must increase"),
        ({"cycles": 0}, "cycles"),
        ({"le_bail_cycles": 0}, "le_bail_cycles"),
        ({"max_passes": 0}, "max_passes"),
        ({"pass_tolerance": 0.0}, "pass_tolerance"),
        ({"broadening": {"LaB6": {"size": 100.0}}}, "size must lie between"),
        ({"stages": [{"cell": True, "celll": True}]}, "unknown keys"),
        ({"scale_start": 0.0}, "scale_start must be positive"),
        ({"sanity": {"max_shift": -0.1}}, "max_shift must be positive"),
        ({"bonds": {"dmax": 0.2}}, "dmax must exceed"),
        ({"bonds": {"cations": ["Sr"]}}, "bonds must be True"),
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


def _touch(folder: Path, *names: str) -> list[Path]:
    """Create empty files ``names`` in ``folder``, returning their paths."""
    paths = [folder / name for name in names]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    return paths


def test_build_refine_job_resolves_relative_paths(tmp_path, monkeypatch) -> None:
    _touch(tmp_path, "data/scan.xrdml", "data/cu.instprm", "cifs/ttb.cif")
    monkeypatch.chdir(tmp_path)

    job = build_refine_job(
        "results/ttb.gpx",
        [{"scale": True}],
        data_file="data/scan.xrdml",
        instprm="data/cu.instprm",
        phases=[{"cif": "cifs/ttb.cif", "name": "TTB"}],
        export_prefix="results/run1",
    )

    root = tmp_path.resolve()
    assert job["gpx"] == str(root / "results" / "ttb.gpx")
    assert job["data_file"] == str(root / "data" / "scan.xrdml")
    assert job["instprm"] == str(root / "data" / "cu.instprm")
    assert job["phases"][0]["cif"] == str(root / "cifs" / "ttb.cif")
    assert job["export_prefix"] == str(root / "results" / "run1")
    for key in ("gpx", "data_file", "instprm", "export_prefix"):
        assert Path(job[key]).is_absolute(), key


@pytest.mark.parametrize("missing", ["cu.instprm", "scan.xrdml", "ttb.cif"])
def test_build_refine_job_names_a_missing_input(tmp_path, missing) -> None:
    names = ["scan.xrdml", "cu.instprm", "ttb.cif"]
    _touch(tmp_path, *(name for name in names if name != missing))

    with pytest.raises(Gsas2Error, match=re.escape(missing)):
        build_refine_job(
            tmp_path / "ttb.gpx",
            [{"scale": True}],
            data_file=tmp_path / "scan.xrdml",
            instprm=tmp_path / "cu.instprm",
            phases=[{"cif": tmp_path / "ttb.cif", "name": "TTB"}],
        )


def test_build_refine_job_keeps_a_dot_in_the_stem(tmp_path) -> None:
    root = tmp_path.resolve()
    bare = build_refine_job(tmp_path / "x0.10", [{"zero": True}])
    assert bare["export_prefix"] == str(root / "x0.10")
    named = build_refine_job(tmp_path / "x0.10.gpx", [{"zero": True}])
    assert named["gpx"] == str(root / "x0.10.gpx")
    assert named["export_prefix"] == str(root / "x0.10")


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
        "esd": None,
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


def test_structure_edits_and_overall_uiso_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 50.0, 1501)
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
    phases = [
        {
            "cif": cif,
            "name": "LaB6",
            "atoms": [
                {"label": "La", "occupancy": 0.9},
                {"label": "Ce", "type": "Ce", "copy": "La", "occupancy": 0.1},
            ],
        }
    ]

    created = run_job(
        {
            **build_refine_job(
                tmp_path / "lab6.gpx",
                [{"scale": True}],
                data_file=data_file,
                instprm=instprm,
                phases=phases,
            ),
            "action": "create",
        },
        tmp_path / "create",
        install,
    )
    atoms = {atom["label"]: atom for atom in created["phases"][0]["atoms"]}
    assert set(atoms) == {"La", "Ce", "B"}
    assert atoms["La"]["occupancy"] == pytest.approx(0.9)
    assert atoms["Ce"]["type"] == "Ce"
    assert atoms["Ce"]["occupancy"] == pytest.approx(0.1)
    assert atoms["Ce"]["xyz"] == atoms["La"]["xyz"]
    assert atoms["Ce"]["multiplicity"] == atoms["La"]["multiplicity"] == 1
    assert atoms["Ce"]["uiso"] == pytest.approx(0.0086)

    job = build_refine_job(
        tmp_path / "lab6.gpx",
        [
            {"name": "scale", "background": {"terms": 3}, "scale": True},
            {"name": "Uiso", "overall_uiso": True},
        ],
        data_file=data_file,
        instprm=instprm,
        phases=phases,
        limits=(20.5, 49.5),
        cycles=5,
        broadening={"LaB6": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
        background_start={"coefficients": [50.0, 0.0, 0.0]},
        overall_uiso_start={"*": 0.005},
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    scale, uiso = result["stages"]
    assert not any("AUiso" in name for name in scale["parameters"])
    refined = [name for name in uiso["parameters"] if "AUiso" in name]
    assert refined, "the one Uiso is refined"
    final = {atom["label"]: atom for atom in result["final"]["phases"][0]["atoms"]}
    values = {atom["uiso"] for atom in final.values()}
    assert len(values) == 1, "every Uiso constrained to one"
    assert values != {0.005}, "and refined from the start"
    assert all(atom["adp"] == "I" for atom in final.values())
    assert all(atom["uiso_esd"] > 0.0 for atom in final.values())
    residuals = uiso["residuals"]
    assert 0.0 < residuals["0:0:Rf^2"] < 100.0


def test_uiso_groups_and_shared_site_coordinates_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 90.0, 3501)
    lines = (
        (21.36, 600.0),
        (30.39, 1000.0),
        (37.44, 400.0),
        (43.51, 250.0),
        (48.96, 500.0),
        (53.99, 300.0),
        (63.22, 150.0),
        (67.55, 350.0),
        (71.75, 200.0),
        (75.84, 150.0),
        (79.83, 120.0),
        (83.75, 100.0),
        (87.62, 90.0),
    )
    counts = 50.0 + sum(
        height * np.exp(-0.5 * ((two_theta - centre) / 0.03) ** 2)
        for centre, height in lines
    )
    data_file = _write_xrdml(
        tmp_path / "lab6.xrdml", 20.0, 90.0, [round(c) for c in counts]
    )
    cif = tmp_path / "lab6.cif"
    cif.write_text(LAB6_CIF, encoding="utf-8")
    instprm = write_instprm(tmp_path / "cu.instprm", CAGLIOTI)
    # B split with C on its site, so that the two must move together.
    phases = [
        {
            "cif": cif,
            "name": "LaB6",
            "atoms": [
                {"label": "B", "occupancy": 0.9},
                {"label": "C", "type": "C", "copy": "B", "occupancy": 0.1},
            ],
        }
    ]
    job = build_refine_job(
        tmp_path / "lab6.gpx",
        [
            {"name": "scale", "background": {"terms": 3}, "scale": True},
            {"name": "overall", "overall_uiso": True},
            {
                # B names its site, and so takes in C, which shares it.
                "name": "groups",
                "overall_uiso": False,
                "uiso_groups": {"LaB6": [["La"], ["B"]]},
            },
            {
                "name": "boron",
                "atom_flags": [{"phase": "LaB6", "labels": ["B", "C"], "flags": "X"}],
            },
        ],
        data_file=data_file,
        instprm=instprm,
        phases=phases,
        limits=(20.5, 89.5),
        cycles=5,
        broadening={"LaB6": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
        background_start={"coefficients": [50.0, 0.0, 0.0]},
        overall_uiso_start={"*": 0.005},
        # The boron stage takes the shared site's Uiso below zero; kept here,
        # as it is its coordinates that are under test.
        on_flagged="accept",
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    assert result["stages"][3]["status"] == "flagged"
    overall, groups, boron = (stage["parameters"] for stage in result["stages"][1:])
    assert sum("AUiso" in name for name in overall) == 1
    # Two groups, La alone and B with C: two independent Uiso.
    assert sum("AUiso" in name for name in groups) == 2
    assert any(name.startswith("0::dA") for name in boron)
    phase = result["final"]["phases"][0]
    atoms = {atom["label"]: atom for atom in phase["atoms"]}
    assert atoms["B"]["uiso"] == pytest.approx(atoms["C"]["uiso"])
    assert atoms["B"]["uiso"] != pytest.approx(atoms["La"]["uiso"])
    assert atoms["B"]["xyz"] == pytest.approx(atoms["C"]["xyz"])
    assert atoms["B"]["xyz"] != pytest.approx([0.5, 0.5, 0.2021], abs=1e-7)
    assert any(esd for esd in atoms["B"]["xyz_esd"])
    assert atoms["La"]["xyz_esd"] == [None, None, None]
    assert len(phase["operators"]) == 48
    # Each stage carries the atoms as it left them: B moves only in the last.
    by_stage = [
        {atom["label"]: atom["xyz"] for atom in stage["atoms"]["LaB6"]}
        for stage in result["stages"]
    ]
    assert by_stage[2]["B"] == pytest.approx([0.5, 0.5, 0.2021])
    assert by_stage[3]["B"] == pytest.approx(atoms["B"]["xyz"])


# A made up P4bm structure with an atom on 4c, x, x + 1/2, z, where the site
# ties y to x, to test the constraints of a shared site against GSAS-II's own.
P4BM_CIF = """data_toy
_cell_length_a  6.0
_cell_length_b  6.0
_cell_length_c  4.0
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 4 b m'
_space_group_IT_number 100
loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_occupancy
   _atom_site_U_iso_or_equiv
Ba  Ba  0.2  0.7  0.5   1.0  0.01
Nb  Nb  0.0  0.0  0.0   1.0  0.01
O   O   0.3  0.1  0.1   1.0  0.01
"""


def test_shared_site_moves_together_under_site_symmetry(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 60.0, 2001)
    counts = 50.0 + sum(
        height * np.exp(-0.5 * ((two_theta - centre) / 0.03) ** 2)
        for centre, height in (
            (21.0, 500.0),
            (29.8, 900.0),
            (36.7, 400.0),
            (42.6, 600.0),
        )
    )
    data_file = _write_xrdml(
        tmp_path / "toy.xrdml", 20.0, 60.0, [round(c) for c in counts]
    )
    cif = tmp_path / "toy.cif"
    cif.write_text(P4BM_CIF, encoding="utf-8")
    job = build_refine_job(
        tmp_path / "toy.gpx",
        [
            {"name": "scale", "background": {"terms": 3}, "scale": True},
            {
                "name": "A site",
                "atom_flags": [{"phase": "toy", "labels": ["Ba", "Sr"], "flags": "X"}],
            },
        ],
        data_file=data_file,
        instprm=write_instprm(tmp_path / "cu.instprm", CAGLIOTI),
        phases=[
            {
                "cif": cif,
                "name": "toy",
                "atoms": [
                    {"label": "Ba", "occupancy": 0.5},
                    {"label": "Sr", "type": "Sr", "copy": "Ba", "occupancy": 0.5},
                ],
            }
        ],
        cycles=3,
        broadening={"toy": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
        # The A site moves far enough to be flagged; kept here, as it is how
        # it moves that is under test.
        on_flagged="accept",
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    # The site's own symmetry, y = x + 1/2 on each atom, makes each atom's x
    # independent there and dependent in the tie, and GSAS-II recasts the
    # equivalences as general constraints, which it notes.
    log = (tmp_path / "refine" / "refine.log").read_text(encoding="utf-8")
    assert "Converting equivalence to constraint" in log
    atoms = {atom["label"]: atom for atom in result["final"]["phases"][0]["atoms"]}
    ba, sr = atoms["Ba"]["xyz"], atoms["Sr"]["xyz"]
    assert ba != pytest.approx([0.2, 0.7, 0.5], abs=1e-7), "the site moved"
    assert sr == pytest.approx(ba)
    # Each keeps the site's y = x + 1/2.
    assert sr[1] - sr[0] == pytest.approx(0.5)
    assert ba[1] - ba[0] == pytest.approx(0.5)


def test_le_bail_with_size_and_mustrain_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 50.0, 1501)
    # LaB6 K alpha 1 and 2 doublets broader than the instrument, with heights
    # unlike the structure's, which a Le Bail fit takes as they come.
    counts = 50.0
    for centre, height in ((21.36, 300.0), (30.39, 1000.0), (37.44, 800.0)):
        sin_theta2 = math.sin(math.radians(centre / 2.0)) * 1.54439 / 1.54056
        alpha2 = 2.0 * math.degrees(math.asin(sin_theta2))
        for line, weight in ((centre, 1.0), (alpha2, 0.5)):
            counts = counts + weight * height * np.exp(
                -0.5 * ((two_theta - line) / 0.05) ** 2
            )
    data_file = _write_xrdml(
        tmp_path / "lab6.xrdml", 20.0, 50.0, [round(c) for c in counts]
    )
    cif = tmp_path / "lab6.cif"
    cif.write_text(LAB6_CIF, encoding="utf-8")
    instprm = write_instprm(tmp_path / "cu.instprm", CAGLIOTI)

    job = build_refine_job(
        tmp_path / "lab6.gpx",
        [
            {
                "name": "background and scale",
                "background": {"terms": 3},
                "scale": True,
                "le_bail": True,
            },
            {"name": "size", "size": True},
            {"name": "mustrain", "mustrain": True},
        ],
        data_file=data_file,
        instprm=instprm,
        phases=[{"cif": cif, "name": "LaB6"}],
        limits=(20.5, 49.5),
        cycles=5,
        broadening={"LaB6": {"size": 1.0, "mustrain": 0.0}},
        le_bail_cycles=3,
        max_passes=4,
    )
    result = run_job(job, tmp_path / "work", install)

    assert result["completed"], result["stages"]
    for stage in result["stages"]:
        assert 1 <= len(stage["passes"]) <= 4
        assert stage["passes"][-1]["gof"] == pytest.approx(stage["gof"])
        assert ":0:Scale" not in {p["parameter"] for p in stage["passes"]}
    first, size, mustrain = result["stages"]
    assert first["le_bail_extraction"]["cycles"] == 3
    assert first["le_bail_extraction"]["rwp"] > 0.0
    assert "le_bail_extraction" not in size
    assert "0:0:Size;i" in size["parameters"]
    assert "0:0:Mustrain;i" not in size["parameters"]
    assert {"0:0:Size;i", "0:0:Mustrain;i"} <= set(mustrain["parameters"])
    # The lines are broader than the instrument, so the size falls from 1 um.
    # (Over three lines the microstrain then trades off against it.)
    assert size["rwp"] < first["rwp"]
    assert 0.001 < size["parameters"]["0:0:Size;i"]["value"] < 1.0
    assert mustrain["rwp"] < size["rwp"]
    phase = result["final"]["phases"][0]
    assert phase["le_bail"]
    assert phase["size"]["esd"] > 0.0
    assert phase["mustrain"]["esd"] > 0.0
    # Le Bail takes the observed heights. With the multiplicities, 8 and 6,
    # and the Lorentz polarisation factor, 3.3 times larger at 100, the
    # heights make F^2(111) / F^2(100) about 6.6, where LaB6 itself gives
    # about 1.6.
    reflections = np.loadtxt(
        result["exports"]["reflections"]["LaB6"], delimiter=",", skiprows=1
    )
    f_obs = {tuple(row[:3].astype(int)): row[7] for row in reflections}
    assert f_obs[(1, 1, 1)] / f_obs[(1, 0, 0)] == pytest.approx(6.6, rel=0.2)

    # Every stage records every value, refined or held, and the job's start.
    assert set(result["start_model"]) == {"LaB6"}
    values = size["values"]
    (lab6,) = values["phases"]
    assert lab6["name"] == "LaB6"
    assert lab6["size"]["held"] is False
    assert lab6["size"]["value"] == pytest.approx(
        size["parameters"]["0:0:Size;i"]["value"]
    )
    assert lab6["mustrain"] == {**lab6["mustrain"], "value": 0.0, "esd": None}
    assert lab6["mustrain"]["held"] is True
    assert lab6["cell_held"] is True and lab6["cell_esd"] is None
    assert lab6["cell"]["length_a"] == pytest.approx(4.15683, abs=1e-4)
    assert lab6["phase_fraction"]["held"] is True
    assert values["instrument"]["Zero"]["held"] is True
    assert values["sample"]["Scale"]["held"] is False
    assert len(values["background"]["coefficients"]) == 3
    assert not any(c["held"] for c in values["background"]["coefficients"])
    # So a later job can start from the size stage, before the microstrain.
    saved = tmp_path / "lab6_result.json"
    saved.write_text(json.dumps(result), encoding="utf-8")
    start = start_from_result(saved, stage="size")
    assert start.microstrain == 0.0
    assert start.size == pytest.approx(lab6["size"]["value"])
    assert start.phases[0].name == "LaB6"
    assert start.start_model == result["start_model"]


# A made up P4bm structure with the two A sites of the tungsten bronze, Sr1
# on 2a (0, 0, z) and Ba2 on 4c (x, x + 1/2, z), Nb1 on 2b (1/2, 0, z) and
# an oxygen on a general position, for the coordinate, origin, occupancy and
# sanity tests.
TOY_TTB_CIF = """data_toy
_cell_length_a  6.0
_cell_length_b  6.0
_cell_length_c  4.0
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 4 b m'
_space_group_IT_number 100
loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_occupancy
   _atom_site_U_iso_or_equiv
Sr1  Sr  0.0  0.0  0.5   0.6  0.01
Ba2  Ba  0.2  0.7  0.5   0.6  0.01
Nb1  Nb  0.5  0.0  0.0   1.0  0.01
O1   O   0.3  0.1  0.1   1.0  0.01
"""


def _toy_ttb(tmp_path: Path) -> dict:
    """The build_refine_job arguments that create the toy tungsten bronze
    project: a pattern with a line at every hkl up to 3 of its cell, the
    CIF and an instrument file."""
    two_theta = np.linspace(20.0, 60.0, 2001)
    counts = np.full(two_theta.size, 50.0)
    for h in range(4):
        for k in range(h + 1):
            for l in range(4):
                if h + k + l == 0:
                    continue
                d = 1.0 / math.sqrt((h * h + k * k) / 36.0 + l * l / 16.0)
                sine = 1.5406 / (2.0 * d)
                if sine >= 1.0:
                    continue
                centre = 2.0 * math.degrees(math.asin(sine))
                height = 2000.0 / (1 + h + k + l)
                counts += height * np.exp(-0.5 * ((two_theta - centre) / 0.03) ** 2)
    cif = tmp_path / "toy.cif"
    cif.write_text(TOY_TTB_CIF, encoding="utf-8")
    return {
        "data_file": _write_xrdml(
            tmp_path / "toy.xrdml", 20.0, 60.0, [round(c) for c in counts]
        ),
        "instprm": write_instprm(tmp_path / "cu.instprm", CAGLIOTI),
        "phases": [{"cif": cif, "name": "toy"}],
        "limits": (20.5, 59.5),
        "cycles": 3,
        "broadening": {"toy": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
    }


def _atoms_by_label(result: dict) -> dict:
    return {atom["label"]: atom for atom in result["final"]["phases"][0]["atoms"]}


def test_coordinates_by_site_and_origin_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    inputs = _toy_ttb(tmp_path)
    scale = {"name": "scale", "background": {"terms": 3}, "scale": True}
    sites = {
        "name": "sites",
        "coordinates": {"toy": {"Sr1": "all", "Ba2": "z", "Nb1": "all", "O1": "xy"}},
        "origin": {"toy": {"site": "Nb1", "axis": "z"}},
    }
    # The sites move far enough to be flagged; kept here, as it is which
    # coordinates move that is under test.
    job = build_refine_job(
        tmp_path / "toy.gpx", [scale, sites], on_flagged="accept", **inputs
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    refined = result["stages"][1]["coordinates"]["toy"]
    assert refined == {
        "Sr1": {"site": "Sr1", "refined": "z", "held": "", "origin": False},
        "Ba2": {"site": "Ba2", "refined": "z", "held": "x", "origin": False},
        "Nb1": {"site": "Nb1", "refined": "", "held": "z", "origin": True},
        "O1": {"site": "O1", "refined": "xy", "held": "z", "origin": False},
    }
    atoms = _atoms_by_label(result)
    # Held: Ba2's x, and y with it, O1's z and Nb1's z, the origin.
    assert atoms["Ba2"]["xyz"][:2] == pytest.approx([0.2, 0.7], abs=1e-9)
    assert atoms["O1"]["xyz"][2] == pytest.approx(0.1, abs=1e-9)
    assert atoms["Nb1"]["xyz"] == pytest.approx([0.5, 0.0, 0.0], abs=1e-9)
    assert atoms["O1"]["xyz_esd"][2] is None
    # Refined: the rest.
    assert atoms["Ba2"]["xyz"][2] != pytest.approx(0.5, abs=1e-7)
    assert atoms["O1"]["xyz"][:2] != pytest.approx([0.3, 0.1], abs=1e-7)
    assert atoms["Sr1"]["xyz_esd"][2] > 0.0

    # Without the origin every site would move along the polar c axis.
    every = {"toy": {label: "all" for label in ("Sr1", "Ba2", "Nb1", "O1")}}
    job = build_refine_job(
        tmp_path / "floating.gpx",
        [scale, {"name": "every", "coordinates": every}],
        on_flagged="accept",
        **inputs,
    )
    result = run_job(job, tmp_path / "floating", install)

    assert not result["completed"]
    assert "polar along z" in result["stages"][1]["error"]


def test_occupancy_exchange_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    inputs = _toy_ttb(tmp_path)
    # Ba put on the 2a site and Sr on the 4c site, for the two to trade.
    inputs["phases"][0]["atoms"] = [
        {"label": "Ba1", "type": "Ba", "copy": "Sr1", "occupancy": 0.1},
        {"label": "Sr2", "type": "Sr", "copy": "Ba2", "occupancy": 0.2},
    ]
    exchange = {"phase": "toy", "sites": ["Sr1", "Ba2"], "elements": ["Sr", "Ba"]}
    job = build_refine_job(
        tmp_path / "toy.gpx",
        [
            {"name": "scale", "background": {"terms": 3}, "scale": True},
            {"name": "exchange", "occupancies": [exchange]},
        ],
        # The exchange takes an occupancy past 1; kept here, as it is the
        # constraints that are under test.
        on_flagged="accept",
        **inputs,
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    stage = result["stages"][1]
    # GSAS-II refines new variables of its own in place of the constrained
    # occupancies, one for each element.
    assert sum(name.startswith("::constr") for name in stage["parameters"]) == 2
    (constraints,) = stage["occupancy_constraints"].values()
    assert [(c["element"], c["atoms"]) for c in constraints] == [
        ("Sr", ["Sr1", "Sr2"]),
        ("Ba", ["Ba1", "Ba2"]),
    ]
    assert [c["total"] for c in constraints] == pytest.approx(
        [2 * 0.6 + 4 * 0.2, 2 * 0.1 + 4 * 0.6]
    )
    atoms = _atoms_by_label(result)

    def content(*labels):
        return sum(atoms[a]["multiplicity"] * atoms[a]["occupancy"] for a in labels)

    # Each element's content over the two sites held, its split moved.
    assert content("Sr1", "Sr2") == pytest.approx(2.0, abs=1e-6)
    assert content("Ba1", "Ba2") == pytest.approx(2.6, abs=1e-6)
    assert atoms["Sr1"]["occupancy"] != pytest.approx(0.6, abs=1e-6)
    assert atoms["Sr1"]["occupancy_esd"] > 0.0


def test_every_stage_rejected_still_exports_the_start_in_gsas2(tmp_path) -> None:
    """A run that keeps no stage leaves the model it started from, computed
    but not refined, with its exports and its bonds."""
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    job = build_refine_job(
        tmp_path / "toy.gpx",
        [{"name": "oxygen", "coordinates": {"toy": {"O1": "all"}}}],
        # Any move at all is flagged, so the only stage is rejected.
        sanity={"max_shift": 1e-9},
        bonds=True,
        export_prefix=tmp_path / "toy",
        **_toy_ttb(tmp_path),
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    (oxygen,) = result["stages"]
    assert oxygen["status"] == "rejected"
    assert result["rejected"] == ["oxygen"]
    assert result["final_from"] == "job start"
    assert "compute_error" not in result
    # The model left behind is the CIF's, not the one the stage refined to.
    atoms = _atoms_by_label(result)
    assert atoms["O1"]["xyz"] == pytest.approx([0.3, 0.1, 0.1], abs=1e-9)
    assert atoms["O1"]["xyz_esd"] == [None, None, None]
    moved = next(a for a in oxygen["atoms"]["toy"] if a["label"] == "O1")
    assert moved["xyz"] != pytest.approx([0.3, 0.1, 0.1], abs=1e-7)
    # It is exported like any other, with a pattern really computed from it,
    # and its bonds are worked out from it.
    assert Path(result["exports"]["instprm"]).is_file()
    pattern = np.genfromtxt(result["exports"]["histogram"], delimiter=",", names=True)
    assert pattern["calculated"].max() > 0.1 * pattern["observed"].max()
    bonds = result["bonds"]["toy"]
    assert bonds and {bond["target"] for bond in bonds} == {"O1"}
    assert {bond["centre"] for bond in bonds} <= {"Sr1", "Ba2", "Nb1"}
    assert all(0.5 <= bond["distance"] <= 3.0 for bond in bonds)
    written = json.loads((tmp_path / "refine" / "refine_result.json").read_text())
    assert written["bonds"] == result["bonds"]


def test_bond_lengths_in_gsas2(tmp_path) -> None:
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    counts = 50.0 + sum(
        height * np.exp(-0.5 * ((np.linspace(20.0, 50.0, 1501) - c) / 0.03) ** 2)
        for c, height in ((21.36, 600.0), (30.39, 1000.0), (37.44, 400.0))
    )
    cif = tmp_path / "lab6.cif"
    cif.write_text(LAB6_CIF, encoding="utf-8")
    job = build_refine_job(
        tmp_path / "lab6.gpx",
        [{"name": "scale", "background": {"terms": 3}, "scale": True}],
        data_file=_write_xrdml(
            tmp_path / "lab6.xrdml", 20.0, 50.0, [round(c) for c in counts]
        ),
        instprm=write_instprm(tmp_path / "cu.instprm", CAGLIOTI),
        phases=[{"cif": cif, "name": "LaB6"}],
        cycles=3,
        bonds={"anions": ["B"], "dmax": 3.1},
    )
    result = run_job(job, tmp_path / "refine", install)

    # La sits among 24 B, each at a sqrt(1/2 + z^2) with B at (1/2, 1/2, z).
    (bond,) = result["bonds"]["LaB6"]
    phase = result["final"]["phases"][0]
    boron = next(atom for atom in phase["atoms"] if atom["label"] == "B")
    z = min(abs(value - round(value)) for value in boron["xyz"])
    a = phase["cell"]["length_a"]
    assert (bond["centre"], bond["target"], bond["count"]) == ("La", "B", 24)
    assert bond["distance"] == pytest.approx(a * math.sqrt(0.5 + z * z))


# Stage statuses and write ups of a sequence of refinements


def _stage(name, status, **rest):
    return {"name": name, "status": status, "rwp": 4.25, "gof": 1.47, **rest}


REJECTED_RUN = {
    "completed": True,
    "on_flagged": "reject",
    "on_unsettled": "reject",
    "rejected": ["oxygen"],
    "stages": [
        _stage("profile", "clean", rwp=4.20, gof=1.45),
        _stage(
            "A site",
            "unsettled",
            passes=[{}] * 60,
            largest_remaining_move={"parameter": "0::dAz:4", "shift_over_esd": 0.42},
        ),
        _stage(
            "oxygen",
            "rejected",
            rejected_because=["sanity check: O3 z moved +0.0890"],
            sanity=[{"kind": "coordinate_shift", "message": "O3 z moved +0.0890"}],
            undetermined=[{"message": "O3 x 0.21000, esd 0.00400 more than its shift"}],
        ),
        _stage("after", "clean", rwp=4.18, gof=1.44),
    ],
}


def test_stage_status_reads_the_status_the_driver_recorded() -> None:
    assert stage_status({"status": "rejected"}) == "rejected"
    # Results written before the driver recorded one.
    assert stage_status({"name": "cell"}) == "clean"
    assert stage_status({"name": "cell", "error": "boom"}) == "failed"


def test_accepted_stages_leaves_out_the_rejected_and_failed() -> None:
    result = {
        "stages": [
            _stage("clean", "clean"),
            _stage("flagged", "flagged"),
            _stage("unsettled", "unsettled"),
            _stage("rejected", "rejected"),
            _stage("failed", "failed"),
        ]
    }

    assert [stage["name"] for stage in accepted_stages(result)] == [
        "clean",
        "flagged",
        "unsettled",
    ]


def test_stage_statuses_say_why_each_stage_is_not_clean() -> None:
    rows = stage_statuses(REJECTED_RUN)

    assert [(row["name"], row["status"]) for row in rows] == [
        ("profile", "clean"),
        ("A site", "unsettled"),
        ("oxygen", "rejected"),
        ("after", "clean"),
    ]
    assert rows[0]["reason"] == ""
    assert rows[1]["reason"] == (
        "not settled in 60 passes; largest remaining move 0.42 esd (0::dAz:4)"
    )
    assert rows[1]["passes"] == 60
    assert rows[2]["reason"] == "sanity check: O3 z moved +0.0890"
    assert rows[2]["undetermined"] == ["O3 x 0.21000, esd 0.00400 more than its shift"]
    assert rows[3]["rwp"] == 4.18


def test_stage_statuses_of_flagged_and_failed_stages() -> None:
    result = {
        "stages": [
            _stage("flagged", "flagged", sanity=[{"message": "La Uiso -0.0040"}]),
            _stage("failed", "failed", error="RefinementError: Ouch #4", rwp=None),
        ]
    }

    rows = stage_statuses(result)

    assert [row["reason"] for row in rows] == [
        "La Uiso -0.0040",
        "RefinementError: Ouch #4",
    ]
    assert rows[1]["passes"] is None


def test_stage_status_table_is_markdown_with_a_row_per_stage() -> None:
    lines = stage_status_table(REJECTED_RUN)

    assert lines[0] == "| Stage | Status | Passes | Rwp (%) | GOF | Why |"
    assert lines[1] == "| --- | --- | --- | --- | --- | --- |"
    assert lines[2] == "| profile | clean | — | 4.200 | 1.450 | — |"
    assert lines[4].startswith("| oxygen | rejected | — | 4.250 | 1.470 | sanity")
    assert len(lines) == 6


def test_stage_status_table_escapes_pipes_and_missing_figures() -> None:
    result = {"stages": [_stage("a | b", "failed", error="x | y", rwp=None, gof=None)]}

    (row,) = stage_status_table(result)[2:]

    assert row == "| a \\| b | failed | — | — | — | x \\| y |"


def test_log_tail_reads_the_last_lines(tmp_path) -> None:
    log = tmp_path / "refine.log"
    log.write_text("\n".join(f"line {i}" for i in range(100)), encoding="utf-8")

    assert log_tail(log, 3) == "line 97\nline 98\nline 99"
    assert log_tail(log).count("\n") == 39
    assert log_tail(tmp_path / "missing.log") == ""


def test_failure_markdown_carries_the_error_the_stages_and_the_log() -> None:
    text = failure_markdown(
        "x = 0.10 powder, coordinates",
        REJECTED_RUN,
        "RefinementError: Ouch #4 stuck",
        "the last of the log",
        intro=["The run stopped."],
    )

    assert text.startswith("# x = 0.10 powder, coordinates\n\nThe run stopped.\n")
    assert "## Error\n\n```\nRefinementError: Ouch #4 stuck\n```" in text
    assert "| oxygen | rejected |" in text
    assert f"## GSAS-II log, last {LOG_TAIL_LINES} lines" in text
    assert "```\nthe last of the log\n```" in text
    assert text.endswith("\n")


def test_failure_markdown_without_a_result_or_a_log() -> None:
    text = failure_markdown("x = 0.12 powder, lebail", None, "Gsas2Error: no", "")

    assert "No stage was refined." in text
    assert "(no log)" in text


def test_summary_markdown_tables_every_mode_and_its_stages() -> None:
    entries = [
        {
            "name": "lebail",
            "run": "2026-09-12 09:00",
            "result": {
                "completed": True,
                "stages": [_stage("size", "clean", rwp=3.93, gof=1.33)],
            },
            "values": {"a (A)": "12.4740(3)", "c (A)": "3.93181(10)"},
        },
        {
            "name": "coordinates",
            "run": "2026-09-12 09:20",
            "result": REJECTED_RUN,
            "values": {"a (A)": "12.4738(3)"},
        },
        {"name": "occupancies", "run": None, "result": None},
    ]

    text = summary_markdown(
        "x = 0.10 powder", entries, columns=["a (A)", "c (A)"], intro=["Every mode."]
    )
    lines = text.splitlines()

    assert lines[0] == "# x = 0.10 powder"
    assert lines[2] == "Every mode."
    header = "| Refinement | Run | Outcome | Rwp (%) | GOF | a (A) | c (A) | Stages |"
    assert lines[4] == header
    assert lines[6] == (
        "| lebail | 2026-09-12 09:00 | completed | 3.930 | 1.330 | 12.4740(3) | "
        "3.93181(10) | size: clean |"
    )
    # The final Rwp and GOF are the last stage kept, not the rejected one.
    assert lines[7].startswith(
        "| coordinates | 2026-09-12 09:20 | completed | 4.180 | 1.440 | "
        "12.4738(3) | — | profile: clean; A site: unsettled; oxygen: rejected; "
        "after: clean |"
    )
    assert lines[8] == "| occupancies | — | not run | — | — | — | — | — |"
    assert "## Stages not clean, and undetermined parameters" in text
    assert "- coordinates, A site: unsettled, not settled in 60 passes" in text
    assert "- coordinates, oxygen: rejected, sanity check: O3 z moved +0.0890" in text
    assert (
        "- coordinates, oxygen (rejected): undetermined O3 x 0.21000, esd 0.00400 "
        "more than its shift" in text
    )


def test_summary_markdown_of_a_failed_mode_and_a_clean_one() -> None:
    failed = {
        "completed": False,
        "error": "Gsas2Error: the driver wrote no result",
        "stages": [
            _stage("profile", "clean"),
            _stage("cell", "failed", rwp=None, error="RefinementError: Ouch #4"),
        ],
    }
    entries = [
        {"name": "fixed-atoms", "run": "now", "result": failed},
        {
            "name": "lebail",
            "run": "now",
            "result": {"completed": True, "stages": [_stage("size", "clean")]},
        },
    ]

    text = summary_markdown("x = 0.12 powder", entries)

    assert "| fixed-atoms | now | failed | 4.250 | 1.470 | profile: clean; " in text
    assert "- fixed-atoms: Gsas2Error: the driver wrote no result" in text
    assert "- fixed-atoms, cell: failed, RefinementError: Ouch #4" in text
    assert "- lebail" not in text


def test_summary_markdown_with_nothing_to_report() -> None:
    entries = [
        {
            "name": "lebail",
            "run": "now",
            "result": {"completed": True, "stages": [_stage("size", "clean")]},
        }
    ]

    text = summary_markdown("x = 0.10 powder", entries)

    assert "Every stage is clean, and no parameter is left undetermined." in text


def test_summary_markdown_of_a_run_with_no_stage_kept() -> None:
    result = {"completed": True, "stages": [_stage("scale", "rejected", rwp=None)]}

    text = summary_markdown(
        "x = 0.10 powder", [{"name": "lebail", "run": "now", "result": result}]
    )

    assert "| lebail | now | completed | — | — | scale: rejected |" in text
    assert (
        "- lebail: no stage was kept, so the model left behind is the one the "
        "refinement started from" in text
    )


def test_build_refine_job_carries_the_policies() -> None:
    job = build_refine_job(
        "out.gpx", [{"scale": True}], on_flagged="accept", on_unsettled="reject"
    )

    assert job["on_flagged"] == "accept"
    assert job["on_unsettled"] == "reject"
    # Left out, so that the driver's defaults stand.
    assert "on_flagged" not in build_refine_job("out.gpx", [{"scale": True}])
    with pytest.raises(ValueError, match="on_unsettled must be one of"):
        build_refine_job("out.gpx", [{"scale": True}], on_unsettled="roll back")


def test_a_flagged_stage_is_rolled_back_in_gsas2(tmp_path) -> None:
    """The rule the x = 0.10 oxygen stage runs into, against GSAS-II: a
    stage whose sanity check raises a flag is rejected, the project goes
    back to the last stage kept, and the run carries on from there."""
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    stages = [
        {"name": "scale", "background": {"terms": 3}, "scale": True},
        {"name": "oxygen", "coordinates": {"toy": {"O1": "all"}}},
        {"name": "cell", "cell": True},
    ]
    inputs = _toy_ttb(tmp_path)
    job = build_refine_job(
        tmp_path / "toy.gpx",
        stages,
        # Any move at all is flagged, so the oxygen stage is rejected.
        sanity={"max_shift": 1e-9},
        export_prefix=tmp_path / "toy",
        **inputs,
    )
    result = run_job(job, tmp_path / "refine", install)

    assert result["completed"], result["stages"]
    assert result["on_flagged"] == "reject"
    assert result["rejected"] == ["oxygen"]
    scale, oxygen, cell = result["stages"]
    assert [stage["status"] for stage in result["stages"]] == [
        "clean",
        "rejected",
        "clean",
    ]
    # The rejected stage is recorded with why, and with the atoms it left.
    assert oxygen["rejected_because"] == [
        f"sanity check: {flag['message']}" for flag in oxygen["sanity"]
    ]
    assert all(
        reason.startswith("sanity check: O1 ") for reason in oxygen["rejected_because"]
    )
    moved = next(a for a in oxygen["atoms"]["toy"] if a["label"] == "O1")
    assert moved["xyz"] != pytest.approx([0.3, 0.1, 0.1], abs=1e-7)
    assert "parameters" in oxygen and oxygen["rwp"] > 0.0

    # The cell stage starts from the atoms of the scale stage, the last kept.
    for stage in (scale, cell):
        clean = next(a for a in stage["atoms"]["toy"] if a["label"] == "O1")
        assert clean["xyz"] == pytest.approx([0.3, 0.1, 0.1], abs=1e-9)
    # What the rejected stage alone refined is held: O1's coordinates are no
    # longer among the parameters, while the cell's are.
    assert not any(name.startswith("0::dA") for name in cell["parameters"])
    assert any(name.startswith("0::A") for name in cell["parameters"])
    assert "coordinates" not in cell

    # The final model and the exports are the cell stage's, not the rejected
    # oxygen stage's.
    atoms = _atoms_by_label(result)
    assert atoms["O1"]["xyz"] == pytest.approx([0.3, 0.1, 0.1], abs=1e-9)
    assert atoms["O1"]["xyz_esd"] == [None, None, None]
    assert result["final"]["phases"][0]["atoms"] == cell["atoms"]["toy"]
    assert Path(result["exports"]["histogram"]).is_file()

    # Kept instead, the stage stands and its move is carried on from.
    kept = build_refine_job(
        tmp_path / "kept.gpx",
        stages,
        sanity={"max_shift": 1e-9},
        on_flagged="accept",
        **inputs,
    )
    kept_result = run_job(kept, tmp_path / "kept", install)

    assert kept_result["rejected"] == []
    assert kept_result["stages"][1]["status"] == "flagged"
    kept_oxygen = next(
        a for a in kept_result["stages"][1]["atoms"]["toy"] if a["label"] == "O1"
    )
    assert kept_oxygen["xyz"] == pytest.approx(moved["xyz"])
    # Its coordinates stay free in the cell stage after it, so O1 goes on
    # moving rather than going back to the CIF.
    final = _atoms_by_label(kept_result)["O1"]
    assert final["xyz"] != pytest.approx([0.3, 0.1, 0.1], abs=1e-3)
    assert all(esd for esd in final["xyz_esd"])


def test_write_instprm_carries_the_goniometer_radius(tmp_path) -> None:
    # GSAS-II's own writer spells the key with the space and its reader
    # strips every space before comparing, so this form is read back.
    path = write_instprm(tmp_path / "cu.instprm", CAGLIOTI, radius_mm=145.0)
    _, values = _parse_instprm(path)

    assert list(values)[-1] == "Gonio. radius"
    assert float(values["Gonio. radius"]) == 145.0


def test_write_instprm_without_a_radius_emits_no_line(tmp_path) -> None:
    _, values = _parse_instprm(write_instprm(tmp_path / "cu.instprm", CAGLIOTI))

    assert "Gonio. radius" not in values


def test_create_from_an_xy_takes_the_radius_from_the_instprm(tmp_path) -> None:
    # The .xy importer carries no geometry, so without the instprm line the
    # histogram would keep GSAS-II's default of 200 mm.
    try:
        install = find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    two_theta = np.linspace(20.0, 50.0, 1501)
    counts = 50.0 + sum(
        height * np.exp(-0.5 * ((two_theta - centre) / 0.03) ** 2)
        for centre, height in ((21.36, 600.0), (30.39, 1000.0), (37.44, 400.0))
    )
    data_file = tmp_path / "lab6.xy"
    data_file.write_text(
        "\n".join(f"{t:.5f} {c:.3f}" for t, c in zip(two_theta, counts)) + "\n",
        encoding="utf-8",
    )
    cif = tmp_path / "lab6.cif"
    cif.write_text(LAB6_CIF, encoding="utf-8")
    instprm = write_instprm(tmp_path / "cu.instprm", CAGLIOTI, radius_mm=145.0)

    result = run_job(
        {
            "action": "create",
            "gpx": str(tmp_path / "lab6.gpx"),
            "data_file": str(data_file),
            "instprm": str(instprm),
            "phases": [{"cif": str(cif), "name": "LaB6"}],
        },
        tmp_path / "work",
        install,
    )

    sample = result["histogram"]["sample"]
    assert sample["gonio_radius"] == pytest.approx(145.0)
    assert sample["gonio_radius"] != pytest.approx(200.0)
    # The Lam1 key of the instprm is what tells GSAS-II the geometry.
    assert sample["type"] == "Bragg-Brentano"
