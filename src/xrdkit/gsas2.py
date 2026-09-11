"""GSAS-II jobs, run under GSAS-II's own Python.

GSAS-II brings its own Python and its own compiled extensions, and cannot be
imported here. Instead a job, a dict naming an action and its inputs, is
written as JSON, and ``gsas2_driver.py`` runs it under the GSAS-II Python
through GSASIIscriptable and writes its result as JSON for :func:`run_job` to
read back. The driver imports nothing from xrdkit.

:func:`find_gsas2` locates the installation, and :func:`write_instprm` writes
the instrument parameter file a powder histogram is read with. This module
imports only the standard library; the one exception, the xrdkit reader used
to write a two column copy of an ``.xrdml`` scan, is imported when needed.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from xrdkit.broadening import Caglioti

__all__ = [
    "GSAS2_HOME_VARIABLE",
    "GSAS2_PYTHON_VARIABLE",
    "Gsas2Error",
    "Gsas2Install",
    "find_gsas2",
    "run_job",
    "write_instprm",
]

GSAS2_PYTHON_VARIABLE = "XRDKIT_GSAS2_PYTHON"
GSAS2_HOME_VARIABLE = "XRDKIT_GSAS2_HOME"

# The default installation, under the home folder: the GSAS-II Python sits at
# its top and the folder holding the GSASII package below it.
DEFAULT_FOLDER = "gsas2main"
DEFAULT_HOME = "GSAS-II"

DRIVER = Path(__file__).with_name("gsas2_driver.py")

# First line of a GSAS-II instrument parameter file, which its reader checks.
INSTPRM_HEADER = "#GSAS-II instrument parameter file; do not add/delete items!"

# Degrees squared to centidegrees squared.
CENTIDEGREES_SQUARED = 1.0e4

# FWHM^2 over the variance of a Gaussian.
EIGHT_LN2 = 8.0 * math.log(2.0)


class Gsas2Error(RuntimeError):
    """A GSAS-II job that exited with an error; ``stderr`` holds its report."""

    def __init__(self, message: str, stdout: str = "", stderr: str = "") -> None:
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr


@dataclass(frozen=True)
class Gsas2Install:
    """A GSAS-II installation.

    ``python`` is the interpreter GSAS-II runs under and ``home`` the folder
    that contains the ``GSASII`` package, which the driver puts on
    ``sys.path``.
    """

    python: Path
    home: Path


def _default_python(folder: Path) -> Path:
    if sys.platform == "win32":
        return folder / "python.exe"
    return folder / "bin" / "python"


def find_gsas2() -> Gsas2Install:
    """Locate GSAS-II.

    ``XRDKIT_GSAS2_PYTHON`` names the GSAS-II Python and ``XRDKIT_GSAS2_HOME``
    the folder containing the ``GSASII`` package. Either left unset falls back
    to the default installation, ``~/gsas2main``, whose Python is
    ``python.exe`` on Windows and ``bin/python`` elsewhere, and whose package
    folder is ``GSAS-II``.

    Raises
    ------
    FileNotFoundError
        If the Python is not a file or the home folder holds no ``GSASII``
        package, naming both environment variables.
    """
    folder = Path.home() / DEFAULT_FOLDER
    python = Path(os.environ.get(GSAS2_PYTHON_VARIABLE) or _default_python(folder))
    home = Path(os.environ.get(GSAS2_HOME_VARIABLE) or folder / DEFAULT_HOME)

    problems = []
    if not python.is_file():
        problems.append(f"no GSAS-II Python at {python}")
    if not (home / "GSASII").is_dir():
        problems.append(f"no GSASII package in {home}")
    if problems:
        raise FileNotFoundError(
            "GSAS-II not found: "
            + "; ".join(problems)
            + f". Set {GSAS2_PYTHON_VARIABLE} to the GSAS-II Python and "
            f"{GSAS2_HOME_VARIABLE} to the folder that contains the GSASII package."
        )
    return Gsas2Install(python=python, home=home)


def write_instprm(
    path: str | Path,
    caglioti: Caglioti,
    zero: float = 0.0,
    x: float = 0.0,
    y: float = 0.0,
    shl: float = 0.002,
    lam1: float = 1.54056,
    lam2: float = 1.54439,
    ratio: float = 0.5,
    polariz: float = 0.7,
) -> Path:
    """Write a GSAS-II instrument parameter file for a Cu K alpha lab pattern.

    The file is the ``key:value`` form GSAS-II reads, one bank of constant
    wavelength X-ray data (``Type:PXC``) with the K alpha 1 and K alpha 2
    wavelengths, which GSAS-II takes to mean Bragg-Brentano geometry.

    The profile widths come from ``caglioti``. xrdkit's U, V and W give the
    squared FWHM in degrees squared, whereas GSAS-II's give the variance of
    the Gaussian component in centidegrees squared, sigma^2 =
    U tan^2 + V tan + W, its FWHM being sqrt(8 ln 2) sigma. Each is therefore
    multiplied by 1e4 / (8 ln 2), about 1803: 1e4 for degrees squared to
    centidegrees squared and 1 / (8 ln 2) for FWHM^2 to sigma^2. This treats
    the whole of the fitted width as Gaussian, which is exact while ``x`` and
    ``y`` are zero and a starting point otherwise.

    Parameters
    ----------
    path
        File to write, conventionally with the ``.instprm`` extension.
    caglioti
        Instrumental resolution, as from :func:`~xrdkit.broadening.fit_caglioti`.
    zero
        Zero shift, in degrees.
    x, y
        Lorentzian widths X / cos(theta) + Y tan(theta), in centidegrees.
    shl
        Axial divergence, (S + H) / L.
    lam1, lam2
        K alpha 1 and K alpha 2 wavelengths, in angstroms.
    ratio
        K alpha 2 over K alpha 1 intensity.
    polariz
        Polarisation fraction.

    Returns
    -------
    Path
        The file written.
    """
    scale = CENTIDEGREES_SQUARED / EIGHT_LN2
    # In the order GSAS-II itself lists the parameters of a dual wavelength
    # pattern.
    values = {
        "Type": "PXC",
        "Bank": 1,
        "Lam1": lam1,
        "Lam2": lam2,
        "Zero": zero,
        "I(L2)/I(L1)": ratio,
        "Polariz.": polariz,
        "U": caglioti.u * scale,
        "V": caglioti.v * scale,
        "W": caglioti.w * scale,
        "X": x,
        "Y": y,
        "Z": 0.0,
        "SH/L": shl,
        "Azimuth": 0.0,
    }
    lines = [INSTPRM_HEADER] + [f"{key}:{value}" for key, value in values.items()]
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _environment(install: Gsas2Install) -> dict[str, str]:
    """The environment the GSAS-II Python runs in.

    GSAS-II opens its data files as text in the locale encoding, which on
    Windows turns the byte order mark that PANalytical writes at the head of
    an ``.xrdml`` file into characters the XML parser rejects. Python's UTF-8
    mode, the default from Python 3.15, reads it correctly, and is turned on
    unless ``PYTHONUTF8`` is already set.

    The GSAS-II Python is a conda environment. On Windows its numpy finds the
    BLAS and LAPACK DLLs only through the folders ``conda activate`` puts on
    PATH, and without them the first matrix inversion kills the process, so
    those of them that exist are put at the front of PATH here.
    """
    environment = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    environment.setdefault("PYTHONUTF8", "1")
    if sys.platform == "win32":
        prefix = install.python.parent
        folders = [
            prefix,
            prefix / "Library" / "mingw-w64" / "bin",
            prefix / "Library" / "usr" / "bin",
            prefix / "Library" / "bin",
            prefix / "Scripts",
        ]
        path = [str(folder) for folder in folders if folder.is_dir()]
        environment["PATH"] = os.pathsep.join([*path, environment.get("PATH", "")])
    return environment


def _write_xy(data_file: Path, workdir: Path) -> Path:
    """A two column copy of an ``.xrdml`` scan in ``workdir``, read by xrdkit."""
    from xrdkit.io import read_xrdml

    scan = read_xrdml(data_file)
    path = workdir / f"{data_file.stem}.xy"
    rows = (f"{t:.6f} {i:.10g}" for t, i in zip(scan.two_theta, scan.intensity))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def run_job(
    job: dict, workdir: str | Path, install: Gsas2Install | None = None
) -> dict:
    """Run a job under the GSAS-II Python and return its result.

    The job, with ``gsas2_home`` added, is written to
    ``workdir/<action>_job.json`` and ``gsas2_driver.py`` run on it as
    ``python -B driver job.json`` in ``workdir``, so that nothing, bytecode
    included, lands in the GSAS-II installation, with the conda environment's
    DLL folders on PATH on Windows. Its output is kept in
    ``workdir/<action>.log``. The driver writes its result to ``job["result"]``,
    by default ``workdir/<action>_result.json``.

    When ``job["data_file"]`` is an ``.xrdml`` scan, a two column ``.xy`` copy
    is written to ``workdir`` first and passed as ``job["data_fallback"]``,
    for the driver to read should GSAS-II's own xrdml importer fail.

    Parameters
    ----------
    job
        The job; ``job["action"]`` selects what the driver does. Not modified.
    workdir
        Folder for the job, its result and its log; created if need be.
    install
        The installation to use, :func:`find_gsas2` by default.

    Raises
    ------
    Gsas2Error
        If the driver exits with an error, carrying its stderr.
    """
    install = install or find_gsas2()
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    action = job["action"]

    job = dict(job)
    job["gsas2_home"] = str(install.home)
    job.setdefault("result", str(workdir / f"{action}_result.json"))
    data_file = job.get("data_file")
    if (
        data_file
        and Path(data_file).suffix.lower() == ".xrdml"
        and "data_fallback" not in job
    ):
        job["data_fallback"] = str(_write_xy(Path(data_file), workdir))

    job_path = workdir / f"{action}_job.json"
    job_path.write_text(json.dumps(job, indent=2), encoding="utf-8")
    result_path = Path(job["result"])
    result_path.unlink(missing_ok=True)

    completed = subprocess.run(
        [str(install.python), "-B", str(DRIVER), str(job_path)],
        check=False,
        cwd=workdir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=_environment(install),
    )
    (workdir / f"{action}.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise Gsas2Error(
            f"GSAS-II job {action!r} failed with exit code "
            f"{completed.returncode}:\n{completed.stderr.strip()}",
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
    if not result_path.is_file():
        raise Gsas2Error(
            f"GSAS-II job {action!r} wrote no result at {result_path}",
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
    return json.loads(result_path.read_text(encoding="utf-8"))
