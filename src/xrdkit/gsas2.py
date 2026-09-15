"""GSAS-II jobs, run under GSAS-II's own Python.

GSAS-II brings its own Python and its own compiled extensions, and cannot be
imported here. Instead a job, a dict naming an action and its inputs, is
written as JSON, and ``gsas2_driver.py`` runs it under the GSAS-II Python
through GSASIIscriptable and writes its result as JSON for :func:`run_job` to
read back. The driver imports nothing from xrdkit.

:func:`find_gsas2` locates the installation, and :func:`write_instprm` writes
the instrument parameter file a powder histogram is read with.
:func:`build_refine_job` assembles a refinement job from a list of stages,
such as the usual sequence :func:`standard_stages` returns, and
:func:`gsas2_fwhm` gives the line width refined instrument parameters imply.
:func:`structure_edits` turns a refined structure into the atom edits that
set it up again in a new project, so that one refinement can start where
another ended.

This module imports only the standard library; the exceptions, the xrdkit
reader used to write a two column copy of an ``.xrdml`` scan, numpy for
widths at many angles at once and the bond lengths a job may ask for, are
imported when needed.
"""

from __future__ import annotations

import copy
import json
import math
import os
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from xrdkit import gsas2_driver

if TYPE_CHECKING:
    import numpy as np

    from xrdkit.broadening import Caglioti

__all__ = [
    "GSAS2_HOME_VARIABLE",
    "GSAS2_PYTHON_VARIABLE",
    "Gsas2Error",
    "Gsas2Install",
    "build_refine_job",
    "find_gsas2",
    "gsas2_fwhm",
    "run_job",
    "standard_stages",
    "structure_edits",
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

# The Gaussian FWHM over sigma, the smallest Gaussian variance in
# centidegrees squared, and the Thompson, Cox and Hastings quintic of
# coefficients of G^5, G^4 L, ... L^5, as GSASIIpwd.getFWHM and getgamFW
# write them.
GSAS2_SQRT_8LN2 = 2.35482
GSAS2_MIN_VARIANCE = 0.001
GSAS2_TCH = (1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0)

# Least squares cycles a refinement stage may take; GSAS-II's own default is 3.
DEFAULT_CYCLES = 10

# What a job's bonds may set: see xrdkit.structure.bond_lengths.
BOND_KEYS = frozenset({"anions", "dmax", "dmin"})

# Atoms closer than this in every fractional coordinate share a site.
SAME_SITE = 1.0e-4


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


def gsas2_fwhm(
    two_theta: float | np.ndarray,
    u: float,
    v: float,
    w: float,
    x: float,
    y: float,
    shl: float = 0.0,
    z: float = 0.0,
) -> float | np.ndarray:
    """The FWHM GSAS-II gives a line at ``two_theta``, in degrees.

    The profile parameters are in GSAS-II's units, as refined: U, V and W
    give the Gaussian variance sigma^2 = U tan^2 + V tan + W in centidegrees
    squared, held at no less than 0.001, and X, Y and Z the Lorentzian FWHM
    gamma = X / cos + Y tan + Z in centidegrees, theta being half of
    ``two_theta``. The two are combined into the pseudo-Voigt FWHM by the
    Thompson, Cox and Hastings quintic, with the Gaussian FWHM taken as
    2.35482 sigma, exactly as GSASIIpwd.getFWHM and getgamFW do.

    ``shl`` is accepted so that the refined parameters can be passed whole,
    but GSAS-II's FWHM does not depend on it: the axial divergence it sets
    makes the line asymmetric without entering the width.
    """
    del shl
    import numpy as np

    theta = np.radians(np.asarray(two_theta, dtype=float) / 2.0)
    tan_theta = np.tan(theta)
    variance = np.maximum(GSAS2_MIN_VARIANCE, u * tan_theta**2 + v * tan_theta + w)
    gaussian = GSAS2_SQRT_8LN2 * np.sqrt(variance)
    lorentzian = z + x / np.cos(theta) + y * tan_theta
    quintic = sum(
        coefficient * gaussian ** (5 - order) * lorentzian**order
        for order, coefficient in enumerate(GSAS2_TCH)
    )
    fwhm = quintic**0.2 / 100.0
    return float(fwhm) if fwhm.ndim == 0 else fwhm


def standard_stages(
    background_type: str = gsas2_driver.DEFAULT_BACKGROUND["type"],
    background_terms: int = gsas2_driver.DEFAULT_BACKGROUND["terms"],
) -> list[dict]:
    """The usual sequence of refinement stages, for a script to edit.

    Background and scale; zero; cell; U, V and W; X and Y; SH/L. Each stage
    adds to those before it. The list is new on every call, so it can be
    changed freely, a stage dropped or a flag added.
    """
    return [
        {
            "name": "background and scale",
            "background": {"type": background_type, "terms": background_terms},
            "scale": True,
        },
        {"name": "zero", "zero": True},
        {"name": "cell", "cell": True},
        {"name": "U V W", "instrument": ["U", "V", "W"]},
        {"name": "X Y", "instrument": ["X", "Y"]},
        {"name": "SH/L", "instrument": ["SH/L"]},
    ]


def build_refine_job(
    gpx: str | Path,
    stages: Sequence[Mapping],
    data_file: str | Path | None = None,
    instprm: str | Path | None = None,
    phases: Sequence[Mapping] | None = None,
    limits: tuple[float, float] | None = None,
    cycles: int = DEFAULT_CYCLES,
    broadening: Mapping[str, Mapping[str, float]] | None = None,
    export_prefix: str | Path | None = None,
    le_bail_cycles: int | None = None,
    max_passes: int | None = None,
    pass_tolerance: float | None = None,
    background_start: Mapping | None = None,
    overall_uiso_start: Mapping[str, float] | None = None,
    scale_start: float | None = None,
    sanity: Mapping | None = None,
    bonds: bool | Mapping | None = None,
    on_flagged: str | None = None,
    on_unsettled: str | None = None,
) -> dict:
    """A job for :func:`run_job` that refines a project stage by stage.

    With ``data_file``, ``instprm`` and ``phases`` the project is created at
    ``gpx`` first, as the create action does; without them ``gpx`` must
    already hold one with a single histogram. The stages are checked here,
    with the driver's own rules, so that a mistake is reported before
    GSAS-II starts.

    Parameters
    ----------
    gpx
        The project file.
    stages
        Refinement stages, as described in ``gsas2_driver``;
        :func:`standard_stages` gives the usual sequence. Copied.
    data_file, instprm
        The pattern and its instrument parameter file.
    phases
        Each a mapping with ``cif`` and ``name``, and optionally ``cell``,
        six lattice parameters that replace those of the CIF, and ``atoms``,
        edits to its atoms made in the project, as
        ``gsas2_driver.check_atom_edits`` describes: ``{"label": ...,
        "occupancy": ...}``, or ``{"label": new, "type": element, "copy":
        existing, "occupancy": ...}`` for another element on an existing
        site, either with ``xyz`` and ``uiso`` if need be;
        :func:`structure_edits` gives those that set up a refined structure.
    limits
        Two theta range to refine over, in degrees; the whole pattern by
        default.
    cycles
        The most least squares cycles a stage may take.
    broadening
        Isotropic sample broadening to start from, by phase name or ``"*"``
        for every phase: ``{"size": microns, "mustrain": microstrain,
        "lgmix": Lorentzian fraction}``, any of them, held fixed until a
        stage refines the size or microstrain. GSAS-II's defaults,
        1 micron and 1000, both Lorentzian, are far from nothing, and the
        size cannot exceed 10 microns; refining the instrument on a standard
        needs ``{"*": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}}``, or the
        instrument terms take up the difference.
    export_prefix
        Path the exported files are named from: ``<prefix>_histogram.csv``,
        ``<prefix>_reflections_<phase>.csv`` and ``<prefix>.instprm``. By
        default ``gpx`` without a ``.gpx`` extension, so that a stem such as
        ``x0.10`` keeps its dot.
    le_bail_cycles
        Le Bail-only cycles to run whenever a stage switches Le Bail
        extraction on, before its least squares; the driver's own default,
        ``gsas2_driver.DEFAULT_LE_BAIL_CYCLES``, if not given.
    max_passes, pass_tolerance
        Refine each stage up to ``max_passes`` times, until no parameter
        moves by more than ``pass_tolerance`` esds (the driver's
        ``DEFAULT_PASS_TOLERANCE`` if not given) from one pass to the next.
        A Le Bail refinement needs this, since GSAS-II stops each one short;
        by default a stage is refined once.
    background_start
        ``{"type": ..., "coefficients": [...]}``, the background to start
        from, such as an earlier refinement's; GSAS-II's own otherwise.
    overall_uiso_start
        The Uiso, by phase name or ``"*"``, that every atom of a phase is
        given when a stage first refines its ``overall_uiso``, or that an
        anisotropic atom is given when a Uiso group takes it in.
    scale_start
        The histogram scale to start from, such as an earlier refinement's;
        GSAS-II's 1 otherwise.
    sanity
        ``{"max_shift": ..., "reference": {phase: [{"label": ..., "xyz":
        [...]}, ...]}}``, either, for the sanity check after every stage: a
        site is flagged when a coordinate has moved more than ``max_shift``
        (fractional) from its reference position. By default
        ``gsas2_driver.DEFAULT_MAX_SHIFT``, from the atoms as the job finds
        them before its first stage.
    bonds
        True, or ``{"anions": [...], "dmax": ..., "dmin": ...}``, any of them,
        for :func:`run_job` to add to the result, as ``bonds``, the cation to
        anion distances of each phase's final structure, as
        :func:`xrdkit.structure.bond_lengths` gives them: to oxygen, from 0.5
        to 3 angstroms, by default.
    on_flagged, on_unsettled
        "accept" or "reject": what becomes of a stage whose sanity check
        raises a flag the last stage kept did not have, and of one not
        settled after ``max_passes``. A rejected stage is recorded and the
        project taken back to the last stage kept, with what the stage
        alone refined held from then on (see the ``gsas2_driver`` notes). By
        default a new flag rejects a stage and an unsettled one is kept.

    Every path, ``gpx``, ``data_file``, ``instprm``, each phase's ``cif`` and
    ``export_prefix``, is made absolute against the current working directory,
    since the driver runs in a working folder of its own.

    Raises
    ------
    Gsas2Error
        If ``data_file``, ``instprm`` or a phase's ``cif`` does not exist.
    ValueError
        If only some of ``data_file``, ``instprm`` and ``phases`` are given, a
        phase lacks ``cif`` or ``name``, the limits are not increasing,
        ``cycles``, ``le_bail_cycles``, ``max_passes``, ``pass_tolerance`` or
        ``scale_start`` is not positive, or a stage, the broadening, an atom
        edit, the starting background, a starting Uiso, the sanity settings or
        the bonds are malformed.
    TypeError
        If a stage or one of its lists is of the wrong type.
    """
    creation = (data_file, instprm, phases)
    if any(item is not None for item in creation) and any(
        item is None for item in creation
    ):
        raise ValueError("give all of data_file, instprm and phases, or none")
    if not isinstance(cycles, int) or cycles < 1:
        raise ValueError(f"cycles must be a positive int, got {cycles!r}")
    if le_bail_cycles is not None and (
        not isinstance(le_bail_cycles, int) or le_bail_cycles < 1
    ):
        raise ValueError(
            f"le_bail_cycles must be a positive int, got {le_bail_cycles!r}"
        )
    if max_passes is not None and (not isinstance(max_passes, int) or max_passes < 1):
        raise ValueError(f"max_passes must be a positive int, got {max_passes!r}")
    if pass_tolerance is not None and not pass_tolerance > 0.0:
        raise ValueError(f"pass_tolerance must be positive, got {pass_tolerance!r}")
    stages = copy.deepcopy([dict(stage) for stage in stages])
    gsas2_driver.accumulate_stages(stages)

    gpx = Path(gpx).resolve()
    if export_prefix is None:
        # Only a .gpx extension is taken off: with_suffix("") would also take
        # the ".10" off a stem such as x0.10.
        name = gpx.name[: -len(".gpx")] if gpx.suffix.lower() == ".gpx" else gpx.name
        export_prefix = gpx.parent / name
    job = {
        "action": "refine",
        "gpx": str(gpx),
        "stages": stages,
        "cycles": cycles,
        "export_prefix": str(Path(export_prefix).resolve()),
    }
    if limits is not None:
        lower, upper = (float(value) for value in limits)
        if not lower < upper:
            raise ValueError(f"limits must increase, got {limits}")
        job["limits"] = [lower, upper]
    if le_bail_cycles is not None:
        job["le_bail_cycles"] = le_bail_cycles
    if max_passes is not None:
        job["max_passes"] = max_passes
    if pass_tolerance is not None:
        job["pass_tolerance"] = float(pass_tolerance)
    if broadening is not None:
        job["broadening"] = gsas2_driver.check_broadening(
            {str(phase): dict(values) for phase, values in broadening.items()}
        )
    if background_start is not None:
        job["background_start"] = gsas2_driver.check_background_start(
            dict(background_start)
        )
    if overall_uiso_start is not None:
        job["overall_uiso_start"] = gsas2_driver.check_uiso_start(
            dict(overall_uiso_start)
        )
    if scale_start is not None:
        job["scale_start"] = gsas2_driver.check_scale_start(scale_start)
    if sanity is not None:
        checked = gsas2_driver.check_sanity_settings(dict(sanity))
        job["sanity"] = {"max_shift": checked["max_shift"]}
        if checked["reference"] is not None:
            job["sanity"]["reference"] = checked["reference"]
    if bonds is not None and bonds is not False:
        job["bonds"] = _check_bonds(bonds)
    policies = {"on_flagged": on_flagged, "on_unsettled": on_unsettled}
    gsas2_driver.check_policies(policies)
    job.update({key: value for key, value in policies.items() if value is not None})
    if data_file is not None:
        entries = []
        for phase in phases:
            if "cif" not in phase or "name" not in phase:
                raise ValueError(f"each phase needs a cif and a name, got {phase}")
            entry = {"cif": str(phase["cif"]), "name": str(phase["name"])}
            if phase.get("cell") is not None:
                entry["cell"] = [float(value) for value in phase["cell"]]
            if phase.get("atoms"):
                entry["atoms"] = gsas2_driver.check_atom_edits(
                    [dict(edit) for edit in phase["atoms"]], entry["name"]
                )
            entry["cif"] = str(_existing_input(phase["cif"], f"CIF of {entry['name']}"))
            entries.append(entry)
        job.update(
            data_file=str(_existing_input(data_file, "data file")),
            instprm=str(_existing_input(instprm, "instrument parameter file")),
            phases=entries,
        )
    return job


def _existing_input(path: str | Path, what: str) -> Path:
    """Return ``path`` made absolute, or raise Gsas2Error if it is not a file."""
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise Gsas2Error(f"{what} not found: {resolved}")
    return resolved


def _check_bonds(bonds: bool | Mapping) -> dict:
    """Check a job's ``bonds``: True for the defaults, or some of
    ``anions``, ``dmax`` and ``dmin``."""
    if bonds is True:
        return {}
    if not isinstance(bonds, Mapping) or set(bonds) - BOND_KEYS:
        raise ValueError(
            f"bonds must be True or give anions, dmax and/or dmin, got {bonds!r}"
        )
    checked: dict = {}
    if "anions" in bonds:
        checked["anions"] = [str(anion) for anion in bonds["anions"]]
        if not checked["anions"]:
            raise ValueError("bonds needs at least one anion")
    for key in ("dmax", "dmin"):
        if key in bonds:
            checked[key] = float(bonds[key])
    if not checked.get("dmax", 3.0) > checked.get("dmin", 0.5):
        raise ValueError(f"bonds dmax must exceed dmin, got {bonds!r}")
    return checked


def _same_site(first: Sequence[float], second: Sequence[float]) -> bool:
    return all(abs((a - b) - round(a - b)) < SAME_SITE for a, b in zip(first, second))


def structure_edits(refined: Sequence[Mapping], base: Sequence[Mapping]) -> list[dict]:
    """Atom edits that set up a refined structure again in a new project
    read from the same CIF, for a phase entry's ``atoms``.

    ``refined`` is a phase's atoms as a refine result reports them, in its
    ``final`` phases or a stage's ``atoms``, and ``base`` the atoms the CIF
    gives, as the create action reports them. Each refined atom is given its
    occupancy, its position and, if it is isotropic, its Uiso. An atom the
    CIF lacks, such as one added on a site by an earlier edit, is added
    beside the CIF atom that shares its site in the refinement, with its
    element, after every edit of the CIF's atoms. An anisotropic atom keeps
    the Uij of the CIF, or of the atom it is added beside.

    Raises
    ------
    ValueError
        If an atom the CIF lacks shares no site with one it has.
    """
    labels = {str(atom["label"]) for atom in base}
    edits, added = [], []
    for atom in refined:
        edit = {
            "label": str(atom["label"]),
            "occupancy": float(atom["occupancy"]),
            "xyz": [float(value) for value in atom["xyz"]],
        }
        if atom.get("adp", "I") == "I" and atom.get("uiso") is not None:
            edit["uiso"] = float(atom["uiso"])
        if edit["label"] in labels:
            edits.append(edit)
            continue
        host = next(
            (
                other
                for other in refined
                if str(other["label"]) in labels
                and _same_site(other["xyz"], atom["xyz"])
            ),
            None,
        )
        if host is None:
            raise ValueError(
                f"atom {edit['label']!r} shares no site with an atom of the CIF"
            )
        added.append({**edit, "type": str(atom["type"]), "copy": str(host["label"])})
    return edits + added


def _add_bonds(result: dict, options: Mapping) -> None:
    """Add the cation to anion distances of each phase's final structure to
    ``result``, as ``bonds`` by phase name, or the error as ``bonds_error``."""
    from xrdkit.structure import bond_lengths

    try:
        result["bonds"] = {
            phase["name"]: bond_lengths(phase, **options)
            for phase in result["final"]["phases"]
            if phase.get("operators")
        }
    except (KeyError, TypeError, ValueError) as error:
        result["bonds_error"] = f"{type(error).__name__}: {error}"


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
    for the driver to read should GSAS-II's own xrdml importer fail. When the
    job has ``bonds`` (see :func:`build_refine_job`) and the result a final
    structure, its bond lengths are added to the result, and the result file
    written again with them.

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
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if "bonds" in job and result.get("final"):
        _add_bonds(result, job["bonds"])
        result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


# Stage statuses and write ups of a sequence of refinements

# The statuses of stages refined and kept; the others are rejected and
# failed.
ACCEPTED_STATUSES = frozenset({"clean", "flagged", "unsettled"})
LOG_TAIL_LINES = 40


def stage_status(stage: Mapping) -> str:
    """A refine result stage's status, as the driver records it: clean,
    flagged, unsettled, rejected or failed (see the ``gsas2_driver``
    notes); worked out from ``error`` for a result written before the driver
    recorded one."""
    if stage.get("status"):
        return str(stage["status"])
    return "failed" if "error" in stage else "clean"


def accepted_stages(result: Mapping) -> list[dict]:
    """The stages of a refine result that were refined and kept."""
    return [
        stage
        for stage in result.get("stages", [])
        if stage_status(stage) in ACCEPTED_STATUSES
    ]


def stage_statuses(result: Mapping) -> list[dict]:
    """Each stage of a refine result as ``{"name", "status", "reason",
    "passes", "rwp", "gof", "undetermined"}``.

    ``reason`` is why a stage was rejected or failed, the new flags of one
    flagged, and for one unsettled the passes it took and its largest
    remaining move in esds; empty for a clean stage. ``undetermined`` lists
    the messages of the parameters the driver found undetermined in it.
    """
    rows = []
    for stage in result.get("stages", []):
        status = stage_status(stage)
        passes = len(stage.get("passes") or []) or None
        if status == "rejected":
            reason = "; ".join(stage.get("rejected_because", []))
        elif status == "failed":
            reason = str(stage.get("error", ""))
        elif status == "unsettled":
            move = stage.get("largest_remaining_move") or {}
            reason = f"not settled in {passes} passes"
            if move.get("shift_over_esd") is not None:
                reason += (
                    f"; largest remaining move {move['shift_over_esd']:.2f} esd "
                    f"({move['parameter']})"
                )
        elif status == "flagged":
            reason = "; ".join(flag["message"] for flag in stage.get("sanity", []))
        else:
            reason = ""
        rows.append(
            {
                "name": stage.get("name"),
                "status": status,
                "reason": reason,
                "passes": passes,
                "rwp": stage.get("rwp"),
                "gof": stage.get("gof"),
                "undetermined": [
                    entry["message"] for entry in stage.get("undetermined", [])
                ],
            }
        )
    return rows


def log_tail(path: str | Path, lines: int = LOG_TAIL_LINES) -> str:
    """The last ``lines`` lines of a log file; empty if there is none."""
    path = Path(path)
    if not path.is_file():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(text[-lines:])


def _figure(value: float | None, digits: int) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


def _cell(text: object) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def stage_status_table(result: Mapping) -> list[str]:
    """The stages of a refine result as markdown table lines: status,
    passes, Rwp, GOF and why."""
    lines = [
        "| Stage | Status | Passes | Rwp (%) | GOF | Why |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in stage_statuses(result):
        lines.append(
            f"| {_cell(row['name'])} | {row['status']} | {row['passes'] or '—'} | "
            f"{_figure(row['rwp'], 3)} | {_figure(row['gof'], 3)} | "
            f"{_cell(row['reason']) or '—'} |"
        )
    return lines


def failure_markdown(
    title: str,
    result: Mapping | None,
    error: str,
    log: str,
    intro: Sequence[str] = (),
) -> str:
    """A write up of a refinement that failed: ``intro``, the error, the
    stages it got through with their statuses, and ``log``, the tail of the
    GSAS-II log."""
    lines = [f"# {title}", ""]
    if intro:
        lines += [*intro, ""]
    lines += ["## Error", "", "```", error.strip(), "```", "", "## Stages", ""]
    if result and result.get("stages"):
        lines += stage_status_table(result)
    else:
        lines.append("No stage was refined.")
    lines += [
        "",
        f"## GSAS-II log, last {LOG_TAIL_LINES} lines",
        "",
        "```",
        log.strip() or "(no log)",
        "```",
    ]
    return "\n".join(lines) + "\n"


def summary_markdown(
    title: str,
    entries: Sequence[Mapping],
    columns: Sequence[str] = (),
    intro: Sequence[str] = (),
) -> str:
    """A summary of a sequence of refinements as markdown: one table with a
    row per refinement, its final Rwp and GOF, ``columns`` of key values and
    the status of each stage, then every stage not clean with why, and the
    parameters left undetermined.

    Each entry is ``{"name", "run", "result", "values"}``: its name, when it
    was run, its refine result (None if it has none) and its key values by
    column, as text. The final Rwp and GOF are those of the last stage kept;
    a refinement that kept no stage is noted as leaving the model it started
    from.
    """
    header = ["Refinement", "Run", "Outcome", "Rwp (%)", "GOF", *columns, "Stages"]
    lines = [f"# {title}", ""]
    if intro:
        lines += [*intro, ""]
    lines += ["| " + " | ".join(header) + " |", "|" + " --- |" * len(header)]
    details = []
    for entry in entries:
        result = entry.get("result")
        if result is None:
            cells = [entry["name"], entry.get("run") or "—", "not run", "—", "—"]
            cells += ["—"] * len(columns) + ["—"]
            lines.append("| " + " | ".join(_cell(cell) for cell in cells) + " |")
            continue
        kept = accepted_stages(result)
        final = kept[-1] if kept else {}
        rows = stage_statuses(result)
        values = entry.get("values") or {}
        cells = [
            entry["name"],
            entry.get("run") or "—",
            "completed" if result.get("completed") else "failed",
            _figure(final.get("rwp"), 3),
            _figure(final.get("gof"), 3),
            *(values.get(column, "—") for column in columns),
            "; ".join(f"{row['name']}: {row['status']}" for row in rows) or "none",
        ]
        lines.append("| " + " | ".join(_cell(cell) for cell in cells) + " |")
        if result.get("error"):
            details.append(f"- {entry['name']}: {result['error']}")
        if not kept and result.get("stages"):
            details.append(
                f"- {entry['name']}: no stage was kept, so the model left behind "
                "is the one the refinement started from"
            )
        for row in rows:
            if row["status"] != "clean":
                details.append(
                    f"- {entry['name']}, {row['name']}: {row['status']}"
                    + (f", {row['reason']}" if row["reason"] else "")
                )
            if row["undetermined"] and (
                row["status"] == "rejected" or row["name"] == final.get("name")
            ):
                where = "rejected" if row["status"] == "rejected" else "final model"
                details.append(
                    f"- {entry['name']}, {row['name']} ({where}): undetermined "
                    + "; ".join(row["undetermined"])
                )
    lines += ["", "## Stages not clean, and undetermined parameters", ""]
    lines += details or ["Every stage is clean, and no parameter is left undetermined."]
    return "\n".join(lines) + "\n"
