"""Run one xrdkit job inside GSAS-II.

This script runs under GSAS-II's own Python, started by
:func:`xrdkit.gsas2.run_job` as ``python -B gsas2_driver.py job.json``. It
must not import xrdkit, which is not installed there. The job JSON names an
``action`` and its inputs, and ``gsas2_home``, the folder containing the
GSASII package; the result is written as JSON to ``job["result"]``.

Each action is a function of GSASIIscriptable and the job returning the result
dict, registered in ``ACTIONS``. GSASII is imported only once the job is read,
so that xrdkit can import this module to check stages without it.

Actions
-------
``create``
    A new project with one histogram and its phases linked to it.
``refine``
    A sequence of refinement stages on a project, created first if the job
    carries the inputs for it, with the pattern, reflections and refined
    instrument parameters exported after the last stage.

Stages
------
A stage is a dict with a ``name`` and the refinement flags to switch on:

``background``
    ``True`` for :data:`DEFAULT_BACKGROUND`, or ``{"type": ..., "terms": n}``.
``scale``, ``zero``, ``displacement``
    Histogram scale, zero shift, and specimen displacement (``Shift`` in
    Bragg-Brentano geometry, ``DisplaceX`` and ``DisplaceY`` in
    Debye-Scherrer).
``instrument``
    A list drawn from :data:`INSTRUMENT_KEYS`.
``cell``, ``le_bail``, ``phase_fractions``
    ``True`` for every phase or a list of phase names.
``atoms``
    A list of ``{"phase": name, "flags": "FXU"}``, any of F (occupancy),
    X (position) and U (displacement).

Flags carry over from stage to stage, a list growing by the new entries, a
``True`` or ``False`` replacing what was there, unless the stage has
``"reset": true``, when it starts again from nothing refined.

Sample broadening
-----------------
GSAS-II gives every phase an isotropic crystallite size of 1 micron and a
microstrain of 1000, both Lorentzian, and adds the width they imply to the
instrumental width of every reflection. Refining instrument parameters on a
standard with those left in place makes the instrument terms absorb the
difference, and they then no longer describe the instrument. A refine job's
``broadening``, ``{phase name or "*": {"size": microns, "mustrain": value,
"lgmix": fraction}}``, any of the three, sets isotropic values, held fixed,
before the first stage.

GSAS-II holds a size to 0.001 to 10 microns, and at 10 microns a Lorentzian
size term still adds 0.088 / cos(theta) centidegrees, the same form as X,
which X then gives up. ``{"size": 10, "mustrain": 0, "lgmix": 0}`` makes the
size term Gaussian instead, adding some 0.0014 / cos^2(theta) centidegrees
squared to a variance of several, so that the sample adds next to nothing.
"""

import copy
import csv
import io
import json
import math
import re
import sys
from pathlib import Path

# Only importers whose format name contains the hint are tried.
XRDML_HINT = "Panalytical"
XY_HINT = "comma/tab/semicolon"
CIF_HINT = "CIF"

DEFAULT_BACKGROUND = {"type": "chebyschev-1", "terms": 6}

# Profile parameters a stage may refine, besides the zero shift.
INSTRUMENT_KEYS = ("U", "V", "W", "X", "Y", "Z", "SH/L")

# Atom refinement flags: occupancy, position and displacement.
ATOM_FLAGS = "FXU"

STAGE_KEYS = frozenset(
    {
        "name",
        "reset",
        "background",
        "scale",
        "zero",
        "displacement",
        "cell",
        "instrument",
        "atoms",
        "le_bail",
        "phase_fractions",
    }
)
BOOLEAN_FLAGS = ("scale", "zero", "displacement")
PHASE_FLAGS = ("cell", "le_bail", "phase_fractions")

# In a phase flag, every phase.
ALL_PHASES = "*"

# Sample parameters the stages switch on and off, and the displacement
# parameters of each geometry.
SAMPLE_KEYS = ("Scale", "Shift", "DisplaceX", "DisplaceY")
DISPLACEMENT = {
    "Bragg-Brentano": ["Shift"],
    "Debye-Scherrer": ["DisplaceX", "DisplaceY"],
}

# Sample broadening a refine job may hold, and the isotropic crystallite
# size range in microns that GSAS-II clamps a size to when it saves one.
BROADENING_KEYS = frozenset({"size", "mustrain", "lgmix"})
MIN_SIZE = 0.001
MAX_SIZE = 10.0

# Instrument parameters reported after a refinement.
REPORTED_INSTRUMENT = ("Zero", *INSTRUMENT_KEYS)

# Lines of refinement output kept as the report of a failed stage.
FAILURE_LINES = 20
FAILURE_WORDS = re.compile(r"error|fail|singular|problem|abort", re.IGNORECASE)


def _plain(value):
    """``value`` with numpy scalars and arrays turned into JSON types."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


# Stages


def empty_flags():
    """Flags with nothing refined."""
    return {
        "background": None,
        "scale": False,
        "zero": False,
        "displacement": False,
        "instrument": [],
        "cell": [],
        "le_bail": [],
        "phase_fractions": [],
        "atoms": {},
    }


def _background(value, name):
    if value is False or value is None:
        return None
    if value is True:
        return dict(DEFAULT_BACKGROUND)
    if isinstance(value, dict):
        unknown = set(value) - {"type", "terms"}
        if unknown:
            raise ValueError(
                f"stage {name!r}: unknown background keys {sorted(unknown)}"
            )
        terms = value.get("terms", DEFAULT_BACKGROUND["terms"])
        if not isinstance(terms, int) or isinstance(terms, bool) or terms < 1:
            raise ValueError(f"stage {name!r}: background terms must be a positive int")
        return {
            "type": str(value.get("type", DEFAULT_BACKGROUND["type"])),
            "terms": terms,
        }
    raise ValueError(f"stage {name!r}: background must be a bool or a dict")


def _phases(value, current, key, name):
    """A phase flag: every phase, some phases by name, or none."""
    if value is True:
        return [ALL_PHASES]
    if value is False or value is None:
        return []
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise ValueError(
            f"stage {name!r}: {key} must be a bool or a list of phase names"
        )
    if ALL_PHASES in current:
        return current
    return current + [phase for phase in value if phase not in current]


def _instrument(value, current, name):
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list):
        raise TypeError(f"stage {name!r}: instrument must be a list")
    unknown = [key for key in value if key not in INSTRUMENT_KEYS]
    if unknown:
        raise ValueError(
            f"stage {name!r}: unknown instrument parameters {unknown}; "
            f"choose from {', '.join(INSTRUMENT_KEYS)}"
        )
    chosen = set(current) | set(value)
    return [key for key in INSTRUMENT_KEYS if key in chosen]


def _atoms(value, current, name):
    if not isinstance(value, list):
        raise TypeError(f"stage {name!r}: atoms must be a list")
    atoms = dict(current)
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != {"phase", "flags"}:
            raise ValueError(
                f'stage {name!r}: each atoms entry must be {{"phase": ..., "flags": ...}}'
            )
        flags = str(entry["flags"]).upper()
        if any(flag not in ATOM_FLAGS for flag in flags):
            raise ValueError(
                f"stage {name!r}: atom flags {entry['flags']!r} must be drawn from "
                f"{ATOM_FLAGS}"
            )
        chosen = set(atoms.get(entry["phase"], "")) | set(flags)
        atoms[entry["phase"]] = "".join(flag for flag in ATOM_FLAGS if flag in chosen)
    return atoms


def accumulate_stages(stages):
    """The flags in force at each stage, as a list of ``{"name", "flags"}``.

    Raises
    ------
    ValueError
        If there are no stages, or a stage has an unknown key or a flag of the
        wrong kind.
    TypeError
        If a stage is not a dict, or its instrument or atoms not a list.
    """
    if not isinstance(stages, list) or not stages:
        raise ValueError("stages must be a non-empty list")
    flags = empty_flags()
    accumulated = []
    for index, stage in enumerate(stages, start=1):
        if not isinstance(stage, dict):
            raise TypeError(f"stage {index} is not a dict")
        name = str(stage.get("name") or f"stage {index}")
        unknown = set(stage) - STAGE_KEYS
        if unknown:
            raise ValueError(f"stage {name!r}: unknown keys {sorted(unknown)}")
        flags = empty_flags() if stage.get("reset") else copy.deepcopy(flags)
        if "background" in stage:
            flags["background"] = _background(stage["background"], name)
        for key in BOOLEAN_FLAGS:
            if key in stage:
                if not isinstance(stage[key], bool):
                    raise ValueError(f"stage {name!r}: {key} must be true or false")
                flags[key] = stage[key]
        if "instrument" in stage:
            flags["instrument"] = _instrument(
                stage["instrument"], flags["instrument"], name
            )
        for key in PHASE_FLAGS:
            if key in stage:
                flags[key] = _phases(stage[key], flags[key], key, name)
        if "atoms" in stage:
            flags["atoms"] = _atoms(stage["atoms"], flags["atoms"], name)
        accumulated.append({"name": name, "flags": flags})
    return accumulated


def check_broadening(broadening):
    """Check a ``broadening`` mapping, returning it with float values.

    Raises
    ------
    ValueError
        If a phase is given anything but some of ``size``, ``mustrain`` and
        ``lgmix``, or a size lies outside the range GSAS-II holds it to, a
        microstrain is negative or ``lgmix`` lies outside 0 to 1.
    TypeError
        If it is not a dict.
    """
    if not isinstance(broadening, dict):
        raise TypeError("broadening must map phase names to size and mustrain")
    checked = {}
    for phase, values in broadening.items():
        if not isinstance(values, dict) or not values or set(values) - BROADENING_KEYS:
            raise ValueError(
                f"broadening of {phase!r} must give size, mustrain and/or lgmix, "
                f"got {values!r}"
            )
        checked[phase] = {key: float(value) for key, value in values.items()}
        size = checked[phase].get("size", MIN_SIZE)
        if not MIN_SIZE <= size <= MAX_SIZE:
            raise ValueError(
                f"broadening of {phase!r}: size must lie between {MIN_SIZE} and "
                f"{MAX_SIZE} microns, the range GSAS-II clamps it to, got {size}"
            )
        if checked[phase].get("mustrain", 0.0) < 0.0:
            raise ValueError(f"broadening of {phase!r}: mustrain cannot be negative")
        if not 0.0 <= checked[phase].get("lgmix", 0.0) <= 1.0:
            raise ValueError(f"broadening of {phase!r}: lgmix must lie between 0 and 1")
    return checked


def set_broadening(project, histogram, broadening):
    """Set isotropic crystallite size and microstrain, held fixed.

    ``broadening`` maps a phase name, or ``"*"`` for every phase, to any of
    ``size`` in microns, ``mustrain``, and ``lgmix``, the Lorentzian fraction
    of both. The values are written into the phase's data for ``histogram``
    directly, since ``set_HAP_refinements`` can set a size but not a
    microstrain.
    """
    broadening = check_broadening(broadening)
    phases = project.phases()
    missing = sorted(set(broadening) - {phase.name for phase in phases} - {ALL_PHASES})
    if missing:
        raise ValueError(f"no phases named {missing} in the project")
    for phase in phases:
        values = broadening.get(phase.name, broadening.get(ALL_PHASES))
        if not values:
            continue
        hap = phase.data["Histograms"][histogram.name]
        for key, entry in (("size", hap["Size"]), ("mustrain", hap["Mustrain"])):
            if key in values:
                entry[0] = "isotropic"
                entry[1][0] = values[key]
                entry[2][0] = False
            if "lgmix" in values:
                entry[1][2] = values["lgmix"]
                entry[2][2] = False


def _broadening(phase, histogram):
    hap = phase.data["Histograms"][histogram.name]
    size, mustrain = hap["Size"], hap["Mustrain"]
    return {
        "size": {
            "type": size[0],
            "value": size[1][0],
            "lorentzian_fraction": size[1][2],
        },
        "mustrain": {
            "type": mustrain[0],
            "value": mustrain[1][0],
            "lorentzian_fraction": mustrain[1][2],
        },
    }


def _selected(phases, name):
    return ALL_PHASES in phases or name in phases


def apply_flags(project, histogram, flags):
    """Switch off every flag the stages control, then switch on ``flags``.

    Each ``set_refinements`` and ``clear_refinements`` call carries a single
    key: GSASIIscriptable returns from either as soon as it has handled a
    ``Background`` given as a bool, dropping any keys after it.
    """
    phases = project.phases()
    names = {phase.name for phase in phases}
    named = [
        *flags["cell"],
        *flags["le_bail"],
        *flags["phase_fractions"],
        *flags["atoms"],
    ]
    missing = sorted({phase for phase in named if phase != ALL_PHASES} - names)
    if missing:
        raise ValueError(
            f"no phases named {missing} in the project; have {sorted(names)}"
        )

    instrument = histogram.InstrumentParameters
    sample = histogram.SampleParameters
    histogram.clear_refinements({"Background": True})
    histogram.clear_refinements(
        {"Instrument Parameters": [k for k in REPORTED_INSTRUMENT if k in instrument]}
    )
    histogram.clear_refinements(
        {"Sample Parameters": [k for k in SAMPLE_KEYS if k in sample]}
    )
    for phase in phases:
        phase.clear_refinements({"Cell": True})
        phase.clear_refinements({"LeBail": True})
        phase.clear_refinements({"Atoms": [atom.label for atom in phase.atoms()]})
        phase.clear_HAP_refinements({"Scale": True}, [histogram])

    background = flags["background"]
    if background:
        histogram.set_refinements(
            {
                "Background": {
                    "type": background["type"],
                    "no. coeffs": background["terms"],
                    "refine": True,
                }
            }
        )
    parameters = (["Zero"] if flags["zero"] else []) + flags["instrument"]
    if parameters:
        histogram.set_refinements({"Instrument Parameters": parameters})
    sample_flags = ["Scale"] if flags["scale"] else []
    if flags["displacement"]:
        geometry = sample.get("Type", "Bragg-Brentano")
        if geometry not in DISPLACEMENT:
            raise ValueError(f"no displacement parameters for {geometry!r} geometry")
        sample_flags += DISPLACEMENT[geometry]
    if sample_flags:
        histogram.set_refinements({"Sample Parameters": sample_flags})
    for phase in phases:
        if _selected(flags["cell"], phase.name):
            phase.set_refinements({"Cell": True})
        if flags["atoms"].get(phase.name):
            phase.set_refinements({"Atoms": {"all": flags["atoms"][phase.name]}})
        if _selected(flags["le_bail"], phase.name):
            phase.set_refinements({"LeBail": True})
        if _selected(flags["phase_fractions"], phase.name):
            phase.set_HAP_refinements({"Scale": True}, [histogram])


# Refinement


class RefinementError(RuntimeError):
    """A refinement that GSAS-II reported as failed without raising."""


class _Tee(io.StringIO):
    """Copies what is printed while it is in force, still printing it."""

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self
        return self

    def __exit__(self, *exc):
        sys.stdout = self._stdout
        return False

    def write(self, text):
        self._stdout.write(text)
        return super().write(text)


def _failure_report(output):
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    flagged = [line for line in lines if FAILURE_WORDS.search(line)]
    return "\n".join((flagged or lines)[-FAILURE_LINES:])


def _refine(project):
    """Refine with the flags as set, and return the covariance data.

    GSAS-II catches the errors of a refinement and returns without saving
    it, so ``G2Project.refine`` returns normally either way. The covariance
    is cleared first so that a refinement that saved nothing shows as one
    without it; the previous covariance is then put back.

    Raises
    ------
    RefinementError
        If the refinement saved no result, with the lines of its output that
        report the problem.
    """
    covariance = project.data["Covariance"]
    previous = covariance["data"]
    covariance["data"] = {}
    with _Tee() as output:
        project.refine()
    covariance = project.data["Covariance"]
    if "Rvals" not in covariance["data"]:
        covariance["data"] = previous
        raise RefinementError(_failure_report(output.getvalue()) or "refinement failed")
    return covariance["data"]


def _esds(covariance):
    """The esd of every refined parameter, and of those derived from them."""
    esds = {}
    for name, (_, esd) in covariance.get("depSigDict", {}).items():
        if esd is not None:
            esds[name] = float(esd)
    matrix = covariance.get("covMatrix")
    for index, name in enumerate(covariance.get("varyList", [])):
        esds[name] = float(math.sqrt(abs(matrix[index][index])))
    return esds


def _stage_record(name, flags, histogram, covariance):
    rvals = covariance["Rvals"]
    esds = _esds(covariance)
    residuals = histogram.residuals
    gof = float(rvals["GOF"])
    return {
        "name": name,
        "flags": flags,
        "rwp": residuals.get("wR"),
        "rp": residuals.get("R"),
        "gof": gof,
        "reduced_chi_squared": gof**2,
        "chi_squared": float(rvals["chisq"]),
        "n_observations": int(rvals["Nobs"]),
        "n_variables": int(rvals["Nvars"]),
        "converged": rvals.get("converged"),
        # GSAS-II's "Max shft/sig": the largest shift over the whole stage,
        # not in its last cycle, over the esd.
        "max_total_shift_over_esd": rvals.get("Max shft/sig"),
        "messages": rvals.get("msg", ""),
        "parameters": {
            name: {"value": float(value), "esd": esds.get(name)}
            for name, value in zip(covariance["varyList"], covariance["variables"])
        },
        "residuals": residuals,
    }


def _final(project, histogram, covariance):
    """Instrument, sample, background and phase values after the last stage."""
    esds = _esds(covariance)
    hfx = f":{histogram.id}:"
    instrument = histogram.InstrumentParameters
    sample = histogram.SampleParameters
    background = histogram.Background[0]
    phases = []
    for phase in project.phases():
        cell, cell_esd = phase.get_cell_and_esd()
        pfx = f"{phase.id}:{histogram.id}:"
        hap = phase.getHAPvalues(histogram.name)
        phases.append(
            {
                "name": phase.name,
                "cell": cell,
                "cell_esd": cell_esd,
                "phase_fraction": {
                    "value": float(hap["Scale"][0]),
                    "esd": esds.get(pfx + "Scale"),
                },
                "weight_fraction": {
                    "value": covariance.get("depSigDict", {}).get(
                        pfx + "WgtFrac", [None]
                    )[0],
                    "esd": esds.get(pfx + "WgtFrac"),
                },
                "le_bail": bool(hap.get("LeBail", False)),
                **_broadening(phase, histogram),
            }
        )
    return {
        "instrument": {
            key: {"value": instrument[key][1], "esd": esds.get(hfx + key)}
            for key in REPORTED_INSTRUMENT
            if key in instrument
        },
        "sample": {
            key: {"value": sample[key][0], "esd": esds.get(hfx + key)}
            for key in SAMPLE_KEYS
            if key in sample
        },
        "background": {"type": background[0], "coefficients": background[3:]},
        "phases": phases,
    }


# Exports


def _export_histogram(histogram, path):
    import numpy as np

    x, observed, _, calculated, background = (
        np.ma.getdata(column) for column in histogram.data["data"][1][:5]
    )
    lower, upper = histogram.Limits("lower"), histogram.Limits("upper")
    inside = (x >= lower) & (x <= upper)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["two_theta", "observed", "calculated", "background", "difference"]
        )
        for row in zip(
            x[inside], observed[inside], calculated[inside], background[inside]
        ):
            t, yo, yc, yb = (float(value) for value in row)
            writer.writerow(
                [f"{t:.6f}", f"{yo:.8g}", f"{yc:.8g}", f"{yb:.8g}", f"{yo - yc:.8g}"]
            )


def _export_reflections(G2sc, reflections, path):
    """One phase's reflections. ``sig`` in the list is the Gaussian variance
    and ``gam`` the Lorentzian FWHM, both in centidegrees."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "h",
                "k",
                "l",
                "multiplicity",
                "d",
                "two_theta",
                "fwhm",
                "f_obs_squared",
                "f_calc_squared",
                "i_corr",
            ]
        )
        for ref in reflections["RefList"]:
            h, k, l, multiplicity, d, two_theta, sig, gam, f_obs, f_calc = (
                float(value) for value in ref[:10]
            )
            fwhm = G2sc.G2pwd.getgamFW(gam, math.sqrt(max(sig, 0.0))) / 100.0
            writer.writerow(
                [
                    round(h),
                    round(k),
                    round(l),
                    round(multiplicity),
                    f"{d:.6f}",
                    f"{two_theta:.6f}",
                    f"{fwhm:.6f}",
                    f"{f_obs:.8g}",
                    f"{f_calc:.8g}",
                    f"{float(ref[11]):.6g}",
                ]
            )


def _safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "phase"


def _export(G2sc, histogram, prefix):
    """The pattern, the reflections of each phase and the instrument file."""
    prefix = Path(prefix)
    exports = {"histogram": str(prefix.with_name(prefix.name + "_histogram.csv"))}
    _export_histogram(histogram, exports["histogram"])
    exports["reflections"] = {}
    for phase, reflections in histogram.reflections().items():
        path = prefix.with_name(f"{prefix.name}_reflections_{_safe_name(phase)}.csv")
        _export_reflections(G2sc, reflections, path)
        exports["reflections"][phase] = str(path)
    exports["instprm"] = str(prefix.with_name(prefix.name + ".instprm"))
    instrument, extra = histogram.data["Instrument Parameters"]
    with open(exports["instprm"], "w", encoding="utf-8") as handle:
        G2sc.G2fil.WriteInstprm(handle, instrument, extra, histogram.SampleParameters)
    return exports


# Actions


def _add_histogram(project, job):
    """Read the pattern, with the xrdml importer or else the .xy fallback.

    Returns the histogram and a note of how it was read.
    """
    data_file = job["data_file"]
    hint = XRDML_HINT if data_file.lower().endswith(".xrdml") else None
    try:
        histogram = project.add_powder_histogram(
            data_file, job["instprm"], fmthint=hint
        )
        return histogram, {"importer": "file", "file": data_file}
    except Exception as error:
        fallback = job.get("data_fallback")
        if not fallback:
            raise
        print(f"Import of {data_file} failed ({error!r}); reading {fallback}")
        histogram = project.add_powder_histogram(
            fallback, job["instprm"], fmthint=XY_HINT
        )
        return histogram, {
            "importer": "xy fallback",
            "file": fallback,
            "error": repr(error),
        }


def _create_project(G2sc, job):
    """A new project at ``gpx`` with the histogram and the phases linked to it.

    A phase entry may carry ``cell``, six lattice parameters that replace
    those read from its CIF.
    """
    project = G2sc.G2Project(newgpx=job["gpx"])
    histogram, source = _add_histogram(project, job)
    for entry in job["phases"]:
        phase = project.add_phase(
            entry["cif"],
            phasename=entry["name"],
            histograms=[histogram],
            fmthint=CIF_HINT,
        )
        if entry.get("cell") is not None:
            cell = [float(value) for value in entry["cell"]]
            if len(cell) != 6:
                raise ValueError(f"phase {entry['name']!r}: cell needs six values")
            general = phase.data["General"]
            general["Cell"][1:7] = cell
            general["Cell"][7] = G2sc.G2lat.calc_V(G2sc.G2lat.cell2A(cell))
    project.save()
    return project, histogram, source


def _histogram_summary(histogram):
    two_theta = histogram.getdata("X")
    instrument = histogram.InstrumentParameters
    sample = histogram.SampleParameters
    return {
        "name": histogram.name,
        "n_points": int(two_theta.size),
        "two_theta_range": [float(two_theta.min()), float(two_theta.max())],
        "limits": [histogram.Limits("lower"), histogram.Limits("upper")],
        # Each parameter is stored as [initial, current, refine]; the current
        # value is the one GSAS-II computes with.
        "instrument": {key: value[1] for key, value in instrument.items()},
        "sample": {
            "type": sample.get("Type"),
            "gonio_radius": sample.get("Gonio. radius"),
        },
    }


def _phase_summary(phase):
    general = phase.data["General"]
    return {
        "name": phase.name,
        "space_group": general["SGData"]["SpGrp"],
        "cell": phase.get_cell(),
        "histograms": list(phase.data["Histograms"]),
    }


def create(G2sc, job):
    """A new project at ``gpx``: one histogram and the phases linked to it.

    ``job`` holds ``gpx``, ``data_file``, ``instprm``, optionally
    ``data_fallback``, and ``phases``, a list of ``{"cif": path, "name": name}``
    with an optional ``cell``.
    """
    project, histogram, source = _create_project(G2sc, job)
    return {
        "gpx": project.filename,
        "source": source,
        "histogram": _histogram_summary(histogram),
        "phases": [_phase_summary(phase) for phase in project.phases()],
    }


def _only_histogram(project, name=None):
    if name is not None:
        return project.histogram(name)
    histograms = project.histograms()
    if len(histograms) != 1:
        raise ValueError(
            f"the project has {len(histograms)} histograms; name one with 'histogram'"
        )
    return histograms[0]


def refine(G2sc, job):
    """Refine a project stage by stage and export the result.

    ``job`` holds ``gpx`` and ``stages``, and optionally ``histogram`` (the
    name of the one to refine, needed only if there are several),
    ``limits`` (two theta minimum and maximum), ``cycles`` (the most least
    squares cycles a stage may take), ``broadening`` (sample size and
    microstrain to hold, see the module notes) and ``export_prefix`` (the
    path the exported files are named from). With ``data_file``,
    ``instprm`` and ``phases`` as for ``create``, the project is created
    first.

    A stage that fails is recorded with its error and ends the run; the
    exports and the final values are those of the last stage that
    succeeded.
    """
    stages = accumulate_stages(job["stages"])
    broadening = check_broadening(job.get("broadening") or {})
    if job.get("data_file"):
        project, histogram, source = _create_project(G2sc, job)
    else:
        project = G2sc.G2Project(job["gpx"])
        histogram, source = _only_histogram(project, job.get("histogram")), None
    if job.get("cycles"):
        project.set_Controls("cycles", int(job["cycles"]))
    if job.get("limits"):
        lower, upper = (float(value) for value in job["limits"])
        histogram.set_refinements({"Limits": [lower, upper]})
    if broadening:
        set_broadening(project, histogram, broadening)

    result = {
        "gpx": project.filename,
        "source": source,
        "histogram": histogram.name,
        "limits": None,
        "stages": [],
        "completed": False,
    }
    covariance = None
    # Whatever a stage raises is recorded in the result, which is written
    # regardless, rather than lost with the run; the same goes for the final
    # values and the exports.
    for stage in stages:
        try:
            apply_flags(project, histogram, stage["flags"])
            covariance = _refine(project)
        except Exception as error:  # noqa: BLE001
            result["stages"].append(
                {**stage, "error": f"{type(error).__name__}: {error}"}
            )
            break
        result["stages"].append(
            _stage_record(stage["name"], stage["flags"], histogram, covariance)
        )
    else:
        result["completed"] = True
    result["limits"] = [histogram.Limits("lower"), histogram.Limits("upper")]

    if covariance is not None:
        try:
            result["final"] = _final(project, histogram, covariance)
        except Exception as error:  # noqa: BLE001
            result["final_error"] = f"{type(error).__name__}: {error}"
        if job.get("export_prefix"):
            try:
                result["exports"] = _export(G2sc, histogram, job["export_prefix"])
            except Exception as error:  # noqa: BLE001
                result["export_error"] = f"{type(error).__name__}: {error}"
    return result


ACTIONS = {"create": create, "refine": refine}


def main(argv):
    if len(argv) != 2:
        sys.exit("usage: gsas2_driver.py job.json")
    job = json.loads(Path(argv[1]).read_text(encoding="utf-8"))
    action = ACTIONS.get(job.get("action"))
    if action is None:
        sys.exit(f"unknown action {job.get('action')!r}; known: {', '.join(ACTIONS)}")

    sys.path.append(job["gsas2_home"])
    from GSASII import GSASIIscriptable as G2sc

    result = action(G2sc, job)
    result["action"] = job["action"]
    Path(job["result"]).write_text(
        json.dumps(_plain(result), indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main(sys.argv)
