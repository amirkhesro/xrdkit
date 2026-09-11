"""Run one xrdkit job inside GSAS-II.

This script runs under GSAS-II's own Python, started by
:func:`xrdkit.gsas2.run_job` as ``python -B gsas2_driver.py job.json``. It
must not import xrdkit, which is not installed there. The job JSON names an
``action`` and its inputs, and ``gsas2_home``, the folder containing the
GSASII package; the result is written as JSON to ``job["result"]``.

Each action is a function of GSASIIscriptable and the job returning the result
dict, registered in ``ACTIONS``.
"""

import json
import sys
from pathlib import Path

# Only importers whose format name contains the hint are tried.
XRDML_HINT = "Panalytical"
XY_HINT = "comma/tab/semicolon"
CIF_HINT = "CIF"


def _plain(value):
    """``value`` with numpy scalars and arrays turned into JSON types."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


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
    ``data_fallback``, and ``phases``, a list of ``{"cif": path, "name": name}``.
    """
    project = G2sc.G2Project(newgpx=job["gpx"])
    histogram, source = _add_histogram(project, job)
    for entry in job["phases"]:
        project.add_phase(
            entry["cif"],
            phasename=entry["name"],
            histograms=[histogram],
            fmthint=CIF_HINT,
        )
    project.save()
    return {
        "gpx": project.filename,
        "source": source,
        "histogram": _histogram_summary(histogram),
        "phases": [_phase_summary(phase) for phase in project.phases()],
    }


ACTIONS = {"create": create}


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
