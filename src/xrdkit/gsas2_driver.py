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
``cell``, ``le_bail``, ``phase_fractions``, ``size``, ``mustrain``, ``overall_uiso``
    ``True`` for every phase or a list of phase names. ``size`` and
    ``mustrain`` refine the isotropic crystallite size and microstrain,
    making either isotropic if it was not. ``overall_uiso`` refines one
    displacement parameter for every atom of the phase: when a stage first
    switches it on, every atom is made isotropic with the Uiso the job's
    ``overall_uiso_start`` gives the phase (by name, or ``"*"``), and all
    the Uiso are constrained equal.
``atoms``
    A list of ``{"phase": name, "flags": "FXU"}``, any of F (occupancy),
    X (position) and U (displacement), for every atom of the phase.
``atom_flags``
    A list of ``{"phase": name, "labels": [...], "flags": "FXU"}``, for the
    named atoms only. Atoms that share a site and have X refined are
    constrained to move together, in the coordinates the site leaves free.
``uiso_groups``
    ``{phase: [[label, ...], ...]}``: groups of sites, each named by the
    label of any atom on it, the Uiso of every atom on a group's sites
    constrained to one value and refined, replacing the phase's groups or
    overall Uiso before (which must be switched off in the same stage,
    ``"overall_uiso": false``). The atoms keep their Uiso as the start.
``coordinates``
    ``{phase: {label: "xz", ...}}``: the coordinates to refine on the site of
    each named atom, any of x, y and z, or ``"all"`` for every one the site
    leaves free. Every atom on the site is refined, the atoms that share it
    constrained to move together; a free coordinate not named is held. A
    coordinate the site ties to another by symmetry follows it, so naming
    either frees both; one the site fixes cannot be named.
``origin``
    ``{phase: {"site": label, "axis": "z"}}``: the site whose coordinate
    along a polar axis is held, whatever else is refined, to fix the origin
    along it; ``{phase: None}`` drops it, and ``None``, for a structure with
    no origin to hold, changes nothing. A stage that would refine that
    coordinate on every site of a phase polar along it fails, since the
    origin would then float.
``occupancies``
    A list of ``{"phase": name, "sites": [label, ...], "elements": [...]}``:
    the occupancies of every atom of each named element on the named sites
    refined, with the element's content over those sites, the sum of
    multiplicity times occupancy, held at its value when the constraint is
    set up. Each element must have an atom on every site named, if need be
    added at occupancy 0 by an atom edit. With every element's content held,
    so is the total over the sites, and with it the vacancies they hold
    between them; each site's own total may change. Other atoms on the sites
    keep their occupancies. ``None``, for a structure with no exchange, adds
    no constraint.
``preferred_orientation``
    ``{phase: [h, k, l]}``, a phase name or ``"*"`` for every phase: the
    March-Dollase ratio of the phase in the histogram refined, for texture
    about the [h, k, l] axis, replacing any axis given before; ``{phase:
    None}`` drops it.
Flags carry over from stage to stage, a list growing by the new entries, a
``True`` or ``False`` replacing what was there, unless the stage has
``"reset": true``, when it starts again from nothing refined. The
constraints a stage needs replace those of the stage before whenever they
differ.

Structure
---------
A phase entry's ``atoms`` edits the structure read from its CIF within the
project, the CIF itself untouched: occupancies, positions and isotropic
displacement parameters changed, and atoms of other elements added on
existing sites with the positions and displacement parameters of the atoms
they share them with (:func:`check_atom_edits`). A refine job's
``background_start`` gives the background coefficients to start from, such
as those of an earlier refinement, and ``scale_start`` the histogram scale.

Sanity check
------------
After every stage each phase's atoms are checked (:func:`check_sanity`):
a negative Uiso or diagonal Uij, an occupancy outside 0 to 1 or a site whose
occupancies add up to more than 1, and a coordinate more than ``max_shift``
(fractional, :data:`DEFAULT_MAX_SHIFT` by default) from its position in the
reference, by default the atoms as the job found them before its first
stage. A job's ``sanity``, ``{"max_shift": ..., "reference": {phase:
[{"label": ..., "xyz": [...]}, ...]}}``, sets either. The flags are recorded
with each stage, and a new one rejects the stage (below).

Rejected, unsettled and failed stages
-------------------------------------
Each stage is recorded with a ``status``, one of :data:`STATUSES`.

``"clean"``
    Settled, and raising no sanity flag the stage before it lacked. Kept.
``"unsettled"``
    Reached ``max_passes`` without settling, and raised no new flag. It is
    recorded with its ``largest_remaining_move``, the parameter that moved
    most in its last pass and by how many esds, and is kept when the job's
    ``on_unsettled`` is ``"accept"``, the default; a job that would rather
    not keep one sets ``on_unsettled`` to ``"reject"``.
``"flagged"``
    Raised a sanity flag the stage before it did not have (for the first
    stage, one the atoms did not have when the job found them), with the
    job's ``on_flagged`` ``"accept"``. Kept.
``"rejected"``
    Raised such a flag with ``on_flagged`` ``"reject"``, the default, or did
    not settle with ``on_unsettled`` ``"reject"``. Rolled back.
``"failed"``
    Raised. Rolled back, and the run ends there.

A rejected stage is recorded, result and all, with ``rejected_because``; the
project then goes back to the state the last stage kept left it in,
parameters, atoms and constraints, and is saved so, and the stage's own
entries are dropped from the stages, so that what it alone refined is held
in every stage after it. A failed stage is recorded with its ``error``.

Each stage records, besides the parameters it refined, ``values``: every
histogram and phase value (zero and the other instrument parameters, scale
and displacement, background coefficients, and each phase's cell, phase
fraction, size and microstrain) as it leaves them, refined or held, a held
one marked ``held`` with no esd, and ``atoms``, every atom's coordinates,
occupancy and displacement parameters, a held one with no esd. The result's
``start_model`` gives each phase's atoms as the job found them.

Every run ends with a final model, whatever became of its stages: the values
and the exports are those of the last stage kept, and where no stage was
kept they are the model the job started from, computed without refining it.
The result's ``final_from`` names the stage they come from, or "job start".

Undetermined parameters
-----------------------
Each stage also records, as ``undetermined``, the refined atom parameters
the data leave undetermined (:func:`find_undetermined`): an occupancy whose
esd is more than half its allowed range, 0 to 1, and a coordinate or an
isotropic Uiso whose esd is larger than its shift from where the job found
it. The result's ``undetermined`` are those of the last stage kept.

Le Bail extraction
------------------
``G2Project.refine`` never resets the Le Bail intensities, so a phase whose
Le Bail flag is switched on in a new project would have no reflection list to
extract into. Whenever a stage switches it on for a phase, the driver first
does what the GSAS-II GUI does when asked to reset the structure factors:
every intensity is set afresh and ``le_bail_cycles`` Le Bail-only cycles are
run, refining nothing, before the least squares. The starting intensities
are random, and the generator is seeded so that a job gives the same result
each time it is run.

GSAS-II extracts the Le Bail intensities afresh within every least squares
cycle, so the function being minimised shifts under the refinement, which
then stops at the first cycle that raises chi squared, well short of the
minimum. A job's ``max_passes`` refines each stage again and again, up to
that many times, until no parameter moves by more than ``pass_tolerance``
esds from one pass to the next. The histogram scale is left out of that test
when every phase is extracted by Le Bail, since the extracted intensities
take up any change in it and it drifts freely. Each stage of a job that
gives ``max_passes``, one included, records its ``passes`` and whether it
settled; without one a stage is refined once and records neither.

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
import random
import re
import sys
from pathlib import Path

# Only importers whose format name contains the hint are tried.
XRDML_HINT = "Panalytical"
XY_HINT = "comma/tab/semicolon"
CIF_HINT = "CIF"

# The version of what this driver writes, recorded by the callers that
# keep its results; raised when the result gains or changes a field.
DRIVER_VERSION = 2

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
        "size",
        "mustrain",
        "overall_uiso",
        "atom_flags",
        "uiso_groups",
        "coordinates",
        "origin",
        "occupancies",
        "preferred_orientation",
    }
)
BOOLEAN_FLAGS = ("scale", "zero", "displacement")

# Fractional coordinates, and the word that frees every one a site allows.
AXES = "xyz"
ALL_COORDINATES = "all"

# What the sanity check after each stage looks for, the fractional shift
# from the reference beyond which a coordinate is flagged, and how far an
# occupancy or a site's total may pass 0 or 1 by rounding alone.
SANITY_KINDS = ("negative_uiso", "occupancy", "coordinate_shift")
DEFAULT_MAX_SHIFT = 0.05
OCCUPANCY_TOLERANCE = 1.0e-6
# What becomes of a stage whose sanity check raises a new flag, and of one
# not settled after max_passes: kept, or rejected and the project taken back
# to the last stage kept.
POLICIES = ("accept", "reject")
DEFAULT_ON_FLAGGED = "reject"
DEFAULT_ON_UNSETTLED = "accept"
# A stage's status: refined and kept, clean, with a new flag let through,
# or not settled; refined and rejected; failed.
STATUSES = ("clean", "flagged", "unsettled", "rejected", "failed")
# The range an occupancy may take; a refined occupancy whose esd is more
# than half of it is undetermined.
OCCUPANCY_RANGE = (0.0, 1.0)

# Atoms closer than this in every fractional coordinate share a site.
SAME_SITE = 1.0e-4
PHASE_FLAGS = (
    "cell",
    "le_bail",
    "phase_fractions",
    "size",
    "mustrain",
    "overall_uiso",
)

# What a phase's atom edits may carry: the label of the atom to change or
# add, for a new atom its element and the label of the atom whose site and
# displacement parameters it shares, and the occupancy, position and
# isotropic displacement parameter to give it.
ATOM_EDIT_KEYS = frozenset({"label", "occupancy", "type", "copy", "xyz", "uiso"})

# Le Bail-only cycles run when a stage switches Le Bail extraction on, and
# the seed of the random starting intensities GSAS-II gives the reflections.
DEFAULT_LE_BAIL_CYCLES = 10
LE_BAIL_SEED = 0

# A stage refined in passes is settled when no parameter moves by more than
# this many esds from one pass to the next.
DEFAULT_PASS_TOLERANCE = 0.1

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
# What GSAS-II prints when a refinement worked. A refinement of no cycles
# leaves no covariance to judge it by, so its output is read instead.
REFINED_MARK = "Refinement successful"


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
        "size": [],
        "mustrain": [],
        "overall_uiso": [],
        "atoms": {},
        "atom_flags": {},
        "uiso_groups": {},
        "coordinates": {},
        "origin": {},
        "occupancies": [],
        "preferred_orientation": {},
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


def _atom_flags(value, current, name):
    """Flags for named atoms: ``[{"phase", "labels", "flags"}]``, growing
    the flags each atom already has."""
    if not isinstance(value, list):
        raise TypeError(f"stage {name!r}: atom_flags must be a list")
    chosen = copy.deepcopy(current)
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != {"phase", "labels", "flags"}:
            raise ValueError(
                f"stage {name!r}: each atom_flags entry must be "
                '{"phase": ..., "labels": [...], "flags": ...}'
            )
        labels = entry["labels"]
        if (
            not isinstance(labels, list)
            or not labels
            or not all(isinstance(label, str) for label in labels)
        ):
            raise ValueError(
                f"stage {name!r}: atom_flags labels must be a list of labels"
            )
        flags = str(entry["flags"]).upper()
        if not flags or any(flag not in ATOM_FLAGS for flag in flags):
            raise ValueError(
                f"stage {name!r}: atom flags {entry['flags']!r} must be drawn from "
                f"{ATOM_FLAGS}"
            )
        phase = chosen.setdefault(str(entry["phase"]), {})
        for label in labels:
            joined = set(phase.get(label, "")) | set(flags)
            phase[label] = "".join(flag for flag in ATOM_FLAGS if flag in joined)
    return chosen


def _uiso_groups(value, current, name):
    """Uiso groups by phase, ``{phase: [[label, ...], ...]}``, each replacing
    the phase's groups before; an empty list drops them."""
    if not isinstance(value, dict):
        raise TypeError(f"stage {name!r}: uiso_groups must map phases to groups")
    groups = copy.deepcopy(current)
    for phase, phase_groups in value.items():
        if not isinstance(phase_groups, list) or not all(
            isinstance(group, list)
            and group
            and all(isinstance(label, str) for label in group)
            for group in phase_groups
        ):
            raise ValueError(
                f"stage {name!r}: the Uiso groups of {phase!r} must be lists of labels"
            )
        labels = [label for group in phase_groups for label in group]
        if len(labels) != len(set(labels)):
            raise ValueError(
                f"stage {name!r}: an atom of {phase!r} is in more than one Uiso group"
            )
        if phase_groups:
            groups[str(phase)] = [list(group) for group in phase_groups]
        else:
            groups.pop(str(phase), None)
    return groups


def _axes(value, name, label):
    """Coordinates named for a site: letters of xyz in order, or "all"."""
    text = value.lower() if isinstance(value, str) else ""
    if text == ALL_COORDINATES:
        return text
    if not text or any(axis not in AXES for axis in text):
        raise ValueError(
            f"stage {name!r}: the coordinates of {label!r} must be drawn from "
            f'xyz, or be "all"; got {value!r}'
        )
    return "".join(axis for axis in AXES if axis in text)


def _coordinates(value, current, name):
    """Coordinates to refine by site, ``{phase: {label: "xz" or "all"}}``,
    each site's growing by the new ones."""
    if not isinstance(value, dict):
        raise TypeError(f"stage {name!r}: coordinates must map phases to sites")
    chosen = copy.deepcopy(current)
    for phase, sites in value.items():
        if not isinstance(sites, dict) or not sites:
            raise ValueError(
                f"stage {name!r}: the coordinates of {phase!r} must map site "
                "labels to coordinates"
            )
        phase_sites = chosen.setdefault(str(phase), {})
        for label, axes in sites.items():
            axes = _axes(axes, name, label)
            before = phase_sites.get(str(label), "")
            if ALL_COORDINATES in (axes, before):
                phase_sites[str(label)] = ALL_COORDINATES
            else:
                phase_sites[str(label)] = "".join(
                    axis for axis in AXES if axis in axes or axis in before
                )
    return chosen


def _origin(value, current, name):
    """The site that holds the origin along a polar axis, by phase,
    ``{phase: {"site": label, "axis": "z"}}``, replacing any before; None
    drops it, and None for the whole leaves them as they were."""
    if value is None:
        return copy.deepcopy(current)
    if not isinstance(value, dict):
        raise TypeError(f"stage {name!r}: origin must map phases to a site")
    origin = copy.deepcopy(current)
    for phase, entry in value.items():
        if entry is None:
            origin.pop(str(phase), None)
            continue
        if (
            not isinstance(entry, dict)
            or set(entry) != {"site", "axis"}
            or str(entry["axis"]).lower() not in tuple(AXES)
        ):
            raise ValueError(
                f"stage {name!r}: the origin of {phase!r} must be "
                '{"site": label, "axis": "x", "y" or "z"}'
            )
        origin[str(phase)] = {
            "site": str(entry["site"]),
            "axis": str(entry["axis"]).lower(),
        }
    return origin


def _preferred_orientation(value, current, name):
    """The March-Dollase axis refined by phase, ``{phase: [h, k, l]}``,
    replacing any before; None for a phase drops it."""
    if not isinstance(value, dict):
        raise TypeError(
            f"stage {name!r}: preferred_orientation must map phases to [h, k, l]"
        )
    chosen = copy.deepcopy(current)
    for phase, hkl in value.items():
        if hkl is None:
            chosen.pop(str(phase), None)
            continue
        if (
            not isinstance(hkl, (list, tuple))
            or len(hkl) != 3
            or not all(isinstance(i, int) and not isinstance(i, bool) for i in hkl)
            or not any(hkl)
        ):
            raise ValueError(
                f"stage {name!r}: the preferred orientation axis of {phase!r} must "
                f"be [h, k, l], three whole numbers not all 0, got {hkl!r}"
            )
        chosen[str(phase)] = [int(i) for i in hkl]
    return chosen


def _names(value):
    """Whether ``value`` is a non-empty list of distinct non-empty strings."""
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(item, str) and item for item in value)
        and len(set(value)) == len(value)
    )


def _occupancies(value, current, name):
    """Occupancy constraints, ``[{"phase", "sites", "elements"}]``, growing
    by the new ones; None adds none."""
    if value is None:
        return copy.deepcopy(current)
    if not isinstance(value, list):
        raise TypeError(f"stage {name!r}: occupancies must be a list")
    chosen = copy.deepcopy(current)
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != {"phase", "sites", "elements"}:
            raise ValueError(
                f"stage {name!r}: each occupancies entry must be "
                '{"phase": ..., "sites": [...], "elements": [...]}'
            )
        if not _names(entry["sites"]) or len(entry["sites"]) < 2:
            raise ValueError(
                f"stage {name!r}: an occupancy constraint needs two or more "
                f"distinct sites, got {entry['sites']!r}"
            )
        if not _names(entry["elements"]):
            raise ValueError(
                f"stage {name!r}: an occupancy constraint needs a list of "
                f"distinct elements, got {entry['elements']!r}"
            )
        clean = {
            "phase": str(entry["phase"]),
            "sites": list(entry["sites"]),
            "elements": list(entry["elements"]),
        }
        if clean not in chosen:
            chosen.append(clean)
    return chosen


def accumulate_stages(stages):
    """The flags in force at each stage, as a list of ``{"name", "flags"}``.

    Raises
    ------
    ValueError
        If there are no stages, or a stage has an unknown key or a flag of the
        wrong kind.
    TypeError
        If a stage is not a dict, or one of its lists or mappings is of the
        wrong type.
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
        if "atom_flags" in stage:
            flags["atom_flags"] = _atom_flags(
                stage["atom_flags"], flags["atom_flags"], name
            )
        if "uiso_groups" in stage:
            flags["uiso_groups"] = _uiso_groups(
                stage["uiso_groups"], flags["uiso_groups"], name
            )
        if "coordinates" in stage:
            flags["coordinates"] = _coordinates(
                stage["coordinates"], flags["coordinates"], name
            )
        if "origin" in stage:
            flags["origin"] = _origin(stage["origin"], flags["origin"], name)
        if "occupancies" in stage:
            flags["occupancies"] = _occupancies(
                stage["occupancies"], flags["occupancies"], name
            )
        if "preferred_orientation" in stage:
            flags["preferred_orientation"] = _preferred_orientation(
                stage["preferred_orientation"], flags["preferred_orientation"], name
            )
        both = [
            phase
            for phase in flags["uiso_groups"]
            if _selected(flags["overall_uiso"], phase)
        ]
        if both:
            raise ValueError(
                f"stage {name!r}: {both} cannot have an overall Uiso and Uiso "
                'groups at once; switch overall_uiso off ("overall_uiso": false)'
            )
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


def check_atom_edits(edits, phase):
    """Check a phase's atom edits, returning them with float occupancies.

    Each edit is ``{"label": ..., "occupancy": ...}`` to change an atom's
    site fraction, or ``{"label": new, "type": element, "copy": existing,
    "occupancy": ...}`` to add an atom of another element on the site of an
    existing one, with its displacement parameters. Either may also carry
    ``xyz``, three fractional coordinates to move the atom to, and ``uiso``,
    to make it isotropic with that Uiso, as when a structure is set up from
    an earlier refinement's.

    Raises
    ------
    ValueError
        If an edit has unknown keys, no label, nothing to change, a new atom
        without its element or occupancy, an element for an existing atom, an
        occupancy outside 0 to 1, or an ``xyz`` or ``uiso`` that is not three
        or one finite numbers.
    TypeError
        If the edits are not a list of dicts.
    """
    if not isinstance(edits, list):
        raise TypeError(f"atom edits of {phase!r} must be a list")
    checked = []
    for edit in edits:
        if not isinstance(edit, dict):
            raise TypeError(f"atom edit of {phase!r} is not a dict: {edit!r}")
        unknown = set(edit) - ATOM_EDIT_KEYS
        if unknown or "label" not in edit or set(edit) == {"label"}:
            raise ValueError(
                f"atom edit of {phase!r} needs a label and something to change "
                "(occupancy, xyz, uiso), or a label, type, copy and occupancy; "
                f"got {edit!r}"
            )
        if "copy" in edit and not {"type", "occupancy"} <= set(edit):
            raise ValueError(
                f"new atom {edit['label']!r} of {phase!r} needs a type and an occupancy"
            )
        if "type" in edit and "copy" not in edit:
            raise ValueError(
                f"atom {edit['label']!r} of {phase!r}: a type is given only for a "
                "new atom, with the label of the atom to copy"
            )
        clean = {
            key: str(value)
            for key, value in edit.items()
            if key in ("label", "type", "copy")
        }
        if "occupancy" in edit:
            clean["occupancy"] = float(edit["occupancy"])
            if not 0.0 <= clean["occupancy"] <= 1.0:
                raise ValueError(
                    f"atom {clean['label']!r} of {phase!r}: occupancy must lie "
                    f"between 0 and 1, got {edit['occupancy']}"
                )
        if "xyz" in edit:
            xyz = edit["xyz"]
            if not isinstance(xyz, (list, tuple)) or len(xyz) != 3:
                raise ValueError(
                    f"atom {clean['label']!r} of {phase!r}: xyz must be three "
                    f"fractional coordinates, got {xyz!r}"
                )
            clean["xyz"] = [float(value) for value in xyz]
        if "uiso" in edit:
            clean["uiso"] = float(edit["uiso"])
        if not all(
            math.isfinite(value)
            for value in [*clean.get("xyz", []), clean.get("uiso", 0.0)]
        ):
            raise ValueError(
                f"atom {clean['label']!r} of {phase!r}: xyz and uiso must be finite"
            )
        checked.append(clean)
    return checked


def edit_atoms(G2sc, phase, edits):
    """Change occupancies, positions and Uiso and add atoms on existing
    sites, as :func:`check_atom_edits` describes, then let GSAS-II set the
    phase up again for the elements it now holds.

    Raises
    ------
    ValueError
        If an atom to change or copy is missing, a new label is taken, or a
        new position would give an atom another site symmetry, which GSAS-II
        does not work out again.
    """
    edits = check_atom_edits(edits, phase.name)
    cx, ct, cs, cia = phase.data["General"]["AtomPtrs"]
    atoms = phase.data["Atoms"]
    by_label = {atom[ct - 1]: atom for atom in atoms}
    for edit in edits:
        label = edit["label"]
        if "copy" in edit:
            if label in by_label:
                raise ValueError(f"phase {phase.name!r} already has an atom {label!r}")
            if edit["copy"] not in by_label:
                raise ValueError(
                    f"phase {phase.name!r} has no atom {edit['copy']!r} to copy"
                )
            atom = copy.deepcopy(by_label[edit["copy"]])
            atom[ct - 1] = label
            atom[ct] = edit["type"]
            atom[cia + 8] = random.randint(0, sys.maxsize)
            atoms.append(atom)
            by_label[label] = atom
        elif label not in by_label:
            raise ValueError(f"phase {phase.name!r} has no atom {label!r}")
        atom = by_label[label]
        if "occupancy" in edit:
            atom[cx + 3] = edit["occupancy"]
        if "xyz" in edit:
            site = G2sc.G2spc.SytSym(edit["xyz"], phase.data["General"]["SGData"])[0]
            if site.strip() != str(atom[cs]).strip():
                raise ValueError(
                    f"atom {label!r} of {phase.name!r} at {edit['xyz']} would sit on "
                    f"a site of symmetry {site.strip()}, not {atom[cs]}"
                )
            atom[cx : cx + 3] = edit["xyz"]
        if "uiso" in edit:
            atom[cia] = "I"
            atom[cia + 1] = edit["uiso"]
    G2sc.SetupGeneral(phase.data, None)


def _xyz_esd(phase, index, site_symmetry, esds, xinel):
    """The esds of an atom's x, y and z: of the refined shifts dAx, dAy and
    dAz, and for a coordinate its site ties to another, as GSAS-II's
    ``GetCSxinel`` codes give it, that other's esd times the multiplier; None
    for one that is fixed or not refined."""
    values = [esds.get(f"{phase.id}::dA{axis}:{index}") for axis in "xyz"]
    if xinel is None:
        return values
    codes, multipliers = xinel(site_symmetry)[:2]
    for axis in range(3):
        if values[axis] is not None or not codes[axis]:
            continue
        for other in range(3):
            if codes[other] == codes[axis] and values[other] is not None:
                ratio = abs(multipliers[axis] / multipliers[other])
                values[axis] = values[other] * ratio
                break
    return values


def _atom_table(phase, esds=None, xinel=None):
    """Every atom of the phase: label, element, position, occupancy, site
    and displacement parameters, with the esds of refined coordinates and
    Uiso. ``xinel`` is GSAS-II's ``GetCSxinel``, to give a coordinate its
    site ties to another the same esd."""
    cx, ct, cs, cia = phase.data["General"]["AtomPtrs"]
    table = []
    for index, atom in enumerate(phase.data["Atoms"]):
        isotropic = atom[cia] == "I"
        row = {
            "label": atom[ct - 1],
            "type": atom[ct],
            "xyz": list(atom[cx : cx + 3]),
            "occupancy": atom[cx + 3],
            "site_symmetry": atom[cs],
            "multiplicity": atom[cs + 1],
            "adp": atom[cia],
            "uiso": atom[cia + 1] if isotropic else None,
            "uij": None if isotropic else list(atom[cia + 2 : cia + 8]),
        }
        if esds is not None:
            row["uiso_esd"] = esds.get(f"{phase.id}::AUiso:{index}")
            row["occupancy_esd"] = esds.get(f"{phase.id}::Afrac:{index}")
            row["xyz_esd"] = _xyz_esd(phase, index, atom[cs], esds, xinel)
        table.append(row)
    return table


def _operators(phase):
    """Every symmetry operator of the phase's space group, centring and
    inversion included, as ``{"rotation", "translation"}``."""
    import numpy as np

    sg = phase.data["General"].get("SGData")
    if not sg:
        return None
    operators = []
    for centring in sg.get("SGCen", [[0, 0, 0]]):
        for rotation, translation in sg["SGOps"]:
            for sign in (1, -1) if sg.get("SGInv") else (1,):
                operators.append(
                    {
                        "rotation": (sign * np.asarray(rotation)).tolist(),
                        "translation": (
                            (sign * np.asarray(translation) + np.asarray(centring))
                            % 1.0
                        ).tolist(),
                    }
                )
    return operators


def check_uiso_start(starts):
    """Check ``overall_uiso_start``: a positive Uiso by phase name or "*"."""
    if not isinstance(starts, dict):
        raise TypeError("overall_uiso_start must map phase names to a Uiso")
    checked = {}
    for phase, value in starts.items():
        checked[str(phase)] = float(value)
        if not checked[str(phase)] > 0.0:
            raise ValueError(f"overall_uiso_start of {phase!r} must be positive")
    return checked


def check_background_start(start):
    """Check ``background_start``, ``{"type": ..., "coefficients": [...]}``,
    returning it with float coefficients, or None if not given.

    Raises
    ------
    ValueError
        If it has other keys or no coefficients.
    """
    if start is None:
        return None
    if not isinstance(start, dict) or set(start) - {"type", "coefficients"}:
        raise ValueError(
            f'background_start must be {{"type": ..., "coefficients": [...]}}, got {start!r}'
        )
    coefficients = [float(value) for value in start.get("coefficients") or []]
    if not coefficients:
        raise ValueError("background_start needs at least one coefficient")
    return {
        "type": str(start.get("type", DEFAULT_BACKGROUND["type"])),
        "coefficients": coefficients,
    }


def check_scale_start(value):
    """Check ``scale_start``, the histogram scale to start from: a positive
    number, or None if not given."""
    if value is None:
        return None
    value = float(value)
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"scale_start must be positive, got {value}")
    return value


def check_displacement_start(value):
    """Check ``displacement_start``, the specimen displacement to start
    from: a finite number, or None if not given."""
    if value is None:
        return None
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"displacement_start must be finite, got {value}")
    return value


def check_phase_fraction_start(starts):
    """Check ``phase_fraction_start``, ``{phase: fraction}``, each positive."""
    if starts is None:
        return {}
    if not isinstance(starts, dict):
        raise TypeError("phase_fraction_start must map phase names to fractions")
    checked = {}
    for phase, value in starts.items():
        value = float(value)
        if not (math.isfinite(value) and value > 0.0):
            raise ValueError(
                f"phase_fraction_start of {phase!r} must be positive, got {value}"
            )
        checked[str(phase)] = value
    return checked


def check_start_model(model):
    """Check ``start_model``, ``{phase: [atom, ...]}``, each atom with a
    ``label`` and ``xyz``: the atoms a refinement is judged against."""
    if model is None:
        return None
    if not isinstance(model, dict):
        raise TypeError("start_model must map phase names to lists of atoms")
    for phase, atoms in model.items():
        if not isinstance(atoms, list) or not all(
            isinstance(atom, dict)
            and "label" in atom
            and isinstance(atom.get("xyz"), (list, tuple))
            and len(atom["xyz"]) == 3
            for atom in atoms
        ):
            raise ValueError(
                f"start_model of {phase!r} must be a list of atoms, each with a "
                "label and xyz"
            )
    return model


def check_sanity_settings(settings):
    """Check a job's ``sanity``, returning ``{"max_shift", "reference"}``
    with the defaults filled in, the reference None if not given.

    Raises
    ------
    ValueError
        If it has other keys, ``max_shift`` is not positive, or a phase's
        reference is not a list of ``{"label": ..., "xyz": [x, y, z]}``.
    TypeError
        If it, or its reference, is not a dict.
    """
    if settings is None:
        settings = {}
    if not isinstance(settings, dict):
        raise TypeError("sanity must be a dict")
    if set(settings) - {"max_shift", "reference"}:
        raise ValueError(
            f'sanity may give only "max_shift" and "reference", got {sorted(settings)}'
        )
    max_shift = float(settings.get("max_shift", DEFAULT_MAX_SHIFT))
    if not max_shift > 0.0:
        raise ValueError(f"sanity max_shift must be positive, got {max_shift}")
    reference = settings.get("reference")
    if reference is not None:
        if not isinstance(reference, dict):
            raise TypeError("the sanity reference must map phases to atoms")
        checked = {}
        for phase, atoms in reference.items():
            if not isinstance(atoms, list) or not all(
                isinstance(atom, dict)
                and "label" in atom
                and isinstance(atom.get("xyz"), (list, tuple))
                and len(atom["xyz"]) == 3
                for atom in atoms
            ):
                raise ValueError(
                    f"the sanity reference of {phase!r} must be a list of "
                    '{"label": ..., "xyz": [x, y, z]}'
                )
            checked[str(phase)] = [
                {"label": str(atom["label"]), "xyz": [float(v) for v in atom["xyz"]]}
                for atom in atoms
            ]
        reference = checked
    return {"max_shift": max_shift, "reference": reference}


def set_overall_uiso(project, phase, start):
    """Make every atom of ``phase`` isotropic with Uiso ``start`` and
    constrain all their Uiso to be one, returning the constraints added as
    lists of variable names."""
    cia = phase.data["General"]["AtomPtrs"][3]
    atoms = phase.data["Atoms"]
    for atom in atoms:
        atom[cia] = "I"
        atom[cia + 1] = float(start)
    names = [f"{phase.id}::AUiso:{index}" for index in range(len(atoms))]
    if len(names) < 2:
        return []
    project.add_EquivConstr(names)
    return [names]


def set_uiso_groups(project, phase, groups, start=None):
    """Constrain the Uiso of every atom on each group's sites to be one,
    returning the constraints added as lists of variable names.

    Each group lists sites by the label of any atom on them, and takes in
    every atom that shares one. The atoms keep their Uiso; one that is
    anisotropic is made isotropic with ``start``.

    Raises
    ------
    ValueError
        If a label is not an atom of the phase, an atom falls in two groups,
        or an anisotropic atom is grouped with no ``start`` given.
    """
    _, ct, _, cia = phase.data["General"]["AtomPtrs"]
    atoms = phase.data["Atoms"]
    index = {atom[ct - 1]: position for position, atom in enumerate(atoms)}
    members = site_members(phase)
    added = []
    seen = set()
    for group in groups:
        missing = [label for label in group if label not in index]
        if missing:
            raise ValueError(f"phase {phase.name!r} has no atoms {missing}")
        group = list(dict.fromkeys(m for label in group for m in members[label]))
        twice = seen & set(group)
        if twice:
            raise ValueError(
                f"atoms {sorted(twice)} of {phase.name!r} are in more than one "
                "Uiso group"
            )
        seen |= set(group)
        for label in group:
            atom = atoms[index[label]]
            if atom[cia] != "I":
                if start is None:
                    raise ValueError(
                        f"atom {label!r} of {phase.name!r} is anisotropic; give "
                        "overall_uiso_start for its group's Uiso"
                    )
                atom[cia] = "I"
                atom[cia + 1] = float(start)
        if len(group) > 1:
            names = [f"{phase.id}::AUiso:{index[label]}" for label in group]
            project.add_EquivConstr(names)
            added.append(names)
    return added


def _set_uiso_scheme(project, phase, flags, uiso_start, schemes):
    """Put the Uiso constraints ``flags`` ask of ``phase`` in place, replacing
    those of another scheme set up before; a stage that asks for none leaves
    them as they are. ``schemes`` holds each phase's scheme and constraints."""
    if _selected(flags["overall_uiso"], phase.name):
        wanted = "overall"
    elif phase.name in flags["uiso_groups"]:
        wanted = tuple(tuple(group) for group in flags["uiso_groups"][phase.name])
    else:
        return
    current = schemes.get(phase.name)
    if current is not None and current[0] == wanted:
        return
    start = uiso_start.get(phase.name, uiso_start.get(ALL_PHASES))
    if current is not None:
        remove_constraints(project, current[1])
    if wanted == "overall":
        if start is None:
            raise ValueError(
                f"overall_uiso_start gives no Uiso for phase {phase.name!r}"
            )
        added = set_overall_uiso(project, phase, start)
    else:
        added = set_uiso_groups(project, phase, wanted, start)
    schemes[phase.name] = (wanted, added)


def remove_constraints(project, added):
    """Remove the phase constraints whose variables are those of each list
    of names in ``added``."""
    wanted = [set(names) for names in added]
    constraints = project.data["Constraints"]["data"]
    constraints["Phase"] = [
        constraint
        for constraint in constraints.get("Phase", [])
        if {str(term[1]) for term in constraint[:-3]} not in wanted
    ]


def atom_flag_strings(phase, flags):
    """The refinement flags of each atom of ``phase``, by label: those given
    every atom; U for an overall Uiso or an atom on a site of a Uiso group;
    X for an atom on a site whose coordinates are refined; F for an atom in
    an occupancy constraint; and those given the atom by name.

    Raises
    ------
    ValueError
        If flags, a Uiso group or a site name an atom the phase does not
        have, or an occupancy constraint cannot be set up.
    """
    name = phase.name
    ct = phase.data["General"]["AtomPtrs"][1]
    labels = [atom[ct - 1] for atom in phase.data["Atoms"]]
    members = site_members(phase)
    chosen = {label: set(flags["atoms"].get(name, "")) for label in labels}
    named = {
        label: set(value) for label, value in flags["atom_flags"].get(name, {}).items()
    }
    missing = set(named) - set(chosen)
    by_site = [(group, "U") for group in flags["uiso_groups"].get(name, [])]
    by_site.append((list(flags["coordinates"].get(name, {})), "X"))
    for group, flag in by_site:
        for label in group:
            if label not in members:
                missing.add(label)
                continue
            for member in members[label]:
                named.setdefault(member, set()).add(flag)
    if missing:
        raise ValueError(f"phase {name!r} has no atoms {sorted(missing)}")
    entries = [entry for entry in flags["occupancies"] if entry["phase"] == name]
    for constraint in occupancy_constraints(phase, entries):
        for label in constraint["atoms"]:
            named.setdefault(label, set()).add("F")
    for label in labels:
        if _selected(flags["overall_uiso"], name):
            chosen[label].add("U")
        chosen[label] |= named.get(label, set())
    return {
        label: "".join(flag for flag in ATOM_FLAGS if flag in value)
        for label, value in chosen.items()
    }


def _wrapped(value):
    """A fractional difference brought into -1/2 to 1/2."""
    return value - round(value)


def _cluster(positions):
    """The indices of ``positions`` grouped by site: each within SAME_SITE
    of the first of its group in every fractional coordinate, whole cell
    translations aside, the groups in the order of their first members."""
    groups = []
    for index, xyz in enumerate(positions):
        for group in groups:
            first = positions[group[0]]
            if all(abs(_wrapped(a - b)) < SAME_SITE for a, b in zip(xyz, first)):
                group.append(index)
                break
        else:
            groups.append([index])
    return groups


def _sites(phase):
    """Every site of the phase, as a list of (index, label, site symmetry)
    of the atoms on it, sites in the order of their first atoms."""
    cx, ct, cs, _ = phase.data["General"]["AtomPtrs"]
    atoms = phase.data["Atoms"]
    groups = _cluster([atom[cx : cx + 3] for atom in atoms])
    return [[(i, atoms[i][ct - 1], atoms[i][cs]) for i in group] for group in groups]


def site_members(phase):
    """The labels of the atoms on each atom's site, by label."""
    members = {}
    for site in _sites(phase):
        labels = [label for _, label, _ in site]
        for label in labels:
            members[label] = labels
    return members


def _element(atom_type):
    """The element of a GSAS-II atom type, such as Sr of "Sr+2"."""
    match = re.match(r"[A-Z][a-z]?", str(atom_type))
    return match.group(0) if match else str(atom_type)


def polar_axes(phase):
    """The axes, as indices, along which the phase's space group lets the
    origin float: those that every symmetry operator leaves alone, so that
    moving every atom along one changes nothing; none for a phase without
    its symmetry operators."""
    operators = _operators(phase)
    if not operators:
        return []
    return [
        axis
        for axis in range(3)
        if all(
            [row[axis] for row in op["rotation"]] == [int(j == axis) for j in range(3)]
            for op in operators
        )
    ]


def coordinate_constraints(xinel, phase, flags, refined):
    """The coordinate constraints ``flags`` ask of ``phase``.

    ``refined`` holds the labels of the atoms whose X flag is set, and
    ``xinel`` is GSAS-II's ``GetCSxinel``, whose codes give the coordinates a
    site leaves free and those it ties together. On each site holding such
    atoms every coordinate the site leaves free is refined if one of them
    has X from ``atoms`` or ``atom_flags``, and otherwise those the site's
    ``coordinates`` name; the origin's coordinate is never refined on its
    site. A coordinate the site ties to another counts as that other.

    Returns the equivalences that make the atoms sharing a site move
    together in each refined coordinate, the holds on each free coordinate
    of their atoms that is not refined, both as lists of variable names, and
    the site's record by the label of its first atom: ``{"site": labels
    joined by "/", "refined": coordinates, "held": coordinates, "origin":
    whether it holds the origin}``, naming each coordinate once, as the
    first of those the site ties together.

    Every atom on a site stays refined, since GSAS-II ignores an
    equivalence, and holds its parameters, unless all of them are varied.
    Where GSAS-II's own symmetry constraints make a parameter independent in
    one equivalence and dependent in another, it recasts them as general
    constraints, with a new variable of its own for the site's shift, and
    notes that it has done so.

    Raises
    ------
    ValueError
        If only some of the atoms on a site are refined, a label is not an
        atom of the phase, a site's coordinates or the origin name one the
        site fixes, a site's coordinates name the origin's, or every site of
        a phase polar along an axis would be refined along it.
    RuntimeError
        If ``xinel`` is None.
    """
    if xinel is None:
        raise RuntimeError("coordinates need GSAS-II's GetCSxinel")
    name = phase.name
    general = set()
    whole = flags["atoms"].get(name, "")
    for label, value in flags["atom_flags"].get(name, {}).items():
        if "X" in value:
            general.add(label)
    listed = flags["coordinates"].get(name, {})
    origin = flags["origin"].get(name)
    sites = _sites(phase)
    site_of = {
        label: position for position, site in enumerate(sites) for _, label, _ in site
    }
    missing = [label for label in listed if label not in site_of]
    if origin is not None and origin["site"] not in site_of:
        missing.append(origin["site"])
    if missing:
        raise ValueError(f"phase {name!r} has no atoms {sorted(missing)}")
    requests = {}
    for label, axes in listed.items():
        request = requests.setdefault(site_of[label], {"all": False, "axes": ""})
        if axes == ALL_COORDINATES:
            request["all"] = True
        else:
            request["axes"] += axes

    equivalences, holds, records, refined_axes = [], [], {}, {}
    for position, site in enumerate(sites):
        labels = [label for _, label, _ in site]
        chosen = [label for label in labels if label in refined]
        if not chosen:
            continue
        if len(chosen) != len(labels):
            raise ValueError(
                f"atoms {sorted(labels)} of {name!r} share a site; refine the "
                f"coordinates of all of them, not only {sorted(chosen)}"
            )
        codes = xinel(site[0][2])[0]
        independent = [
            axis for axis in range(3) if codes[axis] and codes[axis] not in codes[:axis]
        ]

        def first_of(axis, codes=codes, labels=labels, symmetry=site[0][2]):
            if not codes[axis]:
                raise ValueError(
                    f"the site of {labels[0]!r} in {name!r} ({symmetry}) fixes "
                    f"{AXES[axis]}; it cannot be refined or hold the origin"
                )
            return codes.index(codes[axis])

        request = requests.get(position, {"all": False, "axes": ""})
        if "X" in whole or general & set(labels) or request["all"]:
            wanted = set(independent)
        else:
            wanted = set()
        explicit = {first_of(AXES.index(letter)) for letter in request["axes"]}
        wanted |= explicit
        is_origin = origin is not None and origin["site"] in labels
        if is_origin:
            held = first_of(AXES.index(origin["axis"]))
            if held in explicit:
                raise ValueError(
                    f"{origin['axis']} of {origin['site']!r} in {name!r} holds the "
                    "origin and cannot be refined"
                )
            wanted.discard(held)
        for axis in sorted(wanted):
            if len(site) > 1:
                equivalences.append(
                    [f"{phase.id}::dA{AXES[axis]}:{index}" for index, _, _ in site]
                )
        for axis in independent:
            if axis not in wanted:
                holds += [
                    [f"{phase.id}::dA{AXES[axis]}:{index}"] for index, _, _ in site
                ]
        refined_axes[position] = {
            axis
            for axis in range(3)
            if codes[axis] and codes.index(codes[axis]) in wanted
        }
        records[labels[0]] = {
            "site": "/".join(labels),
            "refined": "".join(AXES[axis] for axis in sorted(wanted)),
            "held": "".join(AXES[axis] for axis in independent if axis not in wanted),
            "origin": is_origin,
        }

    for axis in polar_axes(phase):
        if all(
            position in refined_axes and axis in refined_axes[position]
            for position in range(len(sites))
        ):
            raise ValueError(
                f"the space group of {name!r} is polar along {AXES[axis]}, and "
                f"refining {AXES[axis]} on every site leaves the origin free; "
                f"hold one site's {AXES[axis]} with a stage's origin"
            )
    return equivalences, holds, records


def occupancy_constraints(phase, entries):
    """The constraint equations the occupancy constraints ``entries`` of
    ``phase`` stand for, one for each element of each, as ``{"element",
    "sites", "atoms", "names", "multipliers", "total"}``: the atoms of the
    element on the named sites, their occupancy variables and
    multiplicities, and the element's content over the sites now, the sum
    of multiplicity times occupancy, which the constraint holds.

    Raises
    ------
    ValueError
        If a site is not an atom's label, two labels are on one site, an
        element has no atom on a site named, or an atom falls in two
        constraints.
    """
    cx, ct, cs, _ = phase.data["General"]["AtomPtrs"]
    atoms = phase.data["Atoms"]
    index = {atom[ct - 1]: i for i, atom in enumerate(atoms)}
    sites = _sites(phase)
    site_of = {i: position for position, site in enumerate(sites) for i, _, _ in site}
    constraints, used = [], set()
    for entry in entries:
        missing = [label for label in entry["sites"] if label not in index]
        if missing:
            raise ValueError(f"phase {phase.name!r} has no atoms {missing}")
        positions = [site_of[index[label]] for label in entry["sites"]]
        if len(set(positions)) != len(positions):
            raise ValueError(
                f"the occupancy sites {entry['sites']} of {phase.name!r} name "
                "one site twice"
            )
        for element in entry["elements"]:
            chosen = []
            for label, position in zip(entry["sites"], positions):
                on_site = [
                    i
                    for i, _, _ in sites[position]
                    if _element(atoms[i][ct]) == element
                ]
                if not on_site:
                    raise ValueError(
                        f"phase {phase.name!r} has no {element} on the site of "
                        f"{label!r}; add one with an atom edit, at occupancy 0"
                    )
                chosen += on_site
            taken = used & set(chosen)
            if taken:
                raise ValueError(
                    f"atoms {sorted(atoms[i][ct - 1] for i in taken)} of "
                    f"{phase.name!r} are in two occupancy constraints"
                )
            used |= set(chosen)
            multipliers = [float(atoms[i][cs + 1]) for i in chosen]
            constraints.append(
                {
                    "element": element,
                    "sites": list(entry["sites"]),
                    "atoms": [atoms[i][ct - 1] for i in chosen],
                    "names": [f"{phase.id}::Afrac:{i}" for i in chosen],
                    "multipliers": multipliers,
                    "total": sum(
                        m * float(atoms[i][cx + 3]) for m, i in zip(multipliers, chosen)
                    ),
                }
            )
    return constraints


def _replace_constraints(project, managed, key, wanted, add, note):
    """Put in place the constraints ``wanted`` stands for, under ``key`` in
    ``managed``, replacing those set up there before unless they are the
    same, and return the note kept with them. ``add`` adds them, returning
    the lists of names of the constraints it added; ``note`` is kept with
    them when they are new."""
    current = managed.get(key)
    if current is not None and current[0] == wanted:
        return current[2]
    if current is not None:
        remove_constraints(project, current[1])
    managed[key] = (wanted, add() if wanted else [], note)
    return note


def set_structure_constraints(project, phase, flags, xinel, managed):
    """Put the coordinate and occupancy constraints ``flags`` ask of
    ``phase`` in place, replacing those of the stage before where they
    differ, and return what they do: ``{"coordinates": by site,
    "occupancy_constraints": [...]}``, each only if there are any. See
    :func:`coordinate_constraints` and :func:`occupancy_constraints`."""
    notes = {}
    refined = {
        label
        for label, value in atom_flag_strings(phase, flags).items()
        if "X" in value
    }
    equivalences, holds, records = (
        coordinate_constraints(xinel, phase, flags, refined)
        if refined
        else ([], [], {})
    )

    def add_coordinates():
        for names in equivalences:
            project.add_EquivConstr(names)
        for names in holds:
            project.add_HoldConstr(names)
        return [*equivalences, *holds]

    _replace_constraints(
        project,
        managed,
        (phase.name, "coordinates"),
        (tuple(map(tuple, equivalences)), tuple(map(tuple, holds))),
        add_coordinates,
        None,
    )
    if records:
        notes["coordinates"] = records

    entries = [entry for entry in flags["occupancies"] if entry["phase"] == phase.name]
    constraints = occupancy_constraints(phase, entries)

    def add_occupancies():
        for constraint in constraints:
            project.add_EqnConstr(
                constraint["total"], constraint["names"], constraint["multipliers"]
            )
        return [constraint["names"] for constraint in constraints]

    # The totals are left out of what identifies the constraints: once set
    # up, they hold them, and the occupancies they are worked out from move.
    note = _replace_constraints(
        project,
        managed,
        (phase.name, "occupancies"),
        tuple(
            (c["element"], tuple(c["names"]), tuple(c["multipliers"]))
            for c in constraints
        ),
        add_occupancies,
        [{k: v for k, v in c.items() if k != "names"} for c in constraints],
    )
    if note:
        notes["occupancy_constraints"] = note
    return notes


def check_sanity(atoms, reference=None, max_shift=DEFAULT_MAX_SHIFT):
    """What looks wrong with a phase's atoms, given as the driver's atom
    table gives them, as a list of ``{"kind", "atoms", "value",
    "message"}``, ``kind`` one of :data:`SANITY_KINDS`.

    Flagged are a negative Uiso, or diagonal Uij of an anisotropic atom; an
    occupancy outside 0 to 1, and a site whose atoms' occupancies add up to
    more than 1, beyond OCCUPANCY_TOLERANCE; and, with ``reference``, a list
    of ``{"label", "xyz"}``, a site whose coordinate has moved more than
    ``max_shift`` (fractional, whole cells aside) from the position of its
    first atom there.
    """
    flags = []

    def flag(kind, labels, value, message):
        flags.append(
            {"kind": kind, "atoms": labels, "value": value, "message": message}
        )

    for atom in atoms:
        label = atom["label"]
        if atom.get("adp", "I") == "I":
            if atom.get("uiso") is not None and atom["uiso"] < 0.0:
                flag(
                    "negative_uiso",
                    [label],
                    atom["uiso"],
                    f"{label} Uiso {atom['uiso']:.4f}",
                )
        else:
            for key, value in zip(("U11", "U22", "U33"), atom["uij"][:3]):
                if value < 0.0:
                    flag("negative_uiso", [label], value, f"{label} {key} {value:.4f}")
        occupancy = atom["occupancy"]
        if not -OCCUPANCY_TOLERANCE <= occupancy <= 1.0 + OCCUPANCY_TOLERANCE:
            flag("occupancy", [label], occupancy, f"{label} occupancy {occupancy:.4f}")
    start = {entry["label"]: entry["xyz"] for entry in reference or []}
    for group in _cluster([atom["xyz"] for atom in atoms]):
        labels = [atoms[i]["label"] for i in group]
        site = "/".join(labels)
        if len(group) > 1:
            total = sum(atoms[i]["occupancy"] for i in group)
            if total > 1.0 + OCCUPANCY_TOLERANCE:
                flag("occupancy", labels, total, f"site {site} occupancy {total:.4f}")
        first = next((i for i in group if atoms[i]["label"] in start), None)
        if first is None:
            continue
        then = start[atoms[first]["label"]]
        for axis, (now, before) in enumerate(zip(atoms[first]["xyz"], then)):
            shift = _wrapped(now - before)
            if abs(shift) > max_shift:
                flag(
                    "coordinate_shift",
                    labels,
                    shift,
                    f"{site} {AXES[axis]} moved {shift:+.4f}",
                )
    return flags


def find_undetermined(atoms, start=None):
    """The refined atom parameters of a phase that the data leave
    undetermined, as a list of ``{"atom", "parameter", "value", "esd",
    "message"}``, ``parameter`` one of occupancy, x, y, z and uiso.

    ``atoms`` is the driver's atom table with esds, a parameter with no esd
    being one not refined, and ``start`` the atoms where the job found them.
    Undetermined are an occupancy whose esd is more than half of
    :data:`OCCUPANCY_RANGE`, and a coordinate or an isotropic Uiso whose esd
    is larger than its shift from ``start`` (whole cells aside); an atom not
    in ``start`` is judged on its occupancy alone.
    """
    before = {entry["label"]: entry for entry in start or []}
    half = 0.5 * (OCCUPANCY_RANGE[1] - OCCUPANCY_RANGE[0])
    found = []

    def add(label, parameter, value, esd, message):
        found.append(
            {
                "atom": label,
                "parameter": parameter,
                "value": value,
                "esd": esd,
                "message": message,
            }
        )

    for atom in atoms:
        label = atom["label"]
        esd = atom.get("occupancy_esd")
        if esd is not None and esd > half:
            add(
                label,
                "occupancy",
                atom["occupancy"],
                esd,
                f"{label} occupancy {atom['occupancy']:.3f}, esd {esd:.3f} more "
                f"than half its range",
            )
        then = before.get(label)
        if then is None:
            continue
        esds = atom.get("xyz_esd") or [None, None, None]
        for axis, (value, esd, was) in enumerate(zip(atom["xyz"], esds, then["xyz"])):
            shift = abs(_wrapped(value - was))
            if esd is not None and esd > shift:
                add(
                    label,
                    AXES[axis],
                    value,
                    esd,
                    f"{label} {AXES[axis]} {value:.5f}, esd {esd:.5f} more than its "
                    f"shift {shift:.5f}",
                )
        esd = atom.get("uiso_esd")
        if (
            esd is not None
            and atom.get("uiso") is not None
            and then.get("uiso") is not None
        ):
            shift = abs(atom["uiso"] - then["uiso"])
            if esd > shift:
                add(
                    label,
                    "uiso",
                    atom["uiso"],
                    esd,
                    f"{label} Uiso {atom['uiso']:.5f}, esd {esd:.5f} more than its "
                    f"shift {shift:.5f}",
                )
    return found


def check_policies(job):
    """A refine job's ``on_flagged`` and ``on_unsettled``, each "accept" or
    "reject", their defaults filled in.

    Raises
    ------
    ValueError
        If either is anything else.
    """
    found = []
    for key, default in (
        ("on_flagged", DEFAULT_ON_FLAGGED),
        ("on_unsettled", DEFAULT_ON_UNSETTLED),
    ):
        value = job.get(key) or default
        if value not in POLICIES:
            raise ValueError(
                f"{key} must be one of {', '.join(POLICIES)}, not {value!r}"
            )
        found.append(value)
    return tuple(found)


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


def _broadening(phase, histogram, esds):
    """Size and microstrain, with the esd of an isotropic value refined."""
    hap = phase.data["Histograms"][histogram.name]
    pfx = f"{phase.id}:{histogram.id}:"
    broadening = {}
    for key, name in (("size", "Size"), ("mustrain", "Mustrain")):
        entry = hap[name]
        broadening[key] = {
            "type": entry[0],
            "value": entry[1][0],
            "esd": esds.get(f"{pfx}{name};i"),
            "lorentzian_fraction": entry[1][2],
        }
    return broadening


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
        *flags["size"],
        *flags["mustrain"],
        *flags["overall_uiso"],
        *flags["atoms"],
        *flags["atom_flags"],
        *flags["uiso_groups"],
        *flags["coordinates"],
        *flags["origin"],
        *(entry["phase"] for entry in flags["occupancies"]),
        *flags["preferred_orientation"],
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
        phase.clear_HAP_refinements({"Size": True}, [histogram])
        phase.clear_HAP_refinements({"Mustrain": True}, [histogram])
        phase.clear_HAP_refinements({"Pref.Ori.": True}, [histogram])

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
        by_label = atom_flag_strings(phase, flags)
        values = set(by_label.values())
        if len(values) == 1 and "" not in values:
            # Every atom alike, as GSAS-II's "all".
            phase.set_refinements({"Atoms": {"all": values.pop()}})
        elif any(by_label.values()):
            phase.set_refinements(
                {"Atoms": {label: value for label, value in by_label.items() if value}}
            )
        if _selected(flags["le_bail"], phase.name):
            phase.set_refinements({"LeBail": True})
        if _selected(flags["phase_fractions"], phase.name):
            phase.set_HAP_refinements({"Scale": True}, [histogram])
        for flag, key in (("size", "Size"), ("mustrain", "Mustrain")):
            if _selected(flags[flag], phase.name):
                phase.set_HAP_refinements(
                    {key: {"type": "isotropic", "refine": True}}, [histogram]
                )
        orientation = flags["preferred_orientation"]
        hkl = orientation.get(phase.name, orientation.get(ALL_PHASES))
        if hkl is not None:
            # GSAS-II's March-Dollase entry: model, ratio, refine, axis.
            entry = phase.data["Histograms"][histogram.name]["Pref.Ori."]
            entry[0] = "MD"
            entry[3] = list(hkl)
            phase.set_HAP_refinements({"Pref.Ori.": True}, [histogram])


# Refinement


class RefinementError(RuntimeError):
    """A refinement that GSAS-II reported as failed without raising."""


def _copy_into(source, target):
    """Make ``target`` a deep copy of ``source`` in place.

    As GSASIIscriptable's reload does, every dict, and every list of the same
    length, that both have at the same place is kept and filled in rather
    than replaced, so that the objects GSASIIscriptable hands out, which hold
    parts of a project's data, stay valid.
    """
    if isinstance(source, dict):
        for key in [key for key in target if key not in source]:
            del target[key]
        items = source.items()
    else:
        if len(source) != len(target):
            target[:] = copy.deepcopy(source)
            return
        items = enumerate(source)
    for key, value in items:
        current = target[key] if isinstance(source, list) or key in target else None
        if type(value) is type(current) and isinstance(value, (dict, list)):
            _copy_into(value, current)
        else:
            target[key] = copy.deepcopy(value)


def _snapshot(project, covariance, state):
    """The project's data and the driver's own bookkeeping, ``state``, as a
    stage kept leaves them, to go back to."""
    return {
        "data": copy.deepcopy(project.data),
        "covariance": covariance is not None,
        "state": copy.deepcopy(state),
    }


def _restore(project, snapshot):
    """Take the project back to ``snapshot`` and save it so; the covariance
    of the stage kept, if there was one, and a copy of the bookkeeping."""
    _copy_into(snapshot["data"], project.data)
    project.save()
    covariance = project.data["Covariance"]["data"] if snapshot["covariance"] else None
    return covariance, copy.deepcopy(snapshot["state"])


def _unsettled_reason(passes):
    last = passes[-1]
    reason = f"not settled in {len(passes)} passes"
    if last["max_shift_over_esd"] is not None:
        reason += (
            f": {last['parameter']} still moved {last['max_shift_over_esd']:.2f} "
            "esd in the last"
        )
    return reason


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


def _stale(histogram, covariance):
    """Whether the histogram's residuals, and with them its calculated
    pattern, are not those of the refinement's accepted parameters.

    GSAS-II keeps in the histogram the pattern of the last function
    evaluation, which in a refinement of constrained coordinates can be a
    trial step it then refused (its Rwp capped at 100), while the Rwp of the
    covariance comes from the accepted chi squared."""
    accepted = (covariance.get("Rvals") or {}).get("Rwp")
    stored = histogram.residuals.get("wR")
    if accepted is None or stored is None:
        return False
    return abs(float(stored) - float(accepted)) > 1.0e-3 * max(float(accepted), 1.0)


def _refresh_pattern(project, histogram, covariance):
    """The refinement's covariance, the pattern computed again at the
    accepted parameters first when the one the histogram holds is stale.

    Only after the last stage: computing the pattern between stages upsets
    GSAS-II's constraints for the next, and in a stage with Le Bail
    extraction on it would extract the intensities again."""
    if not _stale(histogram, covariance):
        return covariance
    kept = copy.deepcopy(covariance)
    _compute(project)
    project.data["Covariance"]["data"] = kept
    return kept


def _compute(project):
    """Compute the project's pattern without refining it, so that a run
    that kept no stage still has a model to report and export.

    GSAS-II works the calculated pattern out at the start of a refinement,
    so a refinement of no cycles saves it having moved nothing. It refines
    no parameter, so it leaves no covariance either, and an empty one comes
    back then: the values reported are the model's own, none with an esd.
    """
    controls = project.data["Controls"]["data"]
    before = controls.get("cycles")
    previous = project.data["Covariance"]["data"]
    project.set_Controls("cycles", 0)
    try:
        with _Tee() as output:
            project.refine()
        if REFINED_MARK not in output.getvalue():
            raise RefinementError(
                _failure_report(output.getvalue()) or "the pattern was not computed"
            )
        covariance = project.data["Covariance"]["data"]
        if "Rvals" in covariance:
            return covariance
        project.data["Covariance"]["data"] = previous
        return {}
    finally:
        if before is not None:
            project.set_Controls("cycles", before)


def _extract_le_bail(G2sc, project, cycles):
    """Reset the Le Bail intensities and extract them, refining nothing.

    As the GSAS-II GUI does when told to reset the structure factors: the
    project is saved with ``newLeBail`` set, so that every reflection of a Le
    Bail phase is given a fresh random intensity, and GSASIIstrMain.DoLeBail
    runs ``cycles`` extraction cycles on it and saves the project, which is
    then read back.

    Raises
    ------
    RefinementError
        If the extraction fails, with GSAS-II's message.
    """
    random.seed(LE_BAIL_SEED)
    project.data["Controls"]["data"]["newLeBail"] = True
    project.save()
    with _Tee() as output:
        ok, rvals = G2sc.G2strMain.DoLeBail(project.filename, cycles=cycles)
    project.reload()
    if not ok:
        raise RefinementError(
            rvals.get("msg") or _failure_report(output.getvalue()) or "Le Bail failed"
        )
    return {"cycles": cycles, "rwp": rvals.get("Rwp"), "gof": rvals.get("GOF")}


def _all_le_bail(project, le_bail):
    """Whether every phase is extracted by Le Bail."""
    return ALL_PHASES in le_bail or {p.name for p in project.phases()} <= set(le_bail)


def _refine_passes(project, histogram, max_passes, tolerance, ignore=()):
    """Refine again and again until no parameter moves by more than
    ``tolerance`` esds in a pass, or ``max_passes`` have been made.

    Returns the covariance of the last pass and a record of each: Rwp, GOF
    and the largest shift over esd since the pass before, with the parameter
    that made it, leaving out those in ``ignore``.
    """
    passes = []
    before = None
    for _ in range(max_passes):
        covariance = _refine(project)
        values = dict(zip(covariance["varyList"], covariance["variables"]))
        record = {
            "rwp": histogram.residuals.get("wR"),
            "gof": float(covariance["Rvals"]["GOF"]),
            "max_shift_over_esd": None,
            "parameter": None,
            "settled": False,
        }
        if before is not None:
            esds = _esds(covariance)
            shifts = {
                name: abs(value - before[name]) / esds[name]
                for name, value in values.items()
                if name in before and name not in ignore and esds.get(name)
            }
            if shifts:
                record["parameter"] = max(shifts, key=shifts.get)
                record["max_shift_over_esd"] = shifts[record["parameter"]]
            record["settled"] = max(shifts.values(), default=0.0) < tolerance
        passes.append(record)
        if record["settled"]:
            break
        before = values
    return covariance, passes


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


def _stage_values(project, histogram, esds):
    """Every histogram and phase value a later job may start from, as a
    stage leaves them, refined or held: each ``{"value", "esd", "held"}``,
    a held value with no esd; the cell with ``cell_esd`` None and
    ``cell_held`` true when no cell term was refined. The atoms are the
    stage's ``atoms``, a held atom parameter with no esd."""
    hfx = f":{histogram.id}:"

    def entry(value, name):
        esd = esds.get(name)
        return {"value": value, "esd": esd, "held": esd is None}

    instrument = histogram.InstrumentParameters
    sample = histogram.SampleParameters
    background = histogram.Background[0]
    phases = []
    for phase in project.phases():
        pfx = f"{phase.id}:{histogram.id}:"
        hap = phase.getHAPvalues(histogram.name)
        refined = any(f"{phase.id}::A{index}" in esds for index in range(6))
        cell, cell_esd = phase.get_cell_and_esd()
        if not refined:
            cell_esd = None
        broadening = {}
        histogram_data = phase.data["Histograms"][histogram.name]
        for key, name in (("size", "Size"), ("mustrain", "Mustrain")):
            value = histogram_data[name]
            broadening[key] = {
                "type": value[0],
                **entry(value[1][0], f"{pfx}{name};i"),
                "lorentzian_fraction": value[1][2],
            }
        phases.append(
            {
                "name": phase.name,
                "cell": cell,
                "cell_esd": cell_esd,
                "cell_held": not refined,
                "phase_fraction": entry(float(hap["Scale"][0]), pfx + "Scale"),
                **broadening,
            }
        )
    return {
        "instrument": {
            key: entry(instrument[key][1], hfx + key)
            for key in REPORTED_INSTRUMENT
            if key in instrument
        },
        "sample": {
            key: entry(sample[key][0], hfx + key)
            for key in SAMPLE_KEYS
            if key in sample
        },
        "background": {
            "type": background[0],
            "coefficients": [
                entry(value, f"{hfx}Back;{index}")
                for index, value in enumerate(background[3:])
            ],
        },
        "phases": phases,
    }


def _final(project, histogram, covariance, xinel=None):
    """Instrument, sample, background and phase values after the last stage,
    with each phase's atoms and symmetry operators."""
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
                "space_group": phase.data["General"].get("SGData", {}).get("SpGrp"),
                "le_bail": bool(hap.get("LeBail", False)),
                **_broadening(phase, histogram, esds),
                "atoms": _atom_table(phase, esds, xinel),
                "operators": _operators(phase),
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
    those read from its CIF, and ``atoms``, edits to the occupancies and
    atoms read from it (:func:`check_atom_edits`), made in the project only.
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
        if entry.get("atoms"):
            edit_atoms(G2sc, phase, entry["atoms"])
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
        "atoms": _atom_table(phase),
        "operators": _operators(phase),
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
    microstrain to start from, held unless a stage refines them, see the
    module notes), ``background_start`` (``{"type": ..., "coefficients":
    [...]}``, the background to start from), ``overall_uiso_start`` (the
    Uiso an ``overall_uiso`` stage starts every atom of a phase at),
    ``le_bail_cycles`` (the Le Bail-only cycles run when a stage switches Le
    Bail extraction on), ``max_passes`` and ``pass_tolerance``
    (refinements of each stage until it settles, see the module notes; one
    by default), ``scale_start`` (the histogram scale to start from),
    ``displacement_start`` (the specimen displacement to start from),
    ``phase_fraction_start`` (``{phase: fraction}`` to start from),
    ``start_model`` (``{phase: atoms}``, the model the refinement is judged
    against for undetermined parameters and recorded as the result's
    ``start_model``, by default the atoms as the job finds them),
    ``sanity`` (the sanity check's ``max_shift`` and ``reference``, see the
    module notes) and ``export_prefix`` (the path the exported files are
    named from). With ``data_file``, ``instprm`` and ``phases`` as for
    ``create``, the project is created first. ``on_flagged`` and
    ``on_unsettled``, "accept" or "reject", say what becomes of a stage whose
    sanity check raises a new flag and of one that has not settled (see the
    module notes; reject and accept by default).

    Each stage is recorded with its status, residuals, parameters, the atoms
    as it left them, its sanity flags, its undetermined parameters and,
    where it has any, the coordinates it refines and holds by site and its
    occupancy constraints. A stage that fails is recorded with its error and
    ends the run; a rejected stage is recorded with the reasons, named in the
    result's ``rejected``, and the run goes on from the last stage kept.

    Every run ends with ``final`` and, with ``export_prefix``, ``exports``:
    the model of the last stage kept, or, where no stage was kept, the one
    the job started from, computed without refining it. ``final_from`` names
    the stage they come from, or "job start".
    """
    raw = copy.deepcopy(job["stages"])
    stages = accumulate_stages(raw)
    on_flagged, on_unsettled = check_policies(job)
    broadening = check_broadening(job.get("broadening") or {})
    uiso_start = check_uiso_start(job.get("overall_uiso_start") or {})
    background_start = check_background_start(job.get("background_start"))
    scale_start = check_scale_start(job.get("scale_start"))
    displacement_start = check_displacement_start(job.get("displacement_start"))
    fraction_start = check_phase_fraction_start(job.get("phase_fraction_start"))
    given_model = check_start_model(job.get("start_model"))
    sanity = check_sanity_settings(job.get("sanity"))
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
    if background_start:
        coefficients = background_start["coefficients"]
        histogram.set_refinements(
            {
                "Background": {
                    "type": background_start["type"],
                    "no. coeffs": len(coefficients),
                    "coeffs": coefficients,
                    "refine": False,
                }
            }
        )
    if scale_start is not None:
        histogram.SampleParameters["Scale"][0] = scale_start
    if displacement_start is not None:
        for key in ("Shift", "DisplaceX"):
            if key in histogram.SampleParameters:
                histogram.SampleParameters[key][0] = displacement_start
                break
    for phase in project.phases():
        if phase.name in fraction_start:
            hap = phase.data["Histograms"][histogram.name]
            hap["Scale"][0] = fraction_start[phase.name]

    spc = getattr(G2sc, "G2spc", None)
    xinel = spc.GetCSxinel if spc else None
    start_atoms = {phase.name: _atom_table(phase) for phase in project.phases()}
    reference = sanity["reference"] or start_atoms
    # What undetermined parameters are judged against: the model given, such
    # as a structure's CIF start carried from an earlier job, else the atoms
    # as this job found them.
    start_model = given_model or start_atoms
    result = {
        "gpx": project.filename,
        "source": source,
        "histogram": histogram.name,
        "limits": None,
        "stages": [],
        "completed": False,
        "sanity": {
            "max_shift": sanity["max_shift"],
            "reference": "given" if sanity["reference"] else "job start",
        },
        "on_flagged": on_flagged,
        "on_unsettled": on_unsettled,
        "rejected": [],
        # The atoms of each phase as the job found them, before any stage.
        "start_model": start_model,
    }
    le_bail_cycles = int(job.get("le_bail_cycles") or DEFAULT_LE_BAIL_CYCLES)
    max_passes = int(job.get("max_passes") or 1)
    tolerance = float(job.get("pass_tolerance") or DEFAULT_PASS_TOLERANCE)

    def sanity_flags(atoms_by_phase):
        return [
            {**flag, "phase": name}
            for name, atoms in atoms_by_phase.items()
            for flag in check_sanity(atoms, reference.get(name), sanity["max_shift"])
        ]

    def flag_keys(flags):
        return {(flag["phase"], flag["kind"], tuple(flag["atoms"])) for flag in flags}

    covariance = None
    le_bail = []
    # The Uiso constraints in force for each phase, as the scheme that set
    # them ("overall" or the groups) and the constraints added; and the
    # coordinate and occupancy constraints, by phase and kind.
    uiso_schemes = {}
    managed = {}
    # The last stage kept, and the flags a stage's own are judged against,
    # at first those of the atoms as the job found them.
    kept = None
    kept_flags = flag_keys(sanity_flags(start_atoms))
    # The project and the bookkeeping as the last stage kept left them, to go
    # back to when a stage is rejected or fails.
    clean = _snapshot(project, covariance, (le_bail, uiso_schemes, managed))
    # Whatever a stage raises is recorded in the result, which is written
    # regardless, rather than lost with the run; the same goes for the final
    # values and the exports.
    index = 0
    while index < len(stages):
        stage = stages[index]
        index += 1
        extraction = None
        notes = {}
        try:
            apply_flags(project, histogram, stage["flags"])
            for phase in project.phases():
                _set_uiso_scheme(
                    project, phase, stage["flags"], uiso_start, uiso_schemes
                )
                notes[phase.name] = set_structure_constraints(
                    project, phase, stage["flags"], xinel, managed
                )
            if any(phase not in le_bail for phase in stage["flags"]["le_bail"]):
                extraction = _extract_le_bail(G2sc, project, le_bail_cycles)
            le_bail = stage["flags"]["le_bail"]
            if not job.get("max_passes"):
                covariance, passes = _refine(project), None
            else:
                ignore = set()
                if _all_le_bail(project, le_bail):
                    ignore.add(f":{histogram.id}:Scale")
                covariance, passes = _refine_passes(
                    project, histogram, max_passes, tolerance, ignore
                )
        except Exception as error:  # noqa: BLE001
            result["stages"].append(
                {
                    **stage,
                    "status": "failed",
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            covariance, (le_bail, uiso_schemes, managed) = _restore(project, clean)
            break
        record = _stage_record(stage["name"], stage["flags"], histogram, covariance)
        if not le_bail and _stale(histogram, covariance):
            # The histogram holds a trial step's pattern: the Rwp is taken
            # from the accepted chi squared, and the Rp, which that does not
            # give, is left out; the pattern is computed again for the
            # exports at the end.
            record["rwp"] = float(covariance["Rvals"]["Rwp"])
            record["rp"] = None
            record["residuals"] = {}
            record["residuals_stale"] = True
        # The atoms as this stage leaves them: GSAS-II refines coordinate
        # shifts, set back to zero after each refinement, so the parameters
        # alone do not give the positions.
        esds = _esds(covariance)
        record["atoms"] = {
            phase.name: _atom_table(phase, esds, xinel) for phase in project.phases()
        }
        # Every other value as this stage leaves it, refined or held, so that
        # a later job can start from any stage kept.
        record["values"] = _stage_values(project, histogram, esds)
        record["sanity"] = sanity_flags(record["atoms"])
        record["undetermined"] = [
            {**entry, "phase": name}
            for name, atoms in record["atoms"].items()
            for entry in find_undetermined(atoms, start_model.get(name))
        ]
        for key in ("coordinates", "occupancy_constraints"):
            found = {name: note[key] for name, note in notes.items() if key in note}
            if found:
                record[key] = found
        if extraction is not None:
            record["le_bail_extraction"] = extraction
        unsettled = False
        if passes is not None:
            record["passes"] = passes
            record["passes_converged"] = passes[-1]["settled"]
            unsettled = not passes[-1]["settled"]
        if unsettled:
            record["largest_remaining_move"] = {
                "parameter": passes[-1]["parameter"],
                "shift_over_esd": passes[-1]["max_shift_over_esd"],
            }
        new = [flag for flag in record["sanity"] if flag_keys([flag]) - kept_flags]
        reasons = []
        if new and on_flagged == "reject":
            reasons += [f"sanity check: {flag['message']}" for flag in new]
        if unsettled and on_unsettled == "reject":
            reasons.append(_unsettled_reason(passes))
        if reasons:
            record["status"] = "rejected"
            record["rejected_because"] = reasons
            result["stages"].append(record)
            result["rejected"].append(stage["name"])
            covariance, (le_bail, uiso_schemes, managed) = _restore(project, clean)
            # What this stage alone refined is held from here on: its entries
            # are dropped and the flags of the stages after it made again.
            dropped = {"name": stage["name"]}
            if raw[index - 1].get("reset"):
                dropped["reset"] = True
            raw[index - 1] = dropped
            stages = accumulate_stages(raw)
            continue
        record["status"] = "unsettled" if unsettled else "flagged" if new else "clean"
        result["stages"].append(record)
        kept = record
        kept_flags = flag_keys(record["sanity"])
        clean = _snapshot(project, covariance, (le_bail, uiso_schemes, managed))
    else:
        result["completed"] = True
    result["limits"] = [histogram.Limits("lower"), histogram.Limits("upper")]
    result["undetermined"] = kept["undetermined"] if kept else []
    result["final_from"] = kept["name"] if kept else "job start"
    if covariance is None:
        # Nothing was kept, so the project is back where the job found it.
        # Its pattern is computed, without refining anything, so that the
        # model the run leaves behind can still be reported and exported.
        try:
            covariance = _compute(project)
        except Exception as error:  # noqa: BLE001
            result["compute_error"] = f"{type(error).__name__}: {error}"
    elif (
        kept is not None
        and not kept["flags"].get("le_bail")
        and _stale(histogram, covariance)
    ):
        # The pattern the exports and the figure take is the model's own,
        # not a trial step's.
        try:
            covariance = _refresh_pattern(project, histogram, covariance)
        except Exception as error:  # noqa: BLE001
            result["compute_error"] = f"{type(error).__name__}: {error}"
    try:
        result["final"] = _final(project, histogram, covariance or {}, xinel)
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
