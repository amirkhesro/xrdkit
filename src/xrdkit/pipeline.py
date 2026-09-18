"""Stages of a refinement from Le Bail extraction to the occupancies, and
where each mode starts from.

A refinement runs in the modes of :data:`MODES`, each a job of its own that
starts from the saved result of the one before. The builders here give each
mode's stages in the form :func:`xrdkit.gsas2.build_refine_job` and the
GSAS-II driver take them: a list of ``{"name", ...flags}``, each stage
naming the flags it switches, which the driver carries over to the stages
after it (see ``gsas2_driver.accumulate_stages``). So that a stage the
driver rejects holds what it alone refined from then on, a stage names only
what it adds, and the scale, zero and displacement switches are stated in
the first stage of every mode.

:func:`lebail_stages`
    Background and scale, with Le Bail extraction; zero or displacement;
    cell; size; microstrain, as a test.
:func:`fixed_atoms_stages`
    Scale and background, extraction off; zero or displacement, with the
    cell; size, or size and microstrain; one overall Uiso; and optionally a
    March-Dollase preferred orientation.
:func:`coordinates_stages`
    The profile; the Uiso groups; the coordinates of each kind of site in
    turn, the origin held.
:func:`occupancy_stages`
    The profile with the Uiso groups; the occupancies traded within each
    exchange group.

With two phases or more, the stage that would free the histogram scale
frees the phase fractions instead and holds the histogram scale, which the
fractions would otherwise leave undetermined.

:func:`start_from_result` reads a result JSON the driver wrote and gives
the values a mode starts from, as a :class:`StartPoint`.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

from xrdkit import gsas2_driver, writeup
from xrdkit.config import AXIS_COORDINATE, ConfigError, host_elements
from xrdkit.density import parse_formula
from xrdkit.gsas2 import (
    GONIOMETER_RADIUS,
    _same_site,
    accepted_stages,
    build_refine_job,
    failure_markdown,
    log_tail,
    run_job,
    stage_statuses,
    structure_edits,
    summary_markdown,
)
from xrdkit.io import read_scan
from xrdkit.library import DEFAULT_ANIONS, StructureEntry, load_entry
from xrdkit.plotting import plot_rietveld, save_figure
from xrdkit.project import (
    Instrument,
    Project,
    Refine,
    Sample,
    StructureSpec,
    refine_settings,
    resolved_cell,
    results_dir,
)
from xrdkit.structure import composition_edits, site_setup

__all__ = [
    "CYCLES",
    "HELD_INSTRUMENT",
    "LE_BAIL_CYCLES",
    "MODES",
    "PASS_TOLERANCE",
    "START_BROADENING",
    "Inputs",
    "Options",
    "Outcome",
    "PhaseInput",
    "PhaseStart",
    "PipelineError",
    "StartPoint",
    "coordinates_stages",
    "exchange_edits",
    "fixed_atoms_stages",
    "lebail_stages",
    "mode_paths",
    "occupancy_stages",
    "resolve_inputs",
    "run_mode",
    "run_sequence",
    "start_from_result",
]

# The modes, in order, each starting from the saved result of the one before.
MODES = ("lebail", "fixed_atoms", "coordinates", "occupancies")

# The most least squares cycles of one GSAS-II refinement, and the Le Bail
# only cycles run when a stage switches extraction on.
CYCLES = 10
LE_BAIL_CYCLES = 10
# A stage refined in passes is settled when no parameter moves by more than
# this many esds from one pass to the next.
PASS_TOLERANCE = 0.1
# Where a Le Bail refinement starts from and holds until its stages reach
# them: GSAS-II's own 1 micron, with no microstrain, both Lorentzian.
START_BROADENING = {"size": 1.0, "mustrain": 0.0, "lgmix": 1.0}

# Stage names.
BACKGROUND_AND_SCALE = "background and scale"
ZERO = "zero"
DISPLACEMENT = "displacement"
CELL = "cell"
SIZE = "size"
MICROSTRAIN = "microstrain"
SCALE_AND_BACKGROUND = "scale and background"
ZERO_AND_CELL = "zero and cell"
DISPLACEMENT_AND_CELL = "displacement and cell"
SIZE_AND_MICROSTRAIN = "size and microstrain"
OVERALL_UISO = "overall Uiso"
PREFERRED_ORIENTATION = "preferred orientation"
PROFILE = "profile"
UISO_GROUPS = "Uiso groups"
PROFILE_AND_UISO = "profile and Uiso"
# Named from a kind of site, and for occupancies from the elements too when
# two exchange groups share a kind.
KIND_SITES = "{kind} sites"
KIND_OCCUPANCIES = "{kind} site occupancies"

# GSAS-II's names for the parameters a start is read from, in the first
# histogram and the first phase.
_ZERO = ":0:Zero"
_SHIFT = ":0:Shift"
_SCALE = ":0:Scale"
_BACKGROUND = ":0:Back;{index}"
_SIZE = "0:0:Size;i"
_MUSTRAIN = "0:0:Mustrain;i"
_CELL_TERMS = tuple(f"0::A{index}" for index in range(6))
_CELL_KEYS = (
    ("a", "length_a"),
    ("b", "length_b"),
    ("c", "length_c"),
    ("alpha", "angle_alpha"),
    ("beta", "angle_beta"),
    ("gamma", "angle_gamma"),
)
# Two reciprocal metric terms of the final cell this close are taken as tied
# by symmetry, the second following the first.
_TIED = 1.0e-9


class PipelineError(ValueError):
    """Stages that cannot be built as asked, or a result a mode cannot
    start from."""


# Stage builders


def _check_phases(phases: int) -> int:
    if isinstance(phases, bool) or not isinstance(phases, int) or phases < 1:
        raise PipelineError(
            f"phases must be a whole number of at least 1, not {phases!r}"
        )
    return phases


def _scale(phases: int) -> dict:
    """The histogram scale with one phase; the phase fractions, the
    histogram scale held, with more."""
    if _check_phases(phases) > 1:
        return {"scale": False, "phase_fractions": True}
    return {"scale": True}


def _background(background: Mapping) -> dict:
    """``{function, terms}`` as the driver's ``{type, terms}``."""
    if not isinstance(background, Mapping) or set(background) != {"function", "terms"}:
        raise PipelineError(
            f"background must be {{function, terms}}, not {background!r}"
        )
    terms = background["terms"]
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise PipelineError(f"background terms must be at least 1, not {terms!r}")
    return {"type": str(background["function"]), "terms": terms}


def _switches(displacement: bool) -> dict:
    """The zero shift, or the specimen displacement with the zero held."""
    return {"zero": not displacement, "displacement": bool(displacement)}


def lebail_stages(
    background: Mapping,
    phases: int = 1,
    displacement: bool = False,
    mustrain_test: bool = True,
) -> list[dict]:
    """The stages of a Le Bail extraction, each adding to the one before.

    1. ``background and scale``: the background (``{function, terms}``) and
       the histogram scale, or with two phases or more the phase fractions
       with the histogram scale held; Le Bail extraction on for every
       phase; the zero and displacement held.
    2. ``zero``, or ``displacement`` when ``displacement`` is true: that
       one freed, the other held.
    3. ``cell``: every cell parameter the space group leaves free.
    4. ``size``: the isotropic crystallite size.
    5. ``microstrain``: the isotropic microstrain, only when
       ``mustrain_test`` is true.

    No instrument parameter but the zero and no atomic parameter is freed.
    """
    stages = [
        {
            "name": BACKGROUND_AND_SCALE,
            "background": _background(background),
            **_scale(phases),
            "le_bail": True,
            "zero": False,
            "displacement": False,
        },
        {
            "name": DISPLACEMENT if displacement else ZERO,
            **_switches(displacement),
        },
        {"name": CELL, "cell": True},
        {"name": SIZE, "size": True},
    ]
    if mustrain_test:
        stages.append({"name": MICROSTRAIN, "mustrain": True})
    return stages


def _hkl(value: object) -> list[int]:
    if (
        not isinstance(value, (list, tuple))
        or len(value) != 3
        or not all(isinstance(i, int) and not isinstance(i, bool) for i in value)
        or not any(value)
    ):
        raise PipelineError(
            "preferred_orientation must be [h, k, l], three whole numbers not "
            f"all 0, not {value!r}"
        )
    return [int(i) for i in value]


def fixed_atoms_stages(
    background: Mapping,
    phases: int = 1,
    displacement: bool = False,
    mustrain: bool = False,
    preferred_orientation: Sequence[int] | None = None,
) -> list[dict]:
    """The stages of a Rietveld refinement with the atoms fixed, each adding
    to the one before.

    1. ``scale and background``: the background and the histogram scale, or
       with two phases or more the phase fractions with the histogram scale
       held; Le Bail extraction off; zero and displacement held.
    2. ``zero and cell``, or ``displacement and cell`` when
       ``displacement`` is true, the other held; the cell as the space group
       leaves it free.
    3. ``size``, or ``size and microstrain`` when ``mustrain`` is true; the
       microstrain is otherwise held where the job starts it.
    4. ``overall Uiso``: one Uiso for every atom of every phase, every atom
       made isotropic.
    5. ``preferred orientation``, only when ``preferred_orientation`` is an
       [h, k, l] triple: the March-Dollase ratio about that axis, for every
       phase.

    No coordinate or occupancy is freed.
    """
    stages = [
        {
            "name": SCALE_AND_BACKGROUND,
            "background": _background(background),
            **_scale(phases),
            "le_bail": False,
            "zero": False,
            "displacement": False,
        },
        {
            "name": DISPLACEMENT_AND_CELL if displacement else ZERO_AND_CELL,
            **_switches(displacement),
            "cell": True,
        },
        (
            {"name": SIZE_AND_MICROSTRAIN, "size": True, "mustrain": True}
            if mustrain
            else {"name": SIZE, "size": True}
        ),
        {"name": OVERALL_UISO, "overall_uiso": True},
    ]
    if preferred_orientation is not None:
        stages.append(
            {
                "name": PREFERRED_ORIENTATION,
                "preferred_orientation": {"*": _hkl(preferred_orientation)},
            }
        )
    return stages


def _plan_phase(plan: Mapping) -> str:
    phase = plan.get("phase")
    if not isinstance(phase, str) or not phase:
        raise PipelineError(
            "the plan names no phase; set its 'phase' to the GSAS-II phase name"
        )
    return phase


def _profile(
    name: str,
    background: Mapping,
    phases: int,
    displacement: bool,
    mustrain: bool,
) -> dict:
    """Background, scale or phase fractions, zero or displacement, cell and
    size, and the microstrain when ``mustrain`` is true, together."""
    stage = {
        "name": name,
        "background": _background(background),
        **_scale(phases),
        "le_bail": False,
        **_switches(displacement),
        "cell": True,
        "size": True,
    }
    if mustrain:
        stage["mustrain"] = True
    return stage


def _along(free: str, axis: str) -> bool:
    """Whether coordinates ``free``, as the plan gives them, may move a site
    along ``axis``; "all" is taken as able to."""
    return free == "all" or axis in free


def coordinates_stages(
    plan: Mapping,
    background: Mapping,
    phases: int = 1,
    displacement: bool = False,
    mustrain: bool = False,
) -> list[dict]:
    """The stages of a refinement of the coordinates, each adding to the
    one before, for the phase ``plan`` names (``plan["phase"]``), a plan as
    :func:`xrdkit.structure.site_setup` gives it.

    1. ``profile``: the background, the histogram scale (or with two phases
       or more the phase fractions, the histogram scale held), the zero (or
       the displacement, the zero held), the cell and the size, and the
       microstrain when ``mustrain`` is true.
    2. ``Uiso groups``: the plan's Uiso groups, each one Uiso over its
       sites, in place of an overall Uiso.
    3. One stage per kind of site, in the order of the plan's kinds, the
       entry's order, named ``"<kind> sites"``: the free coordinates of
       every site of that kind but the origin's. A kind with nothing to
       free has no stage. The origin site's coordinate along its axis is
       held from the stage of its kind, or from the first of these stages
       when its kind has none.

    No occupancy is freed.

    Raises
    ------
    PipelineError
        If the plan fixes no origin while its entry is polar along an axis
        and a coordinate freed lies along it, since the structure would then
        slide along that axis; or the plan names no phase.
    """
    phase = _plan_phase(plan)
    stages = [
        _profile(PROFILE, background, phases, displacement, mustrain),
        {
            "name": UISO_GROUPS,
            "overall_uiso": False,
            "uiso_groups": {phase: [list(group) for group in plan["uiso_groups"]]},
        },
    ]
    origin = plan.get("origin")
    polar = plan.get("polar_axis")
    kind_stages = []
    for kind in plan["kinds"]:
        free = {
            site: axes
            for site, axes in plan["coordinates"].get(kind, {}).items()
            if axes
        }
        if not free:
            continue
        if origin is None and polar and any(_along(a, polar) for a in free.values()):
            raise PipelineError(
                f"no origin is fixed, but the structure is polar along {polar} and "
                f"the {kind} sites {', '.join(free)} would be freed along it; fix "
                f"the origin on a site along {polar}"
            )
        kind_stages.append(
            {"name": KIND_SITES.format(kind=kind), "coordinates": {phase: free}}
        )
    if origin is not None and kind_stages:
        held = {phase: {"site": origin["name"], "axis": plan["origin_axis"]}}
        names = [stage["name"] for stage in kind_stages]
        own = KIND_SITES.format(kind=origin["kind"])
        kind_stages[names.index(own) if own in names else 0]["origin"] = held
    return stages + kind_stages


def _exchange_groups(plan: Mapping) -> list[Mapping]:
    exchange = plan.get("exchange")
    if exchange is None:
        return []
    return [exchange] if isinstance(exchange, Mapping) else list(exchange)


def occupancy_stages(
    plan: Mapping,
    background: Mapping,
    phases: int = 1,
    displacement: bool = False,
    mustrain: bool = False,
) -> list[dict]:
    """The stages of a refinement of the occupancies, each adding to the
    one before, for the phase ``plan`` names; the coordinates are held
    throughout.

    1. ``profile and Uiso``: the background, the histogram scale (or with
       two phases or more the phase fractions, the histogram scale held),
       the zero (or the displacement), the cell, the size, the microstrain
       when ``mustrain`` is true, and the plan's Uiso groups.
    2. One stage per exchange group of the plan, named ``"<kind> site
       occupancies"`` from the kind of its sites (with its elements too when
       two groups share a kind): the occupancies of its elements traded
       between its sites, each element's content over them held.

    Raises
    ------
    PipelineError
        If the plan names no phase, or an exchange group's sites are of more
        than one kind.
    """
    phase = _plan_phase(plan)
    stages = [
        {
            **_profile(PROFILE_AND_UISO, background, phases, displacement, mustrain),
            "overall_uiso": False,
            "uiso_groups": {phase: [list(group) for group in plan["uiso_groups"]]},
        }
    ]
    kind_of = {site["name"]: site["kind"] for site in plan["sites"]}
    groups = _exchange_groups(plan)
    kinds = []
    for group in groups:
        found = {kind_of[site] for site in group["sites"]}
        if len(found) != 1:
            raise PipelineError(
                f"the exchange sites {', '.join(group['sites'])} are of the kinds "
                f"{', '.join(sorted(found))}, not of one kind"
            )
        kinds.append(found.pop())
    for group, kind in zip(groups, kinds, strict=True):
        name = KIND_OCCUPANCIES.format(kind=kind)
        if kinds.count(kind) > 1:
            name += f" ({', '.join(group['elements'])})"
        stages.append(
            {
                "name": name,
                "occupancies": [
                    {
                        "phase": phase,
                        "sites": list(group["sites"]),
                        "elements": list(group["elements"]),
                    }
                ],
            }
        )
    return stages


# Where a mode starts


@dataclass(frozen=True)
class PhaseStart:
    """A phase as a stage left it: its name, cell (a, b, c in angstroms,
    alpha, beta, gamma in degrees), crystallite size in microns,
    microstrain, phase fraction and atoms."""

    name: str
    cell: dict[str, float]
    size: float
    microstrain: float
    fraction: float
    atoms: list[dict]


@dataclass(frozen=True)
class StartPoint:
    """The values a mode starts from, as a stage of an earlier result left
    them: each phase, in the result's order, and the histogram's zero shift,
    specimen displacement (None where the geometry has no single one), scale,
    background (``{function, terms, coefficients}``) and two theta limits;
    the result's ``start_model``, the atoms of each phase as its job found
    them (None for a result written before it was recorded); and the name
    of the stage taken. ``cell``, ``size``, ``microstrain`` and ``atoms``
    are those of the first phase."""

    phases: tuple[PhaseStart, ...]
    zero: float
    displacement: float | None
    scale: float
    background: dict
    limits: tuple[float, float]
    start_model: object
    stage: str

    @property
    def cell(self) -> dict[str, float]:
        return self.phases[0].cell

    @property
    def size(self) -> float:
        return self.phases[0].size

    @property
    def microstrain(self) -> float:
        return self.phases[0].microstrain

    @property
    def atoms(self) -> list[dict]:
        return self.phases[0].atoms


def _cell_terms(cell: Mapping[str, float]) -> list[float]:
    """GSAS-II's reciprocal metric terms, A0 to A5, of ``cell``."""
    a, b, c, alpha, beta, gamma = (cell[key] for _, key in _CELL_KEYS)
    cosines = [math.cos(math.radians(angle)) for angle in (alpha, beta, gamma)]
    metric = np.array(
        [
            [a * a, a * b * cosines[2], a * c * cosines[1]],
            [a * b * cosines[2], b * b, b * c * cosines[0]],
            [a * c * cosines[1], b * c * cosines[0], c * c],
        ]
    )
    g = np.linalg.inv(metric)
    return [g[0, 0], g[1, 1], g[2, 2], 2 * g[0, 1], 2 * g[0, 2], 2 * g[1, 2]]


def _cell_from_terms(terms: Sequence[float]) -> dict[str, float]:
    reciprocal = np.array(
        [
            [terms[0], terms[3] / 2, terms[4] / 2],
            [terms[3] / 2, terms[1], terms[5] / 2],
            [terms[4] / 2, terms[5] / 2, terms[2]],
        ]
    )
    g = np.linalg.inv(reciprocal)
    a, b, c = (math.sqrt(g[i, i]) for i in range(3))
    return {
        "length_a": a,
        "length_b": b,
        "length_c": c,
        "angle_alpha": math.degrees(math.acos(g[1, 2] / (b * c))),
        "angle_beta": math.degrees(math.acos(g[0, 2] / (a * c))),
        "angle_gamma": math.degrees(math.acos(g[0, 1] / (a * b))),
    }


def _value(entry: object) -> object:
    """A value as the driver records it: ``{"value", ...}`` or the value."""
    return entry["value"] if isinstance(entry, Mapping) else entry


class _Reader:
    """Fields of a result JSON, each missing one a one line PipelineError
    naming the path."""

    def __init__(self, path: Path, data: object) -> None:
        self.path = path
        self.data = data

    def get(self, *keys, within=None):
        value = self.data if within is None else within
        for key in keys:
            try:
                if isinstance(value, Mapping) and key not in value:
                    raise KeyError(key)
                value = value[key]
            except (KeyError, IndexError, TypeError):
                raise PipelineError(
                    f"{self.path}: no {'.'.join(map(str, keys))}"
                ) from None
        return value


def _earlier_stage_model(read, kept, chosen) -> dict:
    """The model of an accepted stage of a result written before the driver
    recorded every value at every stage: what the stage refined from its own
    record, the rest from the final model, which holds it unchanged unless a
    later stage refined it."""
    path, name = read.path, chosen["name"]
    later = kept[kept.index(chosen) + 1 :]
    parameters = chosen.get("parameters", {})

    def at_stage(parameter: str, final_value, label: str):
        if parameter in parameters:
            return parameters[parameter]["value"]
        if any(parameter in record.get("parameters", {}) for record in later):
            raise PipelineError(
                f"{path}: stage {name!r} did not refine {label}, which a later "
                "stage did, and the result, written before every value was "
                "recorded at every stage, does not give it"
            )
        return final_value

    final = read.get("final")
    phases = []
    for index, phase in enumerate(read.get("final", "phases")):
        final_cell = read.get("cell", within=phase)
        final_terms = _cell_terms(final_cell)
        names = [f"{index}::A{term}" for term in range(6)]
        refined = [term for term in range(6) if names[term] in parameters]
        terms = []
        for term in range(6):
            if term in refined:
                terms.append(parameters[names[term]]["value"])
                continue
            tied = next(
                (
                    other
                    for other in refined
                    if final_terms[other] != 0.0
                    and math.isclose(
                        final_terms[term], final_terms[other], rel_tol=_TIED
                    )
                ),
                None,
            )
            if tied is not None:
                terms.append(parameters[names[tied]]["value"])
            else:
                terms.append(
                    at_stage(names[term], final_terms[term], f"cell term A{term}")
                )
        phase_name = read.get("name", within=phase)
        atoms = chosen.get("atoms", {}).get(phase_name)
        if atoms is None:
            raise PipelineError(
                f"{path}: stage {name!r} records no atoms of {phase_name}"
            )
        phases.append(
            {
                "name": phase_name,
                "cell": _cell_from_terms(terms) if later else final_cell,
                "size": at_stage(
                    f"{index}:0:Size;i",
                    read.get("size", "value", within=phase),
                    "the size",
                ),
                "mustrain": at_stage(
                    f"{index}:0:Mustrain;i",
                    read.get("mustrain", "value", within=phase),
                    "the microstrain",
                ),
                "phase_fraction": at_stage(
                    f"{index}:0:Scale",
                    _value(phase.get("phase_fraction", 1.0)),
                    "the phase fraction",
                ),
                "atoms": atoms,
            }
        )
    sample = read.get("sample", within=final)
    shift = sample.get("Shift")
    return {
        "instrument": {
            "Zero": at_stage(
                _ZERO, read.get("instrument", "Zero", "value", within=final), "the zero"
            )
        },
        "sample": {
            "Scale": at_stage(
                _SCALE, read.get("Scale", "value", within=sample), "the scale"
            ),
            **(
                {}
                if shift is None
                else {"Shift": at_stage(_SHIFT, _value(shift), "the displacement")}
            ),
        },
        "background": {
            "type": read.get("background", "type", within=final),
            "coefficients": [
                at_stage(
                    _BACKGROUND.format(index=term), value, f"background term {term}"
                )
                for term, value in enumerate(
                    read.get("background", "coefficients", within=final)
                )
            ],
        },
        "phases": phases,
    }


def start_from_result(path: str | Path, stage: str | None = None) -> StartPoint:
    """The values a mode starts from, read from the result JSON at ``path``
    as its stage ``stage`` left them, by default the last stage accepted.

    The driver records every value at every stage, refined or held
    (``values`` and ``atoms``), so any accepted stage can be started from.
    A result written before it did gives the last stage kept from its final
    model, and an earlier stage from what that stage refined and the final
    model, which it can do only for a value no later stage refined.

    Raises
    ------
    PipelineError
        If there is no such file or it is not JSON, the stage was not
        accepted, or a value is missing (or, in a result written before
        every value was recorded, was refined only by a later stage); the
        message names the path and the stage or field.
    """
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"{path}: no such result file")
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PipelineError(f"{path}: not a result JSON: {error}") from None
    if not isinstance(result, Mapping):
        raise PipelineError(f"{path}: not a result JSON: not an object")
    read = _Reader(path, result)

    kept = accepted_stages(result)
    if not kept:
        raise PipelineError(f"{path}: no stage was accepted")
    names = [record.get("name") for record in kept]
    if stage is None:
        chosen = kept[-1]
    elif stage in names:
        chosen = kept[names.index(stage)]
    else:
        raise PipelineError(
            f"{path}: stage {stage!r} was not accepted; the accepted stages are "
            + ", ".join(map(str, names))
        )
    name = chosen["name"]

    recorded = "values" in chosen
    if recorded:
        model = chosen["values"]
    elif chosen is kept[-1]:
        model = read.get("final")
    else:
        model = _earlier_stage_model(read, kept, chosen)

    phases = []
    for phase in read.get("phases", within=model):
        cell = read.get("cell", within=phase)
        phases.append(
            PhaseStart(
                name=read.get("name", within=phase),
                cell={
                    key: float(read.get(gsas, within=cell)) for key, gsas in _CELL_KEYS
                },
                size=float(_value(read.get("size", within=phase))),
                microstrain=float(_value(read.get("mustrain", within=phase))),
                fraction=float(_value(phase.get("phase_fraction", 1.0))),
                atoms=list(
                    read.get("atoms", phase["name"], within=chosen)
                    if recorded
                    # A result from before the driver recorded atoms has none.
                    else phase.get("atoms", [])
                ),
            )
        )
    if not phases:
        raise PipelineError(f"{path}: no phases")
    shift = read.get("sample", within=model).get("Shift")
    coefficients = [
        float(_value(value))
        for value in read.get("background", "coefficients", within=model)
    ]
    limits = read.get("limits")
    if not isinstance(limits, list) or len(limits) != 2:
        raise PipelineError(f"{path}: no limits")

    return StartPoint(
        phases=tuple(phases),
        zero=float(_value(read.get("instrument", "Zero", within=model))),
        displacement=None if shift is None else float(_value(shift)),
        scale=float(_value(read.get("sample", "Scale", within=model))),
        background={
            "function": read.get("background", "type", within=model),
            "terms": len(coefficients),
            "coefficients": coefficients,
        },
        limits=(float(limits[0]), float(limits[1])),
        start_model=result.get("start_model"),
        stage=name,
    )


# Running the modes


@dataclass(frozen=True)
class Options:
    """What a run takes besides the project file: the start cell of the
    first phase (``{a, b, c, alpha, beta, gamma}``) and the start zero, over
    those the project gives; whether the specimen displacement is refined
    in place of the zero (by default for a pellet and not for a powder);
    whether the microstrain is refined in the Rietveld modes; the [h, k, l]
    axis of a March-Dollase preferred orientation for the fixed atoms mode;
    a cap on the passes of every stage, over the project's by mode; and the
    folder every file is written to in place of ``results/lebail/<key>`` and
    ``results/rietveld/<key>``."""

    cell: Mapping[str, float] | None = None
    zero: float | None = None
    displacement: bool | None = None
    mustrain: bool = False
    preferred_orientation: Sequence[int] | None = None
    max_passes: int | None = None
    out: Path | str | None = None


@dataclass(frozen=True)
class PhaseInput:
    """A phase of a sample as the project gives it: its structure key, which
    names it in GSAS-II, the structure, its library entry if it names one,
    the composition in atoms per formula unit, the formula units per cell,
    and the cell a Le Bail refinement starts from (None for the CIF's), with
    where it came from."""

    key: str
    spec: StructureSpec
    entry: StructureEntry | None
    composition: dict[str, float]
    z: int | None
    cell: dict[str, float] | None
    cell_source: str


@dataclass(frozen=True)
class Inputs:
    """Everything a mode takes from the project file for a sample."""

    project: Project
    sample: Sample
    instrument: Instrument
    instprm: Path
    scan_range: tuple[float, float]
    limits: tuple[float, float]
    refine: Refine
    phases: tuple[PhaseInput, ...]
    displacement: bool
    zero: float
    start_cell_source: str


@dataclass(frozen=True)
class Outcome:
    """What a mode left: the names of its accepted stages, its final model,
    its residuals (Rwp, Rp and chi squared of the last stage kept, and the
    weight fractions with esds when there are two phases or more), the
    parameters left undetermined, the paths it wrote, and the error that
    stopped it, None when it finished."""

    mode: str
    accepted: list[str]
    final: dict | None
    residuals: dict
    undetermined: list
    paths: dict[str, str]
    error: str | None = None


def _noop(text: str) -> None:
    """A reporter that reports nothing."""


def _sample_of(project: Project, sample: Sample | str) -> Sample:
    key = sample.key if isinstance(sample, Sample) else sample
    if key not in project.samples:
        raise PipelineError(
            f"no sample {key!r} in the project; there are "
            + (", ".join(project.samples) or "none")
        )
    return project.samples[key]


def _plural(items: Sequence, one: str, many: str) -> str:
    return one if len(items) == 1 else many


def _instprm_values(path: Path) -> dict[str, str]:
    """The ``key:value`` lines of a GSAS-II instrument parameter file."""
    values = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip()
    return values


def _six(cell: Mapping[str, float]) -> list[float]:
    return [float(cell[key]) for key, _ in _CELL_KEYS]


def _check_cell(cell: Mapping[str, float], where: str) -> dict[str, float]:
    names = [key for key, _ in _CELL_KEYS]
    if not isinstance(cell, Mapping) or set(cell) != set(names):
        raise PipelineError(f"{where} must give exactly {', '.join(names)}")
    checked = {}
    for name in names:
        value = cell[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
            raise PipelineError(
                f"{where}.{name} must be a positive number, not {value!r}"
            )
        checked[name] = float(value)
    return checked


def _lattice_column(row: Mapping[str, str], key: str) -> str:
    """A cell parameter of a lattice results row: under the name the lattice
    command writes, such as a_angstrom or alpha_deg, or where a row lacks
    that column, the plain name of older files, such as a or alpha."""
    unit = "angstrom" if key in ("a", "b", "c") else "deg"
    written = f"{key}_{unit}"
    return row[written] if written in row else row[key]


def _lattice_cell(project: Project, sample: Sample) -> tuple[dict, Path] | None:
    """The cell of the last row of the sample's lattice results, if any."""
    folder = results_dir(project, "lattice", sample.key)
    for path in (folder / f"lattice_{sample.key}.csv", folder / "lattice.csv"):
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        try:
            cell = {key: float(_lattice_column(rows[-1], key)) for key, _ in _CELL_KEYS}
        except (KeyError, TypeError, ValueError):
            raise PipelineError(f"{path}: its last row gives no whole cell") from None
        return cell, path
    return None


def _start_cell(
    project: Project, sample: Sample, spec: StructureSpec, options: Options | None
) -> tuple[dict[str, float] | None, str]:
    """The cell a phase starts from and where it came from: the options',
    then the sample's lattice results (both for the first phase only), then
    the structure's, else the CIF's."""
    if options is not None and options.cell is not None:
        return _check_cell(options.cell, "options.cell"), "options"
    if options is not None:
        found = _lattice_cell(project, sample)
        if found is not None:
            cell, path = found
            return cell, f"lattice results {path}"
    if spec.cell is not None:
        try:
            return resolved_cell(spec), f"structures.{spec.key}.cell"
        except ValueError:
            pass
    if spec.cif is not None:
        return None, f"the CIF {spec.cif}"
    raise PipelineError(
        f"structures.{spec.key}: gives no cell and no CIF to start from"
    )


def _check_unplaced(
    spec: StructureSpec, entry, composition: Mapping[str, float]
) -> None:
    """Every element of the composition must be one the entry's prototype
    carries or one the structure's atoms place."""
    if entry is None:
        return
    held = {element for site in entry.sites for element in site.elements}
    for site in spec.atoms or ():
        held.update(site["atoms"].values())
    unplaced = [element for element in composition if element not in held]
    if unplaced:
        labels = ", ".join(site.label for site in entry.sites)
        raise PipelineError(
            f"structures.{spec.key}: {', '.join(unplaced)} of the composition "
            f"{_plural(unplaced, 'is', 'are')} not in the {entry.name} prototype and "
            f"no atoms place {_plural(unplaced, 'it', 'them')}; place "
            f"{_plural(unplaced, 'it', 'them')} on one of its sites {labels}"
        )


def resolve_inputs(
    project: Project, sample: Sample | str, options: Options | None = None
) -> Inputs:
    """What a run of ``sample`` takes from the project file: its scan and
    the scan's range, its instrument and instrument parameter file, its
    phases in the order of its structures, its refine settings with the two
    theta range clipped to the scan, the start cell of the first phase and
    where it came from, whether the displacement is refined and the start
    zero.

    Raises
    ------
    PipelineError
        If the instrument gives no instrument parameter file, the scan
        cannot be read or the range lies outside it, a cell given is not a
        whole one, a structure has no cell or CIF, an element of a
        composition is on no site, or the start zero is not known.
    """
    options = options or Options()
    sample = _sample_of(project, sample)
    instrument = project.instruments[sample.instrument]
    if instrument.instprm is None:
        raise PipelineError(
            f"instruments.{instrument.key}: gives no instprm, which a GSAS-II "
            "refinement needs"
        )
    try:
        scan = read_scan(sample.file)
    except (OSError, ValueError) as error:
        raise PipelineError(f"{sample.file}: cannot be read: {error}") from None
    scan_range = (float(np.min(scan.two_theta)), float(np.max(scan.two_theta)))
    refine = refine_settings(project, sample)
    if refine.two_theta is None:
        limits = scan_range
    else:
        limits = (
            max(refine.two_theta[0], scan_range[0]),
            min(refine.two_theta[1], scan_range[1]),
        )
        if not limits[0] < limits[1]:
            raise PipelineError(
                f"samples.{sample.key}: two_theta {list(refine.two_theta)} lies "
                f"outside the scan's {scan_range[0]:.2f} to {scan_range[1]:.2f} degrees"
            )

    phases = []
    for index, key in enumerate(sample.structures):
        spec = project.structures[key]
        entry = load_entry(spec.library) if spec.library else None
        composition = parse_formula(spec.composition)
        _check_unplaced(spec, entry, composition)
        cell, source = _start_cell(
            project, sample, spec, options if index == 0 else None
        )
        phases.append(
            PhaseInput(
                key=key,
                spec=spec,
                entry=entry,
                composition=composition,
                z=spec.z if spec.z is not None else (entry.z if entry else None),
                cell=cell,
                cell_source=source,
            )
        )

    if options.zero is not None:
        zero = float(options.zero)
    else:
        try:
            zero = float(_instprm_values(instrument.instprm)["Zero"])
        except (KeyError, ValueError):
            raise PipelineError(
                f"{instrument.instprm}: gives no Zero to start from"
            ) from None
    return Inputs(
        project=project,
        sample=sample,
        instrument=instrument,
        instprm=instrument.instprm,
        scan_range=scan_range,
        limits=(float(limits[0]), float(limits[1])),
        refine=refine,
        phases=tuple(phases),
        displacement=(
            sample.form == "pellet"
            if options.displacement is None
            else bool(options.displacement)
        ),
        zero=zero,
        start_cell_source=phases[0].cell_source,
    )


def mode_paths(
    project: Project, sample: Sample | str, mode: str, options: Options | None = None
) -> dict[str, Path]:
    """Where a mode's files go: ``folder``, ``results/lebail/<key>`` for the
    Le Bail mode and ``results/rietveld/<key>`` for the others, or
    ``options.out``; the ``prefix`` ``<folder>/<key>_<mode>`` the exports
    are named from; the ``gpx`` and ``result`` JSON beside them, named with
    ``with_name`` so that a dot in the key stays; and the ``work`` folder
    ``<folder>/work/<mode>`` holding the jobs, the logs and the start
    instrument parameter file."""
    if mode not in MODES:
        raise PipelineError(f"no mode {mode!r}; the modes are {', '.join(MODES)}")
    options = options or Options()
    sample = _sample_of(project, sample)
    if options.out is not None:
        folder = Path(options.out)
    else:
        command = "lebail" if mode == "lebail" else "rietveld"
        folder = results_dir(project, command, sample)
    prefix = folder / f"{sample.key}_{mode}"
    work = folder / "work" / mode
    return {
        "folder": folder,
        "prefix": prefix,
        "gpx": prefix.with_name(prefix.name + ".gpx"),
        "result": prefix.with_name(prefix.name + "_result.json"),
        "work": work,
        "log": work / "refine.log",
        "start_instprm": work / f"{sample.key}_{mode}_start.instprm",
    }


def _write_start_instprm(
    source: Path, zero: float, path: Path, radius_mm: float | None = None
) -> Path:
    """The instrument parameter file with the zero to start from, and the
    goniometer radius where the instrument gives one.

    GSAS-II's text pattern importer carries no radius and leaves it at its
    default 200 mm, so the file every job reads names the project's instead;
    the .xrdml importer's own value is replaced by the same line. An
    instrument without a radius leaves the file as it found it."""
    lines = Path(source).read_text(encoding="utf-8").splitlines()
    written = False
    for index, line in enumerate(lines):
        if line.startswith("Zero:"):
            lines[index] = f"Zero:{zero!r}"
            written = True
    if not written:
        lines.append(f"Zero:{zero!r}")
    if radius_mm is not None:
        lines = [line for line in lines if not line.startswith(f"{GONIOMETER_RADIUS}:")]
        lines.append(f"{GONIOMETER_RADIUS}:{radius_mm}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


_SYMBOL_PART = re.compile(r"-?\d(?:_\d)?(?:/[a-z])?|[a-z]")


def _spaced_symbol(symbol: str) -> str:
    """A Hermann-Mauguin symbol with its parts spaced, as CIFs write it:
    P4/mbm as P 4/m b m."""
    return " ".join([symbol[0], *_SYMBOL_PART.findall(symbol[1:])])


def _entry_cif(phase: PhaseInput, folder: Path) -> Path:
    """A CIF for a Le Bail refinement of a library entry that has none: its
    space group and cell, and one atom of the first site's element at the
    origin, since Le Bail extraction needs no structure."""
    entry = phase.entry
    cell = phase.cell
    element = next((site.elements[0] for site in entry.sites if site.elements), "C")
    text = "\n".join(
        [
            f"data_{re.sub(r'[^A-Za-z0-9_]', '_', phase.key)}",
            *(f"_cell_{gsas} {cell[key]!r}" for key, gsas in _CELL_KEYS),
            f"_symmetry_space_group_name_H-M '{_spaced_symbol(entry.space_group)}'",
            "loop_",
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
            "_atom_site_U_iso_or_equiv",
            f"{element}1 {element} 0.0 0.0 0.0 1.0 0.01",
            "",
        ]
    )
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{phase.key}_entry.cif"
    path.write_text(text, encoding="utf-8")
    return path


def _need_cifs(inputs: Inputs, mode: str) -> None:
    lacking = [phase.key for phase in inputs.phases if phase.spec.cif is None]
    if lacking:
        raise PipelineError(
            f"the {mode} mode needs the coordinates of a CIF; structures "
            f"{', '.join(lacking)} give none (add cif beside library)"
        )


def _create(inputs: Inputs, instprm: Path, work: Path, edits=None) -> dict:
    """The phases as GSAS-II reads them from their CIFs, with ``edits`` by
    phase name, by name."""
    phases = []
    for phase in inputs.phases:
        entry = {"cif": phase.spec.cif, "name": phase.key}
        if edits and edits.get(phase.key):
            entry["atoms"] = edits[phase.key]
        phases.append(entry)
    job = build_refine_job(
        work / "structure.gpx",
        [{"name": "scale", "scale": True}],
        data_file=inputs.sample.file,
        instprm=instprm,
        phases=phases,
    )
    job["action"] = "create"
    created = run_job(job, work)
    return {phase["name"]: phase for phase in created.get("phases", [])}


def _element_of(atom_type: object) -> str:
    match = re.match(r"[A-Z][a-z]?", str(atom_type))
    return match.group(0) if match else str(atom_type)


def _added_rule(
    phase: PhaseInput, cif_atoms: Sequence[Mapping]
) -> dict[str, dict[str, str]]:
    """The composition rule of a structure's atoms placement, by host atom:
    an element a site names by a label the CIF lacks goes, under that label,
    beside the CIF atom of its host element on each site the placement names
    it on, and on no other. The host element is the one element its sites
    have in common (:func:`xrdkit.config.host_elements`); on a single site
    any CIF atom of it hosts it."""
    labels = {str(atom["label"]): _element_of(atom["type"]) for atom in cif_atoms}
    held = set(labels.values())
    where = f"structures.{phase.key}.atoms"
    sites = {site.get("label") or site["name"]: site for site in phase.spec.atoms or ()}
    absent_on = {
        name: {
            label: element
            for label, element in site["atoms"].items()
            if label not in labels and element not in held
        }
        for name, site in sites.items()
    }
    added = {element for absent in absent_on.values() for element in absent.values()}
    try:
        hosts = host_elements(
            {name: site["atoms"] for name, site in sites.items()}, added, where
        )
    except ConfigError as error:
        raise PipelineError(str(error)) from None
    rule: dict[str, dict[str, str]] = {}
    for name, site in sites.items():
        present = [label for label in site["atoms"] if label in labels]
        if absent_on[name] and not present:
            raise PipelineError(
                f"{where}: site {name} names no atom of the CIF, so "
                f"{', '.join(absent_on[name].values())} has nothing to go beside"
            )
        for label, element in absent_on[name].items():
            host = hosts[element]
            if host is None:
                # Placed on this site alone, the element takes the whole of its
                # content here whichever atom hosts it, so the first will do.
                rule.setdefault(element, {})[present[0]] = label
                continue
            beside = [atom for atom in present if labels[atom] == host]
            if len(beside) != 1:
                raise PipelineError(
                    f"{where}: site {name} must name exactly one {host} atom of "
                    f"the CIF to host {element}, not {len(beside)}"
                )
            rule.setdefault(element, {})[beside[0]] = label
    return rule


def _check_placement(
    phase: PhaseInput, cif_atoms: Sequence[Mapping], edits: Sequence[Mapping]
) -> None:
    """Every atom on a site the structure's atoms placement names, the CIF's
    or one ``edits`` add beside them, must be of an element the placement
    puts on that site."""
    by_label = {str(atom["label"]): atom for atom in cif_atoms}
    on_position = [
        (label, _element_of(atom["type"]), atom.get("xyz"))
        for label, atom in by_label.items()
    ] + [
        (edit["label"], _element_of(edit["type"]), by_label[edit["copy"]].get("xyz"))
        for edit in edits
        if "copy" in edit
    ]
    for site in phase.spec.atoms or ():
        name = site.get("label") or site["name"]
        present = [label for label in site["atoms"] if label in by_label]
        if not present or by_label[present[0]].get("xyz") is None:
            continue
        position = by_label[present[0]]["xyz"]
        placed = set(site["atoms"].values())
        stray = [
            f"{element} ({label})"
            for label, element, xyz in on_position
            if xyz is not None and _same_site(xyz, position) and element not in placed
        ]
        if stray:
            raise PipelineError(
                f"structures.{phase.key}.atoms: site {name} holds "
                f"{', '.join(stray)}, which atoms does not place there; name "
                "every element on the site in its atoms"
            )


def _nominal_edits(phase: PhaseInput, cif_atoms: Sequence[Mapping]) -> list[dict]:
    """Atom edits that put the nominal composition on the CIF's sites."""
    if phase.z is None:
        raise PipelineError(
            f"structures.{phase.key}: gives no z, which putting its composition on "
            "the sites needs; add z"
        )
    structure = {
        "name": phase.key,
        "formula_units": phase.z,
        "composition": {"added": _added_rule(phase, cif_atoms)},
    }
    try:
        edits = composition_edits(cif_atoms, phase.composition, structure)
    except ValueError as error:
        raise PipelineError(f"structures.{phase.key}: {error}") from None
    _check_placement(phase, cif_atoms, edits)
    return edits


def _mean_ueq(atoms: Sequence[Mapping], cell: Mapping[str, float]) -> float:
    """The mean equivalent isotropic displacement parameter of ``atoms``,
    weighted by multiplicity times occupancy, Ueq from the Uij on the
    crystal axes of ``cell`` (GSAS-II's cell keys)."""
    metric = _metric(cell)
    reciprocal = np.linalg.inv(metric)
    lengths = np.sqrt(np.diag(reciprocal))
    total = weight_sum = 0.0
    for atom in atoms:
        weight = float(atom.get("multiplicity", 1)) * float(atom.get("occupancy", 1.0))
        if atom.get("adp", "I") == "I" and atom.get("uiso") is not None:
            ueq = float(atom["uiso"])
        elif atom.get("uij"):
            u11, u22, u33, u12, u13, u23 = (float(value) for value in atom["uij"])
            u = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
            scaled = u * np.outer(lengths, lengths)
            ueq = float(np.sum(scaled * metric)) / 3.0
        else:
            continue
        total += weight * ueq
        weight_sum += weight
    return total / weight_sum if weight_sum else 0.01


def _with_start_uiso(atoms: Sequence[Mapping], start: float) -> list[dict]:
    """``atoms`` with the Uiso the fixed_atoms mode starts each one at: its
    own where it has one, else ``start``, the Uiso the driver gives an atom
    when it makes it isotropic. A CIF with anisotropic Uij carries no Uiso,
    and an atom recorded without one has nothing for the undetermined check
    to judge its refined Uiso against."""
    return [
        {**atom, "uiso": float(start if atom.get("uiso") is None else atom["uiso"])}
        for atom in atoms
    ]


def _metric(cell: Mapping[str, float]) -> np.ndarray:
    a, b, c, alpha, beta, gamma = (float(cell[gsas]) for _, gsas in _CELL_KEYS)
    ca, cb, cg = (math.cos(math.radians(angle)) for angle in (alpha, beta, gamma))
    return np.array(
        [
            [a * a, a * b * cg, a * c * cb],
            [a * b * cg, b * b, b * c * ca],
            [a * c * cb, b * c * ca, c * c],
        ]
    )


def _structure_table(phase: PhaseInput, atoms: Sequence[Mapping]) -> dict:
    """The structure table :func:`xrdkit.structure.site_setup` takes for a
    phase: its sites with the atoms on each, from its library entry and its
    atoms placement or from its CIF sites; its Uiso groups; its origin."""
    spec, entry = phase.spec, phase.entry
    present = {str(atom["label"]): _element_of(atom["type"]) for atom in atoms}
    given = {site.get("label") or site["name"]: site for site in spec.atoms or ()}

    def site_atoms(name: str, placed: Mapping[str, str] | None) -> dict[str, str]:
        if placed:
            found = {label: el for label, el in placed.items() if label in present}
        elif name in present:
            found = {name: present[name]}
        else:
            found = {}
        if not found:
            raise PipelineError(
                f"structures.{spec.key}: site {name} has no atom in the structure; "
                "name its atoms, by their CIF labels, in atoms"
            )
        return found

    if entry is not None:
        sites = []
        for site in entry.sites:
            placed = given[site.label]["atoms"] if site.label in given else None
            found = site_atoms(site.label, placed)
            sites.append(
                {
                    "name": next(iter(found)),
                    "label": site.label,
                    "atoms": found,
                    "wyckoff": site.wyckoff,
                    "kind": site.kind,
                }
            )
        by_label = {site["label"]: site["name"] for site in sites}
        groups: dict[str, list[str]] = {}
        for site in entry.sites:
            groups.setdefault(site.uiso_group or site.label, []).append(
                by_label[site.label]
            )
        label = spec.origin or entry.origin_site
        axis = spec.origin_axis or (
            AXIS_COORDINATE[entry.polar_axis] if entry.polar_axis else None
        )
        library = entry.name
    else:
        if spec.atoms is None:
            raise PipelineError(
                f"structures.{spec.key}: the coordinates and occupancies modes need "
                "its sites; give atoms, the CIF's sites with their kinds"
            )
        sites = []
        for site in spec.atoms:
            found = site_atoms(site["name"], site["atoms"])
            sites.append({**site, "name": next(iter(found)), "atoms": found})
        by_label = {site["name"]: site["name"] for site in sites}
        groups = {}
        for site in sites:
            groups.setdefault(site["kind"], []).append(site["name"])
        label, axis = spec.origin, spec.origin_axis
        library = None
    origin = None
    if spec.origin_fixed and label is not None:
        if label not in by_label:
            raise PipelineError(
                f"structures.{spec.key}: the origin site {label} is not one of its "
                f"sites {', '.join(by_label)}"
            )
        if axis is None:
            raise PipelineError(
                f"structures.{spec.key}: the origin site {label} needs an axis; "
                "give origin = {site, axis}"
            )
        origin = {"site": by_label[label], "axis": axis}
    return {
        "name": spec.key,
        "library": library,
        "sites": sites,
        "free_coordinates": {},
        "uiso_groups": [
            {"name": name, "sites": members} for name, members in groups.items()
        ],
        "origin": origin,
        "exchange": None,
    }


def _plan(phase: PhaseInput, atoms: Sequence[Mapping]) -> dict:
    """The site plan of a phase, with its name and its exchange groups, each
    over the sites of one kind that hold any of its elements."""
    table = _structure_table(phase, atoms)
    try:
        plan = site_setup(table, atoms)
    except ValueError as error:
        raise PipelineError(f"structures.{phase.key}: {error}") from None
    plan["phase"] = phase.key
    groups = []
    for elements in phase.spec.exchange:
        sites = [
            site
            for site in plan["sites"]
            if {_element_of(atom["type"]) for atom in site["atoms"]} & set(elements)
        ]
        kinds = {site["kind"] for site in sites}
        if len(sites) < 2 or len(kinds) != 1:
            raise PipelineError(
                f"structures.{phase.key}: the exchange of {', '.join(elements)} needs "
                "two or more sites of one kind holding them, not "
                + (", ".join(f"{s['name']} ({s['kind']})" for s in sites) or "none")
            )
        groups.append({"elements": list(elements), "sites": [s["name"] for s in sites]})
    plan["exchange"] = groups
    return plan


def exchange_edits(atoms: Sequence[Mapping], plan: Mapping) -> list[dict]:
    """Atoms at occupancy 0 for each exchanged element a site of an exchange
    group lacks, so that the element can move there, each labelled by the
    element and the digits of the site's name."""
    labels = {str(atom["label"]) for atom in atoms}
    by_name = {site["name"]: site for site in plan["sites"]}
    edits = []
    for group in _exchange_groups(plan):
        for name in group["sites"]:
            site = by_name[name]
            present = {_element_of(atom["type"]) for atom in site["atoms"]}
            for element in group["elements"]:
                if element in present:
                    continue
                label = element + re.sub(r"^[A-Za-z]+", "", name)
                while label in labels:
                    label += "x"
                labels.add(label)
                present.add(element)
                edits.append(
                    {"label": label, "type": element, "copy": name, "occupancy": 0.0}
                )
    return edits


def _joined_stages(builder, plans, background, phases, displacement, mustrain):
    """One phase's stages as ``builder`` gives them, or, with several, the
    first phase's profile and every phase's own stages after it, named by
    the phase."""
    lists = [
        builder(
            plan,
            background,
            phases=phases,
            displacement=displacement,
            mustrain=mustrain,
        )
        for plan in plans
    ]
    if len(lists) == 1:
        return lists[0]
    first = dict(lists[0][0])
    if "uiso_groups" in first:
        first["uiso_groups"] = {
            key: groups
            for stage_list in lists
            for key, groups in stage_list[0]["uiso_groups"].items()
        }
    stages = [first]
    for plan, stage_list in zip(plans, lists, strict=True):
        stages += [
            {**stage, "name": f"{plan['phase']}: {stage['name']}"}
            for stage in stage_list[1:]
        ]
    return stages


# The instrument parameters a run never refines, whatever its stages.
HELD_INSTRUMENT = ("U", "V", "W", "X", "Y", "Z", "SH/L")


def _held_changes(start: Path, used: Path, zero_held: bool) -> list[str]:
    before, after = _instprm_values(start), _instprm_values(used)
    changed = []
    for key in (*HELD_INSTRUMENT, *(("Zero",) if zero_held else ())):
        if key not in before or key not in after:
            continue
        try:
            same = math.isclose(
                float(before[key]), float(after[key]), rel_tol=1e-9, abs_tol=1e-12
            )
        except ValueError:
            same = before[key] == after[key]
        if not same:
            changed.append(f"{key} {before[key]} -> {after[key]}")
    return changed


def _xrdkit_version() -> str:
    try:
        return version("xrdkit")
    except PackageNotFoundError:
        return "0.0.0+unknown"


def _plain(value):
    """``value`` with paths as text, for the result JSON."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def _inputs_record(
    inputs: Inputs, options: Options, cells: Mapping, sources: Mapping
) -> dict:
    instrument = inputs.instrument
    return _plain(
        {
            "project": {"name": inputs.project.name, "root": inputs.project.root},
            "sample": {
                "key": inputs.sample.key,
                "file": inputs.sample.file,
                "form": inputs.sample.form,
                "structures": list(inputs.sample.structures),
            },
            "instrument": {
                "key": instrument.key,
                "wavelength": list(instrument.wavelength),
                "ka2": instrument.ka2,
                "radius": instrument.radius,
                "instprm": instrument.instprm,
            },
            "structures": [
                {
                    "key": phase.key,
                    "library": phase.spec.library,
                    "cif": phase.spec.cif,
                    "composition": phase.spec.composition,
                    "z": phase.z,
                    "cell": phase.spec.cell,
                    "exchange": [list(group) for group in phase.spec.exchange],
                    "origin": phase.spec.origin,
                    "origin_axis": phase.spec.origin_axis,
                    "origin_fixed": phase.spec.origin_fixed,
                    "atoms": list(phase.spec.atoms) if phase.spec.atoms else None,
                }
                for phase in inputs.phases
            ],
            "refine": asdict(inputs.refine),
            "scan_range": list(inputs.scan_range),
            "limits": list(inputs.limits),
            "start_cell": cells,
            "start_cell_source": dict(sources),
            "displacement": inputs.displacement,
            "zero": inputs.zero,
            "options": asdict(options),
        }
    )


def _reduced_chi_squared(stage: Mapping) -> float | None:
    """A stage's reduced chi squared, or GOF squared for a record without it."""
    if stage.get("reduced_chi_squared") is not None:
        return stage["reduced_chi_squared"]
    gof = stage.get("gof")
    return None if gof is None else gof**2


def _residuals(result: Mapping, phases: int) -> dict:
    kept = accepted_stages(result)
    last = kept[-1] if kept else {}
    residuals = {
        "rwp": last.get("rwp"),
        "rp": last.get("rp"),
        "chi_squared": last.get("chi_squared"),
        "reduced_chi_squared": _reduced_chi_squared(last),
    }
    if phases > 1:
        residuals["weight_fractions"] = {
            phase["name"]: phase.get("weight_fraction")
            for phase in (result.get("final") or {}).get("phases", [])
        }
    return residuals


def run_mode(
    project: Project,
    sample: Sample | str,
    mode: str,
    options: Options | None = None,
    reporter=None,
) -> Outcome:
    """Run one mode of a sample's refinement, every input from the project
    file, and write its result JSON, GSAS-II project and exports.

    ``lebail`` extracts every phase from the start cell; ``fixed_atoms``
    starts from the Le Bail result's size stage (its microstrain stage with
    ``options.mustrain``), with each phase's CIF atoms at the nominal
    composition; ``coordinates`` and ``occupancies`` start from the result
    of the mode before, its atoms carried over, and ``occupancies`` adds
    every exchanged element a site lacks at occupancy 0. The job's
    ``start_model``, the structure as the Rietveld modes first set it up,
    every atom carrying the Uiso it is started at, is what undetermined
    parameters are judged against. ``reporter``, a
    callable taking a line of text, hears what the run does.

    Raises
    ------
    PipelineError
        If an input is missing or does not fit, the mode before left no
        result, a stage failed, no stage was accepted, the run left no final
        model, or an instrument parameter the run holds came out changed.
        The error carries the run's ``result`` and ``log`` when it got that
        far, else None.
    xrdkit.gsas2.Gsas2Error
        If GSAS-II itself failed, carrying ``result`` None and the ``log``.
    """
    context: dict = {"result": None, "log": None}
    try:
        return _run_mode(
            project, sample, mode, options or Options(), reporter or _noop, context
        )
    except Exception as error:
        error.result = context["result"]
        error.log = context["log"]
        raise


def _run_mode(project, sample, mode, options, report, context) -> Outcome:
    paths = mode_paths(project, sample, mode, options)
    inputs = resolve_inputs(project, sample, options)
    sample = inputs.sample
    work = paths["work"]
    work.mkdir(parents=True, exist_ok=True)
    count = len(inputs.phases)
    max_passes = options.max_passes or inputs.refine.max_passes[mode]
    on_unsettled = inputs.refine.unsettled[mode]
    common = {
        "data_file": sample.file,
        "cycles": CYCLES,
        "export_prefix": paths["prefix"],
        "max_passes": max_passes,
        "pass_tolerance": PASS_TOLERANCE,
        "on_unsettled": on_unsettled,
    }

    plans: list[dict] = []
    sources: dict[str, str] = {}
    if mode == "lebail":
        instprm = _write_start_instprm(
            inputs.instprm,
            inputs.zero,
            paths["start_instprm"],
            inputs.instrument.radius,
        )
        phases = []
        cells = {}
        for phase in inputs.phases:
            if phase.spec.cif is not None:
                cif = phase.spec.cif
            elif phase.entry is not None and phase.cell is not None:
                cif = _entry_cif(phase, work)
            else:
                raise PipelineError(
                    f"structures.{phase.key}: gives no CIF and no cell to write one"
                )
            entry = {"cif": cif, "name": phase.key}
            if phase.cell is not None:
                entry["cell"] = _six(phase.cell)
            cells[phase.key] = phase.cell
            sources[phase.key] = phase.cell_source
            phases.append(entry)
        stages = lebail_stages(
            inputs.refine.background,
            phases=count,
            displacement=inputs.displacement,
            mustrain_test=True,
        )
        job = build_refine_job(
            paths["gpx"],
            stages,
            instprm=instprm,
            phases=phases,
            limits=inputs.limits,
            broadening={"*": START_BROADENING},
            le_bail_cycles=LE_BAIL_CYCLES,
            **common,
        )
    else:
        _need_cifs(inputs, mode)
        before = MODES[MODES.index(mode) - 1]
        before_result = mode_paths(project, sample, before, options)["result"]
        if not before_result.is_file():
            raise PipelineError(
                f"the {mode} mode starts from the {before} result {before_result}, "
                f"which is missing; run {before} first"
            )
        if mode == "fixed_atoms":
            start = start_from_result(
                before_result, MICROSTRAIN if options.mustrain else SIZE
            )
        else:
            start = start_from_result(before_result)
        by_name = {phase.name: phase for phase in start.phases}
        missing = [phase.key for phase in inputs.phases if phase.key not in by_name]
        if missing:
            raise PipelineError(
                f"{before_result}: has no phase {', '.join(missing)} of the sample"
            )
        instprm = _write_start_instprm(
            inputs.instprm,
            start.zero,
            paths["start_instprm"],
            inputs.instrument.radius,
        )
        cells = {key: dict(phase.cell) for key, phase in by_name.items()}
        sources = {
            key: f"the {before} result {before_result.name}, stage {start.stage}"
            for key in cells
        }
        base = _create(inputs, instprm, work / "cif")
        if mode == "fixed_atoms":
            edits = {
                phase.key: _nominal_edits(phase, base[phase.key]["atoms"])
                for phase in inputs.phases
            }
            set_up = _create(inputs, instprm, work / "set_up", edits)
            uiso = {
                key: _mean_ueq(set_up[key]["atoms"], set_up[key]["cell"])
                for key in edits
            }
            start_model = {
                key: _with_start_uiso(set_up[key]["atoms"], uiso[key]) for key in edits
            }
            stages = fixed_atoms_stages(
                {
                    "function": start.background["function"],
                    "terms": start.background["terms"],
                },
                phases=count,
                displacement=inputs.displacement,
                mustrain=options.mustrain,
                preferred_orientation=options.preferred_orientation,
            )
            scale_start = None
        else:
            edits = {}
            for phase in inputs.phases:
                atoms = by_name[phase.key].atoms
                plan = _plan(phase, atoms)
                plans.append(plan)
                try:
                    edits[phase.key] = structure_edits(atoms, base[phase.key]["atoms"])
                except ValueError as error:
                    raise PipelineError(f"structures.{phase.key}: {error}") from None
                if mode == "occupancies":
                    edits[phase.key] += exchange_edits(atoms, plan)
            builder = coordinates_stages if mode == "coordinates" else occupancy_stages
            stages = _joined_stages(
                builder,
                plans,
                {
                    "function": start.background["function"],
                    "terms": start.background["terms"],
                },
                count,
                inputs.displacement,
                options.mustrain,
            )
            earlier = (
                start.start_model if isinstance(start.start_model, Mapping) else {}
            )
            start_model = {key: earlier.get(key) or by_name[key].atoms for key in edits}
            uiso = {
                key: _mean_ueq(by_name[key].atoms, _gsas_cell(by_name[key].cell))
                for key in edits
            }
            scale_start = start.scale
        anions = sorted(
            {
                anion
                for phase in inputs.phases
                for anion in (phase.entry.anions if phase.entry else DEFAULT_ANIONS)
            }
        )
        job = build_refine_job(
            paths["gpx"],
            stages,
            instprm=instprm,
            phases=[
                {
                    "cif": phase.spec.cif,
                    "name": phase.key,
                    "cell": _six(by_name[phase.key].cell),
                    "atoms": edits[phase.key],
                }
                for phase in inputs.phases
            ],
            limits=start.limits,
            broadening={
                key: {
                    "size": phase.size,
                    "mustrain": phase.microstrain if options.mustrain else 0.0,
                    "lgmix": START_BROADENING["lgmix"],
                }
                for key, phase in by_name.items()
                if key in edits
            },
            background_start={
                "type": start.background["function"],
                "coefficients": start.background["coefficients"],
            },
            overall_uiso_start=uiso,
            scale_start=scale_start,
            displacement_start=(
                start.displacement
                if inputs.displacement and start.displacement is not None
                else None
            ),
            phase_fraction_start=(
                {key: by_name[key].fraction for key in edits} if count > 1 else None
            ),
            start_model=start_model,
            bonds={"anions": anions},
            **common,
        )

    job["result"] = str(paths["result"])
    context["log"] = paths["log"]
    report(
        f"{sample.key}: {mode} started, {len(stages)} stages, at most {max_passes} "
        "passes each"
    )
    result = run_job(job, work)
    for row in stage_statuses(result):
        report(
            f"{sample.key}: {mode}: {row['name']} {row['status']}"
            + (f" ({row['reason']})" if row.get("reason") else "")
        )
    result["inputs"] = _inputs_record(inputs, options, cells, sources)
    result["method"] = {
        "mode": mode,
        "stages": _plain(stages),
        "cycles": CYCLES,
        "le_bail_cycles": LE_BAIL_CYCLES if mode == "lebail" else None,
        "max_passes": max_passes,
        "pass_tolerance": PASS_TOLERANCE,
        "driver_version": gsas2_driver.DRIVER_VERSION,
        "xrdkit_version": _xrdkit_version(),
    }
    result["date"] = datetime.now().astimezone().isoformat(timespec="seconds")
    paths["result"].parent.mkdir(parents=True, exist_ok=True)
    paths["result"].write_text(json.dumps(_plain(result), indent=2), encoding="utf-8")
    context["result"] = result

    if not result.get("completed"):
        failed = (result.get("stages") or [{}])[-1]
        raise PipelineError(
            f"{mode}: stage {failed.get('name')!r} failed: {failed.get('error')}"
        )
    kept = accepted_stages(result)
    if not kept:
        raise PipelineError(
            f"{mode}: no stage was accepted; "
            + "; ".join(
                f"{row['name']} {row['status']}" for row in stage_statuses(result)
            )
        )
    if not result.get("final"):
        raise PipelineError(
            f"{mode}: the run left no final model: "
            + str(result.get("final_error") or result.get("compute_error") or "")
        )
    exports = result.get("exports") or {}
    if not exports.get("instprm"):
        raise PipelineError(f"{mode}: the run exported no instrument parameter file")
    changed = _held_changes(instprm, Path(exports["instprm"]), inputs.displacement)
    if changed:
        raise PipelineError(
            f"{mode}: instrument parameters held by the run came out changed: "
            + ", ".join(changed)
        )
    # The figure and the write up, after the result is safely written.
    figures = _figures(result, paths, sample.key, mode)
    earlier = {}
    for before in MODES[: MODES.index(mode)]:
        path = mode_paths(project, sample, before, options)["result"]
        try:
            earlier[before] = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
    markdown = paths["folder"] / f"{mode}.md"
    markdown.write_text(
        writeup.MARKDOWN[mode](
            result,
            result["inputs"],
            {plan["phase"]: plan for plan in plans},
            earlier,
        ),
        encoding="utf-8",
    )
    return Outcome(
        mode=mode,
        accepted=[stage["name"] for stage in kept],
        final=result["final"],
        residuals=_residuals(result, count),
        undetermined=list(result.get("undetermined") or []),
        paths=_plain(
            {
                "folder": paths["folder"],
                "result": paths["result"],
                "gpx": paths["gpx"],
                "work": work,
                "log": paths["log"],
                "start_instprm": instprm,
                "figures": figures,
                "markdown": markdown,
                **{key: value for key, value in exports.items()},
            }
        ),
    )


def _figures(result: Mapping, paths: Mapping, key: str, mode: str) -> list[Path]:
    """The fit drawn with plot_rietveld and saved as PNG and PDF beside the
    result."""
    exports = result.get("exports") or {}
    if not exports.get("histogram"):
        raise PipelineError(f"{mode}: the run exported no pattern to draw")
    fig, _ = plot_rietveld(
        exports["histogram"],
        exports.get("reflections") or None,
        title=f"{key}: {mode}",
        sqrt_scale=True,
        result=result,
    )
    return save_figure(fig, paths["prefix"])


def _gsas_cell(cell: Mapping[str, float]) -> dict[str, float]:
    return {gsas: float(cell[key]) for key, gsas in _CELL_KEYS}


def _summary_values(result: Mapping | None) -> dict[str, str]:
    final = (result or {}).get("final") or {}
    phases = final.get("phases") or []
    if not phases:
        return {}
    cell = phases[0].get("cell") or {}
    zero = (final.get("instrument") or {}).get("Zero") or {}
    values = {
        column: f"{cell[key]:.5f}"
        for column, key in (
            ("a (Å)", "length_a"),
            ("b (Å)", "length_b"),
            ("c (Å)", "length_c"),
        )
        if key in cell
    }
    if "value" in zero:
        values["zero (°)"] = f"{zero['value']:.5f}"
    return values


SUMMARY_COLUMNS = ("a (Å)", "b (Å)", "c (Å)", "zero (°)")


def run_sequence(
    project: Project,
    sample: Sample | str,
    modes: Sequence[str],
    options: Options | None = None,
    reporter=None,
) -> list[Outcome]:
    """Run ``modes`` of a sample's refinement in order, each from the one
    before, and stop at the first that fails.

    A mode that raises, whatever the error, is written up as ``failure.md``
    in its folder from the error, the stages its own run got through (none
    when it wrote no result, never an older result) and the tail of its
    GSAS-II log; the modes after it are not run. ``summary.md`` is written
    in the folder of the last mode run, over the modes run, whatever became
    of them. Returns the outcome of every mode run, the failed one last
    with its ``error``.

    Raises
    ------
    PipelineError
        If a mode is not one of :data:`MODES`, before any is run.
    """
    options = options or Options()
    report = reporter or _noop
    sample = _sample_of(project, sample)
    unknown = [mode for mode in modes if mode not in MODES]
    if unknown:
        raise PipelineError(
            f"no mode {', '.join(map(repr, unknown))}; the modes are {', '.join(MODES)}"
        )
    for mode in modes:
        (mode_paths(project, sample, mode, options)["folder"] / "failure.md").unlink(
            missing_ok=True
        )
    outcomes, entries = [], []
    last_folder = None
    for mode in modes:
        paths = mode_paths(project, sample, mode, options)
        last_folder = paths["folder"]
        run = datetime.now().astimezone().isoformat(timespec="minutes")
        try:
            outcome = run_mode(project, sample, mode, options, report)
        except Exception as error:  # noqa: BLE001
            result = getattr(error, "result", None)
            log = getattr(error, "log", None)
            message = f"{type(error).__name__}: {error}"
            later = list(modes[modes.index(mode) + 1 :])
            intro = [
                f"The {mode} mode of {sample.key} did not finish, so this is what "
                "its run got through rather than a reading of a fit."
                + (f" {', '.join(later)} not run." if later else "")
            ]
            intro.append(
                "The run wrote no result."
                if result is None
                else f"The result was written to `{paths['result']}`, with its stages."
            )
            failure = paths["folder"] / "failure.md"
            failure.parent.mkdir(parents=True, exist_ok=True)
            failure.write_text(
                failure_markdown(
                    f"{sample.key}: {mode}",
                    result,
                    message,
                    log_tail(log) if log else "",
                    intro=intro,
                ),
                encoding="utf-8",
            )
            outcomes.append(
                Outcome(
                    mode=mode,
                    accepted=[s["name"] for s in accepted_stages(result)]
                    if result
                    else [],
                    final=(result or {}).get("final"),
                    residuals=_residuals(result, len(sample.structures))
                    if result
                    else {},
                    undetermined=list((result or {}).get("undetermined") or []),
                    paths=_plain(
                        {
                            "folder": paths["folder"],
                            "failure": failure,
                            **({"result": paths["result"]} if result else {}),
                            **({"log": log} if log else {}),
                        }
                    ),
                    error=message,
                )
            )
            entries.append(
                {
                    "name": mode,
                    "run": run,
                    "result": result,
                    "values": _summary_values(result),
                }
            )
            break
        outcomes.append(outcome)
        result = json.loads(Path(outcome.paths["result"]).read_text(encoding="utf-8"))
        entries.append(
            {
                "name": mode,
                "run": run,
                "result": result,
                "values": _summary_values(result),
            }
        )
    if last_folder is not None:
        summary = last_folder / "summary.md"
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(
            summary_markdown(
                f"{sample.key}: refinement summary",
                entries,
                columns=SUMMARY_COLUMNS,
                intro=[
                    (
                        f"The modes run for {sample.key} in this sequence, "
                        f"{', '.join(entry['name'] for entry in entries)}, as their "
                        "results leave them."
                    )
                ],
            ),
            encoding="utf-8",
        )
        for outcome in outcomes:
            outcome.paths["summary"] = str(summary)
    return outcomes
