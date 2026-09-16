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

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from xrdkit.gsas2 import accepted_stages

__all__ = [
    "CYCLES",
    "LE_BAIL_CYCLES",
    "MODES",
    "PASS_TOLERANCE",
    "START_BROADENING",
    "PhaseStart",
    "PipelineError",
    "StartPoint",
    "coordinates_stages",
    "fixed_atoms_stages",
    "lebail_stages",
    "occupancy_stages",
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
