"""Markdown write ups of the pipeline's modes.

Each generator takes a mode's result JSON as :func:`xrdkit.pipeline.run_mode`
writes it, with the ``inputs`` recorded in it, and optionally the site plan
of each phase (:func:`xrdkit.structure.site_setup`, with ``phase`` set) and
the results of the modes before it, and returns markdown:

:func:`lebail_markdown`
    Settings, the phases extracted, stage outcomes, the cell against the
    start, the zero or displacement, size and the microstrain test,
    background and the reflections followed.
:func:`fixed_atoms_markdown`
    The same for a refinement with the atoms fixed, with the reflection
    misfits, the weight fractions of two phases or more and the sanity check.
:func:`structure_markdown`
    The coordinates and occupancies modes: besides the above, the coordinate
    shifts from the start model in angstroms, the parameters left
    undetermined, the occupancies with each exchange group's held total and
    the bond lengths against the bond limits.

Every stage is looked up by its name, so that a rejected stage has a row of
its own saying so and nothing else depends on where it stood. Every name of
a site kind comes from the plan.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Mapping, Sequence
from itertools import pairwise
from pathlib import Path

import numpy as np

from xrdkit.gsas2 import accepted_stages, stage_status, stage_statuses

__all__ = [
    "fixed_atoms_markdown",
    "lebail_markdown",
    "relative",
    "structure_markdown",
    "with_esd",
]

# A reflection within this many FWHM of another overlaps it, and the worst
# misfits listed.
OVERLAP_WIDTHS = 0.5
WORST_REFLECTIONS = 10
# A microstrain this many esds from zero that also lowers Rwp is taken as
# found by the test.
MICROSTRAIN_SIGNIFICANCE = 3.0
# The stages whose names the Le Bail write up reads.
SIZE_STAGE = "size"
MICROSTRAIN_STAGE = "microstrain"
CELL_NAMES = (
    ("a", "length_a", "a (Å)"),
    ("b", "length_b", "b (Å)"),
    ("c", "length_c", "c (Å)"),
    ("alpha", "angle_alpha", "α (°)"),
    ("beta", "angle_beta", "β (°)"),
    ("gamma", "angle_gamma", "γ (°)"),
)
AXIS_INDEX = {"x": ("h", 0), "y": ("k", 1), "z": ("l", 2)}


# Formatting


def with_esd(value: float, esd: float | None, digits: int = 6) -> str:
    """``value`` with its esd in brackets, the esd to one figure, or two if
    it begins with a 1; a value with no esd, or an esd of 0, is marked
    fixed and given to ``digits`` significant figures."""
    if esd is None or not esd > 0 or not math.isfinite(esd):
        return f"{value:.{digits}g} (fixed)"
    decimals = max(-int(np.floor(np.log10(esd))), 0)
    if round(esd * 10**decimals) == 1:
        decimals += 1
    return f"{value:.{decimals}f}({round(esd * 10**decimals)})"


def relative(value, root: str | Path):
    """``value`` with every absolute path under ``root``, in text or nested
    in lists and mappings, made relative to it, with forward slashes."""
    root = Path(root).resolve()
    if isinstance(value, Mapping):
        return {key: relative(item, root) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [relative(item, root) for item in value]
    if isinstance(value, (str, Path)):
        try:
            path = Path(value)
            if path.is_absolute() and path.resolve().is_relative_to(root):
                return path.resolve().relative_to(root).as_posix()
        except (OSError, ValueError):
            pass
    return value


def _figure(value, digits: int = 3) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


def _cell(text: object) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def _row(cells: Sequence[object]) -> str:
    return "| " + " | ".join(_cell(cell) for cell in cells) + " |"


def _table(header: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    return [_row(header), "|" + " --- |" * len(header), *(_row(row) for row in rows)]


def _joined(items: Sequence[str]) -> str:
    items = list(items)
    if not items:
        return ""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _entry(value) -> tuple[float | None, float | None]:
    """A recorded value and its esd."""
    if isinstance(value, Mapping):
        return value.get("value"), value.get("esd")
    return value, None


def _plans(plan) -> dict[str, Mapping]:
    if not plan:
        return {}
    if "sites" in plan:
        return {plan.get("phase"): plan}
    return dict(plan)


# Sections shared by every mode


def _inputs(result: Mapping, inputs: Mapping | None) -> Mapping:
    return inputs if inputs is not None else (result.get("inputs") or {})


def _stage_index(result: Mapping) -> dict[str, Mapping]:
    return {stage.get("name"): stage for stage in result.get("stages", [])}


def settings_lines(result: Mapping, inputs: Mapping) -> list[str]:
    """The sample, scan, instrument, range, phases with their start cells
    and where they came from, background, and the pass settings."""
    sample = inputs.get("sample") or {}
    instrument = inputs.get("instrument") or {}
    method = result.get("method") or {}
    limits = result.get("limits") or inputs.get("limits")
    root = (inputs.get("project") or {}).get("root")

    def shown(path):
        return relative(path, root) if root and path else path

    lines = ["## Settings", ""]
    if sample:
        lines.append(
            f"- Sample {sample.get('key')}, {sample.get('form')}, scan "
            f"`{shown(sample.get('file'))}`."
        )
    if instrument:
        radius = instrument.get("radius")
        lines.append(
            f"- Instrument {instrument.get('key')}: wavelengths "
            f"{', '.join(f'{w:g}' for w in instrument.get('wavelength', []))} Å"
            + (f", radius {radius:g} mm" if radius else "")
            + f", parameters `{shown(instrument.get('instprm'))}`, held but for the "
            + ("displacement." if inputs.get("displacement") else "zero.")
        )
    if limits:
        lines.append(f"- Range {limits[0]:.2f} to {limits[1]:.2f}° 2θ.")
    background = background_setting(result, inputs)
    if background:
        lines.append(
            f"- Background: {background['function']}, {background['terms']} terms."
        )
    sources = inputs.get("start_cell_source") or {}
    cells = inputs.get("start_cell") or {}
    for structure in inputs.get("structures") or []:
        key = structure.get("key")
        made = (
            f"library entry {structure['library']}"
            if structure.get("library")
            else "its CIF"
        )
        if structure.get("library") and structure.get("cif"):
            made += f" with the coordinates of `{shown(structure['cif'])}`"
        cell = cells.get(key)
        start = (
            ", ".join(
                f"{short} = {cell[short]:.5g}"
                for short, _, _ in CELL_NAMES
                if short in cell
            )
            if cell
            else "the CIF's"
        )
        lines.append(
            f"- Phase {key}, {structure.get('composition')}, from {made}; the cell "
            f"starts from {start} ({sources.get(key, 'as recorded')})."
        )
    if method:
        lines.append(
            f"- At most {method.get('cycles')} least squares cycles a refinement "
            f"and {method.get('max_passes')} passes a stage, until no parameter "
            f"moves by more than {method.get('pass_tolerance')} esd; driver "
            f"version {method.get('driver_version')}, xrdkit "
            f"{method.get('xrdkit_version')}."
        )
    if result.get("date"):
        lines.append(f"- Run {result['date']}.")
    return lines


def background_setting(result: Mapping, inputs: Mapping) -> dict | None:
    """The background function and terms as the job's stages set them, or
    as the refine settings give them."""
    for stage in (result.get("method") or {}).get("stages") or []:
        background = stage.get("background")
        if isinstance(background, Mapping):
            return {
                "function": background.get("type") or background.get("function"),
                "terms": background.get("terms"),
            }
    refine = (inputs.get("refine") or {}).get("background")
    return dict(refine) if refine else None


def stage_outcome_lines(result: Mapping) -> list[str]:
    """Every stage with its status, passes, Rwp, Rp and chi squared, looked
    up by name; a stage rejected or failed says so, and why."""
    reasons = {row["name"]: row for row in stage_statuses(result)}
    rows = []
    for stage in result.get("stages", []):
        name = stage.get("name")
        status = stage_status(stage)
        reason = reasons.get(name, {}).get("reason") or ""
        if status in ("rejected", "failed"):
            rows.append([name, status, "—", "—", "—", "—", reason])
            continue
        passes = stage.get("passes")
        count = "—" if passes is None else str(len(passes))
        rows.append(
            [
                name,
                status,
                count,
                _figure(stage.get("rwp")),
                _figure(stage.get("rp")),
                _figure(stage.get("chi_squared")),
                reason,
            ]
        )
    lines = [
        "## Stage outcomes",
        "",
        *_table(["Stage", "Status", "Passes", "Rwp (%)", "Rp (%)", "χ²", "Why"], rows),
        "",
        final_model_line(result),
    ]
    return lines


def final_model_line(result: Mapping) -> str:
    kept = accepted_stages(result)
    rejected = result.get("rejected") or []
    if not kept:
        return (
            "No stage was kept: the model below is the one the refinement started "
            "from, computed but not refined."
        )
    line = f"The final model is that of stage {kept[-1]['name']}"
    if rejected:
        line += f"; {_joined(rejected)} rejected and rolled back"
    return line + "."


def cell_lines(result: Mapping, inputs: Mapping) -> list[str]:
    """Each phase's cell with esds against the cell it started from."""
    final = result.get("final") or {}
    starts = inputs.get("start_cell") or {}
    lines = ["## Cell", ""]
    for phase in final.get("phases") or []:
        name = phase.get("name")
        cell = phase.get("cell") or {}
        esds = phase.get("cell_esd") or {}
        start = starts.get(name) or {}
        rows = []
        for short, gsas, label in CELL_NAMES:
            if gsas not in cell:
                continue
            value = cell[gsas]
            esd = esds.get(gsas)
            before = start.get(short)
            rows.append(
                [
                    label,
                    with_esd(value, esd),
                    "—" if before is None else f"{before:.6g}",
                    "—" if before is None else f"{value - before:+.5f}",
                ]
            )
        volume = cell.get("volume")
        if volume is not None:
            rows.append(["V (Å³)", with_esd(volume, esds.get("volume")), "—", "—"])
        lines += [
            f"Phase {name}:",
            "",
            *_table(["", "Refined", "Start", "Change"], rows),
            "",
        ]
    return lines


def zero_lines(result: Mapping, inputs: Mapping) -> list[str]:
    """The zero shift, or the specimen displacement, whichever was refined."""
    final = result.get("final") or {}
    if inputs.get("displacement"):
        value, esd = _entry((final.get("sample") or {}).get("Shift"))
        name, unit = "Specimen displacement", "(GSAS-II's Shift, in µm)"
        start = None
    else:
        value, esd = _entry((final.get("instrument") or {}).get("Zero"))
        name, unit = "Zero shift", "(°)"
        start = inputs.get("zero")
    if value is None:
        return []
    line = f"{name} {unit}: {with_esd(value, esd)}"
    if start is not None:
        line += f", from {start:.5g}"
    return ["## Zero or displacement", "", line + "."]


def broadening_lines(result: Mapping, lebail: bool) -> list[str]:
    """Each phase's size and microstrain, and in a Le Bail result the verdict
    of the microstrain test, from the size and microstrain stages by name."""
    final = result.get("final") or {}
    rows = []
    for phase in final.get("phases") or []:
        size, size_esd = _entry(phase.get("size"))
        strain, strain_esd = _entry(phase.get("mustrain"))
        rows.append(
            [
                phase.get("name"),
                "—" if size is None else with_esd(size, size_esd),
                "—" if strain is None else with_esd(strain, strain_esd),
            ]
        )
    lines = [
        "## Size and microstrain",
        "",
        "Isotropic, as GSAS-II refines them: size in µm, microstrain in 10⁻⁶.",
        "",
        *_table(["Phase", "Size (µm)", "Microstrain"], rows),
    ]
    if lebail:
        lines += ["", microstrain_verdict(result)]
    return lines


def microstrain_verdict(result: Mapping) -> str:
    """Whether the microstrain test found a microstrain: its stage kept, the
    microstrain of the first phase MICROSTRAIN_SIGNIFICANCE esds or more from
    zero, and Rwp lowered from the size stage."""
    stages = _stage_index(result)
    strain_stage = stages.get(MICROSTRAIN_STAGE)
    size_stage = stages.get(SIZE_STAGE)
    if strain_stage is None:
        return "The microstrain test was not run."
    status = stage_status(strain_stage)
    if status in ("rejected", "failed"):
        return f"The microstrain test did not finish: its stage was {status}."
    phases = (strain_stage.get("values") or {}).get("phases") or (
        (result.get("final") or {}).get("phases") or []
    )
    if not phases:
        return "The microstrain test left no value to judge."
    value, esd = _entry(phases[0].get("mustrain"))
    if value is None or not esd:
        return "The microstrain test left the microstrain without an esd."
    sigmas = abs(value) / esd
    lowered = (
        size_stage is not None
        and stage_status(size_stage) not in ("rejected", "failed")
        and strain_stage.get("rwp") is not None
        and size_stage.get("rwp") is not None
        and strain_stage["rwp"] < size_stage["rwp"]
    )
    rwp = (
        f", Rwp {size_stage['rwp']:.3f} to {strain_stage['rwp']:.3f}%"
        if size_stage is not None
        and size_stage.get("rwp") is not None
        and strain_stage.get("rwp") is not None
        else ""
    )
    if sigmas >= MICROSTRAIN_SIGNIFICANCE and lowered:
        return (
            f"The microstrain test finds a microstrain: {sigmas:.1f} esds from "
            f"zero{rwp}."
        )
    return (
        f"The microstrain test finds none that the data support: {sigmas:.1f} "
        f"esds from zero{rwp}; the size refined alone is the better estimate."
    )


def weight_fraction_lines(result: Mapping) -> list[str]:
    phases = (result.get("final") or {}).get("phases") or []
    if len(phases) < 2:
        return []
    rows = []
    for phase in phases:
        fraction, fraction_esd = _entry(phase.get("phase_fraction"))
        weight, weight_esd = _entry(phase.get("weight_fraction"))
        rows.append(
            [
                phase.get("name"),
                "—" if fraction is None else with_esd(fraction, fraction_esd),
                "—" if weight is None else with_esd(weight, weight_esd),
            ]
        )
    return [
        "## Phase fractions",
        "",
        *_table(["Phase", "Phase fraction", "Weight fraction"], rows),
    ]


def sanity_lines(result: Mapping) -> list[str]:
    settings = result.get("sanity") or {}
    rows = []
    for stage in result.get("stages", []):
        flags = "; ".join(flag.get("message", "") for flag in stage.get("sanity", []))
        rows.append([stage.get("name"), flags or "none"])
    where = (
        "the reference given"
        if settings.get("reference") == "given"
        else "where the job started"
    )
    return [
        "## Sanity check",
        "",
        (
            "After every stage the atoms are checked for a negative Uiso, an "
            "occupancy outside 0 to 1 or a site whose occupancies add up to more "
            f"than 1, and a coordinate more than {settings.get('max_shift', 0.05):g} "
            f"(fractional) from {where}."
        ),
        "",
        *_table(["Stage", "Flags"], rows),
    ]


def undetermined_lines(result: Mapping) -> list[str]:
    """The parameters the data leave undetermined, as find_undetermined
    judges them against the start model."""
    lines = [
        "## Undetermined parameters",
        "",
        (
            "An occupancy whose esd is more than half its range, 0 to 1, and a "
            "coordinate or an isotropic Uiso whose esd is larger than its shift "
            "from the start model: the data do not determine them, and their "
            "values are not a result."
        ),
        "",
    ]
    found = [
        f"- {entry.get('phase', '')} {entry.get('message', '')}".replace("-  ", "- ")
        for entry in result.get("undetermined") or []
    ]
    return lines + (found or ["No refined parameter is left undetermined."])


# Reflections


def _read_reflections(path: str | Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = []
        for row in csv.DictReader(handle):
            try:
                rows.append(
                    {
                        "h": int(row["h"]),
                        "k": int(row["k"]),
                        "l": int(row["l"]),
                        "multiplicity": int(row["multiplicity"]),
                        "two_theta": float(row["two_theta"]),
                        "fwhm": float(row["fwhm"]),
                        "f_obs_squared": float(row["f_obs_squared"]),
                        "f_calc_squared": float(row["f_calc_squared"]),
                        "i_corr": float(row["i_corr"]),
                    }
                )
            except (KeyError, TypeError, ValueError):
                continue
    return rows


def _hkl(row: Mapping) -> str:
    h, k, l = row["h"], row["k"], row["l"]
    return (
        f"{h}{k}{l}"
        if max(abs(h), abs(k), abs(l)) < 10 and min(h, k, l) >= 0
        else f"{h} {k} {l}"
    )


def reflection_misfits(
    path: str | Path,
    limits: Sequence[float] | None = None,
    polar_axis: str | None = None,
    le_bail_path: str | Path | None = None,
) -> dict:
    """Integrated intensities of a phase's reflections, observed (from
    GSAS-II's partition of the pattern) and calculated, the worst misfits,
    and their ratio by class along the polar axis and by angle.

    With the polar axis, x, y or z, reflections are classed by the absolute
    index along it (0, 1, 2, and 3 or more), each counting only if no
    reflection of another class overlaps it; without one there are no such
    classes. The angle bins split the range in four. With a Le Bail
    reflection list each reflection also carries its Le Bail intensity,
    brought to this scale by the ratio of the observed totals and compared
    over the reflection together with those it overlaps.
    """
    table = _read_reflections(path)
    if not table:
        return {
            "rows": [],
            "strongest": 0.0,
            "worst": [],
            "by_class": {},
            "by_angle": {},
        }
    for row in table:
        row["hkl"] = _hkl(row)
        row["observed"] = row["f_obs_squared"] * row["i_corr"]
        row["calculated"] = row["f_calc_squared"] * row["i_corr"]
    le_bail = {}
    if le_bail_path is not None and Path(le_bail_path).is_file():
        le_bail = {
            _hkl(row): row["f_obs_squared"] * row["i_corr"]
            for row in _read_reflections(le_bail_path)
        }
    common = [row for row in table if row["hkl"] in le_bail]
    to_scale = (
        sum(row["observed"] for row in common)
        / sum(le_bail[row["hkl"]] for row in common)
        if common and sum(le_bail[row["hkl"]] for row in common)
        else None
    )
    for row in table:
        row["overlaps"] = [
            other["hkl"]
            for other in table
            if other is not row
            and abs(other["two_theta"] - row["two_theta"])
            <= OVERLAP_WIDTHS * row["fwhm"]
        ]
        row["le_bail"] = (
            le_bail[row["hkl"]] * to_scale
            if to_scale is not None and row["hkl"] in le_bail
            else None
        )
    by_hkl = {row["hkl"]: row for row in table}
    for row in table:
        group = [row, *(by_hkl[other] for other in row["overlaps"] if other in by_hkl)]
        calculated = sum(member["calculated"] for member in group)
        row["le_bail_ratio"] = (
            None
            if any(member["le_bail"] is None for member in group) or not calculated
            else sum(member["le_bail"] for member in group) / calculated
        )
    strongest = max(row["calculated"] for row in table) or 1.0
    worst = sorted(
        table, key=lambda row: abs(row["observed"] - row["calculated"]), reverse=True
    )[:WORST_REFLECTIONS]
    by_class: dict[int, list] = {}
    if polar_axis in AXIS_INDEX:
        index = AXIS_INDEX[polar_axis][0]
        for row in table:
            level = min(abs(row[index]), 3)
            clean = all(
                min(abs(by_hkl[other][index]), 3) == level
                for other in row["overlaps"]
                if other in by_hkl
            )
            if clean:
                by_class.setdefault(level, []).append(row)
    low = limits[0] if limits else min(row["two_theta"] for row in table)
    high = limits[1] if limits else max(row["two_theta"] for row in table)
    edges = list(np.linspace(low, math.ceil(high), 5))
    by_angle: dict[tuple, list] = {}
    for row in table:
        for lower, upper in pairwise(edges):
            if lower <= row["two_theta"] < upper or (
                upper == edges[-1] and row["two_theta"] == upper
            ):
                by_angle.setdefault((lower, upper), []).append(row)
                break
    return {
        "rows": table,
        "strongest": strongest,
        "worst": worst,
        "by_class": dict(sorted(by_class.items())),
        "by_angle": by_angle,
    }


def _ratio(rows: Sequence[Mapping], key: str = "observed") -> str:
    calculated = sum(row["calculated"] for row in rows)
    if not calculated or (key == "le_bail" and any(row[key] is None for row in rows)):
        return "—"
    return f"{sum(row[key] or 0.0 for row in rows) / calculated:.3f}"


def misfit_lines(
    result: Mapping,
    plans: Mapping[str, Mapping],
    lebail_result: Mapping | None = None,
) -> list[str]:
    """For each phase the reflections with the largest intensity misfit, and
    the misfit by class along the polar axis and by angle."""
    exports = (result.get("exports") or {}).get("reflections") or {}
    le_bail_exports = ((lebail_result or {}).get("exports") or {}).get(
        "reflections"
    ) or {}
    lines = ["## Reflection misfits", ""]
    if not exports:
        return lines + ["No reflection list was exported."]
    lines += [
        (
            "Integrated intensity, F² times GSAS-II's intensity factor, observed from "
            "its partition of the pattern against calculated, the largest "
            f"differences first, in per cent of the strongest reflection. Reflections "
            f"within {OVERLAP_WIDTHS:g} FWHM overlap and share their observed "
            "intensity; the Le Bail ratio, where there is a Le Bail fit, owes nothing "
            "to the structure."
        ),
        "",
    ]
    for name, path in exports.items():
        if not Path(path).is_file():
            lines += [f"Phase {name}: its reflection list `{path}` is missing.", ""]
            continue
        polar = (plans.get(name) or {}).get("polar_axis")
        misfits = reflection_misfits(
            path, result.get("limits"), polar, le_bail_exports.get(name)
        )
        rows = [
            [
                row["hkl"],
                f"{row['two_theta']:.3f}",
                row["multiplicity"],
                f"{row['f_obs_squared']:.4g}",
                f"{row['f_calc_squared']:.4g}",
                _figure(
                    row["f_obs_squared"] / row["f_calc_squared"]
                    if row["f_calc_squared"]
                    else None
                ),
                f"{100 * (row['observed'] - row['calculated']) / misfits['strongest']:+.2f}",
                ", ".join(row["overlaps"]) or "—",
                _figure(row["le_bail_ratio"]),
            ]
            for row in misfits["worst"]
        ]
        lines += [
            f"Phase {name}:",
            "",
            *_table(
                [
                    "hkl",
                    "2θ (°)",
                    "m",
                    "F²obs",
                    "F²calc",
                    "F²obs/F²calc",
                    "Difference (%)",
                    "Overlaps",
                    "Le Bail I / I calc",
                ],
                rows,
            ),
            "",
        ]
        if misfits["by_class"]:
            index = AXIS_INDEX[polar][0]
            lines += [
                (
                    f"By the index along the polar axis, {polar} ({index}), counting "
                    "only reflections no reflection of another class overlaps:"
                ),
                "",
                *_table(
                    [
                        f"|{index}|",
                        "Reflections",
                        "ΣI obs / ΣI calc",
                        "ΣI Le Bail / ΣI calc",
                    ],
                    [
                        [
                            "3 or more" if level == 3 else level,
                            len(group),
                            _ratio(group),
                            _ratio(group, "le_bail"),
                        ]
                        for level, group in misfits["by_class"].items()
                    ],
                ),
                "",
            ]
        lines += [
            *_table(
                ["2θ (°)", "Reflections", "ΣI obs / ΣI calc", "ΣI Le Bail / ΣI calc"],
                [
                    [
                        f"{lower:.1f}–{upper:.1f}",
                        len(group),
                        _ratio(group),
                        _ratio(group, "le_bail"),
                    ]
                    for (lower, upper), group in misfits["by_angle"].items()
                ],
            ),
            "",
        ]
    return lines


def followed_lines(
    result: Mapping, inputs: Mapping, earlier: Mapping[str, Mapping] | None = None
) -> list[str]:
    """F²obs / F²calc of the followed reflections in this mode and in the
    modes before it, phase by phase."""
    followed = [
        tuple(hkl) for hkl in (inputs.get("refine") or {}).get("followed") or []
    ]
    if not followed:
        return []
    runs = [*(earlier or {}).items(), (_mode_name(result), result)]
    header = ["Phase", "hkl", *(mode for mode, _ in runs)]
    rows = []
    phases = list(((result.get("exports") or {}).get("reflections") or {}).keys())
    tables = {}
    for mode, run in runs:
        for name, path in ((run.get("exports") or {}).get("reflections") or {}).items():
            if Path(path).is_file():
                tables[(mode, name)] = {
                    (row["h"], row["k"], row["l"]): row
                    for row in _read_reflections(path)
                }
    for name in phases:
        for hkl in followed:
            cells = [name, " ".join(str(i) for i in hkl)]
            for mode, _ in runs:
                row = tables.get((mode, name), {}).get(hkl)
                cells.append(
                    "—"
                    if row is None or not row["f_calc_squared"]
                    else f"{row['f_obs_squared'] / row['f_calc_squared']:.3f}"
                )
            rows.append(cells)
    return [
        "## Followed reflections",
        "",
        "F²obs / F²calc of each reflection followed, mode by mode.",
        "",
        *_table(header, rows),
    ]


def _mode_name(result: Mapping) -> str:
    return (result.get("method") or {}).get("mode") or "this mode"


# Structure


def _metric(cell: Mapping[str, float]) -> np.ndarray:
    a, b, c = (float(cell[key]) for key in ("length_a", "length_b", "length_c"))
    alpha, beta, gamma = (
        math.radians(float(cell[key]))
        for key in ("angle_alpha", "angle_beta", "angle_gamma")
    )
    return np.array(
        [
            [a * a, a * b * math.cos(gamma), a * c * math.cos(beta)],
            [a * b * math.cos(gamma), b * b, b * c * math.cos(alpha)],
            [a * c * math.cos(beta), b * c * math.cos(alpha), c * c],
        ]
    )


def coordinate_shift(
    atom: Mapping, start: Mapping, cell: Mapping[str, float]
) -> tuple[float, float | None]:
    """How far ``atom`` has moved from ``start``, in angstroms through the
    metric tensor of ``cell``, whole cells aside, and its esd from the
    atom's coordinate esds, the start taken as exact."""
    delta = np.array(
        [
            float(v) - float(s) - round(float(v) - float(s))
            for v, s in zip(atom["xyz"], start["xyz"])
        ]
    )
    metric = _metric(cell)
    distance = float(math.sqrt(max(delta @ metric @ delta, 0.0)))
    esds = [0.0 if e is None else float(e) for e in atom.get("xyz_esd") or [None] * 3]
    if not any(esds):
        return distance, None
    if distance == 0.0:
        return distance, float(
            math.sqrt(sum((metric[i, i] * esds[i] ** 2) for i in range(3)))
        )
    gradient = metric @ delta / distance
    return distance, float(
        math.sqrt(sum((gradient[i] * esds[i]) ** 2 for i in range(3)))
    )


def coordinate_lines(result: Mapping) -> list[str]:
    """Every atom's final position, occupancy and Uiso with esds, and its
    shift from the start model."""
    final = result.get("final") or {}
    model = result.get("start_model") or {}
    lines = ["## Coordinates", ""]
    for phase in final.get("phases") or []:
        name = phase.get("name")
        starts = {atom["label"]: atom for atom in model.get(name) or []}
        rows = []
        for atom in phase.get("atoms") or []:
            esds = atom.get("xyz_esd") or [None, None, None]
            start = starts.get(atom["label"])
            if start is None:
                shift = "not in the start model"
            else:
                distance, esd = coordinate_shift(atom, start, phase["cell"])
                shift = with_esd(distance, esd, 4) if esd else f"{distance:.4f}"
            uiso = atom.get("uiso")
            rows.append(
                [
                    atom["label"],
                    *(with_esd(v, e, 5) for v, e in zip(atom["xyz"], esds)),
                    with_esd(atom.get("occupancy", 1.0), atom.get("occupancy_esd"), 4),
                    "anisotropic"
                    if uiso is None
                    else with_esd(uiso, atom.get("uiso_esd"), 4),
                    shift,
                ]
            )
        lines += [
            f"Phase {name}, the shift from the start model in the refined cell:",
            "",
            *_table(
                ["Atom", "x", "y", "z", "Occupancy", "Uiso (Å²)", "Shift (Å)"], rows
            ),
            "",
        ]
    return lines


def occupancy_lines(result: Mapping) -> list[str]:
    """Each exchange group's constraint: the element, its atoms, the content
    over the sites held, and the content the final atoms give."""
    final = result.get("final") or {}
    atoms = {
        phase.get("name"): {atom["label"]: atom for atom in phase.get("atoms") or []}
        for phase in final.get("phases") or []
    }
    rows = []
    seen = set()
    for stage in accepted_stages(result):
        for name, constraints in (stage.get("occupancy_constraints") or {}).items():
            for constraint in constraints:
                key = (
                    name,
                    constraint.get("element"),
                    tuple(constraint.get("atoms", [])),
                )
                if key in seen:
                    continue
                seen.add(key)
                refined = sum(
                    float(m)
                    * float(atoms.get(name, {}).get(label, {}).get("occupancy", 0.0))
                    for m, label in zip(
                        constraint.get("multipliers", []), constraint.get("atoms", [])
                    )
                )
                rows.append(
                    [
                        name,
                        constraint.get("element"),
                        ", ".join(constraint.get("sites", [])),
                        ", ".join(constraint.get("atoms", [])),
                        f"{float(constraint.get('total', 0.0)):.5f}",
                        f"{refined:.5f}",
                    ]
                )
    lines = ["## Occupancies", ""]
    if not rows:
        return lines + ["No occupancy was traded between sites."]
    return lines + [
        (
            "Each exchanged element's content over its group's sites, the sum of "
            "multiplicity times occupancy, held at its start while the element is "
            "traded between the sites:"
        ),
        "",
        *_table(
            [
                "Phase",
                "Element",
                "Sites",
                "Atoms",
                "Held total per cell",
                "Final total",
            ],
            rows,
        ),
    ]


def bond_lines(result: Mapping, plans: Mapping[str, Mapping]) -> list[str]:
    """Every cation to anion distance of the final model, flagged outside
    the bond limits of its site's kind."""
    bonds = result.get("bonds") or {}
    lines = ["## Bond lengths", ""]
    if not bonds:
        error = result.get("bonds_error")
        return lines + [f"Not worked out{': ' + error if error else ''}."]
    for name, found in bonds.items():
        plan = plans.get(name) or {}
        kind_of = {
            str(atom["label"]): site["kind"]
            for site in plan.get("sites", [])
            for atom in site.get("atoms", [])
        }
        limits = plan.get("bond_limits") or {}
        anions = plan.get("anions")
        rows, flagged = [], []
        for bond in found:
            kind = kind_of.get(bond["centre"])
            low, high = limits.get(kind, (None, None))
            outside = low is not None and not low <= bond["distance"] <= high
            if outside:
                flagged.append(
                    f"{bond['centre_site']}–{bond['target']} {bond['distance']:.3f} Å"
                )
            rows.append(
                [
                    bond["centre_site"],
                    kind or "—",
                    bond["target"],
                    bond["count"],
                    with_esd(bond["distance"], bond.get("esd"), 4),
                    "—" if low is None else f"{low:g} to {high:g}",
                    "**outside**" if outside else "",
                ]
            )
        lines += [
            f"Phase {name}, to "
            + (_joined(anions) if anions else "the anions")
            + ", with the limits of each kind of site:",
            "",
            *_table(
                [
                    "Cation site",
                    "Kind",
                    "Anion",
                    "Count",
                    "Distance (Å)",
                    "Limits (Å)",
                    "Flag",
                ],
                rows,
            ),
            "",
            ("Flagged: " + "; ".join(flagged) + ".")
            if flagged
            else "No bond lies outside its limits.",
            "",
        ]
    return lines


def kind_prose(plans: Mapping[str, Mapping], mode: str) -> list[str]:
    """What the plan's kinds of site are, and in the coordinates mode the
    order they were refined in."""
    lines = []
    for name, plan in plans.items():
        kinds = list(plan.get("kinds") or [])
        if not kinds:
            continue
        described = _joined(
            [
                f"{kind} ({len(plan['kinds'][kind])} site{'' if len(plan['kinds'][kind]) == 1 else 's'})"
                for kind in kinds
            ]
        )
        line = f"Phase {name} has sites of the kinds {described}"
        if mode == "coordinates":
            freed = [
                kind
                for kind in kinds
                if any((plan.get("coordinates") or {}).get(kind, {}).values())
            ]
            line += (
                f"; the coordinates were refined kind by kind, {_joined(freed)}"
                if freed
                else "; none has a coordinate free"
            )
            origin = plan.get("origin")
            if origin:
                line += (
                    f", with {origin['name']} ({origin['kind']}) holding the origin along "
                    f"{plan.get('origin_axis')}"
                )
        groups = plan.get("exchange") or []
        if isinstance(groups, Mapping):
            groups = [groups]
        if mode == "occupancies" and groups:
            line += "; " + "; ".join(
                f"{_joined(group['elements'])} traded between the "
                f"{_kind_of(plan, group['sites'][0])} sites "
                f"{_joined(group['sites'])}"
                for group in groups
            )
        lines.append(line + ".")
    return lines


def _kind_of(plan: Mapping, site: str) -> str:
    for candidate in plan.get("sites", []):
        if candidate["name"] == site:
            return candidate["kind"]
    return "?"


# The write ups


def _title(result: Mapping, inputs: Mapping, heading: str) -> list[str]:
    key = (inputs.get("sample") or {}).get("key") or "sample"
    return [f"# {key}: {heading}", ""]


def lebail_markdown(
    result: Mapping,
    inputs: Mapping | None = None,
    plan=None,
    earlier: Mapping[str, Mapping] | None = None,
) -> str:
    """The write up of a Le Bail mode."""
    inputs = _inputs(result, inputs)
    phases = [structure.get("key") for structure in inputs.get("structures") or []]
    if not phases:
        phases = [
            phase.get("name")
            for phase in (result.get("final") or {}).get("phases") or []
        ]
    lines = [
        *_title(result, inputs, "Le Bail extraction"),
        (
            f"Le Bail extraction of {_joined(phases)}: the intensities of every phase "
            "extracted freely, no atomic parameter refined."
        ),
        "",
        *settings_lines(result, inputs),
        "",
        *stage_outcome_lines(result),
        "",
        *cell_lines(result, inputs),
        *zero_lines(result, inputs),
        "",
        *broadening_lines(result, lebail=True),
        "",
        *weight_fraction_lines(result),
        "",
        *followed_lines(result, inputs, earlier),
    ]
    return _finish(lines)


def fixed_atoms_markdown(
    result: Mapping,
    inputs: Mapping | None = None,
    plan=None,
    earlier: Mapping[str, Mapping] | None = None,
) -> str:
    """The write up of a Rietveld refinement with the atoms fixed."""
    inputs = _inputs(result, inputs)
    plans = _plans(plan)
    lines = [
        *_title(result, inputs, "Rietveld refinement, atoms fixed"),
        (
            "Each phase's structure put in at its nominal composition, the atoms "
            "fixed but for one overall Uiso."
        ),
        "",
        *kind_prose(plans, "fixed_atoms"),
        "",
        *settings_lines(result, inputs),
        "",
        *stage_outcome_lines(result),
        "",
        *cell_lines(result, inputs),
        *zero_lines(result, inputs),
        "",
        *broadening_lines(result, lebail=False),
        "",
        *weight_fraction_lines(result),
        "",
        *misfit_lines(result, plans, (earlier or {}).get("lebail")),
        *followed_lines(result, inputs, earlier),
        "",
        *sanity_lines(result),
        "",
        *undetermined_lines(result),
    ]
    return _finish(lines)


def structure_markdown(
    result: Mapping,
    inputs: Mapping | None = None,
    plan=None,
    earlier: Mapping[str, Mapping] | None = None,
) -> str:
    """The write up of a refinement of the coordinates or the occupancies,
    as the result's recorded mode says."""
    inputs = _inputs(result, inputs)
    plans = _plans(plan)
    mode = _mode_name(result)
    heading = {
        "coordinates": "Rietveld refinement, coordinates",
        "occupancies": "Rietveld refinement, occupancies",
    }.get(mode, "Rietveld refinement of the structure")
    lines = [
        *_title(result, inputs, heading),
        *kind_prose(plans, mode),
        "",
        *settings_lines(result, inputs),
        "",
        *stage_outcome_lines(result),
        "",
        *cell_lines(result, inputs),
        *zero_lines(result, inputs),
        "",
        *broadening_lines(result, lebail=False),
        "",
        *weight_fraction_lines(result),
        "",
        *coordinate_lines(result),
        *undetermined_lines(result),
        "",
        *(occupancy_lines(result) if mode == "occupancies" else []),
        "",
        *bond_lines(result, plans),
        *misfit_lines(result, plans, (earlier or {}).get("lebail")),
        *followed_lines(result, inputs, earlier),
        "",
        *sanity_lines(result),
    ]
    return _finish(lines)


MARKDOWN = {
    "lebail": lebail_markdown,
    "fixed_atoms": fixed_atoms_markdown,
    "coordinates": structure_markdown,
    "occupancies": structure_markdown,
}


def _finish(lines: Sequence[str]) -> str:
    """The lines as markdown, runs of blank lines made one."""
    out: list[str] = []
    for line in lines:
        if line == "" and (not out or out[-1] == ""):
            continue
        out.append(line)
    while out and out[-1] == "":
        out.pop()
    return "\n".join(out) + "\n"
