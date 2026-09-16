"""Refinement settings for samples and their reference structures, from TOML.

A settings file holds two tables of tables, read and checked by
:func:`load_config`:

``samples``
    One table per sample, named as the caller likes (the name of its results
    folder, say), with the keys

    - ``id``: the sample's identifier;
    - ``scan``: its scan file, relative to the data folder;
    - ``composition``: its nominal composition, atoms of each element per
      formula unit, ``{Sr = 0.40, Ba = 0.50, ...}``;
    - ``structure``: the name of its reference structure's table;
    - ``start_cell``: where the cell starts, ``{file, model}`` for a lattice
      refinement result, or the cell parameters given outright, exactly
      those of the structure's crystal system (``{a, c}`` for a tetragonal
      cell, say), each greater than 0;
    - ``two_theta``: the range refined, ``[low, high]`` in degrees;
    - ``background``: ``{function, terms}``;
    - ``refine_microstrain``: whether the microstrain is refined;
    - ``notes``: free text;

    and optionally ``followed_reflections``, a list of hkl labels,
    ``trials``, ``{runs = [{low, terms}, ...], followed = [hkl, ...]}``
    with ``low`` left out for the whole scan, and ``write_up``, a table of
    tables of text, by mode and section.

``structures``
    One table per reference structure, with the keys

    - ``cif``: its CIF, relative to the caller's root folder;
    - ``label``, ``phase_name``: how it is named in write ups and in the
      GSAS-II project;
    - ``space_group``;
    - ``sites``: a list of ``{atoms = {label = element, ...}, wyckoff,
      kind}``, one per site, named by its first atom, its kind a short
      label of the caller's choosing (:data:`SITE_KIND`), such as A, B or O;
    - ``uiso_groups``: a list of ``{name, sites}``, every site in one;
    - ``composition``: ``{added = {element = host element, ...}}``, how a
      nominal composition goes on the sites: every element the sites hold
      is scaled by one factor over them, which keeps its distribution, and
      each added element goes on every site of its host in proportion to
      the host's occupancy there (see
      :func:`xrdkit.structure.composition_edits`);

    and optionally ``formula_units``, formula units per cell (Z);
    ``crystal_system``, one of :data:`xrdkit.library.CRYSTAL_SYSTEMS`;
    ``origin``, ``{site, axis}``, the site coordinate that fixes the
    origin along a polar axis, left out when nothing is to be held;
    ``exchange``, ``{elements, sites}``, the elements whose occupancies are
    traded between the sites, left out when none are;
    ``free_coordinates``, the coordinates to refine by Wyckoff position,
    ``{"8d" = "xyz", "2a" = "z", ...}``, every coordinate for a position not
    named; and ``bond_limits``, ``{kind = {min, max}}`` in angstroms, the
    range outside which a cation to anion bond from a site of that kind is
    flagged. Left out, ``formula_units``, ``crystal_system``, ``origin`` and
    ``exchange`` are None.

    A structure may instead name an entry of the structure library,
    ``library = "ttb/P4bm"``, and give ``atoms``, a table of the CIF's atoms
    on each of the entry's sites by its label, ``{A1 = {Sr1 = "Sr"}, ...}``,
    in place of ``sites``. Its space group, sites, Wyckoff positions, kinds,
    free coordinates, Uiso groups, formula units, crystal system and origin
    (the entry's origin site along its polar axis) then come from the entry,
    and any of those keys given explicitly overrides the entry's value (an
    explicit ``sites`` takes the place of ``atoms``, and then the entry's
    free coordinates, Uiso groups and origin, which name its own sites, do
    not apply). The sites are still named by their first atoms, and every
    reference to a site, in ``uiso_groups``, ``origin`` and ``exchange``,
    may give the entry's label for it instead.

A top level ``unsettled`` table, optional, gives for each mode of the
caller's pipeline, by name, the rule for a stage that has not settled:
``"accept"`` (keep it and go on) or ``"reject"`` (roll it back); a sample's
own ``unsettled`` overrides it mode by mode, and each sample carries the
two merged as its ``unsettled``.

Every sample's composition is checked against its structure's sites: each
element must be held by a site or added by the rule, and, when the
structure gives its formula units, no kind of site may be given more atoms
per cell than it has positions.
"""

from __future__ import annotations

import math
import re
import tomllib
from collections.abc import Mapping
from pathlib import Path

__all__ = [
    "SITE_KIND",
    "ConfigError",
    "check_composition",
    "load_config",
    "read_library_atoms",
    "read_sites",
    "sample_settings",
    "validate_config",
    "wyckoff_multiplicity",
]

# A kind of site: a short label, letters, digits and underscores starting with
# a letter, at most eight. The kinds are whatever a structure's sites, or its
# library entry's, declare.
SITE_KIND = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,7}")

SAMPLE_REQUIRED = (
    "id",
    "scan",
    "composition",
    "structure",
    "start_cell",
    "two_theta",
    "background",
    "refine_microstrain",
    "notes",
)
SAMPLE_OPTIONAL = ("followed_reflections", "trials", "write_up", "unsettled")
# The rules for a stage that has not settled.
UNSETTLED_RULES = ("accept", "reject")
STRUCTURE_REQUIRED = (
    "cif",
    "label",
    "phase_name",
    "space_group",
    "sites",
    "uiso_groups",
    "composition",
)
STRUCTURE_OPTIONAL = (
    "formula_units",
    "crystal_system",
    "origin",
    "exchange",
    "free_coordinates",
    "bond_limits",
    "library",
    "atoms",
)
# Required of a structure that names a library entry; the rest may come from it.
LIBRARY_REQUIRED = ("cif", "label", "phase_name", "composition")
LIBRARY_OPTIONAL = tuple(
    key
    for key in STRUCTURE_REQUIRED + STRUCTURE_OPTIONAL
    if key not in LIBRARY_REQUIRED
)
# The crystal system a start cell's parameters tell, when the structure does
# not say.
INFERRED_SYSTEMS = {
    frozenset({"a"}): "cubic",
    frozenset({"a", "c"}): "tetragonal",
    frozenset({"a", "b", "c"}): "orthorhombic",
}
# The coordinate along each cell axis.
AXIS_COORDINATE = {"a": "x", "b": "y", "c": "z"}

ELEMENT = re.compile(r"[A-Z][a-z]?")
WYCKOFF = re.compile(r"(\d+)([a-z])")
# A composition may put this much more on a kind of site than it has
# positions, for rounding.
CAPACITY_TOLERANCE = 1.0e-9


class ConfigError(ValueError):
    """A settings file with a key missing or unknown, a value of the wrong
    kind, or values that do not fit together."""


def _fail(where: str, message: str) -> ConfigError:
    return ConfigError(f"{where}: {message}")


def _table(value: object, where: str) -> Mapping:
    if not isinstance(value, Mapping):
        raise _fail(where, f"must be a table, not {type(value).__name__}")
    return value


def _keys(table: object, where: str, required: tuple, optional: tuple = ()) -> Mapping:
    table = _table(table, where)
    missing = [key for key in required if key not in table]
    if missing:
        raise _fail(where, "missing key " + ", ".join(repr(key) for key in missing))
    unknown = sorted(set(table) - set(required) - set(optional))
    if unknown:
        raise _fail(
            where,
            "unknown key "
            + ", ".join(repr(key) for key in unknown)
            + "; the keys are "
            + ", ".join(required + optional),
        )
    return table


def _string(value: object, where: str, empty: bool = False) -> str:
    if not isinstance(value, str) or not (empty or value.strip()):
        raise _fail(
            where, "must be a non-empty string" if not empty else "must be text"
        )
    return value


def _number(value: object, where: str, minimum: float | None = None) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise _fail(where, f"must be a number, not {value!r}")
    if minimum is not None and value < minimum:
        raise _fail(where, f"must be at least {minimum:g}, not {value!r}")
    return float(value)


def _positive(value: object, where: str) -> float:
    number = _number(value, where)
    if number <= 0.0:
        raise _fail(where, f"must be greater than 0, not {value!r}")
    return number


def _integer(value: object, where: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise _fail(
            where, f"must be a whole number of at least {minimum}, not {value!r}"
        )
    return value


def _strings(value: object, where: str, minimum: int = 1) -> list[str]:
    if not isinstance(value, list) or len(value) < minimum:
        raise _fail(where, f"must be a list of at least {minimum} strings")
    for index, item in enumerate(value):
        _string(item, f"{where}[{index}]")
    if len(set(value)) != len(value):
        raise _fail(where, f"lists an entry twice: {value}")
    return list(value)


def _element(value: object, where: str) -> str:
    if not isinstance(value, str) or not ELEMENT.fullmatch(value):
        raise _fail(where, f"must be an element symbol such as 'Sr', not {value!r}")
    return value


def wyckoff_multiplicity(symbol: str) -> int:
    """The multiplicity of a Wyckoff position written as ``"4c"``.

    Raises
    ------
    ValueError
        If ``symbol`` is not a number followed by one lower case letter.
    """
    match = WYCKOFF.fullmatch(str(symbol))
    if not match:
        raise ValueError(f"not a Wyckoff position such as '4c': {symbol!r}")
    return int(match.group(1))


# Structures


def read_sites(value: object, where: str) -> list[dict]:
    """A structure's sites, a list of ``{atoms = {label = element, ...},
    wyckoff, kind}`` tables, each as ``{"name", "atoms", "wyckoff", "kind"}``,
    named by its first atom; shared with :mod:`xrdkit.project`.

    Raises
    ------
    ConfigError
        If a site table has a key missing or unknown, an atom is on two
        sites, an element or Wyckoff position is not one, or a kind is not
        a short label (:data:`SITE_KIND`); the message starts with ``where``.
    """
    if not isinstance(value, list) or not value:
        raise _fail(where, "must be a list of site tables")
    labels: set[str] = set()
    return [
        _site(table, f"{where}[{index}]", labels) for index, table in enumerate(value)
    ]


def _site(value: object, at: str, labels: set[str]) -> dict:
    """One site table, its atoms' labels added to ``labels``."""
    table = _keys(value, at, ("atoms", "wyckoff", "kind"))
    atoms = _table(table["atoms"], f"{at}.atoms")
    if not atoms:
        raise _fail(f"{at}.atoms", "must name at least one atom")
    for label, element in atoms.items():
        _element(element, f"{at}.atoms.{label}")
        if label in labels:
            raise _fail(f"{at}.atoms", f"atom {label!r} is on another site too")
        labels.add(label)
    wyckoff = _string(table["wyckoff"], f"{at}.wyckoff")
    try:
        wyckoff_multiplicity(wyckoff)
    except ValueError as error:
        raise _fail(f"{at}.wyckoff", str(error)) from None
    if not isinstance(table["kind"], str) or not SITE_KIND.fullmatch(table["kind"]):
        raise _fail(
            f"{at}.kind",
            "must be a short label, letters, digits and underscores starting with "
            f"a letter, at most eight, not {table['kind']!r}",
        )
    return {
        "name": next(iter(atoms)),
        "atoms": dict(atoms),
        "wyckoff": wyckoff,
        "kind": table["kind"],
    }


def _library_entry(value: object, where: str):
    """The structure library entry ``value`` names."""
    from xrdkit.library import load_entry  # the library imports this module

    name = _string(value, where)
    try:
        return load_entry(name)
    except ValueError as error:
        raise _fail(where, str(error)) from None


def read_library_atoms(
    value: object, where: str, entry, complete: bool = True
) -> tuple[list[dict], dict]:
    """The sites of ``entry``, a structure library entry, holding the atoms
    ``value`` puts on each by the entry's site label, ``{A1 = {Sr1 = "Sr"},
    ...}``, each as :func:`read_sites` gives it with its ``label`` too; and
    the site each label names. Every site of the entry must be given unless
    ``complete`` is false, when only those given are returned. Shared with
    :mod:`xrdkit.project`.

    Raises
    ------
    ConfigError
        If a label is not one of the entry's sites, a site is missing when
        ``complete``, or its atoms do not check out as :func:`read_sites`
        checks them; the message starts with ``where``.
    """
    atoms = _table(value, where)
    known = [site.label for site in entry.sites]
    unknown = [label for label in atoms if label not in known]
    if unknown:
        raise _fail(
            where,
            f"no site {', '.join(map(repr, unknown))} in {entry.name}; its sites "
            f"are {', '.join(known)}",
        )
    missing = [label for label in known if label not in atoms]
    if missing and complete:
        raise _fail(where, f"no atoms for site {', '.join(missing)} of {entry.name}")
    labels: set[str] = set()
    sites = []
    for site in (site for site in entry.sites if site.label in atoms):
        table = {"atoms": atoms[site.label], "wyckoff": site.wyckoff, "kind": site.kind}
        sites.append(
            {**_site(table, f"{where}.{site.label}", labels), "label": site.label}
        )
    aliases = {site["label"]: site["name"] for site in sites}
    names = [site["name"] for site in sites]
    for label, name in aliases.items():
        if label in names and label != name:
            raise _fail(
                where,
                f"{label}, a site label of {entry.name}, is the first atom of "
                "another site; relabel that atom",
            )
    return sites, aliases


def _site_names(
    value: object,
    where: str,
    names: list[str],
    aliases: Mapping[str, str] | None = None,
) -> list[str]:
    """The sites ``value`` lists, each by its name or, through ``aliases``,
    by its library label, resolved to their names."""
    aliases = aliases or {}
    listed = _strings(value, where)
    unknown = [name for name in listed if name not in names and name not in aliases]
    if unknown:
        known = f"sites are named by their first atom: {', '.join(names)}"
        if aliases:
            known += f", or by their library label: {', '.join(aliases)}"
        raise _fail(where, f"no site named {', '.join(map(repr, unknown))}; {known}")
    resolved = [name if name in names else aliases[name] for name in listed]
    if len(set(resolved)) != len(resolved):
        raise _fail(where, f"names a site twice: {', '.join(listed)}")
    return resolved


def _crystal_system(table: Mapping, where: str, entry) -> str | None:
    """The structure's crystal system: its own, else its entry's, else None."""
    from xrdkit.library import CRYSTAL_SYSTEMS  # the library imports this module

    if "crystal_system" not in table:
        return entry.crystal_system if entry else None
    system = table["crystal_system"]
    if system not in CRYSTAL_SYSTEMS:
        raise _fail(
            f"{where}.crystal_system",
            f"must be one of {', '.join(CRYSTAL_SYSTEMS)}, not {system!r}",
        )
    return system


def _origin(
    table: Mapping,
    where: str,
    names: list[str],
    aliases: Mapping[str, str],
    entry,
    from_entry: bool,
) -> dict | None:
    """The site coordinate that fixes the origin: the structure's own, else
    its entry's origin site along its polar axis, else None."""
    if "origin" not in table:
        if from_entry and entry.origin_site and entry.polar_axis:
            return {
                "site": aliases[entry.origin_site],
                "axis": AXIS_COORDINATE[entry.polar_axis],
            }
        return None
    origin = _keys(table["origin"], f"{where}.origin", ("site", "axis"))
    (site,) = _site_names([origin["site"]], f"{where}.origin.site", names, aliases)
    if origin["axis"] not in ("x", "y", "z"):
        raise _fail(
            f"{where}.origin.axis", f"must be x, y or z, not {origin['axis']!r}"
        )
    return {"site": site, "axis": origin["axis"]}


def _exchange(
    value: object,
    where: str,
    sites: list[dict],
    kind_of: Mapping[str, str],
    aliases: Mapping[str, str],
) -> dict:
    """The elements whose occupancies are traded, and the sites they are
    traded between."""
    exchange = _keys(value, where, ("elements", "sites"))
    elements = _strings(exchange["elements"], f"{where}.elements", 2)
    for index, element in enumerate(elements):
        _element(element, f"{where}.elements[{index}]")
    names = list(kind_of)
    between = _site_names(exchange["sites"], f"{where}.sites", names, aliases)
    if len(between) < 2:
        raise _fail(f"{where}.sites", "must name at least two sites")
    if len({kind_of[name] for name in between}) > 1:
        raise _fail(
            f"{where}.sites",
            "must all be of one kind, not "
            + ", ".join(f"{name} ({kind_of[name]})" for name in between),
        )
    on_sites = {
        element
        for site in sites
        if site["name"] in between
        for element in site["atoms"].values()
    }
    absent = [element for element in elements if element not in on_sites]
    if absent:
        raise _fail(
            f"{where}.elements",
            f"{', '.join(absent)} on none of the sites {', '.join(between)}",
        )
    return {"elements": elements, "sites": between}


def _structure(table: object, where: str) -> dict:
    table = _table(table, where)
    entry = None
    if "library" in table:
        _keys(table, where, LIBRARY_REQUIRED, LIBRARY_OPTIONAL)
        entry = _library_entry(table["library"], f"{where}.library")
    else:
        _keys(table, where, STRUCTURE_REQUIRED, STRUCTURE_OPTIONAL)
        if "atoms" in table:
            raise _fail(
                f"{where}.atoms",
                "goes with library, the entry whose sites the atoms are on; "
                "give sites instead",
            )
    structure = {
        key: _string(table[key], f"{where}.{key}")
        for key in ("cif", "label", "phase_name")
    }
    structure["library"] = entry.name if entry else None
    structure["space_group"] = (
        _string(table["space_group"], f"{where}.space_group")
        if "space_group" in table
        else entry.space_group
    )
    structure["formula_units"] = (
        _integer(table["formula_units"], f"{where}.formula_units", 1)
        if "formula_units" in table
        else (entry.z if entry else None)
    )
    structure["crystal_system"] = _crystal_system(table, where, entry)

    # The entry's sites, and what names them, unless sites are given outright.
    aliases: dict[str, str] = {}
    from_entry = entry is not None and "sites" not in table
    if not from_entry:
        sites = read_sites(table["sites"], f"{where}.sites")
    elif "atoms" in table:
        sites, aliases = read_library_atoms(table["atoms"], f"{where}.atoms", entry)
    else:
        raise _fail(
            where,
            f"missing key 'atoms'; give the CIF's atoms on each site of "
            f"{entry.name}, or sites",
        )
    names = [site["name"] for site in sites]
    kind_of = {site["name"]: site["kind"] for site in sites}
    structure["sites"] = sites

    if "free_coordinates" in table:
        free = _table(table["free_coordinates"], f"{where}.free_coordinates")
    else:
        free = {}
    for wyckoff, axes in free.items():
        at = f"{where}.free_coordinates.{wyckoff}"
        try:
            wyckoff_multiplicity(wyckoff)
        except ValueError as error:
            raise _fail(at, str(error)) from None
        if axes != "all" and not (
            isinstance(axes, str)
            and axes
            and set(axes) <= set("xyz")
            and len(set(axes)) == len(axes)
        ):
            raise _fail(at, f"must be 'all' or some of 'xyz', not {axes!r}")
    if "free_coordinates" not in table and from_entry:
        # The entry's own, "" for a site its symmetry fixes.
        free = {site.wyckoff: "".join(site.free) for site in entry.sites}
    structure["free_coordinates"] = dict(free)

    if "uiso_groups" in table:
        groups = table["uiso_groups"]
    elif from_entry:
        by_group: dict[str, list[str]] = {}
        for site in entry.sites:
            if site.uiso_group is not None:
                by_group.setdefault(site.uiso_group, []).append(site.label)
        groups = [
            {"name": name, "sites": members} for name, members in by_group.items()
        ]
    else:
        raise _fail(where, "missing key 'uiso_groups'")
    if not isinstance(groups, list) or not groups:
        raise _fail(f"{where}.uiso_groups", "must be a list of {name, sites} tables")
    placed: dict[str, str] = {}
    structure["uiso_groups"] = []
    for index, group in enumerate(groups):
        at = f"{where}.uiso_groups[{index}]"
        _keys(group, at, ("name", "sites"))
        name = _string(group["name"], f"{at}.name")
        members = _site_names(group["sites"], f"{at}.sites", names, aliases)
        for member in members:
            if member in placed:
                raise _fail(
                    at, f"site {member} is in Uiso group {placed[member]!r} already"
                )
            placed[member] = name
        structure["uiso_groups"].append({"name": name, "sites": members})
    outside = [name for name in names if name not in placed]
    if outside:
        raise _fail(
            f"{where}.uiso_groups", f"sites {', '.join(outside)} are in no Uiso group"
        )

    structure["origin"] = _origin(table, where, names, aliases, entry, from_entry)
    structure["exchange"] = (
        _exchange(table["exchange"], f"{where}.exchange", sites, kind_of, aliases)
        if "exchange" in table
        else None
    )

    rule = _keys(table["composition"], f"{where}.composition", ("added",))
    added = _table(rule["added"], f"{where}.composition.added")
    held = {element for site in sites for element in site["atoms"].values()}
    for element, host in added.items():
        at = f"{where}.composition.added.{element}"
        _element(element, at)
        _element(host, at)
        if element in held:
            raise _fail(
                at,
                f"{element} is on the sites already; only an element "
                "they lack is added",
            )
        if host not in held:
            raise _fail(at, f"its host {host} is on none of the sites")
    structure["composition"] = {"added": dict(added)}

    limits = _table(table.get("bond_limits", {}), f"{where}.bond_limits")
    structure["bond_limits"] = {}
    kinds = list(dict.fromkeys(site["kind"] for site in sites))
    for kind, bounds in limits.items():
        at = f"{where}.bond_limits.{kind}"
        if kind not in kinds:
            raise _fail(
                at,
                f"not a kind of the structure's sites; its kinds are {', '.join(kinds)}",
            )
        _keys(bounds, at, (), ("min", "max"))
        checked = {
            key: _number(value, f"{at}.{key}", 0.0) for key, value in bounds.items()
        }
        if checked.get("min", 0.0) > checked.get("max", math.inf):
            raise _fail(at, "min is above max")
        structure["bond_limits"][kind] = checked
    return structure


# Samples


def check_composition(
    composition: Mapping[str, float], structure: Mapping, where: str = "composition"
) -> None:
    """Check that ``composition``, atoms per formula unit, fits the sites of
    ``structure``, a structure table as :func:`load_config` returns it.

    Every element the sites hold must be in the composition (at 0 if it is
    absent), and every element of the composition must be on a site or added
    by the structure's composition rule. An element takes the kinds of site
    it is on, or its host's; kinds that share an element are counted
    together, and none may be given more atoms per cell than the
    multiplicities of their sites add up to; a structure that gives no
    formula units cannot be checked for that, and is not.

    Raises
    ------
    ConfigError
        If the composition does not fit.
    """
    sites = structure["sites"]
    added = structure["composition"]["added"]
    kinds_of: dict[str, set[str]] = {}
    for site in sites:
        for element in site["atoms"].values():
            kinds_of.setdefault(element, set()).add(site["kind"])
    held = sorted(kinds_of)
    for element, host in added.items():
        kinds_of[element] = set(kinds_of[host])
    missing = [element for element in held if element not in composition]
    if missing:
        raise _fail(
            where,
            f"lacks {', '.join(missing)}, which the sites of the structure hold; "
            "give 0 for an element that is absent",
        )
    foreign = [element for element in composition if element not in kinds_of]
    if foreign:
        raise _fail(
            where,
            f"has {', '.join(foreign)}, which no site of the structure holds and "
            "its composition rule does not add",
        )

    # Kinds of site that share an element are one pool of positions.
    pools: list[set[str]] = []
    for kinds in kinds_of.values():
        joined = set(kinds)
        for pool in [pool for pool in pools if pool & joined]:
            joined |= pool
            pools.remove(pool)
        pools.append(joined)
    units = structure["formula_units"]
    if units is None:
        return
    for pool in pools:
        on = [element for element in composition if kinds_of[element] <= pool]
        content = sum(units * composition[element] for element in on)
        pool_sites = [site for site in sites if site["kind"] in pool]
        capacity = sum(wyckoff_multiplicity(site["wyckoff"]) for site in pool_sites)
        if content > capacity + CAPACITY_TOLERANCE:
            raise _fail(
                where,
                f"puts {content:g} atoms per cell ({units} formula units of "
                + ", ".join(f"{element} {composition[element]:g}" for element in on)
                + f") on the {'/'.join(sorted(pool))} sites, which have "
                f"{capacity} positions ("
                + ", ".join(f"{site['name']} {site['wyckoff']}" for site in pool_sites)
                + ")",
            )


def _start_cell(value: object, where: str, structure: Mapping) -> dict:
    """``{file, model}``, or the cell parameters of the crystal system of
    ``structure``, or, when it has none, of the one the parameters tell."""
    from xrdkit.library import CELL_PARAMETERS  # the library imports this module

    table = _table(value, where)
    if "file" in table or "model" in table:
        _keys(table, where, ("file", "model"))
        return {
            "file": _string(table["file"], f"{where}.file"),
            "model": _string(table["model"], f"{where}.model"),
        }
    if not table:
        raise _fail(
            where,
            "must give file and model (a lattice refinement result) or the cell "
            "parameters of a crystal system",
        )
    every = CELL_PARAMETERS["triclinic"]
    unknown = [key for key in table if key not in every]
    if unknown:
        raise _fail(
            where,
            f"unknown key {', '.join(map(repr, unknown))}; the cell parameters are "
            + ", ".join(every),
        )
    given = [key for key in every if key in table]
    system = structure["crystal_system"]
    if system is None:
        system = INFERRED_SYSTEMS.get(frozenset(given))
        if system is None:
            raise _fail(
                where,
                f"cannot tell the crystal system from {', '.join(given)}; give "
                f"crystal_system in structures.{structure['name']}",
            )
    elif set(given) != set(CELL_PARAMETERS[system]):
        raise _fail(
            where,
            f"a {system} cell (the crystal_system of structures."
            f"{structure['name']}) has {', '.join(CELL_PARAMETERS[system])}, not "
            + ", ".join(given),
        )
    return {key: _positive(table[key], f"{where}.{key}") for key in given}


def _trials(value: object, where: str) -> dict:
    table = _keys(value, where, ("runs",), ("followed",))
    runs = table["runs"]
    if not isinstance(runs, list) or not runs:
        raise _fail(f"{where}.runs", "must be a list of {low, terms} tables")
    checked = []
    for index, run in enumerate(runs):
        at = f"{where}.runs[{index}]"
        _keys(run, at, ("terms",), ("low",))
        checked.append(
            {
                "low": None
                if "low" not in run
                else _number(run["low"], f"{at}.low", 0.0),
                "terms": _integer(run["terms"], f"{at}.terms", 1),
            }
        )
    followed = table.get("followed", [])
    return {
        "runs": checked,
        "followed": _strings(followed, f"{where}.followed") if followed else [],
    }


def _unsettled(value: object, where: str) -> dict[str, str]:
    rules = _table(value, where)
    for mode, rule in rules.items():
        if rule not in UNSETTLED_RULES:
            raise _fail(
                f"{where}.{mode}",
                f"must be one of {', '.join(UNSETTLED_RULES)}, not {rule!r}",
            )
    return dict(rules)


def _sample(table: object, where: str, structures: Mapping[str, dict]) -> dict:
    _keys(table, where, SAMPLE_REQUIRED, SAMPLE_OPTIONAL)
    sample = {
        "id": _string(table["id"], f"{where}.id"),
        "scan": _string(table["scan"], f"{where}.scan"),
        "notes": _string(table["notes"], f"{where}.notes", empty=True),
    }
    composition = _table(table["composition"], f"{where}.composition")
    if not composition:
        raise _fail(f"{where}.composition", "names no element")
    sample["composition"] = {
        _element(element, f"{where}.composition"): _number(
            value, f"{where}.composition.{element}", 0.0
        )
        for element, value in composition.items()
    }
    name = _string(table["structure"], f"{where}.structure")
    if name not in structures:
        raise _fail(
            f"{where}.structure",
            f"no structure {name!r} under structures; there are "
            + (", ".join(structures) or "none"),
        )
    sample["structure"] = name
    sample["start_cell"] = _start_cell(
        table["start_cell"], f"{where}.start_cell", structures[name]
    )

    two_theta = table["two_theta"]
    if not isinstance(two_theta, list) or len(two_theta) != 2:
        raise _fail(f"{where}.two_theta", "must be [low, high] in degrees")
    low, high = (
        _number(value, f"{where}.two_theta[{index}]", 0.0)
        for index, value in enumerate(two_theta)
    )
    if not low < high <= 180.0:
        raise _fail(
            f"{where}.two_theta",
            f"must rise from low to high within 180°, not {two_theta}",
        )
    sample["two_theta"] = (low, high)

    background = _keys(
        table["background"], f"{where}.background", ("function", "terms")
    )
    sample["background"] = {
        "function": _string(background["function"], f"{where}.background.function"),
        "terms": _integer(background["terms"], f"{where}.background.terms", 1),
    }
    if not isinstance(table["refine_microstrain"], bool):
        raise _fail(f"{where}.refine_microstrain", "must be true or false")
    sample["refine_microstrain"] = table["refine_microstrain"]

    followed = table.get("followed_reflections", [])
    sample["followed_reflections"] = (
        _strings(followed, f"{where}.followed_reflections") if followed else []
    )
    sample["trials"] = (
        _trials(table["trials"], f"{where}.trials") if "trials" in table else None
    )
    write_up = _table(table.get("write_up", {}), f"{where}.write_up")
    sample["write_up"] = {}
    for mode, sections in write_up.items():
        sections = _table(sections, f"{where}.write_up.{mode}")
        sample["write_up"][mode] = {
            section: _string(text, f"{where}.write_up.{mode}.{section}", empty=True)
            for section, text in sections.items()
        }
    check_composition(sample["composition"], structures[name], f"{where}.composition")
    return sample


# Files


def validate_config(data: Mapping, source: str = "settings") -> dict:
    """Check settings read from a TOML file and return them with the
    optional keys filled in, each sample and structure carrying its table
    name as ``name`` and each site its first atom's label as ``name``.

    Raises
    ------
    ConfigError
        On the first key missing or unknown, value of the wrong kind, or
        composition that does not fit its structure, naming the file and
        the table.
    """
    try:
        _keys(data, "top level", ("samples", "structures"), ("unsettled",))
        unsettled = _unsettled(data.get("unsettled", {}), "unsettled")
        structures = {}
        for name, table in _table(data["structures"], "structures").items():
            structures[name] = {"name": name, **_structure(table, f"structures.{name}")}
        samples = {}
        ids: dict[str, str] = {}
        for name, table in _table(data["samples"], "samples").items():
            sample = {"name": name, **_sample(table, f"samples.{name}", structures)}
            own = _unsettled(table.get("unsettled", {}), f"samples.{name}.unsettled")
            sample["unsettled"] = {**unsettled, **own}
            if sample["id"] in ids:
                raise _fail(
                    f"samples.{name}.id",
                    f"{sample['id']!r} is the id of samples.{ids[sample['id']]} too",
                )
            ids[sample["id"]] = name
            samples[name] = sample
    except ConfigError as error:
        raise ConfigError(f"{source}: {error}") from None
    return {
        "source": source,
        "unsettled": unsettled,
        "samples": samples,
        "structures": structures,
    }


def load_config(path: str | Path) -> dict:
    """Read the settings file at ``path`` with tomllib and check it (see
    :func:`validate_config` and the module notes for its layout).

    Raises
    ------
    ConfigError
        If the file is not valid TOML or its settings do not check out.
    FileNotFoundError
        If there is no such file.
    """
    path = Path(path)
    with path.open("rb") as handle:
        try:
            data = tomllib.load(handle)
        except tomllib.TOMLDecodeError as error:
            raise ConfigError(f"{path}: not valid TOML: {error}") from None
    return validate_config(data, str(path))


def sample_settings(config: Mapping, identifier: str) -> tuple[dict, dict]:
    """The sample of ``config`` whose id, or table name, is ``identifier``,
    and its structure.

    Raises
    ------
    ConfigError
        If no sample has that id or name.
    """
    samples = config["samples"]
    for name, sample in samples.items():
        if identifier in (sample["id"], name):
            return sample, config["structures"][sample["structure"]]
    raise ConfigError(
        f"{config.get('source', 'settings')}: no sample {identifier!r}; the samples are "
        + ", ".join(f"{sample['id']} ({name})" for name, sample in samples.items())
    )
