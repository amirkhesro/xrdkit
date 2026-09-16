"""The structure library shipped with xrdkit, one TOML file per entry.

An entry describes a structure type apart from any one sample. The files are
package data at ``xrdkit/structures/<family>/<name>.toml``, and each entry is
named ``"<family>/<name>"``, such as ``"ttb/P4bm"``. A file holds two parts,
read and checked by :func:`load_entry`:

``[entry]``
    - ``name``: the entry's name, which must match its path;
    - ``family``: the structure family, as text;
    - ``crystal_system``: one of :data:`CRYSTAL_SYSTEMS`;
    - ``space_group``: its Hermann-Mauguin symbol;
    - ``cell_parameters``: the cell parameters the crystal system leaves
      free, exactly as :data:`CELL_PARAMETERS` lists them;
    - ``z``: formula units per cell;
    - ``reference``: the structure the entry was taken from;

    and optionally ``setting``, the setting or origin choice (``""`` for the
    standard one); ``polar_axis``, ``"a"``, ``"b"`` or ``"c"``, the axis along
    which symmetry leaves the origin free; and ``origin_site``, the label of
    the site whose coordinate along that axis is held to fix the origin.

``[[sites]]``
    One table per site, with ``label``, unique within the entry, ``kind``
    and ``wyckoff``, the Wyckoff position written as ``"4c"``, and
    optionally ``free``, the coordinates the site leaves free, drawn from
    x, y and z, ``uiso_group``, the name of the group whose one Uiso the
    site shares, and ``elements``, the elements the prototype structure puts
    on the site, such as ``["Ba", "Sr"]``.
"""

from __future__ import annotations

import os
import re
import tomllib
from dataclasses import dataclass
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path

from xrdkit.config import wyckoff_multiplicity

__all__ = [
    "CELL_PARAMETERS",
    "CRYSTAL_SYSTEMS",
    "Site",
    "StructureEntry",
    "list_entries",
    "load_entry",
]

# The cell parameters each crystal system leaves free; hexagonal and trigonal
# in the hexagonal setting.
CELL_PARAMETERS = {
    "cubic": ("a",),
    "tetragonal": ("a", "c"),
    "orthorhombic": ("a", "b", "c"),
    "hexagonal": ("a", "c"),
    "trigonal": ("a", "c"),
    "monoclinic": ("a", "b", "c", "beta"),
    "triclinic": ("a", "b", "c", "alpha", "beta", "gamma"),
}
CRYSTAL_SYSTEMS = tuple(CELL_PARAMETERS)

ENTRY_REQUIRED = (
    "name",
    "family",
    "crystal_system",
    "space_group",
    "cell_parameters",
    "z",
    "reference",
)
ENTRY_OPTIONAL = ("setting", "polar_axis", "origin_site")
SITE_REQUIRED = ("label", "kind", "wyckoff")
SITE_OPTIONAL = ("free", "uiso_group", "elements")
AXES = ("a", "b", "c")
COORDINATES = ("x", "y", "z")
ELEMENT = re.compile(r"[A-Z][a-z]?")


@dataclass(frozen=True)
class Site:
    """A site of a structure entry: its label, kind and Wyckoff position,
    the coordinates it leaves free, the Uiso group it is in, if any, and the
    elements the prototype puts on it."""

    label: str
    kind: str
    wyckoff: str
    free: tuple[str, ...] = ()
    uiso_group: str | None = None
    elements: tuple[str, ...] = ()


@dataclass(frozen=True)
class StructureEntry:
    """An entry of the structure library, as :func:`load_entry` reads it."""

    name: str
    family: str
    crystal_system: str
    space_group: str
    cell_parameters: tuple[str, ...]
    z: int
    reference: str
    sites: tuple[Site, ...]
    setting: str = ""
    polar_axis: str | None = None
    origin_site: str | None = None


def _base(root: Traversable | str | os.PathLike | None) -> Traversable:
    if root is None:
        return files("xrdkit") / "structures"
    if isinstance(root, (str, os.PathLike)):
        return Path(root)
    return root


def list_entries(*, root: Traversable | str | os.PathLike | None = None) -> list[str]:
    """The names of the entries in the library, sorted.

    ``root`` is the folder to look in instead of the library shipped with
    xrdkit, laid out the same way.
    """
    base = _base(root)
    if not base.is_dir():
        return []
    names = [
        f"{family.name}/{item.name.removesuffix('.toml')}"
        for family in base.iterdir()
        if family.is_dir()
        for item in family.iterdir()
        if item.is_file() and item.name.endswith(".toml")
    ]
    return sorted(names)


def _fail(name: str, field: str, message: str) -> ValueError:
    return ValueError(f"structure entry {name!r}: {field}: {message}")


def _keys(table: object, name: str, field: str, required: tuple, optional: tuple):
    if not isinstance(table, dict):
        raise _fail(name, field, f"must be a table, not {type(table).__name__}")
    missing = [key for key in required if key not in table]
    if missing:
        raise _fail(name, field, "missing key " + ", ".join(map(repr, missing)))
    unknown = sorted(set(table) - set(required) - set(optional))
    if unknown:
        raise _fail(
            name,
            field,
            "unknown key "
            + ", ".join(map(repr, unknown))
            + "; the keys are "
            + ", ".join(required + optional),
        )
    return table


def _string(value: object, name: str, field: str, empty: bool = False) -> str:
    if not isinstance(value, str) or not (empty or value.strip()):
        raise _fail(
            name, field, "must be text" if empty else "must be a non-empty string"
        )
    return value


def _site(table: object, name: str, field: str) -> Site:
    _keys(table, name, field, SITE_REQUIRED, SITE_OPTIONAL)
    label = _string(table["label"], name, f"{field}.label")
    kind = _string(table["kind"], name, f"{field}.kind")
    wyckoff = _string(table["wyckoff"], name, f"{field}.wyckoff")
    try:
        wyckoff_multiplicity(wyckoff)
    except ValueError as error:
        raise _fail(name, f"{field}.wyckoff", str(error)) from None
    free = table.get("free", [])
    if not isinstance(free, list):
        raise _fail(name, f"{field}.free", "must be a list of coordinates")
    for coordinate in free:
        if coordinate not in COORDINATES:
            raise _fail(
                name,
                f"{field}.free",
                f"must be drawn from x, y and z, not {coordinate!r}",
            )
    if len(set(free)) != len(free):
        raise _fail(name, f"{field}.free", f"names a coordinate twice: {free}")
    group = table.get("uiso_group")
    if group is not None:
        _string(group, name, f"{field}.uiso_group")
    elements = table.get("elements", [])
    if not isinstance(elements, list) or not all(
        isinstance(element, str) and ELEMENT.fullmatch(element) for element in elements
    ):
        raise _fail(
            name,
            f"{field}.elements",
            f"must be a list of element symbols, not {elements!r}",
        )
    if len(set(elements)) != len(elements):
        raise _fail(name, f"{field}.elements", f"names an element twice: {elements}")
    return Site(label, kind, wyckoff, tuple(free), group, tuple(elements))


def _entry(data: dict, name: str) -> StructureEntry:
    unknown = sorted(set(data) - {"entry", "sites"})
    if unknown:
        raise _fail(
            name,
            "top level",
            "unknown key "
            + ", ".join(map(repr, unknown))
            + "; the keys are entry, sites",
        )
    if not isinstance(data.get("entry"), dict):
        raise _fail(name, "entry", "must be one [entry] table")
    table = _keys(data["entry"], name, "entry", ENTRY_REQUIRED, ENTRY_OPTIONAL)
    text = {
        key: _string(table[key], name, f"entry.{key}")
        for key in ("name", "family", "space_group", "reference")
    }
    if text["name"] != name:
        raise _fail(
            name, "entry.name", f"is {text['name']!r}, not the name of its file"
        )
    setting = _string(table.get("setting", ""), name, "entry.setting", empty=True)

    system = table["crystal_system"]
    if system not in CRYSTAL_SYSTEMS:
        raise _fail(
            name,
            "entry.crystal_system",
            f"must be one of {', '.join(CRYSTAL_SYSTEMS)}, not {system!r}",
        )
    parameters = table["cell_parameters"]
    expected = CELL_PARAMETERS[system]
    if not isinstance(parameters, list) or tuple(parameters) != expected:
        raise _fail(
            name,
            "entry.cell_parameters",
            f"a {system} cell has {list(expected)}, not {parameters!r}",
        )
    z = table["z"]
    if isinstance(z, bool) or not isinstance(z, int) or z < 1:
        raise _fail(name, "entry.z", f"must be a whole number of at least 1, not {z!r}")
    polar_axis = table.get("polar_axis")
    if polar_axis is not None and polar_axis not in AXES:
        raise _fail(name, "entry.polar_axis", f"must be a, b or c, not {polar_axis!r}")

    sites = data.get("sites")
    if not isinstance(sites, list) or not sites:
        raise _fail(name, "sites", "must be at least one [[sites]] table")
    checked = tuple(
        _site(site, name, f"sites[{index}]") for index, site in enumerate(sites)
    )
    labels = [site.label for site in checked]
    twice = sorted({label for label in labels if labels.count(label) > 1})
    if twice:
        raise _fail(name, "sites", f"label {', '.join(twice)} is on more than one site")
    origin_site = table.get("origin_site")
    if origin_site is not None and origin_site not in labels:
        raise _fail(
            name,
            "entry.origin_site",
            f"no site labelled {origin_site!r}; the sites are {', '.join(labels)}",
        )

    return StructureEntry(
        name=name,
        family=text["family"],
        crystal_system=system,
        space_group=text["space_group"],
        cell_parameters=expected,
        z=z,
        reference=text["reference"],
        sites=checked,
        setting=setting,
        polar_axis=polar_axis,
        origin_site=origin_site,
    )


def load_entry(
    name: str, *, root: Traversable | str | os.PathLike | None = None
) -> StructureEntry:
    """The library entry ``name``, such as ``"ttb/P4bm"``, read with tomllib
    and checked (see the module notes for its layout).

    ``root`` is the folder to read it from instead of the library shipped
    with xrdkit, laid out the same way.

    Raises
    ------
    ValueError
        If there is no entry of that name, naming those there are, or if
        the entry is not valid TOML or does not check out, naming the entry
        and the field.
    """
    available = list_entries(root=root)
    if name not in available:
        raise ValueError(
            f"no structure entry {name!r}; the entries are "
            + (", ".join(available) or "none")
        )
    family, stem = name.split("/")
    text = (_base(root) / family / f"{stem}.toml").read_text(encoding="utf-8")
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"structure entry {name!r}: not valid TOML: {error}") from None
    return _entry(data, name)
