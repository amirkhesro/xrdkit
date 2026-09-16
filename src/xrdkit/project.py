"""The project file, ``xrdkit.toml``: a project's instruments, structure
models and samples in one place.

The project root is the folder holding ``xrdkit.toml``, and every path in the
file is relative to it. :func:`find_project` finds the file from a folder
inside the project, and :func:`load_project` reads and checks it. The file
holds

``[project]``
    ``name``, and ``version``, the version of this layout, 1.

``[instruments.<key>]``
    ``wavelength``, one or two positive numbers in angstroms (Kα1 and Kα2,
    or Kα1 only), and ``ka2``, true or false; optionally ``radius``, the
    goniometer radius in mm, and ``instprm``, a GSAS-II instrument
    parameter file.

``[structures.<key>]``
    ``library``, the name of a shipped structure entry (see
    :func:`xrdkit.library.list_entries`), or ``cif``, a CIF file, or both, the
    entry giving the sites, their kinds and free coordinates and the CIF the
    coordinates a Rietveld refinement starts from;
    ``composition``, a formula :func:`xrdkit.density.parse_formula` reads;
    ``cell``, ``{a = ..., c = ...}`` in angstroms and degrees, required with
    ``library`` and then with exactly the entry's cell parameters, optional
    with ``cif``; and optionally ``z``, formula units per cell,
    ``exchange``, lists of the elements whose occupancies are traded,
    ``atoms`` and ``origin``.

    ``atoms`` takes the shapes :mod:`xrdkit.config` reads, with its readers.
    With ``library`` it is a table of the entry's sites by label, each the
    atoms on it, ``{A1 = {Sr1 = "Sr", La1 = "La"}, ...}``, by their labels in
    the CIF when there is one, an atom the CIF lacks being put on the site
    beside the CIF's; a site named
    holds exactly those atoms, and one not named the prototype's elements.
    When it is given, every element of the composition must be on a site, so
    an element the prototype does not carry must be placed by it; it may be
    left out when the composition holds only the prototype's elements. With
    ``cif`` alone it is the list of the CIF's sites, ``[{atoms = {Nb1 = "Nb"},
    wyckoff = "2b", kind = "B"}, ...]``, each named by its first atom.

    ``origin`` is the site that fixes the origin: a site label, one of the
    entry's sites with ``library``; ``false``, to fix none; or, with
    ``cif``, ``{site, axis}``, axis x, y or z. Left out, the entry's origin
    site applies, and a CIF structure has none. With ``cif`` and ``atoms``
    the site must be one of the atoms' sites.

``[samples.<key>]``
    ``file``, its scan; ``instrument``, an instrument key; ``structures``,
    one or more structure keys; ``form``, powder or pellet; and optionally
    ``stage``, free text, ``temperature_c``, ``archimedes``, its measured
    density in g/cm3, ``notes``, and ``refine``, a table overriding any key
    of ``[refine]`` for this sample.

``[refine]``
    Optional defaults for refining every sample, over the package's own:
    ``two_theta``, ``[low, high]`` in degrees (the scan's own range when
    left out); ``background``, ``{function, terms}``, by default
    chebyschev-1 with 6 terms; ``max_passes``, the most passes of a stage
    by mode, ``{lebail, fixed_atoms, coordinates, occupancies}``, by default
    60, 60, 100 and 100; ``unsettled``, what becomes of a stage not settled
    by then, accept or reject by mode, by default accept; and ``followed``,
    reflections to follow from mode to mode as ``[h, k, l]`` triples, by
    default none. A table given in part keeps the rest of the table below
    it; :func:`refine_settings` resolves a sample's settings.

Sample and structure keys name folders, so they hold only letters, digits,
dots, underscores and hyphens.
"""

from __future__ import annotations

import json
import math
import os
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from xrdkit.config import (
    UNSETTLED_RULES,
    ConfigError,
    read_library_atoms,
    read_sites,
)
from xrdkit.density import ATOMIC_MASSES, parse_formula
from xrdkit.library import list_entries, load_entry

__all__ = [
    "FORMS",
    "PROJECT_FILE",
    "Instrument",
    "Project",
    "Refine",
    "Sample",
    "StructureSpec",
    "find_project",
    "load_project",
    "load_project_text",
    "project_template",
    "refine_settings",
    "resolved_cell",
    "resolved_z",
    "results_dir",
    "toml_string",
]

PROJECT_FILE = "xrdkit.toml"
VERSION = 1
FORMS = ("powder", "pellet")
CELL_PARAMETERS = ("a", "b", "c", "alpha", "beta", "gamma")

TOP_LEVEL = ("project", "instruments", "structures", "samples", "refine")
PROJECT_KEYS = ("name", "version")
INSTRUMENT_REQUIRED = ("wavelength", "ka2")
INSTRUMENT_OPTIONAL = ("radius", "instprm")
STRUCTURE_REQUIRED = ("composition",)
STRUCTURE_OPTIONAL = ("library", "cif", "cell", "z", "exchange", "origin", "atoms")
SAMPLE_REQUIRED = ("file", "instrument", "structures", "form")
SAMPLE_OPTIONAL = ("stage", "temperature_c", "archimedes", "notes", "refine")
AXES = ("x", "y", "z")

# The refinement modes, in order, the keys of a refine table and the
# package's defaults for them.
REFINE_MODES = ("lebail", "fixed_atoms", "coordinates", "occupancies")
REFINE_KEYS = ("two_theta", "background", "max_passes", "unsettled", "followed")
BACKGROUND_KEYS = ("function", "terms")
DEFAULT_BACKGROUND = {"function": "chebyschev-1", "terms": 6}
DEFAULT_MAX_PASSES = {
    "lebail": 60,
    "fixed_atoms": 60,
    "coordinates": 100,
    "occupancies": 100,
}
DEFAULT_UNSETTLED = dict.fromkeys(REFINE_MODES, "accept")

# A key that names a folder.
FOLDER_KEY = re.compile(r"[A-Za-z0-9._-]+")


@dataclass(frozen=True)
class Instrument:
    """An instrument: its wavelengths in angstroms, Kα1 and Kα2 or Kα1
    alone, whether the scans hold Kα2, and optionally its goniometer radius
    in mm and its GSAS-II instrument parameter file."""

    key: str
    wavelength: tuple[float, ...]
    ka2: bool
    radius: float | None = None
    instprm: Path | None = None


@dataclass(frozen=True)
class StructureSpec:
    """A structure model: a library entry or a CIF, the composition put on
    it, and the cell, formula units, exchanges, origin and atoms given.

    The origin is one of three: the entry's default, ``origin`` None and
    ``origin_fixed`` true; a site given, ``origin`` its label, with
    ``origin_axis`` when the file gives one; or none fixed, ``origin_fixed``
    false. ``atoms`` is None when not given, else its sites as
    :func:`xrdkit.config.read_sites` gives them; those of a library
    structure carry their entry ``label`` too, and only the sites named are
    there."""

    key: str
    composition: str
    library: str | None = None
    cif: Path | None = None
    cell: dict[str, float] | None = None
    z: int | None = None
    exchange: tuple[tuple[str, ...], ...] = ()
    origin: str | None = None
    origin_axis: str | None = None
    origin_fixed: bool = True
    atoms: tuple[dict, ...] | None = None


@dataclass(frozen=True)
class Sample:
    """A sample: its scan, instrument key, structure keys, form, and
    optionally its stage, temperature, measured density, notes and the
    ``[refine]`` keys it overrides, checked, a table given in part."""

    key: str
    file: Path
    instrument: str
    structures: tuple[str, ...]
    form: str
    stage: str | None = None
    temperature_c: float | None = None
    archimedes: float | None = None
    notes: str = ""
    refine: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class Refine:
    """Refinement settings: the range in degrees, None for the scan's own;
    the background function and terms; the most passes of a stage and the
    rule for one not settled, by mode; and the reflections followed, as
    (h, k, l)."""

    two_theta: tuple[float, float] | None = None
    background: dict[str, object] = field(
        default_factory=lambda: dict(DEFAULT_BACKGROUND)
    )
    max_passes: dict[str, int] = field(default_factory=lambda: dict(DEFAULT_MAX_PASSES))
    unsettled: dict[str, str] = field(default_factory=lambda: dict(DEFAULT_UNSETTLED))
    followed: tuple[tuple[int, int, int], ...] = ()


@dataclass(frozen=True)
class Project:
    """A project as :func:`load_project` reads it, its paths made absolute
    under ``root``, the folder holding its ``xrdkit.toml``."""

    root: Path
    name: str
    version: int
    instruments: dict[str, Instrument]
    structures: dict[str, StructureSpec]
    samples: dict[str, Sample]
    refine: Refine = field(default_factory=Refine)


def find_project(start: str | os.PathLike | None = None) -> Path:
    """The ``xrdkit.toml`` of the first folder, from ``start`` (the current
    folder by default) upwards, that holds one.

    Raises
    ------
    FileNotFoundError
        If neither that folder nor any above it holds one.
    """
    folder = Path.cwd() if start is None else Path(start)
    folder = folder.resolve()
    if folder.is_file():
        folder = folder.parent
    for candidate in (folder, *folder.parents):
        path = candidate / PROJECT_FILE
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"no {PROJECT_FILE} in {folder} or any folder above it; "
        "run xrdkit init in the project folder to make one"
    )


class _Checker:
    """Checks for one project file, each failure naming the file, the table
    and the key."""

    def __init__(self, source: Path, root: Path) -> None:
        self.source = source
        self.root = root

    def fail(self, where: str, message: str) -> ValueError:
        return ValueError(f"{self.source}: {where}: {message}")

    def table(self, value: object, where: str) -> Mapping:
        if not isinstance(value, Mapping):
            raise self.fail(where, f"must be a table, not {type(value).__name__}")
        return value

    def keys(
        self, value: object, where: str, required: tuple, optional: tuple = ()
    ) -> Mapping:
        table = self.table(value, where)
        unknown = sorted(set(table) - set(required) - set(optional))
        if unknown:
            raise self.fail(
                where,
                "unknown key "
                + ", ".join(map(repr, unknown))
                + "; the keys are "
                + ", ".join(required + optional),
            )
        missing = [key for key in required if key not in table]
        if missing:
            raise self.fail(where, "missing key " + ", ".join(map(repr, missing)))
        return table

    def string(self, value: object, where: str, empty: bool = False) -> str:
        if not isinstance(value, str) or not (empty or value.strip()):
            raise self.fail(
                where, "must be text" if empty else "must be a non-empty string"
            )
        return value

    def number(self, value: object, where: str, positive: bool = False) -> float:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise self.fail(where, f"must be a number, not {value!r}")
        if positive and value <= 0:
            raise self.fail(where, f"must be positive, not {value!r}")
        return float(value)

    def integer(self, value: object, where: str) -> int:
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise self.fail(where, f"must be a positive whole number, not {value!r}")
        return value

    def config(self, read, *args):
        """``read(*args)``, a reader of :mod:`xrdkit.config`, its error made
        one of this file's."""
        try:
            return read(*args)
        except ConfigError as error:
            raise ValueError(f"{self.source}: {error}") from None

    def path(self, value: object, where: str, exists: bool) -> Path:
        path = self.root / self.string(value, where)
        if exists and not path.is_file():
            raise self.fail(where, f"no such file: {path}")
        return path

    def folder_key(self, key: str, where: str) -> str:
        if not FOLDER_KEY.fullmatch(key) or key in (".", ".."):
            raise self.fail(
                where,
                f"key {key!r} must be letters, digits, dots, underscores and "
                "hyphens only, since it names a folder",
            )
        return key


def _instrument(check: _Checker, key: str, value: object) -> Instrument:
    where = f"instruments.{key}"
    table = check.keys(value, where, INSTRUMENT_REQUIRED, INSTRUMENT_OPTIONAL)
    wavelength = table["wavelength"]
    if not isinstance(wavelength, list) or len(wavelength) not in (1, 2):
        raise check.fail(
            f"{where}.wavelength",
            f"must be one or two positive numbers in a list, not {wavelength!r}",
        )
    ka2 = table["ka2"]
    if not isinstance(ka2, bool):
        raise check.fail(f"{where}.ka2", f"must be true or false, not {ka2!r}")
    if len(wavelength) != (2 if ka2 else 1):
        raise check.fail(
            where,
            "ka2 = true needs two wavelengths, Kα1 and Kα2"
            if ka2
            else "ka2 = false needs one wavelength, Kα1",
        )
    return Instrument(
        key=key,
        wavelength=tuple(
            check.number(item, f"{where}.wavelength[{index}]", positive=True)
            for index, item in enumerate(wavelength)
        ),
        ka2=ka2,
        radius=(
            check.number(table["radius"], f"{where}.radius", positive=True)
            if "radius" in table
            else None
        ),
        instprm=(
            check.path(table["instprm"], f"{where}.instprm", exists=True)
            if "instprm" in table
            else None
        ),
    )


def _structure(check: _Checker, key: str, value: object) -> StructureSpec:
    where = f"structures.{key}"
    check.folder_key(key, where)
    table = check.keys(value, where, STRUCTURE_REQUIRED, STRUCTURE_OPTIONAL)
    if "library" not in table and "cif" not in table:
        raise check.fail(where, "must give library or cif, or both")

    composition = check.string(table["composition"], f"{where}.composition")
    try:
        formula = parse_formula(composition)
    except ValueError as error:
        raise check.fail(f"{where}.composition", str(error)) from None

    library = cif = entry = None
    if "library" in table:
        library = check.string(table["library"], f"{where}.library")
        available = list_entries()
        if library not in available:
            raise check.fail(
                f"{where}.library",
                f"no library entry {library!r}; the entries are "
                + ", ".join(available),
            )
        entry = load_entry(library)
    if "cif" in table:
        cif = check.path(table["cif"], f"{where}.cif", exists=True)

    cell = None
    if "cell" in table:
        allowed = entry.cell_parameters if entry else CELL_PARAMETERS
        given = check.table(table["cell"], f"{where}.cell")
        if entry and set(given) != set(allowed):
            raise check.fail(
                f"{where}.cell",
                f"must give exactly {', '.join(allowed)}, the cell parameters of "
                f"{library}, not {', '.join(given) or 'none'}",
            )
        check.keys(given, f"{where}.cell", (), allowed)
        cell = {
            name: check.number(number, f"{where}.cell.{name}", positive=True)
            for name, number in given.items()
        }
    elif entry:
        raise check.fail(
            f"{where}.cell",
            f"is required with library; give {', '.join(entry.cell_parameters)}",
        )

    z = check.integer(table["z"], f"{where}.z") if "z" in table else None

    exchange = ()
    if "exchange" in table:
        groups = table["exchange"]
        if not isinstance(groups, list) or not all(
            isinstance(group, list) for group in groups
        ):
            raise check.fail(
                f"{where}.exchange",
                f"must be a list of lists of element symbols, not {groups!r}",
            )
        for index, group in enumerate(groups):
            at = f"{where}.exchange[{index}]"
            if len(group) < 2 or len(set(map(str, group))) != len(group):
                raise check.fail(at, f"must name two or more elements, not {group!r}")
            for element in group:
                if element not in ATOMIC_MASSES:
                    raise check.fail(at, f"{element!r} is not an element symbol")
        exchange = tuple(tuple(group) for group in groups)

    atoms = None
    if "atoms" in table:
        at = f"{where}.atoms"
        if entry:
            given, _ = check.config(
                read_library_atoms, table["atoms"], at, entry, False
            )
            _check_placed(check, at, entry, given, formula)
            atoms = tuple(given)
        else:
            atoms = tuple(check.config(read_sites, table["atoms"], at))

    origin, origin_axis, origin_fixed = _origin(check, where, table, entry, atoms)

    return StructureSpec(
        key=key,
        composition=composition,
        library=library,
        cif=cif,
        cell=cell,
        z=z,
        exchange=exchange,
        origin=origin,
        origin_axis=origin_axis,
        origin_fixed=origin_fixed,
        atoms=atoms,
    )


def _plural(elements: list[str], one: str, many: str) -> str:
    return one if len(elements) == 1 else many


def _check_placed(
    check: _Checker,
    where: str,
    entry,
    given: list[dict],
    formula: Mapping[str, float],
) -> None:
    """Check that every element of ``formula`` is on a site of ``entry`` once
    the sites ``given`` hold the atoms given for them and the rest the
    prototype's elements."""
    by_label = {site["label"]: site for site in given}
    held = set()
    for site in entry.sites:
        if site.label in by_label:
            held.update(by_label[site.label]["atoms"].values())
        else:
            held.update(site.elements)
    prototype = {element for site in entry.sites for element in site.elements}
    labels = ", ".join(site.label for site in entry.sites)
    unplaced = [element for element in formula if element not in held]
    new = [element for element in unplaced if element not in prototype]
    if new:
        raise check.fail(
            where,
            f"{', '.join(new)} of the composition {_plural(new, 'is', 'are')} not "
            f"in the {entry.name} prototype; place {_plural(new, 'it', 'them')} "
            f"on one of its sites {labels}",
        )
    if unplaced:
        raise check.fail(
            where,
            f"{', '.join(unplaced)} of the composition "
            f"{_plural(unplaced, 'is', 'are')} on none of the sites once atoms "
            f"replaces the prototype's; keep {_plural(unplaced, 'it', 'them')} on "
            f"one of the sites {labels}",
        )


def _origin(
    check: _Checker,
    where: str,
    table: Mapping,
    entry,
    atoms: tuple[dict, ...] | None,
) -> tuple[str | None, str | None, bool]:
    """The origin site given, its axis, and whether an origin is fixed."""
    if "origin" not in table:
        return None, None, True
    value = table["origin"]
    at = f"{where}.origin"
    if value is False:
        return None, None, False
    axis = None
    if isinstance(value, str):
        site = check.string(value, at)
    elif isinstance(value, Mapping) and entry is None:
        given = check.keys(value, at, ("site", "axis"))
        site = check.string(given["site"], f"{at}.site")
        axis = given["axis"]
        if axis not in AXES:
            raise check.fail(f"{at}.axis", f"must be x, y or z, not {axis!r}")
        at = f"{at}.site"
    elif isinstance(value, Mapping):
        raise check.fail(
            at,
            "the {site, axis} form goes with cif; with library give one of the "
            f"sites of {entry.name}: {', '.join(s.label for s in entry.sites)}",
        )
    else:
        raise check.fail(
            at,
            f"must be a site label, false, or {{site, axis}} with cif, not {value!r}",
        )
    if entry is not None:
        labels = [s.label for s in entry.sites]
        if site not in labels:
            raise check.fail(
                at,
                f"no site {site!r} in {entry.name}; the sites are {', '.join(labels)}",
            )
    elif atoms is not None:
        names = [s["name"] for s in atoms]
        if site not in names:
            raise check.fail(
                at,
                f"no site {site!r} in atoms; sites are named by their first atom: "
                f"{', '.join(names)}",
            )
    return site, axis, True


def _sample(
    check: _Checker,
    key: str,
    value: object,
    instruments: Mapping[str, Instrument],
    structures: Mapping[str, StructureSpec],
) -> Sample:
    where = f"samples.{key}"
    check.folder_key(key, where)
    table = check.keys(value, where, SAMPLE_REQUIRED, SAMPLE_OPTIONAL)
    file = check.path(table["file"], f"{where}.file", exists=True)

    instrument = check.string(table["instrument"], f"{where}.instrument")
    if instrument not in instruments:
        raise check.fail(
            f"{where}.instrument",
            f"no instrument {instrument!r} under instruments; there are "
            + (", ".join(instruments) or "none"),
        )

    named = table["structures"]
    if not isinstance(named, list) or not named:
        raise check.fail(f"{where}.structures", "must list one or more structure keys")
    for index, name in enumerate(named):
        check.string(name, f"{where}.structures[{index}]")
        if name not in structures:
            raise check.fail(
                f"{where}.structures[{index}]",
                f"no structure {name!r} under structures; there are "
                + (", ".join(structures) or "none"),
            )
    if len(set(named)) != len(named):
        raise check.fail(f"{where}.structures", f"names a structure twice: {named}")

    form = table["form"]
    if form not in FORMS:
        raise check.fail(f"{where}.form", f"must be powder or pellet, not {form!r}")

    return Sample(
        key=key,
        file=file,
        instrument=instrument,
        structures=tuple(named),
        form=form,
        stage=(
            check.string(table["stage"], f"{where}.stage", empty=True)
            if "stage" in table
            else None
        ),
        temperature_c=(
            check.number(table["temperature_c"], f"{where}.temperature_c")
            if "temperature_c" in table
            else None
        ),
        archimedes=(
            check.number(table["archimedes"], f"{where}.archimedes", positive=True)
            if "archimedes" in table
            else None
        ),
        notes=check.string(table.get("notes", ""), f"{where}.notes", empty=True),
        refine=(
            _refine(check, table["refine"], f"{where}.refine")
            if "refine" in table
            else {}
        ),
    )


def _refine(check: _Checker, value: object, where: str) -> dict[str, object]:
    """The keys a ``refine`` table gives, checked; a table given in part."""
    table = check.keys(value, where, (), REFINE_KEYS)
    found: dict[str, object] = {}

    if "two_theta" in table:
        at = f"{where}.two_theta"
        pair = table["two_theta"]
        if not isinstance(pair, list) or len(pair) != 2:
            raise check.fail(at, f"must be [low, high] in degrees, not {pair!r}")
        low, high = (
            check.number(item, f"{at}[{index}]") for index, item in enumerate(pair)
        )
        if not 0.0 <= low < high <= 180.0:
            raise check.fail(
                at, f"must rise from low to high within 0 to 180 degrees, not {pair!r}"
            )
        found["two_theta"] = (low, high)

    if "background" in table:
        at = f"{where}.background"
        given = check.keys(table["background"], at, (), BACKGROUND_KEYS)
        background: dict[str, object] = {}
        if "function" in given:
            background["function"] = check.string(given["function"], f"{at}.function")
        if "terms" in given:
            background["terms"] = check.integer(given["terms"], f"{at}.terms")
        found["background"] = background

    if "max_passes" in table:
        at = f"{where}.max_passes"
        given = check.keys(table["max_passes"], at, (), REFINE_MODES)
        found["max_passes"] = {
            mode: check.integer(passes, f"{at}.{mode}")
            for mode, passes in given.items()
        }

    if "unsettled" in table:
        at = f"{where}.unsettled"
        given = check.keys(table["unsettled"], at, (), REFINE_MODES)
        for mode, rule in given.items():
            if rule not in UNSETTLED_RULES:
                raise check.fail(
                    f"{at}.{mode}", f"must be accept or reject, not {rule!r}"
                )
        found["unsettled"] = dict(given)

    if "followed" in table:
        at = f"{where}.followed"
        reflections = table["followed"]
        if not isinstance(reflections, list):
            raise check.fail(
                at, f"must be a list of [h, k, l] triples, not {reflections!r}"
            )
        followed = []
        for index, hkl in enumerate(reflections):
            if (
                not isinstance(hkl, list)
                or len(hkl) != 3
                or any(isinstance(i, bool) or not isinstance(i, int) for i in hkl)
                or not any(hkl)
            ):
                raise check.fail(
                    f"{at}[{index}]",
                    f"must be [h, k, l], three whole numbers not all 0, not {hkl!r}",
                )
            followed.append(tuple(hkl))
        found["followed"] = tuple(followed)
    return found


def _resolved(base: Refine, overrides: Mapping[str, object]) -> Refine:
    """``base`` with ``overrides`` over it, a table given in part keeping the
    rest of ``base``'s."""
    return Refine(
        two_theta=overrides.get("two_theta", base.two_theta),
        background={**base.background, **overrides.get("background", {})},
        max_passes={**base.max_passes, **overrides.get("max_passes", {})},
        unsettled={**base.unsettled, **overrides.get("unsettled", {})},
        followed=overrides.get("followed", base.followed),
    )


def refine_settings(project: Project, sample: Sample | str) -> Refine:
    """The refinement settings of ``sample``, a sample or its key: the
    package defaults, the project's ``[refine]`` over them, and the sample's
    own ``refine`` over that.

    Raises
    ------
    KeyError
        If ``project`` has no sample of that key.
    """
    key = sample.key if isinstance(sample, Sample) else sample
    if key not in project.samples:
        raise KeyError(
            f"no sample {key!r} in the project; there are "
            + (", ".join(project.samples) or "none")
        )
    return _resolved(project.refine, project.samples[key].refine)


def load_project(path: str | os.PathLike | None = None) -> Project:
    """Read the project file at ``path``, an ``xrdkit.toml`` or the folder
    holding one (by default the one :func:`find_project` finds), and check
    it (see the module notes for its layout).

    Raises
    ------
    FileNotFoundError
        If there is no project file there.
    ValueError
        If the file is not valid TOML or does not check out, naming the file,
        the table and the key.
    """
    if path is None:
        source = find_project()
    else:
        source = Path(path)
        if source.is_dir():
            source = source / PROJECT_FILE
    if not source.is_file():
        raise FileNotFoundError(
            f"no project file {source}; run xrdkit init to make one"
        )
    return load_project_text(source.read_text(encoding="utf-8"), source)


def load_project_text(text: str, path: str | os.PathLike) -> Project:
    """Check ``text`` as the project file at ``path`` would be checked, its
    paths taken relative to the folder of ``path``, without reading or
    writing that file; so that a change can be checked before it is made.

    Raises
    ------
    ValueError
        If the text is not valid TOML or does not check out, naming ``path``,
        the table and the key.
    """
    source = Path(path).resolve()
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"{source}: not valid TOML: {error}") from None

    check = _Checker(source, source.parent)
    check.keys(data, "top level", ("project",), TOP_LEVEL[1:])
    project = check.keys(data["project"], "project", PROJECT_KEYS)
    name = check.string(project["name"], "project.name")
    version = project["version"]
    if isinstance(version, bool) or version != VERSION or not isinstance(version, int):
        raise check.fail("project.version", f"must be {VERSION}, not {version!r}")

    instruments = {
        key: _instrument(check, key, value)
        for key, value in check.table(
            data.get("instruments", {}), "instruments"
        ).items()
    }
    structures = {
        key: _structure(check, key, value)
        for key, value in check.table(data.get("structures", {}), "structures").items()
    }
    samples = {
        key: _sample(check, key, value, instruments, structures)
        for key, value in check.table(data.get("samples", {}), "samples").items()
    }
    overrides = _refine(check, data["refine"], "refine") if "refine" in data else {}
    return Project(
        root=source.parent,
        name=name,
        version=version,
        instruments=instruments,
        structures=structures,
        samples=samples,
        refine=_resolved(Refine(), overrides),
    )


def resolved_z(spec: StructureSpec) -> int:
    """The formula units per cell of ``spec``: its own ``z``, or else that of
    its library entry.

    Raises
    ------
    ValueError
        If ``spec`` uses a CIF and gives no ``z``.
    """
    if spec.z is not None:
        return spec.z
    if spec.library is None:
        raise ValueError(
            f"structures.{spec.key}: uses a CIF and gives no z; add z to the structure"
        )
    return load_entry(spec.library).z


# The cell parameters a crystal system leaves out, and their values: equal to
# another parameter, or a fixed angle.
_IMPLIED = {
    "cubic": {"b": "a", "c": "a", "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
    "tetragonal": {"b": "a", "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
    "orthorhombic": {"alpha": 90.0, "beta": 90.0, "gamma": 90.0},
    "hexagonal": {"b": "a", "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
    "trigonal": {"b": "a", "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
    "monoclinic": {"alpha": 90.0, "gamma": 90.0},
    "triclinic": {},
}


def resolved_cell(spec: StructureSpec) -> dict[str, float]:
    """All six cell parameters of ``spec``, a, b, c in angstroms and alpha,
    beta, gamma in degrees, with those its library entry's crystal system
    leaves out filled in (b = a and the angles of a tetragonal cell, say).

    Raises
    ------
    ValueError
        If ``spec`` gives no cell, or uses a CIF and does not give all six,
        since a CIF structure's crystal system is not known here.
    """
    if spec.cell is None:
        raise ValueError(f"structures.{spec.key}: gives no cell")
    cell = dict(spec.cell)
    if spec.library is not None:
        for name, value in _IMPLIED[load_entry(spec.library).crystal_system].items():
            cell[name] = cell[value] if isinstance(value, str) else value
    missing = [name for name in CELL_PARAMETERS if name not in cell]
    if missing:
        raise ValueError(
            f"structures.{spec.key}: uses a CIF, so its cell needs all of "
            f"{', '.join(CELL_PARAMETERS)}; it lacks {', '.join(missing)}"
        )
    return {name: cell[name] for name in CELL_PARAMETERS}


def results_dir(project: Project, command: str, sample: Sample | str) -> Path:
    """The folder a command's results for ``sample`` go in,
    ``<root>/results/<command>/<sample key>``; it is not created."""
    key = sample.key if isinstance(sample, Sample) else sample
    return project.root / "results" / command / key


TEMPLATE = """# xrdkit project file. Every path in it is relative to the folder holding it.
# xrdkit commands find it from this folder or any folder inside it.
# Uncomment and edit the examples below, one table per instrument, structure
# and sample.

[project]
name = {name}
version = 1

# An instrument, named by its key. wavelength is in angstroms: Kα1 and Kα2,
# or Kα1 alone. ka2 says whether the scans hold Kα2. radius, the goniometer
# radius in mm, and instprm, a GSAS-II instrument parameter file, are
# optional.
#
# [instruments.diffractometer]
# wavelength = [1.540598, 1.544426]
# ka2 = true
# radius = 240.0
# instprm = "data/standards/diffractometer.instprm"

# A structure model, named by its key, which names its results folders too.
# Give library, a structure library entry such as ttb/P4bm or perovskite/P4mm
# (xrdkit.list_entries() lists them), or cif, a CIF file, or both: a Rietveld
# refinement needs the CIF's coordinates, Le Bail only the entry.
# composition is the formula put on it. cell is required with library, with
# exactly the entry's cell parameters, and optional with cif. z, exchange
# (groups of elements whose occupancies are traded), origin (the label of the
# site that fixes the origin, or false for none) and atoms are optional. With
# library, atoms places the elements the entry's prototype lacks on its
# sites, A1 = {{ Sr1 = "Sr", La1 = "La" }} say, and must place every such
# element of the composition.
#
# [structures.phase1]
# library = "ttb/P4bm"
# composition = "Sr0.5Ba0.5Nb2O6"
# cell = {{ a = 12.45, c = 3.94 }}
# z = 5
# exchange = [["Sr", "Ba"]]
# origin = "B1"

# A sample, named by its key, which names its results folders too. file is
# its scan, instrument an instruments key and structures one or more
# structures keys. form is powder or pellet. stage, temperature_c, archimedes
# (the measured density in g/cm3) and notes are optional.
#
# [samples.sample1]
# file = "data/raw/sample.xrdml"
# instrument = "diffractometer"
# structures = ["phase1"]
# stage = "calcined"
# form = "powder"
# temperature_c = 1300
# archimedes = 5.21
# notes = ""

# Refinement defaults for every sample, all optional: two_theta, the range in
# degrees, the scan's own when left out; background; max_passes, the most
# passes of a stage, and unsettled, accept or reject for a stage still moving
# after them, by mode; and followed, reflections to follow as [h, k, l]. A
# sample's own refine table overrides any of them.
#
# [refine]
# two_theta = [15.0, 100.0]
# background = {{ function = "chebyschev-1", terms = 6 }}
# max_passes = {{ lebail = 60, fixed_atoms = 60, coordinates = 100, occupancies = 100 }}
# unsettled = {{ lebail = "accept", fixed_atoms = "accept", coordinates = "accept", occupancies = "accept" }}
# followed = [[2, 1, 1], [4, 0, 0]]
#
# [samples.sample1.refine]
# two_theta = [17.0, 100.0]
# background = {{ terms = 8 }}
"""


def toml_string(value: str) -> str:
    """``value`` as a TOML basic string, quoted and escaped."""
    return json.dumps(value, ensure_ascii=False)


def project_template(name: str) -> str:
    """The text of a new project file named ``name``, as ``xrdkit init``
    writes it: ``[project]`` filled in and a commented example of each other
    table."""
    return TEMPLATE.format(name=toml_string(name))
