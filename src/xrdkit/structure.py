"""Interatomic distances in a crystal structure.

:func:`interatomic_distances` finds every distance from chosen atoms to
chosen others, across all symmetry operators and neighbouring cells, in a
cell of any symmetry, with esds from those of the fractional coordinates.
The atoms and operators are in the form the GSAS-II driver reports them:
each atom a mapping with ``label``, ``xyz`` and optionally ``xyz_esd``, each
operator a mapping with a 3 x 3 ``rotation`` and a ``translation``, both
acting on fractional coordinates. :func:`bond_lengths` takes a whole phase
as the driver reports it and gives every cation to anion distance, site by
site.

:func:`site_setup` finds the sites of a structure table of a settings file
(see :mod:`xrdkit.config`) among a phase's atoms, and
:func:`composition_edits` gives the atom edits that put a nominal
composition on them by the table's rule.
"""

from __future__ import annotations

import itertools
import re
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from xrdkit.config import wyckoff_multiplicity

__all__ = [
    "Distance",
    "bond_lengths",
    "cell_contents",
    "composition_edits",
    "interatomic_distances",
    "metric_tensor",
    "site_setup",
]

# Images of an atom closer than this, in angstroms, are the same image
# reached by two operators.
SAME_IMAGE = 1.0e-4

# Atoms closer than this in every fractional coordinate share a site, and
# distances from one site to another equal to within this many angstroms
# are the same bond.
SAME_SITE = 1.0e-4
SAME_DISTANCE = 1.0e-4

# A GSAS-II cell's lengths and angles, in its keys.
CELL_KEYS = (
    "length_a",
    "length_b",
    "length_c",
    "angle_alpha",
    "angle_beta",
    "angle_gamma",
)


@dataclass(frozen=True)
class Distance:
    """A distance from atom ``centre`` to an image of atom ``target``.

    ``operator`` is the index of the symmetry operator that made the image
    and ``translation`` the whole cell shift added to it. ``esd`` is
    propagated from the coordinate esds of the two atoms alone, as if
    uncorrelated, with the cell taken as exact; None if neither atom has
    any.
    """

    centre: str
    target: str
    distance: float
    esd: float | None
    operator: int
    translation: tuple[int, int, int]


def metric_tensor(cell: Sequence[float]) -> np.ndarray:
    """The real space metric tensor G of ``cell``, ``(a, b, c, alpha, beta,
    gamma)`` in angstroms and degrees, so that a fractional vector v has
    length sqrt(v . G v)."""
    a, b, c = (float(value) for value in cell[:3])
    alpha, beta, gamma = np.radians([float(value) for value in cell[3:6]])
    return np.array(
        [
            [a * a, a * b * np.cos(gamma), a * c * np.cos(beta)],
            [a * b * np.cos(gamma), b * b, b * c * np.cos(alpha)],
            [a * c * np.cos(beta), b * c * np.cos(alpha), c * c],
        ]
    )


def _esd(atom: Mapping) -> np.ndarray:
    values = atom.get("xyz_esd") or [None, None, None]
    return np.array([0.0 if value is None else float(value) for value in values])


def interatomic_distances(
    cell: Sequence[float],
    atoms: Sequence[Mapping],
    operators: Sequence[Mapping],
    centres: Iterable[str],
    targets: Iterable[str],
    dmax: float,
    dmin: float = 0.5,
) -> list[Distance]:
    """Every distance from each atom in ``centres`` to images of the atoms in
    ``targets`` between ``dmin`` and ``dmax`` angstroms, shortest first for
    each centre, centres in the order given.

    The images are those of every operator, shifted into the cells around
    the centre; the same image reached by two operators, as on a special
    position, is counted once.

    The esd of a distance d = |x_t' - x_c|, with x_t' = R x_t + t, is
    propagated from the esds of the fractional coordinates of both atoms,
    through the derivatives -G(x_t' - x_c)/d for the centre and
    R^T G(x_t' - x_c)/d for the target; a distance from an atom to an image
    of itself moves with (R^T - I) G(x_t' - x_c)/d. Correlations between
    coordinates, and the esds of the cell, are not included.

    Raises
    ------
    ValueError
        If a centre or target is not an atom's label, or ``dmax`` is not
        above ``dmin``.
    """
    if not dmax > dmin:
        raise ValueError(f"dmax must exceed dmin, got {dmax} and {dmin}")
    by_label = {str(atom["label"]): atom for atom in atoms}
    centres, targets = list(centres), list(targets)
    missing = sorted(set(centres + targets) - set(by_label))
    if missing:
        raise ValueError(f"no atoms labelled {missing}")
    metric = metric_tensor(cell)
    rotations = [np.asarray(op["rotation"], dtype=float) for op in operators]
    shifts = [np.asarray(op["translation"], dtype=float) for op in operators]
    cells = [np.array(n) for n in itertools.product((-1, 0, 1), repeat=3)]

    found = []
    for centre in centres:
        atom = by_label[centre]
        x_c, esd_c = np.asarray(atom["xyz"], dtype=float), _esd(atom)
        near = []
        for target in targets:
            other = by_label[target]
            x_t, esd_t = np.asarray(other["xyz"], dtype=float), _esd(other)
            images = []
            for index, (rotation, shift) in enumerate(zip(rotations, shifts)):
                image = rotation @ x_t + shift
                nearest = np.round(image - x_c)
                for extra in cells:
                    offset = extra - nearest
                    delta = image + offset - x_c
                    distance = float(np.sqrt(delta @ metric @ delta))
                    if not dmin <= distance <= dmax:
                        continue
                    position = x_c + delta
                    if any(
                        np.sqrt((position - seen) @ metric @ (position - seen))
                        < SAME_IMAGE
                        for seen in images
                    ):
                        continue
                    images.append(position)
                    gradient = metric @ delta / distance
                    if target == centre:
                        variance = np.sum(
                            ((rotation.T - np.eye(3)) @ gradient * esd_c) ** 2
                        )
                    else:
                        variance = np.sum((gradient * esd_c) ** 2) + np.sum(
                            (rotation.T @ gradient * esd_t) ** 2
                        )
                    has_esd = esd_c.any() or esd_t.any()
                    near.append(
                        Distance(
                            centre=centre,
                            target=target,
                            distance=distance,
                            esd=float(np.sqrt(variance)) if has_esd else None,
                            operator=index,
                            translation=tuple(int(n) for n in offset),
                        )
                    )
        found += sorted(near, key=lambda item: item.distance)
    return found


def _element(atom_type: str) -> str:
    match = re.match(r"[A-Z][a-z]?", str(atom_type))
    return match.group(0) if match else str(atom_type)


def _sites(atoms: Sequence[Mapping]) -> list[list[int]]:
    """The indices of ``atoms`` grouped by the site they share."""
    groups: list[list[int]] = []
    for index, atom in enumerate(atoms):
        for group in groups:
            delta = np.asarray(atom["xyz"], float) - np.asarray(
                atoms[group[0]]["xyz"], float
            )
            if np.all(np.abs(delta - np.round(delta)) < SAME_SITE):
                group.append(index)
                break
        else:
            groups.append([index])
    return groups


def bond_lengths(
    phase: Mapping,
    anions: Iterable[str] = ("O",),
    dmax: float = 3.0,
    dmin: float = 0.5,
) -> list[dict]:
    """Every distance from a cation site of ``phase`` to an anion site
    between ``dmin`` and ``dmax`` angstroms.

    ``phase`` is a phase as the GSAS-II driver reports it: ``cell``, a
    mapping with GSAS-II's ``length_a`` to ``angle_gamma`` or the six values
    themselves, ``atoms``, each with ``label``, ``type`` and ``xyz`` and
    optionally ``xyz_esd``, and ``operators``. Atoms that share a site count
    once, as the first of them; a site is an anion site if its first atom
    is of an element in ``anions``, and a cation site otherwise. Distances
    from one cation site to images of one anion site that are equal to
    within SAME_DISTANCE are the same bond, counted. A phase with no anion
    site has no such bonds: the result is empty, with a warning saying so.

    Returns a list of ``{"centre", "centre_site", "target", "target_site",
    "distance", "esd", "count"}``: the labels of the first atoms and of
    every atom on each site, joined by "/", the distance and its esd as
    :func:`interatomic_distances` gives them, and how many such bonds the
    cation site has; by cation site in the phase's order, then by distance.

    Raises
    ------
    ValueError
        If ``dmax`` is not above ``dmin``.
    """
    cell = phase["cell"]
    if isinstance(cell, Mapping):
        cell = [cell[key] for key in CELL_KEYS]
    atoms = list(phase["atoms"])
    anions = set(anions)
    site_labels, centres, targets = {}, [], []
    for group in _sites(atoms):
        first = atoms[group[0]]
        site_labels[first["label"]] = "/".join(atoms[i]["label"] for i in group)
        chosen = targets if _element(first["type"]) in anions else centres
        chosen.append(first["label"])
    if not targets:
        warnings.warn(
            f"no site of an anion, {', '.join(sorted(anions))}, among the atoms of "
            "the phase, so no cation to anion bond lengths",
            stacklevel=2,
        )
        return []
    bonds: dict[tuple, dict] = {}
    for found in interatomic_distances(
        cell, atoms, phase["operators"], centres, targets, dmax, dmin
    ):
        key = (
            found.centre,
            found.target,
            round(found.distance / SAME_DISTANCE),
        )
        if key in bonds:
            bonds[key]["count"] += 1
            continue
        bonds[key] = {
            "centre": found.centre,
            "centre_site": site_labels[found.centre],
            "target": found.target,
            "target_site": site_labels[found.target],
            "distance": found.distance,
            "esd": found.esd,
            "count": 1,
        }
    order = {label: position for position, label in enumerate(centres)}
    return sorted(
        bonds.values(), key=lambda bond: (order[bond["centre"]], bond["distance"])
    )


def cell_contents(atoms: Iterable[Mapping]) -> dict[str, float]:
    """Atoms of each element per cell, multiplicity times occupancy, from
    atoms with ``type``, ``multiplicity`` and ``occupancy``."""
    contents: dict[str, float] = {}
    for atom in atoms:
        kind = _element(atom["type"])
        contents[kind] = (
            contents.get(kind, 0.0) + atom["multiplicity"] * atom["occupancy"]
        )
    return contents


def composition_edits(
    atoms: Sequence[Mapping], composition: Mapping[str, float], structure: Mapping
) -> list[dict]:
    """Atom edits that give a cell of ``atoms`` the nominal ``composition``,
    atoms of each element per formula unit, by the rule of ``structure``, a
    structure table as :func:`xrdkit.config.load_config` returns it.

    ``atoms`` are the structure's as read from its CIF, each with
    ``label``, ``type``, ``multiplicity`` and ``occupancy``, and the cell
    holds ``structure["formula_units"]`` formula units. Every element the
    atoms hold is scaled by one factor on every site it occupies, which
    keeps the CIF's distribution of it over the sites; an atom whose
    occupancy that leaves unchanged gets no edit. Each element of
    ``structure["composition"]["added"]``, which maps it to a host element,
    is put on every atom of its host at one fraction of that atom's CIF
    occupancy, so that it goes where the host is, as an atom labelled by the
    added element and the rest of the host's label (La1 on Sr1); an added
    element whose content is 0 is not added. The edits are in the form the
    GSAS-II driver takes: ``{"label", "occupancy"}``, and ``{"label",
    "type", "copy", "occupancy"}`` for an added atom.

    Raises
    ------
    ValueError
        If the structure gives no formula units, the atoms hold an element
        the composition lacks, the composition has one the atoms lack that
        the rule does not add, an added element is on the atoms already, or
        a label for an added atom is taken.
    """
    if structure.get("formula_units") is None:
        name = structure.get("name")
        raise ValueError(
            f"structure {name} gives no formula_units (Z), which putting a "
            "composition on its sites needs; add formula_units to it"
            if name
            else "the structure gives no formula_units (Z), which putting a "
            "composition on its sites needs"
        )
    added_on = structure["composition"]["added"]
    cif = cell_contents(atoms)
    per_cell = {
        kind: structure["formula_units"] * value for kind, value in composition.items()
    }
    missing = sorted(set(cif) - set(per_cell))
    if missing:
        raise ValueError(f"the atoms hold {missing}, which the composition lacks")
    lacking = sorted(
        kind for kind in set(per_cell) - set(cif) - set(added_on) if per_cell[kind]
    )
    if lacking:
        raise ValueError(
            f"the composition has {lacking}, which the atoms lack and the rule "
            "does not add"
        )
    for kind, host in added_on.items():
        if kind in cif:
            raise ValueError(f"{kind} is to be added on {host} but is there already")
        if per_cell.get(kind) and not cif.get(host):
            raise ValueError(f"{kind} is to be added on {host}, which the atoms lack")
    labels = {atom["label"] for atom in atoms}
    edits = []
    for atom in atoms:
        kind = _element(atom["type"])
        if cif[kind]:
            occupancy = atom["occupancy"] * per_cell[kind] / cif[kind]
            if occupancy != atom["occupancy"]:
                edits.append({"label": atom["label"], "occupancy": occupancy})
        for added, host in added_on.items():
            if host != kind or not per_cell.get(added):
                continue
            label = added + atom["label"][len(kind) :]
            if label in labels:
                raise ValueError(
                    f"label {label!r} for {added} on {atom['label']} is taken"
                )
            labels.add(label)
            edits.append(
                {
                    "label": label,
                    "type": added,
                    "copy": atom["label"],
                    "occupancy": atom["occupancy"] * per_cell[added] / cif[kind],
                }
            )
    return edits


def site_setup(structure: Mapping, atoms: Sequence[Mapping]) -> dict:
    """The sites of ``structure``, a structure table as
    :func:`xrdkit.config.load_config` returns it, found among ``atoms``, a
    phase's atoms as the GSAS-II driver reports them (``label``, ``type``,
    ``xyz`` and optionally ``multiplicity``).

    Each configured site is the position of its atoms; any other atom on
    that position, such as one added by :func:`composition_edits`, is on
    the site too. Returns a dict with

    - ``sites``: in the structure's order, each ``{"name", "kind",
      "wyckoff", "multiplicity", "free", "atoms"}``, ``free`` the
      coordinates to refine from the table's ``free_coordinates`` ("all"
      for a position it does not name) and ``atoms`` those on the position,
      in the order of ``atoms``;
    - ``kinds``: the sites of each kind, in that order;
    - ``uiso_groups`` and ``group_names``: the site names of each Uiso
      group, and its name;
    - ``origin``: the site that fixes the origin, and ``origin_axis``, both
      None for a structure with no origin;
    - ``coordinates``: by kind, the coordinates to refine on each site but
      the origin's (on every site when there is no origin), ``{kind: {site
      name: free}}``;
    - ``exchange``: ``{"elements", "sites"}``, or None for a structure with
      no exchange.

    Raises
    ------
    ValueError
        If a configured atom is not among ``atoms``, a site's atoms are not
        on one position or two sites on the same one, an atom is not of the
        element given for it, a site's multiplicity is not that of its
        Wyckoff position, or an atom is on no configured site.
    """
    groups = _sites(atoms)
    group_of = {index: number for number, group in enumerate(groups) for index in group}
    index_of = {str(atom["label"]): index for index, atom in enumerate(atoms)}
    free = structure.get("free_coordinates", {})
    sites, used = [], {}
    for site in structure["sites"]:
        name = site["name"]
        absent = [label for label in site["atoms"] if label not in index_of]
        if absent:
            raise ValueError(f"site {name}: no atom labelled {absent} in the structure")
        numbers = {group_of[index_of[label]] for label in site["atoms"]}
        if len(numbers) > 1:
            raise ValueError(
                f"site {name}: atoms {list(site['atoms'])} are not on one position"
            )
        (number,) = numbers
        if number in used:
            raise ValueError(f"sites {used[number]} and {name} are on one position")
        used[number] = name
        for label, element in site["atoms"].items():
            found = _element(atoms[index_of[label]]["type"])
            if found != element:
                raise ValueError(
                    f"site {name}: atom {label} is {found}, not {element} as given"
                )
        members = [atoms[index] for index in groups[number]]
        multiplicity = wyckoff_multiplicity(site["wyckoff"])
        if "multiplicity" in members[0] and members[0]["multiplicity"] != multiplicity:
            raise ValueError(
                f"site {name}: multiplicity {members[0]['multiplicity']}, not "
                f"{multiplicity} as on {site['wyckoff']}"
            )
        sites.append(
            {
                "name": name,
                "kind": site["kind"],
                "wyckoff": site["wyckoff"],
                "multiplicity": multiplicity,
                "free": free.get(site["wyckoff"], "all"),
                "atoms": members,
            }
        )
    unplaced = [
        str(atoms[index]["label"])
        for number, group in enumerate(groups)
        if number not in used
        for index in group
    ]
    if unplaced:
        raise ValueError(f"atoms {unplaced} are on no configured site")
    by_name = {site["name"]: site for site in sites}
    kinds: dict[str, list[dict]] = {}
    for site in sites:
        kinds.setdefault(site["kind"], []).append(site)
    held = structure.get("origin")
    origin = by_name[held["site"]] if held else None
    exchange = structure.get("exchange")
    return {
        "sites": sites,
        "kinds": kinds,
        "uiso_groups": [list(group["sites"]) for group in structure["uiso_groups"]],
        "group_names": [group["name"] for group in structure["uiso_groups"]],
        "origin": origin,
        "origin_axis": held["axis"] if held else None,
        "coordinates": {
            kind: {site["name"]: site["free"] for site in found if site is not origin}
            for kind, found in kinds.items()
        },
        "exchange": (
            {
                "elements": list(exchange["elements"]),
                "sites": list(exchange["sites"]),
            }
            if exchange
            else None
        ),
    }
