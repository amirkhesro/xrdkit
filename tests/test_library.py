"""Tests for xrdkit.library."""

import re
from fractions import Fraction
from importlib.resources import files

import pytest

import xrdkit
from xrdkit.library import (
    DEFAULT_BOND_LIMITS,
    Site,
    StructureEntry,
    bond_limits,
    list_entries,
    load_entry,
)
from xrdkit.symmetry import space_group_operations

# A valid entry, for a library of one under tmp_path; each validation test
# breaks one thing in it.
ENTRY = """
[entry]
name = "test/cubic"
family = "test"
crystal_system = "cubic"
space_group = "Pm-3m"
cell_parameters = ["a"]
z = 1
reference = "none"

[[sites]]
label = "A1"
kind = "A"
wyckoff = "1a"
free = []
uiso_group = "A"

[[sites]]
label = "O1"
kind = "O"
wyckoff = "3c"
"""


def write_entry(root, text: str, name: str = "test/cubic") -> None:
    family, stem = name.split("/")
    folder = root / family
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{stem}.toml").write_text(text, encoding="utf-8")


def multiplicity(wyckoff: str) -> int:
    """The multiplicity of a Wyckoff position written as "18b"."""
    match = re.fullmatch(r"(\d+)[a-z]", wyckoff)
    assert match, f"not a Wyckoff position: {wyckoff!r}"
    return int(match.group(1))


def cell_contents_by_kind(entry: StructureEntry) -> dict[str, int]:
    """Positions per cell of each kind of site."""
    contents: dict[str, int] = {}
    for site in entry.sites:
        contents[site.kind] = contents.get(site.kind, 0) + multiplicity(site.wyckoff)
    return contents


# Shipped entries


def test_ttb_entry() -> None:
    entry = load_entry("ttb/P4bm")

    assert entry.name == "ttb/P4bm"
    assert entry.family == "tetragonal tungsten bronze"
    assert entry.crystal_system == "tetragonal"
    assert entry.space_group == "P4bm"
    assert entry.setting == ""
    assert entry.cell_parameters == ("a", "c")
    assert entry.z == 5
    assert entry.polar_axis == "c"
    assert entry.origin_site == "B1"
    assert entry.reference == "COD 2100720"
    assert len(entry.sites) == 9
    assert [(site.label, site.wyckoff) for site in entry.sites] == [
        ("A1", "2a"),
        ("A2", "4c"),
        ("B1", "2b"),
        ("B2", "8d"),
        ("O1", "4c"),
        ("O2", "8d"),
        ("O3", "8d"),
        ("O4", "2b"),
        ("O5", "8d"),
    ]
    assert entry.sites[1] == Site("A2", "A", "4c", ("x", "z"), "A", ("Ba", "Sr"))
    assert {site.uiso_group for site in entry.sites} == {"A", "B", "O"}


def test_p4mbm_bronze_entry() -> None:
    entry = load_entry("ttb/P4mbm")

    assert entry.family == "tetragonal tungsten bronze"
    assert entry.crystal_system == "tetragonal"
    assert entry.space_group == "P4/mbm"
    assert entry.z == 5
    assert entry.cell_parameters == ("a", "c")
    assert len(entry.sites) == 9
    assert entry.polar_axis is None
    assert entry.origin_site is None
    assert entry.reference == "COD 1563700"
    assert [(site.label, site.wyckoff, site.free) for site in entry.sites] == [
        ("A1", "2a", ()),
        ("A2", "4g", ("x",)),
        ("B1", "2c", ()),
        ("B2", "8j", ("x", "y")),
        ("O1", "4h", ("x",)),
        ("O2", "8j", ("x", "y")),
        ("O3", "8j", ("x", "y")),
        ("O4", "2d", ()),
        ("O5", "8i", ("x", "y")),
    ]


def test_tetragonal_perovskite_entry() -> None:
    entry = load_entry("perovskite/P4mm")

    assert entry.crystal_system == "tetragonal"
    assert entry.space_group == "P4mm"
    assert entry.z == 1
    assert entry.cell_parameters == ("a", "c")
    assert len(entry.sites) == 4
    assert entry.polar_axis == "c"
    assert entry.origin_site == "A1"
    assert [site.wyckoff for site in entry.sites] == ["1a", "1b", "1b", "2c"]


def test_rhombohedral_perovskite_entry() -> None:
    entry = load_entry("perovskite/R3c")

    assert entry.crystal_system == "trigonal"
    assert entry.space_group == "R3c"
    assert entry.setting == "hexagonal"
    assert entry.z == 6
    assert entry.cell_parameters == ("a", "c")
    assert len(entry.sites) == 3
    assert entry.polar_axis == "c"
    assert entry.origin_site == "A1"
    assert [site.wyckoff for site in entry.sites] == ["6a", "6a", "18b"]


def test_orthorhombic_pbnm_perovskite_entry() -> None:
    entry = load_entry("perovskite/Pbnm")

    assert entry.crystal_system == "orthorhombic"
    assert entry.space_group == "Pbnm"
    assert entry.setting == ""
    assert entry.z == 4
    assert entry.cell_parameters == ("a", "b", "c")
    assert len(entry.sites) == 4
    assert entry.polar_axis is None
    assert entry.origin_site is None
    assert [(site.wyckoff, site.free) for site in entry.sites] == [
        ("4c", ("x", "y")),
        ("4b", ()),
        ("4c", ("x", "y")),
        ("8d", ("x", "y", "z")),
    ]


def test_orthorhombic_amm2_perovskite_entry() -> None:
    entry = load_entry("perovskite/Amm2")

    assert entry.crystal_system == "orthorhombic"
    assert entry.space_group == "Amm2"
    assert entry.z == 2
    assert entry.cell_parameters == ("a", "b", "c")
    assert len(entry.sites) == 4
    assert entry.polar_axis == "c"
    assert entry.origin_site == "A1"
    assert [(site.wyckoff, site.free) for site in entry.sites] == [
        ("2a", ("z",)),
        ("2b", ("z",)),
        ("2a", ("z",)),
        ("4e", ("y", "z")),
    ]


@pytest.mark.parametrize(
    ("name", "contents"),
    [
        ("perovskite/Pm-3m", {"A": 1, "B": 1, "O": 3}),
        ("perovskite/P4mm", {"A": 1, "B": 1, "O": 3}),
        ("perovskite/R3c", {"A": 6, "B": 6, "O": 18}),
        ("perovskite/Pbnm", {"A": 4, "B": 4, "O": 12}),
        ("perovskite/Amm2", {"A": 2, "B": 2, "O": 6}),
        # Five AB2O6 on six A positions, a sixth of them empty.
        ("ttb/P4bm", {"A": 6, "B": 10, "O": 30}),
        ("ttb/P4mbm", {"A": 6, "B": 10, "O": 30}),
    ],
)
def test_multiplicities_give_the_cell_contents(
    name: str, contents: dict[str, int]
) -> None:
    assert cell_contents_by_kind(load_entry(name)) == contents


@pytest.mark.parametrize("name", list_entries())
def test_only_r3c_carries_a_setting(name: str) -> None:
    expected = "hexagonal" if name == "perovskite/R3c" else ""

    assert load_entry(name).setting == expected


@pytest.mark.parametrize("name", list_entries())
def test_origin_site_is_free_along_the_polar_axis(name: str) -> None:
    entry = load_entry(name)
    if entry.polar_axis is None:
        assert entry.origin_site is None
        return
    origin = next(site for site in entry.sites if site.label == entry.origin_site)

    assert {"a": "x", "b": "y", "c": "z"}[entry.polar_axis] in origin.free


def test_perovskite_entry() -> None:
    entry = load_entry("perovskite/Pm-3m")

    assert entry.name == "perovskite/Pm-3m"
    assert entry.family == "perovskite"
    assert entry.crystal_system == "cubic"
    assert entry.space_group == "Pm-3m"
    assert entry.cell_parameters == ("a",)
    assert entry.z == 1
    assert entry.polar_axis is None
    assert entry.origin_site is None
    assert "1577245" in entry.reference
    assert entry.sites == (
        Site("A1", "A", "1a", (), "A", ("Sr",)),
        Site("B1", "B", "1b", (), "B", ("Ti",)),
        Site("O1", "O", "3c", (), "O", ("O",)),
    )


def test_list_entries_has_all_seven() -> None:
    names = list_entries()

    assert names == sorted(names)
    assert set(names) >= {
        "perovskite/Amm2",
        "perovskite/P4mm",
        "perovskite/Pbnm",
        "perovskite/Pm-3m",
        "perovskite/R3c",
        "ttb/P4bm",
        "ttb/P4mbm",
    }


def test_unknown_entry_names_the_available_ones() -> None:
    with pytest.raises(ValueError, match="perovskite/Pm-3m"):
        load_entry("ttb/P4/mbm")


@pytest.mark.parametrize("name", list_entries())
def test_every_shipped_entry_loads(name: str) -> None:
    entry = load_entry(name)

    assert entry.name == name
    assert entry.sites
    assert entry.anions == ("O",)
    assert set(entry.bond_limits) <= set(entry.kinds)
    for kind in entry.kinds:
        low, high = bond_limits(entry, kind)
        assert 0 < low < high
    # The tungsten bronzes carry their own limits for bonds to O from the A
    # and B sites; the rest take the package default for every kind.
    expected = {"A": (2.45, 3.0), "B": (1.8, 2.2)} if name.startswith("ttb/") else {}
    assert entry.bond_limits == expected


def test_entry_files_are_package_data() -> None:
    assert (files("xrdkit") / "structures" / "ttb" / "P4bm.toml").is_file()


def test_library_is_exported() -> None:
    for name in ("Site", "StructureEntry", "list_entries", "load_entry"):
        assert name in xrdkit.__all__
    assert xrdkit.load_entry is load_entry


# Validation, on entries under tmp_path


def test_valid_entry_under_a_root(tmp_path) -> None:
    write_entry(tmp_path, ENTRY)

    assert list_entries(root=tmp_path) == ["test/cubic"]
    entry = load_entry("test/cubic", root=tmp_path)
    assert entry.sites[1] == Site("O1", "O", "3c")
    assert entry.setting == ""


def test_missing_entry_table(tmp_path) -> None:
    write_entry(tmp_path, ENTRY[ENTRY.index("[[sites]]") :])

    with pytest.raises(ValueError, match=r"'test/cubic': entry: must be one \[entry\]"):
        load_entry("test/cubic", root=tmp_path)


def test_no_sites(tmp_path) -> None:
    write_entry(tmp_path, ENTRY[: ENTRY.index("[[sites]]")])

    with pytest.raises(ValueError, match=r"'test/cubic': sites: must be at least one"):
        load_entry("test/cubic", root=tmp_path)


def test_duplicate_site_label(tmp_path) -> None:
    write_entry(tmp_path, ENTRY.replace('label = "O1"', 'label = "A1"'))

    with pytest.raises(ValueError, match=r"'test/cubic': sites: label A1"):
        load_entry("test/cubic", root=tmp_path)


def test_free_coordinate_not_xyz(tmp_path) -> None:
    write_entry(tmp_path, ENTRY.replace("free = []", 'free = ["w"]'))

    with pytest.raises(ValueError, match=r"'test/cubic': sites\[0\]\.free: .*'w'"):
        load_entry("test/cubic", root=tmp_path)


def test_site_elements(tmp_path) -> None:
    write_entry(tmp_path, ENTRY.replace("free = []", 'free = []\nelements = ["Sr"]'))

    entry = load_entry("test/cubic", root=tmp_path)
    assert entry.sites[0].elements == ("Sr",)
    assert entry.sites[1].elements == ()


@pytest.mark.parametrize(
    ("elements", "message"),
    [
        ('["sr"]', r"sites\[0\]\.elements: must be a list of element symbols"),
        ('"Sr"', r"sites\[0\]\.elements: must be a list of element symbols"),
        ('["Sr", "Sr"]', r"sites\[0\]\.elements: names an element twice"),
    ],
)
def test_site_elements_not_element_symbols(tmp_path, elements, message) -> None:
    write_entry(
        tmp_path, ENTRY.replace("free = []", f"free = []\nelements = {elements}")
    )

    with pytest.raises(ValueError, match=message):
        load_entry("test/cubic", root=tmp_path)


def test_anions_and_bond_limits(tmp_path) -> None:
    write_entry(
        tmp_path,
        ENTRY.replace(
            'reference = "none"\n',
            'reference = "none"\nanions = ["F", "Cl"]\n\n'
            "[entry.bond_limits]\nA = [2.2, 2.9]\n",
        ),
    )

    entry = load_entry("test/cubic", root=tmp_path)

    assert entry.anions == ("F", "Cl")
    assert entry.bond_limits == {"A": (2.2, 2.9)}
    assert bond_limits(entry, "A") == (2.2, 2.9)
    assert bond_limits(entry, "O") == DEFAULT_BOND_LIMITS == (1.6, 3.0)
    with pytest.raises(
        ValueError,
        match=r"^'B' is not a kind of site of test/cubic; its kinds are A, O$",
    ):
        bond_limits(entry, "B")


def test_anions_and_bond_limits_default(tmp_path) -> None:
    write_entry(tmp_path, ENTRY)

    entry = load_entry("test/cubic", root=tmp_path)

    assert entry.anions == ("O",)
    assert entry.bond_limits == {}
    assert bond_limits(entry, "A") == (1.6, 3.0)


def test_kinds_are_the_entry_own(tmp_path) -> None:
    text = (
        ENTRY.replace('kind = "A"', 'kind = "M"')
        .replace('kind = "O"', 'kind = "X_2"')
        .replace(
            'reference = "none"\n',
            'reference = "none"\n\n[entry.bond_limits]\nM = [1.9, 2.4]\n',
        )
    )
    write_entry(tmp_path, text)

    entry = load_entry("test/cubic", root=tmp_path)

    assert entry.kinds == ("M", "X_2")
    assert bond_limits(entry, "M") == (1.9, 2.4)
    assert bond_limits(entry, "X_2") == (1.6, 3.0)


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        (
            'reference = "none"\n',
            'reference = "none"\nanions = "O"\n',
            r"entry\.anions: must be a list of one or more element symbols",
        ),
        (
            'reference = "none"\n',
            'reference = "none"\nanions = []\n',
            r"entry\.anions: must be a list of one or more element symbols",
        ),
        (
            'reference = "none"\n',
            'reference = "none"\nanions = ["O", "O"]\n',
            r"entry\.anions: names an element twice",
        ),
        (
            'reference = "none"\n',
            'reference = "none"\n\n[entry.bond_limits]\nB = [1.8, 2.2]\n',
            (
                r"entry\.bond_limits\.B: not a kind of the entry's sites; its kinds "
                r"are A, O$"
            ),
        ),
        (
            'reference = "none"\n',
            'reference = "none"\n\n[entry.bond_limits]\nA = [2.2, 1.8]\n',
            r"entry\.bond_limits\.A: must be \[min, max\] in angstroms",
        ),
        (
            'reference = "none"\n',
            'reference = "none"\n\n[entry.bond_limits]\nA = [2.2]\n',
            r"entry\.bond_limits\.A: must be \[min, max\] in angstroms",
        ),
        (
            'kind = "A"',
            'kind = "A site"',
            r"sites\[0\]\.kind: must be a short label",
        ),
    ],
)
def test_anions_bond_limits_and_kinds_rejected(tmp_path, old, new, message) -> None:
    write_entry(tmp_path, ENTRY.replace(old, new))

    with pytest.raises(ValueError, match=message) as raised:
        load_entry("test/cubic", root=tmp_path)
    assert "\n" not in str(raised.value)


# A point on each Wyckoff position the entries use, as the position's own
# form (International Tables) with x, y and z given values of no symmetry.
_X, _Y, _Z = Fraction(13, 100), Fraction(29, 100), Fraction(37, 100)
_H, _Q = Fraction(1, 2), Fraction(1, 4)
WYCKOFF_POINTS = {
    ("P4bm", "2a"): (0, 0, _Z),
    ("P4bm", "2b"): (0, _H, _Z),
    ("P4bm", "4c"): (_X, _X + _H, _Z),
    ("P4bm", "8d"): (_X, _Y, _Z),
    ("P4/mbm", "2a"): (0, 0, 0),
    ("P4/mbm", "2c"): (0, _H, _H),
    ("P4/mbm", "2d"): (0, _H, 0),
    ("P4/mbm", "4g"): (_X, _X + _H, 0),
    ("P4/mbm", "4h"): (_X, _X + _H, _H),
    ("P4/mbm", "8i"): (_X, _Y, 0),
    ("P4/mbm", "8j"): (_X, _Y, _H),
    ("Amm2", "2a"): (0, 0, _Z),
    ("Amm2", "2b"): (_H, 0, _Z),
    ("Amm2", "4e"): (_H, _Y, _Z),
    ("P4mm", "1a"): (0, 0, _Z),
    ("P4mm", "1b"): (_H, _H, _Z),
    ("P4mm", "2c"): (_H, 0, _Z),
    ("Pbnm", "4b"): (_H, 0, 0),
    ("Pbnm", "4c"): (_X, _Y, _Q),
    ("Pbnm", "8d"): (_X, _Y, _Z),
    ("Pm-3m", "1a"): (0, 0, 0),
    ("Pm-3m", "1b"): (_H, _H, _H),
    ("Pm-3m", "3c"): (0, _H, _H),
    ("R3c", "6a"): (0, 0, _Z),
    ("R3c", "18b"): (_X, _Y, _Z),
}


def site_symmetry_free(space_group: str, point) -> tuple[tuple[str, ...], int]:
    """The coordinates a site at ``point`` leaves free, a coordinate tied to
    an earlier one following it, and the order of its site symmetry group."""
    point = tuple(Fraction(value) for value in point)
    stabiliser = []
    for rotation, translation in space_group_operations(space_group):
        image = [
            sum(rotation[i][j] * point[j] for j in range(3)) + translation[i]
            for i in range(3)
        ]
        if all((image[i] - point[i]) % 1 == 0 for i in range(3)):
            stabiliser.append(rotation)
    # Row reduce R - I over the site symmetry, the columns taken z, y, x, so
    # that a coordinate tied to an earlier one is the one determined.
    order = (2, 1, 0)
    rows = [
        [Fraction(rotation[i][c] - (1 if i == c else 0)) for c in order]
        for rotation in stabiliser
        for i in range(3)
    ]
    determined, rank = [], 0
    for column in range(3):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column] != 0), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i][column] != 0:
                factor = rows[i][column] / rows[rank][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[rank])]
        determined.append(order[column])
        rank += 1
    free = tuple("xyz"[axis] for axis in range(3) if axis not in determined)
    return free, len(stabiliser)


@pytest.mark.parametrize("name", list_entries())
def test_every_shipped_site_frees_what_its_wyckoff_position_allows(name: str) -> None:
    entry = load_entry(name)
    group_order = len(space_group_operations(entry.space_group))

    for site in entry.sites:
        point = WYCKOFF_POINTS[(entry.space_group, site.wyckoff)]
        free, site_order = site_symmetry_free(entry.space_group, point)
        # The point is on the position: its site symmetry has the right order.
        assert site_order * int(site.wyckoff[:-1]) == group_order, site.label
        assert site.free == free, site.label


@pytest.mark.parametrize("name", list_entries())
def test_every_shipped_site_has_its_prototype_elements(name: str) -> None:
    assert all(site.elements for site in load_entry(name).sites)


@pytest.mark.parametrize(
    ("system", "parameters"),
    [
        ("cubic", '["a", "c"]'),
        ("tetragonal", '["a"]'),
        ("hexagonal", '["a", "b", "c"]'),
        ("monoclinic", '["a", "b", "c"]'),
        ("triclinic", '["a", "b", "c", "beta"]'),
    ],
)
def test_cell_parameters_not_those_of_the_crystal_system(
    tmp_path, system: str, parameters: str
) -> None:
    text = ENTRY.replace('"cubic"', f'"{system}"').replace('["a"]', parameters)
    write_entry(tmp_path, text)

    with pytest.raises(ValueError, match=r"'test/cubic': entry\.cell_parameters"):
        load_entry("test/cubic", root=tmp_path)


@pytest.mark.parametrize(
    ("system", "parameters"),
    [
        ("orthorhombic", '["a", "b", "c"]'),
        ("trigonal", '["a", "c"]'),
        ("monoclinic", '["a", "b", "c", "beta"]'),
        ("triclinic", '["a", "b", "c", "alpha", "beta", "gamma"]'),
    ],
)
def test_cell_parameters_of_each_crystal_system(
    tmp_path, system: str, parameters: str
) -> None:
    text = ENTRY.replace('"cubic"', f'"{system}"').replace('["a"]', parameters)
    write_entry(tmp_path, text)

    assert load_entry("test/cubic", root=tmp_path).crystal_system == system


def test_unknown_crystal_system(tmp_path) -> None:
    write_entry(tmp_path, ENTRY.replace('"cubic"', '"rhombic"'))

    with pytest.raises(ValueError, match=r"'test/cubic': entry\.crystal_system"):
        load_entry("test/cubic", root=tmp_path)


@pytest.mark.parametrize("z", ["0", "-5", "2.5", "true"])
def test_z_not_a_positive_integer(tmp_path, z: str) -> None:
    write_entry(tmp_path, ENTRY.replace("z = 1", f"z = {z}"))

    with pytest.raises(ValueError, match=r"'test/cubic': entry\.z"):
        load_entry("test/cubic", root=tmp_path)


def test_name_not_that_of_the_file(tmp_path) -> None:
    write_entry(tmp_path, ENTRY, name="test/other")

    with pytest.raises(ValueError, match=r"'test/other': entry\.name"):
        load_entry("test/other", root=tmp_path)


def test_origin_site_not_a_site(tmp_path) -> None:
    text = ENTRY.replace('reference = "none"', 'reference = "none"\norigin_site = "B1"')
    write_entry(tmp_path, text)

    with pytest.raises(ValueError, match=r"'test/cubic': entry\.origin_site"):
        load_entry("test/cubic", root=tmp_path)
