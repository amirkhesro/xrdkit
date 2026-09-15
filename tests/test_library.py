"""Tests for xrdkit.library."""

from importlib.resources import files

import pytest

import xrdkit
from xrdkit.config import wyckoff_multiplicity
from xrdkit.library import Site, StructureEntry, list_entries, load_entry

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


def cell_contents_by_kind(entry: StructureEntry) -> dict[str, int]:
    contents: dict[str, int] = {}
    for site in entry.sites:
        contents[site.kind] = contents.get(site.kind, 0) + wyckoff_multiplicity(
            site.wyckoff
        )
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
    assert entry.sites[1] == Site("A2", "A", "4c", ("x", "z"), "A")
    assert {site.uiso_group for site in entry.sites} == {"A", "B", "O"}


def test_ttb_entry_holds_five_ab2o6_with_a_sixth_of_a_empty() -> None:
    assert cell_contents_by_kind(load_entry("ttb/P4bm")) == {"A": 6, "B": 10, "O": 30}


def test_ttb_origin_site_is_free_along_the_polar_axis() -> None:
    entry = load_entry("ttb/P4bm")
    origin = next(site for site in entry.sites if site.label == entry.origin_site)

    assert "z" in origin.free


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
        Site("A1", "A", "1a", (), "A"),
        Site("B1", "B", "1b", (), "B"),
        Site("O1", "O", "3c", (), "O"),
    )
    assert cell_contents_by_kind(entry) == {"A": 1, "B": 1, "O": 3}


def test_list_entries_has_both() -> None:
    names = list_entries()

    assert "ttb/P4bm" in names
    assert "perovskite/Pm-3m" in names
    assert names == sorted(names)


def test_unknown_entry_names_the_available_ones() -> None:
    with pytest.raises(ValueError, match="perovskite/Pm-3m"):
        load_entry("ttb/P4/mbm")


@pytest.mark.parametrize("name", list_entries())
def test_every_shipped_entry_loads(name: str) -> None:
    entry = load_entry(name)

    assert entry.name == name
    assert entry.sites


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
