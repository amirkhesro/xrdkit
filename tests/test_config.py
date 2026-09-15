"""Tests for xrdkit.config."""

import copy
import inspect
import re
import tomllib

import pytest

from xrdkit import gsas2_driver as driver
from xrdkit.config import (
    ConfigError,
    check_composition,
    load_config,
    sample_settings,
    validate_config,
    wyckoff_multiplicity,
)
from xrdkit.density import theoretical_density
from xrdkit.structure import bond_lengths, composition_edits, site_setup

# A perovskite-like structure, A on two sites, and one sample on it.
SETTINGS = """
[samples.x10_powder]
id = "10"
scan = "raw/10.xrdml"
composition = { Sr = 0.4, Ba = 0.5, La = 0.1, Ti = 1.0, O = 3.0 }
structure = "abo"
start_cell = { file = "results/lattice/lattice_fits.csv", model = "zero" }
two_theta = [17.0, 100.0]
background = { function = "chebyschev-1", terms = 6 }
refine_microstrain = false
notes = "decided"

[samples.x10_powder.trials]
runs = [{ terms = 6 }, { low = 17.0, terms = 8 }]
followed = ["100"]

[samples.x10_powder.write_up.lebail]
range = "why 17"

[samples.x12_powder]
id = "12"
scan = "raw/12.xrdml"
composition = { Sr = 0.38, Ba = 0.5, La = 0.12, Ti = 1.0, O = 3.0 }
structure = "abo"
start_cell = { a = 4.0, c = 4.1 }
two_theta = [10.0, 90.0]
background = { function = "chebyschev-1", terms = 6 }
refine_microstrain = true
notes = ""

[structures.abo]
cif = "cifs/abo.cif"
label = "COD 1"
phase_name = "P"
space_group = "P4mm"
formula_units = 1
sites = [
    { atoms = { Ba1 = "Ba", Sr1 = "Sr" }, wyckoff = "1a", kind = "A" },
    { atoms = { Ti1 = "Ti" }, wyckoff = "1b", kind = "B" },
    { atoms = { O1 = "O" }, wyckoff = "1b", kind = "O" },
    { atoms = { O2 = "O" }, wyckoff = "2c", kind = "O" },
]
free_coordinates = { "1a" = "z", "2c" = "xz" }
uiso_groups = [
    { name = "A", sites = ["Ba1"] },
    { name = "B", sites = ["Ti1"] },
    { name = "O", sites = ["O1", "O2"] },
]
origin = { site = "Ti1", axis = "z" }
exchange = { elements = ["Sr", "Ba"], sites = ["Ba1", "Ba1x"] }
bond_limits = { B = { min = 1.8, max = 2.2 } }
composition = { added = { La = "Sr" } }
"""


def settings() -> dict:
    data = tomllib.loads(SETTINGS)
    # The exchange needs two sites; add a second A site for it.
    abo = data["structures"]["abo"]
    abo["sites"].insert(1, {"atoms": {"Ba1x": "Ba"}, "wyckoff": "1a", "kind": "A"})
    abo["uiso_groups"][0]["sites"].append("Ba1x")
    return data


def test_load_config_reads_and_fills_in(tmp_path) -> None:
    path = tmp_path / "samples.toml"
    path.write_text(SETTINGS, encoding="utf-8")
    with pytest.raises(ConfigError, match="exchange.sites.*no site named 'Ba1x'"):
        load_config(path)

    config = validate_config(settings(), "samples.toml")

    sample, structure = sample_settings(config, "10")
    assert sample["name"] == "x10_powder" and structure["name"] == "abo"
    assert sample["two_theta"] == (17.0, 100.0)
    assert sample["trials"]["runs"] == [
        {"low": None, "terms": 6},
        {"low": 17.0, "terms": 8},
    ]
    assert sample["write_up"] == {"lebail": {"range": "why 17"}}
    assert [site["name"] for site in structure["sites"]] == [
        "Ba1",
        "Ba1x",
        "Ti1",
        "O1",
        "O2",
    ]
    other, _ = sample_settings(config, "x12_powder")
    assert other["start_cell"] == {"a": 4.0, "c": 4.1}
    # Optional keys come back filled in.
    assert other["trials"] is None
    assert other["followed_reflections"] == [] and other["write_up"] == {}
    with pytest.raises(ConfigError, match="no sample '13'.*10 \\(x10_powder\\)"):
        sample_settings(config, "13")


def test_load_config_rejects_bad_toml(tmp_path) -> None:
    path = tmp_path / "bad.toml"
    path.write_text("[samples\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="not valid TOML"):
        load_config(path)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (
            lambda d: d["samples"]["x10_powder"].pop("scan"),
            "samples.x10_powder: missing key 'scan'",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(backgroud=1),
            "unknown key 'backgroud'",
        ),
        (
            lambda d: d["structures"]["abo"].pop("sites"),
            "structures.abo: missing key 'sites'",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(structure="sbn"),
            "no structure 'sbn' under structures; there are abo",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(
                start_cell={"a": 4.0, "b": 4.1}
            ),
            "start_cell: cannot tell the crystal system from a, b",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(start_cell={}),
            "must give file and model",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(two_theta=[50.0, 20.0]),
            "must rise from low to high",
        ),
        (
            lambda d: d["samples"]["x10_powder"]["background"].update(terms=0),
            "background.terms: must be a whole number of at least 1",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(refine_microstrain="no"),
            "must be true or false",
        ),
        (
            lambda d: d["samples"]["x12_powder"].update(id="10"),
            "'10' is the id of samples.x10_powder too",
        ),
        (
            lambda d: d["structures"]["abo"]["sites"][0].update(kind="C"),
            "kind: must be one of A, B, O",
        ),
        (
            lambda d: d["structures"]["abo"]["sites"][0].update(wyckoff="a1"),
            "not a Wyckoff position",
        ),
        (
            lambda d: d["structures"]["abo"]["uiso_groups"][2]["sites"].pop(),
            "sites O2 are in no Uiso group",
        ),
        (
            lambda d: d["structures"]["abo"]["uiso_groups"][1]["sites"].append("O1"),
            "site O1 is in Uiso group 'O' already|site O1 is in Uiso group 'B' already",
        ),
        (
            lambda d: d["structures"]["abo"]["origin"].update(site="Nb1"),
            "origin.site: no site named 'Nb1'",
        ),
        (
            lambda d: d["structures"]["abo"]["exchange"].update(sites=["Ba1", "Ti1"]),
            "must all be of one kind",
        ),
        (
            lambda d: d["structures"]["abo"]["composition"]["added"].update(Nb="Ta"),
            "its host Ta is on none of the sites",
        ),
        (
            lambda d: d["structures"]["abo"]["free_coordinates"].update({"2c": "xw"}),
            "must be 'all' or some of 'xyz'",
        ),
    ],
)
def test_validate_config_names_the_problem(change, message) -> None:
    data = settings()
    change(data)
    with pytest.raises(ConfigError, match=message) as error:
        validate_config(data, "samples.toml")
    assert str(error.value).startswith("samples.toml: ")


def test_composition_must_match_the_sites() -> None:
    config = validate_config(settings())
    structure = config["structures"]["abo"]
    good = {"Sr": 0.4, "Ba": 0.5, "La": 0.1, "Ti": 1.0, "O": 3.0}
    check_composition(good, structure)

    with pytest.raises(ConfigError, match="lacks Ba, which the sites"):
        check_composition({k: v for k, v in good.items() if k != "Ba"}, structure)
    with pytest.raises(ConfigError, match="has Nb, which no site"):
        check_composition({**good, "Nb": 0.1}, structure)
    # Two A sites of multiplicity 1 hold at most 2 atoms per cell; the three
    # O positions exactly 3.
    check_composition({**good, "Ba": 1.5}, structure)
    with pytest.raises(
        ConfigError,
        match=r"puts 2\.1 atoms per cell .* on the A sites, which have 2 positions",
    ):
        check_composition({**good, "Ba": 1.6}, structure)
    with pytest.raises(ConfigError, match="on the O sites, which have 3 positions"):
        check_composition({**good, "O": 3.5}, structure)

    # The same checks run when the file is read.
    data = settings()
    data["samples"]["x12_powder"]["composition"]["O"] = 4.0
    with pytest.raises(ConfigError, match="samples.x12_powder.composition: puts 4"):
        validate_config(data)


def test_kinds_sharing_an_element_count_together() -> None:
    data = settings()
    # Ti on the A sites too: A and B are then one pool of 3 positions.
    data["structures"]["abo"]["sites"][1]["atoms"] = {"Ti2": "Ti"}
    data["structures"]["abo"]["exchange"]["sites"] = ["Ba1", "Ba1x"]
    data["structures"]["abo"]["sites"].append(
        {"atoms": {"Ba1x": "Ba"}, "wyckoff": "1a", "kind": "A"}
    )
    data["structures"]["abo"]["uiso_groups"][1]["sites"].append("Ti2")
    config = validate_config(copy.deepcopy(data))
    structure = config["structures"]["abo"]
    composition = {"Sr": 0.4, "Ba": 0.5, "La": 0.1, "Ti": 2.0, "O": 3.0}
    check_composition(composition, structure)
    with pytest.raises(ConfigError, match="on the A/B sites, which have 4 positions"):
        check_composition({**composition, "Ti": 3.1}, structure)


def test_wyckoff_multiplicity() -> None:
    assert wyckoff_multiplicity("8d") == 8
    assert wyckoff_multiplicity("16k") == 16
    with pytest.raises(ValueError, match="not a Wyckoff position"):
        wyckoff_multiplicity("d8")


# The rule for a stage that has not settled


def test_unsettled_rules_default_at_the_top_level() -> None:
    data = settings()
    data["unsettled"] = {"lebail": "accept", "coordinates": "reject"}

    config = validate_config(data, "samples.toml")

    assert config["unsettled"] == {"lebail": "accept", "coordinates": "reject"}
    for name in ("x10_powder", "x12_powder"):
        sample, _ = sample_settings(config, name)
        assert sample["unsettled"] == {"lebail": "accept", "coordinates": "reject"}


def test_a_sample_overrides_the_default_mode_by_mode() -> None:
    data = settings()
    data["unsettled"] = {"lebail": "accept", "coordinates": "reject"}
    data["samples"]["x10_powder"]["unsettled"] = {
        "coordinates": "accept",
        "occupancies": "reject",
    }

    config = validate_config(data, "samples.toml")

    sample, _ = sample_settings(config, "10")
    assert sample["unsettled"] == {
        "lebail": "accept",
        "coordinates": "accept",
        "occupancies": "reject",
    }
    # The samples that do not override it keep the default.
    other, _ = sample_settings(config, "12")
    assert other["unsettled"] == {"lebail": "accept", "coordinates": "reject"}


def test_unsettled_left_out_altogether() -> None:
    config = validate_config(settings(), "samples.toml")

    assert config["unsettled"] == {}
    assert sample_settings(config, "10")[0]["unsettled"] == {}


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (
            lambda d: d.update(unsettled={"lebail": "roll back"}),
            "unsettled.lebail: must be one of accept, reject, not 'roll back'",
        ),
        (
            lambda d: d["samples"]["x10_powder"].update(unsettled={"lebail": "hold"}),
            "samples.x10_powder.unsettled.lebail: must be one of accept, reject",
        ),
        (
            lambda d: d.update(unsettled="accept"),
            "unsettled: must be a table",
        ),
        (
            lambda d: d.update(unsettle={"lebail": "accept"}),
            "top level: unknown key 'unsettle'",
        ),
    ],
)
def test_bad_unsettled_rules_name_the_problem(change, message) -> None:
    data = settings()
    change(data)
    with pytest.raises(ConfigError, match=message):
        validate_config(data, "samples.toml")


# Structures other than the one the settings were first written for


def abo_atom(label, kind, xyz, multiplicity=1, occupancy=1.0) -> dict:
    return {
        "label": label,
        "type": kind,
        "xyz": xyz,
        "multiplicity": multiplicity,
        "occupancy": occupancy,
    }


# The atoms of the settings' abo structure, as the GSAS-II driver reports them.
ABO_ATOMS = [
    abo_atom("Ba1", "Ba", [0.0, 0.0, 0.0], occupancy=0.5),
    abo_atom("Sr1", "Sr", [0.0, 0.0, 0.0], occupancy=0.5),
    abo_atom("Ba1x", "Ba", [0.0, 0.0, 0.5], occupancy=0.0),
    abo_atom("Ti1", "Ti", [0.5, 0.5, 0.5]),
    abo_atom("O1", "O", [0.5, 0.5, 0.0]),
    abo_atom("O2", "O", [0.5, 0.0, 0.5], multiplicity=2),
]


def test_exchange_is_optional() -> None:
    data = settings()
    del data["structures"]["abo"]["exchange"]

    structure = validate_config(data)["structures"]["abo"]

    assert structure["exchange"] is None
    assert site_setup(structure, ABO_ATOMS)["exchange"] is None
    # A stage builder with no exchange passes None, and no constraint is set.
    (stage,) = driver.accumulate_stages([{"name": "none", "occupancies": None}])
    assert stage["flags"]["occupancies"] == []


def test_origin_is_optional() -> None:
    data = settings()
    del data["structures"]["abo"]["origin"]

    structure = validate_config(data)["structures"]["abo"]

    assert structure["origin"] is None
    setup = site_setup(structure, ABO_ATOMS)
    assert setup["origin"] is None and setup["origin_axis"] is None
    # No coordinate is held for the origin: Ti1's free coordinates are refined.
    assert setup["coordinates"] == {
        "A": {"Ba1": "z", "Ba1x": "z"},
        "B": {"Ti1": "all"},
        "O": {"O1": "all", "O2": "xz"},
    }
    (stage,) = driver.accumulate_stages([{"name": "none", "origin": None}])
    assert stage["flags"]["origin"] == {}


def test_the_sites_need_not_hold_every_kind() -> None:
    data = settings()
    abo = data["structures"]["abo"]
    abo["sites"] = [site for site in abo["sites"] if site["kind"] != "B"]
    abo["uiso_groups"] = [group for group in abo["uiso_groups"] if group["name"] != "B"]
    abo["origin"]["site"] = "Ba1"
    abo["bond_limits"] = {"A": {"min": 2.0}}
    for sample in data["samples"].values():
        del sample["composition"]["Ti"]

    structure = validate_config(data)["structures"]["abo"]

    assert sorted({site["kind"] for site in structure["sites"]}) == ["A", "O"]
    # Bond lengths of a phase with no O: none, and a note saying why.
    phase = {
        "cell": [4.0, 4.0, 4.0, 90.0, 90.0, 90.0],
        "atoms": [
            {"label": "Ba1", "type": "Ba", "xyz": [0.0, 0.0, 0.0]},
            {"label": "Ti1", "type": "Ti", "xyz": [0.5, 0.5, 0.5]},
        ],
        "operators": [
            {"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "translation": [0, 0, 0]}
        ],
    }
    with pytest.warns(UserWarning, match="no site of an anion, O"):
        assert bond_lengths(phase) == []


def test_start_cell_by_crystal_system() -> None:
    def start_cell(cell: dict, crystal_system: str | None = None) -> dict:
        data = settings()
        data["samples"]["x12_powder"]["start_cell"] = cell
        if crystal_system is not None:
            data["structures"]["abo"]["crystal_system"] = crystal_system
        config = validate_config(data, "samples.toml")
        return sample_settings(config, "12")[0]["start_cell"]

    # Without crystal_system the keys tell cubic, tetragonal and orthorhombic.
    assert start_cell({"a": 4.0}) == {"a": 4.0}
    assert start_cell({"a": 4.0, "c": 4.1}) == {"a": 4.0, "c": 4.1}
    assert start_cell({"c": 4.2, "b": 4.1, "a": 4.0}) == {"a": 4.0, "b": 4.1, "c": 4.2}
    # With it, the keys must be the crystal system's.
    assert start_cell({"a": 5.6, "c": 13.9}, "hexagonal") == {"a": 5.6, "c": 13.9}
    monoclinic = {"a": 5.0, "b": 6.0, "c": 7.0, "beta": 95.0}
    assert start_cell(monoclinic, "monoclinic") == monoclinic

    for cell, crystal_system, message in [
        (
            {"a": 0.0, "c": 4.1},
            None,
            "samples.x12_powder.start_cell.a: must be greater than 0, not 0.0",
        ),
        (
            {"a": 4.0, "b": 4.1},
            None,
            (
                "samples.x12_powder.start_cell: cannot tell the crystal system from "
                "a, b; give crystal_system in structures.abo"
            ),
        ),
        (
            monoclinic,
            None,
            "cannot tell the crystal system from a, b, c, beta; give crystal_system",
        ),
        (
            {"a": 4.0, "b": 4.0, "c": 4.1},
            "tetragonal",
            (
                "samples.x12_powder.start_cell: a tetragonal cell (the crystal_system "
                "of structures.abo) has a, c, not a, b, c"
            ),
        ),
        (
            {"a": 4.0, "q": 1.0},
            None,
            (
                "samples.x12_powder.start_cell: unknown key 'q'; the cell parameters "
                "are a, b, c, alpha, beta, gamma"
            ),
        ),
        (
            {"a": 4.0},
            "rhombic",
            (
                "structures.abo.crystal_system: must be one of cubic, tetragonal, "
                "orthorhombic, hexagonal, trigonal, monoclinic, triclinic, not "
                "'rhombic'"
            ),
        ),
    ]:
        with pytest.raises(ConfigError, match=re.escape(message)):
            start_cell(cell, crystal_system)


def test_formula_units_is_optional() -> None:
    data = settings()
    del data["structures"]["abo"]["formula_units"]

    structure = validate_config(data)["structures"]["abo"]

    assert structure["formula_units"] is None
    composition = {"Sr": 0.5, "Ba": 0.5, "Ti": 1.0, "O": 3.0}
    with pytest.raises(ValueError, match="structure abo gives no formula_units"):
        composition_edits(ABO_ATOMS, composition, structure)
    # Nothing assumes a Z: theoretical_density takes one.
    z = inspect.signature(theoretical_density).parameters["z"]
    assert z.default is inspect.Parameter.empty


# A structure named by a library entry, its atoms mapped onto the entry's
# sites, and one sample on it.
LIBRARY_SETTINGS = """
[samples.bronze]
id = "b1"
scan = "raw/b1.xrdml"
composition = { Sr = 0.5, Ba = 0.5, Nb = 2.0, O = 6.0 }
structure = "bronze"
start_cell = { a = 12.45, c = 3.94 }
two_theta = [15.0, 100.0]
background = { function = "chebyschev-1", terms = 6 }
refine_microstrain = false
notes = ""

[structures.bronze]
library = "ttb/P4bm"
cif = "cifs/bronze.cif"
label = "a bronze"
phase_name = "TTB"
exchange = { elements = ["Sr", "Ba"], sites = ["A1", "A2"] }
composition = { added = {} }

[structures.bronze.atoms]
A1 = { Sr1 = "Sr" }
A2 = { Ba2 = "Ba", Sr2 = "Sr" }
B1 = { Nb1 = "Nb" }
B2 = { Nb2 = "Nb" }
O1 = { O1 = "O" }
O2 = { O2 = "O" }
O3 = { O3 = "O" }
O4 = { O4 = "O" }
O5 = { O5 = "O" }
"""


def test_a_structure_from_the_library() -> None:
    structure = validate_config(tomllib.loads(LIBRARY_SETTINGS))["structures"]["bronze"]

    assert structure["library"] == "ttb/P4bm"
    assert structure["space_group"] == "P4bm"
    assert structure["formula_units"] == 5
    assert structure["crystal_system"] == "tetragonal"
    assert [
        (site["name"], site["label"], site["wyckoff"], site["kind"])
        for site in structure["sites"]
    ] == [
        ("Sr1", "A1", "2a", "A"),
        ("Ba2", "A2", "4c", "A"),
        ("Nb1", "B1", "2b", "B"),
        ("Nb2", "B2", "8d", "B"),
        ("O1", "O1", "4c", "O"),
        ("O2", "O2", "8d", "O"),
        ("O3", "O3", "8d", "O"),
        ("O4", "O4", "2b", "O"),
        ("O5", "O5", "8d", "O"),
    ]
    assert structure["sites"][1]["atoms"] == {"Ba2": "Ba", "Sr2": "Sr"}
    assert structure["uiso_groups"] == [
        {"name": "A", "sites": ["Sr1", "Ba2"]},
        {"name": "B", "sites": ["Nb1", "Nb2"]},
        {"name": "O", "sites": ["O1", "O2", "O3", "O4", "O5"]},
    ]
    assert structure["origin"] == {"site": "Nb1", "axis": "z"}
    assert structure["free_coordinates"] == {
        "2a": "z",
        "4c": "xz",
        "2b": "z",
        "8d": "xyz",
    }
    # The entry's site labels name the sites, as their first atoms do.
    assert structure["exchange"] == {"elements": ["Sr", "Ba"], "sites": ["Sr1", "Ba2"]}

    # Any key given explicitly overrides the entry.
    data = tomllib.loads(LIBRARY_SETTINGS)
    bronze = data["structures"]["bronze"]
    bronze["origin"] = {"site": "A2", "axis": "z"}
    bronze["formula_units"] = 4
    bronze["uiso_groups"] = [
        {"name": "all", "sites": ["A1", "A2", "B1", "B2", "O1", "O2", "O3", "O4", "O5"]}
    ]
    overridden = validate_config(data)["structures"]["bronze"]
    assert overridden["origin"] == {"site": "Ba2", "axis": "z"}
    assert overridden["formula_units"] == 4
    assert overridden["uiso_groups"] == [
        {"name": "all", "sites": [site["name"] for site in structure["sites"]]}
    ]

    # Every site of the entry needs its atoms, and the entry must exist.
    del bronze["atoms"]["O5"]
    with pytest.raises(
        ConfigError, match=re.escape("structures.bronze.atoms: no atoms for site O5")
    ):
        validate_config(data)
    bronze["library"] = "ttb/P4mm"
    with pytest.raises(
        ConfigError,
        match=re.escape("structures.bronze.library: no structure entry 'ttb/P4mm'"),
    ):
        validate_config(data)
