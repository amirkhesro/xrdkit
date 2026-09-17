"""Tests for xrdkit.project and xrdkit init."""

import re
from dataclasses import replace
from pathlib import Path

import pytest

import xrdkit
from xrdkit.cli import main
from xrdkit.config import ConfigError
from xrdkit.project import (
    PROJECT_FILE,
    Instrument,
    Project,
    Refine,
    Sample,
    StructureSpec,
    find_project,
    load_project,
    load_project_text,
    refine_settings,
    resolved_cell,
    resolved_z,
    results_dir,
)

# A complete project; the files it names are made empty by project_dir.
VALID = """
[project]
name = "SBLNT"
version = 1

[instruments.aeris]
wavelength = [1.540598, 1.544426]
ka2 = true
radius = 240.0
instprm = "data/standards/aeris.instprm"

[instruments.mono]
wavelength = [1.540598]
ka2 = false

[structures.ttb_x010]
library = "ttb/P4bm"
composition = "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
cell = { a = 12.45, c = 3.94 }
z = 5
exchange = [["Sr", "Ba"], ["Nb", "Ti"]]
origin = "B1"

[structures.cod-2100720]
cif = "cifs/2100720.cif"
composition = "Sr0.48Ba0.52Nb2O6"
cell = { a = 12.4844, c = 3.9572 }

[samples.x010_calcined]
file = "data/raw/10c.xrdml"
instrument = "aeris"
structures = ["ttb_x010", "cod-2100720"]
stage = "calcined"
form = "powder"
temperature_c = 1300
archimedes = 5.21
notes = "first batch"

[samples."x0.10.pellet"]
file = "data/raw/10s.xrdml"
instrument = "mono"
structures = ["ttb_x010"]
form = "pellet"
"""

PLACEHOLDERS = (
    "data/raw/10c.xrdml",
    "data/raw/10s.xrdml",
    "cifs/2100720.cif",
    "data/standards/aeris.instprm",
)


@pytest.fixture
def project_dir(tmp_path: Path) -> Path:
    for name in PLACEHOLDERS:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_bytes(b"")
    (tmp_path / PROJECT_FILE).write_text(VALID, encoding="utf-8")
    return tmp_path.resolve()


def rewrite(root: Path, old: str, new: str) -> Path:
    """Replace every ``old`` in the project file with ``new``."""
    path = root / PROJECT_FILE
    text = path.read_text(encoding="utf-8")
    assert old in text, old
    path.write_text(text.replace(old, new), encoding="utf-8")
    return path


# Loading


def test_valid_project_loads_as_written(project_dir: Path) -> None:
    project = load_project(project_dir / PROJECT_FILE)

    assert isinstance(project, Project)
    assert project.root == project_dir
    assert project.name == "SBLNT"
    assert project.version == 1

    assert project.instruments == {
        "aeris": Instrument(
            key="aeris",
            wavelength=(1.540598, 1.544426),
            ka2=True,
            radius=240.0,
            instprm=project_dir / "data/standards/aeris.instprm",
        ),
        "mono": Instrument(key="mono", wavelength=(1.540598,), ka2=False),
    }
    assert project.structures == {
        "ttb_x010": StructureSpec(
            key="ttb_x010",
            composition="Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6",
            library="ttb/P4bm",
            cell={"a": 12.45, "c": 3.94},
            z=5,
            exchange=(("Sr", "Ba"), ("Nb", "Ti")),
            origin="B1",
        ),
        "cod-2100720": StructureSpec(
            key="cod-2100720",
            composition="Sr0.48Ba0.52Nb2O6",
            cif=project_dir / "cifs/2100720.cif",
            cell={"a": 12.4844, "c": 3.9572},
        ),
    }
    assert project.samples == {
        "x010_calcined": Sample(
            key="x010_calcined",
            file=project_dir / "data/raw/10c.xrdml",
            instrument="aeris",
            structures=("ttb_x010", "cod-2100720"),
            form="powder",
            stage="calcined",
            temperature_c=1300.0,
            archimedes=5.21,
            notes="first batch",
        ),
        "x0.10.pellet": Sample(
            key="x0.10.pellet",
            file=project_dir / "data/raw/10s.xrdml",
            instrument="mono",
            structures=("ttb_x010",),
            form="pellet",
        ),
    }


def test_load_project_takes_the_folder(project_dir: Path) -> None:
    assert load_project(project_dir).name == "SBLNT"


def test_load_project_finds_the_file_from_the_current_folder(
    project_dir: Path, monkeypatch
) -> None:
    nested = project_dir / "results" / "rietveld"
    nested.mkdir(parents=True)
    monkeypatch.chdir(nested)

    assert load_project().root == project_dir


def test_a_project_with_only_its_project_table_loads(tmp_path: Path) -> None:
    (tmp_path / PROJECT_FILE).write_text(
        '[project]\nname = "empty"\nversion = 1\n', encoding="utf-8"
    )

    project = load_project(tmp_path)
    assert (project.instruments, project.structures, project.samples) == ({}, {}, {})


def test_load_project_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="xrdkit init"):
        load_project(tmp_path)


def test_load_project_invalid_toml(project_dir: Path) -> None:
    rewrite(project_dir, "version = 1", "version = ")

    with pytest.raises(ValueError, match="not valid TOML"):
        load_project(project_dir)


def test_project_is_exported() -> None:
    for name in ("Project", "find_project", "load_project"):
        assert name in xrdkit.__all__
    assert xrdkit.load_project is load_project


# find_project


def test_find_project_from_a_nested_subfolder(project_dir: Path) -> None:
    nested = project_dir / "data" / "raw" / "deeper"
    nested.mkdir(parents=True)

    assert find_project(nested) == project_dir / PROJECT_FILE


def test_find_project_from_the_current_folder(project_dir: Path, monkeypatch) -> None:
    monkeypatch.chdir(project_dir / "cifs")

    assert find_project() == project_dir / PROJECT_FILE


def test_find_project_fails_outside_a_project(tmp_path: Path) -> None:
    nested = tmp_path / "not" / "a" / "project"
    nested.mkdir(parents=True)

    with pytest.raises(FileNotFoundError, match="run xrdkit init"):
        find_project(nested)


# Validation: each case breaks one rule of the valid project; the error names
# the table and the key.

BROKEN = [
    # Unknown keys anywhere.
    (
        "unknown top level",
        "[project]",
        "[extra]\n[project]",
        r"top level: unknown key 'extra'",
    ),
    (
        "unknown project key",
        "version = 1",
        "version = 1\nauthor = 'A'",
        r"project: unknown key 'author'",
    ),
    (
        "unknown instrument key",
        "ka2 = false",
        "ka2 = false\nsource = 'Cu'",
        r"instruments\.mono: unknown key 'source'",
    ),
    (
        "unknown structure key",
        "z = 5",
        "z = 5\nphase = 'TTB'",
        r"structures\.ttb_x010: unknown key 'phase'",
    ),
    (
        "unknown sample key",
        'form = "pellet"',
        'form = "pellet"\nbatch = 2',
        r"samples\.x0\.10\.pellet: unknown key 'batch'",
    ),
    (
        "unknown cif cell key",
        "a = 12.4844, c",
        "a = 12.4844, d = 1.0, c",
        r"structures\.cod-2100720\.cell: unknown key 'd'",
    ),
    # version equal to 1.
    ("version 2", "version = 1", "version = 2", r"project\.version: must be 1"),
    ("version as text", "version = 1", 'version = "1"', r"project\.version: must be 1"),
    # Exactly one of library or cif.
    (
        "neither library nor cif",
        'cif = "cifs/2100720.cif"\n',
        "",
        r"structures\.cod-2100720: must give library or cif, or both",
    ),
    # Library name in list_entries().
    (
        "unknown library entry",
        '"ttb/P4bm"',
        '"ttb/P42bm"',
        r"structures\.ttb_x010\.library: no library entry 'ttb/P42bm'.*perovskite/Pm-3m",
    ),
    # cif path exists.
    (
        "missing cif",
        "cifs/2100720.cif",
        "cifs/missing.cif",
        r"structures\.cod-2100720\.cif: no such file",
    ),
    # Cell with a library.
    (
        "no cell with library",
        "cell = { a = 12.45, c = 3.94 }\n",
        "",
        r"structures\.ttb_x010\.cell: is required with library",
    ),
    (
        "cell keys not the entry's",
        "{ a = 12.45, c = 3.94 }",
        "{ a = 12.45, b = 12.45, c = 3.94 }",
        r"structures\.ttb_x010\.cell: must give exactly a, c",
    ),
    (
        "cell not positive",
        "{ a = 12.45, c = 3.94 }",
        "{ a = 12.45, c = 0.0 }",
        r"structures\.ttb_x010\.cell\.c: must be positive",
    ),
    (
        "cif cell not positive",
        "a = 12.4844,",
        "a = -12.4844,",
        r"structures\.cod-2100720\.cell\.a: must be positive",
    ),
    # z a positive integer.
    (
        "z zero",
        "z = 5",
        "z = 0",
        r"structures\.ttb_x010\.z: must be a positive whole number",
    ),
    (
        "z fractional",
        "z = 5",
        "z = 5.5",
        r"structures\.ttb_x010\.z: must be a positive whole number",
    ),
    # exchange a list of lists of element symbols.
    (
        "exchange not lists",
        '[["Sr", "Ba"], ["Nb", "Ti"]]',
        '["Sr", "Ba"]',
        r"structures\.ttb_x010\.exchange: must be a list of lists",
    ),
    (
        "exchange unknown element",
        '["Nb", "Ti"]',
        '["Nb", "Tx"]',
        r"structures\.ttb_x010\.exchange\[1\]: 'Tx' is not an element",
    ),
    # origin a site label of the entry.
    (
        "origin not a site",
        'origin = "B1"',
        'origin = "Nb1"',
        r"structures\.ttb_x010\.origin: no site 'Nb1' in ttb/P4bm",
    ),
    # composition parses.
    (
        "composition does not parse",
        '"Sr0.48Ba0.52Nb2O6"',
        '"Sr0.48Bq0.52Nb2O6"',
        r"structures\.cod-2100720\.composition: .*Bq",
    ),
    # Sample file exists.
    (
        "missing sample file",
        "data/raw/10s.xrdml",
        "data/raw/missing.xrdml",
        r"samples\.x0\.10\.pellet\.file: no such file",
    ),
    # Instrument key exists.
    (
        "unknown instrument",
        'instrument = "mono"',
        'instrument = "lab"',
        r"samples\.x0\.10\.pellet\.instrument: no instrument 'lab'",
    ),
    # Structures name structure keys.
    (
        "unknown structure",
        '["ttb_x010"]',
        '["ttb_x012"]',
        r"samples\.x0\.10\.pellet\.structures\[0\]: no structure 'ttb_x012'",
    ),
    (
        "no structures",
        '["ttb_x010"]',
        "[]",
        r"samples\.x0\.10\.pellet\.structures: must list one or more",
    ),
    # form powder or pellet.
    (
        "form not powder or pellet",
        'form = "pellet"',
        'form = "film"',
        r"samples\.x0\.10\.pellet\.form: must be powder or pellet",
    ),
    # archimedes positive.
    (
        "archimedes not positive",
        "archimedes = 5.21",
        "archimedes = 0",
        r"samples\.x010_calcined\.archimedes: must be positive",
    ),
    # wavelength one or two positive numbers.
    (
        "no wavelength value",
        "[1.540598]",
        "[]",
        r"instruments\.mono\.wavelength: must be one or two",
    ),
    (
        "three wavelengths",
        "[1.540598]",
        "[1.540598, 1.544426, 1.39222]",
        r"instruments\.mono\.wavelength: must be one or two",
    ),
    (
        "wavelength not positive",
        "[1.540598]",
        "[-1.540598]",
        r"instruments\.mono\.wavelength\[0\]: must be positive",
    ),
    # ka2 true needs two wavelengths, false one.
    (
        "ka2 with one wavelength",
        "[1.540598, 1.544426]",
        "[1.540598]",
        r"instruments\.aeris: ka2 = true needs two wavelengths",
    ),
    (
        "no ka2 with two wavelengths",
        "wavelength = [1.540598]\nka2 = false",
        "wavelength = [1.540598, 1.544426]\nka2 = false",
        r"instruments\.mono: ka2 = false needs one wavelength",
    ),
    # instprm exists when given.
    (
        "missing instprm",
        "data/standards/aeris.instprm",
        "data/standards/missing.instprm",
        r"instruments\.aeris\.instprm: no such file",
    ),
    # Sample and structure keys name folders.
    (
        "sample key with a space",
        '"x0.10.pellet"',
        '"x0.10 pellet"',
        r"samples\.x0\.10 pellet: key 'x0\.10 pellet' must be letters",
    ),
    (
        "structure key with a slash",
        "[structures.cod-2100720]",
        '[structures."cod/2100720"]',
        r"structures\.cod/2100720: key 'cod/2100720' must be letters",
    ),
    # atoms: a library structure's by the entry's site labels, and every
    # element of the composition on a site.
    (
        "atoms on a site the entry lacks",
        'origin = "B1"',
        'origin = "B1"\natoms = { A9 = { Sr1 = "Sr" } }',
        r"structures\.ttb_x010\.atoms: no site 'A9' in ttb/P4bm; its sites are A1,",
    ),
    (
        "atoms with an element that is not one",
        'origin = "B1"',
        'origin = "B1"\natoms = { A1 = { Sr1 = "sr" } }',
        r"structures\.ttb_x010\.atoms\.A1\.atoms\.Sr1: must be an element symbol",
    ),
    (
        "atoms leaving an element the prototype lacks unplaced",
        'origin = "B1"',
        'origin = "B1"\natoms = { A1 = { Sr1 = "Sr", La1 = "La" } }',
        (
            r"structures\.ttb_x010\.atoms: Ti of the composition is not in the "
            r"ttb/P4bm prototype; place it on one of its sites A1, A2, B1, B2, O1, O2, "
            r"O3, O4, O5$"
        ),
    ),
    (
        "atoms taking a prototype element off every site",
        'origin = "B1"',
        (
            'origin = "B1"\natoms = { A1 = { La1 = "La" }, A2 = { Ba2 = "Ba" }, '
            'B1 = { Nb1 = "Nb", Ti1 = "Ti" } }'
        ),
        (
            r"structures\.ttb_x010\.atoms: Sr of the composition is on none of the "
            r"sites once atoms replaces the prototype's"
        ),
    ),
    (
        "cif atoms with a kind that is not a label",
        "c = 3.9572 }",
        'c = 3.9572 }\natoms = [{ atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "Q site" }]',
        r"structures\.cod-2100720\.atoms\[0\]\.kind: must be a short label",
    ),
    (
        "cif atoms not a list",
        "c = 3.9572 }",
        'c = 3.9572 }\natoms = { Nb1 = "Nb" }',
        r"structures\.cod-2100720\.atoms: must be a list of site tables",
    ),
    # origin: a site label, false, or {site, axis} with cif.
    (
        "origin true",
        'origin = "B1"',
        "origin = true",
        (
            r"structures\.ttb_x010\.origin: must be a site label, false, or \{site, "
            r"axis\} with cif, not True"
        ),
    ),
    (
        "origin table with library",
        'origin = "B1"',
        'origin = { site = "B1", axis = "z" }',
        (
            r"structures\.ttb_x010\.origin: the \{site, axis\} form goes with cif; "
            r"with library give one of the sites of ttb/P4bm: A1, A2"
        ),
    ),
    (
        "cif origin axis not x, y or z",
        "c = 3.9572 }",
        'c = 3.9572 }\norigin = { site = "Nb1", axis = "c" }',
        r"structures\.cod-2100720\.origin\.axis: must be x, y or z, not 'c'",
    ),
    (
        "cif origin without an axis",
        "c = 3.9572 }",
        'c = 3.9572 }\norigin = { site = "Nb1" }',
        r"structures\.cod-2100720\.origin: missing key 'axis'",
    ),
    (
        "cif origin not a site of its atoms",
        "c = 3.9572 }",
        (
            'c = 3.9572 }\natoms = [{ atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "B" }]'
            '\norigin = { site = "Nb2", axis = "z" }'
        ),
        (
            r"structures\.cod-2100720\.origin\.site: no site 'Nb2' in atoms; sites are "
            r"named by their first atom: Nb1"
        ),
    ),
    # refine, the project's and a sample's.
    (
        "unknown refine key",
        "[project]",
        "[refine]\ncycles = 10\n\n[project]",
        (
            r"refine: unknown key 'cycles'; the keys are two_theta, background, "
            r"max_passes, unsettled, followed"
        ),
    ),
    (
        "unknown sample refine key",
        'form = "pellet"',
        'form = "pellet"\nrefine = { range = [10.0, 90.0] }',
        r"samples\.x0\.10\.pellet\.refine: unknown key 'range'",
    ),
    (
        "two_theta descending",
        "[project]",
        "[refine]\ntwo_theta = [100.0, 15.0]\n\n[project]",
        r"refine\.two_theta: must rise from low to high within 0 to 180 degrees",
    ),
    (
        "two_theta not a pair",
        "[project]",
        "[refine]\ntwo_theta = [15.0]\n\n[project]",
        r"refine\.two_theta: must be \[low, high\] in degrees",
    ),
    (
        "two_theta not numbers",
        'form = "pellet"',
        'form = "pellet"\nrefine = { two_theta = ["15", 100.0] }',
        r"samples\.x0\.10\.pellet\.refine\.two_theta\[0\]: must be a number",
    ),
    (
        "unknown background key",
        "[project]",
        '[refine]\nbackground = { function = "chebyschev-1", order = 6 }\n\n[project]',
        r"refine\.background: unknown key 'order'",
    ),
    (
        "background terms zero",
        "[project]",
        "[refine]\nbackground = { terms = 0 }\n\n[project]",
        r"refine\.background\.terms: must be a positive whole number, not 0",
    ),
    (
        "background function not text",
        "[project]",
        "[refine]\nbackground = { function = 1 }\n\n[project]",
        r"refine\.background\.function: must be a non-empty string",
    ),
    (
        "max_passes zero",
        'form = "pellet"',
        'form = "pellet"\nrefine = { max_passes = { coordinates = 0 } }',
        (
            r"samples\.x0\.10\.pellet\.refine\.max_passes\.coordinates: must be a "
            r"positive whole number, not 0"
        ),
    ),
    (
        "max_passes unknown mode",
        "[project]",
        '[refine]\nmax_passes = { "fixed-atoms" = 60 }\n\n[project]',
        (
            r"refine\.max_passes: unknown key 'fixed-atoms'; the keys are lebail, "
            r"fixed_atoms, coordinates, occupancies"
        ),
    ),
    (
        "unsettled rule not accept or reject",
        "[project]",
        '[refine]\nunsettled = { lebail = "keep" }\n\n[project]',
        r"refine\.unsettled\.lebail: must be accept or reject, not 'keep'",
    ),
    (
        "followed not triples",
        "[project]",
        '[refine]\nfollowed = ["211"]\n\n[project]',
        r"refine\.followed\[0\]: must be \[h, k, l\], three whole numbers",
    ),
    (
        "followed 000",
        "[project]",
        "[refine]\nfollowed = [[2, 1, 1], [0, 0, 0]]\n\n[project]",
        r"refine\.followed\[1\]: must be \[h, k, l\], three whole numbers not all 0",
    ),
]


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [case[1:] for case in BROKEN],
    ids=[case[0] for case in BROKEN],
)
def test_validation(project_dir: Path, old: str, new: str, message: str) -> None:
    path = rewrite(project_dir, old, new)

    with pytest.raises(ValueError, match=message) as raised:
        load_project(path)
    assert str(raised.value).startswith(str(path))
    assert "\n" not in str(raised.value)


# atoms, origin and refine


def test_the_valid_project_has_no_additions(project_dir: Path) -> None:
    project = load_project(project_dir)
    ttb = project.structures["ttb_x010"]

    # La and Ti are not in the prototype, but with atoms left out nothing is
    # checked, so a file written before atoms loads as it did.
    assert ttb.atoms is None
    assert (ttb.origin, ttb.origin_axis, ttb.origin_fixed) == ("B1", None, True)
    cod = project.structures["cod-2100720"]
    assert (cod.origin, cod.origin_axis, cod.origin_fixed, cod.atoms) == (
        None,
        None,
        True,
        None,
    )
    assert project.refine == Refine()
    assert all(sample.refine == {} for sample in project.samples.values())


def test_library_atoms_place_the_elements_the_prototype_lacks(
    project_dir: Path,
) -> None:
    rewrite(
        project_dir,
        'origin = "B1"',
        (
            "origin = false\n"
            'atoms = { A1 = { Sr1 = "Sr", La1 = "La" }, '
            'B1 = { Nb1 = "Nb", Ti1 = "Ti" }, B2 = { Nb2 = "Nb", Ti2 = "Ti" } }'
        ),
    )

    ttb = load_project(project_dir).structures["ttb_x010"]

    assert [(site["label"], site["name"], site["atoms"]) for site in ttb.atoms] == [
        ("A1", "Sr1", {"Sr1": "Sr", "La1": "La"}),
        ("B1", "Nb1", {"Nb1": "Nb", "Ti1": "Ti"}),
        ("B2", "Nb2", {"Nb2": "Nb", "Ti2": "Ti"}),
    ]
    assert [(site["wyckoff"], site["kind"]) for site in ttb.atoms] == [
        ("2a", "A"),
        ("2b", "B"),
        ("8d", "B"),
    ]
    # origin = false switches the entry's origin off.
    assert (ttb.origin, ttb.origin_fixed) == (None, False)


def test_library_atoms_with_two_common_elements_to_host_are_refused(
    project_dir: Path,
) -> None:
    # A1 and A2 both hold Sr and Ba, so either could host La on them.
    path = rewrite(
        project_dir,
        'origin = "B1"',
        (
            'origin = "B1"\n'
            'atoms = { A1 = { Sr1 = "Sr", Ba1 = "Ba", La1 = "La" }, '
            'A2 = { Ba2 = "Ba", Sr2 = "Sr", La2 = "La" }, B1 = { Nb1 = "Nb", Ti1 = "Ti" } }'
        ),
    )

    with pytest.raises(
        ConfigError,
        match=(
            r"structures\.ttb_x010\.atoms: La is placed on sites A1, A2, which have "
            r"more than one element in common to host it: Ba, Sr; place it on "
            r"sites that share exactly one$"
        ),
    ) as raised:
        load_project(path)
    assert str(raised.value).startswith(str(path))


def test_library_atoms_with_no_common_element_to_host_are_refused(
    project_dir: Path,
) -> None:
    # A1 holds Sr and B1 Nb, so nothing could host La on both.
    path = rewrite(
        project_dir,
        'origin = "B1"',
        (
            'origin = "B1"\n'
            'atoms = { A1 = { Sr1 = "Sr", La1 = "La" }, '
            'B1 = { Nb1 = "Nb", La2 = "La", Ti1 = "Ti" } }'
        ),
    )

    with pytest.raises(
        ConfigError,
        match=(
            r"structures\.ttb_x010\.atoms: La is placed on sites A1, B1, which have "
            r"no element in common to host it; place it on sites that share "
            r"exactly one$"
        ),
    ) as raised:
        load_project(path)
    assert str(raised.value).startswith(str(path))


def test_library_atoms_may_be_left_out_for_the_prototype_elements(
    project_dir: Path,
) -> None:
    rewrite(project_dir, '"Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"', '"Sr0.5Ba0.5Nb2O6"')
    rewrite(project_dir, 'origin = "B1"\n', "")

    ttb = load_project(project_dir).structures["ttb_x010"]

    assert ttb.atoms is None
    # The entry's own origin, neither given nor switched off.
    assert (ttb.origin, ttb.origin_axis, ttb.origin_fixed) == (None, None, True)


def test_cif_atoms_and_origin_site_and_axis(project_dir: Path) -> None:
    rewrite(
        project_dir,
        "c = 3.9572 }",
        (
            "c = 3.9572 }\n"
            "atoms = [\n"
            '    { atoms = { Ba2 = "Ba", Sr2 = "Sr" }, wyckoff = "4c", kind = "A" },\n'
            '    { atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "B" },\n'
            '    { atoms = { O1 = "O" }, wyckoff = "4c", kind = "O" },\n'
            "]\n"
            'origin = { site = "Nb1", axis = "z" }'
        ),
    )

    cod = load_project(project_dir).structures["cod-2100720"]

    assert cod.atoms == (
        {
            "name": "Ba2",
            "atoms": {"Ba2": "Ba", "Sr2": "Sr"},
            "wyckoff": "4c",
            "kind": "A",
        },
        {"name": "Nb1", "atoms": {"Nb1": "Nb"}, "wyckoff": "2b", "kind": "B"},
        {"name": "O1", "atoms": {"O1": "O"}, "wyckoff": "4c", "kind": "O"},
    )
    assert (cod.origin, cod.origin_axis, cod.origin_fixed) == ("Nb1", "z", True)


def test_library_and_cif_together(project_dir: Path) -> None:
    rewrite(
        project_dir,
        'library = "ttb/P4bm"',
        'library = "ttb/P4bm"\ncif = "cifs/2100720.cif"',
    )

    ttb = load_project(project_dir).structures["ttb_x010"]

    assert ttb.library == "ttb/P4bm"
    assert ttb.cif == project_dir / "cifs/2100720.cif"


def test_cif_atoms_take_any_kinds(project_dir: Path) -> None:
    rewrite(
        project_dir,
        "c = 3.9572 }",
        "c = 3.9572 }\n"
        "atoms = [\n"
        '    { atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "M" },\n'
        '    { atoms = { F1 = "F" }, wyckoff = "4c", kind = "X_1" },\n'
        "]",
    )

    cod = load_project(project_dir).structures["cod-2100720"]

    assert [(site["name"], site["kind"]) for site in cod.atoms] == [
        ("Nb1", "M"),
        ("F1", "X_1"),
    ]


def test_cif_origin_label_without_atoms_is_unchecked(project_dir: Path) -> None:
    rewrite(project_dir, "c = 3.9572 }", 'c = 3.9572 }\norigin = "Nb1"')

    cod = load_project(project_dir).structures["cod-2100720"]

    assert (cod.origin, cod.origin_axis, cod.origin_fixed) == ("Nb1", None, True)


REFINE = """
[refine]
two_theta = [15.0, 100.0]
background = { terms = 8 }
max_passes = { lebail = 30, coordinates = 150 }
unsettled = { occupancies = "reject" }
followed = [[2, 1, 1], [4, 0, 0]]
"""


def test_refine_settings_without_a_sample_override(project_dir: Path) -> None:
    path = project_dir / PROJECT_FILE
    path.write_text(path.read_text(encoding="utf-8") + REFINE, encoding="utf-8")

    project = load_project(project_dir)
    expected = Refine(
        two_theta=(15.0, 100.0),
        background={"function": "chebyschev-1", "terms": 8},
        max_passes={
            "lebail": 30,
            "fixed_atoms": 60,
            "coordinates": 150,
            "occupancies": 100,
        },
        unsettled={
            "lebail": "accept",
            "fixed_atoms": "accept",
            "coordinates": "accept",
            "occupancies": "reject",
        },
        followed=((2, 1, 1), (4, 0, 0)),
    )

    assert project.refine == expected
    assert refine_settings(project, "x010_calcined") == expected
    assert refine_settings(project, project.samples["x0.10.pellet"]) == expected


def test_refine_settings_with_a_sample_override(project_dir: Path) -> None:
    path = project_dir / PROJECT_FILE
    path.write_text(
        path.read_text(encoding="utf-8")
        + "\n[samples.x010_calcined.refine]\n"
        + "two_theta = [17.0, 90.0]\n"
        + 'background = { function = "cosine" }\n'
        + "max_passes = { lebail = 5 }\n"
        + "followed = []\n"
        + REFINE,
        encoding="utf-8",
    )

    project = load_project(project_dir)
    calcined = refine_settings(project, "x010_calcined")

    assert project.samples["x010_calcined"].refine == {
        "two_theta": (17.0, 90.0),
        "background": {"function": "cosine"},
        "max_passes": {"lebail": 5},
        "followed": (),
    }
    assert calcined.two_theta == (17.0, 90.0)
    assert calcined.background == {"function": "cosine", "terms": 8}
    assert calcined.max_passes == {
        "lebail": 5,
        "fixed_atoms": 60,
        "coordinates": 150,
        "occupancies": 100,
    }
    assert calcined.unsettled == project.refine.unsettled
    assert calcined.followed == ()
    # The other sample keeps the project's settings, and the project's are
    # untouched by the override.
    assert refine_settings(project, "x0.10.pellet") == project.refine
    assert project.refine.max_passes["lebail"] == 30


def test_refine_settings_package_defaults(project_dir: Path) -> None:
    settings = refine_settings(load_project(project_dir), "x010_calcined")

    assert settings.two_theta is None
    assert settings.background == {"function": "chebyschev-1", "terms": 6}
    assert settings.max_passes == {
        "lebail": 60,
        "fixed_atoms": 60,
        "coordinates": 100,
        "occupancies": 100,
    }
    assert settings.unsettled == dict.fromkeys(
        ("lebail", "fixed_atoms", "coordinates", "occupancies"), "accept"
    )
    assert settings.followed == ()


def test_refine_settings_unknown_sample(project_dir: Path) -> None:
    with pytest.raises(KeyError, match="no sample 'x9' in the project"):
        refine_settings(load_project(project_dir), "x9")


# load_project_text, resolved_z and resolved_cell


def test_load_project_text_checks_text_as_the_file_at_path(project_dir: Path) -> None:
    path = project_dir / PROJECT_FILE
    text = path.read_text(encoding="utf-8")

    assert load_project_text(text, path) == load_project(path)
    with pytest.raises(ValueError, match=r"samples\.x0\.10\.pellet\.form"):
        load_project_text(text.replace('"pellet"', '"film"'), path)
    assert path.read_text(encoding="utf-8") == text


def test_resolved_z(project_dir: Path) -> None:
    project = load_project(project_dir)
    ttb = project.structures["ttb_x010"]
    cod = project.structures["cod-2100720"]

    assert resolved_z(ttb) == 5
    assert resolved_z(replace(ttb, z=None)) == 5
    assert resolved_z(replace(cod, z=10)) == 10
    with pytest.raises(ValueError, match=r"structures\.cod-2100720: uses a CIF"):
        resolved_z(cod)


def test_resolved_cell(project_dir: Path) -> None:
    project = load_project(project_dir)

    assert resolved_cell(project.structures["ttb_x010"]) == {
        "a": 12.45,
        "b": 12.45,
        "c": 3.94,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    with pytest.raises(ValueError, match="uses a CIF, so its cell needs all of"):
        resolved_cell(project.structures["cod-2100720"])


# results_dir


def test_results_dir(project_dir: Path) -> None:
    project = load_project(project_dir)
    sample = project.samples["x010_calcined"]

    expected = project_dir / "results" / "rietveld" / "x010_calcined"
    assert results_dir(project, "rietveld", sample) == expected
    assert results_dir(project, "rietveld", "x010_calcined") == expected
    assert not (project_dir / "results").exists()


# xrdkit init


def test_init_writes_a_loadable_file(tmp_path: Path, monkeypatch, capsys) -> None:
    folder = tmp_path / "my-project"
    folder.mkdir()
    monkeypatch.chdir(folder)

    assert main(["init"]) == 0

    path = folder / PROJECT_FILE
    assert capsys.readouterr().out.strip() == str(path)
    project = load_project(path)
    assert project.name == "my-project"
    assert project.version == 1
    assert (project.instruments, project.structures, project.samples) == ({}, {}, {})
    for name in ("data/raw", "cifs", "results"):
        assert (folder / name).is_dir()


def test_init_refuses_a_second_time(tmp_path: Path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)
    assert main(["init"]) == 0
    written = (tmp_path / PROJECT_FILE).read_text(encoding="utf-8")
    capsys.readouterr()

    assert main(["init", "--name", "other"]) == 1

    error = capsys.readouterr().err
    assert "xrdkit init:" in error
    assert "already exists" in error
    assert (tmp_path / PROJECT_FILE).read_text(encoding="utf-8") == written


@pytest.mark.parametrize("name", ["SBLNT", 'Amir\'s "bronzes" \\ 2026 Kα'])
def test_init_name(tmp_path: Path, monkeypatch, name: str) -> None:
    monkeypatch.chdir(tmp_path)

    assert main(["init", "--name", name]) == 0
    assert load_project(tmp_path).name == name


def test_init_refuses_a_blank_name(tmp_path: Path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)

    assert main(["init", "--name", "  "]) == 1
    assert "blank" in capsys.readouterr().err
    assert not (tmp_path / PROJECT_FILE).exists()


def test_init_examples_load_once_uncommented(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert main(["init"]) == 0
    path = tmp_path / PROJECT_FILE
    text = path.read_text(encoding="utf-8")
    uncommented = re.sub(r"(?m)^# (\[|\w+ = )", r"\1", text)
    path.write_text(uncommented, encoding="utf-8")
    (tmp_path / "data/raw/sample.xrdml").write_bytes(b"")
    (tmp_path / "data/standards").mkdir()
    (tmp_path / "data/standards/diffractometer.instprm").write_bytes(b"")

    project = load_project(path)

    assert list(project.instruments) == ["diffractometer"]
    assert project.structures["phase1"].library == "ttb/P4bm"
    assert project.samples["sample1"].structures == ("phase1",)
    assert project.refine.two_theta == (15.0, 100.0)
    assert project.refine.followed == ((2, 1, 1), (4, 0, 0))
    assert refine_settings(project, "sample1").background == {
        "function": "chebyschev-1",
        "terms": 8,
    }
