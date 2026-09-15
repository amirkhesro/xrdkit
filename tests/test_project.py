"""Tests for xrdkit.project and xrdkit init."""

import re
from pathlib import Path

import pytest

import xrdkit
from xrdkit.cli import main
from xrdkit.project import (
    PROJECT_FILE,
    Instrument,
    Project,
    Sample,
    StructureSpec,
    find_project,
    load_project,
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

PLACEHOLDERS = ("data/raw/10c.xrdml", "data/raw/10s.xrdml", "cifs/2100720.cif")


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
        "library and cif",
        'library = "ttb/P4bm"',
        'library = "ttb/P4bm"\ncif = "cifs/2100720.cif"',
        r"structures\.ttb_x010: must give exactly one of library and cif",
    ),
    (
        "neither library nor cif",
        'cif = "cifs/2100720.cif"\n',
        "",
        r"structures\.cod-2100720: must give exactly one of library and cif",
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
    (tmp_path / "data/raw/10c.xrdml").write_bytes(b"")

    project = load_project(path)

    assert list(project.instruments) == ["aeris"]
    assert project.structures["ttb_x010"].library == "ttb/P4bm"
    assert project.samples["x010_calcined"].structures == ("ttb_x010",)
