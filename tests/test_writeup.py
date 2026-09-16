"""Tests for xrdkit.writeup: the markdown write ups on synthetic results in
tests/data/writeup, and run_mode's figure and write up with GSAS-II faked."""

import json
from pathlib import Path

import pytest
from test_pipeline import FakeGsas2, fake_project

from xrdkit import pipeline, writeup
from xrdkit.pipeline import run_mode, run_sequence
from xrdkit.writeup import (
    coordinate_shift,
    fixed_atoms_markdown,
    lebail_markdown,
    reflection_misfits,
    relative,
    structure_markdown,
    with_esd,
)

DATA = Path(__file__).parent / "data" / "writeup"


def load(name: str) -> dict:
    """A synthetic result, its exported reflection lists found in DATA."""
    result = json.loads((DATA / name).read_text(encoding="utf-8"))
    exports = result.get("exports") or {}
    exports["reflections"] = {
        phase: str(DATA / path)
        for phase, path in exports.get("reflections", {}).items()
    }
    return result


def plan() -> dict:
    return json.loads((DATA / "plan.json").read_text(encoding="utf-8"))


# Formatting


def test_with_esd() -> None:
    assert with_esd(1.23456, None) == "1.23456 (fixed)"
    assert with_esd(1.23456, 0.0) == "1.23456 (fixed)"
    assert with_esd(12.4765, 0.0003) == "12.4765(3)"
    # An esd beginning with a 1 keeps two figures.
    assert with_esd(1.23456, 0.0012) == "1.2346(12)"
    # An esd above 1 rounds the value to whole units.
    assert with_esd(123.4, 3.4) == "123(3)"
    assert with_esd(1234.5, 12.0) == "1234(12)"
    assert with_esd(0.5, None, 3) == "0.5 (fixed)"


def test_relative(tmp_path) -> None:
    inside = tmp_path / "results" / "x.json"
    value = {"a": str(inside), "b": [str(tmp_path / "c.csv"), 3], "c": "plain"}

    assert relative(value, tmp_path) == {
        "a": "results/x.json",
        "b": ["c.csv", 3],
        "c": "plain",
    }


# The write ups


def test_lebail_markdown() -> None:
    result = load("lebail_result.json")

    text = lebail_markdown(result)

    assert text.startswith("# s1: Le Bail extraction\n")
    assert "Le Bail extraction of P:" in text
    assert "structures.P.cell" in text
    assert "| microstrain | clean | 2 | 7.500 | 6.000 | 1.440 |  |" in text
    assert "| a (Å) | 6.0100(10) | 6 | +0.01000 |" in text
    assert "Zero shift (°): 0.012(2), from 0." in text
    assert "The microstrain test finds a microstrain: 9.0 esds from zero" in text
    assert "- Background: cosine, 4 terms." in text
    assert "## Followed reflections" in text
    for name in ("Phase 8", "Phase 10", "Chebyshev"):
        assert name not in text


def test_a_rejected_stage_has_a_row_saying_so() -> None:
    result = load("fixed_atoms_two_phases_result.json")

    text = fixed_atoms_markdown(result)

    assert "| overall Uiso | rejected | — | — | — | — |" in text
    assert "sanity check: M1 Uiso -0.004" in text
    assert "The final model is that of stage size; overall Uiso rejected" in text


def test_two_phases_show_the_weight_fractions() -> None:
    text = fixed_atoms_markdown(load("fixed_atoms_two_phases_result.json"))

    assert "## Phase fractions" in text
    assert "| P | 0.620(10) | 0.62(2) |" in text
    assert "| Q | 0.380(10) | 0.38(2) |" in text
    assert "Phase Q:" in text


def test_structure_markdown_coordinates() -> None:
    result = load("coordinates_result.json")

    text = structure_markdown(
        result, plan=plan(), earlier={"lebail": load("lebail_result.json")}
    )

    assert text.startswith("# s1: Rietveld refinement, coordinates\n")
    # The kinds named from the plan, none of them A, B or O.
    assert "Phase P has sites of the kinds M (2 sites) and X (1 site)" in text
    assert "refined kind by kind, M and X" in text
    assert "with M1 (M) holding the origin along z" in text
    assert "from the start model" in text
    assert "M1 occupancy 0.620, esd 0.600 more than half its range" in text
    # Bonds against the plan's limits, N1-X1 flagged under the M limits.
    assert "| N1 | M | X1 | 2 | 1.520(3) | 2.4 to 3 | **outside** |" in text
    assert "to F, with the limits" in text
    # Misfits classed along the polar axis, z, by l.
    assert "By the index along the polar axis, z (l)" in text
    assert "| lebail | coordinates |" in text


def test_structure_markdown_occupancies() -> None:
    text = structure_markdown(load("occupancies_result.json"), plan=plan())

    assert text.startswith("# s1: Rietveld refinement, occupancies\n")
    assert "Sr and Ba traded between the M sites M1 and N1" in text
    assert (
        "| M site occupancies | clean | P | Sr | M1, N1 | M1 0.6000, Sr2 0.0000 | "
        "1.20000 | 1.20000 |"
    ) in text


def test_coordinate_shift_in_angstroms() -> None:
    cell = {
        "length_a": 6.0,
        "length_b": 6.0,
        "length_c": 4.0,
        "angle_alpha": 90.0,
        "angle_beta": 90.0,
        "angle_gamma": 90.0,
    }
    moved = {"xyz": [0.0, 0.0, 0.52], "xyz_esd": [None, None, 0.002]}

    distance, esd = coordinate_shift(moved, {"xyz": [0.0, 0.0, 1.5]}, cell)

    # 0.02 of c, a whole cell aside, with the z esd of 0.002 as 0.008 A.
    assert distance == pytest.approx(0.08)
    assert esd == pytest.approx(0.008)


def test_reflection_misfits_without_a_polar_axis() -> None:
    misfits = reflection_misfits(DATA / "fit_reflections_P.csv", [20.5, 59.5])

    assert misfits["by_class"] == {}
    assert [row["hkl"] for row in misfits["worst"]][:2] == ["110", "001"]
    assert sum(len(group) for group in misfits["by_angle"].values()) == 4


# run_mode's figure and write up


def test_run_mode_writes_the_figure_and_the_write_up(tmp_path, monkeypatch) -> None:
    project = fake_project(tmp_path)
    monkeypatch.setattr(pipeline, "run_job", FakeGsas2())

    outcome = run_mode(project, "chain", "lebail")

    folder = project.root / "results" / "lebail" / "chain"
    assert (folder / "chain_lebail.png").is_file()
    assert (folder / "chain_lebail.pdf").is_file()
    markdown = (folder / "lebail.md").read_text(encoding="utf-8")
    assert markdown.startswith("# chain: Le Bail extraction")
    assert outcome.paths["markdown"] == str(folder / "lebail.md")
    assert [Path(path).suffix for path in outcome.paths["figures"]] == [".png", ".pdf"]


def test_a_write_up_that_fails_keeps_the_result(tmp_path, monkeypatch) -> None:
    project = fake_project(tmp_path)
    monkeypatch.setattr(pipeline, "run_job", FakeGsas2())

    def broken(*arguments):
        raise RuntimeError("the write up broke")

    monkeypatch.setitem(writeup.MARKDOWN, "lebail", broken)

    (outcome,) = run_sequence(project, "chain", ["lebail", "fixed_atoms"])

    folder = project.root / "results" / "lebail" / "chain"
    assert outcome.error == "RuntimeError: the write up broke"
    assert json.loads((folder / "chain_lebail_result.json").read_text("utf-8"))[
        "stages"
    ]
    failure = (folder / "failure.md").read_text(encoding="utf-8")
    assert "The result was written to" in failure
    assert "RuntimeError: the write up broke" in failure
    assert (folder / "summary.md").is_file()


def test_a_negative_microstrain_is_not_physical() -> None:
    result = load("lebail_result.json")
    for stage in result["stages"]:
        if stage["name"] == "microstrain":
            stage["values"]["phases"][0]["mustrain"] = {
                "value": -180.0,
                "esd": 15.0,
                "held": False,
            }

    text = lebail_markdown(result)

    assert "The microstrain test gives a negative microstrain, 12.0 esds" in text
    assert "which is not physical" in text


def test_a_held_bond_prints_its_distance_to_four_places() -> None:
    result = load("coordinates_result.json")
    result["bonds"]["P"][0]["esd"] = None

    text = structure_markdown(result, plan=plan())

    assert "| M1 | M | X1 | 4 | 2.7100 (held) | 2.4 to 3 |  |" in text
