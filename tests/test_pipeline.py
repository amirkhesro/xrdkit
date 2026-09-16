"""Tests for xrdkit.pipeline: the stage builders and start_from_result, on
synthetic plans and result files; no GSAS-II run."""

import json
import math
import re
from pathlib import Path

import pytest

import xrdkit
from xrdkit import gsas2_driver as driver
from xrdkit.library import load_entry
from xrdkit.pipeline import (
    CYCLES,
    LE_BAIL_CYCLES,
    MODES,
    PASS_TOLERANCE,
    START_BROADENING,
    PipelineError,
    StartPoint,
    coordinates_stages,
    fixed_atoms_stages,
    lebail_stages,
    occupancy_stages,
    start_from_result,
)
from xrdkit.structure import site_setup

BACKGROUND = {"function": "chebyschev-1", "terms": 8}
DRIVER_BACKGROUND = {"type": "chebyschev-1", "terms": 8}


def flags(stages: list[dict]) -> list[dict]:
    """The flags in force at each stage, as the driver accumulates them."""
    return [stage["flags"] for stage in driver.accumulate_stages(stages)]


def names(stages: list[dict]) -> list[str]:
    return [stage["name"] for stage in stages]


def test_constants() -> None:
    assert MODES == ("lebail", "fixed_atoms", "coordinates", "occupancies")
    assert (CYCLES, LE_BAIL_CYCLES, PASS_TOLERANCE) == (10, 10, 0.1)
    assert START_BROADENING == {"size": 1.0, "mustrain": 0.0, "lgmix": 1.0}
    assert issubclass(PipelineError, ValueError)


def test_pipeline_is_exported() -> None:
    for name in (
        "PipelineError",
        "StartPoint",
        "coordinates_stages",
        "fixed_atoms_stages",
        "lebail_stages",
        "occupancy_stages",
        "start_from_result",
    ):
        assert name in xrdkit.__all__
        assert getattr(xrdkit, name) is getattr(xrdkit.pipeline, name)


# Le Bail


def test_lebail_stages_one_phase() -> None:
    stages = lebail_stages(BACKGROUND)

    assert names(stages) == [
        "background and scale",
        "zero",
        "cell",
        "size",
        "microstrain",
    ]
    first, zero, cell, size, strain = flags(stages)
    assert first == {
        **driver.empty_flags(),
        "background": DRIVER_BACKGROUND,
        "scale": True,
        "le_bail": ["*"],
    }
    assert zero == {**first, "zero": True}
    assert cell == {**zero, "cell": ["*"]}
    assert size == {**cell, "size": ["*"]}
    assert strain == {**size, "mustrain": ["*"]}
    # Nothing but the flags above: no instrument parameter, no atom.
    assert strain["instrument"] == [] and strain["phase_fractions"] == []


def test_lebail_stages_two_phases_free_the_fractions_not_the_scale() -> None:
    first, *_, last = flags(lebail_stages(BACKGROUND, phases=2))

    assert first["scale"] is False and first["phase_fractions"] == ["*"]
    assert last["scale"] is False and last["phase_fractions"] == ["*"]
    assert first["le_bail"] == ["*"]


def test_lebail_stages_displacement_and_no_microstrain_test() -> None:
    stages = lebail_stages(BACKGROUND, displacement=True, mustrain_test=False)

    assert names(stages) == ["background and scale", "displacement", "cell", "size"]
    last = flags(stages)[-1]
    assert (last["zero"], last["displacement"]) == (False, True)
    assert last["mustrain"] == []


@pytest.mark.parametrize(
    ("background", "phases", "message"),
    [
        (
            {"type": "chebyschev-1", "terms": 6},
            1,
            "background must be {function, terms}",
        ),
        (
            {"function": "chebyschev-1", "terms": 0},
            1,
            "background terms must be at least 1",
        ),
        (BACKGROUND, 0, "phases must be a whole number of at least 1, not 0"),
        (BACKGROUND, True, "phases must be a whole number of at least 1, not True"),
    ],
)
def test_lebail_stages_rejects(background, phases, message) -> None:
    with pytest.raises(PipelineError, match=re.escape(message)):
        lebail_stages(background, phases=phases)


# Fixed atoms


def test_fixed_atoms_stages_one_phase() -> None:
    stages = fixed_atoms_stages(BACKGROUND)

    assert names(stages) == [
        "scale and background",
        "zero and cell",
        "size",
        "overall Uiso",
    ]
    first, cell, size, uiso = flags(stages)
    assert first == {
        **driver.empty_flags(),
        "background": DRIVER_BACKGROUND,
        "scale": True,
    }
    assert cell == {**first, "zero": True, "cell": ["*"]}
    assert size == {**cell, "size": ["*"]}
    assert uiso == {**size, "overall_uiso": ["*"]}
    assert uiso["le_bail"] == [] and uiso["mustrain"] == []


def test_fixed_atoms_stages_two_phases_displacement_and_microstrain() -> None:
    stages = fixed_atoms_stages(BACKGROUND, phases=2, displacement=True, mustrain=True)

    assert names(stages) == [
        "scale and background",
        "displacement and cell",
        "size and microstrain",
        "overall Uiso",
    ]
    first, cell, size, _ = flags(stages)
    assert (first["scale"], first["phase_fractions"]) == (False, ["*"])
    assert (cell["zero"], cell["displacement"], cell["cell"]) == (False, True, ["*"])
    assert (size["size"], size["mustrain"]) == (["*"], ["*"])


def test_fixed_atoms_stages_preferred_orientation() -> None:
    stages = fixed_atoms_stages(BACKGROUND, preferred_orientation=(0, 0, 1))

    assert names(stages)[-1] == "preferred orientation"
    assert stages[-1] == {
        "name": "preferred orientation",
        "preferred_orientation": {"*": [0, 0, 1]},
    }
    before, after = flags(stages)[-2:]
    assert before["preferred_orientation"] == {}
    assert after == {**before, "preferred_orientation": {"*": [0, 0, 1]}}


@pytest.mark.parametrize("hkl", [(0, 0, 0), (1, 0), (1.0, 0, 0), "001", (True, 0, 0)])
def test_fixed_atoms_stages_rejects_a_preferred_orientation_that_is_not_hkl(
    hkl,
) -> None:
    with pytest.raises(
        PipelineError, match=r"preferred_orientation must be \[h, k, l\]"
    ):
        fixed_atoms_stages(BACKGROUND, preferred_orientation=hkl)


def test_driver_preferred_orientation_flag() -> None:
    first, dropped = flags(
        [
            {"name": "on", "preferred_orientation": {"P": [1, 1, 0]}},
            {"name": "off", "preferred_orientation": {"P": None}},
        ]
    )
    assert first["preferred_orientation"] == {"P": [1, 1, 0]}
    assert dropped["preferred_orientation"] == {}
    with pytest.raises(ValueError, match=r"axis of 'P' must be \[h, k, l\]"):
        driver.accumulate_stages([{"preferred_orientation": {"P": [0, 0, 0]}}])


# Coordinates and occupancies, on a plan with kinds of its own


def plan(**changes) -> dict:
    """A plan as site_setup gives it, with the phase named, for a structure
    whose kinds are X, M and Q in that order: M1 fixes the origin along z,
    the Q site has nothing free."""
    sites = [
        {"name": "X1", "kind": "X"},
        {"name": "M1", "kind": "M"},
        {"name": "M2", "kind": "M"},
        {"name": "Q1", "kind": "Q"},
        {"name": "X2", "kind": "X"},
    ]
    kinds = {}
    for site in sites:
        kinds.setdefault(site["kind"], []).append(site)
    base = {
        "phase": "P",
        "sites": sites,
        "kinds": kinds,
        "uiso_groups": [["M1", "M2"], ["X1", "X2"], ["Q1"]],
        "group_names": ["M", "X", "Q"],
        "origin": sites[1],
        "origin_axis": "z",
        "polar_axis": "z",
        "coordinates": {
            "X": {"X1": "xz", "X2": "all"},
            "M": {"M2": "z"},
            "Q": {"Q1": ""},
        },
        "exchange": {"elements": ["Sr", "Ba"], "sites": ["M1", "M2"]},
    }
    return {**base, **changes}


def test_coordinates_stages_follow_the_kinds_in_order() -> None:
    stages = coordinates_stages(plan(), BACKGROUND)

    # Q has nothing to free, so no stage; X comes first, as in the plan.
    assert names(stages) == ["profile", "Uiso groups", "X sites", "M sites"]
    profile, groups, x_sites, m_sites = flags(stages)
    assert profile == {
        **driver.empty_flags(),
        "background": DRIVER_BACKGROUND,
        "scale": True,
        "zero": True,
        "cell": ["*"],
        "size": ["*"],
    }
    assert groups == {
        **profile,
        "uiso_groups": {"P": [["M1", "M2"], ["X1", "X2"], ["Q1"]]},
    }
    assert x_sites["coordinates"] == {"P": {"X1": "xz", "X2": "all"}}
    assert x_sites["origin"] == {}
    # The origin is held from the stage of its own kind.
    assert m_sites["coordinates"] == {"P": {"X1": "xz", "X2": "all", "M2": "z"}}
    assert m_sites["origin"] == {"P": {"site": "M1", "axis": "z"}}
    assert m_sites["occupancies"] == [] and m_sites["overall_uiso"] == []


def test_coordinates_stages_origin_kind_without_a_stage() -> None:
    changed = plan(
        coordinates={"X": {"X1": "xz"}, "M": {}, "Q": {"Q1": "z"}},
        origin={"name": "M1", "kind": "M"},
    )

    stages = coordinates_stages(changed, BACKGROUND)

    assert names(stages) == ["profile", "Uiso groups", "X sites", "Q sites"]
    # M has no stage, so the origin is held from the first coordinate stage.
    assert stages[2]["origin"] == {"P": {"site": "M1", "axis": "z"}}
    assert "origin" not in stages[3]


def test_coordinates_stages_two_phases_displacement_microstrain() -> None:
    profile = flags(
        coordinates_stages(
            plan(), BACKGROUND, phases=2, displacement=True, mustrain=True
        )
    )[0]

    assert (profile["scale"], profile["phase_fractions"]) == (False, ["*"])
    assert (profile["zero"], profile["displacement"]) == (False, True)
    assert profile["mustrain"] == ["*"]


def test_coordinates_stages_one_phase_never_flags_fractions() -> None:
    for stage in coordinates_stages(plan(), BACKGROUND):
        assert "phase_fractions" not in stage


def test_coordinates_stages_refuse_a_polar_structure_with_no_origin() -> None:
    with pytest.raises(
        PipelineError,
        match=(
            r"^no origin is fixed, but the structure is polar along z and the X "
            r"sites X1, X2 would be freed along it"
        ),
    ):
        coordinates_stages(plan(origin=None, origin_axis=None), BACKGROUND)


def test_coordinates_stages_without_origin_off_the_polar_axis() -> None:
    changed = plan(
        origin=None,
        origin_axis=None,
        coordinates={"X": {"X1": "x"}, "M": {"M1": "xy"}, "Q": {}},
    )

    stages = coordinates_stages(changed, BACKGROUND)

    assert names(stages) == ["profile", "Uiso groups", "X sites", "M sites"]
    assert all("origin" not in stage for stage in stages)
    # A structure that is not polar needs no origin, whatever is freed.
    stages = coordinates_stages(plan(origin=None, polar_axis=None), BACKGROUND)
    assert all("origin" not in stage for stage in stages)


# A polar library entry that fixes no origin, whose F site its symmetry
# fixes and whose X sites are free only across the polar axis.
POLAR_ENTRY = """
[entry]
name = "test/polar"
family = "test"
crystal_system = "tetragonal"
space_group = "P4mm"
cell_parameters = ["a", "c"]
z = 1
reference = "none"
polar_axis = "c"

[[sites]]
label = "F1"
kind = "F"
wyckoff = "1a"
free = []

[[sites]]
label = "X1"
kind = "X"
wyckoff = "4d"
free = ["x", "y"]
"""


def polar_plan(tmp_path, monkeypatch) -> dict:
    """The plan site_setup gives for the entry above, the entry read from
    under tmp_path."""
    import xrdkit.structure

    folder = tmp_path / "test"
    folder.mkdir()
    (folder / "polar.toml").write_text(POLAR_ENTRY, encoding="utf-8")
    monkeypatch.setattr(
        xrdkit.structure,
        "load_entry",
        lambda name: load_entry(name, root=tmp_path),
    )
    structure = {
        "library": "test/polar",
        "sites": [
            {
                "name": "M1",
                "label": "F1",
                "atoms": {"M1": "Mo"},
                "wyckoff": "1a",
                "kind": "F",
            },
            {
                "name": "N1",
                "label": "X1",
                "atoms": {"N1": "N"},
                "wyckoff": "4d",
                "kind": "X",
            },
        ],
        "uiso_groups": [{"name": "all", "sites": ["M1", "N1"]}],
    }
    atoms = [
        {"label": "M1", "type": "Mo", "xyz": [0.0, 0.0, 0.0]},
        {"label": "N1", "type": "N", "xyz": [0.2, 0.3, 0.0]},
    ]
    return {**site_setup(structure, atoms), "phase": "P"}


def test_a_site_its_symmetry_fixes_has_no_coordinate_stage(
    tmp_path, monkeypatch
) -> None:
    plan = polar_plan(tmp_path, monkeypatch)

    # The entry's own free lists, "" for the fixed site, not "all".
    assert [site["free"] for site in plan["sites"]] == ["", "xy"]
    assert plan["polar_axis"] == "z" and plan["origin"] is None

    stages = coordinates_stages(plan, BACKGROUND)

    # No F stage, and no refusal: nothing freed lies along the polar axis.
    assert names(stages) == ["profile", "Uiso groups", "X sites"]
    assert stages[-1]["coordinates"] == {"P": {"N1": "xy"}}
    assert all("origin" not in stage for stage in stages)


def test_coordinates_stages_need_the_phase_name() -> None:
    changed = plan()
    del changed["phase"]
    with pytest.raises(PipelineError, match="the plan names no phase"):
        coordinates_stages(changed, BACKGROUND)


def test_occupancy_stages() -> None:
    stages = occupancy_stages(plan(), BACKGROUND)

    assert names(stages) == ["profile and Uiso", "M site occupancies"]
    profile, exchange = flags(stages)
    assert profile == {
        **driver.empty_flags(),
        "background": DRIVER_BACKGROUND,
        "scale": True,
        "zero": True,
        "cell": ["*"],
        "size": ["*"],
        "uiso_groups": {"P": [["M1", "M2"], ["X1", "X2"], ["Q1"]]},
    }
    # The coordinates are held throughout.
    assert profile["coordinates"] == {} and exchange["coordinates"] == {}
    assert exchange["occupancies"] == [
        {"phase": "P", "sites": ["M1", "M2"], "elements": ["Sr", "Ba"]}
    ]


def test_occupancy_stages_two_phases_and_groups_of_one_kind() -> None:
    changed = plan(
        exchange=[
            {"elements": ["Sr", "Ba"], "sites": ["M1", "M2"]},
            {"elements": ["La", "Ca"], "sites": ["M2", "M1"]},
            {"elements": ["F", "O"], "sites": ["X1", "X2"]},
        ]
    )

    stages = occupancy_stages(changed, BACKGROUND, phases=2)

    assert names(stages) == [
        "profile and Uiso",
        "M site occupancies (Sr, Ba)",
        "M site occupancies (La, Ca)",
        "X site occupancies",
    ]
    assert (stages[0]["scale"], stages[0]["phase_fractions"]) == (False, True)
    assert len(flags(stages)[-1]["occupancies"]) == 3


def test_occupancy_stages_without_exchange_and_across_kinds() -> None:
    assert names(occupancy_stages(plan(exchange=None), BACKGROUND)) == [
        "profile and Uiso"
    ]
    with pytest.raises(
        PipelineError,
        match=r"^the exchange sites M1, X1 are of the kinds M, X, not of one kind$",
    ):
        occupancy_stages(
            plan(exchange={"elements": ["Sr", "Ba"], "sites": ["M1", "X1"]}),
            BACKGROUND,
        )


# start_from_result


def reciprocal_terms(a: float, c: float) -> tuple[float, float]:
    """A0 and A2 of a tetragonal cell."""
    return 1.0 / a**2, 1.0 / c**2


def parameter(value: float) -> dict:
    return {"value": value, "esd": 0.001}


def result(**changes) -> dict:
    """A driver result of four stages, the second rejected: a cell of
    a = 5.0, c = 7.0 after "cell", refined on to 5.1 and 7.1 by "strain"."""
    a0, a2 = reciprocal_terms(5.0, 7.0)
    atoms_then = [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.1]}]
    atoms_final = [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.2]}]
    base_params = {
        ":0:Scale": parameter(2.0),
        ":0:Back;0": parameter(100.0),
        ":0:Back;1": parameter(-3.0),
    }
    stages = [
        {"name": "scale", "status": "clean", "parameters": base_params, "atoms": {}},
        {"name": "wild", "status": "rejected", "parameters": {}, "atoms": {}},
        {
            "name": "cell",
            "status": "unsettled",
            "parameters": {
                **base_params,
                ":0:Zero": parameter(0.02),
                "0::A0": parameter(a0),
                "0::A2": parameter(a2),
                "0:0:Size;i": parameter(0.3),
            },
            "atoms": {"P": atoms_then},
        },
        {
            "name": "strain",
            "status": "clean",
            "parameters": {
                ":0:Scale": parameter(2.5),
                ":0:Back;0": parameter(110.0),
                ":0:Back;1": parameter(-4.0),
                ":0:Zero": parameter(0.03),
                "0::A0": parameter(reciprocal_terms(5.1, 7.1)[0]),
                "0::A2": parameter(reciprocal_terms(5.1, 7.1)[1]),
                "0:0:Size;i": parameter(0.35),
                "0:0:Mustrain;i": parameter(900.0),
            },
            "atoms": {"P": atoms_final},
        },
    ]
    base = {
        "limits": [15.0, 100.0],
        "stages": stages,
        "final_from": "strain",
        "final": {
            "instrument": {"Zero": {"value": 0.03, "esd": 0.001}},
            "sample": {
                "Scale": {"value": 2.5, "esd": 0.01},
                "Shift": {"value": 0.0, "esd": None},
            },
            "background": {"type": "chebyschev-1", "coefficients": [110.0, -4.0]},
            "phases": [
                {
                    "name": "P",
                    "cell": {
                        "length_a": 5.1,
                        "length_b": 5.1,
                        "length_c": 7.1,
                        "angle_alpha": 90.0,
                        "angle_beta": 90.0,
                        "angle_gamma": 90.0,
                        "volume": 5.1 * 5.1 * 7.1,
                    },
                    "size": {"value": 0.35},
                    "mustrain": {"value": 900.0},
                    "atoms": atoms_final,
                }
            ],
        },
    }
    return {**base, **changes}


def write(tmp_path: Path, data: object, name: str = "mode_result.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_start_from_result_last_accepted_stage(tmp_path) -> None:
    start = start_from_result(write(tmp_path, result()))

    assert isinstance(start, StartPoint)
    assert start.stage == "strain"
    assert start.cell == {
        "a": 5.1,
        "b": 5.1,
        "c": 7.1,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    assert (start.zero, start.displacement, start.size) == (0.03, 0.0, 0.35)
    assert (start.microstrain, start.scale) == (900.0, 2.5)
    assert start.background == {
        "function": "chebyschev-1",
        "terms": 2,
        "coefficients": [110.0, -4.0],
    }
    assert start.limits == (15.0, 100.0)
    assert start.atoms == [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.2]}]
    assert start.start_model is None
    with pytest.raises(AttributeError):
        start.zero = 1.0


def test_start_from_result_an_earlier_stage_of_an_old_result(tmp_path) -> None:
    # Written before every value was recorded at every stage, and with no
    # start_model; nothing after "cell" refines what "cell" held.
    data = result()
    del data["stages"][3]["parameters"]["0:0:Mustrain;i"]
    data["final"]["phases"][0]["mustrain"]["value"] = 0.0

    start = start_from_result(write(tmp_path, data), stage="cell")

    assert start.stage == "cell"
    # a and c as that stage refined them, b following a by symmetry.
    for key, value in {"a": 5.0, "b": 5.0, "c": 7.0}.items():
        assert start.cell[key] == pytest.approx(value)
    for key in ("alpha", "beta", "gamma"):
        assert start.cell[key] == pytest.approx(90.0)
    assert (start.zero, start.size, start.scale) == (0.02, 0.3, 2.0)
    assert start.background["coefficients"] == [100.0, -3.0]
    assert start.atoms == [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.1]}]
    assert start.microstrain == 0.0
    # Neither this stage nor a later one refined the displacement: the final.
    assert start.displacement == 0.0
    assert start.start_model is None


def test_start_from_result_hexagonal_ties(tmp_path) -> None:
    data = result()
    del data["stages"][3]["parameters"]["0:0:Mustrain;i"]
    final_cell = data["final"]["phases"][0]["cell"]
    final_cell.update(length_a=4.0, length_b=4.0, length_c=9.0, angle_gamma=120.0)
    a0 = 4.0 / (3.0 * 3.9**2)
    for stage in data["stages"][2:]:
        stage["parameters"]["0::A0"] = parameter(a0)
        stage["parameters"]["0::A2"] = parameter(1.0 / 8.8**2)
    data["stages"][3]["parameters"]["0::A0"] = parameter(4.0 / (3.0 * 16.0))
    data["stages"][3]["parameters"]["0::A2"] = parameter(1.0 / 81.0)

    start = start_from_result(write(tmp_path, data), stage="cell")

    # A1 and A3 follow A0 in a hexagonal cell.
    assert start.cell["a"] == pytest.approx(3.9)
    assert start.cell["b"] == pytest.approx(3.9)
    assert start.cell["c"] == pytest.approx(8.8)
    assert start.cell["gamma"] == pytest.approx(120.0)
    assert math.isclose(start.cell["alpha"], 90.0, abs_tol=1e-9)


def recorded(value: float, held: bool = False) -> dict:
    """A value as the driver records it at a stage."""
    return {"value": value, "esd": None if held else 0.001, "held": held}


def stage_values(
    cell: tuple[float, float],
    zero: float,
    size: float,
    microstrain: float | None,
    coefficients: list[float],
    phases: tuple[str, ...] = ("P",),
    fractions: tuple[float, ...] = (1.0,),
) -> dict:
    """A stage's values as the driver records them: a held microstrain when
    ``microstrain`` is None, at 0; the phases alike but for name, fraction
    and a cell 1 A longer in a for each after the first."""
    return {
        "instrument": {"Zero": recorded(zero), "U": recorded(2.0, held=True)},
        "sample": {
            "Scale": recorded(1.0, held=len(phases) > 1),
            "Shift": recorded(0.0, held=True),
        },
        "background": {
            "type": "chebyschev-1",
            "coefficients": [recorded(value) for value in coefficients],
        },
        "phases": [
            {
                "name": name,
                "cell": {
                    "length_a": cell[0] + index,
                    "length_b": cell[0] + index,
                    "length_c": cell[1],
                    "angle_alpha": 90.0,
                    "angle_beta": 90.0,
                    "angle_gamma": 90.0,
                    "volume": (cell[0] + index) ** 2 * cell[1],
                },
                "cell_esd": None,
                "cell_held": False,
                "phase_fraction": recorded(fraction, held=len(phases) == 1),
                "size": {"type": "isotropic", **recorded(size)},
                "mustrain": {
                    "type": "isotropic",
                    **(
                        recorded(0.0, held=True)
                        if microstrain is None
                        else recorded(microstrain)
                    ),
                },
            }
            for index, (name, fraction) in enumerate(zip(phases, fractions))
        ],
    }


def test_start_from_result_every_value_at_every_stage(tmp_path) -> None:
    data = result(start_model={"P": [{"label": "M1", "xyz": [0.0, 0.0, 0.0]}]})
    scale, _, cell, strain = data["stages"]
    scale["atoms"] = {"P": [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.0]}]}
    scale["values"] = stage_values((5.5, 7.5), 0.0, 1.0, None, [100.0, -3.0])
    cell["values"] = stage_values((5.0, 7.0), 0.02, 0.3, None, [100.0, -3.0])
    strain["values"] = stage_values((5.1, 7.1), 0.03, 0.35, 900.0, [110.0, -4.0])

    start = start_from_result(write(tmp_path, data), stage="cell")

    # The microstrain was held at "cell" and refined only by "strain", and is
    # read as held, where an old result could not give it.
    assert start.microstrain == 0.0
    assert start.cell["a"] == 5.0 and start.cell["c"] == 7.0
    assert (start.zero, start.size, start.scale) == (0.02, 0.3, 1.0)
    assert start.background["coefficients"] == [100.0, -3.0]
    assert start.displacement == 0.0
    assert start.atoms == [{"label": "M1", "type": "Sr", "xyz": [0.0, 0.0, 0.1]}]
    assert start.start_model == {"P": [{"label": "M1", "xyz": [0.0, 0.0, 0.0]}]}
    # The first stage too, which held all but the scale and background.
    first = start_from_result(write(tmp_path, data), stage="scale")
    assert (first.zero, first.size, first.cell["a"]) == (0.0, 1.0, 5.5)
    # And the last, from its own values rather than the final model.
    last = start_from_result(write(tmp_path, data))
    assert (last.stage, last.microstrain, last.zero) == ("strain", 900.0, 0.03)


def test_start_from_result_two_phases(tmp_path) -> None:
    data = result()
    atoms_q = [{"label": "Q1", "type": "Ti", "xyz": [0.5, 0.5, 0.5]}]
    for stage in data["stages"]:
        stage["atoms"] = {**stage["atoms"], "Q": atoms_q}
    scale, _, cell, strain = data["stages"]
    for record, values in (
        (scale, ((5.5, 7.5), 0.0, 1.0, None, [100.0, -3.0])),
        (cell, ((5.0, 7.0), 0.02, 0.3, None, [100.0, -3.0])),
        (strain, ((5.1, 7.1), 0.03, 0.35, 900.0, [110.0, -4.0])),
    ):
        record["atoms"] = {"P": [{"label": "M1"}], "Q": atoms_q}
        record["values"] = stage_values(
            *values, phases=("P", "Q"), fractions=(0.7, 0.3)
        )

    start = start_from_result(write(tmp_path, data), stage="cell")

    assert [phase.name for phase in start.phases] == ["P", "Q"]
    p, q = start.phases
    assert (p.fraction, q.fraction) == (0.7, 0.3)
    assert (p.cell["a"], q.cell["a"]) == (5.0, 6.0)
    assert q.atoms == atoms_q and p.atoms == [{"label": "M1"}]
    assert (q.size, q.microstrain) == (0.3, 0.0)
    # The histogram values once, and the first phase's as the single fields.
    assert start.scale == 1.0 and start.zero == 0.02
    assert (start.cell, start.size, start.atoms) == (p.cell, p.size, p.atoms)
    with pytest.raises(AttributeError):
        start.phases = ()


def test_start_from_result_failures(tmp_path) -> None:
    missing = tmp_path / "missing.json"
    with pytest.raises(
        PipelineError, match=rf"^{re.escape(str(missing))}: no such result file$"
    ):
        start_from_result(missing)

    garbled = tmp_path / "garbled.json"
    garbled.write_text("{", encoding="utf-8")
    with pytest.raises(PipelineError, match="garbled.json: not a result JSON"):
        start_from_result(garbled)

    path = write(tmp_path, result())
    with pytest.raises(
        PipelineError,
        match=(
            r"mode_result\.json: stage 'wild' was not accepted; the accepted stages "
            r"are scale, cell, strain$"
        ),
    ):
        start_from_result(path, stage="wild")
    with pytest.raises(PipelineError, match=r"stage 'nothing' was not accepted"):
        start_from_result(path, stage="nothing")
    # In a result written before every value was recorded, the microstrain
    # held at "cell" and refined later cannot be read.
    with pytest.raises(
        PipelineError,
        match=(
            r"mode_result\.json: stage 'cell' did not refine the microstrain, which "
            r"a later stage did, and the result, written before every value was "
            r"recorded at every stage, does not give it$"
        ),
    ):
        start_from_result(path, stage="cell")

    none_kept = result()
    for stage in none_kept["stages"]:
        stage["status"] = "rejected"
    with pytest.raises(PipelineError, match=r"none\.json: no stage was accepted$"):
        start_from_result(write(tmp_path, none_kept, "none.json"))

    for name, change in (
        ("no_final", lambda data: data.pop("final")),
        ("no_zero", lambda data: data["final"]["instrument"].pop("Zero")),
        ("no_phase", lambda data: data["final"]["phases"].clear()),
        ("no_limits", lambda data: data.pop("limits")),
    ):
        data = result()
        change(data)
        with pytest.raises(PipelineError, match=rf"{name}\.json: no ") as raised:
            start_from_result(write(tmp_path, data, f"{name}.json"))
        assert "\n" not in str(raised.value)
