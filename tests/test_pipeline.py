"""Tests for xrdkit.pipeline: the stage builders and start_from_result, on
synthetic plans and result files; no GSAS-II run."""

import copy
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import xrdkit
from xrdkit import gsas2_driver as driver
from xrdkit import pipeline
from xrdkit.broadening import Caglioti
from xrdkit.gsas2 import Gsas2Error, build_refine_job, find_gsas2, write_instprm
from xrdkit.library import load_entry
from xrdkit.pipeline import (
    CYCLES,
    LE_BAIL_CYCLES,
    MODES,
    PASS_TOLERANCE,
    START_BROADENING,
    Options,
    PipelineError,
    StartPoint,
    coordinates_stages,
    fixed_atoms_stages,
    lebail_stages,
    occupancy_stages,
    resolve_inputs,
    run_mode,
    run_sequence,
    start_from_result,
)
from xrdkit.project import load_project
from xrdkit.structure import composition_edits, site_setup
from xrdkit.symmetry import is_absent, space_group_operations

BACKGROUND = {"function": "chebyschev-1", "terms": 8}
# The instrument the GSAS-II tests use.
CAGLIOTI = Caglioti(
    u=1.047e-02,
    v=-1.184e-02,
    w=8.359e-03,
    esd_u=1.714e-03,
    esd_v=2.049e-03,
    esd_w=5.663e-04,
    n_peaks=14,
    rms=0.00255,
    covariance=np.zeros((3, 3)),
)
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


# Running the modes, with GSAS-II played by a fake


TOY_CIF = """data_toy
_cell_length_a  6.0
_cell_length_b  6.0
_cell_length_c  4.0
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 4 b m'
loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_occupancy
   _atom_site_U_iso_or_equiv
Sr1  Sr  0.0  0.0  0.5   0.6  0.01
Ba2  Ba  0.2  0.7  0.5   0.6  0.01
Nb1  Nb  0.5  0.0  0.0   1.0  0.01
O1   O   0.3  0.1  0.1   1.0  0.01
"""

# The toy's atoms as the driver's create action reports them.
TOY_ATOMS = [
    {
        "label": "Sr1",
        "type": "Sr",
        "xyz": [0.0, 0.0, 0.5],
        "occupancy": 0.6,
        "multiplicity": 2,
        "adp": "I",
        "uiso": 0.01,
    },
    {
        "label": "Ba2",
        "type": "Ba",
        "xyz": [0.2, 0.7, 0.5],
        "occupancy": 0.6,
        "multiplicity": 4,
        "adp": "I",
        "uiso": 0.01,
    },
    {
        "label": "Nb1",
        "type": "Nb",
        "xyz": [0.5, 0.0, 0.0],
        "occupancy": 1.0,
        "multiplicity": 2,
        "adp": "I",
        "uiso": 0.01,
    },
    {
        "label": "O1",
        "type": "O",
        "xyz": [0.3, 0.1, 0.1],
        "occupancy": 1.0,
        "multiplicity": 8,
        "adp": "I",
        "uiso": 0.01,
    },
]

INSTPRM = """#GSAS-II instrument parameter file
Type:PXC
Lam1:1.540598
Lam2:1.544426
I(L2)/I(L1):0.5
Zero:0.0
Polariz.:0.7
Azimuth:0.0
U:2.0
V:-1.0
W:5.0
X:0.0
Y:1.0
Z:0.0
SH/L:0.002
"""

XRDML_SCAN = """<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.0">
  <sample><id>toy</id></sample>
  <xrdMeasurement>
    <usedWavelength>
      <kAlpha1>1.5405980</kAlpha1>
      <kAlpha2>1.5444260</kAlpha2>
      <ratioKAlpha2KAlpha1>0.5</ratioKAlpha2KAlpha1>
    </usedWavelength>
    <scan>
      <dataPoints>
        <positions axis="2Theta" unit="deg">
          <startPosition>{start}</startPosition>
          <endPosition>{end}</endPosition>
        </positions>
        <commonCountingTime>1.0</commonCountingTime>
        <counts>{counts}</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
"""

PROJECT = """
[project]
name = "fake"
version = 1

[refine]
two_theta = [20.5, 59.5]
background = { function = "chebyschev-1", terms = 3 }

[instruments.lab]
wavelength = [1.540598, 1.544426]
ka2 = true
radius = 240.0
instprm = "lab.instprm"

[structures.bronze]
library = "ttb/P4bm"
cif = "toy.cif"
composition = "Sr0.24Ba0.48Nb0.4O1.6"
cell = { a = 6.0, c = 4.0 }

[structures.toy]
cif = "toy.cif"
composition = "Sr0.24Ba0.48Nb0.4O1.6"
z = 5
exchange = [["Sr", "Ba"]]
origin = { site = "Nb1", axis = "z" }
atoms = [
    { atoms = { Sr1 = "Sr" }, wyckoff = "2a", kind = "A" },
    { atoms = { Ba2 = "Ba" }, wyckoff = "4c", kind = "A" },
    { atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "B" },
    { atoms = { O1 = "O" }, wyckoff = "8d", kind = "O" },
]

[structures.second]
cif = "toy.cif"
composition = "Sr0.24Ba0.48Nb0.4O1.6"
z = 5

[structures.lanthanum]
library = "ttb/P4bm"
cif = "toy.cif"
composition = "Sr0.2Ba0.48La0.04Nb0.4O1.6"
cell = { a = 6.0, c = 4.0 }

[samples."x0.10.powder"]
file = "toy.xrdml"
instrument = "lab"
structures = ["bronze"]
form = "powder"

[samples.pellet]
file = "toy.xrdml"
instrument = "lab"
structures = ["bronze"]
form = "pellet"

[samples.chain]
file = "toy.xrdml"
instrument = "lab"
structures = ["toy"]
form = "powder"

[samples.two]
file = "toy.xrdml"
instrument = "lab"
structures = ["toy", "second"]
form = "powder"

[samples.unplaced]
file = "toy.xrdml"
instrument = "lab"
structures = ["lanthanum"]
form = "powder"
"""


def fake_project(tmp_path: Path, extra: str = ""):
    root = tmp_path / "project"
    root.mkdir()
    (root / "toy.cif").write_text(TOY_CIF, encoding="utf-8")
    (root / "lab.instprm").write_text(INSTPRM, encoding="utf-8")
    counts = " ".join(["100"] * 401)
    (root / "toy.xrdml").write_text(
        XRDML_SCAN.format(start=20.0, end=60.0, counts=counts), encoding="utf-8"
    )
    (root / "xrdkit.toml").write_text(PROJECT + extra, encoding="utf-8")
    return load_project(root)


def _edited(atoms: list[dict], edits: list[dict] | None) -> list[dict]:
    """``atoms`` with the driver's atom edits applied."""
    atoms = [dict(atom) for atom in atoms]
    by_label = {atom["label"]: atom for atom in atoms}
    for edit in edits or ():
        if edit["label"] in by_label:
            by_label[edit["label"]]["occupancy"] = edit["occupancy"]
            for key in ("xyz", "uiso"):
                if key in edit:
                    by_label[edit["label"]][key] = edit[key]
        else:
            host = by_label[edit["copy"]]
            atom = {
                **host,
                "label": edit["label"],
                "type": edit["type"],
                "occupancy": edit["occupancy"],
            }
            atoms.append(atom)
            by_label[atom["label"]] = atom
    return atoms


FAKE_HISTOGRAM = "two_theta,observed,calculated,background,difference\n" + "".join(
    f"{20.5 + 0.5 * i:.2f},{100 + 50 * (i % 7 == 0)},{98 + 50 * (i % 7 == 0)},90,2\n"
    for i in range(78)
)
FAKE_REFLECTIONS = (
    "h,k,l,multiplicity,d,two_theta,fwhm,f_obs_squared,f_calc_squared,i_corr\n"
    "1,1,0,4,4.24,20.92,0.1,80,100,1.0\n"
    "0,0,1,2,4.0,22.20,0.1,210,200,1.0\n"
    "2,0,1,8,2.4,37.4,0.1,50,60,1.0\n"
)


class FakeGsas2:
    """Plays run_job: a create job reports the toy's atoms, edits applied; a
    refine job writes a result of its stages, every one clean unless told
    otherwise, with exports and a log, and keeps every job it was given."""

    def __init__(
        self, status="clean", completed=True, final=True, change=None, raises=None
    ):
        self.jobs = []
        self.status = status
        self.completed = completed
        self.final = final
        self.change = change
        self.raises = raises

    def __call__(self, job, workdir, install=None):
        self.jobs.append(copy.deepcopy(job))
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        phases = job["phases"]
        if job["action"] == "create":
            return {
                "phases": [
                    {
                        "name": phase["name"],
                        "cell": {
                            "length_a": 6.0,
                            "length_b": 6.0,
                            "length_c": 4.0,
                            "angle_alpha": 90.0,
                            "angle_beta": 90.0,
                            "angle_gamma": 90.0,
                        },
                        "atoms": _edited(TOY_ATOMS, phase.get("atoms")),
                    }
                    for phase in phases
                ]
            }
        (workdir / "refine.log").write_text("fake GSAS-II log line\n", encoding="utf-8")
        if self.raises is not None:
            raise self.raises
        terms = job["stages"][0]["background"]["terms"]
        atoms = {
            phase["name"]: _edited(TOY_ATOMS, phase.get("atoms")) for phase in phases
        }
        cells = {
            phase["name"]: dict(
                zip(
                    (
                        "length_a",
                        "length_b",
                        "length_c",
                        "angle_alpha",
                        "angle_beta",
                        "angle_gamma",
                    ),
                    phase.get("cell") or [6.0, 6.0, 4.0, 90.0, 90.0, 90.0],
                )
            )
            for phase in phases
        }

        def value(number, held=False):
            return {"value": number, "esd": None if held else 0.001, "held": held}

        values = {
            "instrument": {"Zero": value(0.01)},
            "sample": {"Scale": value(2.0), "Shift": value(0.02)},
            "background": {
                "type": "chebyschev-1",
                "coefficients": [value(10.0)] * terms,
            },
            "phases": [
                {
                    "name": name,
                    "cell": cells[name],
                    "cell_esd": None,
                    "cell_held": False,
                    "phase_fraction": value(0.5),
                    "size": value(0.4),
                    "mustrain": value(0.0, held=True),
                }
                for name in atoms
            ],
        }
        stages = [
            {
                "name": stage["name"],
                "status": self.status,
                "parameters": {},
                "rwp": 10.0,
                "rp": 8.0,
                "gof": 1.5,
                "chi_squared": 2.25,
                "atoms": atoms,
                "values": values,
            }
            for stage in job["stages"]
        ]
        if not self.completed:
            stages[-1] = {
                "name": stages[-1]["name"],
                "status": "failed",
                "error": "RefinementError: singular matrix",
            }
        prefix = Path(job["export_prefix"])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        instprm = Path(job["instprm"]).read_text(encoding="utf-8")
        if self.change:
            instprm = instprm.replace(*self.change)
        exports = {
            "histogram": str(prefix.with_name(prefix.name + "_histogram.csv")),
            "reflections": {
                name: str(prefix.with_name(f"{prefix.name}_reflections_{name}.csv"))
                for name in atoms
            },
            "instprm": str(prefix.with_name(prefix.name + ".instprm")),
        }
        Path(exports["histogram"]).write_text(FAKE_HISTOGRAM, encoding="utf-8")
        for path in exports["reflections"].values():
            Path(path).write_text(FAKE_REFLECTIONS, encoding="utf-8")
        Path(exports["instprm"]).write_text(instprm, encoding="utf-8")
        result = {
            "limits": job["limits"],
            "stages": stages,
            "completed": self.completed,
            "rejected": [s["name"] for s in stages if s.get("status") == "rejected"],
            "undetermined": [],
            "start_model": job.get("start_model") or atoms,
            "final_from": stages[-1]["name"],
            "exports": exports,
        }
        if self.final:
            result["final"] = {
                "instrument": {"Zero": {"value": 0.01, "esd": 0.001}},
                "sample": {
                    "Scale": {"value": 2.0, "esd": 0.01},
                    "Shift": {"value": 0.02, "esd": 0.001},
                },
                "background": {"type": "chebyschev-1", "coefficients": [10.0] * terms},
                "phases": [
                    {
                        "name": name,
                        "cell": cells[name],
                        "size": {"value": 0.4},
                        "mustrain": {"value": 0.0},
                        "phase_fraction": {"value": 0.5},
                        "weight_fraction": {"value": 0.5, "esd": 0.02},
                        "atoms": atoms[name],
                    }
                    for name in atoms
                ],
            }
        Path(job["result"]).write_text(json.dumps(result), encoding="utf-8")
        return result


@pytest.fixture
def fake(monkeypatch):
    gsas2 = FakeGsas2()
    monkeypatch.setattr(pipeline, "run_job", gsas2)
    return gsas2


def refine_jobs(gsas2: FakeGsas2) -> list[dict]:
    return [job for job in gsas2.jobs if job["action"] == "refine"]


def test_run_mode_file_layout_and_a_key_with_a_dot(tmp_path, fake) -> None:
    project = fake_project(tmp_path)
    lines = []

    outcome = run_mode(project, "x0.10.powder", "lebail", reporter=lines.append)

    folder = project.root / "results" / "lebail" / "x0.10.powder"
    (job,) = refine_jobs(fake)
    assert job["gpx"] == str(folder / "x0.10.powder_lebail.gpx")
    assert job["export_prefix"] == str(folder / "x0.10.powder_lebail")
    assert job["result"] == str(folder / "x0.10.powder_lebail_result.json")
    assert (folder / "work" / "lebail" / "refine.log").is_file()
    assert (folder / "work" / "lebail" / "x0.10.powder_lebail_start.instprm").is_file()
    assert Path(outcome.paths["instprm"]) == folder / "x0.10.powder_lebail.instprm"
    assert outcome.error is None and outcome.accepted == [
        "background and scale",
        "zero",
        "cell",
        "size",
        "microstrain",
    ]
    assert outcome.residuals == {
        "rwp": 10.0,
        "rp": 8.0,
        "chi_squared": 2.25,
        "reduced_chi_squared": 2.25,
    }
    saved = json.loads((folder / "x0.10.powder_lebail_result.json").read_text("utf-8"))
    assert saved["inputs"]["sample"]["key"] == "x0.10.powder"
    assert saved["inputs"]["instrument"]["radius"] == 240.0
    assert saved["inputs"]["limits"] == [20.5, 59.5]
    assert saved["method"]["stages"][0]["name"] == "background and scale"
    assert saved["method"]["driver_version"] == driver.DRIVER_VERSION
    assert "xrdkit_version" in saved["method"] and saved["date"]
    assert lines and all(isinstance(line, str) for line in lines)

    # The Rietveld modes go to results/rietveld/<key>, from the Le Bail result.
    run_mode(project, "x0.10.powder", "fixed_atoms")
    rietveld = project.root / "results" / "rietveld" / "x0.10.powder"
    assert refine_jobs(fake)[-1]["gpx"] == str(
        rietveld / "x0.10.powder_fixed_atoms.gpx"
    )
    assert (rietveld / "x0.10.powder_fixed_atoms_result.json").is_file()

    # Or all to one folder.
    out = tmp_path / "elsewhere"
    run_mode(project, "x0.10.powder", "lebail", Options(out=out))
    assert (out / "x0.10.powder_lebail_result.json").is_file()


def test_run_mode_start_cell_order(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    run_mode(project, "x0.10.powder", "lebail")
    assert refine_jobs(fake)[-1]["phases"][0]["cell"] == [
        6.0,
        6.0,
        4.0,
        90.0,
        90.0,
        90.0,
    ]
    assert resolve_inputs(project, "x0.10.powder").start_cell_source == (
        "structures.bronze.cell"
    )

    lattice = project.root / "results" / "lattice" / "x0.10.powder"
    lattice.mkdir(parents=True)
    (lattice / "lattice_x0.10.powder.csv").write_text(
        "a,esd_a,b,esd_b,c,esd_c,alpha,esd_alpha,beta,esd_beta,gamma,esd_gamma\n"
        "6.2,0.01,6.2,0.01,4.1,0.01,90,0,90,0,90,0\n"
        "6.1,0.01,6.1,0.01,4.05,0.01,90,0,90,0,90,0\n",
        encoding="utf-8",
    )
    run_mode(project, "x0.10.powder", "lebail")
    assert refine_jobs(fake)[-1]["phases"][0]["cell"] == [
        6.1,
        6.1,
        4.05,
        90.0,
        90.0,
        90.0,
    ]
    saved = json.loads(
        (
            project.root / "results/lebail/x0.10.powder/x0.10.powder_lebail_result.json"
        ).read_text("utf-8")
    )
    assert saved["inputs"]["start_cell_source"]["bronze"].startswith("lattice results")

    cell = {"a": 5.9, "b": 5.9, "c": 3.9, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    run_mode(project, "x0.10.powder", "lebail", Options(cell=cell))
    assert refine_jobs(fake)[-1]["phases"][0]["cell"] == [
        5.9,
        5.9,
        3.9,
        90.0,
        90.0,
        90.0,
    ]
    assert resolve_inputs(
        project, "x0.10.powder", Options(cell=cell)
    ).start_cell_source == ("options")
    with pytest.raises(PipelineError, match="options.cell must give exactly a, b, c"):
        resolve_inputs(project, "x0.10.powder", Options(cell={"a": 5.9}))


def test_run_mode_displacement_by_form_and_option(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    run_mode(project, "x0.10.powder", "lebail")
    run_mode(project, "pellet", "lebail")
    run_mode(project, "x0.10.powder", "lebail", Options(displacement=True))
    run_mode(project, "pellet", "lebail", Options(displacement=False))

    second = [job["stages"][1]["name"] for job in refine_jobs(fake)]
    assert second == ["zero", "displacement", "displacement", "zero"]


def test_run_mode_two_phases_free_the_fractions(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    outcome = run_mode(project, "two", "lebail")

    (job,) = refine_jobs(fake)
    assert [phase["name"] for phase in job["phases"]] == ["toy", "second"]
    assert (job["stages"][0]["scale"], job["stages"][0]["phase_fractions"]) == (
        False,
        True,
    )
    assert outcome.residuals["weight_fractions"] == {
        "toy": {"value": 0.5, "esd": 0.02},
        "second": {"value": 0.5, "esd": 0.02},
    }


def test_run_mode_passes_the_start_model(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    for mode in ("lebail", "fixed_atoms", "coordinates", "occupancies"):
        run_mode(project, "chain", mode)

    _, fixed, coordinates, occupancies = refine_jobs(fake)
    # The structure as fixed atoms set it up, nominal occupancies in place,
    # is the reference every later mode judges against.
    set_up = fake.jobs[[j["action"] for j in fake.jobs].index("create") + 1]
    assert fixed["start_model"] == {
        "toy": _edited(TOY_ATOMS, set_up["phases"][0].get("atoms"))
    }
    assert coordinates["start_model"] == fixed["start_model"]
    assert occupancies["start_model"] == fixed["start_model"]
    # Coordinates by kind; the B site, Nb1, fixes the origin and is the only
    # B site, so B has no stage and the origin is held from the first.
    assert [stage["name"] for stage in coordinates["stages"]] == [
        "profile",
        "Uiso groups",
        "A sites",
        "O sites",
    ]
    assert coordinates["stages"][2]["origin"] == {"toy": {"site": "Nb1", "axis": "z"}}
    assert coordinates["stages"][2]["coordinates"] == {
        "toy": {"Sr1": "all", "Ba2": "all"}
    }
    assert [stage["name"] for stage in occupancies["stages"]] == [
        "profile and Uiso",
        "A site occupancies",
    ]
    added = [edit for edit in occupancies["phases"][0]["atoms"] if "copy" in edit]
    assert {(edit["label"], edit["copy"]) for edit in added} == {
        ("Ba1", "Sr1"),
        ("Sr2", "Ba2"),
    }


def test_start_model_records_the_starting_uiso_of_anisotropic_atoms(
    tmp_path, fake, monkeypatch
) -> None:
    # A CIF with anisotropic Uij carries no Uiso, but fixed_atoms makes every
    # atom isotropic before it refines. Unless the start model says what Uiso
    # each atom starts from, the undetermined check has nothing to judge a
    # refined Uiso against and never fires.
    anisotropic = [
        {**atom, "adp": "A", "uiso": None, "uij": [0.012, 0.012, 0.012, 0.0, 0.0, 0.0]}
        for atom in TOY_ATOMS
    ]
    monkeypatch.setitem(globals(), "TOY_ATOMS", anisotropic)
    project = fake_project(tmp_path)

    run_mode(project, "chain", "lebail")
    run_mode(project, "chain", "fixed_atoms")

    start = refine_jobs(fake)[1]["start_model"]["toy"]
    # Ueq of the diagonal Uij above, which is what the driver is told to make
    # the atoms isotropic with.
    assert [atom["uiso"] for atom in start] == pytest.approx([0.012] * len(start))
    # Nb1 barely moved and its esd is eight times the shift; every other atom
    # moved by ten times its esd.
    refined = [
        {
            **atom,
            "adp": "I",
            "uij": None,
            "uiso": 0.0125 if atom["label"] == "Nb1" else 0.013,
            "uiso_esd": 0.004 if atom["label"] == "Nb1" else 0.0001,
            "occupancy_esd": None,
            "xyz_esd": [None, None, None],
        }
        for atom in start
    ]

    found = driver.find_undetermined(refined, start)

    assert [(entry["atom"], entry["parameter"]) for entry in found] == [("Nb1", "uiso")]
    assert found[0]["message"] == (
        "Nb1 Uiso 0.01250, esd 0.00400 more than its shift 0.00050"
    )


def test_run_mode_unplaced_element(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    with pytest.raises(
        PipelineError,
        match=(
            r"^structures\.lanthanum: La of the composition is not in the ttb/P4bm "
            r"prototype and no atoms place it; place it on one of its sites A1, A2, "
            r"B1, B2, O1, O2, O3, O4, O5$"
        ),
    ):
        run_mode(project, "unplaced", "lebail")
    assert fake.jobs == []


def test_run_mode_places_an_element_beside_the_cif_atom(tmp_path, fake) -> None:
    project = fake_project(
        tmp_path,
        '\n[structures.lanthanum.atoms]\nA1 = { Sr1 = "Sr", La1 = "La" }\n',
    )

    run_mode(project, "unplaced", "lebail")
    run_mode(project, "unplaced", "fixed_atoms")

    fixed = refine_jobs(fake)[-1]
    (la,) = [edit for edit in fixed["phases"][0]["atoms"] if edit.get("type") == "La"]
    assert la["copy"] == "Sr1"


# A synthetic cell in which Sr is on two sites, Sr1 (multiplicity 1) and Sr2
# (multiplicity 2), for placing La by an atoms table.
TWO_SITE_ATOMS = [
    {
        "label": "Sr1",
        "type": "Sr",
        "xyz": [0.0, 0.0, 0.0],
        "multiplicity": 1,
        "occupancy": 0.6,
    },
    {
        "label": "Sr2",
        "type": "Sr",
        "xyz": [0.5, 0.0, 0.5],
        "multiplicity": 2,
        "occupancy": 0.2,
    },
    {
        "label": "Nb1",
        "type": "Nb",
        "xyz": [0.5, 0.5, 0.5],
        "multiplicity": 1,
        "occupancy": 1.0,
    },
    {
        "label": "O1",
        "type": "O",
        "xyz": [0.5, 0.5, 0.0],
        "multiplicity": 3,
        "occupancy": 1.0,
    },
]
TWO_SITE_COMPOSITION = {"Sr": 0.8, "La": 0.2, "Nb": 1.0, "O": 3.0}


def two_site_phase(atoms):
    """A phase of the synthetic cell, one formula unit, with ``atoms`` its
    atoms table by site label, or None."""
    sites = None
    if atoms is not None:
        sites = [
            {"name": next(iter(held)), "atoms": held, "label": label}
            for label, held in atoms.items()
        ]
    return SimpleNamespace(
        key="toy",
        spec=SimpleNamespace(atoms=sites),
        z=1,
        composition=TWO_SITE_COMPOSITION,
    )


def test_nominal_edits_place_an_added_element_on_the_named_site_only() -> None:
    phase = two_site_phase({"A1": {"Sr1": "Sr", "La1": "La"}, "A2": {"Sr2": "Sr"}})

    edits = pipeline._nominal_edits(phase, TWO_SITE_ATOMS)

    added = [edit for edit in edits if "copy" in edit]
    # All 0.2 La per cell on the multiplicity 1 site beside Sr1; Sr scaled
    # by 0.8 on both its sites as without the table.
    assert added == [
        {"label": "La1", "type": "La", "copy": "Sr1", "occupancy": pytest.approx(0.2)}
    ]
    by_label = {edit["label"]: edit["occupancy"] for edit in edits}
    assert by_label["Sr1"] == pytest.approx(0.48)
    assert by_label["Sr2"] == pytest.approx(0.16)

    # Placed on both sites, La is shared in proportion to the Sr of the CIF.
    phase = two_site_phase(
        {"A1": {"Sr1": "Sr", "La1": "La"}, "A2": {"Sr2": "Sr", "La2": "La"}}
    )
    added = {
        edit["label"]: (edit["copy"], edit["occupancy"])
        for edit in pipeline._nominal_edits(phase, TWO_SITE_ATOMS)
        if "copy" in edit
    }
    assert added == {
        "La1": ("Sr1", pytest.approx(0.12)),
        "La2": ("Sr2", pytest.approx(0.04)),
    }

    # With Ba2 beside Sr2 on the second site, listed first, the host is still
    # Sr, the one element both sites have: La goes beside Sr1 and Sr2 in
    # proportion to their 0.6 and 0.4 of Sr per cell, not beside Ba2.
    atoms = TWO_SITE_ATOMS + [
        {
            "label": "Ba2",
            "type": "Ba",
            "xyz": [0.5, 0.0, 0.5],
            "multiplicity": 2,
            "occupancy": 0.5,
        },
    ]
    phase = two_site_phase(
        {
            "A1": {"Sr1": "Sr", "La1": "La"},
            "A2": {"Ba2": "Ba", "Sr2": "Sr", "La2": "La"},
        }
    )
    phase.composition = {**TWO_SITE_COMPOSITION, "Ba": 1.0, "La": 0.25}
    added = {
        edit["label"]: (edit["copy"], edit["occupancy"])
        for edit in pipeline._nominal_edits(phase, atoms)
        if "copy" in edit
    }
    assert added == {
        "La1": ("Sr1", pytest.approx(0.15)),
        "La2": ("Sr2", pytest.approx(0.05)),
    }


def test_nominal_edits_without_atoms_keep_the_rule_by_element() -> None:
    # Without atoms the composition takes only the CIF's elements, scaled.
    phase = two_site_phase(None)
    phase.composition = {**TWO_SITE_COMPOSITION, "La": 0.0}
    edits = pipeline._nominal_edits(phase, TWO_SITE_ATOMS)

    assert edits == [
        {"label": "Sr1", "occupancy": pytest.approx(0.48)},
        {"label": "Sr2", "occupancy": pytest.approx(0.16)},
    ]

    # A rule by element puts La beside every Sr atom, in proportion.
    structure = {"formula_units": 1, "composition": {"added": {"La": "Sr"}}}
    added = {
        edit["label"]: (edit["copy"], edit["occupancy"])
        for edit in composition_edits(TWO_SITE_ATOMS, TWO_SITE_COMPOSITION, structure)
        if "copy" in edit
    }
    assert added == {
        "La1": ("Sr1", pytest.approx(0.12)),
        "La2": ("Sr2", pytest.approx(0.04)),
    }


def test_nominal_edits_reject_an_element_the_atoms_do_not_place() -> None:
    # Ba shares Sr1's position in the CIF, but the table puts only Sr and La
    # on A1.
    atoms = TWO_SITE_ATOMS + [
        {
            "label": "Ba1",
            "type": "Ba",
            "xyz": [0.0, 0.0, 0.0],
            "multiplicity": 1,
            "occupancy": 0.1,
        },
    ]
    phase = two_site_phase({"A1": {"Sr1": "Sr", "La1": "La"}})
    phase.composition = {**TWO_SITE_COMPOSITION, "Ba": 0.1}

    with pytest.raises(
        PipelineError,
        match=r"^structures\.toy\.atoms: site A1 holds Ba \(Ba1\), which atoms "
        r"does not place there",
    ):
        pipeline._nominal_edits(phase, atoms)


def test_run_mode_held_instrument_check(tmp_path, monkeypatch) -> None:
    project = fake_project(tmp_path)
    monkeypatch.setattr(pipeline, "run_job", FakeGsas2(change=("U:2.0", "U:2.5")))

    with pytest.raises(
        PipelineError,
        match=r"^lebail: instrument parameters held by the run came out changed: "
        r"U 2\.0 -> 2\.5$",
    ):
        run_mode(project, "x0.10.powder", "lebail")


def old_result(folder: Path) -> Path:
    """An older result of the Le Bail mode, which no failure may report."""
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "x0.10.powder_lebail_result.json"
    path.write_text(
        json.dumps(
            {"completed": True, "stages": [{"name": "OLD STAGE", "status": "clean"}]}
        ),
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize(
    ("gsas2", "error", "wrote"),
    [
        (
            FakeGsas2(raises=Gsas2Error("GSAS-II crashed")),
            "Gsas2Error: GSAS-II crashed",
            False,
        ),
        (
            FakeGsas2(status="rejected"),
            "PipelineError: lebail: no stage was accepted",
            True,
        ),
        (
            FakeGsas2(final=False),
            "PipelineError: lebail: the run left no final model",
            True,
        ),
        (
            FakeGsas2(completed=False),
            "PipelineError: lebail: stage 'microstrain' failed",
            True,
        ),
    ],
    ids=["gsas2 error", "every stage rejected", "no final model", "a stage failed"],
)
def test_run_sequence_failures(tmp_path, monkeypatch, gsas2, error, wrote) -> None:
    project = fake_project(tmp_path)
    monkeypatch.setattr(pipeline, "run_job", gsas2)
    folder = project.root / "results" / "lebail" / "x0.10.powder"
    old_result(folder)

    outcomes = run_sequence(project, "x0.10.powder", ["lebail", "fixed_atoms"])

    (outcome,) = outcomes
    assert outcome.error.startswith(error)
    failure = (folder / "failure.md").read_text(encoding="utf-8")
    assert error in failure
    assert "OLD STAGE" not in failure
    assert "fixed_atoms not run" in failure
    assert "fake GSAS-II log line" in failure
    assert ("The run wrote no result." in failure) is not wrote
    summary = (folder / "summary.md").read_text(encoding="utf-8")
    assert "lebail" in summary and "OLD STAGE" not in summary
    assert not (project.root / "results" / "rietveld").exists()


def test_run_sequence_failure_at_set_up(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    outcomes = run_sequence(project, "unplaced", ["lebail"])

    assert outcomes[0].error.startswith("PipelineError: structures.lanthanum: La")
    folder = project.root / "results" / "lebail" / "unplaced"
    failure = (folder / "failure.md").read_text(encoding="utf-8")
    assert "The run wrote no result." in failure
    assert "(no log)" in failure
    assert (folder / "summary.md").is_file()


def test_run_sequence_runs_the_modes_in_order(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    outcomes = run_sequence(
        project, "chain", ["lebail", "fixed_atoms", "coordinates", "occupancies"]
    )

    assert [outcome.mode for outcome in outcomes] == list(MODES)
    assert all(outcome.error is None for outcome in outcomes)
    summary = project.root / "results" / "rietveld" / "chain" / "summary.md"
    assert summary.is_file() and outcomes[-1].paths["summary"] == str(summary)
    assert not (summary.parent / "failure.md").exists()
    with pytest.raises(PipelineError, match="no mode 'rietveld'"):
        run_sequence(project, "chain", ["rietveld"])


def test_build_refine_job_takes_the_start_model_and_starts() -> None:
    model = {"toy": [{"label": "Sr1", "xyz": [0.0, 0.0, 0.5]}]}

    job = build_refine_job(
        "toy.gpx",
        [{"name": "scale", "scale": True}],
        start_model=model,
        displacement_start=0.01,
        phase_fraction_start={"toy": 0.3},
    )

    assert job["start_model"] == model
    assert (job["displacement_start"], job["phase_fraction_start"]) == (
        0.01,
        {"toy": 0.3},
    )
    with pytest.raises(
        ValueError, match="start_model of 'toy' must be a list of atoms"
    ):
        build_refine_job(
            "toy.gpx", [{"scale": True}], start_model={"toy": [{"label": "X"}]}
        )


def test_run_mode_needs_the_mode_before(tmp_path, fake) -> None:
    project = fake_project(tmp_path)

    with pytest.raises(PipelineError, match="run lebail first$") as raised:
        run_mode(project, "chain", "fixed_atoms")
    assert raised.value.result is None and raised.value.log is None


def test_lebail_then_fixed_atoms_in_gsas2(tmp_path) -> None:
    try:
        find_gsas2()
    except FileNotFoundError as error:
        pytest.skip(str(error))

    root = tmp_path / "project"
    root.mkdir()
    two_theta = np.linspace(20.0, 60.0, 2001)
    counts = np.full(two_theta.size, 50.0)
    operations = space_group_operations("P4bm")
    for h in range(5):
        for k in range(h + 1):
            for l in range(4):
                if h + k + l == 0 or is_absent((h, k, l), operations):
                    continue
                d = 1.0 / math.sqrt((h * h + k * k) / 36.0 + l * l / 16.0)
                if 1.544426 / (2.0 * d) >= 1.0:
                    continue
                # K alpha 1 and 2, broader than the instrument alone.
                for wavelength, weight in ((1.540598, 1.0), (1.544426, 0.5)):
                    sine = wavelength / (2.0 * d)
                    centre = 2.0 * math.degrees(math.asin(sine))
                    height = weight * 2000.0 / (1 + h + k + l)
                    counts += height / (1.0 + ((two_theta - centre) / 0.06) ** 2)
    (root / "toy.xrdml").write_text(
        XRDML_SCAN.format(
            start=20.0, end=60.0, counts=" ".join(str(round(c)) for c in counts)
        ),
        encoding="utf-8",
    )
    (root / "toy.cif").write_text(TOY_CIF, encoding="utf-8")
    write_instprm(root / "cu.instprm", CAGLIOTI)
    (root / "xrdkit.toml").write_text(
        """
[project]
name = "toy"
version = 1

[refine]
two_theta = [20.5, 59.5]
background = { function = "chebyschev-1", terms = 3 }

[instruments.cu]
wavelength = [1.540598, 1.544426]
ka2 = true
instprm = "cu.instprm"

[structures.bronze]
library = "ttb/P4bm"
cif = "toy.cif"
composition = "Sr0.24Ba0.48Nb0.4O1.6"
cell = { a = 6.0, c = 4.0 }

[samples.toy]
file = "toy.xrdml"
instrument = "cu"
structures = ["bronze"]
form = "powder"
""",
        encoding="utf-8",
    )
    project = load_project(root)

    outcomes = run_sequence(
        project, "toy", ["lebail", "fixed_atoms"], Options(max_passes=2)
    )

    assert [outcome.error for outcome in outcomes] == [None, None], outcomes
    lebail, fixed = outcomes
    assert lebail.accepted[:3] == ["background and scale", "zero", "cell"]
    cell = lebail.final["phases"][0]["cell"]
    assert cell["length_a"] == pytest.approx(6.0, abs=0.01)
    assert cell["length_c"] == pytest.approx(4.0, abs=0.01)
    assert fixed.accepted[:3] == ["scale and background", "zero and cell", "size"]
    saved = json.loads(Path(fixed.paths["result"]).read_text(encoding="utf-8"))
    # The line heights are not the structure's, so one Uiso may be driven
    # negative and that stage rejected; it is recorded either way.
    assert [stage["name"] for stage in saved["stages"]][-1] == "overall Uiso"
    assert set(saved["start_model"]) == {"bronze"}
    assert saved["inputs"]["structures"][0]["library"] == "ttb/P4bm"
    assert (root / "results" / "rietveld" / "toy" / "summary.md").is_file()
