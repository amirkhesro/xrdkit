"""Tests for xrdkit.gsas2_driver, with GSASIIscriptable replaced by fakes.

The driver imports GSASII only once it runs a job, so its stage handling can
be tested here against small stand-ins for a project, a histogram and phases
that record what is asked of them.
"""

from types import SimpleNamespace

import pytest

from xrdkit import gsas2_driver as driver

INSTRUMENT = ("Type", "Lam1", "Lam2", "Zero", "U", "V", "W", "X", "Y", "Z", "SH/L")


class FakeHistogram:
    def __init__(self):
        self.name = "PWDR fake"
        self.id = 0
        self.InstrumentParameters = {key: ["PXC", "PXC", False] for key in INSTRUMENT}
        self.SampleParameters = {
            "Type": "Bragg-Brentano",
            "Scale": [1.0, False],
            "Shift": [0.0, False],
            "Transparency": [0.0, False],
        }
        self.Background = [["chebyschev-1", True, 3, 10.0, 1.0, 0.5], {}]
        self.residuals = {"wR": 12.5, "R": 9.25}
        self.calls = []

    def set_refinements(self, refs):
        self.calls.append(("set", refs))

    def clear_refinements(self, refs):
        self.calls.append(("clear", refs))

    def Limits(self, which):
        return {"lower": 11.0, "upper": 98.0}[which]


class FakePhase:
    def __init__(self, name, pid):
        self.name = name
        self.id = pid
        self.calls = []
        # GSAS-II's default sample broadening.
        self.data = {
            "Histograms": {
                "PWDR fake": {
                    "Size": ["isotropic", [1.0, 1.0, 1.0], [False, False, False]],
                    "Mustrain": [
                        "isotropic",
                        [1000.0, 1000.0, 1.0],
                        [True, False, False],
                    ],
                }
            }
        }

    def atoms(self):
        return [SimpleNamespace(label="La"), SimpleNamespace(label="B")]

    def set_refinements(self, refs):
        self.calls.append(("set", refs))

    def clear_refinements(self, refs):
        self.calls.append(("clear", refs))

    def set_HAP_refinements(self, refs, histograms):
        self.calls.append(("set HAP", refs, [h.name for h in histograms]))

    def clear_HAP_refinements(self, refs, histograms):
        self.calls.append(("clear HAP", refs, [h.name for h in histograms]))

    def get_cell_and_esd(self):
        return {"length_a": 4.156826}, {"length_a": 0.0}

    def getHAPvalues(self, histogram):
        return {"Scale": [1.0, False], "LeBail": False}


class FakeProject:
    """Refines by filling in the covariance, except on the stages in ``fail``,
    where, as GSAS-II does, it prints an error and saves nothing."""

    def __init__(self, fail=()):
        self.filename = "fake.gpx"
        self.histogram_ = FakeHistogram()
        self.phases_ = [FakePhase("LaB6", 0)]
        self.data = {"Covariance": {"data": {}}}
        self.controls = {}
        self.fail = set(fail)
        self.refinements = 0

    def histograms(self):
        return [self.histogram_]

    def histogram(self, name):
        return self.histogram_

    def phases(self):
        return self.phases_

    def set_Controls(self, key, value):
        self.controls[key] = value

    def refine(self):
        self.refinements += 1
        if self.refinements in self.fail:
            print(" ***** Refinement error *****")
            print("**** ERROR - Refinement failed")
            return
        self.data["Covariance"]["data"] = {
            "Rvals": {
                "GOF": 1.5,
                "chisq": 2250.0,
                "Nobs": 1004,
                "Nvars": 2,
                "converged": True,
                "Max shft/sig": 0.01,
            },
            "varyList": [":0:Scale", ":0:Zero"],
            "variables": [2.0 * self.refinements, -0.01],
            "covMatrix": [[0.04, 0.0], [0.0, 1.0e-6]],
            "depSigDict": {},
        }


def _set_calls(calls):
    return [call[1:] for call in calls if call[0].startswith("set")]


# accumulate_stages


def test_stages_accumulate() -> None:
    stages = driver.accumulate_stages(
        [
            {"name": "background and scale", "background": True, "scale": True},
            {"name": "zero", "zero": True},
            {"name": "U V W", "instrument": ["W", "U", "V"]},
            {"name": "X Y", "instrument": ["X", "Y"], "cell": True},
        ]
    )

    assert [stage["name"] for stage in stages] == [
        "background and scale",
        "zero",
        "U V W",
        "X Y",
    ]
    first, second, third, last = (stage["flags"] for stage in stages)
    assert first["background"] == {"type": "chebyschev-1", "terms": 6}
    assert first["scale"] and not first["zero"]
    assert second["scale"] and second["zero"]
    assert third["instrument"] == ["U", "V", "W"]
    assert last["instrument"] == ["U", "V", "W", "X", "Y"]
    assert last["cell"] == [driver.ALL_PHASES]
    assert last["background"] == first["background"]
    # Earlier stages are not changed by later ones.
    assert third["cell"] == []


def test_stage_reset_and_override() -> None:
    stages = driver.accumulate_stages(
        [
            {"background": {"terms": 3}, "scale": True, "zero": True},
            {"zero": False, "atoms": [{"phase": "LaB6", "flags": "x"}]},
            {"atoms": [{"phase": "LaB6", "flags": "U"}]},
            {"name": "only cell", "reset": True, "cell": ["LaB6"]},
        ]
    )

    flags = [stage["flags"] for stage in stages]
    assert stages[0]["name"] == "stage 1"
    assert flags[0]["background"] == {"type": "chebyschev-1", "terms": 3}
    assert flags[1]["scale"] and not flags[1]["zero"]
    assert flags[2]["atoms"] == {"LaB6": "XU"}
    assert flags[3] == {**driver.empty_flags(), "cell": ["LaB6"]}


@pytest.mark.parametrize(
    ("stages", "message"),
    [
        ([], "non-empty"),
        ([{"name": "a", "displacment": True}], "unknown keys"),
        ([{"instrument": ["U", "Zero"]}], "unknown instrument"),
        ([{"atoms": [{"phase": "LaB6", "flags": "XB"}]}], "atom flags"),
        ([{"scale": "yes"}], "scale must be"),
        ([{"background": {"terms": 0}}], "terms"),
        ([{"cell": 3}], "cell must be"),
    ],
)
def test_bad_stages_rejected(stages, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.accumulate_stages(stages)


# apply_flags


def test_apply_flags_clears_then_sets() -> None:
    project = FakeProject()
    histogram = project.histogram_
    (flags,) = (
        stage["flags"]
        for stage in driver.accumulate_stages(
            [
                {
                    "background": {"type": "cosine", "terms": 4},
                    "scale": True,
                    "zero": True,
                    "displacement": True,
                    "instrument": ["U", "V", "W"],
                    "cell": True,
                    "atoms": [{"phase": "LaB6", "flags": "U"}],
                    "phase_fractions": ["LaB6"],
                }
            ]
        )
    )

    driver.apply_flags(project, histogram, flags)

    clears = [call for call in histogram.calls if call[0] == "clear"]
    assert histogram.calls[: len(clears)] == clears, "every clear comes first"
    assert ("clear", {"Background": True}) in clears
    assert (
        "clear",
        {"Instrument Parameters": ["Zero", "U", "V", "W", "X", "Y", "Z", "SH/L"]},
    ) in clears
    assert ("clear", {"Sample Parameters": ["Scale", "Shift"]}) in clears
    assert _set_calls(histogram.calls) == [
        ({"Background": {"type": "cosine", "no. coeffs": 4, "refine": True}},),
        ({"Instrument Parameters": ["Zero", "U", "V", "W"]},),
        ({"Sample Parameters": ["Scale", "Shift"]},),
    ]
    (phase,) = project.phases_
    assert ("clear", {"Atoms": ["La", "B"]}) in phase.calls
    assert ("clear HAP", {"Scale": True}, ["PWDR fake"]) in phase.calls
    assert _set_calls(phase.calls) == [
        ({"Cell": True},),
        ({"Atoms": {"all": "U"}},),
        ({"Scale": True}, ["PWDR fake"]),
    ]


def test_apply_flags_unknown_phase() -> None:
    project = FakeProject()
    flags = {**driver.empty_flags(), "cell": ["Si"]}

    with pytest.raises(ValueError, match="Si"):
        driver.apply_flags(project, project.histogram_, flags)


# refine


def test_refine_records_each_stage() -> None:
    project = FakeProject()
    G2sc = SimpleNamespace(G2Project=lambda gpx: project)
    job = {
        "gpx": "fake.gpx",
        "limits": [11, 98],
        "cycles": 5,
        "stages": [
            {"name": "scale", "scale": True},
            {"name": "zero", "zero": True},
        ],
    }

    result = driver.refine(G2sc, job)

    assert result["completed"]
    assert project.controls == {"cycles": 5}
    assert ("set", {"Limits": [11.0, 98.0]}) in project.histogram_.calls
    assert result["limits"] == [11.0, 98.0]
    first, second = result["stages"]
    assert first["name"] == "scale"
    assert first["rwp"] == 12.5
    assert first["rp"] == 9.25
    assert first["gof"] == 1.5
    assert first["reduced_chi_squared"] == pytest.approx(2.25)
    assert first["n_variables"] == 2
    assert first["n_observations"] == 1004
    assert first["parameters"][":0:Scale"] == {"value": 2.0, "esd": pytest.approx(0.2)}
    assert second["parameters"][":0:Scale"]["value"] == 4.0
    assert second["flags"]["scale"] and second["flags"]["zero"]
    final = result["final"]
    assert final["instrument"]["Zero"]["esd"] == pytest.approx(1.0e-3)
    assert final["instrument"]["U"]["esd"] is None
    assert final["sample"]["Scale"]["esd"] == pytest.approx(0.2)
    assert final["background"] == {
        "type": "chebyschev-1",
        "coefficients": [10.0, 1.0, 0.5],
    }
    assert final["phases"][0]["cell"] == {"length_a": 4.156826}
    assert "exports" not in result


def test_refine_stops_at_a_failed_stage() -> None:
    project = FakeProject(fail={2})
    G2sc = SimpleNamespace(G2Project=lambda gpx: project)
    job = {
        "gpx": "fake.gpx",
        "stages": [{"scale": True}, {"zero": True}, {"instrument": ["U"]}],
    }

    result = driver.refine(G2sc, job)

    assert not result["completed"]
    assert project.refinements == 2
    first, failed = result["stages"]
    assert "error" not in first
    assert failed["name"] == "stage 2"
    assert failed["error"].startswith("RefinementError")
    assert "Refinement failed" in failed["error"]
    # The covariance of the last good stage is kept for the final values.
    assert project.data["Covariance"]["data"]["variables"][0] == 2.0
    assert result["final"]["sample"]["Scale"]["value"] == 1.0


def test_refine_sets_broadening_before_the_stages() -> None:
    project = FakeProject()
    G2sc = SimpleNamespace(G2Project=lambda gpx: project)
    job = {
        "gpx": "fake.gpx",
        "broadening": {"*": {"size": 10, "mustrain": 0, "lgmix": 0}},
        "stages": [{"scale": True}],
    }

    result = driver.refine(G2sc, job)

    hap = project.phases_[0].data["Histograms"]["PWDR fake"]
    assert hap["Size"][1] == [10.0, 1.0, 0.0]
    assert hap["Mustrain"][1] == [0.0, 1000.0, 0.0]
    assert hap["Mustrain"][2][0] is False
    phase = result["final"]["phases"][0]
    assert phase["size"] == {
        "type": "isotropic",
        "value": 10.0,
        "lorentzian_fraction": 0.0,
    }
    assert phase["mustrain"]["value"] == 0.0


def test_broadening_of_one_part_leaves_the_rest() -> None:
    project = FakeProject()
    driver.set_broadening(project, project.histogram_, {"LaB6": {"mustrain": 50}})

    hap = project.phases_[0].data["Histograms"]["PWDR fake"]
    assert hap["Size"][1] == [1.0, 1.0, 1.0]
    assert hap["Mustrain"][1] == [50.0, 1000.0, 1.0]


def test_refine_reports_default_broadening() -> None:
    project = FakeProject()
    result = driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": [{"scale": True}]},
    )

    phase = result["final"]["phases"][0]
    assert phase["size"]["value"] == 1.0
    assert phase["mustrain"]["value"] == 1000.0


@pytest.mark.parametrize(
    ("broadening", "message"),
    [
        ({"LaB6": {"size": 0.0}}, "size must lie between"),
        ({"LaB6": {"size": 100.0}}, "clamps"),
        ({"LaB6": {"mustrain": -1.0}}, "cannot be negative"),
        ({"LaB6": {"lgmix": 1.5}}, "lgmix must lie"),
        ({"LaB6": {"strain": 0.0}}, "size, mustrain and/or lgmix"),
        ({"LaB6": {}}, "size, mustrain and/or lgmix"),
    ],
)
def test_bad_broadening_rejected(broadening, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.check_broadening(broadening)


def test_broadening_must_be_a_dict() -> None:
    with pytest.raises(TypeError, match="must map"):
        driver.check_broadening(["LaB6"])


def test_stage_of_the_wrong_type() -> None:
    with pytest.raises(TypeError, match="stage 2 is not a dict"):
        driver.accumulate_stages([{"zero": True}, ["scale"]])
    with pytest.raises(TypeError, match="instrument must be a list"):
        driver.accumulate_stages([{"instrument": {"U": True}}])


def test_broadening_for_unknown_phase() -> None:
    project = FakeProject()
    with pytest.raises(ValueError, match="Si"):
        driver.set_broadening(project, project.histogram_, {"Si": {"size": 10.0}})


def test_refine_rejects_bad_stages_before_opening_the_project() -> None:
    def no_project(gpx):
        raise AssertionError("project opened")

    with pytest.raises(ValueError, match="unknown keys"):
        driver.refine(
            SimpleNamespace(G2Project=no_project),
            {"gpx": "fake.gpx", "stages": [{"zero": True, "shift": True}]},
        )


def test_plain_makes_json_types() -> None:
    assert driver._plain({"a": (1, float("nan")), 2: {3.0}}) == {
        "a": [1, None],
        "2": [3.0],
    }
