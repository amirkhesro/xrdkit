"""Tests for xrdkit.gsas2_driver, with GSASIIscriptable replaced by fakes.

The driver imports GSASII only once it runs a job, so its stage handling can
be tested here against small stand-ins for a project, a histogram and phases
that record what is asked of them.
"""

from types import SimpleNamespace

import pytest

from xrdkit import gsas2_driver as driver

INSTRUMENT = ("Type", "Lam1", "Lam2", "Zero", "U", "V", "W", "X", "Y", "Z", "SH/L")

# GSAS-II's atom record of a nuclear phase: label, type, refinement flags,
# x, y, z, occupancy, site symmetry, multiplicity, I or A, Uiso, six Uij and
# a random id; so cx = 3, ct = 1, cs = 7 and cia = 9.
ATOM_POINTERS = [3, 1, 7, 9]


def _atom(label, element, xyz, uiso=None, uij=None):
    adp = ["I", uiso, 0, 0, 0, 0, 0, 0] if uij is None else ["A", 0.0, *uij]
    return [label, element, "", *xyz, 1.0, "m-3m", 1, *adp, hash(label)]


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
            "General": {"AtomPtrs": ATOM_POINTERS},
            "Atoms": [
                _atom("La", "La", [0.0, 0.0, 0.0], uiso=0.0086),
                _atom("B", "B", [0.5, 0.5, 0.2], uij=[0.01, 0.01, 0.006, 0, 0, 0]),
            ],
            "Histograms": {
                "PWDR fake": {
                    "Size": ["isotropic", [1.0, 1.0, 1.0], [False, False, False]],
                    "Mustrain": [
                        "isotropic",
                        [1000.0, 1000.0, 1.0],
                        [True, False, False],
                    ],
                }
            },
        }

    def atoms(self):
        return [SimpleNamespace(label=atom[0]) for atom in self.data["Atoms"]]

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
        self.data = {
            "Covariance": {"data": {}},
            "Controls": {"data": {}},
            "Constraints": {"data": {"Phase": []}},
            # As in a GSAS-II project, each phase's data is part of the
            # project's, so that going back to a snapshot of it takes the
            # atoms back too.
            "Phases": {phase.name: phase.data for phase in self.phases_},
        }
        self.controls = {}
        self.fail = set(fail)
        self.refinements = 0
        self.saved = []
        self.reloads = 0
        self.constraints = []
        self.equations = []

    def add_EquivConstr(self, varlist):
        # As GSAS-II stores an equivalence: multiplier and variable pairs,
        # then two Nones and "e".
        self.constraints.append(list(varlist))
        self.data["Constraints"]["data"]["Phase"].append(
            [[1.0, name] for name in varlist] + [None, None, "e"]
        )

    def add_HoldConstr(self, varlist):
        # GSAS-II holds each variable with a constraint of its own.
        for name in varlist:
            self.constraints.append([name])
            self.data["Constraints"]["data"]["Phase"].append(
                [[1.0, name], None, None, "h"]
            )

    def add_EqnConstr(self, total, varlist, multlist):
        self.equations.append((total, list(varlist), list(multlist)))
        self.data["Constraints"]["data"]["Phase"].append(
            [[m, name] for m, name in zip(multlist, varlist)] + [total, None, "c"]
        )

    def phase_constraints(self):
        return [
            [term[1] for term in constraint[:-3]]
            for constraint in self.data["Constraints"]["data"]["Phase"]
        ]

    def histograms(self):
        return [self.histogram_]

    def histogram(self, name):
        return self.histogram_

    def phases(self):
        return self.phases_

    def set_Controls(self, key, value):
        self.controls[key] = value
        self.data["Controls"]["data"][key] = value

    def save(self):
        self.saved.append(dict(self.data["Controls"]["data"]))

    def reload(self):
        self.reloads += 1

    def refine(self):
        self.refinements += 1
        if self.refinements in self.fail:
            print(" ***** Refinement error *****")
            print("**** ERROR - Refinement failed")
            return
        print(" ***** Refinement successful *****")
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


def test_size_and_mustrain_flags() -> None:
    stages = driver.accumulate_stages(
        [{"name": "size", "size": True}, {"name": "strain", "mustrain": ["LaB6"]}]
    )
    assert stages[0]["flags"]["size"] == [driver.ALL_PHASES]
    assert stages[0]["flags"]["mustrain"] == []
    assert stages[1]["flags"]["mustrain"] == ["LaB6"]

    project = FakeProject()
    driver.apply_flags(project, project.histogram_, stages[1]["flags"])

    (phase,) = project.phases_
    assert ("clear HAP", {"Size": True}, ["PWDR fake"]) in phase.calls
    assert ("clear HAP", {"Mustrain": True}, ["PWDR fake"]) in phase.calls
    assert _set_calls(phase.calls) == [
        ({"Size": {"type": "isotropic", "refine": True}}, ["PWDR fake"]),
        ({"Mustrain": {"type": "isotropic", "refine": True}}, ["PWDR fake"]),
    ]
    with pytest.raises(ValueError, match="Si"):
        driver.apply_flags(
            project, project.histogram_, {**driver.empty_flags(), "size": ["Si"]}
        )


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
        "esd": None,
        "lorentzian_fraction": 0.0,
    }
    assert phase["mustrain"]["value"] == 0.0


def _le_bail_g2sc(project, calls, ok=True):
    """A GSASIIscriptable whose DoLeBail records the controls it was saved
    with and the cycles it was asked for."""

    def do_le_bail(gpx, cycles):
        calls.append({"cycles": cycles, "controls": project.saved[-1]})
        if not ok:
            print(" ***** LeBail fit error *****")
            return False, {"msg": "Ouch #8: no reflections in data range."}
        return True, {"Rwp": 20.0, "GOF": 3.0}

    return SimpleNamespace(
        G2Project=lambda gpx: project,
        G2strMain=SimpleNamespace(DoLeBail=do_le_bail),
    )


def test_le_bail_extracted_when_switched_on() -> None:
    project = FakeProject()
    calls = []
    job = {
        "gpx": "fake.gpx",
        "le_bail_cycles": 4,
        "stages": [
            {"name": "background and scale", "scale": True, "le_bail": True},
            {"name": "zero", "zero": True},
            {"name": "no Le Bail", "reset": True, "scale": True},
            {"name": "Le Bail again", "le_bail": ["LaB6"]},
        ],
    }

    result = driver.refine(_le_bail_g2sc(project, calls), job)

    assert result["completed"]
    assert [call["cycles"] for call in calls] == [4, 4]
    assert all(call["controls"]["newLeBail"] for call in calls)
    assert project.reloads == 2
    first, zero, off, again = result["stages"]
    assert first["le_bail_extraction"] == {"cycles": 4, "rwp": 20.0, "gof": 3.0}
    assert "le_bail_extraction" not in zero
    assert "le_bail_extraction" not in off
    assert again["le_bail_extraction"]["cycles"] == 4


def test_le_bail_extraction_defaults_and_failure() -> None:
    project = FakeProject()
    calls = []
    job = {"gpx": "fake.gpx", "stages": [{"scale": True}, {"le_bail": True}]}

    result = driver.refine(_le_bail_g2sc(project, calls, ok=False), job)

    assert calls[0]["cycles"] == driver.DEFAULT_LE_BAIL_CYCLES
    assert not result["completed"]
    assert project.refinements == 1, "the failed stage is not refined"
    assert result["stages"][1]["error"].startswith("RefinementError: Ouch #8")


def _settling_project():
    """A project whose zero halves its distance to -0.01 with each refinement,
    esd 0.001, while its scale grows by 10 esds each time."""
    project = FakeProject()
    original = project.refine

    def refine():
        original()
        covariance = project.data["Covariance"]["data"]
        covariance["variables"][1] = -0.01 + 0.004 * 0.5**project.refinements

    project.refine = refine
    return project


def test_passes_until_settled_leaving_out_a_le_bail_scale() -> None:
    project = _settling_project()
    job = {
        "gpx": "fake.gpx",
        "max_passes": 20,
        "stages": [{"name": "zero", "scale": True, "zero": True, "le_bail": True}],
    }

    result = driver.refine(_le_bail_g2sc(project, []), job)

    (stage,) = result["stages"]
    # Zero shifts of 1, 0.5, 0.25, 0.125 and 0.0625 esds: settled at the sixth.
    assert project.refinements == 6
    assert stage["passes_converged"]
    passes = stage["passes"]
    assert passes[0]["max_shift_over_esd"] is None
    assert [p["parameter"] for p in passes[1:]] == [":0:Zero"] * 5
    assert [p["max_shift_over_esd"] for p in passes[1:]] == pytest.approx(
        [1.0, 0.5, 0.25, 0.125, 0.0625]
    )
    assert [p["settled"] for p in passes] == [False] * 5 + [True]
    assert stage["parameters"][":0:Scale"]["value"] == 12.0


def test_passes_stop_at_the_most_without_settling() -> None:
    project = _settling_project()
    job = {
        "gpx": "fake.gpx",
        "max_passes": 3,
        "pass_tolerance": 0.5,
        "stages": [{"scale": True, "zero": True}],
    }

    result = driver.refine(SimpleNamespace(G2Project=lambda gpx: project), job)

    (stage,) = result["stages"]
    # Without Le Bail the scale counts, and it moves 10 esds every pass.
    assert project.refinements == 3
    assert not stage["passes_converged"]
    assert stage["passes"][-1]["parameter"] == ":0:Scale"
    assert stage["passes"][-1]["max_shift_over_esd"] == pytest.approx(10.0)


def test_one_pass_by_default() -> None:
    project = _settling_project()
    result = driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": [{"zero": True}, {"scale": True}]},
    )

    assert project.refinements == 2
    assert all("passes" not in stage for stage in result["stages"])


def test_refine_reports_size_and_mustrain_esds() -> None:
    project = FakeProject()
    original = project.refine

    def refine_size():
        original()
        covariance = project.data["Covariance"]["data"]
        covariance["varyList"].append("0:0:Size;i")
        covariance["variables"].append(0.15)
        covariance["covMatrix"] = [
            [0.04, 0.0, 0.0],
            [0.0, 1.0e-6, 0.0],
            [0.0, 0.0, 1.0e-4],
        ]

    project.refine = refine_size
    result = driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": [{"size": True}]},
    )

    (stage,) = result["stages"]
    assert stage["parameters"]["0:0:Size;i"] == {
        "value": 0.15,
        "esd": pytest.approx(0.01),
    }
    phase = result["final"]["phases"][0]
    assert phase["size"]["esd"] == pytest.approx(0.01)
    assert phase["mustrain"]["esd"] is None


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


# Structure


def test_edit_atoms_changes_occupancies_and_adds_atoms_on_sites() -> None:
    phase = FakePhase("LaB6", 0)
    setups = []
    G2sc = SimpleNamespace(SetupGeneral=lambda data, dirname: setups.append(data))

    driver.edit_atoms(
        G2sc,
        phase,
        [
            {"label": "La", "occupancy": 0.9},
            {"label": "Ce", "type": "Ce", "copy": "La", "occupancy": 0.1},
            {"label": "C", "type": "C", "copy": "B", "occupancy": 0.05},
        ],
    )

    la, b, ce, c = phase.data["Atoms"]
    assert la[6] == 0.9
    assert ce[:3] == ["Ce", "Ce", ""]
    assert ce[3:6] == la[3:6] and ce[7:17] == la[7:17]
    assert ce[6] == 0.1
    assert ce[17] != la[17], "a new random id"
    # The copy carries the anisotropic displacement parameters with it.
    assert c[9:17] == b[9:17] and c[9] == "A"
    assert c[6] == 0.05 and b[6] == 1.0
    assert setups == [phase.data]


@pytest.mark.parametrize(
    ("edits", "message"),
    [
        ([{"label": "Xx", "occupancy": 0.5}], "no atom 'Xx'"),
        (
            [{"label": "Ce", "type": "Ce", "copy": "Xx", "occupancy": 0.1}],
            "no atom 'Xx' to copy",
        ),
        ([{"label": "B", "type": "C", "copy": "La", "occupancy": 0.1}], "already has"),
    ],
)
def test_edit_atoms_rejects_missing_and_taken_labels(edits, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.edit_atoms(
            SimpleNamespace(SetupGeneral=None), FakePhase("LaB6", 0), edits
        )


@pytest.mark.parametrize(
    ("edits", "message"),
    [
        ([{"label": "La"}], "needs a label"),
        ([{"occupancy": 0.5}], "needs a label"),
        ([{"label": "La", "occ": 0.5}], "needs a label"),
        ([{"label": "Ce", "copy": "La", "occupancy": 0.1}], "needs a type"),
        ([{"label": "Ce", "type": "Ce", "copy": "La"}], "needs a type"),
        ([{"label": "La", "type": "Ce", "occupancy": 0.1}], "only for a new atom"),
        ([{"label": "La", "occupancy": 1.2}], "between 0 and 1"),
        ([{"label": "La", "xyz": [0.0, 0.5]}], "three fractional"),
        ([{"label": "La", "uiso": float("nan")}], "must be finite"),
    ],
)
def test_bad_atom_edits_rejected(edits, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.check_atom_edits(edits, "LaB6")


def test_atom_edits_must_be_a_list_of_dicts() -> None:
    with pytest.raises(TypeError, match="must be a list"):
        driver.check_atom_edits({"label": "La"}, "LaB6")
    with pytest.raises(TypeError, match="not a dict"):
        driver.check_atom_edits(["La"], "LaB6")


def test_overall_uiso_makes_every_atom_isotropic_and_constrains_them() -> None:
    project = FakeProject()
    (phase,) = project.phases_
    stages = [{"scale": True}, {"name": "Uiso", "overall_uiso": True}]

    result = driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": stages, "overall_uiso_start": {"*": 0.007}},
    )

    assert result["completed"]
    assert [atom[9:11] for atom in phase.data["Atoms"]] == [["I", 0.007]] * 2
    assert project.constraints == [["0::AUiso:0", "0::AUiso:1"]]
    first, second = (_set_calls(calls) for calls in _stage_calls(phase.calls))
    assert ({"Atoms": {"all": "U"}},) not in first
    assert ({"Atoms": {"all": "U"}},) in second
    atoms = result["final"]["phases"][0]["atoms"]
    assert [(atom["label"], atom["adp"], atom["uiso"]) for atom in atoms] == [
        ("La", "I", 0.007),
        ("B", "I", 0.007),
    ]


def _stage_calls(calls):
    """The phase calls of each apply_flags, split where the clears begin."""
    stages, current = [], []
    for call in calls:
        if call == ("clear", {"Cell": True}) and current:
            stages.append(current)
            current = []
        current.append(call)
    return [*stages, current]


def test_overall_uiso_joins_the_other_atom_flags_and_is_set_up_once() -> None:
    project = FakeProject()
    (phase,) = project.phases_
    stages = [
        {"overall_uiso": ["LaB6"]},
        {"atoms": [{"phase": "LaB6", "flags": "X"}]},
    ]

    driver.refine(
        _g2sc(project),
        {"gpx": "fake.gpx", "stages": stages, "overall_uiso_start": {"LaB6": 0.01}},
    )

    assert len(project.constraints) == 1
    assert ({"Atoms": {"all": "XU"}},) in _set_calls(phase.calls)


def test_overall_uiso_without_a_start_fails_the_stage() -> None:
    project = FakeProject()
    result = driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": [{"overall_uiso": True}]},
    )

    assert not result["completed"]
    assert "no Uiso for phase 'LaB6'" in result["stages"][0]["error"]
    assert project.refinements == 1, "the starting model is computed, not refined"
    assert result["final_from"] == "job start"


def test_background_start_is_set_before_the_stages() -> None:
    project = FakeProject()
    job = {
        "gpx": "fake.gpx",
        "background_start": {"type": "chebyschev-1", "coefficients": [900, -40, 7]},
        "stages": [{"background": {"terms": 3}, "scale": True}],
    }

    driver.refine(SimpleNamespace(G2Project=lambda gpx: project), job)

    calls = _set_calls(project.histogram_.calls)
    start = {
        "Background": {
            "type": "chebyschev-1",
            "no. coeffs": 3,
            "coeffs": [900.0, -40.0, 7.0],
            "refine": False,
        }
    }
    assert calls.index((start,)) < calls.index(
        ({"Background": {"type": "chebyschev-1", "no. coeffs": 3, "refine": True}},)
    )


@pytest.mark.parametrize(
    ("start", "message"),
    [
        ({"coefficients": []}, "at least one"),
        ({"type": "cosine"}, "at least one"),
        ({"coeffs": [1.0]}, "background_start must be"),
    ],
)
def test_bad_background_start_rejected(start, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.check_background_start(start)


def test_bad_uiso_start_rejected() -> None:
    with pytest.raises(ValueError, match="must be positive"):
        driver.check_uiso_start({"LaB6": 0.0})
    with pytest.raises(TypeError, match="must map"):
        driver.check_uiso_start(0.01)


# Named atom flags, Uiso groups and shared sites


def _xinel(site_symmetry):
    """GetCSxinel for the fakes: x and y tied on "m", z alone on "4",
    nothing free on "m-3m"."""
    return {
        "m": [[1, 1, 2], [1.0, 1.0, 1.0]],
        "4": [[0, 0, 1], [0.0, 0.0, 1.0]],
        "1": [[1, 2, 3], [1.0, 1.0, 1.0]],
        "m-3m": [[0, 0, 0], [0.0, 0.0, 0.0]],
    }[site_symmetry]


def _g2sc(project):
    return SimpleNamespace(
        G2Project=lambda gpx: project, G2spc=SimpleNamespace(GetCSxinel=_xinel)
    )


def _shared_site_project():
    """La and Ce sharing a site of symmetry m, B alone on 1."""
    project = FakeProject()
    (phase,) = project.phases_
    la, b = phase.data["Atoms"]
    la[7], b[7] = "m", "1"
    ce = list(la)
    ce[:2] = ["Ce", "Ce"]
    phase.data["Atoms"].append(ce)
    return project, phase


def test_atom_flags_and_uiso_groups_accumulate() -> None:
    stages = driver.accumulate_stages(
        [
            {"overall_uiso": True},
            {
                "overall_uiso": False,
                "uiso_groups": {"LaB6": [["La", "Ce"], ["B"]]},
            },
            {"atom_flags": [{"phase": "LaB6", "labels": ["La", "Ce"], "flags": "x"}]},
            {"atom_flags": [{"phase": "LaB6", "labels": ["B"], "flags": "X"}]},
            {"uiso_groups": {"LaB6": []}},
        ]
    )

    flags = [stage["flags"] for stage in stages]
    assert flags[0]["uiso_groups"] == {}
    assert flags[1]["overall_uiso"] == []
    assert flags[1]["uiso_groups"] == {"LaB6": [["La", "Ce"], ["B"]]}
    assert flags[2]["atom_flags"] == {"LaB6": {"La": "X", "Ce": "X"}}
    assert flags[3]["atom_flags"] == {"LaB6": {"La": "X", "Ce": "X", "B": "X"}}
    assert flags[3]["uiso_groups"] == flags[1]["uiso_groups"]
    assert flags[4]["uiso_groups"] == {}


@pytest.mark.parametrize(
    ("stages", "message"),
    [
        (
            [{"overall_uiso": True, "uiso_groups": {"LaB6": [["La"]]}}],
            "overall Uiso and Uiso groups",
        ),
        ([{"uiso_groups": {"LaB6": [["La"], ["La", "B"]]}}], "more than one"),
        ([{"uiso_groups": {"LaB6": [[]]}}], "lists of labels"),
        ([{"atom_flags": [{"phase": "LaB6", "flags": "X"}]}], "each atom_flags"),
        (
            [{"atom_flags": [{"phase": "LaB6", "labels": [], "flags": "X"}]}],
            "list of labels",
        ),
        (
            [{"atom_flags": [{"phase": "LaB6", "labels": ["La"], "flags": "B"}]}],
            "atom flags",
        ),
    ],
)
def test_bad_atom_flags_and_uiso_groups_rejected(stages, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.accumulate_stages(stages)


def test_named_atom_flags_set_atom_by_atom() -> None:
    project = FakeProject()
    (phase,) = project.phases_
    (flags,) = (
        stage["flags"]
        for stage in driver.accumulate_stages(
            [
                {
                    "uiso_groups": {"LaB6": [["La"]]},
                    "atom_flags": [{"phase": "LaB6", "labels": ["B"], "flags": "X"}],
                }
            ]
        )
    )

    driver.apply_flags(project, project.histogram_, flags)

    assert ({"Atoms": {"La": "U", "B": "X"}},) in _set_calls(phase.calls)
    with pytest.raises(ValueError, match="no atoms"):
        driver.atom_flag_strings(phase, {**flags, "atom_flags": {"LaB6": {"Xx": "X"}}})


def test_uiso_groups_replace_the_overall_constraint() -> None:
    project, phase = _shared_site_project()
    stages = [
        {"name": "overall", "overall_uiso": True},
        {
            "name": "groups",
            "overall_uiso": False,
            "uiso_groups": {"LaB6": [["La", "Ce"], ["B"]]},
        },
        {"name": "again"},
    ]

    result = driver.refine(
        _g2sc(project),
        {"gpx": "fake.gpx", "stages": stages, "overall_uiso_start": {"*": 0.01}},
    )

    assert result["completed"], result["stages"]
    # The overall constraint, on all three atoms, went; the La-Ce group came.
    assert project.constraints[0] == ["0::AUiso:0", "0::AUiso:1", "0::AUiso:2"]
    assert project.phase_constraints() == [["0::AUiso:0", "0::AUiso:2"]]
    assert [atom[9:11] for atom in phase.data["Atoms"]] == [["I", 0.01]] * 3


def test_atoms_sharing_a_site_move_together() -> None:
    project, _ = _shared_site_project()
    stages = [
        {"atom_flags": [{"phase": "LaB6", "labels": ["La", "Ce"], "flags": "X"}]},
        {"atom_flags": [{"phase": "LaB6", "labels": ["B"], "flags": "X"}]},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert result["completed"], result["stages"]
    # Each stage records the atoms as it left them.
    first, second = result["stages"]
    assert [atom["label"] for atom in first["atoms"]["LaB6"]] == ["La", "B", "Ce"]
    assert second["atoms"]["LaB6"][0]["xyz"] == [0.0, 0.0, 0.0]
    # On m, x and y share a code, so x and z are tied, once, y following x by
    # the site's symmetry; B shares no site.
    assert project.phase_constraints() == [
        ["0::dAx:0", "0::dAx:2"],
        ["0::dAz:0", "0::dAz:2"],
    ]
    # Both atoms on the site stay refined, since GSAS-II ignores an
    # equivalence with parameters that are not.
    (phase,) = project.phases_
    assert ({"Atoms": {"La": "X", "Ce": "X"}},) in _set_calls(phase.calls)
    assert ({"Atoms": {"all": "X"}},) in _set_calls(phase.calls)


def test_part_of_a_shared_site_cannot_be_refined() -> None:
    project, _ = _shared_site_project()
    stages = [{"atom_flags": [{"phase": "LaB6", "labels": ["La"], "flags": "X"}]}]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert "share a site" in result["stages"][0]["error"]


def test_coordinate_esds_fill_in_what_the_site_ties() -> None:
    phase = FakePhase("LaB6", 0)
    esds = {"0::dAx:0": 0.002, "0::dAz:0": 0.004}
    phase.data["Atoms"][0][7] = "m"

    assert driver._xyz_esd(phase, 0, "m", esds, _xinel) == [0.002, 0.002, 0.004]
    assert driver._xyz_esd(phase, 0, "m", esds, None) == [0.002, None, 0.004]
    assert driver._xyz_esd(phase, 1, "1", esds, _xinel) == [None, None, None]


def test_operators_expand_centring_and_inversion() -> None:
    phase = FakePhase("LaB6", 0)
    phase.data["General"]["SGData"] = {
        "SGOps": [[[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.5]]],
        "SGCen": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "SGInv": True,
    }

    operators = driver._operators(phase)

    assert len(operators) == 4
    assert operators[1] == {
        "rotation": [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
        "translation": [0.0, 0.0, 0.5],
    }
    assert operators[2]["translation"] == [0.5, 0.5, 0.0]
    assert driver._operators(FakePhase("LaB6", 0)) is None


def test_plain_makes_json_types() -> None:
    assert driver._plain({"a": (1, float("nan")), 2: {3.0}}) == {
        "a": [1, None],
        "2": [3.0],
    }


def test_edit_atoms_moves_atoms_and_sets_their_uiso() -> None:
    phase = FakePhase("LaB6", 0)
    phase.data["General"]["SGData"] = {}
    asked = []

    def site_symmetry(xyz, sg_data):
        asked.append(list(xyz))
        return "m-3m", 1, 1, {}

    G2sc = SimpleNamespace(
        SetupGeneral=lambda data, dirname: None,
        G2spc=SimpleNamespace(SytSym=site_symmetry),
    )
    driver.edit_atoms(
        G2sc,
        phase,
        [
            {"label": "B", "xyz": [0.5, 0.5, 0.25], "uiso": 0.012},
            {"label": "Ce", "type": "Ce", "copy": "B", "occupancy": 0.1, "uiso": 0.02},
        ],
    )

    _, b, ce = phase.data["Atoms"]
    assert b[3:6] == [0.5, 0.5, 0.25]
    assert b[9:11] == ["I", 0.012]
    # The new atom takes the moved position, and a Uiso of its own.
    assert ce[3:6] == [0.5, 0.5, 0.25]
    assert ce[9:11] == ["I", 0.02]
    assert asked == [[0.5, 0.5, 0.25]]

    G2sc.G2spc.SytSym = lambda xyz, sg_data: ("4mm", 6, 1, {})
    with pytest.raises(ValueError, match="site of symmetry 4mm"):
        driver.edit_atoms(G2sc, phase, [{"label": "La", "xyz": [0.1, 0.0, 0.0]}])


def test_scale_start_is_set_before_the_stages() -> None:
    project = FakeProject()
    driver.refine(
        SimpleNamespace(G2Project=lambda gpx: project),
        {"gpx": "fake.gpx", "stages": [{"zero": True}], "scale_start": 14.2},
    )

    assert project.histogram_.SampleParameters["Scale"][0] == 14.2
    with pytest.raises(ValueError, match="scale_start must be positive"):
        driver.check_scale_start(0)


# Coordinates by site, the origin, occupancies and the sanity check


def test_coordinates_origin_and_occupancies_accumulate() -> None:
    entry = {"phase": "LaB6", "sites": ["La", "B"], "elements": ["La"]}
    stages = driver.accumulate_stages(
        [
            {
                "coordinates": {"LaB6": {"La": "Zx", "B": "y"}},
                "origin": {"LaB6": {"site": "B", "axis": "Z"}},
            },
            {"coordinates": {"LaB6": {"La": "y", "B": "all"}}},
            {"occupancies": [entry], "origin": {"LaB6": None}},
            {"occupancies": [entry]},
        ]
    )

    flags = [stage["flags"] for stage in stages]
    assert flags[0]["coordinates"] == {"LaB6": {"La": "xz", "B": "y"}}
    assert flags[0]["origin"] == {"LaB6": {"site": "B", "axis": "z"}}
    assert flags[1]["coordinates"] == {"LaB6": {"La": "xyz", "B": "all"}}
    assert flags[1]["origin"] == flags[0]["origin"]
    assert flags[2]["origin"] == {}
    assert flags[3]["occupancies"] == [entry], "the same entry is not added twice"
    assert all(set(stage) == {"name", "flags"} for stage in stages)


@pytest.mark.parametrize(
    ("stages", "message"),
    [
        ([{"coordinates": {"LaB6": {"La": "xw"}}}], "drawn from xyz"),
        ([{"coordinates": {"LaB6": {}}}], "map site labels"),
        ([{"origin": {"LaB6": {"site": "B", "axis": "w"}}}], "origin of 'LaB6'"),
        ([{"origin": {"LaB6": {"site": "B"}}}], "origin of 'LaB6'"),
        (
            [{"occupancies": [{"phase": "LaB6", "sites": ["La"], "elements": ["La"]}]}],
            "two or more",
        ),
        (
            [
                {
                    "occupancies": [
                        {"phase": "LaB6", "sites": ["La", "B"], "elements": []}
                    ]
                }
            ],
            "distinct elements",
        ),
        (
            [{"occupancies": [{"phase": "LaB6", "sites": ["La", "B"]}]}],
            "each occupancies entry",
        ),
        ([{"hold_if": ["uiso"]}], r"unknown keys \['hold_if'\]"),
    ],
)
def test_bad_coordinates_origin_occupancies_rejected(stages, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.accumulate_stages(stages)


def test_coordinates_by_site_hold_what_is_not_named() -> None:
    project, phase = _shared_site_project()
    stages = [
        {"name": "some", "coordinates": {"LaB6": {"Ce": "y", "B": "xz"}}},
        {"name": "more", "coordinates": {"LaB6": {"La": "z"}}},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert result["completed"], result["stages"]
    first, second = result["stages"]
    # On m, y follows x: naming y frees x, the first of the two, on both
    # atoms of the site, and z is held on each; B alone is held in y.
    assert first["coordinates"] == {
        "LaB6": {
            "La": {"site": "La/Ce", "refined": "x", "held": "z", "origin": False},
            "B": {"site": "B", "refined": "xz", "held": "y", "origin": False},
        }
    }
    assert project.constraints[:4] == [
        ["0::dAx:0", "0::dAx:2"],
        ["0::dAz:0"],
        ["0::dAz:2"],
        ["0::dAy:1"],
    ]
    assert ({"Atoms": {"all": "X"}},) in _set_calls(phase.calls)
    # The next stage frees z on the site: its holds go, a tie comes.
    assert second["coordinates"]["LaB6"]["La"]["refined"] == "xz"
    assert project.phase_constraints() == [
        ["0::dAx:0", "0::dAx:2"],
        ["0::dAz:0", "0::dAz:2"],
        ["0::dAy:1"],
    ]
    kinds = [c[-1] for c in project.data["Constraints"]["data"]["Phase"]]
    assert kinds == ["e", "e", "h"]


def test_a_coordinate_the_site_fixes_cannot_be_refined() -> None:
    project, phase = _shared_site_project()
    for atom in phase.data["Atoms"]:
        atom[7] = "4"
    stages = [{"coordinates": {"LaB6": {"La": "x"}}}]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert "fixes x" in result["stages"][0]["error"]
    # The stage itself refined nothing; the one refinement is the no cycle
    # one that computes the starting model the run falls back to.
    assert project.refinements == 1
    assert result["final_from"] == "job start"


# P4: the identity and the fourfold about c, which leaves z alone, so that
# the origin floats along c.
P4_SYMMETRY = {
    "SGOps": [
        [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0]],
        [[[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0.0, 0.0, 0.0]],
    ],
    "SGCen": [[0, 0, 0]],
    "SGInv": False,
}


def _polar_project():
    project, phase = _shared_site_project()
    phase.data["General"]["SGData"] = P4_SYMMETRY
    return project, phase


def test_origin_is_held_and_a_polar_axis_needs_one() -> None:
    every = {"LaB6": {"La": "all", "B": "all"}}
    project, phase = _polar_project()
    assert driver.polar_axes(phase) == [2]

    result = driver.refine(
        _g2sc(project), {"gpx": "fake.gpx", "stages": [{"coordinates": every}]}
    )
    assert "polar along z" in result["stages"][0]["error"]

    project, _ = _polar_project()
    stages = [{"coordinates": every, "origin": {"LaB6": {"site": "B", "axis": "z"}}}]
    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert result["completed"], result["stages"]
    assert result["stages"][0]["coordinates"]["LaB6"]["B"] == {
        "site": "B",
        "refined": "xy",
        "held": "z",
        "origin": True,
    }
    assert ["0::dAz:1"] in project.phase_constraints()


@pytest.mark.parametrize(
    ("stage", "message"),
    [
        (
            {
                "coordinates": {"LaB6": {"B": "z"}},
                "origin": {"LaB6": {"site": "B", "axis": "z"}},
            },
            "holds the origin",
        ),
        (
            {
                "coordinates": {"LaB6": {"B": "x"}},
                "origin": {"LaB6": {"site": "Xx", "axis": "z"}},
            },
            "no atoms",
        ),
        ({"coordinates": {"LaB6": {"Xx": "x"}}}, "no atoms"),
    ],
)
def test_bad_coordinates_fail_the_stage(stage, message) -> None:
    project, _ = _polar_project()
    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": [stage]})

    assert message in result["stages"][0]["error"]


def test_uiso_groups_take_in_every_atom_on_their_sites() -> None:
    project, phase = _shared_site_project()
    job = {
        "gpx": "fake.gpx",
        "stages": [{"uiso_groups": {"LaB6": [["Ce"], ["B"]]}}],
        "overall_uiso_start": {"*": 0.01},
    }

    result = driver.refine(_g2sc(project), job)

    assert result["completed"], result["stages"]
    # Ce names its site, which La shares.
    assert project.phase_constraints() == [["0::AUiso:0", "0::AUiso:2"]]
    assert ({"Atoms": {"all": "U"}},) in _set_calls(phase.calls)

    project, _ = _shared_site_project()
    job["stages"] = [{"uiso_groups": {"LaB6": [["La"], ["Ce", "B"]]}}]
    result = driver.refine(_g2sc(project), job)
    assert "more than one Uiso group" in result["stages"][0]["error"]


def _two_site_project():
    """Sr and Ba on two sites, 2a and 4c as in P4bm, with La on the second."""
    project = FakeProject()
    (phase,) = project.phases_

    def atom(label, element, xyz, occupancy, site, multiplicity):
        record = _atom(label, element, xyz, uiso=0.01)
        record[6:9] = [occupancy, site, multiplicity]
        return record

    phase.data["Atoms"] = [
        atom("Sr1", "Sr", [0.0, 0.0, 0.5], 0.6, "4", 2),
        atom("Ba1", "Ba", [0.0, 0.0, 0.5], 0.1, "4", 2),
        atom("Ba2", "Ba+2", [0.17, 0.67, 0.5], 0.6, "m", 4),
        atom("Sr2", "Sr", [0.17, 0.67, 0.5], 0.2, "m", 4),
        atom("La2", "La", [0.17, 0.67, 0.5], 0.05, "m", 4),
    ]
    return project, phase


def test_occupancies_hold_each_element_content_over_the_sites() -> None:
    project, phase = _two_site_project()
    entry = {"phase": "LaB6", "sites": ["Sr1", "Ba2"], "elements": ["Sr", "Ba"]}
    stages = [
        {"name": "exchange", "occupancies": [entry]},
        {"name": "again", "scale": True},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert result["completed"], result["stages"]
    # One equation for each element, over its atoms on both sites, weighted
    # by multiplicity; set up once, the second stage asking the same.
    assert project.equations == [
        (pytest.approx(2 * 0.6 + 4 * 0.2), ["0::Afrac:0", "0::Afrac:3"], [2.0, 4.0]),
        (pytest.approx(2 * 0.1 + 4 * 0.6), ["0::Afrac:1", "0::Afrac:2"], [2.0, 4.0]),
    ]
    # La, on one of the sites but not named, keeps its occupancy.
    assert ({"Atoms": {"Sr1": "F", "Ba1": "F", "Ba2": "F", "Sr2": "F"}},) in _set_calls(
        phase.calls
    )
    exchange, again = result["stages"]
    (sr, ba) = exchange["occupancy_constraints"]["LaB6"]
    assert sr["element"] == "Sr" and sr["atoms"] == ["Sr1", "Sr2"]
    assert sr["total"] == pytest.approx(2.0)
    assert ba["atoms"] == ["Ba1", "Ba2"]
    assert again["occupancy_constraints"] == exchange["occupancy_constraints"]


@pytest.mark.parametrize(
    ("entry", "message"),
    [
        ({"sites": ["Sr1", "Ba2"], "elements": ["Sr", "Ca"]}, "no Ca on the site of"),
        ({"sites": ["Sr1", "Ba1"], "elements": ["Sr"]}, "one site twice"),
        ({"sites": ["Sr1", "Xx"], "elements": ["Sr"]}, "no atoms"),
    ],
)
def test_bad_occupancy_constraints_fail_the_stage(entry, message) -> None:
    project, _ = _two_site_project()
    stages = [{"occupancies": [{"phase": "LaB6", **entry}]}]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert message in result["stages"][0]["error"]
    assert project.equations == []


def test_check_sanity_flags_uiso_occupancy_and_shifts() -> None:
    def atom(label, xyz, occupancy, uiso=0.01, uij=None):
        return {
            "label": label,
            "xyz": xyz,
            "occupancy": occupancy,
            "adp": "I" if uij is None else "A",
            "uiso": uiso if uij is None else None,
            "uij": uij,
        }

    atoms = [
        atom("A", [0.0, 0.0, 0.52], 0.7, uiso=-0.002),
        atom("A2", [0.0, 0.0, 0.52], 0.4),
        atom("B", [0.98, 0.5, 0.0], 1.2, uij=[0.01, -0.001, 0.02, 0, 0, 0]),
        atom("C", [0.25, 0.25, 0.25], 1.0),
    ]
    reference = [
        {"label": "A", "xyz": [0.0, 0.0, 0.5]},
        {"label": "B", "xyz": [0.03, 0.5, 0.0]},
        {"label": "C", "xyz": [0.25, 0.25, 0.25]},
    ]

    flags = driver.check_sanity(atoms, reference, max_shift=0.01)

    # B has moved across the cell edge by 0.05, not by 0.95.
    assert [(flag["kind"], flag["message"]) for flag in flags] == [
        ("negative_uiso", "A Uiso -0.0020"),
        ("negative_uiso", "B U22 -0.0010"),
        ("occupancy", "B occupancy 1.2000"),
        ("occupancy", "site A/A2 occupancy 1.1000"),
        ("coordinate_shift", "A/A2 z moved +0.0200"),
        ("coordinate_shift", "B x moved -0.0500"),
    ]
    assert flags[4]["atoms"] == ["A", "A2"]
    assert len(driver.check_sanity(atoms)) == 4, "no shifts without a reference"


def _uiso_going_negative():
    """A project whose refinement leaves La with a negative Uiso; one of no
    cycles, which computes the pattern rather than refining, leaves it be."""
    project = FakeProject()
    (phase,) = project.phases_
    original = project.refine

    def refine():
        original()
        if project.controls.get("cycles") != 0:
            phase.data["Atoms"][0][10] = -0.004

    project.refine = refine
    return project


def test_a_flagged_stage_kept_carries_its_flag_on() -> None:
    """With on_flagged accept, the flag stands and the run carries on."""
    project = _uiso_going_negative()
    stages = [
        {"name": "Uiso", "scale": True},
        {"name": "after", "cell": True},
    ]

    job = {"gpx": "fake.gpx", "stages": stages, "on_flagged": "accept"}
    result = driver.refine(_g2sc(project), job)

    assert result["completed"]
    assert project.refinements == 2
    first, after = result["stages"]
    assert first["status"] == "flagged"
    assert first["sanity"] == [
        {
            "kind": "negative_uiso",
            "atoms": ["La"],
            "value": -0.004,
            "message": "La Uiso -0.0040",
            "phase": "LaB6",
        }
    ]
    # The same flag is not a new one, so the stage after it is clean.
    assert after["status"] == "clean"
    assert result["final"]["phases"][0]["atoms"][0]["uiso"] == -0.004
    assert result["final_from"] == "after"
    assert result["sanity"] == {
        "max_shift": driver.DEFAULT_MAX_SHIFT,
        "reference": "job start",
    }


def test_sanity_reference_and_shift_from_the_job() -> None:
    project = FakeProject()
    job = {
        "gpx": "fake.gpx",
        "stages": [{"scale": True}],
        "sanity": {
            "max_shift": 0.01,
            "reference": {"LaB6": [{"label": "B", "xyz": [0.5, 0.5, 0.25]}]},
        },
    }

    result = driver.refine(_g2sc(project), job)

    (stage,) = result["stages"]
    assert [flag["message"] for flag in stage["sanity"]] == ["B z moved -0.0500"]
    assert result["sanity"] == {"max_shift": 0.01, "reference": "given"}


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"max_shift": 0.0}, "must be positive"),
        ({"shift": 0.1}, 'only "max_shift"'),
        ({"reference": {"LaB6": [{"label": "B"}]}}, "must be a list"),
    ],
)
def test_bad_sanity_settings_rejected(settings, message) -> None:
    with pytest.raises(ValueError, match=message):
        driver.check_sanity_settings(settings)


# Rejecting a stage and going back to the last one kept


def _oxygen_moving_project():
    """A project whose second refinement moves B's z by 0.089, as the
    x = 0.10 oxygen stage moved O3's, and whose later ones leave it be."""
    project = FakeProject()
    (phase,) = project.phases_
    original = project.refine

    def refine():
        original()
        if project.refinements == 2:
            phase.data["Atoms"][1][5] = 0.289

    project.refine = refine
    return project


def _b_z(atoms):
    return next(atom["xyz"][2] for atom in atoms if atom["label"] == "B")


def test_a_flagged_stage_is_rejected_and_rolled_back() -> None:
    project = _oxygen_moving_project()
    stages = [
        {"name": "A site", "scale": True},
        {"name": "oxygen", "zero": True},
        {"name": "occupancies", "cell": True},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert result["completed"]
    assert result["on_flagged"] == "reject"
    assert result["rejected"] == ["oxygen"]
    a_site, oxygen, occupancies = result["stages"]
    assert [stage["status"] for stage in result["stages"]] == [
        "clean",
        "rejected",
        "clean",
    ]
    # The rejected stage is recorded in full, with why and the atoms as it
    # left them.
    assert oxygen["rejected_because"] == ["sanity check: B z moved +0.0890"]
    assert [flag["kind"] for flag in oxygen["sanity"]] == ["coordinate_shift"]
    assert _b_z(oxygen["atoms"]["LaB6"]) == pytest.approx(0.289)
    assert "parameters" in oxygen and "rwp" in oxygen
    # The stage after it starts from where the A site stage left the atoms.
    assert _b_z(a_site["atoms"]["LaB6"]) == pytest.approx(0.2)
    assert _b_z(occupancies["atoms"]["LaB6"]) == pytest.approx(0.2)
    assert _b_z(result["final"]["phases"][0]["atoms"]) == pytest.approx(0.2)
    # What the rejected stage alone refined, the zero, is held from then on,
    # while the stages before and after it keep theirs.
    assert occupancies["flags"]["zero"] is False
    assert occupancies["flags"]["scale"] is True
    assert occupancies["flags"]["cell"]
    assert project.refinements == 3


def test_a_flagged_stage_is_kept_when_the_policy_accepts() -> None:
    project = _oxygen_moving_project()
    stages = [
        {"name": "A site", "scale": True},
        {"name": "oxygen", "zero": True},
        {"name": "occupancies", "cell": True},
    ]
    job = {"gpx": "fake.gpx", "stages": stages, "on_flagged": "accept"}

    result = driver.refine(_g2sc(project), job)

    assert result["rejected"] == []
    # The last stage carries the same flag, not a new one, so it is clean.
    assert [stage["status"] for stage in result["stages"]] == [
        "clean",
        "flagged",
        "clean",
    ]
    # Nothing is rolled back: the move is carried into the stages after it.
    assert _b_z(result["stages"][2]["atoms"]["LaB6"]) == pytest.approx(0.289)
    assert result["stages"][2]["flags"]["zero"] is True


def test_a_flag_the_job_started_with_rejects_nothing() -> None:
    project = FakeProject()
    (phase,) = project.phases_
    phase.data["Atoms"][0][10] = -0.004

    result = driver.refine(
        _g2sc(project),
        {"gpx": "fake.gpx", "stages": [{"name": "scale", "scale": True}]},
    )

    (stage,) = result["stages"]
    assert stage["status"] == "clean"
    assert [flag["message"] for flag in stage["sanity"]] == ["La Uiso -0.0040"]
    assert result["rejected"] == []


def test_every_stage_rejected_leaves_the_run_where_it_started() -> None:
    project = _uiso_going_negative()

    result = driver.refine(
        _g2sc(project),
        {"gpx": "fake.gpx", "stages": [{"name": "Uiso", "scale": True}]},
    )

    assert result["completed"]
    assert result["rejected"] == ["Uiso"]
    assert result["stages"][0]["status"] == "rejected"
    assert result["undetermined"] == [], "no stage was kept"
    # The project is back where the job found it, and the final model is
    # that starting one, computed by a refinement of no cycles.
    assert project.phases_[0].data["Atoms"][0][10] == pytest.approx(0.0086)
    assert result["final_from"] == "job start"
    assert result["final"]["phases"][0]["atoms"][0]["uiso"] == pytest.approx(0.0086)
    assert project.controls["cycles"] == 0
    assert project.refinements == 2


def test_a_failed_stage_goes_back_to_the_last_one_kept() -> None:
    project = FakeProject(fail=(2,))
    (phase,) = project.phases_
    original = project.refine

    def refine():
        original()
        if project.refinements == 1:
            phase.data["Atoms"][1][5] = 0.21

    project.refine = refine
    stages = [
        {"name": "first", "scale": True},
        {"name": "second", "zero": True},
        {"name": "third", "cell": True},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    assert not result["completed"]
    assert [stage["status"] for stage in result["stages"]] == ["clean", "failed"]
    assert result["stages"][1]["error"].startswith("RefinementError:")
    assert project.refinements == 2, "the run stops at the failed stage"
    # The final values are those of the last stage kept.
    assert _b_z(result["final"]["phases"][0]["atoms"]) == pytest.approx(0.21)


# Unsettled stages


def test_an_unsettled_stage_is_kept_by_default_with_its_largest_move() -> None:
    project = _settling_project()
    job = {
        "gpx": "fake.gpx",
        "max_passes": 3,
        "stages": [{"name": "A site", "scale": True, "zero": True}],
    }

    result = driver.refine(SimpleNamespace(G2Project=lambda gpx: project), job)

    assert result["on_unsettled"] == "accept"
    (stage,) = result["stages"]
    assert stage["status"] == "unsettled"
    assert not stage["passes_converged"]
    assert stage["largest_remaining_move"] == {
        "parameter": ":0:Scale",
        "shift_over_esd": pytest.approx(10.0),
    }
    assert result["rejected"] == []


def _settling_after(refinements):
    """A project whose scale moves 10 esds a pass until it has been refined
    ``refinements`` times, and stays put from then on."""
    project = FakeProject()
    original = project.refine

    def refine():
        original()
        covariance = project.data["Covariance"]["data"]
        covariance["variables"][0] = 2.0 * min(project.refinements, refinements)

    project.refine = refine
    return project


def test_an_unsettled_stage_is_rejected_when_the_policy_says_so() -> None:
    project = _settling_after(3)
    job = {
        "gpx": "fake.gpx",
        "max_passes": 3,
        "on_unsettled": "reject",
        "stages": [
            {"name": "A site", "scale": True, "zero": True},
            {"name": "oxygen", "cell": True},
        ],
    }

    result = driver.refine(SimpleNamespace(G2Project=lambda gpx: project), job)

    assert result["rejected"] == ["A site"]
    assert [stage["status"] for stage in result["stages"]] == ["rejected", "clean"]
    a_site, oxygen = result["stages"]
    assert a_site["status"] == "rejected"
    assert a_site["rejected_because"] == [
        "not settled in 3 passes: :0:Scale still moved 10.00 esd in the last"
    ]
    # Its largest remaining move is recorded all the same.
    assert a_site["largest_remaining_move"]["shift_over_esd"] == pytest.approx(10.0)
    # What it alone refined is held in the stage after it.
    assert oxygen["flags"]["scale"] is False
    assert oxygen["flags"]["zero"] is False
    assert oxygen["flags"]["cell"]


def test_a_settled_stage_is_clean() -> None:
    project = _settling_project()
    job = {
        "gpx": "fake.gpx",
        "max_passes": 20,
        "on_unsettled": "reject",
        "stages": [{"name": "zero", "scale": True, "zero": True, "le_bail": True}],
    }

    result = driver.refine(_le_bail_g2sc(project, []), job)

    (stage,) = result["stages"]
    assert stage["status"] == "clean"
    assert "largest_remaining_move" not in stage


def test_policies_default_and_anything_else_is_rejected() -> None:
    assert driver.check_policies({}) == ("reject", "accept")
    assert driver.check_policies({"on_flagged": "accept"}) == ("accept", "accept")
    assert driver.check_policies({"on_unsettled": "reject"}) == ("reject", "reject")
    with pytest.raises(ValueError, match="on_unsettled must be one of accept, reject"):
        driver.check_policies({"on_unsettled": "roll back"})
    with pytest.raises(ValueError, match="on_flagged must be one of accept, reject"):
        driver.check_policies({"on_flagged": "hold"})


# Parameters the data leave undetermined


def _atoms_with_esds(**changes):
    """An atom table of one atom, as _atom_table gives it, with esds."""
    atom = {
        "label": "Sr1",
        "type": "Sr",
        "xyz": [0.0, 0.0, 0.5],
        "occupancy": 0.4,
        "site_symmetry": "4",
        "multiplicity": 2,
        "adp": "I",
        "uiso": 0.01,
        "uij": None,
        "uiso_esd": None,
        "occupancy_esd": None,
        "xyz_esd": [None, None, None],
    }
    return [{**atom, **changes}]


def test_find_undetermined_occupancy_over_half_its_range() -> None:
    # The Sr and Ba split of 0.1(6): an esd larger than half of 0 to 1.
    atoms = _atoms_with_esds(occupancy=0.1, occupancy_esd=0.6)

    (found,) = driver.find_undetermined(atoms, _atoms_with_esds())

    assert found == {
        "atom": "Sr1",
        "parameter": "occupancy",
        "value": 0.1,
        "esd": 0.6,
        "message": "Sr1 occupancy 0.100, esd 0.600 more than half its range",
    }
    # Half the range exactly is not more than half of it.
    assert driver.find_undetermined(_atoms_with_esds(occupancy_esd=0.5)) == []
    # An occupancy not refined has no esd and is never undetermined.
    assert driver.find_undetermined(_atoms_with_esds()) == []


def test_find_undetermined_coordinate_and_uiso_against_their_shift() -> None:
    start = _atoms_with_esds()
    # z moved 0.002 with an esd of 0.005, and Uiso 0.004 with an esd of 0.001.
    atoms = _atoms_with_esds(
        xyz=[0.0, 0.0, 0.502],
        xyz_esd=[None, None, 0.005],
        uiso=0.014,
        uiso_esd=0.001,
    )

    found = driver.find_undetermined(atoms, start)

    assert [entry["parameter"] for entry in found] == ["z"]
    assert found[0]["message"] == (
        "Sr1 z 0.50200, esd 0.00500 more than its shift 0.00200"
    )
    # Moved further than its esd, it is determined.
    moved = _atoms_with_esds(xyz=[0.0, 0.0, 0.51], xyz_esd=[None, None, 0.005])
    assert driver.find_undetermined(moved, start) == []
    # A Uiso whose esd is larger than its shift is undetermined too.
    barely = _atoms_with_esds(uiso=0.0105, uiso_esd=0.001)
    (entry,) = driver.find_undetermined(barely, start)
    assert entry["parameter"] == "uiso"


def test_find_undetermined_wraps_a_coordinate_over_a_cell() -> None:
    start = _atoms_with_esds(xyz=[0.0, 0.0, 0.999])
    # 0.001 the other side of the cell edge: a shift of 0.002, not 0.998.
    atoms = _atoms_with_esds(xyz=[0.0, 0.0, 0.001], xyz_esd=[None, None, 0.005])

    (entry,) = driver.find_undetermined(atoms, start)

    assert entry["parameter"] == "z"
    assert "shift 0.00200" in entry["message"]


def test_an_atom_not_in_the_start_is_judged_on_its_occupancy_alone() -> None:
    atoms = _atoms_with_esds(
        label="Ba1",
        occupancy_esd=0.6,
        xyz_esd=[None, None, 0.5],
        uiso_esd=0.5,
    )

    found = driver.find_undetermined(atoms, _atoms_with_esds())

    assert [entry["parameter"] for entry in found] == ["occupancy"]


def test_undetermined_parameters_are_recorded_with_each_stage() -> None:
    """A stage that refines an occupancy the data do not determine."""
    project = FakeProject()
    (phase,) = project.phases_
    original = project.refine

    def refine():
        original()
        covariance = project.data["Covariance"]["data"]
        covariance["varyList"] = [":0:Scale", "0::Afrac:0"]
        covariance["variables"] = [2.0, 0.1]
        covariance["covMatrix"] = [[0.04, 0.0], [0.0, 0.36]]
        phase.data["Atoms"][0][6] = 0.1

    project.refine = refine
    stages = [{"name": "A site", "scale": True}]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    (stage,) = result["stages"]
    assert stage["status"] == "clean"
    assert stage["undetermined"] == [
        {
            "atom": "La",
            "parameter": "occupancy",
            "value": 0.1,
            "esd": 0.6,
            "message": "La occupancy 0.100, esd 0.600 more than half its range",
            "phase": "LaB6",
        }
    ]
    # The result's own are those of the last stage kept.
    assert result["undetermined"] == stage["undetermined"]


def test_a_rejected_stage_keeps_its_undetermined_out_of_the_result() -> None:
    project = FakeProject()
    (phase,) = project.phases_
    original = project.refine

    def refine():
        original()
        covariance = project.data["Covariance"]["data"]
        if project.refinements == 2:
            covariance["varyList"] = [":0:Scale", "0::Afrac:0"]
            covariance["variables"] = [2.0, 1.6]
            covariance["covMatrix"] = [[0.04, 0.0], [0.0, 0.36]]
            # Past 1, so the stage is flagged and rejected.
            phase.data["Atoms"][0][6] = 1.6

    project.refine = refine
    stages = [
        {"name": "A site", "scale": True},
        {"name": "occupancies", "zero": True},
    ]

    result = driver.refine(_g2sc(project), {"gpx": "fake.gpx", "stages": stages})

    a_site, occupancies = result["stages"]
    assert occupancies["status"] == "rejected"
    assert [entry["message"] for entry in occupancies["undetermined"]] == [
        "La occupancy 1.600, esd 0.600 more than half its range"
    ]
    # The result's are the A site stage's, the last kept: none.
    assert a_site["undetermined"] == []
    assert result["undetermined"] == []


def test_the_starting_model_is_computed_not_refined() -> None:
    """The no cycle refinement that gives a run with nothing kept a model."""
    project = _uiso_going_negative()
    job = {"gpx": "fake.gpx", "cycles": 7, "stages": [{"name": "Uiso", "scale": True}]}

    result = driver.refine(_g2sc(project), job)

    assert result["final_from"] == "job start"
    assert "compute_error" not in result and "final_error" not in result
    # The cycles the job asked for are put back after it.
    assert project.controls["cycles"] == 7
    # The model has no esd, nothing having been refined into it.
    atom = result["final"]["phases"][0]["atoms"][0]
    assert atom["uiso"] == pytest.approx(0.0086)
    assert atom["uiso_esd"] is None and atom["xyz_esd"] == [None, None, None]


def test_a_starting_model_that_cannot_be_computed_is_recorded() -> None:
    project = _uiso_going_negative()
    # The first refinement is the rejected stage's, the second the compute.
    project.fail = {2}

    result = driver.refine(
        _g2sc(project),
        {"gpx": "fake.gpx", "stages": [{"name": "Uiso", "scale": True}]},
    )

    assert result["completed"], "the stages all ran"
    assert result["rejected"] == ["Uiso"]
    assert result["compute_error"].startswith("RefinementError:")
    # The values of the starting model are still reported, without the
    # residuals the computation would have given.
    assert result["final_from"] == "job start"
    assert result["final"]["phases"][0]["atoms"][0]["uiso"] == pytest.approx(0.0086)


def test_a_cap_of_one_pass_still_records_it() -> None:
    """A run capped at one pass records that pass and that it did not settle,
    so that a result has the same shape whatever the cap."""
    project = _settling_project()
    job = {
        "gpx": "fake.gpx",
        "max_passes": 1,
        "stages": [{"name": "zero", "zero": True}],
    }

    result = driver.refine(SimpleNamespace(G2Project=lambda gpx: project), job)

    (stage,) = result["stages"]
    assert project.refinements == 1
    assert len(stage["passes"]) == 1
    assert stage["passes"][0]["max_shift_over_esd"] is None
    assert not stage["passes_converged"]
    assert stage["status"] == "unsettled"
    assert stage["largest_remaining_move"] == {
        "parameter": None,
        "shift_over_esd": None,
    }


class StalePatternProject(FakeProject):
    """Refines as GSAS-II does a stage of constrained coordinates whose last
    function evaluation was a trial step it refused: the covariance's Rwp is
    the model's, the histogram's residuals the trial's, until the pattern is
    computed again with no cycles."""

    def refine(self):
        super().refine()
        covariance = self.data["Covariance"]["data"]
        if self.controls.get("cycles") == 0:
            self.histogram_.residuals = {"wR": 4.0, "R": 3.0}
        elif "Rvals" in covariance:
            covariance["Rvals"]["Rwp"] = 4.0
            self.histogram_.residuals = {"wR": 100.0, "R": 56.0}


def test_a_stale_pattern_is_recorded_from_the_accepted_fit_and_computed_again() -> None:
    project = StalePatternProject()
    G2sc = SimpleNamespace(G2Project=lambda gpx: project)
    job = {
        "gpx": "fake.gpx",
        "cycles": 5,
        "stages": [{"name": "zero", "zero": True}],
    }

    result = driver.refine(G2sc, job)

    (stage,) = result["stages"]
    # The Rwp of the accepted chi squared, and no Rp, which that does not give.
    assert (stage["rwp"], stage["rp"]) == (4.0, None)
    assert stage["residuals_stale"] is True
    # The pattern was computed again before the exports, and the refinement's
    # covariance kept.
    assert project.histogram_.residuals == {"wR": 4.0, "R": 3.0}
    assert result["final"]["instrument"]["Zero"]["esd"] == pytest.approx(1.0e-3)


def test_stale_compares_the_histogram_with_the_accepted_rwp() -> None:
    histogram = FakeHistogram()

    assert not driver._stale(histogram, {"Rvals": {"GOF": 1.5}})
    assert not driver._stale(histogram, {"Rvals": {"Rwp": 12.5001}})
    assert driver._stale(histogram, {"Rvals": {"Rwp": 4.0}})
