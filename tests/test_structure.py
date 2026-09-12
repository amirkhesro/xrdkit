"""Tests for xrdkit.structure."""

import math

import numpy as np
import pytest

from xrdkit.structure import (
    bond_lengths,
    cell_contents,
    composition_edits,
    interatomic_distances,
    metric_tensor,
    site_setup,
)

IDENTITY = {"rotation": np.eye(3).tolist(), "translation": [0.0, 0.0, 0.0]}
# P4: the fourfold about c, x, y, z -> -y, x, z, and its powers.
P4 = [
    {
        "rotation": np.linalg.matrix_power(
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]], n
        ).tolist(),
        "translation": [0.0, 0.0, 0.0],
    }
    for n in range(4)
]


def test_metric_tensor_hexagonal() -> None:
    metric = metric_tensor([3.0, 3.0, 5.0, 90.0, 90.0, 120.0])
    for vector, length in (([1, 0, 0], 3.0), ([1, 1, 0], 3.0), ([0, 0, 1], 5.0)):
        v = np.array(vector, dtype=float)
        assert math.sqrt(v @ metric @ v) == pytest.approx(length)


def test_distances_to_both_neighbouring_images_with_esd() -> None:
    atoms = [
        {"label": "A", "xyz": [0.0, 0.0, 0.0]},
        {"label": "B", "xyz": [0.5, 0.0, 0.0], "xyz_esd": [0.001, None, None]},
    ]

    found = interatomic_distances(
        [4.0, 4.0, 4.0, 90, 90, 90], atoms, [IDENTITY], ["A"], ["B"], dmax=2.5
    )

    # B at +x and at -x, a / 2 away; an x esd of 0.001 is 0.004 A along x.
    assert [d.distance for d in found] == pytest.approx([2.0, 2.0])
    assert [d.esd for d in found] == pytest.approx([0.004, 0.004])
    assert sorted(d.translation for d in found) == [(-1, 0, 0), (0, 0, 0)]
    assert all(d.centre == "A" and d.target == "B" for d in found)


def test_symmetry_images_counted_once() -> None:
    atoms = [
        {"label": "M", "xyz": [0.0, 0.0, 0.2], "xyz_esd": [None, None, 0.01]},
        {"label": "O", "xyz": [0.1, 0.0, 0.2]},
    ]
    cell = [5.0, 5.0, 3.0, 90, 90, 90]

    to_oxygen = interatomic_distances(cell, atoms, P4, ["M"], ["O"], dmax=1.0, dmin=0.1)
    # The fourfold puts O at four places 0.5 A round M.
    assert [d.distance for d in to_oxygen] == pytest.approx([0.5] * 4)
    assert len({d.operator for d in to_oxygen}) == 4

    to_itself = interatomic_distances(cell, atoms, P4, ["M"], ["M"], dmax=3.5)
    # M lies on the axis, where every operator gives the same image: the two
    # along c, 3 A away, each once, and their length does not depend on z.
    assert [d.distance for d in to_itself] == pytest.approx([3.0, 3.0])
    assert [d.esd for d in to_itself] == pytest.approx([0.0, 0.0])


def test_no_esd_without_coordinate_esds() -> None:
    atoms = [{"label": "A", "xyz": [0, 0, 0]}, {"label": "B", "xyz": [0.25, 0.25, 0]}]
    (found,) = interatomic_distances(
        [4.0, 4.0, 10.0, 90, 90, 90], atoms, [IDENTITY], ["A"], ["B"], dmax=1.5
    )
    assert found.distance == pytest.approx(math.sqrt(2.0))
    assert found.esd is None


def test_distances_rejects_unknown_labels_and_bad_range() -> None:
    atoms = [{"label": "A", "xyz": [0, 0, 0]}]
    with pytest.raises(ValueError, match="no atoms labelled"):
        interatomic_distances(
            [4, 4, 4, 90, 90, 90], atoms, [IDENTITY], ["A"], ["X"], 3.0
        )
    with pytest.raises(ValueError, match="dmax must exceed"):
        interatomic_distances(
            [4, 4, 4, 90, 90, 90], atoms, [IDENTITY], ["A"], ["A"], 0.4
        )


# A cubic perovskite written out in P1, with Sr sharing the A site with Ba,
# as the GSAS-II driver reports a phase.
PEROVSKITE = {
    "cell": {
        "length_a": 4.0,
        "length_b": 4.0,
        "length_c": 4.0,
        "angle_alpha": 90.0,
        "angle_beta": 90.0,
        "angle_gamma": 90.0,
        "volume": 64.0,
    },
    "atoms": [
        {"label": "Ba", "type": "Ba+2", "xyz": [0.0, 0.0, 0.0]},
        {"label": "Sr", "type": "Sr", "xyz": [0.0, 0.0, 0.0]},
        {"label": "Ti", "type": "Ti", "xyz": [0.5, 0.5, 0.5]},
        {
            "label": "O1",
            "type": "O-2",
            "xyz": [0.5, 0.5, 0.0],
            "xyz_esd": [None, None, 0.001],
        },
        {"label": "O2", "type": "O", "xyz": [0.5, 0.0, 0.5]},
        {"label": "O3", "type": "O", "xyz": [0.0, 0.5, 0.5]},
    ],
    "operators": [IDENTITY],
}


def test_bond_lengths_by_site_to_the_anions() -> None:
    bonds = bond_lengths(PEROVSKITE)

    # Twelve A-O bonds of a / sqrt(2) and six B-O of a / 2, each oxygen
    # counted by its images; Sr shares the A site and counts once, and no
    # O-O distance is taken.
    assert [
        (bond["centre"], bond["centre_site"], bond["target"], bond["count"])
        for bond in bonds
    ] == [
        ("Ba", "Ba/Sr", "O1", 4),
        ("Ba", "Ba/Sr", "O2", 4),
        ("Ba", "Ba/Sr", "O3", 4),
        ("Ti", "Ti", "O1", 2),
        ("Ti", "Ti", "O2", 2),
        ("Ti", "Ti", "O3", 2),
    ]
    assert [bond["distance"] for bond in bonds] == pytest.approx(
        [2.0 * math.sqrt(2.0)] * 3 + [2.0] * 3
    )
    # Ti-O1 lies along c, where O1's z esd of 0.001 is 0.004 A.
    assert bonds[3]["esd"] == pytest.approx(0.004)
    assert bonds[4]["esd"] is None
    assert bonds[0]["target_site"] == "O1"


def test_bond_lengths_range_anions_and_a_plain_cell() -> None:
    plain = {**PEROVSKITE, "cell": [4.0, 4.0, 4.0, 90.0, 90.0, 90.0]}

    assert [bond["centre"] for bond in bond_lengths(plain, dmax=2.5)] == ["Ti"] * 3
    # With Ti as the anion, the oxygens count as cations, and the A site's
    # nearest Ti are a sqrt(3) / 2 away.
    to_ti = bond_lengths(plain, anions=["Ti"], dmax=3.5)
    assert {bond["centre"] for bond in to_ti} == {"Ba", "O1", "O2", "O3"}
    (to_ti,) = (bond for bond in to_ti if bond["centre"] == "Ba")
    assert to_ti["centre"] == "Ba" and to_ti["count"] == 8
    assert to_ti["distance"] == pytest.approx(2.0 * math.sqrt(3.0))
    with pytest.raises(ValueError, match="dmax must exceed"):
        bond_lengths(plain, dmax=0.4)


# A structure table as xrdkit.config.load_config gives it: two A sites, Sr
# and Ba sharing the first, one B site and two O sites.
STRUCTURE = {
    "formula_units": 1,
    "sites": [
        {
            "name": "Ba1",
            "atoms": {"Ba1": "Ba", "Sr1": "Sr"},
            "wyckoff": "1a",
            "kind": "A",
        },
        {"name": "Sr2", "atoms": {"Sr2": "Sr"}, "wyckoff": "1a", "kind": "A"},
        {"name": "Ti1", "atoms": {"Ti1": "Ti"}, "wyckoff": "1b", "kind": "B"},
        {"name": "O1", "atoms": {"O1": "O"}, "wyckoff": "1b", "kind": "O"},
        {"name": "O2", "atoms": {"O2": "O"}, "wyckoff": "2c", "kind": "O"},
    ],
    "free_coordinates": {"1a": "z", "2c": "xz"},
    "uiso_groups": [
        {"name": "A cations", "sites": ["Ba1", "Sr2"]},
        {"name": "B", "sites": ["Ti1"]},
        {"name": "O", "sites": ["O1", "O2"]},
    ],
    "origin": {"site": "Ti1", "axis": "z"},
    "exchange": {"elements": ["Sr", "Ba"], "sites": ["Ba1", "Sr2"]},
    "composition": {"added": {"La": "Sr", "Nb": "Ti"}},
}


def atom(label, kind, xyz, multiplicity, occupancy):
    return {
        "label": label,
        "type": kind,
        "xyz": xyz,
        "multiplicity": multiplicity,
        "occupancy": occupancy,
    }


ATOMS = [
    atom("Ba1", "Ba+2", [0.0, 0.0, 0.0], 1, 0.6),
    atom("Sr1", "Sr", [0.0, 0.0, 0.0], 1, 0.2),
    atom("Sr2", "Sr", [0.0, 0.0, 0.5], 1, 0.4),
    atom("Ti1", "Ti", [0.5, 0.5, 0.5], 1, 1.0),
    atom("O1", "O-2", [0.5, 0.5, 0.0], 1, 1.0),
    atom("O2", "O", [0.5, 0.0, 0.5], 2, 1.0),
]


def test_cell_contents() -> None:
    assert cell_contents(ATOMS) == pytest.approx(
        {"Ba": 0.6, "Sr": 0.6, "Ti": 1.0, "O": 3.0}
    )


def test_composition_edits_keep_the_distribution_and_add_by_host() -> None:
    target = {"Sr": 0.3, "Ba": 0.5, "La": 0.15, "Ti": 0.9, "Nb": 0.1, "O": 3.0}

    edits = composition_edits(ATOMS, target, STRUCTURE)

    # Sr halved on both its sites and La put where Sr was, a quarter of its
    # CIF occupancy; Ba scaled by 5/6; Ti at 0.9 and Nb at 0.1 on its site;
    # the oxygens, unchanged, get no edit.
    by_label = {edit["label"]: edit for edit in edits}
    assert [edit["label"] for edit in edits] == [
        "Ba1",
        "Sr1",
        "La1",
        "Sr2",
        "La2",
        "Ti1",
        "Nb1",
    ]
    assert by_label["Ba1"]["occupancy"] == pytest.approx(0.5)
    assert by_label["Sr1"]["occupancy"] == pytest.approx(0.1)
    assert by_label["Sr2"]["occupancy"] == pytest.approx(0.2)
    assert by_label["La1"] == {
        "label": "La1",
        "type": "La",
        "copy": "Sr1",
        "occupancy": pytest.approx(0.05),
    }
    assert by_label["La2"]["occupancy"] == pytest.approx(0.1)
    assert by_label["Nb1"]["occupancy"] == pytest.approx(0.1)

    # Applied, the edits give the target contents.
    edited = [dict(a) for a in ATOMS]
    for edit in edits:
        if "copy" in edit:
            source = next(a for a in edited if a["label"] == edit["copy"])
            edited.append({**source, "label": edit["label"], "type": edit["type"]})
        next(a for a in edited if a["label"] == edit["label"])["occupancy"] = edit[
            "occupancy"
        ]
    assert cell_contents(edited) == pytest.approx(target)

    # With no La, none is added.
    labels = [
        e["label"] for e in composition_edits(ATOMS, {**target, "La": 0.0}, STRUCTURE)
    ]
    assert "La1" not in labels and "La2" not in labels


def test_composition_edits_reject_a_mismatch() -> None:
    target = {"Sr": 0.3, "Ba": 0.5, "Ti": 1.0, "O": 3.0}
    with pytest.raises(ValueError, match="which the composition lacks"):
        composition_edits(
            ATOMS, {k: v for k, v in target.items() if k != "Ba"}, STRUCTURE
        )
    with pytest.raises(ValueError, match="the rule does not add"):
        composition_edits(ATOMS, {**target, "Ca": 0.1}, STRUCTURE)
    taken = ATOMS + [atom("La1", "Ba", [0.0, 0.0, 0.0], 1, 0.0)]
    with pytest.raises(ValueError, match="is there already|is taken"):
        composition_edits(taken, {**target, "La": 0.1}, STRUCTURE)


def test_site_setup_finds_the_sites() -> None:
    # La added on the Sr2 position joins that site.
    atoms = ATOMS + [atom("La2", "La", [0.0, 0.0, 0.5], 1, 0.1)]

    setup = site_setup(STRUCTURE, atoms)

    assert [site["name"] for site in setup["sites"]] == [
        "Ba1",
        "Sr2",
        "Ti1",
        "O1",
        "O2",
    ]
    assert [a["label"] for a in setup["sites"][1]["atoms"]] == ["Sr2", "La2"]
    assert [a["label"] for a in setup["sites"][0]["atoms"]] == ["Ba1", "Sr1"]
    assert [site["name"] for site in setup["kinds"]["A"]] == ["Ba1", "Sr2"]
    assert setup["uiso_groups"] == [["Ba1", "Sr2"], ["Ti1"], ["O1", "O2"]]
    assert setup["group_names"] == ["A cations", "B", "O"]
    assert setup["origin"]["name"] == "Ti1" and setup["origin_axis"] == "z"
    # The origin site is not refined; 1b is not named, so O1 takes "all".
    assert setup["coordinates"] == {
        "A": {"Ba1": "z", "Sr2": "z"},
        "B": {},
        "O": {"O1": "all", "O2": "xz"},
    }
    assert setup["sites"][4]["multiplicity"] == 2
    assert setup["exchange"] == {"elements": ["Sr", "Ba"], "sites": ["Ba1", "Sr2"]}


@pytest.mark.parametrize(
    ("atoms", "message"),
    [
        ([a for a in ATOMS if a["label"] != "O2"], "site O2: no atom labelled"),
        (
            [
                {**a, "xyz": [0.0, 0.0, 0.25]} if a["label"] == "Sr1" else a
                for a in ATOMS
            ],
            "site Ba1: atoms .* are not on one position",
        ),
        (
            [{**a, "type": "Ca"} if a["label"] == "Sr2" else a for a in ATOMS],
            "site Sr2: atom Sr2 is Ca, not Sr",
        ),
        (
            [{**a, "multiplicity": 4} if a["label"] == "O2" else a for a in ATOMS],
            "site O2: multiplicity 4, not 2 as on 2c",
        ),
        (
            ATOMS + [atom("O3", "O", [0.0, 0.5, 0.5], 2, 1.0)],
            r"atoms \['O3'\] are on no configured site",
        ),
        (
            [
                {**a, "xyz": [0.0, 0.0, 0.0]} if a["label"] == "Sr2" else a
                for a in ATOMS
            ],
            "sites Ba1 and Sr2 are on one position",
        ),
    ],
)
def test_site_setup_rejects_atoms_that_do_not_fit(atoms, message) -> None:
    with pytest.raises(ValueError, match=message):
        site_setup(STRUCTURE, atoms)
