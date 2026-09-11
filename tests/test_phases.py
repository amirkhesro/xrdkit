"""Tests for xrdkit.phases. The COD is never contacted: urlopen is replaced."""

import csv
import io
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlsplit

import numpy as np
import pytest

from xrdkit import (
    CodRecord,
    SimulatedReflection,
    cod_fetch,
    cod_search,
    match_candidate,
    phases,
    simulate_pattern,
    write_cif_index,
)
from xrdkit.peaks import KALPHA1_WAVELENGTH
from xrdkit.phases import CIF_INDEX_COLUMNS


def cod_entry(**fields: object) -> dict[str, object]:
    """One object of a COD JSON result, trimmed to the fields read, as strings."""
    entry: dict[str, object] = {
        "file": "2100720",
        "a": "12.4844",
        "b": "12.4844",
        "c": "3.9572",
        "alpha": "90",
        "beta": "90",
        "gamma": "90",
        "vol": "616.77",
        "sg": "P 4 b m",
        "sgNumber": "100",
        "formula": "- Ba0.52 Nb2 O6 Sr0.48 -",
        "authors": "Podlozhenov, S.; Graetsch, H. A.",
        "journal": "Acta Crystallographica Section B",
        "year": "2006",
        "volume": "62",
        "doi": "10.1107/S0108768106038869",
    }
    entry.update(fields)
    return entry


SEARCH_RESULT = [
    cod_entry(),
    cod_entry(file="2311739", sg="P 4 b m (a,b,2*c)", c="7.9006"),
    cod_entry(file="1537507", sg="P 63/m m c", sgNumber="194", c="15.38"),
    cod_entry(file="2103856", sg="X4bm", sgNumber=None, c="7.8698"),
]

CIF_TEXT = b"#---\n# COD entry\ndata_2100720\n_cell_length_a 12.4844\n"


class FakeResponse(io.BytesIO):
    """What urlopen hands back: a readable body usable in a with block."""


# What the fake COD answers, by URL without its query: a body or an error to raise.
ANSWERS: dict[str, bytes | Exception] = {}


@pytest.fixture
def requests(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Replace urlopen, answer from ANSWERS and record every URL asked for."""
    asked: list[str] = []

    def fake_urlopen(request, timeout=None):
        url = request.full_url
        asked.append(url)
        body = ANSWERS[url.split("?")[0]]
        if isinstance(body, Exception):
            raise body
        return FakeResponse(body)

    ANSWERS.clear()
    monkeypatch.setattr(phases, "urlopen", fake_urlopen)
    return asked


def answer(url: str, body: bytes | Exception) -> None:
    ANSWERS[url] = body


def query_of(url: str) -> dict[str, str]:
    return {key: values[0] for key, values in parse_qs(urlsplit(url).query).items()}


def search_answer(result: list[dict[str, object]]) -> None:
    answer(f"{phases.COD_URL}/result", json.dumps(result).encode())


def test_search_asks_for_an_exact_element_set(requests: list[str]) -> None:
    search_answer(SEARCH_RESULT)

    cod_search(["Sr", "Ba", "Nb", "O"])

    assert query_of(requests[0]) == {
        "format": "json",
        "el1": "Sr",
        "el2": "Ba",
        "el3": "Nb",
        "el4": "O",
        "strictmin": "4",
        "strictmax": "4",
    }


def test_search_without_exact_sets_no_element_count(requests: list[str]) -> None:
    search_answer([])

    cod_search(["Nb", "O"], exact=False)

    query = query_of(requests[0])
    assert "strictmin" not in query
    assert "strictmax" not in query


def test_search_passes_a_space_group_number_and_extras(requests: list[str]) -> None:
    search_answer([])

    cod_search(["Sr", "Nb", "O"], space_group=100, extra={"year": 2006, "el1": "Ba"})

    query = query_of(requests[0])
    assert query["space_group_number"] == "100"
    assert query["year"] == "2006"
    assert query["el1"] == "Ba"


def test_search_reads_every_field(requests: list[str]) -> None:
    search_answer(SEARCH_RESULT[:1])

    (record,) = cod_search(["Sr", "Ba", "Nb", "O"])

    assert record == CodRecord(
        cod_id="2100720",
        formula="Ba0.52 Nb2 O6 Sr0.48",
        space_group="P 4 b m",
        space_group_number=100,
        a=12.4844,
        b=12.4844,
        c=3.9572,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        volume=616.77,
        authors="Podlozhenov, S.; Graetsch, H. A.",
        journal="Acta Crystallographica Section B",
        year=2006,
        doi="10.1107/S0108768106038869",
    )


def test_search_takes_blank_fields_as_none(requests: list[str]) -> None:
    search_answer([cod_entry(year=None, doi=None, vol=None, sgNumber="", a=None)])

    (record,) = cod_search(["Sr", "Ba", "Nb", "O"])

    assert record.year is None
    assert record.doi is None
    assert record.volume is None
    assert record.space_group_number is None
    assert record.a is None


@pytest.mark.parametrize("symbol", ["P4bm", "P 4 b m", " P4bm "])
def test_search_matches_a_symbol_ignoring_spaces_and_setting(
    requests: list[str], symbol: str
) -> None:
    search_answer(SEARCH_RESULT)

    records = cod_search(["Sr", "Ba", "Nb", "O"], space_group=symbol)

    assert [record.cod_id for record in records] == ["2100720", "2311739"]
    assert "sg" not in query_of(requests[0])


def test_search_symbol_drops_an_origin_choice(requests: list[str]) -> None:
    search_answer([cod_entry(sg="R -3 :H")])

    assert len(cod_search(["Sr", "Nb", "O"], space_group="R-3")) == 1


@pytest.mark.parametrize(
    "elements",
    [[], ["H", "He", "Li", "Be", "B", "C", "N", "O", "F"], ["Sr", "O", "Sr"]],
    ids=["none", "nine", "repeated"],
)
def test_search_rejects_bad_element_lists(elements: list[str]) -> None:
    with pytest.raises(ValueError):
        cod_search(elements)


def test_record_reference_reads_as_a_citation() -> None:
    record = phases._record(cod_entry())

    assert record.reference() == (
        "Podlozhenov, S.; Graetsch, H. A. (2006) Acta Crystallographica Section B,"
        " doi:10.1107/S0108768106038869"
    )


def test_fetch_writes_the_cif_as_served(requests: list[str], tmp_path: Path) -> None:
    answer(f"{phases.COD_URL}/2100720.cif", CIF_TEXT)

    path = cod_fetch(2100720, tmp_path / "cifs")

    assert requests == [f"{phases.COD_URL}/2100720.cif"]
    assert path == tmp_path / "cifs" / "2100720.cif"
    assert path.read_bytes() == CIF_TEXT


@pytest.mark.parametrize("cod_id", ["210072", "21007200", "abc1234", "../1234567"])
def test_fetch_rejects_malformed_ids(cod_id: str, tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        cod_fetch(cod_id, tmp_path)


def test_fetch_rejects_a_body_that_is_not_a_cif(
    requests: list[str], tmp_path: Path
) -> None:
    answer(f"{phases.COD_URL}/2100720.cif", b"<html>maintenance</html>")

    with pytest.raises(ValueError):
        cod_fetch("2100720", tmp_path)
    assert not (tmp_path / "2100720.cif").exists()


def test_fetch_lets_a_missing_entry_raise(requests: list[str], tmp_path: Path) -> None:
    url = f"{phases.COD_URL}/9999999.cif"
    answer(url, HTTPError(url, 404, "Not Found", None, None))

    with pytest.raises(HTTPError):
        cod_fetch("9999999", tmp_path)


def read_index(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert tuple(reader.fieldnames) == CIF_INDEX_COLUMNS
        return list(reader)


def records_from(result: list[dict[str, object]]) -> list[CodRecord]:
    return [phases._record(entry) for entry in result]


def test_index_rows_of_cod_records(tmp_path: Path) -> None:
    path = write_cif_index(records_from(SEARCH_RESULT[:2]), tmp_path / "index.csv")

    rows = read_index(path)
    assert [row["file"] for row in rows] == ["2100720.cif", "2311739.cif"]
    assert rows[0] == {
        "file": "2100720.cif",
        "source": "COD",
        "identifier": "2100720",
        "formula": "Ba0.52 Nb2 O6 Sr0.48",
        "space_group": "P 4 b m",
        "a": "12.4844",
        "b": "12.4844",
        "c": "3.9572",
        "alpha": "90.0",
        "beta": "90.0",
        "gamma": "90.0",
        "reference": "Podlozhenov, S.; Graetsch, H. A. (2006) Acta Crystallographica"
        " Section B, doi:10.1107/S0108768106038869",
        "notes": "",
    }


def test_index_update_replaces_in_place_and_keeps_notes(tmp_path: Path) -> None:
    path = tmp_path / "index.csv"
    write_cif_index(records_from(SEARCH_RESULT[:2]), path)
    rows = read_index(path)
    rows[0]["notes"] = "closest to SBN50"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CIF_INDEX_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    updated = records_from([cod_entry(a="12.49"), SEARCH_RESULT[2]])
    write_cif_index(updated, path)

    rows = read_index(path)
    assert [row["identifier"] for row in rows] == ["2100720", "2311739", "1537507"]
    assert rows[0]["a"] == "12.49"
    assert rows[0]["notes"] == "closest to SBN50"


def test_index_takes_rows_from_other_sources(tmp_path: Path) -> None:
    path = tmp_path / "index.csv"
    write_cif_index(records_from(SEARCH_RESULT[:1]), path)

    write_cif_index(
        [
            {
                "file": "sbn50_own.cif",
                "source": "own refinement",
                "identifier": "sbn50-2026",
                "space_group": "P4bm",
                "notes": "Rietveld, x = 0",
            }
        ],
        path,
    )

    rows = read_index(path)
    assert [row["source"] for row in rows] == ["COD", "own refinement"]
    assert rows[1]["a"] == ""
    assert rows[1]["notes"] == "Rietveld, x = 0"


@pytest.mark.parametrize(
    "row",
    [{"source": "ICSD"}, {"source": "ICSD", "identifier": "1", "colour": "red"}],
    ids=["no identifier", "unknown column"],
)
def test_index_rejects_bad_rows(row: dict[str, str], tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        write_cif_index([row], tmp_path / "index.csv")


# Rock salt NaCl, a = 5.640 A, written by hand so the simulation needs no download.
NACL_CIF = """\
data_NaCl
_symmetry_space_group_name_H-M   'F m -3 m'
_symmetry_Int_Tables_number      225
_cell_length_a                   5.640
_cell_length_b                   5.640
_cell_length_c                   5.640
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1 Na 0.0 0.0 0.0 1.0
Cl1 Cl 0.5 0.5 0.5 1.0
"""
NACL_A = 5.640


@pytest.fixture
def nacl_cif(tmp_path: Path) -> Path:
    path = tmp_path / "nacl.cif"
    path.write_text(NACL_CIF, encoding="utf-8")
    return path


def test_simulated_nacl_lines_are_the_face_centred_ones(nacl_cif: Path) -> None:
    pytest.importorskip("pymatgen")

    pattern = simulate_pattern(nacl_cif, two_theta_range=(20, 80))

    families = [tuple(sorted(abs(i) for i in r.hkl)) for r in pattern]
    assert families[:5] == [(1, 1, 1), (0, 0, 2), (0, 2, 2), (1, 1, 3), (2, 2, 2)]
    for reflection in pattern:
        h, k, l = reflection.hkl
        assert len({h % 2, k % 2, l % 2}) == 1
        d = NACL_A / np.sqrt(h * h + k * k + l * l)
        expected = 2 * np.degrees(np.arcsin(KALPHA1_WAVELENGTH / (2 * d)))
        assert reflection.two_theta == pytest.approx(expected, abs=1e-3)


def test_simulated_nacl_is_normalised_to_its_200_line(nacl_cif: Path) -> None:
    pytest.importorskip("pymatgen")

    pattern = simulate_pattern(nacl_cif, two_theta_range=(20, 80))

    strongest = max(pattern, key=lambda r: r.intensity)
    assert strongest.intensity == pytest.approx(100)
    assert sorted(abs(i) for i in strongest.hkl) == [0, 0, 2]
    assert strongest.two_theta == pytest.approx(31.70, abs=0.01)
    # Na and Cl scatter out of phase for all odd hkl, so 111 is weak.
    assert pattern[0].intensity < 15


def test_simulate_pattern_follows_the_wavelength(nacl_cif: Path) -> None:
    pytest.importorskip("pymatgen")

    copper = simulate_pattern(nacl_cif, two_theta_range=(20, 80))
    cobalt = simulate_pattern(nacl_cif, wavelength=1.78897, two_theta_range=(20, 80))

    assert cobalt[1].two_theta > copper[1].two_theta
    assert all(20 <= r.two_theta <= 80 for r in copper)


def sim(*lines: tuple[float, float]) -> list[SimulatedReflection]:
    return [SimulatedReflection(t, i, (n, 0, 0)) for n, (t, i) in enumerate(lines)]


def test_match_explains_observed_peaks_within_tolerance() -> None:
    simulated = sim((24.83, 100), (26.70, 40), (30.00, 5))

    match = match_candidate([24.85, 26.73, 40.0], simulated)

    assert [peak.observed for peak in match.explained] == [24.85, 26.73]
    assert [peak.reflection.two_theta for peak in match.explained] == [24.83, 26.70]
    assert match.explained[0].offset == pytest.approx(0.02)
    assert match.missing == []
    assert match.score == 2


def test_match_takes_the_strongest_reflection_within_tolerance() -> None:
    simulated = sim((24.84, 3), (24.87, 60))

    (peak,) = match_candidate([24.85], simulated).explained

    assert peak.reflection.intensity == 60


def test_match_counts_strong_unseen_reflections_as_missing() -> None:
    simulated = sim((24.83, 100), (28.00, 50), (31.00, 8), (35.00, 10))

    match = match_candidate([24.85], simulated, min_intensity=10)

    assert [r.two_theta for r in match.missing] == [28.00]
    assert match.score == 0


def test_match_forgives_reflections_under_the_known_phase() -> None:
    simulated = sim((24.83, 100), (28.00, 50), (32.10, 70))

    match = match_candidate([24.85], simulated, exclude_two_theta=[27.97, 45.0])

    assert [r.two_theta for r in match.missing] == [32.10]
    assert match.score == 0


def test_match_with_nothing_near_scores_minus_the_missing() -> None:
    match = match_candidate([24.85, 26.73], sim((20.0, 100), (22.0, 30)))

    assert match.explained == []
    assert match.score == -2


def test_match_rejects_a_negative_tolerance() -> None:
    with pytest.raises(ValueError):
        match_candidate([24.85], sim((24.85, 100)), tolerance=-0.1)


def test_index_notes_apply_to_cod_records_and_replace_old_ones(tmp_path: Path) -> None:
    path = tmp_path / "index.csv"
    write_cif_index(records_from(SEARCH_RESULT[:1]), path, notes="old")

    write_cif_index(records_from(SEARCH_RESULT[:2]), path, notes="candidate")
    write_cif_index(
        [{"source": "own", "identifier": "1", "notes": "mine"}], path, notes="x"
    )

    assert [row["notes"] for row in read_index(path)] == [
        "candidate",
        "candidate",
        "mine",
    ]
