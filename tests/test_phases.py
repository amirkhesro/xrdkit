"""Tests for xrdkit.phases. The COD is never contacted: urlopen is replaced."""

import csv
import io
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlsplit

import pytest

from xrdkit import CodRecord, cod_fetch, cod_search, phases, write_cif_index
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
