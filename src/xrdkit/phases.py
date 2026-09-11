"""Reference structures for phase identification from the Crystallography Open Database.

The COD is searched through its REST interface, which returns one JSON object
per entry, and each entry's CIF is downloaded by its seven digit COD id. The
CIFs gathered for a project are described by a CSV index, one row per file,
written by :func:`write_cif_index`.

Only the standard library is used, so fetching references adds no dependency.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

__all__ = [
    "CIF_INDEX_COLUMNS",
    "COD_URL",
    "CodRecord",
    "cod_fetch",
    "cod_search",
    "write_cif_index",
]

COD_URL = "https://www.crystallography.net/cod"

# Columns of the CIF index, in order.
CIF_INDEX_COLUMNS = (
    "file",
    "source",
    "identifier",
    "formula",
    "space_group",
    "a",
    "b",
    "c",
    "alpha",
    "beta",
    "gamma",
    "reference",
    "notes",
)

# The COD search form takes at most eight elements, el1 to el8.
MAX_ELEMENTS = 8

# Seconds to wait for the COD before giving up on a request.
TIMEOUT_S = 60.0

USER_AGENT = "xrdkit (https://github.com/amirkhesro/xrdkit)"


@dataclass(frozen=True)
class CodRecord:
    """One COD entry as returned by a search.

    Cell lengths are in angstroms, angles in degrees and the volume in cubic
    angstroms. Any value the COD leaves blank is ``None``.
    """

    cod_id: str
    formula: str
    space_group: str | None
    space_group_number: int | None
    a: float | None
    b: float | None
    c: float | None
    alpha: float | None
    beta: float | None
    gamma: float | None
    volume: float | None
    authors: str
    journal: str
    year: int | None
    doi: str | None = None

    @property
    def filename(self) -> str:
        """Name :func:`cod_fetch` gives this entry's CIF."""
        return f"{self.cod_id}.cif"

    def reference(self) -> str:
        """Return a one line citation, such as 'Smith, J. (2006) Acta Cryst. B'."""
        text = self.authors
        if self.year is not None:
            text += f" ({self.year})"
        if self.journal:
            text += f" {self.journal}"
        if self.doi:
            text += f", doi:{self.doi}"
        return text.strip()


def cod_search(
    elements: Iterable[str],
    exact: bool = True,
    space_group: str | int | None = None,
    extra: Mapping[str, str | int | float] | None = None,
) -> list[CodRecord]:
    """Search the COD for entries containing ``elements``.

    Parameters
    ----------
    elements
        Element symbols, one to eight of them, such as ``["Sr", "Ba", "Nb", "O"]``.
    exact
        If true, only entries made of exactly these elements are returned;
        otherwise entries containing these elements and any others.
    space_group
        Keep only entries in this space group. An integer is the International
        Tables number and is passed to the COD's ``space_group_number`` search.
        A string is a Hermann-Mauguin symbol such as ``"P4bm"``, compared with
        each entry's symbol after spaces and any trailing setting, such as
        ``" (a,b,2*c)"`` or ``":H"``, are removed. The COD ignores a symbol in
        the query itself, so the symbol is matched here and must be written as
        the COD writes it (``"P121/c1"`` rather than ``"P21/c"``); entries
        whose symbol the COD gives in a nonstandard centring such as ``X4bm``
        are not matched by either form.
    extra
        Further COD search parameters, such as ``{"year": 2006}``. They are
        added to the query as given and override those set from the arguments.

    Raises
    ------
    ValueError
        If no element or more than eight are given, or one is repeated.
    urllib.error.URLError
        If the COD cannot be reached or answers with an error.
    """
    symbols = [element.strip() for element in elements]
    if not symbols:
        raise ValueError("no elements given")
    if len(symbols) > MAX_ELEMENTS:
        raise ValueError(
            f"the COD searches at most {MAX_ELEMENTS} elements, got {len(symbols)}"
        )
    if len(set(symbols)) != len(symbols):
        raise ValueError(f"repeated element in {symbols}")

    params: dict[str, str | int | float] = {"format": "json"}
    for number, element in enumerate(symbols, start=1):
        params[f"el{number}"] = element
    if exact:
        params["strictmin"] = len(symbols)
        params["strictmax"] = len(symbols)
    if isinstance(space_group, int):
        params["space_group_number"] = space_group
    if extra:
        params.update(extra)

    body = _get(f"{COD_URL}/result?{urlencode(params)}")
    records = [_record(entry) for entry in json.loads(body.decode("utf-8"))]
    if isinstance(space_group, str):
        wanted = _normalise_symbol(space_group)
        records = [
            record
            for record in records
            if record.space_group is not None
            and _normalise_symbol(record.space_group) == wanted
        ]
    return records


def cod_fetch(cod_id: str | int, folder: str | Path) -> Path:
    """Download the CIF of COD entry ``cod_id`` into ``folder`` as ``<id>.cif``.

    The folder is created if needed and an existing file of the same name is
    replaced. The file is written exactly as the COD serves it.

    Raises
    ------
    ValueError
        If ``cod_id`` is not a seven digit COD id, or what comes back is not a CIF.
    urllib.error.URLError
        If the COD cannot be reached or has no such entry.
    """
    identifier = str(cod_id).strip()
    if not re.fullmatch(r"\d{7}", identifier):
        raise ValueError(f"not a COD id: {cod_id!r}")
    body = _get(f"{COD_URL}/{identifier}.cif")
    if b"data_" not in body:
        raise ValueError(f"COD {identifier} did not return a CIF")
    path = Path(folder) / f"{identifier}.cif"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(body)
    return path


def write_cif_index(
    records: Iterable[CodRecord | Mapping[str, object]], path: str | Path
) -> Path:
    """Write or update the CSV index of reference CIFs at ``path``.

    Each record is either a :class:`CodRecord`, entered with source ``COD``,
    its COD id as the identifier and ``<id>.cif`` as the file, or a mapping
    keyed by :data:`CIF_INDEX_COLUMNS` for a CIF from elsewhere. A row is
    identified by its source and identifier: a record matching a row already in
    the index replaces it in place, and one that matches none is added at the
    end. Rows not named by any record are kept, and so are the notes of a
    replaced row when the new record brings none, so notes written by hand
    survive a fresh download.

    Raises
    ------
    ValueError
        If a mapping has a column not in the index, or no source or identifier.
    """
    path = Path(path)
    rows: dict[tuple[str, str], dict[str, str]] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                entry = {column: row.get(column) or "" for column in CIF_INDEX_COLUMNS}
                rows[entry["source"], entry["identifier"]] = entry

    for record in records:
        entry = _index_row(record)
        key = (entry["source"], entry["identifier"])
        previous = rows.get(key)
        if previous is not None and not entry["notes"]:
            entry["notes"] = previous["notes"]
        rows[key] = entry

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CIF_INDEX_COLUMNS)
        writer.writeheader()
        writer.writerows(rows.values())
    return path


def _get(url: str) -> bytes:
    """Return the body of ``url``, raising on an HTTP error."""
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=TIMEOUT_S) as response:
        return response.read()


def _record(entry: Mapping[str, object]) -> CodRecord:
    """Build a :class:`CodRecord` from one object of the COD's JSON result."""
    return CodRecord(
        cod_id=str(entry["file"]),
        formula=_formula(entry.get("formula")),
        space_group=_text(entry.get("sg")),
        space_group_number=_integer(entry.get("sgNumber")),
        a=_number(entry.get("a")),
        b=_number(entry.get("b")),
        c=_number(entry.get("c")),
        alpha=_number(entry.get("alpha")),
        beta=_number(entry.get("beta")),
        gamma=_number(entry.get("gamma")),
        # "volume" in the COD result is the journal volume; the cell's is "vol".
        volume=_number(entry.get("vol")),
        authors=_text(entry.get("authors")) or "",
        journal=_text(entry.get("journal")) or "",
        year=_integer(entry.get("year")),
        doi=_text(entry.get("doi")),
    )


def _index_row(record: CodRecord | Mapping[str, object]) -> dict[str, str]:
    """Return the index row of one record, every column as text."""
    if isinstance(record, CodRecord):
        values: Mapping[str, object] = {
            "file": record.filename,
            "source": "COD",
            "identifier": record.cod_id,
            "formula": record.formula,
            "space_group": record.space_group,
            "a": record.a,
            "b": record.b,
            "c": record.c,
            "alpha": record.alpha,
            "beta": record.beta,
            "gamma": record.gamma,
            "reference": record.reference(),
            "notes": "",
        }
    else:
        unknown = sorted(set(record) - set(CIF_INDEX_COLUMNS))
        if unknown:
            raise ValueError(f"not index columns: {', '.join(unknown)}")
        values = record
    row = {
        column: "" if values.get(column) is None else str(values[column])
        for column in CIF_INDEX_COLUMNS
    }
    if not row["source"] or not row["identifier"]:
        raise ValueError(f"index row needs a source and an identifier: {row}")
    return row


def _normalise_symbol(symbol: str) -> str:
    """Strip spaces and a trailing setting from a Hermann-Mauguin symbol."""
    symbol = re.sub(r"\(.*\)\s*$", "", symbol)
    symbol = symbol.split(":", 1)[0]
    return re.sub(r"\s+", "", symbol)


def _formula(value: object) -> str:
    """Drop the dashes the COD puts round a formula: '- Nb O3 Sr -' to 'Nb O3 Sr'."""
    return (_text(value) or "").strip("- ").strip()


def _text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _number(value: object) -> float | None:
    # The COD gives esds in separate fields, but drop a '(3)' in case one slips in.
    text = _text(value)
    return None if text is None else float(re.sub(r"\(\d+\)$", "", text))


def _integer(value: object) -> int | None:
    text = _text(value)
    return None if text is None else int(text)
