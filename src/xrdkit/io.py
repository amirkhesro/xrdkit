"""Readers for X-ray diffraction data files."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

__all__ = ["XRDScan", "read_scan", "read_xrdml", "read_xy"]

# The scan files the package reads, and the suffixes of a two or three column
# text pattern: two theta, intensity and, in an .xye, the esd of the intensity.
XY_SUFFIXES = (".xy", ".xye")
SCAN_SUFFIXES = (".xrdml", *XY_SUFFIXES)

# A line of a text pattern begins the data when its first character could
# begin a number. Anything else before the data is a header or a comment.
_NUMBER_START = tuple("0123456789+-.")


@dataclass
class XRDScan:
    """A single powder diffraction scan.

    ``wavelength`` and ``time_per_step`` are None when the file carries
    neither, as a two column pattern does not: the wavelength is then
    supplied from the project file or the command line. ``esd`` holds the
    esds of the intensities when the file gives them, as an ``.xye`` does.
    """

    two_theta: np.ndarray
    intensity: np.ndarray
    wavelength: float | None
    start_angle: float
    end_angle: float
    step_size: float
    time_per_step: float | None
    sample_id: str
    source_path: str
    esd: np.ndarray | None = None


def _namespace(root: ET.Element) -> str:
    """Return the default namespace of ``root`` as a ``{uri}`` prefix."""
    if root.tag.startswith("{"):
        return root.tag[: root.tag.index("}") + 1]
    return ""


def _text(element: ET.Element | None) -> str:
    if element is None or element.text is None:
        return ""
    return element.text.strip()


def read_xrdml(path: str | Path) -> XRDScan:
    """Read a PANalytical ``.xrdml`` file into an :class:`XRDScan`.

    Parameters
    ----------
    path
        Path to the ``.xrdml`` file.

    Raises
    ------
    ValueError
        If the file contains no intensity data.
    """
    path = Path(path)
    root = ET.parse(path).getroot()
    ns = _namespace(root)

    scan = root.find(f".//{ns}scan")
    if scan is None:
        raise ValueError(f"No scan found in {path}")

    data_points = scan.find(f"{ns}dataPoints")
    if data_points is None:
        raise ValueError(f"No dataPoints element found in {path}")

    # Intensities live in <counts> or <intensities> depending on the writer.
    counts_element = data_points.find(f"{ns}counts")
    if counts_element is None:
        counts_element = data_points.find(f"{ns}intensities")
    raw_counts = _text(counts_element)
    if not raw_counts:
        raise ValueError(f"No intensity data found in {path}")

    intensity = np.array(raw_counts.split(), dtype=float)
    if intensity.size == 0:
        raise ValueError(f"No intensity data found in {path}")

    two_theta_positions = None
    for positions in data_points.findall(f"{ns}positions"):
        if positions.get("axis") == "2Theta":
            two_theta_positions = positions
            break
    if two_theta_positions is None:
        raise ValueError(f"No 2Theta axis found in {path}")

    start_angle = float(_text(two_theta_positions.find(f"{ns}startPosition")))
    end_angle = float(_text(two_theta_positions.find(f"{ns}endPosition")))

    n_points = intensity.size
    two_theta = np.linspace(start_angle, end_angle, n_points)
    step_size = (end_angle - start_angle) / (n_points - 1) if n_points > 1 else 0.0

    common_time = _text(data_points.find(f"{ns}commonCountingTime"))
    if common_time:
        time_per_step = float(common_time)
    else:
        per_point = _text(data_points.find(f"{ns}countingTimes"))
        time_per_step = float(per_point.split()[0]) if per_point else 0.0

    wavelength = float(_text(root.find(f".//{ns}usedWavelength/{ns}kAlpha1")))

    sample_id = _text(root.find(f".//{ns}sample/{ns}id"))

    return XRDScan(
        two_theta=two_theta,
        intensity=intensity,
        wavelength=wavelength,
        start_angle=start_angle,
        end_angle=end_angle,
        step_size=step_size,
        time_per_step=time_per_step,
        sample_id=sample_id,
        source_path=str(path),
    )


def _fields(line: str) -> list[str]:
    """The fields of a data line, separated by runs of spaces or tabs, by a
    comma, or by a comma and spaces."""
    return line.replace(",", " ").split()


def read_xy(path: str | Path, wavelength: float | None = None) -> XRDScan:
    """Read a two or three column text pattern into an :class:`XRDScan`.

    The first two columns are the two theta and the intensity; a third, as an
    ``.xye`` carries, is the esd of the intensity. Any number of header or
    comment lines may come first, a line counting as one until the data
    begins: the data begins at the first line whose first character could
    begin a number. Fields are separated by spaces, tabs or a comma, blank
    lines are ignored wherever they fall, and either line ending is read.

    ``wavelength`` is the K alpha 1 wavelength in angstroms to record, since
    the file carries none. None is allowed here; the command line is what
    insists on one where the wavelength is needed.

    Raises
    ------
    ValueError
        If the file holds no data, or a data line has fewer than two numbers,
        a field that is not a number, or a two theta that does not increase
        on the line before it. The message names the file and the line.
    """
    path = Path(path)
    two_theta: list[float] = []
    intensity: list[float] = []
    esds: list[float] = []
    # utf-8-sig so that a byte order mark, which some converters write, does
    # not hide the first line of data behind an unreadable character.
    text = path.read_text(encoding="utf-8-sig")
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if not two_theta and not line.startswith(_NUMBER_START):
            continue
        fields = _fields(line)
        if len(fields) < 2:
            raise ValueError(
                f"{path}, line {number}: needs two numbers, a two theta and an "
                "intensity"
            )
        try:
            values = [float(field) for field in fields]
        except ValueError:
            bad = next(field for field in fields if not _is_number(field))
            raise ValueError(
                f"{path}, line {number}: {bad!r} is not a number"
            ) from None
        if two_theta and values[0] <= two_theta[-1]:
            raise ValueError(
                f"{path}, line {number}: two theta {values[0]:g} does not increase "
                f"on {two_theta[-1]:g} before it"
            )
        two_theta.append(values[0])
        intensity.append(values[1])
        if len(values) > 2:
            esds.append(values[2])
    if not two_theta:
        raise ValueError(f"{path} holds no data")

    points = len(two_theta)
    start_angle, end_angle = two_theta[0], two_theta[-1]
    return XRDScan(
        two_theta=np.array(two_theta, dtype=float),
        intensity=np.array(intensity, dtype=float),
        wavelength=wavelength,
        start_angle=start_angle,
        end_angle=end_angle,
        step_size=(end_angle - start_angle) / (points - 1) if points > 1 else 0.0,
        time_per_step=None,
        sample_id=path.stem,
        source_path=str(path),
        # Only when every row gives one, so that the esds line up with the
        # intensities point for point.
        esd=np.array(esds, dtype=float) if len(esds) == points else None,
    )


def _is_number(field: str) -> bool:
    try:
        float(field)
    except ValueError:
        return False
    return True


def read_scan(path: str | Path, wavelength: float | None = None) -> XRDScan:
    """Read a scan of any format the package knows, by its suffix: ``.xrdml``
    with :func:`read_xrdml`, ``.xy`` and ``.xye`` with :func:`read_xy`.

    ``wavelength`` in angstroms replaces the scan's own where it is given,
    and supplies the one a text pattern carries no trace of.

    Raises
    ------
    ValueError
        If the suffix is not one of those, or the reader refuses the file.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".xrdml":
        scan = read_xrdml(path)
        return scan if wavelength is None else replace(scan, wavelength=wavelength)
    if suffix in XY_SUFFIXES:
        return read_xy(path, wavelength)
    raise ValueError(
        f"{path}: cannot read a {path.suffix or '(none)'!r} scan; xrdkit reads "
        + ", ".join(SCAN_SUFFIXES[:-1])
        + f" and {SCAN_SUFFIXES[-1]}"
    )
