"""Readers for X-ray diffraction data files."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["XRDScan", "read_xrdml"]


@dataclass
class XRDScan:
    """A single powder diffraction scan."""

    two_theta: np.ndarray
    intensity: np.ndarray
    wavelength: float
    start_angle: float
    end_angle: float
    step_size: float
    time_per_step: float
    sample_id: str
    source_path: str


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
