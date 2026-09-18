"""The instrument parameter file of a diffractometer, from a standard scan.

A standard of known cell, LaB6 (SRM 660) or silicon (SRM 640) measured on the
instrument with the optics a sample will be measured with, carries no sample
broadening worth the name, so the widths of its reflections are the
instrument's own. :func:`fit_instrument_widths` fits each reflection of such a
scan as a K alpha doublet and the widths to the Caglioti relation, which needs
no GSAS-II. :func:`refine_instrument` then takes that fit as the starting
point of a GSAS-II refinement of the profile terms, and gives back the file
GSAS-II exports, which is what every later refinement reads.

The two are separate so that the width fit, which is quick and answerable on
its own, can be looked at before the refinement is run.
"""

from __future__ import annotations

import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from xrdkit.broadening import Caglioti, fit_caglioti, fit_profile
from xrdkit.gsas2 import (
    Gsas2Install,
    build_refine_job,
    run_job,
    standard_stages,
    write_instprm,
)
from xrdkit.io import XRDScan
from xrdkit.peaks import exclude_kalpha2, find_peaks

__all__ = [
    "REFINED_KEYS",
    "WIDTH_WINDOW",
    "InstrumentRefinement",
    "WidthFit",
    "fit_instrument_widths",
    "kalpha2_wavelength",
    "refine_instrument",
]

# The two theta range the widths are fitted over. Below the lower limit the
# reflections of a standard are few and asymmetric; above the upper one the
# K alpha doublet is wide enough that the fit is about the splitting rather
# than the instrument.
WIDTH_WINDOW = (10.0, 98.0)

# The instrument parameters the refinement frees, in the order it frees them:
# the zero, the Gaussian width terms, the Lorentzian terms and the axial
# divergence. The cell is held at the standard's certified value, which is the
# point of using a standard.
REFINED_KEYS = ("Zero", "U", "V", "W", "X", "Y", "SH/L")

# A standard is ground and sieved to leave no size or strain broadening, so
# the phase is given a large crystallite size and no microstrain, held.
STANDARD_BROADENING = {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}

# Refinement cycles per stage, as the standard sequence uses.
CYCLES = 10


@dataclass(frozen=True)
class WidthFit:
    """The Caglioti width fit of a standard scan.

    ``u``, ``v`` and ``w`` are in degrees squared of two theta and ``rms`` in
    degrees, the root mean square of the FWHM residuals. ``two_theta`` and
    ``fwhm`` are the reflections the fit used, those whose profile fit
    converged; ``n_peaks`` counts them. ``wavelength`` is the scan's, in
    angstroms.
    """

    n_peaks: int
    u: float
    v: float
    w: float
    esd_u: float
    esd_v: float
    esd_w: float
    rms: float
    wavelength: float
    window: tuple[float, float]
    two_theta: tuple[float, ...]
    fwhm: tuple[float, ...]
    covariance: tuple[tuple[float, float, float], ...]

    @property
    def caglioti(self) -> Caglioti:
        """The fit as a :class:`~xrdkit.broadening.Caglioti`, for
        :func:`~xrdkit.gsas2.write_instprm`."""
        return Caglioti(
            u=self.u,
            v=self.v,
            w=self.w,
            esd_u=self.esd_u,
            esd_v=self.esd_v,
            esd_w=self.esd_w,
            n_peaks=self.n_peaks,
            rms=self.rms,
            covariance=np.array(self.covariance, dtype=float),
        )


@dataclass(frozen=True)
class InstrumentRefinement:
    """What the GSAS-II refinement of a standard came to.

    ``parameters`` maps each of :data:`REFINED_KEYS` to its refined ``value``
    and ``esd``. ``instprm`` is the file GSAS-II exported, copied into place.
    ``result`` is the whole job result, for a caller that wants more.
    """

    parameters: dict[str, dict[str, float]]
    rwp: float
    gof: float
    stages: tuple[str, ...]
    instprm: Path
    start_instprm: Path
    result: dict


def fit_instrument_widths(
    scan: XRDScan, window: tuple[float, float] = WIDTH_WINDOW
) -> WidthFit:
    """Fit the widths of a standard's reflections to the Caglioti relation.

    Every peak found in ``window``, its K alpha 2 satellites excluded, is
    fitted as a split pseudo-Voigt doublet, and the widths of the fits that
    converged go into a weighted fit of FWHM squared against tan theta, each
    weighted by one over the square of its own esd.

    Raises
    ------
    ValueError
        If the scan carries no wavelength, or fewer than three reflections
        converged, which is the least a three parameter relation can be fitted
        to.
    """
    if scan.wavelength is None:
        raise ValueError("the scan carries no wavelength")
    peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=window))
    angles: list[float] = []
    widths: list[float] = []
    weights: list[float] = []
    for peak in peaks:
        profile = fit_profile(scan.two_theta, scan.intensity, peak.two_theta, peak.fwhm)
        if profile.converged:
            angles.append(profile.two_theta)
            widths.append(profile.fwhm)
            weights.append(1.0 / profile.esd_fwhm**2)
    if len(angles) < 3:
        raise ValueError(
            f"{len(angles)} reflections of the standard converged between "
            f"{window[0]:g} and {window[1]:g} degrees; the width fit needs three"
        )
    caglioti = fit_caglioti(np.array(angles), np.array(widths), np.array(weights))
    return WidthFit(
        n_peaks=caglioti.n_peaks,
        u=caglioti.u,
        v=caglioti.v,
        w=caglioti.w,
        esd_u=caglioti.esd_u,
        esd_v=caglioti.esd_v,
        esd_w=caglioti.esd_w,
        rms=caglioti.rms,
        wavelength=float(scan.wavelength),
        window=(float(window[0]), float(window[1])),
        two_theta=tuple(angles),
        fwhm=tuple(widths),
        covariance=tuple(tuple(float(v) for v in r) for r in caglioti.covariance),
    )


def refine_instrument(
    fit: WidthFit,
    scan_file: str | Path,
    cif: str | Path,
    phase: str,
    cell: tuple[float, float, float, float, float, float],
    folder: str | Path,
    stem: str,
    work: str | Path,
    install: Gsas2Install | None = None,
    radius_mm: float | None = None,
) -> InstrumentRefinement:
    """Refine the instrument parameters of ``scan_file`` in GSAS-II from ``fit``.

    The width fit is written as ``<stem>_start.instprm`` in ``folder`` and is
    what the refinement starts from; the standard's phase is held at ``cell``,
    with no size or microstrain broadening, and the stages of
    :func:`~xrdkit.gsas2.standard_stages` are run without the cell stage. The
    file GSAS-II exports is copied to ``<stem>.instprm`` in ``folder``, and
    that is the file every later refinement reads.

    Raises
    ------
    xrdkit.gsas2.Gsas2Error
        If GSAS-II failed.
    KeyError
        If the run left no final instrument parameters.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    start = write_instprm(
        folder / f"{stem}_start.instprm", fit.caglioti, radius_mm=radius_mm
    )
    stages = [stage for stage in standard_stages() if stage["name"] != "cell"]
    work = Path(work)
    job = build_refine_job(
        (work / f"{stem}.gpx").resolve(),
        stages,
        data_file=Path(scan_file).resolve(),
        instprm=start.resolve(),
        phases=[{"cif": Path(cif).resolve(), "name": phase, "cell": list(cell)}],
        limits=fit.window,
        cycles=CYCLES,
        broadening={phase: dict(STANDARD_BROADENING)},
        export_prefix=(work / stem).resolve(),
    )
    result = run_job(job, work, install)
    last = result["stages"][-1]
    final = result["final"]["instrument"]
    instprm = folder / f"{stem}.instprm"
    shutil.copyfile(result["exports"]["instprm"], instprm)
    return InstrumentRefinement(
        parameters={
            key: {
                "value": float(final[key]["value"]),
                "esd": None
                if final[key].get("esd") is None
                else float(final[key]["esd"]),
            }
            for key in REFINED_KEYS
            if key in final
        },
        rwp=float(last["rwp"]),
        gof=float(last["gof"]),
        stages=tuple(stage["name"] for stage in stages),
        instprm=instprm,
        start_instprm=start,
        result=result,
    )


def kalpha2_wavelength(path: str | Path) -> float | None:
    """The K alpha 2 wavelength an ``.xrdml`` records, or None.

    A two or three column text pattern carries no wavelength at all, so the
    answer for one is None: such a scan is taken as K alpha 1 only unless the
    project file says otherwise.
    """
    path = Path(path)
    if path.suffix.lower() != ".xrdml":
        return None
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    namespace = root.tag[: root.tag.index("}") + 1] if root.tag.startswith("{") else ""
    element = root.find(f".//{namespace}usedWavelength/{namespace}kAlpha2")
    if element is None or not (element.text or "").strip():
        return None
    try:
        return float(element.text.strip())
    except ValueError:
        return None
