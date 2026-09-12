# xrdkit user guide

1. [Introduction](#1-introduction)
2. [Data quality, in numbers](#2-data-quality-in-numbers)
3. [Files and formats](#3-files-and-formats)
4. [Workflow 1: plotting a pattern with hkl indices](#4-workflow-1-plotting-a-pattern-with-hkl-indices)
5. [Workflow 2: lattice parameters and theoretical density](#5-workflow-2-lattice-parameters-and-theoretical-density)
6. [Workflow 3: Rietveld refinement](#6-workflow-3-rietveld-refinement)
7. [Known limitations of version 0.1.0](#7-known-limitations-of-version-010)

If Python is not yet installed on your machine, start with
[GETTING_STARTED.md](GETTING_STARTED.md), which installs Python and xrdkit
on Windows or macOS, sets up the folder layout used below, and shows how to
run the code blocks in this guide.

## 1. Introduction

xrdkit is a library for the analysis of laboratory powder X-ray diffraction
patterns from electroceramic samples. It reads a scan from the diffractometer,
finds the peaks in it, indexes them against a cell, refines lattice parameters,
calculates a theoretical density, draws publication quality figures, and drives
a GSAS-II Rietveld refinement from a description of the refinement written as
data. It is a library rather than a program: every routine is imported and
called from a script of your own, so that the same analysis can be run over a
whole series of samples without anyone clicking through a dialogue.

The kit supports three workflows. The first is plotting a pattern with hkl
indices on it, which is what a phase check and most figures in a paper need.
The second is lattice parameters and a theoretical density from a Le Bail fit,
which is what a composition series or a solid solution study needs. The third
is a full Rietveld refinement of coordinates and occupancies, which is
expensive in beam time and in effort and is normally worth it only for selected
samples. The usual order is to plot every sample first, identify the phases
present, index the pattern of each sample against the expected cell, refine the
cell and the density from that, and only then take the few samples that matter
through a Rietveld refinement.

The kit never modifies raw data. The reader opens a scan file for reading and
returns the numbers it holds; nothing is written back into the file, and no
routine in the kit writes anything into the raw data directory. Output goes
where you send it: peak lists and indexing tables to the CSV path you give,
figures to the path you give, refinement projects and exports to the working
directory you name in the job. A raw scan can therefore be reprocessed any
number of times and remains the record of what the instrument measured.

## 2. Data quality, in numbers

What the data must be depends on what you intend to get out of them. A scan
that is ample for a phase check is not good enough for a lattice parameter, and
a scan that gives a good lattice parameter is still not good enough for refined
occupancies. The criteria below are the ones the three workflows need.

| | Plotting and phase check | Lattice parameters and density, by Le Bail | Rietveld refinement |
| --- | --- | --- | --- |
| Purpose | Show what the sample is, identify the phases, label the reflections | Refine a and c, and the cell volume, well enough for a theoretical density | Refine the structure: coordinates, occupancies, displacement parameters |
| Angular range | 10 to 80 degrees two theta, extended to 90 for a figure meant for publication | 10 to 120 degrees. The reflections above 90 degrees are what pin the cell down, because a given error in d shifts them furthest in two theta | 5 to 130 degrees, or wider where the instrument allows |
| Step size | 0.01 to 0.03 degrees | 0.013 to 0.026 degrees | 0.01 to 0.02 degrees |
| Counting statistics | Strongest peak at least 20 times the background and above about 2000 counts. For a secondary phase at 1 to 2 weight per cent to be visible on a square root or logarithmic scale, the background itself needs at least 100 counts per step, which means three to five times the count time of a quick scan | Strongest peak above 10000 counts, and the high angle peaks at least 10 times the background, since those are the ones the cell is refined on | Strongest peak above 20000 counts, and background above 200 counts per step at high angle. On a laboratory instrument this normally means several hours per scan |
| Sample preparation | Flat sample, flush with the surface of the holder | Crushed pellet powder preferred. A pellet surface is acceptable if a displacement term is refined | Powder crushed and sieved below about 45 micrometres, then back loaded, side loaded or mounted in a capillary, to limit preferred orientation |
| Standards needed | None | Either an internal standard, silicon SRM 640 at 10 to 20 weight per cent mixed into the powder, or a refined displacement term. Also an instrument parameter file from a LaB6 (SRM 660) scan on the same instrument and the same optics | An instrument parameter file from a LaB6 scan on the same configuration, a CIF of the expected structure, and a known composition |
| Radiation | K alpha 2 present is acceptable | K alpha 2 present is acceptable, since the fit follows the K alpha 1 positions | Monochromatic Cu K alpha 1 strongly preferred |
| What goes wrong otherwise | Weak phases hide in the noise and are missed, so the sample is reported as single phase when it is not. Peaks are too noisy for the peak finder to separate, and labels go on the wrong reflections | The cell is refined on low angle reflections alone, where a displacement error and a cell error look alike, so a and c come out precise and wrong, and the density with them | The refinement runs, but the parameters are not determined by the data. Occupancies and coordinates drift to whatever fits the noise, and the esds do not say so unless they are examined |

Two further numbers govern the first two workflows. For hkl labelling, peak
positions good to about 0.02 degrees are enough, and a zero error of up to 0.05
degrees does not change which reflection a peak is assigned to. That same zero
error does change the lattice parameters, which is why the second workflow
either uses an internal standard or refines a displacement or zero term: an
uncorrected specimen displacement of 0.1 mm shifts a by about 0.001 angstrom.
For a useful density, a and c are needed to a relative precision of 0.01 per
cent, which is about 0.001 angstrom on a 12 angstrom axis, because the density
goes as the reciprocal of the cell volume and the relative error in the density
is therefore about three times the relative error in a lattice parameter.

For a Rietveld refinement on good data, expect an Rwp of a few per cent and a
goodness of fit between 1 and 2. Coordinates and occupancies should be refined
only where the data justify it. The kit reports this rather than leaving it to
be guessed at: a stage that reached its pass cap while still moving is recorded
as unsettled, and a refined parameter whose esd is larger than its shift, or in
the case of an occupancy larger than half its allowed range, is recorded as
undetermined. An undetermined value is not a result, whatever the fit looks
like.

### Checking a scan against these criteria

Read the scan and print what it is: the range covered, the step, the count
time, the strongest intensity and the median intensity. The median stands in
for the background, since most of the points in a powder pattern are
background.

```python
import numpy as np

from xrdkit import read_xrdml

scan = read_xrdml("data/raw/10s.xrdml")
print(f"range   {scan.start_angle:.2f} to {scan.end_angle:.2f} degrees")
print(f"step    {scan.step_size:.4f} degrees")
print(f"points  {scan.intensity.size}")
print(f"time    {scan.time_per_step:.1f} s per step")
print(f"maximum {scan.intensity.max():.0f} counts")
print(f"median  {np.median(scan.intensity):.0f} counts")
print(f"peak over median {scan.intensity.max() / np.median(scan.intensity):.0f}")
```

```
range   10.01 to 99.98 degrees
step    0.0217 degrees
points  4141
time    34.2 s per step
maximum 13964 counts
median  998 counts
peak over median 14
```

Read against the table, that scan is comfortable for plotting and for phase
identification, and it reaches nearly to 100 degrees, but its strongest peak is
below the 10000 counts a Le Bail cell refinement wants and its peak to
background ratio of 14 is short of 20, so weak secondary phases may not be
visible in it.

## 3. Files and formats

| What you supply | Format | One per | Where it comes from | Required |
| --- | --- | --- | --- | --- |
| Raw scan | `.xrdml` | Scan | The diffractometer. PANalytical Aeris, X'Pert3 and Empyrean all write this form, and the reader takes them unchanged | Required |
| Standard scan | `.xrdml` | Instrument and optical configuration | A LaB6 (SRM 660) or silicon (SRM 640) powder measured on the same instrument with the same optics. It gives the instrumental peak widths and, through them, the instrument parameter file | Required for the Le Bail and Rietveld workflows |
| CIF | `.cif` | Phase | The Crystallography Open Database, the ICSD, or the supporting information of a paper. Record which, and the entry number, alongside the file | Required for the Rietveld workflow, and for phase matching |
| Instrument parameter file | `.instprm` | Instrument and optical configuration | Written by the kit from the LaB6 scan, with `xrdkit.gsas2.write_instprm`, from a Caglioti fit to the standard's peak widths | Required for the Rietveld workflow |
| Sample metadata | Your own table, or a TOML settings file | Sample | You. The sample name, the nominal composition and the full formula, the formula units per cell Z, whether the sample is calcined or sintered, the instrument, and the scan settings | Required |
| Archimedes density | A number, in g/cm3 | Sample | Your own measurement | Optional, and needed only to quote a relative density against the theoretical one |
| GSAS-II installation | A separate Python installation | Machine | Installed separately from xrdkit, and pointed at with the environment variables `XRDKIT_GSAS2_PYTHON` and `XRDKIT_GSAS2_HOME`, or left at `~/gsas2main` | Required for the Rietveld workflow only |

### What the reader accepts

The reader in `xrdkit.io` is `read_xrdml`, and `.xrdml` is the only format it
accepts. It parses the XML, takes the intensities from the `counts` element or,
where the writer used that name instead, from `intensities`, takes the start
and the end of the two theta axis and derives the step from the number of
points, takes the counting time per step, takes the K alpha 1 wavelength from
the used wavelength block, and takes the sample identifier. It returns an
`XRDScan` whose fields are `two_theta`, `intensity`, `wavelength`,
`start_angle`, `end_angle`, `step_size`, `time_per_step`, `sample_id` and
`source_path`.

Plain text column formats are not supported: there is no reader for two column
or three column `.xy` or `.xye` files, none for Bruker `.raw` or `.brml`, and
none for `.gsas` or `.fxye`. A scan in one of those forms has to be converted
to `.xrdml`, or an `XRDScan` has to be built directly from the columns, since
it is an ordinary dataclass and every routine downstream of the reader takes an
`XRDScan` rather than a file path. Section 7 shows how.

### Folder layout

The recommended layout separates what the instrument produced from everything
derived from it.

```
data/raw/          scans as the instrument wrote them, never edited
data/standards/    LaB6 and silicon scans, and the instrument parameter files
cifs/              one CIF per phase, with a note of its source
config/            the TOML settings the Rietveld workflow reads
results/           peak lists, indexing tables, refinement projects and exports
figures/           the png and pdf files the plotting routines write
```

Code lives in its own repository and data in another, so that the analysis can
be versioned without the scans going with it.

## 4. Workflow 1: plotting a pattern with hkl indices

This workflow takes a scan and produces a figure with the reflections labelled.
There are two routes through it, a single pattern and a stack of several
patterns, and they are alternatives rather than steps: the peak finding and
the indexing in the middle are the same either way, and only the plotting
differs.

### Step 1. Read the scan

Reading gives an `XRDScan`, which carries the two theta and intensity arrays
and the wavelength that every later step needs.

```python
from xrdkit import apply_style, read_xrdml

apply_style()
scan = read_xrdml("data/raw/10s.xrdml")
```

`apply_style` sets the matplotlib rcParams to a clean single column journal
style. Call it once, before any plotting.

### Step 2a. Plot a single pattern

`plot_pattern` draws one scan and returns the figure and the axes, so that a
title, the limits, or anything else can be set before the figure is saved. A
linear scale shows the strong reflections in proportion; a square root scale
brings up the weak ones, which is what a phase check wants.

```python
from xrdkit import plot_pattern, save_figure

fig, ax = plot_pattern(scan, scale="linear", colour="black")
ax.set_title("Sample 10, as measured")
save_figure(fig, "figures/pattern_10")

fig, ax = plot_pattern(scan, scale="sqrt")
save_figure(fig, "figures/pattern_10_sqrt", formats=("png", "pdf"), dpi=600)
```

`save_figure` treats the path as a stem and writes one file per format, png and
pdf by default, at 300 dpi by default. It returns the list of paths written.

### Step 2b. Plot several scans stacked

`plot_stacked` draws a list of scans one above another with a constant offset,
each labelled at its upper right. The traces are normalised by default, so that
samples measured for different times can be compared. It returns four things:
the figure, the axes, the vertical base of each slot, and the line drawn for
each scan. The bases and the lines are what the annotation step needs in order
to put labels against a chosen trace.

```python
from xrdkit import plot_stacked, save_figure

scans = [read_xrdml(path) for path in ("data/raw/10s.xrdml", "data/raw/12s.xrdml")]
fig, ax, bases, lines = plot_stacked(
    scans, labels=["x = 0.10", "x = 0.12"], scale="sqrt", offset=None
)
save_figure(fig, "figures/stack")
```

With `offset=None` the spacing is 1.2 times the tallest scaled trace, which
keeps the tallest peak of one pattern clear of the pattern above it. Pass a
number to set the spacing yourself. For the two scans above, the bases came
back as `[0.0, 1.2]`.

### Step 3. Find the peaks

`find_peaks` locates the reflections and returns them in two theta order, each
with its position, intensity, prominence, full width at half maximum, d spacing
and relative intensity as a percentage of the strongest peak found. Restrict
the search to the range you intend to plot, so that the relative intensities
refer to that range.

Peaks that look like K alpha 2 satellites are flagged rather than dropped, so
that they can be removed when it suits. Remove them before indexing: a
satellite sits at the position its parent's d spacing gives at the longer
wavelength, so indexing it against the cell is meaningless.

```python
from xrdkit import exclude_kalpha2, find_peaks, peaks_to_csv

peaks = find_peaks(scan, min_prominence=0.02, two_theta_range=(10.0, 80.0))
clean = exclude_kalpha2(peaks)
peaks_to_csv(clean, "results/peaks_10.csv")
print(f"{len(peaks)} peaks found, {len(clean)} after removing K alpha 2 satellites")
```

```
41 peaks found, 35 after removing K alpha 2 satellites
```

`exclude_kalpha2` returns a new list and recomputes the relative intensities
against the strongest of the survivors, so the strongest is again 100. Below
about 50 degrees the doublet is not resolved and there is nothing to flag.

### Step 4. Index against a known cell

`index_and_refine` assigns a reflection of the cell to each peak and refines
the cell as it goes. It works in cycles: the first indexes only the low angle
peaks with a loose tolerance, where a cell that is still some way off can be
trusted to put reflections near the right peaks, and each later cycle indexes
the whole list against the cell the previous cycle gave. It also searches for
the zero offset first, unless you pass one.

The indexing in `xrdkit.indexing` is tetragonal, and `P4bm` is the only space
group whose reflection conditions it knows; pass `space_group=None` to apply
none. `TTB_CELL` is the tetragonal tungsten bronze starting cell, a = 12.45 and
c = 3.94 angstrom.

```python
from xrdkit import TTB_CELL, index_and_refine, indexed_to_csv, indexing_summary

indexed, fit = index_and_refine(
    clean, start_cell=TTB_CELL, wavelength=scan.wavelength, space_group="P4bm"
)
indexed_to_csv(indexed, "results/indexed_10.csv")
print(f"a = {fit.cell.a:.4f}, c = {fit.cell.c:.4f} angstrom")
print(f"{fit.n_peaks} peaks used, rms {fit.rms_two_theta:.4f} degrees")
print(f"zero offset {fit.zero_offset:.3f} degrees")
print(indexing_summary(indexed))
```

```
a = 12.4799, c = 3.9323 angstrom
28 peaks used, rms 0.0076 degrees
zero offset 0.170 degrees
{'n_peaks': 35, 'n_indexed': 35, 'n_unindexed': 0, 'n_ambiguous': 7, 'rms_difference': 0.007559815493900485}
```

A cell refined this way is good enough to label reflections with. It is not a
lattice parameter for publication, which is what Workflow 2 is for: the fit
here is a linear least squares on the indexed peaks alone, and the zero offset
it reports is the one that indexed the most peaks rather than one refined
alongside the cell.

### Step 5a. Label a single pattern

`annotate_hkl` writes an hkl label above each indexed peak. Given the line of
the trace, each label rides on top of its own peak, a couple of points above
the trace, which keeps the label and the reflection it names together however
the pattern rises and falls. Set the figure size and the x limits before
calling it, because it measures the rendered labels against the geometry as it
stands.

```python
from xrdkit import annotate_hkl, plot_pattern, save_figure

fig, ax = plot_pattern(scan, scale="sqrt")
ax.set_xlim(10.0, 80.0)
labels = annotate_hkl(ax, indexed, min_relative_intensity=5.0, line=ax.lines[0])
save_figure(fig, "figures/pattern_hkl")
print(f"{len(labels)} labels written")
```

```
20 labels written
```

The alternative is a single row of labels at a fixed height, which suits a
normalised trace with room left above it. Labels are placed strongest first,
and a weaker peak that would collide with one already placed keeps no label, so
in a crowded stretch it is the weak reflections that lose theirs.

```python
fig, ax = plot_pattern(scan, scale="sqrt", normalise=True)
ax.set_xlim(10.0, 80.0)
ax.set_ylim(0.0, 1.35)
labels = annotate_hkl(ax, indexed, y=1.05, min_relative_intensity=10.0, rotation=90)
save_figure(fig, "figures/pattern_hkl_row")
```

### Step 5b. Label the top trace of a stack

In a stack the labels normally go on the top trace only, which is enough to
tell the reader what every trace below shows. Pass the line of that trace as
`line`, and index the peaks of that same scan, not of another one: a label
stands over an observed position, so it belongs to the pattern it was found in.

```python
from xrdkit import annotate_hkl, plot_stacked, save_figure

top_scan = scans[-1]
top_peaks = exclude_kalpha2(find_peaks(top_scan, two_theta_range=(10.0, 80.0)))
top_indexed, top_fit = index_and_refine(
    top_peaks, start_cell=TTB_CELL, wavelength=top_scan.wavelength
)

fig, ax, bases, lines = plot_stacked(
    scans, labels=["x = 0.10", "x = 0.12"], scale="sqrt"
)
ax.set_xlim(10.0, 80.0)
top = annotate_hkl(ax, top_indexed, line=lines[-1], min_relative_intensity=10.0)
save_figure(fig, "figures/stack_hkl")
print(f"a = {top_fit.cell.a:.4f}, c = {top_fit.cell.c:.4f} angstrom, {len(top)} labels")
```

```
a = 12.4707, c = 3.9268 angstrom, 19 labels
```

To put the labels in a row above a chosen trace instead of on it, pass that
trace's base from `bases` as `y`, raised by enough to clear its tallest peak,
and leave `line` out.

### Step 6. Mark the peaks the cell does not account for

A peak that matched no reflection carries no hkl, and is exactly the peak a
reader should be looking at. `mark_peaks` writes an asterisk above each
position given.

A peak counts as indexed if any reflection falls within the tolerance, and the
tolerance `index_and_refine` works at is generous enough that a pattern of a
single well behaved phase often leaves nothing unindexed at all.

```python
unindexed = [entry.peak.two_theta for entry in indexed if not entry.is_indexed]
print(f"{len(unindexed)} peaks unaccounted for: {[round(v, 3) for v in unindexed]}")
```

```
0 peaks unaccounted for: []
```

`is_indexed` is a property of each `IndexedPeak`, so the list comprehension
above is the whole of the selection. To see which peaks the cell really
accounts for, index again with `index_peaks` against the refined cell at a
tighter tolerance, passing the zero offset the refinement found, and mark what
is left over.

```python
from xrdkit import index_peaks, mark_peaks

tight = index_peaks(
    clean, fit.cell, scan.wavelength, tolerance=0.02, zero_offset=fit.zero_offset
)
unindexed = [entry.peak.two_theta for entry in tight if not entry.is_indexed]

fig, ax = plot_pattern(scan, scale="sqrt", normalise=True)
ax.set_xlim(10.0, 80.0)
ax.set_ylim(0.0, 1.35)
annotate_hkl(ax, tight, y=1.02, min_relative_intensity=10.0)
markers = mark_peaks(ax, unindexed, y=1.02)
save_figure(fig, "figures/pattern_star")
print(f"{len(markers)} peaks unaccounted for: {[round(v, 3) for v in unindexed]}")
```

```
1 peaks unaccounted for: [26.921]
```

### The tunables

| Setting | Where | What it does |
| --- | --- | --- |
| `scale` | `plot_pattern`, `plot_stacked` | `"linear"`, `"sqrt"` or `"log"`. Use `"sqrt"` for anything where weak reflections matter |
| `normalise` | `plot_pattern`, `plot_stacked` | Scale each trace to a maximum of 1. Off by default for a single pattern, on for a stack |
| `offset` | `plot_stacked` | Vertical spacing between traces. `None` gives 1.2 times the tallest scaled trace |
| `figsize` | `plot_stacked` | Figure size in inches, 3.5 by 5.0 by default. A single pattern is 3.5 by 2.6 |
| `min_relative_intensity` | `annotate_hkl` | Skip peaks below this percentage of the strongest. 5 by default; raise it to thin out a crowded figure |
| `rotation` | `annotate_hkl` | Label angle in degrees anticlockwise, 90 by default, because upright labels take the least horizontal room |
| `max_levels` | `annotate_hkl` | How many rows of labels to try. 1 by default, and ignored when a `line` is given |
| `ambiguous` | `annotate_hkl` | What to do with a peak that matched more than one reflection: `"first"` labels the assigned one, `"all"` joins every candidate with a solidus, `"skip"` leaves it unlabelled |
| `formats`, `dpi` | `save_figure` | Which files to write, and at what resolution. `("png", "pdf")` at 300 dpi by default |

### Reading the unindexed peaks

An asterisk is a question, not an answer. There are four common explanations,
and they are told apart by where the peak sits and how strong it is.

A secondary phase gives peaks that do not move with composition across a series
and that keep a fixed ratio to one another. Look for the strongest two or three
reflections of the phase you suspect, not just one, and check them all.

K beta gives a faint copy of a strong reflection at lower angle. Where the
filter or the monochromator is imperfect, each strong peak has a companion at
the position its d spacing gives at the K beta wavelength, about 1.392 angstrom
for copper, which is some ten per cent shorter than K alpha 1 and so puts the
companion well below its parent.

Tungsten L lines come from a contaminated or an aged tube. They appear at fixed
angles that do not move with the sample, so the same stray peaks turn up in
every pattern measured on that instrument, which is what identifies them.

K alpha 2 satellites sit just above their parents and carry roughly half the
intensity. These should have been removed by `exclude_kalpha2` before indexing;
one that survives is a sign that the flagging tolerance or the intensity range
did not suit the pattern, rather than a sign of anything wrong with the sample.

## 5. Workflow 2: lattice parameters and theoretical density

### 5.1 Which route

Route A refines the cell from the peak positions alone, by least squares, with
`xrdkit.lattice.refine_lattice`. It needs nothing but the indexed peaks that
Workflow 1 already produced, so it runs in under a second and without GSAS-II.
Because it sees only the positions of the peaks the finder resolved, it cannot
use an overlapped reflection, and its esds describe the scatter of those
positions about the model rather than the agreement of a calculated pattern
with the measurement. That is enough to track a trend across a composition
series, where every sample is treated alike and what matters is how a and c
move from one sample to the next, and it is enough to decide whether a sample
is worth a longer scan.

Route B extracts an intensity for every reflection from the whole measured
profile, a Le Bail fit in GSAS-II, and refines the cell against that. It uses
the counts between the peaks as well as the peaks themselves, so overlapped
reflections and weak high angle reflections both contribute, and it refines the
zero point alongside the cell within one model rather than as a correction
applied beforehand. Route B gives the values to publish and the values to
compute a density from. Take Route A when the answer is a trend; take Route B
when the answer is a number that goes in a table, and always when a density
follows from it.

### 5.2 One time setup for Route B

Route B needs GSAS-II installed and an instrument parameter file for the
diffractometer. Both are done once, not once per sample.

Install GSAS-II from its own home page,
<https://advancedphotonsource.github.io/GSAS-II-tutorials>, which is a separate
Python installation and not a dependency of xrdkit. Point the kit at
it with the environment variables `XRDKIT_GSAS2_PYTHON` and `XRDKIT_GSAS2_HOME`,
or install it at `~/gsas2main`, which is where the kit looks by default. Check
that it is found before anything else.

```python
from xrdkit import find_gsas2

install = find_gsas2()
print(f"GSAS-II Python: {install.python}")
print(f"GSAS-II home:   {install.home}")
```

```
GSAS-II Python: C:\Users\amirk\gsas2main\python.exe
GSAS-II home:   C:\Users\amirk\gsas2main\GSAS-II
```

The instrument parameter file can be made with xrdkit alone, and no GSAS-II
session of your own is needed. It is made in two steps. First a starting file
is written from a Caglioti fit to the measured widths of the standard, which
`write_instprm` converts into the units GSAS-II uses. Fit each reflection of
the LaB6 scan as a K alpha doublet with `fit_profile`, weight the widths by
their esds, and fit U, V and W with `fit_caglioti`.

```python
import numpy as np

from xrdkit import (
    exclude_kalpha2,
    find_peaks,
    fit_caglioti,
    fit_profile,
    read_xrdml,
    write_instprm,
)

standard = read_xrdml("data/standards/lab6.xrdml")
peaks = exclude_kalpha2(find_peaks(standard, two_theta_range=(10.0, 98.0)))

angles, widths, weights = [], [], []
for peak in peaks:
    profile = fit_profile(
        standard.two_theta, standard.intensity, peak.two_theta, peak.fwhm
    )
    if profile.converged:
        angles.append(profile.two_theta)
        widths.append(profile.fwhm)
        weights.append(1.0 / profile.esd_fwhm**2)

caglioti = fit_caglioti(np.array(angles), np.array(widths), np.array(weights))
write_instprm("data/standards/start.instprm", caglioti)
print(f"{caglioti.n_peaks} reflections fitted")
print(
    f"U = {caglioti.u:.4f}, V = {caglioti.v:.4f}, W = {caglioti.w:.4f} degrees squared"
)
print(f"rms of the width fit {caglioti.rms:.5f} degrees")
```

```
14 reflections fitted
U = 0.0107, V = -0.0122, W = 0.0085 degrees squared
rms of the width fit 0.00254 degrees
```

Second, that starting file is refined against the standard's own structure, so
that the Lorentzian terms X and Y and the asymmetry SH/L are measured rather
than guessed. `standard_stages` gives the usual sequence, and the cell stage is
dropped because the certified lattice parameter of the standard is held: NIST
SRM 660c LaB6 is a = 4.156826 angstrom. The standard is taken to contribute no
broadening of its own, which is what the size of 10 micrometres and the
microstrain of zero say; left at the GSAS-II defaults of 1 micrometre and 1000
microstrain, the instrument terms would absorb the difference.

```python
import shutil
from pathlib import Path

from xrdkit import build_refine_job, run_job, standard_stages

stages = [stage for stage in standard_stages() if stage["name"] != "cell"]
print([stage["name"] for stage in stages])

job = build_refine_job(
    Path("results/instrument/lab6.gpx").resolve(),
    stages,
    data_file=Path("data/standards/lab6.xrdml").resolve(),
    instprm=Path("data/standards/start.instprm").resolve(),
    phases=[
        {
            "cif": Path("cifs/lab6.cif").resolve(),
            "name": "LaB6",
            "cell": [4.156826] * 3 + [90.0] * 3,
        }
    ],
    limits=(10.0, 98.0),
    cycles=10,
    broadening={"LaB6": {"size": 10.0, "mustrain": 0.0, "lgmix": 0.0}},
    export_prefix=Path("results/instrument/lab6").resolve(),
)
result = run_job(job, "results/instrument/gsas2_work")

last = result["stages"][-1]
print(f"Rwp {last['rwp']:.3f} per cent, GOF {last['gof']:.3f}")
for key in ("Zero", "U", "V", "W", "X", "Y", "SH/L"):
    entry = result["final"]["instrument"][key]
    print(f"  {key:5s} {entry['value']:11.6g}  esd {entry['esd']:.3g}")
shutil.copyfile(result["exports"]["instprm"], "data/standards/aeris.instprm")
```

```
['background and scale', 'zero', 'U V W', 'X Y', 'SH/L']
Rwp 7.160 per cent, GOF 1.909
  Zero   -0.0234803  esd 0.000283
  U         21.8437  esd 1.93
  V        -26.2779  esd 2.25
  W         10.8037  esd 0.651
  X         4.56529  esd 0.151
  Y         -1.5982  esd 0.35
  SH/L    0.0220585  esd 0.000396
```

That run takes a few seconds. The refined file is the one the run exports,
`result["exports"]["instprm"]`, not the starting file, so copy it to where the
sample refinements will read it from. Note every path handed to
`build_refine_job` above is made absolute with `resolve`. `run_job` runs the
driver inside the working directory it is given, so a relative path in the job
would be looked for under that directory and not found.

A negative Y, as here, means the Lorentzian width the data want is smaller than
the terms can describe together. It is not fatal, and the file is usable, but it
is worth refining a second variant with Y held at zero and adopting that
instead when its Rwp is no worse by more than a few tenths of a percentage
point, since one fewer parameter and a fixed rather than a negative Y is the
simpler description of the same widths.

The instrument parameter file holds the wavelengths and their intensity ratio,
the polarisation, the zero point, the Gaussian width terms U, V and W, the
Lorentzian terms X and Y, and the asymmetry SH/L: that is, everything about the
peak shape that belongs to the diffractometer rather than to the sample. It is
held fixed in every sample refinement afterwards. The reason is that a sample's
own size and microstrain broadening and the instrument's resolution have the
same kind of angular dependence, so refining both on one pattern lets each take
up whatever the other leaves. Measuring the instrument once on a standard that
is known to be sharp fixes the instrument half, and whatever width is left over
in a sample pattern is then the sample's.

### 5.3 Route A step by step

Route A starts from the indexed peaks of Workflow 1. The function is
`refine_lattice` in `xrdkit.lattice`, which refines a and c together with the
zero point and, if asked, a specimen displacement, and reports an esd for each.
It is not `refine_cell`, which belongs to `xrdkit.indexing`: that one fits a
cell to indexed peaks by linear least squares and reports no esds, and is what
the indexing cycles use internally.

Refine twice, once with the zero point held at zero and once with it free, and
compare.

```python
from xrdkit import (
    TTB_CELL,
    exclude_kalpha2,
    find_peaks,
    index_and_refine,
    read_xrdml,
    refine_lattice,
)

scan = read_xrdml("data/raw/sample.xrdml")
clean = exclude_kalpha2(find_peaks(scan, two_theta_range=(10.0, 80.0)))
indexed, cell_fit = index_and_refine(
    clean, start_cell=TTB_CELL, wavelength=scan.wavelength
)

held = refine_lattice(indexed, scan.wavelength, TTB_CELL, fit_zero=False)
refined = refine_lattice(indexed, scan.wavelength, TTB_CELL, fit_zero=True)

for name, fit in (("zero held at zero", held), ("zero refined", refined)):
    print(f"{name}:")
    print(f"  a = {fit.a:.4f} +/- {fit.esd_a:.4f} angstrom")
    print(f"  c = {fit.c:.4f} +/- {fit.esd_c:.4f} angstrom")
    print(f"  rms {fit.rms_two_theta:.4f} degrees on {fit.n_peaks} peaks")
print(
    f"zero {refined.zero:.4f} +/- {refined.esd_zero:.4f} degrees, "
    f"{refined.zero / refined.esd_zero:.0f} times its esd"
)
```

```
zero held at zero:
  a = 12.4400 +/- 0.0040 angstrom
  c = 3.9224 +/- 0.0017 angstrom
  rms 0.0554 degrees on 28 peaks
zero refined:
  a = 12.4777 +/- 0.0011 angstrom
  c = 3.9320 +/- 0.0003 angstrom
  rms 0.0069 degrees on 28 peaks
zero 0.1611 +/- 0.0041 degrees, 40 times its esd
```

Two things say that this zero point is real rather than a parameter absorbing
noise. The residual falls by a factor of eight, from 0.0554 to 0.0069 degrees,
which is far more than one extra parameter can buy on twenty eight peaks; and
the offset itself is forty times its own esd, so it is nowhere near zero. Where
a zero point comes out at one or two times its esd and the residual barely
moves, it is not measuring anything and is better held.

Notice also how far the cell moves: holding the zero at zero puts a at 12.4400
angstrom, and refining it puts a at 12.4777, a shift of 0.038 angstrom, which
is thirty times the esd of the better fit. This is why the data quality table
asks for either an internal standard or a refined displacement term. An
uncorrected shift of this size does not change a single hkl label, and it ruins
every lattice parameter.

`refine_lattice` can also refine a specimen displacement, with
`fit_displacement=True` and a goniometer radius in millimetres. Do not refine a
displacement and a zero point together on one scan unless the peaks span a wide
range of two theta: a zero point is a constant and a displacement follows
cos(theta), and over a short range the two cannot be told apart.

### 5.4 Route B step by step

A Le Bail refinement is described as a list of stages, each adding flags to the
ones before it. The stage that switches `le_bail` on extracts an intensity for
every reflection before its least squares, and the later stages keep the
extraction going while they refine the zero point, the cell and the crystallite
size. Microstrain is left out here, which is the safer default: size and
microstrain both broaden the peaks and are separated only by how that broadening
grows with angle, so refining both on one laboratory scan often gives two
numbers that trade against each other. Add `{"name": "microstrain", "mustrain":
True}` at the end when the data are good enough to want it.

The instrument parameters are not in the list, so they stay as the instrument
parameter file has them.

```python
from pathlib import Path

from xrdkit import build_refine_job, run_job

stages = [
    {
        "name": "background and scale",
        "background": {"type": "chebyschev-1", "terms": 6},
        "scale": True,
        "le_bail": True,
    },
    {"name": "zero", "zero": True},
    {"name": "cell", "cell": True},
    {"name": "size", "size": True},
]

job = build_refine_job(
    Path("results/lebail/sample.gpx").resolve(),
    stages,
    data_file=Path("data/raw/sample.xrdml").resolve(),
    instprm=Path("data/standards/aeris.instprm").resolve(),
    phases=[
        {
            "cif": Path("cifs/ttb.cif").resolve(),
            "name": "TTB",
            "cell": [refined.a, refined.a, refined.c, 90.0, 90.0, 90.0],
        }
    ],
    limits=(10.0, 98.0),
    cycles=10,
    le_bail_cycles=10,
    max_passes=60,
    pass_tolerance=0.1,
    broadening={"TTB": {"size": 1.0, "mustrain": 0.0, "lgmix": 1.0}},
    export_prefix=Path("results/lebail/sample").resolve(),
)
result = run_job(job, "results/lebail/gsas2_work")
print("completed:", result["completed"], " model from:", result["final_from"])
```

```
completed: True  model from: size
```

The cell of the phase is replaced with the Route A cell, so the refinement
starts from the best answer already available rather than from the cell of the
CIF, which came from a different composition. The limits stop at 98 degrees
because the scan ends at 99.98 and the last doublet is cut by the end of it.
`max_passes` and `pass_tolerance` matter for a Le Bail refinement in particular:
GSAS-II stops each one at the first cycle that raises chi squared, well short of
the minimum, so each stage is refined again and again until no parameter moves
by more than a tenth of an esd in a pass. That is why the zero stage below needed
thirty nine passes where the first stage needed four, and it is most of the run
time. On a laboratory scan of four thousand points the whole job takes about two
minutes.

Read the stages back to see where the fit improved and whether each one settled.

```python
for stage in result["stages"]:
    print(
        f"{stage['name']:22s} Rwp {stage['rwp']:6.3f}  GOF {stage['gof']:5.3f}  "
        f"variables {stage['n_variables']:2d}  passes {len(stage['passes']):2d}  "
        f"{stage['status']}"
    )
```

```
background and scale   Rwp 37.827  GOF 7.865  variables  7  passes  4  clean
zero                   Rwp  9.966  GOF 3.053  variables  8  passes 39  clean
cell                   Rwp  9.841  GOF 3.019  variables 10  passes  8  clean
size                   Rwp  7.945  GOF 2.577  variables 11  passes 10  clean
```

Every stage here is clean, meaning it settled and raised no sanity flag. A stage
recorded as unsettled reached the pass cap while a parameter was still moving by
more than the tolerance. That is not in itself an error, and the stage is kept,
but the esds of a stage that was still moving understate the real uncertainty,
so a cell taken from an unsettled stage should be treated as provisional.
Raising `max_passes` is the first thing to try; a stage that will not settle
however many passes it is given is usually refining something the data do not
determine.

The cell and the crystallite size come from the final model.

```python
phase = result["final"]["phases"][0]
a, c = phase["cell"]["length_a"], phase["cell"]["length_c"]
esd_a, esd_c = phase["cell_esd"]["length_a"], phase["cell_esd"]["length_c"]
size = phase["size"]
print(f"a = {a:.4f} +/- {esd_a:.4f} angstrom")
print(f"c = {c:.4f} +/- {esd_c:.4f} angstrom")
print(
    f"V = {phase['cell']['volume']:.3f} +/- {phase['cell_esd']['volume']:.3f} cubic angstrom"
)
print(f"crystallite size {size['value'] * 1000:.0f} +/- {size['esd'] * 1000:.0f} nm")
```

```
a = 12.4743 +/- 0.0004 angstrom
c = 3.9301 +/- 0.0002 angstrom
V = 611.551 +/- 0.058 cubic angstrom
crystallite size 181 +/- 4 nm
```

The size GSAS-II reports is in micrometres, which is why it is multiplied by a
thousand above. It is an apparent size, the whole of the sample broadening
attributed to small crystallites because microstrain was not refined, so read it
as an upper bound on the broadening rather than as a measurement of the
crystallites.

Look at the fit before believing any of it. The job exports the fitted pattern
and the reflection positions as CSV files, and `plot_rietveld` draws them with
the observed points, the calculated curve, the background, the reflection ticks
and the difference. Passing the result as well writes Rwp, the goodness of fit
and the refined cell into the figure.

```python
from xrdkit import apply_style, plot_rietveld, save_figure

apply_style()
fig, ax = plot_rietveld(
    result["exports"]["histogram"],
    result["exports"]["reflections"],
    sqrt_scale=True,
    result=result,
)
print(save_figure(fig, "figures/lebail_fit"))
```

```
[WindowsPath('figures/lebail_fit.png'), WindowsPath('figures/lebail_fit.pdf')]
```

What to look for is the shape of the difference curve rather than its size. In
this fit the largest excursions sit on the two strongest reflections, at 29.66
and 32.26 degrees, and each is a swing from negative to positive across the
peak. A difference that changes sign across a peak like that is a profile or
position mismatch, the calculated peak being a little too wide or a little
displaced, and it is what limits Rwp here. A difference that is one sided, a
lump of unexplained intensity where the ticks show no reflection, would mean
something else entirely: a second phase. The two are easy to tell apart in the
figure and impossible to tell apart from Rwp alone, which is why the figure is
drawn before the numbers are quoted, not after.

### 5.5 Theoretical density

The density module takes the refined cell, a composition and the number of
formula units per cell, and gives the cell volume, the mass of one formula unit
and the X-ray density. The composition is a dictionary of element symbol to the
number of atoms of that element in one formula unit, and not a formula string:
there is no formula parser in the kit, and passing a string such as
`"Sr0.4Ba0.5"` raises a `ValueError` about the individual characters. Coefficients
need not be integers. A composition that is empty, that names an element the kit
has no atomic mass for, or that has a negative coefficient raises a `ValueError`
naming the problem.

For the tetragonal tungsten bronze Sr0.40Ba0.50La0.10Nb1.90Ti0.10O6 there are
five formula units per cell.

```python
from xrdkit import TetragonalCell, cell_volume, formula_mass, theoretical_density

composition = {"Sr": 0.40, "Ba": 0.50, "La": 0.10, "Nb": 1.90, "Ti": 0.10, "O": 6.0}
cell = TetragonalCell(a=a, c=c)
volume, esd_volume = cell_volume(cell, esd_a, esd_c)
mass = formula_mass(composition)
density, esd_density = theoretical_density(composition, 5, volume, esd_volume)
print(f"M = {mass:.3f} g/mol per formula unit")
print(f"V = {volume:.3f} +/- {esd_volume:.3f} cubic angstrom")
print(f"theoretical density {density:.4f} +/- {esd_density:.4f} g/cm3")
```

```
M = 394.905 g/mol per formula unit
V = 611.551 +/- 0.050 cubic angstrom
theoretical density 5.3614 +/- 0.0004 g/cm3
```

The volume `cell_volume` returns agrees with the volume GSAS-II reports, 611.551
cubic angstrom, but its esd is a little smaller, 0.050 against 0.058, because
`cell_volume` propagates the esds of a and c as though the two were
uncorrelated, which the fits here give no covariance to do better with. Use
GSAS-II's own volume esd where one is available and the difference matters.

A relative density compares the measured density of a pellet with the
theoretical one. Suppose the Archimedes measurement gave 5.15 plus or minus 0.02
g/cm3. The relative uncertainty of the ratio is the quadrature sum of three
terms: twice the relative error in a, once the relative error in c, and the
relative error in the Archimedes measurement. The density goes as the reciprocal
of the volume and the volume is a squared times c, so a relative error in a
enters twice over and one in c once. Where a and c carry the same relative
error, those two terms come to about three times it, which is the rule of thumb
quoted in Section 2.

```python
import math

archimedes, esd_archimedes = 5.15, 0.02
relative = 100.0 * archimedes / density
uncertainty = relative * math.sqrt(
    (2.0 * esd_a / a) ** 2 + (esd_c / c) ** 2 + (esd_archimedes / archimedes) ** 2
)
print(f"relative density {relative:.2f} +/- {uncertainty:.2f} per cent")
print(f"  2 esd_a / a      {2.0 * esd_a / a:.2e}")
print(f"  esd_c / c        {esd_c / c:.2e}")
print(f"  esd_rho / rho    {esd_archimedes / archimedes:.2e}")
```

```
relative density 96.06 +/- 0.37 per cent
  2 esd_a / a      7.18e-05
  esd_c / c        3.93e-05
  esd_rho / rho    3.88e-03
```

The three terms are worth reading. The lattice parameters contribute about
eight parts in a hundred thousand between them and the Archimedes measurement
about four parts in a thousand, so the Archimedes measurement dominates by a
factor of about forty seven and the uncertainty on the relative density is
essentially its uncertainty alone. Once the cell is refined to this precision, more diffraction
time buys nothing; a better balance, or more repeats of the weighing, is what
improves the relative density.

### 5.6 Judging the result

Take the following in order before quoting a cell or a density.

Check Rwp and the goodness of fit against the band the data deserve. A Le Bail
fit has an intensity free for every reflection, so it should fit better than a
Rietveld refinement of the same pattern: an Rwp of a few per cent and a goodness
of fit approaching 1 on a long scan. A goodness of fit well above 2, as in the
example above, usually means the counting statistics are better than the model,
which is a fair description of a scan whose peaks are sharper than the
instrument file says or whose background has structure the polynomial cannot
follow.

Check that the zero point and the specimen displacement were not both refined on
one scan. They are nearly the same parameter over a short angular range, and a
fit that refines both will report small esds for two numbers that are trading
against each other.

Check the cell esds. For a density, a and c are wanted to about 0.001 angstrom,
and a Le Bail fit on a scan that reaches past 90 degrees should give esds at or
below that. An esd of several thousandths is a sign that the high angle
reflections carried too little weight, which is a counting problem rather than a
refinement problem.

Check that Route A and Route B agree. They are different methods on the same
peaks, so they should land within a few esds of one another once both have a
zero point. In the example, Route A gave a = 12.4777(11) and Route B a =
12.4743(4), a difference of 0.0034 angstrom, which is three times the Route A
esd and eight times the Route B esd. A gap of that size is expected and is
itself informative: Route A weights the resolved peaks it found, Route B weights
the whole profile, and the two disagree by rather more than either esd suggests.
Quote Route B, and treat the gap, not the esd, as the honest measure of how well
the cell is known. A gap of a whole per cent means one of the two indexed
something wrongly.

Check that no secondary phase has been swept into the background. A Chebyshev
polynomial with six terms is flexible enough to follow a broad hump, and a weak
second phase under the main pattern will be partly absorbed by it while the
residuals stay respectable. Look along the tick row for intensity that no tick
accounts for, on the square root scale where a weak phase is visible at all.

Finally, know when to remeasure rather than refine. A scan that stops at 80
degrees cannot give a publishable cell, because the reflections that pin it down
were never measured, and no refinement recovers them. A scan whose high angle
peaks are only two or three times the background gives esds that are honest and
too large, and the answer is a longer count and not more parameters. A sample
mounted proud of the holder, with no internal standard and no displacement term,
gives a cell that is precise and wrong, and remounting it flush takes less time
than arguing with the numbers.

## 6. Workflow 3: Rietveld refinement

### 6.1 What Rietveld adds, and when to do it

A Le Bail fit lets every reflection take whatever intensity fits best. That is
why it is so good at cells: the positions and the widths of the peaks are
modelled and the intensities are simply granted, so nothing about the structure
can bias the cell. It is also why a Le Bail fit tests nothing. A pattern from
the wrong structure in the right cell fits a Le Bail model exactly as well as
the right one does.

A Rietveld refinement computes every intensity from a structural model: which
atoms sit on which sites, at what coordinates, with what occupancies and what
thermal displacement. The intensities are then predictions, and the fit is a
test of the model. That is what makes it worth doing, and it is also what makes
it demanding. Getting an occupancy or a coordinate out of a powder pattern means
the calculated intensities have to be wrong in a way the data can see, which
needs the counting statistics, the angular range and the profile description set
out in the Rietveld column of Section 2. Refine one or two samples per series
this way, chosen because the question needs a structure, and use Workflow 2 for
everything else.

### 6.2 What you need

Four things, beyond the raw scan.

The Le Bail result of Workflow 2 supplies the starting cell and confirms that
the instrument parameter file describes the peak shapes. Start a Rietveld
refinement from a cell that has already been refined against the whole pattern,
not from the cell in the CIF, which came from somebody else's composition.

A CIF of the reference structure supplies the sites, their Wyckoff positions and
the starting coordinates. It does not have to be the same composition as the
sample; it has to be the same structure type in the same space group.

The nominal composition, as atoms per formula unit, says what the sample is
meant to be. The kit puts it on the sites for you rather than making you edit
occupancies by hand.

A decision about which sites share what. This is the part that cannot be
automated, because it is the chemistry. Four decisions are needed, and the kit
takes them as a structure description.

Site kinds group the sites into the families that are refined together. For the
tetragonal tungsten bronze here there are three: A for the two cation channel
sites, B for the two niobium sites, and O for the five oxygens. Coordinates are
freed one kind at a time.

Uiso groups say which sites share one thermal displacement parameter. A powder
pattern cannot support an independent Uiso on nine sites, so the A site cations
share one, the two niobium sites share another, and the five oxygens share a
third. Three parameters instead of nine.

The origin site says which coordinate is held to stop the structure sliding.
P4bm is polar along c, so nothing in the symmetry fixes where the origin sits on
that axis; every atom can move together along c without changing the pattern at
all. Holding the z of one site fixes it. Here that is Nb1, on the special 2b
position.

The composition rule says how the nominal composition is distributed. Elements
the CIF already holds keep the CIF's distribution over their sites and are
scaled by one factor each. Elements the CIF lacks are added onto a host element's
sites in proportion to the host's occupancy: lanthanum goes where strontium is,
and titanium where niobium is.

The kit reads all of this from a TOML file with `load_config`, which checks it as
it reads and raises `ConfigError` naming the table and the key on the first
thing that is wrong, including a composition that will not fit on the sites. A
minimal `config/samples.toml` for this sample and structure is below.

```toml
[samples.x10]
id = "10"
scan = "data/raw/sample.xrdml"
composition = { Sr = 0.40, Ba = 0.50, La = 0.10, Nb = 1.90, Ti = 0.10, O = 6.0 }
structure = "ttb"
start_cell = { a = 12.4743, c = 3.9301 }
two_theta = [17.0, 98.0]
background = { function = "chebyschev-1", terms = 6 }
refine_microstrain = false
notes = "x = 0.10 calcined powder, Rietveld through the coordinates."

[structures.ttb]
cif = "cifs/ttb.cif"
label = "COD 2100720"
phase_name = "TTB"
space_group = "P4bm"
formula_units = 5
sites = [
    { atoms = { Ba2 = "Ba", Sr2 = "Sr" }, wyckoff = "4c", kind = "A" },
    { atoms = { Sr1 = "Sr" }, wyckoff = "2a", kind = "A" },
    { atoms = { Nb1 = "Nb" }, wyckoff = "2b", kind = "B" },
    { atoms = { Nb2 = "Nb" }, wyckoff = "8d", kind = "B" },
    { atoms = { O1 = "O" }, wyckoff = "4c", kind = "O" },
    { atoms = { O2 = "O" }, wyckoff = "8d", kind = "O" },
    { atoms = { O3 = "O" }, wyckoff = "8d", kind = "O" },
    { atoms = { O4 = "O" }, wyckoff = "2b", kind = "O" },
    { atoms = { O5 = "O" }, wyckoff = "8d", kind = "O" },
]
free_coordinates = { "8d" = "xyz", "4c" = "xz", "2a" = "z", "2b" = "z" }
uiso_groups = [
    { name = "A site cations", sites = ["Ba2", "Sr1"] },
    { name = "Nb/Ti", sites = ["Nb1", "Nb2"] },
    { name = "O", sites = ["O1", "O2", "O3", "O4", "O5"] },
]
origin = { site = "Nb1", axis = "z" }
exchange = { elements = ["Sr", "Ba"], sites = ["Ba2", "Sr1"] }

[structures.ttb.composition]
added = { La = "Sr", Ti = "Nb" }
```

`free_coordinates` says what each Wyckoff position of this space group leaves
free, so that the kit never refines a coordinate that symmetry fixes.
`exchange` names the elements the occupancy stage trades and the sites it trades
them between. `start_cell` here is the cell Workflow 2 refined.

```python
from xrdkit import load_config, sample_settings

config = load_config("config/samples.toml")
sample, structure = sample_settings(config, "10")
print(f"sample {sample['name']}, scan {sample['scan']}")
print(
    f"phase {structure['phase_name']} in {structure['space_group']}, "
    f"Z = {structure['formula_units']}"
)
print(
    "sites:",
    ", ".join(f"{s['name']} ({s['wyckoff']}, {s['kind']})" for s in structure["sites"]),
)
print("Uiso groups:", ", ".join(g["name"] for g in structure["uiso_groups"]))
print(f"origin held: {structure['origin']['site']} {structure['origin']['axis']}")
```

```
sample x10, scan data/raw/sample.xrdml
phase TTB in P4bm, Z = 5
sites: Ba2 (4c, A), Sr1 (2a, A), Nb1 (2b, B), Nb2 (8d, B), O1 (4c, O), O2 (8d, O), O3 (8d, O), O4 (2b, O), O5 (8d, O)
Uiso groups: A site cations, Nb/Ti, O
origin held: Nb1 z
```

### 6.3 The stage sequence, and why it is ordered this way

A refinement is a list of stages, and each stage adds flags to the ones before
it, so by the last stage everything named along the way is refining together.
The order is not a matter of taste: each stage frees parameters that only make
sense once the ones before them are near right.

The profile comes first: background, scale, zero point, cell and crystallite
size. This puts the calculated peaks in the right places with roughly the right
heights and widths. Nothing structural can be judged until they are, because a
peak in the wrong place produces intensity errors that look exactly like
occupancy errors. What can go wrong here is the background taking up intensity
that belongs to the peaks, which shows as a scale that keeps falling.

One overall Uiso comes next. Every atom gets the same thermal parameter. It is a
single number, it is always determined by the data, and it soaks up the overall
falling off of intensity with angle that would otherwise be pushed into the
occupancies. What can go wrong is that it comes out negative, which never means
the atoms are colder than still; it means intensity is missing at high angle, so
look at absorption, at the sample height, or at the instrument file.

Grouped Uiso follows: the three groups here instead of one overall value. Note
that the overall Uiso has to be switched off explicitly in the stage that
introduces the groups, with `"overall_uiso": False`, because the stages
accumulate and a phase cannot have both at once. `build_refine_job` checks the
whole sequence before GSAS-II starts and refuses it with a message saying so, so
this costs a second rather than a refinement.

Coordinates come next, one site kind at a time, heaviest scatterers first: the
niobium sites, then the A site cations, then the oxygens. Heavy atoms dominate
the intensities, so their positions are the best determined and moving them
first stops the light atoms from chasing errors that are not theirs. The stage
that frees the kind holding the origin also carries the origin hold. What goes
wrong here is that oxygen positions in a laboratory pattern are weakly
determined and strongly correlated with one another, so they wander.

Occupancies come last, under the composition constraint: the exchange elements
are traded between the exchange sites with each element's total content over
them held, so the refinement decides how strontium and barium are distributed
between the two A sites without being free to change how much of either the
sample contains. This is the stage that usually takes longest, and the one most
often left undetermined. The run below stops before it, but it is built the same
way, as `{"name": "A site occupancies", "occupancies": [exchange]}` with
`exchange` naming the phase, the sites and the elements.

#### Passes

GSAS-II stops a refinement at the first least squares cycle that raises chi
squared, which for a structural refinement is often well short of the minimum.
The kit therefore refines each stage repeatedly, a pass at a time, until nothing
but the scale moves by more than `pass_tolerance` esds from one pass to the next,
and at most `max_passes` times. A tolerance of 0.1 esd is the usual setting. A
stage that reaches the cap is not an error, but it is a fact about the result and
it is recorded.

#### Sanity rules and rollback

After every stage the atoms are checked. Three things raise a flag: a negative
Uiso, an occupancy outside 0 to 1 or a site whose occupancies add up to more than
1, and a coordinate that has moved more than `max_shift`, 0.05 fractional by
default, from where the job found it. A stage that raises a flag the stage before
it did not have is rejected by default: it is recorded in full, with
`rejected_because`, and then the project is put back to the state the last kept
stage left, with whatever that stage alone was refining held from then on. The
run carries on from there. Set `on_flagged` to `"accept"` to keep such a stage
instead.

#### The labels

A stage is clean when it settled and raised no new flag. Quote it.

A stage is rejected when it raised a sanity flag the stage before it did not
have. It is rolled back, as above, so it is not part of the model at all, and
`rejected` on the result names every stage this happened to. Its residuals are
still recorded, and they are often better than those of the stage that was kept,
which is the point: a rejected stage fitted the data by moving the structure
somewhere it should not be.

A stage is unsettled when it reached the pass cap while something was still
moving. It is kept, and `largest_remaining_move` names the parameter and how many
esds it shifted in the last pass. The values are usable but the esds understate
the truth, because they describe a minimum the refinement had not reached. Say so
if you quote them.

A parameter is undetermined when the data do not determine it: an occupancy whose
esd is more than half its allowed range, or a coordinate or Uiso whose esd is
larger than the distance it moved from where it started. This is independent of
the stage's status, and a clean stage can still leave undetermined parameters.
An undetermined value is not a result and should not go in a table of refined
parameters. Report it as held, or report the refinement without it.

### 6.4 Step by step

The structure has to be set up before it can be refined, and setting it up needs
the atoms as GSAS-II reads them from the CIF, which means creating the project
first. `create_phase` below does that with the `create` action, which builds the
project and reports the phase without refining anything. It is called twice:
once to see the CIF as it stands, and again with the occupancy edits that put the
nominal composition on the sites.

```python
from pathlib import Path

from xrdkit import (
    build_refine_job,
    cell_contents,
    composition_edits,
    run_job,
    site_setup,
)

phase_name = structure["phase_name"]


def create_phase(atoms=None):
    """The phase as GSAS-II reads it from the CIF, with ``atoms`` edits applied."""
    phase = {"cif": Path(structure["cif"]).resolve(), "name": phase_name}
    if atoms:
        phase["atoms"] = atoms
    job = build_refine_job(
        Path("results/rietveld/structure.gpx").resolve(),
        [{"scale": True}],
        data_file=Path(sample["scan"]).resolve(),
        instprm=Path("data/standards/aeris.instprm").resolve(),
        phases=[phase],
    )
    job["action"] = "create"
    (created,) = run_job(job, "results/rietveld/gsas2_work/create")["phases"]
    return created


from_cif = create_phase()
edits = composition_edits(from_cif["atoms"], sample["composition"], structure)
phase = create_phase(edits)
plan = site_setup(structure, phase["atoms"])

print("per cell, from the CIF:", cell_contents(from_cif["atoms"]))
print("per cell, as set up:   ", cell_contents(phase["atoms"]))
print("sites by kind:", {k: [s["name"] for s in v] for k, v in plan["kinds"].items()})
print("coordinates freed:", plan["coordinates"])
```

```
per cell, from the CIF: {'Ba': 2.6, 'Sr': 2.4008, 'Nb': 10.0, 'O': 30.0}
per cell, as set up:    {'Ba': 2.5, 'Sr': 2.0, 'Nb': 9.5, 'O': 30.0, 'La': 0.5, 'Ti': 0.5}
sites by kind: {'A': ['Ba2', 'Sr1'], 'B': ['Nb1', 'Nb2'], 'O': ['O1', 'O2', 'O3', 'O4', 'O5']}
coordinates freed: {'A': {'Ba2': 'xz', 'Sr1': 'z'}, 'B': {'Nb2': 'xyz'}, 'O': {'O1': 'xz', 'O2': 'xyz', 'O3': 'xyz', 'O4': 'z', 'O5': 'xyz'}}
```

Check the second line against the composition before going on. Five formula
units of Sr0.40Ba0.50La0.10Nb1.90Ti0.10O6 are Sr 2.0, Ba 2.5, La 0.5, Nb 9.5,
Ti 0.5 and O 30.0 per cell, which is what the edits produced, so the model is
the sample the furnace was asked for. Note also that `coordinates` lists no
entry for Nb1: it is the origin site, so its z is held rather than refined, and
`site_setup` has already left it out.

The stages are assembled from the plan.

```python
origin = {phase_name: {"site": plan["origin"]["name"], "axis": plan["origin_axis"]}}

stages = [
    {
        "name": "profile",
        "background": {
            "type": sample["background"]["function"],
            "terms": sample["background"]["terms"],
        },
        "scale": True,
        "zero": True,
        "cell": True,
        "size": True,
    },
    {"name": "overall Uiso", "overall_uiso": True},
    {
        "name": "Uiso groups",
        "overall_uiso": False,
        "uiso_groups": {phase_name: plan["uiso_groups"]},
    },
]
for kind in ("B", "A", "O"):
    stage = {"name": f"{kind} site coordinates"}
    if plan["coordinates"].get(kind):
        stage["coordinates"] = {phase_name: plan["coordinates"][kind]}
    if kind == plan["origin"]["kind"]:
        stage["origin"] = origin
    stages.append(stage)

print([stage["name"] for stage in stages])
```

```
['profile', 'overall Uiso', 'Uiso groups', 'B site coordinates', 'A site coordinates', 'O site coordinates']
```

Now build the job and run it. As in Section 5, every path is made absolute,
because the driver runs inside the working directory it is given. `bonds=True`
asks the run to add the final model's cation to anion distances to the result.
This refinement takes about three to four minutes on a laptop.

The scan refined here is a longer powder scan of the same composition as the one
Section 5 used, over the narrower range 17 to 98 degrees rather than 10 to 98, so
the residuals below are not comparable with the Le Bail ones in Section 5.4.
Compare a Rietveld refinement with a Le Bail fit only when both were run on the
same scan over the same range.

```python
start = sample["start_cell"]

job = build_refine_job(
    Path("results/rietveld/x10.gpx").resolve(),
    stages,
    data_file=Path(sample["scan"]).resolve(),
    instprm=Path("data/standards/aeris.instprm").resolve(),
    phases=[
        {
            "cif": Path(structure["cif"]).resolve(),
            "name": phase_name,
            "cell": [start["a"], start["a"], start["c"], 90.0, 90.0, 90.0],
            "atoms": edits,
        }
    ],
    limits=tuple(sample["two_theta"]),
    cycles=10,
    max_passes=30,
    pass_tolerance=0.1,
    overall_uiso_start={phase_name: 0.01},
    broadening={phase_name: {"size": 1.0, "mustrain": 0.0, "lgmix": 1.0}},
    export_prefix=Path("results/rietveld/x10").resolve(),
    bonds=True,
)
result = run_job(job, "results/rietveld/gsas2_work/refine")
print("completed:", result["completed"])
print("final model from:", result["final_from"])
print("rolled back:", result["rejected"])
```

```
completed: True
final model from: A site coordinates
rolled back: ['O site coordinates']
```

The run finished, but the last stage was rolled back, so the model it leaves is
the one the A site stage produced. Read the stages to see what happened.

```python
for stage in result["stages"]:
    print(
        f"{stage['name']:20s} Rwp {stage['rwp']:5.3f}  Rp {stage['rp']:5.3f}  "
        f"GOF {stage['gof']:5.3f}  RF2 {stage['residuals']['0:0:Rf^2']:5.2f}  "
        f"variables {stage['n_variables']:2d}  passes {len(stage['passes']):2d}  "
        f"{stage['status']}"
    )
```

```
profile              Rwp 4.594  Rp 3.552  GOF 1.589  RF2  5.70  variables 11  passes  6  clean
overall Uiso         Rwp 4.249  Rp 3.328  GOF 1.469  RF2  5.91  variables 12  passes  2  clean
Uiso groups          Rwp 4.173  Rp 3.271  GOF 1.444  RF2  6.92  variables 14  passes  2  clean
B site coordinates   Rwp 4.151  Rp 3.263  GOF 1.436  RF2  6.83  variables 17  passes  7  clean
A site coordinates   Rwp 4.155  Rp 3.262  GOF 1.438  RF2  6.48  variables 20  passes 30  unsettled
O site coordinates   Rwp 4.127  Rp 3.248  GOF 1.431  RF2  7.75  variables 32  passes 30  rejected
```

`stage_statuses` gives the same rows as data, with a reason on every stage that
is not clean, and `stage_status_table` gives them as markdown ready for a write
up. The reasons and the undetermined parameters are what matter.

```python
from xrdkit import stage_statuses

for row in stage_statuses(result):
    if row["reason"]:
        print(f"{row['name']}: {row['reason']}")
for entry in result["undetermined"]:
    print(f"undetermined: {entry['message']}")
```

```
A site coordinates: not settled in 30 passes; largest remaining move 9.50 esd (0::dAz:2)
O site coordinates: sanity check: O3 z moved +0.1019
undetermined: Ba2 x 0.17207, esd 0.00019 more than its shift 0.00004
undetermined: Ba2 y 0.67207, esd 0.00019 more than its shift 0.00004
undetermined: Sr2 x 0.17207, esd 0.00019 more than its shift 0.00004
undetermined: Sr2 y 0.67207, esd 0.00019 more than its shift 0.00004
undetermined: La2 x 0.17207, esd 0.00019 more than its shift 0.00004
undetermined: La2 y 0.67207, esd 0.00019 more than its shift 0.00004
```

The final model carries every atom as the last kept stage left it, with an esd on
each parameter that was refined and none on the ones that were not.

```python
final = result["final"]["phases"][0]
print(
    f"a = {final['cell']['length_a']:.4f} +/- {final['cell_esd']['length_a']:.4f} angstrom"
)
print(
    f"c = {final['cell']['length_c']:.4f} +/- {final['cell_esd']['length_c']:.4f} angstrom"
)
for atom in final["atoms"]:
    esds = atom["xyz_esd"] or [None, None, None]
    coordinates = ", ".join(
        f"{value:.5f}" if esd is None else f"{value:.5f}({esd * 1e5:.0f})"
        for value, esd in zip(atom["xyz"], esds)
    )
    print(
        f"  {atom['label']:4s} {atom['type']:3s} occ {atom['occupancy']:.4f}  "
        f"Uiso {atom['uiso']:.4f}({atom['uiso_esd'] * 1e4:.0f})  {coordinates}"
    )
```

```
a = 12.4738 +/- 0.0003 angstrom
c = 3.9320 +/- 0.0001 angstrom
  Ba2  Ba  occ 0.6250  Uiso 0.0382(12)  0.17207(19), 0.67207(19), 0.46429(506)
  Sr2  Sr  occ 0.2051  Uiso 0.0382(12)  0.17207(19), 0.67207(19), 0.46429(506)
  Sr1  Sr  occ 0.5898  Uiso 0.0382(12)  0.00000, 0.00000, 0.46303(620)
  Nb1  Nb  occ 0.9500  Uiso 0.0208(8)  0.50000, 0.00000, 0.01582
  Nb2  Nb  occ 0.9500  Uiso 0.0208(8)  0.07459(19), 0.21232(21), -0.02946(384)
  O1   O   occ 1.0000  Uiso 0.0516(29)  0.28318, 0.78318, 0.96860
  O2   O   occ 1.0000  Uiso 0.0516(29)  0.13908, 0.06924, 0.95850
  O3   O   occ 1.0000  Uiso 0.0516(29)  0.99357, 0.34369, 0.95880
  O4   O   occ 1.0000  Uiso 0.0516(29)  0.50000, 0.00000, 0.47880
  O5   O   occ 1.0000  Uiso 0.0516(29)  0.07542, 0.20469, 0.46600
  La2  La  occ 0.0513  Uiso 0.0382(12)  0.17207(19), 0.67207(19), 0.46429(506)
  La1  La  occ 0.1475  Uiso 0.0382(12)  0.00000, 0.00000, 0.46303(620)
  Ti1  Ti  occ 0.0500  Uiso 0.0208(8)  0.50000, 0.00000, 0.01582
  Ti2  Ti  occ 0.0500  Uiso 0.0208(8)  0.07459(19), 0.21232(21), -0.02946(384)
```

Every oxygen carries its CIF coordinates with no esd, because the stage that
would have refined them was rolled back. Nb1 has no esd on any coordinate: two of
them are fixed by the 2b site symmetry and the third is the origin hold. Atoms
that share a site share their coordinates and their Uiso exactly, as they should.

The bond lengths are the first physical check on a refined structure. A niobium
to oxygen distance in a tungsten bronze should lie between about 1.8 and 2.2
angstrom; anything well outside that means the model has gone somewhere it
should not be, whatever the residuals say.

```python
for bond in result["bonds"][phase_name]:
    if not bond["centre"].startswith("Nb"):
        continue
    esd = "held" if bond["esd"] is None else f"+/- {bond['esd']:.4f}"
    print(
        f"{bond['centre_site']:12s} to {bond['target_site']:4s}  "
        f"{bond['distance']:.4f} {esd:16s} angstrom  x{bond['count']}"
    )
```

```
Nb1/Ti1      to O4    1.8204 held             angstrom  x1
Nb1/Ti1      to O3    1.9643 held             angstrom  x4
Nb1/Ti1      to O4    2.1115 held             angstrom  x1
Nb2/Ti2      to O3    1.9259 +/- 0.0026       angstrom  x1
Nb2/Ti2      to O5    1.9505 +/- 0.0151       angstrom  x1
Nb2/Ti2      to O2    1.9582 +/- 0.0026       angstrom  x1
Nb2/Ti2      to O1    1.9821 +/- 0.0025       angstrom  x1
Nb2/Ti2      to O5    1.9861 +/- 0.0151       angstrom  x1
Nb2/Ti2      to O2    2.0139 +/- 0.0025       angstrom  x1
```

Every distance is inside the window, so the octahedra are still octahedra. The
distances from Nb1 are marked held because neither it nor the oxygens around it
were refined in the kept stages, so the numbers are those of the CIF.

Plot the fit, as always before quoting anything.

```python
from xrdkit import apply_style, plot_rietveld, save_figure

apply_style()
fig, ax = plot_rietveld(
    result["exports"]["histogram"],
    result["exports"]["reflections"],
    sqrt_scale=True,
    result=result,
)
print(save_figure(fig, "figures/rietveld_x10"))
```

```
[WindowsPath('figures/rietveld_x10.png'), WindowsPath('figures/rietveld_x10.pdf')]
```

One caution about that figure when a stage has been rolled back. The cell written
into it is the final model's, but the Rwp and the goodness of fit are those of the
last stage that did not raise an error, and a rejected stage did not raise an
error: it was rolled back. So this figure reports Rwp 4.13 per cent, from the O
site stage, while the curve it draws is the A site model whose Rwp is 4.155. Take
the residuals from the stage table, not from the figure, whenever `rejected` is
not empty.

When a stage fails outright rather than being rejected, `run_job` raises
`Gsas2Error` and there is no result to read. `failure_markdown` turns what there
is into a write up: the error, the stages that did run, and the tail of the
GSAS-II log. The snippet below asks for an instrument file that does not exist,
so that the failure is real.

```python
from xrdkit import Gsas2Error, failure_markdown, log_tail

broken = build_refine_job(
    Path("results/broken/x10.gpx").resolve(),
    [{"name": "profile", "scale": True}],
    data_file=Path(sample["scan"]).resolve(),
    instprm=Path("data/standards/missing.instprm").resolve(),
    phases=[{"cif": Path(structure["cif"]).resolve(), "name": phase_name}],
)
work = Path("results/broken/gsas2_work")
try:
    run_job(broken, work)
except Gsas2Error as error:
    report = failure_markdown(
        "x = 0.10, coordinates", None, str(error), log_tail(work / "refine.log")
    )
    Path("results/broken/failure.md").write_text(report, encoding="utf-8")
    headings = [line for line in report.splitlines() if line.startswith("## ")]
    print(f"{len(report.splitlines())} lines written, sections {headings}")
    print(next(line for line in report.splitlines() if "failed with exit code" in line))
```

```
104 lines written, sections ['## Error', '## Stages', '## GSAS-II log, last 40 lines']
GSAS-II job 'refine' failed with exit code 1:
```

### 6.5 Reading the outcome

Take the run above as it stands. Four stages came out clean: the profile, the
overall Uiso, the Uiso groups and the niobium coordinates. Rwp fell from 4.594 to
4.151 per cent across them and the goodness of fit from 1.589 to 1.436, so each
stage earned the parameters it added.

The A site coordinates stage came out unsettled. It ran the full thirty passes
with one parameter still moving by nine and a half esds a pass, and the six
undetermined entries say what that parameter was doing: the x and y of the A site
atoms have esds larger than the distance they moved. The stage is kept, and its
cell and Uiso are perfectly usable, but those two coordinates are not a result.
Raising `max_passes` is worth one attempt; the honest reading is that this
pattern does not determine where the A site cations sit in the plane.

The O site coordinates stage was rejected, and the reason names exactly what
happened: O3 moved 0.1019 fractional along z, twice the 0.05 the sanity check
allows. Note that Rwp fell when it did so, to 4.127. That is the whole argument
for having a sanity check at all: the residual improved and the structure got
worse, because twelve weakly determined oxygen coordinates can always find
somewhere to go that fits the noise a little better. The kit rolled the stage
back, held what it alone refined, and left the A site model as the answer.

Telling a data limitation from a model error is the next question, and the
signatures differ. A data limitation shows as parameters that are individually
undetermined while the fit as a whole is good: correlated coordinates that trade
against one another, mixed occupancies on a shared site that the scattering
contrast cannot separate, esds larger than shifts. That is this run. A model
error shows as a fit that is bad in a structured way: a whole class of
reflections fitted badly points at the wrong space group, intensity under peaks
with no tick mark at a missing phase, and a residual that is good at low angle
and bad at high angle, or the reverse, points at the composition or at the
thermal parameters rather than at the counting.

For a paper, report Rwp, Rp and the goodness of fit of the stage the model came
from, the refined cell with its esds, the wavelength and the range refined over,
which parameters were refined and which were held and why, and the status labels.
Say that a stage was unsettled if it was. Do not tabulate an undetermined
coordinate as a refined value; either hold it and say so, or leave that part of
the structure out of the claim. A referee can tell the difference, and the labels
are in the result precisely so that the decision is not left to memory.

Go back for better data when the limitation is the data. The Rietveld row of
Section 2 is the standard: monochromatic K alpha 1, 5 to 130 degrees, a strongest
peak above 20000 counts and a background above 200 counts per step at high angle,
powder sieved below 45 micrometres and loaded to limit preferred orientation. The
scan used here is a 34 second per step powder scan reaching 98 degrees with a
strongest peak near 10000 counts, which is a good Le Bail scan and a marginal
Rietveld one. That, and not the refinement strategy, is why the A site
coordinates would not settle.

### 6.6 Troubleshooting

| Symptom | Cause | What to do |
| --- | --- | --- |
| `Gsas2Error` saying the GSAS-II Python was not found | The kit looks at `XRDKIT_GSAS2_PYTHON`, then `XRDKIT_GSAS2_HOME`, then `~/gsas2main` | Run `find_gsas2()` on its own to see what it resolves to, and set the two environment variables to the GSAS-II Python executable and to the installation folder |
| GSAS-II fails to read the instrument parameter file | A blank line, or a trailing blank line, in the `.instprm`; GSAS-II reads it key by key and does not tolerate one | Rewrite it with `write_instprm`, or take the one a refinement exported, which is always well formed |
| `FileNotFoundError` naming a path under the working directory | A relative path in the job. `run_job` runs the driver inside the working directory it is given, so the job's paths are resolved from there | Pass every path through `Path(...).resolve()`, as in every snippet here |
| Every stage comes out unsettled | The pass cap is too low for the problem, or parameters are trading against one another | Raise `max_passes` once. If nothing changes, look at which parameter `largest_remaining_move` names on each stage; if it is the same one every time, that parameter is not determined and should be held |
| A negative Uiso on the first structural stage | Intensity missing at high angle, not cold atoms | Check the instrument file against a fresh standard scan, check the sample sat flush in the holder, and check the composition: too much heavy scattering in the model shows up this way |
| Rwp at fixed atoms far above the Le Bail value | Expected, but only by so much. A Le Bail fit has a free intensity per reflection, so it always fits better | A gap of one to two percentage points is normal. A gap of five or more means the structure is wrong, not merely imperfect: check the space group, the composition on the sites, and whether a second phase is present |

## 7. Known limitations of version 0.1.0

Six things the kit does not do yet, all of them met somewhere in this guide.
Each is on the list for the next release.

The reader accepts `.xrdml` and nothing else. There is no reader for two or
three column `.xy` or `.xye`, for Bruker `.raw` or `.brml`, or for `.gsas` or
`.fxye`. A scan in another format can still be used, because `XRDScan` is an
ordinary dataclass and every routine downstream of the reader takes one of those
rather than a file path. Load the columns yourself and fill the fields in.

```python
import numpy as np

from xrdkit import XRDScan, find_peaks

two_theta, intensity = np.loadtxt("data/raw/sample.xy", unpack=True)
scan = XRDScan(
    two_theta=two_theta,
    intensity=intensity,
    wavelength=1.540598,
    start_angle=float(two_theta[0]),
    end_angle=float(two_theta[-1]),
    step_size=float(two_theta[1] - two_theta[0]),
    time_per_step=0.0,
    sample_id="x10",
    source_path="data/raw/sample.xy",
)
print(
    f"{scan.sample_id}: {scan.start_angle:.2f} to {scan.end_angle:.2f} degrees, "
    f"{len(find_peaks(scan))} peaks"
)
```

```
x10: 10.01 to 99.98 degrees, 49 peaks
```

Indexing is tetragonal only, and `P4bm` is the only space group whose reflection
conditions `xrdkit.indexing` knows. Pass `space_group=None` to apply none, which
works for any tetragonal cell but labels reflections the conditions would have
forbidden. Nothing in the module handles a lower symmetry cell.

`standard_stages` is the instrument calibration sequence, background and scale,
zero, cell, U V W, X Y, SH/L, and not a sample refinement sequence. There is no
helper that builds a Le Bail or a Rietveld stage list, so those are written out
in full, as in Sections 5.4 and 6.4.

The caption `plot_rietveld` writes takes its Rwp and goodness of fit from the
last stage that has no error entry, and a stage that was rejected has none: it
was rolled back, not failed. The cell in the same caption comes from the final
model. So when a run has rejected a stage, the two halves of that caption come
from different stages, and the residuals should be read from the stage table
instead, as Section 6.4 says.

`formula_mass`, and through it `theoretical_density`, takes a dictionary of
element symbol to atoms per formula unit. It does not parse a formula string:
passing one raises a `ValueError` naming the individual characters, which is a
confusing way to be told that the argument was of the wrong kind.

Every path handed to `build_refine_job` must be absolute. `run_job` runs the
GSAS-II driver inside the working directory it is given, so a relative path in
the job is resolved from there and the file is not found. Passing each one
through `Path(...).resolve()`, as every snippet in Sections 5 and 6 does, is the
whole of the workaround.
