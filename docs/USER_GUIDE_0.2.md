# xrdkit user guide

This guide is being rewritten for 0.2.0 around the commands; Sections beyond 3 are still in [docs/USER_GUIDE.md](USER_GUIDE.md).

## 1. Introduction

xrdkit analyses laboratory powder X-ray diffraction patterns. It reads a scan
from the diffractometer, finds the peaks in it, indexes them against a cell,
refines lattice parameters, calculates a theoretical density, draws
publication quality figures, and drives a GSAS-II Rietveld refinement from a
description of the refinement written as data.

In 0.2.0 each piece of that is a command you type. There is a library
underneath, and Part II documents it for anyone who wants to build something
of their own, but nothing in Part I needs a line of Python. A command reads
your scan, writes its results beside it, and prints what it wrote.

The kit never modifies raw data. A scan file is opened for reading and the
numbers in it are returned; nothing is written back into it, and no command
writes anything into the folder your scans live in. Output goes where you send
it, so a raw scan can be reprocessed any number of times and remains the
record of what the instrument measured.

### 1.1 The workflows and their commands

The workflows are meant to be taken in order, and each is one command.

Before anything else, ask whether the scan is good enough for what you want
from it. `xrdkit check` reads a scan, prints its range, step, counting time
and intensities, and says whether each workflow is within reach, with the
criteria it failed where one is not. A scan that is ample for a phase check is
not good enough for a lattice parameter, and the check says so before you
spend an afternoon on it.

Plotting a pattern with hkl indices on it is what a phase check and most
figures in a paper need. `xrdkit plot` plots one scan and writes its peak
list, and given a start cell it indexes the peaks, writes the indexing and
labels a second figure with hkl. `xrdkit stack` draws several scans one above
the other, which is how a composition series is shown.

Identifying the phases the pattern is made of says whose pattern the rest of
the analysis is describing. `xrdkit phases` searches the Crystallography Open
Database for entries made of exactly the elements you name, fetches their
CIFs, simulates each pattern, weighs it against your peaks, and reports what
explains every peak left over. It needs the `phases` extra and an internet
connection.

Measuring the diffractometer itself comes next, and is done once rather than
once per sample. `xrdkit instrument` takes a scan of a standard of known cell
and gives back the instrument parameter file that every later refinement
reads. Section 3 is about this command.

Lattice parameters and a theoretical density are what a composition series or
a solid solution study needs. `xrdkit lattice` refines the cell of a scan from
a start cell, with the zero point or a specimen displacement, and reports the
volume and the density; `xrdkit density` calculates a theoretical density from
a formula and a cell on its own, and a relative density against a measured
one.

A full Rietveld refinement of coordinates and occupancies is expensive in beam
time and in effort, and is normally worth it only for selected samples.
`xrdkit lebail` extracts the cell of a sample in GSAS-II, and `xrdkit rietveld`
carries that result through the fixed atoms, coordinates and occupancies modes
in turn. Both need GSAS-II.

The usual order is to check every scan, plot it, identify the phases present,
measure the instrument once, refine the cell and the density of each sample,
and only then take the few samples that matter through a Rietveld refinement.

### 1.2 What you need

Python 3.11 or later, and the package, which is installed with
`pip install xrdkit`. Phase identification needs pymatgen, which comes with
the `phases` extra, `pip install "xrdkit[phases]"`, and is worth installing
only when you intend to use `xrdkit phases`.

GSAS-II is needed by three commands, `xrdkit instrument`, `xrdkit lebail` and
`xrdkit rietveld`, and by nothing else. It is a separate Python installation
rather than a dependency of xrdkit: install it from its own home page,
<https://advancedphotonsource.github.io/GSAS-II-tutorials>, and either point
the kit at it with the environment variables `XRDKIT_GSAS2_PYTHON` and
`XRDKIT_GSAS2_HOME` or install it at `~/gsas2main`, which is where the kit
looks by default. Everything else works without it.

Of your own you need the scans of your samples; a scan of a standard measured
on the same instrument with the same optics, which is what Section 3 turns
into an instrument parameter file; and a CIF for each phase you intend to
refine, from the Crystallography Open Database, the ICSD, or the supporting
information of a paper. Keep a note of where each CIF came from.

### 1.3 How the commands and their output are shown

A sentence says what a block is, then the command sits alone in a block of its
own, written as you would type it:

```console
xrdkit --version
```

and everything it printed follows in a block of its own:

```text
xrdkit 0.1.0
```

There is no prompt character to delete, so a command can be copied straight
out. Every command is run from the project folder, the one holding
`xrdkit.toml`, unless the text says otherwise. Where a command prints an
absolute path, the path shown is that of the folder these examples were run
in, and yours will be your own.

A few blocks are not the output of a command but a folder layout or an
extract from a file you edit yourself. The sentence before the block always
says which it is.

### 1.4 The example files behind every number

Every number and every output block in this guide came from running the
command shown, in a folder laid out as Section 2 describes, on these files:

| Path in the examples | What it is |
| --- | --- |
| `data/raw/pellet_a.xrdml` | a sintered pellet of a tetragonal tungsten bronze, the composition Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6 |
| `data/raw/pellet_b.xrdml` | a sintered pellet of the same series at a different composition |
| `data/raw/powder_a.xrdml` | the calcined powder of the first composition, measured with the same range, step and counting time |
| `data/standards/lab6.xrdml` | a scan of NIST SRM 660c lanthanum hexaboride on the same instrument |
| `cifs/lab6.cif` | the structure of that standard |
| `cifs/2100720.cif` | the tungsten bronze of COD entry 2100720, the reference structure of the samples |

None of these ships with the package, and your own numbers will differ from
the ones printed here. What should not differ is the shape of the output: the
same lines, in the same order, with your values in them. Where a number in the
text is worth comparing with your own, the text says what a good value looks
like.

## 2. Setting up a folder

Everything to do with one project lives in one folder, and every command is
run from it. The folder holds the scans as the instrument wrote them, the
CIFs, a project file that names them, and the results the commands write.

### 2.1 Starting a project with xrdkit init

Make an empty folder, change into it, and run the one command that starts a
project.

```console
xrdkit init --name ttb-example
```

It prints the one file it wrote.

```text
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\xrdkit.toml
```

`--name` is the project name recorded in the file, and defaults to the name of
the folder. An existing `xrdkit.toml` is never overwritten, so running the
command twice is safe.

Besides the file, `init` makes three folders if they are missing, and the
layout after it has run is this.

```text
cifs/              one CIF per phase, with a note of its source
data/raw/          scans as the instrument wrote them, never edited
results/           peak lists, indexing tables, refinement projects and exports
xrdkit.toml        the project file
```

Two more folders are worth making by hand, because the commands use them
without creating them in advance. `data/standards/` holds the scan of the
standard and the instrument parameter file made from it, which describe the
diffractometer rather than any sample. `figures/` holds the png and pdf files
the plotting commands write when you give them `--out`.

The package ships no data. Your scans, your CIFs and your project file live in
a folder of your own, laid out as above, and the commands write their results
beside them.

### 2.2 The project file, table by table

`xrdkit.toml` names the things a command would otherwise need on the command
line every time. Every path in it is relative to the folder holding the file,
and the commands find it from that folder or any folder inside it. `init`
writes it with `[project]` filled in and everything else as commented
examples, so setting up a project is a matter of uncommenting and editing.

`[project]` carries the name and the file format version, and is the only
table `init` leaves live.

`[instruments.KEY]` describes a diffractometer. `wavelength` is in angstroms,
either the K alpha 1 and K alpha 2 pair or K alpha 1 alone, and `ka2` says
whether the scans hold the K alpha 2 satellites. `radius`, the goniometer
radius in millimetres, and `instprm`, the instrument parameter file, are
optional, but a Le Bail or Rietveld refinement needs the `instprm` and a
refined specimen displacement needs the `radius`. You do not have to write
this table by hand: Section 3 shows `xrdkit instrument` appending it.

`[structures.KEY]` describes a phase, and its key names the results folders
of anything refined against it. A structure gives `library`, an entry of the
structure library, or `cif`, a CIF file, or both. The entry gives what a
refinement needs to know about the structure type: its space group and crystal
system, its sites and their kinds, the coordinates each site's Wyckoff
position leaves free, the Uiso groups, the site that holds the origin along a
polar axis, the anions and the bond limits. The CIF gives the coordinates. A
Le Bail extraction needs no coordinates, so it runs from a library entry
alone; the Rietveld modes put atoms in, so they need a CIF beside the entry.
`composition` is the formula put on the structure, `cell` its cell parameters,
and `z` the formula units per cell. Nothing in the kit writes a structure
table, so this one is written by hand. The example used throughout this guide
is the tungsten bronze of the samples.

```text
[structures.ttb_p4bm]
library = "ttb/P4bm"
cif = "cifs/2100720.cif"
composition = "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
cell = { a = 12.45, c = 3.94 }
z = 5
exchange = [["Sr", "Ba"]]
atoms = { A1 = { Sr1 = "Sr", La1 = "La" }, A2 = { Sr2 = "Sr", Ba2 = "Ba", La2 = "La" }, B1 = { Nb1 = "Nb", Ti1 = "Ti" }, B2 = { Nb2 = "Nb", Ti2 = "Ti" } }
```

Three of those keys need a word of explanation.

`atoms`, beside a library entry, says which atoms of the CIF sit on which site
of the entry, as a table of the entry's site labels. It is needed in two
cases. The first is an element of the composition that the entry's prototype
does not carry: the tungsten bronze prototype carries neither lanthanum nor
titanium, so the table above places them. An element the CIF lacks is placed
only on the sites the table names it on, beside the atom of its host element
on each, the host being the element present on all those sites, and its amount
is split among them in proportion to the multiplicity and occupancy of the
host. Naming one site alone would put all of it there, and the fit would be
poorer if the real structure has it on both. Leaving such an element unplaced
stops the run with a message naming the element and the entry's sites. The
second case is a CIF whose labels differ from the entry's, since a refinement
matches each site the entry names to the CIF atom of the same label.

`origin` names the site whose coordinate along a polar axis is held, to stop
the whole structure sliding along it. Left out, the entry's own applies:
`ttb/P4bm` is polar along c and holds its B1 site. A site label of the entry
chooses another, and `origin = false` holds none. With no origin held, the
coordinates mode refuses to free a coordinate along the polar axis, naming the
axis, since the refinement would have nothing to fix the structure in place
along it.

`exchange` lists groups of elements whose occupancies the occupancies mode
trades, `[["Sr", "Ba"]]` here. The sites a group is traded between are every
site of one kind that holds any of its elements.

`[samples.KEY]` names a scan and ties it to an instrument and one or more
structures. `file` is the scan, relative to the project folder; `instrument`
is an instruments key; `structures` is a list of structures keys; and `form`
is `powder` or `pellet`, which decides whether a refinement frees the zero
point or a specimen displacement. `stage`, `temperature_c`, `archimedes`, the
measured density in grams per cubic centimetre, and `notes` are optional and
are carried into the results for the record. Section 3.4 shows
`xrdkit add-sample` writing this table for you.

`[refine]` holds what the refinement commands refine with, for every sample,
and a sample's own `[samples.KEY.refine]` table overrides any key of it. Every
key is optional, and a table given in part, such as `background = { terms = 8 }`,
keeps the rest of its values.

| Key | Default | What it sets |
| --- | --- | --- |
| `two_theta` | the scan's range | the range refined, `[low, high]` in degrees, clipped to the scan |
| `background` | `{ function = "chebyschev-1", terms = 6 }` | the background function and its number of terms |
| `max_passes` | `{ lebail = 60, fixed_atoms = 60, coordinates = 100, occupancies = 100 }` | the most passes of a stage, by mode |
| `unsettled` | `accept` for every mode | whether a stage still moving after those passes is kept (`accept`) or rolled back (`reject`) |
| `followed` | none | reflections to follow from mode to mode, as `[h, k, l]` triples |

The range is the one key worth setting per sample, because of what sits at the
bottom of it. Reflections the structure expects at low angle but the pattern
does not show give a Le Bail extraction nothing to fit, so it drags the
background up to cover them and pushes the microstrain negative to narrow
them, after which the size and the microstrain correlate strongly and neither
is worth reading. Starting the range above those reflections removes the
correlation, which is why a sample whose low angle reflections are absent or
very weak is given a range of its own rather than the scan's.

### 2.3 The formats the commands read, and the wavelength rule

The commands read three formats, chosen by the suffix: PANalytical `.xrdml`,
and the two and three column text patterns `.xy` and `.xye`. A sample's `file`
key in the project file goes through the same reader, so a text pattern is
accepted wherever an `.xrdml` is, the Le Bail and Rietveld commands included.

An `.xrdml` carries its own wavelength, its counting time and a sample
identifier, and needs nothing further. A text pattern carries none of those.
The first two columns are read as the two theta and the intensity, and a
third, which an `.xye` carries, as the esd of the intensity. Any number of
header or comment lines may come first: a line counts as one until the data
begins, and the data begins at the first line whose first character could
begin a number, so column headings and comments opened with a hash or a
semicolon are passed over alike. Fields may be separated by one or more
spaces, by tabs, or by a comma with or without spaces around it; intensities
may be written as integers or as decimals; either line ending is read; and
blank lines are ignored wherever they fall. A line of fewer than two numbers,
a field that is not a number once the data has begun, or a two theta that does
not increase on the line before it is refused, with a message naming the file
and the line. Because the counting time is one of the things such a file does
not carry, `xrdkit check` reports the time per step of a `.xy` or `.xye` scan
as unknown rather than as a number.

A text pattern therefore has to be told what radiation measured it. A sample
key takes its wavelength from the instrument table of the project file and
needs nothing further, which is the way to run the Le Bail and Rietveld
workflows on one. A file named directly on the command line needs
`--wavelength ANGSTROM` on any command that uses the wavelength, which is
`xrdkit lattice` always and `xrdkit plot` when it is labelling hkl. `xrdkit
check` and `xrdkit stack` never use the wavelength and need nothing. A command
that needs it and is given neither refuses before it writes anything, with one
line naming the file and saying to give `--wavelength` or a sample key.

Most diffractometer software exports a two column `.xy` directly, and where it
does not, any converter that writes one will do; PowDLL and ConvX are two free
examples. The reader was written to take what such converters produce, which
is why the rules above are as loose as they are. There is no reader for Bruker
`.raw` or `.brml`, or for `.gsas` or `.fxye`; convert one of those to `.xy`
first.

## 3. The instrument file

### 3.1 What the file is, and why a standard is needed

The instrument parameter file holds the wavelengths and their intensity ratio,
the polarisation, the zero point, the Gaussian width terms U, V and W, the
Lorentzian terms X and Y, and the asymmetry SH/L: everything about the peak
shape that belongs to the diffractometer rather than to the sample. It is held
fixed in every sample refinement afterwards.

The reason is that a sample's own size and microstrain broadening and the
instrument's resolution have the same kind of angular dependence, so refining
both on one pattern lets each take up whatever the other leaves. Measuring the
instrument once on a standard that is known to be sharp fixes the instrument
half, and whatever width is left over in a sample pattern is then the
sample's.

That is what a standard is for. A line position and line shape standard, NIST
SRM 660c lanthanum hexaboride or SRM 640 silicon, has a certified lattice
parameter and is ground to contribute no broadening worth the name, so the
widths of its reflections are the instrument's own. Measure it on the same
instrument with the same optics as the samples, because a change of slits or
of the monochromator changes the answer. One scan serves every sample measured
in that configuration, which is why this is done once rather than once per
sample.

The command does the work in two parts. It fits every reflection of the
standard as a K alpha doublet and fits the widths to the Caglioti relation,
which gives the Gaussian terms and needs no GSAS-II. It then refines that
starting point against the standard's own structure in GSAS-II, so that the
Lorentzian terms and the asymmetry are measured rather than guessed, with the
cell held at the certified value and the standard given no broadening of its
own.

### 3.2 Running xrdkit instrument

The standard's scan and CIF are named, the certified cell is given as the one
free parameter of a cubic cell, and the goniometer radius is recorded so that
a later refinement can free a specimen displacement in the right geometry.
The certified lattice parameter of SRM 660c is 0.415 682 6 nm with an expanded
uncertainty of 0.000 008 nm at 22.5 degrees Celsius, from the NIST certificate
of 10 March 2015, which is 4.156826 angstrom.

```console
xrdkit instrument data/standards/lab6.xrdml --cif cifs/lab6.cif --cell 4.156826 --radius 145 --name diffractometer
```

It prints the width fit, the stages it ran, the residuals and every refined
parameter, then the instrument table it appended and the files it wrote.

```text
peaks   14 fitted from 10 to 98 degrees
U = 0.0107, V = -0.0122, W = 0.0085 degrees squared
rms of the width fit 0.00254 degrees
stages  background and scale, zero, U V W, X Y, SH/L
Rwp 7.160 per cent, GOF 1.909
  Zero   -0.0234803  esd 0.000283
  U         21.8437  esd 1.93
  V        -26.2779  esd 2.25
  W         10.8037  esd 0.651
  X         4.56529  esd 0.151
  Y         -1.5982  esd 0.35
  SH/L    0.0220585  esd 0.000396
[instruments.diffractometer]
wavelength = [1.540598, 1.544426]
ka2 = true
radius = 145
instprm = "data/standards/lab6.instprm"
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\data\standards\lab6.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\data\standards\lab6_instrument.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\xrdkit.toml
```

The run took about twenty seconds. Read the first three lines and the
residuals before you use the file. The U, V and W of the width fit are in
degrees squared and are the Caglioti relation as the peak widths alone give
it; the rms is how far the fitted widths sit from that relation, and a few
thousandths of a degree says the relation describes them well. The U, V and W
printed lower down are the refined values in the centidegrees squared GSAS-II
works in, which is why they look nothing like the first three. An Rwp of a few
per cent and a goodness of fit between one and two on a standard is what a
sound file looks like.

A negative Y, as here, means the Lorentzian width the data want is smaller
than the terms can describe together. It is not fatal and the file is usable,
but it is worth refining a second variant with Y held at zero and adopting
that instead when its Rwp is no worse by more than a few tenths of a
percentage point, since one fewer parameter and a fixed rather than a negative
Y is the simpler description of the same widths.

Two files are written, both in `data/standards` when the command is run in a
project. `lab6.instprm` is the file itself, the one GSAS-II exported
after the refinement rather than the starting point, and is what the
`instprm` key of the instrument table points at. The starting file keeps the
same stem with `_start` on the end, so the two are never confused, and the
GSAS-II project and its logs go under `work/` beside them.

`lab6_instrument.csv`, named from the same stem, is the record of how the
file was made, so that a file found months later can be traced back to the
scan it came from. Its first columns say what went in: the scan and the CIF,
the phase name, the wavelength, the six cell parameters held, and the two
theta window the widths were fitted over. Then come the results of the width
fit, the number of reflections it used and its U, V, W and rms. Then every
refined parameter with its esd, in the order the command prints them, the zero
first and the asymmetry last. Then the Rwp and the goodness of fit, the
goniometer radius, and the path of the file written. The last three columns
are the standing record every results file in the kit carries: the method in
one sentence, the date, and the version of xrdkit that wrote the row.

### 3.3 The options

`SCAN` is the standard's scan, a `.xrdml`, `.xy` or `.xye` file or a sample
key of the project file, read by the rule of Section 2.3.

`--cif` is required and names the standard's CIF.

`--cell` takes the certified cell, as the free parameters of its crystal
system in angstroms and degrees: one number for a cubic standard as here, two
for a tetragonal, hexagonal or trigonal one, three for an orthorhombic one,
four with `--system monoclinic`, and six for a triclinic one. `--system`
settles the crystal system where the count of numbers leaves a choice. The
cell is held throughout, which is the point of using a standard.

`--phase` is the phase name GSAS-II gives the standard, and defaults to the
stem of the CIF file, `lab6` in the run above. It matters only in the GSAS-II
project, so the default is usually right.

`--window` is the two theta range the widths are fitted over, ten to
ninety-eight degrees by default. Below the lower limit the reflections of a
standard are few and asymmetric; above the upper one the K alpha doublet is
wide enough that the fit is about the splitting rather than the instrument.
Narrow it if your scan stops short of ninety-eight degrees or if the high
angle reflections are too weak to fit.

`--radius` is the goniometer radius in millimetres, and is written into the
instrument file. A text pattern importer gives GSAS-II no radius, and its
default of two hundred millimetres would make a refined specimen displacement
wrong by the ratio of the two, so give this option if any sample will be
refined with a displacement. It is also recorded in the instrument table.

`--wavelength` overrides the K alpha 1 wavelength, which is otherwise the
scan's own or, for a sample key, the instrument's.

`--name` appends an `[instruments.KEY]` table to `xrdkit.toml` naming the
wavelengths, whether the scans carry K alpha 2, the radius where one was
given, and the path of the file written, relative to the project. The table is
checked with the rest of the file before anything is written, and a key that
is already there is refused so that an existing instrument is never
overwritten.

`--stem` names the output files, and defaults to the stem of the scan, which
is what the run above takes, so the files are named after the standard they
were measured from. Give it where one folder holds the files of more than one
optical configuration and the standard's name alone would not tell them
apart.

`--out` is the folder the two files go in, and defaults to `data/standards`
under the project root when an `xrdkit.toml` is found and to the current
folder when none is.

`--json` prints the files written and the numbers as JSON on standard output,
with the progress lines on standard error so that the JSON stands alone.

Every input and option is checked, and GSAS-II located, before anything is
written, so a mistyped path or a missing GSAS-II costs you the message and
nothing else.

### 3.4 Adding a sample

`xrdkit add-sample` writes a `[samples.KEY]` table for you rather than leaving
you to write one by hand. It comes after Section 3 rather than in Section 2
because the table it writes names an instrument, and unless you give
`--instrument` it expects the project to hold exactly one, which is true only
once `xrdkit instrument --name` has appended it. It also needs a structures
key, so the structure table of Section 2.2 must be in the file first. The
scan added here is a sintered pellet, so `--form pellet` is given rather than
left at the default.

```console
xrdkit add-sample data/raw/pellet_a.xrdml --structure ttb_p4bm --form pellet
```

It prints the table it appended.

```text
[samples.pellet_a]
file = "data/raw/pellet_a.xrdml"
instrument = "diffractometer"
structures = ["ttb_p4bm"]
form = "pellet"
```

`FILE` is the sample's scan and must lie inside the project folder, since the
table stores it relative to the project. `--structure KEY` is required and is
repeated once for each phase the sample holds. `--name` is the sample key and
defaults to the stem of the file, `pellet_a` here. `--instrument` names the
instrument and defaults to the only one when the project has exactly one.
`--form` is `powder` or `pellet` and defaults to `powder`, and is given as
`pellet` above because that is what the scan is. It is worth setting
correctly, because a pellet is the case where a refinement frees a specimen
displacement rather than the zero point: the surface of a pellet sits where
the press left it rather than flush with the holder, and that shift looks
exactly like a cell error at low angle. `--stage` is free text for the
processing stage, such as calcined or sintered. `--temperature-c` and
`--archimedes` record the processing temperature in degrees Celsius and the
measured density in grams per cubic centimetre, and `--notes` is free text.

The rest of the project file is left exactly as it was, byte for byte, and the
whole file is checked with the new table added before anything is written, so
a table that would not load is refused rather than saved. A key that is
already in the file is refused too; give `--name` to add the same scan under
another key.

The other two example samples are added the same way. The second pellet is a
pellet too.

```console
xrdkit add-sample data/raw/pellet_b.xrdml --structure ttb_p4bm --form pellet
```

```text
[samples.pellet_b]
file = "data/raw/pellet_b.xrdml"
instrument = "diffractometer"
structures = ["ttb_p4bm"]
form = "pellet"
```

The calcined powder takes the default form, so `--form` is left out.

```console
xrdkit add-sample data/raw/powder_a.xrdml --structure ttb_p4bm
```

```text
[samples.powder_a]
file = "data/raw/powder_a.xrdml"
instrument = "diffractometer"
structures = ["ttb_p4bm"]
form = "powder"
```

The project now names three samples, and every command from here on can be
given a key rather than a path.

## 4. Checking a scan

What the data must be depends on what you intend to get out of them. A scan
that is ample for a phase check is not good enough for a lattice parameter,
and a scan that gives a good lattice parameter is still not good enough for
refined occupancies. `xrdkit check` measures a scan against the criteria of
each workflow and says which of them it is within reach of, which is worth
knowing before you spend an afternoon on a refinement the counting statistics
will not support.

### 4.1 Running the command

Given a sample key, the command reads the scan the project file names for it.

```console
xrdkit check pellet_a
```

It prints the sample and its composition, the numbers it measured, and a
verdict for each workflow, then the file it wrote.

```text
sample  pellet_a, Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6
range   10.01 to 99.98 degrees
step    0.0217 degrees
points  4141
time    34.2 s per step
maximum 13964 counts
median  998 counts
peak over median 14
high angle maximum 1833 counts, from 69.99 degrees
high angle median  1012 counts

plotting: suitable
phase_identification: not suitable (peak over median 14.0, below 20, weak phases may not be visible)
le_bail: not suitable (range 10.01 to 99.98 degrees, short of 10 to 120; high angle peak over median 1.8, below 10)
rietveld: not suitable (range 10.01 to 99.98 degrees, short of 5 to 130; step 0.0217 degrees, outside 0.01 to 0.02; strongest peak 13964 counts, below 20000)
```

The numbers are the measured range and step, the number of points, the
counting time per step, and the intensities. The median stands in for the
background, since most of the points in a powder pattern are background, and
the peak over median is the strongest peak divided by it. The high angle
numbers are the same two over the last third of the scanned range, which for
this scan begins at 69.99 degrees, because a cell is refined on the high angle
reflections and it is their contrast with the background that decides whether
that can be done.

A workflow that is out of reach is followed by every criterion it failed, each
with the value measured and the threshold it fell short of, so there is no
guessing about which number to improve.

For a sample of the project file the report is also kept, under
`results/check/KEY`, so that the state of a scan at the time it was analysed
stays on the record.

```text
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\check\pellet_a\check_pellet_a.txt
```

The same command takes a path instead of a key, which is the way to look at a
scan that is not in the project at all.

```console
xrdkit check data/raw/powder_a.xrdml
```

There is no sample line, because a bare file has no key and no composition,
and nothing is written: the report goes to the terminal alone.

```text
range   10.01 to 99.98 degrees
step    0.0217 degrees
points  4141
time    34.2 s per step
maximum 10032 counts
median  1079 counts
peak over median 9
high angle maximum 1505 counts, from 69.99 degrees
high angle median  1008 counts

plotting: suitable
phase_identification: not suitable (peak over median 9.3, below 20, weak phases may not be visible)
le_bail: not suitable (range 10.01 to 99.98 degrees, short of 10 to 120; high angle peak over median 1.5, below 10)
rietveld: not suitable (range 10.01 to 99.98 degrees, short of 5 to 130; step 0.0217 degrees, outside 0.01 to 0.02; strongest peak 10032 counts, below 20000)
```

### 4.2 The criteria, in numbers

These are the criteria the command applies. Plotting and phase identification
share a column, because they ask nearly the same of a scan; where they differ,
the counting statistics row says so.

| | Plotting and phase check | Lattice parameters and density, by Le Bail | Rietveld refinement |
| --- | --- | --- | --- |
| Purpose | Show what the sample is, identify the phases, label the reflections | Refine a and c, and the cell volume, well enough for a theoretical density | Refine the structure: coordinates, occupancies, displacement parameters |
| Angular range | 10 to 80 degrees two theta, extended to 90 for a figure meant for publication | 10 to 120 degrees. The reflections above 90 degrees are what pin the cell down, because a given error in d shifts them furthest in two theta | 5 to 130 degrees, or wider where the instrument allows |
| Step size | 0.01 to 0.03 degrees | 0.013 to 0.026 degrees | 0.01 to 0.02 degrees |
| Counting statistics | Strongest peak above about 2000 counts, which is all plotting asks of a scan. The phase check adds a background of at least 100 counts per step and a strongest peak at least 20 times it. For a secondary phase at 1 to 2 weight per cent to be visible on a square root or logarithmic scale, the background itself needs at least 100 counts per step, which means three to five times the count time of a quick scan | Strongest peak above 10000 counts, and the high angle peaks at least 10 times the background, since those are the ones the cell is refined on | Strongest peak above 20000 counts, and background above 200 counts per step at high angle. On a laboratory instrument this normally means several hours per scan |
| Sample preparation | Flat sample, flush with the surface of the holder | Crushed pellet powder preferred. A pellet surface is acceptable if a displacement term is refined | Powder crushed and sieved below about 45 micrometres, then back loaded, side loaded or mounted in a capillary, to limit preferred orientation |
| Standards needed | None | Either an internal standard, silicon SRM 640 at 10 to 20 weight per cent mixed into the powder, or a refined displacement term. Also an instrument parameter file from a LaB6 (SRM 660) scan on the same instrument and the same optics | An instrument parameter file from a LaB6 scan on the same configuration, a CIF of the expected structure, and a known composition |
| Radiation | K alpha 2 present is acceptable | K alpha 2 present is acceptable, since the fit follows the K alpha 1 positions | Monochromatic Cu K alpha 1 strongly preferred |
| What goes wrong otherwise | Weak phases hide in the noise and are missed, so the sample is reported as single phase when it is not. Peaks are too noisy for the peak finder to separate, and labels go on the wrong reflections | The cell is refined on low angle reflections alone, where a displacement error and a cell error look alike, so a and c come out precise and wrong, and the density with them | The refinement runs, but the parameters are not determined by the data. Occupancies and coordinates drift to whatever fits the noise, and the esds do not say so unless they are examined |

Only the angular range, the step size and the counting statistics are checked
by the command, because only those can be read off a scan. The other rows are
for you: nothing in a data file says how the powder was mounted or whether a
standard was measured, and those decide as much as the numbers do.

Two further numbers govern the first three workflows. For hkl labelling, peak
positions good to about 0.02 degrees are enough, and a zero error of up to
0.05 degrees does not change which reflection a peak is assigned to. That same
zero error does change the lattice parameters, which is why the lattice
workflow either uses an internal standard or refines a displacement or zero
term: an uncorrected specimen displacement of 0.1 mm shifts a by about 0.001
angstrom. For a useful density, a and c are needed to a relative precision of
0.01 per cent, which is about 0.001 angstrom on a 12 angstrom axis, because
the density goes as the reciprocal of the cell volume and the relative error
in the density is therefore about three times the relative error in a lattice
parameter.

### 4.3 What the verdicts said about the example scans

Read against the table, `pellet_a` is comfortable for plotting: it runs from
10.01 to 99.98 degrees at 0.0217 degrees per step, both inside the plotting
row, and its strongest peak of 13964 counts is well above the 2000 that row
asks for. Its peak over median of 14 is short of the 20 the phase check wants,
so a secondary phase at a per cent or two could be sitting in the background
unseen, and the phase check is reported as out of reach for that reason alone.

For a Le Bail fit it fails twice over, and the second failure is the
instructive one. The range stops at 100 degrees rather than 120, which loses
the reflections that pin the cell down. More telling, the high angle peak over
median is 1.8 against the 10 asked for: above 70 degrees the strongest peak
reaches 1833 counts on a background of 1012, so the reflections the cell would
be refined on are barely clear of the noise. A longer count at high angle is
what that number is asking for, not a wider range alone.

The Rietveld verdict adds a step size outside 0.01 to 0.02 degrees and a
strongest peak of 13964 against the 20000 wanted. Taken together the three
verdicts say what these scans are: quick survey scans, good for seeing what
the sample is and for labelling its reflections, and not counted long enough
for a refined cell or a refined structure.

`powder_a` tells the same story a little more sharply. Its strongest peak is
10032 counts against the pellet's 13964, and its peak over median is 9.3, so
it is the weaker of the two for a phase check even though it is the calcined
powder rather than a pellet surface. The two scans were measured with the same
range, step and counting time, so the difference is in the samples and not in
the measurement.

### 4.4 The options

`SCAN` is the only argument: a path to a `.xrdml`, `.xy` or `.xye` file, or a
sample key of the project file. A key brings the composition into the heading
and sends the report to `results/check/KEY`; a path prints the report and
writes nothing.

`--json` prints the same numbers and verdicts as JSON instead of the report,
for a script that wants to read them rather than a person.

The command never uses the wavelength, so a `.xy` or `.xye` scan needs no
`--wavelength` here. What such a scan does not carry is the counting time, and
the report says so rather than inventing a number: the time line reads
`unknown` instead of a value in seconds. Every other line is measured from the
two theta and intensity columns and reads the same as for an `.xrdml`.

## 5. Plotting and indexing

`xrdkit plot` draws one pattern and writes its peak list, and where it can
index the peaks it writes the indexing too and labels a second figure with
hkl. `xrdkit stack` draws several patterns one above another, which is how a
composition series is shown. Between them they cover the first thing anyone
does with a new scan.

### 5.1 Plotting a sample

Given a sample key, the command takes the wavelength from the instrument
table and the start cell, crystal system and space group from the sample's
first structure, so there is nothing to type but the key.

```console
xrdkit plot pellet_a
```

It prints the refined cell and how the indexing went, then every file it
wrote.

```text
a = 12.4803, c = 3.9324 angstrom, zero 0.170 degrees, 45 of 49 peaks indexed, rms 0.0109 degrees
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\peaks_pellet_a.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\pattern_pellet_a.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\pattern_pellet_a.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\indexed_pellet_a.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\pattern_hkl_pellet_a.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\plot\pellet_a\pattern_hkl_pellet_a.pdf
```

Six files: the peak list, the pattern as a png and a pdf, the indexing, and
the hkl labelled pattern as a png and a pdf. The peak list gives each
reflection its position, intensity, prominence, full width at half maximum, d
spacing and relative intensity as a percentage of the strongest peak found.
The indexing file gives each peak the reflection assigned to it, how far the
two differ, and whether more than one reflection was within reach. Figures are
written as a png to look at and a pdf to put in a document, at 300 dots per
inch.

The zero offset of 0.170 degrees on this scan is worth noticing. It is a
pellet, and a pellet surface that sits a little proud of the holder shifts
every peak in the same direction, which is exactly what an apparent zero
offset looks like. It does not change which reflection a peak is assigned to,
which is why the labels are sound, but it would change a lattice parameter,
which is what the lattice workflow deals with.

### 5.2 Plotting a scan that is not in the project

Without a sample key there is no structure to take a cell from, so the cell is
given on the command line. The reflection conditions of the space group are
worth giving too, to keep reflections the group forbids out of the labels.

```console
xrdkit plot data/raw/powder_a.xrdml --cell 12.45 3.94 --space-group P4bm
```

The same line is printed, and the files go to `results` and `figures` under
the folder the command was run in rather than to a sample's own folder.

```text
a = 12.4761, c = 3.9323 angstrom, zero -0.030 degrees, 48 of 49 peaks indexed, rms 0.0159 degrees
results\peaks_powder_a.csv
figures\pattern_powder_a.png
figures\pattern_powder_a.pdf
results\indexed_powder_a.csv
figures\pattern_hkl_powder_a.png
figures\pattern_hkl_powder_a.pdf
```

This is the same material as `pellet_a` in powder rather than pellet form, and
the zero offset of minus 0.030 degrees against the pellet's plus 0.170 is the
difference between a flat powder bed and a pellet surface. The cell comes out
within a few thousandths of an angstrom of the pellet's, as it should.

### 5.3 What the indexing does

The cell is given as the free parameters of its crystal system, and the count
of numbers chooses the system: one number is cubic, two are tetragonal, three
are orthorhombic, four are monoclinic and six are triclinic. Two numbers are
ambiguous, since tetragonal, hexagonal and trigonal all take a and c, so
`--system` settles it; tetragonal is taken when nothing is said. The cell
above, 12.45 and 3.94 angstrom, is the tetragonal tungsten bronze.

The indexing assigns a reflection of the cell to each peak and refines the
cell as it goes. It works in cycles: the first indexes only the low angle
peaks with a loose tolerance, where a cell that is still some way off can be
trusted to put reflections near the right peaks, and each later cycle indexes
the whole list against the cell the cycle before gave. It searches for the
zero offset first, unless `--zero` gives one.

Reflection conditions are known for eight space groups: Pm-3m, P4mm, P4bm,
P4/mbm, R3c, R3m, Pbnm and Amm2. Given one of those, the reflections the group
forbids are never offered to a peak. Given any other symbol the command says
so in a note and indexes without conditions, which means a peak may be
labelled with a reflection its group does not allow, so the labels are
provisional until the group is checked. Given no space group at all, no
conditions are applied, which is the same situation without the note.

A cell refined this way is good enough to label reflections with. It is not a
lattice parameter for publication, which is what `xrdkit lattice` is for: the
fit here is a linear least squares on the indexed peaks alone, and the zero
offset it reports is the one that indexed the most peaks rather than one
refined alongside the cell.

Peaks that look like K alpha 2 satellites are flagged rather than dropped, and
are removed before the indexing, because a satellite sits at the position its
parent's d spacing gives at the longer wavelength and indexing it against the
cell is meaningless. Below about 50 degrees the doublet is not resolved and
there is nothing to flag. `--no-satellites` drops them from the peak list as
well, for a list meant to be read rather than to be indexed.

The count of peaks indexed is the number to look at. Of the 49 peaks found in
`pellet_a`, 45 were assigned a reflection and four were not. A peak that
matched no reflection is exactly the peak a reader should be looking at, and
there are four common explanations. A secondary phase gives peaks that do not
move with composition across a series and that keep a fixed ratio to one
another; naming that phase is what `xrdkit phases` is for. K beta gives a
faint copy of a strong reflection at lower angle, where the filter or the
monochromator is imperfect, at the position the parent's d spacing gives at
about 1.392 angstrom for copper. Tungsten L lines come from a contaminated or
an aged tube and appear at fixed angles that do not move with the sample, so
the same stray peaks turn up in every pattern measured on that instrument. A
surviving K alpha 2 satellite sits just above its parent and carries roughly
half the intensity, and is a sign that the flagging did not suit the pattern
rather than a sign of anything wrong with the sample.

### 5.4 The hkl labels

Each label rides on top of its own peak, a couple of points above the trace,
which keeps the label and the reflection it names together however the pattern
rises and falls. The top of the axes is raised above the tallest point of the
trace before the labels are placed, so that the label on the strongest peak
stays inside the figure.

Labels are placed strongest first and each is tested against the ones already
placed, so a weaker peak whose label would collide with one already there
keeps no label at all. Labels are never stacked or shifted sideways onto a
neighbour: in a crowded stretch it is the weak reflections that lose theirs,
which is the behaviour that keeps every label that is drawn pointing at the
peak it belongs to. Peaks below a few per cent of the strongest are skipped
from the start, for the same reason.

### 5.5 Stacking several scans

A stack takes two or more scans, files and sample keys mixed as you like, and
one label for each in the same order.

```console
xrdkit stack pellet_a pellet_b --labels "x = 0.10" "x = 0.12"
```

It writes one figure and prints it.

```text
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\stack\pellet_a_pellet_b\stack_pellet_a_pellet_b.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\stack\pellet_a_pellet_b\stack_pellet_a_pellet_b.pdf
```

The traces are drawn bottom to top in the order given, each labelled at its
upper right, and are normalised by default so that scans counted for different
times can be compared. The spacing is 1.2 times the tallest scaled trace,
which keeps the tallest peak of one pattern clear of the pattern above it.

### 5.6 The options

`--scale` takes `linear`, `sqrt` or `log` and belongs to both commands. A
linear scale shows the strong reflections in proportion, which suits a figure
about the main phase. A square root scale brings up the weak ones, which is
what a phase check wants and what a stack is usually drawn on. Normalisation,
where it applies, happens after the transform, so a normalised square root
trace runs from 0 to 1 on the square root scale.

`--label` gives `xrdkit plot` the legend label for its trace, and defaults to
the scan's sample identifier, or to a sample's key and composition.
`--labels` gives `xrdkit stack` one label per scan in the same order as the
scans, and is required: a stack with unlabelled traces says nothing.

`--offset` sets the vertical spacing of a stack in place of the default of 1.2
times the tallest scaled trace. Raise it where labels crowd the trace above,
lower it to fit more traces in one figure. `--no-normalise` keeps the measured
intensities instead of scaling each trace to 1, which is what to use when the
point of the figure is that one sample counted higher than another.

`--cell` and `--system` are the start cell and its crystal system, as Section
5.3 describes. `--space-group` names the group whose reflection conditions
apply, and defaults to the structure's for a sample and to none for a file.
`--zero` gives the zero offset in degrees rather than letting the indexing
search for one.

`--no-satellites` drops the peaks flagged as K alpha 2 satellites from the
peak list and the indexing. `--wavelength` overrides the K alpha 1 wavelength
used for the d spacings and the indexing, and is needed for a `.xy` or `.xye`
file that is being indexed, since such a file carries no wavelength of its
own.

`--stem` names the output files, and defaults to the scan's file stem, or for
a stack to the scan stems joined by an underscore. `--out` is the output root,
under which `results` and `figures` are created as needed; without it a sample
writes to `results/COMMAND/KEY` under the project root and a bare file writes
to the current folder. `--json` prints the files written and the results as
JSON.
