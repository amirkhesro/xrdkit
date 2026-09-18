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
atoms = { A1 = { Sr1 = "Sr", La1 = "La" }, A2 = { Ba2 = "Ba", Sr2 = "Sr" }, B1 = { Nb1 = "Nb", Ti1 = "Ti" }, B2 = { Nb2 = "Nb", Ti2 = "Ti" } }
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
very weak is given a range of its own rather than the scan's. The example
project gives `powder_a` exactly that, an extract from `xrdkit.toml` written by
hand beneath its sample table, and Section 8.4 shows what happens without it.

```text
[samples.powder_a.refine]
two_theta = [17.0, 99.98]
```

Two keys of the structure table above earn their place only once a refinement
is run, and are explained where they are used. `cif` beside `library` supplies
the coordinates that the Rietveld modes put atoms at, which a Le Bail
extraction does not need and a Rietveld refinement cannot do without. `z`, the
formula units per cell, is what turns a cell into a density in Section 7 and
what the nominal composition is spread over when the Rietveld modes set up the
site occupancies.

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

## 6. Identifying phases

`xrdkit phases` searches the Crystallography Open Database for structures made
of the elements your sample is made of, simulates the pattern of each, weighs
it against the peaks of your scan, and ranks them. It is what says whose
pattern the rest of the analysis is describing.

The command needs the `phases` extra, `pip install "xrdkit[phases]"`, because
the simulation is done with pymatgen, and it needs an internet connection,
because the database is fetched rather than shipped. Both are checked before
anything is written.

### 6.1 When it is needed

Three occasions, and none of them is optional.

The first scan of a new composition. Until the phases are known, a cell
refinement is fitting a model to a pattern that may not be that model's, and
the numbers it gives will look no different for being wrong. Identify the
phases once per composition, not once per sample: the next batch of the same
nominal composition, made the same way, is a phase purity check against the
first rather than an identification from nothing.

Any peak the main phase leaves unindexed. Section 5 reports these as the
difference between the peaks found and the peaks indexed, and Section 5.3 gives
the four things they are usually made of. Two of those, K beta and a tungsten L
line, belong to the instrument and are settled by where they sit. The other
two, a secondary phase and a mis-set cell, are settled here.

Every phase purity scan of a sintered pellet. Sintering is where second phases
appear that the calcined powder did not have, from a reaction with the
crucible or the setter, from a volatile component leaving, or from a melt at a
grain boundary. A pellet that fired well and looks right can still carry two
weight per cent of something else, which is enough to change a dielectric
measurement and not enough to be obvious in the pattern unless it is looked
for.

### 6.2 Running the command

Given a sample key and nothing else, the elements searched on are those of the
compositions of the sample's structures, which for `pellet_a` is the whole
doped formula.

```console
xrdkit phases pellet_a
```

That finds nothing at all.

```text
35 peaks observed from 10 to 80 degrees
0 COD entries made of exactly Sr, Ba, La, Nb, Ti, O
```

The command then stops, with one line on the error stream.

```text
xrdkit phases: the COD holds no entry of exactly those elements; widen --elements or drop --space-group
```

This is worth dwelling on rather than working around, because it is the normal
result for a doped sample. The search asks for entries made of exactly those
six elements, and nobody has deposited a structure of this particular
substitution. What you are trying to identify is not the doped composition but
the structure type it is built on, and the elements of the parent are what to
search on. Lanthanum and titanium are the dopants, so they come out, and the
search is run on the four elements of the parent tungsten bronze.

The other thing the pattern needs is its zero offset. A simulated pattern is
calculated at true angles and a measured one is not, and `pellet_a` was shown
in Section 5 to carry an offset of 0.170 degrees, which is more than the
matching tolerance and so would push good candidates out of reach. The offset
is given with `--zero` and is subtracted from every observed position before
anything is compared.

```console
xrdkit phases pellet_a --elements Sr Ba Nb O --zero 0.170
```

The run prints the peaks it found, the entries the search returned, where the
index went, the ranking, and every file it wrote. It took about fifteen
seconds, most of it fetching eleven CIFs.

```text
35 peaks observed from 10 to 80 degrees
11 COD entries made of exactly Sr, Ba, Nb, O
  1520953  Ba0.247 Nb2 O6 Sr0.744 P 4 b m
  1537507  Ba3 Nb2 O9 Sr          P 63/m m c
  2020161  Ba3 Nb2 O9 Sr          P 63/m
  2100719  Ba0.67 Nb2 O6 Sr0.33   P 4 b m
  2100720  Ba0.52 Nb2 O6 Sr0.48   P 4 b m
  2100721  Ba0.39 Nb2 O6 Sr0.61   P 4 b m
  2100722  Ba0.14 Nb2 O6 Sr0.86   P 4 b m
  2103856  Ba0.39 Nb2 O6 Sr0.61   X4bm
  2238610  Ba0.4 Nb2 O6 Sr0.6     P 4 b m
  2311739  Ba0.476 Nb2 O6 Sr0.524 P 4 b m (a,b,2*c)
  2311740  Ba0.47 Nb2 O6 Sr0.53   P 4 b m (a,b,2*c)
index written to C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\index.csv
  2100721  explained 34/35  missing  0  score  34
  2311739  explained 33/35  missing  2  score  31
  2311740  explained 33/35  missing  2  score  31
  1520953  explained 27/35  missing  3  score  24
  2100720  explained 27/35  missing  5  score  22
  2238610  explained 26/35  missing  7  score  19
  2100719  explained 26/35  missing  9  score  17
  2103856  explained 35/35  missing 19  score  16  rejected, strongest line absent
  1537507  explained 14/35  missing  6  score   8
  2020161  explained 16/35  missing  8  score   8
  2100722  explained 16/35  missing 19  score  -3  rejected, strongest line absent
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\1520953.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\1537507.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2020161.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2100719.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2100720.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2100721.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2100722.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2103856.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2238610.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2311739.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\2311740.cif
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\cifs\cod\index.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\phases\pellet_a\phases_pellet_a.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\phases\pellet_a\phases_pellet_a_unexplained.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\phases\pellet_a\phases_pellet_a_record.csv
```

That block is the run of 18 September 2026. The database is live and is added
to continually, so a rerun may return an entry that was not there before, and
the ranking may shift when it does. This is a feature rather than a nuisance,
but it means the numbers here are a record of one run rather than a constant,
and it is why every run writes a record of its own.

The CIFs are kept under `cifs/cod`, one file per entry named by its COD id,
and are not fetched again on a later run. Beside them `index.csv` is the
register of reference CIFs, with the source, the identifier, the formula, the
space group, the cell and the citation of each, and it is updated rather than
rewritten, so notes you add to a row by hand survive a fresh download. Under
`results/phases/pellet_a` go three files: the candidates with their counts and
scores, the peaks no candidate explained, and a record of the run carrying the
scan, the elements, the zero, the window, the tolerance, how many entries were
searched and fetched, the main phase, the method, the date and the version of
xrdkit.

### 6.3 Reading the ranking

Read the table in this order.

Take the rejections first, because they are the only column that is a verdict
rather than a measure. A candidate whose strongest reflection falls where
nothing at all was observed is not present in the sample, whatever it scores:
the strongest line of a phase is the last one to disappear, so if it is not
there the phase is not there. That rule rejects 2100722 here, whose cell at
Ba 0.14 is too far from the sample's for its lines to land anywhere near the
right places, and 2103856, which is the same tungsten bronze but in the
nonstandard X4bm centring, which the simulator expands into a structure with
nineteen strong lines the sample does not have. Note that 2103856 explains all
thirty-five peaks and is still rejected, which is the rule doing its job: a
structure that predicts everything, including a great deal that is not there,
has explained nothing. A rejection is a statement about the CIF as it was
read, not only about the sample, and an entry rejected for a reason like that
one is worth looking at rather than deleting.

Then read the missing count. It is the column that discriminates, because a
strong line predicted where nothing was seen is hard evidence against a phase,
while an explained peak may be a coincidence in a crowded pattern. The two
hexagonal perovskites, 1537507 and 2020161, explain fewer than half the peaks
and miss six and eight strong lines; the tungsten bronzes near the sample's
composition miss none, two or three. That gap is the identification: the
sample is the tetragonal tungsten bronze, and it is not the hexagonal
perovskite.

Read the explained count and the score last, and do not read them as a ranking
of composition. 2100721 scores 34 and 2100720 scores 22, and both are the same
phase; what separates them is how close each entry's cell is to the sample's,
which this command is not measuring and `xrdkit lattice` is. The entry to
carry forward is the one whose composition is nearest the nominal one, which
for Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6 is 2100720 at Ba 0.52 to Sr 0.48, and that is
why the project file of Section 2.2 names 2100720 rather than the entry at the
top of this table. Its cell is replaced by a refined one in the later
workflows in any case.

Where two candidates score the same, as 2311739 and 2311740 do here, they are
ordered by COD id, which is arbitrary and is meant to be: the command does not
pretend to separate what it cannot.

### 6.4 The peaks nothing explained

The last file is the one to open when the ranking has been read. It lists every
observed peak that no candidate accounted for, each with its two theta and its
d spacing, which is the form a spacing is quoted in and searched on, and beside
it the reflection of the main phase that falls within the tolerance of it, or
the word unidentified where none does.

For this run the file holds its header and nothing else: every one of the
thirty-five peaks was explained by at least one candidate, which is what a
single phase sample should give.

`--main` names the phase those peaks are attributed to, and defaults to the
top ranked candidate that was not rejected. Give it the COD id you have
decided on, which here would be `--main 2100720`, so that the attribution is
made against the structure you intend to carry forward rather than against
whichever entry happened to score highest. It may also name an entry outside
the candidate set, which is fetched like any other, for the case where you
suspect a particular impurity and want the leftover peaks checked against it.

A weak reflection of the phase you already have is the commonest explanation
of a peak the indexing passed over, and it is worth ruling out before anything
else. Where a peak does come back as unidentified, the next step is to widen
the search: add the elements a crucible or an incompletely decomposed
carbonate could contribute and run again. What is not the next step is
deciding that one peak does not matter.

### 6.5 Where a licensed database fits

xrdkit identifies phases from the Crystallography Open Database, which is free
and requires no licence. Where you have access to a licensed search and match
system with the PDF database behind it, run that search first, because it is
the authoritative identification. It is a better index than any free database:
the PDF is larger and covers the inorganic literature far more completely, it
is curated and quality graded, its cards carry measured rather than calculated
intensities for many phases, and the matching is tuned for it. A card number
with its quality mark and citation is also the conventional record of phase
identification in publications. Nothing in xrdkit replaces that search, and
nothing in xrdkit reads its output.

The two are complementary. The command here serves routine screening on your
own computer, and the licensed route, where it is available, serves the
definitive check and the citable record.

Record what the licensed search found, for the main phase and for every
secondary phase it named: the card number you chose, its quality mark and
citation, the formula, the space group, the cell, and the release of the
database the card came from. That last one matters because cards are revised
and withdrawn, and a card number alone does not say which version was seen.

A card is a pattern, and the later workflows need a structure. A CIF exported
from the licensed database seeds a Le Bail extraction or a Rietveld refinement
in exactly the same way as one fetched from the COD: put it in `cifs` and name
it in a structure table. Where the database exports no CIF, come back to this
command for one.

### 6.6 What the data quality means for this answer

Section 4 reported `pellet_a` as not suitable for phase identification, on a
peak over median of 14 against the 20 the criteria ask for. That verdict
stands, and it should be read against the result above rather than forgotten
because the result looks tidy.

What the ranking here can be trusted to say is which structure type the main
phase is. That conclusion rests on strong reflections and on a gap between the
tungsten bronzes and the hexagonal perovskites that is far larger than the
noise. What it cannot be trusted to say is that there is no second phase. A
phase at one or two weight per cent shows itself in the background, and on this
scan the background is not counted deeply enough for such a phase to rise
clearly out of it; the empty unexplained list is therefore weak evidence of
purity rather than strong evidence. If phase purity is the question being
asked, the scan to answer it with is a longer one, counted three to five times
as long, and the verdict in Section 4 is the command telling you so in advance.

### 6.7 The options

`SCAN` is a path to a `.xrdml`, `.xy` or `.xye` file, or a sample key.

`--elements` gives the element symbols the entries are to be made of, at most
eight, which is what the database's search form takes. For a sample key it
defaults to the elements of the compositions of the sample's structures, and
as Section 6.2 shows that default is the right starting point only for an
undoped composition. For a scan named as a path there is no structure to take
elements from, so `--elements` is required and the command refuses without it.
The search is for entries made of exactly those elements and no others.

`--space-group` restricts the search to one space group symbol, and by default
no restriction is applied. It is worth using when the structure type is already
known and the list is long, and worth leaving off on a first look, since a
symbol given in a nonstandard setting will exclude the very entries you want.

`--zero` is the zero offset in degrees, subtracted from every observed peak
before anything is compared, and defaults to 0. Take it from the plot or
lattice run on the same scan. On a scan that has never been indexed, leave it
at 0 and raise `--tolerance` to about 0.3 degrees instead, so that an
uncorrected shift does not throw the matching out.

`--window` is the two theta range the peaks are taken and the patterns
simulated over, ten to eighty degrees by default. Compare only over the range
that was measured: a reflection simulated outside the scan would be counted as
missing and would tell against a phase that is present.

`--tolerance` is how far an observed peak may lie from a simulated reflection
and still count as explained, 0.15 degrees by default.

`--max-candidates` weighs only the first N entries the search returned, in
COD id order, which is a way to keep a first look quick when the search comes
back with dozens of entries. By default every entry found is fetched and
weighed.

`--main` names the phase the peaks left over are attributed to, as Section 6.4
describes.

`--wavelength` overrides the K alpha 1 wavelength, which is otherwise the
scan's own or the sample's instrument's, and is needed for a `.xy` or `.xye`
file named as a path.

`--stem` names the output files and defaults to the sample key or the scan's
file stem. `--out` is the output root, under which the CIFs and the results
are placed. `--json` prints the candidates, the peaks left over and the files
written as JSON on standard output, with the progress lines on standard error
so that the JSON stands alone.

## 7. Lattice parameters and density

`xrdkit lattice` refines the cell of a scan and, given a composition and the
formula units per cell, turns it into a theoretical density. `xrdkit density`
does the density half on its own, from a cell you already have. This is the
workflow a composition series or a solid solution study rests on, and it is
also the one where a careless answer looks exactly like a careful one, so the
section spends as much space on judging the result as on getting it.

### 7.1 Refining a pellet

Given a sample key, the command takes the start cell, crystal system and space
group from the sample's first structure, the wavelength and the goniometer
radius from its instrument, the formula and the formula units from the
structure's composition and `z`, and the form from the sample table. For
`pellet_a` the form is `pellet`, so the specimen displacement is refined and
the zero is held.

```console
xrdkit lattice pellet_a
```

```text
sample  pellet_a, Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6
peaks   49 found, 43 refitted, 1 fits rejected, 5 satellites excluded, 1 recovered; 43 of 44 indexed (97.7 per cent), 33 used in the refinement
coarse window to 35.00 degrees
tetragonal cell a = 12.4740 +/- 0.0011, c = 3.9295 +/- 0.0004 angstrom
V = 611.426 +/- 0.136 cubic angstrom
zero 0.0000 degrees (held)
displacement -0.2033 +/- 0.0069 mm, radius 145 mm
rms 0.0128 degrees; a low indexed fraction or a high rms means the cell should not be trusted
M = 394.905 g/mol per formula unit, Z = 5
theoretical density 5.3625 +/- 0.0012 g/cm3
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\pellet_a\peaks_pellet_a.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\pellet_a\lattice_pellet_a.csv
```

### 7.2 The five steps behind that line of counts

The peaks line is the whole method in one row, and it is worth reading slowly,
because every later judgement rests on it.

First the peaks are found, and each is examined for the signature of a K alpha
2 satellite: a peak sitting where its neighbour's d spacing would put the
longer wavelength, at no more than a set fraction of that neighbour's height.
Forty-nine peaks were found here and five were flagged.

Second, every peak that is not flagged is refitted, as the K alpha 1 line of a
doublet rather than as a single line, so that its position is the K alpha 1
position rather than the centre of a blend. A fit that does not converge, or
that lands somewhere the peak is not, is rejected and the peak is dropped:
forty-three were refitted and one fit was rejected.

Third, those positions are indexed against the start cell and the cell is
refined by least squares, with the zero or the displacement.

Fourth comes the recovery pass. A flagged peak takes no part in any of that,
because a satellite indexed against the cell is meaningless. Some flagged
peaks, though, are not satellites at all: they are genuine reflections that
happened to sit where a satellite would. The refined cell is used to test
each one, and a flagged peak that falls within the indexing tolerance of a
real reflection, and whose doublet refit is accepted, is brought back in. One
peak was recovered here.

Fifth, with the recovered peaks added, the indexing and the refinement are run
once more, and that second result is what is reported.

The counts relate as follows: the peaks found split into those refitted,
those whose fit was rejected and those excluded as satellites, and the indexed
fraction is quoted against the peaks that were never flagged, forty-four here
and forty-nine for a scan with nothing flagged, of which a subset carries
enough weight to be used in the refinement.

A recovered peak is worth treating with some suspicion. To have been flagged
at all it had to sit under the ceiling on its height over its parent's, which
means a good part of its intensity really does come from the parent's K alpha
2 line, and the refit models it as a lone doublet all the same. Its fitted
position therefore carries a bias. The count is printed so that a reader can
see how much of a refinement rests on such peaks: one out of thirty-three is
nothing to worry about, and ten would be worth a second look with
`--no-satellites`, which excludes every flagged peak and recovers none.

The coarse window is the two theta limit of the first indexing cycle, chosen
from the data rather than fixed. The first cycle works on the low angle peaks
alone with a loose tolerance, where a cell that is still some way off can be
trusted to put reflections near the right peaks; later cycles use the whole
list against the cell the cycle before gave.

### 7.3 Refining a powder, and why the form matters

`powder_a` is the calcined powder of the same composition. Its form is
`powder`, so the zero is refined and no displacement is fitted.

```console
xrdkit lattice powder_a
```

```text
sample  powder_a, Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6
peaks   49 found, 46 refitted, 3 fits rejected, 0 satellites excluded, 0 recovered; 46 of 49 indexed (93.9 per cent), 35 used in the refinement
coarse window to 35.00 degrees
tetragonal cell a = 12.4728 +/- 0.0016, c = 3.9312 +/- 0.0006 angstrom
V = 611.581 +/- 0.218 cubic angstrom
zero -0.0388 +/- 0.0079 degrees
displacement not refined
rms 0.0185 degrees; a low indexed fraction or a high rms means the cell should not be trusted
M = 394.905 g/mol per formula unit, Z = 5
theoretical density 5.3611 +/- 0.0019 g/cm3
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\powder_a\peaks_powder_a.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\powder_a\lattice_powder_a.csv
```

The powder refines a zero of minus 0.039 degrees, which is a real instrument
zero and is small. The pellet, in Section 7.1, refines a displacement of minus
0.203 millimetres with its zero held. The two cells agree: a of 12.4740
against 12.4728 and c of 3.9295 against 3.9312, within about one and two esds,
and volumes of 611.426 against 611.581 cubic angstrom.

That agreement is the point. A zero point is a constant offset and a specimen
displacement follows the cosine of theta, and over a short angular range the
two are nearly the same parameter, so a fit that refines both will report
small esds for two numbers trading against each other. What decides which to
refine is not the data but the mounting: a powder bed is flush with the holder
and its offset is the instrument's, while a pellet surface sits where the
press and the polishing left it and its offset is geometric. Choosing the
wrong one does not fail. It gives a precise and wrong answer.

### 7.4 The same pellet done the wrong way

The command chooses by the sample's form, so the way to see the error is to
run the same scan as a bare file, where there is no form to consult and the
zero is refined instead.

```console
xrdkit lattice data/raw/pellet_a.xrdml --cell 12.45 3.94 --space-group P4bm --formula "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6" --z 5
```

```text
peaks   49 found, 43 refitted, 1 fits rejected, 5 satellites excluded, 1 recovered; 43 of 44 indexed (97.7 per cent), 33 used in the refinement
coarse window to 35.00 degrees
tetragonal cell a = 12.4801 +/- 0.0013, c = 3.9314 +/- 0.0005 angstrom
V = 612.324 +/- 0.167 cubic angstrom
zero 0.1707 +/- 0.0059 degrees
displacement not refined
rms 0.0131 degrees; a low indexed fraction or a high rms means the cell should not be trusted
M = 394.905 g/mol per formula unit, Z = 5
theoretical density 5.3546 +/- 0.0015 g/cm3
results\lattice\peaks_pellet_a.csv
results\lattice\lattice.csv
```

Read the two runs side by side. The displacement of the pellet has been
absorbed into a false zero of plus 0.171 degrees, which is the same number
Section 5 reported as the indexing zero of this scan. The cell has moved with
it: a from 12.4740 to 12.4801, an increase of 0.0061 angstrom, and the volume
from 611.426 to 612.324, an increase of 0.898 cubic angstrom, which carries
the density down from 5.3625 to 5.3546 grams per cubic centimetre. Against the
powder of the same composition the wrong run is out by 0.0073 angstrom in a
and 0.743 cubic angstrom in volume, while the right one sits within one esd of
it.

Now look at the residuals. The rms is 0.0131 degrees for the wrong run and
0.0128 for the right one. The fit is not worse. Nothing in the goodness of the
fit says which of the two answers to believe, and that is exactly why the form
is a property of the sample recorded in the project file rather than something
chosen after seeing the numbers.

Where a zero has been measured on a standard, give it with `--zero` and it is
held at that value whatever the form, which is the honest arrangement: the
zero comes from the instrument, and whatever offset is left over is the
sample's own and is fitted as a displacement.

### 7.5 The indexed fraction and the rms

The last line of the report carries its own warning, and it is not decoration.

The indexed fraction says how much of the pattern the cell accounts for. A
cell that is wrong in a way that matters usually leaves peaks unindexed, so
97.7 per cent for the pellet and 93.9 per cent for the powder are both
reassuring. The rms says how far the indexed peaks sit from where the cell
puts them, and a hundredth or two of a degree is what a good scan on a lab
instrument gives.

Neither number on its own proves the cell is right. The case to be wary of is
a pseudo-cubic pattern, where the sub cell reflections very nearly coincide
and a wrong cell with a spurious zero can be made to index every peak in the
list with a respectable residual. What that combination means is that the
model has enough freedom to follow the data wherever they go, and the fit has
stopped being a test. A low fraction or a high rms is therefore evidence
against a cell, but a high fraction and a low rms are not, on their own,
evidence for one. Check the cell against the phase identification of Section
6, against another sample of the same composition, and against a scan that
reaches to high angle.

Two more checks are worth making before a cell is quoted. The esds of a and c
should be about a thousandth of an angstrom for a density to be worth
calculating, and several thousandths is a sign that the high angle
reflections carried too little weight, which is a counting problem and not a
refinement problem. And a scan that stops at eighty degrees cannot give a
publishable cell at all, because the reflections that pin it down were never
measured and no refinement recovers them. Know when to remeasure rather than
refine.

### 7.6 The theoretical density, and the measured one

The density follows from the cell, the composition and the number of formula
units per cell: the mass of the contents of one cell divided by its volume.
The report gives the formula mass and Z before it so that both can be checked
at a glance, since a wrong Z is the commonest way to get a density wrong by a
clean factor.

`xrdkit density` does that calculation from a cell you already have, which is
the way to redo a density after a better cell arrives without rerunning the
refinement. The cell here is the one the pellet run gave, with its esds.

```console
xrdkit density --formula "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6" --z 5 --cell 12.4740 3.9295 --esd-cell 0.0011 0.0004
```

```text
M = 394.905 g/mol per formula unit
tetragonal cell a = 12.4740 +/- 0.0011, c = 3.9295 +/- 0.0004 angstrom
V = 611.433 +/- 0.125 cubic angstrom
theoretical density 5.3624 +/- 0.0011 g/cm3
results\density.csv
```

The volume differs from the refinement's in the third decimal place, and the
density in the fourth, because the refinement propagates the full covariance
of the fit while this calculation has only the two esds you typed. Quote the
refinement's number where you have it.

The relative density is the measured density as a percentage of this one, and
it is what says how well a pellet sintered. Give the measurement with
`--archimedes`, followed optionally by its own esd, and the report adds the
relative density with the two relative errors added in quadrature. For a
sample of the project file the measurement is better recorded once in the
sample table, as the `archimedes` key that `xrdkit add-sample --archimedes`
writes, and then `xrdkit lattice` picks it up without being told. Neither
example sample carries one here, which is why no run above prints a relative
density.

A word on what the theoretical density is worth. It is the density the
structure would have with no porosity and with the composition exactly as
weighed out, so it is an upper bound and not a prediction. Its relative error
is about three times the relative error in a lattice parameter, because the
volume goes as the cube of a length, which is why a cell good to 0.001
angstrom on a twelve angstrom axis is the target.

### 7.7 The files written

Two files come out of every lattice run, both under `results/lattice/KEY` for
a sample and under `results/lattice` for a bare file.

`peaks_STEM.csv` has one row per peak found and seventeen columns. The first
is `found_two_theta`, the position the peak finder gave. Then
`fitted_two_theta` and `esd_fitted_two_theta`, the position the doublet refit
gave and its esd, and `corrected_two_theta`, that position after the zero or
the displacement correction has been applied, which is the number the indexing
actually used. Then `d_spacing`, from the corrected position and the
wavelength. Then three flags: `fit_rejected`, true where the refit was thrown
out; `kalpha2_satellite`, true where the peak was flagged; and `recovered`,
true for a flagged peak the refined cell brought back in. Then `intensity`,
`relative_intensity` as a percentage of the strongest, and `fwhm`. Then the
indexing: `h`, `k` and `l` of the reflection assigned, `calculated_two_theta`
where the refined cell puts that reflection, `difference` between the
corrected and the calculated positions, and `n_candidates`, how many
reflections were within tolerance, which is worth reading because a peak with
several candidates is assigned on intensity and could have been assigned
otherwise.

`lattice_KEY.csv` is appended rather than rewritten, one row per run, so a
series of samples or a series of attempts on one sample builds up in a single
table ready to be plotted. Its columns carry, in order, what identified the
run and what went in: the sample, structure and scan file, the wavelength, the
form, the space group and crystal system, and the six start cell parameters.
Then what came out: the six refined cell parameters with their esds, the
volume with its esd, the zero and the displacement each with an esd and a flag
saying whether it was refined, and the radius. Then the seven counts of the
peaks line, the rms and the coarse window. Then the density block: the
formula, Z, the formula mass, the theoretical density and its esd, the
Archimedes density and its esd, and the relative density and its esd. The last
three columns are the standing record: the method in one sentence, the date
and the version of xrdkit that wrote the row. `xrdkit density` appends to
`results/density.csv` in the same way, with the inputs, the results and the
same three closing columns.

### 7.8 The options

Both commands share `--stem`, which names the output files and defaults to the
sample key or the scan's file stem, `--out`, the output root, and `--json`,
which prints the numbers and the files as JSON.

For `xrdkit lattice`, `SCAN` is a path or a sample key. `--cell` and
`--system` give the start cell and settle its crystal system by the count of
numbers, as Section 5.3 describes; for a sample they default to the first
structure's. `--space-group` names the group whose reflection conditions
apply, and defaults to the structure's; a symbol outside the eight groups
whose conditions are known is indexed without them, with a note.
`--wavelength` overrides the K alpha 1 wavelength.

`--zero` holds the zero at a value in degrees instead of refining it, and is
how a zero measured on a standard is imposed. `--displacement` refines the
specimen displacement for a scan that is not a pellet, and is unnecessary for
one that is, since a pellet always refines it. `--radius` is the goniometer
radius in millimetres and defaults to the instrument's; a displacement cannot
be refined without one.

`--no-satellites` excludes every peak flagged as a K alpha 2 satellite and
recovers none, which is the conservative setting where the recovered count is
high enough to worry about. `--formula` and `--z` give the composition and the
formula units per cell for the density, and default to the structure's for a
sample. `--archimedes` takes the measured density in grams per cubic
centimetre and optionally its esd, and adds the relative density to the
report.

For `xrdkit density`, the cell is given either as `--cell` with `--system`, in
which case `--esd-cell` takes the esds of those numbers in the same order, or
as `--volume` with `--esd-volume`, which is the way in for a crystal system
the kit does not otherwise handle or for a volume from elsewhere. `--formula`
and `--z` are required unless the optional `SAMPLE` argument is given, in
which case the formula, Z, cell and Archimedes density not given as options
come from the sample and its first structure, and the row goes to
`results/density/KEY` instead.

## 8. Le Bail and Rietveld

The two commands of this section drive GSAS-II. `xrdkit lebail` extracts the
cell from the whole profile rather than from peak positions alone, and
`xrdkit rietveld` carries that result on into a refinement of the structure
itself.

### 8.1 What they are for, and what each needs

A Le Bail extraction gives every reflection a free intensity and fits the
whole pattern with them. Nothing structural is refined and nothing structural
is learned, which is the point: the cell, the zero, the size and the
microstrain are separated from the intensities rather than fought over with
them. That makes it the better cell than the one Section 7 gives, because it
uses the whole profile and not the positions of the peaks a peak finder
happened to resolve, and it is the step that prepares a Rietveld refinement.

A Rietveld refinement calculates the intensities from a structure and refines
the structure until they match. It answers what the atoms are doing, and it
costs beam time and judgement in proportion.

Both need GSAS-II, installed as Section 1.2 describes, and both need the
instrument parameter file of Section 3, named in the instrument table, so that
the peak shape belonging to the diffractometer is held rather than refined.

Where they differ is the structure. A Le Bail extraction needs no coordinates,
so a structure given as a `library` entry alone is enough: the entry's space
group and cell decide where the reflections are, and the command writes a CIF
of that entry for GSAS-II to read. The Rietveld modes put atoms in, so they
need a `cif` beside the entry, and stop with a message saying so when there is
none.

### 8.2 The refine table

Both commands take what they refine with from the `[refine]` table of the
project file, and a sample's own `[samples.KEY.refine]` table overrides any
key of it. Every key is optional, and a table given in part keeps the rest of
its values, so `background = { terms = 8 }` changes the number of terms and
leaves the function alone.

| Key | Default | What it sets |
| --- | --- | --- |
| `two_theta` | the scan's range | the range refined, `[low, high]` in degrees, clipped to the scan |
| `background` | `{ function = "chebyschev-1", terms = 6 }` | the background function and its number of terms |
| `max_passes` | `{ lebail = 60, fixed_atoms = 60, coordinates = 100, occupancies = 100 }` | the most passes of a stage, by mode |
| `unsettled` | `accept` for every mode | whether a stage still moving after those passes is kept or rolled back |
| `followed` | none | reflections to follow from mode to mode, as `[h, k, l]` triples |

`two_theta` is the key worth setting per sample, and Section 8.4 shows why.
`background` should be as few terms as will follow the background and no more,
since a flexible polynomial will absorb intensity that belongs to the peaks.
`max_passes` and `unsettled` are about how hard the refinement tries and what
becomes of a stage that is still moving when it stops, both of which Section
8.5 explains. `followed` names reflections to report in every mode, which is
how a particular superstructure line is watched across a refinement.

### 8.3 Running a Le Bail extraction

The command takes a sample key and nothing else. The start cell is the
sample's lattice result from Section 7 where there is one, its structure's
cell otherwise, or `--cell` where that is given.

```console
xrdkit lebail powder_a
```

It reports each stage as it finishes, then the files, then the numbers. The
run took about seven minutes.

```text
powder_a: lebail started, 5 stages, at most 60 passes each
powder_a: lebail: background and scale clean
powder_a: lebail: zero clean
powder_a: lebail: cell clean
powder_a: lebail: size clean
powder_a: lebail: microstrain clean
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail_result.json
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail.gpx
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\lebail.md
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail_histogram.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail_reflections_ttb_p4bm.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\powder_a_lebail.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a\summary.md
ttb_p4bm: tetragonal cell a = 12.4740 +/- 0.0003, c = 3.9318 +/- 0.0001 angstrom
V = 611.789 +/- 0.042 cubic angstrom
size 0.2051 +/- 0.0088 micron, microstrain 1011 +/- 98
zero -0.0357 +/- 0.0008 degrees
Rwp 3.888 per cent, reduced chi squared 1.799
start cell of ttb_p4bm from lattice results C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\powder_a\lattice_powder_a.csv
```

Read the last line first: the start cell came from the lattice result of
Section 7, so the two workflows are joined without either being told about the
other, and the record says which file was used.

Then compare the cell with the one that file holds. The lattice refinement
gave a = 12.4728(16) and c = 3.9312(6); the extraction gives a = 12.4740(3)
and c = 3.9318(1). The two agree within the lattice esds, and the extraction's
esds are five times smaller, which is what using the whole profile buys. Quote
this cell rather than the other one.

The zero of minus 0.0357 degrees is the same instrument zero the lattice run
found, to within its esd. An Rwp of a few per cent on a Le Bail fit and a
reduced chi squared approaching 2 is a sound result for a scan of this length:
a Le Bail fit has an intensity free for every reflection, so it should fit
better than a Rietveld refinement of the same pattern, and a reduced chi
squared well above 2 usually means the counting statistics are better than the
model rather than that the model is wrong.

The size and the microstrain are the two numbers to treat carefully. Both
broaden the peaks, and they are separated only by how that broadening grows
with angle, so on a laboratory scan they can trade against each other. Here
the microstrain comes out at 1011 with an esd of 98, which is ten esds from
zero and therefore a real measurement rather than a parameter absorbing noise.
The next section shows what the same sample gives when that is not true.

### 8.4 The same sample over the whole scan

The range for `powder_a` is set to 17 to 99.98 degrees by the per sample table
of Section 2.2, rather than the scan's own 10.01 to 99.98. Comment that table
out and the extraction runs over everything measured.

```console
xrdkit lebail powder_a --out results/lebail/powder_a_full
```

It writes the same nine files, to the folder `--out` names. The stage lines
and the cell are almost unchanged, and everything else is worse.

```text
powder_a: lebail started, 5 stages, at most 60 passes each
powder_a: lebail: background and scale clean
powder_a: lebail: zero clean
powder_a: lebail: cell clean
powder_a: lebail: size clean
powder_a: lebail: microstrain clean
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail_result.json
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail.gpx
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\lebail.md
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail_histogram.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail_reflections_ttb_p4bm.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\powder_a_lebail.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lebail\powder_a_full\summary.md
ttb_p4bm: tetragonal cell a = 12.4739 +/- 0.0003, c = 3.9318 +/- 0.0001 angstrom
V = 611.787 +/- 0.045 cubic angstrom
size 0.1266 +/- 0.0037 micron, microstrain -378 +/- 104
zero -0.0358 +/- 0.0009 degrees
Rwp 4.829 per cent, reduced chi squared 2.666
start cell of ttb_p4bm from lattice results C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\lattice\powder_a\lattice_powder_a.csv
```

The microstrain has gone negative, minus 378 with an esd of 104. A microstrain
below zero is not a physical quantity; it is the refinement narrowing the
calculated peaks because something else is making them too wide. The size has
moved with it, from 0.205 to 0.127 microns, which is the correlation between
the two showing itself: each has taken up what the other left. And the fit is
worse on both measures, Rwp 4.829 against 3.888 and reduced chi squared 2.666
against 1.799.

What causes it sits at the bottom of the range. The structure expects two
reflections at low angle that this pattern does not show, and a Le Bail
extraction cannot take an intensity to zero, so it has nothing to fit them
with. It drags the background up to cover them and pushes the microstrain
negative to narrow them, after which the size and the microstrain correlate
strongly and neither is worth reading. Starting the range above those
reflections removes the correlation, which is why a sample whose low angle
reflections are absent or very weak is given a range of its own rather than
the scan's.

Notice what did not change. The cell is the same to the fourth decimal place
and the zero to the fourth, so nothing in the numbers you were after says that
anything went wrong. The evidence is in the microstrain sign, in the size
moving, and in the residuals, and a reader who looked only at the cell would
have taken the second run as readily as the first.

### 8.5 Stages and passes

A refinement is a list of stages, and each stage adds flags to the ones before
it, so by the last stage everything named along the way is refining together.
The order is not a matter of taste: each stage frees parameters that only make
sense once the ones before them are near right. A Le Bail extraction runs five
of them, and the printed lines name each as it finishes. Background and scale
come first, then the zero, then the cell, then the size, and last the
microstrain.

One new kind of parameter at a time is the rule. Freeing the cell before the
background has settled lets the background take up intensity that belongs to
the peaks, and freeing the size and the microstrain together, as the last two
stages deliberately do not, gives two numbers that trade against each other.

Within a stage the refinement is repeated in passes. GSAS-II stops at the
first least squares cycle that raises chi squared, which is often well short
of the minimum, so the kit runs the stage again, and again, until nothing but
the scale moves by more than a tenth of an esd from one pass to the next, or
until the pass cap for that mode is reached. The cap is sixty passes for a Le
Bail extraction, and the write up records how many each stage took: eleven,
seventeen, ten, twenty and twelve for the run above.

Each stage is checked when it finishes, and its status says what became of it.
A clean stage settled within the cap and raised nothing against itself, as all
five did here. A stage still moving at the cap is recorded as unsettled and is
kept or rolled back according to the `unsettled` key. A stage that raises a
sanity check is rolled back, the refinement returns to the state the last
stage kept left it in, and what that stage alone refined stays held for the
rest of the run. The point of rolling back rather than stopping is that a
sequence can be left to run and still leave a usable record.

### 8.6 The files, and reading the write up

Nine files go to `results/lebail/KEY`, and the logs go under `work` beside
them.

The result JSON is the machine readable record of the whole run: the inputs,
every stage with its residuals and parameters and esds, the final model, the
method and the date. Everything else is derived from it. The `.gpx` is the
GSAS-II project, which opens in the GSAS-II interface if you want to look at
the refinement yourself, with a `.bak0.gpx` beside it and a `.lst` listing.
The histogram CSV holds the observed, calculated, background and difference
curves point by point, for a figure of your own; the reflections CSV holds
every reflection with its position, its indices and its extracted intensity.
The `.instprm` is the instrument file exactly as the run used it, kept so that
the refinement can be reproduced even if the original is edited later. The png
and the pdf are the fitted pattern with its difference curve. `summary.md`
gathers a run of several modes into one table, and matters more in Section 8.8
than here. Where a mode fails, a `failure.md` is written in place of the write
up, carrying the error, the stages that did finish and the tail of the GSAS-II
log.

`lebail.md` is the one to read. It opens with the settings, every one of them:
the sample and its scan, the instrument and its parameter file, the range, the
background, the phase and where its cell started, the pass rule, the driver
and xrdkit versions, and the time of the run. That header is what makes the
result reproducible a year later.

Then comes a table of stage outcomes, one row per stage with its status, its
passes, Rwp, Rp and reduced chi squared, and a column saying why where a
status is not clean. Reading down the residual columns tells you which stage
earned its place: here Rwp fell from 8.182 to 7.313 when the zero was freed
and from 7.248 to 3.940 when the size was, and the microstrain stage moved it
only from 3.940 to 3.888.

Then the cell, refined against start with the change in each parameter, so
that a cell that has walked a long way from where it started is obvious. Then
the zero or the displacement, with the value it started from.

Last comes the size and microstrain section, and it ends with a verdict rather
than a number. For this run it reads that the microstrain test finds a
microstrain, 10.3 esds from zero, with Rwp falling from 3.940 to 3.888 per
cent. That is the test the fifth stage exists to perform: the microstrain is
refined last and on its own so that its effect on the residual can be seen,
and the verdict says whether the number it produced is worth keeping. A
microstrain a fraction of an esd from zero, or one that improves the residual
not at all, should be read as no microstrain rather than as a small one.

Where the `followed` key names reflections, a further section reports each of
them in every mode, with its observed and calculated structure factors, which
is how a superstructure line is watched from one refinement to the next.

### 8.7 The options

`SAMPLE` is a sample key of the project file, and is the only argument. There
is no way to name a scan file here: a refinement needs the instrument, the
structure and the form, and those live in the project file.

`--cell` and `--system` give the start cell of the first phase, in place of
the sample's lattice results or its structure's cell, as the free parameters
of its crystal system. Use it where the lattice result is one you have decided
against.

`--zero` starts the zero at a value in degrees rather than at the
instrument's. `--displacement` refines the specimen displacement with the zero
held, which is what a pellet does in any case, and is the flag for a powder
you have reason to think sits proud of the holder.

`--out` sends every file to a folder of your choosing rather than to
`results/lebail/KEY`, which is how Section 8.4 keeps two runs of one sample
side by side. `--json` prints the files written and the outcome as JSON, with
the progress lines on the error stream so that the JSON stands alone.

Every phase of the sample is extracted, not only the first; the first is the
one whose cell `--cell` sets and whose numbers are printed first.

### 8.8 Rietveld

The Rietveld half of this section, covering `xrdkit rietveld`, its modes and
stages, the write ups it produces and the undetermined parameters it reports,
follows here.
