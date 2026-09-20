# xrdkit user guide

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

### 8.8 Running a Rietveld refinement

`xrdkit rietveld` takes the Le Bail result of Section 8.3 and carries it
through three modes in turn. With no options it runs all three.

```console
xrdkit rietveld powder_a
```

It reports each stage as it finishes, then every file, then one line per mode.

```text
powder_a: fixed_atoms started, 4 stages, at most 60 passes each
powder_a: fixed_atoms: scale and background clean
powder_a: fixed_atoms: zero and cell clean
powder_a: fixed_atoms: size clean
powder_a: fixed_atoms: overall Uiso clean
powder_a: coordinates started, 5 stages, at most 100 passes each
powder_a: coordinates: profile clean
powder_a: coordinates: Uiso groups clean
powder_a: coordinates: A sites clean
powder_a: coordinates: B sites clean
powder_a: coordinates: O sites rejected (sanity check: O2 z moved +0.0872; sanity check: O5 z moved +0.0743)
powder_a: occupancies started, 2 stages, at most 100 passes each
powder_a: occupancies: profile and Uiso clean
powder_a: occupancies: A site occupancies clean
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms_result.json
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms.gpx
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\fixed_atoms.md
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms_histogram.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms_reflections_ttb_p4bm.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_fixed_atoms.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates_result.json
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates.gpx
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\coordinates.md
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates_histogram.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates_reflections_ttb_p4bm.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_coordinates.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies_result.json
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies.gpx
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\occupancies.md
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies.png
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies.pdf
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies_histogram.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies_reflections_ttb_p4bm.csv
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\powder_a_occupancies.instprm
C:\Users\amirk\AppData\Local\Temp\claude\C--Users-amirk-Source-repos-xrdkit\742c29b6-8095-4509-8f93-500b972aa5a3\scratchpad\guide_project\results\rietveld\powder_a\summary.md
fixed_atoms: scale and background, zero and cell, size, overall Uiso; Rwp 4.290 per cent, reduced chi squared 2.190
coordinates: profile, Uiso groups, A sites, B sites; Rwp 4.158 per cent, reduced chi squared 2.062
occupancies: profile and Uiso, A site occupancies; Rwp 4.088 per cent, reduced chi squared 1.990
```

The whole run took under four minutes, and it was not spent evenly. The fixed
atoms mode took about twenty seconds and the occupancies mode about seventeen,
while the coordinates mode took a little over three minutes, almost all of it
in one stage that needed sixty-one passes to settle. That is the usual shape:
the stages that free many weakly determined parameters at once are the slow
ones.

Read the three closing lines first. Each names the stages that were accepted,
so a stage missing from that list was rolled back, and gives the residuals of
the model the mode ended with. Rwp fell from 4.290 to 4.158 to 4.088 per cent
across the three modes, and the reduced chi squared from 2.190 to 1.990.

### 8.9 The three modes

Each mode starts from the saved result of the one before it, and the first
starts from the Le Bail result, so the sequence is a chain of files rather
than one long refinement held in memory. That is what makes it possible to
rerun one mode without repeating the others.

| Mode | Starts from | Frees | Holds |
| --- | --- | --- | --- |
| `fixed_atoms` | the Le Bail size stage, or its microstrain stage with `--mustrain` | the scale, background, zero or displacement, cell, size, and one Uiso for every atom | the coordinates and occupancies of the CIF, at the nominal composition |
| `coordinates` | the fixed_atoms result | the profile, the entry's Uiso groups, then the free coordinates of each kind of site in turn | the origin site along its axis, and the occupancies |
| `occupancies` | the coordinates result | the profile and Uiso groups, then the occupancies of each exchange group over its sites | the coordinates, and each element's content over the group's sites |

The fixed atoms mode puts the structure in at its nominal composition and asks
only whether the profile can be made to fit with the atoms where the CIF has
them. Nothing structural is judged until it can, because a peak in the wrong
place produces intensity errors that look exactly like occupancy errors. The
one Uiso every atom shares is always determined by the data and soaks up the
overall falling off of intensity with angle that would otherwise be pushed
into the occupancies. If it comes out negative, that never means the atoms are
colder than still; it means intensity is missing at high angle, so look at
absorption, at the sample height, or at the instrument file.

The coordinates mode frees one kind of site at a time, heaviest scatterers
first. Heavy atoms dominate the intensities, so their positions are the best
determined, and moving them first stops the light atoms from chasing errors
that are not theirs. A kind of site with nothing free, every site of it fixed
by symmetry or holding the origin, has no stage of its own.

The occupancies mode trades the exchanged elements between their sites under
the composition constraint, so the refinement decides how the elements are
distributed without being free to change how much of each the sample contains.
It is usually the longest and the one most often left undetermined.

The stage order and the pass rule are those of Section 8.5: one new kind of
parameter at a time, each stage repeated until nothing but the scale moves by
more than a tenth of an esd, and a stage that raises a sanity check rolled
back. The pass cap differs by mode, sixty for fixed atoms and a hundred for
the other two, and the printed lines say which cap is in force.

### 8.10 Choosing modes with --from and --through

`--from` and `--through` run part of the sequence, which is what to use after
changing something that affects only the later modes.

Because each mode reads the saved result of the one before, the result it
needs has to be there. `--from coordinates` needs the fixed atoms result in
the folder, and where it is missing the command stops before writing anything,
naming the file and the command that writes it. `--from` after `--through` is
refused outright.

### 8.11 How the structure is set up from the project file

The mode that puts the atoms in has to reconcile three things: the sites the
library entry names, the atoms the CIF actually has, and the composition the
sample was weighed out to. The `atoms` table of the structure is what joins
them, and Section 2.2 gives the one this project uses.

An element of the composition that the CIF does not carry is placed only on
the sites the `atoms` table names it on, and on each of those it goes beside
the atom of its host element, the host being the one element present on every
one of those sites. Its amount is split among them in proportion to
multiplicity times the host's occupancy, so a site with more room for it, or
more of its host, takes more of it. Where no single element is present on all
the named sites the placement is ambiguous, and the file is refused when it
loads rather than at the end of a refinement.

This project file shows both cases. Lanthanum is named on A1 alone, so all of
it goes there, beside `Sr1`; had the table named A1 and A2, it would have been
split between them in proportion to the strontium on each. Titanium is named
on B1 and B2, and niobium is the element present on both, so the titanium is
split between the two B sites in proportion to the niobium. Leaving such an
element out of the table altogether stops the run with a message naming the
element and the entry's sites.

The other half of the table is `library` beside `cif`. The entry supplies the
space group, the sites and their kinds, the coordinates each Wyckoff position
leaves free, the Uiso groups, the anions and the bond limits; the CIF supplies
the coordinates. Where the CIF labels its atoms differently from the entry's
site names, the `atoms` table says which atom sits on which site, since a
refinement matches each site to the CIF atom of the same label.

`origin` names the site whose coordinate along a polar axis is held, to stop
the whole structure sliding along it. Left out, the entry's own applies, which
for this tungsten bronze is its B1 site. A site label of the entry chooses
another and `origin = false` holds none, and in that case the coordinates mode
refuses to free a coordinate along the polar axis, naming the axis, because
the refinement would have nothing to fix the structure in place along it.

`exchange` lists the groups of elements whose occupancies are traded, one
occupancy stage per group, named from the kind of site it works on. The sites
a group is traded between are every site of one kind that holds any of its
elements. Where a refinement frees phase fractions as well, the histogram
scale is held whenever it does, since the scale and the fractions are the same
quantity counted twice.

### 8.12 The options

`SAMPLE` is a sample key, and as with `xrdkit lebail` there is no way to name
a scan file.

`--from` and `--through` choose the modes, as Section 8.10 describes.

`--mustrain` refines the microstrain with the size. It is held at zero in
every Rietveld mode by default, because size and microstrain trade against
each other on a laboratory scan, and the Le Bail run is where that question is
settled; give it only when that run found a microstrain worth carrying.

`--preferred-orientation H K L` adds a last fixed atoms stage freeing a
March-Dollase ratio about that axis. Say plainly what its state is: preferred
orientation has not yet been exercised on a real refinement. The stage is
built and its flags are checked, but no measured scan has been carried through
it, so check any result it gives against a scan of the same sample loaded to
limit texture before relying on it.

`--out DIR` writes every file into that folder itself, not under a `results`
folder inside it, so a Rietveld run given the same `--out` as the Le Bail run
before it finds the Le Bail result waiting there. A relative path is taken
from the folder the command is run in.

`--overwrite` is needed to replace a result that is already there. The command
refuses by default: where the result JSON of any mode it is to run is already
in the folder it stops before refining or writing anything, names the files,
and returns 1. A mode it is not running is left alone, so `--from coordinates`
after a fixed atoms run needs no `--overwrite`. There is no dated subfolder,
because each mode has to read the result of the one before from the same
place.

`--json` prints the files written and the outcome of each mode as JSON, with
the progress lines on the error stream.

### 8.13 The files written

Every mode writes the same eight files under `results/rietveld/KEY`, named
from the sample key and the mode, and the run as a whole adds `summary.md`.

The result JSON is the record of the mode: every stage with its status,
residuals, parameters, atoms and every value refined or held; the final model;
and three keys that make it reproducible. `inputs` holds the project, the
sample, the instrument, every structure as the project file gives it, the
refine settings, the scan's range and the range refined, the start cell of
every phase and where it came from, whether the displacement was refined, the
start zero and the command's options. `method` holds the mode, the stages as
the job gave them, the cycles, the pass cap and tolerance, the driver's
version and the version of xrdkit. `date` is when the mode finished.

Beside it are the GSAS-II project, the fitted pattern and the reflections as
CSVs, the instrument parameters exactly as that mode used them, the figure as
a png and a pdf, and `MODE.md`, the write up. `summary.md` gathers the modes
into one table with each mode's outcome, residuals, cell and the status of
every stage, which is the file to open first when a run is a few days old.
Where a mode fails, `failure.md` is written in its place, carrying the error,
the stages that mode did get through and the tail of its GSAS-II log, and the
modes after it are not run. The jobs, the logs and the start instrument file
of each mode are under `work/MODE`.

A failure in drawing the figure or writing the write up happens after the
result JSON is written, so the result is kept and `failure.md` says so.

### 8.14 Reading the write ups

Each write up opens with the settings, in the same form as the Le Bail one,
then a stage outcomes table with one row per stage: its status, its passes,
Rwp, Rp and reduced chi squared, and a column saying why where the status is
not clean. A rejected stage keeps its row and its reason but carries no
figures, since its numbers are not part of the model. The line under the table
names the stage the final model came from.

Then the cell with its esds beside the cell it started from, and the zero or
displacement beside its start.

The structure write ups add four things. Every atom's coordinates, occupancy
and Uiso with esds, and its shift from the start model in angstroms in the
refined cell, so that a shift can be read as a distance rather than as a
fraction. The bond lengths from each cation site to the anions, with the
limits the library entry gives for that kind of site and a flag on any bond
outside them, which is the quickest check that a structure has not gone
somewhere impossible. The occupancies table, giving for each element of an
exchange group the sites it was traded between, the occupancy of each of its
atoms after the stage, the content per cell held while it was traded, and the
content the stage's own atoms give, which differs from the held total only if
the constraint failed. And the reflection misfits, the reflections whose
observed and calculated intensities differ most, beside the ratio the Le Bail
fit gives, which owes nothing to the structure: where the two agree, the
misfit is the structure's.

Last comes the undetermined list, and it is the section to read before quoting
anything. A parameter is undetermined when the data do not determine it: an
occupancy whose esd is more than half its range of 0 to 1, or a coordinate or
Uiso whose esd is larger than its shift from the start model. The start model
is the structure as the fixed atoms mode set it up, the CIF at the nominal
composition, and it is carried to every later mode, so a shift is always
measured from the same place.

This is independent of a stage's status, and a clean stage can leave
undetermined parameters. What it means for a number destined for a paper is
simple: an undetermined value is not a result, because the data cannot tell it
from where it started. Do not tabulate it as refined. Either hold it and say
so, or leave that part of the structure out of the claim. A referee can tell
the difference, and the labels are in the result precisely so that the
decision is not left to memory.

### 8.15 What this refinement came to

Take the run above as it stands, mode by mode.

The fixed atoms mode came out clean on all four stages, with Rwp falling from
5.034 to 4.290 per cent and the reduced chi squared from 3.012 to 2.190.
Almost all of that came from the last stage: freeing the scale, the background,
the zero, the cell and the size together moved Rwp only from 5.034 to 4.997,
and the one overall Uiso took it to 4.290. That is the Uiso doing exactly what
Section 8.9 says it does, taking up the fall of intensity with angle.

The coordinates mode came out clean on four stages and rolled the fifth back.
The profile and Uiso groups stages settled in two passes each. The A sites
stage needed sixty-one passes to settle and moved Rwp from 4.222 to 4.198, and
the B sites stage settled in four and took it to 4.158.

The O sites stage was rejected, and the reason names what happened: O2 moved
0.0872 along z and O5 moved 0.0743, against the 0.05 fractional the sanity
check allows. Note what the residual was doing while that went on. Inside the
stage Rwp had fallen to 4.036, better than the 4.158 the mode ended with, and
the stage was rejected all the same. That is the whole argument for having a
sanity check: the residual improved and the structure got worse, because
weakly determined oxygen coordinates can always find somewhere to go that fits
the noise a little better. The stage was rolled back, what it alone refined
was held, and the B sites model is the answer.

Worth recording, since it is easy to assume otherwise: the O sites stage did
not fail by running out of passes. It settled, within the hundred pass cap,
and was rejected on the shift check afterwards. A stage that reaches the cap
is a different label, unsettled, and is kept rather than rolled back.

The occupancies mode came out clean on both stages, taking Rwp to 4.088 per
cent and the reduced chi squared to 1.990, with the occupancy stage settling
in six passes. And here the labels matter more than the residual. The
refinement put strontium at 0.2(9) on A1 and 0.4(4) on A2, and barium at
0.2(6) on A1 and 0.5(3) on A2, with lanthanum held at its nominal 0.25 on A1.
The composition constraint held: strontium keeps 2.0 atoms per cell and barium
2.5, and the stage's own totals match.

Two of those four are named in the undetermined list, the strontium and the
barium on A1, each with an esd larger than half the range an occupancy can
take. The other two escape the test but should not be quoted either: an
occupancy of 0.4 with an esd of 0.4 is not a measurement. What this pattern
says about the A site distribution is nothing at all, and the honest report is
that the split was refined, came out undetermined, and is therefore left as
the nominal composition.

That is a data limitation rather than a model error, and the two have
different signatures. A data limitation shows as parameters that are
individually undetermined while the fit as a whole is good: correlated
coordinates that trade against one another, mixed occupancies on a shared site
that the scattering contrast cannot separate, esds larger than shifts. That is
this run, and it is not surprising: strontium and barium differ by eighteen
electrons on sites that also carry lanthanum, on a scan Section 4 reported as
unsuitable for a Rietveld refinement. A model error shows instead as a fit
that is bad in a structured way: a whole class of reflections fitted badly
points at the wrong space group, intensity under peaks with no tick mark at a
missing phase, and a residual good at low angle and bad at high angle, or the
reverse, points at the composition or at the thermal parameters.

What can be quoted from this run is the cell, the profile terms and the
residuals, with the coordinates of the A and B sites, and the statement that
the oxygen positions and the A site occupancies were not determined. That is a
real result and a modest one.

### 8.16 When a stage is rolled back, and when a mode fails

A rolled back stage is not an error and does not stop the run. The refinement
returns to the state the last kept stage left it in, whatever that stage alone
was refining stays held, and the run carries on, so a sequence can be left
unattended and still leave a usable record.

What to do next depends on which label appeared. A rejected stage means the
model went somewhere the sanity check forbids, and the first question is
whether the check was right. Look at the shift it names in the write up's
coordinates table, in angstroms rather than fractions, and at the bond lengths
for that site: a cation to anion distance outside the entry's limits says the
structure really had gone wrong. If the shift is large but the bonds are still
sound, the start model may simply have been poor, and the answer is a better
CIF rather than a looser check.

An unsettled stage means the pass cap arrived while something was still
moving. Raising the cap in the `max_passes` key is worth one attempt. If it
still will not settle, the parameter it names is being traded against another
and the refinement is not going to resolve them.

A failed mode is different: GSAS-II itself raised an error, `failure.md`
carries it with the tail of the log, and the modes after it are not run. Read
the error first, then the last stage that did finish, because a failure
usually follows a stage that had already gone somewhere unreasonable.

In every case the question to ask before changing the strategy is whether the
data can answer what is being asked of them. The Rietveld row of the table in
Section 4.2 is the standard: 5 to 130 degrees, a strongest peak above 20000
counts, a background above 200 counts per step at high angle, powder sieved
below 45 micrometres and loaded to limit preferred orientation. The scan used
here reaches 98 degrees with a strongest peak near 10000 counts, which is a
good Le Bail scan and a marginal Rietveld one. That, and not the refinement
strategy, is why the A site occupancies would not resolve, and no rearrangement
of the stages will change it. Go back for better data when the limitation is
the data.

## 9. What the commands write

Every command prints one line per file it writes, so the terminal is always
the first answer to where something went. This section is the second: what
each file holds, and where the record of how it was made is kept.

The rule the kit follows is that any file carrying a number you might put in a
paper also carries the inputs it came from, the method in one sentence, the
date and the version of xrdkit that wrote it. For a CSV those are columns; for
a result JSON they are the `inputs`, `method` and `date` keys; for a write up
they are the settings block at the top. Nothing that matters is left to be
remembered.

Running every command of Sections 2 to 8 once, as this guide did, adds 144
files to the project folder, counting the project file itself and the eleven
CIFs fetched from the database. Most of them are the GSAS-II working files of
the refinements, which are kept so that a run can be reproduced or inspected
and which nobody reads day to day.

### 9.1 init

`xrdkit.toml`, the project file, is the only file written, and `data/raw`,
`cifs` and `results` are created if missing. Nothing is recorded in it but the
project name and the format version, because nothing has been measured yet. An
existing project file is never overwritten.

### 9.2 add-sample

`xrdkit.toml` is appended to, with a `[samples.KEY]` table, and the rest of the
file is left byte for byte as it was. Nothing else is written. The whole file
is checked with the new table in place before anything is saved, and a key
already present is refused.

### 9.3 instrument

Two files you will use and several you will not, all under `data/standards`
when the command is run in a project.

`data/standards/STEM.instprm` is the instrument parameter file itself, the one
GSAS-II exported, and is what every later refinement reads.
`data/standards/STEM_instrument.csv` is its record, appended one row per run,
with the scan, the CIF, the phase, the wavelength, the six cell parameters
held, the fitting window, the width fit's peak count, U, V, W and rms, every
refined parameter with its esd, the Rwp and goodness of fit, the radius, the
path of the file written, and then the method, date and version.

`data/standards/STEM_start.instprm` is the starting file the refinement began
from, kept so that the two are never confused, and
`data/standards/work/STEM/` holds the GSAS-II project, its backup and listing,
the two column copy of the scan, the exported histogram and reflections, the
job as JSON, the result as JSON and the log.

With `--name`, `xrdkit.toml` is appended to as well, with an
`[instruments.KEY]` table.

### 9.4 check

For a sample key, `results/check/KEY/check_KEY.txt`, the report exactly as it
was printed, headed by the key and composition. For a scan named as a path,
nothing is written at all. With `--json` the report goes to the terminal as
JSON instead and, for a sample, `check_KEY.json` is written in place of the
text.

The report is a snapshot of the scan rather than a derived number, so it
carries no method or version line; what makes it reproducible is that it is
computed from the scan alone.

### 9.5 plot

Six files for a sample with a cell to index against, under
`results/plot/KEY`: `peaks_KEY.csv`, `pattern_KEY.png` and `.pdf`,
`indexed_KEY.csv`, and `pattern_hkl_KEY.png` and `.pdf`. Without a cell the
indexing file and the hkl figure are not written. For a scan named as a path
the peak list and indexing go to `results` and the figures to `figures` under
the working folder, which is why this guide's second plot run left
`results/peaks_powder_a.csv`, `results/indexed_powder_a.csv` and four files in
`figures`.

`peaks_KEY.csv` gives each peak its position, intensity, prominence, width, d
spacing and relative intensity. `indexed_KEY.csv` gives each peak the
reflection assigned to it, how far the two differ and how many candidates were
within tolerance. Both are replaced on a rerun, and neither carries a method
line: the refined cell they came from is printed and belongs in the lattice
record rather than here.

### 9.6 stack

Two files, `stack_STEM.png` and `.pdf`, under `results/stack/STEM` when a
sample is among the scans and under `figures` otherwise. Nothing else, and
both are replaced on a rerun.

### 9.7 phases

The CIFs go to `cifs/cod/<id>.cif`, one per candidate, and are not fetched
again once they are there. `cifs/cod/index.csv` is the register of reference
CIFs, with the file, source, identifier, formula, space group, the six cell
parameters, the reference and your notes. It is updated rather than rewritten,
so a row you have annotated survives a fresh download.

Under `results/phases/KEY` go three files. `phases_KEY.csv` is the ranking,
with the rank, COD id, formula, space group, the explained and missing counts,
the score, the rejection reason and the path of the CIF.
`phases_KEY_unexplained.csv` is one row per peak no candidate explained, with
its two theta, d spacing, and the main phase and reflection that account for
it or nothing where none does. `phases_KEY_record.csv` is the record of the
run, appended one row per run: the scan, the sample, the wavelength, the
elements, the space group filter, the zero, the window, the tolerance, the
peaks observed, the candidates searched and fetched, the peaks left over, the
main phase, the index path, and the method, date and version. The ranking and
the unexplained list are replaced on a rerun; the record grows.

### 9.8 lattice

Two files, under `results/lattice/KEY` for a sample and directly under
`results/lattice` for a scan named as a path.

`peaks_STEM.csv` is one row per peak with the seventeen columns Section 7.7
describes, from the found position through the flags to the assigned
reflection, and is replaced on a rerun. `lattice_KEY.csv` is appended one row
per run, so a series of samples or a series of attempts builds into one table,
and carries the inputs, the refined cell and volume with esds, the zero and
displacement with their refined flags, the radius, the seven peak counts, the
rms and coarse window, the density block, and the method, date and version.

### 9.9 density

One file, `results/density.csv`, appended one row per run, or
`results/density_STEM.csv` with `--stem`, or
`results/density/KEY/density_KEY.csv` for a sample key. Its columns are the
formula and Z, the crystal system and the six cell parameters with esds, the
volume with its esd, the formula mass, the theoretical density with its esd,
the Archimedes density with its esd, the relative density with its esd, and
the method, date and version.

### 9.10 lebail

Nine files are printed, under `results/lebail/KEY` or under the folder
`--out` names, and two more sit beside them unannounced: the project backup
and the GSAS-II listing. With the four working files that is fifteen in all.

`KEY_lebail_result.json` is the full record: every stage with its status,
residuals, parameters and atoms, the final model, and the `inputs`, `method`
and `date` keys that make the run reproducible. `KEY_lebail.gpx` is the GSAS-II
project, with a `.bak0.gpx` backup and a `.lst` listing beside it.
`KEY_lebail_histogram.csv` holds the observed, calculated, background and
difference curves point by point, and `KEY_lebail_reflections_PHASE.csv` every
reflection with its position, indices and extracted intensity.
`KEY_lebail.instprm` is the instrument file exactly as the run used it.
`KEY_lebail.png` and `.pdf` are the fitted pattern. `lebail.md` is the write
up, whose settings block carries the same record in prose. `summary.md`
gathers the run into one table.

`work/lebail/` holds the job as JSON, the two column copy of the scan, the
start instrument file and the GSAS-II log.

All of these are replaced on a rerun of the same mode into the same folder,
which is why `--out` exists.

### 9.11 rietveld

Eight files per mode are printed, under `results/rietveld/KEY`, named from the sample key
and the mode, and the same shape as the Le Bail set: the result JSON with its
`inputs`, `method` and `date` keys, the GSAS-II project with its backup and
listing, the histogram and reflections CSVs, the instrument file as used, the
figure as png and pdf, and `MODE.md`. One `summary.md` covers the modes
together. `work/MODE/` holds the jobs, the logs, the start instrument file and
the CIF the mode built for GSAS-II.

Counting the project backup and the listing that go unprinted, a full three
mode run leaves sixty-three files: thirty-one in the folder itself and
thirty-two under `work`. That is what the `results/rietveld/powder_a` folder
of this guide holds.

Unlike every other command, `xrdkit rietveld` refuses to replace a result it
would write. Where the result JSON of any mode it is to run is already in the
folder it stops before refining or writing anything, names the files, and
returns 1. `--overwrite` replaces them and `--out` writes the new run
elsewhere. The reason is that each mode reads the result of the one before
from the same folder, so a half replaced folder would be a mixture of two
runs.

### 9.12 Appended, replaced and refused

Four files grow a row at a time and are never rewritten: the instrument record
`STEM_instrument.csv`, the lattice record `lattice_KEY.csv`, the density
record `results/density.csv` and the phases record
`phases_KEY_record.csv`. Each is a history of every run of that command on
that sample, in the order they were made, and each carries its own method,
date and version, so a row read a year later says what made it. A file with
columns other than the ones this version writes is refused rather than
appended to, with a message saying to move it aside.

`cifs/cod/index.csv` is updated in place, row by row, keeping notes you have
added.

Everything else is replaced by a rerun, except the Rietveld results, which are
refused without `--overwrite`.

## 10. Known limitations

Eleven things the kit does not do yet, or does in a way worth knowing about
before you rely on it. Each says what the limitation is, when you meet it and
what to do meanwhile.

### 10.1 Space group coverage is partial

Indexing and refinement take a cell of any of the seven crystal systems, but
reflection conditions, the equivalence of reflections and their multiplicities
come from the symmetry operations of eight space groups only: Pm-3m, P4mm,
P4bm, P4/mbm, R3c, R3m, Pbnm and Amm2. Given any other symbol `xrdkit plot`
and `xrdkit lattice` print a note and index without conditions, so a peak may
be labelled with a reflection the group forbids and the labels are provisional
until checked by hand. With no space group at all no conditions are applied.
Trigonal and rhombohedral groups are handled on hexagonal axes only, so a cell
in the rhombohedral setting has to be converted first.

### 10.2 A recovered satellite carries a bias

A peak `xrdkit lattice` recovers from the K alpha 2 satellites is necessarily
a blend. To have been flagged at all, a good part of its height must come from
its parent's K alpha 2 line, and the refit models it as a lone doublet, so its
fitted position carries a bias. The report prints the number recovered so that
you can see how much of a refinement rests on such peaks; where that number is
more than a few, run again with `--no-satellites` and compare.

### 10.3 The Rietveld modes need a CIF beside a library entry

A library entry carries the sites, their kinds and what their Wyckoff
positions leave free, but no coordinates. A structure given as `library` alone
runs `xrdkit lebail` and then stops at the fixed atoms mode with a message
asking for `cif`. Find a CIF of the structure type, from the COD as Section 6
describes or from a licensed database, and name it beside the entry.

### 10.4 A CIF whose labels differ from the entry's must be mapped by hand

The coordinates and occupancies modes match each site the entry names to the
CIF atom of the same label and to nothing else. Where a CIF labels its atoms
otherwise, the `atoms` table has to name them site by site, as Section 8.11
describes, and the run stops naming the site when there is none.

### 10.5 The amount of an added element follows its host

An element the CIF lacks is placed only on the sites the `atoms` table names
it on, and its amount is split among them in proportion to multiplicity times
the host element's occupancy. There is no way to say how much of it goes on
each site. Where the real distribution is known to be different, the way round
it at present is to edit a CIF that already has the element where you want it
and name that CIF instead.

### 10.6 Preferred orientation has not been exercised

`--preferred-orientation H K L` adds a stage freeing a March-Dollase ratio
about that axis, and its flags are checked, but no refinement of a measured
scan has yet been run with it. Check any result it gives against a scan of the
same sample loaded to limit texture before relying on it.

### 10.7 One instrument, one histogram, one wavelength

The commands assume one instrument with a constant wavelength: one histogram,
with K alpha 1 and 2 or K alpha 1 alone, as the instrument parameter file
describes it. Time of flight and energy dispersive data, and several
histograms of one sample refined together, are not handled.

### 10.8 The instrument file cannot be made without GSAS-II

`xrdkit instrument` does the Caglioti width fit itself and needs nothing but
the scan for it, but the fit is not written out on its own: only the refined
file is, and the refinement is GSAS-II's. So there is no way to get an
instrument parameter file, even a starting one, without GSAS-II installed. If
you only want to look at the widths, the library function behind the command
returns them, and Part II says how to call it.

### 10.9 What xrdkit phases can and cannot find

Three limits meet here, and Section 6 shows all of them.

The search is for entries made of exactly the elements given. A doped
composition therefore finds nothing, and has to be searched on the elements of
its parent, with the dopants left out. This is a property of the search rather
than a defect, but it catches everyone once.

The figure of merit counts lines explained and lines missing without weighting
either by intensity, so a structure with many weak lines in the right places
can score well, and the rejection rule, which throws out any candidate whose
strongest line is absent, is doing most of the real work. Read the missing
count before the score.

The zero offset is not searched for. A displaced specimen shifts every peak by
more than the matching tolerance, and the ranking is then meaningless unless
`--zero` is given, which is why Section 6 takes the offset from the plot run
of Section 5.

### 10.10 The refined range comes only from the project file

`xrdkit lebail` and `xrdkit rietveld` take the two theta range from the
`[refine]` table or a sample's own, and there is no command line option for
it. Trying a different range means editing the project file, as Section 8.4
does, which is awkward when comparing two ranges on one sample. Use `--out` to
keep the two runs apart while you do it.

### 10.11 A text pattern carries no counting time

A two or three column `.xy` or `.xye` file records positions and intensities
and nothing else. `xrdkit check` therefore reports its time per step as
unknown rather than as a number, and the counting statistics have to be judged
from the intensities alone. Where the counting time matters to the record,
keep the instrument's own file beside the converted one.

### 10.12 Where to go next

Part I ends here. Everything above is done with the commands, and for most
work that is all that is needed.

Part II is the library reference: one section per module, with the functions a
command calls and the ones it does not, for anyone who wants to do something
the commands do not cover, to process a series in a loop, or to build a figure
of their own. It follows below.

## 11. Using the library

Part I is the whole of xrdkit as a set of commands, and for most work that is
all that is needed. Part II is the library underneath them: the modules the
commands call, the public functions and classes of each, and what every
argument does. Nothing here replaces a command. It is for the work that falls
outside one.

### 11.1 Why a script is worth writing

Four kinds of work ask for a script rather than a command.

A figure of your own. `xrdkit plot` writes one pattern and one hkl labelled
pattern in a fixed layout, and that is all it writes. The plotting functions
hand back the matplotlib figure and axes, so a script can set a title, choose
the window, write the refined cell into a corner, draw two patterns on one
axes, or anything else matplotlib can do, and then save. Section 15 does
exactly that.

A series in a loop. Twelve compositions checked, indexed and tabulated is
twelve commands and twelve results folders, or one script over a list of
sample keys that writes a single table at the end. The project file is
readable from the library as well, so the loop can take its samples from
`xrdkit.toml` rather than from a list typed out by hand, which Section 27
covers.

A stage sequence the commands do not offer. `xrdkit rietveld` runs three
modes in a fixed order, each with a fixed set of stages. A refinement that
needs a stage the modes do not have, or the same mode run twice from
different starting values, is built from a refinement job directly, which
Sections 28 and 29 cover.

Size and strain. A Scherrer size, a Williamson-Hall plot and the separation of
size broadening from strain broadening are in the library and behind no
command at all, because they rest on a judgement about which reflections to
use that a command cannot make for you. Section 23 is the whole of that.

### 11.2 Importing

Every name in Part II is imported from the module that defines it, as in
`from xrdkit.io import read_scan` or
`from xrdkit.peaks import exclude_kalpha2, find_peaks`. The module is named in
the heading of each section, so the import line follows from the section you
are reading.

The package re-exports most of these names at the top level as well, so
`from xrdkit import read_scan` reaches the same function and is what the older
guide used. Importing from the module is used throughout Part II because it
says where a function lives, which is what you want when you go looking for
the rest of its module. Two names are worth knowing at the top level:
`xrdkit.TTB_CELL`, the tetragonal tungsten bronze start cell of 12.45 and 3.94
angstrom that the examples use, and `xrdkit.list_entries`, which lists the
structure library entries of Section 25.

### 11.3 How the code blocks work

Every block of Python in Part II is part of a script. Save it in a file whose
name ends in `.py`, in the project folder or in a `scripts` folder inside it,
and run it from the project folder, the one holding `xrdkit.toml`, with
`py scripts/name.py` on Windows or `python3 scripts/name.py` on macOS. Nothing
here is meant to be typed at a Python prompt.

One line of prose above each block says which file it belongs to and what to
do with it: the whole of a file of that name, a new file to start, or lines to
add to the end of a file already started. A block that continues a script is
run with everything before it in the same file, never on its own. Each script
begins with its settings, marked by the comment line `Edit these lines for
each new sample. Nothing below needs changing.`, and those lines are the only
ones that change from one sample to the next: the scan file, the labels, and
the stem the output files are named from. A block introduced as what the
script printed is its output, which is what your own run should print in its
place.

The settings are this guide's examples and are meant to be replaced. Change
those lines and nothing else. In particular, do not use Replace All on a
number: `10.0` is the start of a two theta window in one line of a script and
part of something else further down, and replacing every one of them silently
changes the analysis.

### 11.4 Where the numbers come from

Every printed block in Part II came from running the script shown, in the same
project folder as Part I and on the same example files, the ones Section 1.4
lists. The two parts can therefore be read across. The cell `plot_pattern.py`
refines in Section 15 is the cell `xrdkit plot pellet_a` printed in Section
5.1, to the last figure, because the command and the script do the same
arithmetic on the same scan. The report `scan_quality.py` prints in Section 13
is the report `xrdkit check pellet_a` printed in Section 4.1. Where a number
differs from Part I, the text says which argument made it differ.

The scripts of Part II name their output files apart from Part I's, so that
running them leaves everything Part I wrote where it was.

### 11.5 What is in Part II

One section per module, in the order a workflow meets them. Section 12 is
`io`, the readers and the scan they return, and Section 13 is `quality`, the
assessment behind `xrdkit check`. Section 14 is `peaks`, the peak finder and
the K alpha 2 rule, and Section 15 is `plotting`, every figure the kit draws.
Sections 16 to 18 are the cell and what is done with it: `indexing` assigns
reflections to peaks and refines a cell as it goes, `cell` is the `Cell` class
and its geometry, and `lattice` is the least squares refinement behind
`xrdkit lattice`. Section 19 is `density`, formula masses and theoretical
densities, and Section 20 is `symmetry`, the space group operations,
reflection conditions and multiplicities that Section 10.1 says are known for
eight groups. Sections 21 to 23 are the widths: `broadening` fits peak shapes
and the Caglioti function, `instrument` turns a standard scan into
instrumental widths, and `sizestrain` takes those to a crystallite size and a
microstrain. Section 24 is `phases`, the search, the simulation and the
matching behind `xrdkit phases`, and Section 25 is `library`, the structure
library entries. Section 26 is `structure`, cell contents, bond lengths and
the site set up a refinement needs, and Section 27 is `project`, reading and
writing `xrdkit.toml`. Sections 28 and 29 are the refinement itself: `gsas2`
builds and runs a job under the GSAS-II Python, and `pipeline` is the mode
sequence behind `xrdkit lebail` and `xrdkit rietveld`. Section 30 lists the
modules you are not expected to call directly and says what calls them
instead, and Section 31 indexes every public name in the package with the
section that documents it.

## 12. io

`xrdkit.io` reads a scan file into an `XRDScan`, which is what every other
module in the kit takes. It has four public names: the dataclass, one reader
for each format, and one reader that chooses between them by the suffix.
Nothing in this module writes to a scan file.

### 12.1 XRDScan

`XRDScan` is a plain dataclass. The two arrays are the pattern and everything
else is what the file said about it.

| Field | What it holds |
| --- | --- |
| `two_theta` | positions in degrees, as a numpy array, rising |
| `intensity` | counts at each position, as a numpy array of the same length |
| `wavelength` | the K alpha 1 wavelength in angstroms, or `None` when the file carries none |
| `start_angle`, `end_angle` | the first and last two theta of the scan, in degrees |
| `step_size` | the step in degrees, from the range and the number of points, or `0.0` for a scan of one point |
| `time_per_step` | the counting time in seconds, or `None` when the file carries none |
| `sample_id` | whatever identifier the file carries, which for a text pattern is the file stem |
| `source_path` | the path it was read from, as a string, as it was given |
| `esd` | esds of the intensities when every row gave one, as an `.xye` does, and `None` otherwise |

Two of those are `None` more often than a reader expects. A two or three
column text pattern carries no wavelength and no counting time, so both come
back as `None` and have to be supplied from the project file or the command
line. Section 10.11 says what the missing time means for `xrdkit check`.

`sample_id` is what the diffractometer wrote into the file, not the sample key
of your project. The example scans keep the identifiers they were measured
under, so `pellet_a.xrdml` reports a `sample_id` of `10s`. Use the key for
naming results, and the identifier only for telling you what the instrument
thought it was measuring.

### 12.2 The readers

`read_xrdml(path)` reads a PANalytical `.xrdml` file. The two theta axis is
rebuilt from the start position, the end position and the number of
intensities rather than read point by point, so the step is constant by
construction. Intensities come from the `counts` element or from the
`intensities` element, whichever the writer used. The counting time comes from
the common counting time where the file gives one and from the first per point
time otherwise. It raises `ValueError` when the file holds no scan, no data
points, no intensities or no two theta axis.

`read_xy(path, wavelength=None)` reads a two or three column text pattern: two
theta, intensity, and in an `.xye` the esd of the intensity. Any number of
header or comment lines may come first, and the data begin at the first line
whose first character could begin a number. Fields are separated by spaces, by
tabs or by a comma, blank lines are ignored wherever they fall, either line
ending is read, and a byte order mark is stripped. `wavelength` in angstroms
is recorded on the scan, since the file carries none; `None` is allowed here,
and it is the command line that insists on one where the wavelength is needed.
It raises `ValueError` when the file holds no data, or a data line has fewer
than two numbers, a field that is not a number, or a two theta that does not
increase on the line before it, and the message names the file and the line.

`read_scan(path, wavelength=None)` chooses the reader by the suffix: `.xrdml`
goes to `read_xrdml`, and `.xy` and `.xye` go to `read_xy`. A `wavelength`
given here replaces the one an `.xrdml` carries and supplies the one a text
pattern lacks. Any other suffix raises `ValueError` naming the formats the kit
reads. This is the reader to call in a script that should accept whatever the
user has, and it is the one the commands call.

### 12.3 Reading a scan and looking at it

The whole of `scan_fields.py`:

```python
from xrdkit.io import read_xrdml

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"

scan = read_xrdml(SCAN_FILE)
print(f"sample id   {scan.sample_id}")
print(f"source      {scan.source_path}")
print(f"wavelength  {scan.wavelength} angstrom")
print(f"points      {scan.two_theta.size}")
print(f"range       {scan.start_angle:.2f} to {scan.end_angle:.2f} degrees")
print(f"step        {scan.step_size:.4f} degrees")
print(f"time        {scan.time_per_step:.1f} s per step")
print(f"esds        {scan.esd is not None}")
print(f"first point {scan.two_theta[0]:.4f} degrees, {scan.intensity[0]:.0f} counts")
print(f"last point  {scan.two_theta[-1]:.4f} degrees, {scan.intensity[-1]:.0f} counts")
```

It prints the fields a later step needs.

```text
sample id   10s
source      data\raw\pellet_a.xrdml
wavelength  1.540598 angstrom
points      4141
range       10.01 to 99.98 degrees
step        0.0217 degrees
time        34.2 s per step
esds        False
first point 10.0099 degrees, 652 counts
last point  99.9840 degrees, 913 counts
```

The range, the step, the point count and the counting time are the numbers
`xrdkit check pellet_a` reported in Section 4.1, because the command reads
them off this same object. The wavelength of 1.540598 angstrom is the one the
diffractometer recorded, and it is what every d spacing in Part II is
calculated from unless a script overrides it.

## 13. quality

`xrdkit.quality` is the whole of `xrdkit check`. It measures a scan against
the criteria of Section 4.2 and says which workflows the scan is within reach
of. Only the criteria that can be read off the scan itself are applied:
angular range, step size and counting statistics. Sample preparation,
standards and radiation are left to you.

### 13.1 The functions and the dataclasses

`assess_scan(scan)` takes an `XRDScan` and returns a `ScanQuality`. It raises
`ValueError` when the scan holds no intensities.

`format_report(quality)` takes that `ScanQuality` and returns the plain text
report as one string, with no trailing newline.

`ScanQuality` carries the measured numbers and one verdict per workflow.

| Field | What it holds |
| --- | --- |
| `start_angle`, `end_angle`, `step_size`, `points` | the scanned range in degrees, the step, and how many points |
| `time_per_step` | the counting time in seconds, or `None` |
| `maximum`, `median` | the strongest intensity and the median, in counts |
| `peak_over_median` | the first divided by the second, and infinite when the median is zero and the maximum is not |
| `high_angle_start` | where the high angle third begins, at the start plus two thirds of the range |
| `high_angle_maximum`, `high_angle_median` | the same two intensities over that third |
| `verdicts` | a dict of workflow name to `Verdict` |

`Verdict` has two fields: `suitable`, a boolean, and `reasons`, a list of
strings, one per criterion failed, each naming the value measured and the
threshold it fell short of. A suitable verdict has an empty list.

The median stands in for the background, because most of the points in a
powder pattern are background. A scan counts as covering an angular range when
it starts and ends within one step of the limits, and a step size passes when
it lies inside the range the table gives, limits included, with a tolerance of
a millionth of a degree so that a step computed as 0.020000000001 still counts
as 0.02.

Two module level names are worth knowing. `CRITERIA` is the thresholds
themselves, a dict of workflow name to a dict of criterion name to threshold,
and `WORKFLOWS` is the workflow names in the order the report prints them:
`plotting`, `phase_identification`, `le_bail`, `rietveld`. Reading `CRITERIA`
is how to put the numbers of Section 4.2 into a table of your own rather than
copying them out.

### 13.2 Assessing a scan

The whole of `scan_quality.py`:

```python
from xrdkit.io import read_scan
from xrdkit.quality import assess_scan, format_report

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"

quality = assess_scan(read_scan(SCAN_FILE))
print(format_report(quality))
```

It prints this.

```text
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

That is the report `xrdkit check pellet_a` printed in Section 4.1, line for
line, less the sample line the command adds from the project file and the path
of the file it wrote. Section 4.3 reads it. The command is `assess_scan` and
`format_report` with the project file around them, so a script that wants the
numbers rather than the report takes them off the `ScanQuality` instead of
parsing the text: `quality.peak_over_median` is the 14.0 of the second
verdict, and `quality.verdicts["le_bail"].reasons` is the two reasons of the
third, already separated.

## 14. peaks

`xrdkit.peaks` finds the reflections in a pattern and decides which of them
are K alpha 2 satellites. It has five public names: the `Peak` dataclass, the
finder, the flagging rule, the filter that applies it, and a CSV writer.

### 14.1 Peak

| Field | What it holds |
| --- | --- |
| `two_theta` | the position in degrees, refined off the grid by a parabola through the peak and its two neighbours |
| `intensity` | the height in counts at that point, background included |
| `prominence` | the height above the higher of the two saddles either side, in counts |
| `fwhm` | the full width at half the prominence, in degrees |
| `d_spacing` | from Bragg's law at the scan's wavelength, in angstroms |
| `relative_intensity` | the height as a percentage of the strongest peak in the list, so the strongest is always 100 |
| `kalpha2_of` | the index in the list of the K alpha 1 parent this peak is a satellite of, or `None` |
| `background` | the counts under the peak, the lowest intensity within one degree of it, or `0.0` when it was not estimated |

The position is refined rather than taken from the grid because the step is
0.0217 degrees on the example scans and the positions that matter are wanted
to about a tenth of that. The refinement falls back to the grid position when
the peak sits on the first or the last point, when the three points are
collinear, or when the vertex lands outside the neighbouring points.

### 14.2 find_peaks

`find_peaks(scan, min_prominence=0.02, min_distance=0.15,
two_theta_range=None, flag_satellites=True)` returns a list of `Peak` ordered
by two theta.

| Argument | What it does |
| --- | --- |
| `scan` | the `XRDScan` to search |
| `min_prominence` | the least prominence a peak must have, as a fraction of the strongest intensity inside the searched range. 0.02 by default, so a peak must stand 2 per cent of the tallest peak clear of its surroundings. Lower it to catch weak reflections along with more noise, raise it for a list of the major peaks alone |
| `min_distance` | the least separation between peaks, in degrees. 0.15 by default, which is wide enough to keep one reflection from being reported twice and narrow enough to keep a resolved K alpha 2 satellite as a peak of its own |
| `two_theta_range` | a `(low, high)` window in degrees to search inside. `None`, the default, searches the whole scan. Both the relative intensities and the prominence threshold are taken from the strongest intensity inside the window, so narrowing the window changes both |
| `flag_satellites` | whether to run `flag_kalpha2` over the result before returning it. `True` by default. `False` leaves every `kalpha2_of` as `None`, which is what to pass for radiation with no K alpha 2 |

It raises `ValueError` when `two_theta_range` selects no points, or when the
step cannot be determined because the scan has one point and no step. It
returns an empty list, rather than raising, when the window holds points but
no peak clears the threshold.

Satellites are returned rather than dropped. Flagging and filtering are
separate steps so that a peak list meant to be read can keep them with their
parents named, while a peak list meant to be indexed has them removed.

### 14.3 flag_kalpha2

`flag_kalpha2(peaks, position_tolerance=0.5, intensity_ratio=(0.2, 0.8),
wavelength_ratio=1.002486)` sets `kalpha2_of` on the peaks that look like
satellites, in place, and returns the same list for chaining.

A satellite sits where the parent's d spacing puts it at the longer K alpha 2
wavelength, always to high angle, and carries roughly half the parent's
intensity. A peak is flagged only when both hold. Its parent is the nearest
peak below it in two theta whose predicted satellite position falls within
`position_tolerance` times that peak's own full width at half maximum, and the
peak is then flagged only if the ratio of the two intensities lies inside
`intensity_ratio`. A peak whose nearest candidate parent fails the intensity
test is not handed on to a farther one.

| Argument | What it does |
| --- | --- |
| `peaks` | the list to flag, modified in place |
| `position_tolerance` | the largest gap allowed between a peak and its parent's predicted satellite position, as a fraction of the parent's full width at half maximum. 0.5 by default |
| `intensity_ratio` | the allowed `(low, high)` range of the peak's height above background over its parent's. `(0.2, 0.8)` by default, around the theoretical 0.5. Raising the ceiling towards 1 starts flagging genuine reflections that happen to sit at the K alpha 2 spacing above a neighbour of comparable height |
| `wavelength_ratio` | K alpha 2 over K alpha 1, 1.002486 for copper. `None` flags nothing, which is the setting for monochromated radiation |

It raises `ValueError` when `intensity_ratio` is not a rising pair of positive
numbers.

Both intensities are heights above the background, `intensity` less
`background`, and not the raw heights: raw heights on a high background give a
ratio pulled towards one. Prominences are no substitute either, because a
satellite on its parent's tail has its prominence measured down to the saddle
between the two, which pulls the ratio the other way.

Only resolved satellites can be caught this way. Below about 50 degrees the
pair is not separated enough for the finder to report two peaks, so there is
nothing to flag, which is what Section 5.3 means when it says the doublet is
not resolved there.

### 14.4 exclude_kalpha2 and peaks_to_csv

`exclude_kalpha2(peaks)` returns a new list holding the peaks whose
`kalpha2_of` is `None`. The survivors are copies rather than the same objects,
so the original list is untouched, and their `relative_intensity` is
recomputed against the strongest of them, so it again reaches 100. An empty
list comes back when nothing survives.

`peaks_to_csv(peaks, path)` writes the list as CSV with a header and returns
the path, creating the parent folder if it does not exist. The columns are
`two_theta`, `intensity`, `prominence`, `fwhm`, `d_spacing`,
`relative_intensity` and `kalpha2_of`. Positions and widths keep three decimal
places, d spacings four, counts one and relative intensities two, and
`kalpha2_of` is the parent's index or blank. `background` is not written: it is
an intermediate of the flagging rather than a measurement to report.

### 14.5 Finding the peaks of a scan

The whole of `peak_list.py`:

```python
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks, peaks_to_csv

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
STEM = "pellet_a"

HEADING = "  i  two_theta  intensity   fwhm  d_spacing  relative  kalpha2_of"


def show(peaks, indices):
    """Print the given peaks of ``peaks`` one to a line."""
    for i in indices:
        peak = peaks[i]
        print(
            f"{i:3d}  {peak.two_theta:9.3f}  {peak.intensity:9.1f}"
            f"  {peak.fwhm:5.3f}  {peak.d_spacing:9.4f}"
            f"  {peak.relative_intensity:8.2f}"
            f"  {'' if peak.kalpha2_of is None else peak.kalpha2_of:>10}"
        )


scan = read_scan(SCAN_FILE)
peaks = find_peaks(scan, min_prominence=0.02, two_theta_range=(10.0, 80.0))
print(f"{len(peaks)} peaks found")
print(HEADING)
show(peaks, range(4))

flagged = [i for i, peak in enumerate(peaks) if peak.kalpha2_of is not None]
print(f"{len(flagged)} flagged as K alpha 2 satellites")
print(HEADING)
show(peaks, [peaks[flagged[0]].kalpha2_of, flagged[0]])

clean = exclude_kalpha2(peaks)
print(f"{len(clean)} peaks after exclude_kalpha2")
print(HEADING)
show(clean, range(4))
print(peaks_to_csv(clean, f"results/library/peaks_{STEM}.csv"))
```

It prints the first four peaks, then the first flagged satellite with the
parent it was matched to, then the first four of what is left, then the file
it wrote.

```text
41 peaks found
  i  two_theta  intensity   fwhm  d_spacing  relative  kalpha2_of
  0     22.748     2383.0  0.196     3.9059     17.07            
  1     25.882     6968.0  0.119     3.4397     49.90            
  2     26.921     1397.0  0.188     3.3092     10.00            
  3     27.895     7179.0  0.140     3.1958     51.41            
6 flagged as K alpha 2 satellites
  i  two_theta  intensity   fwhm  d_spacing  relative  kalpha2_of
 17     51.920     2904.0  0.153     1.7597     20.80            
 18     52.061     2157.0  0.033     1.7553     15.45          17
35 peaks after exclude_kalpha2
  i  two_theta  intensity   fwhm  d_spacing  relative  kalpha2_of
  0     22.748     2383.0  0.196     3.9059     17.07            
  1     25.882     6968.0  0.119     3.4397     49.90            
  2     26.921     1397.0  0.188     3.3092     10.00            
  3     27.895     7179.0  0.140     3.1958     51.41            
results\library\peaks_pellet_a.csv
```

Peak 18 is the satellite of peak 17. It sits 0.141 degrees above it, which is
where the longer wavelength puts the same d spacing, and its height is about
three quarters of its parent's over a background the two share. All six
flagged peaks lie above 50 degrees, as Section 14.3 says they must.

The first four peaks read the same before and after the filter, relative
intensities included, and that is the expected result rather than a sign that
nothing happened. `exclude_kalpha2` recomputes the relative intensities
against the strongest survivor, and the strongest peak of this scan is not a
satellite, so the divisor does not change. Had the tallest peak in the list
been one, every relative intensity would have risen.

The window of 10 to 80 degrees is why 41 peaks are found here where Section
5.1 reports 49: the command searches the whole scan, which runs to 99.98
degrees, and the six satellites are still in its count.

## 15. plotting

`xrdkit.plotting` draws every figure the kit produces. There are eight public
names: one that sets the style, four that draw something, two that write text
onto axes already drawn, and one that saves. All of them work on matplotlib
figures and axes and hand them back, so anything matplotlib can do to a figure
can be done to one of these before it is saved. Nothing in this module calls
`pyplot`, and nor should a script that uses it: take the figure and the axes
that come back and work on those.

### 15.1 apply_style

`apply_style()` sets the global matplotlib rcParams to a clean single column
journal style: serif text, thin lines, ticks on all four sides pointing in, no
grid, and a figure size of 3.5 by 2.6 inches. Call it once, before anything is
drawn. It changes the rcParams of the process, so a script that draws figures
of its own alongside the kit's should call it first and then override what it
wants.

### 15.2 plot_pattern

`plot_pattern(scan, ax=None, scale="linear", normalise=False, label=None,
colour="black", linewidth=0.7)` draws one pattern and returns the figure and
the axes as a tuple.

| Argument | What it does |
| --- | --- |
| `scan` | the `XRDScan` to draw |
| `ax` | axes to draw on. `None`, the default, makes a new single column figure |
| `scale` | `"linear"`, `"sqrt"` or `"log"`. Linear shows the strong reflections in proportion; square root brings up the weak ones, which is what a phase check wants. Anything else raises `ValueError` |
| `normalise` | scale the trace to a maximum of 1, after the transform, so a normalised square root trace runs from 0 to 1 on the square root scale. Off by default for a single pattern |
| `label` | the legend label for the line. No legend is drawn unless you ask for one |
| `colour`, `linewidth` | the appearance of the line |

The axes come back with the two theta label, the intensity label and the tick
locators already set, and the limits as matplotlib chose them. Set the limits
yourself before annotating, because the annotation reads them as they stand.

### 15.3 plot_stacked

`plot_stacked(scans, labels=None, offset=None, scale="linear",
normalise=True, colour="black", linewidth=0.7, figsize=(3.5, 5.0))` draws
several patterns one above another and returns four things: the figure, the
axes, the vertical base of each slot, and the line drawn for each scan, the
last two in the order the scans were given.

| Argument | What it does |
| --- | --- |
| `scans` | the scans to stack, drawn bottom to top in the order given |
| `labels` | one label per scan, written at the upper right of each trace. Each scan's `sample_id` by default, which for the example files is the instrument's identifier rather than the project key, so give the labels |
| `offset` | the vertical spacing between traces. `None`, the default, gives 1.2 times the tallest scaled trace, which keeps the tallest peak of one pattern clear of the pattern above it. Raise it where labels crowd the trace above, lower it to fit more traces in |
| `scale`, `normalise`, `colour`, `linewidth` | as for `plot_pattern`, except that `normalise` is on by default, so that scans counted for different times can be compared |
| `figsize` | the figure size in inches, 3.5 by 5.0 by default, against the 3.5 by 2.6 of a single pattern |

It raises `ValueError` when `scans` is empty, when `scale` is unknown, or when
`labels` is a different length from `scans`.

The bases and the lines are what the annotation needs. A label placed against
`lines[-1]` rides on the top trace; a label placed at `bases[0]` plus enough
room to clear the tallest peak of the bottom trace sits in a row above it.

### 15.4 annotate_hkl

`annotate_hkl(ax, indexed, y=0.0, min_relative_intensity=5.0, fontsize=7,
rotation=90, min_separation=None, label_height=None, ambiguous="first",
max_levels=1, line=None)` writes an hkl label above the indexed peaks of one
trace and returns the labels written, in two theta order, as a list of
matplotlib `Text`.

| Argument | What it does |
| --- | --- |
| `ax` | the axes to write on |
| `indexed` | the list of `IndexedPeak` for that trace, from the indexing of Section 16. Peaks with no assignment are skipped; `mark_peaks` is for those |
| `y` | the base of the level 0 labels, in data coordinates. Ignored when `line` is given |
| `min_relative_intensity` | skip peaks below this percentage of the strongest peak of the scan. 5 by default; raise it to thin out a crowded figure |
| `fontsize` | the label size in points, 7 by default |
| `rotation` | the label angle in degrees anticlockwise, 90 by default, because upright labels take the least horizontal room |
| `min_separation` | how far apart in degrees two labels on the same level must be. `None`, the default, measures the rendered labels instead and keeps whichever do not overlap, which packs them as tightly as the text really allows |
| `label_height` | the height of one level in data coordinates, by default 0.04 of the y range of the axes |
| `ambiguous` | what to do with a peak that matched more than one reflection: `"first"` labels the assigned one, `"all"` joins every candidate with a solidus, `"skip"` leaves the peak unlabelled. Anything else raises `ValueError` |
| `max_levels` | how many rows of labels to try. 1 by default, which keeps them in a single row. Ignored when `line` is given, and below 1 raises `ValueError` |
| `line` | the trace the peaks belong to, one of the lines `plot_stacked` returned or `ax.lines[0]` after `plot_pattern`. Each label then rides on top of its own peak rather than sitting in a row. Giving one always measures, so `min_separation` has no effect |

Room is allotted by priority and not by position. The peaks are taken
strongest first, and a peak that finds no free level keeps no label at all, so
in a crowded stretch it is the weak reflections that lose theirs and the
strong ones a reader is looking for stay labelled. Labels are never stacked
into a column or shifted sideways onto a neighbour, which is what keeps every
label that is drawn pointing at the peak it belongs to.

The figure size and the x limits are read as they stand. Set both before
calling this: annotating and then resizing the figure or changing the limits
leaves the labels where the old geometry put them.

### 15.5 mark_peaks

`mark_peaks(ax, positions, y, marker="*", fontsize=9)` writes `marker` above
each two theta in `positions`, with its base at `y` in data coordinates, and
returns one `Text` per position in the order given. It is for the peaks a cell
does not account for, whether unindexed or from a second phase, which are
worth pointing at even though they carry no hkl. It places nothing and checks
nothing: every position given gets a marker, wherever it falls.

### 15.6 plot_caglioti and plot_rietveld

These two draw the figures of the later workflows, and the sections that own
those workflows show them in use. They are listed here because they live in
this module.

`plot_caglioti(fit, two_theta, fwhm, esd=None, included=None, ax=None,
two_theta_range=None, labels=("Caglioti fit", "Used in fit", "Excluded"))`
plots measured peak widths against two theta with a fitted Caglioti curve
through them. Widths the fit used are filled circles and the rest are open
ones, so the two are told apart without colour; `esd` draws error bars,
`included` says which widths the fit used, and `two_theta_range` sets the x
limits, which otherwise span the data widened out to whole major ticks. With
other `labels` the same figure serves for sample widths against an
instrumental curve. Section 22 uses it.

`plot_rietveld(pattern, reflections=None, title=None, phase_labels=None,
ax=None, sqrt_scale=False, difference_offset=None, result=None)` plots a
Rietveld fit and returns the figure and the axes. The observed points are
small open grey circles, the calculated pattern a black line through them and
the background a thin grey line; below the pattern each phase has a row of
tick marks at its reflections, first phase at the top, and below those the
difference is drawn about a zero of its own. `pattern` is the table the
GSAS-II driver exports, as the path of its CSV file or as anything indexed by
column name, and `reflections` is the reflection lists in the same form, a
mapping of phase name to list or a bare sequence. `sqrt_scale` takes the
square root of the counts, and the difference is then that of the square roots
so that it stays on the scale of the curves above it. `difference_offset` sets
the height of the difference zero, which by default goes below the lowest tick
row. `result` is the driver's refine result, as a dict or the path of its JSON
file, and writes Rwp, the goodness of fit and each phase's refined cell in a
small block below the legend. It raises `ValueError` when a table lacks a
column or when `phase_labels` does not have one name per reflection list.
Section 29 uses it.

The caption `plot_rietveld` writes takes its Rwp and goodness of fit from the
last stage that has no error entry, and a stage that was rejected has none: it
was rolled back, not failed. The cell in the same caption comes from the final
model. So when a run has rejected a stage, the two halves of that caption come
from different stages, and the residuals should be read from the stage table
instead, as Section 8.16 says.

### 15.7 save_figure

`save_figure(fig, path, formats=("png", "pdf"), dpi=300)` saves the figure
once per format beside `path` and returns the paths written, as a list. The
path is treated as a stem: an extension that already names one of the formats
is replaced, so `"pattern.png"` and `"pattern"` behave the same, while other
dots are kept, which matters for a name like `"x0.10_calcined"`. The parent
folder is created if it does not exist. The default pair is a png to look at
and a pdf to put in a document, at 300 dots per inch; raise `dpi` for a
journal that asks for 600, and pass a single format where only one is wanted.

### 15.8 A pattern with a figure of your own

`xrdkit plot` writes the pattern and the hkl labelled pattern in a fixed
layout. Everything up to that layout is the same work, so the script below
does what the command does and then keeps the axes instead of saving at once.
The customisation is modest and is the point of the exercise: a title of your
own, and the refined cell written into the corner of the axes, so that the
figure carries the numbers its labels were placed from. Both are ordinary
matplotlib on the `Axes` object that came back.

It also marks the peaks that a tighter tolerance leaves unexplained. The
indexing is run twice: once as the command runs it, which refines the cell,
and once more at a tolerance of 0.02 degrees against that refined cell, which
is a question about the fit rather than a way of getting a better one.
Section 16 is where `index_peaks` and `index_and_refine` are documented.

The whole of `plot_pattern.py`:

```python
from xrdkit.cell import Cell
from xrdkit.indexing import index_and_refine, index_peaks, indexing_summary
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks
from xrdkit.plotting import (
    annotate_hkl,
    apply_style,
    mark_peaks,
    plot_pattern,
    save_figure,
)

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
TITLE = "pellet_a, Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
STEM = "library_pellet_a"
START_CELL = Cell.tetragonal(12.45, 3.94)
SPACE_GROUP = "P4bm"
WINDOW = (10.0, 100.0)

apply_style()
scan = read_scan(SCAN_FILE)
peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=WINDOW))
indexed, fit = index_and_refine(
    peaks,
    start_cell=START_CELL,
    wavelength=scan.wavelength,
    space_group=SPACE_GROUP,
)
summary = indexing_summary(indexed)
print(f"a = {fit.cell.a:.4f}, c = {fit.cell.c:.4f} angstrom")
print(f"zero offset {fit.zero_offset:.3f} degrees")
print(f"{summary['n_indexed']} of {summary['n_peaks']} peaks indexed")
print(f"rms {summary['rms_difference']:.4f} degrees")

tight = index_peaks(
    peaks,
    fit.cell,
    scan.wavelength,
    tolerance=0.02,
    zero_offset=fit.zero_offset,
    space_group=SPACE_GROUP,
)
unexplained = [entry.peak.two_theta for entry in tight if not entry.is_indexed]
positions = [round(value, 3) for value in unexplained]
print(f"{len(unexplained)} peaks outside 0.02 degrees: {positions}")

fig, ax = plot_pattern(scan, scale="sqrt")
ax.set_xlim(*WINDOW)
bottom, top = ax.get_ylim()
highest = float(max(ax.lines[0].get_ydata()))
ax.set_ylim(bottom, highest + 0.30 * (top - bottom))

labels = annotate_hkl(ax, indexed, min_relative_intensity=5.0, line=ax.lines[0])
markers = mark_peaks(ax, unexplained, y=highest)
print(f"{len(labels)} labels, {len(markers)} markers")

# The customisation the plot command does not offer: a title of your own, and
# the refined cell written into the corner of the axes, so that the figure
# carries the numbers its labels were placed from.
ax.set_title(TITLE)
ax.text(
    0.99,
    0.03,
    f"a = {fit.cell.a:.4f} angstrom\nc = {fit.cell.c:.4f} angstrom\n"
    f"rms {fit.rms_two_theta:.4f} degrees",
    transform=ax.transAxes,
    horizontalalignment="right",
    verticalalignment="bottom",
    fontsize=5,
)
print(save_figure(fig, f"figures/pattern_{STEM}"))
```

It prints the refinement, then the peaks the tighter tolerance left, then what
went onto the figure and where the figure went.

```text
a = 12.4803, c = 3.9324 angstrom
zero offset 0.170 degrees
43 of 43 peaks indexed
rms 0.0099 degrees
2 peaks outside 0.02 degrees: [26.921, 94.318]
21 labels, 2 markers
[WindowsPath('figures/pattern_library_pellet_a.png'), WindowsPath('figures/pattern_library_pellet_a.pdf')]
```

The cell of 12.4803 and 3.9324 angstrom, the zero offset of 0.170 degrees and
the rms of 0.0099 degrees are what `xrdkit plot pellet_a` printed in Section
5.1, figure for figure, because this is the same arithmetic on the same scan
with the same start cell. What differs is the count: the command reported 45
of 49 peaks indexed and the script reports 43 of 43, because the script hands
the indexing the list `exclude_kalpha2` returned and the command counts the
full list, satellites included.

Twenty one of the forty three peaks carry a label. The rest were either below
5 per cent of the strongest peak or lost their label to a stronger neighbour,
which is the rule Section 5.4 describes. Raising `min_relative_intensity`
thins the figure further; there is no setting that labels everything, because
there is no room.

The two markers sit over 26.921 and 94.318 degrees, the peaks whose nearest
reflection is more than 0.02 degrees away once the cell is fixed. They are
indexed at the default tolerance and so are not a second phase on this
evidence. The one at 26.921 degrees is the weak peak of Section 14.5, at 10
per cent of the strongest, where a small error in position is easy to come by;
the one at 94.318 degrees is at the far end of the scan, where the same error
in d spacing shows up as a larger error in two theta.

### 15.9 Labelling the top trace of a stack

`xrdkit stack` draws the stack but writes no hkl labels on it, so a stack with
its top trace labelled is a script. In a stack the labels normally go on the
top trace only, which is enough to tell a reader what every trace below shows.
Pass the line of that trace as `line`, and index the peaks of that same scan
and not of another one: a label stands over an observed position, so it
belongs to the pattern it was found in.

The whole of `stack_patterns.py`:

```python
from xrdkit.cell import Cell
from xrdkit.indexing import index_and_refine
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks
from xrdkit.plotting import annotate_hkl, apply_style, plot_stacked, save_figure

# Edit these lines for each new comparison. Nothing below needs changing.
SCAN_FILES = ["data/raw/pellet_a.xrdml", "data/raw/pellet_b.xrdml"]
LABELS = ["x = 0.10", "x = 0.12"]
STEM = "library_pellet_a_pellet_b"
START_CELL = Cell.tetragonal(12.45, 3.94)
SPACE_GROUP = "P4bm"
WINDOW = (10.0, 100.0)

apply_style()
scans = [read_scan(path) for path in SCAN_FILES]
fig, ax, bases, lines = plot_stacked(scans, labels=LABELS, scale="sqrt")
ax.set_xlim(*WINDOW)
print(f"{len(lines)} traces, bases {[round(base, 3) for base in bases]}")

top_scan = scans[-1]
top_peaks = exclude_kalpha2(find_peaks(top_scan, two_theta_range=WINDOW))
top_indexed, top_fit = index_and_refine(
    top_peaks,
    start_cell=START_CELL,
    wavelength=top_scan.wavelength,
    space_group=SPACE_GROUP,
)
print(f"top trace: a = {top_fit.cell.a:.4f}, c = {top_fit.cell.c:.4f} angstrom")

top_labels = annotate_hkl(ax, top_indexed, line=lines[-1], min_relative_intensity=10.0)
print(f"{len(top_labels)} labels on the top trace")
print(save_figure(fig, f"figures/stack_{STEM}"))
```

It prints the stack, the cell of the top trace and the labels it carried.

```text
2 traces, bases [0.0, 1.2]
top trace: a = 12.4709, c = 3.9267 angstrom
20 labels on the top trace
[WindowsPath('figures/stack_library_pellet_a_pellet_b.png'), WindowsPath('figures/stack_library_pellet_a_pellet_b.pdf')]
```

The bases of 0.0 and 1.2 are the default spacing at work: every trace is
normalised to a maximum of 1, so 1.2 times the tallest scaled trace is 1.2,
and the second slot starts there. Passing `offset` replaces that number.

The cell of the top trace, `pellet_b`, comes out at 12.4709 and 3.9267
angstrom against `pellet_a`'s 12.4803 and 3.9324, which is the composition
difference the stack was drawn to show. Twenty labels are written rather than
the twenty one of Section 15.8, because `min_relative_intensity` is 10 here
rather than 5: a stack has less vertical room per trace, so fewer labels fit.

To put the labels in a row above a chosen trace instead of on it, pass that
trace's base from `bases` as `y`, raised by enough to clear its tallest peak,
and leave `line` out.

## 16. indexing

`xrdkit.indexing` assigns reflections of a cell to observed peaks and refines
the cell as it goes. Section 5.3 says what that is for and Section 7 runs it
from the command line; this section is the functions themselves. Eleven of the
module's thirteen public names are here: two dataclasses for what the indexing
produces, one for the cell fit and one for the zero offset search, six
functions and a CSV writer. The other two are covered elsewhere,
`TetragonalCell` in Section 17 and `TTB_CELL` in Section 11.2.

### 16.1 Reflection and IndexedPeak

`Reflection` is one calculated reflection of a cell at one wavelength.

| Field | What it holds |
| --- | --- |
| `h`, `k`, `l` | the Miller indices of the family's representative |
| `hkl` | those three as a tuple |
| `d_spacing` | from the cell's reciprocal metric, in angstroms |
| `two_theta` | from Bragg's law at the wavelength, in degrees |
| `multiplicity` | how many reflections the family holds under the Laue group, or 1 for one built by hand |

`IndexedPeak` is one observed peak with its assignment.

| Field | What it holds |
| --- | --- |
| `peak` | the `Peak` of Section 14.1, as it was given |
| `corrected_two_theta` | that peak's position with the zero offset subtracted, in degrees |
| `reflection` | the `Reflection` assigned to it, or `None` |
| `difference` | corrected less calculated, in degrees, or `nan` when nothing was assigned |
| `candidates` | every reflection within the tolerance, nearest first |
| `is_indexed` | whether a reflection was assigned |

`candidates` is the field to read before a cell is believed. A peak with one
candidate has been identified. A peak with three has been labelled with
whichever of the three is nearest, which is a different thing, and the label
may change the next time the cell moves. That is why every refinement in the
kit, here and in Section 18, is run on the peaks with exactly one candidate.

### 16.2 generate_reflections

`generate_reflections(cell, wavelength, two_theta_max, two_theta_min=0.0,
space_group=None)` returns the allowed reflections of `cell` inside a two
theta window, as a list of `Reflection`.

| Argument | What it does |
| --- | --- |
| `cell` | the `Cell` of Section 17 to calculate d spacings from |
| `wavelength` | the radiation wavelength in angstroms |
| `two_theta_max` | the upper limit of the window in degrees, above `two_theta_min` and below 180 |
| `two_theta_min` | the lower limit in degrees, 0 by default |
| `space_group` | a symbol whose systematic absences to remove, one of the eight of Section 20, or `None` for none |

How far h, k and l run is worked out from the cell and the window rather than
fixed. The smallest d spacing that can diffract inside the window sets the
limits, since h can reach a over that d spacing before falling below it, and
likewise k with b and l with c, so a large cell is enumerated further than a
small one without anyone having to choose a cut-off.

Reflections equivalent under the Laue group are merged into one entry,
labelled by the representative of Section 20.4 and carrying the family's
multiplicity, and a space group drops the ones it forbids. Without a space
group the Laue group of the holohedry of the cell's crystal system is used,
which merges the families but removes nothing.

The order is by two theta, and two reflections whose angles agree to within
`COINCIDENCE_TOLERANCE`, a module constant of 1e-9 degrees, count as
coincident and are ordered by hkl instead, compared as a tuple, lowest first.
Exact coincidences are common in a tetragonal cell, where (550) and (710) both
give an h squared plus k squared of 50, and the metric returns such a pair
about 1e-14 degrees apart because it sums their terms in a different order.
Without the tolerance the order of that pair, and so the label a peak between
them takes, would be settled by rounding.

It raises `ValueError` when the window is empty or falls outside 0 to 180
degrees, or the space group is not one of the eight.

### 16.3 index_peaks

`index_peaks(peaks, cell, wavelength, tolerance=DEFAULT_TOLERANCE,
zero_offset=DEFAULT_ZERO_OFFSET, space_group=None)` returns one `IndexedPeak`
per peak, in the order the peaks were given.

| Argument | What it does |
| --- | --- |
| `peaks` | the observed peaks to index, satellites already removed |
| `cell` | the cell to index against, which is not refined here |
| `wavelength` | the radiation wavelength in angstroms |
| `tolerance` | the largest difference allowed between a corrected and a calculated position, in degrees. `DEFAULT_TOLERANCE` is 0.05, a few times the width of a well resolved peak on a lab instrument |
| `zero_offset` | a zero point correction in degrees, subtracted from every observed position before the comparison. `DEFAULT_ZERO_OFFSET` is 0 |
| `space_group` | as for `generate_reflections` |

The calculated window is the corrected peak positions widened by the tolerance
at both ends, so a peak at either end of the list can still find a partner
just outside the observed range. Each peak then takes the nearest reflection
within the tolerance; distances within `COINCIDENCE_TOLERANCE` of each other
are ties, and a tie goes to the lowest hkl, so a peak sitting on two
coincident reflections always takes the same label. A peak with nothing inside
the tolerance keeps `reflection=None` and a difference of `nan`, and still
appears in the list. An empty list of peaks gives an empty list back.

Nothing is refined here and no zero offset is looked for. Both of those are
`index_and_refine`, which calls this function in a loop.

### 16.4 refine_cell and CellFit

`refine_cell(indexed, wavelength=None, start_cell=None)` refines a cell from
indexed peaks by linear least squares and returns a `CellFit`. This is the fit
the indexing cycles use internally. It reports no esds; the refinement that
does is `refine_lattice` in Section 18.

One over d squared is linear in the components of the reciprocal metric, and
the crystal system of `start_cell` decides which of them are free: one for a
cubic cell, two for tetragonal, hexagonal and trigonal, three for
orthorhombic, four for monoclinic and six for triclinic. Tetragonal is assumed
when no start cell is given. A component that no indexed peak carries any
information on, its column being all zero, is held at the value `start_cell`
gives it and named in `held`. The fitted reciprocal metric is inverted to the
direct metric and the six parameters are read off that, so no crystal system
needs a formula of its own.

The observed d spacings are recomputed from `corrected_two_theta` rather than
taken off the peaks, which is what carries a zero point correction already
applied into the fit.

`CellFit` carries the result.

| Field | What it holds |
| --- | --- |
| `cell` | the fitted `Cell` |
| `n_peaks` | how many indexed peaks it was fitted on |
| `rms_two_theta` | the root mean square difference in degrees between their corrected positions and the positions the fitted cell puts them at |
| `held` | the names of the parameters held at the start cell's values, empty when all were fitted |
| `zero_offset` | the offset the indexing ran with, in degrees; 0 from `refine_cell` on its own |
| `coarse_two_theta_max` | the upper limit of `index_and_refine`'s first cycle, in corrected degrees; `None` from `refine_cell` on its own |

Those last two fields are filled in by `index_and_refine` and left at their
defaults by `refine_cell`, which is worth knowing before a script reads a zero
offset off a fit that never searched for one.

It raises `ValueError` when no wavelength is given, when there are not more
indexed peaks than free components, when a component has to be held and no
start cell was given, when the peaks cannot separate the free components, and
when the fitted reciprocal metric is not positive definite, which usually
means the assignments were wrong.

### 16.5 index_and_refine

`index_and_refine(peaks, start_cell, wavelength, zero_offset=0.0,
coarse_tolerance=0.4, coarse_two_theta_max=None,
fine_tolerance=DEFAULT_TOLERANCE, space_group=None, n_cycles=2,
search_zero=True)` indexes and refines in cycles and returns two things: the
whole peak list indexed against the final cell, and the `CellFit` of the last
cycle.

| Argument | What it does |
| --- | --- |
| `peaks` | the observed peaks to index |
| `start_cell` | the cell the first cycle indexes against, which also sets the crystal system |
| `wavelength` | the radiation wavelength in angstroms |
| `zero_offset` | a zero point correction in degrees. Giving one turns the search off, since the offset is then already known |
| `coarse_tolerance` | the tolerance of the first cycle in degrees, 0.4 by default, wide enough for a cell that is still some way off |
| `coarse_two_theta_max` | the upper limit of the first cycle in corrected degrees, or `None` to choose it from the data |
| `fine_tolerance` | the tolerance of every later cycle and of the indexing returned, `DEFAULT_TOLERANCE` by default |
| `space_group` | as for `generate_reflections` |
| `n_cycles` | how many cycles to run, at least one. The first is the coarse one |
| `search_zero` | look for the zero offset before indexing. Only done when `zero_offset` is left at zero |

The first cycle indexes the low angle peaks alone at the coarse tolerance,
where a cell that is still some way off can be trusted to put reflections near
the right peaks, and every later cycle indexes the whole list at the fine
tolerance against the cell the cycle before gave. Each cycle refines on the
peaks that matched exactly one reflection, so an ambiguous peak never chooses
between two candidates on the strength of a cell that has not converged.

The coarse window is chosen from the data and not fixed. Limits of 35, 50 and
70 degrees are tried in turn, and the first that holds enough lone peaks is
used: enough being the free parameters of the start cell's crystal system plus
two, and never fewer than three. If none of the three does, the whole pattern
is used and the refinement says how many peaks it lacks. One fixed limit does
not suit every cell, because a large cell has plenty of reflections below 35
degrees while a small or pseudo-cubic one has few there and those crowd
together, so it can be left with too few lone peaks to refine on at all.

It raises `ValueError` when `n_cycles` is below one, when there are no peaks
or none falls below a `coarse_two_theta_max` that was given, and from
`refine_cell` when a cycle has too few peaks to refine on.

### 16.6 estimate_zero_offset and ZeroSearch

`estimate_zero_offset(peaks, cell, wavelength, search=(-0.4, 0.4), step=0.01,
tolerance=DEFAULT_TOLERANCE, two_theta_max=None, space_group=None)` is the
search `index_and_refine` runs first, and it returns a `ZeroSearch` with the
best `offset` in degrees, the `n_indexed` peaks that offset left with exactly
one candidate, the `rms` of their differences, and the whole profile tried as
`offsets` and `counts`.

A specimen sitting proud of its holder moves every reflection by close to a
constant, which no cell can absorb, so a pattern indexed against an
uncorrected cell loses most of its peaks. Every offset across `search` is
tried at `step`, and the one leaving the most peaks with exactly one candidate
wins; a tie goes to the offset whose matched peaks sit closest to their
calculated positions.

Every trial gets a cell refined for it, from the low angle peaks, rather than
sharing one. That matters. A constant offset and a cell that is a per cent out
shift the low angle peaks by much the same amount, so a search holding one
cell returns whichever offset best hides the error in that cell rather than
the one that is right.

It raises `ValueError` when `search` is not a rising pair, when `step` is not
positive, or when `two_theta_max` leaves no peak to search on.

### 16.7 indexing_summary and indexed_to_csv

`indexing_summary(indexed)` returns a dict of five numbers: `n_peaks`,
`n_indexed`, `n_unindexed`, `n_ambiguous`, the peaks with more than one
candidate within the tolerance, and `rms_difference`, the root mean square of
the differences of the indexed peaks in degrees, which is `nan` when nothing
was indexed.

`indexed_to_csv(indexed, path)` writes the list as CSV with a header and
returns the path, creating the parent folder if it does not exist. The columns
are `two_theta`, `corrected_two_theta`, `d_spacing`, `relative_intensity`,
`h`, `k`, `l`, `calculated_two_theta`, `difference` and `n_candidates`. A peak
with no assignment leaves the hkl, the calculated position and the difference
blank but keeps its row, so the file is the whole peak list and not only the
part of it that worked.

### 16.8 The reflections of a cell, and the peaks they index

The script below takes the peak list of Section 14.5, the same scan through
the same finder with the same window, and does both halves of the module on
it: the reflections the start cell of Section 2.2 calculates, and then the
indexing and the refinement those reflections feed.

The whole of `reflections.py`:

```python
from xrdkit.cell import Cell
from xrdkit.indexing import (
    generate_reflections,
    index_and_refine,
    index_peaks,
    indexed_to_csv,
    indexing_summary,
    refine_cell,
)
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
STEM = "pellet_a"
START_CELL = Cell.tetragonal(12.45, 3.94)
SPACE_GROUP = "P4bm"
WINDOW = (10.0, 80.0)
REFLECTION_LIMIT = 30.0

scan = read_scan(SCAN_FILE)
every = generate_reflections(START_CELL, scan.wavelength, REFLECTION_LIMIT)
allowed = generate_reflections(
    START_CELL, scan.wavelength, REFLECTION_LIMIT, space_group=SPACE_GROUP
)
allowed_hkl = {reflection.hkl for reflection in allowed}
absent = sorted({reflection.hkl for reflection in every} - allowed_hkl)
print(f"{len(every)} reflections below {REFLECTION_LIMIT:.0f} degrees")
print(f"{len(allowed)} of them allowed in {SPACE_GROUP}, absent {absent}")
print("    h   k   l  d_spacing  two_theta  multiplicity")
for reflection in allowed[:8]:
    print(
        f"  {reflection.h:3d} {reflection.k:3d} {reflection.l:3d}"
        f"  {reflection.d_spacing:9.4f}  {reflection.two_theta:9.3f}"
        f"  {reflection.multiplicity:12d}"
    )

peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=WINDOW))
plain = index_peaks(peaks, START_CELL, scan.wavelength, space_group=SPACE_GROUP)
matched = sum(1 for entry in plain if entry.is_indexed)
print(f"{matched} of {len(peaks)} peaks indexed against the start cell, zero offset 0")

indexed, fit = index_and_refine(
    peaks,
    start_cell=START_CELL,
    wavelength=scan.wavelength,
    space_group=SPACE_GROUP,
)
summary = indexing_summary(indexed)
print(f"coarse window to {fit.coarse_two_theta_max:.2f} degrees")
print(f"zero offset {fit.zero_offset:.3f} degrees")
print(f"a = {fit.cell.a:.4f}, c = {fit.cell.c:.4f} angstrom, held {fit.held}")
print(
    f"{summary['n_indexed']} of {summary['n_peaks']} peaks indexed, "
    f"{summary['n_ambiguous']} with more than one candidate"
)
print(f"rms {summary['rms_difference']:.4f} degrees")

ambiguous = [entry for entry in indexed if len(entry.candidates) > 1]
first = ambiguous[0]
print(
    f"first ambiguous peak at {first.peak.two_theta:.3f} degrees, "
    f"corrected {first.corrected_two_theta:.3f}"
)
for candidate in first.candidates:
    print(f"  {candidate.hkl} at {candidate.two_theta:.3f} degrees")

lone = [entry for entry in indexed if len(entry.candidates) == 1]
again = refine_cell(lone, scan.wavelength, START_CELL)
print(
    f"refine_cell on {again.n_peaks} lone peaks: "
    f"a = {again.cell.a:.4f}, c = {again.cell.c:.4f} angstrom"
)
print(
    f"rms {again.rms_two_theta:.4f} degrees, "
    f"zero_offset {again.zero_offset:.1f}, "
    f"coarse_two_theta_max {again.coarse_two_theta_max}"
)
print(indexed_to_csv(indexed, f"results/library/indexed_{STEM}.csv"))
```

It prints the head of the reflection list, then what the start cell alone can
index, then the refinement, then a peak with more than one candidate, then the
same cell again from `refine_cell` on its own, then the file it wrote.

```text
15 reflections below 30 degrees
12 of them allowed in P4bm, absent [(1, 0, 0), (1, 0, 1), (3, 0, 0)]
    h   k   l  d_spacing  two_theta  multiplicity
    1   1   0     8.8035     10.040             4
    2   0   0     6.2250     14.216             4
    2   1   0     5.5678     15.905             8
    2   2   0     4.4017     20.157             4
    0   0   1     3.9400     22.549             2
    3   1   0     3.9370     22.566             8
    1   1   1     3.5963     24.737             8
    3   2   0     3.4530     25.780             8
8 of 35 peaks indexed against the start cell, zero offset 0
coarse window to 35.00 degrees
zero offset 0.170 degrees
a = 12.4799, c = 3.9323 angstrom, held ()
35 of 35 peaks indexed, 7 with more than one candidate
rms 0.0076 degrees
first ambiguous peak at 43.644 degrees, corrected 43.474
  (6, 0, 0) at 43.473 degrees
  (5, 1, 1) at 43.519 degrees
refine_cell on 28 lone peaks: a = 12.4799, c = 3.9323 angstrom
rms 0.0076 degrees, zero_offset 0.0, coarse_two_theta_max None
results\library\indexed_pellet_a.csv
```

Read the reflection list first. Fifteen families reach 30 degrees for this
cell and twelve of them are allowed in P4bm; the three that are not are (100),
(101) and (300), which Section 20.3 gets from the operations of the group. The
multiplicities are the families and not the reflections: (001) stands for two
reflections and (211) for sixteen, and the list holds one row per family
rather than one row per reflection.

The two rows at 22.549 and 22.566 degrees are why the tolerances of this
module are what they are. (001) and (310) lie 0.017 degrees apart, which is
less than the step of the scan, so no peak finder will separate them and any
peak there has two candidates whatever the tolerance.

Next the indexing. Against the start cell with no zero offset, 8 of the 35
peaks find a reflection inside the default tolerance of 0.05 degrees. That is
not a bad cell; it is an uncorrected one. `index_and_refine` searches for the
offset first, finds 0.170 degrees, chooses a coarse window of 35 degrees from
the data, and comes back with all 35 peaks indexed at an rms of 0.0076
degrees. The 0.170 degrees is `estimate_zero_offset`'s answer, the same number
Section 5.1 reported for this scan and the same one Section 7.4 calls a false
zero, for the reason Section 18 returns to.

Seven of the thirty five peaks carry more than one candidate. The first, at
43.644 degrees, could be (600) or (511), which sit 0.046 degrees apart: the
corrected position is 0.001 degrees from the first and 0.045 from the second,
so (600) is the label, but the peak is not evidence for it. Those seven take
no part in the refinement, which is why `refine_cell` on the twenty eight lone
peaks returns the same cell and the same rms. It is the same arithmetic: the
last cycle of `index_and_refine` is exactly that call. Its `CellFit` comes
back with `zero_offset` at 0 and `coarse_two_theta_max` at `None`, because
`refine_cell` on its own neither searched nor windowed.

The cell of 12.4799 and 3.9323 angstrom is all but the cell of Section 15.8,
which came to 12.4803 and 3.9324. The arguments are the same and the scan is
the same; what differs is the window, 10 to 80 degrees here against 10 to 100
there, which adds eight peaks at the top of the pattern and moves a by four
ten thousandths of an angstrom. Neither is the cell of Section 7.1, and
Section 18 says why.

## 17. cell

`xrdkit.cell` holds one public name, `Cell`, and every cell anywhere in the
kit is one. It carries all six parameters and a crystal system, in angstroms
and degrees, and it is frozen, so two cells of the same parameters and system
are equal and neither can be altered in place.

Every quantity a cell can be asked for comes from its metric tensor G, whose
entries are the dot products of the three axes. The volume is the square root
of its determinant and a d spacing comes from its inverse, one over d squared
being h dotted into G star h. There is therefore no formula per crystal
system in this module and none in this section: a monoclinic cell and a
tetragonal one are the same three lines of arithmetic on different numbers.

`TetragonalCell(a, c)`, in `xrdkit.indexing`, is the older name for a
tetragonal cell and returns `Cell.tetragonal(a, c)`; it is kept so that
existing scripts keep working, and new ones should use `Cell.tetragonal`.

### 17.1 The parameters, the constructors and what is refused

The six parameters are `a`, `b`, `c`, `alpha`, `beta` and `gamma`, and
`crystal_system` is one of the seven. A cell always holds all six. The system
decides which of them are free and the rest follow: b equals a in a tetragonal
cell, every angle is 90 in an orthorhombic one, gamma is 120 in a hexagonal
one. Monoclinic cells take b as the unique axis, so alpha and gamma are 90 and
beta is free. Hexagonal and trigonal cells are on hexagonal axes, so a cell in
the rhombohedral setting has to be converted to hexagonal axes before it can
be one of these.

There is a named constructor per system, each taking the free parameters in
the order above: `Cell.cubic(a)`, `Cell.tetragonal(a, c)`,
`Cell.orthorhombic(a, b, c)`, `Cell.hexagonal(a, c)`, `Cell.trigonal(a, c)`,
`Cell.monoclinic(a, b, c, beta)` and `Cell.triclinic(a, b, c, alpha, beta,
gamma)`. Each fills in the dependent parameters and hands the result through
the same checks as the constructor of the class itself, which takes all six
and the system.

Four things are refused, each with a message naming the parameter. A length or
an angle that is not a number or is not above zero, an angle that is not below
180, a parameter that disagrees with what the crystal system fixes by more
than a millionth, and a crystal system that is not one of the seven. Then the
angles are checked together: three angles that individually look reasonable
can still describe no cell at all, and a set whose metric is not positive
definite is refused as well.

Start a new file named `cell_geometry.py`.

```python
from xrdkit.cell import Cell

# Edit these lines for each new cell. Nothing below needs changing.
A, C = 12.4740, 3.9305
HKL = (3, 1, 0)

REFUSED = (
    ("a zero length", lambda: Cell.tetragonal(A, 0.0)),
    ("an angle of 180", lambda: Cell.monoclinic(A, A, C, 180.0)),
    ("b set apart from a", lambda: Cell(A, 12.0, C, 90.0, 90.0, 90.0, "tetragonal")),
    ("a crystal system of its own", lambda: Cell.from_parameters("rhombic", {"a": A})),
    (
        "angles that make no cell",
        lambda: Cell.triclinic(5.0, 5.0, 5.0, 20.0, 20.0, 150.0),
    ),
)

cell = Cell.tetragonal(A, C)
print(cell)
print(f"crystal system {cell.crystal_system}")
print(f"a {cell.a}, b {cell.b}, c {cell.c}")
print(f"alpha {cell.alpha}, beta {cell.beta}, gamma {cell.gamma}")
for description, build in REFUSED:
    try:
        build()
    except ValueError as error:
        print(f"{description}: {error}")
```

It prints the cell, its parameters and the five refusals.

```text
Cell(a=12.474, b=12.474, c=3.9305, alpha=90.0, beta=90.0, gamma=90.0, crystal_system='tetragonal')
crystal system tetragonal
a 12.474, b 12.474, c 3.9305
alpha 90.0, beta 90.0, gamma 90.0
a zero length: cell parameter c must be greater than 0, not 0.0
an angle of 180: cell parameter beta must be below 180 degrees, not 180.0
b set apart from a: cell parameter b must be 12.474 in a tetragonal cell, not 12
a crystal system of its own: unknown crystal_system 'rhombic'; the crystal systems are cubic, tetragonal, orthorhombic, hexagonal, trigonal, monoclinic, triclinic
angles that make no cell: cell parameters alpha, beta and gamma (20, 20, 150) make no cell: the metric is not positive definite
```

The dependent parameters are filled in and not merely tolerated: `b` comes
back as 12.474 without being asked for. The tolerance of a millionth is what
lets a value computed rather than typed, a b of 12.473999999999998 from a
refinement, pass as equal to a, while a b of 12.0 beside an a of 12.474 is
refused as what it is, a cell that is not tetragonal. The last refusal is the
one that cannot be caught parameter by parameter: 5, 5 and 5 angstrom with
angles of 20, 20 and 150 degrees has nothing wrong with any one of its six
numbers and is not a cell.

### 17.2 from_parameters, from_entry and CELL_PARAMETERS

`Cell.from_parameters(crystal_system, parameters)` builds a cell from a
mapping. It is the constructor for a cell whose system is not known until the
program runs, which is what a refinement reading a project file or a library
entry is doing. `parameters` gives either exactly the parameters the system
leaves free or all six, and all six are checked against the system as usual.

Which parameters a system leaves free is `CELL_PARAMETERS`, a dict in
`xrdkit.library` of crystal system to a tuple of names, and it is the order
everything in the kit uses: the free vector of Section 18, the esds of a fit,
`Cell.parameters`. Reading it is how to write a loop over the systems rather
than a chain of conditions, and `cell.parameter_names` is the entry for a
cell's own system.

`Cell.from_entry(entry, parameters)` is the same thing with the crystal system
taken from a structure library entry, the `StructureEntry` of Section 25. An
entry names which parameters its cell leaves free but carries no values for
them, since the values are a property of the sample and not of the structure
type, so they come separately.

Both raise `ValueError` naming the key when a parameter is unknown, is missing
or is not free for the system, and both raise on an unknown crystal system or
a value out of range.

Continues `cell_geometry.py`. Add these lines at the end of the file.

```python
from xrdkit.library import CELL_PARAMETERS, load_entry

REFUSED_PARAMETERS = (
    ("one the system fixes", {"a": A, "b": A, "c": C}),
    ("one left out", {"a": A}),
    ("one that is not a parameter", {"a": A, "d": C}),
)

for system, names in CELL_PARAMETERS.items():
    print(f"{system:14s} {', '.join(names)}")
hexagonal = Cell.from_parameters("hexagonal", {"a": 5.0, "c": 13.0})
print(f"hexagonal a {hexagonal.a}, b {hexagonal.b}, gamma {hexagonal.gamma}")
six = Cell.from_parameters(
    "tetragonal",
    {"a": A, "b": A, "c": C, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
)
print(f"the same cell from all six: {six == cell}")
entry = load_entry("ttb/P4bm")
print(f"entry {entry.name}, {entry.crystal_system}, {entry.space_group}")
print(f"its cell parameters {entry.cell_parameters}")
from_entry = Cell.from_entry(entry, {"a": A, "c": C})
print(f"the same cell from the entry: {from_entry == cell}")
for description, parameters in REFUSED_PARAMETERS:
    try:
        Cell.from_parameters("tetragonal", parameters)
    except ValueError as error:
        print(f"{description}: {error}")
```

It prints the free parameters of the seven systems, three ways of reaching the
same cell, and the three refusals.

```text
cubic          a
tetragonal     a, c
orthorhombic   a, b, c
hexagonal      a, c
trigonal       a, c
monoclinic     a, b, c, beta
triclinic      a, b, c, alpha, beta, gamma
hexagonal a 5.0, b 5.0, gamma 120.0
the same cell from all six: True
entry ttb/P4bm, tetragonal, P4bm
its cell parameters ('a', 'c')
the same cell from the entry: True
one the system fixes: cell parameter 'b' is not free in a tetragonal cell; give a, c or all six
one left out: missing cell parameter 'c'; a tetragonal cell needs a, c
one that is not a parameter: unknown cell parameter 'd'; the cell parameters are a, b, c, alpha, beta, gamma
```

The entry `ttb/P4bm` is the one the `[structures.ttb_p4bm]` table of Section
2.2 names, and `from_entry` builds the same cell that `Cell.tetragonal` did,
because the entry's crystal system is tetragonal and its cell parameters are a
and c. Giving all six works as well and gives an equal cell, which is what a
script reading six numbers out of a CIF wants.

The three refusals are worth reading as a group, because they are the three
ways a mapping can be wrong and the messages tell them apart. Giving `b` is
giving a parameter the system fixes, and the message says what to give
instead. Leaving `c` out is giving too few. Giving `d` is not a cell parameter
at all.

### 17.3 The geometry: volume, d spacings and the metric

`cell.metric_tensor` is G, so a fractional vector v has a length of the square
root of v dotted into G v. It is a plain property, computed when asked for.
`cell.reciprocal_metric` is G star, the inverse of G, cached on first use and
returned read only, since everything that indexes a cell asks for it over and
over.

`cell.volume` is the square root of the determinant of G, in cubic angstroms.

`cell.d_spacing(h, k, l)` is the d spacing of one reflection in angstroms,
from one over d squared being h dotted into G star h. It raises `ValueError`
for (000), which has no d spacing. `cell.d_spacings(hkl)` is the same for an
array of shape (n, 3) and returns an array of n, which is the form the
indexing of Section 16 and the refinement of Section 18 use, both of them
asking for hundreds of d spacings at a time. It raises `ValueError` when the
array is not of that shape or a row is (000).

Continues `cell_geometry.py`. Add these lines at the end of the file.

```python
import numpy as np


def show_matrix(matrix, decimals):
    """Print a 3 by 3 matrix, rounded so that no minus zero is printed."""
    for row in matrix:
        values = [round(value, decimals) + 0.0 for value in row]
        width = decimals + 4
        print("  " + "  ".join(f"{value:{width}.{decimals}f}" for value in values))


print(f"V = {cell.volume:.3f} cubic angstrom")
print(f"sqrt(det G) = {np.sqrt(np.linalg.det(cell.metric_tensor)):.3f}")
print(f"d{HKL} = {cell.d_spacing(*HKL):.4f} angstrom")
spacings = cell.d_spacings(np.array([(1, 1, 0), (0, 0, 1), HKL]))
print("d_spacings " + "  ".join(f"{value:.4f}" for value in spacings))
print("metric tensor G")
show_matrix(cell.metric_tensor, 6)
print("reciprocal metric G*")
show_matrix(cell.reciprocal_metric, 8)
identity = cell.metric_tensor @ cell.reciprocal_metric
print(f"G G* is the identity: {np.allclose(identity, np.eye(3))}")
try:
    cell.d_spacing(0, 0, 0)
except ValueError as error:
    print(f"(000): {error}")
```

It prints the volume two ways, one d spacing and three, the two metric
tensors, and the refusal.

```text
V = 611.588 cubic angstrom
sqrt(det G) = 611.588
d(3, 1, 0) = 3.9446 angstrom
d_spacings 8.8204  3.9305  3.9446
metric tensor G
  155.600676    0.000000    0.000000
    0.000000  155.600676    0.000000
    0.000000    0.000000   15.448830
reciprocal metric G*
    0.00642671    0.00000000    0.00000000
    0.00000000    0.00642671    0.00000000
    0.00000000    0.00000000    0.06472982
G G* is the identity: True
(000): (000) has no d spacing
```

The two volumes agree because they are the same calculation: `cell.volume` is
the square root of that determinant. In a tetragonal cell G is diagonal, its
first two entries a squared and its third c squared, and G star is diagonal
with the reciprocals, which is the whole of the tetragonal d spacing formula
without anyone writing it down. Read across the diagonal of G star and the
first entry is one over 155.60, the second the same and the third one over
15.45; a reflection's one over d squared is h squared times the first plus k
squared times the second plus l squared times the third. A monoclinic cell
fills in one off-diagonal pair and a triclinic one fills in all three, and
nothing else in the arithmetic changes.

The off-diagonal entries print as zeros here after the rounding the script
does. They are not exactly zero in floating point, being of the order of 1e-15
in G and 1e-19 in G star, which is the inversion of a matrix whose diagonal
spans an order of magnitude and is far below anything that matters.

### 17.4 parameters, replace and to_dict

Four names read a cell out again.

`cell.parameter_names` is the tuple of names the crystal system leaves free,
which is `CELL_PARAMETERS` for that system. `cell.parameters` is those names
with their values, as a dict in the same order.

`cell.replace(**changes)` returns a new cell of the same crystal system with
the free parameters named in `changes` set to new values, the dependent ones
following. A cell is frozen, so this is how a cell is moved: a refinement step
that widens a by a thousandth, or the central differences the volume esd of
Section 18 is propagated through. It raises `ValueError` naming the parameter
when one is not free for the system, which is the check that stops a script
setting `b` on a tetragonal cell and silently getting something that is not a
cell of that system.

`cell.to_dict()` is all six parameters by name with `crystal_system` beside
them, which is the form for a CSV row, a JSON file or a project file table.

Continues `cell_geometry.py`. Add these lines at the end of the file.

```python
print(f"parameter_names {cell.parameter_names}")
print(f"parameters {cell.parameters}")
wider = cell.replace(a=12.50)
print(f"replace(a=12.50) gives a {wider.a}, b {wider.b}, c {wider.c}")
print(f"its volume {wider.volume:.3f} against {cell.volume:.3f} cubic angstrom")
print(f"to_dict {cell.to_dict()}")
try:
    cell.replace(b=12.50)
except ValueError as error:
    print(f"replace(b=12.50): {error}")
```

It prints the free parameters, a cell moved, and the dict.

```text
parameter_names ('a', 'c')
parameters {'a': 12.474, 'c': 3.9305}
replace(a=12.50) gives a 12.5, b 12.5, c 3.9305
its volume 614.141 against 611.588 cubic angstrom
to_dict {'a': 12.474, 'b': 12.474, 'c': 3.9305, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0, 'crystal_system': 'tetragonal'}
replace(b=12.50): cell parameter 'b' is not free in a tetragonal cell; the free parameters are a, c
```

`replace(a=12.50)` moved b with a, because b is not free to be left behind,
and the volume followed. `replace(b=12.50)` is refused for the same reason.
Note the difference between `parameters`, which gives the two free numbers and
is what a fit reports, and `to_dict`, which gives all six and the system and
is what a file records.

## 18. lattice

`xrdkit.lattice` is the refinement behind `xrdkit lattice`. It refines the
free parameters of a cell together with the systematic errors that move every
peak of a pattern, and unlike `refine_cell` in Section 16 it reports an
estimated standard deviation for each. There are three public names: the
refinement, the dataclass it returns, and a flattener for a CSV row.

Everything the command does around this function is Section 7 and is not
repeated here: choosing the start cell and the space group from the project
file, refitting each peak as a K alpha 1 line of a doublet, recovering flagged
satellites, and freeing the displacement for a pellet and the zero for a
powder. A script that wants that work done should run the command.

### 18.1 refine_lattice

`refine_lattice(indexed, wavelength, start_cell, fit_zero=True,
fit_displacement=False, radius_mm=None, start_zero=0.0)` refines by least
squares and returns a `LatticeFit`.

| Argument | What it does |
| --- | --- |
| `indexed` | the indexed peaks, from `index_peaks` or the first element of `index_and_refine` |
| `wavelength` | the radiation wavelength in angstroms |
| `start_cell` | the cell to start from, which sets the crystal system and so which parameters are free |
| `fit_zero` | refine the zero point error. `False` holds it at `start_zero`, which is the thing to do when the instrument zero is known from a standard |
| `fit_displacement` | refine the specimen displacement. Held at zero when `False` |
| `radius_mm` | the goniometer radius in millimetres, needed only to refine a displacement |
| `start_zero` | the zero point error to start from in degrees: refined from there when `fit_zero`, held there when not |

The free vector is the cell parameters of the crystal system, in the order of
`start_cell.parameter_names`, then the zero and then the displacement, and the
last two appear only when they are refined. Each peak is modelled as its Bragg
position for the trial cell, plus the zero, plus a displacement term of minus
two s cos(theta) over R radians, the usual flat plate term: a specimen below
the focusing circle moves every peak to low angle, most of all at low angle.

The observed positions are `peak.two_theta` and not `corrected_two_theta`,
because the zero point is what is being refined rather than something already
assumed. An indexing run with a zero offset is still the right input; its
assignments are used and its correction is not.

Only peaks that matched exactly one reflection are used, so an assignment that
was a choice between rivals cannot pull the cell. The fit needs at least three
more such peaks than it has parameters, and raises `ValueError` when it has
fewer, or when a displacement is asked for without a radius.

A zero point error and a specimen displacement are not independent. One is a
constant and the other follows the cosine of theta, and over a short angular
range the two are very nearly the same parameter, so a fit that frees both on
a narrow pattern reports small esds for two numbers trading against each
other. Section 7.3 says why the choice between them belongs to the mounting
and not to the data.

### 18.2 LatticeFit

| Field | What it holds |
| --- | --- |
| `cell` | the refined `Cell` |
| `a`, `b`, `c`, `alpha`, `beta`, `gamma` | all six parameters of that cell, in angstroms and degrees |
| `esd_a` to `esd_gamma` | the esd of each, or `None` for one that was not refined, whether held or fixed by the crystal system |
| `volume`, `esd_volume` | the cell volume in cubic angstroms and its esd |
| `zero`, `esd_zero` | the zero point error in degrees, added to every calculated position, and its esd or `None` |
| `displacement`, `esd_displacement` | the specimen displacement in millimetres, positive below the focusing circle, or `None` when it was not refined, and its esd |
| `radius_mm` | the radius the displacement was refined with, as it was given |
| `n_peaks` | how many lone peaks the fit used |
| `rms_two_theta` | the root mean square residual in degrees |
| `residuals` | observed less calculated for each of those peaks, in degrees, as an array |
| `hkl` | the reflection of each, as a list of tuples in the same order |
| `converged` | whether the optimiser reported success |
| `parameter_names` | the free cell parameters, in the order of the vector |
| `covariance` | the covariance matrix of the refined parameters |

`covariance` is in the order of the free vector: the cell parameters first,
then the zero and the displacement where they were refined. It is the inverse
curvature at the solution scaled by the reduced chi squared, so the error bars
reflect how well the model actually fits rather than only how sharply the
surface curves, and its off-diagonal entries are where the correlation between
a parameter and a systematic error can be read.

`esd_volume` is propagated through the cell block of that matrix, correlations
included, the derivative of the volume with respect to each free parameter
being taken by central differences. It is therefore not the same number as
adding the parameters' own esds in quadrature, which is what `cell_volume` in
Section 19 does with esds typed in by hand.

`lattice_fit_to_dict(fit)` flattens the scalar fields into a dict ready for a
CSV row: all six parameters with their esds, `c_over_a` for a tetragonal,
hexagonal or trigonal cell, the volume and its esd, the zero, the displacement
and the radius, the peak count, the rms and `converged`. The residuals, the
hkl and the covariance are left out, being more than one value each. It is
what writes the `lattice_STEM.csv` of Section 7.7.

### 18.3 Refining a cell on a peak list you have edited

A command cannot be handed an edited peak list. `xrdkit lattice` finds the
peaks, refits them, indexes them and refines, and there is no point in the
middle where a peak can be struck out. That is the reason this script exists,
and it is the only reason: everything else it does the command does better,
because the command refits every peak as a doublet first.

Striking peaks out is worth doing when a pattern holds a second phase, or a
tube line, or a reflection whose position is known to be poor. The honest way
to record such an edit is to say which rows went and why. The script below
makes the edit in code, by position in the peak list, so that the run is
reproducible; deleting those rows from `results/library/peaks_pellet_a.csv` by
hand and reading the file back is the same edit and gives the same answer.

The two rows here are the ones Section 15.8 marked, the peaks whose nearest
reflection is more than 0.02 degrees away once the cell is fixed. The sample
is `pellet_a` fitted as the command fits a pellet: the displacement free, the
zero held at zero, and the goniometer radius of 145 millimetres that its
instrument table of Section 3.2 records.

The whole of `lattice_fit.py`:

```python
import numpy as np

from xrdkit.cell import Cell
from xrdkit.indexing import index_and_refine
from xrdkit.io import read_scan
from xrdkit.lattice import refine_lattice
from xrdkit.peaks import exclude_kalpha2, find_peaks

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
START_CELL = Cell.tetragonal(12.45, 3.94)
SPACE_GROUP = "P4bm"
WINDOW = (10.0, 100.0)
RADIUS_MM = 145.0
DROPPED_ROWS = (2, 41)

scan = read_scan(SCAN_FILE)
found = exclude_kalpha2(find_peaks(scan, two_theta_range=WINDOW))
dropped = ", ".join(f"{found[row].two_theta:.3f}" for row in DROPPED_ROWS)
peaks = [peak for row, peak in enumerate(found) if row not in DROPPED_ROWS]
print(f"{len(found)} peaks found, rows {list(DROPPED_ROWS)} deleted: {dropped} degrees")

indexed, _ = index_and_refine(
    peaks,
    start_cell=START_CELL,
    wavelength=scan.wavelength,
    space_group=SPACE_GROUP,
)
fit = refine_lattice(
    indexed,
    scan.wavelength,
    START_CELL,
    fit_zero=False,
    fit_displacement=True,
    radius_mm=RADIUS_MM,
)

print(f"parameter_names {fit.parameter_names}, converged {fit.converged}")
for name in ("a", "b", "c", "alpha", "beta", "gamma"):
    esd = getattr(fit, f"esd_{name}")
    shown = "not refined" if esd is None else f"+/- {esd:.4f}"
    print(f"  {name:5s} {getattr(fit, name):8.4f}  {shown}")
print(f"V = {fit.volume:.3f} +/- {fit.esd_volume:.4f} cubic angstrom")
print(f"zero {fit.zero:.4f} degrees, esd {fit.esd_zero}")
print(
    f"displacement {fit.displacement:.4f} +/- {fit.esd_displacement:.4f} mm, "
    f"radius {fit.radius_mm:.0f} mm"
)
print(f"rms {fit.rms_two_theta:.4f} degrees on {fit.n_peaks} peaks")

worst = int(np.argmax(np.abs(fit.residuals)))
print(f"largest residual {fit.residuals[worst]:+.4f} degrees at {fit.hkl[worst]}")
print(f"first three hkl used {fit.hkl[:3]}")
deviations = np.sqrt(np.diag(fit.covariance))
correlation = fit.covariance / np.outer(deviations, deviations)
print("correlation of a, c and the displacement")
for row in correlation:
    print("  " + "  ".join(f"{value:+.3f}" for value in row))
```

It prints the edit, then every field of the fit that is one number, then the
worst residual and its reflection, then the correlations.

```text
43 peaks found, rows [2, 41] deleted: 26.921, 94.318 degrees
parameter_names ('a', 'c'), converged True
  a      12.4740  +/- 0.0006
  b      12.4740  not refined
  c       3.9305  +/- 0.0002
  alpha  90.0000  not refined
  beta   90.0000  not refined
  gamma  90.0000  not refined
V = 611.593 +/- 0.0674 cubic angstrom
zero 0.0000 degrees, esd None
displacement -0.2020 +/- 0.0034 mm, radius 145 mm
rms 0.0063 degrees on 32 peaks
largest residual -0.0121 degrees at (0, 0, 1)
first three hkl used [(0, 0, 1), (3, 2, 0), (2, 1, 1)]
correlation of a, c and the displacement
  +1.000  +0.120  -0.756
  +0.120  +1.000  -0.527
  -0.756  -0.527  +1.000
```

Take the parameters first. Only a and c are refined, so only those two carry
an esd; b is equal to a and the three angles are 90 because the crystal system
says so, and all four come back with an esd of `None` rather than of zero. The
zero was held, so `esd_zero` is `None` as well while `zero` is the 0.0000
degrees it was held at. `parameter_names` is `('a', 'c')`, which is the order
of the first two rows and columns of the covariance.

The displacement comes out at minus 0.2020 millimetres, about sixty times its
own esd, and the rms is 0.0063 degrees on 32 peaks. Thirty two and not forty
one, because every one of the forty one was indexed but nine of them carried
more than one candidate and were left out.

The correlations are the reason the section keeps saying the zero and the
displacement are the same parameter twice over. With the zero held there are
three free numbers and a and the displacement correlate at minus 0.756: a
shift of the specimen and a change of the a axis do much the same thing to
this pattern, and the fit can only tell them apart because they do it in
different proportions at different angles. Had the zero been free as well
there would have been a fourth row, correlating with the displacement more
strongly still.

Now compare the cell with Section 7.1, which is the same scan, the same start
cell, the same space group and the same choice of what to refine, run through
`xrdkit lattice pellet_a`.

| | Section 7.1 | this script |
| --- | --- | --- |
| a, angstrom | 12.4740 +/- 0.0011 | 12.4740 +/- 0.0006 |
| c, angstrom | 3.9295 +/- 0.0004 | 3.9305 +/- 0.0002 |
| V, cubic angstrom | 611.426 +/- 0.136 | 611.593 +/- 0.067 |
| displacement, mm | -0.2033 +/- 0.0069 | -0.2020 +/- 0.0034 |
| rms, degrees | 0.0128 | 0.0063 |
| peaks used | 33 | 32 |

The a axis agrees to four decimal places and the displacement agrees well
within one esd, which is the check that the script is doing the command's
arithmetic. The c axis is a thousandth of an angstrom higher here, which is
between two and five times the esds involved, and the volume follows it up by
0.167 cubic angstrom.

Two differences in the input account for that. The command refits every
unflagged peak as the K alpha 1 line of a doublet before it indexes anything,
which moves a blended position to where the K alpha 1 line really is, and it
recovers one flagged peak that turned out to be a genuine reflection. This
script uses the peak finder's positions as they stand, which are parabola
vertices through a blend wherever the doublet is unresolved, and it never sees
the recovered peak at all. It also has two rows struck out that the command
kept.

Do not read the smaller esds here as the better answer. They are smaller
because the rms is half the command's, and the rms is half the command's
partly because two of the worst-fitting peaks were deleted by hand. An esd
measures the spread of the peaks that were kept. It says nothing about the
peaks that were not, and nothing at all about the doublet bias that the
command removes and this script leaves in. Where both can be run, quote the
command.

## 19. density

`xrdkit.density` turns a cell and a composition into the density the material
would have with no porosity, and compares a measured density with it. It is
what `xrdkit density` and the density half of `xrdkit lattice` are, and
Section 7.6 says what the answer is worth. There are six public names: the
table of masses, two functions over a formula, one over a cell and two over a
density.

### 19.1 parse_formula, formula_mass and ATOMIC_MASSES

`parse_formula(text)` reads a chemical formula and returns a dict of element
symbol to atoms per formula unit, in order of first appearance. A formula is a
run of element symbols, each followed by an optional count that may be
fractional, as in `"Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"` or `"LaB6"`. Parentheses
group symbols under one multiplier, as in `"Ca(OH)2"`, and may nest, and an
element named more than once has its counts added.

It raises `ValueError` when the formula is empty, when it holds anything other
than element symbols, counts and balanced parentheses, when it names an
element the mass table lacks, or when a pair of parentheses is empty, and the
message quotes the text it could not read.

`formula_mass(composition)` returns the mass of one formula unit in grams per
mole. `composition` is either a mapping of element symbol to coefficient or a
formula string, which it passes through `parse_formula` first, so the two ways
of saying the same thing give the same number. It raises `ValueError` on an
empty composition, an element with no mass in the table, or a negative
coefficient, and from `parse_formula` on a formula it cannot read.

`ATOMIC_MASSES` is the table behind both: a dict of element symbol to standard
atomic weight in grams per mole, every element from hydrogen to uranium, from
the IUPAC 2021 table. The fourteen elements whose standard atomic weight is an
interval rather than a number, hydrogen, lithium, boron, carbon, nitrogen,
oxygen, magnesium, silicon, sulphur, chlorine, argon, bromine, thallium and
lead, take the conventional value IUPAC gives for them. The eight with no
stable isotope and so no standard atomic weight, technetium, promethium,
polonium, astatine, radon, francium, radium and actinium, take the mass number
of their longest lived isotope. Read the dict to put a mass into a table of
your own rather than copying one out.

### 19.2 cell_volume

`cell_volume(cell, esd=None)` returns the volume of `cell` in cubic angstroms
and its esd, as a pair. `esd` maps free parameters of the cell, the names of
`cell.parameter_names`, to their esds. The esd of the volume is propagated to
first order, the derivative with respect to each free parameter taken by
central differences with the dependent parameters moving along, and the
parameters are treated as uncorrelated, since no covariance is given here. A
free parameter left out of the mapping counts as zero, and the esd comes back
as `None` when no mapping was given at all or when every esd in it is zero. It
raises `ValueError` when `esd` names a parameter that is not free in the cell.

That treatment is the difference between this function and
`LatticeFit.esd_volume` of Section 18, which carries the full covariance of
the fit and so the correlations between the parameters as well. Section 19.4
puts the two side by side.

### 19.3 theoretical_density and relative_density

`theoretical_density(composition, z, volume_a3, esd_volume_a3=None)` returns
the X-ray density in grams per cubic centimetre and its esd, as a pair. The
density is Z M over N sub A V, with M the formula mass, Z the formula units
per cell and V the cell volume. `composition` is whatever `formula_mass`
takes. Only the volume carries an error here, so the density has the same
relative esd as the volume, and the esd is `None` when no volume esd was
given. It raises `ValueError` when Z or the volume is not positive, and from
`formula_mass` on a composition it cannot read.

`relative_density(measured, theoretical, esd_measured=None,
esd_theoretical=None)` returns the measured density as a percentage of the
theoretical one, and its esd, as a pair. The esd adds the two relative errors
in quadrature. An esd left as `None` counts as zero, and the esd comes back as
`None` only when neither was given. It raises `ValueError` when either density
is not positive.

### 19.4 A density from a refined cell

The cell below is the one Section 18 refined, with its esds, typed into the
settings the way `xrdkit density` is given a cell on the command line. The
composition and the Z are the `composition` and `z` of the
`[structures.ttb_p4bm]` table of Section 2.2. The Archimedes measurement is an
example rather than one of the example files, which carry none, and it is
there to show what `relative_density` does with one.

The whole of `density_from_cell.py`:

```python
from xrdkit.cell import Cell
from xrdkit.density import (
    ATOMIC_MASSES,
    cell_volume,
    formula_mass,
    parse_formula,
    relative_density,
    theoretical_density,
)

# Edit these lines for each new sample. Nothing below needs changing.
FORMULA = "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
Z = 5
A, C = 12.4740, 3.9305
ESD_A, ESD_C = 0.0006, 0.0002
ARCHIMEDES, ESD_ARCHIMEDES = 5.15, 0.02

composition = parse_formula(FORMULA)
print(composition)
print(f"{sum(composition.values()):.1f} atoms per formula unit")
for element, n in composition.items():
    print(
        f"  {element:2s} {n:5.2f} x {ATOMIC_MASSES[element]:9.4f} "
        f"= {n * ATOMIC_MASSES[element]:8.4f} g/mol"
    )
mass = formula_mass(FORMULA)
print(f"M = {mass:.3f} g/mol per formula unit")
print(f"from the dict: {formula_mass(composition) == mass}")

cell = Cell.tetragonal(A, C)
volume, esd_volume = cell_volume(cell, esd={"a": ESD_A, "c": ESD_C})
print(f"V = {volume:.3f} +/- {esd_volume:.4f} cubic angstrom")
print(f"without esds: {cell_volume(cell)}")

density, esd_density = theoretical_density(FORMULA, Z, volume, esd_volume)
print(f"theoretical density {density:.4f} +/- {esd_density:.4f} g/cm3")
print(f"relative esd {esd_density / density:.2e}, volume {esd_volume / volume:.2e}")

ratio, esd_ratio = relative_density(ARCHIMEDES, density, ESD_ARCHIMEDES, esd_density)
print(f"relative density {ratio:.2f} +/- {esd_ratio:.2f} per cent")

for description, call in (
    ("an element the table lacks", lambda: formula_mass("SrNqO3")),
    ("a formula it cannot read", lambda: parse_formula("Sr0.4-Ba0.5")),
    ("an empty group", lambda: parse_formula("Ca()2")),
):
    try:
        call()
    except ValueError as error:
        print(f"{description}: {error}")
```

It prints the composition and the mass term by term, then the volume, then the
density, then the relative density, then three formulas it refuses.

```text
{'Sr': 0.4, 'Ba': 0.5, 'La': 0.1, 'Nb': 1.9, 'Ti': 0.1, 'O': 6.0}
9.0 atoms per formula unit
  Sr  0.40 x   87.6200 =  35.0480 g/mol
  Ba  0.50 x  137.3270 =  68.6635 g/mol
  La  0.10 x  138.9055 =  13.8905 g/mol
  Nb  1.90 x   92.9064 = 176.5221 g/mol
  Ti  0.10 x   47.8670 =   4.7867 g/mol
  O   6.00 x   15.9990 =  95.9940 g/mol
M = 394.905 g/mol per formula unit
from the dict: True
V = 611.588 +/- 0.0666 cubic angstrom
without esds: (611.5884570180001, None)
theoretical density 5.3611 +/- 0.0006 g/cm3
relative esd 1.09e-04, volume 1.09e-04
relative density 96.06 +/- 0.37 per cent
an element the table lacks: unknown element 'Nq' in formula 'SrNqO3'
a formula it cannot read: cannot read '-Ba0.5' in formula 'Sr0.4-Ba0.5'
an empty group: empty group '()' in formula 'Ca()2'
```

The mass of 394.905 grams per mole is the M that `xrdkit lattice pellet_a`
reported in Section 7.1, which is the same composition through the same table.
Printing it term by term is worth doing once for a composition of your own:
the commonest way to get a density badly wrong is a coefficient that does not
mean what it was meant to, and a line of the table showing oxygen at six atoms
and 95.994 grams per mole says at a glance whether the formula was read as
intended.

Now the volume. Section 18 refined a of 12.474044 and c of 3.930501 and got a
volume of 611.593 cubic angstrom; this script types those in to four decimal
places and gets 611.588. The five thousandths are the rounding, and they are a
fair warning: a volume goes as a squared times c, so a thousandth of an
angstrom of rounding in a shows up multiplied by two in the volume and again
in the density.

The esd differs too, and for a different reason: 0.0666 cubic angstrom here
against 0.0674 in Section 18. Both are the same first order propagation
through the same derivatives. What differs is what they are propagated
through. Section 18 has the covariance of the fit, in which a and c correlate
at plus 0.120, and pushes the esd through all of it. This function has two
numbers typed in and no way of knowing they were ever correlated, so it adds
them in quadrature as though they were independent, which here understates the
esd by about one per cent. Quote the refinement's number where you have it,
as Section 7.6 says.

The theoretical density of 5.3611 grams per cubic centimetre carries a
relative esd of 1.09 times ten to the minus four, exactly the volume's,
because the formula mass and Z are taken as exact. Against them the measured
density dominates the relative density entirely: the Archimedes value of 5.15
plus or minus 0.02 is a relative error of four parts in a thousand, some
thirty five times the cell's, so the 96.06 plus or minus 0.37 per cent is very
nearly the weighing error alone. This is the usual state of affairs, and it is
why Section 7.5 asks for a cell good to a thousandth of an angstrom and then
stops asking for more.

## 20. symmetry

`xrdkit.symmetry` is where the reflection conditions, the equivalence of
reflections and the multiplicities of Section 16 come from. It knows eight
space groups, which is the limitation Section 10.1 records, and it works
entirely in exact arithmetic: rotations are integer matrices and translations
are fractions, so nothing here is ever decided by a tolerance.

There are nine public names: the tuple of symbols, two functions that produce
operations, one that tests a reflection against them, and five over a Laue
group.

### 20.1 What an operation is, and parse_xyz

An operation is a pair (R, t): R a three by three integer rotation matrix, as
a tuple of rows, and t a translation of three `Fraction` reduced modulo 1, so
that a fractional position x goes to R x plus t. A reflection h, a row of
three integers, transforms the other way round, as h times R.

`parse_xyz(text)` returns the operation of one line of the xyz notation the
`_symmetry_equiv_pos_as_xyz` field of a CIF uses, such as `"-y,x,z"` or
`"x+1/2,-y+1/2,z"`. Spaces are ignored, case does not matter, and the
translation is reduced modulo 1. It raises `ValueError` when the text does not
have three components, when a component holds something it cannot read or has
no x, y or z in it, when a term after the first carries no sign, or when a
translation is not a multiple of one twelfth, which every space group's
translations are.

### 20.2 SUPPORTED_SPACE_GROUPS and space_group_operations

`SUPPORTED_SPACE_GROUPS` is the tuple of symbols the module knows: `Pm-3m`,
`P4mm`, `P4bm`, `P4/mbm`, `R3c`, `R3m`, `Pbnm` and `Amm2`. Each is stored as a
few generators from International Tables Volume A in the setting its symbol
names, with the centring translations among the generators. `R3c` and `R3m`
are on hexagonal axes in the obverse setting, so a cell in the rhombohedral
setting has to be converted first, and `Pbnm` is the cab setting of Pnma.

`space_group_operations(symbol)` returns every operation of the group,
centring included, with the identity first. The generators are expanded by
closure once and the result is cached, so a zero offset search asking for the
same group two hundred times pays for it once. It raises `ValueError` when the
symbol is not one of the eight, listing the eight, and when the generators do
not close to the order the group should have, which is a check on the table
rather than on the caller.

`xrdkit plot` and `xrdkit lattice` catch the first of those: given a symbol
outside the eight they print a note and index without conditions, so a peak
may be labelled with a reflection its group forbids and the labels are
provisional until they are checked by hand. Section 10.1 is the whole of that
limitation.

### 20.3 is_absent

`is_absent(hkl, operations)` returns whether a reflection is systematically
absent under those operations. A reflection is absent when some operation
leaves it fixed, h times R equal to h, while h dotted into t is not an
integer. That one test gives the integral, zonal and serial conditions of a
group together, without any of them being written out as a rule: the centring
conditions come from the operations whose R is the identity, and a glide or a
screw gives the zone or the row it forbids.

Equivalent reflections share their absence, so `generate_reflections` tests
one member of each family and drops the family on its answer.

### 20.4 The Laue group and its orbits

`laue_group(operations)` returns the Laue group of a set of operations: their
distinct rotations together with the inversion, closed under multiplication
and sorted. The translations play no part, which is why two space groups
differing only in their glides share a Laue group.

`holohedry(crystal_system)` returns the Laue group of the holohedral class of
a crystal system, for when no space group is given: m-3m, 4/mmm, mmm, 6/mmm on
hexagonal axes for both hexagonal and trigonal, 2/m with b unique, or -1. It
raises `ValueError` on anything that is not one of the seven systems.

`laue_orbit(hkl, laue)` returns the distinct reflections equivalent to a given
one under that group, largest first, and `multiplicity(hkl, laue)` is how many
there are.

`representative(hkl, laue)` is the label of the orbit: of the members with no
negative index, or of all of them if none qualifies, the lexicographically
largest, h first, then k, then l. That gives h greater than or equal to k
greater than or equal to 0 and l greater than or equal to 0 for a tetragonal
group, (100) for the cubic family of that name, (110) rather than (2 -1 0) on
hexagonal axes, and it keeps the sign of l in a monoclinic (h 0 -l). On
obverse hexagonal axes -3m keeps (h k l) and (k h l) with l not zero in
separate orbits, so a label never swaps h and k: (101) rather than (011),
which is absent in R3m, and (104) and (012) in R3c, where (014) and (102) are
absent.

### 20.5 The group of the example structure

P4bm is the space group of the `[structures.ttb_p4bm]` table of Section 2.2
and of every indexing run in this guide. The script below takes it apart: the
operations, the conditions they imply, the one absent reflection Section 16.8
met, and the orbit of a reflection with its label and its multiplicity.

The whole of `space_group.py`:

```python
from xrdkit.symmetry import (
    SUPPORTED_SPACE_GROUPS,
    holohedry,
    is_absent,
    laue_group,
    laue_orbit,
    multiplicity,
    parse_xyz,
    representative,
    space_group_operations,
)

# Edit these lines for each new group. Nothing below needs changing.
SPACE_GROUP = "P4bm"
CRYSTAL_SYSTEM = "tetragonal"
HKL = (3, 1, 0)
ZONES = {
    "0kl": [(0, k, 1) for k in range(1, 5)],
    "h0l": [(h, 0, 1) for h in range(1, 5)],
    "h00": [(h, 0, 0) for h in range(1, 5)],
    "hk0": [(h, k, 0) for h in range(1, 4) for k in range(1, 4)],
    "hhl": [(h, h, 1) for h in range(1, 5)],
}


def as_xyz(operation):
    """Write an operation back in the notation parse_xyz reads."""
    rotation, translation = operation
    components = []
    for row, shift in zip(rotation, translation):
        terms = [
            f"{'-' if value < 0 else '+'}{axis}"
            for axis, value in zip("xyz", row)
            if value
        ]
        if shift:
            terms.append(f"+{shift}")
        components.append("".join(terms).lstrip("+"))
    return ",".join(components)


print(
    f"{len(SUPPORTED_SPACE_GROUPS)} space groups: {', '.join(SUPPORTED_SPACE_GROUPS)}"
)
operations = space_group_operations(SPACE_GROUP)
print(f"{SPACE_GROUP} has {len(operations)} operations")
for start in range(0, len(operations), 4):
    row = operations[start : start + 4]
    print("  " + "  ".join(f"{as_xyz(operation):16s}" for operation in row).rstrip())
rotation, translation = parse_xyz("x+1/2,-y+1/2,z")
print(f"parse_xyz rotation {rotation}")
print(f"parse_xyz translation {translation}")

for zone, family in ZONES.items():
    absent = [hkl for hkl in family if is_absent(hkl, operations)]
    print(f"{zone:4s} absent {absent}" if absent else f"{zone:4s} none absent")
print(
    f"(100) absent {is_absent((1, 0, 0), operations)}, "
    f"(200) absent {is_absent((2, 0, 0), operations)}"
)

laue = laue_group(operations)
print(f"Laue group of {SPACE_GROUP}: {len(laue)} rotations")
print(f"same as the {CRYSTAL_SYSTEM} holohedry: {laue == holohedry(CRYSTAL_SYSTEM)}")
orbit = laue_orbit(HKL, laue)
print(f"orbit of {HKL}, largest first")
for start in range(0, len(orbit), 4):
    row = orbit[start : start + 4]
    print("  " + "  ".join(f"{member!s:12s}" for member in row).rstrip())
print(
    f"representative {representative(HKL, laue)}, "
    f"multiplicity {multiplicity(HKL, laue)}"
)

try:
    space_group_operations("P4/mmm")
except ValueError as error:
    print(f"an unsupported symbol: {error}")
```

It prints the eight groups, the operations of this one written back in the
notation they were read from, one operation as it is really held, the zones,
the Laue group, an orbit, and the message an unsupported symbol gets.

```text
8 space groups: Pm-3m, P4mm, P4bm, P4/mbm, R3c, R3m, Pbnm, Amm2
P4bm has 8 operations
  x,y,z             -y,x,z            x+1/2,-y+1/2,z    -x,-y,z
  y+1/2,x+1/2,z     -y+1/2,-x+1/2,z   y,-x,z            -x+1/2,y+1/2,z
parse_xyz rotation ((1, 0, 0), (0, -1, 0), (0, 0, 1))
parse_xyz translation (Fraction(1, 2), Fraction(1, 2), Fraction(0, 1))
0kl  absent [(0, 1, 1), (0, 3, 1)]
h0l  absent [(1, 0, 1), (3, 0, 1)]
h00  absent [(1, 0, 0), (3, 0, 0)]
hk0  none absent
hhl  none absent
(100) absent True, (200) absent False
Laue group of P4bm: 16 rotations
same as the tetragonal holohedry: True
orbit of (3, 1, 0), largest first
  (3, 1, 0)     (3, -1, 0)    (1, 3, 0)     (1, -3, 0)
  (-1, 3, 0)    (-1, -3, 0)   (-3, 1, 0)    (-3, -1, 0)
representative (3, 1, 0), multiplicity 8
an unsupported symbol: unsupported space group 'P4/mmm'; the supported space groups are Pm-3m, P4mm, P4bm, P4/mbm, R3c, R3m, Pbnm, Amm2
```

The eight operations are the four rotations of a four fold axis and four
mirrors, four of them carrying a translation of a half along a and b. That
half is the b glide of the symbol, and everything the group forbids follows
from it. No condition applies to hk0 or to hhl; a 0kl is absent when k is odd,
an h0l when h is odd, and an h00 when h is odd, which is the h0l condition at
l equal to zero. Those are the conditions International Tables lists for P4bm,
and nothing in the module wrote them down: the script asked `is_absent` about
a few reflections and read them off the answers.

That is where the (100), (101) and (300) of Section 16.8 come from. (100) is
absent and (200) is not, and any tungsten bronze pattern indexed in P4bm will
have no reflection labelled (100), which is a thing worth knowing before a
peak near 7 degrees is explained.

The Laue group of P4bm has sixteen rotations and is the tetragonal holohedry,
4/mmm. It must be: the Laue group throws the translations away, so the b glide
that distinguishes P4bm from P4mm cannot survive into it, and the inversion is
added whether the group has one or not. This is why the multiplicities in
Section 16.8 are the same with a space group and without one, while the
absences are not.

The orbit of (310) holds eight reflections and its representative is (310)
itself, the largest of the two members with no negative index. A multiplicity
of eight is what `generate_reflections` puts on that row, and it is the number
an intensity calculation needs: eight reflections of the same d spacing
contribute to one peak.

## 21. broadening

`xrdkit.broadening` measures peak widths properly and takes the instrument out
of them. It has seventeen public names: four line shape functions, two that
convert between a pseudo-Voigt and the Voigt it approximates, two dataclasses
and two functions for the two fits, the correction and its dataclass, a model
of specimen height spread, and the three names of the breadth model
comparison.

The widths `find_peaks` reports are read straight off the counts at half the
prominence. They carry no esd, and they include whatever of the K alpha 2
satellite falls inside the half maximum, which is all of it at low angle and a
changing fraction of it once the pair begins to separate. That is good enough
to find a peak and not good enough to measure a broadening, which is what this
module exists for.

Nothing here makes the instrument parameter file. `xrdkit instrument` does
that, and Section 3 is the whole of it; Section 22 is the two functions the
command is built from.

### 21.1 The line shapes

`pseudo_voigt(two_theta, centre, fwhm, eta)` is a line of unit area whose
Lorentzian and Gaussian components share one FWHM and are mixed as `eta` times
the Lorentzian plus one minus `eta` times the Gaussian, so `eta` runs from 0
for a pure Gaussian to 1 for a pure Lorentzian.

`split_pseudo_voigt(two_theta, centre, fwhm, eta, asymmetry)` is the same line
with different half widths either side. The half width below `centre` is
`fwhm` times one minus `asymmetry`, over two, and the half width above it
`fwhm` times one plus `asymmetry`, over two, so `fwhm` is still the full width
at half maximum and the two halves meet at the same height. A negative
`asymmetry` puts a tail on the low angle side, which is what axial divergence
does to a low angle reflection, and an asymmetry of zero is `pseudo_voigt`
again. The refinement holds it short of plus or minus `MAX_ASYMMETRY`, 0.95,
at which one half width would vanish.

A symmetric line fitted to an asymmetric reflection misses its top and so
misreads its width, which is why the split shape rather than the plain one is
what `fit_profile` uses.

### 21.2 Where the satellite falls, and how much room a doublet has

`kalpha2_position(two_theta, wavelength_ratio=KALPHA2_RATIO)` returns where
the K alpha 2 line of a K alpha 1 line at `two_theta` falls. The satellite
diffracts at the same d spacing at the longer wavelength, so it always lies to
high angle, and by more the higher the angle.
`KALPHA2_INTENSITY_RATIO`, 0.5, is the integrated intensity of the satellite
over its parent that `fit_profile` assumes.

`doublet_gaps(two_theta, others, wavelength_ratio=KALPHA2_RATIO)` returns the
clear space either side of a whole doublet, in degrees. `others` are the K
alpha 1 positions of every other line that might lie nearby, and each brings
its own satellite along. The first number runs down from the K alpha 1 line to
the nearest line below it and the second up from the K alpha 2 line to the
nearest line above. A line falling between the two members closes both gaps to
zero, and a side with no line at all is infinitely clear.

That is the function to select reflections with before any width is measured.
A width is only the sample's if nothing else is under it, and what counts as
nearby is not the peaks the finder happened to report but every reflection the
cell calculates, which Section 23.6 uses it on.

### 21.3 ProfileFit and fit_profile

`fit_profile(two_theta, intensity, centre, fwhm_guess, window=None,
fit_asymmetry=True, wavelength_ratio=KALPHA2_RATIO,
intensity_ratio=KALPHA2_INTENSITY_RATIO)` fits one reflection as a K alpha 1
and K alpha 2 doublet of split pseudo-Voigt lines on a linear background.

The satellite is tied to its parent rather than fitted: it sits where the
parent's d spacing puts the longer wavelength, carries `intensity_ratio` of
its area, and shares its width, mixing parameter and asymmetry. The free
parameters are therefore seven, the K alpha 1 position, FWHM, mixing
parameter, asymmetry and area and the background level and slope, however many
lines are drawn. Counts are weighted as Poisson, one over the counts or one,
whichever is larger.

| Argument | What it does |
| --- | --- |
| `two_theta`, `intensity` | the pattern, in degrees and counts |
| `centre` | the starting K alpha 1 position in degrees |
| `fwhm_guess` | the starting FWHM in degrees; a width read off the raw doublet, such as `Peak.fwhm`, is close enough |
| `window` | the `(low, high)` range to fit, in degrees. By default it runs `WINDOW_FWHM`, ten, starting widths below `centre` and the same beyond the satellite, cut at the ends of the data |
| `fit_asymmetry` | refine the asymmetry. `False` holds it at zero and the lines are symmetric |
| `wavelength_ratio` | K alpha 2 over K alpha 1 wavelength |
| `intensity_ratio` | K alpha 2 over K alpha 1 integrated intensity. Zero fits a single line, which is what monochromated radiation wants |

It raises `ValueError` when `fwhm_guess` is not positive or the window holds
no more points than there are free parameters.

`ProfileFit` carries the result. `two_theta` and `fwhm` belong to the K alpha
1 line and the satellite sits at `kalpha2_two_theta` with the same width and
shape.

| Field | What it holds |
| --- | --- |
| `two_theta`, `esd_two_theta` | the K alpha 1 position and its esd, in degrees |
| `fwhm`, `esd_fwhm` | its full width at half maximum and its esd, in degrees |
| `eta`, `esd_eta` | the mixing parameter and its esd |
| `asymmetry`, `esd_asymmetry` | the asymmetry and its esd, the esd `None` when it was held |
| `area`, `esd_area` | the integrated intensity of the K alpha 1 line above background, in counts times degrees |
| `background`, `slope` | the linear background at the middle of the window and its slope per degree |
| `kalpha2_two_theta` | where the satellite was placed, in degrees |
| `window`, `n_points` | the range actually fitted and how many points it held |
| `reduced_chi_squared`, `r_wp` | the residual per degree of freedom and the weighted profile R factor, as a fraction |
| `converged` | whether the optimiser reported success |
| `truncated` | whether the data stop within `MIN_MARGIN_FWHM`, three, fitted widths of either end of the doublet, so that a tail and the background under it went unseen |
| `fitted` | the fitted curve over the window, for drawing or inspecting |

Esds are scaled by the reduced chi squared, as every other fit in the kit
scales its own.

`converged` is not by itself the test of a good fit, and Section 21.8 shows
why. A caller has to decide what to reject. `xrdkit lattice` rejects a fit
that failed, did not converge, or moved the position by more than the found
peak's own FWHM, and keeps the found position instead; that is the rule
Section 7.2 describes and the one Section 23.6 applies.

### 21.4 Caglioti and fit_caglioti

`fit_caglioti(two_theta, fwhm, weights=None)` fits the Caglioti relation,
FWHM squared equals U tan squared theta plus V tan theta plus W, by weighted
least squares. The relation is linear in U, V and W, so this is a single
linear solve, and the covariance is scaled by the reduced chi squared so that
the esds reflect how well the relation really describes the widths. It needs
at least `MIN_CAGLIOTI_POINTS`, four, widths, one more than the three
parameters so that a residual is left to scale by.

It raises `ValueError` when the arrays differ in length, there are fewer than
four points, a position lies outside 0 to 180 degrees, a width is not
positive, or a weight is negative or not finite.

`Caglioti` carries U, V and W in degrees squared of two theta with their esds,
`n_peaks`, `rms`, the root mean square of the FWHM residuals in degrees, and
the 3 by 3 `covariance`. Two methods evaluate it: `fwhm(two_theta)` is the
width at any angle, taken as zero where a fitted negative V drives the
quadratic below zero outside the range of the data, and `fwhm_esd(two_theta)`
is its esd from the covariance, `nan` where the width is zero.

That pair is what makes the relation useful. A sample reflection almost never
sits at the angle of a standard reflection, and `fwhm` and `fwhm_esd` give the
instrumental width and its esd wherever it is wanted.

A note on the weights. The docstring of `fit_caglioti` says that for a width
with esd s the weight in the sum of squared FWHM squared residuals is
1 / (2 FWHM s) squared, which is the weight that carries the esd of a width
through to the square of that width. `fit_instrument_widths` in Section 22,
which is the kit's own caller and what `xrdkit instrument` runs, passes
1 / s squared instead, which weights the width rather than its square. Both
are defensible and the fit is judged on the FWHM residuals either way, but the
two are not the same weighting and the example standard gives U of 0.0107
against 0.0105 depending on which is used. Pass the weights you mean.

### 21.5 Between a pseudo-Voigt and the Voigt it stands for

A pseudo-Voigt is a sum of a Gaussian and a Lorentzian; a Voigt is their
convolution. Broadening adds in the convolution, so a correction has to go
through the Voigt. Thompson, Cox and Hastings give the pair of approximations
that pass between them, and two functions apply them.

`pseudo_voigt_from_components(fwhm_gaussian, fwhm_lorentzian)` returns the
FWHM and mixing parameter of the pseudo-Voigt matching the Voigt of those two
components. It raises `ValueError` when either width is negative or both are
zero.

`pseudo_voigt_components(fwhm, eta)` is the inverse: the Gaussian and
Lorentzian FWHM of the Voigt a pseudo-Voigt matches. The cubic is solved for
the Lorentzian fraction of the width and the quintic then for the Gaussian
width that makes up the rest, both by bracketed root finding, so the round
trip closes. It raises `ValueError` when `fwhm` is not positive or `eta` lies
outside 0 to 1.

`integral_breadth(fwhm_gaussian, fwhm_lorentzian)` is the integral breadth of
that Voigt, area over height, exact rather than approximated. It is the
breadth for which the Scherrer constant is 1 whatever the crystallite shape,
which is why Section 23 has a function that takes it.

### 21.6 correct_broadening and BroadeningCorrection

`correct_broadening(fwhm_obs, eta_obs, fwhm_inst, eta_inst, esd_fwhm_obs=0.0,
esd_eta_obs=0.0, esd_fwhm_inst=0.0, esd_eta_inst=0.0,
significance=DEFAULT_SIGNIFICANCE)` takes the instrumental broadening out of
an observed profile.

Both profiles are split into their Gaussian and Lorentzian FWHM. Lorentzians
convolve by adding widths and Gaussians by adding squares, so the sample
Lorentzian is the observed one less the instrumental one and the sample
Gaussian is the root of the difference of squares, and the two are recombined
into one pseudo-Voigt. Esds are propagated linearly through numerical
derivatives, treating the four inputs as independent; the observed width and
mixing parameter of a fit are in fact correlated, so the esds are indicative
rather than exact, and the docstring says so.

The important argument is `significance`, `DEFAULT_SIGNIFICANCE` of 2.0 by
default. The observed width must exceed the instrumental one by that many
combined esds before any sample breadth is reported at all. Below it the
correction returns with `unresolved` true and every sample quantity `None`,
which is the module refusing to turn noise into a crystallite size. `excess`
is the observed less the instrumental FWHM over their combined esd, reported
either way, so a caller can see how close a reflection came.

`BroadeningCorrection` is what comes back.

| Field | What it holds |
| --- | --- |
| `fwhm`, `esd_fwhm` | the sample FWHM and its esd, in the units given |
| `eta`, `esd_eta` | its mixing parameter and esd |
| `fwhm_gaussian`, `fwhm_lorentzian` and their esds | the two components of the sample profile |
| `integral_breadth`, `esd_integral_breadth` | the integral breadth of the corrected Voigt |
| `fwhm_linear`, `fwhm_quadrature` and their esds | the two simple estimates for comparison: the plain difference of the widths, as if both were Lorentzian, and the root of the difference of their squares, as if both were Gaussian |
| `gaussian_clipped`, `lorentzian_clipped` | whether that component came out below zero and was set to zero |
| `unresolved` | whether the reflection failed the significance test, in which case every field above is `None` |
| `excess` | the observed less the instrumental FWHM over their combined esd |

A component that comes out below zero is set to zero and flagged rather than
carried as a negative width, which noise can easily do to the Gaussian part of
a size broadened peak; its esd is then that of the difference it was clipped
from.

### 21.7 What a breadth follows with angle

Three things broaden a line and each has its own angular dependence, which is
what lets them be told apart. Size broadening goes as one over cos theta,
strain as tan theta, and a spread of specimen heights across the irradiated
surface as cos theta, so a height spread narrows as the angle rises where the
other two widen.

`height_spread_breadth(two_theta, delta_s_mm, radius_mm)` is that third one. A
specimen displaced by s shifts a line by minus two s cos theta over R radians,
so a surface whose heights spread uniformly over `delta_s_mm` spreads the
shift over twice that, and that spread is both the FWHM and the integral
breadth of the broadening it adds. It raises `ValueError` when `delta_s_mm` is
negative, `radius_mm` is not positive, or a position lies outside 0 to 180
degrees. It is worth knowing about because a pressed pellet can have one and
the ground standard of Section 3.1 cannot, so it is broadening the correction
of Section 21.6 will not have removed.

`fit_breadth_models(two_theta, breadth, esd, wavelength, radius_mm,
k=BREADTH_MODEL_K)` fits the corrected breadths with four models and returns a
`BreadthModels`. Three have one parameter each, `size`, `strain` and
`height`, each a breadth proportional to its own dependence; the fourth,
`size+height`, adds the size and height breadths in quadrature and is fitted
in the squares of the two terms, each held at zero or above, which keeps the
derivatives finite when either vanishes. Every fit weights by the breadth esds
and is judged by its reduced chi squared on the breadths themselves, so the
four compare directly. It needs at least `MIN_BREADTH_MODEL_POINTS`, three,
breadths.

`BreadthModelFit` carries one such fit: the `model` name, the `size` in the
units of the wavelength, the `strain`, the `delta_s` in millimetres, each with
its esd and each `None` when the model lacks it or it fits at zero, the
`chi_squared`, `degrees_of_freedom` and `reduced_chi_squared`, the
`normalised_residuals`, the three terms in radians, and `at_bound`, which
marks a two parameter fit that has fallen back onto one term. Its
`breadth(two_theta)` method evaluates the fitted model. `BreadthModels` holds
the four by name in `fits` and offers `preferred`, the one with the lowest
reduced chi squared.

Read `preferred` with the caution the docstring asks for. Size and height
breadths go as one over cos theta and cos theta, both close to flat over a
narrow range of angle, so the two parameter fit is strongly correlated and its
esds are large unless the reflections span widely. A model winning on reduced
chi squared over a short range of angle is not the same as a model being
right.

### 21.8 Fitting one reflection of a standard

The script fits the strongest reflection of the example standard, prints every
field the fit carries, and then does the same thing at an angle where there is
no reflection at all, which is what a caller has to be able to tell apart.

The whole of `profile_fit.py`:

```python
from xrdkit.broadening import (
    doublet_gaps,
    fit_profile,
    integral_breadth,
    kalpha2_position,
    pseudo_voigt_components,
    pseudo_voigt_from_components,
)
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks

# Edit these lines for each new standard. Nothing below needs changing.
SCAN_FILE = "data/standards/lab6.xrdml"
WINDOW = (10.0, 98.0)
REFLECTION = 1
EMPTY_AT = 58.0

scan = read_scan(SCAN_FILE)
peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=WINDOW))
peak = peaks[REFLECTION]
print(
    f"{len(peaks)} reflections; fitting the one found at {peak.two_theta:.3f} degrees"
)
print(f"its found width is {peak.fwhm:.4f} degrees, which is only the guess")

fit = fit_profile(scan.two_theta, scan.intensity, peak.two_theta, peak.fwhm)
print(f"two_theta    {fit.two_theta:.4f} +/- {fit.esd_two_theta:.4f} degrees")
print(f"fwhm         {fit.fwhm:.4f} +/- {fit.esd_fwhm:.4f} degrees")
print(f"eta          {fit.eta:.3f} +/- {fit.esd_eta:.3f}")
print(f"asymmetry    {fit.asymmetry:.3f} +/- {fit.esd_asymmetry:.3f}")
print(f"area         {fit.area:.1f} +/- {fit.esd_area:.1f} counts degrees")
print(f"background   {fit.background:.1f} counts, slope {fit.slope:.1f} per degree")
print(f"kalpha2      {fit.kalpha2_two_theta:.4f} degrees")
print(f"window       {fit.window[0]:.3f} to {fit.window[1]:.3f}, {fit.n_points} points")
print(
    f"reduced chi squared {fit.reduced_chi_squared:.3f}, Rwp {100 * fit.r_wp:.2f} per cent"
)
print(
    f"converged {fit.converged}, truncated {fit.truncated}, fitted {fit.fitted.shape}"
)

print(f"kalpha2_position of the found centre {kalpha2_position(peak.two_theta):.4f}")
gaussian, lorentzian = pseudo_voigt_components(fit.fwhm, fit.eta)
print(f"components   G {gaussian:.4f}, L {lorentzian:.4f} degrees")
width, mixing = pseudo_voigt_from_components(gaussian, lorentzian)
print(f"back again   fwhm {width:.4f} degrees, eta {mixing:.3f}")
print(f"integral breadth {integral_breadth(gaussian, lorentzian):.4f} degrees")
others = [other.two_theta for other in peaks if other is not peak]
below, above = doublet_gaps(peak.two_theta, others)
print(f"clear either side {below:.3f} and {above:.3f} degrees")

empty = fit_profile(scan.two_theta, scan.intensity, EMPTY_AT, peak.fwhm)
print(f"at {EMPTY_AT:.0f} degrees, where there is no reflection:")
print(f"  converged {empty.converged}, Rwp {100 * empty.r_wp:.2f} per cent")
print(f"  fwhm {empty.fwhm:.4f} +/- {empty.esd_fwhm:.4f} degrees")
print(f"  area {empty.area:.1f} +/- {empty.esd_area:.1f} counts degrees")
```

It prints the found peak, the fit, the conversions, the clear space around the
doublet, and then the fit to nothing.

```text
14 reflections; fitting the one found at 30.365 degrees
its found width is 0.0971 degrees, which is only the guess
two_theta    30.3678 +/- 0.0008 degrees
fwhm         0.0781 +/- 0.0014 degrees
eta          0.626 +/- 0.026
asymmetry    -0.364 +/- 0.015
area         1929.7 +/- 17.3 counts degrees
background   612.3 counts, slope 34.7 per degree
kalpha2      30.4451 degrees
window       29.396 to 31.395, 93 points
reduced chi squared 7.041, Rwp 5.67 per cent
converged True, truncated False, fitted (93,)
kalpha2_position of the found centre 30.4425
components   G 0.0516, L 0.0430 degrees
back again   fwhm 0.0781 degrees, eta 0.626
integral breadth 0.1039 degrees
clear either side 8.971 and 6.976 degrees
at 58 degrees, where there is no reflection:
  converged True, Rwp 3.95 per cent
  fwhm 0.0115 +/- 0.2145 degrees
  area 1.2 +/- 9.4 counts degrees
```

Start with the width. The finder gave 0.0971 degrees and the fit gives 0.0781
plus or minus 0.0014. The difference is the satellite: at 30 degrees the K
alpha 2 line sits 0.077 degrees above its parent, well inside the half
maximum, so the raw width is the width of the pair and the fitted width is the
width of the K alpha 1 line alone. That is the whole reason this module
exists: the raw width is a quarter larger than the real one, in a quantity
that later gets turned into a crystallite size.

The asymmetry is minus 0.364, a genuinely lopsided line with the tail to low
angle, which is axial divergence and is what the split shape is for. The
mixing parameter of 0.626 says the line is more Lorentzian than Gaussian, and
the components bear that out: a Gaussian FWHM of 0.0516 and a Lorentzian of
0.0430 degrees, which `pseudo_voigt_from_components` turns straight back into
the 0.0781 and 0.626 they came from. The integral breadth of the same Voigt is
0.1039 degrees, a third as much again as the FWHM, as it is for any profile
with appreciable Lorentzian content.

The reduced chi squared of 7.0 is worth a word, because it looks alarming and
is not a reason to reject this fit. Poisson weights on a strong reflection of
a standard measure the counting statistics exactly, and everything the split
pseudo-Voigt does not describe, which on a sharp line is most of the true peak
shape, lands in the residual. An Rwp of 5.7 per cent on a reflection of this
height is a good fit of a shape that is not quite the right shape. What the
number does mean is that the esds are scaled by it, so the esd of 0.0014
degrees on the width is already widened to account for the misfit.

Now the last three lines. At 58 degrees, where the standard has no reflection,
the fit still converges, and it reports a width of 0.0115 plus or minus 0.2145
degrees and an area of 1.2 plus or minus 9.4 counts times degrees. `converged`
is `True` and the answer is nothing at all: both the width and the area are
consistent with zero and their esds are larger than their values. A script
that filters on `converged` alone will collect results like this one and
average them into whatever comes next. Filter on the esds as well, and on how
far the position moved from where the peak was found, which Section 23.6 does.

## 22. instrument

`xrdkit.instrument` is the two halves of `xrdkit instrument`. Section 3 is the
command, what the instrument parameter file holds and why a standard is
needed; this section is the functions, and it points back rather than
repeating any of it.

There are seven public names: two constants, the two functions, the frozen
dataclass each returns, and a reader for the K alpha 2 wavelength of an
`.xrdml`.

The two functions are separate on purpose. `fit_instrument_widths` fits the
reflections of a standard and the Caglioti relation through their widths,
which takes a second or two and needs no GSAS-II at all. `refine_instrument`
takes that fit as the starting point of a GSAS-II refinement, which measures
the Lorentzian terms and the asymmetry rather than guessing them and writes
the file every later refinement reads. Splitting them means the width fit,
which is quick and answerable on its own, can be looked at before the
refinement is run, and it means a script that only wants the instrumental
resolution function never has to have GSAS-II installed.

### 22.1 fit_instrument_widths and WidthFit

`fit_instrument_widths(scan, window=WIDTH_WINDOW)` fits the widths of a
standard's reflections to the Caglioti relation and returns a `WidthFit`.

Every peak found inside `window` has its K alpha 2 satellites excluded and is
then fitted as a split pseudo-Voigt doublet by `fit_profile`, and the widths
of the fits that converged go into `fit_caglioti` weighted by one over the
square of each width's own esd. `WIDTH_WINDOW` is ten to ninety-eight degrees:
below the lower limit the reflections of a standard are few and asymmetric,
and above the upper one the doublet is wide enough that the fit is about the
splitting rather than the instrument. Section 3.3 is the `--window` option
that changes it.

It raises `ValueError` when the scan carries no wavelength, or when fewer than
three reflections converged, three being the least a three parameter relation
can be fitted to at all.

`WidthFit` is frozen, so nothing downstream can alter a measured instrument by
accident.

| Field | What it holds |
| --- | --- |
| `n_peaks` | how many reflections the fit used |
| `u`, `v`, `w` and their esds | the Caglioti parameters in degrees squared of two theta |
| `rms` | the root mean square of the FWHM residuals, in degrees |
| `wavelength` | the scan's K alpha 1 wavelength, in angstroms |
| `window` | the two theta range the widths were fitted over |
| `two_theta`, `fwhm` | the reflections used and their fitted widths, as tuples |
| `covariance` | the 3 by 3 covariance of U, V and W, as a tuple of tuples |

`fit.caglioti` is the same fit as a `Caglioti` of Section 21.4, which is what
carries the `fwhm` and `fwhm_esd` methods, and what
`xrdkit.gsas2.write_instprm` takes. A script wanting the instrumental width
at an arbitrary angle asks for `fit.caglioti.fwhm(angle)`.

One thing this fit does not give is the instrumental mixing parameter.
`correct_broadening` needs `eta_inst` beside `fwhm_inst`, and nothing in
`WidthFit` or `Caglioti` carries one: the Caglioti relation is about widths
only. A script has to get the instrumental `eta` from the standard's own
profile fits, which is what Section 23.5 does and says.

### 22.2 refine_instrument and InstrumentRefinement

`refine_instrument(fit, scan_file, cif, phase, cell, folder, stem, work,
install=None, radius_mm=None)` refines the instrument parameters in GSAS-II,
starting from a `WidthFit`.

The width fit is written as `<stem>_start.instprm` in `folder` and is what the
refinement starts from. The standard's phase is held at `cell`, with the
crystallite size set large and the microstrain at zero, `STANDARD_BROADENING`,
because a standard is ground and sieved to contribute no broadening of its own
and leaving those at the GSAS-II defaults would let the instrument terms
absorb the difference. The stages are `standard_stages` without the cell
stage, since holding the certified cell is the point of using a standard, and
each runs `CYCLES`, ten. The file GSAS-II exports is copied to
`<stem>.instprm` in `folder`, and that is the file every later refinement
reads.

`REFINED_KEYS` is the order the parameters are freed in: `Zero`, `U`, `V`,
`W`, `X`, `Y` and `SH/L`, the zero, then the Gaussian width terms, then the
Lorentzian terms, then the axial divergence.

`InstrumentRefinement`, also frozen, carries `parameters`, a dict of each of
those keys to its refined `value` and `esd`; `rwp` and `gof` from the last
stage; `stages`, the names in order; `instprm` and `start_instprm`, the two
paths; and `result`, the whole job result for a caller that wants more. It
raises `Gsas2Error` when GSAS-II failed and `KeyError` when the run left no
final instrument parameters.

This guide does not run it. The refinement writes into `data/standards` and
under `work/`, which belong to the command, and the record of what it came to
on the example standard is already in Section 3.2: the five stages, an Rwp of
7.160 per cent, a goodness of fit of 1.909, and every refined parameter with
its esd, including the negative `Y` that section says what to do about. The
refined `U`, `V` and `W` printed there are in the centidegrees squared GSAS-II
works in, which is why they look nothing like the ones the width fit gives.

### 22.3 kalpha2_wavelength

`kalpha2_wavelength(path)` returns the K alpha 2 wavelength an `.xrdml`
records, or `None`. Any other suffix gives `None`, as does a file that cannot
be parsed or that carries no such element, so a two or three column text
pattern always answers `None` and is taken as K alpha 1 only unless the
project file says otherwise. It is how `xrdkit instrument` fills in the
`wavelength` pair of the instrument table it appends, the one Section 3.2
prints.

### 22.4 The width fit of the example standard

The script is the first half of `xrdkit instrument` and nothing else, run on
the same standard scan with the same default window, so the numbers can be
read straight against Section 3.2.

The whole of `instrument_widths.py`:

```python
from xrdkit.instrument import WIDTH_WINDOW, fit_instrument_widths
from xrdkit.io import read_scan

# Edit these lines for each new standard. Nothing below needs changing.
STANDARD_FILE = "data/standards/lab6.xrdml"
AT = (25.0, 40.0, 60.0, 80.0)

scan = read_scan(STANDARD_FILE)
fit = fit_instrument_widths(scan)
print(
    f"peaks   {fit.n_peaks} fitted from {fit.window[0]:.0f} to {fit.window[1]:.0f} degrees"
)
print(f"U = {fit.u:.4f}, V = {fit.v:.4f}, W = {fit.w:.4f} degrees squared")
print(f"rms of the width fit {fit.rms:.5f} degrees")
print(f"esds    {fit.esd_u:.4f}, {fit.esd_v:.4f}, {fit.esd_w:.4f}")
print(f"wavelength {fit.wavelength} angstrom, default window {WIDTH_WINDOW}")

print("  two_theta  measured  Caglioti  difference")
for angle, width in zip(fit.two_theta, fit.fwhm):
    modelled = fit.caglioti.fwhm(angle)
    print(f"  {angle:9.3f}  {width:8.4f}  {modelled:8.4f}  {width - modelled:+10.4f}")

print("the resolution function away from the measured reflections")
for angle in AT:
    print(
        f"  {angle:5.1f} degrees  {fit.caglioti.fwhm(angle):.4f}"
        f" +/- {fit.caglioti.fwhm_esd(angle):.4f} degrees"
    )
```

It prints the fit, then every reflection with the width the relation puts
there, then the resolution function at four angles no reflection sits at.

```text
peaks   14 fitted from 10 to 98 degrees
U = 0.0107, V = -0.0122, W = 0.0085 degrees squared
rms of the width fit 0.00254 degrees
esds    0.0017, 0.0020, 0.0006
wavelength 1.540598 angstrom, default window (10.0, 98.0)
  two_theta  measured  Caglioti  difference
     21.346    0.0852    0.0811     +0.0041
     30.368    0.0781    0.0772     +0.0009
     37.424    0.0728    0.0747     -0.0019
     43.488    0.0720    0.0730     -0.0009
     48.938    0.0710    0.0718     -0.0008
     53.966    0.0709    0.0710     -0.0001
     63.196    0.0721    0.0709     +0.0012
     67.525    0.0728    0.0715     +0.0014
     71.722    0.0732    0.0725     +0.0008
     75.820    0.0775    0.0740     +0.0035
     79.846    0.0823    0.0759     +0.0063
     83.819    0.0793    0.0784     +0.0009
     87.769    0.0786    0.0815     -0.0029
     95.648    0.0884    0.0897     -0.0014
the resolution function away from the measured reflections
   25.0 degrees  0.0794 +/- 0.0013 degrees
   40.0 degrees  0.0739 +/- 0.0007 degrees
   60.0 degrees  0.0708 +/- 0.0008 degrees
   80.0 degrees  0.0760 +/- 0.0010 degrees
```

The first three lines are the first three lines of Section 3.2, figure for
figure: fourteen reflections from ten to ninety-eight degrees, U of 0.0107, V
of minus 0.0122, W of 0.0085 degrees squared and an rms of 0.00254 degrees.
They have to be, because the command calls this function and prints what it
returns.

The table underneath is what the rms of 0.00254 degrees is made of. The widths
run from 0.0709 degrees at 54 degrees up to 0.0884 at 95, a shallow bowl with
its minimum near 60 degrees, which is the shape a negative V gives. Most
reflections sit within a thousandth of a degree of the relation. Two do not:
the first reflection at 21.3 degrees is 0.0041 high and the one at 79.8
degrees is 0.0063 high. The first is the most asymmetric reflection in the
scan, where the split shape is working hardest; the second is the weakest, at
four and a half per cent of the strongest, and the same reflection whose
mixing parameter comes out at 0.31 against 0.6 everywhere else, which is a fit
with too little signal rather than a real feature of the instrument.

The last four lines are the point of fitting a relation rather than tabulating
widths. No sample reflection will sit at 21.346 or 30.368 degrees, and the
instrumental width at 25, 40, 60 and 80 degrees is what a correction actually
needs. The esds there, one to two thousandths of a degree, are what
`correct_broadening` combines with the sample's own esd, which on a good
sample reflection is three to seven thousandths, to decide whether a
reflection is broadened at all. Two of that combination is the threshold, so a
reflection on this instrument has to be the best part of a hundredth of a
degree wider than the standard before anything can be said about it at all.
Section 23 is what follows from that.

## 23. sizestrain

`xrdkit.sizestrain` turns sample breadths into a crystallite size and a
microstrain. It is the one module in the kit that no command touches, and that
is deliberate rather than an omission: every step from a peak width to a size
rests on a judgement about which reflections to use that a command cannot make
for you. This section is therefore the longest of the three and carries a
whole analysis rather than a demonstration.

Every function here takes breadths that already have the instrument taken out,
in degrees of two theta, such as `correct_broadening` of Section 21.6 gives,
and works in radians internally. Sizes come out in the units of the
wavelength, angstroms in practice. There are eight public names: three
dataclasses, four functions and one alias, together with `SCHERRER_K`, 0.9,
the constant for a FWHM, and `INTEGRAL_BREADTH_K`, 1.0, the constant for an
integral breadth.

### 23.1 scherrer_size and scherrer_size_integral

`scherrer_size(fwhm_deg, two_theta, wavelength, k=SCHERRER_K,
esd_fwhm=0.0)` returns the Scherrer size and its esd, as a pair, from one
breadth: D equals K lambda over beta cos theta. Scalars give floats and arrays
give arrays.

The whole breadth is put down to size, so D is an apparent size along the
normal to the reflecting planes, and a lower bound on the real size if the
line is also strain broadened. The esd comes from that of the breadth alone,
the position being taken as exact. It raises `ValueError` when a breadth, the
wavelength or `k` is not positive, a position lies outside 0 to 180 degrees,
or an esd is negative.

`scherrer_size_integral(breadth_deg, two_theta, wavelength,
k=INTEGRAL_BREADTH_K, esd_breadth=0.0)` is the same function with the constant
for an integral breadth. Area over height is the breadth for which K equal to
1 gives the volume weighted column length whatever the crystallite shape,
which is why it has a name of its own rather than being a constant a caller
has to remember.

Match the constant to the breadth. A FWHM with K of 1, or an integral breadth
with K of 0.9, is an error of ten per cent in the size and nothing warns about
it.

### 23.2 williamson_hall and WilliamsonHall

`williamson_hall(two_theta, breadth_deg, esd, wavelength, k=SCHERRER_K)` fits
beta cos theta equals K lambda over D plus epsilon times 4 sin theta, by
weighted linear least squares, and returns a `WilliamsonHall`.

Size broadening goes as one over cos theta and strain broadening as tan theta,
so plotting beta cos theta against 4 sin theta puts the size in the intercept
and the strain in the slope. Both are free: the line is never forced through
the origin, so a crystallite too large to broaden anything shows up as an
intercept consistent with zero rather than being assumed away. Each point is
weighted by one over the square of its own esd in the plotted quantity. It
needs at least `MIN_WILLIAMSON_HALL_POINTS`, three, so that a two parameter
line is left a residual.

It raises `ValueError` when the arrays differ in length, there are fewer than
three points, a position lies outside 0 to 180 degrees, a breadth is not
finite, an esd is not positive, or the wavelength or `k` is not positive.

| Field | What it holds |
| --- | --- |
| `size`, `esd_size` | K lambda over the intercept, in the units of the wavelength, and its esd. `None` when the intercept is not positive, since no size then fits |
| `strain`, `esd_strain` | the slope, the apparent strain of beta equals 4 epsilon tan theta, and its esd |
| `intercept`, `esd_intercept`, `slope`, `esd_slope` | the line itself, the intercept in radians |
| `covariance` | the 2 by 2 covariance of intercept and slope, scaled by the reduced chi squared |
| `reduced_chi_squared` | the weighted residual per degree of freedom |
| `n_reflections`, `k`, `wavelength` | what went in |
| `x`, `y`, `esd_y` | 4 sin theta, beta cos theta in radians, and its esd, per reflection |
| `residuals`, `normalised_residuals` | y less the line, in radians, and the same over `esd_y` |

`slope_significance` is the slope over its esd, and `line(x)` and
`line_esd(x)` evaluate the fitted line and its esd, which is what draws the
band in Section 23.8.

The strain this returns is an upper limit on the local lattice strain, not a
measurement of it, because everything that broadens like tan theta lands in
the slope. The size is the more fragile of the two. The data almost never
reach x of zero, so the intercept is an extrapolation across a gap, and its
esd is what says whether the extrapolation was worth making. Read
`esd_intercept` against `intercept` before quoting a size at all.

### 23.3 component_size_strain and ComponentSizeStrain

`component_size_strain(two_theta, lorentzian_deg, esd_lorentzian,
gaussian_deg, esd_gaussian, wavelength, k=INTEGRAL_BREADTH_K)` is the cross
check of de Keijser, Langford, Mittemeijer and Vogels on Williamson-Hall. Size
broadening is taken to be Lorentzian and strain broadening Gaussian, so the
size comes from the weighted mean of beta L cos theta and the strain from the
weighted slope of beta G through 4 tan theta, each component fitted against
its own dependence alone with no constant term.

The breadths it wants are integral breadths of the two components, which are
pi over 2 times the Lorentzian FWHM and the square root of pi over log 2, over
2, times the Gaussian FWHM; `BroadeningCorrection` reports the FWHM of each,
so a caller converts. A component the correction clipped to zero is a
measurement of zero, with its esd, and is used as such. It needs at least
`MIN_COMPONENT_POINTS`, two.

`ComponentSizeStrain` carries `size` and `esd_size`, `size_term` and its esd,
`strain` and its esd, a `size_reduced_chi_squared` and a
`strain_reduced_chi_squared` saying how well each dependence held, the counts
and constants that went in, and the two sets of residuals.

The worked example below does not use it. With six resolved reflections over
thirty degrees of two theta, a second separation of the same two effects from
the same six numbers adds a second answer rather than a check on the first,
and the honest thing is to say so. It earns its place where the reflections
are many and span widely, which a benchtop scan of a weakly broadened sample
is not.

### 23.4 resolution_limit and ResolutionLimit

`resolution_limit(two_theta, fwhm_inst, esd_inst, esd_obs, wavelength,
k=SCHERRER_K, significance=DEFAULT_SIGNIFICANCE, convolution="quadrature")`
gives the largest crystallite size a reflection could have told apart from the
instrument. This is the function that keeps an analysis honest, and it should
be run before the sizes are, not after.

`correct_broadening` counts a reflection as resolved when its FWHM exceeds the
instrumental one by `significance` combined esds. The sample breadth that just
reaches that threshold depends on how the two profiles combine. Under
`"quadrature"` they add as Gaussians, and the breadth is the root of the
difference of squares, which is the larger breadth and so the smaller size:
the conservative bound, since a sample breadth of any shape below it could
have gone unseen. Under `"linear"` they add as Lorentzians and the breadth is
the threshold itself, the smaller breadth and the larger size, which is the
bound to use when the sample broadening is known to be Lorentzian. The breadth
is then turned into a size by `scherrer_size`.

`ResolutionLimit` carries `threshold`, the smallest excess that counts, in
degrees; `breadth`, the sample FWHM that would produce it; `size`, the
Scherrer size of that breadth; and `convolution`, `significance` and `k`, so
that a row written to a file records the assumption it was made under. Each is
a scalar or an array, following the input.

Read `size` as a lower bound. A reflection that shows no measurable broadening
has crystallites larger than its own limit, and says nothing at all about how
much larger. It raises `ValueError` when `convolution` is unknown, an
instrumental width or `significance` is not positive, both esds of a
reflection are zero, or an esd is negative.

### 23.5 The instrument side of the analysis

The rest of this section is one analysis of the example pellet against the
example standard, in four blocks of one script. This first block measures the
instrument, and it is the width fit of Section 22.4 with one thing added.

`correct_broadening` needs an instrumental mixing parameter as well as an
instrumental width, and Section 22.1 says that `WidthFit` carries no such
thing: the Caglioti relation describes widths and says nothing about shape. So
the standard's reflections are fitted a second time here, with `fit_profile`
directly, and their mixing parameters are kept alongside their angles. A
sample reflection then takes the instrumental `eta` interpolated between the
two standard reflections either side of it. That interpolation is the one step
in this analysis that no xrdkit function provides, and it is the weakest link
in the chain; it is also why `ESD_ETA_INST` is set to a frank 0.02 in the
settings rather than to something the data justify.

Start a new file named `size_strain.py`.

```python
import numpy as np

from xrdkit.broadening import correct_broadening, doublet_gaps, fit_profile
from xrdkit.cell import Cell
from xrdkit.indexing import generate_reflections, index_and_refine
from xrdkit.instrument import fit_instrument_widths
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
STANDARD_FILE = "data/standards/lab6.xrdml"
STEM = "pellet_a"
START_CELL = Cell.tetragonal(12.45, 3.94)
SPACE_GROUP = "P4bm"
WINDOW = (20.0, 90.0)
STANDARD_WINDOW = (10.0, 98.0)
MIN_CLEAR = 0.3
ESD_ETA_INST = 0.02

standard = read_scan(STANDARD_FILE)
widths = fit_instrument_widths(standard, window=STANDARD_WINDOW)
resolution = widths.caglioti
standard_angles, standard_etas = [], []
for peak in exclude_kalpha2(find_peaks(standard, two_theta_range=STANDARD_WINDOW)):
    profile = fit_profile(
        standard.two_theta, standard.intensity, peak.two_theta, peak.fwhm
    )
    if profile.converged:
        standard_angles.append(profile.two_theta)
        standard_etas.append(profile.eta)
print(f"instrument: {widths.n_peaks} reflections, rms {widths.rms:.5f} degrees")
print(f"U = {widths.u:.4f}, V = {widths.v:.4f}, W = {widths.w:.4f} degrees squared")
print(
    f"instrumental eta runs {min(standard_etas):.2f} to {max(standard_etas):.2f}"
    f" over {standard_angles[0]:.1f} to {standard_angles[-1]:.1f} degrees"
)
```

It prints the instrument.

```text
instrument: 14 reflections, rms 0.00254 degrees
U = 0.0107, V = -0.0122, W = 0.0085 degrees squared
instrumental eta runs 0.31 to 0.75 over 21.3 to 95.6 degrees
```

The width fit is the one of Section 22.4 and Section 3.2. The mixing parameter
runs from 0.31 to 0.75 across the scan, and Section 22.4 has already said that
the 0.31 belongs to the weakest reflection in the standard and is a fit with
too little signal rather than the instrument changing shape. Every sample
reflection used below sits between 25 and 60 degrees, and the instrumental
mixing parameters the interpolation hands them run from 0.59 to 0.69, so it is
working across a gently sloping stretch of the standard and the outlier at
79.8 degrees never enters it at all.

### 23.6 The reflections worth measuring

A width belongs to one reflection or to nothing. The selection here is three
tests, in order, and each of them throws away reflections that a less careful
script would have measured.

First, overlap. The peaks are indexed against the refined cell as Section 16
indexes them, and every calculated reflection of that cell, not merely the
peaks the finder reported, is offered to `doublet_gaps`. A peak keeps its
place only if the nearest other calculated line is `MIN_CLEAR` degrees clear
of its whole doublet. This is the test that catches a peak which is really two
reflections the instrument never resolved, which no amount of profile fitting
will separate and which would otherwise be reported as a broad line and turned
into a small crystallite.

Second, the refit. Each surviving peak is fitted with `fit_profile`, and a fit
that did not converge or that moved the position by more than the found peak's
own FWHM is rejected, which is the rule of Section 7.2 and the warning of
Section 21.8.

Third, significance. `correct_broadening` decides whether what is left is
broadened at all, at the default two combined esds.

Continues `size_strain.py`. Add these lines at the end of the file.

```python
scan = read_scan(SCAN_FILE)
peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=WINDOW))
indexed, cell_fit = index_and_refine(
    peaks,
    start_cell=START_CELL,
    wavelength=scan.wavelength,
    space_group=SPACE_GROUP,
)
calculated = [
    reflection.two_theta + cell_fit.zero_offset
    for reflection in generate_reflections(
        cell_fit.cell,
        scan.wavelength,
        WINDOW[1] + 2.0,
        WINDOW[0] - 2.0,
        space_group=SPACE_GROUP,
    )
]
print(f"{len(peaks)} peaks, {len(calculated)} calculated reflections")

rows = []
for entry in indexed:
    if not entry.is_indexed:
        continue
    own = entry.reflection.two_theta + cell_fit.zero_offset
    others = [angle for angle in calculated if abs(angle - own) > 1e-9]
    if min(doublet_gaps(entry.peak.two_theta, others)) < MIN_CLEAR:
        continue
    profile = fit_profile(
        scan.two_theta, scan.intensity, entry.peak.two_theta, entry.peak.fwhm
    )
    moved = abs(profile.two_theta - entry.peak.two_theta)
    if not profile.converged or moved > entry.peak.fwhm:
        print(f"  refit at {entry.peak.two_theta:.3f} rejected, moved {moved:.3f}")
        continue
    fwhm_inst = float(resolution.fwhm(profile.two_theta))
    esd_inst = float(resolution.fwhm_esd(profile.two_theta))
    eta_inst = float(np.interp(profile.two_theta, standard_angles, standard_etas))
    correction = correct_broadening(
        profile.fwhm,
        profile.eta,
        fwhm_inst,
        eta_inst,
        profile.esd_fwhm,
        profile.esd_eta,
        esd_inst,
        ESD_ETA_INST,
    )
    rows.append(
        (entry.reflection.hkl, profile, fwhm_inst, esd_inst, eta_inst, correction)
    )

print(f"{len(rows)} reflections clear of their neighbours by {MIN_CLEAR} degrees")
print("        hkl  two_theta  observed  instrument  excess  sample breadth")
for hkl, profile, fwhm_inst, _, _, correction in rows:
    sample = (
        "unresolved"
        if correction.unresolved
        else f"{correction.fwhm:.4f} +/- {correction.esd_fwhm:.4f}"
    )
    print(
        f"  {hkl!s:>9}  {profile.two_theta:9.3f}  {profile.fwhm:8.4f}"
        f"  {fwhm_inst:10.4f}  {correction.excess:6.2f}  {sample}"
    )
```

It prints what each test left.

```text
39 peaks, 168 calculated reflections
  refit at 26.921 rejected, moved 1.030
10 reflections clear of their neighbours by 0.3 degrees
        hkl  two_theta  observed  instrument  excess  sample breadth
  (3, 2, 0)     25.890    0.0926      0.0791    2.01  0.0436 +/- 0.0125
  (2, 1, 1)     27.904    0.1094      0.0782    2.88  0.0509 +/- 0.0212
  (4, 0, 0)     28.764    0.0631      0.0779   -0.13  unresolved
  (4, 1, 0)     29.670    0.0889      0.0775    1.09  unresolved
  (3, 2, 1)     34.666    0.0868      0.0756    3.08  0.0450 +/- 0.0077
  (5, 2, 0)     39.013    0.0907      0.0742    1.54  unresolved
  (5, 3, 0)     42.359    0.0881      0.0732    2.38  0.0490 +/- 0.0150
  (6, 0, 1)     49.722    0.0972      0.0716    2.02  0.0842 +/- 0.0335
  (6, 3, 1)     54.701    0.0878      0.0710    2.63  0.0729 +/- 0.0175
  (8, 0, 0)     59.336    0.0541      0.0707   -0.86  unresolved
```

Thirty nine peaks go in and ten come out of the first two tests. That is the
first honest number in this analysis: three quarters of the pattern of a
tungsten bronze with a twelve angstrom axis is too crowded to measure a width
on, and the reflections that survive are all below sixty degrees, because the
higher the angle the more reflections there are per degree.

The rejected refit is the peak found at 26.921 degrees, which moved 1.030
degrees when it was fitted. That is the peak Section 15.8 marked as
unexplained and Section 18 deleted from its peak list, and here it fails for a
third reason: the window `fit_profile` chooses is ten of its starting widths
wide, and inside that window there is a much stronger neighbour for the fit to
walk to. Without the shift test this analysis would have measured that
neighbour twice and called the second measurement (201).

Of the ten that remain, six clear the significance test and four do not. Look
at the four. Two, (400) at 28.764 and (800) at 59.336 degrees, are actually
narrower than the instrument, with an excess of minus 0.13 and minus 0.86
esds. A sample cannot be sharper than the instrument; what this means is that
those two fitted widths are a little low and their esds are honest. The other
two, (410) and (520), sit at 1.09 and 1.54 esds, broader than the instrument
but not by enough to be sure of. The module reports them as `unresolved`
rather than handing back a small positive breadth with an esd of nearly the
same size, and that refusal is the whole point of the `significance` argument.

### 23.7 Scherrer, Williamson-Hall and what the scan supports

Continues `size_strain.py`. Add these lines at the end of the file.

```python
from xrdkit.sizestrain import (
    SCHERRER_K,
    resolution_limit,
    scherrer_size,
    williamson_hall,
)

resolved = [row for row in rows if not row[5].unresolved]
angles = np.array([row[1].two_theta for row in resolved])
breadths = np.array([row[5].fwhm for row in resolved])
esd_breadths = np.array([row[5].esd_fwhm for row in resolved])
sizes, esd_sizes = scherrer_size(
    breadths, angles, scan.wavelength, esd_fwhm=esd_breadths
)
limits = resolution_limit(
    np.array([row[1].two_theta for row in rows]),
    np.array([row[2] for row in rows]),
    np.array([row[3] for row in rows]),
    np.array([row[1].esd_fwhm for row in rows]),
    scan.wavelength,
)
print(f"{len(resolved)} of {len(rows)} reflections resolved at {SCHERRER_K} Scherrer K")
print("        hkl  two_theta   Scherrer size       resolution limit")
for (hkl, profile, *_), size, esd in zip(resolved, sizes, esd_sizes):
    index = [row[1].two_theta for row in rows].index(profile.two_theta)
    print(
        f"  {hkl!s:>9}  {profile.two_theta:9.3f}  {size:6.0f} +/- {esd:5.0f} A"
        f"  {float(limits.size[index]):10.0f} A"
    )
print(f"Scherrer sizes run {sizes.min():.0f} to {sizes.max():.0f} angstrom")
print(
    f"resolution limits run {limits.size.min():.0f} to {limits.size.max():.0f}"
    f" angstrom, {limits.convolution}, {limits.significance:.0f} esds"
)

fit = williamson_hall(angles, breadths, esd_breadths, scan.wavelength)
print(f"Williamson-Hall on {fit.n_reflections} reflections")
print(f"  intercept {fit.intercept:.3e} +/- {fit.esd_intercept:.3e} radians")
print(f"  size      {fit.size:.0f} +/- {fit.esd_size:.0f} angstrom")
print(f"  strain    {fit.strain:.2e} +/- {fit.esd_strain:.2e}")
print(f"  slope over its esd {fit.slope_significance:.2f}")
print(f"  reduced chi squared {fit.reduced_chi_squared:.3f}")
bound = SCHERRER_K * scan.wavelength / (fit.intercept + 2.0 * fit.esd_intercept)
print(f"  intercept is {fit.intercept / fit.esd_intercept:.2f} esds from zero")
print(f"  so the size is above {bound:.0f} angstrom at two esds, with no upper bound")
```

It prints the sizes reflection by reflection against the resolution limits,
then the line through them.

```text
6 of 10 reflections resolved at 0.9 Scherrer K
        hkl  two_theta   Scherrer size       resolution limit
  (3, 2, 0)     25.890    1871 +/-   539 A        1697 A
  (2, 1, 1)     27.904    1608 +/-   670 A        1318 A
  (3, 2, 1)     34.666    1850 +/-   318 A        2462 A
  (5, 3, 0)     42.359    1738 +/-   533 A        1908 A
  (6, 0, 1)     49.722    1039 +/-   414 A        1339 A
  (6, 3, 1)     54.701    1227 +/-   295 A        2006 A
Scherrer sizes run 1039 to 1871 angstrom
resolution limits run 279 to 2462 angstrom, quadrature, 2 esds
Williamson-Hall on 6 reflections
  intercept 3.064e-04 +/- 2.126e-04 radians
  size      4526 +/- 3140 angstrom
  strain    4.13e-04 +/- 1.67e-04
  slope over its esd 2.48
  reduced chi squared 0.266
  intercept is 1.44 esds from zero
  so the size is above 1896 angstrom at two esds, with no upper bound
```

Read the two columns of the table against each other and not down the page.
The Scherrer sizes run from 1039 to 1871 angstrom, which looks like a
measurement until the resolution limits beside them are read: 1697, 1318,
2462, 1908, 1339 and 2006 angstrom. Four of the six sizes lie below their own
reflection's limit and two lie above it. A size below its limit is a size the
reflection could only just have detected, which is exactly what a breadth at
two or three esds means. These are not six measurements of a crystallite size.
They are six reflections sitting on the edge of what this instrument can see.

The Williamson-Hall fit says the same thing more sharply. The slope, the
strain, comes out at 4.13 times ten to the minus four plus or minus 1.67, two
and a half esds from zero, which is a real if marginal signal. The intercept
is 3.064 times ten to the minus four radians plus or minus 2.126, one and a
half esds from zero, so it is consistent with there being no size broadening
in this pattern at all. The size of 4526 plus or minus 3140 angstrom that
follows from it is not a crystallite size and must not be quoted as one: its
esd is seven tenths of its value, and a quantity known to seventy per cent is
a quantity that has not been measured.

What the scan does support is a bound. Taking the intercept two esds high
gives a size above 1896 angstrom, about 190 nanometres, with no upper bound at
all. That is the sentence to write down. The figure of Section 23.8 shows why
it is the only one available: the six points span 4 sin theta from 0.90 to
1.84 and the intercept sits at zero, so the line is extrapolated back across
almost as much x again as it was measured over, and the band opens out as it
goes.

The reduced chi squared of 0.266 is worth one more word. It is well below one,
which means the six breadths scatter about the line by less than their esds
say they should. That is not a good fit; it is a sign that the esds are
generous, which is what one expects after `correct_broadening` has propagated
four inputs through numerical derivatives and treated them as independent when
two of them are correlated. It is a reason to trust the bound above and not to
sharpen it.

The honest summary of this sample, then, is one line: the crystallites are
larger than about 190 nanometres, and an apparent strain of four parts in ten
thousand is present at two and a half esds. A benchtop scan of a well
crystallised ceramic usually ends here, and a script that reports a crystallite
size from it has not measured one.

### 23.8 The record and the figure

Every number above rests on choices, and a file of results that does not carry
them cannot be checked six months later. The last block writes both: a CSV
with one row per reflection, carrying the inputs beside the outputs, and the
Williamson-Hall figure.

The CSV follows the standing rule every results file in the kit follows, which
Section 3.2 describes for `lab6_instrument.csv`: the last three columns are
the method in one sentence, the date, and the version of xrdkit that wrote the
row. Beside those go the scan and the standard, the hkl, the observed and
instrumental widths and mixing parameters with their esds, the excess in esds,
whether the reflection resolved, the sample breadth, the Scherrer size, the
resolution limit, and the wavelength, Scherrer constant, significance and
clearance the run used, then the Williamson-Hall result on every row. Both the
unresolved reflections and the resolved ones get a row, because which
reflections did not resolve is part of the answer.

The figure is drawn with matplotlib directly rather than through
`plot_pattern` or `plot_stacked`, since no function in `xrdkit.plotting` draws
a Williamson-Hall plot. `apply_style` is called first so the figure matches
the rest of the kit's, a `Figure` is made and an `Axes` added to it, and
`save_figure` writes the pair of files; `pyplot` is not used anywhere, as
Section 15 says it should not be.

Continues `size_strain.py`. Add these lines at the end of the file.

```python
import csv
import datetime
from pathlib import Path

from matplotlib.figure import Figure

from xrdkit import __version__
from xrdkit.plotting import apply_style, save_figure

METHOD = (
    "profile fit of each well separated reflection as a K alpha doublet; "
    "instrumental width from a Caglioti fit to a standard; "
    "Thompson-Cox-Hastings deconvolution; Scherrer and Williamson-Hall"
)
COLUMNS = (
    "scan",
    "standard",
    "h",
    "k",
    "l",
    "two_theta_deg",
    "fwhm_obs_deg",
    "esd_fwhm_obs_deg",
    "eta_obs",
    "fwhm_inst_deg",
    "esd_fwhm_inst_deg",
    "eta_inst",
    "excess_esds",
    "resolved",
    "fwhm_sample_deg",
    "esd_fwhm_sample_deg",
    "scherrer_size_a",
    "esd_scherrer_size_a",
    "resolution_limit_a",
    "wavelength_a",
    "scherrer_k",
    "significance",
    "min_clear_deg",
    "wh_size_a",
    "esd_wh_size_a",
    "wh_strain",
    "esd_wh_strain",
    "wh_reduced_chi_squared",
    "wh_size_lower_bound_a",
    "method",
    "date",
    "xrdkit_version",
)

today = datetime.datetime.now(datetime.UTC).astimezone().date().isoformat()
path = Path(f"results/library/sizestrain_{STEM}.csv")
path.parent.mkdir(parents=True, exist_ok=True)
with path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, COLUMNS)
    writer.writeheader()
    for index, (hkl, profile, fwhm_inst, esd_inst, eta_inst, correction) in enumerate(
        rows
    ):
        row = dict(zip(COLUMNS, [""] * len(COLUMNS)))
        size = esd_size = ""
        if not correction.unresolved:
            place = resolved.index(rows[index])
            size, esd_size = f"{sizes[place]:.0f}", f"{esd_sizes[place]:.0f}"
        row.update(
            scan=SCAN_FILE,
            standard=STANDARD_FILE,
            h=hkl[0],
            k=hkl[1],
            l=hkl[2],
            two_theta_deg=f"{profile.two_theta:.4f}",
            fwhm_obs_deg=f"{profile.fwhm:.5f}",
            esd_fwhm_obs_deg=f"{profile.esd_fwhm:.5f}",
            eta_obs=f"{profile.eta:.4f}",
            fwhm_inst_deg=f"{fwhm_inst:.5f}",
            esd_fwhm_inst_deg=f"{esd_inst:.5f}",
            eta_inst=f"{eta_inst:.4f}",
            excess_esds=f"{correction.excess:.3f}",
            resolved=not correction.unresolved,
            fwhm_sample_deg="" if correction.unresolved else f"{correction.fwhm:.5f}",
            esd_fwhm_sample_deg=(
                "" if correction.unresolved else f"{correction.esd_fwhm:.5f}"
            ),
            scherrer_size_a=size,
            esd_scherrer_size_a=esd_size,
            resolution_limit_a=f"{float(limits.size[index]):.0f}",
            wavelength_a=scan.wavelength,
            scherrer_k=SCHERRER_K,
            significance=limits.significance,
            min_clear_deg=MIN_CLEAR,
            wh_size_a=f"{fit.size:.0f}",
            esd_wh_size_a=f"{fit.esd_size:.0f}",
            wh_strain=f"{fit.strain:.3e}",
            esd_wh_strain=f"{fit.esd_strain:.3e}",
            wh_reduced_chi_squared=f"{fit.reduced_chi_squared:.4f}",
            wh_size_lower_bound_a=f"{bound:.0f}",
            method=METHOD,
            date=today,
            xrdkit_version=__version__,
        )
        writer.writerow(row)
print(f"{path} ({len(rows)} rows, {len(COLUMNS)} columns)")


apply_style()
figure = Figure(figsize=(3.5, 2.6))
axes = figure.add_subplot()
axes.errorbar(
    fit.x,
    1000.0 * fit.y,
    yerr=1000.0 * fit.esd_y,
    fmt="o",
    markersize=3.0,
    color="black",
    linewidth=0.7,
    capsize=2.0,
)
grid = np.linspace(0.0, 1.05 * fit.x.max(), 100)
axes.plot(grid, 1000.0 * fit.line(grid), color="black", linewidth=0.7)
band = 1000.0 * fit.line_esd(grid)
axes.fill_between(
    grid,
    1000.0 * fit.line(grid) - band,
    1000.0 * fit.line(grid) + band,
    color="black",
    alpha=0.12,
    linewidth=0.0,
)
axes.set_xlim(0.0, 1.05 * fit.x.max())
axes.set_xlabel(r"$4\sin\theta$")
axes.set_ylabel(r"$\beta\cos\theta$ / mrad")
axes.set_title(f"{STEM}, Williamson-Hall")
axes.text(
    0.03,
    0.95,
    f"intercept {1000 * fit.intercept:.2f} +/- {1000 * fit.esd_intercept:.2f} mrad\n"
    f"strain {fit.strain:.2e} +/- {fit.esd_strain:.2e}\n"
    f"size > {bound:.0f} angstrom",
    transform=axes.transAxes,
    verticalalignment="top",
    fontsize=5,
)
figure.tight_layout()
print(save_figure(figure, f"figures/library_williamson_hall_{STEM}"))
```

It prints the two files it wrote.

```text
results\library\sizestrain_pellet_a.csv (10 rows, 32 columns)
[WindowsPath('figures/library_williamson_hall_pellet_a.png'), WindowsPath('figures/library_williamson_hall_pellet_a.pdf')]
```

The figure puts the intercept, the strain and the bound in the corner, and
draws the esd band of the fitted line across the whole x range so that the
extrapolation to zero is visible rather than implied. Read it as the picture
of the last two paragraphs: six points with honest error bars, a line with a
real slope, and a band at x equal to zero that reaches from the bound down to
very nearly the axis.

## 24. phases

`xrdkit.phases` is what `xrdkit phases` is built from: the Crystallography
Open Database search, the CIF downloads, the index of what was gathered, the
simulation of a candidate's powder pattern and the scoring of that pattern
against the peaks that were measured. Section 6 is the command, what it found
on the example sample and how to read the ranking, and this section points
there rather than repeating any of it.

The module splits in two along a dependency. Everything that talks to the COD
and everything that keeps a record uses the standard library alone, so it
works in any installation. Everything that simulates a pattern needs pymatgen,
which is the optional `phases` extra, `xrdkit[phases]`, and is imported only
at the moment a pattern is simulated.

### 24.1 The COD: cod_search and CodRecord

`cod_search(elements, exact=True, space_group=None, extra=None)` searches the
COD through its REST interface and returns a list of `CodRecord`.

| Argument | What it does |
| --- | --- |
| `elements` | one to eight element symbols, as `["Sr", "Ba", "Nb", "O"]`. The COD's search form takes at most eight |
| `exact` | `True` returns only entries made of exactly these elements; `False` returns entries holding these and any others |
| `space_group` | keep only entries in this group. An integer is the International Tables number and goes to the COD's own search; a string is a Hermann-Mauguin symbol, matched here after spaces and any trailing setting are stripped |
| `extra` | further COD parameters, such as `{"year": 2006}`, added to the query as given and overriding the rest |

It raises `ValueError` when no element or more than eight are given, or one is
repeated, and `urllib.error.URLError` when the COD cannot be reached.

The symbol form of `space_group` needs care, and the docstring says why: the
COD ignores a symbol in the query itself, so the matching is done here and the
symbol has to be written as the COD writes it, `"P121/c1"` rather than
`"P21/c"`. An entry whose symbol the COD gives in a nonstandard centring, such
as the `X4bm` of Section 6.3, is matched by neither form.

`CodRecord` is one entry as the search returns it: `cod_id`, `formula`,
`space_group` and `space_group_number`, the six cell parameters and the
`volume`, `authors`, `journal`, `year` and `doi`, and `temperature` and
`pressure` of the cell measurement. Anything the COD leaves blank is `None`,
and old entries often carry no temperature at all. Its `filename` property is
the name `cod_fetch` gives its CIF, `<id>.cif`, and `reference()` is a one
line citation built from the authors, year, journal and doi.

### 24.2 cod_fetch, fetch_candidates and write_cif_index

`cod_fetch(cod_id, folder)` downloads one entry's CIF into `folder` as
`<id>.cif`, creating the folder if needed, and returns the path. The file is
written exactly as the COD serves it, and an existing file of the same name is
replaced. It raises `ValueError` when `cod_id` is not a seven digit COD id or
what comes back is not a CIF.

`fetch_candidates(records, folder, pause=time.sleep, pause_s=PHASES_PAUSE_S)`
is the loop around it: every record whose CIF is not already in `folder` is
downloaded, with a pause of `PHASES_PAUSE_S`, half a second, after each one so
that the database is not hammered. It returns only the files it fetched, so a
second run over the same records returns an empty list. That is the behaviour
Section 6.2 describes as the CIFs not being fetched again: a project accrues a
local set of reference CIFs and the network is used only for what is new.
`pause` is there so that a test can pass a callable that does nothing.

`write_cif_index(records, path, notes="")` writes or updates the CSV index of
reference CIFs and returns the path. Its columns are `CIF_INDEX_COLUMNS`: the
file, the source, the identifier, the formula, the space group, the six cell
parameters, the reference and any notes. A `CodRecord` is entered with a
source of `COD` and its COD id as the identifier; a CIF from anywhere else is
given as a mapping keyed by those same columns and brings its own notes.

The updating is the part worth knowing. A row is identified by its source and
identifier, so a record matching a row already in the file replaces it in
place, one matching none is added at the end, and rows no record names are
left alone. The notes of a replaced row survive when the new record brings
none, which is what lets a note written by hand outlive a fresh download. It
raises `ValueError` when a mapping has a column that is not an index column,
or no source or identifier.

None of these three is run in this guide. They reach the network, and the
record of what they returned on the example sample is Section 6.2: the eleven
entries the search found, the index they were written to, and the eleven CIFs
under `cifs/cod`. That block is a record of one run rather than a constant,
because the database is added to continually.

### 24.3 simulate_pattern and SimulatedReflection

`simulate_pattern(cif_path, wavelength=KALPHA1_WAVELENGTH,
two_theta_range=(10, 100))` simulates the powder pattern of the structure in a
CIF and returns a list of `SimulatedReflection`.

The pattern is pymatgen's `XRDCalculator` on the first structure in the file,
in its conventional cell, at a single wavelength, so it carries no K alpha 2
lines at all. Partial occupancies are kept, so a disordered site scatters as
its average, which is what a tungsten bronze needs. It raises `ImportError`
when pymatgen is not installed, and `ValueError` when no structure can be read
from the file.

`SimulatedReflection` is a `NamedTuple` of `two_theta`, `intensity` and `hkl`.
The intensity is relative to the strongest line inside the range simulated,
taken as 100, so it depends on the window: narrowing the range can change
every intensity in it. The `hkl` is one representative of the family, with
four indices for a hexagonal cell, and where several families fall at one
angle the first of them is given.

That single wavelength is the reason a measured pattern has to be prepared
before it is compared. A simulated pattern has no satellites, so the observed
peaks must have theirs excluded; and a simulated pattern is at true angles, so
a measured zero offset must be taken off the observed positions. Section 24.6
does both.

### 24.4 match_candidate, CandidateMatch and ExplainedPeak

`match_candidate(observed_two_theta, simulated, tolerance=0.05,
min_intensity=10, exclude_two_theta=None)` weighs one simulated pattern
against a list of observed positions and returns a `CandidateMatch`.

An observed position is explained when a simulated reflection of any intensity
lies within `tolerance` degrees of it; where several do, the strongest is
taken as its cause, and a tie on intensity goes to the nearer. A simulated
reflection counts as missing when it is stronger than `min_intensity` and lies
more than `tolerance` from every observed position and from every position in
`exclude_two_theta`. That last argument is for the reflections of phases
already known to be present: a line hidden under one of those is not evidence
against the candidate.

Compare only over the range that was measured. A reflection simulated outside
the observed scan has nothing to be near and would be counted as missing,
which is why both the simulation and the peak search use one window.

`CandidateMatch` holds `explained`, a list of `ExplainedPeak`, and `missing`,
a list of `SimulatedReflection`, and its `score` is the count explained less
the count missing. `ExplainedPeak` is a `NamedTuple` of the `observed`
position and the `reflection` that accounts for it, with an `offset` property,
observed less simulated. It raises `ValueError` when `tolerance` is negative.

Section 6.3 reads that score properly and this section will not repeat it. The
short form: the missing count is the column that discriminates, a strong line
predicted where nothing was seen being hard evidence against a phase, while an
explained peak in a crowded pattern may be a coincidence.

### 24.5 The whole identification: rank_candidates and attribute_unexplained

Four more functions are the command's steps, each separate so that a caller
can stop after any of them.

`observed_peaks(scan, window=PHASES_WINDOW, zero=0.0)` is the peak positions
of a scan inside `window`, satellites excluded and `zero` subtracted, which is
the preparation Section 24.3 asks for in one call. It raises `ValueError` when
the scan carries no wavelength.

`rank_candidates(records, folder, observed, wavelength, window=PHASES_WINDOW,
tolerance=PHASES_TOLERANCE, simulate=None)` simulates every record's CIF from
`folder`, matches it, and returns a list of `Candidate` best first with ties
broken by COD id and `rank` numbering it from 1. A candidate whose strongest
simulated line is not observed is marked `rejected`, whatever it scores,
because the strongest line of a phase that is present is always there; that is
the rule Section 6.3 watches do its job on two of the eleven entries.
`Candidate` carries the `rank`, `cod_id`, `formula` and `space_group`, the
`explained`, `missing` and `score` counts, the `rejected` reason or an empty
string, the `cif` the pattern came from, and `explained_positions`.

`attribute_unexplained(observed, explained, cif, wavelength, phase,
window=PHASES_WINDOW, tolerance=PHASES_TOLERANCE, simulate=None)` takes the
peaks no candidate accounted for and tries the main phase's own pattern on
each. It returns a list of `UnexplainedPeak`, each with its `two_theta` and
`d_spacing` and, where the main phase has a reflection within tolerance, that
reflection's `phase`, `hkl`, `reflection_two_theta` and `intensity`; its
`identified` property says whether there was one. A peak with none is
unidentified and is the one worth chasing, which Section 6.4 discusses.

Both take a `simulate` callable that replaces `simulate_pattern`. That is the
seam for a caller who has patterns from elsewhere, or a test with no pymatgen,
and it is why the two functions can be exercised without the extra installed.

`require_phases_extra()` checks that pymatgen is importable and raises
`MissingPhasesExtra` if not. `MissingPhasesExtra` subclasses `ImportError`, so
a caller can tell a missing optional dependency from any other import failure,
and `rank_candidates` and `attribute_unexplained` convert the plain
`ImportError` that `simulate_pattern` raises into it. Call
`require_phases_extra` at the top of a script that will later simulate, so
that a missing extra costs the message and not a minute of downloads first.

Three constants set the defaults: `PHASES_WINDOW`, 10 to 80 degrees, the range
peaks are taken and patterns simulated over; `PHASES_TOLERANCE`, 0.15 degrees;
and `PHASES_PAUSE_S`, half a second. `COD_URL` is the database's address.

### 24.6 Matching one candidate against a measured pattern

The script below is the middle of `xrdkit phases` and nothing else: one CIF
already on disk, simulated and matched against the peaks of `pellet_a`. The
search, the downloads, the ranking of eleven candidates and the file of
unexplained peaks are all Section 6, and none of them is here, so the script
needs no network at all.

The CIF is `cifs/2100720.cif`, the tungsten bronze the project file of Section
2.2 names, which Section 6.3 explains is chosen for having the composition
nearest the sample's rather than for topping the ranking. The zero of 0.170
degrees is the offset Section 5.1 measured on this scan and Section 6.2 passed
to the command with `--zero`.

One thing about what you will see on screen. pymatgen writes its CIF parser
warnings to standard error, and several of the entries Section 6.2 fetched
raise them: 2311739 and 2311740 warn about stoichiometry, and 2103856 warns
that it found no symmetry operators and is defaulting to P1, which is the same
nonstandard centring Section 6.3 rejects. This particular CIF raises none, but
the printed block below is standard output only either way, so a warning would
never appear in it.

The whole of `phase_match.py`:

```python
from xrdkit.io import read_scan
from xrdkit.peaks import exclude_kalpha2, find_peaks
from xrdkit.phases import (
    PHASES_TOLERANCE,
    PHASES_WINDOW,
    match_candidate,
    simulate_pattern,
)

# Edit these lines for each new sample. Nothing below needs changing.
SCAN_FILE = "data/raw/pellet_a.xrdml"
CIF_FILE = "cifs/2100720.cif"
ZERO = 0.170
SHOWN = 6

scan = read_scan(SCAN_FILE)
peaks = exclude_kalpha2(find_peaks(scan, two_theta_range=PHASES_WINDOW))
observed = [peak.two_theta - ZERO for peak in peaks]
print(
    f"{len(observed)} peaks from {PHASES_WINDOW[0]:.0f} to {PHASES_WINDOW[1]:.0f}"
    f" degrees, zero {ZERO} degrees taken off"
)

simulated = simulate_pattern(
    CIF_FILE, wavelength=scan.wavelength, two_theta_range=PHASES_WINDOW
)
strongest = max(simulated, key=lambda reflection: reflection.intensity)
print(f"{len(simulated)} reflections simulated from {CIF_FILE}")
print("  two_theta  intensity  hkl")
for reflection in simulated[:SHOWN]:
    print(
        f"  {reflection.two_theta:9.3f}  {reflection.intensity:9.2f}  {reflection.hkl}"
    )
print(f"strongest {strongest.hkl} at {strongest.two_theta:.3f} degrees")

match = match_candidate(observed, simulated, tolerance=PHASES_TOLERANCE)
print(
    f"explained {len(match.explained)}/{len(observed)},"
    f" missing {len(match.missing)}, score {match.score}"
)
print(f"strongest line observed: {strongest not in match.missing}")
print("  observed  simulated  offset  hkl")
for peak in match.explained[:SHOWN]:
    print(
        f"  {peak.observed:8.3f}  {peak.reflection.two_theta:9.3f}"
        f"  {peak.offset:+6.3f}  {peak.reflection.hkl}"
    )
print("the strong reflections nothing was observed near")
for reflection in match.missing:
    print(
        f"  {reflection.two_theta:9.3f}  {reflection.intensity:9.2f}  {reflection.hkl}"
    )
```

It prints the peaks it prepared, the head of the simulated pattern, the match,
the first explained peaks and every strong reflection that was missed.

```text
35 peaks from 10 to 80 degrees, zero 0.17 degrees taken off
119 reflections simulated from cifs/2100720.cif
  two_theta  intensity  hkl
     10.012       0.02  (1, 1, 0)
     15.860       1.18  (2, 1, 0)
     20.101       0.00  (2, 2, 0)
     22.449      24.70  (0, 0, 1)
     22.503       8.17  (3, 1, 0)
     24.634       0.72  (1, 1, 1)
strongest (3, 1, 1) at 31.997 degrees
explained 27/35, missing 5, score 22
strongest line observed: True
  observed  simulated  offset  hkl
    22.578     22.449  +0.129  (0, 0, 1)
    25.712     25.708  +0.004  (3, 2, 0)
    26.751     26.650  +0.101  (2, 0, 1)
    27.725     27.607  +0.118  (2, 1, 1)
    28.593     28.577  +0.016  (4, 0, 0)
    29.489     29.476  +0.013  (4, 1, 0)
the strong reflections nothing was observed near
     45.824      38.77  (0, 0, 2)
     45.938      12.09  (6, 2, 0)
     55.429      32.40  (4, 1, 2)
     57.004      14.21  (4, 2, 2)
     71.569      17.11  (5, 5, 2)
```

The match is `explained 27/35, missing 5, score 22`, which is the 2100720 line
of the ranking in Section 6.3, figure for figure. It has to be: the command
calls these two functions with these arguments on this file.

The head of the simulated pattern is worth reading against Section 16.8, which
calculated the reflections of the same structure type from a cell alone. The
positions are the same reflections in the same order, and what is new is the
intensity column, which is what a CIF adds to a cell: (110) at 10.012 degrees
is allowed by the space group and has an intensity of 0.02, so it is a
reflection that exists and will never be seen. (220) at 20.101 is 0.00. A
simulated pattern is mostly lines like these, which is why `min_intensity`
exists: only the lines above 10 per cent count as missing when they are
absent.

The five missing reflections are the evidence against this particular entry,
and they are worth naming: (002) at 45.824 degrees carries 38.8 per cent of
the strongest line, (412) at 55.429 carries 32.4, and three more between 12
and 18 per cent. A phase whose second strongest line is unobserved is not
quite the phase in front of you. Section 6.3 makes the same point from the
other side, that 2100721 misses none: both entries are the same structure
type, and the missing count is measuring how far each one's cell sits from the
sample's, which is `xrdkit lattice`'s question and not this one's.

The offsets in the explained table say the same thing again. The first few
peaks sit 0.004 to 0.129 degrees above where this entry puts them, every one
in the same direction, which is a cell that is out rather than a scatter about
the right answer. The direction is the one to expect: the entry's cell is
12.4844 and 3.9572 angstrom and Section 7.1 refined this sample to 12.4740 and
3.9295, so the sample's spacings are the smaller and its reflections are the
higher in angle. The tolerance of 0.15 degrees is wide enough to absorb that,
which is exactly what it is for: this step identifies a structure type, and
the cell is refined afterwards.

## 25. library

`xrdkit.library` is the structure library shipped with the package: one TOML
file per entry, each describing a structure type apart from any one sample.
An entry is what a refinement needs to know before it has any coordinates at
all, and Section 2.2 is where a project file names one, as the `library` key
of a structure table.

There are eight public names: two dataclasses, two functions over the library,
one over an entry, and three constants.

### 25.1 What an entry is, and what it is not

An entry carries the site plan of a structure type. It names the space group
and crystal system, which cell parameters the system leaves free, the formula
units per cell, and one table per site giving that site's label, its kind, its
Wyckoff position, the coordinates that position leaves free, the Uiso group it
shares and the elements the prototype puts on it. It also names the anions
bonds are measured to and the range a bond from each kind of site may have.

It carries no coordinates. Not one number saying where an atom actually sits
is in the file, because those belong to a particular structure and the entry
describes a type. That is the whole reason the Rietveld modes need a CIF
beside a library entry, which is the limitation Section 10.3 records: the
entry says there is a site on 8d with x, y and z free, and only the CIF says
what x, y and z are.

It carries no reflection conditions either. Nothing in an entry lists what the
space group forbids; the symbol alone is stored, and `xrdkit.symmetry` derives
the operations, the absences and the multiplicities from it when they are
wanted, which Section 20 is the whole of. An entry and the symmetry module
therefore never disagree, because only one of them holds the fact.

### 25.2 list_entries, load_entry and the file layout

`list_entries(root=None)` returns the names of the entries in the library,
sorted. A name is `"<family>/<name>"`, from the file at
`xrdkit/structures/<family>/<name>.toml`, so `ttb/P4bm` is `P4bm.toml` in the
`ttb` folder. `root` is a folder to look in instead of the library shipped
with xrdkit, laid out the same way, which is how a project keeps entries of
its own.

`load_entry(name, root=None)` reads one entry with `tomllib`, checks it and
returns a `StructureEntry`. It raises `ValueError` when there is no entry of
that name, naming the ones there are, and when the file is not valid TOML or
does not check out, naming the entry and the field.

A file holds an `[entry]` table and one `[[sites]]` table per site. `[entry]`
requires `name`, which must match the file's path, `family`,
`crystal_system`, `space_group`, `cell_parameters`, `z` and `reference`, and
allows `setting`, `polar_axis`, `origin_site`, `anions` and `bond_limits`.
Each `[[sites]]` requires `label`, unique within the entry, `kind`, and
`wyckoff`, and allows `free`, `uiso_group` and `elements`.

### 25.3 StructureEntry, Site and the constants

`StructureEntry` is frozen and holds exactly those fields.

| Field | What it holds |
| --- | --- |
| `name`, `family`, `reference` | the entry's name, the structure family as text, and the structure it was taken from |
| `crystal_system`, `space_group`, `setting` | one of `CRYSTAL_SYSTEMS`, the Hermann-Mauguin symbol, and the setting or origin choice, empty for the standard one |
| `cell_parameters` | the parameters the crystal system leaves free, exactly as `CELL_PARAMETERS` lists them |
| `z` | formula units per cell |
| `sites` | the `Site` tables in file order |
| `polar_axis`, `origin_site` | the axis along which symmetry leaves the origin free, and the label of the site whose coordinate along it is held |
| `anions` | the elements bonds are measured to, `DEFAULT_ANIONS` when the file names none |
| `bond_limits` | by kind of site, `(min, max)` in angstroms, only for the kinds the file lists |

Its `kinds` property is the kinds of site the entry declares, in the order of
its sites, each once.

`Site` is frozen too: `label`, `kind`, `wyckoff`, `free`, the coordinates the
Wyckoff position leaves free drawn from x, y and z, `uiso_group`, the name of
the group whose one Uiso the site shares or `None`, and `elements`, the
elements the prototype puts on it.

`CELL_PARAMETERS` is the dict Section 17.2 documents, crystal system to the
names of its free parameters, and `CRYSTAL_SYSTEMS` is its keys. They live
here rather than in `xrdkit.cell` because an entry has to be checked against
them before any cell exists.

`bond_limits(entry, kind)` returns the range in angstroms of a bond from a
site of that kind to an anion: the entry's own where it gives one, else
`DEFAULT_BOND_LIMITS`, 1.6 to 3.0 angstrom. It raises `ValueError` when
`kind` is not a kind of the entry's sites. `DEFAULT_ANIONS` is `("O",)`.

### 25.4 The entry the examples use

The script prints the library, then the entry of the tungsten bronze, then
what its validator refuses.

Start a new file named `library_entry.py`.

```python
from xrdkit.library import CELL_PARAMETERS, bond_limits, list_entries, load_entry

# Edit these lines for each new entry. Nothing below needs changing.
ENTRY = "ttb/P4bm"

print(f"{len(list_entries())} entries in the library")
for name in list_entries():
    print(f"  {name}")

entry = load_entry(ENTRY)
print(f"name           {entry.name}")
print(f"family         {entry.family}")
print(f"crystal system {entry.crystal_system}, space group {entry.space_group}")
print(f"setting        {entry.setting!r}")
print(
    f"cell parameters {entry.cell_parameters}, "
    f"the {entry.crystal_system} pair of CELL_PARAMETERS: "
    f"{CELL_PARAMETERS[entry.crystal_system]}"
)
print(f"z              {entry.z} formula units per cell")
print(f"reference      {entry.reference}")
print(f"polar axis     {entry.polar_axis}, origin held on {entry.origin_site}")
print(f"anions         {entry.anions}")
print(f"kinds          {entry.kinds}")
for kind in entry.kinds:
    low, high = bond_limits(entry, kind)
    own = "the entry's own" if kind in entry.bond_limits else "the package default"
    print(f"  {kind} to an anion  {low} to {high} angstrom, {own}")

print(f"{len(entry.sites)} sites")
print("  label  kind  wyckoff  free   uiso_group  elements")
for site in entry.sites:
    free = "".join(site.free) or "-"
    print(
        f"  {site.label:5s}  {site.kind:4s}  {site.wyckoff:7s}  {free:5s}"
        f"  {site.uiso_group or '-':10s}  {', '.join(site.elements)}"
    )
```

It prints the entries, the entry, its bond limits and its sites.

```text
7 entries in the library
  perovskite/Amm2
  perovskite/P4mm
  perovskite/Pbnm
  perovskite/Pm-3m
  perovskite/R3c
  ttb/P4bm
  ttb/P4mbm
name           ttb/P4bm
family         tetragonal tungsten bronze
crystal system tetragonal, space group P4bm
setting        ''
cell parameters ('a', 'c'), the tetragonal pair of CELL_PARAMETERS: ('a', 'c')
z              5 formula units per cell
reference      COD 2100720
polar axis     c, origin held on B1
anions         ('O',)
kinds          ('A', 'B', 'O')
  A to an anion  2.45 to 3.0 angstrom, the entry's own
  B to an anion  1.8 to 2.2 angstrom, the entry's own
  O to an anion  1.6 to 3.0 angstrom, the package default
9 sites
  label  kind  wyckoff  free   uiso_group  elements
  A1     A     2a       z      A           Sr
  A2     A     4c       xz     A           Ba, Sr
  B1     B     2b       z      B           Nb
  B2     B     8d       xyz    B           Nb
  O1     O     4c       xz     O           O
  O2     O     8d       xyz    O           O
  O3     O     8d       xyz    O           O
  O4     O     2b       z      O           O
  O5     O     8d       xyz    O           O
```

Seven entries in two families, and `ttb/P4bm` is the one the
`[structures.ttb_p4bm]` table of Section 2.2 names. Read its nine sites as the
plan the later sections work to. Two A sites, two B sites and five O sites, on
Wyckoff positions 2a, 4c, 2b and 8d; a cell holds two A1, four A2, two B1,
eight B2 and thirty O, which at Z of 5 and a formula of AB2O6 is five A
cations in six A positions, ten B and thirty O. One A position in six is
empty, which is what makes this a bronze rather than a perovskite.

The `free` column is the Wyckoff position speaking. A1 on 2a is at (0, 0, z),
so only z is free; A2 on 4c is at (x, x + 1/2, z), so x and z are free and y
follows x; B2 on 8d is the general position and has all three. Those are the
coordinates the Rietveld coordinates mode of Section 8 is allowed to refine,
and they are read off the entry rather than guessed from the CIF.

The Uiso groups are the other economy. Nine sites share three displacement
parameters, one for the A sites, one for the B and one for the O, which is
what makes a refinement of this structure possible on a laboratory pattern at
all. The elements column is the prototype's, strontium on A1 and barium and
strontium on A2; the lanthanum and titanium of the sample's composition are
not there, which is exactly what the `atoms` table of Section 2.2 is for and
what Section 26.4 does.

The bond limits are the entry's judgement about its own structure type. A to
an anion runs 2.45 to 3.0 angstrom and B to an anion 1.8 to 2.2, both narrower
than the package default of 1.6 to 3.0 that the O kind falls back on, since
the entry lists no limits for a site that is itself an anion.

Continues `library_entry.py`. Add these lines at the end of the file.

```python
import tempfile
from pathlib import Path

GOOD = """
[entry]
name = "demo/One"
family = "demonstration"
crystal_system = "tetragonal"
space_group = "P4bm"
cell_parameters = ["a", "c"]
z = 1
reference = "written for this section"

[[sites]]
label = "A1"
kind = "A"
wyckoff = "2a"
"""
BROKEN = {
    "a cell parameter the system does not leave free": (
        'cell_parameters = ["a", "c"]',
        'cell_parameters = ["a", "b", "c"]',
    ),
    "an origin site that is not a site": ("z = 1", 'z = 1\norigin_site = "A9"'),
}

with tempfile.TemporaryDirectory() as folder:
    root = Path(folder) / "demo"
    root.mkdir()
    (root / "One.toml").write_text(GOOD, encoding="utf-8")
    print(f"a library of my own holds {list_entries(root=root.parent)}")
    print(f"and it loads: {load_entry('demo/One', root=root.parent).name}")
    for description, (before, after) in BROKEN.items():
        (root / "One.toml").write_text(GOOD.replace(before, after), encoding="utf-8")
        try:
            load_entry("demo/One", root=root.parent)
        except ValueError as error:
            print(f"{description}:")
            print(f"  {error}")

try:
    load_entry("ttb/P4bmm")
except ValueError as error:
    print(f"an entry that is not there:\n  {error}")
try:
    bond_limits(entry, "C")
except ValueError as error:
    print(f"a kind the entry has no site of:\n  {error}")
```

It builds a small library of its own in a temporary folder, breaks its one
entry twice, and then asks the shipped library two questions it cannot answer.

```text
a library of my own holds ['demo/One']
and it loads: demo/One
a cell parameter the system does not leave free:
  structure entry 'demo/One': entry.cell_parameters: a tetragonal cell has ['a', 'c'], not ['a', 'b', 'c']
an origin site that is not a site:
  structure entry 'demo/One': entry.origin_site: no site labelled 'A9'; the sites are A1
an entry that is not there:
  no structure entry 'ttb/P4bmm'; the entries are perovskite/Amm2, perovskite/P4mm, perovskite/Pbnm, perovskite/Pm-3m, perovskite/R3c, ttb/P4bm, ttb/P4mbm
a kind the entry has no site of:
  'C' is not a kind of site of ttb/P4bm; its kinds are A, B, O
```

The `root` argument is what makes this possible and is worth knowing about for
its own sake: an entry for a structure type the shipped library lacks goes in
a folder of your own laid out the same way, and `list_entries` and
`load_entry` read it with no change to the package.

Both refusals are cross checks rather than type checks, which is the character
of this validator. The first catches `cell_parameters` that do not match the
crystal system in the same table, and says what a tetragonal cell should have.
The second catches an `origin_site` naming a site the file does not define,
and lists the sites it does. A file is checked as a whole when it loads, so a
refinement never gets as far as discovering that its origin site does not
exist.

The last two are the ordinary errors of the two functions, and both name what
was available: `load_entry` lists every entry there is, and `bond_limits`
lists the kinds the entry declares.

## 26. structure

`xrdkit.structure` is the geometry and the bookkeeping a Rietveld refinement
needs before it starts: which atoms sit on which site, what the distances
between them are, and what occupancies put a nominal composition on those
sites. Section 8.11 is how the pipeline uses all of it from a project file,
and this section is the four functions themselves.

There are seven public names. `metric_tensor`, `interatomic_distances`,
`Distance` and `bond_lengths` are the geometry; `cell_contents`,
`composition_edits` and `site_setup` are the bookkeeping.

Every function here takes a structure in the form the GSAS-II driver reports
one: a `cell` of six numbers or GSAS-II's own keys, `atoms` as mappings with
`label`, `type`, `xyz` and optionally `xyz_esd`, `multiplicity` and
`occupancy`, and `operators` as mappings with a 3 by 3 `rotation` and a
`translation`. Nothing in the package turns a CIF into that form without
GSAS-II, so a script that wants to work on a CIF alone reads its atom loop
itself, which is what Section 26.1 does and says.

### 26.1 metric_tensor, interatomic_distances and bond_lengths

`metric_tensor(cell)` is the real space metric tensor G of a cell given as
`(a, b, c, alpha, beta, gamma)` in angstroms and degrees, so a fractional
vector v has length the square root of v dotted into G v. It is the same
quantity `Cell.metric_tensor` of Section 17.3 gives, reached from six numbers
rather than from a `Cell`, and it is what `Cell` itself calls.

`interatomic_distances(cell, atoms, operators, centres, targets, dmax,
dmin=0.5)` returns every distance from each atom in `centres` to images of the
atoms in `targets` between `dmin` and `dmax` angstroms, shortest first for
each centre. The images are those of every operator, shifted into the
twenty-seven cells around the centre; the same image reached by two operators,
as happens on a special position, is counted once.

Each comes back as a `Distance`: the `centre` and `target` labels, the
`distance`, its `esd`, the index of the `operator` that made the image and the
whole cell `translation` added to it. The esd is propagated from the esds of
the fractional coordinates of both atoms alone, treated as uncorrelated, with
the cell taken as exact, and is `None` when neither atom has any. It raises
`ValueError` when a centre or target is not an atom's label, or `dmax` is not
above `dmin`.

`bond_lengths(phase, anions=None, dmax=3.0, dmin=0.5)` is that function
wrapped for a whole phase. Atoms that share a position count as one site,
taken as the first of them; a site is an anion site when its first atom is of
an element in `anions` and a cation site otherwise, and every cation site to
anion site distance is found. Distances from one cation site to images of one
anion site that agree to within a ten thousandth of an angstrom are the same
bond and are counted rather than listed twice.

It returns a list of dicts with `centre` and `centre_site`, `target` and
`target_site`, the first atom's label and every label on that site joined by a
solidus, the `distance`, its `esd` and the `count`, ordered by cation site in
the phase's order and then by distance.

`anions` is required in practice: left out it warns that it is deprecated and
falls back to `DEFAULT_ANIONS`, for callers written before it was required.
Pass the anions of the structure's library entry, which `site_setup` reports.
A phase with no anion site at all gets an empty list and a warning saying so.

### 26.2 Reading a CIF and finding its sites

`site_setup(structure, atoms)` matches the sites a structure table names
against the atoms a phase actually has, and returns the plan everything
downstream works from. `structure` is a structure table as
`xrdkit.config.load_config` returns it: a `name`, a `library` entry name or
`None`, a list of `sites` each with a `name`, an optional entry `label`, the
`atoms` on it as label to element, a `wyckoff` and a `kind`, plus
`free_coordinates`, `uiso_groups`, `origin` and `exchange`.

The returned plan is a dict of `sites`, each with its name, kind, Wyckoff
position, multiplicity, free coordinates and the atoms on it; `kinds`, the
sites of each kind; `uiso_groups` and `group_names`; `origin` and
`origin_axis`; `coordinates`, by kind, the coordinates to refine on every site
but the origin's; `exchange`; `polar_axis`; `anions`; and `bond_limits` for
every kind, the structure's own over the entry's over the package default.

It raises `ValueError` when a named atom is not among the atoms, when a site's
atoms are not all on one position or two sites land on the same one, when an
atom is not of the element given for it, when a site's multiplicity is not
that of its Wyckoff position, or when an atom is on no configured site. That
last one matters: a CIF with an atom nobody accounted for is refused rather
than silently ignored.

The script below does the whole of Section 26 on the example CIF and the entry
of Section 25. It begins by reading the CIF, because it has to: the functions
take atoms in the driver's form and the driver needs GSAS-II, so a script
working from a CIF alone parses the atom loop itself. The multiplicity of each
atom is not read from the CIF at all but taken from the Wyckoff position the
library entry gives for its site, which is the entry and the CIF doing exactly
the jobs Section 25.1 divides between them.

Start a new file named `site_plan.py`.

```python
import re
from pathlib import Path

from xrdkit.config import wyckoff_multiplicity
from xrdkit.library import load_entry

# Edit these lines for each new structure. Nothing below needs changing.
CIF_FILE = "cifs/2100720.cif"
ENTRY = "ttb/P4bm"
# The atoms table of Section 2.2: which CIF atom sits on which site of the
# entry, and which element each is.
ATOMS = {
    "A1": {"Sr1": "Sr", "La1": "La"},
    "A2": {"Ba2": "Ba", "Sr2": "Sr"},
    "B1": {"Nb1": "Nb", "Ti1": "Ti"},
    "B2": {"Nb2": "Nb", "Ti2": "Ti"},
}

ESD = re.compile(r"\(\d+\)$")


def number(text):
    """A CIF number without its esd in brackets."""
    return float(ESD.sub("", text))


def read_cell(path):
    """The six cell parameters a CIF gives, esds stripped."""
    text = Path(path).read_text(encoding="utf-8")
    keys = [f"_cell_length_{axis}" for axis in ("a", "b", "c")]
    keys += [f"_cell_angle_{angle}" for angle in ("alpha", "beta", "gamma")]
    found = []
    for key in keys:
        match = re.search(rf"^{key}\s+(\S+)", text, re.MULTILINE)
        found.append(number(match.group(1)))
    return tuple(found)


def read_atom_loop(path):
    """The _atom_site loop of a CIF, one dict of tag to text per atom."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    start = next(
        index
        for index, line in enumerate(lines)
        if line.strip() == "loop_"
        and lines[index + 1].strip() == "_atom_site_type_symbol"
    )
    tags, row = [], start + 1
    while lines[row].strip().startswith("_"):
        tags.append(lines[row].strip())
        row += 1
    rows = []
    while row < len(lines) and lines[row].strip() and lines[row].strip() != "loop_":
        rows.append(dict(zip(tags, lines[row].split())))
        row += 1
    return rows


entry = load_entry(ENTRY)
CELL = read_cell(CIF_FILE)
wyckoff = {site.label: site.wyckoff for site in entry.sites}
site_of = {label: name for name, placed in ATOMS.items() for label in placed}
atoms = []
for row in read_atom_loop(CIF_FILE):
    label = row["_atom_site_label"]
    name = site_of.get(label, label)
    atoms.append(
        {
            "label": label,
            "type": row["_atom_site_type_symbol"],
            "xyz": [
                number(row[f"_atom_site_fract_{axis}"]) for axis in ("x", "y", "z")
            ],
            "multiplicity": wyckoff_multiplicity(wyckoff[name]),
            "occupancy": number(row["_atom_site_occupancy"]),
        }
    )
print(f"{len(atoms)} atoms read from {CIF_FILE}")
print(f"its cell is {CELL}")
print("  label  type  site  wyckoff  mult  occupancy")
for atom in atoms:
    name = site_of.get(atom["label"], atom["label"])
    print(
        f"  {atom['label']:5s}  {atom['type']:4s}  {name:4s}  {wyckoff[name]:7s}"
        f"  {atom['multiplicity']:4d}  {atom['occupancy']:9.4f}"
    )
```

It prints the cell and the atoms, each against the site it belongs to.

```text
10 atoms read from cifs/2100720.cif
its cell is (12.4844, 12.4844, 3.9572, 90.0, 90.0, 90.0)
  label  type  site  wyckoff  mult  occupancy
  Ba2    Ba    A2    4c          4     0.6500
  Sr2    Sr    A2    4c          4     0.2462
  Sr1    Sr    A1    2a          2     0.7080
  Nb1    Nb    B1    2b          2     1.0000
  Nb2    Nb    B2    8d          8     1.0000
  O1     O     O1    4c          4     1.0000
  O2     O     O2    8d          8     1.0000
  O3     O     O3    8d          8     1.0000
  O4     O     O4    2b          2     1.0000
  O5     O     O5    8d          8     1.0000
```

Ten atoms on nine sites: `Ba2` and `Sr2` share the 4c A2 position, which is
the disorder the structure type is built on, and every other site has one atom
in this CIF. The `ATOMS` table in the settings is the `atoms` table of Section
2.2 read as a dictionary, and it names `La1`, `Ti1` and `Ti2` as well, which
this CIF does not have; those are the composition's business and Section 26.4
puts them in.

Continues `site_plan.py`. Add these lines at the end of the file.

```python
from xrdkit.structure import site_setup

present = {atom["label"] for atom in atoms}
sites = []
for site in entry.sites:
    placed = ATOMS.get(site.label)
    on_site = (
        {label: element for label, element in placed.items() if label in present}
        if placed
        else {site.label: site.elements[0]}
    )
    sites.append(
        {
            "name": next(iter(on_site)),
            "label": site.label,
            "atoms": on_site,
            "wyckoff": site.wyckoff,
            "kind": site.kind,
        }
    )
groups = {}
for site in entry.sites:
    groups.setdefault(site.uiso_group or site.label, []).append(
        next(s["name"] for s in sites if s["label"] == site.label)
    )
structure = {
    "name": "ttb_p4bm",
    "library": entry.name,
    "sites": sites,
    "free_coordinates": {},
    "uiso_groups": [
        {"name": name, "sites": members} for name, members in groups.items()
    ],
    "origin": {"site": "Nb1", "axis": "z"},
    "exchange": {"elements": ["Sr", "Ba"], "sites": ["Sr1", "Ba2"]},
}

plan = site_setup(structure, atoms)
print(f"{len(plan['sites'])} sites planned")
print("  name  kind  wyckoff  mult  free  atoms")
for site in plan["sites"]:
    labels = ", ".join(atom["label"] for atom in site["atoms"])
    print(
        f"  {site['name']:4s}  {site['kind']:4s}  {site['wyckoff']:7s}"
        f"  {site['multiplicity']:4d}  {site['free'] or '-':4s}  {labels}"
    )
print(f"kinds          { {k: len(v) for k, v in plan['kinds'].items()} }")
print(f"uiso groups    {plan['group_names']} over {plan['uiso_groups']}")
print(
    f"origin         {plan['origin']['name']} along {plan['origin_axis']},"
    f" polar axis {plan['polar_axis']}"
)
print(f"anions         {plan['anions']}")
print(f"bond limits    {plan['bond_limits']}")
print(f"exchange       {plan['exchange']}")
print(f"coordinates    {plan['coordinates']}")
```

It prints the plan.

```text
9 sites planned
  name  kind  wyckoff  mult  free  atoms
  Sr1   A     2a          2  z     Sr1
  Ba2   A     4c          4  xz    Ba2, Sr2
  Nb1   B     2b          2  z     Nb1
  Nb2   B     8d          8  xyz   Nb2
  O1    O     4c          4  xz    O1
  O2    O     8d          8  xyz   O2
  O3    O     8d          8  xyz   O3
  O4    O     2b          2  z     O4
  O5    O     8d          8  xyz   O5
kinds          {'A': 2, 'B': 2, 'O': 5}
uiso groups    ['A', 'B', 'O'] over [['Sr1', 'Ba2'], ['Nb1', 'Nb2'], ['O1', 'O2', 'O3', 'O4', 'O5']]
origin         Nb1 along z, polar axis z
anions         ['O']
bond limits    {'A': (2.45, 3.0), 'B': (1.8, 2.2), 'O': (1.6, 3.0)}
exchange       {'elements': ['Sr', 'Ba'], 'sites': ['Sr1', 'Ba2']}
coordinates    {'A': {'Sr1': 'z', 'Ba2': 'xz'}, 'B': {'Nb2': 'xyz'}, 'O': {'O1': 'xz', 'O2': 'xyz', 'O3': 'xyz', 'O4': 'z', 'O5': 'xyz'}}
```

Each site is named by the first CIF atom on it rather than by the entry's
label, which is the convention the pipeline uses: `A2` becomes `Ba2` because
that is the atom a refinement will address. The free coordinates come from the
entry, the multiplicities from the Wyckoff positions, the Uiso groups gather
the nine sites into three, and the bond limits are the entry's for A and B and
the package default for O, exactly as Section 25.4 printed them.

`coordinates` is the one to look at twice. It lists eight sites, not nine:
`Nb1` is missing from the B kind because it is the origin site, and its z is
held to stop the whole structure sliding along the polar axis. That is the
rule Section 8.11 states, and it is applied here rather than remembered later.

Continues `site_plan.py`. Add these lines at the end of the file.

```python
from xrdkit.structure import bond_lengths, metric_tensor
from xrdkit.symmetry import space_group_operations

metric = metric_tensor(CELL)
print("metric tensor G, from the cell alone")
for row in metric:
    print("  " + "  ".join(f"{round(value, 6) + 0.0:11.6f}" for value in row))

operators = [
    {
        "rotation": [list(row) for row in rotation],
        "translation": [float(shift) for shift in translation],
    }
    for rotation, translation in space_group_operations(entry.space_group)
]
phase = {"cell": list(CELL), "atoms": atoms, "operators": operators}
bonds = bond_lengths(phase, anions=plan["anions"], dmax=3.0)
print(f"{len(bonds)} distinct cation to anion distances under 3.0 angstrom")
print("  site     kind  to     distance  count  within the entry's limits")
by_name = {site["name"]: site for site in plan["sites"]}
for bond in bonds:
    site = by_name[bond["centre"]]
    low, high = plan["bond_limits"][site["kind"]]
    inside = "yes" if low <= bond["distance"] <= high else "no"
    print(
        f"  {bond['centre_site']:8s} {site['kind']:4s}  {bond['target_site']:5s}"
        f"  {bond['distance']:8.4f}  {bond['count']:5d}  {inside}"
    )
```

It prints the metric tensor and every cation to anion distance under three
angstroms.

```text
metric tensor G, from the cell alone
   155.860243     0.000000     0.000000
     0.000000   155.860243     0.000000
     0.000000     0.000000    15.659432
16 distinct cation to anion distances under 3.0 angstrom
  site     kind  to     distance  count  within the entry's limits
  Ba2/Sr2  A     O1       2.7157      1  yes
  Ba2/Sr2  A     O3       2.7754      2  yes
  Ba2/Sr2  A     O1       2.8575      1  yes
  Ba2/Sr2  A     O3       2.9665      2  yes
  Sr1      A     O2       2.6853      4  yes
  Sr1      A     O5       2.7249      4  yes
  Sr1      A     O2       2.8587      4  yes
  Nb1      B     O4       1.8321      1  yes
  Nb1      B     O3       1.9661      4  yes
  Nb1      B     O4       2.1251      1  yes
  Nb2      B     O5       1.8409      1  yes
  Nb2      B     O3       1.9424      1  yes
  Nb2      B     O2       1.9573      1  yes
  Nb2      B     O1       1.9971      1  yes
  Nb2      B     O2       2.0130      1  yes
  Nb2      B     O5       2.1199      1  yes
```

The symmetry operators come from `space_group_operations` of Section 20,
converted into the mappings this module takes, so no GSAS-II is needed to
generate the images either. The metric tensor is diagonal because the cell is
tetragonal, its first two entries a squared and its third c squared, which
Section 17.3 says is the whole of the tetragonal geometry.

Sixteen distinct distances, and every one of them falls inside the limits its
kind's entry gives. That is the check this function exists for. The two B
sites are octahedra of niobium and oxygen, six bonds each, 1.83 to 2.13
angstrom and comfortably inside the entry's 1.8 to 2.2; the A sites run 2.69
to 2.97, inside 2.45 to 3.0. A refined structure whose bonds have wandered
outside those ranges has gone wrong somewhere, and this is how a script finds
out.

The counts are the multiplicity of each bond about its centre. `Sr1` has four
of each of its distances because it sits on 2a, where the four fold axis
repeats every neighbour four times; `Nb2` on the general position has one of
each, since nothing repeats a bond from a site with no symmetry of its own.

### 26.3 cell_contents

`cell_contents(atoms)` is atoms of each element per cell, multiplicity times
occupancy summed over the atoms, from atoms carrying `type`, `multiplicity`
and `occupancy`. It is two lines of arithmetic and it is the thing to print
whenever a composition is in doubt, because it is what a refinement is
actually working with.

### 26.4 composition_edits, and how an added element is shared

`composition_edits(atoms, composition, structure)` returns the atom edits that
give a cell of `atoms` the nominal `composition`, in atoms of each element per
formula unit, by the rule the `structure` table carries. The cell holds
`structure["formula_units"]` formula units.

Two rules are at work, and they are different.

Every element the atoms already hold is scaled by one factor on every site it
occupies. That keeps the CIF's distribution of it over the sites, whatever
that distribution was, and an atom whose occupancy the scaling leaves
unchanged gets no edit at all.

Every element in `structure["composition"]["added"]` is placed against a host.
The rule maps the added element either to a host element, when it goes on
every atom of that element and is labelled by the added element and the rest
of the host's label, `La1` on `Sr1`; or to a mapping of host label to added
label, when it goes beside those atoms only, labelled as given. Either way its
whole content goes onto its host atoms, each at one fraction of that atom's
own occupancy, so it is shared among them in proportion to multiplicity times
the host's occupancy. An added element whose content is zero is not added.

That proportional sharing is the rule Section 8.11 states and Section 10.5
records as a limitation: the amount on each site follows its host, and there
is no way to say how much goes on each site. Where that is not what the
structure does, the answer is a refinement of the occupancies afterwards, not
a different starting split.

It raises `ValueError` when the structure gives no formula units, when the
atoms hold an element the composition lacks, when the composition has one the
atoms lack that the rule does not add, when an added element is on the atoms
already, when its host is not among them, or when a label for an added atom is
taken.

The edits come back in the form the GSAS-II driver takes: `{"label",
"occupancy"}` for a change, and `{"label", "type", "copy", "occupancy"}` for
an atom to add, `copy` naming the atom whose position it takes.

Continues `site_plan.py`. Add these lines at the end of the file.

```python
from xrdkit.density import parse_formula
from xrdkit.structure import cell_contents, composition_edits

FORMULA = "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
# Section 2.2's atoms table read as a placement rule: each added element
# against the host atoms it goes beside, and the label it takes on each.
ADDED = {
    "La": {"Sr1": "La1"},
    "Ti": {"Nb1": "Ti1", "Nb2": "Ti2"},
}

composition = parse_formula(FORMULA)
print(f"nominal {composition}")
print(f"the CIF cell holds {cell_contents(atoms)}")
edits = composition_edits(
    atoms,
    composition,
    {"name": "ttb_p4bm", "formula_units": entry.z, "composition": {"added": ADDED}},
)
print(f"{len(edits)} edits for Z = {entry.z}")
print("  label  type  copy  occupancy  was")
was = {atom["label"]: atom["occupancy"] for atom in atoms}
for edit in edits:
    print(
        f"  {edit['label']:5s}  {edit.get('type', '-'):4s}  {edit.get('copy', '-'):4s}"
        f"  {edit['occupancy']:9.4f}  {was.get(edit['label'], 0.0):9.4f}"
    )
edited = []
for atom in atoms:
    change = next((e for e in edits if e["label"] == atom["label"]), None)
    edited.append(atom if change is None else {**atom, **change})
for edit in edits:
    if "type" in edit:
        host = next(atom for atom in atoms if atom["label"] == edit["copy"])
        edited.append({**host, **edit})
print(f"the edited cell holds {cell_contents(edited)}")
for element in ("La", "Ti"):
    share = {
        edit["label"]: next(
            atom["multiplicity"] for atom in atoms if atom["label"] == edit["copy"]
        )
        * edit["occupancy"]
        for edit in edits
        if edit.get("type") == element
    }
    print(f"{element} per cell by site: {share}, {sum(share.values()):.1f} in all")
```

It prints the nominal composition, what the CIF holds, the edits, and what the
cell holds once they are applied.

```text
nominal {'Sr': 0.4, 'Ba': 0.5, 'La': 0.1, 'Nb': 1.9, 'Ti': 0.1, 'O': 6.0}
the CIF cell holds {'Ba': 2.6, 'Sr': 2.4008, 'Nb': 10.0, 'O': 30.0}
8 edits for Z = 5
  label  type  copy  occupancy  was
  Ba2    -     -        0.6250     0.6500
  Sr2    -     -        0.2051     0.2462
  Sr1    -     -        0.5898     0.7080
  La1    La    Sr1      0.2500     0.0000
  Nb1    -     -        0.9500     1.0000
  Ti1    Ti    Nb1      0.0500     0.0000
  Nb2    -     -        0.9500     1.0000
  Ti2    Ti    Nb2      0.0500     0.0000
the edited cell holds {'Ba': 2.5, 'Sr': 2.0, 'Nb': 9.5, 'O': 30.0, 'La': 0.5, 'Ti': 0.5}
La per cell by site: {'La1': 0.5}, 0.5 in all
Ti per cell by site: {'Ti1': 0.1, 'Ti2': 0.4}, 0.5 in all
```

Read the edits in two groups. The three scalings first. The CIF holds 2.4008
strontium per cell and the composition wants 2.0, so every strontium atom is
multiplied by 0.8331: `Sr1` falls from 0.7080 to 0.5898 and `Sr2` from 0.2462
to 0.2051, and the ratio between the two sites is untouched. Barium goes from
2.6 to 2.5 the same way, and niobium from 10 to 9.5. Oxygen is already right
and gets no edit.

Then the additions. Lanthanum is named against `Sr1` alone, so the whole half
an atom per cell goes there: `La1` copies `Sr1`'s position with an occupancy
of 0.25, which on a site of multiplicity 2 is 0.5 atoms. Titanium is named
against both niobium atoms, and the last two lines of the block show what the
sharing rule does with it. Both `Ti1` and `Ti2` get an occupancy of 0.05,
which looks like an even split and is not one: `Nb1` has multiplicity 2 and
`Nb2` has 8, so the titanium lands 0.1 on B1 and 0.4 on B2, four times as much
on the site with four times the room. Had the `atoms` table named B2 alone,
all 0.5 would have gone there.

The last line is the check worth making in any script that does this. The
edited cell holds 2.0 strontium, 2.5 barium, 0.5 lanthanum, 9.5 niobium, 0.5
titanium and 30 oxygen, which is five times the nominal formula unit exactly,
as it must be at Z of 5. If that line does not come out right the rule has
been given something it cannot do, and the place to find out is here rather
than at the end of a refinement.

## 27. project

`xrdkit.project` reads `xrdkit.toml`. Section 2.2 is the file itself, table by
table, and what every key means; this section is the loader, the five frozen
dataclasses it returns and the four resolvers that answer the questions a
command asks of a project.

There are twelve public names: two constants, five dataclasses, three
functions that read a project and four that resolve something from one.

The reason to use them rather than `tomllib` is that the file is checked as a
whole. A sample naming an instrument that is not there, a structure whose cell
is not its entry's, a path that does not exist: all of them are caught when
the file loads, so a command fails before it starts rather than in the middle
of a refinement.

### 27.1 find_project, load_project and load_project_text

`PROJECT_FILE` is `"xrdkit.toml"`, and the project root is the folder holding
it. Every path in the file is relative to that root, and every path the loader
returns has been made absolute against it.

`find_project(start=None)` returns the `xrdkit.toml` of the first folder, from
`start` or the current folder upwards, that holds one. It raises
`FileNotFoundError` naming the folder it began at when neither that folder nor
any above it has one. That upward walk is why a command works from anywhere
inside a project.

`load_project(path=None)` reads and checks the file at `path`, an
`xrdkit.toml` or the folder holding one, and by default the one
`find_project` finds. It returns a `Project`.

`load_project_text(text, path)` is the same check on text you already have,
with `path` standing in for where it came from, so relative paths resolve
against that file's folder and the messages name it. It is what a test uses,
and what Section 27.4 uses to show a refusal without writing anything.

Both raise `ValueError` naming the file and the exact key at fault.

### 27.2 The five dataclasses

All five are frozen.

`Project` is the whole file: `root`, the folder holding it; `name` and
`version`; `instruments`, `structures` and `samples`, each a dict by key; and
`refine`, the project's `[refine]` defaults as a `Refine`.

`Instrument` carries `key`, `wavelength` as a tuple of one or two floats in
angstroms, `ka2`, and optionally `radius` in millimetres and `instprm`, the
path of the GSAS-II instrument parameter file.

`StructureSpec` carries `key`, `composition` as the formula text,
`library` and `cif`, `cell` as the parameters the file gave, `z`, `exchange`
as a tuple of tuples, and the origin as three fields: `origin`, the site
label or `None`; `origin_axis`; and `origin_fixed`, false only when the file
said `origin = false`. `atoms` is `None` when the table gives none, else its
sites as `xrdkit.config.read_sites` gives them, each carrying its entry
`label` for a library structure.

`Sample` carries `key`, `file`, `instrument`, `structures` as a tuple of keys,
`form`, which is one of `FORMS`, powder or pellet, and the optional `stage`,
`temperature_c`, `archimedes` and `notes`, with `refine`, the keys of
`[refine]` this sample overrides, checked but not yet merged.

`Refine` carries `two_theta`, `None` for the scan's own range; `background`;
`max_passes` and `unsettled`, both by mode; and `followed`, the reflections
carried from mode to mode as (h, k, l) triples.

### 27.3 The four resolvers

Each answers a question the file does not answer literally.

`resolved_z(spec)` is the formula units per cell: the structure's own `z`, or
its library entry's when it gives none. It raises `ValueError` when the
structure uses a CIF and gives no `z`, since a CIF carries no Z the loader
trusts.

`resolved_cell(spec)` is all six cell parameters, with those the entry's
crystal system leaves out filled in: b equal to a and the three right angles
of a tetragonal cell, gamma of 120 for a hexagonal one. It raises `ValueError`
when the structure gives no cell, and when it uses a CIF without a library
entry and does not give all six, since the crystal system is then unknown
here.

`refine_settings(project, sample)` is the settings a sample is refined with:
the package defaults, the project's `[refine]` over them, and the sample's own
`refine` over that. A table given in part keeps the rest of the table below
it, so `background = { terms = 8 }` changes the terms and keeps the function.
It raises `KeyError` naming the samples there are when there is no such
sample.

`results_dir(project, command, sample)` is
`<root>/results/<command>/<sample key>`, the folder a command's results go in.
It is not created, and it is why Section 9 can say where every file lands
without each command being asked.

`toml_string(value)` quotes and escapes a string as TOML, which is what
`xrdkit init` and `xrdkit instrument` write names and paths with.
`project_template(name)` is the text of a new project file, `[project]`
filled in and every other table commented out, which is what `xrdkit init`
writes and Section 2.1 shows.

### 27.4 The example project read back

The script reads the project of this guide and prints it. Every path is
printed relative to the root, since the absolute ones say only which machine
the guide was built on.

Start a new file named `project_file.py`.

```python
from xrdkit.project import (
    PROJECT_FILE,
    find_project,
    load_project,
    refine_settings,
    resolved_cell,
    resolved_z,
    results_dir,
    toml_string,
)

# Edit these lines for each new project. Nothing below needs changing.
STRUCTURE = "ttb_p4bm"
SAMPLES = ("powder_a", "pellet_a")
COMMAND = "lattice"

path = find_project()
project = load_project(path)
root = project.root


def under_root(value):
    """A path as it reads from the project root, so no machine shows in it."""
    return value.relative_to(root)


print(f"{PROJECT_FILE} found at {under_root(path)}")
print(f"project {toml_string(project.name)}, version {project.version}")
print(
    f"{len(project.instruments)} instruments, {len(project.structures)} structures,"
    f" {len(project.samples)} samples"
)

for key, instrument in project.instruments.items():
    print(f"instrument {key}")
    print(f"  wavelength {instrument.wavelength}, ka2 {instrument.ka2}")
    print(f"  radius {instrument.radius} mm, instprm {under_root(instrument.instprm)}")

for key, spec in project.structures.items():
    print(f"structure {key}")
    print(f"  library {spec.library}, cif {under_root(spec.cif)}")
    print(f"  composition {spec.composition}, z {spec.z}")
    print(
        f"  exchange {spec.exchange}, origin {spec.origin}, fixed {spec.origin_fixed}"
    )
    print(f"  atoms placed on {[site['label'] for site in spec.atoms or ()]}")

for key, sample in project.samples.items():
    print(f"sample {key}")
    print(f"  file {under_root(sample.file)}, instrument {sample.instrument}")
    print(f"  structures {sample.structures}, form {sample.form}")
    print(
        f"  stage {sample.stage!r}, archimedes {sample.archimedes},"
        f" refine {sample.refine}"
    )
```

It prints where the file was found and then every table in it.

```text
xrdkit.toml found at xrdkit.toml
project "ttb-example", version 1
1 instruments, 1 structures, 3 samples
instrument diffractometer
  wavelength (1.540598, 1.544426), ka2 True
  radius 145.0 mm, instprm data\standards\lab6.instprm
structure ttb_p4bm
  library ttb/P4bm, cif cifs\2100720.cif
  composition Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6, z 5
  exchange (('Sr', 'Ba'),), origin None, fixed True
  atoms placed on ['A1', 'A2', 'B1', 'B2']
sample pellet_a
  file data\raw\pellet_a.xrdml, instrument diffractometer
  structures ('ttb_p4bm',), form pellet
  stage None, archimedes None, refine {}
sample pellet_b
  file data\raw\pellet_b.xrdml, instrument diffractometer
  structures ('ttb_p4bm',), form pellet
  stage None, archimedes None, refine {}
sample powder_a
  file data\raw\powder_a.xrdml, instrument diffractometer
  structures ('ttb_p4bm',), form powder
  stage None, archimedes None, refine {'two_theta': (17.0, 99.98)}
```

That is the file of Section 2.2 as the loader sees it. Three things are worth
noticing. The wavelengths are a tuple of two, so the instrument is a K alpha
doublet; `radius` and `instprm` are both there, which is what a refined
displacement and a GSAS-II refinement respectively need. The structure's
`atoms` has become a list of sites carrying their entry labels, `A1`, `A2`,
`B1` and `B2`, which is the placement rule of Section 26.4 in the form the
pipeline takes it. And `powder_a` is the only sample with anything in
`refine`, the per sample range the extraction of Section 8.3 runs over.

Continues `project_file.py`. Add these lines at the end of the file.

```python
spec = project.structures[STRUCTURE]
print(f"resolved_z({STRUCTURE}) = {resolved_z(spec)}")
print(f"resolved_cell({STRUCTURE}) = {resolved_cell(spec)}")
print(f"the table itself gives cell = {spec.cell}")

for key in SAMPLES:
    settings = refine_settings(project, key)
    print(f"refine_settings({key})")
    print(f"  two_theta  {settings.two_theta}")
    print(f"  background {settings.background}")
    print(f"  max_passes {settings.max_passes}")
    print(f"  unsettled  {settings.unsettled}")
    print(f"  followed   {settings.followed}")
print(f"the project's own [refine] two_theta is {project.refine.two_theta}")
for key in SAMPLES:
    print(
        f"results_dir({COMMAND}, {key}) = {under_root(results_dir(project, COMMAND, key))}"
    )
```

It prints what each resolver makes of the file.

```text
resolved_z(ttb_p4bm) = 5
resolved_cell(ttb_p4bm) = {'a': 12.45, 'b': 12.45, 'c': 3.94, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
the table itself gives cell = {'a': 12.45, 'c': 3.94}
refine_settings(powder_a)
  two_theta  (17.0, 99.98)
  background {'function': 'chebyschev-1', 'terms': 6}
  max_passes {'lebail': 60, 'fixed_atoms': 60, 'coordinates': 100, 'occupancies': 100}
  unsettled  {'lebail': 'accept', 'fixed_atoms': 'accept', 'coordinates': 'accept', 'occupancies': 'accept'}
  followed   ()
refine_settings(pellet_a)
  two_theta  None
  background {'function': 'chebyschev-1', 'terms': 6}
  max_passes {'lebail': 60, 'fixed_atoms': 60, 'coordinates': 100, 'occupancies': 100}
  unsettled  {'lebail': 'accept', 'fixed_atoms': 'accept', 'coordinates': 'accept', 'occupancies': 'accept'}
  followed   ()
the project's own [refine] two_theta is None
results_dir(lattice, powder_a) = results\lattice\powder_a
results_dir(lattice, pellet_a) = results\lattice\pellet_a
```

`resolved_cell` is the clearest of the four. The table gives two numbers,
because a tetragonal cell has two free parameters; the resolver returns six,
because that is what GSAS-II and every geometry function want. Nothing in the
file says b equals a; the library entry's crystal system does, and the
resolver reads it from there.

`refine_settings` shows the three layers at work. `powder_a` comes back with
`two_theta` of 17.0 to 99.98 from its own table, while `pellet_a` comes back
with `None`, meaning the scan's own range, and the project file sets no
`[refine]` range at all. Everything else, the background, the passes and the
unsettled rules, is identical between the two because both fall through to the
package defaults. A sample's table is a patch over those layers, not a
replacement of them.

Continues `project_file.py`. Add these lines at the end of the file.

```python
import os

from xrdkit.project import load_project_text

BROKEN = {
    "a sample naming an instrument that is not there": (
        'instrument = "diffractometer"',
        'instrument = "other"',
    ),
    "a structure whose cell is not its entry's": (
        "cell = { a = 12.45, c = 3.94 }",
        "cell = { a = 12.45, b = 12.45, c = 3.94 }",
    ),
}


def without_root(message):
    """A message with the project root taken off, so no machine shows in it."""
    return str(message).replace(f"{root}{os.sep}", "").replace(f"{root}", ".")


text = path.read_text(encoding="utf-8")
print(f"the text loads as {load_project_text(text, path).name}")
for description, (before, after) in BROKEN.items():
    try:
        load_project_text(text.replace(before, after), path)
    except ValueError as error:
        print(f"{description}:")
        print(f"  {without_root(error)}")
```

It loads the file's own text, then breaks it twice.

```text
the text loads as ttb-example
a sample naming an instrument that is not there:
  xrdkit.toml: samples.pellet_a.instrument: no instrument 'other' under instruments; there are diffractometer
a structure whose cell is not its entry's:
  xrdkit.toml: structures.ttb_p4bm.cell: must give exactly a, c, the cell parameters of ttb/P4bm, not a, b, c
```

Both refusals are cross checks between tables, which is what the loader is
for. A sample may name only an instrument the file defines, and the message
lists the ones it does. A structure's `cell` must give exactly the parameters
its library entry's crystal system leaves free, and the message names the
entry. Neither is a TOML error; the file parses perfectly in both cases, and
it is the checking on top of the parse that catches them.

## 28. gsas2

`xrdkit.gsas2` runs GSAS-II. Section 8 is what the two refinement commands do
with it, and this section is the layer underneath: how a job is built, how it
is run, and the handful of things about GSAS-II that a caller has to know
because the module cannot hide them.

There are eleven public names: two environment variable names, the error, the
installation, the finder, the instrument parameter writer, the width formula,
the stage list, the job builder, the runner and the structure edit helper.

### 28.1 Why there is a driver at all

GSAS-II brings its own Python and its own compiled extensions, and cannot be
imported into yours. So nothing here imports GSAS-II. A job, a dict naming an
action and its inputs, is written as JSON; `gsas2_driver.py` is run under the
GSAS-II Python as a subprocess, imports GSASIIscriptable, does the work and
writes its result as JSON; and `run_job` reads that back. The driver imports
nothing from xrdkit, so the two installations never have to agree about
anything but the shape of the JSON.

Four things about that arrangement are worth knowing, because each is a
GSAS-II quirk the module works around and each shows up in a caller's life
sooner or later.

The import path. The driver puts the folder holding the `GSASII` package on
`sys.path` itself, from the `gsas2_home` the job carries, which is why
`Gsas2Install` has a `home` as well as a `python`.

The conda folders on PATH. The GSAS-II Python is a conda environment, and on
Windows its numpy finds the BLAS and LAPACK libraries only through the folders
`conda activate` would put on PATH. Without them the first matrix inversion
kills the process outright, with no Python traceback at all. `run_job` puts
those folders at the front of PATH for the subprocess, which is why you do not
have to activate anything.

The encoding. GSAS-II opens its data files in the locale encoding, and on
Windows that turns the byte order mark at the head of a `.xrdml` file into
characters the XML parser rejects. `run_job` sets `PYTHONUTF8` for the
subprocess unless it is already set. It also writes a two column `.xy` copy of
an `.xrdml` scan into the working folder and passes it as a fallback, for the
driver to read should GSAS-II's own importer fail anyway.

Refinement does not raise. GSAS-II's least squares reports trouble by what it
leaves behind rather than by an exception, so the driver judges each stage
after it runs and records a status, and a stage that fails does not stop the
job. Section 8.16 is how those statuses are read.

### 28.2 find_gsas2 and Gsas2Install

`find_gsas2()` locates the installation and returns a frozen `Gsas2Install`
with `python`, the interpreter, and `home`, the folder containing the `GSASII`
package. `XRDKIT_GSAS2_PYTHON` and `XRDKIT_GSAS2_HOME`, the two names
`GSAS2_PYTHON_VARIABLE` and `GSAS2_HOME_VARIABLE` hold, name them; either left
unset falls back to `~/gsas2main`, whose Python is `python.exe` on Windows and
`bin/python` elsewhere and whose package folder is `GSAS-II`.

It raises `FileNotFoundError` naming both variables and both problems it
found. Call it first in a script that will later refine, so that a missing
installation costs the message rather than the minutes before it.

`Gsas2Error` is a `RuntimeError` carrying the subprocess's `stdout` and
`stderr`, raised when the driver exits non-zero or writes no result.

### 28.3 write_instprm and gsas2_fwhm

`write_instprm(path, caglioti, zero=0.0, x=0.0, y=0.0, shl=0.002,
lam1=1.54056, lam2=1.54439, ratio=0.5, polariz=0.7, radius_mm=None)` writes
the instrument parameter file a powder histogram is read with, and returns the
path. It is what `xrdkit instrument` writes its starting file with from the
Caglioti fit of Section 22.

The unit conversion in it is the thing to know. xrdkit's U, V and W give the
squared FWHM in degrees squared; GSAS-II's give the variance of the Gaussian
component in centidegrees squared, its FWHM being the square root of eight ln
two times sigma. Each is therefore multiplied by ten thousand over eight ln
two, about 1803, which is why the numbers Section 3.2 prints from the width
fit and from the refinement look nothing like each other. The whole fitted
width is taken as Gaussian, which is exact while `x` and `y` are zero and a
starting point otherwise.

`radius_mm` writes the goniometer radius line GSAS-II writes itself. Without
it a text pattern leaves GSAS-II at its default of two hundred millimetres,
and a refined specimen displacement comes out wrong by the ratio of the two.

`gsas2_fwhm(two_theta, u, v, w, x, y, shl=0.0, z=0.0)` is the line width
GSAS-II's own refined parameters imply, in degrees, taking them in GSAS-II's
units: U, V and W the Gaussian variance in centidegrees squared, held at no
less than 0.001, and X, Y and Z the Lorentzian FWHM in centidegrees. The two
are combined by the Thompson, Cox and Hastings quintic exactly as GSAS-II
does. `shl` is accepted so that a refined set can be passed whole, and is
ignored: axial divergence makes a line asymmetric without entering its width.

### 28.4 standard_stages, build_refine_job and the stage language

A stage is a dict with a `name` and the refinement flags it switches on. The
driver carries every flag forward to the stages after it, so a stage names
only what it adds, and a stage the driver rejects holds what it alone refined
from then on. The flags are the whole vocabulary of a refinement here:
`background`, `scale`, `zero`, `displacement`, `instrument` (a list of
GSAS-II's parameter names), `cell`, `le_bail`, `phase_fractions`, `size`,
`mustrain`, `overall_uiso`, `atoms`, `atom_flags`, `uiso_groups`,
`coordinates`, `origin` and `occupancies`. Section 29 builds the four
sequences the commands use out of them.

`standard_stages(background_type, background_terms)` is the usual sequence for
a script to edit: background and scale; zero; cell; U, V and W; X and Y;
SH/L. The list is new on every call, so a stage can be dropped or a flag added
freely. It is what `xrdkit instrument` runs, less its cell stage.

`build_refine_job(gpx, stages, ...)` assembles the job. With `data_file`,
`instprm` and `phases` the project is created at `gpx` first; without all
three, `gpx` must already hold one. The stages are checked here, by the
driver's own rules, so a mistake is reported before GSAS-II starts. Every
path is made absolute against the current working directory, because the
driver runs in a working folder of its own.

Its arguments are many and Section 8 exercises most of them; the ones a hand
built refinement reaches for first are `limits`, the two theta range;
`cycles`, the most least squares cycles a stage may take; `broadening`, the
sample size and microstrain to start from and hold until a stage refines
them; `export_prefix`, what the exported files are named from; and the pass
pair below.

`max_passes` and `pass_tolerance` are the pair worth understanding. GSAS-II
stops a refinement when its own convergence test is met, which for a Le Bail
extraction is well short of settled. Given `max_passes`, the driver refines a
stage again and again, up to that many times, until no parameter moves by more
than `pass_tolerance` esds from one pass to the next. Without it a stage is
refined once and carries no `passes` in the result. Section 8.5 is what the
passes mean for a run and Section 29 gives the pipeline's own values.

It raises `Gsas2Error` when an input file does not exist, and `ValueError` or
`TypeError` when a stage, the broadening, an atom edit or any of the starting
values is malformed.

### 28.5 run_job and structure_edits

`run_job(job, workdir, install=None)` writes the job to
`workdir/<action>_job.json`, runs the driver on it as `python -B driver
job.json` in `workdir`, keeps the output in `workdir/<action>.log` and returns
the result the driver wrote. The `-B` keeps bytecode out of the GSAS-II
installation. It raises `Gsas2Error` carrying the driver's stderr when the
driver exits non-zero or leaves no result.

The result is a dict. `completed` says whether every stage ran, `final_from`
names the stage the final model was taken from, `stages` is one record per
stage with its `rwp`, `rp`, `gof`, `n_variables`, `passes` and `status`,
`final` holds the refined `instrument` and `phases`, `rejected` lists the
stages rolled back, `undetermined` the parameters the driver could not
determine, and `exports` the files written.

Two readers go with it. `stage_status(stage)` is one stage's status, and
`accepted_stages(result)` the stages refined and kept, which is what Section
29.4 counts. `stage_statuses(result)` gives each stage as a row with its
`name`, `status`, `reason`, `passes`, `rwp`, `gof` and `undetermined`, the
`reason` being why a stage was rejected or failed, the new flags of one
flagged, or the passes and largest remaining move of one unsettled.

`structure_edits(refined, base)` turns a refined structure back into the atom
edits that set it up again in a new project read from the same CIF, so that
one refinement can start where another ended. It is how the Rietveld modes of
Section 29 carry their atoms from mode to mode. It raises `ValueError` when an
atom the CIF lacks shares no site with one it has.

Four more names in this module render a run as markdown rather than read it.
`stage_status_table(result)` is the stage table of Section 8.14 as markdown
lines; `log_tail(path, lines=40)` is the last lines of a GSAS-II log, empty
when there is none; `failure_markdown(title, result, error, log, intro=())`
is the `failure.md` a mode that did not finish leaves, from the error, the
stages its run got through and that log tail; and `summary_markdown(title,
entries, columns=(), intro=())` is the `summary.md` written over a sequence
of modes. `run_sequence` of Section 29.4 writes both files, Section 8.14
reads them and Section 8.16 is what a failure page says. They render the
frame of a page; the per mode pages themselves are the `writeup` module of
Section 30.

### 28.6 A refinement of your own, stage by stage

The script below is a hand built refinement: three stages on `powder_a`
against the CIF of Section 24, starting from the cell the Le Bail extraction
of Section 8.3 gave. It is the shape to copy when a refinement needs a stage
sequence the modes do not offer.

Everything the commands do around this is Section 8 and is not here: choosing
the modes, carrying each one's result into the next, putting the nominal
composition on the atoms, the write ups, the sanity checks and the rollback.
In particular this refinement uses the CIF's own composition rather than the
sample's, so its atoms are the published strontium barium niobate rather than
the doped formula of Section 2.2, and its R factors are not comparable with
Section 8.8's. What it demonstrates is the machinery, not a result to quote.

Everything it writes goes under `results/library/gsas2_powder_a`, the GSAS-II
project, the logs, the job and result JSON and the exports together, so that
nothing lands in the folders the commands own.

Start a new file named `refine_stages.py`.

```python
from pathlib import Path

from xrdkit.gsas2 import build_refine_job, find_gsas2, run_job, stage_statuses
from xrdkit.io import read_scan
from xrdkit.project import find_project, load_project, refine_settings

# Edit these lines for each new sample. Nothing below needs changing.
SAMPLE = "powder_a"
PHASE = "ttb_p4bm"
CIF_FILE = "cifs/2100720.cif"
START_CELL = (12.4740, 12.4740, 3.9318, 90.0, 90.0, 90.0)
OUT = Path("results/library/gsas2_powder_a")
MAX_PASSES = 20
STAGES = [
    {"name": "scale and background", "background": True, "scale": True},
    {"name": "zero and cell", "zero": True, "cell": True},
    {"name": "profile", "instrument": ["U", "V", "W"]},
]

install = find_gsas2()
print(f"GSAS-II found: {install.python.name} beside {install.home.name}")

project = load_project(find_project())
root = project.root
sample = project.samples[SAMPLE]
instrument = project.instruments[sample.instrument]
settings = refine_settings(project, sample)
scan = read_scan(sample.file)
limits = settings.two_theta or (scan.start_angle, scan.end_angle)
print(
    f"{SAMPLE}: {scan.two_theta.size} points from {scan.start_angle:.2f}"
    f" to {scan.end_angle:.2f} degrees"
)
print(f"refining {limits[0]} to {limits[1]} degrees, background {settings.background}")
print(f"instrument file {instrument.instprm.relative_to(root)}")
print(f"stages {[stage['name'] for stage in STAGES]}")

job = build_refine_job(
    OUT / f"{SAMPLE}.gpx",
    STAGES,
    data_file=sample.file,
    instprm=instrument.instprm,
    phases=[{"cif": CIF_FILE, "name": PHASE, "cell": list(START_CELL)}],
    limits=limits,
    cycles=10,
    max_passes=MAX_PASSES,
    pass_tolerance=0.1,
    broadening={PHASE: {"size": 1.0, "mustrain": 0.0, "lgmix": 1.0}},
    export_prefix=OUT / SAMPLE,
)
print(
    f"job action {job['action']}, {len(job['stages'])} stages, {job['cycles']} cycles"
)
result = run_job(job, OUT / "work")
print(f"completed {result['completed']}, final model from {result['final_from']}")
```

It prints the installation, the inputs it took from the project file, the job
it built and whether GSAS-II finished it.

```text
GSAS-II found: python.exe beside GSAS-II
powder_a: 4141 points from 10.01 to 99.98 degrees
refining 17.0 to 99.98 degrees, background {'function': 'chebyschev-1', 'terms': 6}
instrument file data\standards\lab6.instprm
stages ['scale and background', 'zero and cell', 'profile']
job action refine, 3 stages, 10 cycles
completed True, final model from profile
```

Note where each input came from. The scan, the instrument parameter file, the
two theta range and the background all came out of the project file through
Section 27's loader and resolvers, not out of the settings at the top of the
script; only the phase, the CIF and the start cell are the script's own. That
is the division worth keeping in a script of your own, because it means a
sample is described in one place.

The run took about forty seconds.

Continues `refine_stages.py`. Add these lines at the end of the file.

```python
print("  stage                 Rwp     Rp    GOF  variables  passes  status")
for stage in result["stages"]:
    print(
        f"  {stage['name']:20s} {stage['rwp']:5.3f}  {stage['rp']:5.3f}"
        f"  {stage['gof']:5.3f}  {stage['n_variables']:9d}"
        f"  {len(stage['passes']):6d}  {stage['status']}"
    )
for row in stage_statuses(result):
    if row["reason"]:
        print(f"  {row['name']}: {row['reason']}")

phase = result["final"]["phases"][0]
cell, esd = phase["cell"], phase["cell_esd"]
print(f"a = {cell['length_a']:.4f} +/- {esd['length_a']:.4f} angstrom")
print(f"c = {cell['length_c']:.4f} +/- {esd['length_c']:.4f} angstrom")
print(f"V = {cell['volume']:.3f} +/- {esd['volume']:.3f} cubic angstrom")
zero = result["final"]["instrument"]["Zero"]
print(f"zero {zero['value']:.4f} +/- {zero['esd']:.4f} degrees")
for key in ("U", "V", "W"):
    entry = result["final"]["instrument"][key]
    print(f"{key} {entry['value']:9.4f} +/- {entry['esd']:.4f} centidegrees squared")
print("exports")
for name, written in sorted(result["exports"].items()):
    files = written.values() if isinstance(written, dict) else [written]
    for one in files:
        print(f"  {name}: {Path(one).relative_to(root)}")
```

It prints every stage, the refined cell and instrument, and the files GSAS-II
exported.

```text
  stage                 Rwp     Rp    GOF  variables  passes  status
  scale and background 8.383  5.787  2.890          7       2  clean
  zero and cell        7.791  5.588  2.687         10       6  clean
  profile              5.439  4.222  1.876         13      14  clean
a = 12.4755 +/- 0.0006 angstrom
c = 3.9322 +/- 0.0002 angstrom
V = 612.003 +/- 0.073 cubic angstrom
zero -0.0343 +/- 0.0013 degrees
U  444.9634 +/- 54.6354 centidegrees squared
V -209.2662 +/- 42.2414 centidegrees squared
W   49.3307 +/- 7.4550 centidegrees squared
exports
  histogram: results\library\gsas2_powder_a\powder_a_histogram.csv
  instprm: results\library\gsas2_powder_a\powder_a.instprm
  reflections: results\library\gsas2_powder_a\powder_a_reflections_ttb_p4bm.csv
```

Read the stage table downwards. Rwp falls from 8.383 to 7.791 to 5.439 per
cent as the three stages free seven, then ten, then thirteen parameters, and
the goodness of fit falls from 2.890 to 1.876 with it. The largest single
gain is the profile stage, which frees U, V and W: the instrument file was
measured on a standard in Section 3, and this sample's lines are broader than
the standard's, so until those three are free the calculated pattern has the
wrong widths everywhere. That is a demonstration of the machinery and not good
practice, incidentally. Refining the instrument terms against a sample throws
away the separation Section 3.1 exists to make; the right way to account for
sample broadening is the size and microstrain of Section 29, which is what the
commands refine.

The passes column is `max_passes` at work. The first stage settled in two
passes and the last needed fourteen, which is what freeing three correlated
width terms at once costs. Had the cap of twenty been reached, the stage would
have come back `unsettled` rather than `clean`.

The cell comes out at a = 12.4755(6) and c = 3.9322(2) angstrom, against the
12.4740(3) and 3.9318(1) the Le Bail extraction of Section 8.3 gave and this
script started from. The two agree in c and differ in a by 0.0015 angstrom,
five of the extraction's esds, which is the difference between fitting
intensities freely and calculating them from a structure whose composition is
not quite the sample's.

The exports are the three files every refinement writes: the observed and
calculated pattern, the reflection list per phase, and the refined instrument
parameter file. They are what `plot_rietveld` of Section 15.6 draws from, and
Section 9.10 lists them as the commands write them.

## 29. pipeline

`xrdkit.pipeline` is `xrdkit lebail` and `xrdkit rietveld`: the four modes,
the stages each one runs, what each takes from the result of the one before,
and the run itself. Section 8 is the commands and how to read what they
report, and this section is the functions they are made of.

The public surface divides in three. Four stage builders turn a project's
settings into a list of stages, and need no GSAS-II at all. `start_from_result`
and `StartPoint` read what a previous mode left. `run_mode` and `run_sequence`
do the run, with `Options` going in and `Outcome` coming out.

`MODES` is the four in order, `lebail`, `fixed_atoms`, `coordinates` and
`occupancies`. `CYCLES` is 10, the most least squares cycles of one GSAS-II
refinement; `LE_BAIL_CYCLES` is 10, the extraction-only cycles run whenever a
stage switches Le Bail extraction on; and `PASS_TOLERANCE` is 0.1, the esds a
parameter may move between passes and still count as settled, which is the
`pass_tolerance` of Section 28.4. `START_BROADENING` is where a Le Bail
refinement starts its size and microstrain and holds them until its own stages
reach them.

### 29.1 The stage language, and why a stage names only what it adds

The builders give stages in the form Section 28.4 describes, and they rely on
the driver carrying every flag forward. A stage therefore names only what it
switches on, and the scale, zero and displacement switches are stated
explicitly in the first stage of every mode so that the mode does not inherit
a guess.

That accumulation is what makes rollback work. A stage the driver rejects
holds what that stage alone refined from then on, so the run continues with
the parameters of the stages before it and nothing the rejected stage touched.
Section 8.16 is how a rejected stage is read, and Section 8.8 shows one: the
`O sites` stage of the coordinates mode.

With two phases or more, the stage that would free the histogram scale frees
the phase fractions instead and holds the scale, since the two are the same
quantity counted twice.

### 29.2 The four stage builders

`lebail_stages(background, phases=1, displacement=False, mustrain_test=True)`
gives five stages: `background and scale`, with extraction on for every phase
and the zero and displacement held; `zero`, or `displacement` when
`displacement` is true; `cell`; `size`; and `microstrain`, which
`mustrain_test` drops. No instrument parameter but the zero and no atomic
parameter is freed anywhere in it.

`fixed_atoms_stages(background, phases=1, displacement=False, mustrain=False,
preferred_orientation=None)` gives `scale and background` with extraction off,
`zero and cell` (or `displacement and cell`), `size` or `size and
microstrain`, `overall Uiso`, and a `preferred orientation` stage when an
[h, k, l] axis is given.

`coordinates_stages(plan, background, phases=1, displacement=False,
mustrain=False)` gives `profile`, then `Uiso groups`, the plan's groups in
place of one overall Uiso, then one stage per kind of site named `<kind>
sites`, in the order of the plan's kinds. A kind with nothing to free has no
stage, and the origin site's coordinate along its axis is held throughout.

`occupancy_stages(plan, background, phases=1, displacement=False,
mustrain=False)` gives `profile and Uiso` and then one stage per exchange
group, named `<kind> site occupancies`, trading the occupancies of that
group's elements between its sites with each element's total held. The
coordinates are held throughout. It raises `PipelineError` when the plan names
no phase or an exchange group's sites are of more than one kind.

The last two take a `plan`, the site plan `site_setup` of Section 26.2 gives,
with the phase name added as `plan["phase"]`.

### 29.3 The stages of the example sample

The script builds all four sequences for the sample and structure of Section
2.2. It needs no GSAS-II, and it is worth running before a refinement rather
than after, because it is the cheapest way to see what a mode is going to do.

The plan comes from `site_setup` as Section 26.2 built it, with one change:
the atoms are read out of a result an earlier refinement wrote rather than
parsed from the CIF. The driver records every atom at every stage, so a result
JSON carries its atoms in exactly the form `site_setup` takes, which saves the
CIF parsing of Section 26.1 whenever a refinement has already been run.

Start a new file named `refine_pipeline.py`.

```python
import json
from pathlib import Path

from xrdkit.library import load_entry
from xrdkit.pipeline import (
    CYCLES,
    MODES,
    PASS_TOLERANCE,
    coordinates_stages,
    fixed_atoms_stages,
    lebail_stages,
    occupancy_stages,
)
from xrdkit.project import find_project, load_project, refine_settings
from xrdkit.structure import site_setup

# Edit these lines for each new sample. Nothing below needs changing.
SAMPLE = "powder_a"
STRUCTURE = "ttb_p4bm"
ATOMS_FROM = "results/rietveld/powder_a/powder_a_fixed_atoms_result.json"

project = load_project(find_project())
root = project.root
sample = project.samples[SAMPLE]
spec = project.structures[STRUCTURE]
entry = load_entry(spec.library)
settings = refine_settings(project, sample)
background = settings.background
print(f"modes {MODES}, {CYCLES} cycles, pass tolerance {PASS_TOLERANCE} esds")
print(f"background {background}, form {sample.form}")

# The plan of Section 26.2, from the atoms an earlier refinement recorded.
atoms = json.loads(Path(ATOMS_FROM).read_text(encoding="utf-8"))
atoms = atoms["final"]["phases"][0]["atoms"]
placed = {site["label"]: site["atoms"] for site in spec.atoms or ()}
present = {atom["label"] for atom in atoms}
sites = []
for site in entry.sites:
    on_site = {
        label: element
        for label, element in (placed.get(site.label) or {site.label: "O"}).items()
        if label in present
    }
    sites.append(
        {
            "name": next(iter(on_site)),
            "label": site.label,
            "atoms": on_site,
            "wyckoff": site.wyckoff,
            "kind": site.kind,
        }
    )
by_label = {site["label"]: site["name"] for site in sites}
groups = {}
for site in entry.sites:
    groups.setdefault(site.uiso_group or site.label, []).append(by_label[site.label])
plan = site_setup(
    {
        "name": STRUCTURE,
        "library": entry.name,
        "sites": sites,
        "free_coordinates": {},
        "uiso_groups": [
            {"name": name, "sites": members} for name, members in groups.items()
        ],
        "origin": {"site": by_label[entry.origin_site], "axis": "z"},
        "exchange": {
            "elements": list(spec.exchange[0]),
            "sites": [by_label["A1"], by_label["A2"]],
        },
    },
    atoms,
)
plan["phase"] = STRUCTURE

displacement = sample.form == "pellet"
for name, stages in (
    ("lebail", lebail_stages(background, displacement=displacement)),
    ("fixed_atoms", fixed_atoms_stages(background, displacement=displacement)),
    ("coordinates", coordinates_stages(plan, background, displacement=displacement)),
    ("occupancies", occupancy_stages(plan, background, displacement=displacement)),
):
    print(f"{name}: {len(stages)} stages")
    for stage in stages:
        switches = ", ".join(key for key in stage if key != "name")
        print(f"  {stage['name']:22s} {switches}")
```

It prints the constants, then each mode with its stages and the flags each
stage switches on.

```text
modes ('lebail', 'fixed_atoms', 'coordinates', 'occupancies'), 10 cycles, pass tolerance 0.1 esds
background {'function': 'chebyschev-1', 'terms': 6}, form powder
lebail: 5 stages
  background and scale   background, scale, le_bail, zero, displacement
  zero                   zero, displacement
  cell                   cell
  size                   size
  microstrain            mustrain
fixed_atoms: 4 stages
  scale and background   background, scale, le_bail, zero, displacement
  zero and cell          zero, displacement, cell
  size                   size
  overall Uiso           overall_uiso
coordinates: 5 stages
  profile                background, scale, le_bail, zero, displacement, cell, size
  Uiso groups            overall_uiso, uiso_groups
  A sites                coordinates
  B sites                coordinates, origin
  O sites                coordinates
occupancies: 2 stages
  profile and Uiso       background, scale, le_bail, zero, displacement, cell, size, overall_uiso, uiso_groups
  A site occupancies     occupancies
```

Those stage names are the ones Section 8.3 and Section 8.8 report, in the same
order, because the commands print what these builders return. Read the flag
columns across and the accumulation is visible: `background`, `scale` and
`le_bail` appear in the first stage of the Le Bail mode and again in the first
stage of every Rietveld mode, because each mode is a fresh GSAS-II project and
has to state them again, while `cell` appears once in the Le Bail mode and is
carried by the driver into the three stages after it.

Two details in the Rietveld modes repay attention. The `B sites` stage of the
coordinates mode carries an `origin` flag that the `A sites` and `O sites`
stages do not: the origin site is `B1`, and its z is held as soon as that
kind's coordinates are freed. And the occupancies mode has one exchange stage,
`A site occupancies`, because the project file's `exchange` names one group,
strontium and barium, and both sit on A sites.

### 29.4 start_from_result, run_mode and one mode run

`start_from_result(path, stage=None)` reads a result JSON and returns a
`StartPoint`: each phase as a `PhaseStart` with its name, cell, size,
microstrain, fraction and atoms; the histogram's `zero`, `displacement`,
`scale`, `background` and `limits`; the run's `start_model`; and the `stage`
taken, the last accepted by default. `cell`, `size`, `microstrain` and `atoms`
are properties reaching the first phase. It raises `PipelineError` when the
file is missing or not a result, when the stage was not accepted, or when a
value it needs is not recorded.

`Options` is what a run takes besides the project file: `cell` and `zero` to
start from, over the project's; `displacement`, whether the specimen
displacement is refined in place of the zero, by default from the sample's
form; `mustrain`; `preferred_orientation`; `max_passes`, a cap over the
project's; and `out`, the folder every file goes to in place of
`results/lebail/<key>` and `results/rietveld/<key>`.

`out` has to be an absolute path. `mode_paths` takes it as given and the
GSAS-II driver runs in a working folder of its own, so a relative `out` sends
the result JSON somewhere neither of them expects and the run fails when the
driver tries to write it. The docstring does not say so; the script below
makes its own absolute against the project root and says why in a comment.

`run_mode(project, sample, mode, options=None, reporter=None)` runs one mode
and returns an `Outcome`: the `mode`, the names of its `accepted` stages, its
`final` model, its `residuals`, the `undetermined` parameters, the `paths` it
wrote and the `error` that stopped it, `None` when it finished. `reporter` is
a callable taking a line of text, which is how the commands print their
progress. It raises `PipelineError` when an input is missing or a run leaves
no usable model, and `Gsas2Error` when GSAS-II itself failed; either carries
the run's `result` and `log`.

`run_sequence(project, sample, modes, options=None, reporter=None)` runs
several modes in order, each from the one before, stopping at the first that
fails and writing that mode a `failure.md` from the error and the stages its
own run got through, then a `summary.md` over the modes run. That is
`xrdkit rietveld` entire, and Section 8.8 is the run of it on this sample;
this guide does not call it, because it would write over the results Part I
recorded.

Five more public names are the parts a run is assembled from, and a script
reaches for them only to look at what a run would do before doing it.
`resolve_inputs(project, sample, options=None)` is everything a run takes
from the project file, as an `Inputs`: the sample, its instrument and
instrument parameter file, the scan's range and the range to refine over, the
resolved `Refine` settings, the phases as a tuple of `PhaseInput`, and
whether the displacement is refined and the zero to start from. A
`PhaseInput` is one phase as the project gives it, its structure key, spec,
library entry, composition, formula units and start cell with where that cell
came from. `mode_paths(project, sample, mode, options=None)`, which Section
29.4's block uses, says where a mode's files go without running anything.
`exchange_edits(atoms, plan)` adds an atom at occupancy zero for each
exchanged element a site of an exchange group lacks, so that the element has
somewhere to move to; the occupancies mode calls it before its own stages.
And `HELD_INSTRUMENT` is the instrument parameters a run never refines, U, V,
W, X, Y, Z and SH/L: they belong to the instrument file of Section 3, and
`run_mode` checks after every run that none of them came back changed.

The block below runs one mode instead: `fixed_atoms` on `powder_a`, starting
from the Le Bail result of Section 8.3, with everything written to
`results/library/pipeline_powder_a`. A mode looks for the previous mode's
result in its own output folder, so the Le Bail result is copied there first;
that copy is the only reason the block touches `results/lebail` at all, and it
reads it rather than writing to it.

Continues `refine_pipeline.py`. Add these lines at the end of the file.

```python
import shutil

from xrdkit.pipeline import Options, mode_paths, run_mode, start_from_result

LEBAIL_RESULT = "results/lebail/powder_a/powder_a_lebail_result.json"
# Options.out must be absolute: the GSAS-II driver runs in a working folder
# of its own, so a relative path would be written there instead.
OUT = root / "results/library/pipeline_powder_a"

start = start_from_result(LEBAIL_RESULT)
print(f"start_from_result: stage {start.stage!r}, zero {start.zero:.4f} degrees")
print(f"  cell a {start.cell['a']:.4f}, c {start.cell['c']:.4f} angstrom")
print(f"  size {start.size:.4f} micron, microstrain {start.microstrain:.0f}")
print(f"  scale {start.scale:.4f}, limits {start.limits}")
print(f"  background {start.background['function']}, {start.background['terms']} terms")

options = Options(out=OUT)
paths = mode_paths(project, sample, "fixed_atoms", options)
before = mode_paths(project, sample, "lebail", options)["result"]
before.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(LEBAIL_RESULT, before)
print(f"the mode starts from {before.relative_to(root)}")

outcome = run_mode(project, sample, "fixed_atoms", options, reporter=print)
print(f"mode {outcome.mode}, error {outcome.error}")
print(f"accepted {outcome.accepted}")
print(
    f"Rwp {outcome.residuals['rwp']:.3f} per cent,"
    f" reduced chi squared {outcome.residuals['reduced_chi_squared']:.3f}"
)
result = json.loads(paths["result"].read_text(encoding="utf-8"))
print("  stage                 passes  status")
for stage in result["stages"]:
    print(f"  {stage['name']:22s} {len(stage['passes']):6d}  {stage['status']}")
cell = outcome.final["phases"][0]["cell"]
print(f"a = {cell['length_a']:.4f}, c = {cell['length_c']:.4f} angstrom")
```

It prints the start point, the run as the reporter hears it, and what the mode
came to.

```text
start_from_result: stage 'microstrain', zero -0.0357 degrees
  cell a 12.4740, c 3.9318 angstrom
  size 0.2051 micron, microstrain 1011
  scale 0.0004, limits (17.0, 99.98)
  background chebyschev-1, 6 terms
the mode starts from results\library\pipeline_powder_a\powder_a_lebail_result.json
powder_a: fixed_atoms started, 4 stages, at most 60 passes each
powder_a: fixed_atoms: scale and background clean
powder_a: fixed_atoms: zero and cell clean
powder_a: fixed_atoms: size clean
powder_a: fixed_atoms: overall Uiso clean
mode fixed_atoms, error None
accepted ['scale and background', 'zero and cell', 'size', 'overall Uiso']
Rwp 4.290 per cent, reduced chi squared 2.190
  stage                 passes  status
  scale and background        2  clean
  zero and cell               3  clean
  size                        2  clean
  overall Uiso                2  clean
a = 12.4735, c = 3.9318 angstrom
```

The run took about twenty seconds and reproduces the fixed atoms line of
Section 8.8 exactly: four stages accepted, in two, three, two and two passes,
an Rwp of 4.290 per cent and a reduced chi squared of 2.190. It has to, since
`xrdkit rietveld` calls this function with these arguments; what changes is
only where the files land.

The start point is worth reading beside the run. `start_from_result` with no
stage named gives the last accepted stage, `microstrain`, and the cell, size
and microstrain that stage left. `run_mode` does not take that one: its
docstring says the fixed atoms mode starts from the Le Bail result's `size`
stage, because the microstrain stage of a Le Bail run is a test rather than a
measurement, as Section 8.3 explains. Naming a stage is therefore not an
unusual thing to do, and the default is not always the right start.

The zero of minus 0.0357 degrees, the cell of 12.4740 and 3.9318 angstrom and
the limits of 17 to 99.98 all came out of the Le Bail result rather than the
project file. That is the chain Section 8.9 describes: each mode is a separate
GSAS-II project, and the only thing joining them is the result JSON.

### 29.5 Reading a project's refinements back

The last block does no refining at all. It walks the two results folders,
reads every result JSON that belongs to a sample of the project, and prints
one line per sample and mode. This is the shape of a script that tabulates a
composition series: the refinements are run once by the commands, and
everything afterwards is reading files.

Continues `refine_pipeline.py`. Add these lines at the end of the file.

```python
from xrdkit.gsas2 import accepted_stages

FOLDERS = ("results/lebail", "results/rietveld")

found_results = []
for folder in FOLDERS:
    for found in Path(folder).glob("*/*_result.json"):
        key = found.parent.name
        if key not in project.samples:
            print(f"  {found.parent} is not a sample of the project, skipped")
            continue
        mode = found.name[len(f"{key}_") : -len("_result.json")]
        found_results.append((key, MODES.index(mode), mode, found))

print("  sample    mode          stages  Rwp    chi2   a         c")
for key, _, mode, found in sorted(found_results):
    record = json.loads(found.read_text(encoding="utf-8"))
    kept = accepted_stages(record)
    last = kept[-1]
    cell = record["final"]["phases"][0]["cell"]
    print(
        f"  {key:9s} {mode:13s} {len(kept):5d}  {last['rwp']:5.3f}"
        f"  {last['gof'] ** 2:5.3f}  {cell['length_a']:.4f}"
        f"  {cell['length_c']:.4f}"
    )
```

It prints one line per result found.

```text
  results\lebail\powder_a_full is not a sample of the project, skipped
  sample    mode          stages  Rwp    chi2   a         c
  powder_a  lebail            5  3.888  1.799  12.4740  3.9318
  powder_a  fixed_atoms       4  4.290  2.190  12.4735  3.9318
  powder_a  coordinates       4  4.158  2.062  12.4734  3.9318
  powder_a  occupancies       2  4.088  1.990  12.4736  3.9319
```

Four results for one sample, in mode order, and the numbers are those of
Sections 8.3 and 8.8. Rwp rises from the Le Bail extraction's 3.888 to the
fixed atoms mode's 4.290 and then falls through 4.158 to 4.088, which is the
shape to expect: an extraction has an intensity free for every reflection and
should always fit better than a structure, so the step from 3.888 to 4.290 is
the cost of calculating the intensities instead, and the fall after it is the
structure being improved.

The skipped folder is the other half of the lesson.
`results/lebail/powder_a_full` is the run of Section 8.4 over the whole scan,
made with `--out`, and its name
is not a sample key, so the loop passes over it rather than reporting a sample
that does not exist. A results folder accrues runs that are not the current
ones, and a script reading them back has to say which it is using.

`accepted_stages` from Section 28.5 is what counts the stages, so a mode with
a rejected stage shows fewer here than its builder produced: the coordinates
mode shows four of the five stages Section 29.3 built, the `O sites` stage
having been rolled back, which Section 8.8 reports and Section 8.16 reads.

## 30. The modules you are not expected to call

Four modules of the package are not documented as a library, because nothing
you would write needs them. Each exists so that something else in Sections 12
to 29 can be simple, and each is named here with what it is for, what calls
it, and where to look if you turn out to need one of its functions after all.

Nothing in this section has a worked example, and none of it is a stable
promise in the way the rest of Part II is: these are the seams of the package,
and they may move.

### 30.1 cli

`xrdkit.cli` is argparse over the library. Every `xrdkit` command in Part I is
one function here that reads the arguments, calls the project loader of
Section 27 and then the library, and prints what came back. There is no
analysis in it; anything a command computes, it computes by calling something
this guide has already documented, which is why Part II never has to say
"unlike the command".

Three names are public: `main(argv=None)`, the entry point the `xrdkit`
executable runs; `build_parser()`, the whole argparse parser, which is what
`--help` prints and what the tests read the options off; and `CommandError`,
raised for anything a user can fix, which `main` catches and prints as one
line on the error stream before exiting non-zero.

The reason not to call it is that a command's arguments are a user interface
and the library's are an API. If you want a command's behaviour in a script,
call the functions it calls; if you want the command itself, run it. Section 9
lists every file each command writes, which is what a script wrapping one
actually needs to know.

### 30.2 config

`xrdkit.config` reads the settings tables that describe a structure's sites.
It is the older half of what is now the project file: `xrdkit.toml` is read by
`project.py` of Section 27, and `project.py` calls this module for the parts
of a structure table that describe sites and compositions.

Ten names are public. `load_config(path)`, `validate_config(config)` and
`sample_settings(config, sample_id)` read and check a standalone settings
file, which is the arrangement the older guide used and which the project file
has replaced. `read_sites(value, where)` and `read_library_atoms(value,
where)` read the two shapes a structure's `atoms` may take, the list of a
CIF's sites and the table of a library entry's, and they are what gives
`StructureSpec.atoms` of Section 27.2 its form. `check_composition` and
`host_elements` settle which element hosts an added one, the rule Section 26.4
applies and Section 8.11 explains. `wyckoff_multiplicity(wyckoff)` is the
multiplicity of a Wyckoff position such as `"8d"`, which Section 26.1 uses to
give each atom its multiplicity. `SITE_KIND` is the pattern a kind of site
must match, and `ConfigError` is what any of them raises.

Of those, `wyckoff_multiplicity` is the one worth knowing: it is a plain
function of a string and Section 26 calls it directly. The rest are reached
through `load_project` of Section 27.1, which is the supported way in.

### 30.3 gsas2_driver

`xrdkit.gsas2_driver` is the script that runs inside GSAS-II. Section 28.1
says why it exists: GSAS-II brings its own Python and its own compiled
extensions, so nothing in xrdkit imports it, and instead `run_job` of Section
28.5 writes a job as JSON and runs this file under the GSAS-II Python as a
subprocess. It imports nothing from xrdkit, because xrdkit is not installed
there.

It declares no `__all__` and is not re-exported from the package, so it has no
public names in the sense Section 31 uses. `gsas2.py` does import a handful of
its checkers, `accumulate_stages`, `check_broadening`, `check_atom_edits` and
their neighbours, so that a job's stages are checked in your Python before
GSAS-II is started at all, which is what lets `build_refine_job` of Section
28.4 report a malformed stage in a tenth of a second rather than after a
minute of loading.

Everything a refinement can ask for is defined here rather than in `gsas2.py`:
the two actions, `create` and `refine`; the stage flags Section 28.4 lists;
the pass loop, the sanity check and the rollback of Section 8.16; and the
shape of the result JSON that Sections 28.5 and 29.4 read back. The module
docstring is the reference for all of it, and it is the place to look when a
stage does not do what you expected.

### 30.4 writeup

`xrdkit.writeup` turns a result JSON into the markdown pages Section 8.14
reads. It is called by `pipeline.py` of Section 29 after every mode, never by
a user: its input is the result of a run that has just finished, which a
script would have to have run to have.

Five names are public. `lebail_markdown`, `fixed_atoms_markdown` and
`structure_markdown` render a mode's page, the last one serving both the
coordinates and the occupancies modes, each taking the result, the run's
inputs and its site plan. `with_esd(value, esd, digits=6)` writes a value with
its esd in brackets, the esd to one figure or two when it begins with a one,
and marks a value with no esd as fixed; it is why every number in those pages
is written the same way. `relative(value, root)` rewrites every absolute path
in a result, however deeply nested, as a path under the project root with
forward slashes, which is what keeps a machine's folder layout out of the
pages a run leaves behind.

The four markdown helpers `pipeline.py` uses around these live in `gsas2.py`
instead and are in Section 28.5, which is the seam to be aware of: this module
renders a mode, and those render the frame of a run.

## 31. Every public name

The table lists every public name in the package, alphabetically, with the
module that defines it and the section of Part II that documents it. Public
means a name in a module's `__all__`, and every name `xrdkit/__init__.py`
re-exports, so `from xrdkit import find_peaks` and
`from xrdkit.peaks import find_peaks` reach the same function and the table
lists it once, under the module that defines it. The fields of a dataclass are
not listed separately; they are in the table of the section that documents the
class.

Two hundred and one names are in it. `__version__`, the package version
string, is the only one with no section of its own.

Three things are worth knowing before reading it.

A name in the table is one this guide documents, not one the package promises
forever. The modules of Section 30 are in it too, since they have `__all__` of
their own, and their rows point at Section 30 rather than at a worked example.

`TetragonalCell` is the one deprecated name. It is kept so that older scripts
keep working and returns `Cell.tetragonal(a, c)`, which is what to write
instead, and Section 17 says so.

Seven names are re-exported from `xrdkit/__init__.py` without being in the
`__all__` of the module that defines them: `accepted_stages`,
`failure_markdown`, `log_tail`, `stage_status`, `stage_status_table`,
`stage_statuses` and `summary_markdown`, all of `gsas2`. They are public by
every other measure, Sections 28.5 and 29.5 use three of them, and they are
listed here under `gsas2`.

| Name | Module | Documented in |
| --- | --- | --- |
| `accepted_stages` | `gsas2` | Section 28.5 |
| `annotate_hkl` | `plotting` | Section 15.4 |
| `apply_style` | `plotting` | Section 15.1 |
| `assess_scan` | `quality` | Section 13.1 |
| `ATOMIC_MASSES` | `density` | Section 19.1 |
| `attribute_unexplained` | `phases` | Section 24.5 |
| `bond_lengths` | `structure` | Section 26.1 |
| `bond_limits` | `library` | Section 25.3 |
| `BreadthModelFit` | `broadening` | Section 21.7 |
| `BreadthModels` | `broadening` | Section 21.7 |
| `BroadeningCorrection` | `broadening` | Section 21.6 |
| `build_parser` | `cli` | Section 30.1 |
| `build_refine_job` | `gsas2` | Section 28.4 |
| `Caglioti` | `broadening` | Section 21.4 |
| `Candidate` | `phases` | Section 24.5 |
| `CandidateMatch` | `phases` | Section 24.4 |
| `Cell` | `cell` | Section 17.1 |
| `cell_contents` | `structure` | Section 26.3 |
| `CELL_PARAMETERS` | `library` | Section 25.3 |
| `cell_volume` | `density` | Section 19.2 |
| `CellFit` | `indexing` | Section 16.4 |
| `check_composition` | `config` | Section 30.2 |
| `CIF_INDEX_COLUMNS` | `phases` | Section 24.2 |
| `cod_fetch` | `phases` | Section 24.2 |
| `cod_search` | `phases` | Section 24.1 |
| `COD_URL` | `phases` | Section 24.5 |
| `CodRecord` | `phases` | Section 24.1 |
| `CommandError` | `cli` | Section 30.1 |
| `component_size_strain` | `sizestrain` | Section 23.3 |
| `ComponentSizeStrain` | `sizestrain` | Section 23.3 |
| `composition_edits` | `structure` | Section 26.4 |
| `ConfigError` | `config` | Section 30.2 |
| `coordinates_stages` | `pipeline` | Section 29.2 |
| `correct_broadening` | `broadening` | Section 21.6 |
| `CRITERIA` | `quality` | Section 13.1 |
| `CRYSTAL_SYSTEMS` | `library` | Section 25.3 |
| `CYCLES` | `pipeline` | Section 29 |
| `DEFAULT_ANIONS` | `library` | Section 25.3 |
| `DEFAULT_BOND_LIMITS` | `library` | Section 25.3 |
| `Distance` | `structure` | Section 26.1 |
| `doublet_gaps` | `broadening` | Section 21.2 |
| `estimate_zero_offset` | `indexing` | Section 16.6 |
| `exchange_edits` | `pipeline` | Section 29.4 |
| `exclude_kalpha2` | `peaks` | Section 14.4 |
| `ExplainedPeak` | `phases` | Section 24.4 |
| `failure_markdown` | `gsas2` | Section 28.5 |
| `fetch_candidates` | `phases` | Section 24.2 |
| `find_gsas2` | `gsas2` | Section 28.2 |
| `find_peaks` | `peaks` | Section 14.2 |
| `find_project` | `project` | Section 27.1 |
| `fit_breadth_models` | `broadening` | Section 21.7 |
| `fit_caglioti` | `broadening` | Section 21.4 |
| `fit_instrument_widths` | `instrument` | Section 22.1 |
| `fit_profile` | `broadening` | Section 21.3 |
| `fixed_atoms_markdown` | `writeup` | Section 30.4 |
| `fixed_atoms_stages` | `pipeline` | Section 29.2 |
| `flag_kalpha2` | `peaks` | Section 14.3 |
| `format_report` | `quality` | Section 13.1 |
| `FORMS` | `project` | Section 27.2 |
| `formula_mass` | `density` | Section 19.1 |
| `generate_reflections` | `indexing` | Section 16.2 |
| `gsas2_fwhm` | `gsas2` | Section 28.3 |
| `GSAS2_HOME_VARIABLE` | `gsas2` | Section 28.2 |
| `GSAS2_PYTHON_VARIABLE` | `gsas2` | Section 28.2 |
| `Gsas2Error` | `gsas2` | Section 28.2 |
| `Gsas2Install` | `gsas2` | Section 28.2 |
| `height_spread_breadth` | `broadening` | Section 21.7 |
| `HELD_INSTRUMENT` | `pipeline` | Section 29.4 |
| `holohedry` | `symmetry` | Section 20.4 |
| `host_elements` | `config` | Section 30.2 |
| `index_and_refine` | `indexing` | Section 16.5 |
| `index_peaks` | `indexing` | Section 16.3 |
| `indexed_to_csv` | `indexing` | Section 16.7 |
| `IndexedPeak` | `indexing` | Section 16.1 |
| `indexing_summary` | `indexing` | Section 16.7 |
| `Inputs` | `pipeline` | Section 29.4 |
| `Instrument` | `project` | Section 27.2 |
| `InstrumentRefinement` | `instrument` | Section 22.2 |
| `integral_breadth` | `broadening` | Section 21.5 |
| `interatomic_distances` | `structure` | Section 26.1 |
| `is_absent` | `symmetry` | Section 20.3 |
| `kalpha2_position` | `broadening` | Section 21.2 |
| `kalpha2_wavelength` | `instrument` | Section 22.3 |
| `lattice_fit_to_dict` | `lattice` | Section 18.2 |
| `LatticeFit` | `lattice` | Section 18.1 |
| `laue_group` | `symmetry` | Section 20.4 |
| `laue_orbit` | `symmetry` | Section 20.4 |
| `LE_BAIL_CYCLES` | `pipeline` | Section 29 |
| `lebail_markdown` | `writeup` | Section 30.4 |
| `lebail_stages` | `pipeline` | Section 29.2 |
| `list_entries` | `library` | Section 25.2 |
| `load_config` | `config` | Section 30.2 |
| `load_entry` | `library` | Section 25.2 |
| `load_project` | `project` | Section 27.1 |
| `load_project_text` | `project` | Section 27.1 |
| `log_tail` | `gsas2` | Section 28.5 |
| `main` | `cli` | Section 30.1 |
| `mark_peaks` | `plotting` | Section 15.5 |
| `match_candidate` | `phases` | Section 24.4 |
| `metric_tensor` | `structure` | Section 26.1 |
| `MissingPhasesExtra` | `phases` | Section 24.5 |
| `mode_paths` | `pipeline` | Section 29.4 |
| `MODES` | `pipeline` | Section 29 |
| `multiplicity` | `symmetry` | Section 20.4 |
| `observed_peaks` | `phases` | Section 24.5 |
| `occupancy_stages` | `pipeline` | Section 29.2 |
| `Options` | `pipeline` | Section 29.4 |
| `Outcome` | `pipeline` | Section 29.4 |
| `parse_formula` | `density` | Section 19.1 |
| `parse_xyz` | `symmetry` | Section 20.1 |
| `PASS_TOLERANCE` | `pipeline` | Section 29 |
| `Peak` | `peaks` | Section 14.2 |
| `peaks_to_csv` | `peaks` | Section 14.4 |
| `PhaseInput` | `pipeline` | Section 29.4 |
| `PHASES_PAUSE_S` | `phases` | Section 24.2 |
| `PHASES_TOLERANCE` | `phases` | Section 24.5 |
| `PHASES_WINDOW` | `phases` | Section 24.5 |
| `PhaseStart` | `pipeline` | Section 29.4 |
| `PipelineError` | `pipeline` | Section 29.2 |
| `plot_caglioti` | `plotting` | Section 15.6 |
| `plot_pattern` | `plotting` | Section 15.2 |
| `plot_rietveld` | `plotting` | Section 15.6 |
| `plot_stacked` | `plotting` | Section 15.3 |
| `ProfileFit` | `broadening` | Section 21.3 |
| `Project` | `project` | Section 27.1 |
| `PROJECT_FILE` | `project` | Section 27.1 |
| `project_template` | `project` | Section 27.3 |
| `pseudo_voigt` | `broadening` | Section 21.1 |
| `pseudo_voigt_components` | `broadening` | Section 21.5 |
| `pseudo_voigt_from_components` | `broadening` | Section 21.5 |
| `rank_candidates` | `phases` | Section 24.5 |
| `read_library_atoms` | `config` | Section 30.2 |
| `read_scan` | `io` | Section 12.2 |
| `read_sites` | `config` | Section 30.2 |
| `read_xrdml` | `io` | Section 12.2 |
| `read_xy` | `io` | Section 12.2 |
| `Refine` | `project` | Section 27.2 |
| `refine_cell` | `indexing` | Section 16.4 |
| `refine_instrument` | `instrument` | Section 22.2 |
| `refine_lattice` | `lattice` | Section 18.1 |
| `refine_settings` | `project` | Section 27.3 |
| `REFINED_KEYS` | `instrument` | Section 22.2 |
| `Reflection` | `indexing` | Section 16.1 |
| `relative` | `writeup` | Section 30.4 |
| `relative_density` | `density` | Section 19.3 |
| `representative` | `symmetry` | Section 20.4 |
| `require_phases_extra` | `phases` | Section 24.5 |
| `resolution_limit` | `sizestrain` | Section 23.4 |
| `ResolutionLimit` | `sizestrain` | Section 23.4 |
| `resolve_inputs` | `pipeline` | Section 29.4 |
| `resolved_cell` | `project` | Section 27.3 |
| `resolved_z` | `project` | Section 27.3 |
| `results_dir` | `project` | Section 27.3 |
| `run_job` | `gsas2` | Section 28.5 |
| `run_mode` | `pipeline` | Section 29.4 |
| `run_sequence` | `pipeline` | Section 29.4 |
| `Sample` | `project` | Section 27.2 |
| `sample_settings` | `config` | Section 30.2 |
| `save_figure` | `plotting` | Section 15.7 |
| `ScanQuality` | `quality` | Section 13.1 |
| `scherrer_size` | `sizestrain` | Section 23.1 |
| `scherrer_size_integral` | `sizestrain` | Section 23.1 |
| `simulate_pattern` | `phases` | Section 24.3 |
| `SimulatedReflection` | `phases` | Section 24.3 |
| `Site` | `library` | Section 25.3 |
| `SITE_KIND` | `config` | Section 30.2 |
| `site_setup` | `structure` | Section 26.2 |
| `space_group_operations` | `symmetry` | Section 20.2 |
| `split_pseudo_voigt` | `broadening` | Section 21.1 |
| `stage_status` | `gsas2` | Section 28.5 |
| `stage_status_table` | `gsas2` | Section 28.5 |
| `stage_statuses` | `gsas2` | Section 28.5 |
| `standard_stages` | `gsas2` | Section 28.4 |
| `START_BROADENING` | `pipeline` | Section 29 |
| `start_from_result` | `pipeline` | Section 29.4 |
| `StartPoint` | `pipeline` | Section 29.4 |
| `structure_edits` | `gsas2` | Section 28.5 |
| `structure_markdown` | `writeup` | Section 30.4 |
| `StructureEntry` | `library` | Section 25.3 |
| `StructureSpec` | `project` | Section 27.2 |
| `summary_markdown` | `gsas2` | Section 28.5 |
| `SUPPORTED_SPACE_GROUPS` | `symmetry` | Section 20.2 |
| `TetragonalCell` | `indexing` | Section 17, deprecated; use `Cell.tetragonal` |
| `theoretical_density` | `density` | Section 19.3 |
| `toml_string` | `project` | Section 27.3 |
| `TTB_CELL` | `indexing` | Section 16 |
| `UnexplainedPeak` | `phases` | Section 24.5 |
| `validate_config` | `config` | Section 30.2 |
| `Verdict` | `quality` | Section 13.1 |
| `__version__` | `xrdkit` | Section 31 |
| `WIDTH_WINDOW` | `instrument` | Section 22.1 |
| `WidthFit` | `instrument` | Section 22.1 |
| `williamson_hall` | `sizestrain` | Section 23.2 |
| `WilliamsonHall` | `sizestrain` | Section 23.2 |
| `with_esd` | `writeup` | Section 30.4 |
| `WORKFLOWS` | `quality` | Section 13.1 |
| `write_cif_index` | `phases` | Section 24.2 |
| `write_instprm` | `gsas2` | Section 28.3 |
| `wyckoff_multiplicity` | `config` | Section 30.2 |
| `XRDScan` | `io` | Section 12.1 |
| `ZeroSearch` | `indexing` | Section 16.6 |
