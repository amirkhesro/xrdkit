# Getting started with xrdkit

This guide sets up a laptop that has nothing installed on it, so that the
workflows in [USER_GUIDE.md](USER_GUIDE.md) can be run. It assumes no previous
experience of Python and no previous experience of the command line. Everything
below is either a thing to click or a line to type.

Nothing here asks you to write any code. xrdkit is a set of commands typed into
a terminal, and what is installed is Python, the package itself, and GSAS-II
for the three commands that need it. There is no editor to set up and no script
to keep.

Part A covers Windows 11 and Part B covers macOS. Do one of them, not both.
Part C is the same on either system, and where a command differs the Windows
form and the macOS form are given side by side.

Allow about half an hour for the whole of Parts A to C, and longer if GSAS-II
is needed, which is Step 6 and is optional.

1. [Part A: Windows 11](#part-a-windows-11)
2. [Part B: macOS](#part-b-macos)
3. [Part C: on either system](#part-c-on-either-system)

## Part A: Windows 11

### Step A1. Install Python 3.11 or later

Open a web browser and go to <https://www.python.org/downloads/windows/>. Under
the most recent stable release, download the file described as the Windows
installer (64-bit). Any release numbered 3.11 or later is suitable.

Run the file that was downloaded. The first screen of the installer carries two
buttons and two tick boxes at the bottom.

1. Tick the box labelled Add python.exe to PATH. This is the important one. It
   is what allows Python to be started by name from a terminal window, and
   without it none of the commands below will be found.
2. Click Install Now.
3. Wait for the progress bar to finish, then click Close. If a final screen
   offers to disable the path length limit, accept it.

Now confirm the installation. Press the Windows key, type PowerShell, and open
Windows PowerShell. At the prompt type the following and press Enter.

```
py --version
```

The reply should name a version of 3.11 or later, for example:

```
Python 3.13.14
```

If instead the reply is that the term py is not recognised, the tick box in
point 1 was missed. Run the installer again, choose Modify, and make sure the
option to add Python to the environment variables is selected.

### Step A2. Open the terminal and move to a folder

The terminal is the window commands are typed into. On Windows it is
PowerShell. Press the Windows key, type PowerShell, and open Windows
PowerShell. A window opens with a prompt that looks like this, where `<name>`
is your Windows user name.

```
PS C:\Users\<name>>
```

The part before the angle bracket is the folder the terminal is currently
working in, which matters because every file name typed into a command is read
relative to it. Three commands manage this.

```
pwd                        show the folder currently worked in
ls                         list what is in it
cd C:\Users\<name>\xrd     change to another folder
```

`cd ..` moves one folder up. A folder name containing a space has to be
surrounded by double quotes, as in `cd "C:\Users\<name>\My Work"`. Typing part
of a name and pressing the Tab key completes it.

To close the terminal, type `exit` or close the window.

### Step A3. Create the project folder

Everything to do with one project lives in one folder, and every command in the
user guide is run from it. Create it under your user folder. In PowerShell,
type the following lines one at a time, pressing Enter after each.

```
cd C:\Users\<name>
mkdir xrd
cd xrd
mkdir data\standards, figures
```

Replace `<name>` with your own Windows user name, which is the name that
already appeared in the prompt in Step A2. The project folder can be called
something other than `xrd` and can live somewhere other than the user folder,
but it should not be inside OneDrive, because the synchronisation can move
files while a command is writing them.

Those two are the folders the commands use without making them first:
`data\standards` for the scan of a standard and the instrument parameter file
made from it, and `figures` for the png and pdf files the plotting commands
write. The rest of the layout is made in Step 5 by `xrdkit init`.

### Step A4. Install xrdkit

With the terminal still open, type:

```
py -m pip install xrdkit
```

This downloads xrdkit and the three libraries it needs, which are numpy, scipy
and matplotlib, and takes a minute or two. The last line printed should begin
with the word Successfully. pip installs xrdkit into Python itself and not into
the folder the terminal happens to be in, so it does not matter which folder
this is run from, and it is done once for the whole computer rather than once
for each project.

One command needs more than that. `xrdkit phases`, which identifies the phases
of a pattern against the Crystallography Open Database, needs the library
pymatgen, which comes with the optional extra named phases. It is a large
download; install it now if the phases are to be identified, and skip it
otherwise. That command also needs the network while it runs, because the
database is searched rather than shipped. To install the extra, type:

```
py -m pip install "xrdkit[phases]"
```

Confirm the installation with one line, which imports xrdkit and prints its
version.

```
py -c "import xrdkit; print(xrdkit.__version__)"
```

The reply should be a version number:

```
0.2.0
```

Part A is finished. Continue at Part C.

## Part B: macOS

### Step B1. Install Python 3.11 or later

macOS may already carry a Python, but it is the one the system itself uses and
it should be left alone. Install a Python of your own instead.

Open a web browser and go to <https://www.python.org/downloads/macos/>. Under
the most recent stable release, download the file described as the macOS
64-bit universal2 installer. Any release numbered 3.11 or later is suitable.

Open the file that was downloaded, which is a package ending in `.pkg`, and
work through the installer.

1. Click Continue on the introduction, the read me and the licence screens, and
   Agree when asked to accept the licence.
2. Click Install. Enter your Mac password when it is asked for, and wait for
   the installation to finish.
3. When it finishes, a Finder window opens showing the installed folder, which
   is Python 3.x inside Applications. Double click the file in it named Install
   Certificates.command, let the small terminal window that appears finish, and
   close it. This sets up the certificates Python uses for secure downloads,
   which `xrdkit phases` needs in order to reach the Crystallography Open
   Database.

Now confirm the installation. Open Terminal, which Step B2 explains how to
find, then type the following and press Return.

```
python3 --version
```

The reply should name a version of 3.11 or later, for example:

```
Python 3.13.7
```

If the version named is older than 3.11, the terminal is still finding the
system Python. Close the Terminal window, open a new one, and try again; the
installer puts its own Python ahead of the system one only in terminal windows
opened after it ran.

### Step B2. Open the terminal and move to a folder

The terminal is the window commands are typed into. On macOS it is Terminal.
Press Command and the space bar together, type Terminal, and press Return.
Terminal can also be found in the Utilities folder inside Applications.

A window opens with a prompt ending in a percent sign. The name just before it
is the folder the terminal is currently working in, which matters because every
file name typed into a command is read relative to it. Three commands manage
this.

```
pwd                   show the folder currently worked in
ls                    list what is in it
cd /Users/<name>/xrd  change to another folder
```

The tilde character is shorthand for your home folder, so `cd ~/xrd` and
`cd /Users/<name>/xrd` mean the same thing. `cd ..` moves one folder up. A
folder name containing a space has to be surrounded by double quotes, as in
`cd "/Users/<name>/My Work"`, and the quotation marks have to go around the
full path rather than around a tilde, which is not read as the home folder
inside them. Typing part of a name and pressing the Tab key completes it.
Dragging a folder from the Finder onto the Terminal window types its full path.

To close the terminal, type `exit` or close the window.

### Step B3. Create the project folder

Everything to do with one project lives in one folder, and every command in the
user guide is run from it. Create it in your home folder. In Terminal, type the
following lines one at a time, pressing Return after each.

```
cd ~
mkdir xrd
cd xrd
mkdir -p data/standards figures
```

The project folder can be called something other than `xrd` and can live
somewhere other than the home folder, but it should not be inside iCloud Drive,
because the synchronisation can move files while a command is writing them.

Those two are the folders the commands use without making them first:
`data/standards` for the scan of a standard and the instrument parameter file
made from it, and `figures` for the png and pdf files the plotting commands
write. The rest of the layout is made in Step 5 by `xrdkit init`.

### Step B4. Install xrdkit

With the terminal still open, type:

```
python3 -m pip install xrdkit
```

This downloads xrdkit and the three libraries it needs, which are numpy, scipy
and matplotlib, and takes a minute or two. The last line printed should begin
with the word Successfully. pip installs xrdkit into Python itself and not into
the folder the terminal happens to be in, so it does not matter which folder
this is run from, and it is done once for the whole computer rather than once
for each project.

One command needs more than that. `xrdkit phases`, which identifies the phases
of a pattern against the Crystallography Open Database, needs the library
pymatgen, which comes with the optional extra named phases. It is a large
download; install it now if the phases are to be identified, and skip it
otherwise. That command also needs the network while it runs, because the
database is searched rather than shipped. To install the extra, type:

```
python3 -m pip install "xrdkit[phases]"
```

The quotation marks matter here, because without them the shell tries to read
the square brackets as a file name pattern.

Confirm the installation with one line, which imports xrdkit and prints its
version.

```
python3 -c "import xrdkit; print(xrdkit.__version__)"
```

The reply should be a version number:

```
0.2.0
```

Part B is finished. Continue at Part C.

## Part C: on either system

From here the two systems differ in only two ways. Python is started with `py`
on Windows and `python3` on macOS, and the terminal is PowerShell on Windows
and Terminal on macOS. The `xrdkit` command itself is typed the same way on
both. The steps of the two parts above are numbered A1 to A4 and B1 to B4, and
a reference below to Step 3 or Step 4 means the one in whichever part was
followed.

### Step 5. Start the project

Everything from here is typed in the terminal, in the project folder made in
Step 3, which is reached with `cd C:\Users\<name>\xrd` on Windows or `cd ~/xrd`
on macOS. One command starts a project.

```
xrdkit init --name xrd
```

It prints the one file it wrote, which on Windows is:

```
C:\Users\<name>\xrd\xrdkit.toml
```

`--name` is the project name recorded in the file and defaults to the name of
the folder, so it can be left out. An existing `xrdkit.toml` is never
overwritten, so running the command twice is safe.

Besides the file, `init` makes three folders if they are missing, and with the
two made by hand in Step 3 the layout is this.

```
xrd\
    cifs\              one CIF per phase, with a note of its source
    data\raw\          scans as the instrument wrote them, never edited
    data\standards\    the standard scan and the instrument parameter file
    figures\           the png and pdf files the plotting commands write
    results\           peak lists, refinement projects and exports
    xrdkit.toml        the project file
```

Check it with `ls`. Section 2.1 of the user guide is the same command in full,
and Section 2.2 goes through `xrdkit.toml` table by table: the instrument, the
structures and the samples are filled in there, and every command after this
one can then take a sample key in place of a file name.

### Step 6. GSAS-II, if it is needed

GSAS-II is needed by three commands and by nothing else: `xrdkit instrument`,
which makes the instrument parameter file; `xrdkit lebail`, which extracts a
cell by Le Bail fitting; and `xrdkit rietveld`, which is the full Rietveld
refinement. Sections 3 and 8 of the user guide are those commands. Everything
else in the kit works without it, so skip this step until one of the three is
actually wanted.

GSAS-II is not a Python library and cannot be installed with pip. It brings its
own Python with it, and xrdkit runs its refinement jobs under that Python
rather than importing it. Install it from its own home page,
<https://advancedphotonsource.github.io/GSAS-II-tutorials>, following the
installation instructions there for your system. Nothing about the installation
is specific to xrdkit, and xrdkit writes nothing into it.

The one thing xrdkit has to know is where GSAS-II ended up. It looks in three
places, in this order: the environment variable `XRDKIT_GSAS2_PYTHON`, for the
GSAS-II Python itself, and `XRDKIT_GSAS2_HOME`, for the folder that contains
the `GSASII` package; failing those, a folder named `gsas2main` in your home
folder, which is `C:\Users\<name>\gsas2main` on Windows and
`/Users/<name>/gsas2main` on macOS. In that default installation the GSAS-II
Python sits at the top of `gsas2main` and the package folder, named `GSAS-II`,
sits below it. Install GSAS-II there and there is nothing to set.

If it is installed anywhere else, set the two variables. On Windows, for the
current terminal window only:

```
$env:XRDKIT_GSAS2_PYTHON = "C:\path\to\gsas2\python.exe"
$env:XRDKIT_GSAS2_HOME = "C:\path\to\gsas2\GSAS-II"
```

and for every window from now on, after which a new terminal window has to be
opened before the setting takes effect:

```
setx XRDKIT_GSAS2_PYTHON "C:\path\to\gsas2\python.exe"
setx XRDKIT_GSAS2_HOME "C:\path\to\gsas2\GSAS-II"
```

On macOS, for the current terminal window only:

```
export XRDKIT_GSAS2_PYTHON=/path/to/gsas2/bin/python
export XRDKIT_GSAS2_HOME=/path/to/gsas2/GSAS-II
```

and for every window from now on, by adding those two lines to the file
`~/.zshrc`, which the shell reads when it starts.

There is nothing separate to run to check this. Each of the three commands
looks GSAS-II up before it does anything else, and one that cannot find it
stops with a single line naming what was missing and both variables, having
written nothing.

### Step 7. Copying a scan from the diffractometer

The diffractometer writes the scan to a file. The reader accepts `.xrdml` and
the two and three column text patterns `.xy` and `.xye`: most diffractometer
software exports one of these directly, and for any other format a converter
that writes `.xy` will do. Copy the file, by USB stick or over the network
share the instrument writes to, into the `data/raw` folder of the project. A
`.xy` or `.xye` carries no wavelength, so it is given on the command line with
`--wavelength`, or by the instrument table of the project file when the scan is
named as a sample. Section 2.3 of the user guide is the formats and that
wavelength rule in full.

Copy it rather than moving it, and do not open it in anything that might write
it back. The raw file is the record of what the instrument measured, the kit
never writes into it, and every result can be produced again from it. Section 9
of the user guide lists what each command writes and where.

Two points of practice save trouble later. Give the file a name you will
recognise, such as the sample identifier, and note the original instrument file
name alongside it, in a notebook or in a table of your own. Keep the standard
scans apart from the samples: a lanthanum hexaboride or silicon scan measured
on the same instrument with the same optics belongs in `data/standards`, not in
`data/raw`, because it describes the diffractometer rather than a sample.
`xrdkit instrument` is what turns it into the instrument parameter file that
`xrdkit lebail` and `xrdkit rietveld` both need.

Check that the file arrived, with `ls data/raw` on either system.

### Step 8. Check that everything works

Two commands exercise the reader, the quality criteria, the plotting style and
the figure writer in one go. Put a scan of your own in `data/raw` first, as
Step 7 describes, and use its name in place of `powder_a.xrdml` below.

The first reports what the scan is and whether it is good enough for each
workflow.

```
xrdkit check data/raw/powder_a.xrdml
```

It prints the numbers it measured, then a verdict per workflow.

```
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

Your numbers will be your own, and so will the verdicts. Not suitable is not a
failure of the installation: it is the command saying what that scan will and
will not support, and Section 4.2 of the user guide is the table of criteria
behind every line of it.

The second draws the pattern.

```
xrdkit plot data/raw/powder_a.xrdml --scale sqrt
```

It prints the files it wrote, which on Windows are:

```
results\peaks_powder_a.csv
figures\pattern_powder_a.png
figures\pattern_powder_a.pdf
```

On macOS the same paths are written with forward slashes. Open
`figures/pattern_powder_a.png` and you should see the diffraction pattern, on a
square root intensity scale, with two theta along the bottom.

If that worked, the installation is sound and everything in the user guide will
run. Three things go wrong at this point more often than anything else. A reply
that the term xrdkit is not recognised means Step 4 of Part A or Part B did not
finish, or finished into a different Python from the one the terminal finds. A
message that the file does not exist means the terminal is not in the project
folder, or the scan is not in `data/raw` under the name that was typed. An
error from the reader means the file is not one of the formats it accepts;
Section 2.3 of the user guide lists them and says how to supply the wavelength
a text pattern leaves out.

### Step 9. VESTA, if structures are to be looked at

VESTA draws crystal structures. Nothing in the kit calls it and no command
needs it, but it is the quickest way to see what a file actually holds. A CIF
downloaded from the Crystallography Open Database opens in it in a second,
which is how to tell at a glance that the entry is the phase it was taken for
rather than something else with a similar name. It also opens the structure a
GSAS-II refinement exports, so a refined model can be looked at beside the one
it started from, and it draws the structure figures that go in a paper.

VESTA is free for academic use and there is a version for Windows and a version
for macOS. The download page is
<https://jp-minerals.org/vesta/en/download.html>.

### Where to go next

Steps 5, 7 and 8 come first, and none of them can be skipped. The project file
has to exist, the scan has to be in `data/raw` under a name you know, and the
two commands of Step 8 have to have run, before the user guide can be followed
at all.

Then read [USER_GUIDE.md](USER_GUIDE.md). Part I, Sections 1 to 10, is the
whole kit as commands, in the order a project runs them: `init` and the project
file in Section 2, `instrument` and `add-sample` in Section 3, `check` in
Section 4, `plot` and `stack` in Section 5, `phases` in Section 6, `lattice`
and `density` in Section 7, and `lebail` and `rietveld` in Section 8. Section 9
is what each command writes and where, and Section 10 is the known limitations,
which is worth reading before relying on a result. Part II, Sections 11 to 31,
is the library reference, for anyone who would rather call the same functions
from code of their own than type the commands.

#### Doing it again for the next sample

None of it has to be set up twice. Once the project file names the instrument
and the structures, the next dataset is one copy into `data/raw` and one
`xrdkit add-sample`, after which every command takes the new sample key and
writes its results under that key. Nothing is edited but the project file, and
nothing is repeated but the commands themselves.
