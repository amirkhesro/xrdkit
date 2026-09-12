# Getting started with xrdkit

This guide sets up a laptop that has nothing installed on it, so that the
workflows in [USER_GUIDE.md](USER_GUIDE.md) can be run. It assumes no previous
experience of Python and no previous experience of the command line. Everything
below is either a thing to click or a line to type.

Part A covers Windows 11 and Part B covers macOS. Do one of them, not both.
Part C is the same on either system, and where a command differs the Windows
form and the macOS form are given side by side.

Allow about half an hour for the whole of Parts A to C, and longer if GSAS-II
is needed, which is Step 7 and is optional.

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

Everything to do with one project lives in one folder, and the workflows in the
user guide are written on the assumption that they are run from it. Create it
under your user folder. In PowerShell, type the following lines one at a time,
pressing Enter after each.

```
cd C:\Users\<name>
mkdir xrd
cd xrd
mkdir data\raw, data\standards, cifs, config, results, figures
```

Replace `<name>` with your own Windows user name, which is the name that
already appeared in the prompt in Step A2. The project folder can be called
something other than `xrd` and can live somewhere other than the user folder,
but it should not be inside OneDrive, because the synchronisation can move
files while a script is writing them.

The result is the layout that Section 3 of the user guide describes.

```
C:\Users\<name>\xrd\
    data\raw\          scans as the instrument wrote them, never edited
    data\standards\    LaB6 and silicon scans, and the instrument parameter files
    cifs\              one CIF per phase, with a note of its source
    config\            the TOML settings the Rietveld workflow reads
    results\           peak lists, indexing tables, refinement projects and exports
    figures\           the png and pdf files the plotting routines write
```

Check it with `ls`, which should list the six folders just made.

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

Phase matching against the Crystallography Open Database needs one further
library, pymatgen, which comes with the optional extra named phases. It is a
large download and is worth installing only if the search and match routines in
`xrdkit.phases` are going to be used. If they are, type:

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
0.1.0
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
   which the phase matching routines need in order to reach the Crystallography
   Open Database.

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

Everything to do with one project lives in one folder, and the workflows in the
user guide are written on the assumption that they are run from it. Create it
in your home folder. In Terminal, type the following lines one at a time,
pressing Return after each.

```
cd ~
mkdir xrd
cd xrd
mkdir -p data/raw data/standards cifs config results figures
```

The project folder can be called something other than `xrd` and can live
somewhere other than the home folder, but it should not be inside iCloud Drive,
because the synchronisation can move files while a script is writing them.

The result is the layout that Section 3 of the user guide describes.

```
/Users/<name>/xrd/
    data/raw/          scans as the instrument wrote them, never edited
    data/standards/    LaB6 and silicon scans, and the instrument parameter files
    cifs/              one CIF per phase, with a note of its source
    config/            the TOML settings the Rietveld workflow reads
    results/           peak lists, indexing tables, refinement projects and exports
    figures/           the png and pdf files the plotting routines write
```

Check it with `ls`, which should list the six folders just made.

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

Phase matching against the Crystallography Open Database needs one further
library, pymatgen, which comes with the optional extra named phases. It is a
large download and is worth installing only if the search and match routines in
`xrdkit.phases` are going to be used. If they are, type:

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
0.1.0
```

Part B is finished. Continue at Part C.

## Part C: on either system

From here the two systems differ in only two ways. Python is started with `py`
on Windows and `python3` on macOS, and the terminal is PowerShell on Windows
and Terminal on macOS. Both forms are given wherever it matters. The steps of
the two parts above are numbered A1 to A4 and B1 to B4, and a reference below
to Step 3 or Step 4 means the one in whichever part was followed.

### Step 5. Install an editor

An editor is where the scripts are written. Visual Studio Code is recommended
but is not required, and anything that saves plain text will do: Notepad on
Windows, or TextEdit on macOS with Format set to Make Plain Text. A word
processor such as Word will not do, because it saves formatting along with the
words and Python cannot read the result.

Visual Studio Code is free. Download it from <https://code.visualstudio.com>,
which offers the right version for the system it is visited from.

On Windows, run the installer that was downloaded, accept the licence, and
click Next through to Install. On the screen offering additional tasks, leaving
every box at its default is fine.

On macOS, open the file that was downloaded, which unpacks to an application
named Visual Studio Code, and drag that application into the Applications
folder. Open it from there.

Then add the Python extension, which colours the code, points out mistakes as
they are typed, and offers a button that runs the open file.

1. Open Visual Studio Code.
2. Click the Extensions icon in the bar down the left side, which is the one
   made of four small squares. The keyboard equivalent is Ctrl and Shift and X
   together on Windows, or Command and Shift and X together on macOS.
3. Type Python into the search box.
4. Install the extension named Python that is published by Microsoft, which is
   normally the first result.

To open the project folder in the editor, choose File, then Open Folder on
Windows or Open on macOS, and select the `xrd` folder made in Step 3 of Part A
or Part B.

### Step 6. Running the code blocks in the user guide

Every code block in [USER_GUIDE.md](USER_GUIDE.md) is a piece of a script.
There is no menu and no dialogue to click through: the way to use the kit is to
put the lines of a block into a file, save the file, and run it. The procedure
is the same every time.

1. In the editor, make a new file and paste the block into it.
2. Save it in the project folder made in Step 3, under a name ending in `.py`,
   for example `plot.py`. The name is yours to choose, but it should not be
   `xrdkit.py`, because a file of that name would be found instead of the
   library.
3. In the terminal, change to the project folder, which is
   `cd C:\Users\<name>\xrd` on Windows or `cd ~/xrd` on macOS.
4. Run it, with `py plot.py` on Windows or `python3 plot.py` on macOS.

Four things follow from running a script that way.

Anything the script prints appears in the terminal, underneath the command,
which is where the printed output shown in the user guide comes from.

Anything the script writes, such as a figure or a CSV file, appears in the
project folder, in the sub folder the path names. Figures are not displayed in
a window; they are saved as files, and `save_figure` writes a png and a pdf of
each by default.

Every path in the code blocks, such as `data/raw/10s.xrdml` or
`figures/pattern_10`, is relative to the folder the terminal is working in when
the script is run. Running from anywhere else gives an error saying the file
was not found. If that happens, check with `pwd` that the terminal really is in
the project folder. Forward slashes in a path work on Windows as well as on
macOS, which is why the guide uses them throughout.

A block part way through a workflow usually depends on the blocks before it in
the same section, which is where its `scan` or its `peaks` came from. Put the
blocks of a section into one file, in the order they appear, rather than
running each alone.

### Step 7. GSAS-II, if it is needed

GSAS-II is needed for two things only: Workflow 3, which is Rietveld
refinement, and Route B of Workflow 2, which extracts the lattice parameters
from a Le Bail fit. Workflow 1 and Route A of Workflow 2 do not use it, and
neither does anything else in the kit. Skip this step until one of those two is
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

Either way, check what the kit resolves to before going any further. Save this
as `check_gsas2.py` and run it as Step 6 describes.

```python
from xrdkit import find_gsas2

install = find_gsas2()
print(f"GSAS-II Python: {install.python}")
print(f"GSAS-II home:   {install.home}")
```

It prints the two paths it found, as in the example in Section 5.2 of the user
guide. If instead it raises `FileNotFoundError`, the message names what was
missing and both variables, and nothing in Workflow 3 will run until it is
resolved.

### Step 8. Copying a scan from the diffractometer

The diffractometer writes a file ending in `.xrdml`, which is the only format
the reader accepts. Copy it, by USB stick or over the network share the
instrument writes to, into the `data/raw` folder of the project.

Copy it rather than moving it, and do not open it in anything that might write
it back. The raw file is the record of what the instrument measured, the kit
never writes into it, and every result can be produced again from it. Section 1
of the user guide says what the kit does and does not touch.

Two points of practice save trouble later. Give the file a name you will
recognise, such as the sample identifier, and note the original instrument file
name alongside it, in a notebook or in a table of your own. Keep the standard
scans apart from the samples: a LaB6 or silicon scan measured on the same
instrument with the same optics belongs in `data/standards`, not in `data/raw`,
because it describes the instrument rather than a sample. Workflow 2 and
Workflow 3 both need one.

Check that the file arrived, with `ls data/raw` on either system.

### Step 9. Check that everything works

This plots one scan, which exercises the reader, the plotting style and the
figure writer in one go. Put a scan of your own in `data/raw` first, as Step 8
describes.

Save the following as `check.py` in the project folder, changing
`data/raw/10s.xrdml` to the name of the file that was copied there.

```python
from xrdkit import apply_style, plot_pattern, read_xrdml, save_figure

apply_style()
scan = read_xrdml("data/raw/10s.xrdml")
print(f"{scan.sample_id}: {len(scan.two_theta)} points, wavelength {scan.wavelength} A")

fig, ax = plot_pattern(scan, scale="sqrt")
paths = save_figure(fig, "figures/check")
print("wrote", ", ".join(str(path) for path in paths))
```

Run it with `py check.py` on Windows or `python3 check.py` on macOS, from the
project folder. The output names the sample the file carries, the number of
points in the scan and the wavelength it was measured at, then the files
written. On Windows it looks like this:

```
10s: 4141 points, wavelength 1.540598 A
wrote figures\check.png, figures\check.pdf
```

On macOS the two paths are written with forward slashes instead. Open
`figures/check.png` and you should see the diffraction pattern, on a square
root intensity scale, with two theta along the bottom.

If that worked, the installation is sound and everything in the user guide will
run. Three things go wrong at this point more often than anything else. A
`ModuleNotFoundError` naming xrdkit means Step 4 of Part A or Part B did not
finish, or was run under a different Python from the one running the script. A `FileNotFoundError`
means the terminal is not in the project folder, or the scan is not in
`data/raw` under the name the script uses. An error from the reader means the
file is not an `.xrdml`; Section 3 of the user guide gives what the reader
accepts and Section 7 gives what to do with a scan in another format.

### Step 10. VESTA, if structures are to be looked at

VESTA draws crystal structures. Nothing in the kit calls it and no workflow
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

Steps 8 and 9 come first, and neither can be skipped. Section 2 of the user
guide reads a scan of your own and prints what it is, so the scan has to be in
`data/raw` under a name you know, which is Step 8, and the check in Step 9 has
to have run, before Section 2 can be followed at all.

Read Section 2 of [USER_GUIDE.md](USER_GUIDE.md) first, which says in numbers
what a scan has to be for each of the three workflows, and how to check a scan
against it. Then take Workflow 1, in Section 4, which plots a pattern and
labels its reflections and needs nothing beyond what was installed above.
Workflow 2, in Section 5, refines the lattice parameters and calculates a
theoretical density. Workflow 3, in Section 6, is the full Rietveld refinement,
and needs GSAS-II from Step 7.

#### Reusing the scripts

None of it has to be written twice. Once a workflow has been saved as a script
and has run on one scan, the next dataset needs no new code at all: open the
script, change the name of the scan file and the sample labels at the top of
it, change the figure and results names as well if the first set is worth
keeping, and run the script again. That is the whole of the work for each
further sample, and it is why the guide is written as scripts rather than as
commands typed one at a time.
