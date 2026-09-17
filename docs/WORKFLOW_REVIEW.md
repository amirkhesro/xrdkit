# xrdkit workflow review

Review of xrdkit 0.1.0 and the XRD-Analysis project against the group's routine XRD work: phase purity of calcined powders and sintered pellets, lattice parameters and theoretical density for comparison with Archimedes density, Rietveld refinement of selected compositions, and indexed single and stacked figures for papers and talks. Prepared 15 September 2026 from the source in `xrdkit/src/xrdkit`, the two guides in `xrdkit/docs`, and the scripts, config and results in `XRD-Analysis`.

## 1. Summary

The kit has a sound scientific core and an unusually careful Rietveld driver, but it is a library with no command line, and it was grown around one composition series on one instrument. The two consequences for a user in the group are these. First, every routine task still means copying a script of 30 to 180 lines from the user guide and editing it, and the full refinement pipeline lives in a 3,773 line script in the data repository rather than in the package. Second, nothing outside the GSAS-II driver knows any crystal system but tetragonal, so a cubic, rhombohedral or orthorhombic perovskite cannot be indexed, fitted for lattice parameters, or described in the refinement config at all, and a two phase sample cannot be described either.

The plan in Section 4 turns the kit into a command driven tool (`xrdkit check`, `xrdkit plot`, `xrdkit purity`, `xrdkit lattice`, `xrdkit rietveld`, `xrdkit report`) over a generalised core, with a structure library for perovskites and TTBs, an instrument registry covering the Aeris, Empyrean, X'Pert3, D2 Phaser and STADI P, a licence free COD route beside the PDF-4 route, and a report writer that produces the group's Table 1 directly. The order is chosen so that each step is usable on its own and no command has to be written twice.

## 2. What is already good

These should be kept as they are. The science modules are correct where they apply: the density formula and atomic masses, the Caglioti and instrumental broadening treatment, the size and strain analysis, and the bond length calculation from the metric tensor. The GSAS-II driver is the strongest part of the kit: it runs unattended under GSAS-II's own Python, refines each stage in passes until it settles, records every residual and parameter with its esd, rolls back a stage that produces a negative Uiso or an impossible occupancy, and labels parameters the data do not determine. Very few groups have anything like it. The figure style is right for publication (Arial 9 pt, inward ticks, all four spines, 300 dpi, png and pdf together), the hkl labels sit in one rotated row with collision handling, and the stacked figure matches what the group publishes. The separation of code from data, the read only rule on `data/raw`, the TOML config checked on load, the test suite (14 files) and the two guides are all better than the norm for an academic package.

## 3. Shortcomings

### 3.1 There is no command line, so users still write scripts

`pyproject.toml` has no `[project.scripts]` entry and nothing under `src/xrdkit` uses argparse, click or typer. The user guide is honest about this: "It is a library rather than a program: every routine is imported and called from a script of your own" (USER_GUIDE.md line 33). To do the routine work a new user copies five scripts from the guide, `plot_pattern.py` (74 lines), `stack_patterns.py` (34), `identify_phases.py` (73), `make_instprm.py` (78) and `lattice_density.py` (130), about 320 lines of code before Rietveld, then `refine_rietveld.py` (176 lines) plus a 40 line TOML for Workflow 4. Each script carries three to ten settings lines that must be edited per sample, and the guide has to warn against Replace All because the stem `10` also appears as a threshold in the same file (lines 91 to 98).

The XRD-Analysis project shows where this leads. It holds 14 scripts, including `refine_sample.py` at 3,773 lines, of which roughly 85 to 90 per cent is generic pipeline (job building, stage sequences, fit diagnostics, 1,600 lines of markdown write up) that belongs in the package. Sample constants are baked into scripts (`plot_indexed.py` lines 23 to 35, `screen_cod.py` 44 to 47, `common.py` line 50), the Aeris instprm path is hard coded in `common.py` line 32, and `refine_sample.py` imports from `instrument_gsas2.py` so the two must sit together. The refinement script also needs a row in `data/samples.csv` and a row in `cifs/index.csv` in addition to the `samples.toml` table, so a sample is described in three places.

### 3.2 Everything outside the GSAS-II driver is tetragonal and TTB only

`TetragonalCell(a, c)` is the only cell class (`indexing.py` 78 to 96); `refine_cell`, `refine_lattice` and `density.cell_volume` all assume it (`indexing.py` 520 to 615, `lattice.py` 57 to 70, `density.py` 88 to 105); `config._start_cell` accepts only `{a, c}` (`config.py` 459 to 472). `SUPPORTED_SPACE_GROUPS = ("P4bm",)` (`indexing.py` 52) and P4bm is the default everywhere; P4/mbm, which shares its conditions, is not even listed. `TTB_CELL = (12.45, 3.94)` is exported at package level and is the start cell in `scripts/common.py` line 112.

A cubic perovskite cannot be indexed today even by setting c = a: `generate_reflections` enumerates h >= k >= 0 and l >= 0 separately, so (100) and (001) become two candidates for one peak, `index_and_refine` uses only single candidate peaks, and `refine_cell` then refines c independently of a. With the P4bm default the cubic (100) would be labelled (001) because (100) is forbidden. Rhombohedral, orthorhombic and monoclinic cells have no class, no d spacing formula and no reflection conditions. `structure.metric_tensor` is general but is used only for bond lengths.

On the Rietveld side the driver is general (symmetry and polar axes from GSAS-II's SGData, coordinate freedom from `GetCSxinel`, any six parameter cell), but `config.py` prevents a perovskite from being described: `exchange` is a required key and must name two or more elements on two or more sites of one kind (`config.py` 100 to 111, 327 to 350); `origin` is required even for centrosymmetric groups; every structure must have a site of each kind A, B and O (247 to 249); a sample has exactly one `structure` (526 to 533). Pm-3m, R3c and Pbnm perovskites have one A site and one B site, so no valid `[structures]` table exists for them, and a two phase sample such as R3c plus Pm-3m in BiFeO3-BaTiO3 cannot be expressed. `refine_sample.py` then reads the cell as `[a, a, c, 90, 90, 90]` (lines 503, 2041), reads A0 and A2 as 1/a² and 1/c² (527 to 538, wrong for a hexagonal setting), and uses parameter names that assume one phase and one histogram (`":0:Zero"`, `"0::A0"`, `"0:0:Size;i"`).

### 3.3 One file format and one instrument

`read_xrdml` is the only reader (`io.py` 42 to 115). There is nothing for Bruker `.raw` or `.brml` (D2 Phaser), STOE STADI P, or plain `.xy`/`.xye`, and the guide tells users to build an `XRDScan` by hand (Section 8.1). The reader takes only the first `<scan>` in the file, ignores `<listPositions>` and reads Kα1 but not Kα2, the anode or the goniometer radius, which `common.py` lines 52 to 57 hand codes as 145 mm with a note that the file contains it. There is no registry of instruments: the instprm is a bare file path, `write_instprm` always writes a Cu Kα1 plus Kα2 file (`gsas2.py` 164 to 165) so a Kα1 only configuration such as the STADI P with its Ge monochromator cannot be written without overriding the wavelengths by hand, and the Kα2 flagging in `find_peaks` assumes Cu regardless of the scan's wavelength (`peaks.py` 212 to 233), which would misplace the satellites on a Co tube.

### 3.4 Phase identification has no PDF-4 route and a weak COD match

`cod_search` and `cod_fetch` work and record provenance properly. The match, however, scores a candidate as explained lines minus missing lines (`phases.py` 342 to 344) with no intensity weighting, no normalisation to line count and no cell or zero refinement, so a large cell candidate with hundreds of weak lines explains everything. Simulation needs pymatgen (an optional extra that installs hundreds of megabytes) and `screen_cod.py` had to add a 180 s per CIF timeout and a subprocess pool to make it usable. There is no way to bring an ICDD PDF-4 or HighScore result in: no reader for card numbers, exported stick patterns or reference `.xy` files (the guide says so at lines 777 to 779 and 828 to 829), and no reference tick rows can be drawn under a pattern or a stack, so `match_foreign_peaks.py` draws its own with `vlines`. The peak finder's threshold of two per cent of the strongest peak (`peaks.py` 23, 181) is above the level at which a one to two weight per cent second phase shows, so the routine purity check will not list exactly the peaks it exists to find, and there is no noise based or absolute counts option.

### 3.5 Density and lattice parameters stop short of what the group reports

`theoretical_density` takes Z as an explicit argument and nothing derives it from the CIF or space group, although `config.py` already carries `formula_units`. `formula_mass` takes only a mapping, not a formula string. There is no relative density function, so the Archimedes over theoretical comparison, which is the whole purpose of the exercise, is eight lines the user writes (guide 1530 to 1537) and `scripts/density.py` does not compute it at all. `ATOMIC_MASSES` has 33 elements and lacks B (so the LaB6 standard's density cannot be computed), C and H (carbonate and hydroxide precursors), and common dopants such as Pr, Eu, Dy, Er, Yb, Sc, V, Ga, In and Sb. `refine_lattice(fit_zero=False)` holds the zero at `start_zero`, not at zero, contrary to its own docstring (`lattice.py` 124, 154, 168), and `scripts/refine_lattice.py` therefore compares models that are not the ones its comments describe.

### 3.6 Rietveld results are recorded but not reported

The driver records Rwp, Rp, GOF, chi² and reduced chi² for every stage (`gsas2_driver.py` 2010 to 2035) and per phase cell, volume, space group, phase fraction and weight fraction in `final` (2047 to 2072), so every column of the group's Table 1 exists in the result JSON. Nothing tabulates it: `summary_markdown` prints Rwp and GOF only, `refine_sample.py` prints a and c only, and there is no function that writes a per composition table of space group, a, b, c, V, phase fraction, Rwp, Rp and chi². `plot_rietveld` already draws observed circles, calculated line, difference and one tick row per phase, which matches the group's figure, but its caption takes Rwp and GOF from the last stage without an error entry, which is the wrong stage after a rollback (`plotting.py` 821 to 826, listed by the author), shows neither Rp, chi² nor fractions, cannot draw the difference in a separate panel, and cannot stack several compositions' fits.

### 3.7 Scientific and robustness issues

These are ordered by how much they affect routine numbers.

Peak positions for lattice parameters are biased by Kα2. `find_peaks` takes a three point parabola vertex on the raw doublet (`peaks.py` 102 to 123), which below about 50° sits between Kα1 and the centroid while resolved high angle peaks give the true Kα1 position; the mixed bias looks like a cell or zero error in `refine_lattice`. The proper doublet fit exists in `broadening.fit_profile` but the indexing path never calls it.

Displacement is not the default for pellets. The zero is refined per sample and displacement never (`refine_sample.py` 435, 1328, 2299; `common.py` line 43); in Bragg-Brentano geometry the sample dependent aberration is displacement and the zero belongs to the instrument. The sample zero of minus 0.0356° against the LaB6 zero of minus 0.0235° is a displacement signature.

Microstrain is refined on every Le Bail run although the author found it unstable (`samples.toml` 99 to 115), with `lgmix` fixed at 1.0 so the Gaussian and Lorentzian parts are never used to separate size from strain. No preferred orientation, absorption, surface roughness or transparency term exists (`SAMPLE_KEYS` in the driver line 288 has scale and displacement only), which matters most for the polished pellets; texture is only tested after the fact with a tetragonal c axis March-Dollase.

`find_undetermined` judges esd against the shift from where the current job started, so a mode that starts from converged values flags every Uiso (summary.md lines 12 to 14 show all fifteen). The occupancy constraint holds each element's content but not each site's total, which is why the Sr and Ba exchange is 100 per cent correlated and Ba1 ran to minus 0.22 before rollback. A stage with scale and every phase fraction free has no guard against the resulting singularity.

Smaller items: every path to `build_refine_job` must be absolute (author listed); `prefix.with_suffix(".gpx")` turns a stem such as `x0.10` into `x0.gpx`; `build_refine_job` never checks that the scan, instprm or CIF exist, so a missing file surfaces as a GSAS-II traceback; `write_instprm` omits `Gonio. radius` and `Diff-type`; GSAS-II is looked for only in `~/gsas2main` or the two environment variables; `_format_axes` forces 10° major ticks so a 20 to 35° detail figure gets two ticks; `plot_stacked` uses one colour for all traces and has no legend or reference ticks.

## 4. Steps to improve it

The steps are ordered so that each leaves the kit usable and so that the command line is built once, on the generalised core, rather than written over the tetragonal code and rewritten later. Steps 1 to 4 are the ones that change daily work; 5 to 8 finish the job. Each step is a natural boundary for one Claude Code session or a short series.

### Step 1. A thin command line over what exists

Add `[project.scripts] xrdkit = "xrdkit.cli:main"` and a `cli` module (argparse is enough and adds no dependency; typer is nicer if a dependency is acceptable). Start with the commands that need no new science: `xrdkit check scan.xrdml` (the Section 2 quality report against the table in the guide, with a verdict per workflow), `xrdkit plot scan.xrdml --label "x = 0.10"`, `xrdkit stack a.xrdml b.xrdml --labels ...`, and `xrdkit density --formula ... --z 5 --cell a c --archimedes 5.21`. Every command writes its outputs under `results/` by sample name and prints the paths, exactly as the guide's scripts do now, and takes `--json` so that Claude Code and other tools can consume the result. This step forces the API cleanup that the later steps need and gives users something to run on day one. It also fixes the two items the author listed (absolute paths, `.gpx` suffix) and adds existence checks for inputs.

### Step 2. A project file and a structure library

Replace the three sample descriptions (`samples.csv`, `samples.toml`, `cifs/index.csv`) with one `xrdkit.toml` at the project root, created by `xrdkit init`, holding the instrument, the samples (file, name, composition, stage, form, temperature, Archimedes density) and the structures. Ship a structure library in the package so that a user writes `structure = "perovskite/cubic"`, `"perovskite/R3c"`, `"perovskite/Pbnm"`, `"ttb/P4bm"` and so on, each entry carrying the crystal system, the reflection conditions, the site kinds, Z, the free coordinates by Wyckoff position and sensible defaults for Uiso groups and origin, and a `--cif` override for anything else. `exchange` and `origin` become optional and are derived from the structure where they apply. A sample may list several phases. Add `xrdkit add-sample` so that dropping a file into `data/raw` and one command is all that registering a sample takes.

### Step 3. Generalise the cell

Introduce a `Cell` with six parameters and a crystal system, a general d spacing from the metric tensor already in `structure.py`, reflection generation with proper equivalence (so cubic (100) is one reflection, not two), and reflection conditions from a small table of the space groups the group uses (Pm-3m, P4mm, R3c and R3m in the hexagonal setting, Pbnm and Amm2, P4bm and P4/mbm), with pymatgen as an optional route for anything else. Make `refine_cell`, `refine_lattice`, `cell_volume` and `theoretical_density` accept it, derive Z from the structure entry, accept a formula string, complete `ATOMIC_MASSES` to the full periodic table, and add `relative_density`. This is the step that makes the kit useful for perovskites, and it is what `xrdkit lattice sample --archimedes 5.21` runs: peaks (fitted with `fit_profile` so that Kα2 no longer biases the positions), index against the structure's cell, refine with displacement for pellets and zero for powders by default, then print a, b, c, V, theoretical density and relative density with esds and write them to `results/lattice/lattice.csv`.

### Step 4. Move the refinement pipeline into the package

Take the generic 85 to 90 per cent of `refine_sample.py` into `xrdkit.pipeline`: Le Bail and Rietveld stage list helpers (`lebail_stages`, `rietveld_stages` by mode) alongside the existing `standard_stages`; multi phase jobs with phase fractions and a guard that holds the histogram scale when fractions are refined; the write up, summary and failure reporting; and an instrument aware start (displacement, radius from the instrument entry). Then `xrdkit lebail sample` gives lattice parameters by Le Bail for the compositions that warrant it, and `xrdkit rietveld sample --through occupancies` runs the full sequence. Change `find_undetermined` to judge against the CIF or original start, hold each site's total in the occupancy constraint, hold microstrain by default with `--mustrain` to free it, and add a preferred orientation flag for pellets. Keep `refine_sample.py` in XRD-Analysis only as a thin wrapper until the command replaces it.

### Step 5. Reporting

Add `xrdkit report sample1 sample2 ...` writing the group's Table 1 (composition, space group, a, b, c, volume, phase fraction, Rwp, Rp, chi²) from the result JSONs as CSV, markdown and a LaTeX or Word table, using the residuals the driver already records. Fix the `plot_rietveld` caption to use `accepted_stages`, add Rp, chi² and weight fractions to it, add a separate difference panel option and a `plot_rietveld_stack` for several compositions, and give `plot_stacked` per trace colours, a legend, an `xlim` argument and reference tick rows from a CIF or a stick list. Add `xrdkit figure rietveld sample` and `xrdkit figure stack ...` so that no figure needs a script.

### Step 6. Readers and an instrument registry

Add readers for `.xy`/`.xye`, Bruker `.brml` (a zip of XML, straightforward) and `.raw` (through the GSAS-II importer or xylib), and STOE exports, and make `read_xrdml` read every scan, `listPositions`, Kα2, the anode and the goniometer radius. Add `[instruments]` to the project file with wavelengths, radius, whether Kα2 is present and the instprm path, and `xrdkit instrument add aeris --standard data/standards/aeris/*.xrdml` to fit the LaB6 scan and write the instprm in one command, with a report of how the widths compare with previous standards. Let `write_instprm` write a Kα1 only file and the radius. Every command then takes `--instrument` or reads it from the sample entry, and Kα2 handling follows the instrument rather than assuming Cu.

### Step 7. Phase identification with and without a licence

Make `xrdkit purity sample` the routine check: peaks with a noise based threshold (a multiple of the square root of the background rather than two per cent of the strongest peak), indexed against the sample's structure, the unexplained peaks listed with their relative intensities, and a plain verdict. For the unexplained peaks offer two routes. `xrdkit phases sample --cod` searches the COD by the elements in the composition and precursors, with an intensity weighted figure of merit normalised to line count and a Kα2 aware tolerance, and caches simulated patterns so pymatgen runs once per CIF. `xrdkit phases sample --pdf` takes what a licensed user can export from HighScore or PDF-4 (a card list with d and I, or an exported reference pattern) and records the card number as the citable identification. Both routes draw reference tick rows under the pattern and write `results/phases/<sample>.md` with provenance.

### Step 8. Documentation and tests

Rewrite the user guide around the commands, one page per routine task, and keep the present library sections as a reference for people who script. Add synthetic pattern tests for each crystal system and each reader, a two phase Rietveld test, and a test that `xrdkit report` reproduces Table 1 from a stored result. Version this as 0.2.0.

## 5. What the routine day then looks like

After Step 7 the group's routine is six commands, none of which requires editing a script:

```
xrdkit add-sample data/raw/10_2026.xrdml --name x0.10_calcined --structure ttb/P4bm --composition "Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"
xrdkit check x0.10_calcined
xrdkit purity x0.10_calcined
xrdkit lattice x0.10_sintered --archimedes 5.21
xrdkit rietveld x0.10_sintered --through occupancies
xrdkit report x0.10_sintered x0.12_sintered --format latex
```

with `xrdkit stack` and `xrdkit figure` for the figures, and `xrdkit phases --cod` or `--pdf` when a purity check turns something up.

## 6. Notes on order and effort

Steps 1 and 2 are small and can be done first without touching the science. Step 3 is the largest single piece of new code and the one the perovskite work depends on, so it should follow immediately. Step 4 is mostly moving code that already works. Steps 5 to 7 can be taken in any order; the report writer in Step 5 is the quickest win for papers. Step 6 is worth doing before anyone but the Aeris users adopts the kit. The scientific fixes in Section 3.7 are folded into the step that touches the relevant code rather than done separately.
