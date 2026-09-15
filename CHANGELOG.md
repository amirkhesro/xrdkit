# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `docs/USER_GUIDE.md` — a user guide covering the data each workflow needs in
  numbers, the files and formats to supply, and the three workflows end to end:
  plotting with hkl indices, lattice parameters and theoretical density by Le
  Bail, and Rietveld refinement. Every snippet in it was executed against real
  scans before it was written down.
- The `xrdkit` command, with its first subcommand, `xrdkit check SCAN
  [--json]`. It reads a `.xrdml` scan and prints the numbers `check_scan.py`
  prints in the user guide, followed by a verdict for each workflow (plotting,
  phase identification, Le Bail, Rietveld). A workflow the scan is not good
  enough for gets a list of the criteria it failed, each with its measured
  value and threshold. A missing or unreadable file exits with status 1.
- `xrdkit.quality`: `assess_scan` and `format_report`, which do the work behind
  `xrdkit check`, and `CRITERIA`, the thresholds of the user guide's Section 2
  table held in one place.
- `xrdkit plot SCAN`, which writes the peak list to `results/peaks_{stem}.csv`
  and the pattern to `figures/pattern_{stem}` (png and pdf). With `--cell A C`
  it also indexes the peaks from that tetragonal start cell, writes
  `results/indexed_{stem}.csv` and an hkl labelled `figures/pattern_hkl_{stem}`,
  and prints the refined cell and zero. `--space-group`, `--zero`, `--label` and
  `--no-satellites` tune it. The hkl figure's y axis runs 20 per cent of its
  span above the tallest point of the trace, so the label on the strongest peak
  stays inside the axes.
- `xrdkit stack SCAN [SCAN ...] --labels TEXT [TEXT ...]`, which writes the
  scans stacked to `figures/stack_{stem}`, with `--offset` and
  `--no-normalise`.
- Both figure commands take `--stem`, `--out DIR` (the output root, under which
  `results/` and `figures/` are created as needed), `--scale
  {linear,sqrt,log}` and `--json`. Each prints the files it wrote, and exits
  with status 1 when an input file is missing.
- `xrdkit density --formula TEXT --z N`, with either `--cell A C [--esd-cell
  EA EC]` for a tetragonal cell or `--volume V [--esd-volume EV]` for any
  symmetry, and optionally `--archimedes RHO [ESD]`. It prints the formula
  mass, the cell volume, the theoretical density and, given a measured one, the
  relative density, each with its esd where one exists, and writes one row
  holding every input, every result, the method and the date to
  `results/density.csv`, or `results/density_{stem}.csv` with `--stem`. `--json`
  prints the same fields.
- `xrdkit.density.parse_formula`, which reads a formula such as
  `"Sr0.4Ba0.5La0.1Nb1.9Ti0.1O6"` or `"Ca(OH)2"` into element counts, and
  `relative_density`, a measured density as a percentage of the theoretical one
  with the relative errors of the two added in quadrature.
- `ATOMIC_MASSES` covers every element from H to U: IUPAC 2021 conventional
  standard atomic weights, and for the elements with no stable isotope the
  mass number of the longest lived one.
- `xrdkit.library`, a structure library shipped as package data, one TOML file
  per entry under `xrdkit/structures/<family>/<name>.toml`. `list_entries()`
  names the entries, and `load_entry(name)` reads and validates one into a
  `StructureEntry` of `Site`s. It holds seven entries: the perovskites
  `perovskite/Pm-3m` (ideal cubic), `perovskite/P4mm` (tetragonal BaTiO3 type),
  `perovskite/R3c` (BiFeO3 type, hexagonal axes), `perovskite/Pbnm` (GdFeO3
  type, in the Pbnm setting) and `perovskite/Amm2` (orthorhombic BaTiO3 type),
  and the tetragonal tungsten bronzes `ttb/P4bm` (polar, from COD 2100720) and
  `ttb/P4mbm` (centrosymmetric, P4/mbm).
- `xrdkit.project`, the project file `xrdkit.toml`: `[project]`, and tables of
  instruments, structures (a structure library entry or a CIF, with composition,
  cell, z, exchange and origin) and samples (scan, instrument, structures, form,
  stage, temperature, Archimedes density, notes), every path relative to the
  folder holding the file. `find_project` finds the file from any folder inside
  the project, `load_project` reads and validates it into a `Project`, naming
  the table and key of anything wrong, and `results_dir` gives
  `results/<command>/<sample>`. No command reads the project file yet.
- `xrdkit init [--name TEXT]`, which writes `xrdkit.toml` in the current folder
  with `[project]` filled in and a commented example instrument, structure and
  sample, and makes `data/raw`, `cifs` and `results`. It refuses to overwrite
  an existing `xrdkit.toml`.

### Changed

- The peak to background criterion of 20 applies to phase identification
  only, not to plotting: a low ratio does not spoil a figure, but weak phases
  may not be visible. `xrdkit check` says so in its reason.
- The tests that read a measured scan take the folder from the
  `XRDKIT_TEST_RAW_DIR` environment variable, and skip with a message naming
  it when it is unset or the folder holds no `.xrdml` files. No path into a
  private data folder is written down in the suite any more.
- `.gitignore` also covers `*.gpx`, `*.xy`, `*.xye` and `*.cif`, and the
  `results/`, `output/`, `figures/` and `data/` folders, so a refinement run
  in a working copy cannot leave data or output staged by accident.

### Documentation

- `docs/USER_GUIDE.md` and `docs/GETTING_STARTED.md` follow a hand trial of both
  by a new user with no Python experience. Every code block in the user guide
  now says which script file it belongs to and whether it starts that file or
  continues it, Section 2 and each workflow end with the complete script in one
  piece, and every `save_figure` call prints the paths it wrote. The setup guide
  says that pip installs into Python rather than into a folder, adds VESTA as an
  optional viewer, and says what has to be done before the user guide is opened
  and what to change when a script is reused on the next dataset.
- `docs/USER_GUIDE.md` Section 5.2, Route A, is expanded to say that a licensed
  search and match is the authoritative phase identification and the citable
  record, that Route B serves routine screening, and that a CIF exported from
  the licensed database can seed Workflows 3 and 4.

### Fixed

- `formula_mass` and `theoretical_density` accept a formula string as well as a
  mapping of element to atoms per formula unit.
- `build_refine_job` makes every path it is given absolute, against the
  caller's working directory rather than the driver's, so relative paths work,
  and raises `Gsas2Error` naming the file when the data file, the instrument
  parameter file or a phase's CIF does not exist.
- A `.gpx` stem containing a dot, such as `x0.10`, keeps it in the default
  export prefix; only a `.gpx` extension is taken off.

### Planned

- A reader for plain text column formats (`.xy`, `.xye`) and for Bruker `.raw`
  and `.brml`; only `.xrdml` is read today.
- Indexing beyond tetragonal `P4bm`, which is the only space group whose
  reflection conditions `xrdkit.indexing` applies.
- Stage list helpers for Le Bail and Rietveld sequences; `standard_stages` is
  the instrument calibration sequence only, so the others are written out by
  hand.
- `plot_rietveld` to take its Rwp and GOF from the stage the final model came
  from, rather than from the last stage without an error entry, which is the
  wrong one whenever a stage was rejected.

## [0.1.0] - 2026-09-12

### Added

- `xrdkit.io` — reads PANalytical `.xrdml` scans into an `XRDScan`.
- `xrdkit.plotting` — publication figures for patterns, stacks, Caglioti fits
  and Rietveld fits, with peak marks and hkl labels.
- `xrdkit.peaks` — peak finding with Ka2 flagging and exclusion.
- `xrdkit.indexing` — reflection generation, peak indexing, cell refinement
  and a zero-offset search.
- `xrdkit.lattice` — least-squares lattice parameter refinement from indexed
  peaks.
- `xrdkit.density` — theoretical density from cell volume and formula mass.
- `xrdkit.sizestrain` — Scherrer sizes, Williamson-Hall analysis and
  instrumental resolution limits.
- `xrdkit.broadening` — profile and integral-breadth fitting, Caglioti fits
  and instrumental broadening corrections.
- `xrdkit.phases` — COD search and fetch, pattern simulation and candidate
  phase matching.
- `xrdkit.structure` — metric tensors, cell contents, site setup and
  interatomic distances.
- `xrdkit.config` — TOML settings files for samples and refinements, checked
  on load.
- `xrdkit.gsas2` and `xrdkit.gsas2_driver` — an unattended, stage-by-stage
  Rietveld refinement driver that runs GSAS-II under its own Python.

[0.1.0]: https://github.com/amirkhesro/xrdkit/releases/tag/v0.1.0
