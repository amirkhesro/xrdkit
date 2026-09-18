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
  and the pattern to `figures/pattern_{stem}` (png and pdf). With `--cell`
  it also indexes the peaks from that start cell, writes
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
- `xrdkit density --formula TEXT --z N`, with either `--cell` (see below)
  or `--volume V [--esd-volume EV]` for any symmetry, and optionally `--archimedes RHO [ESD]`. It prints the formula
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
- In the project file, `ka2 = true` needs exactly two wavelengths and
  `ka2 = false` exactly one, and an `instprm` that is given must exist.
  `xrdkit.project` also gains `resolved_z` (a structure's z, or its library
  entry's), `resolved_cell` (all six cell parameters from a library entry's
  crystal system), `load_project_text` and `toml_string`.
- `xrdkit add-sample FILE --structure KEY [--structure KEY ...] [--name KEY]
  [--instrument KEY] [--stage TEXT] [--form powder|pellet] [--temperature-c N]
  [--archimedes N] [--notes TEXT]`, which appends a `[samples.KEY]` table to
  the project file found from the current folder, leaving the rest of the file
  byte for byte, and prints it. FILE must lie inside the project folder and is
  stored relative to it; the key defaults to the file's stem and the instrument
  to the only one. The file with the table added is validated before anything
  is written, and an existing key is refused.
- `xrdkit check`, `plot` and `stack` take a sample key of the project file
  wherever they take a scan path, and `xrdkit density` takes one as an optional
  `SAMPLE` argument. An existing file is used exactly as before; any other
  argument without a folder separator is looked up in `[samples]` of the
  project found from the current folder. A sample is read at its instrument's
  K alpha 1 wavelength (its K alpha 2 over K alpha 1 ratio flags the
  satellites, and an instrument without K alpha 2 flags none), its outputs go
  to `results/<command>/<key>` under the project root unless `--out` is given,
  `check` heads its report with the key and composition and keeps it there,
  `plot` labels the trace with them, and `density` takes the formula, z, cell
  and Archimedes density of the sample and its first structure. Options given
  on the command line win. `stack` takes keys and paths mixed. `xrdkit plot`
  gains `--wavelength`. Given a sample and no `--cell`, `plot` indexes from
  the cell, crystal system and space group of the sample's first structure; if
  that indexing fails it prints a one line note and plots without hkl labels.
- `xrdkit.cell.Cell`, a frozen unit cell of any crystal system: six
  parameters and a crystal system, checked against each other on
  construction. Named constructors (`Cell.cubic`, `tetragonal`,
  `orthorhombic`, `hexagonal`, `trigonal` on hexagonal axes, `monoclinic` with
  b unique, `triclinic`), `from_parameters` and `from_entry`; the direct and
  reciprocal metric tensors, `volume`, `d_spacing` and the vectorised
  `d_spacings`, `parameters` and `parameter_names` (the free ones), `replace`
  and `to_dict`.
- `xrdkit.symmetry`, which works out systematic absences, Laue equivalence and
  multiplicity from the operations of a space group rather than from written
  out conditions: `parse_xyz`, `space_group_operations` (Pm-3m, P4mm, P4bm,
  P4/mbm, R3c and R3m on hexagonal axes, Pbnm and Amm2, each checked against
  its order, and against GSAS-II for every hkl from -6 to 6), `is_absent`,
  `laue_group`, `holohedry`, `laue_orbit`, `multiplicity` and
  `representative`.
- `generate_reflections` works for any crystal system and supported space
  group: equivalent reflections are merged into one, labelled with a
  conventional representative (h >= k >= 0, l >= 0 for a tetragonal cell) and
  carrying the family's `multiplicity`, a new field of `Reflection`.
- `refine_cell` fits the reciprocal metric of any crystal system by linear
  least squares, holding a component no peak carries information on at the
  start cell's value and naming the components when the peaks cannot separate
  them. `refine_lattice` refines the free parameters of any crystal system with
  the zero and the specimen displacement.
- `LatticeFit` gives all six cell parameters with their esds, the volume and
  its esd, the free parameter names and the scaled covariance matrix. The
  volume esd is propagated through the covariance, so it now includes the
  correlation between the cell parameters, which the earlier propagation from
  a and c alone left out.
- `xrdkit plot --cell` and `xrdkit density --cell` take one to six numbers,
  the free parameters of the crystal system the count implies: a (cubic), a c
  (tetragonal), a b c (orthorhombic), a b c alpha beta gamma (triclinic).
  `--system` settles the counts that leave a choice: two numbers are a and c of
  a hexagonal or trigonal cell with `--system hexagonal` or `trigonal`, and
  four are a b c beta with `--system monoclinic`. `density --esd-cell` takes as
  many esds, in the same order. `plot --space-group` takes any symbol; one
  whose reflection conditions are not known is indexed without them, with a
  one line note. `density` prints the crystal system and cell it used, and its
  CSV gains `crystal_system` and all six cell parameters with their esds.
- `xrdkit lattice SCAN`, the lattice parameter workflow. It finds the peaks,
  refits each position with `fit_profile` as the K alpha 1 line of its
  doublet (a single line for an instrument without K alpha 2), keeping the
  found position where a fit fails or moves it by more than the peak's FWHM,
  indexes them against the start cell from `--cell` or the sample's first
  structure with the adaptive coarse window, and refines the cell with
  `refine_lattice`. A pellet refines the specimen displacement with the zero
  held at `--zero` (default 0); a powder, or a scan that is not a sample,
  refines the zero, and `--displacement` frees the displacement too. The
  goniometer radius comes from the instrument or `--radius`. With a formula
  and Z, from the options or the structure, it gives the theoretical density
  from the refined volume and its correlated esd, and with an Archimedes
  density the relative density. It prints the peak counts, the coarse window,
  the cell, volume, zero and displacement with their esds, the rms and the
  densities, writes `results/lattice/peaks_STEM.csv` (found and fitted
  positions, the fit, hkl and difference for every peak) and appends one row
  per run, with every input, result, method, date and version, to
  `results/lattice/lattice.csv`, or `results/lattice/KEY/lattice_KEY.csv`
  for a sample. Takes `--cell`, `--system`, `--space-group`,
  `--wavelength`, `--formula`, `--z`, `--archimedes`, `--zero`,
  `--displacement`, `--radius`, `--no-satellites`, `--stem`, `--out` and
  `--json`.

### Changed

- `xrdkit rietveld` refuses to replace an existing result: when the result
  JSON of a mode it is to run is already in the output folder it stops before
  refining, names the files and returns 1. `--overwrite` replaces them. It used
  to overwrite `results/rietveld/KEY` without a word.
- Breaking: `flag_kalpha2` flags a peak only when it lies within half its
  nearest lower parent's FWHM (`KALPHA2_POSITION_TOLERANCE`) of that parent's
  K alpha 2 position and its height above the background over the parent's
  lies in 0.2 to 0.8 (`KALPHA2_INTENSITY_BAND`). `tolerance` in degrees is
  replaced by `position_tolerance`, a fraction of the FWHM,
  `wavelength_ratio=None` flags nothing, and `Peak` gains `background`, the
  lowest intensity within `BACKGROUND_HALF_WIDTH` of the peak.
- Breaking: `xrdkit lattice` leaves peaks flagged as K alpha 2 satellites out
  of the indexing and refinement, where they used to enter both at their raw
  found positions. A flagged peak that the refined cell puts within the
  indexing tolerance of a reflection is refitted, and recovered only if the
  refit is accepted, when the indexing and refinement run once more with it; a
  rejected refit leaves it a satellite. `--no-satellites` now skips that
  recovery, and no longer drops the flagged peaks from the peaks file and the
  found count. The report prints the satellites excluded, the peaks recovered
  and the indexed fraction; the peaks file gains `recovered` and the results
  row `n_satellites` and `n_peaks_recovered`, so results files written before
  must be regenerated rather than appended to. `n_peaks_indexed` now counts
  only the peaks that took part.
- Breaking: `cell_volume(cell, esd=None)` takes a `Cell` and a mapping of free
  parameter name to esd, in place of a tetragonal cell and `esd_a` and `esd_c`.
- Breaking: `CellFit.c_fitted` is replaced by `CellFit.held`, the names of the
  parameters held at the start cell's value, empty when every one was fitted.
- Breaking: the `space_group` argument of `generate_reflections`,
  `index_peaks`, `estimate_zero_offset` and `index_and_refine` defaults to
  `None`, no reflection conditions, instead of `"P4bm"`. Code that relied on
  the old default must pass `space_group="P4bm"`. `xrdkit plot --space-group`
  no longer defaults to P4bm either, and `DEFAULT_SPACE_GROUP` is gone from
  `xrdkit.indexing`.
- Breaking: `LatticeFit` has new fields, and `cell`, now a stored `Cell`
  rather than a property, is its first field, so building one by position
  changes. `a`, `c`, `esd_a` and `esd_c` remain; an esd of a parameter the
  crystal system fixes is `None`. `lattice_fit_to_dict` writes the new fields,
  and `c_over_a` only for tetragonal, hexagonal and trigonal cells.
- `TetragonalCell` is deprecated in favour of `Cell.tetragonal`; it is now a
  function that returns one, and `TTB_CELL` is a `Cell`.
- `xrdkit plot` on a sample whose structure cannot be indexed plots the
  pattern and returns 0 with a note, where it used to fail; with `--cell` given
  a failure to index is still an error.

- The Rietveld settings of `xrdkit.config` no longer assume a tetragonal P4bm
  bronze. A configuration that validated before gives the same result; each
  structure also carries `library` and `crystal_system`, None unless given.
  - `exchange` is optional. Without it the structure's `exchange` and
    `site_setup`'s `exchange` are None, and the GSAS-II driver takes
    `"occupancies": None` in a stage as no constraint.
  - `origin` is optional. Without it no coordinate is held for the origin:
    `site_setup` gives `origin` and `origin_axis` as None and lists the free
    coordinates of every site. The driver takes `"origin": None` in a stage as
    no change.
  - A structure's sites need not include one of each kind A, B and O.
    `bond_lengths` of a phase with no anion site returns no bonds, with a
    warning, instead of passing an empty target list on.
  - A sample's `start_cell` is `{file, model}` or exactly the cell parameters
    of the structure's crystal system, each greater than 0. `a = 0` is now
    refused. The crystal system comes from the structure's new optional
    `crystal_system`, or its library entry. Without either it is inferred
    from `a` (cubic), `a, c` (tetragonal) or `a, b, c` (orthorhombic), and any
    other set of parameters asks for `crystal_system`.
  - `formula_units` is optional. Without it the composition is not checked
    against the capacity of the sites, and `composition_edits` refuses, naming
    the structure, instead of assuming a Z.
  - A structure may give `library = "<entry>"` and `atoms`, the CIF's atoms on
    each of the entry's sites by label, instead of `space_group`, `sites` and
    `uiso_groups`. The sites, Wyckoff positions, kinds, free coordinates, Uiso
    groups, Z, crystal system and origin then come from the entry, and a key
    given explicitly overrides the entry's value. Sites keep the names of
    their first atoms, and `uiso_groups`, `origin` and `exchange` may name a
    site by its library label instead.
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
- The peaks file of `xrdkit lattice` gains `corrected_two_theta` and
  `d_spacing` after `esd_fitted_two_theta`, `relative_intensity` after
  `intensity`, and `n_candidates` last, so that it is a complete indexing
  table; the results file `lattice.csv` keeps its columns.

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
- `docs/USER_GUIDE.md`: the TTB indexing examples pass `space_group` explicitly,
  and `lattice_density.py` is updated to the general `cell_volume` signature.

### Fixed

- The undetermined list never named a Uiso for a structure whose CIF gives
  anisotropic Uij: `fixed_atoms` makes every atom isotropic before it refines
  but recorded the CIF atoms, which carry no Uiso, as the start model. The
  start model now records the Uiso every atom is started at.
- `xrdkit lebail` and `xrdkit rietveld` given a relative `--out` refined the
  sample and then failed writing the result JSON with `FileNotFoundError`:
  the GSAS-II driver runs in a work folder of its own and took the path from
  there. `--out` is now made absolute, against the folder the command is run
  in, before anything else.
- The Rietveld pipeline's start cell reads the lattice results the lattice
  command writes. It looked for columns `a` to `gamma`, but the command writes
  `a_angstrom`, `b_angstrom`, `c_angstrom`, `alpha_deg`, `beta_deg` and
  `gamma_deg`, so `resolve_inputs` raised "its last row gives no whole cell"
  on every sample with lattice results. The plain names of older files are
  still read.
- The Rietveld set up places an element that a structure's `atoms` table adds
  only beside the CIF atoms the table names on that element's sites, with the
  label the table gives it. It used to reduce the table to a host element and
  add the element beside every atom of that element, so La placed on A1 beside
  Sr1 also landed beside Sr2 on A2. The element's whole content now goes on
  the named sites, shared in proportion to the hosts' CIF occupancies when it is
  placed on more than one. `composition_edits` takes, besides a host element,
  host atoms as `{host label: added label}`. A site of the table that holds an
  element the table does not place there, in the CIF or after the edits, stops
  the set up with a `PipelineError` naming the site and element. Structures
  without `atoms` are unchanged.
- The host of an element an `atoms` table places on more than one site is
  chosen by element, not by the order of the table: it is the one element all
  those sites hold, beside whose CIF atom on each site the element goes, so Ca
  on A1 and A2 goes beside Sr1 and Sr2 even when A2 lists Ba2 first. Sites with
  more than one element in common, or none, make `load_project` raise
  `ConfigError` naming the element, the sites and any candidate hosts. An
  element on a single site goes beside any CIF atom of it, the host deciding
  only its label. `xrdkit.config.host_elements` holds the rule.
- `formula_mass` and `theoretical_density` accept a formula string as well as a
  mapping of element to atoms per formula unit.
- `build_refine_job` makes every path it is given absolute, against the
  caller's working directory rather than the driver's, so relative paths work,
  and raises `Gsas2Error` naming the file when the data file, the instrument
  parameter file or a phase's CIF does not exist.
- A `.gpx` stem containing a dot, such as `x0.10`, keeps it in the default
  export prefix; only a `.gpx` extension is taken off.
- Exactly coincident reflections from different Laue orbits, such as (553) and
  (713) or (860) and (10 0 0) in a tetragonal cell, come out of
  `generate_reflections` in hkl order, lowest first, and a peak on them is
  indexed as the lowest. Their angles differ by about 1e-14 degrees through the
  metric, which used to decide the order and so the label. Angles, and distances
  from a peak, within the new `COINCIDENCE_TOLERANCE` of 1e-9 degrees now count
  as equal and are ordered by hkl. Which reflections are generated is unchanged.

### Planned

- A reader for plain text column formats (`.xy`, `.xye`) and for Bruker `.raw`
  and `.brml`; only `.xrdml` is read today.
- Reflection conditions for space groups beyond the eight in
  `xrdkit.symmetry`, through pymatgen as an optional dependency.
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
