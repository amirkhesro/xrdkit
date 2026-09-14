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

### Changed

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
- `docs/USER_GUIDE.md` Section 5.1 notes a licensed search and match as the
  authoritative phase identification route, alongside the COD route.

### Fixed

- Nothing yet.

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
- `formula_mass` to accept a formula string as well as a dictionary of element
  to atoms per formula unit.
- `build_refine_job` to accept relative paths, resolving them against the
  caller's working directory rather than the driver's.

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
