# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
