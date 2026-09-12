# xrdkit

A reusable X-ray diffraction analysis toolkit for electroceramics research.

`xrdkit` provides shared, importable routines for loading, processing and
analysing XRD patterns, so that the same analysis code can be reused across
projects instead of being copied between one-off scripts.

New here? [docs/USER_GUIDE.md](docs/USER_GUIDE.md) says what data you need, in
what form, and how to get from a raw scan to each result.

Starting from a machine with nothing installed?
[docs/GETTING_STARTED.md](docs/GETTING_STARTED.md) sets up Python, xrdkit and
the folder layout on Windows and macOS.

## Installation

```bash
pip install xrdkit
```

Phase matching against the COD needs pymatgen, which comes with the `phases`
extra:

```bash
pip install "xrdkit[phases]"
```

GSAS-II is optional and is needed only for the Rietveld driver in
`xrdkit.gsas2`. It is not a Python dependency: install it separately and
point `XRDKIT_GSAS2_PYTHON` and `XRDKIT_GSAS2_HOME` at it, or leave it at
`~/gsas2main`. Every other module works without it.

## Data

This repository contains **code only**. Raw and processed diffraction data
(including `.xrdml` files) live in a separate data repository and are never
committed here.

## GSAS-II refinement

GSAS-II runs under its own Python, so `xrdkit.gsas2` writes a job as JSON and
runs `gsas2_driver.py` on it there (`run_job`), with `-B`; xrdkit writes
nothing into the GSAS-II installation. `find_gsas2` locates it
(`XRDKIT_GSAS2_PYTHON`, `XRDKIT_GSAS2_HOME`, or `~/gsas2main`).

GSAS-II's own version and update check, run whenever its scripting module is
imported, writes one bytecode file inside its installation
(`Lib/site-packages/scipy/__pycache__/__config__.cpython-313.pyc` in the
`~/gsas2main` install); it is written once, not on every run, and does not
come from xrdkit.

`build_refine_job` assembles a refinement from a list of stages, each adding
flags to the ones before (see the `gsas2_driver` module notes for every key):

| Stage key | Refines |
| --- | --- |
| `background`, `scale`, `zero`, `displacement`, `instrument` | histogram and instrument |
| `cell`, `le_bail`, `phase_fractions`, `size`, `mustrain` | phase, by name or all |
| `overall_uiso` | one Uiso for every atom of a phase |
| `uiso_groups` | `{phase: [[site, ...], ...]}`: one Uiso for every atom on each group's sites |
| `atoms`, `atom_flags` | F, X, U for every atom, or for named atoms |
| `coordinates` | `{phase: {site: "xz" or "all"}}`: the coordinates to free on each site; the rest are held, atoms sharing a site move together |
| `origin` | `{phase: {"site": label, "axis": "z"}}`: a site's coordinate held to fix the origin; a stage that would leave the origin floating along a polar axis fails |
| `occupancies` | `[{"phase", "sites": [...], "elements": [...]}]`: the named elements' occupancies on the named sites, each element's content over them held (so the vacancies over them too) |

Job options besides the stages: `limits`, `cycles`, `broadening`,
`background_start`, `scale_start`, `overall_uiso_start`, `le_bail_cycles`,
`max_passes` and `pass_tolerance` (refine each stage until it settles),
`on_flagged` and `on_unsettled` (below), `sanity` (`{"max_shift": 0.05,
"reference": {phase: [{"label", "xyz"}]}}`) and `bonds` (`True` or
`{"anions", "dmax", "dmin"}`). A phase's `atoms` edits the structure read
from its CIF in the project only: occupancies, positions and Uiso, and atoms
added on existing sites; `structure_edits` gives the edits that set up a
refined structure again, so that one refinement starts where another ended.

Every stage of the result records its residuals and parameters with esds, the
atoms as it left them, the coordinates it refined and held by site, its
occupancy constraints, and a sanity check: negative Uiso, occupancies outside
0 to 1 (and site totals above 1), and sites moved more than `max_shift`
(fractional) from the reference. With `bonds`, `run_job` adds the final
model's cation to anion distances (`xrdkit.structure.bond_lengths`, to oxygen
up to 3 Å by default) to the result.

### What becomes of a stage that goes wrong

So that a sequence can be run unattended and still leave a usable record,
every stage carries a `status` and the run goes on from the last stage kept:

| Status | The stage | The run |
| --- | --- | --- |
| `clean` | settled, and raised no sanity flag the stage before it lacked | kept |
| `unsettled` | reached `max_passes` still moving, and raised no new flag | kept |
| `flagged` | raised a new sanity flag, with `on_flagged` `"accept"` | kept |
| `rejected` | raised a new flag under `on_flagged` `"reject"` (the default), or did not settle under `on_unsettled` `"reject"` | rolled back |
| `failed` | raised | rolled back, and the run ends |

**Rejected.** A sanity flag means the model has gone somewhere it should not
be — a Uiso below zero, an occupancy outside 0 to 1, a site further from the
reference than `max_shift` — so the stage that raised one is rejected by
default. It is recorded in full, residuals, parameters and atoms, with
`rejected_because` and its name in the result's `rejected`. The project then
goes back to the state the last stage kept left it in — parameters, atoms
and constraints — and is saved so, and the stage's own entries are dropped
from the stages, so that what it alone refined is held in every stage after
it. `on_flagged: "accept"` keeps it instead, as `"flagged"`.

**Unsettled.** Reaching the pass cap only means the stage was still moving
when the passes ran out, which is not in itself wrong, so it is kept by
default and `largest_remaining_move` records the parameter that moved most
in its last pass and by how many esds. `on_unsettled: "reject"` rolls one
back like a flagged stage, for a caller that wants nothing unsettled in its
model.

**Undetermined.** Each stage also records, as `undetermined`, the refined
atom parameters the data do not determine
(`gsas2_driver.find_undetermined`): an occupancy whose esd is more than half
its allowed range, 0 to 1, and a coordinate or an isotropic Uiso whose esd
is larger than its shift from where the job found it. Such a value is not a
result, whatever its stage's status. The result's `undetermined` are those
of the last stage kept.

**The model a run leaves.** Every run ends with `final`, and with
`export_prefix` an `exports`, whatever became of its stages: the model of
the last stage kept, or, where no stage was kept, the one the job started
from, computed by a refinement of no cycles that moves nothing. The result's
`final_from` names the stage they come from, or `"job start"`; a starting
model has no refined parameter and so no esd.

### Writing a sequence up

For a pipeline that runs several refinements in turn, `xrdkit.gsas2` writes
what they came to:

| Function | Gives |
| --- | --- |
| `stage_status` | one stage's status, worked out for a result written before the driver recorded one |
| `accepted_stages` | the stages of a result that were refined and kept |
| `stage_statuses` | every stage as `{name, status, reason, passes, rwp, gof, undetermined}`, `reason` saying why it is not clean |
| `stage_status_table` | those as markdown table lines |
| `log_tail` | the last lines of a GSAS-II log |
| `failure_markdown` | a write up of a refinement that could not be finished: the error, the stages it got through and the log tail |
| `summary_markdown` | a summary of a sequence: a row per refinement with its final Rwp, GOF and key values, the status of each stage, then every stage not clean with why and the parameters left undetermined |

## Refinement settings

`xrdkit.config` reads the settings of a refinement pipeline from a TOML file
with the standard library's `tomllib` (`load_config`): one table per sample
under `samples` and one per reference structure under `structures`. It
checks them as it reads, and a key missing or unknown, a value of the wrong
kind, or a sample composition that does not fit its structure's sites (an
element on no site and not added by the rule, or more atoms per cell on a
kind of site than it has positions) raises `ConfigError` naming the file and
the table. `sample_settings` finds a sample by its id or table name, with its
structure.

| Table | Keys |
| --- | --- |
| `samples.<name>` | `id`, `scan`, `composition` (atoms per formula unit), `structure`, `start_cell` (`{file, model}` or `{a, c}`), `two_theta`, `background` (`{function, terms}`), `refine_microstrain`, `notes`; optionally `followed_reflections`, `trials` (`{runs = [{low, terms}], followed}`), `write_up` (text by mode and section), `unsettled` (the rule for each mode, over the top level one) |
| `unsettled` | at the top level, the default rule for a stage of each mode that has not settled, `"accept"` (the pipeline's default) or `"reject"`, by the caller's own mode names; each sample carries it merged with its own as its `unsettled` |
| `structures.<name>` | `cif`, `label`, `phase_name`, `space_group`, `formula_units`, `sites` (`{atoms = {label = element}, wyckoff, kind}`, kind A, B or O), `uiso_groups` (`{name, sites}`), `origin` (`{site, axis}`), `exchange` (`{elements, sites}`), `composition` (`{added = {element = host}}`); optionally `free_coordinates` (by Wyckoff position), `bond_limits` (`{kind = {min, max}}`) |

From a structure table, `xrdkit.structure.site_setup` finds the sites among a
phase's atoms as the driver reports them, checking that each site's atoms
share one position, are of the elements given and have the multiplicity of
its Wyckoff position, and that no atom is left over. It returns the sites by
kind, the Uiso groups, the coordinates to refine on each site, the origin site
and the exchange, in the forms the stage keys above take.
`composition_edits` gives the atom edits that put a nominal composition on
the sites by the table's rule: each element the CIF holds is scaled by one
factor over its sites, which keeps its distribution, and each added element
goes on its host's sites in proportion to the host's occupancy.
`cell_contents` gives the atoms of each element per cell.

For example, the sample pipeline of the XRD analysis repository reads every
sample and structure setting from its `config/samples.toml` and runs Le Bail,
fixed atoms, coordinates and occupancies in turn, each from the saved result
of the one before:

```
python -B scripts/refine_sample.py 10                          # Le Bail only
python -B scripts/refine_sample.py 10 --trials                 # Le Bail with the configured trial ranges and background terms
python -B scripts/refine_sample.py 10 --through occupancies    # all four modes
python -B scripts/refine_sample.py 10 --from coordinates       # one mode, from the saved fixed atoms result
python -B scripts/refine_sample.py 10 --from coordinates --through occupancies
```

It passes each mode's `unsettled` rule and pass cap to its job, writes the
mode's write up with `failure_markdown` and stops with a non-zero status
where a mode cannot be finished, and writes `summary.md` beside the write
ups with `summary_markdown` at the end of every run.

## Status

Under active development. The API is not yet stable and may change without
notice.

## Citing

If xrdkit contributes to work you publish, please cite it. The metadata is in
[CITATION.cff](CITATION.cff), which GitHub renders as a ready-made citation
under **Cite this repository**.

## License

MIT — see [LICENSE](LICENSE).
