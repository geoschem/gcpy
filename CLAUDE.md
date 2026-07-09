# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What GCPy is

GCPy (`geoschem-gcpy`) is a Python toolkit for working with output from the
[GEOS-Chem](https://geos-chem.readthedocs.io) atmospheric chemistry model. It is scoped
narrowly and deliberately:

- **In scope**: plots/tables from GEOS-Chem output, benchmark simulation comparisons,
  horizontal/vertical grid utilities, GCHP cubed-sphere regridding, example/community scripts.
- **Out of scope**: general NetCDF manipulation (use xarray/NCO/CDO instead), statistical
  analysis (use scipy/scikit-learn/R), machine learning (use pytorch/tensorflow/julia). Don't
  add these kinds of features to GCPy itself.

## Environment setup

GCPy depends on a pinned conda/mamba environment (cartopy, xesmf, esmf, xarray, etc. — many
of these are not pip-installable in a compatible way, so `pip install` alone will not work).

```bash
# environment.yml is a symlink to docs/environment_files/gcpy_environment_py313.yml
mamba env create -n gcpy_env --file=environment.yml
conda activate gcpy_env
pip install -e .
```

Sanity check that the install worked:

```bash
python -c "import gcpy"
```

There is also a Python 3.12 environment file at `docs/environment_files/gcpy_environment_py312.yml`.
CI (`.github/workflows/build-gcpy-environment-py31{2,3}.yml`) builds these environments and runs
the same `import gcpy` check on Python 3.10–3.14 — there is currently no automated unit test
suite in this repo, so "does it build and import" is the practical correctness bar enforced by CI.

## Linting

Style follows PEP 8, enforced via the repo-root `.pylintrc`. Run `pylint` on any file you modify
(it picks up `.pylintrc` automatically from the current directory):

```bash
pylint gcpy/<path_to_file>.py
```

`.pylintrc` deliberately relaxes checks that don't match GCPy's conventions rather than being a
generic default config — don't "fix" these if pylint stays quiet about them:

- Structural-size metrics are disabled: line length, module length, and argument/local/branch/
  statement/return counts. GCPy modules and functions are large and numerically dense by nature;
  these checks were flagging normal code, not real problems.
- `wildcard-import`/`unused-wildcard-import` are disabled because every `gcpy/**/__init__.py`
  does `from .X import *` by design, to re-export submodule contents at the package level.
- `good-names` in `.pylintrc` allowlists short scientific variable names used throughout the
  codebase (`ds`, `lat`, `lon`, `ax`, `nx`, `ny`, loop var `i`, etc.) plus `AP`/`BP` (hybrid-grid
  coefficient names in `gcpy/grid.py` that mirror standard GEOS-Chem terminology).
- Private (`_`-prefixed) helper functions are exempt from the missing-docstring check; public
  functions/classes still need numpydoc-style docstrings, which is the norm across the codebase.

## Running benchmarks

The benchmark driver is invoked as a module, not a plain script, because it relies on relative
imports from `gcpy.benchmark.modules`:

```bash
conda activate gcpy_env
python -m gcpy.benchmark.run_benchmark <path-to-benchmark-config.yaml>
```

Config templates live in `gcpy/benchmark/config/` (e.g. `1yr_fullchem_benchmark.yml`). If run
without an X display (e.g. in a queue), set `QT_QPA_PLATFORM=offscreen` before running, or
matplotlib will try to open a window. `gcpy/benchmark/benchmark_slurm.sh` submits the driver via
SLURM.

## Building docs

Docs are Sphinx-based, built from `docs/source/conf.py`. ReadTheDocs builds them using the conda
env at `docs/environment_files/gcpy_environment_py313.yml` (see `.readthedocs.yaml`). Locally:

```bash
cd docs
make html
```

## Architecture

### `gcpy/` top-level modules

These are the core, standalone utility modules, all re-exported at the package level via
`gcpy/__init__.py` (`from .X import *`), so `gcpy.<function>` works directly without knowing
which submodule defines it:

- `grid.py` / `vgrid_defs.py` — GEOS-Chem horizontal and vertical grid definitions.
- `cstools.py` — cubed-sphere grid detection/manipulation (originally from Liam Bindle /
  Sebastian Eastham).
- `regrid.py` — builds xESMF regridder objects; `file_regrid.py` is the higher-level driver
  that regrids a whole file between lat/lon and/or cubed-sphere (including stretched) grids,
  usable as a CLI via `argparse`.
- `regrid_restart_file.py` — combines a restart file, ESMF regrid weights, and a template file
  (plus optional stretched-grid params) to produce a regridded GCHP restart file.
- `grid_stretching_transforms.py` — low-level vector-rotation math used for stretched-grid work.
- `units.py` — unit conversions, primarily to support benchmarking.
- `util.py` — general xarray/numpy helper utilities (largest module; includes variable-type
  verification helpers used throughout the other modules, e.g. `verify_variable_type`).
- `date_time.py` — datetime/string handling utilities.
- `raveller_1D.py`, `append_grid_corners.py` — cubed-sphere satellite-track and grid-corner
  utilities.
- `constants.py` — physical/chemical constants and shared global variables.

Since everything is star-imported into the top-level namespace, when adding a new top-level
module remember to add its import to `gcpy/__init__.py`, and watch for name collisions across
modules.

### `gcpy/plot/` — plotting subsystem

`core.py` holds shared state/helpers used by the rest of the subpackage (colormaps loaded from
`colormaps/`, the `gcpy_style` matplotlib stylesheet, panel-naming helpers). `single_panel.py`,
`six_plot.py`, `compare_single_level.py`, and `compare_zonal_mean.py` build on top of `core.py`
to produce the standard GEOS-Chem comparison plots (these are what the benchmark modules call
into to generate figures).

### `gcpy/benchmark/` — benchmark report generation

This is the most complex subsystem and spans multiple files that must be read together:

- `run_benchmark.py` — the single entry point/driver (`python -m gcpy.benchmark.run_benchmark
  <config.yml>`). It reads a YAML config and orchestrates 1-hour/1-month/1-year benchmarks for
  GCC-vs-GCC, GCHP-vs-GCC, GCHP-vs-GCHP, and GCHP-vs-GCC-diff-of-diffs comparisons.
- `modules/` — the actual benchmarking logic, imported by `run_benchmark.py` rather than run
  standalone: `benchmark_funcs.py` (core comparison logic), `benchmark_utils.py`,
  `benchmark_mass_cons_table.py`, `benchmark_species_changes.py`, `benchmark_drydep.py`,
  `benchmark_models_vs_obs.py`, `benchmark_models_vs_sondes.py`, `benchmark_gcclassic_stats.py`,
  `benchmark_scrape_gc{classic,hp}_timers.py`, `oh_metrics.py`, `ste_flux.py`, `budget_*.py`.
  Several `.yml`/`.csv` files here (e.g. `benchmark_categories.yml`, `emission_species.yml`,
  `lumped_species.yml`, `aod_species.yml`, `GC_72_vertical_levels.csv`) are data/config read by
  those modules, not code.
- `config/` — YAML configs (per benchmark type/duration) that parameterize `run_benchmark.py`
  runs; `cloud/` — AWS-specific template configs; `benchmark_slurm.sh` — SLURM submission script.

### `gcpy/kpp/`, `gcpy/profile/`, `gcpy/community/`

Independent, narrowly-scoped tool collections, not part of the core `gcpy` import surface:

- `kpp/` — utilities for KPP solver analysis output (`kppsa_*.py`).
- `profile/` — parses/plots profiling output from gprofng and Intel VTune
  (`vtune_*.py`, `gprofng_functions.py`).
- `community/` — user-submitted scripts; each has its own author of record noted in
  `gcpy/community/README.md` — contact that author with questions rather than assuming GCST
  ownership.

### `gcpy/examples/`

Standalone example scripts grouped by topic (`diagnostics/`, `grids/`, `hemco/`, `plotting/`,
`timeseries/`, `working_with_files/`, `xarray_examples/`, `dry_run/`). These demonstrate API
usage for end users; they are not exercised by CI and aren't part of the package's public
interface.

## Contribution conventions worth knowing

- Development happens off `main` via PRs (GitHub Flow); `dev` is the integration branch.
- Any user-facing change should get a one-line entry in `CHANGELOG.md`.
- Source modules generally start with a docstring giving attribution/citation context (several
  modules were contributed by named authors, e.g. `cstools.py`, `community/*`) — preserve that
  when editing.
