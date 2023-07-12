# README.md (gcpy/benchmark)

This directory contains development materials for GEOS-Chem benchmarking.

| File or folder     | Description |
| --------------     | ----------- |
| `cloud/`           | Contains template config files for benchmarks on the AWS  cloud platform. |
| `config/`          | Contains configuration files with options for 1-month and 1-year benchmarks. |
| `modules/`         | Contains scripts for creating 1-year benchmarks. These are imported by the `run.benchmark.py` script. |
| `plot_driver.sh`   | Script to submit the `run_benchmark.py` script to a computational queue using the SLURM scheduler. |
| `run_benchmark.py` | Driver script for GEOS-Chem benchmarks (1-hour, 1-month, 1-year). |
