#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH -p sapphire,huce_cascade,seas_compute,shared
#SBATCH --mem=180000
#SBATCH --mail-type=END

#============================================================================
# This us a sample SLURM script that you can use to run the GCPy
# benchmark plotting code as a SLURM batch job.
#
# You can modify the SLURM parameters above for your setup.
#
# Tips:
# -----
# (1) Use fewer cores to reduce the memory footprint. This may prevent
#     your job from running out of memory.  Python under Linux seems
#     to have an issue where not all memory is released back to the OS.
#
# (2) We recommend that you generate only one benchmark comparison
#     (GCC vs GCC, GCHP vs GCC, GCHP vs GCC, or diff of diffs)
#     at a time.  Otherwise your job will probaly run out of memory.
#
# (3) For diff-of-diffs plots, we recommend using 6 cores.
#============================================================================

# Apply all bash initialization settings
. ~/.bashrc

# Make sure to set multiple threads; Joblib will use multiple
# cores to parallelize certain plotting operations.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=500m

# Use a non-interactive backend for matplotlib (we're printing to file)
export MPLBACKEND=agg

# Turn on Python environment (edit for your setup)
conda activate gcpy_env

# Specify a YAML file with benchmark options
# Uncomment the file that you wish:
config="1mo_benchmark.yml"
#config="1yr_fullchem_benchmark.yml"
#config="1yr_tt_benchmark.yml"

# Call the run_benchmark script to make the plots
python -m gcpy.benchmark.run_benchmark "${config}" > "${config/.yml/.log}" 2>&1

# Turn off python environment
conda deactivate

exit 0

