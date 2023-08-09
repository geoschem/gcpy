#!/bin/bash

#SBATCH -c 12
#SBATCH -N 1
#SBATCH -t 0-4:00
#SBATCH -p huce_intel
#SBATCH --mem=50000
#SBATCH --mail-type=END

#============================================================================
# This us a sample SLURM script that you can use to submit
# the run_1mo_benchmark.py or the run_1yr_benchmark.py
# script to a computational queue.
#
# You can modify the SLURM parameters above for your setup.
#============================================================================

# Apply all bash initialization settings
. ~/.bashrc

# Make sure to set multiple threads; Joblib will use multiple
# cores to parallelize certain plotting operations.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=500m

# Turn on Python environment (edit for your setup)
conda activate gcpy

# Uncomment this line to make 1-month benchmark plots & tables
./run_benchmark.py 1mo_benchmark.yml > plots_1mo.log

# Uncomment this line to make 1-year benchmark plots & tables
# ./run_benchmark.py 1yr_fullchem_benchmark.yml > plots_1yr_fullchem.log

# Uncomment this line to make 1-year TransportTracers plots & tables
# ./run_benchmark.py 1yr_tt_benchmark.yml > plots_1yr_tt.log

# Turn off python environment
conda deactivate

exit 0

