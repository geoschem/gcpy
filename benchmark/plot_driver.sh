#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH -p huce_intel
#SBATCH --mem=12000
#SBATCH --mail-type=END

#============================================================
# This us a sample SLURM script that you can use to submit
# the run_1mo_benchmark.py or the run_1yr_benchmark.py
# script to a computational queue.
#
# You can modify the SLURM parameters above for your setup.
#============================================================

# Turn on Python environment (edit for your setup)
source activate gcpy

# Uncomment this line to make 1-month benchmark plots
./run_1mo_benchmark.py > plots_1mo.log

# Uncomment this line to make 1-year benchmark plots
#./run_1yr_benchmark.py > plots_1yr.log

# Turn off python environment
source deactivate

exit 0

