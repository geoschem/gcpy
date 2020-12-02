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

# Load environment file
. ~/.bashrc

# Turn on Python environment (edit for your setup)
conda activate gcpy

# Uncomment this line to make 1-month benchmark plots & tables
./run_1mo_benchmark.py > plots_1mo.log

# Uncomment this line to make 1-year benchmark plots & tables
#./run_1yr_fullchem_benchmark.py > plots_1yr_fullchem.log

# Uncomment this line to make 1-year TransportTracers plots & tables
#./run_1yr_tt_benchmark.py > plots_1yr_tt.log

# Turn off python environment
conda deactivate

exit 0

