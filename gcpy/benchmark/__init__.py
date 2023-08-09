"""
GCPY initialization script: gcpy/src/gcpy/benchmark
"""

from .modules.benchmark_models_vs_obs import *
from .modules.run_1yr_fullchem_benchmark import run_benchmark as run_1yr_benchmark
from .modules.run_1yr_tt_benchmark import run_benchmark as run_1yr_tt_benchmark
from .run_benchmark import *
