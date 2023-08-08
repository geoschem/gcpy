"""
GCPY initialization script: gcpy/src/gcpy
Imports nested packages for convenience.
"""

# Figure this out later
try:
    from ._version import __version__
except ImportError as exc:
    MSG = "Could not find the GCPy version number in _version.py! "
    MSG += "Your GCPy installation may be corrupted.  Please install again."
    raise ImportError(MSG) from exc

from . import benchmark
from . import examples

from .append_grid_corners import *
from .benchmark_funcs import *
from .budget_ox import *
from .budget_tt import *
from .constants import *
from .cstools import *
from .date_time import *
from .grid import *
from .file_regrid import *
from .grid_stretching_transforms import *
from .mean_oh_from_logs import *
from .oh_metrics import *
from .plot import *
from .raveller_1D import *
from .regrid import *
from .ste_flux import *
from .units import *
from .util import *
