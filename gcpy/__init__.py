'''
GCPY initialization script.  Imports nested packages for convenience.
'''

# Figure this out later
try:
    from ._version import __version__
except ImportError:
    raise ImportError('gcpy was not properly installed; some functionality '
                      'may be not work. If installing from source code, '
                      'please re-install in place by running\n'
                      '$ pip install -e .'
                      '\nElse, please reinstall using your package manager.')

from .util import *
from .date_time import *
from .units import *
from .ste_flux import *
from .regrid import *
from .plot import *
from .oh_metrics import *
from .mean_oh_from_logs import *
from .grid import *
from .constants import *
from .budget_tt import *
from .benchmark import *
from .file_regrid import *
from .grid_stretching_transforms import *
