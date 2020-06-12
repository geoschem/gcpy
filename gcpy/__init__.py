'''
GCPY initialization script.  Imports nested packages for convenience.
'''

from .util import *
from .grid import *
from .regrid import *
from .benchmark import *
from .units import *
from . import plot
from . import grid
from .constants import *
# Figure this out later
#try:
#    from .version import __version__
#except ImportError:
#    raise ImportError('gcpy was not properly installed; some functionality '
#                      'may be not work. If installing from source code, '
#                      'please re-install in place by running\n'
#                      '$ pip install -e .'
#                      '\nElse, please reinstall using your package manager.')
