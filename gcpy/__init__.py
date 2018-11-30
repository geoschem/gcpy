from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core import open_dataset, open_mfdataset, get_collection_data
from .core import convert_bpch_names_to_netcdf_names, add_species_to_dataset
from .core import check_paths, compare_varnames
from .benchmark import compare_single_level, compare_zonal_mean, add_bookmarks_to_pdf
from . import plot


#try:
#    from . version import __version__
#except ImportError:
#    raise ImportError('gcpy was not properly installed; some functionality '
#                      'may be not work. If installing from source code, '
#                      'please re-install in place by running\n'
#                      '$ pip install -e .'
#                      '\nElse, please reinstall using your package manager.')
