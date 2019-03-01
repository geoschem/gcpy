import os
from gcpy import benchmark


# Directories (edit)
maindir  = '/path/to/main/directory'
refdir = os.path.join(maindir, 'ref_directory_name') 
devdir = os.path.join(maindir, 'dev_directory_name')
plotsdir =os.path.join(maindir,'plots_directory_name')

# Files (edit)
spcfile = 'SpeciesConc_diagnostics_filename'
hcofile = 'HEMCO_diagnostics_filename'

# Paths to data
refspc = os.path.join(maindir, refdir, spcfile)
devspc = os.path.join(maindir, devdir, spcfile)
refhco = os.path.join(maindir, refdir, hcofile)
devhco = os.path.join(maindir, devdir, hcofile)

# Plot strings (edit)
refstr = 'ref_version_string'
devstr = 'dev_version_string'

# Concentration plots (includes lumped species and separates by category)
benchmark.make_gcc_1mo_benchmark_conc_plots(refspc, refstr, devspc, devstr, dst=plotsdir, overwrite=True)

# Emissions plots (one file for now)
benchmark.make_gcc_1mo_benchmark_emis_plots(refhco, refstr, devhco, devstr, dst=plotsdir, overwrite=True)
