import os
from gcpy import benchmark

# Benchmark information (*MUST EDIT*)
maindir  = '/path/to/main/directory'
refstr   = 'ref_version_string'
devstr   = 'dev_version_string'
datestr  = '20160701'

# Directories (edit as needed)
refdir   = os.path.join(maindir, refstr) 
devdir   = os.path.join(maindir, devstr)
plotsdir = os.path.join(maindir, devstr, 'output')

# Files (edit as needed)
spcfile  = 'GEOSChem.SpeciesConc.' + datestr + '_0000z.nc4'
hcofile  = 'HEMCO_diagnostics.'    + datestr + '0000.nc'
aerfile  = 'GEOSChem.Aerosols.'    + datestr + '_0000z.nc4'
jvfile   = 'GEOSChem.JValues.'     + datestr + '_0000z.nc4'

# Paths to data
refspc = os.path.join(maindir, refdir, spcfile)
devspc = os.path.join(maindir, devdir, spcfile)
refhco = os.path.join(maindir, refdir, hcofile)
devhco = os.path.join(maindir, devdir, hcofile)

# Concentration plots (includes lumped species and separates by category)
benchmark.make_gcc_1mo_benchmark_conc_plots(refspc, refstr, devspc, devstr, dst=plotsdir, overwrite=True)

# Emissions plots (one file for now)
benchmark.make_gcc_1mo_benchmark_emis_plots(refhco, refstr, devhco, devstr, dst=plotsdir, overwrite=True)

# Emissions tables
benchmark.make_gcc_1mo_benchmark_emis_tables(refhco, refstr, devhco, devstr, dst=plotsdir, overwrite=True)
