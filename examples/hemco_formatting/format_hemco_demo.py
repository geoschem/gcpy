import xarray as xr
from copy import deepcopy as dc

# ----------------------------------------------------------------- #
# Preparing the file for the demonstration
# ----------------------------------------------------------------- #

# Load and print the demo data. This is a restart file, which we will
# modify to violate HEMCO standards for the sake of the demo.
# TO DO: Fix to absolute path
data = xr.open_dataset('./examples/hemco_formatting/GEOSChem.Restart.20200101_0000z.nc4')

# First, we will drop all but two of the restart, chem, and met
# variables to make this a manageable example.
species_labels = [label for label in data.data_vars.keys()
                  if label[:10] == 'SpeciesRst' or
                  label[:4] == 'Chem' or
                  label[:3] == 'Met']
species_labels = species_labels[:-1]
data = data.drop(species_labels)

# If we were to run ncdump on this file, we would get the output saved
# in demo_original.txt. We will compare this file to equivalent output
# for our final file.
data_mod = dc(data)

# We will now intentionally change the file to be HEMCO incompatible.

# TO DO: Levels question: do we need the formula terms to be always a, b,
# p0, and ps? Or is a, b, and ps adequate, as in the example shown
# on the readdocs

# First, remove one of the attributes from latitude and longitude.
# These changes should all be handled with no manual edits from
# the user.
data_mod['lat'].attrs.pop('units')
data_mod['lon'].attrs.pop('axis')

# We will also reverse latitude so that it's monotonically decreasing.
# This should throw an error.
data_mod['lat'] = data_mod['lat'][::-1]

# Second, remove the time variable. Often, files without an explicit
# time dimension will exclude time from the netcdf. This is bad for
# HEMCO, and we want to make sure that the functions can deal with it.
data_mod = data_mod.drop('time').squeeze()

# Third, mess with the level attributes. We'll add an extra formula
# term that doesn't exist in the dataset. This change should throw an
# error.
data_mod['lev'].attrs['formula_terms'] += ' xs: xx'

# We also change the positive direction. So long as gchp=False is
# passed to the function, this should be handled by the functions.
data_mod['lev'].attrs['positive'] = 'down'

# Finally, we'll change some things in the variable SpeciesRst_ACET.
# We'll add a fourth attribute, which we hope won't be clobbered.
# This should be the only difference between demo_original.txt 
# and the updated demo_post_formatting.txt.
data_mod['SpeciesRst_ACET'].attrs['test'] = (
    'Testing that additional attributes are not clobbered')

# We also delete the units on data_mod SpeciesRst_ACET
del(data_mod['SpeciesRst_ACET'].attrs['units'])

# ----------------------------------------------------------------- #
# Using format_hemco_data to save a HEMCO-compatible file
# ----------------------------------------------------------------- #

# Using format_hemco_data.py is easy and requires only four steps.
data_fix = dc(data_mod)

# 1. Import the module.
from gcpy import format_hemco_data as hemco

# 2. Format the required dimensions (time, lat, lon, and lev) for 
# HEMCO.
# We have to provide the file start time because there is no time
# dimension in this file. If there was, we could still provide a 
# start time, but it would be overwritten (with a warning) with
# the first time value in the dataset.
def test_format_hemco_dimensions(data):
    try:
        data = hemco.format_hemco_dimensions(
            data, 
            start_time='2020-01-01 00:00:00'
        )
    except Exception as error:
        print(f"format_hemco_dimensions_failed: {error}")
    return data

# Let's test this!
data_fix = test_format_hemco_dimensions(data_fix)
print('-'*70)

# We return  an error that "lat is not monotonically increasing."
# Good! We changed that intentionally. Let's undo that and 
# try again.
data_fix['lat'] = data_fix['lat'][::-1]
data_fix = test_format_hemco_dimensions(data_fix)
print('-'*70)

# We also get a warning message that it is assigning the time coordinate 
# from the provided start_time. This is needed for HEMCO compliance, but
# the user should be aware of the specification of the time dimension.

# We find that 'PS' and 'xx' are included in lev_formula_terms but not in
# data_fix. This is a warning, so we don't need to do anything. Onto the 
# next step!

# 3. Format any variables in the netcdf
# Run the checking function.
def test_check_variables(data):
    try:
        hemco.check_hemco_variables(data_fix)
    except Exception as error:
        print(f"check_hemco_variables failed: {error}")

test_check_variables(data_fix)
print('-'*70)

# We get the following error:
# Checking dataset variables for HEMCO compliance.
#   SpeciesRst_ACET missing ['units']
# check_hemco_variables failed: Required units missing from dataset variables.

# We add units back in using the convenience function from the package so
# that we avoid clobbering anything important.
data_fix = hemco.format_hemco_variable(
    data_fix, 
    "SpeciesRst_ACET", 
    units=data["SpeciesRst_ACET"].attrs["units"]
)

# Test one more time
test_check_variables(data_fix)
print('-'*70)

# 4. Save out.
hemco.save_hemco_netcdf(data_fix, './examples/hemco_formatting/', 
                        'HEMCO_demonstration_file_post_fixes.nc')
