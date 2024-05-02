#!/usr/bin/env python3
"""
Example script using gcpy.community.format_hemco_data.py

NOTE: Before starting this demo, please download the file:

https://gcgrid.s3.amazonaws.com/HEMCO/GCClassic_Output/14.0.0/2019/GEOSChem.ProdLoss.20190101_0000z.nc4

to this folder and rename it to HEMCO_demonstration_file.nc.
"""
import xarray as xr
from copy import deepcopy as dc

# ----------------------------------------------------------------- #
# Preparing the file for the demonstration
# ----------------------------------------------------------------- #

# Load the data file

# NOTE: You can copy any data from the HEMCO data path to this folder:

data = xr.open_dataset("./HEMCO_demonstration_file.nc")


# We will now intentionally change the file to be HEMCO incompatible.

# First, remove one of the attributes from latitude and longitude.
# These changes should all be handled with no manual edits from
# the user.
data["lat"].attrs.pop("units")
data["lon"].attrs.pop("axis")

# We will also reverse latitude so that it"s monotonically decreasing.
# This should throw an error.
data["lat"] = data["lat"][::-1]

# Second, remove the time variable. Often, files without an explicit
# time dimension will exclude time from the netcdf. This is bad for
# HEMCO, and we want to make sure that the functions can deal with it.
data = data.drop("time").squeeze()

# Third, mess with the level attributes. We"ll add an extra formula
# term that doesn"t exist in the dataset. This change should throw an
# error.
data["lev"].attrs["formula_terms"] += " xs: xx"

# We also change the positive direction. So long as gchp=False is
# passed to the function, this should be handled by the functions.
data["lev"].attrs["positive"] = "down"

# Finally, we"ll change some things in the variable SpeciesRst_ACET.
# We"ll add a fourth attribute, which we hope won"t be clobbered.
# This should be the only difference between demo_original.txt
# and the updated demo_post_formatting.txt.
data["Loss_Ox"].attrs["test"] = (
    "Testing that additional attributes are not clobbered"
)

# Save long name and units strings so we can restore it later
save_long_name = dc(data["Loss_Ox"].attrs["long_name"])
save_units = dc(data["Loss_Ox"].attrs["units"])

# We also delete the units on data SpeciesRst_ACET
del(data["Loss_Ox"].attrs["units"])

# ----------------------------------------------------------------- #
# Using format_hemco_data to save a HEMCO-compatible file
# ----------------------------------------------------------------- #
# Using format_hemco_data.py is easy and requires only four steps.
data_fix = dc(data)

# 1. Import the module.
from gcpy.community import format_hemco_data as hemco

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
            start_time="2019-01-01 00:00:00"
        )
    except Exception as error:
        print(f"format_hemco_dimensions_failed: {error}")
    return data

# Let"s test this!
data_fix = test_format_hemco_dimensions(data_fix)
print("-"*70)

# We return  an error that "lat is not monotonically increasing."
# Good! We changed that intentionally. Let"s undo that and
# try again.
data_fix["lat"] = data_fix["lat"][::-1]
data_fix = test_format_hemco_dimensions(data_fix)
print("-"*70)

# We also get a warning message that it is assigning the time coordinate
# from the provided start_time. This is needed for HEMCO compliance, but
# the user should be aware of the specification of the time dimension.

# We find that "PS" and "xx" are included in lev_formula_terms but not in
# data_fix. This is a warning, so we don"t need to do anything. Onto the
# next step!

# 3. Format any variables in the netcdf
# Run the checking function.
def test_check_variables(data):
    try:
        hemco.check_hemco_variables(data_fix)
    except Exception as error:
        print(f"check_hemco_variables failed: {error}")

test_check_variables(data_fix)
print("-"*70)

# We get the following error:
# Checking dataset variables for HEMCO compliance.
#   Loss_Ox missing ["units"]
# check_hemco_variables failed: Required units missing from dataset variables.

# We add units back in using the convenience function from the package so
# that we avoid clobbering anything important.
data_fix = hemco.format_hemco_variable(
    data_fix,
    "Loss_Ox",
    long_name=save_long_name,
    units=save_units,
)

# Test one more time
test_check_variables(data_fix)
print("-"*70)

# 4. Save out.
hemco.save_hemco_netcdf(
    data_fix,
    save_dir=".",
    save_name="./HEMCO_demonstration_file_post_fixes.nc"
)
