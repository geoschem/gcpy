#!/usr/bin/env python
"""
Regrids a 4x5 GEOS-Chem Classic restart file to cubed-sphere resolutions.
"""

# Imports
from os.path import join
import numpy as np
import xarray as xr
import sparselt.esmf
import sparselt.xr

def main():

    # Path to regridding weights (EDIT AS NEEDED)
    weights_dir="/path/to/regridding/weights/"

    # List of simulation types (EDIT AS NEEDED)
    simulation_list = ["carboncycle"]

    # List of months (EDIT AS NEEDED)
    month_list = ["01", "07"]

    # List of cubed-sphere grids (EDIT AS NEEDED)
    cubed_sphere_grid_list = ["c24", "c48", "c90", "c180", "c360"]

    # Preserves all global and variable attributes
    with xr.set_options(keep_attrs=True):

        # Loop over simulation types
        for sim in simulation_list:

            # Loop over months
            for mm in month_list:

                # Read input data
                infile = f"GEOSChem.Restart.{sim}.2019{mm}01_0000z.nc4"
                print(f"Reading {infile}")
                ds_in = xr.open_dataset(infile)

                # Rename GCClassic "SpeciesRst_" prefix to GCHP "SPC_" prefix
                old_to_new_names = {}
                for v in ds_in.data_vars.keys():
                    if "SpeciesRst_" in v:
                        new_name = v.replace("SpeciesRst_", "SPC_")
                        old_to_new_names[v] = new_name
                ds_in = ds_in.rename(old_to_new_names)

                # Loop over cubed-sphere grids
                for cs in cubed_sphere_grid_list:

                    # Number of grid points per side
                    cs_res = int(cs[1:])

                    # Regridding transform file
                    regrid_file = f"regrid_weights_latlon46x72_to_{cs}.nc"
                    weights_file = join(weights_dir, regrid_file)

                    # Create a linear transform object from the regridding
                    # weights file for the combination of source and target
                    # horizontal resolutions.  NOTE: GCHP restart files use
                    # a grid where lat = 6*cs_res.
                    transform = sparselt.esmf.load_weights(
                        weights_file,
                        input_dims=[('lat', 'lon'), (46, 72)],
                        output_dims=[('lat', 'lon'), (6*cs_res, cs_res)]
                    )

                    # Regrid to cubed-sphere
                    ds_out = sparselt.xr.apply(transform, ds_in)

                    # Redefine coordinate arrays to be consistent
                    # with GCHP restart file expectations
                    coords_dict = {
                        "lon": np.arange(1, cs_res+1, dtype=np.float64),
                        "lat": np.arange(1, 6*cs_res+1, dtype=np.float64),
                        "lev": np.arange(1, 73, dtype=np.float64),
                    }
                    ds_out = ds_out.assign_coords(coords_dict)

                    # Write to output resolution
                    outfile = f"GEOSChem.Restart.{sim}.2015{mm}01_0000z.{cs}.nc4"
                    print(f"Writing {outfile}")
                    ds_out.to_netcdf(outfile)

                    # Cleanup
                    del transform
                    del ds_out

                # Cleanup
                del ds_in


# Only execute when we run as a standalone script
if __name__ == '__main__':
    main()
