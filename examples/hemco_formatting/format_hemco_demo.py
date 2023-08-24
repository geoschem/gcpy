import xarray as xr
import format_hemco_data as hemco

# Load and print the demo data. This is a restart file, which we will
# modify to violate HEMCO standards for the sake of the demo.
data = xr.open_dataset('GEOSChem.Restart.20200101_0000z.nc4')
print(data)