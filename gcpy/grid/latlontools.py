import numpy as np
import xarray as xr

def latlon_extract(nc_file):
    # Attempt to extract lat and lon data from a netCDF4 dataset
    lon = nc_file['lon'][:]
    lat = nc_file['lat'][:]
    lon_b = latlon_est_bnds(lon)
    lat_b = latlon_est_bnds(lat,force_poles=True)
    return lon_b, lat_b, lon, lat

def latlon_gridarea(lon_b, lat_b, r_earth=6.375e6):
    # Calculate grid areas (m2) for a rectilinear grid
    lon_abs = []
    lastlon = lon_b[0]
    for i,lon in enumerate(lon_b):
        while lon < lastlon:
            lon += 360.0
        lon_abs.append(lon)
        lastlon = lon
    
    n_lat = lat_b.size - 1
    n_lon = lon_b.size - 1
    # Total surface area in each meridional band (allows for a regional domain)
    merid_area = 2*np.pi*r_earth*r_earth*(lon_abs[-1]-lon_abs[0])/(360.0*n_lon)
    grid_area = np.empty([n_lon,n_lat])
    lat_b_rad = np.pi * lat_b / 180.0
    for i_lat in range(n_lat):
        # Fraction of meridional area which applies
        sin_diff = np.sin(lat_b_rad[i_lat+1])-np.sin(lat_b_rad[i_lat])
        grid_area[:,i_lat] = sin_diff * merid_area

    # Transpose this - convention is [lat, lon]
    grid_area = np.transpose(grid_area)
    return grid_area

def latlon_est_bnds(indata,force_poles=False):
    # Estimate lat/lon edges based on a vector of mid-points
    dx = np.median(np.diff(indata))
    x0 = indata.data[0] - (dx/2.0)
    outdata = np.array([x0 + i*dx for i in range(0,indata.size + 1)])
    if force_poles:
        outdata[outdata<-90] = -90.0
        outdata[outdata>90] = 90.0
    return outdata

def latlon_est_mid(indata):
    # Calculate midpoints from edges
    return np.array([0.5*(indata[i] + indata[i+1]) for i in range(len(indata)-1)])
    #return (indata[1:] + indata[:-1]) / 2.0

def make_llvec(lbnd,dl):
    n_mid = np.int(np.round((lbnd[1]-lbnd[0])/dl))
    ledge = np.linspace(start=lbnd[0],stop=lbnd[1],num=n_mid+1)
    lmid  = latlon_est_mid(ledge)
    return lmid,ledge

def gen_hrz_grid(lon_stride,lat_stride,half_polar=False,center_180=False,lon_range=None,lat_range=None):
    # Define a simple rectilinear grid
    
    # Generate longitude edge vector
    n_lon = np.int(np.round(360.0 / lon_stride))
    if center_180:
        start_lon = (-180.0) - (lon_stride/2.0)
    else:
        start_lon = -180.0
    lon_b = np.linspace(start_lon,start_lon+360.0,n_lon + 1)
    
    # Generate latitude edge vector
    # If half-polar, first and last cell are centered on the poles
    n_lat = np.int(np.round(180.0/lat_stride))
    if half_polar:
        n_lat += 1
        lat_max = 90.0 + (lat_stride/2.0)
    else:
        lat_max = 90.0
    lat_b = np.linspace(-lat_max,lat_max,n_lat + 1)
    
    # This will deal with the half-polar issue
    lat_b[0]  = -90.0
    lat_b[-1] =  90.0
    
    lat = latlon_est_mid(lat_b)
    lon = latlon_est_mid(lon_b)
    
    if lon_range is not None:
        lon = lon[lon>=(lon_range[0] - lon_stride*0.01)]
        lon = lon[lon<=(lon_range[-1] + lon_stride*0.01)]
        lon_b = latlon_est_bnds(lon)
        
        lat = lat[lat>=(lat_range[0] - lat_stride*0.01)]
        lat = lat[lat<=(lat_range[-1] + lat_stride*0.01)]
        lat_b = latlon_est_bnds(lat)
    
    return xr.Dataset({'lat': (['lat'], lat),'lon': (['lon'], lon),
                     'lat_b': (['lat_b'], lat_b),'lon_b': (['lon_b'], lon_b)})
