from . import latlontools

# Define all the major GEOS-Chem grids
# Horizontal global grids first
gmao_4x5_global       = latlontools.gen_hrz_grid(lat_stride=4.0, lon_stride=5.0, half_polar=True,center_180=True)
gmao_2x25_global      = latlontools.gen_hrz_grid(lat_stride=2.0, lon_stride=2.5, half_polar=True,center_180=True)
gmao_05x0666_global   = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=2/3, half_polar=True,center_180=True)
gmao_05x0625_global   = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=5/8, half_polar=True,center_180=True)
gmao_025x03125_global = latlontools.gen_hrz_grid(lat_stride=0.25,lon_stride=5/16,half_polar=True,center_180=True)

# Horizontal nested grids
gmao_05x0666_us       = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=2/3, half_polar=True,center_180=True,lon_range=[-140, -40],lat_range=[ 10, 70])
gmao_05x0666_ch       = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=2/3, half_polar=True,center_180=True,lon_range=[  70, 150],lat_range=[-11, 55])

gmao_05x0625_as       = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=5/8, half_polar=True,center_180=True,lon_range=[  60, 150],lat_range=[-11, 55])
gmao_05x0625_eu       = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=5/8, half_polar=True,center_180=True,lon_range=[ -30,  50],lat_range=[ 30, 70])
gmao_05x0625_us       = latlontools.gen_hrz_grid(lat_stride=0.5, lon_stride=5/8, half_polar=True,center_180=True,lon_range=[-140, -40],lat_range=[ 10, 70])

# All grids
global_grid_inventory = [gmao_4x5_global,
                         gmao_2x25_global,
                         gmao_05x0666_global,
                         gmao_05x0625_global,
                         gmao_025x03125_global]

nested_grid_inventory = [gmao_05x0666_us,
                         gmao_05x0666_ch,
                         gmao_05x0625_as,
                         gmao_05x0625_eu,
                         gmao_05x0625_us]

def get_grid(grid_shape,is_nested=None,first_call=True):
    # Try to match a grid based only on its size
    # Target not yet found
    is_target = False
    out_grid = None
    # First try global
    if is_nested is None or (not is_nested):
        for grid in global_grid_inventory:
            is_target = grid.lon.size == grid_shape[0] and grid.lat.size == grid_shape[1]
            if is_target:
                out_grid = grid
                break
    if not is_target and (is_nested is None or is_nested):
        for grid in nested_grid_inventory:
            is_target = grid.lon.size == grid_shape[0] and grid.lat.size == grid_shape[1]
            if is_target:
                out_grid = grid
                break
    if not is_target and first_call:
        # Try transposing but prevent recursion
        out_grid = get_grid([grid_shape[1],grid_shape[0]],is_nested,False)
    
    # Return result
    return out_grid
