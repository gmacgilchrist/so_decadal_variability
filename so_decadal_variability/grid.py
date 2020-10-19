import xarray as xr
import numpy as np
from xgcm import Grid
from xgcm.autogenerate import generate_grid_ds

# Creating an xgcm grid
def _degrees_to_meters(dlon, dlat, lon, lat):
        """Converts lat/lon differentials into distances in meters

        PARAMETERS
        ----------
        dlon : xarray.DataArray longitude differentials
        dlat : xarray.DataArray latitude differentials
        lon  : xarray.DataArray longitude values
        lat  : xarray.DataArray latitude values

        RETURNS
        -------
        dx  : xarray.DataArray distance inferred from dlon
        dy  : xarray.DataArray distance inferred from dlat
        """

        distance_1deg_equator = 111000.0
        dx = dlon * xr.ufuncs.cos(xr.ufuncs.deg2rad(lat)) * distance_1deg_equator
        dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
        return dx, dy

def get_xgcm(ds):
    
    ds = generate_grid_ds(ds, {'X':'lon','Y':'lat'})
    xgrid = Grid(ds, periodic=['X'])
    
    # Get horizontal distances
    dlonG = xgrid.diff(ds['lon'], 'X', boundary_discontinuity=360)
    dlonC = xgrid.diff(ds['lon_left'], 'X', boundary_discontinuity=360)

    dlatG = xgrid.diff(ds['lat'], 'Y', boundary='fill', fill_value=np.nan)
    dlatC = xgrid.diff(ds['lat_left'], 'Y', boundary='fill', fill_value=np.nan)
    
    ds['dxG'], ds['dyG'] = _degrees_to_meters(dlonG, dlatG, ds['lon'], ds['lat'])
    ds['dxC'], ds['dyC'] = _degrees_to_meters(dlonC, dlatC, ds['lon'], ds['lat'])
    
    # Find depths of outer points
    depth_i_vals = np.append(ds['depth_bnds'].isel(bnds=0), ds['depth_bnds'].isel(bnds=1)[-1])
    ds['depth_i'] = xr.DataArray(depth_i_vals,coords={'depth_i':depth_i_vals},dims={'depth_i'})
    # Get depth distance
    ds['dz'] = ds['depth_i'].diff('depth_i').rename({'depth_i':'depth'}).assign_coords({'depth':ds['depth']})
    
    # Regerate grid
    coords = {
        'X':{'center':'lon','left':'lon_left'},
        'Y':{'center':'lat','left':'lat_left'},
        'Z':{'center':'depth','outer':'depth_i'}
    }
    metrics = {
        'X':['dxC','dxG'],
        'Y':['dyC','dyG'],
        'Z':['dz']
    }
    xgrid = Grid(ds,coords=coords,metrics=metrics,periodic=['X'])
    
    return ds,xgrid