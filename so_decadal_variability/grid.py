import xarray as xr
import numpy as np
from xgcm import Grid
from xgcm.autogenerate import generate_grid_ds
from so_decadal_variability.process import _get_specifics

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
    
def _calc_outerdepths(depth):
    '''Given the cell-centre depth points, determine the outer points.
    Assumption that the first outer point is at the ocean surface, z=0.'''
    
    nk = len(depth)
    depth_i_vals = np.zeros(nk+1)
    for k in range(nk):
        if k>0:
            depth_i_vals[k] = (depth[k-1]+depth[k])/2
    depth_i_vals[nk] = depth[nk-1]+(depth[nk-1]-depth_i_vals[nk-1])
    return depth_i_vals

def get_xgcm(ds,gridlon,gridlat):

    ds = generate_grid_ds(ds, {'X':gridlon,'Y':gridlat})
    xgrid = Grid(ds, periodic=['X'])
    
    # Get horizontal distances
    dlonG = xgrid.diff(ds[gridlon], 'X', boundary_discontinuity=360)
    dlonC = xgrid.diff(ds[gridlon+'_left'], 'X', boundary_discontinuity=360)

    dlatG = xgrid.diff(ds[gridlat], 'Y', boundary='fill', fill_value=np.nan)
    dlatC = xgrid.diff(ds[gridlat+'_left'], 'Y', boundary='fill', fill_value=np.nan)
    
    ds['dxG'], ds['dyG'] = _degrees_to_meters(dlonG, dlatG, ds[gridlon], ds[gridlat])
    ds['dxC'], ds['dyC'] = _degrees_to_meters(dlonC, dlatC, ds[gridlon], ds[gridlat])
    
    # Find depths of outer points
    if 'depth_bnds' in ds.data_vars:
        depth_i_vals = np.append(ds['depth_bnds'].isel(bnds=0), ds['depth_bnds'].isel(bnds=1)[-1])
    else:
        depth_i_vals = _calc_outerdepths(ds['depth'])
    # Put into a dataarray
    ds['depth_i'] = xr.DataArray(depth_i_vals,coords={'depth_i':depth_i_vals},dims={'depth_i'})
    
    # Get depth distance
    ds['dz'] = ds['depth_i'].diff('depth_i').rename({'depth_i':'depth'}).assign_coords({'depth':ds['depth']})
    
    # Regenerate grid
    coords = {
        'X':{'center':gridlon,'left':gridlon+'_left'},
        'Y':{'center':gridlat,'left':gridlat+'_left'},
        'Z':{'center':'depth','outer':'depth_i'}
    }
    metrics = {
        'X':['dxC','dxG'],
        'Y':['dyC','dyG'],
        'Z':['dz']
    }
    xgrid = Grid(ds,coords=coords,metrics=metrics,periodic=['X'])
    
    return ds,xgrid